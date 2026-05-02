# ========================================
# Step 3: Residual Calculation and vQTL Analysis (Z-Transform)
# ========================================

# Protein phenotype adjustment and variance QTL discovery
# Tools: R (residualization), OSCA (vQTL)

library(data.table)
library(dplyr)
library(furrr)
library(here)
library(purrr)

# Enable parallel processing
plan(multisession, workers = 6)
options(future.globals.maxSize = 2000 * 1024^2)

# -------- 3.1 Load Data --------

# Load protein phenotypes (discovery cohort)
phenotypes <- fread("discovery_cohort_proteins.csv")

# Load covariates
# Standard covariates: sex, age, PC1-40
std_cov <- fread("path/to/std_covariates.csv") %>%
  mutate(
    agesex = age * sex,
    age2 = age^2,
    age2sex = age^2 * sex
  ) %>%
  filter(FID %in% phenotypes$FID) %>%
  drop_na()

# Specific covariates: assessment center, genotype array, batch
spec_cov <- fread("path/to/spec_covariates.csv") %>%
  mutate(
    f.54.0.0 = as.factor(f.54.0.0),        # Assessment center
    f.22000.0.0 = as.factor(f.22000.0.0),  # Genotype array
    Batch = as.factor(Batch)                # Olink batch
  ) %>%
  filter(FID %in% phenotypes$FID) %>%
  drop_na()

# -------- 3.2 Calculate Protein-Specific Sample Dates --------

# Load Olink metadata
batch <- fread("path/to/olink_batch_number.dat")
plate <- fread("path/to/olink_assay.dat")
processing_date <- fread("path/to/olink_processing_start_date.dat")

# Merge protein measurements with dates
protein_with_date <- phenotypes %>%
  left_join(batch, by = c("PlateID" = "PlateID")) %>%
  left_join(processing_date, by = c("PlateID" = "PlateID"))

# Create date matrix: rows = samples, columns = proteins
# Each cell contains the processing date for that protein-sample combination
date_matrix <- data.frame(FID = unique(protein_with_date$FID))

for (panel_type in unique(plate$Panel)) {
  # Get proteins in this panel
  assays <- plate$Assay[plate$Panel == panel_type]
  
  # Get processing dates for this panel
  dates <- protein_with_date %>%
    filter(Panel == panel_type) %>%
    select(FID, Processing_StartDate)
  
  # Fill dates for all proteins in panel
  for (assay in assays) {
    date_matrix[[paste0(assay, "_date")]] <- 
      dates$Processing_StartDate[match(date_matrix$FID, dates$FID)]
  }
}

# Calculate days from baseline assessment
baseline_date <- fread("path/to/ukb_baseline_dates.txt")  # Field 21842

date_matrix_days <- date_matrix %>%
  left_join(baseline_date, by = "FID") %>%
  mutate(across(
    ends_with("_date"),
    ~as.numeric(difftime(as.Date(.), f.21842.0.0, units = "days"))
  )) %>%
  select(-f.21842.0.0)

# -------- 3.3 Merge All Data --------

# Combine phenotypes, covariates, and dates
merged_data <- phenotypes %>%
  inner_join(std_cov, by = c("FID", "IID")) %>%
  inner_join(spec_cov, by = c("FID", "IID")) %>%
  inner_join(date_matrix_days, by = "FID")

cat("Samples with complete data:", nrow(merged_data), "\n")

# -------- 3.4 Calculate Residuals --------

# Function: regress out covariates and dates
calculate_residuals <- function(df, protein_name) {
  # Define regression formula including protein-specific date
  date_var <- paste0(protein_name, "_date")
  
  formula <- as.formula(paste(
    protein_name, "~ sex + age + agesex + age2 + age2sex +",
    paste0("PC", 1:40, collapse = " + "),
    "+ f.22000.0.0 + f.54.0.0 + Batch +",
    date_var
  ))
  
  # Fit linear model (only complete cases)
  cols_needed <- c("FID", "IID", protein_name, "sex", "age", "agesex", 
                   "age2", "age2sex", paste0("PC", 1:40),
                   "f.22000.0.0", "f.54.0.0", "Batch", date_var)
  
  complete_idx <- complete.cases(df[, cols_needed])
  model <- lm(formula, data = df[complete_idx, ])
  
  # Extract residuals
  result <- data.frame(
    FID = df$FID[complete_idx],
    IID = df$IID[complete_idx],
    Residual = resid(model)
  )
  
  return(result)
}

# Calculate residuals for all proteins (parallel)
protein_names <- setdiff(names(phenotypes), c("FID", "IID"))

residuals_list <- future_map(
  protein_names,
  ~calculate_residuals(merged_data, .x),
  .options = furrr_options(seed = TRUE)
) %>%
  set_names(protein_names)

cat("Residuals calculated for", length(residuals_list), "proteins\n")

# -------- 3.5 Z-Transformation (Standardization) --------

# Apply Z-score standardization to residuals instead of INT
cat("\nApplying Z-transformation (Standardization) to residuals...\n")

residuals_z <- map(residuals_list, function(df) {
  df %>%
    # Using scale() for Z-transformation: (x - mean(x)) / sd(x)
    mutate(Residual_Z = as.numeric(scale(Residual))) %>% 
    select(FID, IID, Residual_Z)
})

# -------- 3.6 Save Phenotype Files for OSCA --------

# OSCA requires: FID IID phenotype
output_dir <- "path/to/osca_phenotypes_Z"
dir.create(output_dir, showWarnings = FALSE)

walk2(residuals_z, names(residuals_z), function(df, prot_name) {
  output_file <- file.path(output_dir, paste0(prot_name, "_Z.txt"))
  fwrite(df, output_file, sep = "\t", col.names = TRUE)
})

cat("Z-transformed phenotype files saved to:", output_dir, "\n")

# -------- 3.7 Generate OSCA vQTL Scripts --------

# vQTL model:
# Var(Y|X) = exp(β₀ + β₁X)
# H₀: β₁ = 0 (no variance heterogeneity)
# Test statistic: Levene's test or Bartlett's test

# OSCA command template
osca_template <- "#!/bin/bash

phenoname=$1
phenofile=$2
chr=$3
OutPath=$4

# OSCA vQTL analysis
/path/to/osca_Linux \\
  --vqtl \\
  --bfile /path/to/qc_chr${chr} \\
  --pheno \"$phenofile\" \\
  --vqtl-mtd 2 \\
  --thread-num 10 \\
  --out \"$OutPath/vQTL_${phenoname}_chr${chr}\"

# Filter significant results (P < 1.71e-11)
awk 'BEGIN{FS=\" \"}{if ($7<1.71e-11){print}}' \\
  \"$OutPath/vQTL_${phenoname}_chr${chr}.vqtl\" \\
  > \"$OutPath/vQTL_${phenoname}_chr${chr}_sig.txt\"

echo \"vQTL analysis completed: $phenoname, chr$chr\"
"

writeLines(osca_template, "osca_vqtl_template.sh")

# Generate scripts for each protein × chromosome combination
script_dir <- "path/to/vqtl_scripts"
dir.create(script_dir, showWarnings = FALSE)

for (prot in names(residuals_z)) {
  for (chr in 1:22) {
    script_content <- sprintf(
      "#!/bin/bash
bash osca_vqtl_template.sh %s %s %d %s",
      prot,
      file.path(output_dir, paste0(prot, "_Z.txt")),
      chr,
      "path/to/vqtl_results"
    )
    
    script_file <- file.path(script_dir, 
                             sprintf("vqtl_%s_chr%d.sh", prot, chr))
    writeLines(script_content, script_file)
    Sys.chmod(script_file, mode = "0755")
  }
}

cat("vQTL scripts generated in:", script_dir, "\n")

# -------- 3.8 Summary --------

cat("\n=== Residual Calculation and vQTL Setup Complete ===\n")
cat("Proteins processed:", length(residuals_z), "\n")
cat("Samples with residuals:", nrow(residuals_z[[1]]), "\n")
cat("\nNext steps:\n")
cat("1. Run: bash osca_vqtl_template.sh [protein] [phenofile] [chr] [outpath]\n")
cat("2. Combine results across chromosomes\n")
cat("3. Apply multiple testing correction\n")

# Close parallel processing
plan(sequential)