# ========================================
# Step 1: Protein Data Quality Control and Cohort Definition
# ========================================
# Quality Control Pipeline for Olink Protein Data
# UK Biobank Pharma Proteomics Project (UKB-PPP)
# Olink Explore 1536/3072 platforms measuring 2,923 plasma proteins

library(data.table)
library(dplyr)
library(ukbtools)

# -------- 1.1 Load Raw Protein Data --------

# Load Olink NPX (Normalized Protein eXpression) values
# Data: 53,014 participants across batches 0-7
protein_raw <- fread("path/to/olink_npx_data.csv")
colnames(protein_raw)[1] <- "f.eid"

cat("Initial data:\n")
cat("Participants:", nrow(protein_raw), "\n")
cat("Proteins measured:", ncol(protein_raw) - 1, "\n")

# -------- 1.2 Protein-Level Quality Control --------

# Exclude proteins with missing rate > 50%
missing_rates <- colSums(is.na(protein_raw[, -1])) / nrow(protein_raw)
proteins_exclude <- names(missing_rates[missing_rates > 0.5])

protein_qc <- protein_raw %>%
  select(-all_of(proteins_exclude))

cat("\nProtein QC:\n")
cat("Proteins excluded (missing >50%):", length(proteins_exclude), "\n")
cat("Proteins retained:", ncol(protein_qc) - 1, "\n")

# -------- 1.3 Sample-Level Quality Control --------

# 1.3.1 Load QC filters
withdrawn <- fread("path/to/w75556_withdrawal.csv")
genetic_qc <- fread("path/to/ukb_genetic_qc.txt") %>%
  select(
    f.eid,
    ancestry = f.22006.0.0,          # Genetic ethnic grouping
    sex_aneuploidy = f.22019.0.0,    # Sex chromosome aneuploidy
    genetic_sex = f.22001.0.0,       # Genetic sex
    reported_sex = f.31.0.0,         # Self-reported sex
    het_missing_outlier = f.22027.0.0 # Heterozygosity/missing outlier
  )

# 1.3.2 Apply sample filters
samples_pass <- genetic_qc %>%
  filter(
    !f.eid %in% withdrawn$f.eid,     # Not withdrawn
    !is.na(ancestry),                # Has ancestry info
    sex_aneuploidy == 0,             # No sex chr aneuploidy
    genetic_sex == reported_sex,     # No sex mismatch
    het_missing_outlier == 0         # No het/missing outlier
  )

cat("Samples after QC:", nrow(samples_pass), "\n")

# 1.3.3 Relatedness exclusion (KING coefficient > 0.0884)
kinship_data <- fread("path/to/ukb_rel.dat")

samples_to_remove <- ukb_gen_samples_to_remove(
  ukbdata = kinship_data,
  ukb.with.data = samples_pass$f.eid,
  cutoff = 0.0884
)

samples_final <- samples_pass %>%
  filter(!f.eid %in% samples_to_remove)

cat("Samples after relatedness exclusion:", nrow(samples_final), "\n")

# -------- 1.4 Define Discovery and Validation Cohorts --------

# Load batch information
batch_info <- fread("path/to/olink_batch_number.dat")
protein_samples <- protein_qc %>%
  filter(f.eid %in% samples_final$f.eid) %>%
  left_join(batch_info, by = c("f.eid" = "sample_id"))

# 1.4.1 Discovery cohort: European ancestry + batches 0-6
discovery_cohort <- protein_samples %>%
  filter(
    ancestry == 1,                   # Field 22006: Caucasian
    Batch %in% 0:6                   # Batches 0-6
  ) %>%
  mutate(FID = f.eid, IID = f.eid) %>%
  select(FID, IID, everything(), -f.eid, -ancestry, -Batch)

cat("\nDiscovery cohort (European, batches 0-6):", nrow(discovery_cohort), "\n")

# 1.4.2 Validation cohort: All others with complete data
# Load self-reported ethnicity (Field 21000)
ethnicity <- fread("path/to/ukb_ethnicity.txt") %>%
  select(f.eid, ethnicity = f.21000.0.0)

validation_cohort <- protein_samples %>%
  filter(!f.eid %in% discovery_cohort$FID) %>%
  left_join(ethnicity, by = "f.eid") %>%
  mutate(
    ethnic_group = case_when(
      ethnicity %in% 1001:1003 ~ "White",
      ethnicity %in% 3001:3004 ~ "Asian",
      ethnicity %in% 4001:4003 ~ "Black",
      ethnicity %in% 2001:2004 ~ "Mixed",
      TRUE ~ "Other"
    ),
    FID = f.eid,
    IID = f.eid
  ) %>%
  select(FID, IID, ethnic_group, everything(), -f.eid)

cat("Validation cohort:", nrow(validation_cohort), "\n")
cat("  - White:", sum(validation_cohort$ethnic_group == "White"), "\n")
cat("  - Asian:", sum(validation_cohort$ethnic_group == "Asian"), "\n")
cat("  - Black:", sum(validation_cohort$ethnic_group == "Black"), "\n")
cat("  - Mixed:", sum(validation_cohort$ethnic_group == "Mixed"), "\n")
cat("  - Other:", sum(validation_cohort$ethnic_group == "Other"), "\n")

# -------- 1.5 Save Cohort Data --------

fwrite(discovery_cohort, "discovery_cohort_proteins.csv")
fwrite(discovery_cohort[, c("FID", "IID")], 
       "discovery_samples.txt", sep = "\t", col.names = FALSE)

fwrite(validation_cohort, "validation_cohort_proteins.csv")
fwrite(validation_cohort[, c("FID", "IID")], 
       "validation_samples.txt", sep = "\t", col.names = FALSE)

cat("\n=== Protein QC and Cohort Definition Complete ===\n")