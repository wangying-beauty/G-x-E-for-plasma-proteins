# ========================================
# Step 2: Genotype Quality Control
# ========================================

# Quality Control for UK Biobank Genetic Data
# Tools: PLINK2, LiftOver

library(data.table)
library(dplyr)

# -------- 2.1 Variant-Level Quality Control --------

# Extract variants with INFO > 0.5 from .mfi files
cat("Extracting variants with INFO > 0.5...\n")

info_variants <- data.frame()
for (chr in 1:22) {
  mfi_file <- sprintf("path/to/ukb_mfi_chr%d_v3.txt", chr)
  if (file.exists(mfi_file)) {
    chr_mfi <- fread(mfi_file) %>%
      filter(V8 > 0.5) %>%        # INFO score column
      select(V2)                   # Variant ID
    info_variants <- rbind(info_variants, chr_mfi)
  }
}

fwrite(info_variants, "variants_info0.5.txt", col.names = FALSE)
cat("Variants with INFO > 0.5:", nrow(info_variants), "\n")


# -------- 2.2 Generate PLINK2 QC Commands --------

# PLINK2 parameters:
# - MAF >= 0.05
# - INFO >= 0.5
# - HWE P >= 1e-6
# - Remove duplicate variants

plink2_commands <- character()

for (chr in 1:22) {
  cmd <- sprintf(
    "plink2 \\
    --bgen path/to/ukb_imp_chr%d_v3.bgen ref-first \\
    --sample path/to/ukb_imp_chr%d_v3.sample \\
    --keep samples_passes_qc.txt \\
    --extract variants_info0.5.txt \\
    --maf 0.05 \\
    --hwe 1e-6 \\
    --rm-dup force-first \\
    --make-bed \\
    --out qc_chr%d \\
    --threads 10 \\
    --memory 32000",
    chr, chr, chr
  )
  plink2_commands <- c(plink2_commands, cmd)
}

writeLines(
  c("#!/bin/bash", "", plink2_commands),
  "run_genotype_qc.sh"
)

cat("\nPLINK2 commands written to: run_genotype_qc.sh\n")

# -------- 2.4 Merge Chromosomes --------

merge_list <- data.frame(file = sprintf("qc_chr%d", 2:22))
fwrite(merge_list, "merge_list.txt", col.names = FALSE)

merge_cmd <- "plink2 \\
  --bfile qc_chr1 \\
  --pmerge-list merge_list.txt bfile \\
  --make-bed \\
  --out ukb_qc_merged \\
  --threads 10"

writeLines(c("#!/bin/bash", "", merge_cmd), "merge_chromosomes.sh")

# -------- 2.5 Liftover Instructions --------

liftover_notes <- "
# Liftover to hg38:
# 1. Convert .bim to UCSC BED format
# 2. Run liftOver with hg19ToHg38 chain
# 3. Update .bim coordinates
# See: https://genome.ucsc.edu/cgi-bin/hgLiftOver
"
writeLines(liftover_notes, "liftover_instructions.txt")

cat("\n=== Genotype QC Setup Complete ===\n")
cat("Expected final variants: ~6,134,546\n")