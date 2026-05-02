# -------- 4.1 Apply Bonferroni Correction --------

cat("\n=== Applying Multiple Testing Correction ===\n")

# Calculate Bonferroni threshold
# Threshold = 0.05 / (N_SNPs × N_proteins)
bonferroni_threshold <- 1.71e-11  # Adjust based on your analysis

# Filter significant results
significant_vqtls <- all_vqtl_results %>%
  filter(P < bonferroni_threshold) %>%
  mutate(log10P = -log10(P)) %>%
  arrange(desc(log10P))

cat(sprintf("Significant vQTLs at P < %.2e: %d\n", 
            bonferroni_threshold, nrow(significant_vqtls)))

# Save significant results
fwrite(significant_vqtls, 
       file.path(combined_results_dir, "significant_vQTLs_bonferroni.txt"),
       sep = "\t")

# -------- 4.2 LD-based Clumping --------

cat("\n=== Step 5: LD-based Clumping ===\n")

# Add genomic position information if not already present
# Assuming SNP column format: "rs123" or "chr:pos"
significant_vqtls <- significant_vqtls %>%
  mutate(
    Chr = as.character(Chr),
    # Add bp position from SNP name if needed
    # bp = as.numeric(gsub(".*:", "", SNP))
  )

# Define clumping windows
# Standard: ±1 Mb window
# MHC/HLA region (Chr 6: 25.5-34.0 Mb): treat as single region
significant_vqtls <- significant_vqtls %>%
  mutate(
    HLA = ifelse(Chr == "6" & bp > 25.5e6 & bp < 34.0e6, 1, 0),
    position_start = ifelse(HLA == 1, 25.5e6, bp - 1e6),
    position_end = ifelse(HLA == 1, 34.0e6, bp + 1e6),
    position_start = ifelse(position_start < 0, 0, position_start)
  ) %>%
  arrange(desc(log10P))

# Get unique protein-chromosome combinations
protein_chr_pairs <- significant_vqtls %>%
  distinct(Protein, Chr)

cat(sprintf("Clumping %d protein-chromosome pairs...\n", 
            nrow(protein_chr_pairs)))

# Perform clumping for each protein-chromosome pair
clumped_regions <- list()

for (i in seq(nrow(protein_chr_pairs))) {
  
  if (i %% 50 == 0) {
    cat(sprintf("Processing pair %d/%d...\n", i, nrow(protein_chr_pairs)))
  }
  
  # Extract current protein-chromosome data
  current_data <- significant_vqtls %>%
    filter(Protein == protein_chr_pairs$Protein[i] & 
             Chr == protein_chr_pairs$Chr[i]) %>%
    filter(!is.na(position_start) & !is.na(position_end))
  
  if (nrow(current_data) == 0) next
  
  # Use GenomicRanges to merge overlapping windows
  clumped_result <- current_data %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      seqnames.field = "Chr",
      start.field = "position_start",
      end.field = "position_end"
    ) %>%
    GenomicRanges::reduce() %>%
    as_tibble()
  
  if (nrow(clumped_result) == 0) next
  
  # Add region information
  clumped_result <- clumped_result %>%
    select(Chr = seqnames, region_start = start, region_end = end) %>%
    mutate(
      region_size = region_end - region_start,
      region_index = paste0(i, "_", 1:n()),
      Protein = protein_chr_pairs$Protein[i]
    )
  
  clumped_regions[[i]] <- clumped_result
}

# Combine all clumped regions
clumped_regions_merged <- bind_rows(clumped_regions)

cat(sprintf("Total clumped regions identified: %d\n", 
            nrow(clumped_regions_merged)))

# -------- 4.3 Select Lead vQTLs --------

cat("\n=== Selecting Lead vQTLs within Each Region ===\n")

# Assign vQTLs to clumped regions
vqtls_in_regions <- significant_vqtls %>%
  mutate(Chr = as.character(Chr)) %>%
  left_join(clumped_regions_merged %>% mutate(Chr = as.character(Chr)), 
            by = c("Protein", "Chr")) %>%
  mutate(in_region = ifelse(bp >= region_start & bp <= region_end, 1, 0)) %>%
  filter(in_region == 1) %>%
  mutate(distance_to_midreg = abs(bp - (region_start + region_end) / 2))

# Select lead vQTL in each region
# Selection criteria (in order of priority):
# 1. Highest significance (max log10P)
# 2. Largest sample size (max NMISS)
# 3. Closest to region midpoint (min distance_to_midreg)
# 4. Smallest genomic position (min bp)
lead_vqtls <- vqtls_in_regions %>%
  group_by(Protein, Chr, region_index) %>%
  filter(log10P == max(log10P)) %>%
  filter(NMISS == max(NMISS)) %>%
  filter(distance_to_midreg == min(distance_to_midreg)) %>%
  filter(bp == min(bp)) %>%
  ungroup()

cat(sprintf("Lead vQTLs identified: %d\n", nrow(lead_vqtls)))

# Save clumped results
fwrite(lead_vqtls, 
       file.path(combined_results_dir, "lead_vQTLs_clumped.txt"),
       sep = "\t")