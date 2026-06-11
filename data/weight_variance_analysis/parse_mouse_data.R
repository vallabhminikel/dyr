#!/usr/bin/env Rscript
# Parse and clean mouse and hamster weight data
# This script processes two TSV files with longitudinal weight data
# - SG_SH_weights_individual.tsv: Hamster data
# - 2020_aso_weights_saline.tsv: Mouse data

library(tidyverse)

# Load the data files
cat("Loading data files...\n")
blind <- readxl::read_xlsx("DYR_experiments_blind_match.xlsx") %>% as_tibble %>%
  filter(expt == "RG")
rat_data <- read_tsv("RG_weights.tsv",) %>%
  left_join(blind,by = "animal") %>%
  mutate (
    animal = as.character(animal),
    species = "rat",
    inoculated = if_else(inoculum == "uninoculated",FALSE,TRUE),
    strain = cohort_disp,
    dpi = wtdpi
  )


# File 1: SG_SH_weights_individual.tsv (HAMSTER DATA)
# - Rename wtdpi to dpi for consistency
# - Split cohort_disp into genotype and inoculation status
sg_sh_data <- read_tsv("SG_SH_weights_individual.tsv", show_col_types = FALSE) %>%
  rename(dpi = wtdpi) %>%
  mutate(
    animal = as.character(animal),  # Convert animal ID to character
    species = "hamster"  # Add species identifier
  ) %>%
  separate(cohort_disp,
           into = c("genotype", "inoculation"),
           sep = " (?=\\d|uninoculated)",  # Split before strain code or "uninoculated"
           remove = FALSE) %>%
  mutate(
    # Clean up genotype (remove trailing spaces)
    genotype = str_trim(genotype),
    # Classify inoculation status as binary
    inoculated = if_else(inoculation == "uninoculated", FALSE, TRUE),
    # Keep the strain information
    strain = if_else(inoculation == "uninoculated", NA_character_, inoculation)
  )

cat("SG_SH HAMSTER data processed:\n")
cat(sprintf("  - %d observations\n", nrow(sg_sh_data)))
cat(sprintf("  - Genotypes: %s\n", paste(unique(sg_sh_data$genotype), collapse = ", ")))
cat(sprintf("  - Strains: %s\n", paste(na.omit(unique(sg_sh_data$strain)), collapse = ", ")))

# File 2: 2020_aso_weights_saline.tsv (MOUSE DATA)
# - Already has dpi column
# - Create consistent inoculation status based on treatment
aso_data <- read_tsv("2020_aso_weights_saline.tsv", show_col_types = FALSE) %>%
  mutate(
    # Convert animal ID to character for consistency
    animal = as.character(animal),
    # Add species identifier
    species = "mouse",
    # Classify inoculation based on treatment
    inoculated = if_else(treatment == "saline", FALSE, TRUE),
    # Use strain column directly but set to NA for saline
    strain = if_else(treatment == "saline", NA_character_, strain),
    # Create genotype field (this dataset doesn't have explicit genotype in column)
    # Assuming all are WT based on the data structure, adjust if needed
    genotype = "WT"  # May need adjustment based on your experimental design
  )

cat("\nASO MOUSE data processed:\n")
cat(sprintf("  - %d observations\n", nrow(aso_data)))
cat(sprintf("  - Treatments: %s\n", paste(unique(aso_data$treatment), collapse = ", ")))
cat(sprintf("  - Strains: %s\n", paste(na.omit(unique(aso_data$strain)), collapse = ", ")))


# Verify dpi consistency
cat("\nDPI ranges:\n")
cat(sprintf("  - SG_SH: %d to %d\n", min(sg_sh_data$dpi), max(sg_sh_data$dpi)))
cat(sprintf("  - ASO: %d to %d\n", min(aso_data$dpi), max(aso_data$dpi)))

# Create a combined dataset with common columns
combined_data <- bind_rows(
  sg_sh_data %>%
    select(animal, dpi, weight, species, genotype, inoculated, strain) %>%
    mutate(dataset = "SG_SH"),
  aso_data %>%
    select(animal, dpi, weight, species, genotype, inoculated, strain) %>%
    mutate(dataset = "ASO"),
  rat_data %>%
    select(animal,dpi,weight,species,genotype,inoculated,strain)
)

cat("\nCombined dataset:\n")
cat(sprintf("  - %d total observations\n", nrow(combined_data)))
cat(sprintf("  - %d unique animals\n", n_distinct(combined_data$animal)))
cat(sprintf("  - Species: %s\n", paste(unique(combined_data$species), collapse = ", ")))

# Summary statistics by species and group
summary_stats <- combined_data %>%
  group_by(species, genotype, inoculated, strain) %>%
  summarise(
    n_animals = n_distinct(animal),
    n_observations = n(),
    mean_weight = mean(weight, na.rm = TRUE),
    sd_weight = sd(weight, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nSummary by species and group:\n")
print(summary_stats, n = Inf)

# Save processed data
write_tsv(sg_sh_data, "sg_sh_data_processed.tsv")
write_tsv(aso_data, "aso_data_processed.tsv")
write_tsv(rat_data,"rat_data_processed.tsv")
write_tsv(combined_data, "combined_data_processed.tsv")

cat("\nProcessed data saved to:\n")
cat("  - sg_sh_data_processed.tsv\n")
cat("  - aso_data_processed.tsv\n")
cat("  - combined_data_processed.tsv\n")

