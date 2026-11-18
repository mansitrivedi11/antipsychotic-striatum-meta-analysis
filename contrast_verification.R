# Statistical provenance and contrast verification
# Project: Antipsychotic striatum transcriptomics meta-analysis

pkgs <- c("dplyr", "readr", "tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

library(dplyr)
library(readr)

# 1. Paths & input validation 

primary_file <- "data/processed/MASTER_primary_MGI_v2_harmonized.rds"
report_dir <- "docs/statistical_validation"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(primary_file)) {
  stop("Primary harmonized dataset not found: ", primary_file)
}

primary <- readRDS(primary_file)

cat("\n--- Script 06C: Statistical Provenance ---\n")
cat("Observations:", nrow(primary), "| Datasets:", n_distinct(primary$dataset), "\n\n")

# 2. Required columns

req_cols <- c(
  "dataset", "species", "drug", "contrast_id", "baseline",
  "vehicle_verified", "treatment_duration", "platform",
  "gene_symbol", "harmonized_mouse_gene_id",
  "log2FC", "t_stat", "pvalue", "SE", "SV", "variance_source"
)

missing_cols <- setdiff(req_cols, names(primary))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# Helper to collapse unique metadata values cleanly
collapse_vals <- function(x) {
  paste(sort(unique(as.character(x))), collapse = "; ")
}

# 3. Contrast-level inventory 
# One row per dataset, drug, and contrast

contrast_inventory <- primary %>%
  group_by(dataset, species, drug, contrast_id) %>%
  summarise(
    n_genes = n_distinct(harmonized_mouse_gene_id, na.rm = TRUE),
    n_rows = n(),
    platform = collapse_vals(platform),
    baseline = collapse_vals(baseline),
    vehicle_verified = collapse_vals(vehicle_verified),
    treatment_duration = collapse_vals(treatment_duration),
    variance_source = collapse_vals(variance_source),
    .groups = "drop"
  ) %>%
  arrange(dataset, drug)

cat("Contrast inventory:\n")
print(contrast_inventory, n = Inf, width = Inf)
write_csv(contrast_inventory, file.path(report_dir, "06C_contrast_inventory.csv"))

# 4. Statistical provenance classification 


provenance_inventory <- contrast_inventory %>%
  mutate(
    statistical_provenance = case_when(
      grepl("limma moderated", variance_source, ignore.case = TRUE) ~ "limma_moderated_t_derived",
      grepl("Gemma t-derived", variance_source, ignore.case = TRUE) ~ "Gemma_t_derived_method_unverified",
      TRUE ~ "unknown"
    ),
    ordinary_sampling_SE_verified = FALSE,
    inverse_variance_readiness = "REQUIRES_REVIEW"
  )

cat("\nStatistical provenance:\n")
provenance_inventory %>%
  dplyr::select(dataset, drug, statistical_provenance, ordinary_sampling_SE_verified, inverse_variance_readiness) %>%
  print(n = Inf, width = Inf)

write_csv(provenance_inventory, file.path(report_dir, "06C_statistical_provenance.csv"))

# 5. Effect-size direction audit 
# Baseline labels alone don't guarantee positive log2FC = treated > control

direction_audit <- contrast_inventory %>%
  mutate(
    baseline_status = case_when(
      grepl("reference substance role", baseline, ignore.case = TRUE) ~ "generic_reference_annotation",
      grepl("verify", baseline, ignore.case = TRUE) ~ "vehicle_formulation_unverified",
      grepl("tween", baseline, ignore.case = TRUE) ~ "named_vehicle_annotation",
      TRUE ~ "other_or_unknown"
    ),
    direction_verified = FALSE,
    direction_review_required = TRUE,
    direction_review_note = "Confirm numerator and denominator from original DE contrast."
  )

cat("\nEffect-size direction audit:\n")
direction_audit %>%
  dplyr::select(dataset, drug, baseline, baseline_status, direction_verified) %>%
  print(n = Inf, width = Inf)

write_csv(direction_audit, file.path(report_dir, "06C_effect_direction_audit.csv"))

# 6. Within-study contrast dependence 
# Flags multi-arm studies sharing controls/batches (covariance not estimated here)

dependence_inventory <- primary %>%
  distinct(dataset, drug, contrast_id) %>%
  group_by(dataset) %>%
  summarise(
    n_drugs = n_distinct(drug),
    n_contrasts = n(),
    drugs = paste(sort(unique(drug)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    possible_dependence = n_contrasts > 1,
    shared_controls_verified = FALSE,
    covariance_available = FALSE,
    review_required = possible_dependence
  )

cat("\nWithin-study dependence:\n")
print(dependence_inventory, n = Inf, width = Inf)
write_csv(dependence_inventory, file.path(report_dir, "06C_within_study_dependence.csv"))

# 7. Gene-level independent replication 
# Counts independent GEO datasets rather than treating multiple drug arms as independent

gene_replication <- primary %>%
  dplyr::filter(!is.na(harmonized_mouse_gene_id), nzchar(harmonized_mouse_gene_id)) %>%
  group_by(harmonized_mouse_gene_id) %>%
  summarise(
    n_datasets = n_distinct(dataset),
    n_drugs = n_distinct(drug),
    n_species = n_distinct(species),
    n_observations = n(),
    n_valid_effects = sum(is.finite(log2FC) & is.finite(SV) & SV > 0),
    .groups = "drop"
  ) %>%
  mutate(
    independent_replication = n_datasets >= 2,
    preliminary_primary_eligibility = n_datasets >= 2 & n_valid_effects >= 2
  )

replication_summary <- tibble::tibble(
  metric = c(
    "Total harmonized genes",
    "Genes represented in >=2 datasets",
    "Genes represented in >=3 datasets",
    "Genes represented in both species",
    "Genes preliminarily eligible for primary analysis"
  ),
  n_genes = c(
    nrow(gene_replication),
    sum(gene_replication$n_datasets >= 2),
    sum(gene_replication$n_datasets >= 3),
    sum(gene_replication$n_species >= 2),
    sum(gene_replication$preliminary_primary_eligibility)
  )
)

cat("\nIndependent replication:\n")
print(replication_summary)

write_csv(gene_replication, file.path(report_dir, "06C_gene_replication.csv"))
write_csv(replication_summary, file.path(report_dir, "06C_gene_replication_summary.csv"))

# 8. Cross-study effect-size distributions


effect_distribution <- primary %>%
  group_by(dataset, drug, species) %>%
  summarise(
    n_rows = n(),
    median_log2FC = median(log2FC, na.rm = TRUE),
    median_abs_log2FC = median(abs(log2FC), na.rm = TRUE),
    median_SE = median(SE, na.rm = TRUE),
    median_SV = median(SV, na.rm = TRUE),
    fraction_positive = mean(log2FC > 0, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nEffect-size distributions:\n")
print(effect_distribution, n = Inf, width = Inf)
write_csv(effect_distribution, file.path(report_dir, "06C_effect_size_distributions.csv"))

# 9. Study-level heterogeneity inventory

heterogeneity_inventory <- primary %>%
  distinct(dataset, species, platform, treatment_duration, drug) %>%
  group_by(dataset) %>%
  summarise(
    species = collapse_vals(species),
    platform = collapse_vals(platform),
    treatment_duration = collapse_vals(treatment_duration),
    drugs = collapse_vals(drug),
    .groups = "drop"
  )

cat("\nHeterogeneity inventory:\n")
print(heterogeneity_inventory, n = Inf, width = Inf)
write_csv(heterogeneity_inventory, file.path(report_dir, "06C_heterogeneity_inventory.csv"))

# 10. Manual verification checklist 

verification_checklist <- contrast_inventory %>%
  transmute(
    dataset,
    drug,
    contrast_id,
    baseline_recorded = baseline,
    variance_source_recorded = variance_source,
    original_DE_method = NA_character_,
    t_statistic_type = NA_character_,
    original_SE_available = NA_character_,
    contrast_numerator = NA_character_,
    contrast_denominator = NA_character_,
    control_group_identity = NA_character_,
    shared_control_samples = NA_character_,
    biological_replicates_confirmed = NA_character_,
    source_reference = NA_character_,
    verification_status = "PENDING"
  )

write_csv(verification_checklist, file.path(report_dir, "06C_manual_verification_checklist.csv"))

cat("\nManual verification checklist:\n")
verification_checklist %>%
  dplyr::select(dataset, drug, contrast_id, verification_status) %>%
  print(n = Inf)

# 11. Statistical readiness summary 
# PASS_FROM_06B refers to the prior audit step and is not re-verified here

readiness <- tibble::tibble(
  requirement = c(
    "Numerical SE/SV consistency",
    "Independent dataset replication evaluated",
    "Original DE method verified",
    "Ordinary sampling variance verified",
    "Treatment-control direction verified",
    "Shared-control dependence evaluated",
    "Appropriate dependence-aware model selected"
  ),
  status = c(
    "PASS_FROM_06B",
    "PASS",
    "PENDING",
    "PENDING",
    "PENDING",
    "PENDING",
    "PENDING"
  )
)

cat("\nStatistical readiness:\n")
print(readiness, n = Inf, width = Inf)
write_csv(readiness, file.path(report_dir, "06C_statistical_readiness.csv"))

# 12. Session info & wrap-up 

writeLines(capture.output(sessionInfo()), file.path(report_dir, "06C_session_info.txt"))

cat("\n--- Script 06C Completed ---\n")
cat("Primary observations:", nrow(primary), "\n")
cat("Independent datasets:", n_distinct(primary$dataset), "\n")
cat("Drug contrasts:", nrow(contrast_inventory), "\n")
cat("Reports saved in:", report_dir, "\n")
cat("No master datasets were modified.\n")
