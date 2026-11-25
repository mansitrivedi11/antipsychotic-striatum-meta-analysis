# Statistical validation and meta-analysis readiness


library(dplyr)
library(readr)

processed_dir <- "data/processed"
report_dir <- "docs/statistical_validation"

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

primary_file <- file.path(
  processed_dir, "MASTER_primary_MGI_v2_harmonized.rds"
)

qc_sensitivity_file <- file.path(
  processed_dir, "MASTER_QC_sensitivity_MGI_v2_harmonized.rds"
)

identifier_sensitivity_file <- file.path(
  processed_dir, "MASTER_primary_MGI_v2_no_identifier_conflicts.rds"
)

if (!file.exists(primary_file)) {
  stop("Primary harmonized dataset not found: ", primary_file)
}

primary <- readRDS(primary_file)

cat("Primary observations:", nrow(primary), "\n")

required_columns <- c(
  "dataset", "species", "drug", "drug_class", "contrast_id",
  "gene_symbol", "harmonized_mouse_gene_id", "platform",
  "log2FC", "t_stat", "pvalue", "SE", "SV", "variance_source"
)

missing_columns <- setdiff(required_columns, names(primary))

if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

save_report <- function(data, filename) {
  write_csv(data, file.path(report_dir, filename))
}

# Summarise finite values without returning Inf for empty groups.

finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 0) max(x) else NA_real_
}

finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 0) median(x) else NA_real_
}


# Check statistical column types

numeric_columns <- c("log2FC", "t_stat", "pvalue", "SE", "SV")

column_types <- tibble(
  column = numeric_columns,
  class = vapply(
    primary[numeric_columns],
    function(x) paste(class(x), collapse = "/"),
    character(1)
  )
)

print(column_types)

non_numeric <- numeric_columns[
  !vapply(primary[numeric_columns], is.numeric, logical(1))
]

if (length(non_numeric) > 0) {
  stop("Non-numeric columns: ", paste(non_numeric, collapse = ", "))
}


# Numerical validity of effect sizes and variance estimates

numerical_validation <- primary %>%
  group_by(dataset, drug, species) %>%
  summarise(
    n_rows = n(),

    n_missing_log2FC = sum(is.na(log2FC)),
    n_nonfinite_log2FC = sum(!is.na(log2FC) & !is.finite(log2FC)),

    n_missing_t = sum(is.na(t_stat)),
    n_nonfinite_t = sum(!is.na(t_stat) & !is.finite(t_stat)),

    n_missing_SE = sum(is.na(SE)),
    n_nonfinite_SE = sum(!is.na(SE) & !is.finite(SE)),
    n_nonpositive_SE = sum(is.finite(SE) & SE <= 0),

    n_missing_SV = sum(is.na(SV)),
    n_nonfinite_SV = sum(!is.na(SV) & !is.finite(SV)),
    n_nonpositive_SV = sum(is.finite(SV) & SV <= 0),

    n_missing_pvalue = sum(is.na(pvalue)),
    n_invalid_pvalue = sum(
      !is.na(pvalue) &
        (!is.finite(pvalue) | pvalue < 0 | pvalue > 1)
    ),

    .groups = "drop"
  )

cat("\nNumerical validation\n")
print(numerical_validation, n = Inf, width = Inf)

save_report(numerical_validation, "06B_numerical_validation.csv")


# Check whether sampling variance equals squared standard error.
# Arithmetic agreement does not establish statistical validity.

variance_consistency <- primary %>%
  mutate(
    se_squared = SE^2,
    variance_difference = abs(SV - se_squared),
    variance_relative_difference =
      variance_difference / pmax(abs(SV), abs(se_squared), 1e-12)
  ) %>%
  group_by(dataset, drug, variance_source) %>%
  summarise(
    n_rows = n(),

    n_comparable = sum(
      is.finite(SV) & is.finite(se_squared)
    ),

    n_inconsistent = sum(
      is.finite(variance_relative_difference) &
        variance_relative_difference > 1e-6
    ),

    max_relative_difference = finite_max(
      variance_relative_difference
    ),

    .groups = "drop"
  )

cat("\nSE and SV consistency\n")
print(variance_consistency, n = Inf, width = Inf)

save_report(variance_consistency, "06B_SE_SV_consistency.csv")


# Check standard errors implied by the reported t-statistics.
#
# For a coefficient t-statistic, SE = abs(log2FC / t).
# Moderated t-statistics may imply moderated uncertainty rather
# than ordinary sampling uncertainty, so this is a numerical check.

t_validation <- primary %>%
  mutate(
    implied_SE = if_else(
      is.finite(log2FC) & is.finite(t_stat) & t_stat != 0,
      abs(log2FC / t_stat),
      NA_real_
    ),

    implied_SE_relative_difference =
      abs(SE - implied_SE) /
      pmax(abs(SE), abs(implied_SE), 1e-12)
  ) %>%
  group_by(dataset, drug, variance_source) %>%
  summarise(
    n_rows = n(),

    n_testable = sum(
      is.finite(implied_SE) & is.finite(SE)
    ),

    n_matching = sum(
      is.finite(implied_SE_relative_difference) &
        implied_SE_relative_difference <= 1e-4
    ),

    n_mismatching = sum(
      is.finite(implied_SE_relative_difference) &
        implied_SE_relative_difference > 1e-4
    ),

    median_relative_difference = finite_median(
      implied_SE_relative_difference
    ),

    .groups = "drop"
  )

cat("\nT-statistic-derived standard errors\n")
print(t_validation, n = Inf, width = Inf)

save_report(t_validation, "06B_t_stat_SE_validation.csv")


# Review the source of variance estimates across datasets

variance_inventory <- primary %>%
  group_by(dataset, drug, platform, variance_source) %>%
  summarise(
    n_rows = n(),
    median_SE = finite_median(SE),
    median_SV = finite_median(SV),
    .groups = "drop"
  )

cat("\nVariance source inventory\n")
print(variance_inventory, n = Inf, width = Inf)

save_report(
  variance_inventory,
  "06B_variance_source_inventory.csv"
)


# Review treatment contrasts and available baseline metadata.
# These fields still require verification against the original studies.

baseline_columns <- intersect(
  c("baseline", "vehicle_verified", "treatment_duration", "subregion"),
  names(primary)
)

direction_inventory <- primary %>%
  distinct(
    across(
      all_of(c("dataset", "drug", "contrast_id", baseline_columns))
    )
  )

cat("\nEffect direction and baseline metadata\n")
print(direction_inventory, n = Inf, width = Inf)

save_report(
  direction_inventory,
  "06B_effect_direction_inventory.csv"
)


# Identify potentially dependent contrasts.
#
# Contrasts from the same dataset may share control samples,
# animals, experimental batches or model-fitting procedures.
# Dataset membership alone does not establish covariance.

contrast_inventory <- primary %>%
  group_by(dataset, species) %>%
  summarise(
    n_drugs = n_distinct(drug),
    n_contrast_ids = n_distinct(contrast_id),
    n_rows = n(),
    drugs = paste(sort(unique(drug)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    possible_within_study_dependence =
      n_drugs > 1 | n_contrast_ids > 1
  )

cat("\nPotential contrast dependence\n")
print(contrast_inventory, n = Inf, width = Inf)

save_report(
  contrast_inventory,
  "06B_contrast_dependence_inventory.csv"
)


# Gene-level coverage across datasets, drugs and species

gene_coverage <- primary %>%
  group_by(harmonized_mouse_gene_id) %>%
  summarise(
    n_observations = n(),
    n_datasets = n_distinct(dataset),
    n_drugs = n_distinct(drug),
    n_species = n_distinct(species),

    n_contrast_ids = n_distinct(
      paste(dataset, drug, contrast_id, sep = "::")
    ),

    n_valid_variances = sum(is.finite(SV) & SV > 0),

    .groups = "drop"
  )

gene_coverage_summary <- tibble(
  metric = c(
    "All harmonized genes",
    "Genes in at least 2 datasets",
    "Genes in at least 3 datasets",
    "Genes in both species",
    "Genes with at least 2 valid variances"
  ),

  n_genes = c(
    nrow(gene_coverage),
    sum(gene_coverage$n_datasets >= 2),
    sum(gene_coverage$n_datasets >= 3),
    sum(gene_coverage$n_species == 2),
    sum(gene_coverage$n_valid_variances >= 2)
  )
)

cat("\nGene-level coverage\n")
print(gene_coverage_summary, n = Inf)

save_report(gene_coverage, "06B_gene_coverage.csv")

save_report(
  gene_coverage_summary,
  "06B_gene_coverage_summary.csv"
)


# Identify genes contributing multiple effects within a dataset.
# These observations should not automatically be treated as independent.

within_study_gene_counts <- primary %>%
  count(harmonized_mouse_gene_id, dataset, name = "n_effects") %>%
  filter(n_effects > 1)

cat(
  "\nGene-dataset combinations with multiple effects:",
  nrow(within_study_gene_counts), "\n"
)

save_report(
  within_study_gene_counts,
  "06B_within_study_multiple_effects.csv"
)


# Check availability of the sensitivity datasets.
# These files are not analysed or modified here.

sensitivity_inventory <- tibble(
  analysis = c(
    "Primary",
    "GSE66275 QC sensitivity",
    "Identifier conflict sensitivity"
  ),

  file = c(
    primary_file,
    qc_sensitivity_file,
    identifier_sensitivity_file
  )
) %>%
  mutate(exists = file.exists(file))

cat("\nSensitivity dataset availability\n")
print(sensitivity_inventory, n = Inf, width = Inf)

save_report(
  sensitivity_inventory,
  "06B_sensitivity_file_inventory.csv"
)


# Preliminary meta-analysis readiness
#
# Numerical checks and methodological verification are kept separate.
# Passing the numerical checks does not establish meta-analysis readiness.

duplicate_effects <- primary %>%
  count(
    dataset, drug, contrast_id,
    harmonized_mouse_gene_id
  ) %>%
  filter(n > 1)

variance_difference <- abs(primary$SV - primary$SE^2) /
  pmax(abs(primary$SV), abs(primary$SE^2), 1e-12)

valid_log2FC <- all(is.finite(primary$log2FC))

valid_SE <- all(
  is.finite(primary$SE) & primary$SE > 0
)

valid_SV <- all(
  is.finite(primary$SV) & primary$SV > 0
)

valid_variance <- all(
  is.finite(primary$SE) &
    is.finite(primary$SV) &
    variance_difference <= 1e-6
)

unique_effects <- nrow(duplicate_effects) == 0

readiness <- tibble(
  check = c(
    "All log2FC finite",
    "All SE positive and finite",
    "All SV positive and finite",
    "SE and SV numerically consistent",
    "No repeated gene-contrast identifiers",
    "Variance derivation method independently verified",
    "Within-study covariance independently established",
    "Contrast direction independently verified"
  ),

  status = c(
    if (valid_log2FC) "PASS" else "REVIEW",
    if (valid_SE) "PASS" else "REVIEW",
    if (valid_SV) "PASS" else "REVIEW",
    if (valid_variance) "PASS" else "REVIEW",
    if (unique_effects) "PASS" else "REVIEW",
    "NOT VERIFIED",
    "NOT VERIFIED",
    "NOT VERIFIED"
  )
)

cat("\nPreliminary meta-analysis readiness\n")
print(readiness, n = Inf, width = Inf)

save_report(readiness, "06B_preliminary_readiness.csv")


# Record the R environment for reproducibility

writeLines(
  capture.output(sessionInfo()),
  file.path(report_dir, "06B_session_info.txt")
)

cat("\nStatistical validation completed.\n")
cat("Reports saved in:", report_dir, "\n")
cat("No input datasets were modified.\n")

cat(
  "Review variance derivation, effect direction and contrast",
  "dependence before further statistical analysis.\n"
)
