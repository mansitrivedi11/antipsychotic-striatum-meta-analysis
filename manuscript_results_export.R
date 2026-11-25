# Manuscript results export
# Comparative transcriptomic analysis of haloperidol and
# risperidone responses in the rodent striatum.
#
# Uses existing analytical outputs to prepare manuscript text,
# figure legends and results tables. No new statistical tests.

library(readr)
library(dplyr)
library(tibble)

results_dir <- "results_v2"
output_dir <- file.path(results_dir, "13_manuscript_export")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

src_dir <- c(
  s06 = file.path(results_dir, "06_exploratory_analysis", "tables"),
  s07 = file.path(results_dir, "07_QC_sensitivity", "tables"),
  s08 = file.path(results_dir, "08_statistical_convergence", "tables"),
  s09 = file.path(results_dir, "09_biological_interpretation", "tables")
)


# Helpers

read_source <- function(script, filename) {
  path <- file.path(src_dir[[script]], filename)

  if (!file.exists(path))
    stop("Required file missing: ", path)

  dat <- read_csv(path, show_col_types = FALSE, name_repair = "minimal")

  if (nrow(dat) == 0)
    stop("Empty source file: ", path)

  if (anyDuplicated(names(dat)) || anyNA(names(dat)) ||
      any(names(dat) == "")) {
    stop("Invalid column names in: ", path)
  }

  dat
}

find_column <- function(dat, options, label) {
  matches <- intersect(options, names(dat))

  if (length(matches) == 0) {
    stop(
      "Cannot identify ", label,
      "\nAvailable columns: ", paste(names(dat), collapse = ", ")
    )
  }

  matches[1]
}

numeric_checked <- function(x, label) {
  values <- suppressWarnings(as.numeric(as.character(x)))

  if (any(!is.finite(values)))
    stop("Invalid numeric values in ", label)

  values
}

write_text <- function(lines, filename) {
  writeLines(lines, file.path(output_dir, filename), useBytes = TRUE)
}

format_count <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

verification <- tibble(
  check = character(),
  observed = character(),
  expected = character(),
  status = character()
)

add_check <- function(name, observed, expected, passed) {
  verification[nrow(verification) + 1, ] <<- list(
    name,
    as.character(observed),
    as.character(expected),
    if (isTRUE(passed)) "PASS" else "FAIL"
  )
}


# 1. Dataset and contrast inventory

study_summary <- read_source("s06", "06_study_summary.csv")
contrasts <- read_source("s08", "08_contrast_inventory.csv")

dataset_col <- find_column(
  contrasts,
  c("dataset", "dataset_id", "study_id", "geo_accession", "gse"),
  "dataset identifier"
)

dataset_ids <- trimws(as.character(contrasts[[dataset_col]]))

if (anyNA(dataset_ids) || any(dataset_ids == ""))
  stop("Missing dataset identifiers in the contrast inventory.")

n_datasets <- length(unique(dataset_ids))
n_contrasts <- nrow(contrasts)

# Reference counts from the archived primary analysis.

add_check(
  "Number of drug-study contrasts",
  n_contrasts, 7, n_contrasts == 7
)

add_check(
  "Number of datasets",
  n_datasets, 5, n_datasets == 5
)


# 2. Eligible genes and shared candidates

eligible <- read_source("s09", "09_eligible_genes.csv")
candidates <- read_source("s06", "06_shared_candidate_genes.csv")

n_eligible <- nrow(eligible)
n_candidates <- nrow(candidates)

add_check(
  "Eligible genes",
  n_eligible, 16102, n_eligible == 16102
)

add_check(
  "Shared directional candidates",
  n_candidates, 926, n_candidates == 926
)

if (!"gene_id" %in% names(candidates))
  stop("Candidate table must contain gene_id.")

candidate_ids <- trimws(as.character(candidates$gene_id))

valid_ids <- all(
  !is.na(candidate_ids) &
    grepl("^ENSMUSG[0-9]{11}$", candidate_ids)
)

unique_ids <- !anyDuplicated(candidate_ids)

add_check(
  "Candidate identifier validity",
  valid_ids, TRUE, valid_ids
)

add_check(
  "Candidate identifier uniqueness",
  unique_ids, TRUE, unique_ids
)


# 3. Within-drug directional consistency

consistency <- read_source("s06", "06_consistency_summary.csv")

required_cols <- c(
  "drug", "n_studies", "n_genes",
  "n_consistent", "percentage_consistent"
)

missing_cols <- setdiff(required_cols, names(consistency))

if (length(missing_cols) > 0) {
  stop(
    "Missing consistency columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

consistency <- consistency %>%
  mutate(
    calculated_percentage = 100 * n_consistent / n_genes,
    percentage_difference = abs(
      calculated_percentage - percentage_consistent
    )
  )

consistency_valid <- all(
  is.finite(consistency$percentage_difference) &
    consistency$percentage_difference <= 0.15
)

add_check(
  "Directional consistency percentages",
  consistency_valid, TRUE, consistency_valid
)


# 4. QC sensitivity

qc <- read_source("s07", "07_candidate_comparison.csv")

status_col <- find_column(
  qc,
  c("status", "candidate_status", "comparison_status", "membership"),
  "QC membership"
)

qc_status <- tolower(trimws(as.character(qc[[status_col]])))

# Earlier versions of the QC analysis used different membership labels.

qc_status <- case_when(
  qc_status %in% c(
    "retained", "both", "shared", "stable",
    "in_both", "both analyses"
  ) ~ "retained",

  qc_status %in% c(
    "primary only", "primary_only", "primary-only", "primary"
  ) ~ "primary_only",

  qc_status %in% c(
    "sensitivity only", "sensitivity_only",
    "sensitivity-only", "sensitivity"
  ) ~ "sensitivity_only",

  TRUE ~ NA_character_
)

if (anyNA(qc_status))
  stop("Unrecognized QC membership values.")

qc_counts <- table(factor(
  qc_status,
  levels = c("retained", "primary_only", "sensitivity_only")
))

n_retained <- as.integer(qc_counts["retained"])
n_primary_only <- as.integer(qc_counts["primary_only"])
n_sensitivity_only <- as.integer(qc_counts["sensitivity_only"])

n_sensitivity <- n_retained + n_sensitivity_only
retention_percentage <- 100 * n_retained / n_candidates

add_check(
  "QC retained candidates",
  n_retained, 887, n_retained == 887
)

add_check(
  "QC primary-only candidates",
  n_primary_only, 39, n_primary_only == 39
)

add_check(
  "QC sensitivity-only candidates",
  n_sensitivity_only, 27, n_sensitivity_only == 27
)

add_check(
  "QC sensitivity candidate count",
  n_sensitivity, 914, n_sensitivity == 914
)

add_check(
  "QC primary candidate reconciliation",
  n_retained + n_primary_only,
  n_candidates,
  n_retained + n_primary_only == n_candidates
)


# 5. Permutation analysis

perms <- read_source("s08", "08_permutation_results.csv")
tests <- read_source("s08", "08_statistical_test_summary.csv")

stat_col <- find_column(
  tests, c("statistic", "test", "metric", "measure"),
  "statistic identifier"
)

observed_col <- find_column(
  tests, c("observed", "observed_count", "observed_value"),
  "observed statistic"
)

null_col <- find_column(
  tests, c("null_mean", "permutation_mean", "mean_null"),
  "null mean"
)

p_col <- find_column(
  tests, c("empirical_p", "empirical_p_value", "p_value"),
  "empirical p-value"
)

target <- which(as.character(tests[[stat_col]]) %in% c(
  "shared_direction_candidates",
  "shared_candidates",
  "shared_directional_candidates"
))

if (length(target) != 1)
  stop("Cannot uniquely identify the shared-candidate test.")

observed <- numeric_checked(
  tests[[observed_col]][target], "observed candidates"
)

null_mean <- numeric_checked(
  tests[[null_col]][target], "permutation null mean"
)

empirical_p <- numeric_checked(
  tests[[p_col]][target], "empirical p-value"
)

if (empirical_p < 0 || empirical_p > 1)
  stop("Empirical p-value is outside the valid range.")

n_permutations <- nrow(perms)

add_check(
  "Permutation count",
  n_permutations, 1000, n_permutations == 1000
)

add_check(
  "Observed permutation statistic",
  observed, n_candidates, observed == n_candidates
)

# Describe the reported result without assuming the archived outcome.

null_comparison <- if (observed > null_mean) {
  "The observed count was higher than the mean null count."
} else {
  "The observed count did not exceed the mean null count."
}

permutation_interpretation <- if (empirical_p < 0.05) {
  paste(
    "The upper-tail permutation test provided evidence of",
    "excess aggregate directional convergence under the specified null model."
  )
} else {
  paste(
    "The permutation test did not provide evidence of excess",
    "aggregate directional convergence under the specified null model."
  )
}


# 6. GO enrichment

go_all <- read_source("s09", "09_GO_BP_all_results.csv")

go_sig <- read_source(
  "s09", "09_GO_BP_FDR_significant.csv"
)

n_go_tested <- nrow(go_all)
n_go_significant <- nrow(go_sig)

add_check(
  "GO terms tested",
  n_go_tested, 6215, n_go_tested == 6215
)

add_check(
  "FDR-significant GO terms",
  n_go_significant, 4, n_go_significant == 4
)

go_id_col <- find_column(
  go_sig, c("go_id", "GO_ID", "GOID"), "GO identifier"
)

go_term_col <- find_column(
  go_sig, c("term", "go_term", "description", "TERM"),
  "GO term"
)

go_fdr_col <- find_column(
  go_sig, c("fdr", "FDR", "p_adjusted", "p_adj", "padj"),
  "GO FDR"
)

go_fold_col <- find_column(
  go_sig,
  c("fold_enrichment", "enrichment", "fold_change", "enrichment_ratio"),
  "GO fold enrichment"
)

go_hits_col <- find_column(
  go_sig,
  c("candidate_hits", "n_candidate_hits", "overlap", "n_hits"),
  "GO candidate hits"
)

go_export <- tibble(
  go_id = as.character(go_sig[[go_id_col]]),
  term = as.character(go_sig[[go_term_col]]),
  candidate_hits = numeric_checked(
    go_sig[[go_hits_col]], "GO candidate hits"
  ),
  fold_enrichment = numeric_checked(
    go_sig[[go_fold_col]], "GO fold enrichment"
  ),
  fdr = numeric_checked(go_sig[[go_fdr_col]], "GO FDR")
) %>%
  arrange(fdr)

add_check(
  "GO significance threshold",
  all(go_export$fdr >= 0 & go_export$fdr < 0.05),
  TRUE,
  all(go_export$fdr >= 0 & go_export$fdr < 0.05)
)

write_csv(
  go_export,
  file.path(output_dir, "13_manuscript_GO_results.csv")
)


# 7. Consolidated results

summary_table <- tribble(
  ~metric, ~value, ~source,

  "Datasets",
  n_datasets, "Script 08",

  "Drug-study contrasts",
  n_contrasts, "Script 08",

  "Eligible genes",
  n_eligible, "Script 09",

  "Shared directional candidates",
  n_candidates, "Script 06",

  "QC retained candidates",
  n_retained, "Script 07",

  "QC primary-only candidates",
  n_primary_only, "Script 07",

  "QC sensitivity-only candidates",
  n_sensitivity_only, "Script 07",

  "QC sensitivity candidates",
  n_sensitivity, "Script 07",

  "QC primary retention percentage",
  round(retention_percentage, 4), "Calculated from Scripts 06 and 07",

  "Permutations",
  n_permutations, "Script 08",

  "Observed shared candidates",
  observed, "Script 08",

  "Permutation null mean",
  round(null_mean, 5), "Script 08",

  "Empirical upper-tail p-value",
  empirical_p, "Script 08",

  "GO biological-process terms tested",
  n_go_tested, "Script 09",

  "FDR-significant GO terms",
  n_go_significant, "Script 09"
) %>%
  mutate(value = as.character(value))

write_csv(
  summary_table,
  file.path(output_dir, "13_manuscript_results_summary.csv")
)


# 8. Manuscript Results section

consistency_lines <- with(
  consistency,
  sprintf(
    "%s: %s of %s genes (%.1f%%) across %s studies.",
    drug,
    format_count(n_consistent),
    format_count(n_genes),
    percentage_consistent,
    n_studies
  )
)

go_lines <- with(
  go_export,
  sprintf(
    "%s (%s): %s candidate hits; fold enrichment %.2f; FDR %.4g.",
    term, go_id, candidate_hits, fold_enrichment, fdr
  )
)

results_lines <- c(
  "RESULTS",
  "",
  "Study integration",
  "",
  sprintf(
    paste(
      "The exploratory analysis integrated %d publicly available",
      "rodent striatal transcriptomic datasets, representing %d",
      "haloperidol- or risperidone-related drug-study contrasts."
    ),
    n_datasets, n_contrasts
  ),
  "",
  "Directional consistency",
  "",
  sprintf(
    paste(
      "A total of %s genes met the study-coverage criteria",
      "for the cross-drug exploratory analysis. Of these, %s",
      "satisfied the prespecified within-drug and cross-drug",
      "directional consistency criteria."
    ),
    format_count(n_eligible),
    format_count(n_candidates)
  ),
  "",
  "Within-drug directional consistency by study coverage:",
  consistency_lines,
  "",
  paste(
    "These percentages are descriptive and should not be",
    "interpreted as evidence of statistically significant convergence."
  ),
  "",
  "QC sensitivity",
  "",
  sprintf(
    paste(
      "Of the %d primary exploratory candidates, %d (%.2f%%)",
      "were retained under the alternative GSE66275 sample-selection",
      "analysis. A further %d candidates were primary-only and %d",
      "were sensitivity-only, yielding %d sensitivity candidates."
    ),
    n_candidates, n_retained, retention_percentage,
    n_primary_only, n_sensitivity_only, n_sensitivity
  ),
  "",
  "Permutation assessment",
  "",
  sprintf(
    paste(
      "The observed shared directional candidate count was %d.",
      "Across %d study-aware permutations, the mean null count",
      "was %.2f. The empirical one-sided upper-tail p-value was %.4f."
    ),
    observed, n_permutations, null_mean, empirical_p
  ),
  null_comparison,
  "",
  permutation_interpretation,
  "",
  "Exploratory GO biological-process enrichment",
  "",
  sprintf(
    paste(
      "Of %s GO biological-process terms tested,",
      "%d met the specified BH-FDR threshold."
    ),
    format_count(n_go_tested), n_go_significant
  ),
  "",
  go_lines,
  "",
  paste(
    "These GO terms were identified using an overrepresentation",
    "model applied to a directionally selected exploratory gene set.",
    "They do not independently establish statistically significant",
    "cross-drug convergence."
  )
)

write_text(results_lines, "13_manuscript_results.txt")


# 9. Figure legends

figure_legends <- c(
  "FIGURE LEGENDS",
  "",
  "Figure 1. Study integration and analytical workflow.",
  paste(
    "Overview of public dataset integration, gene harmonization,",
    "exploratory candidate selection, QC sensitivity, permutation",
    "testing and GO biological-process enrichment."
  ),
  "",
  "Figure 2. Within-drug directional consistency.",
  paste(
    "Proportions of directionally consistent genes for haloperidol",
    "and risperidone across different numbers of contributing studies.",
    "Denominators vary by study coverage. Results are descriptive."
  ),
  "",
  "Figure 3. Cross-drug expression changes.",
  sprintf(
    paste(
      "Unweighted mean log2 fold changes for the %d shared",
      "exploratory candidates. As genes were selected for matching",
      "directions, the plotted relationship is not independent",
      "evidence of cross-drug convergence."
    ),
    n_candidates
  ),
  "",
  "Figure 4. QC sensitivity of exploratory candidates.",
  sprintf(
    paste(
      "%d primary candidates were retained under the alternative",
      "GSE66275 sample-selection analysis; %d were primary-only",
      "and %d were sensitivity-only."
    ),
    n_retained, n_primary_only, n_sensitivity_only
  ),
  "",
  "Figure 5. Permutation assessment of convergence.",
  sprintf(
    paste(
      "Distribution of shared candidate counts from %d study-aware",
      "permutations. The observed count was %d and the null mean",
      "was %.2f (one-sided empirical p = %.4f)."
    ),
    n_permutations, observed, null_mean, empirical_p
  ),
  permutation_interpretation,
  "",
  "Figure 6. Exploratory GO enrichment.",
  sprintf(
    paste(
      "Fold enrichment of the %d GO biological-process terms",
      "meeting the specified BH-FDR criterion. Overlapping terms",
      "are not independent biological discoveries."
    ),
    n_go_significant
  ),
  "",
  "Supplementary Figure S1. Selected exploratory genes.",
  paste(
    "Unweighted mean expression changes for selected genes across",
    "haloperidol and risperidone comparisons. These are descriptive",
    "estimates rather than FDR-significant gene-level discoveries."
  )
)

write_text(figure_legends, "13_figure_legends.txt")


# 10. Supplementary table legends

supplementary_legends <- c(
  "SUPPLEMENTARY TABLE LEGENDS",
  "",
  "Table S1A. Study-level dataset summary.",
  "Table S1B. Drug-study contrast inventory.",
  "Table S2A. GO annotation coverage.",
  "Table S2B. Genes eligible for exploratory selection.",
  "Table S3. Within-drug directional consistency.",
  "Table S4. Shared exploratory directional candidates.",
  "Table S5A. QC candidate membership comparison.",
  "Table S5B. QC gene-stability output.",
  "Table S6A. Permutation-test summary.",
  "Table S6B. Individual permutation results.",
  "Table S7A. Complete GO biological-process results.",
  "Table S7B. FDR-significant GO biological-process terms.",
  "Table S8. Candidate genes contributing to enriched GO terms.",
  "Table S9. Selected genes of biological interest.",
  "",
  paste(
    "All tables preserve the original analytical results.",
    "Candidate membership is exploratory and does not imply",
    "gene-level statistical significance."
  )
)

write_text(
  supplementary_legends,
  "13_supplementary_table_legends.txt"
)


# 11. Methodological limitations

limitations <- c(
  "METHODOLOGICAL LIMITATIONS",
  "",
  "1. Candidate selection",
  sprintf(
    paste(
      "The %d candidates were identified using directional",
      "consistency criteria rather than gene-level FDR-significant",
      "meta-analysis."
    ),
    n_candidates
  ),
  "",
  "2. Permutation inference",
  permutation_interpretation,
  paste(
    "Inference depends on the assumptions of the implemented",
    "gene-label permutation null model."
  ),
  "",
  "3. GO enrichment",
  paste(
    "GO overrepresentation was performed on an exploratory",
    "directionally selected gene set. Gene dependence and",
    "selection effects may influence enrichment inference."
  ),
  "",
  "4. Source contrast provenance",
  paste(
    "Original GEO treatment assignments, contrast definitions,",
    "vehicle metadata and treatment durations require independent",
    "source verification."
  ),
  "",
  "5. Ortholog mapping",
  paste(
    "Some rat-to-mouse mappings were supported by gene-symbol-based",
    "annotation rather than independent rat NCBI ortholog verification."
  ),
  "",
  "6. Effect estimates and variances",
  paste(
    "The provenance of all source log2FC and variance estimates",
    "has not been independently authenticated."
  ),
  "",
  "7. Reproducibility",
  paste(
    "Original analysis scripts were not consistently saved during",
    "execution. Consolidated outputs are available, but complete",
    "clean-session reproduction remains outstanding."
  ),
  "",
  "8. Biological generalizability",
  paste(
    "Rodent striatal transcriptional responses cannot be assumed",
    "to generalize directly to human clinical drug response."
  )
)

write_text(limitations, "13_manuscript_limitations.txt")


# 12. Final verification

write_csv(
  verification,
  file.path(output_dir, "13_results_verification.csv")
)

write_text(
  capture.output(sessionInfo()),
  "13_session_info.txt"
)

n_passed <- sum(verification$status == "PASS")
n_failed <- sum(verification$status == "FAIL")

cat("\nManuscript results export\n")
print(verification, n = Inf)

cat("\nPassed:", n_passed, "\n")
cat("Failed:", n_failed, "\n")
cat("Output directory:", output_dir, "\n")

if (n_failed > 0) {
  warning(
    "Some verification checks failed. Review ",
    "13_results_verification.csv before using the manuscript text."
  )
} else {
  cat("\nManuscript export completed without failed checks.\n")
  cat("Scientific provenance limitations still apply.\n")
}
