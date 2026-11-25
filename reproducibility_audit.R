# Reproducibility and provenance audit

# Checks existing results, supplementary tables and source files.
# No differential expression or statistical analyses are rerun.

project_dir <- getwd()
results_dir <- file.path(project_dir, "results_v2")
output_dir <- file.path(results_dir, "12_reproducibility_audit")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source_dirs <- c(
  s06 = file.path(results_dir, "06_exploratory_analysis", "tables"),
  s07 = file.path(results_dir, "07_QC_sensitivity", "tables"),
  s08 = file.path(results_dir, "08_statistical_convergence", "tables"),
  s09 = file.path(results_dir, "09_biological_interpretation", "tables"),
  s11 = file.path(results_dir, "11_supplementary_tables", "tables")
)


# Helpers

read_csv_checked <- function(path) {
  if (!file.exists(path))
    stop("Required file missing: ", path)

  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if (nrow(dat) == 0)
    stop("Empty input table: ", path)

  if (anyDuplicated(names(dat)) || anyNA(names(dat)) ||
      any(names(dat) == "")) {
    stop("Invalid column names: ", path)
  }

  dat
}

source_file <- function(script, filename) {
  file.path(source_dirs[[script]], filename)
}

read_source <- function(script, filename) {
  read_csv_checked(source_file(script, filename))
}

get_column <- function(dat, options) {
  found <- intersect(options, names(dat))
  if (length(found) == 0) return(NULL)
  found[1]
}

numeric_column <- function(dat, column) {
  if (is.null(column)) return(NULL)

  x <- dat[[column]]
  if (is.factor(x)) x <- as.character(x)

  suppressWarnings(as.numeric(x))
}

audit <- data.frame(
  category = character(),
  check = character(),
  status = character(),
  observed = character(),
  expected = character(),
  explanation = character(),
  stringsAsFactors = FALSE
)

add_check <- function(category, check, status,
                      observed, expected, explanation = "") {

  stopifnot(status %in% c(
    "PASS", "FAIL", "NOT VERIFIED", "NOT APPLICABLE"
  ))

  audit[nrow(audit) + 1, ] <<- list(
    category,
    check,
    status,
    paste(observed, collapse = "; "),
    paste(expected, collapse = "; "),
    explanation
  )
}

check_condition <- function(category, check, condition,
                            observed, expected, explanation = "") {

  add_check(
    category, check,
    if (isTRUE(condition)) "PASS" else "FAIL",
    observed, expected, explanation
  )
}


# 1. Primary candidate identifiers

candidates <- read_source("s06", "06_shared_candidate_genes.csv")

check_condition(
  "Candidates", "Primary candidate count",
  nrow(candidates) == 926,
  nrow(candidates), 926
)

if ("gene_id" %in% names(candidates)) {

  ids <- trimws(as.character(candidates$gene_id))

  missing_ids <- sum(is.na(ids) | ids == "")
  valid_ids <- ids[!is.na(ids) & ids != ""]

  duplicate_ids <- sum(duplicated(valid_ids))
  invalid_ids <- sum(!grepl("^ENSMUSG[0-9]{11}$", valid_ids))

  check_condition(
    "Candidates", "Missing gene identifiers",
    missing_ids == 0, missing_ids, 0
  )

  check_condition(
    "Candidates", "Duplicate gene identifiers",
    duplicate_ids == 0, duplicate_ids, 0
  )

  check_condition(
    "Candidates", "Mouse Ensembl identifier format",
    invalid_ids == 0, invalid_ids, 0,
    "Checks identifier syntax, not biological mapping."
  )

} else {

  add_check(
    "Candidates", "Gene identifier verification", "FAIL",
    "gene_id column unavailable", "gene_id column",
    "Expected identifier column was not found."
  )
}


# 2. Directional consistency

consistency <- read_source("s06", "06_consistency_summary.csv")

required_cols <- c(
  "drug", "n_studies", "n_genes",
  "n_consistent", "percentage_consistent"
)

if (all(required_cols %in% names(consistency))) {

  calculated <- with(consistency, 100 * n_consistent / n_genes)
  difference <- abs(calculated - consistency$percentage_consistent)

  valid_counts <- with(
    consistency,
    all(
      is.finite(n_genes) &
      is.finite(n_consistent) &
      n_genes > 0 &
      n_consistent >= 0 &
      n_consistent <= n_genes
    )
  )

  check_condition(
    "Directional consistency", "Valid gene counts",
    valid_counts, nrow(consistency),
    "All rows contain valid counts"
  )

  check_condition(
    "Directional consistency", "Reported percentages",
    all(is.finite(difference)) && all(difference <= 0.15),
    round(max(difference), 5),
    "Maximum discrepancy <= 0.15 percentage points",
    "Allows rounding of reported percentages."
  )

  check_condition(
    "Directional consistency", "Expected summary rows",
    nrow(consistency) == 5, nrow(consistency), 5
  )

} else {

  add_check(
    "Directional consistency", "Summary schema", "FAIL",
    names(consistency), required_cols,
    "Required columns unavailable."
  )
}


# 3. QC sensitivity

qc <- read_source("s07", "07_candidate_comparison.csv")

status_col <- get_column(
  qc, c("status", "candidate_status", "comparison_status", "membership")
)

if (!is.null(status_col)) {

  statuses <- tolower(trimws(as.character(qc[[status_col]])))

  # Standardize labels used in earlier versions of the QC analysis.

  membership_labels <- list(
    retained = c(
      "retained", "both", "shared", "stable",
      "in_both", "both analyses"
    ),
    primary_only = c(
      "primary only", "primary_only", "primary-only", "primary"
    ),
    sensitivity_only = c(
      "sensitivity only", "sensitivity_only",
      "sensitivity-only", "sensitivity"
    )
  )

  membership <- rep(NA_character_, length(statuses))

  for (label in names(membership_labels)) {
    membership[statuses %in% membership_labels[[label]]] <- label
  }

  if (anyNA(membership)) {

    add_check(
      "QC sensitivity", "Candidate membership categories",
      "NOT VERIFIED", unique(statuses),
      "Recognized primary/sensitivity categories",
      "Unrecognized labels require inspection."
    )

  } else {

    counts <- table(factor(membership, levels = names(membership_labels)))
    expected <- c(887L, 39L, 27L)

    check_condition(
      "QC sensitivity", "Candidate membership counts",
      identical(as.integer(counts), expected),
      as.integer(counts), expected
    )

    primary_total <- sum(counts[c("retained", "primary_only")])

    sensitivity_total <- sum(
      counts[c("retained", "sensitivity_only")]
    )

    check_condition(
      "QC sensitivity", "Primary candidate total",
      primary_total == 926, primary_total, 926
    )

    check_condition(
      "QC sensitivity", "Sensitivity candidate total",
      sensitivity_total == 914, sensitivity_total, 914
    )

    retention <- as.numeric(counts["retained"]) / primary_total

    check_condition(
      "QC sensitivity", "Primary candidate retention",
      abs(retention - 887 / 926) < 1e-10,
      round(100 * retention, 4),
      round(100 * 887 / 926, 4)
    )
  }

} else {

  add_check(
    "QC sensitivity", "Candidate membership", "NOT VERIFIED",
    names(qc), "Recognized membership column",
    "Column mapping requires inspection."
  )
}


# 4. Permutation results

perms <- read_source("s08", "08_permutation_results.csv")
tests <- read_source("s08", "08_statistical_test_summary.csv")

check_condition(
  "Permutation", "Number of permutations",
  nrow(perms) == 1000, nrow(perms), 1000
)

perm_col <- get_column(
  perms,
  c(
    "shared_direction_candidates", "shared_candidates",
    "n_shared_candidates", "shared_direction_consistent",
    "shared_directional_candidates"
  )
)

stat_col <- get_column(tests, c("statistic", "test", "metric", "measure"))

observed_col <- get_column(
  tests, c("observed", "observed_count", "observed_value")
)

null_col <- get_column(
  tests, c("null_mean", "permutation_mean", "mean_null")
)

p_col <- get_column(
  tests, c("empirical_p", "empirical_p_value", "p_value")
)

perm_columns <- list(
  perm_col, stat_col, observed_col, null_col, p_col
)

if (all(!vapply(perm_columns, is.null, logical(1)))) {

  test_names <- as.character(tests[[stat_col]])

  target <- which(test_names %in% c(
    "shared_direction_candidates",
    "shared_candidates",
    "shared_directional_candidates"
  ))

  if (length(target) == 1) {

    perm_values <- numeric_column(perms, perm_col)
    observed <- numeric_column(tests, observed_col)[target]
    reported_mean <- numeric_column(tests, null_col)[target]
    reported_p <- numeric_column(tests, p_col)[target]

    valid_perms <- all(is.finite(perm_values))

    check_condition(
      "Permutation", "Valid permutation counts",
      valid_perms &&
        all(perm_values >= 0) &&
        all(perm_values == floor(perm_values)),
      sum(is.finite(perm_values)), nrow(perms)
    )

    check_condition(
      "Permutation", "Observed shared candidates",
      is.finite(observed) && observed == 926,
      observed, 926
    )

    if (valid_perms) {

      calculated_mean <- mean(perm_values)

      check_condition(
        "Permutation", "Reported null mean",
        is.finite(reported_mean) &&
          abs(calculated_mean - reported_mean) <= 0.01,
        round(calculated_mean, 5), reported_mean
      )

      # Recalculate the upper-tail empirical p-value with +1 correction.

      calculated_p <- (1 + sum(perm_values >= observed)) /
        (length(perm_values) + 1)

      check_condition(
        "Permutation", "Empirical upper-tail p-value",
        is.finite(reported_p) &&
          abs(calculated_p - reported_p) <= 0.001,
        round(calculated_p, 6), reported_p,
        "Assumes Script 08 used the +1 correction."
      )
    }

  } else {

    add_check(
      "Permutation", "Shared-candidate statistic", "NOT VERIFIED",
      test_names, "One matching shared-candidate statistic",
      "Statistic identifier requires inspection."
    )
  }

} else {

  add_check(
    "Permutation", "Statistical output schema", "NOT VERIFIED",
    "One or more expected columns unavailable",
    "Recognized permutation and summary columns",
    "Requires comparison with Script 08."
  )
}


# 5. GO enrichment

go_all <- read_source("s09", "09_GO_BP_all_results.csv")
go_sig <- read_source("s09", "09_GO_BP_FDR_significant.csv")
eligible <- read_source("s09", "09_eligible_genes.csv")

check_condition(
  "GO enrichment", "Eligible genes",
  nrow(eligible) == 16102, nrow(eligible), 16102
)

check_condition(
  "GO enrichment", "GO terms tested",
  nrow(go_all) == 6215, nrow(go_all), 6215
)

check_condition(
  "GO enrichment", "Reported significant GO terms",
  nrow(go_sig) == 4, nrow(go_sig), 4
)

id_options <- c("go_id", "GO_ID", "GOID")
fdr_options <- c("fdr", "FDR", "p_adjusted", "p_adj", "padj")

go_id_col <- get_column(go_all, id_options)
sig_id_col <- get_column(go_sig, id_options)

fdr_col <- get_column(go_all, fdr_options)
sig_fdr_col <- get_column(go_sig, fdr_options)

hits_col <- get_column(
  go_all,
  c("candidate_hits", "n_candidate_hits", "overlap", "n_hits")
)

go_columns <- list(go_id_col, sig_id_col, fdr_col, sig_fdr_col)

if (all(!vapply(go_columns, is.null, logical(1)))) {

  all_ids <- as.character(go_all[[go_id_col]])
  sig_ids <- as.character(go_sig[[sig_id_col]])

  check_condition(
    "GO enrichment", "Significant terms represented in full results",
    all(sig_ids %in% all_ids),
    sum(sig_ids %in% all_ids), length(sig_ids)
  )

  check_condition(
    "GO enrichment", "Unique GO identifiers",
    !anyDuplicated(all_ids) && !anyDuplicated(sig_ids),
    paste(sum(duplicated(all_ids)), sum(duplicated(sig_ids)), sep = ", "),
    "0, 0"
  )

  sig_fdr <- numeric_column(go_sig, sig_fdr_col)

  check_condition(
    "GO enrichment", "Significant-term FDR threshold",
    all(is.finite(sig_fdr)) && all(sig_fdr < 0.05),
    signif(sig_fdr, 4), "All FDR < 0.05"
  )

  if (!is.null(hits_col)) {

    full_fdr <- numeric_column(go_all, fdr_col)
    full_hits <- numeric_column(go_all, hits_col)

    expected_ids <- all_ids[
      !is.na(full_fdr) &
      !is.na(full_hits) &
      full_fdr < 0.05 &
      full_hits >= 3
    ]

    check_condition(
      "GO enrichment", "Significant-term selection",
      setequal(expected_ids, sig_ids) &&
        length(expected_ids) == length(sig_ids),
      length(expected_ids), length(sig_ids),
      "Checks FDR < 0.05 and at least three candidate hits."
    )

  } else {

    add_check(
      "GO enrichment", "Significant-term selection", "NOT VERIFIED",
      "Candidate-hit column unavailable",
      "FDR < 0.05 and at least three hits"
    )
  }

} else {

  add_check(
    "GO enrichment", "GO result matching", "NOT VERIFIED",
    "Required identifier or FDR column unavailable",
    "Matching GO identifiers and FDR values"
  )
}


# 6. Supplementary table consistency

# Compare selected supplementary tables with their analytical sources.

supplement_pairs <- data.frame(
  script = c("s06", "s08", "s09", "s09"),
  original = c(
    "06_shared_candidate_genes.csv",
    "08_permutation_results.csv",
    "09_GO_BP_all_results.csv",
    "09_GO_BP_FDR_significant.csv"
  ),
  supplement = c(
    "Table_S4_shared_candidates.csv",
    "Table_S6B_permutation_results.csv",
    "Table_S7A_GO_BP_all_results.csv",
    "Table_S7B_GO_BP_FDR_significant.csv"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(supplement_pairs))) {

  pair <- supplement_pairs[i, ]

  original_path <- source_file(pair$script, pair$original)
  supplement_path <- source_file("s11", pair$supplement)

  if (file.exists(original_path) && file.exists(supplement_path)) {

    original <- read_csv_checked(original_path)
    supplement <- read_csv_checked(supplement_path)

    matching <- identical(names(original), names(supplement)) &&
      isTRUE(all.equal(original, supplement, check.attributes = FALSE))

    check_condition(
      "Supplementary tables", pair$supplement,
      matching,
      if (matching) "Matching content" else "Content differs",
      "Matching source table",
      "Compares column names, parsed values and row order."
    )

  } else {

    add_check(
      "Supplementary tables", pair$supplement, "NOT VERIFIED",
      "Source or supplementary file unavailable",
      "Both files available"
    )
  }
}


# 7. Master dataset availability

master_files <- c(
  "data/processed/MASTER_primary_MGI_v2_harmonized.rds",
  "data/processed/MASTER_QC_sensitivity_MGI_v2_harmonized.rds"
)

for (path in master_files) {

  present <- file.exists(file.path(project_dir, path))

  check_condition(
    "Input availability", basename(path),
    present,
    if (present) "Present" else "Missing",
    "Present",
    "Availability does not establish source provenance."
  )
}


# 8. Script availability

script_dir <- file.path(project_dir, "scripts_v2")

scripts_found <- if (dir.exists(script_dir)) {
  list.files(script_dir, pattern = "\\.[Rr]$", full.names = TRUE)
} else {
  character()
}

script_inventory <- data.frame(
  filename = basename(scripts_found),
  path = scripts_found,
  size_bytes = if (length(scripts_found)) {
    file.info(scripts_found)$size
  } else {
    numeric()
  },
  stringsAsFactors = FALSE
)

write.csv(
  script_inventory,
  file.path(output_dir, "12_script_inventory.csv"),
  row.names = FALSE
)

for (script_number in 6:11) {

  pattern <- paste0("^0?", script_number, "(_|\\.)")
  matches <- grepl(pattern, basename(scripts_found))

  check_condition(
    "Script availability", paste("Script", script_number),
    any(matches), sum(matches),
    "At least one matching R script",
    "Checks presence, not independent executability."
  )
}


# 9. Package environment

packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "metafor",
  "AnnotationDbi", "org.Mm.eg.db", "GO.db"
)

installed <- vapply(
  packages, requireNamespace, logical(1), quietly = TRUE
)

versions <- vapply(
  packages,
  function(pkg) {
    if (installed[[pkg]]) {
      as.character(packageVersion(pkg))
    } else {
      NA_character_
    }
  },
  character(1)
)

package_inventory <- data.frame(
  package = packages,
  installed = unname(installed),
  version = unname(versions),
  stringsAsFactors = FALSE
)

write.csv(
  package_inventory,
  file.path(output_dir, "12_package_inventory.csv"),
  row.names = FALSE
)

check_condition(
  "Environment", "Package dependency availability",
  all(installed), sum(installed), length(packages),
  "Checks selected analysis packages, not every dependency."
)


# 10. Scientific provenance

# These questions cannot be resolved using consolidated outputs alone.

provenance <- data.frame(
  item = c(
    "Original GEO treatment and control assignments",
    "Biological contrast definitions",
    "Treatment duration and vehicle metadata",
    "Source log2FC provenance",
    "Source variance and standard-error provenance",
    "Rat-to-mouse ortholog verification",
    "Cross-study independence",
    "Full workflow reproducibility"
  ),
  status = rep("NOT VERIFIED", 8),
  explanation = c(
    "Requires checking original GEO sample metadata.",
    "Requires checking source contrast definitions.",
    "Requires source-level experimental metadata review.",
    "Requires tracing estimates to original processed results.",
    "Requires tracing variance calculations to source data.",
    "Some rat mappings remain symbol-based rather than independently verified.",
    "Requires review of study and sample overlap.",
    "Requires a clean-session run from documented inputs."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  provenance,
  file.path(output_dir, "12_provenance_register.csv"),
  row.names = FALSE
)

for (i in seq_len(nrow(provenance))) {

  add_check(
    "Scientific provenance",
    provenance$item[i],
    provenance$status[i],
    "Not independently authenticated",
    "Verified source documentation",
    provenance$explanation[i]
  )
}


# 11. Source file inventory and checksums

source_paths <- unlist(
  lapply(source_dirs, function(folder) {
    if (dir.exists(folder)) {
      list.files(folder, pattern = "\\.csv$", full.names = TRUE)
    } else {
      character()
    }
  }),
  use.names = FALSE
)

source_paths <- sort(unique(source_paths))

file_inventory <- data.frame(
  path = source_paths,
  size_bytes = file.info(source_paths)$size,
  md5 = unname(tools::md5sum(source_paths)),
  stringsAsFactors = FALSE
)

write.csv(
  file_inventory,
  file.path(output_dir, "12_source_file_inventory.csv"),
  row.names = FALSE
)


# 12. Audit report

write.csv(
  audit,
  file.path(output_dir, "12_audit_checks.csv"),
  row.names = FALSE
)

status_counts <- table(factor(
  audit$status,
  levels = c("PASS", "FAIL", "NOT VERIFIED", "NOT APPLICABLE")
))

failed <- audit[audit$status == "FAIL", , drop = FALSE]

unverified <- audit[
  audit$status == "NOT VERIFIED", , drop = FALSE
]

report <- c(
  "Reproducibility and provenance audit",
  paste("Date:", Sys.Date()),
  paste("Project directory:", project_dir),
  "",
  "Audit summary",
  paste("Total checks:", nrow(audit)),
  paste("Passed:", status_counts["PASS"]),
  paste("Failed:", status_counts["FAIL"]),
  paste("Not verified:", status_counts["NOT VERIFIED"]),
  paste("Not applicable:", status_counts["NOT APPLICABLE"]),
  "",
  "Interpretation",
  "PASS: the specified computational check succeeded.",
  "FAIL: a discrepancy requires investigation.",
  "NOT VERIFIED: insufficient evidence to complete the check.",
  "",
  "Limitations",
  "This audit does not independently authenticate GEO metadata,",
  "verify every ortholog mapping or validate source variance estimates.",
  "Full reproducibility in a clean environment remains unverified.",
  "Exploratory candidates are not statistically significant gene-level discoveries.",
  "",
  "Historical results used as reference benchmarks",
  "Primary exploratory candidates: 926",
  "QC-retained primary candidates: 887",
  "Permutation analysis: no excess aggregate convergence",
  "GO enrichment: four overlapping exploratory terms",
  "These are archived findings, not new statistical conclusions.",
  "",
  "Next steps",
  "Review failed and unverified checks before manuscript preparation.",
  "Reconcile discrepancies with the original analytical outputs."
)

writeLines(report, file.path(output_dir, "12_audit_summary.txt"))

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "12_session_info.txt")
)

cat("\nReproducibility audit\n")
print(status_counts)

cat("\nFailed checks\n")

if (nrow(failed) > 0) {
  print(
    failed[, c("category", "check", "observed", "expected")],
    row.names = FALSE
  )
} else {
  cat("No failed checks.\n")
}

cat("\nUnverified checks\n")

if (nrow(unverified) > 0) {
  print(unverified[, c("category", "check")], row.names = FALSE)
} else {
  cat("No unverified checks.\n")
}

cat("\nOutput directory:", output_dir, "\n")

if (nrow(failed) > 0) {
  warning("Audit completed with failed checks. Review 12_audit_checks.csv.")
} else {
  cat("\nAudit completed without failed checks.\n")
  cat("Unverified provenance items may still remain.\n")
}
