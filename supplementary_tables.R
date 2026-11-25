# Script 11: Supplementary tables and results check
# Antipsychotic transcriptional responses in the striatum
#
# Gathers the outputs from scripts 06-09 into supplementary tables and
# runs a few sanity checks on them. No tests are rerun here.

library(dplyr)
library(readr)
library(tibble)

res_dir <- "results_v2"
out_dir <- file.path(res_dir, "11_supplementary_tables")
table_dir <- file.path(out_dir, "tables")
report_dir <- file.path(out_dir, "reports")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

# where each script keeps its tables
src_dir <- c(
  s06 = file.path(res_dir, "06_exploratory_analysis", "tables"),
  s07 = file.path(res_dir, "07_QC_sensitivity", "tables"),
  s08 = file.path(res_dir, "08_statistical_convergence", "tables"),
  s09 = file.path(res_dir, "09_biological_interpretation", "tables")
)


# Helpers -----------------------------------------------------------------

read_source <- function(script, file) {
  path <- file.path(src_dir[[script]], file)

  if (!file.exists(path)) stop("Can't find source file: ", path)

  dat <- read_csv(path, show_col_types = FALSE, name_repair = "minimal")

  if (nrow(dat) == 0) stop("No rows in: ", path)
  if (anyDuplicated(names(dat)) > 0 || anyNA(names(dat)) || any(names(dat) == "")) {
    stop("Bad or duplicated column names in: ", path)
  }

  dat
}

check_columns <- function(dat, cols, label) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) {
    stop(label, " is missing: ", paste(missing, collapse = ", "),
         "\nColumns present: ", paste(names(dat), collapse = ", "))
  }
}

make_check <- function(name, passed, observed, expected, note = "") {
  tibble(
    check = name,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected),
    note = note
  )
}


# Tables to export --------------------------------------------------------

supp <- tribble(
  ~out_file,                                   ~script, ~in_file,                                   ~description,
  "Table_S1A_study_summary.csv",               "s06",   "06_study_summary.csv",                     "Study-level dataset summary",
  "Table_S1B_contrast_inventory.csv",          "s08",   "08_contrast_inventory.csv",                "Drug-study contrast inventory",
  "Table_S2A_annotation_coverage.csv",         "s09",   "09_annotation_coverage.csv",               "GO annotation coverage",
  "Table_S2B_eligible_genes.csv",              "s09",   "09_eligible_genes.csv",                    "Genes eligible for candidate selection",
  "Table_S3_directional_consistency.csv",      "s06",   "06_consistency_summary.csv",               "Within-drug directional consistency by study coverage",
  "Table_S4_shared_candidates.csv",            "s06",   "06_shared_candidate_genes.csv",            "Exploratory shared directional candidates",
  "Table_S5A_QC_candidate_comparison.csv",     "s07",   "07_candidate_comparison.csv",              "Candidate membership under primary and sensitivity analyses",
  "Table_S5B_QC_gene_stability.csv",           "s07",   "07_gene_stability.csv",                    "Gene-level QC sensitivity information",
  "Table_S6A_permutation_test_summary.csv",    "s08",   "08_statistical_test_summary.csv",          "Observed statistics and permutation-based comparisons",
  "Table_S6B_permutation_results.csv",         "s08",   "08_permutation_results.csv",               "Individual permutation results",
  "Table_S7A_GO_BP_all_results.csv",           "s09",   "09_GO_BP_all_results.csv",                 "Complete exploratory GO biological-process enrichment results",
  "Table_S7B_GO_BP_FDR_significant.csv",       "s09",   "09_GO_BP_FDR_significant.csv",             "GO terms meeting the BH-FDR criterion",
  "Table_S8_significant_GO_genes.csv",         "s09",   "09_significant_GO_contributing_genes.csv", "Candidate genes contributing to enriched GO terms",
  "Table_S9_genes_of_interest.csv",            "s09",   "09_genes_of_interest.csv",                 "Selected exploratory genes with descriptive expression estimates"
)

# read everything in and write out the copies
tabs <- list()
manifest <- vector("list", nrow(supp))

for (i in seq_len(nrow(supp))) {
  row <- supp[i, ]
  dat <- read_source(row$script, row$in_file)

  # the consistency table gets used again in the checks below
  if (row$in_file == "06_consistency_summary.csv") {
    check_columns(
      dat,
      c("drug", "n_studies", "n_genes", "n_consistent", "percentage_consistent"),
      "Script 06 consistency summary"
    )
  }

  out_path <- file.path(table_dir, row$out_file)
  write_csv(dat, out_path, na = "")

  tabs[[row$in_file]] <- dat
  manifest[[i]] <- tibble(
    filename = row$out_file,
    description = row$description,
    source = paste("Script", sub("s", "", row$script)),
    rows = nrow(dat),
    columns = ncol(dat),
    output_path = out_path
  )
}

manifest <- bind_rows(manifest)


# Integrity checks --------------------------------------------------------

candidates <- tabs[["06_shared_candidate_genes.csv"]]
perms <- tabs[["08_permutation_results.csv"]]
go_all <- tabs[["09_GO_BP_all_results.csv"]]
go_sig <- tabs[["09_GO_BP_FDR_significant.csv"]]
eligible <- tabs[["09_eligible_genes.csv"]]
consistency <- tabs[["06_consistency_summary.csv"]]

checks <- list(
  make_check("Primary shared candidate count", nrow(candidates) == 926,
             nrow(candidates), 926, "Expected count from Script 06"),
  make_check("Permutation count", nrow(perms) == 1000,
             nrow(perms), 1000, "Expected number of permutations"),
  make_check("Number of GO terms tested", nrow(go_all) == 6215,
             nrow(go_all), 6215, "Expected from Script 09"),
  make_check("Number of FDR-significant GO terms", nrow(go_sig) == 4,
             nrow(go_sig), 4, "Expected from Script 09"),
  make_check("Eligible gene count", nrow(eligible) == 16102,
             nrow(eligible), 16102, "Expected eligible gene count")
)

# duplicate gene IDs in the candidate list (column name has varied)
id_col <- intersect(c("gene_ensembl_id", "harmonized_mouse_gene_id"), names(candidates))[1]

if (!is.na(id_col)) {
  n_dup <- sum(duplicated(candidates[[id_col]]))
  checks[[length(checks) + 1]] <- make_check(
    "Duplicate candidate gene identifiers", n_dup == 0, n_dup, 0,
    paste("Uses", id_col)
  )
} else {
  checks[[length(checks) + 1]] <- tibble(
    check = "Duplicate candidate gene identifiers",
    status = "NOT CHECKED",
    observed = "Identifier column unavailable",
    expected = "Unique gene identifiers",
    note = "Review candidate table schema"
  )
}

# recompute the consistency percentages and compare to what was reported
max_diff <- consistency |>
  mutate(diff = abs(100 * n_consistent / n_genes - percentage_consistent)) |>
  pull(diff) |>
  max(na.rm = TRUE)

checks[[length(checks) + 1]] <- make_check(
  "Directional consistency percentages",
  is.finite(max_diff) && max_diff <= 0.15,
  round(max_diff, 5),
  "Maximum difference <= 0.15 percentage points",
  "Allows for rounding in the original table"
)

check_table <- bind_rows(checks)


# Save manifest, checks and checksums -------------------------------------

write_csv(manifest, file.path(out_dir, "11_supplementary_manifest.csv"))
write_csv(check_table, file.path(report_dir, "11_integrity_checks.csv"))

source_files <- file.path(src_dir[supp$script], supp$in_file)

checksums <- tibble(
  file = source_files,
  md5 = unname(tools::md5sum(source_files))
)
write_csv(checksums, file.path(report_dir, "11_source_checksums.csv"))


# Summary report ----------------------------------------------------------

n_pass <- sum(check_table$status == "PASS")
n_fail <- sum(check_table$status == "FAIL")
n_unchecked <- sum(check_table$status == "NOT CHECKED")

report <- c(
  "SCRIPT 11: SUPPLEMENTARY TABLES",
  paste("Date:", Sys.Date()),
  "",
  "OUTPUT SUMMARY",
  paste("Supplementary files:", nrow(manifest)),
  paste("Integrity checks passed:", n_pass),
  paste("Integrity checks failed:", n_fail),
  paste("Integrity checks not completed:", n_unchecked),
  "",
  "SCIENTIFIC INTERPRETATION",
  "The 926 shared candidates were selected descriptively.",
  "They are not FDR-significant gene-level discoveries.",
  "The permutation analysis did not identify excess",
  "aggregate directional convergence.",
  "GO enrichment is exploratory.",
  "",
  "LIMITATIONS",
  "This script only checks that the outputs are consistent.",
  "It does not independently validate the original GEO",
  "contrasts, variance estimates or ortholog assignments."
)

writeLines(report, file.path(report_dir, "11_supplementary_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(report_dir, "11_session_info.txt"))

cat("\nSUPPLEMENTARY TABLE MANIFEST\n")
print(manifest, n = Inf)

cat("\nINTEGRITY CHECK RESULTS\n")
print(check_table, n = Inf)

cat("\nSCRIPT 11 DONE\n")
cat("Supplementary files:", nrow(manifest), "\n")
cat("Checks passed:", n_pass, "\n")
cat("Checks failed:", n_fail, "\n")
cat("Checks not completed:", n_unchecked, "\n")
cat("Output directory:", out_dir, "\n")

if (n_fail > 0) {
  warning("Some integrity checks failed. Look at the report before using these tables.")
}
