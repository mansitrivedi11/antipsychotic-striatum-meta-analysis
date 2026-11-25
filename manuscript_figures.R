library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Set up directories

results_dir <- "results_v2"

input_dirs <- list(
  s06 = file.path(results_dir, "06_exploratory_analysis", "tables"),
  s07 = file.path(results_dir, "07_QC_sensitivity", "tables"),
  s08 = file.path(results_dir, "08_statistical_convergence", "tables"),
  s09 = file.path(results_dir, "09_biological_interpretation", "tables")
)

output_dir <- file.path(results_dir, "10_final_figures")
png_dir <- file.path(output_dir, "png")
pdf_dir <- file.path(output_dir, "pdf")

dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)


# Functions used throughout the script

read_table <- function(folder, filename) {
  path <- file.path(folder, filename)

  if (!file.exists(path))
    stop("Missing input file: ", path)

  data <- read_csv(path, show_col_types = FALSE)

  if (nrow(data) == 0)
    stop("Empty input file: ", path)

  data
}

# Column names vary slightly between earlier versions of the scripts.

find_column <- function(data, options, label) {
  matches <- intersect(options, names(data))

  if (length(matches) == 0) {
    stop(
      "Cannot find ", label,
      "\nExpected: ", paste(options, collapse = ", "),
      "\nAvailable: ", paste(names(data), collapse = ", ")
    )
  }

  matches[1]
}

check_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))

  if (length(missing) > 0)
    stop(label, " is missing: ", paste(missing, collapse = ", "))
}

as_numeric_checked <- function(x, label) {
  if (is.factor(x)) x <- as.character(x)

  values <- suppressWarnings(as.numeric(x))

  if (any(is.na(values) & !is.na(x)))
    stop("Non-numeric values in ", label)

  values
}

save_figure <- function(plot, name, width = 9, height = 6) {

  ggsave(
    file.path(png_dir, paste0(name, ".png")),
    plot, width = width, height = height,
    units = "in", dpi = 300, bg = "white"
  )

  ggsave(
    file.path(pdf_dir, paste0(name, ".pdf")),
    plot, width = width, height = height,
    units = "in", device = "pdf", bg = "white"
  )

  message("Saved: ", name)
}

theme_manuscript <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(12, 15, 12, 12)
    )
}


# Figure 1. Study integration and analytical workflow

workflow <- data.frame(
  stage = c(
    "Public GEO datasets",
    "Gene identifier harmonization",
    "Directional consistency",
    "QC sensitivity",
    "Permutation analysis",
    "GO biological interpretation"
  ),
  detail = c(
    "5 datasets; 7 haloperidol/risperidone comparisons",
    "Harmonized mouse gene identifiers",
    "926 exploratory shared candidates",
    "887 primary candidates retained",
    "No excess aggregate convergence; p = 0.885",
    "4 overlapping enriched biological-process terms"
  ),
  y = 6:1
)

p1 <- ggplot(workflow, aes(x = 1, y = y)) +
  geom_tile(
    width = 0.92, height = 0.76,
    fill = "#E8F0F7", colour = "#4B7294",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = stage),
    fontface = "bold", size = 4, nudge_y = 0.14
  ) +
  geom_text(
    aes(label = detail),
    size = 3.1, nudge_y = -0.15
  ) +
  scale_x_continuous(limits = c(0.4, 1.6), breaks = NULL) +
  scale_y_continuous(limits = c(0.4, 6.6), breaks = NULL) +
  labs(title = "Study integration and analytical workflow") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  )

save_figure(p1, "Figure_1_study_workflow", 10, 8)


# Figure 2. Within-drug directional consistency

consistency <- read_table(
  input_dirs$s06, "06_consistency_summary.csv"
)

check_columns(
  consistency,
  c("drug", "n_studies", "n_genes",
    "n_consistent", "percentage_consistent"),
  "Script 06 consistency summary"
)

consistency <- consistency %>%
  mutate(
    drug = tools::toTitleCase(as.character(drug)),
    n_studies = as.integer(n_studies),
    n_genes = as_numeric_checked(n_genes, "n_genes"),
    n_consistent = as_numeric_checked(n_consistent, "n_consistent"),
    percent_consistent = as_numeric_checked(
      percentage_consistent, "percentage_consistent"
    )
  )

if (any(!is.finite(consistency$percent_consistent) |
        consistency$percent_consistent < 0 |
        consistency$percent_consistent > 100)) {
  stop("Invalid consistency percentages.")
}

# Check the reported percentages against the original gene counts.

calculated_pct <- with(
  consistency, 100 * n_consistent / n_genes
)

if (any(abs(calculated_pct - consistency$percent_consistent) > 0.15))
  stop("Consistency percentages do not match the gene counts.")

p2 <- ggplot(
  consistency,
  aes(factor(n_studies), percent_consistent, fill = drug)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent_consistent)),
    position = position_dodge(width = 0.75),
    vjust = -0.4, size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, max(consistency$percent_consistent) * 1.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Within-drug directional consistency",
    subtitle = "Observed agreement by number of contributing studies",
    x = "Number of contributing studies",
    y = "Directionally consistent genes (%)",
    fill = "Drug",
    caption = paste(
      "Denominators vary across study-coverage groups.",
      "Directional consistency does not establish",
      "statistically significant convergence."
    )
  ) +
  theme_manuscript()

save_figure(p2, "Figure_2_directional_consistency", 9, 6)


# Figure 3. Mean expression changes across the two drugs

candidates <- read_table(
  input_dirs$s06, "06_shared_candidate_genes.csv"
)

hal_col <- find_column(
  candidates,
  c("haloperidol_mean", "haloperidol_mean_log2FC",
    "mean_haloperidol_log2FC", "haloperidol_mean_log2fc",
    "haloperidol_mean_fc"),
  "haloperidol mean log2FC"
)

ris_col <- find_column(
  candidates,
  c("risperidone_mean", "risperidone_mean_log2FC",
    "mean_risperidone_log2FC", "risperidone_mean_log2fc",
    "risperidone_mean_fc"),
  "risperidone mean log2FC"
)

candidates <- candidates %>%
  mutate(
    haloperidol_plot = as_numeric_checked(.data[[hal_col]], hal_col),
    risperidone_plot = as_numeric_checked(.data[[ris_col]], ris_col)
  )

# Confirm that the table matches the final exploratory analysis.

if (nrow(candidates) != 926)
  stop("Expected 926 shared candidates; found ", nrow(candidates))

if (any(!is.finite(candidates$haloperidol_plot) |
        !is.finite(candidates$risperidone_plot))) {
  stop("Non-finite candidate effect estimates.")
}

p3 <- ggplot(
  candidates,
  aes(haloperidol_plot, risperidone_plot)
) +
  geom_hline(yintercept = 0, colour = "grey65") +
  geom_vline(xintercept = 0, colour = "grey65") +
  geom_point(alpha = 0.5, size = 1.6, colour = "#326B91") +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", colour = "#B2182B"
  ) +
  coord_equal() +
  labs(
    title = "Mean expression changes among shared candidates",
    subtitle = paste(nrow(candidates), "directionally selected genes"),
    x = "Haloperidol mean log2FC",
    y = "Risperidone mean log2FC",
    caption = paste(
      "Unweighted study means, not pooled effect estimates.",
      "Genes were selected for matching expression directions;",
      "the plot does not independently demonstrate convergence."
    )
  ) +
  theme_manuscript()

save_figure(p3, "Figure_3_cross_drug_effects", 8, 7)


# Figure 4. Candidate retention after QC sensitivity analysis

qc <- read_table(
  input_dirs$s07, "07_candidate_comparison.csv"
)

status_col <- find_column(
  qc,
  c("status", "candidate_status", "comparison_status", "membership"),
  "QC candidate membership"
)

status_raw <- tolower(trimws(as.character(qc[[status_col]])))

# Standardize membership labels from earlier output versions.

status_clean <- case_when(
  status_raw %in% c(
    "retained", "both", "shared", "stable",
    "in_both", "both analyses"
  ) ~ "Retained",

  status_raw %in% c(
    "primary only", "primary_only", "primary-only", "primary"
  ) ~ "Primary only",

  status_raw %in% c(
    "sensitivity only", "sensitivity_only",
    "sensitivity-only", "sensitivity"
  ) ~ "Sensitivity only",

  TRUE ~ NA_character_
)

if (anyNA(status_clean)) {
  stop(
    "Unrecognized QC labels: ",
    paste(unique(status_raw[is.na(status_clean)]), collapse = ", ")
  )
}

qc_levels <- c("Retained", "Primary only", "Sensitivity only")

qc_summary <- data.frame(
  status = factor(qc_levels, levels = qc_levels),
  n_genes = as.integer(
    table(factor(status_clean, levels = qc_levels))
  )
)

# These counts should agree with the final Script 07 results.

expected_qc <- c(887L, 39L, 27L)

if (!identical(qc_summary$n_genes, expected_qc)) {
  stop(
    "QC counts differ from the archived results.\n",
    "Observed: ", paste(qc_summary$n_genes, collapse = ", "),
    "\nExpected: ", paste(expected_qc, collapse = ", ")
  )
}

n_retained <- qc_summary$n_genes[1]
n_primary <- sum(qc_summary$n_genes[1:2])
retention_pct <- 100 * n_retained / n_primary

p4 <- ggplot(qc_summary, aes(status, n_genes, fill = status)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_genes), vjust = -0.4, size = 4) +
  scale_y_continuous(
    limits = c(0, max(qc_summary$n_genes) * 1.15)
  ) +
  labs(
    title = "QC sensitivity of exploratory candidates",
    subtitle = sprintf(
      "%d of %d primary candidates retained (%.2f%%)",
      n_retained, n_primary, retention_pct
    ),
    x = NULL,
    y = "Number of genes",
    caption = paste(
      "Primary-only: 39 genes; sensitivity-only: 27 genes.",
      "This assesses one sample-selection decision."
    )
  ) +
  theme_manuscript() +
  theme(legend.position = "none")

save_figure(p4, "Figure_4_QC_candidate_retention", 8, 6)


# Figure 5. Permutation null distribution

permutations <- read_table(
  input_dirs$s08, "08_permutation_results.csv"
)

tests <- read_table(
  input_dirs$s08, "08_statistical_test_summary.csv"
)

permutation_col <- find_column(
  permutations,
  c("shared_direction_candidates", "shared_candidates",
    "n_shared_candidates", "shared_direction_consistent",
    "shared_directional_candidates"),
  "permutation candidate count"
)

stat_col <- find_column(
  tests,
  c("statistic", "test", "metric", "measure"),
  "statistical test identifier"
)

observed_col <- find_column(
  tests,
  c("observed", "observed_count", "observed_value"),
  "observed statistic"
)

null_col <- find_column(
  tests,
  c("null_mean", "permutation_mean", "mean_null"),
  "null mean"
)

p_col <- find_column(
  tests,
  c("empirical_p", "empirical_p_value", "p_value"),
  "empirical p-value"
)

test_names <- as.character(tests[[stat_col]])

shared_index <- which(
  test_names %in% c(
    "shared_direction_candidates",
    "shared_candidates",
    "shared_directional_candidates"
  )
)

if (length(shared_index) != 1) {
  stop(
    "Cannot uniquely identify the shared-candidate test.\n",
    "Available tests: ", paste(test_names, collapse = ", ")
  )
}

observed <- as_numeric_checked(
  tests[[observed_col]][shared_index], observed_col
)

null_mean <- as_numeric_checked(
  tests[[null_col]][shared_index], null_col
)

empirical_p <- as_numeric_checked(
  tests[[p_col]][shared_index], p_col
)

permutation_counts <- as_numeric_checked(
  permutations[[permutation_col]], permutation_col
)

if (any(!is.finite(permutation_counts)))
  stop("Invalid permutation counts.")

if (length(permutation_counts) != 1000)
  stop("Expected 1,000 permutations; found ", length(permutation_counts))

if (observed != 926)
  stop("Unexpected observed candidate count: ", observed)

if (abs(mean(permutation_counts) - null_mean) > 0.1)
  stop("Permutation mean does not match the statistical summary.")

permutation_plot <- data.frame(
  shared_candidates = permutation_counts
)

p5 <- ggplot(permutation_plot, aes(shared_candidates)) +
  geom_histogram(
    bins = 35, fill = "#A9C7DD", colour = "white"
  ) +
  geom_vline(
    xintercept = observed,
    colour = "#B2182B", linewidth = 1
  ) +
  geom_vline(
    xintercept = null_mean,
    colour = "#234F70", linetype = "dashed", linewidth = 0.9
  ) +
  labs(
    title = "Permutation assessment of directional convergence",
    subtitle = sprintf(
      "Observed = %.0f | Null mean = %.2f | One-sided p = %.3f",
      observed, null_mean, empirical_p
    ),
    x = "Shared directional candidates",
    y = "Permutation frequency",
    caption = paste(
      "1,000 study-aware permutations.",
      "Red line: observed count; dashed blue line: null mean.",
      "No excess aggregate convergence detected."
    )
  ) +
  theme_manuscript()

save_figure(p5, "Figure_5_permutation_null", 9, 6)


# Figure 6. GO biological-process enrichment

go <- read_table(
  input_dirs$s09, "09_GO_BP_FDR_significant.csv"
)

go_id_col <- find_column(
  go, c("go_id", "GO_ID", "GOID"), "GO identifier"
)

term_col <- find_column(
  go, c("term", "go_term", "description", "TERM"), "GO term"
)

fold_col <- find_column(
  go,
  c("fold_enrichment", "enrichment", "fold_change", "enrichment_ratio"),
  "GO fold enrichment"
)

fdr_col <- find_column(
  go,
  c("fdr", "FDR", "p_adjusted", "p_adj", "padj"),
  "GO FDR"
)

hits_col <- find_column(
  go,
  c("candidate_hits", "n_candidate_hits", "overlap", "n_hits"),
  "GO candidate hits"
)

go_plot <- data.frame(
  go_id = as.character(go[[go_id_col]]),
  term = as.character(go[[term_col]]),
  fold_enrichment = as_numeric_checked(go[[fold_col]], fold_col),
  fdr = as_numeric_checked(go[[fdr_col]], fdr_col),
  candidate_hits = as_numeric_checked(go[[hits_col]], hits_col)
)

if (nrow(go_plot) != 4)
  stop("Expected four significant GO terms; found ", nrow(go_plot))

if (any(!is.finite(go_plot$fdr) |
        go_plot$fdr >= 0.05 |
        go_plot$candidate_hits < 3)) {
  stop("GO results do not meet the reported selection criteria.")
}

go_plot <- go_plot %>%
  arrange(fdr, desc(fold_enrichment)) %>%
  mutate(
    label = paste0(term, "\n(", go_id, ")"),
    label = factor(label, levels = rev(unique(label)))
  )

p6 <- ggplot(go_plot, aes(fold_enrichment, label)) +
  geom_col(width = 0.65, fill = "#4B7EA8") +
  geom_text(
    aes(
      label = sprintf(
        "FDR = %.3g; n = %d", fdr, as.integer(candidate_hits)
      )
    ),
    hjust = -0.1, size = 3.4
  ) +
  scale_x_continuous(
    limits = c(0, max(go_plot$fold_enrichment) * 1.5)
  ) +
  labs(
    title = "Exploratory GO biological-process enrichment",
    subtitle = paste(
      "Four overlapping terms driven predominantly",
      "by ribosomal protein genes"
    ),
    x = "Fold enrichment",
    y = NULL,
    caption = paste(
      "BH-FDR < 0.05 under a hypergeometric model.",
      "Terms overlap and do not represent",
      "independent pathway discoveries."
    )
  ) +
  theme_manuscript()

save_figure(p6, "Figure_6_GO_enrichment", 10, 6)


# Supplementary Figure S1. Selected genes of interest

genes <- read_table(
  input_dirs$s09, "09_genes_of_interest.csv"
)

gene_col <- find_column(
  genes, c("gene_symbol", "symbol", "gene"), "gene symbol"
)

hal_gene_col <- find_column(
  genes,
  c("haloperidol_mean", "haloperidol_mean_log2FC",
    "mean_haloperidol_log2FC"),
  "selected-gene haloperidol mean"
)

ris_gene_col <- find_column(
  genes,
  c("risperidone_mean", "risperidone_mean_log2FC",
    "mean_risperidone_log2FC"),
  "selected-gene risperidone mean"
)

gene_plot <- data.frame(
  gene_symbol = as.character(genes[[gene_col]]),
  haloperidol_mean = as_numeric_checked(
    genes[[hal_gene_col]], hal_gene_col
  ),
  risperidone_mean = as_numeric_checked(
    genes[[ris_gene_col]], ris_gene_col
  )
)

gene_plot <- gene_plot %>%
  pivot_longer(
    cols = c(haloperidol_mean, risperidone_mean),
    names_to = "drug",
    values_to = "mean_log2FC"
  ) %>%
  mutate(
    drug = recode(
      drug,
      haloperidol_mean = "Haloperidol",
      risperidone_mean = "Risperidone"
    )
  )

p_s1 <- ggplot(gene_plot, aes(gene_symbol, mean_log2FC, fill = drug)) +
  geom_hline(yintercept = 0, colour = "grey45") +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  labs(
    title = "Selected exploratory candidate genes",
    subtitle = "Unweighted mean expression changes",
    x = NULL,
    y = "Mean log2FC",
    fill = "Drug",
    caption = paste(
      "Descriptive estimates across available studies.",
      "Not FDR-significant gene-level discoveries."
    )
  ) +
  theme_manuscript() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(p_s1, "Supplementary_Figure_S1_selected_genes", 9, 6)


# Check that all manuscript figures were saved.

manifest <- data.frame(
  figure = c(
    paste("Figure", 1:6),
    "Supplementary Figure S1"
  ),
  filename = c(
    "Figure_1_study_workflow",
    "Figure_2_directional_consistency",
    "Figure_3_cross_drug_effects",
    "Figure_4_QC_candidate_retention",
    "Figure_5_permutation_null",
    "Figure_6_GO_enrichment",
    "Supplementary_Figure_S1_selected_genes"
  )
)

manifest$png_exists <- file.exists(
  file.path(png_dir, paste0(manifest$filename, ".png"))
)

manifest$pdf_exists <- file.exists(
  file.path(pdf_dir, paste0(manifest$filename, ".pdf"))
)

write_csv(
  manifest, file.path(output_dir, "10_figure_manifest.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "10_session_info.txt")
)

print(manifest, row.names = FALSE)

if (!all(manifest$png_exists & manifest$pdf_exists))
  stop("One or more figure files are missing.")

cat("\nAll manuscript figures saved successfully.\n")
cat("PNG figures:", png_dir, "\n")
cat("PDF figures:", pdf_dir, "\n")
