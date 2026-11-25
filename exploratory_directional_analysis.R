# Exploratory striatum analysis

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

input_file <- "data/processed/MASTER_primary_MGI_v2_harmonized.rds"
output_dir <- "results_v2/06_exploratory_analysis"

figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Master dataset not found: ", input_file,
       "\nCheck your working directory.")
}

master <- as_tibble(readRDS(input_file))

required_columns <- c(
  "dataset", "drug", "species", "harmonized_mouse_gene_id",
  "gene_symbol", "log2FC", "pvalue"
)

missing_columns <- setdiff(required_columns, names(master))

if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

save_table <- function(data, filename) {
  write_csv(data, file.path(table_dir, filename))
}

save_figure <- function(plot, filename, width, height) {
  ggsave(
    file.path(figure_dir, filename),
    plot = plot, width = width, height = height, dpi = 300
  )
}

cat("Master observations:", nrow(master), "\n")
cat("Unique genes:", n_distinct(master$harmonized_mouse_gene_id), "\n")


# Select the primary drug comparisons

analysis_data <- master %>%
  mutate(
    drug = tolower(trimws(as.character(drug))),
    dataset = as.character(dataset),
    gene_id = as.character(harmonized_mouse_gene_id),
    gene_symbol = as.character(gene_symbol),
    log2FC = as.numeric(log2FC),
    pvalue = as.numeric(pvalue)
  ) %>%
  filter(
    drug %in% c("haloperidol", "risperidone"),
    !is.na(gene_id),
    gene_id != "",
    is.finite(log2FC)
  )

# Each gene should occur only once per dataset and drug.

duplicates <- analysis_data %>%
  count(dataset, drug, gene_id) %>%
  filter(n > 1)

save_table(duplicates, "06_duplicate_gene_checks.csv")

if (nrow(duplicates) > 0) {
  stop("Duplicate gene-dataset-drug records found. ",
       "Review 06_duplicate_gene_checks.csv.")
}

cat("Selected observations:", nrow(analysis_data), "\n")

print(analysis_data %>% count(dataset, species, drug))


# Study-level summary

study_summary <- analysis_data %>%
  group_by(dataset, species, drug) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    n_positive = sum(log2FC > 0),
    n_negative = sum(log2FC < 0),
    n_nominal_p05 = sum(is.finite(pvalue) & pvalue < 0.05),
    .groups = "drop"
  ) %>%
  arrange(drug, dataset)

save_table(study_summary, "06_study_summary.csv")

print(study_summary, n = Inf)


# Within-drug directional consistency
#
# A gene is considered consistent if it occurs in at least two
# studies and all observed effects have the same non-zero sign.
# This is a descriptive criterion, not a significance test.

gene_consistency <- analysis_data %>%
  group_by(drug, gene_id) %>%
  summarise(
    gene_symbol = first(
      gene_symbol[!is.na(gene_symbol)],
      default = NA_character_
    ),
    n_studies = n_distinct(dataset),
    mean_log2FC = mean(log2FC),
    median_log2FC = median(log2FC),
    min_log2FC = min(log2FC),
    max_log2FC = max(log2FC),
    n_positive = sum(log2FC > 0),
    n_negative = sum(log2FC < 0),
    n_zero = sum(log2FC == 0),
    .groups = "drop"
  ) %>%
  mutate(
    all_same_direction = n_studies >= 2 &
      n_zero == 0 &
      (n_positive == n_studies | n_negative == n_studies),

    direction = case_when(
      n_positive == n_studies ~ "Upregulated",
      n_negative == n_studies ~ "Downregulated",
      TRUE ~ "Mixed"
    ),

    min_abs_effect = pmin(abs(min_log2FC), abs(max_log2FC))
  ) %>%
  arrange(drug, desc(n_studies), gene_id)

save_table(gene_consistency, "06_gene_consistency_all.csv")

multi_study_genes <- gene_consistency %>%
  filter(n_studies >= 2)

consistent_genes <- multi_study_genes %>%
  filter(all_same_direction)

save_table(multi_study_genes, "06_multi_study_genes.csv")
save_table(consistent_genes, "06_directionally_consistent_genes.csv")


# Summarise consistency by drug and study coverage

consistency_summary <- multi_study_genes %>%
  group_by(drug, n_studies) %>%
  summarise(
    n_genes = n(),
    n_consistent = sum(all_same_direction),
    percentage_consistent = 100 * n_consistent / n_genes,
    .groups = "drop"
  ) %>%
  arrange(drug, n_studies)

save_table(consistency_summary, "06_consistency_summary.csv")

cat("\nDirectional consistency by study coverage\n")
print(consistency_summary, n = Inf)


# Compare transcriptional responses between the two drugs.
# Mean log2FC values are unweighted descriptive summaries.

haloperidol <- gene_consistency %>%
  filter(drug == "haloperidol") %>%
  transmute(
    gene_id,
    gene_symbol_haloperidol = gene_symbol,
    haloperidol_mean = mean_log2FC,
    haloperidol_studies = n_studies,
    haloperidol_consistent = all_same_direction
  )

risperidone <- gene_consistency %>%
  filter(drug == "risperidone") %>%
  transmute(
    gene_id,
    gene_symbol_risperidone = gene_symbol,
    risperidone_mean = mean_log2FC,
    risperidone_studies = n_studies,
    risperidone_consistent = all_same_direction
  )

cross_drug <- inner_join(haloperidol, risperidone, by = "gene_id") %>%
  mutate(
    gene_symbol = coalesce(
      gene_symbol_haloperidol,
      gene_symbol_risperidone
    ),

    shared_direction = case_when(
      haloperidol_mean > 0 & risperidone_mean > 0 ~ "Both up",
      haloperidol_mean < 0 & risperidone_mean < 0 ~ "Both down",
      haloperidol_mean == 0 | risperidone_mean == 0 ~ "Zero effect",
      TRUE ~ "Opposite directions"
    )
  )

save_table(cross_drug, "06_cross_drug_all_genes.csv")

# Require at least two contributing studies for each drug.

cross_drug_supported <- cross_drug %>%
  filter(
    haloperidol_studies >= 2,
    risperidone_studies >= 2
  )

save_table(cross_drug_supported, "06_cross_drug_multi_study.csv")

# Shared candidates must agree within and between drugs.
# Rank candidates by their average absolute expression change.

shared_candidates <- cross_drug_supported %>%
  filter(
    haloperidol_consistent,
    risperidone_consistent,
    shared_direction %in% c("Both up", "Both down")
  ) %>%
  mutate(
    combined_abs_effect =
      (abs(haloperidol_mean) + abs(risperidone_mean)) / 2
  ) %>%
  arrange(desc(combined_abs_effect), gene_id)

save_table(shared_candidates, "06_shared_candidate_genes.csv")

cat("\nCross-drug comparison\n")
cat("Genes represented in both drugs:", nrow(cross_drug), "\n")
cat("Genes with at least two studies per drug:",
    nrow(cross_drug_supported), "\n")
cat("Shared directional candidates:", nrow(shared_candidates), "\n")


# Figure 1: Gene coverage across studies

p1 <- ggplot(
  study_summary,
  aes(x = dataset, y = n_genes, fill = drug)
) +
  geom_col(position = "dodge") +
  labs(
    title = "Gene coverage across antipsychotic studies",
    x = "GEO dataset",
    y = "Number of genes",
    fill = "Drug"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_figure(p1, "Figure_1_study_overview.png", 9, 5)


# Figure 2: Directional consistency

p2 <- ggplot(
  consistency_summary,
  aes(
    x = factor(n_studies),
    y = percentage_consistent,
    fill = drug
  )
) +
  geom_col(position = "dodge") +
  facet_wrap(~drug) +
  labs(
    title = "Directional consistency across studies",
    subtitle = "Descriptive agreement; not statistical significance",
    x = "Number of contributing studies",
    y = "Genes with concordant direction (%)"
  ) +
  ylim(0, 100) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

save_figure(p2, "Figure_2_directional_consistency.png", 9, 5)


# Figure 3: Haloperidol versus risperidone

p3 <- ggplot(
  cross_drug_supported,
  aes(x = haloperidol_mean, y = risperidone_mean)
) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Haloperidol and risperidone transcriptional responses",
    subtitle = "Mean log2FC across available independent studies",
    x = "Haloperidol: mean log2FC",
    y = "Risperidone: mean log2FC"
  ) +
  theme_bw(base_size = 12)

save_figure(p3, "Figure_3_cross_drug_comparison.png", 7, 6)


# Figure 4: Heatmap of the top 20 shared candidates
# Ranking reflects effect magnitude rather than statistical significance.

top_candidates <- shared_candidates %>%
  slice_head(n = 20) %>%
  pull(gene_id)

if (length(top_candidates) > 0) {

  heatmap_data <- analysis_data %>%
    filter(gene_id %in% top_candidates) %>%
    mutate(
      study_drug = paste(dataset, drug, sep = " | "),
      gene_label = if_else(
        is.na(gene_symbol) | gene_symbol == "",
        gene_id,
        gene_symbol
      )
    )

  p4 <- ggplot(
    heatmap_data,
    aes(x = study_drug, y = gene_label, fill = log2FC)
  ) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0
    ) +
    labs(
      title = "Exploratory shared candidate genes",
      subtitle = "Log2 fold changes across individual contrasts",
      x = "Dataset and drug",
      y = "Gene",
      fill = "log2FC"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1),
      panel.grid = element_blank()
    )

  save_figure(p4, "Figure_4_candidate_heatmap.png", 11, 8)

} else {
  message("No shared candidates found. Heatmap not generated.")
}


# Pairwise study correlations
#
# Compare effect estimates for genes shared between pairs of studies.
# Correlation describes similarity, not gene-level significance.

study_pairs <- analysis_data %>%
  select(dataset, drug, gene_id, log2FC)

pairwise_results <- list()

for (current_drug in c("haloperidol", "risperidone")) {

  drug_data <- study_pairs %>%
    filter(drug == current_drug)

  datasets <- sort(unique(drug_data$dataset))

  if (length(datasets) < 2) next

  comparisons <- combn(datasets, 2, simplify = FALSE)

  for (pair in comparisons) {

    comparison_data <- drug_data %>%
      filter(dataset %in% pair) %>%
      select(gene_id, dataset, log2FC) %>%
      pivot_wider(
        names_from = dataset,
        values_from = log2FC
      )

    x <- comparison_data[[pair[1]]]
    y <- comparison_data[[pair[2]]]

    valid <- is.finite(x) & is.finite(y)
    n_shared <- sum(valid)

    correlation <- NA_real_

    if (n_shared >= 3 &&
        sd(x[valid]) > 0 &&
        sd(y[valid]) > 0) {

      correlation <- cor(x[valid], y[valid])
    }

    same_direction <- if (n_shared > 0) {
      100 * mean(sign(x[valid]) == sign(y[valid]))
    } else {
      NA_real_
    }

    pairwise_results[[length(pairwise_results) + 1]] <- tibble(
      drug = current_drug,
      dataset_1 = pair[1],
      dataset_2 = pair[2],
      n_shared_genes = n_shared,
      pearson_correlation = correlation,
      percentage_same_direction = same_direction
    )
  }
}

pairwise_summary <- bind_rows(pairwise_results)

save_table(pairwise_summary, "06_pairwise_study_correlations.csv")

cat("\nPairwise study comparisons\n")
print(pairwise_summary, n = Inf)


# Save analysis summary and session information

n_haloperidol <- n_distinct(
  analysis_data$dataset[analysis_data$drug == "haloperidol"]
)

n_risperidone <- n_distinct(
  analysis_data$dataset[analysis_data$drug == "risperidone"]
)

summary_lines <- c(
  "EXPLORATORY STRIATUM ANALYSIS",
  paste("Analysis date:", Sys.Date()),
  "",
  paste("Master observations:", nrow(master)),
  paste("Selected observations:", nrow(analysis_data)),
  paste("Independent datasets:", n_distinct(analysis_data$dataset)),
  "",
  paste("Haloperidol datasets:", n_haloperidol),
  paste("Risperidone datasets:", n_risperidone),
  "",
  paste("Genes represented in both drugs:", nrow(cross_drug)),
  paste(
    "Genes represented in >=2 studies per drug:",
    nrow(cross_drug_supported)
  ),
  paste(
    "Directionally consistent shared candidates:",
    nrow(shared_candidates)
  ),
  "",
  "INTERPRETATION NOTES",
  "",
  "1. This analysis is exploratory.",
  "2. Directional consistency is not statistical significance.",
  "3. Mean log2FC values are unweighted descriptive summaries.",
  "4. Studies differ in species, treatment duration and design.",
  "5. Original contrast and variance provenance require caution.",
  "6. Candidate genes require independent biological validation.",
  "7. No causal or clinical conclusions can be drawn."
)

writeLines(
  summary_lines,
  file.path(output_dir, "06_analysis_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "06_session_info.txt")
)

cat("\nExploratory analysis completed.\n")
cat("Results saved to:", output_dir, "\n")
