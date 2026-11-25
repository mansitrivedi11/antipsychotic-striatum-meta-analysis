# Biological interpretation and GO enrichment

packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2",
  "AnnotationDbi", "org.Mm.eg.db", "GO.db"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


# Input and output paths

input_file <- "data/processed/MASTER_primary_MGI_v2_harmonized.rds"
out_dir <- "results_v2/09_biological_interpretation"

table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Primary harmonized dataset not found: ", input_file)
}

master <- as_tibble(readRDS(input_file))

required_cols <- c(
  "dataset", "drug", "harmonized_mouse_gene_id",
  "gene_symbol", "log2FC"
)

missing_cols <- setdiff(required_cols, names(master))

if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

save_table <- function(dat, filename) {
  write_csv(dat, file.path(table_dir, filename))
}


# Prepare the primary drug comparisons

dat <- master %>%
  transmute(
    dataset = as.character(dataset),
    drug = tolower(trimws(as.character(drug))),
    gene_id = as.character(harmonized_mouse_gene_id),
    gene_symbol = as.character(gene_symbol),
    log2FC = as.numeric(log2FC)
  ) %>%
  filter(
    drug %in% c("haloperidol", "risperidone"),
    !is.na(gene_id),
    gene_id != "",
    is.finite(log2FC)
  )

duplicates <- dat %>%
  count(dataset, drug, gene_id) %>%
  filter(n > 1)

if (nrow(duplicates) > 0) {
  stop("Duplicate gene-dataset-drug records detected.")
}

contrast_inventory <- dat %>%
  distinct(dataset, drug) %>%
  arrange(dataset, drug)

save_table(contrast_inventory, "09_contrast_inventory.csv")

cat("Selected observations:", nrow(dat), "\n")
cat("Datasets:", n_distinct(dat$dataset), "\n")

print(contrast_inventory, n = Inf)


# Reconstruct the directional candidates from Script 06
#
# Genes must occur in at least two studies per drug.
# All available effects must agree in direction within each drug,
# and both drugs must show the same direction.
# Zero effects do not qualify as directionally consistent.

gene_summary <- dat %>%
  group_by(drug, gene_id) %>%
  summarise(
    n_studies = n_distinct(dataset),
    mean_log2FC = mean(log2FC),
    n_positive = sum(log2FC > 0),
    n_negative = sum(log2FC < 0),
    n_zero = sum(log2FC == 0),
    .groups = "drop"
  ) %>%
  mutate(
    consistent_up = n_studies >= 2 & n_positive == n_studies,
    consistent_down = n_studies >= 2 & n_negative == n_studies
  )

haloperidol <- gene_summary %>%
  filter(drug == "haloperidol") %>%
  transmute(
    gene_id,
    haloperidol_studies = n_studies,
    haloperidol_mean = mean_log2FC,
    haloperidol_up = consistent_up,
    haloperidol_down = consistent_down
  )

risperidone <- gene_summary %>%
  filter(drug == "risperidone") %>%
  transmute(
    gene_id,
    risperidone_studies = n_studies,
    risperidone_mean = mean_log2FC,
    risperidone_up = consistent_up,
    risperidone_down = consistent_down
  )

eligible <- inner_join(haloperidol, risperidone, by = "gene_id") %>%
  filter(
    haloperidol_studies >= 2,
    risperidone_studies >= 2
  ) %>%
  mutate(
    candidate = (haloperidol_up & risperidone_up) |
      (haloperidol_down & risperidone_down),

    shared_direction = case_when(
      haloperidol_up & risperidone_up ~ "Both up",
      haloperidol_down & risperidone_down ~ "Both down",
      TRUE ~ "Not consistent"
    )
  )

universe_ids <- unique(eligible$gene_id)
candidate_ids <- unique(eligible$gene_id[eligible$candidate])

cat("\nEligible genes:", length(universe_ids), "\n")
cat("Shared candidates:", length(candidate_ids), "\n")

# Check against the archived primary analysis.

if (length(universe_ids) != 16102) {
  stop(
    "Eligible universe differs from Script 06. Expected 16102; found ",
    length(universe_ids)
  )
}

if (length(candidate_ids) != 926) {
  stop(
    "Candidate count differs from Script 06. Expected 926; found ",
    length(candidate_ids)
  )
}

save_table(eligible, "09_eligible_genes.csv")


# Retrieve mouse GO biological-process annotations
#
# GOALL includes inherited annotations. Only biological-process
# terms are used for the enrichment analysis.

go_raw <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = universe_ids,
  keytype = "ENSEMBL",
  columns = c("GOALL", "ONTOLOGYALL")
)

go_map <- as_tibble(go_raw) %>%
  transmute(
    gene_id = as.character(ENSEMBL),
    go_id = as.character(GOALL),
    ontology = as.character(ONTOLOGYALL)
  ) %>%
  filter(
    !is.na(gene_id),
    !is.na(go_id),
    ontology == "BP"
  ) %>%
  distinct(gene_id, go_id)

if (nrow(go_map) == 0) {
  stop("No GO biological-process annotations were retrieved.")
}

annotated_universe <- intersect(universe_ids, unique(go_map$gene_id))
annotated_candidates <- intersect(candidate_ids, annotated_universe)

go_map <- go_map %>%
  filter(gene_id %in% annotated_universe)

N <- length(annotated_universe)
n <- length(annotated_candidates)

annotation_summary <- tibble(
  category = c(
    "Eligible genes",
    "Candidate genes",
    "GO-annotated eligible genes",
    "GO-annotated candidate genes"
  ),
  n_genes = c(
    length(universe_ids),
    length(candidate_ids),
    N,
    n
  )
)

print(annotation_summary, n = Inf)

save_table(annotation_summary, "09_annotation_coverage.csv")

if (N == 0 || n == 0 || n >= N) {
  stop("Insufficient annotated genes for enrichment analysis.")
}

save_table(go_map, "09_GO_annotation_mapping.csv")


# GO biological-process enrichment
#
# Test terms containing 10-500 annotated background genes.
# Use a one-sided hypergeometric test and apply BH correction
# across all tested terms.

term_sizes <- go_map %>%
  count(go_id, name = "term_size")

tested_terms <- term_sizes %>%
  filter(term_size >= 10, term_size <= 500)

candidate_hits <- go_map %>%
  filter(gene_id %in% annotated_candidates) %>%
  count(go_id, name = "candidate_hits")

enrichment <- tested_terms %>%
  left_join(candidate_hits, by = "go_id") %>%
  mutate(
    candidate_hits = coalesce(candidate_hits, 0L),
    expected_hits = n * term_size / N,
    fold_enrichment = if_else(
      expected_hits > 0,
      candidate_hits / expected_hits,
      NA_real_
    ),
    pvalue = phyper(
      q = candidate_hits - 1,
      m = term_size,
      n = N - term_size,
      k = n,
      lower.tail = FALSE
    ),
    fdr = p.adjust(pvalue, method = "BH")
  ) %>%
  arrange(pvalue, go_id)

if (nrow(enrichment) == 0) {
  stop("No GO terms passed the term-size filters.")
}


# Add GO term descriptions

go_terms <- AnnotationDbi::select(
  GO.db::GO.db,
  keys = unique(enrichment$go_id),
  keytype = "GOID",
  columns = "TERM"
) %>%
  as_tibble() %>%
  transmute(
    go_id = as.character(GOID),
    term = as.character(TERM)
  ) %>%
  distinct(go_id, .keep_all = TRUE)

enrichment <- enrichment %>%
  left_join(go_terms, by = "go_id") %>%
  relocate(term, .after = go_id)

save_table(enrichment, "09_GO_BP_all_results.csv")


# Report terms with FDR < 0.05 and at least three candidate hits.
# The hit threshold is applied after multiple-testing correction.

significant_terms <- enrichment %>%
  filter(fdr < 0.05, candidate_hits >= 3) %>%
  arrange(fdr, pvalue)

save_table(significant_terms, "09_GO_BP_FDR_significant.csv")

cat("\nGO terms tested:", nrow(enrichment), "\n")
cat("FDR-significant terms:", nrow(significant_terms), "\n")

cat("\nTop 15 GO terms by unadjusted p-value\n")

print(
  enrichment %>% slice_head(n = 15),
  n = 15,
  width = Inf
)

cat("\nFDR-significant GO terms\n")

if (nrow(significant_terms) > 0) {
  print(significant_terms, n = Inf, width = Inf)
} else {
  cat("No GO terms met the reporting criteria.\n")
}


# Identify genes contributing to enriched terms

symbol_lookup <- dat %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_id) %>%
  summarise(
    gene_symbol = first(gene_symbol),
    .groups = "drop"
  )

candidate_go <- go_map %>%
  filter(gene_id %in% annotated_candidates) %>%
  inner_join(
    enrichment %>% select(go_id, term, fdr),
    by = "go_id"
  ) %>%
  left_join(symbol_lookup, by = "gene_id") %>%
  relocate(gene_symbol, .after = gene_id)

save_table(candidate_go, "09_candidate_GO_membership.csv")

significant_gene_membership <- candidate_go %>%
  filter(go_id %in% significant_terms$go_id) %>%
  arrange(fdr, go_id, gene_symbol)

save_table(
  significant_gene_membership,
  "09_significant_GO_contributing_genes.csv"
)

cat("\nGenes contributing to significant GO terms\n")

if (nrow(significant_gene_membership) > 0) {
  print(
    significant_gene_membership,
    n = min(nrow(significant_gene_membership), 100),
    width = Inf
  )
} else {
  cat("No significant GO gene memberships.\n")
}


# Selected genes for biological discussion
#
# These genes are of exploratory interest, not validated biomarkers.

genes_of_interest <- c(
  "Penk", "Neat1", "Ecrg4", "Ddit4l",
  "Atp5mk", "Itgb6", "Frmd7"
)

selected_gene_ids <- dat %>%
  filter(gene_symbol %in% genes_of_interest) %>%
  distinct(gene_id, gene_symbol)

gene_details <- eligible %>%
  inner_join(selected_gene_ids, by = "gene_id") %>%
  select(
    gene_id,
    gene_symbol,
    candidate,
    shared_direction,
    haloperidol_studies,
    haloperidol_mean,
    risperidone_studies,
    risperidone_mean
  ) %>%
  arrange(gene_symbol)

save_table(gene_details, "09_genes_of_interest.csv")

cat("\nSelected genes\n")
print(gene_details, n = Inf, width = Inf)


# Screen GO terms related to neuronal function
# Keyword screening is descriptive, not a separate statistical test.

neuronal_keywords <- c(
  "dopamin", "synap", "neuron", "neuropeptide",
  "neurotransmitter", "axon", "dendrit",
  "learning", "memory", "behavior", "behaviour", "striatal"
)

neuronal_pattern <- paste(neuronal_keywords, collapse = "|")

neuronal_terms <- enrichment %>%
  filter(
    !is.na(term),
    grepl(neuronal_pattern, term, ignore.case = TRUE)
  ) %>%
  arrange(pvalue)

save_table(neuronal_terms, "09_neuronal_related_GO_terms.csv")

cat("\nTop neuronal-related GO terms\n")

if (nrow(neuronal_terms) > 0) {
  print(
    neuronal_terms %>%
      select(
        go_id, term, candidate_hits,
        fold_enrichment, pvalue, fdr
      ) %>%
      slice_head(n = 20),
    n = 20,
    width = Inf
  )
} else {
  cat("No neuronal-related terms identified.\n")
}


# Figure 10: GO biological-process enrichment
#
# Show the 15 lowest-FDR terms with at least three candidate hits,
# including terms that do not reach statistical significance.

plot_terms <- enrichment %>%
  filter(
    candidate_hits >= 3,
    is.finite(fdr),
    !is.na(term)
  ) %>%
  slice_min(fdr, n = 15, with_ties = FALSE) %>%
  mutate(
    term_label = paste0(term, " (", go_id, ")"),
    term_label = factor(term_label, levels = rev(term_label)),
    minus_log10_fdr = -log10(pmax(fdr, .Machine$double.xmin))
  )

if (nrow(plot_terms) > 0) {

  p1 <- ggplot(plot_terms, aes(x = minus_log10_fdr, y = term_label)) +
    geom_col(fill = "#4B7EA8") +
    geom_vline(
      xintercept = -log10(0.05),
      linetype = "dashed",
      color = "#B2182B"
    ) +
    labs(
      title = "Exploratory GO biological-process enrichment",
      subtitle = "BH-adjusted enrichment among eligible genes",
      x = "-log10(BH-adjusted p-value)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(
    file.path(figure_dir, "Figure_10_GO_BP_enrichment.png"),
    plot = p1,
    width = 11,
    height = 7,
    dpi = 300
  )

} else {
  message("No GO terms available for plotting.")
}


# Figure 11: Expression effects for selected genes
# Values are unweighted study-level means, not pooled meta-analysis estimates.

gene_plot <- gene_details %>%
  select(gene_symbol, haloperidol_mean, risperidone_mean) %>%
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

if (nrow(gene_plot) > 0) {

  p2 <- ggplot(
    gene_plot,
    aes(x = gene_symbol, y = mean_log2FC, fill = drug)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(
      title = "Exploratory effect summaries for selected genes",
      subtitle = "Unweighted means across available studies",
      x = NULL,
      y = "Mean log2 fold change",
      fill = "Drug"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    file.path(figure_dir, "Figure_11_selected_gene_effects.png"),
    plot = p2,
    width = 9,
    height = 5,
    dpi = 300
  )

} else {
  message("No selected genes available for plotting.")
}


# Save the analysis summary

report <- c(
  "Biological interpretation and GO enrichment",
  paste("Date:", Sys.Date()),
  "",
  "Dataset",
  paste("Selected observations:", nrow(dat)),
  paste("Datasets:", n_distinct(dat$dataset)),
  "",
  "Gene selection",
  paste("Eligible genes:", length(universe_ids)),
  paste("Exploratory shared candidates:", length(candidate_ids)),
  "",
  "Annotation coverage",
  paste("GO BP-annotated eligible genes:", N),
  paste("GO BP-annotated candidates:", n),
  "",
  "GO enrichment",
  paste("GO terms tested:", nrow(enrichment)),
  paste(
    "Terms with BH-FDR < 0.05 and >=3 hits:",
    nrow(significant_terms)
  ),
  "",
  "Method",
  "Mouse GO biological-process annotations were obtained",
  "from org.Mm.eg.db using harmonized mouse Ensembl IDs.",
  "Inherited GO annotations were included.",
  "The reference universe consisted of eligible genes with",
  "GO biological-process annotations.",
  "Terms containing 10 to 500 background genes were tested.",
  "One-sided hypergeometric tests were used.",
  "BH correction was applied across all tested terms.",
  "",
  "Interpretation",
  "Candidates were selected using descriptive cross-study",
  "directional consistency.",
  "The archived Script 08 results did not identify excess",
  "aggregate directional convergence under the permutation null.",
  "GO enrichment is therefore exploratory.",
  "",
  "Limitations",
  "GO terms overlap and are not independent.",
  "Annotation coverage and gene-selection biases remain.",
  "The hypergeometric model assumes exchangeability of genes",
  "within the annotated eligible universe.",
  "Gene-gene correlation is not modeled.",
  "Candidates are not FDR-significant differential-expression discoveries.",
  "GO enrichment does not validate the underlying convergence hypothesis.",
  "No causal or clinical conclusions are supported."
)

write_lines(
  report,
  file.path(out_dir, "09_biological_interpretation_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "09_session_info.txt")
)


# Check that all expected result tables were exported

expected_files <- c(
  "09_contrast_inventory.csv",
  "09_eligible_genes.csv",
  "09_annotation_coverage.csv",
  "09_GO_annotation_mapping.csv",
  "09_GO_BP_all_results.csv",
  "09_GO_BP_FDR_significant.csv",
  "09_candidate_GO_membership.csv",
  "09_significant_GO_contributing_genes.csv",
  "09_genes_of_interest.csv",
  "09_neuronal_related_GO_terms.csv"
)

missing_outputs <- expected_files[
  !file.exists(file.path(table_dir, expected_files))
]

if (length(missing_outputs) > 0) {
  warning(
    "Missing output files: ",
    paste(missing_outputs, collapse = ", ")
  )
} else {
  cat("\nAll expected result tables were saved.\n")
}

cat("\nBiological interpretation completed.\n")
cat("Results saved to:", out_dir, "\n")
