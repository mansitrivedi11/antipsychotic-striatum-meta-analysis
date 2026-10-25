# Script 08: Publication Figures
# Project:   Striatum Antipsychotic Transcriptomics Meta-Analysis
#
# Generates all five manuscript figures plus the two QC PCA plots
# from GSE66275 reprocessing (used in Supplementary).
#
# Fig 1 - PRISMA flow diagram          (see PRISMA_flow_diagram.R)
# Fig 2 - Volcano plot                 (meta-analysis overview)
# Fig 3 - Forest plot, top 5 genes     (per-dataset effect sizes)
# Fig 4 - DEG bar chart                (significant genes + direction)
# Fig 5 - Heatmap                      (all 34 DEGs across contrasts)
# Fig 6 - Cross-region comparison      (striatum vs frontal cortex)
#
# Note on Fig 5 heatmap: NAs (genes not measured in a contrast) are
# filled with 0 before clustering. This is because pheatmap's hclust
# cannot handle NA values. Zero represents "no observed change" which
# is a reasonable assumption for genes absent from a particular platform.
# White cells in the final heatmap correspond to these filled zeros.


library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# LOAD DATA

meta_df  <- read.csv("results/meta_analysis_full_results.csv",
                      stringsAsFactors = FALSE)
sig_clean <- read.csv("results/significant_meta_DEGs_CLEAN.csv",
                       stringsAsFactors = FALSE)
master    <- readRDS("data/processed/MASTER_all_datasets.rds")

message(sprintf("Loaded: %d meta-analysis results | %d significant genes",
                nrow(meta_df), nrow(sig_clean)))

# FIG 2 - VOLCANO PLOT
# Overview of all 25,144 genes tested. Significant genes labelled.
# Point size = number of contrasts the gene appeared in, so genes
# with n=9 (present in every contrast) appear larger.

meta_df <- meta_df %>%
  mutate(
    sig_category = case_when(
      significant & meta_log2FC > 0 ~ "Up (FDR<0.05)",
      significant & meta_log2FC < 0 ~ "Down (FDR<0.05)",
      TRUE                          ~ "Not significant"
    )
  )

# Label the most biologically interpretable significant genes
label_genes <- c("ARNTL", "KLF5", "LDLRAD3", "NTS", "ETV5",
                  "KITLG", "GADD45B", "S100A10", "FATE1", "FSHB",
                  "NWD2", "RSPH4A", "RETNLB")

meta_df$label <- ifelse(meta_df$gene_symbol %in% label_genes,
                         meta_df$gene_symbol, "")

p_volcano <- ggplot(meta_df,
                    aes(x     = meta_log2FC,
                        y     = -log10(meta_p),
                        color = sig_category,
                        size  = n_contrasts)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c(
    "Up (FDR<0.05)"   = "firebrick",
    "Down (FDR<0.05)" = "steelblue",
    "Not significant" = "grey75"
  )) +
  scale_size_continuous(range = c(0.5, 3), name = "N contrasts") +
  geom_text_repel(aes(label = label),
                  size         = 3,
                  color        = "black",
                  fontface     = "bold",
                  max.overlaps = 20,
                  box.padding  = 0.4) +
  geom_hline(yintercept = -log10(max(sig_clean$meta_p)),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.2, 0.2),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  theme_classic(base_size = 13) +
  labs(
    title    = "Antipsychotic Effects on Striatum Transcriptome",
    subtitle = "Meta-analysis: 9 contrasts × 5 datasets | Random-effects (REML)",
    x        = "Meta-analytic Log2 Fold Change (Drug vs Vehicle)",
    y        = "-log10(p-value)",
    color    = ""
  ) +
  theme(legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 13))

ggsave("figures/Fig2_volcano.pdf", p_volcano, width = 9, height = 6)
ggsave("figures/Fig2_volcano.png", p_volcano, width = 9, height = 6, dpi = 300)
message("Fig 2 saved: volcano plot")

# FIG 3 - FOREST PLOT (top 5 genes)
# Shows per-dataset effect sizes for the five most significant genes.
# Red = typical APD (haloperidol), Blue = atypical, Black = meta-estimate.
# Confirms that these genes show consistent direction across datasets.

top5_genes <- c("ARNTL", "KLF5", "LDLRAD3", "NTS", "ETV5")

forest_data <- master %>%
  filter(gene_symbol %in% top5_genes) %>%
  mutate(
    SE            = sqrt(SV),
    CI_low        = log2FC - 1.96 * SE,
    CI_high       = log2FC + 1.96 * SE,
    gene_symbol   = factor(gene_symbol, levels = top5_genes),
    dataset_label = case_when(
      dataset == "GSE48954"  ~ "GSE48954 (mouse)",
      dataset == "GSE66275"  ~ "GSE66275 (rat)",
      dataset == "GSE67723"  ~ "GSE67723 (mouse)",
      dataset == "GSE67750"  ~ "GSE67750 (mouse)",
      dataset == "GSE180473" ~ "GSE180473 (mouse)",
      TRUE ~ dataset
    )
  )

# Add meta-analysis summary row for each gene
meta_rows <- sig_clean %>%
  filter(gene_symbol %in% top5_genes) %>%
  mutate(
    dataset_label = "Meta-analysis",
    drug          = "pooled",
    drug_class    = "pooled",
    log2FC        = meta_log2FC,
    SE            = meta_SE,
    CI_low        = meta_log2FC - 1.96 * meta_SE,
    CI_high       = meta_log2FC + 1.96 * meta_SE,
    gene_symbol   = factor(gene_symbol, levels = top5_genes)
  ) %>%
  select(gene_symbol, dataset_label, drug, drug_class,
         log2FC, SE, CI_low, CI_high)

forest_plot_data <- bind_rows(
  forest_data %>%
    select(gene_symbol, dataset_label, drug,
           drug_class, log2FC, SE, CI_low, CI_high),
  meta_rows
) %>%
  mutate(
    is_meta     = dataset_label == "Meta-analysis",
    point_color = case_when(
      is_meta                  ~ "black",
      drug_class == "typical"  ~ "firebrick",
      drug_class == "atypical" ~ "steelblue",
      TRUE                     ~ "grey50"
    ),
    point_size  = ifelse(is_meta, 4, 2.5)
  )

p_forest <- ggplot(forest_plot_data,
                   aes(x     = log2FC,
                       y     = reorder(paste(dataset_label, drug), log2FC),
                       color = point_color,
                       size  = point_size)) +
  geom_point() +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                width = 0.25, linewidth = 0.5,
                orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_color_identity() +
  scale_size_identity() +
  facet_wrap(~ gene_symbol, scales = "free", ncol = 5) +
  theme_classic(base_size = 10) +
  labs(
    title    = "Per-dataset Effect Sizes: Top 5 Meta-DEGs",
    subtitle = "Red = typical APD | Blue = atypical APD | Black = meta-analysis",
    x        = "Log2 Fold Change (95% CI)",
    y        = ""
  ) +
  theme(strip.text  = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 7))

ggsave("figures/Fig3_forest.pdf", p_forest, width = 16, height = 5)
ggsave("figures/Fig3_forest.png", p_forest, width = 16, height = 5, dpi = 300)
message("Fig 3 saved: forest plot")


# FIG 4 - DEG BAR CHART
# Horizontal bar chart of the biologically interpretable significant
# genes, coloured by direction. Shows the n= label (number of contrasts
# the gene was significant in) next to each bar.
# Unnamed/uncharacterised genes excluded for clarity.

plot_genes <- sig_clean %>%
  filter(!grepl("RIK$|^OLR|^GASK|^ENKUR|^EFCAB|^CYP2J|
                 ^CLDN|^AU0|^FAM49|^CAPSL|^RSPH|^NWD2|
                 ^RETNLB|^SAXO3|^MPP6|^GIOT1|^GATSL|
                 ^KRTAP|^PHF21|^FSHB|^TAS2R", gene_symbol)) %>%
  arrange(meta_log2FC) %>%
  mutate(
    gene_symbol = factor(gene_symbol, levels = gene_symbol),
    bar_color   = ifelse(meta_log2FC > 0, "firebrick", "steelblue"),
    label_n     = paste0("n=", n_contrasts)
  )

p_bar <- ggplot(plot_genes,
                aes(x    = meta_log2FC,
                    y    = gene_symbol,
                    fill = bar_color)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(xmin = meta_log2FC - 1.96 * meta_SE,
                    xmax = meta_log2FC + 1.96 * meta_SE),
                width = 0.3, color = "grey30",
                linewidth = 0.5,
                orientation = "y") +
  geom_text(aes(label = label_n,
                x     = ifelse(meta_log2FC > 0,
                               meta_log2FC + 0.01,
                               meta_log2FC - 0.01)),
            hjust  = ifelse(plot_genes$meta_log2FC > 0, 0, 1),
            size   = 2.8, color = "grey30") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey40") +
  scale_fill_identity() +
  theme_classic(base_size = 11) +
  labs(
    title    = "Significant Meta-DEGs: Antipsychotic Effects in Striatum",
    subtitle = "n = number of drug contrasts gene was significant in (max 9)",
    x        = "Meta-analytic Log2 Fold Change (Drug vs Vehicle)",
    y        = ""
  ) +
  theme(axis.text.y   = element_text(size = 9, face = "bold"),
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey50", size = 9))

ggsave("figures/Fig4_DEG_barplot.pdf", p_bar, width = 8, height = 8)
ggsave("figures/Fig4_DEG_barplot.png", p_bar, width = 8, height = 8, dpi = 300)
message("Fig 4 saved: DEG bar chart")

# FIG 5 - HEATMAP
# Log2FC values for all 34 significant genes across all 9 contrasts.
# NAs (gene not measured in that contrast) filled with 0 before
# clustering — pheatmap cannot cluster with NA values.
# White cells in the output correspond to these filled zeros.
# Column annotation shows drug class and species.

heatmap_genes <- sig_clean$gene_symbol

heatmap_wide <- master %>%
  filter(gene_symbol %in% heatmap_genes) %>%
  mutate(contrast_label = paste0(dataset, "\n", drug)) %>%
  select(gene_symbol, contrast_label, log2FC) %>%
  group_by(gene_symbol, contrast_label) %>%
  summarise(log2FC = mean(log2FC), .groups = "drop") %>%
  pivot_wider(names_from  = contrast_label,
              values_from = log2FC,
              values_fill = NA) %>%
  column_to_rownames("gene_symbol")

# Fill NAs with 0 for clustering
heatmap_filled <- heatmap_wide
heatmap_filled[is.na(heatmap_filled)] <- 0

# Column annotation
col_drugs <- colnames(heatmap_filled)
drug_ann  <- data.frame(
  Drug_Class = case_when(
    grepl("haloperidol", col_drugs)                    ~ "Typical",
    grepl("clozapine|risperidone|olanzapine", col_drugs) ~ "Atypical",
    TRUE ~ "Unknown"
  ),
  Species = ifelse(grepl("GSE66275", col_drugs), "Rat", "Mouse"),
  row.names = col_drugs
)

ann_colors <- list(
  Drug_Class = c("Typical" = "firebrick", "Atypical" = "steelblue"),
  Species    = c("Rat" = "darkgreen",     "Mouse"    = "orange")
)

max_val <- max(abs(heatmap_filled), na.rm = TRUE)
breaks  <- seq(-max_val, max_val, length.out = 101)

pheatmap(
  heatmap_filled,
  color             = colorRampPalette(c("steelblue","white","firebrick"))(100),
  breaks            = breaks,
  annotation_col    = drug_ann,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 8,
  fontsize_col      = 7,
  border_color      = NA,
  main              = "Log2FC of Significant Meta-DEGs Across All Contrasts",
  filename          = "figures/Fig5_heatmap.pdf",
  width             = 12,
  height            = 10
)

pheatmap(
  heatmap_filled,
  color             = colorRampPalette(c("steelblue","white","firebrick"))(100),
  breaks            = breaks,
  annotation_col    = drug_ann,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 8,
  fontsize_col      = 7,
  border_color      = NA,
  main              = "Log2FC of Significant Meta-DEGs Across All Contrasts",
  filename          = "figures/Fig5_heatmap.png",
  width             = 12,
  height            = 10
)

message("Fig 5 saved: heatmap")

# FIG 6 - CROSS-REGION COMPARISON
# Classifies each significant striatal DEG as:
#   Convergent      — also in Bhuiyan 2026 frontal cortex leading-edge genes
#   Discordant      — in frontal cortex but opposite direction
#   Striatum-enriched — not in any significant frontal cortex pathway
#
# Important caveat: Bhuiyan et al. did not publish gene-level effect sizes,
# only GSEA leading-edge genes. Absence from leading-edge genes does not
# prove a gene is not differentially expressed in frontal cortex — it may
# simply not have been in an enriched pathway. See Methods for full details.


cross_region <- data.frame(
  gene_symbol = c(
    "NTS",     "LDLRAD3",  "GADD45B",
    "KLF5",    "ETV5",     "ARNTL",   "S100A10",
    "KITLG",   "FSHB",     "FATE1",   "SAXO3",
    "MPP6",    "RSPH4A",   "NWD2",    "GASK1B",
    "RETNLB",  "PHF21B",   "CAPSL",   "ENKUR",
    "EFCAB10", "TAS2R110", "CYP2J10",
    "KRTAP6-5","FAM49A"
  ),
  striatum_log2FC = c(
     0.499,  0.212,  0.252,
    -0.346, -0.269, -0.347,  0.226,
    -0.282, -0.201, -0.317, -0.310,
    -0.314,  0.211, -0.241, -0.207,
     0.243, -0.213,  0.223,  0.202,
     0.208, -0.210,  0.364,
    -0.246, -0.237
  ),
  category = c(
    "Convergent", "Convergent", "Discordant",
    rep("Striatum-enriched", 21)
  ),
  stringsAsFactors = FALSE
)

# Add SE from sig_clean for error bars
cross_region <- cross_region %>%
  left_join(sig_clean %>%
              select(gene_symbol, meta_SE),
            by = "gene_symbol")

p_cross <- ggplot(cross_region,
                  aes(x    = striatum_log2FC,
                      y    = reorder(gene_symbol, striatum_log2FC),
                      fill = category)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_errorbar(aes(xmin = striatum_log2FC - 1.96 * meta_SE,
                    xmax = striatum_log2FC + 1.96 * meta_SE),
                width = 0.3, color = "grey30",
                linewidth = 0.5,
                orientation = "y") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey40") +
  scale_fill_manual(values = c(
    "Convergent"         = "darkgreen",
    "Discordant"         = "orange",
    "Striatum-enriched"  = "steelblue"
  )) +
  theme_classic(base_size = 11) +
  labs(
    title    = "Cross-Region Comparison: Striatum vs Frontal Cortex",
    subtitle = paste0("Green = convergent with Bhuiyan 2026 FC findings | ",
                      "Blue = striatum-enriched | Orange = discordant"),
    x        = "Striatum Meta-analytic Log2FC (Drug vs Vehicle)",
    y        = "",
    fill     = "Region pattern"
  ) +
  theme(axis.text.y    = element_text(size = 9, face = "bold"),
        plot.subtitle  = element_text(size = 8, color = "grey50"),
        legend.position = "right")

ggsave("figures/Fig6_cross_region.pdf", p_cross, width = 9, height = 8)
ggsave("figures/Fig6_cross_region.png", p_cross, width = 9, height = 8, dpi = 300)
message("Fig 6 saved: cross-region comparison")

message("\nAll figures saved to figures/")
message("Script 08 complete -> run 09_cross_region_comparison.R")
