# Script 09: Cross-Region Comparison with Bhuiyan 2026 Frontal Cortex
# Project:   Striatum Antipsychotic Transcriptomics Meta-Analysis
#
# NOTE: Bhuiyan 2026 (media-2.xlsx) only provides GSEA pathway data, not gene-level 
# logFC. As a proxy for FC regulation, we extract leading-edge genes from sig pathways 
# (FDR < 0.05) and estimate direction by averaging the pathway NES.
# Caveat: Absent genes are classified as "striatum-enriched" rather than "striatum-specific" 
# since lack of pathway enrichment doesn't definitively mean lack of DE.


library(tidyverse)
library(readxl)

if(!dir.exists("results")) dir.create("results", recursive = TRUE)


# Load Bhuiyan 2026 data ---------------------------------------------------------------------


supp_path <- "data/raw/Bhuiyan2026_supplementary.xlsx"

if (!file.exists(supp_path)) {
  stop("Missing Bhuiyan 2026 supp file! Grab media-2.xlsx from bioRxiv:\n",
       "https://www.biorxiv.org/content/10.64898/2026.06.24.734301")
}

gsea_raw <- read_excel(supp_path, sheet = "Exploratory_RodentMacaque_GSEA_")
cat("Loaded GSEA data:", nrow(gsea_raw), "pathways found.\n")


# Parse out leading edge genes to map to their mean NES
gene_nes <- gsea_raw %>%
  filter(padj < 0.05, !is.na(leadingEdge)) %>%
  select(pathway, NES, leadingEdge) %>%
  mutate(gene_symbol = strsplit(as.character(leadingEdge), ",\\s*")) %>%
  unnest(gene_symbol) %>%
  mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
  group_by(gene_symbol) %>%
  summarise(
    avg_NES = mean(as.numeric(NES), na.rm = TRUE),
    n_pathways = n(),
    .groups = "drop"
  ) %>%
  mutate(fc_direction = if_else(avg_NES > 0, "Up in APD-treated", "Down in APD-treated"))

cat("Extracted", nrow(gene_nes), "unique leading-edge genes from significant FC pathways.\n")

write.csv(gene_nes,
          "results/frontal_cortex_leading_edge_genes.csv",
          row.names = FALSE)


# Cross-reference with our Striatal Meta-DEGs --------------------------


sig_striatum <- read.csv("results/significant_meta_DEGs_CLEAN.csv", stringsAsFactors = FALSE)
cat("\nLoaded our data:", nrow(sig_striatum), "striatal meta-DEGs.\n")

comparison <- sig_striatum %>%
  select(gene_symbol, meta_log2FC, meta_SE, FDR,
         n_contrasts, direction) %>%
  left_join(gene_nes %>%
              select(gene_symbol, fc_direction, avg_NES, n_pathways),
            by = "gene_symbol") %>%
  mutate(
    in_frontal_cortex = !is.na(fc_direction),

    region_category = case_when(
      !in_frontal_cortex                       ~ "Striatum-enriched",
      direction == fc_direction                ~ "Convergent",
      direction != fc_direction                ~ "Discordant",
      TRUE                                     ~ "Striatum-enriched"
    )
  )


# SUMMARY


cat("\n=== CROSS-REGION CLASSIFICATION ===\n")
print(table(comparison$region_category))

cat("\nConvergent genes (same direction in both regions):\n")
comparison %>%
  filter(region_category == "Convergent") %>%
  select(gene_symbol, meta_log2FC, direction, fc_direction, FDR) %>%
  print()

cat("\nDiscordant genes (opposite direction between regions):\n")
comparison %>%
  filter(region_category == "Discordant") %>%
  select(gene_symbol, meta_log2FC, direction, fc_direction, FDR) %>%
  print()

cat("\nStriatum-enriched genes (not in FC leading-edge):\n")
comparison %>%
  filter(region_category == "Striatum-enriched") %>%
  select(gene_symbol, meta_log2FC, direction, n_contrasts) %>%
  arrange(FDR) %>%
  head(20) %>% 
  print()


# Check top candidate genes targeted in paper narrative


key_genes <- c("KLF5", "NTS", "LDLRAD3", "ETV5", "ARNTL",
                "S100A10", "GADD45B", "KITLG")

cat("\n=== KEY GENE CLASSIFICATION ===\n")
print(comparison %>%
        filter(gene_symbol %in% key_genes) %>%
        select(gene_symbol, meta_log2FC, direction,
               in_frontal_cortex, fc_direction,
               region_category, n_contrasts) %>%
        arrange(match(gene_symbol, key_genes)))


# SAVE

write.csv(comparison,
          "results/cross_region_comparison.csv",
          row.names = FALSE)

cat("\nSaved: results/cross_region_comparison.csv\n")


# Safely extract counts for manuscript text without risking NA generation
n_conv <- sum(comparison$region_category == "Convergent")
n_disc <- sum(comparison$region_category == "Discordant")
n_enr  <- sum(comparison$region_category == "Striatum-enriched")

cat(sprintf(
  "\n[MANUSCRIPT TEXT]: Out of %d total striatal meta-DEGs evaluated, %d were convergent with the frontal cortex, %d were discordant, and %d were striatum-enriched.\n",
  nrow(comparison), n_conv, n_disc, n_enr
))

cat("\nDone.\n")
