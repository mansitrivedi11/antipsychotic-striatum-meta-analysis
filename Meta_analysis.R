# Script 07: Random-Effects Meta-Analysis
# Project:   Striatum Antipsychotic Transcriptomics Meta-Analysis
#
# Runs a gene-level random-effects meta-analysis across all 9 drug
# contrasts using the metafor package (Viechtbauer 2010).
# Method: REML (restricted maximum likelihood) — standard in the
# Brain Data Alchemy pipeline (Hagenauer et al. 2024).
#
# Each gene gets its own rma() model fitted on the log2FC values and
# sampling variances (SV) across whichever datasets measured that gene.
# We require a gene to appear in at least 2 contrasts before including
# it — anything with only one observation is just a single study result,
# not a meta-analysis.

library(tidyverse)
library(metafor)

dir.create("results", recursive = TRUE, showWarnings = FALSE)

# LOAD MASTER DATASET
master <- readRDS("data/processed/MASTER_all_datasets.rds")

message(sprintf("Loaded: %d rows | %d genes | %d contrasts",
                nrow(master),
                n_distinct(master$gene_symbol),
                n_distinct(paste(master$dataset, master$drug))))

# FILTER: GENES IN AT LEAST 2 CONTRASTS
gene_counts <- master %>%
  group_by(gene_symbol) %>%
  summarise(n_contrasts = n_distinct(paste(dataset, drug)),
            .groups = "drop")

eligible_genes <- gene_counts %>%
  filter(n_contrasts >= 2) %>%
  pull(gene_symbol)

message(sprintf("Genes in >=2 contrasts (eligible): %d", length(eligible_genes)))
message(sprintf("Genes in exactly 1 contrast (excluded): %d",
                sum(gene_counts$n_contrasts == 1)))

# META-ANALYSIS
# One rma() model per gene. 
message(format(Sys.time()))
run_meta_one_gene <- function(gene) {

  d <- master %>%
    filter(gene_symbol == gene,
           !is.na(log2FC), !is.na(SV),
           !is.infinite(log2FC), !is.infinite(SV),
           SV > 0) %>%
    # One row per unique contrast — avoids double-counting
    distinct(dataset, drug, .keep_all = TRUE)

  if (nrow(d) < 2) return(NULL)

  tryCatch({
    fit <- rma(yi      = d$log2FC,
               vi      = d$SV,
               method  = "REML",
               control = list(maxiter = 1000))

    data.frame(
      gene_symbol  = gene,
      meta_log2FC  = as.numeric(fit$beta),
      meta_SE      = fit$se,
      meta_z       = fit$zval,
      meta_p       = fit$pval,
      I2           = fit$I2,
      tau2         = fit$tau2,
      n_contrasts  = nrow(d),
      datasets     = paste(unique(d$dataset), collapse = ";"),
      drugs        = paste(unique(d$drug),    collapse = ";"),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

# Run across all eligible genes
# Progress printed every 2500 genes so you know it's still running
results_list <- vector("list", length(eligible_genes))

for (i in seq_along(eligible_genes)) {
  results_list[[i]] <- run_meta_one_gene(eligible_genes[i])

  if (i %% 2500 == 0) {
    message(sprintf("  %d / %d genes done (%.0f%%)",
                    i, length(eligible_genes),
                    100 * i / length(eligible_genes)))
  }
}

message(sprintf("Meta-analysis complete: %s", format(Sys.time())))

# COMBINE + FDR CORRECTION
meta_df <- bind_rows(Filter(Negate(is.null), results_list))

meta_df$FDR <- p.adjust(meta_df$meta_p, method = "BH")

meta_df <- meta_df %>%
  mutate(
    direction   = ifelse(meta_log2FC > 0,
                         "Up in APD-treated",
                         "Down in APD-treated"),
    significant = FDR < 0.05 & abs(meta_log2FC) > 0.2
  ) %>%
  arrange(FDR)

message(sprintf("\nTotal genes tested:    %d", nrow(meta_df)))
message(sprintf("Significant (FDR<0.05, |log2FC|>0.2): %d",
                sum(meta_df$significant, na.rm = TRUE)))
message(sprintf("  Upregulated:   %d",
                sum(meta_df$significant & meta_df$meta_log2FC > 0,
                    na.rm = TRUE)))
message(sprintf("  Downregulated: %d",
                sum(meta_df$significant & meta_df$meta_log2FC < 0,
                    na.rm = TRUE)))

message("\nTop 20 significant genes:")
print(meta_df %>%
        filter(significant) %>%
        select(gene_symbol, meta_log2FC, meta_p, FDR,
               I2, n_contrasts) %>%
        head(20))

# CLEAN SIGNIFICANT GENE LIST
#
# Remove probe annotation artefacts - probes that mapped to multiple
# genes get concatenated with | or /// separators and slip through
# the gene symbol filter. We also remove olfactory receptors (OR genes),
# vomeronasal receptors (VMN), and other families unlikely to represent
# genuine striatal antipsychotic targets.
# The 34 genes that remain are the core finding of this paper.

sig_clean <- meta_df %>%
  filter(significant) %>%
  filter(
    !grepl("\\|",              gene_symbol),   # multi-gene probes
    !grepl("^LOC|^RGD",       gene_symbol),   # uncharacterised loci
    !grepl("^OR[0-9]",        gene_symbol),   # olfactory receptors
    !grepl("^VMN",            gene_symbol),   # vomeronasal receptors
    !grepl("^DEFB|^DEFA",     gene_symbol),   # defensins
    !grepl("^IFNA|^IFNB",     gene_symbol),   # interferons
    !grepl("^GSDMC",          gene_symbol)    # gasdermin artefacts
  ) %>%
  arrange(FDR)

message(sprintf("\nAfter removing probe artefacts: %d significant genes",
                nrow(sig_clean)))
message(sprintf("  Upregulated:   %d",
                sum(sig_clean$meta_log2FC > 0)))
message(sprintf("  Downregulated: %d",
                sum(sig_clean$meta_log2FC < 0)))

message("\nGenes present in all 9 contrasts:")
print(sig_clean %>%
        filter(n_contrasts == 9) %>%
        select(gene_symbol, meta_log2FC, FDR, direction))

# SAVE

write.csv(meta_df,
          "results/meta_analysis_full_results.csv",
          row.names = FALSE)

write.csv(sig_clean,
          "results/significant_meta_DEGs_CLEAN.csv",
          row.names = FALSE)

message("\nSaved: results/meta_analysis_full_results.csv")
message("Saved: results/significant_meta_DEGs_CLEAN.csv")
message("\nScript 07 complete -> run 08_figures.R")
