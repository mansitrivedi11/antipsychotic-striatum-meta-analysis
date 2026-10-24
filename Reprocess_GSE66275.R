# Script 03: GEO2R Reprocessing — GSE66275 Striatum Subset
# Project:   Striatum Antipsychotic Transcriptomics Meta-Analysis


library(GEOquery)
library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# DOWNLOAD GSE66275 FROM GEO
gse <- getGEO("GSE66275", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse[[1]]

cat(sprintf("Downloaded: %d probes × %d samples\n",
            nrow(exprs(gse)), ncol(exprs(gse))))

# Inspect the 4 key columns
pheno <- pData(gse)
cat("=== tissue:ch1 (brain regions) ===\n")
print(table(pheno$`tissue:ch1`))

cat("\n=== treated with:ch1 (drug groups) ===\n")
print(table(pheno$`treated with:ch1`))

cat("\n=== tissue subtype:ch1 ===\n")
print(table(pheno$`tissue subtype:ch1`))

cat("\n=== strain:ch1 ===\n")
print(table(pheno$`strain:ch1`))

cat("\n=== Sample cross-table: tissue x treatment ===\n")
print(table(pheno$`tissue:ch1`, pheno$`treated with:ch1`))
# -------------------------------------------------------------------------------------------------
# Filter Striatum + Check Expression Range
# -------------------------------------------------------------------------------------------------

striatum_mask <- pheno$`tissue subtype:ch1` == "Striatum"

cat(sprintf("\nStriatum samples: %d of %d total\n",
            sum(striatum_mask), nrow(pheno)))

pheno_str <- pheno[striatum_mask, ]
expr_str  <- exprs(gse)[, striatum_mask]

cat(sprintf("Expression matrix: %d probes × %d samples\n",
            nrow(expr_str), ncol(expr_str)))

cat("\nTreatment groups in striatum samples:\n")
print(table(pheno_str$`treated with:ch1`))

# Ensure correct data scaling (log2 transformation fallback)
expr_min <- min(expr_str, na.rm=TRUE)
expr_max <- max(expr_str, na.rm=TRUE)
cat(sprintf("Expression range: %.2f to %.2f\n", expr_min, expr_max))

if (expr_min >= 0 && expr_max <= 20) {
  cat("Range looks like log2-normalised microarray data (4-15)\n")
} else if (expr_max > 100) {
  cat("log2 transform needed\n")
  expr_str <- log2(expr_str + 1)
  cat(sprintf("   After log2 transform: %.2f to %.2f\n",
              min(expr_str), max(expr_str)))
} else {
  cat("Unexpected range, check data manually\n")
}
# Quality Control, Normalization & Filtering

expr_norm <- normalizeBetweenArrays(expr_str, method="quantile")

probe_var <- apply(expr_norm, 1, var)

cat(sprintf("Zero variance probes: %d of %d total\n", sum(probe_var == 0), nrow(expr_norm)))
cat(sprintf("Near-zero variance (<0.01): %d\n", sum(probe_var < 0.01)))

expr_norm_filtered <- expr_norm[probe_var > 0.01, ]
cat(sprintf("Probes retained: %d\n", nrow(expr_norm_filtered)))

treatment_raw <- pheno_str$`treated with:ch1`

treatment_clean <- case_when(
  grepl("Haloperidol", treatment_raw, ignore.case=TRUE) ~ "haloperidol",
  grepl("Risperidone",  treatment_raw, ignore.case=TRUE) ~ "risperidone",
  grepl("vehicle",      treatment_raw, ignore.case=TRUE) ~ "vehicle",
  TRUE ~ treatment_raw
)

cat("\nTreatment groups in striatum:\n")
print(table(treatment_clean))

top_var_idx  <- order(probe_var[probe_var > 0.01],
                       decreasing=TRUE)[1:5000]
expr_for_pca <- expr_norm_filtered[top_var_idx, ]

pca     <- prcomp(t(expr_for_pca), scale.=TRUE)
var_exp <- summary(pca)$importance[2, 1:2] * 100

pca_df <- data.frame(
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  treatment = treatment_clean,
  sample    = colnames(expr_norm_filtered)
)

# Diagnostic plot generation
p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2,
                              color=treatment, label=sample)) +
  geom_point(size=4, alpha=0.9) +
  geom_text_repel(size=2.5, max.overlaps=20) +
  scale_color_manual(values=c("haloperidol"="firebrick",
                               "risperidone"="steelblue",
                               "vehicle"="grey40")) +
  theme_classic(base_size=12) +
  labs(title    = "PCA - GSE66275 Striatum Only (n=15)",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       color = "Treatment")

ggsave("figures/QC_GSE66275_striatum_PCA.pdf",
       p_pca, width=8, height=6)
cat("PCA plot saved: figures/QC_GSE66275_striatum_PCA.pdf\n")

# OUTLIER REMOVAL + REBUILD PCA + LIMMA
# Outliers identified from PCA:
#  GSM1618555 - risperidone, far right on PC1 (~105)
#  GSM1618554 - risperidone, bottom on PC2 (~-65)
# Both are risperidone samples clearly separated from main cluster.

outlier_samples <- c("GSM1618555", "GSM1618554")

keep_mask <- !colnames(expr_norm_filtered) %in% outlier_samples

expr_clean       <- expr_norm_filtered[, keep_mask]
pheno_clean      <- pheno_str[keep_mask, ]
treatment_final  <- treatment_clean[keep_mask]

cat("Samples after outlier removal:\n")
print(table(treatment_final))
cat(sprintf("Matrix: %d probes × %d samples\n",
            nrow(expr_clean), ncol(expr_clean)))

# Re-run PCA on clean matrix to generate structural verification plot
probe_var_clean <- apply(expr_clean, 1, var)
top_idx_clean   <- order(probe_var_clean, decreasing = TRUE)[1:5000]
pca_clean       <- prcomp(t(expr_clean[top_idx_clean, ]), scale. = TRUE)
var_exp_clean   <- summary(pca_clean)$importance[2, 1:2] * 100

pca_df_clean <- data.frame(
  PC1 = pca_clean$x[, 1],
  PC2 = pca_clean$x[, 2],
  treatment = treatment_final,
  sample = colnames(expr_clean)
)

# Plot clean PCA (using cleaner, low-profile minimalist styling)
p_pca_clean <- ggplot(pca_df_clean, aes(x = PC1, y = PC2, color = treatment, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("haloperidol" = "firebrick", "risperidone" = "steelblue", "vehicle" = "grey40")) +
  theme_minimal(base_size = 11) +
  labs(
    title = "PCA Diagnostics: GSE66275 Clean Subset (n=13)",
    x = paste0("PC1 (", round(var_exp_clean[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp_clean[2], 1), "%)")
  )

ggsave("figures/QC_GSE66275_striatum_PCA_clean.pdf", p_pca_clean, width = 7, height = 5)

# Linear Modeling Matrix Setup (limma)

treatment_factor <- factor(treatment_final, levels = c("vehicle", "haloperidol", "risperidone"))
design <- model.matrix(~ treatment_factor)
colnames(design) <- c("Intercept", "haloperidol_vs_vehicle", "risperidone_vs_vehicle")

# Verify design allocation in console
print(colnames(design))
print(table(treatment_factor))

# Run limma
message("Fitting limma linear model...")
fit  <- lmFit(expr_clean, design)
fit2 <- eBayes(fit)

results_hal <- topTable(fit2, coef="haloperidol_vs_vehicle",
                         number=Inf, sort.by="P")
results_ris <- topTable(fit2, coef="risperidone_vs_vehicle",
                         number=Inf, sort.by="P")
# Log summary performance baseline
cat(sprintf("Haloperidol DEGs (FDR<0.05): %d\n",
            sum(results_hal$adj.P.Val < 0.05, na.rm=TRUE)))
cat(sprintf("Risperidone  DEGs (FDR<0.05): %d\n",
            sum(results_ris$adj.P.Val < 0.05, na.rm=TRUE)))
cat(sprintf("Haloperidol nominal (p<0.05): %d\n",
            sum(results_hal$P.Value < 0.05, na.rm=TRUE)))
cat(sprintf("Risperidone  nominal (p<0.05): %d\n",
            sum(results_ris$P.Value < 0.05, na.rm=TRUE)))

# Probe Annotation & Cross-Gene Aggregation
ann       <- fData(gse)
ann_clean <- ann[, c("ID", "Gene symbol")]

# Add probe IDs
results_hal$probe_id <- rownames(results_hal)
results_ris$probe_id <- rownames(results_ris)

hal_ann <- results_hal %>% left_join(ann_clean, by = c("probe_id" = "ID"))
ris_ann <- results_ris %>% left_join(ann_clean, by = c("probe_id" = "ID"))

# Collapse to one value per gene
collapse_to_gene <- function(df, drug_name, drug_class) {

  df %>%
    rename(gene_symbol = `Gene symbol`) %>%
    filter(!is.na(gene_symbol),
           gene_symbol != "",
           gene_symbol != "---") %>%
    # Take first gene if multiple probes mapped (e.g. "Ube2d2///Ube2d3")
    mutate(gene_symbol = gsub("///.*", "", gene_symbol),
           gene_symbol = trimws(gene_symbol),
           # Uppercase for cross-species alignment
           gene_symbol = toupper(gene_symbol)) %>%
    # Sampling variance: SV = (logFC / t)^2
    mutate(SE = logFC / t,
           SV = SE^2) %>%
    filter(!is.na(SV), !is.infinite(SV),
           SV > 0, !is.na(logFC), t != 0) %>%
    # Collapse multiple probes per gene > mean
    group_by(gene_symbol) %>%
    summarise(
      log2FC   = mean(logFC,   na.rm=TRUE),
      SV          = mean(SV,      na.rm=TRUE),
      pvalue   = mean(P.Value, na.rm=TRUE),
      n_probes = n(),
      .groups  = "drop"
    ) %>%
    mutate(
      drug       = drug_name,
      drug_class = drug_class,
      dataset    = "GSE66275",
      species    = "rat",
      subregion  = "striatum",
      note       = "GEO2R reprocessed — striatum subset, 2 outliers removed"
    )
}

gse66275_hal <- collapse_to_gene(hal_ann, "haloperidol", "typical")
gse66275_ris <- collapse_to_gene(ris_ann, "risperidone",  "atypical")

cat(sprintf("\nHaloperidol: %d genes retained\n", nrow(gse66275_hal)))
cat(sprintf("Risperidone:  %d genes retained\n", nrow(gse66275_ris)))

# Sanity checks
cat(sprintf("Haloperidol log2FC range: %.3f to %.3f\n",
            min(gse66275_hal$log2FC), max(gse66275_hal$log2FC)))
cat(sprintf("Risperidone  log2FC range: %.3f to %.3f\n",
            min(gse66275_ris$log2FC), max(gse66275_ris$log2FC)))

cat("\nTop 10 genes (haloperidol, by |log2FC|):\n")
print(gse66275_hal %>%
        arrange(desc(abs(log2FC))) %>%
        select(gene_symbol, log2FC, pvalue) %>%
        head(10))

# Combine and save
GSE66275_clean <- bind_rows(gse66275_hal, gse66275_ris)

dir.create("figures", recursive=TRUE, showWarnings=FALSE)
dir.create("data/processed", recursive=TRUE, showWarnings=FALSE)
saveRDS(GSE66275_clean, "data/processed/GSE66275_reprocessed_striatum.rds")
write.csv(GSE66275_clean,
          "data/processed/GSE66275_reprocessed_striatum.csv",
          row.names=FALSE)

message(paste("Reprocessing pipeline complete. Total rows exported:", nrow(GSE66275_clean)))
