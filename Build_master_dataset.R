# Master Dataset Assembly

library(tidyverse)

gse48954  <- read.csv("data/processed/GSE48954_APD_striatum_clean.csv",
                       stringsAsFactors = FALSE)

gse66275  <- read.csv("data/processed/GSE66275_reprocessed_striatum.csv",
                       stringsAsFactors = FALSE)

gse67723  <- read.csv("data/processed/GSE67723_gemma_striatum.csv",
                       stringsAsFactors = FALSE)

gse67750  <- read.csv("data/processed/GSE67750_gemma_striatum.csv",
                       stringsAsFactors = FALSE)

gse180473 <- read.csv("data/processed/GSE180473_gemma_striatum.csv",
                       stringsAsFactors = FALSE)


stopifnot("GSE180473 missing olanzapine — re-run Script 05" =
            "olanzapine" %in% gse180473$drug)
stopifnot("GSE180473 missing risperidone — re-run Script 05" =
            "risperidone" %in% gse180473$drug)


message(sprintf("GSE48954:  %d rows | %d genes per drug",
                nrow(gse48954),  nrow(gse48954) / 3))
message(sprintf("GSE66275:  %d rows | %d genes per drug",
                nrow(gse66275),  nrow(gse66275) / 2))
message(sprintf("GSE67723:  %d rows", nrow(gse67723)))
message(sprintf("GSE67750:  %d rows", nrow(gse67750)))
message(sprintf("GSE180473: %d rows | drugs: %s",
                nrow(gse180473),
                paste(unique(gse180473$drug), collapse = ", ")))


core_cols <- c("gene_symbol", "log2FC", "SV", "pvalue",
                "n_probes", "drug", "drug_class",
                "dataset", "species", "subregion")

# COMBINE

master <- bind_rows(
  gse48954  [, intersect(core_cols, colnames(gse48954))],
  gse66275  [, intersect(core_cols, colnames(gse66275))],
  gse67723  [, intersect(core_cols, colnames(gse67723))],
  gse67750  [, intersect(core_cols, colnames(gse67750))],
  gse180473 [, intersect(core_cols, colnames(gse180473))]
)

# VALIDATE

message("\n=== MASTER DATASET SUMMARY ===")
message(sprintf("Total rows:      %d", nrow(master)))
message(sprintf("Unique genes:    %d", n_distinct(master$gene_symbol)))
message(sprintf("Unique datasets: %d", n_distinct(master$dataset)))
message(sprintf("Drug contrasts:  %d",
                n_distinct(paste(master$dataset, master$drug))))

message("\nRows per drug contrast:")
print(master %>%
        count(dataset, drug, drug_class, species) %>%
        arrange(dataset, drug))


n_na_gene  <- sum(is.na(master$gene_symbol))
n_na_log2fc <- sum(is.na(master$log2FC))
n_na_sv    <- sum(is.na(master$SV))
n_bad_sv   <- sum(master$SV <= 0, na.rm = TRUE)

message(sprintf("\nNA in gene_symbol: %d", n_na_gene))
message(sprintf("NA in log2FC:      %d", n_na_log2fc))
message(sprintf("NA in SV:          %d", n_na_sv))
message(sprintf("SV <= 0:           %d", n_bad_sv))

if (any(c(n_na_gene, n_na_log2fc, n_na_sv, n_bad_sv) > 0)) {
  warning("Data quality issues found - check above before proceeding")
} else {
  message("All key columns clean ")
}


n_contrasts <- n_distinct(paste(master$dataset, master$drug))
stopifnot("Expected 9 drug contrasts" = n_contrasts == 9)
message("9 drug contrasts confirmed ")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write.csv(master,
          "data/processed/MASTER_all_datasets.csv",
          row.names = FALSE)

saveRDS(master,
        "data/processed/MASTER_all_datasets.rds")

message("\nSaved: data/processed/MASTER_all_datasets.csv")
message("Saved: data/processed/MASTER_all_datasets.rds")
