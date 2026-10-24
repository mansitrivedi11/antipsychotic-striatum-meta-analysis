# Script 04: Drug Contrast Extraction — GSE48954
# Project:   Striatum Antipsychotic Transcriptomics Meta-Analysis

library(tidyverse)

tsv_path <- "data/raw/GSE48954_ResultSet553384.tsv"

if (!file.exists(tsv_path)) {
  stop(paste(
    "File not found:", tsv_path,
    "\nDownload from Gemma: GSE48954 → Differential Expression → TSV #553384"
  ))
}

# Count comment rows to skip
con        <- file(tsv_path, "r")
all_lines  <- readLines(con, n = 20)
close(con)

n_skip <- sum(startsWith(all_lines, "#"))

# Read the data
df_raw <- read.table(
  tsv_path,
  header       = TRUE,
  sep          = "\t",
  skip         = n_skip,
  quote        = "",
  fill         = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)

# Verify required columns exist
required_cols <- c(
  "gene_official_symbol",
  "contrast_96278_log2fc", "contrast_96278_tstat",   # haloperidol
  "contrast_96264_log2fc", "contrast_96264_tstat",   # risperidone
  "contrast_96265_log2fc", "contrast_96265_tstat"    # clozapine
)

missing <- setdiff(required_cols, colnames(df_raw))
if (length(missing) > 0) {
  stop(paste("Missing required columns:", paste(missing, collapse=", ")))
}
# Drug contrast mapping
apd_map <- list(
  haloperidol = list(
    fc_col     = "contrast_96278_log2fc",
    t_col      = "contrast_96278_tstat",
    drug_class = "typical"
  ),
  risperidone = list(
    fc_col     = "contrast_96264_log2fc",
    t_col      = "contrast_96264_tstat",
    drug_class = "atypical"
  ),
  clozapine = list(
    fc_col     = "contrast_96265_log2fc",
    t_col      = "contrast_96265_tstat",
    drug_class = "atypical"
  )
)

extract_one_drug <- function(df, drug_name, fc_col, t_col, drug_class) {
  df %>%
    select(
      gene_symbol = gene_official_symbol,
      log2FC      = all_of(fc_col),
      tstat       = all_of(t_col)
    ) %>%
    
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    mutate(
      log2FC = as.numeric(log2FC),
      tstat  = as.numeric(tstat)
    ) %>%
    filter(!is.na(log2FC), !is.na(tstat), tstat != 0) %>%
    mutate(gene_symbol = toupper(gene_symbol)) %>%
    mutate(
      SE = log2FC / tstat,
      SV = SE^2
    ) %>%
    filter(!is.na(SV), !is.infinite(SV), SV > 0) %>%
    # Condense cross-probe artifacts down to unique target representation
    group_by(gene_symbol) %>%
    summarise(
      log2FC   = mean(log2FC, na.rm = TRUE),
      SV       = mean(SV,     na.rm = TRUE),
      n_probes = n(),
      .groups  = "drop"
    ) %>%
    mutate(
      drug       = drug_name,
      drug_class = drug_class,
      dataset    = "GSE48954",
      species    = "mouse",
      subregion  = "striatum",
      pvalue     = NA_real_   
    )

}

# Run extraction loop dynamically across all specified drug parameters
results_list <- lapply(names(apd_map), function(drug_name) {
  params <- apd_map[[drug_name]]
  extract_one_drug(
    df         = df_raw,
    drug_name  = drug_name,
    fc_col     = params$fc_col,
    t_col      = params$t_col,
    drug_class = params$drug_class
  )
})
GSE48954_clean <- bind_rows(results_list)

# Downstream sanity checks and validation assertions
gene_counts <- GSE48954_clean %>%
  group_by(drug) %>%
  summarise(n = n(), .groups = "drop")

if (!all(gene_counts$n == gene_counts$n[1])) {
  warning("Asymmetric gene records detected across target contrasts.")
}

stopifnot(all(GSE48954_clean$SV > 0))

# Export clean artifacts
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write.csv(GSE48954_clean, "data/processed/GSE48954_APD_striatum_clean.csv", row.names = FALSE)

message(paste("Extraction pipeline complete. Target data records saved:", nrow(GSE48954_clean)))
