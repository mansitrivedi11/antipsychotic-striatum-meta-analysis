#  ------------------------------------------------------------------------------------------------
# Script 05: Gemma API Extraction — GSE67750, GSE67723, GSE180473
# Project:   Striatum Antipsychotic Transcriptomics Meta-Analysis
# ---------------------------------------------------------------------------------------------------

library(gemma.R)
library(tidyverse)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

extract_from_gemma <- function(rs_id, dataset, drug,
                                drug_class, species, subregion) {

  message(sprintf("Extracting %s | %s | ResultSet %s", dataset, drug, rs_id))

  raw <- tryCatch(
    get_differential_expression_values(resultSet = rs_id),
    error = function(e) {
      message("  Error: ", e$message)
      return(NULL)
    }
  )

  if (is.null(raw)) return(NULL)

  # Data lives in raw[[1]] — Gemma returns a named list
  de <- raw[[1]]
  message(sprintf("  Retrieved %d rows", nrow(de)))

  # Find the log2fc and tstat columns
  # Pattern: contrast_XXXXXX_log2fc / _tstat / _pvalue
  fc_col <- grep("_log2fc$", colnames(de), value = TRUE)[1]
  t_col  <- grep("_tstat$",  colnames(de), value = TRUE)[1]
  p_col  <- grep("^pvalue$", colnames(de), value = TRUE)[1]

  message(sprintf("  Using: %s | %s | %s", fc_col, t_col, p_col))

  # Build clean data frame, compute SV, collapse probes to genes
  de %>%
    transmute(
      gene_symbol = toupper(GeneSymbol),
      log2FC      = as.numeric(.data[[fc_col]]),
      t_stat      = as.numeric(.data[[t_col]]),
      pvalue      = as.numeric(.data[[p_col]])
    ) %>%
    mutate(
      SE = log2FC / t_stat,
      SV = SE^2
    ) %>%
    filter(
      !is.na(gene_symbol), gene_symbol != "",
      !is.na(log2FC),      !is.infinite(log2FC),
      !is.na(SV),          !is.infinite(SV),
      SV > 0,              t_stat != 0
    ) %>%
    group_by(gene_symbol) %>%
    summarise(
      log2FC   = mean(log2FC, na.rm = TRUE),
      SV       = mean(SV,     na.rm = TRUE),
      pvalue   = mean(pvalue, na.rm = TRUE),
      n_probes = n(),
      .groups  = "drop"
    ) %>%
    mutate(
      drug       = drug,
      drug_class = drug_class,
      dataset    = dataset,
      species    = species,
      subregion  = subregion
    )
}

result_GSE67750 <- extract_from_gemma(
  rs_id      = 530218,
  dataset    = "GSE67750",
  drug       = "haloperidol",
  drug_class = "typical",
  species    = "mouse",
  subregion  = "striatum"
)
Sys.sleep(0.5)

result_GSE67723 <- extract_from_gemma(
  rs_id      = 543220,
  dataset    = "GSE67723",
  drug       = "risperidone",
  drug_class = "atypical",
  species    = "mouse",
  subregion  = "striatum"
)


Sys.sleep(0.5)

message("\nExtracting GSE180473 | ResultSet 506295 | Both drugs")
raw_180473 <- tryCatch(
  get_differential_expression_values(resultSet = 506295),
  error = function(e) { message("Error: ", e$message); NULL }
)

if (!is.null(raw_180473)) {

  de_180473 <- raw_180473[[1]]
  message(sprintf("Retrieved %d rows", nrow(de_180473)))
  message("Columns: ", paste(colnames(de_180473), collapse = ", "))

  extract_180473_contrast <- function(de, fc_col, t_col,
                                       drug, drug_class) {
    de %>%
      transmute(
        gene_symbol = toupper(GeneSymbol),
        log2FC      = as.numeric(.data[[fc_col]]),
        t_stat      = as.numeric(.data[[t_col]]),
        pvalue      = as.numeric(pvalue)
      ) %>%
      mutate(SE = log2FC / t_stat, SV = SE^2) %>%
      filter(
        !is.na(gene_symbol), gene_symbol != "",
        !is.na(log2FC),      !is.infinite(log2FC),
        !is.na(SV),          !is.infinite(SV),
        SV > 0,              t_stat != 0
      ) %>%
      group_by(gene_symbol) %>%
      summarise(
        log2FC   = mean(log2FC, na.rm = TRUE),
        SV       = mean(SV,     na.rm = TRUE),
        pvalue   = mean(pvalue, na.rm = TRUE),
        n_probes = n(),
        .groups  = "drop"
      ) %>%
      mutate(
        drug       = drug,
        drug_class = drug_class,
        dataset    = "GSE180473",
        species    = "mouse",
        subregion  = "striatum"
      )
  }

  result_GSE180473_ris <- extract_180473_contrast(
    de         = de_180473,
    fc_col     = "contrast_181985_log2fc",
    t_col      = "contrast_181985_tstat",
    drug       = "risperidone",
    drug_class = "atypical"
  )

  result_GSE180473_ola <- extract_180473_contrast(
    de         = de_180473,
    fc_col     = "contrast_181986_log2fc",
    t_col      = "contrast_181986_tstat",
    drug       = "olanzapine",
    drug_class = "atypical"
  )

  message(sprintf("  Risperidone: %d genes", nrow(result_GSE180473_ris)))
  message(sprintf("  Olanzapine:  %d genes", nrow(result_GSE180473_ola)))

} else {
  result_GSE180473_ris <- NULL
  result_GSE180473_ola <- NULL
}

write.csv(result_GSE67750,
          "data/processed/GSE67750_gemma_striatum.csv",
          row.names = FALSE)

write.csv(result_GSE67723,
          "data/processed/GSE67723_gemma_striatum.csv",
          row.names = FALSE)

# GSE180473: combine both drug contrasts into one file
if (!is.null(result_GSE180473_ris) && !is.null(result_GSE180473_ola)) {
  result_GSE180473 <- bind_rows(result_GSE180473_ris, result_GSE180473_ola)
  write.csv(result_GSE180473,
            "data/processed/GSE180473_gemma_striatum.csv",
            row.names = FALSE)
  message(sprintf("GSE180473 saved: %d rows (2 drugs)", nrow(result_GSE180473)))
}

message("\n=== SUMMARY ===")
message(sprintf("GSE67750  (haloperidol): %d genes", nrow(result_GSE67750)))
message(sprintf("GSE67723  (risperidone): %d genes", nrow(result_GSE67723)))
message(sprintf("GSE180473 (risperidone): %d genes", nrow(result_GSE180473_ris)))
message(sprintf("GSE180473 (olanzapine):  %d genes", nrow(result_GSE180473_ola)))

message("\nScript 05 complete -> run 06_build_master_dataset.R")
