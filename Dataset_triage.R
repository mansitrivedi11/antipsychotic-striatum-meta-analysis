# --------------------------------------------------------------------------------------------
# Script 02: Dataset Screening and Triage
# Project: Striatum Antipsychotic Transcriptomics Meta-Analysis
# --------------------------------------------------------------------------------------------

library(gemma.R)
library(tidyverse)

# Load raw triage logs from script 01

triage <- read.csv("docs/DATASET_TRIAGE_LOG.csv",
                   stringsAsFactors=FALSE)

included <- triage[triage$decision == "INCLUDE", ]
message(paste("Datasets included after triage:", nrow(included)))

# ------------------------------------------------------------------------------------------------
# Map Gemma Result Set IDs to Include Datasets
# ------------------------------------------------------------------------------------------------

result_set_map <- list(

  GSE48954 = list(
    rs_id        = "553384",
    contrast     = "All drugs vs saline/tween80 — per-drug columns extracted",
    note         = "Pooled contrast — individual drug columns extracted from TSV",
    haloperidol  = "contrast_96278",
    risperidone  = "contrast_96264",
    clozapine    = "contrast_96265",
    reprocessed  = FALSE,
    n_final      = "156 (striatum)"
  ),

  GSE66275 = list(
    rs_id        = "551546",
    contrast     = "haloperidol + risperidone vs baseline (mixed tissue)",
    note         = "Multi-tissue design required GEO2R reprocessing for striatum subset",
    reprocessed  = TRUE,
    outliers_removed = c("GSM1618554", "GSM1618555"),
    outlier_reason   = "PCA outliers — both risperidone group",
    n_final      = "13 (haloperidol=5, risperidone=3, vehicle=5)"
  ),

  GSE67723 = list(
    rs_id        = "543220",
    contrast     = "risperidone vs reference substance role",
    note         = "Used treatment contrast only; strain contrast (543218) excluded",
    reprocessed  = FALSE,
    n_final      = "16"
  ),

  GSE67750 = list(
    rs_id        = "530218",
    contrast     = "haloperidol vs reference substance role",
    note         = "Clean single-drug single-tissue design",
    reprocessed  = FALSE,
    n_final      = "28"
  ),

  GSE180473 = list(
    rs_id        = "506295",
    contrast     = "risperidone + olanzapine vs reference (striatum-specific)",
    note         = "Gemma provided tissue-specific analyses; used organism part: striatum",
    contrast_risperidone = "contrast_181985",
    contrast_olanzapine  = "contrast_181986",
    reprocessed  = FALSE,
    n_final      = "~24 (striatum subset)"
  )
)

# -------------------------------------------------------------------------------------------------
# FINAL INCLUDED DATASET TABLE (for Table 1 in manuscript)
# ------------------------------------------------------------------------------------------------

# Extract key mapping parameters dynamically into a data frame
map_df <- data.frame(
  Dataset = names(result_set_map),
  Result_Set_ID = sapply(result_set_map, function(x) x$rs_id),
  N_subset = sapply(result_set_map, function(x) x$n_final),
  stringsAsFactors = FALSE
)
# Merge back with baseline tracking metadata from script 01
table1 <- included %>%
  rename(Dataset = dataset, Species = species, Platform = platform) %>%
  left_join(map_df, by = "Dataset") %>%
  mutate(
    Brain_Region = ifelse(Dataset == "GSE66275", "Caudate-putamen", "Striatum"),
    Drugs = case_when(
      Dataset == "GSE48954"  ~ "Haloperidol, Risperidone, Clozapine",
      Dataset == "GSE66275"  ~ "Haloperidol, Risperidone",
      Dataset == "GSE67723"  ~ "Risperidone",
      Dataset == "GSE67750"  ~ "Haloperidol",
      Dataset == "GSE180473" ~ "Risperidone, Olanzapine"
    ),
   Drug_Class = case_when(
      Dataset %in% c("GSE48954", "GSE66275") ~ "Typical + Atypical",
      Dataset == "GSE67750"                                   ~ "Typical",
      TRUE                                                                    ~ "Atypical"
    ),
   Duration = case_when(
      Dataset == "GSE48954"  ~ "Acute/subacute (multiple timepoints)",
      Dataset == "GSE66275"  ~ "21 days chronic",
      Dataset == "GSE67723"  ~ "Subacute",
      Dataset == "GSE67750"  ~ "Chronic",
      Dataset == "GSE180473" ~ "28 days chronic"
    ),
   N_contrasts = case_when(
      Dataset == "GSE48954"  ~ 3L,
      Dataset == "GSE66275"  ~ 2L,
      Dataset == "GSE67723"  ~ 1L,
      Dataset == "GSE67750"  ~ 1L,
      Dataset == "GSE180473" ~ 2L
    ),

    Source = ifelse(reprocessing_needed == "none", 
                    "Gemma API direct", 
                    "GEO2R reprocessed")
) %>%
select(Dataset, Species, Platform, Brain_Region, Drugs, Drug_Class, Duration, N_subset, Result_Set_ID, Source)

# Console Verification
print(table1)

# Write output files
write.csv(table1, "docs/INCLUDED_DATASETS_FINAL.csv", row.names = FALSE)
