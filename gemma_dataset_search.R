# -------------------------------------------------------------------------------------------------------------
# Script 01: Gemma Database Search & Dataset Identification
# Project:   Converging Transcriptional Effects of Antipsychotic Drugs
#            in the Rodent Striatum: A Gemma-Based Meta-Analysis
# ---------------------------------------------------------------------------------------------------------------
# Pipeline:  Brain Data Alchemy 2024 (Hagenauer et al.)
#            
# Gemma:     Lim et al. 2021, Database, baab006
#
# NOTE: The initial dataset search was conducted manually via the
#       Gemma web interface (gemma.msl.ubc.ca) in October 2026,
#       using the annotation filters "Striatum" (Organism Part) and
#       Haloperidol, Risperidone, Clozapine, Quetiapine (Treatments).
#       This script reproduces that search programmatically via
#       the gemma.R API and documents the search strategy.
#       Search date: October 2026


library(gemma.R)
library(tidyverse)

# Setup search arrays
drug_terms <- c(
  "haloperidol",    # typical D2 antagonist
  "clozapine",      # atypical — prototype
  "risperidone",    # atypical D2/5-HT2A
  "olanzapine",     # atypical multi-receptor
  "quetiapine",     # atypical
  "antipsychotic"   # catch-all
)

region_terms <- c(
  "striatum",
  "caudate putamen",
  "caudate",
  "nucleus accumbens",
  "neostriatum",
  "basal ganglia"
)

target_taxa <- c("rat", "mouse")

cat(sprintf("Search terms: %d drug × %d region combinations\n",
            length(drug_terms), length(region_terms)))
cat(sprintf("Taxa: %s\n\n", paste(target_taxa, collapse=", ")))

results_all <- list()

# Query Gemma API across combinations

for (drug in drug_terms) {
  for (region in region_terms) {

    query <- paste(drug, region)

    tryCatch({
      res <- get_datasets(
        query = query,
        taxa  = target_taxa,
        limit = 100
      )

      if (!is.null(res) && nrow(res) > 0) {
        res$query_drug   <- drug
        res$query_region <- region
        results_all[[paste(drug, region)]] <- res
        message(paste0("Found: ", query, " -> ", nrow(res), " hits"))
      }

      Sys.sleep(0.5)  

    }, error = function(e) {
      warning(paste0("Failed query: ", query, " - ", e$message))
    })
  }
}

# Deduplicate hits and drop flagged/troubled data

all_df <- bind_rows(results_all) %>%
  distinct(experiment.shortName, .keep_all = TRUE) %>%
  filter(is.na(experiment.troubled) | !experiment.troubled)

cat(sprintf("\nUnique non-troubled datasets found: %d\n\n", nrow(all_df)))

print(all_df %>%
  select(experiment.shortName, taxon.name,
         experiment.sampleCount, experiment.name) %>%
  arrange(desc(experiment.sampleCount)))


# SAVE SEARCH RESULTS


dir.create("data/raw", recursive=TRUE, showWarnings=FALSE)
dir.create("docs",     recursive=TRUE, showWarnings=FALSE)

write.csv(all_df, "data/raw/gemma_search_results.csv", row.names=FALSE)

# ------------------------------------------------------------------------------------------------
# Manual Triage Log (Verified against manual web UI search)
# -------------------------------------------------------------------------------------------------

triage_log <- data.frame(
  dataset = c("GSE37755","GSE67750","GSE66275","GSE227464",
              "GSE153791","GSE143475","GSE180473",
              "GSE67723","GSE48954","GSE45229"),
  species   = c("Mouse","Mouse","Rat","Rat",
                "Mouse","Mouse","Mouse","Mouse","Mouse","Mouse"),
  n_samples = c(275, 28, 45, 12, 16, 96, 96, 16, 156, 20),
  platform  = c("Microarray","RNA-Seq","Microarray","RNA-Seq",
                "RNA-Seq","RNA-Seq","Microarray",
                "Microarray","Microarray","Microarray"),

  decision = c("EXCLUDE","INCLUDE","INCLUDE","EXCLUDE",
               "EXCLUDE","EXCLUDE","INCLUDE","INCLUDE","INCLUDE","EXCLUDE"),

  exclusion_reason = c(
    "Genetic responder vs non-responder comparison — no drug vs vehicle contrast",
    NA,
    NA,
    "n<3 per group; 0 DEGs at FDR<0.05 in Gemma treatment contrast",
    "Striatal astrocyte study — Huntington's disease model, not antipsychotic treatment",
    "Striatal astrocyte study — Huntington's disease model, not antipsychotic treatment",
    NA,
    NA,
    NA,
    "Insufficient n per group after stratification by tissue subregion and dose level"
  ),

  gemma_de_ready = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),

  reprocessing_needed = c(NA, "none", "GEO2R subsetting", NA, NA, NA, "none", "none", "TSV extraction", NA),
  stringsAsFactors = FALSE
)

write.csv(triage_log, "docs/DATASET_TRIAGE_LOG.csv", row.names=FALSE)
