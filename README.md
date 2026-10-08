# antipsychotic-striatum-meta-analysis

A systematic, genome-wide transcriptomic meta-analysis mapping the molecular signatures of typical and atypical antipsychotic medications within the mammalian striatum. 

## 🧬 Project Overview
This project investigates the long-term transcriptional effects of antipsychotic drug administration within the striatum—a critical subcortical node highly enriched in Dopamine D2 receptors and heavily implicated in both the therapeutic efficacy and extrapyramidal motor side effects of antipsychotic therapy. 

By pulling raw functional genomics datasets (RNA-Seq and microarrays) from public repositories, this study seeks to capture a robust molecular signature that complements existing cortical findings and addresses the classical dopaminergic pathways fundamental to the dopamine hypothesis of psychosis.

## 🛠️ Pipeline Architecture
This repository implements and adapts the standardized, open-science processing workflow developed by the **Brain Data Alchemy Project**. 

The sequential execution flow is structured as follows:
1. **Data Discovery (`BrainDataAlchemy.R`):** Mining metadata and applying curated inclusion/exclusion filters to target striatal tissue samples.
2. **Automated Ingestion :** Interface with the Gemma database API to download standardized differential expression statistics.
3. **Data Normalization :** Standalone quality control fallback utilizing the Bioconductor `limma` pipeline for manual array curation.
4. **Matrix Alignment :** Parsing experimental contrasts, handling ortholog mappings, and isolating target drug dosages.
5. **Statistical Modeling :** Execution of a random-effects meta-analysis model via the `metafor` package to identify significant Differentially Expressed Genes (DEGs).
6. **Functional Annotation & Visualization:** Functional enrichment analyses using `fgsea` alongside publication-ready data plots (Volcano, Hierarchical Clustered Heatmaps).

## ⚖️ Acknowledgments & License
The core analytical logic and infrastructure used in this repository are heavily inspired by the public framework created for the study: *“A Meta-Analysis of the Converging Effects of Different Classes of Antipsychotics on the Frontal Cortex Transcriptome...”* ([raeedbhuiyan/Antipsychotic_Meta-Analysis](https://github.com)).

This project is open-source and distributed under the **MIT License**.
# antipsychotic-striatum-meta-analysis
