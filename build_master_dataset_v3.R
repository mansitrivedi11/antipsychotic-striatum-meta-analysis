# Assemble study-level masters (NO cross-species pooling here)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

dir.create('data/processed', recursive=TRUE, showWarnings=FALSE)
dir.create('docs', recursive=TRUE, showWarnings=FALSE)

input_files <- c(
  GSE48954='data/processed/GSE48954_APD_striatum_clean.csv',
  GSE66275='data/processed/GSE66275_reprocessed_striatum.csv',
  GSE67723='data/processed/GSE67723_gemma_striatum.csv',
  GSE67750='data/processed/GSE67750_gemma_striatum.csv',
  GSE180473='data/processed/GSE180473_gemma_striatum.csv'
)
sensitivity_file <- 'data/processed/GSE66275_reprocessed_striatum_QC_sensitivity.csv'
needed_files <- c(input_files, GSE66275_QC_sensitivity=sensitivity_file)
missing_files <- needed_files[!file.exists(needed_files)]
if (length(missing_files)) stop('Missing files:\n',paste(missing_files,collapse='\n'))

expected <- tibble::tribble(
  ~dataset, ~drug,
  'GSE48954','haloperidol', 'GSE48954','risperidone',
  'GSE48954','clozapine', 'GSE66275','haloperidol',
  'GSE66275','risperidone', 'GSE67723','risperidone',
  'GSE67750','haloperidol', 'GSE180473','risperidone',
  'GSE180473','olanzapine'
)
class_map <- c(haloperidol='typical',risperidone='atypical',
               clozapine='atypical',olanzapine='atypical')

study_metadata <- tibble::tribble(
  ~dataset, ~expected_species, ~platform_default, ~duration_default, ~baseline_default,
  'GSE48954','mouse','microarray','acute time course','Tween 80',
  'GSE66275','rat','microarray','21 days','vehicle (verify formulation)',
  'GSE67723','mouse','microarray',NA_character_,'reference substance role',
  'GSE67750','mouse','RNA-seq',NA_character_,'reference substance role',
  'GSE180473','mouse','microarray','28 days','reference substance role'
)

optional_columns <- c('gene_ncbi_id','gene_ensembl_id','probe_id','probe_name',
                      'n_probes','t_stat','fdr','platform','treatment_duration',
                      'baseline','contrast_id','result_set','vehicle_verified',
                      'SE','SV')

load_dataset <- function(path, dataset_label) {
  df <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE,
                 colClasses='character')
  essential <- c('gene_symbol','log2FC','pvalue','drug','drug_class',
                 'dataset','species','subregion')
  absent <- setdiff(essential,names(df))
  if(length(absent)) stop(path,': missing required columns: ',paste(absent,collapse=', '))
  if(!nrow(df)) stop(path,': empty input')
  if(anyNA(df$dataset) || !all(df$dataset==dataset_label))
    stop(path,': dataset label mismatch')
  for(nm in setdiff(optional_columns,names(df))) df[[nm]] <- NA_character_
  for(nm in c('log2FC','pvalue','SV','SE','t_stat','fdr','n_probes')) {
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }
  
  missing_sv <- !is.finite(df$SV) | df$SV<=0
  can_derive <- missing_sv & is.finite(df$log2FC) &
    is.finite(df$t_stat) & df$t_stat!=0
  df$SV[can_derive] <- (df$log2FC[can_derive]/df$t_stat[can_derive])^2
  df$SE[can_derive] <- sqrt(df$SV[can_derive])
  df$variance_source <- ifelse(dataset_label=='GSE66275',
                              'limma moderated t-derived',
                              'Gemma t-derived (verify moderation)')

  invalid <- !is.finite(df$log2FC) | !is.finite(df$pvalue) |
    df$pvalue<0 | df$pvalue>1 | !is.finite(df$SV) | df$SV<=0
  if(any(invalid)) stop(path,': ',sum(invalid),
                        ' invalid effect/variance/p-value rows. Review extraction.')
  if(anyNA(df$gene_symbol) || any(!nzchar(trimws(df$gene_symbol))) ||
     any(grepl('[|;,]|///',df$gene_symbol)))
    stop(path,': missing/ambiguous gene symbols')
  df$gene_symbol_original <- trimws(df$gene_symbol)

  if(anyDuplicated(df[,c('dataset','drug','gene_symbol_original')]))
    stop(path,': duplicated dataset/drug/gene_symbol rows')
  meta <- study_metadata[study_metadata$dataset==dataset_label,,drop=FALSE]
  if(nrow(meta)!=1) stop('Unknown dataset: ',dataset_label)
  if(any(df$species!=meta$expected_species,na.rm=TRUE) || anyNA(df$species))
    stop(path,': species mismatch')
  if(anyNA(df$subregion) || any(df$subregion!='striatum'))
    stop(path,': non-striatum labels')
  if(anyNA(df$drug) || anyNA(df$drug_class) ||
     any(!df$drug %in% names(class_map)) ||
     any(df$drug_class!=unname(class_map[df$drug])))
    stop(path,': unexpected drug or drug class')
  fill_if_blank <- function(x, default) {
    missing <- is.na(x) | !nzchar(trimws(x))
    x[missing] <- default
    x
  }
  df$platform <- fill_if_blank(df$platform,meta$platform_default)
  df$treatment_duration <- fill_if_blank(df$treatment_duration,meta$duration_default)
  df$baseline <- fill_if_blank(df$baseline,meta$baseline_default)

  df$study_id <- dataset_label
  df$ortholog_status <- ifelse(df$species=='mouse',
                               'mouse reference (not cross-species mapped)',
                               'rat; ortholog mapping pending')
  df$harmonized_mouse_gene_id <- NA_character_
  as_tibble(df)
}

build_master <- function(gse66275_path,label) {
  paths <- input_files
  paths[['GSE66275']] <- gse66275_path
  master <- bind_rows(lapply(names(paths),function(ds) load_dataset(paths[[ds]],ds)))
  actual <- distinct(master,dataset,drug)
  missing <- anti_join(expected,actual,by=c('dataset','drug'))
  extra <- anti_join(actual,expected,by=c('dataset','drug'))
  if(nrow(missing) || nrow(extra)) {
    print(missing);print(extra)
    stop(label,': unexpected/missing contrasts')
  }
  counts <- count(master,dataset,drug,name='n_genes')
  if(any(counts$n_genes<5000))
    stop(label,': a contrast has fewer than 5,000 genes')
  message('\n',label,': ',nrow(master),' rows; ',
          n_distinct(master$gene_symbol_original),' distinct symbols (NOT orthologs)')
  print(counts,n=Inf)
  master
}

primary <- build_master(input_files[['GSE66275']],'PRIMARY (GSE66275 all 15)')
sensitivity <- build_master(sensitivity_file,'QC SENSITIVITY (GSE66275 exclusions)')

non_rat_primary <- filter(primary,dataset!='GSE66275')
non_rat_sensitivity <- filter(sensitivity,dataset!='GSE66275')
if(!isTRUE(all.equal(non_rat_primary,non_rat_sensitivity,
                     check.attributes=FALSE)))
  stop('Primary and sensitivity masters differ outside GSE66275')

summary <- primary %>%
  group_by(dataset,drug,species,platform) %>%
  summarise(n_genes=n(),median_log2FC=median(log2FC),
            median_SV=median(SV),nominal_p05=sum(pvalue<0.05),
            .groups='drop')
print(summary,n=Inf)
write.csv(summary,'docs/MASTER_contrast_summary.csv',row.names=FALSE)
write.csv(primary,'data/processed/MASTER_all_datasets.csv',row.names=FALSE)
saveRDS(primary,'data/processed/MASTER_all_datasets.rds')
write.csv(sensitivity,'data/processed/MASTER_all_datasets_GSE66275_QC_sensitivity.csv',
          row.names=FALSE)
saveRDS(sensitivity,'data/processed/MASTER_all_datasets_GSE66275_QC_sensitivity.rds')
writeLines(capture.output(sessionInfo()),'docs/session_info_06.txt')
