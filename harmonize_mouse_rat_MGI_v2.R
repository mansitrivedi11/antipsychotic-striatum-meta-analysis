# MGI mouse-rat ortholog harmonization
# Project: Antipsychotic striatum transcriptomics

pkgs <- c("dplyr", "readr", "tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install missing packages: ", paste(missing, collapse = ", "))
}

library(dplyr)
library(readr)

options(timeout = 600)

# 1. Directories & input checks -------------------------------------------

reference_dir <- "data/reference"
report_dir <- "docs/ortholog_mapping"
processed_dir <- "data/processed"

for (d in c(reference_dir, report_dir, processed_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

primary_file <- file.path(processed_dir, "MASTER_all_datasets.rds")
sensitivity_file <- file.path(processed_dir, "MASTER_all_datasets_GSE66275_QC_sensitivity.rds")

if (!file.exists(primary_file)) {
  stop("Primary master dataset not found: ", primary_file)
}

master <- readRDS(primary_file)
message("Primary master loaded: ", nrow(master), " rows")

req_cols <- c(
  "dataset", "drug", "species", "gene_symbol",
  "gene_ncbi_id", "gene_ensembl_id", "log2FC", "SE", "SV"
)

missing_cols <- setdiff(req_cols, names(master))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

if (any(is.na(master$species)) || !all(master$species %in% c("mouse", "rat"))) {
  stop("Unexpected species values in master dataset.")
}

message("Datasets: ", n_distinct(master$dataset))
message("Drug contrasts: ", nrow(distinct(master, dataset, drug)))

# 2. Helper functions -----------------------------------------------------

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "-", "NA", "N/A", "NULL", "null")] <- NA_character_
  x
}

clean_ensembl <- function(x) {
  sub("\\.[0-9]+$", "", clean_id(x))
}

download_reference <- function(url, dest) {
  if (file.exists(dest) && is.finite(file.info(dest)$size) && file.info(dest)$size > 0) {
    message("Using cached reference: ", dest)
    return(invisible(dest))
  }

  message("Downloading: ", basename(dest))
  tmp <- paste0(dest, ".partial")
  if (file.exists(tmp)) unlink(tmp)

  status <- tryCatch(
    utils::download.file(url, destfile = tmp, method = "libcurl", mode = "wb", quiet = FALSE),
    error = function(e) stop("Download failed: ", conditionMessage(e))
  )

  if (!identical(status, 0L) || !file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Reference download produced an invalid file.")
  }

  if (!file.rename(tmp, dest)) {
    stop("Could not move downloaded reference.")
  }

  invisible(dest)
}

# 3. Download & read MGI reference files ----------------------------------

mgi_base <- "https://jax.org"
homology_url <- paste0(mgi_base, "HOM_AllOrganism.rpt")
annotation_url <- paste0(mgi_base, "HGNC_AllianceHomology.rpt")

homology_file <- file.path(reference_dir, "HOM_AllOrganism.rpt")
annotation_file <- file.path(reference_dir, "HGNC_AllianceHomology.rpt")

download_reference(homology_url, homology_file)
download_reference(annotation_url, annotation_file)

message("Reading MGI homology reference...")
homology_raw <- read_tsv(
  homology_file,
  col_names = FALSE,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE,
  progress = FALSE
)

if (ncol(homology_raw) < 5) {
  stop("Unexpected MGI homology report structure.")
}
message("Homology report rows: ", nrow(homology_raw))

# Read annotation file line-by-line because some MGI rows have a trailing empty 17th tab
message("Reading MGI mouse annotation reference...")
annotation_lines <- readLines(annotation_file, warn = FALSE, encoding = "UTF-8")
annotation_fields <- strsplit(annotation_lines, "\t", fixed = TRUE)
field_counts <- lengths(annotation_fields)

if (!all(field_counts %in% c(16L, 17L))) {
  stop("Unexpected field count in MGI annotation report.")
}

has_extra <- field_counts == 17L
if (any(has_extra)) {
  extra_vals <- vapply(annotation_fields[has_extra], function(x) x[17], character(1))
  if (any(nzchar(trimws(extra_vals)))) {
    stop("Unexpected nonempty 17th field in MGI annotation report.")
  }
}

annotation_mat <- do.call(rbind, lapply(annotation_fields, `[`, 1:16))
annotation_raw <- tibble::as_tibble(annotation_mat, .name_repair = "minimal")
names(annotation_raw) <- paste0("X", seq_len(16))

if (
  annotation_raw$X1[1] != "MGI Accession ID" ||
  annotation_raw$X5[1] != "EntrezGene ID" ||
  annotation_raw$X10[1] != "Ensembl Gene ID"
) {
  stop("MGI annotation header does not match expected structure.")
}

annotation_raw <- annotation_raw[-1, , drop = FALSE]
message("Mouse annotation rows: ", nrow(annotation_raw))

# 4. Extract unambiguous 1:1 mouse-rat homology pairs ---------------------
# HOM_AllOrganism.rpt cols: X1 = homology class, X3 = tax ID, X4 = symbol, X5 = NCBI Gene ID

homology <- homology_raw %>%
  transmute(
    homology_class = clean_id(X1),
    tax_id = clean_id(X3),
    gene_symbol = clean_id(X4),
    ncbi_gene_id = clean_id(X5)
  ) %>%
  dplyr::filter(
    tax_id %in% c("10090", "10116"),
    !is.na(homology_class),
    !is.na(ncbi_gene_id)
  ) %>%
  distinct()

message("Mouse/rat homology records: ", nrow(homology))

eligible_classes <- homology %>%
  group_by(homology_class) %>%
  summarise(
    n_mouse = n_distinct(ncbi_gene_id[tax_id == "10090"]),
    n_rat = n_distinct(ncbi_gene_id[tax_id == "10116"]),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_mouse == 1, n_rat == 1) %>%
  pull(homology_class)

mouse_homology <- homology %>%
  dplyr::filter(tax_id == "10090", homology_class %in% eligible_classes) %>%
  transmute(homology_class, mouse_ncbi_id = ncbi_gene_id, mouse_symbol_mgi = gene_symbol) %>%
  distinct()

rat_homology <- homology %>%
  dplyr::filter(tax_id == "10116", homology_class %in% eligible_classes) %>%
  transmute(homology_class, rat_ncbi_id = ncbi_gene_id, rat_symbol_mgi = gene_symbol) %>%
  distinct()

# Enforce strict 1:1 partner uniqueness across both species
ortholog_pairs <- mouse_homology %>%
  inner_join(rat_homology, by = "homology_class") %>%
  distinct() %>%
  group_by(mouse_ncbi_id) %>%
  dplyr::filter(n_distinct(rat_ncbi_id) == 1) %>%
  ungroup() %>%
  group_by(rat_ncbi_id) %>%
  dplyr::filter(n_distinct(mouse_ncbi_id) == 1) %>%
  ungroup() %>%
  distinct()

message("Unambiguous mouse-rat homology pairs: ", nrow(ortholog_pairs))
if (nrow(ortholog_pairs) < 10000) {
  stop("Unexpectedly low mouse-rat homology coverage.")
}

write_csv(ortholog_pairs, file.path(report_dir, "MGI_v2_mouse_rat_ortholog_pairs.csv"))

# 5. Mouse NCBI-to-Ensembl mappings ---------------------------------------
# MGI columns: X1 = MGI ID, X2 = Symbol, X5 = EntrezGene ID, X10 = Ensembl Gene ID

mouse_reference_raw <- annotation_raw %>%
  transmute(
    mouse_symbol = clean_id(X2),
    mouse_ncbi_id = clean_id(X5),
    mouse_ensembl_id = clean_ensembl(X10)
  ) %>%
  dplyr::filter(
    !is.na(mouse_ncbi_id),
    !is.na(mouse_ensembl_id),
    grepl("^ENSMUSG[0-9]+$", mouse_ensembl_id)
  ) %>%
  distinct()

message("Raw mouse NCBI-Ensembl mappings: ", nrow(mouse_reference_raw))

ambiguous_mouse_ncbi <- mouse_reference_raw %>%
  group_by(mouse_ncbi_id) %>%
  dplyr::filter(n_distinct(mouse_ensembl_id) > 1) %>%
  ungroup()

ambiguous_mouse_ensembl <- mouse_reference_raw %>%
  group_by(mouse_ensembl_id) %>%
  dplyr::filter(n_distinct(mouse_ncbi_id) > 1) %>%
  ungroup()

write_csv(ambiguous_mouse_ncbi, file.path(report_dir, "MGI_v2_ambiguous_mouse_NCBI.csv"))
write_csv(ambiguous_mouse_ensembl, file.path(report_dir, "MGI_v2_ambiguous_mouse_Ensembl.csv"))

mouse_reference <- mouse_reference_raw %>%
  group_by(mouse_ncbi_id) %>%
  dplyr::filter(n_distinct(mouse_ensembl_id) == 1) %>%
  ungroup() %>%
  group_by(mouse_ensembl_id) %>%
  dplyr::filter(n_distinct(mouse_ncbi_id) == 1) %>%
  ungroup() %>%
  distinct(mouse_ncbi_id, mouse_ensembl_id, .keep_all = TRUE)

message("Validated MGI mouse NCBI-Ensembl mappings: ", nrow(mouse_reference))
if (nrow(mouse_reference) < 10000) {
  stop("Unexpectedly low MGI mouse annotation coverage.")
}

write_csv(mouse_reference, file.path(report_dir, "MGI_v2_mouse_NCBI_Ensembl_reference.csv"))

# 6. Validate mouse identifiers in master ---------------------------------

mouse_master_ids <- master %>%
  dplyr::filter(species == "mouse") %>%
  transmute(
    gene_symbol = clean_id(gene_symbol),
    mouse_ncbi_id = clean_id(gene_ncbi_id),
    existing_ensembl_id = clean_ensembl(gene_ensembl_id)
  ) %>%
  distinct()

mouse_validation <- mouse_master_ids %>%
  left_join(
    mouse_reference %>% dplyr::select(mouse_ncbi_id, reference_ensembl_id = mouse_ensembl_id),
    by = "mouse_ncbi_id"
  ) %>%
  mutate(
    annotation_status = case_when(
      !is.na(existing_ensembl_id) & !is.na(reference_ensembl_id) &
        existing_ensembl_id == reference_ensembl_id ~ "reference_confirmed",
      !is.na(existing_ensembl_id) & !is.na(reference_ensembl_id) &
        existing_ensembl_id != reference_ensembl_id ~ "reference_conflict",
      is.na(existing_ensembl_id) & !is.na(reference_ensembl_id) ~ "reference_id_available",
      !is.na(existing_ensembl_id) & is.na(reference_ensembl_id) ~ "not_confirmed_by_reference",
      TRUE ~ "unmapped"
    )
  )

write_csv(mouse_validation, file.path(report_dir, "MGI_v2_mouse_identifier_validation.csv"))

cat("\nMouse identifier validation:\n")
print(count(mouse_validation, annotation_status, name = "n"), n = Inf)
