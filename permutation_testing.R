# Statistical Assessment of Directional Convergence

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

set.seed(20261116)

n_permutations <- 1000

input_file <- file.path(
  "data", "processed",
  "MASTER_primary_MGI_v2_harmonized.rds"
)

output_dir <- file.path(
  "results_v2", "08_statistical_convergence"
)

table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Primary harmonized master dataset not found.")
}

master <- as_tibble(readRDS(input_file))

required <- c(
  "dataset",
  "drug",
  "harmonized_mouse_gene_id",
  "log2FC"
)

missing <- setdiff(required, names(master))

if (length(missing) > 0) {
  stop(
    "Missing columns: ",
    paste(missing, collapse = ", ")
  )
}


data <- master %>%
  transmute(
    dataset = as.character(dataset),
    drug = tolower(trimws(as.character(drug))),
    gene_id = as.character(harmonized_mouse_gene_id),
    log2FC = as.numeric(log2FC)
  ) %>%
  filter(
    drug %in% c("haloperidol", "risperidone"),
    !is.na(gene_id),
    gene_id != "",
    is.finite(log2FC)
  ) %>%
  mutate(
    contrast = paste(dataset, drug, sep = "__")
  )

duplicates <- data %>%
  count(gene_id, contrast) %>%
  filter(n > 1)

if (nrow(duplicates) > 0) {
  stop(
    "Duplicate gene-contrast observations detected. ",
    "Resolve these before continuing."
  )
}

contrast_info <- data %>%
  distinct(dataset, drug, contrast) %>%
  arrange(dataset, drug)

write_csv(
  contrast_info,
  file.path(table_dir, "08_contrast_inventory.csv")
)

cat("\nCONTRASTS INCLUDED\n")
print(contrast_info, n = Inf)

# 2. Construct the gene-by-contrast direction matrix

wide <- data %>%
  mutate(
    direction = sign(log2FC)
  ) %>%
  select(gene_id, contrast, direction) %>%
  pivot_wider(
    names_from = contrast,
    values_from = direction
  ) %>%
  arrange(gene_id)

gene_ids <- wide$gene_id

direction_matrix <- as.matrix(
  wide[, setdiff(names(wide), "gene_id")]
)

storage.mode(direction_matrix) <- "numeric"

contrast_names <- colnames(direction_matrix)

halo_cols <- which(
  grepl("__haloperidol$", contrast_names)
)

risp_cols <- which(
  grepl("__risperidone$", contrast_names)
)

if (length(halo_cols) != 3 ||
    length(risp_cols) != 4) {
  stop(
    "Expected three haloperidol and four risperidone ",
    "contrasts. Check the input dataset."
  )
}

# 3. Define the convergence statistics

calculate_statistics <- function(mat) {

  halo <- mat[, halo_cols, drop = FALSE]
  risp <- mat[, risp_cols, drop = FALSE]

  halo_n <- rowSums(!is.na(halo))
  risp_n <- rowSums(!is.na(risp))

  halo_up <- rowSums(halo == 1, na.rm = TRUE)
  halo_down <- rowSums(halo == -1, na.rm = TRUE)

  risp_up <- rowSums(risp == 1, na.rm = TRUE)
  risp_down <- rowSums(risp == -1, na.rm = TRUE)

  halo_consistent_up <- (
    halo_n >= 2 &
    halo_up == halo_n
  )

  halo_consistent_down <- (
    halo_n >= 2 &
    halo_down == halo_n
  )

  risp_consistent_up <- (
    risp_n >= 2 &
    risp_up == risp_n
  )

  risp_consistent_down <- (
    risp_n >= 2 &
    risp_down == risp_n
  )

  halo_consistent <- (
    halo_consistent_up |
    halo_consistent_down
  )

  risp_consistent <- (
    risp_consistent_up |
    risp_consistent_down
  )

  shared <- (
    (halo_consistent_up & risp_consistent_up) |
    (halo_consistent_down & risp_consistent_down)
  )

  c(
    haloperidol_consistent = sum(halo_consistent),
    risperidone_consistent = sum(risp_consistent),
    both_drugs_consistent = sum(
      halo_consistent & risp_consistent
    ),
    shared_direction_candidates = sum(shared)
  )
}

observed <- calculate_statistics(direction_matrix)

cat("\nOBSERVED STATISTICS\n")
print(observed)

if (observed["shared_direction_candidates"] != 926) {
  stop(
    "Observed shared candidate count does not match ",
    "Script 06 (926). Do not interpret permutation results ",
    "until the discrepancy is resolved."
  )
}

# 4. Identify contrasts belonging to the same study

study_groups <- split(
  seq_along(contrast_names),
  sub("__.*$", "", contrast_names)
)

# 5. Permute gene directions while preserving coverage
#
# Within each study, genes observed for both drugs have
# their paired directions shuffled together.
#
# Genes observed in only one drug are shuffled separately.
#
# Missing positions are never changed.
#
# This preserves within-study pairing for shared genes,
# but it does not preserve cross-study gene dependence.

permute_study <- function(mat, cols) {

  block <- mat[, cols, drop = FALSE]

  if (ncol(block) == 1) {

    observed_rows <- which(!is.na(block[, 1]))

    if (length(observed_rows) > 1) {
      block[observed_rows, 1] <- sample(
        block[observed_rows, 1],
        length(observed_rows),
        replace = FALSE
      )
    }

  } else if (ncol(block) == 2) {

    both <- which(
      !is.na(block[, 1]) &
      !is.na(block[, 2])
    )

    first_only <- which(
      !is.na(block[, 1]) &
      is.na(block[, 2])
    )

    second_only <- which(
      is.na(block[, 1]) &
      !is.na(block[, 2])
    )

    if (length(both) > 1) {

      shuffled <- sample(
        both,
        length(both),
        replace = FALSE
      )

      block[both, ] <- block[
        shuffled,
        ,
        drop = FALSE
      ]
    }

    if (length(first_only) > 1) {
      block[first_only, 1] <- sample(
        block[first_only, 1],
        length(first_only),
        replace = FALSE
      )
    }

    if (length(second_only) > 1) {
      block[second_only, 2] <- sample(
        block[second_only, 2],
        length(second_only),
        replace = FALSE
      )
    }

  } else {
    stop("Unexpected number of drug contrasts per study.")
  }

  block
}

permute_matrix <- function(mat) {

  result <- mat

  for (cols in study_groups) {
    result[, cols] <- permute_study(mat, cols)
  }

  result
}

# 6. Run the permutation test

statistic_names <- names(observed)

permutation_results <- matrix(
  NA_real_,
  nrow = n_permutations,
  ncol = length(observed),
  dimnames = list(NULL, statistic_names)
)

cat("\nRUNNING", n_permutations, "PERMUTATIONS\n")

for (i in seq_len(n_permutations)) {

  permuted <- permute_matrix(direction_matrix)

  permutation_results[i, ] <- calculate_statistics(
    permuted
  )

  if (i %% 100 == 0) {
    cat("Completed:", i, "/", n_permutations, "\n")
  }
}

permutation_table <- as_tibble(
  permutation_results
) %>%
  mutate(
    permutation = row_number(),
    .before = 1
  )

write_csv(
  permutation_table,
  file.path(table_dir, "08_permutation_results.csv")
)

# 7. Calculate empirical p-values
#
# One-sided test:
# Is observed agreement greater than expected under
# the permutation null?
#
# The +1 correction avoids reporting p = 0.

test_summary <- tibble(
  statistic = statistic_names,
  observed = as.numeric(observed),
  null_mean = colMeans(permutation_results),
  null_sd = apply(permutation_results, 2, sd),
  null_2.5_percent = apply(
    permutation_results,
    2,
    quantile,
    probs = 0.025
  ),
  null_97.5_percent = apply(
    permutation_results,
    2,
    quantile,
    probs = 0.975
  )
) %>%
  mutate(
    empirical_p = vapply(
      seq_along(statistic),
      function(j) {
        (
          1 + sum(
            permutation_results[, j] >= observed[j]
          )
        ) / (n_permutations + 1)
      },
      numeric(1)
    ),
    observed_minus_null = observed - null_mean,
    observed_to_null_ratio = if_else(
      null_mean > 0,
      observed / null_mean,
      NA_real_
    )
  )

write_csv(
  test_summary,
  file.path(table_dir, "08_statistical_test_summary.csv")
)

cat("\nSTATISTICAL TEST SUMMARY\n")
print(test_summary, n = Inf, width = Inf)

# 8. Figure: Observed versus null distribution

plot_data <- permutation_table %>%
  select(
    permutation,
    all_of(statistic_names)
  ) %>%
  pivot_longer(
    cols = -permutation,
    names_to = "statistic",
    values_to = "permuted_count"
  )

observed_plot <- tibble(
  statistic = statistic_names,
  observed = as.numeric(observed)
)

p1 <- ggplot(
  plot_data,
  aes(x = permuted_count)
) +
  geom_histogram(
    bins = 35,
    fill = "#A9C7DD",
    color = "white"
  ) +
  geom_vline(
    data = observed_plot,
    aes(xintercept = observed),
    color = "#B2182B",
    linewidth = 0.9
  ) +
  facet_wrap(
    ~statistic,
    scales = "free",
    ncol = 2
  ) +
  labs(
    title = "Observed directional agreement versus permutation null",
    subtitle = paste(
      n_permutations,
      "study-aware gene-label permutations"
    ),
    x = "Number of genes",
    y = "Permutation frequency"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(
    figure_dir,
    "Figure_8_permutation_null_distributions.png"
  ),
  p1,
  width = 11,
  height = 8,
  dpi = 300
)

# 9. Figure: Shared candidate null distribution

shared_null <- permutation_table %>%
  select(
    permutation,
    shared_direction_candidates
  )

shared_observed <- as.numeric(
  observed["shared_direction_candidates"]
)

shared_p <- test_summary %>%
  filter(
    statistic == "shared_direction_candidates"
  ) %>%
  pull(empirical_p)

p2 <- ggplot(
  shared_null,
  aes(x = shared_direction_candidates)
) +
  geom_histogram(
    bins = 40,
    fill = "#A9C7DD",
    color = "white"
  ) +
  geom_vline(
    xintercept = shared_observed,
    color = "#B2182B",
    linewidth = 1
  ) +
  labs(
    title = "Are shared antipsychotic responses greater than chance?",
    subtitle = paste0(
      "Observed = ",
      shared_observed,
      " | Empirical p = ",
      signif(shared_p, 3)
    ),
    x = "Shared directionally consistent genes",
    y = "Permutation frequency"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(
    figure_dir,
    "Figure_9_shared_convergence_test.png"
  ),
  p2,
  width = 9,
  height = 6,
  dpi = 300
)

# 10. Save an interpretation summary

shared_test <- test_summary %>%
  filter(
    statistic == "shared_direction_candidates"
  )

interpretation <- if (
  shared_test$empirical_p < 0.05 &&
  shared_test$observed > shared_test$null_mean
) {

  paste(
    "Observed shared directional agreement exceeded",
    "the study-aware permutation null at the",
    "nominal 0.05 threshold."
  )

} else {

  paste(
    "The analysis did not provide evidence that",
    "shared directional agreement exceeded",
    "the study-aware permutation null at the",
    "nominal 0.05 threshold."
  )

}

summary_lines <- c(
  "STATISTICAL ASSESSMENT OF CONVERGENCE",
  paste("Date:", Sys.Date()),
  "",
  paste("Permutations:", n_permutations),
  paste("Observed shared candidates:", shared_observed),
  paste(
    "Expected shared candidates under null:",
    round(shared_test$null_mean, 2)
  ),
  paste(
    "Null 95% interval:",
    paste0(
      round(shared_test$null_2.5_percent, 1),
      " to ",
      round(shared_test$null_97.5_percent, 1)
    )
  ),
  paste(
    "Observed / expected ratio:",
    round(shared_test$observed_to_null_ratio, 3)
  ),
  paste(
    "Empirical one-sided p-value:",
    signif(shared_test$empirical_p, 4)
  ),
  "",
  "INTERPRETATION",
  interpretation,
  "",
  "LIMITATIONS",
  "This is an exploratory aggregate directional test.",
  "Gene-specific statistical significance is not established.",
  "Permutation assumes exchangeability within study groups.",
  "Gene dependence and shared biological pathways may affect",
  "the calibration of the empirical null.",
  "The procedure preserves measured-gene coverage, sign",
  "frequencies and paired drug directions within studies.",
  "It does not preserve cross-study gene dependence.",
  "Cross-species and treatment-design differences remain.",
  "The four aggregate tests are exploratory and unadjusted.",
  "A small p-value does not establish a shared drug mechanism."
)

writeLines(
  summary_lines,
  file.path(
    output_dir,
    "08_convergence_summary.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "08_session_info.txt"
  )
)
cat("Results saved to:", output_dir, "\n")
