# QC Sensitivity Analysis
# Comparing the primary harmonized dataset with the alternative
# GSE66275 QC sensitivity dataset.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

primary_file <- file.path( "data", "processed","MASTER_primary_MGI_v2_harmonized.rds")
sensitivity_file <- file.path("data", "processed", "MASTER_QC_sensitivity_MGI_v2_harmonized.rds")

output_dir <- file.path("results_v2", "07_QC_sensitivity")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(primary_file)) {
    stop("Primary master dataset not found.")
}

if (!file.exists(sensitivity_file)) {
    stop("QC sensitivity master dataset not found.")
}

primary <- as_tibble(readRDS(primary_file))
sensitivity <- as_tibble(readRDS(sensitivity_file))

cat("\nPRIMARY MASTER:", nrow(primary), "rows\n")
cat("SENSITIVITY MASTER:", nrow(sensitivity), "rows\n")

prepare_data <- function(data) {
    required <- c("dataset", "drug", "species", "harmonized_mouse_gene_id", "gene_symbol",
 "log2FC"
    )
    missing <- setdiff(required, names(data))
    if (length(missing) > 0) {
        stop(
            "Missing columns: ",
            paste(missing, collapse = ", ")
        )
    }

    data <- data %>%
        mutate(
            dataset = as.character(dataset),
            drug = tolower(trimws(as.character(drug))),
            gene_id = as.character(harmonized_mouse_gene_id),
            gene_symbol = as.character(gene_symbol),
            log2FC = as.numeric(log2FC)
        ) %>%
        filter(
            drug %in% c("haloperidol", "risperidone"),
            !is.na(gene_id),
            gene_id != "",
            is.finite(log2FC)
        )

    duplicates <- data %>%
        count(dataset, drug, gene_id) %>%
        filter(n > 1)

    if (nrow(duplicates) > 0) {
        stop(
            "Duplicate gene-dataset-drug records detected."
        )
    }

    data
}

primary_data <- prepare_data(primary)
sensitivity_data <- prepare_data(sensitivity)

comparison <- full_join(
    primary_data %>%
        select(dataset, drug, gene_id, log2FC) %>%
        rename(primary_log2FC = log2FC),

    sensitivity_data %>%
        select(dataset, drug, gene_id, log2FC) %>%
        rename(sensitivity_log2FC = log2FC),

    by = c("dataset", "drug", "gene_id")
) %>%
    mutate(
        difference = sensitivity_log2FC - primary_log2FC,
        direction_changed = case_when(
            is.na(primary_log2FC) |
                is.na(sensitivity_log2FC) ~ NA,
            TRUE ~ sign(primary_log2FC) !=
                sign(sensitivity_log2FC)
        )
    )

write_csv(
    comparison,
    file.path(table_dir, "07_all_effect_comparisons.csv")
)

dataset_changes <- comparison %>%
    group_by(dataset, drug) %>%
    summarise(
        n_primary = sum(!is.na(primary_log2FC)),
        n_sensitivity = sum(!is.na(sensitivity_log2FC)),
        n_shared = sum(
            !is.na(primary_log2FC) &
                !is.na(sensitivity_log2FC)
        ),
        n_effects_changed = sum(
            !is.na(difference) &
                abs(difference) > 1e-10
        ),
        n_directions_changed = sum(
            direction_changed,
            na.rm = TRUE
        ),
        median_absolute_difference = if (
            any(is.finite(difference))
        ) {
            median(abs(difference), na.rm = TRUE)
        } else {
            NA_real_
        },
        .groups = "drop"
    )

write_csv(
    dataset_changes,
    file.path(table_dir, "07_dataset_changes.csv")
)

cat("\nDATASET CHANGES\n")
print(dataset_changes, n = Inf)

calculate_consistency <- function(data) {

    data %>%
        group_by(drug, gene_id) %>%
        summarise(
            gene_symbol = first(
                gene_symbol[
                    !is.na(gene_symbol) &
                        gene_symbol != ""
                ],
                default = NA_character_
            ),
            n_studies = n_distinct(dataset),
            mean_log2FC = mean(log2FC),
            n_positive = sum(log2FC > 0),
            n_negative = sum(log2FC < 0),
            n_zero = sum(log2FC == 0),
            .groups = "drop"
        ) %>%
        mutate(
            consistent = (
                n_studies >= 2 &
                    n_zero == 0 &
                    (
                        n_positive == n_studies |
                            n_negative == n_studies
                    )
            ),
            direction = case_when(
                n_positive == n_studies ~ "Up",
                n_negative == n_studies ~ "Down",
                TRUE ~ "Mixed"
            )
        )
}

primary_consistency <- calculate_consistency(primary_data)
sensitivity_consistency <- calculate_consistency(sensitivity_data)

identify_candidates <- function(consistency) {

    haloperidol <- consistency %>%
        filter(drug == "haloperidol") %>%
        transmute(
            gene_id,
            gene_symbol_haloperidol = gene_symbol,
            haloperidol_mean = mean_log2FC,
            haloperidol_studies = n_studies,
            haloperidol_consistent = consistent
        )

    risperidone <- consistency %>%
        filter(drug == "risperidone") %>%
        transmute(
            gene_id,
            gene_symbol_risperidone = gene_symbol,
            risperidone_mean = mean_log2FC,
            risperidone_studies = n_studies,
            risperidone_consistent = consistent
        )

    inner_join(
        haloperidol,
        risperidone,
        by = "gene_id"
    ) %>%
        mutate(
            gene_symbol = coalesce(
                gene_symbol_haloperidol,
                gene_symbol_risperidone
            ),
            shared_direction = case_when(
                haloperidol_mean > 0 &
                    risperidone_mean > 0 ~ "Both up",

                haloperidol_mean < 0 &
                    risperidone_mean < 0 ~ "Both down",

                haloperidol_mean == 0 |
                    risperidone_mean == 0 ~ "Zero effect",

                TRUE ~ "Opposite"
            )
        ) %>%
        filter(
            haloperidol_studies >= 2,
            risperidone_studies >= 2,
            haloperidol_consistent,
            risperidone_consistent,
            shared_direction %in% c(
                "Both up", "Both down"
            )
        )
}

primary_candidates <- identify_candidates(primary_consistency)
sensitivity_candidates <- identify_candidates(sensitivity_consistency)

write_csv(
    primary_candidates,
    file.path(table_dir, "07_primary_candidates.csv")
)

write_csv(
    sensitivity_candidates,
    file.path(table_dir, "07_sensitivity_candidates.csv")
)

candidate_comparison <- full_join(
    primary_candidates %>%
        transmute(
            gene_id,
            gene_symbol,
            primary_haloperidol = haloperidol_mean,
            primary_risperidone = risperidone_mean,
            primary_direction = shared_direction,
            in_primary = TRUE
        ),

    sensitivity_candidates %>%
        transmute(
            gene_id,
            sensitivity_haloperidol = haloperidol_mean,
            sensitivity_risperidone = risperidone_mean,
            sensitivity_direction = shared_direction,
            in_sensitivity = TRUE
        ),

    by = "gene_id"
) %>%
    mutate(
        in_primary = coalesce(in_primary, FALSE),
        in_sensitivity = coalesce(in_sensitivity, FALSE),

        status = case_when(
            in_primary & in_sensitivity ~ "Retained",
            in_primary & !in_sensitivity ~ "Primary only",
            !in_primary & in_sensitivity ~ "Sensitivity only"
        )
    ) %>%
    arrange(status, gene_id)

write_csv(
    candidate_comparison,
    file.path(table_dir, "07_candidate_comparison.csv")
)

n_primary <- nrow(primary_candidates)
n_sensitivity <- nrow(sensitivity_candidates)

n_retained <- sum(
    candidate_comparison$status == "Retained"
)

n_primary_only <- sum(
    candidate_comparison$status == "Primary only"
)

n_sensitivity_only <- sum(
    candidate_comparison$status == "Sensitivity only"
)

retention_percentage <- if (n_primary > 0) {
    100 * n_retained / n_primary
} else {
    NA_real_
}

jaccard <- if (nrow(candidate_comparison) > 0) {
    n_retained / nrow(candidate_comparison)
} else {
    NA_real_
}

cat("\nCANDIDATE COMPARISON\n")
cat("Primary candidates:", n_primary, "\n")
cat("Sensitivity candidates:", n_sensitivity, "\n")
cat("Retained candidates:", n_retained, "\n")
cat("Primary only:", n_primary_only, "\n")
cat("Sensitivity only:", n_sensitivity_only, "\n")
cat("Retention (%):", round(retention_percentage, 2), "\n")
cat("Jaccard similarity:", round(jaccard, 3), "\n")

genes_of_interest <- c(
    "Penk",
    "Neat1",
    "Ecrg4",
    "Ddit4l",
    "Atp5mk",
    "Itgb6",
    "Frmd7"
)

gene_lookup <- bind_rows(
    primary_data %>%
        select(gene_id, gene_symbol),

    sensitivity_data %>%
        select(gene_id, gene_symbol)
) %>%
    distinct() %>%
    filter(gene_symbol %in% genes_of_interest)

gene_stability <- full_join(
    primary_consistency %>%
        filter(gene_id %in% gene_lookup$gene_id) %>%
        select(
            drug,
            gene_id,
            gene_symbol,
            n_studies,
            mean_log2FC,
            consistent,
            direction
        ) %>%
        rename(
            primary_n_studies = n_studies,
            primary_mean = mean_log2FC,
            primary_consistent = consistent,
            primary_direction = direction
        ),

    sensitivity_consistency %>%
        filter(gene_id %in% gene_lookup$gene_id) %>%
        select(
            drug,
            gene_id,
            n_studies,
            mean_log2FC,
            consistent,
            direction
        ) %>%
        rename(
            sensitivity_n_studies = n_studies,
            sensitivity_mean = mean_log2FC,
            sensitivity_consistent = consistent,
            sensitivity_direction = direction
        ),

    by = c("drug", "gene_id")
)

write_csv(
    gene_stability,
    file.path(table_dir, "07_gene_stability.csv")
)

cat("\nGENES OF INTEREST\n")
print(gene_stability, n = Inf)

retention_plot <- tibble(
    category = c(
        "Retained",
        "Primary only",
        "Sensitivity only"
    ),
    n_genes = c(
        n_retained,
        n_primary_only,
        n_sensitivity_only
    )
)

p1 <- ggplot(
    retention_plot,
    aes(
        x = category,
        y = n_genes,
        fill = category
    )
) +
    geom_col(width = 0.65) +
    geom_text(
        aes(label = n_genes),
        vjust = -0.3
    ) +
    labs(
        title = "Candidate-gene stability after QC sensitivity analysis",
        x = NULL,
        y = "Number of genes"
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none"
    )

ggsave(
    file.path(
        figure_dir,
        "Figure_5_candidate_retention.png"
    ),
    p1,
    width = 8,
    height = 5,
    dpi = 300
)

gse66275 <- comparison %>%
    filter(
        dataset == "GSE66275",
        is.finite(primary_log2FC),
        is.finite(sensitivity_log2FC)
    )

if (nrow(gse66275) > 0) {

    p2 <- ggplot(
        gse66275,
        aes(
            x = primary_log2FC,
            y = sensitivity_log2FC
        )
    ) +
        geom_point(
            alpha = 0.2,
            size = 0.7
        ) +
        geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed"
        ) +
        facet_wrap(~drug, scales = "free") +
        labs(
            title = "GSE66275: Primary versus QC sensitivity",
            x = "Primary log2FC",
            y = "Sensitivity log2FC"
        ) +
        theme_bw(base_size = 12)

    ggsave(
        file.path(
            figure_dir,
            "Figure_6_GSE66275_effect_comparison.png"
        ),
        p2,
        width = 9,
        height = 5,
        dpi = 300
    )

}

retained_genes <- candidate_comparison %>%
    filter(status == "Retained") %>%
    pull(gene_id)

effect_stability <- bind_rows(
    primary_consistency %>%
        filter(gene_id %in% retained_genes) %>%
        mutate(analysis = "Primary"),

    sensitivity_consistency %>%
        filter(gene_id %in% retained_genes) %>%
        mutate(analysis = "Sensitivity")
)

if (nrow(effect_stability) > 0) {

    p3 <- ggplot(
        effect_stability,
        aes(
            x = analysis,
            y = mean_log2FC,
            group = gene_id
        )
    ) +
        geom_line(alpha = 0.08) +
        geom_point(alpha = 0.15, size = 0.5) +
        facet_wrap(~drug, scales = "free_y") +
        labs(
            title = "Effect stability among retained candidates",
            x = NULL,
            y = "Mean log2FC"
        ) +
        theme_bw(base_size = 12)

    ggsave(
        file.path(
            figure_dir,
            "Figure_7_retained_gene_effect_stability.png"
        ),
        p3,
        width = 9,
        height = 5,
        dpi = 300
    )

}

unexpected_changes <- dataset_changes %>%
    filter(
        dataset != "GSE66275",
        n_primary != n_sensitivity |
            n_effects_changed > 0
    )

if (nrow(unexpected_changes) > 0) {

    warning(
        "Differences were detected outside GSE66275. ",
        "Review 07_dataset_changes.csv."
    )

}

if (all(dataset_changes$n_effects_changed == 0) &&
    all(dataset_changes$n_primary ==
        dataset_changes$n_sensitivity)) {

    warning(
        "No differences detected between the two masters. ",
        "Verify that the sensitivity file contains the intended ",
        "alternative GSE66275 analysis."
    )

}

summary_lines <- c(
    "SCRIPT 07: QC SENSITIVITY ANALYSIS",
    paste("Date:", Sys.Date()),
    "",
    paste("Primary observations:", nrow(primary_data)),
    paste("Sensitivity observations:", nrow(sensitivity_data)),
    "",
    paste("Primary shared candidates:", n_primary),
    paste("Sensitivity shared candidates:", n_sensitivity),
    paste("Candidates retained:", n_retained),
    paste("Primary-only candidates:", n_primary_only),
    paste("Sensitivity-only candidates:", n_sensitivity_only),
    paste(
        "Primary candidate retention (%):",
        round(retention_percentage, 2)
    ),
    paste(
        "Jaccard similarity:",
        round(jaccard, 3)
    ),
    "",
    "INTERPRETATION",
    "",
    "Candidate retention measures stability of exploratory",
    "directional-consistency criteria under an alternative",
    "GSE66275 sample-selection analysis.",
    "",
    "Candidate retention is not statistical validation.",
    "The alternative sample selection is a sensitivity",
    "analysis, not proof that excluded samples were invalid.",
    "",
    "No causal or clinical conclusions are supported."
)

writeLines(
    summary_lines,
    file.path(
        output_dir,
        "07_sensitivity_summary.txt"
    )
)
writeLines(
    capture.output(sessionInfo()),
    file.path(
        output_dir,
        "07_session_info.txt"
    )
)
cat("Results saved to:", output_dir, "\n")
