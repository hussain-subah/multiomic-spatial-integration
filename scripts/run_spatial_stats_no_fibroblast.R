# ============================================================
# Sensitivity analysis: re-run spatial contrasts with Fibroblast excluded
#
# `Fibroblast` claims ~45-64% of the entire cell-type composition (see
# results/cell_proportions sanity checks) -- implausible for a rare
# perivascular cell type in brain tissue -- and its swing between AD-CAA
# and Control is large enough to mechanically inflate the apparent share
# of nearly every other cell type through the sum-to-one constraint alone.
# This script drops Fibroblast, renormalizes the remaining cell types to
# sum to 1 per ROI, and re-runs the same contrast framework as
# run_spatial_stats.R so the two result sets can be compared directly.
# ============================================================

library(glmmTMB)
library(emmeans)
library(dplyr)
library(readr)

source("R/contrast_utils.R")
source("R/exclude_fibroblast_utils.R")

setwd("D:/repo/multiomic-spatial-integration-branchv2")
input_file <- "results/cell_proportions/spatial_celltype_proportions_for_R.csv"
output_dir <- "results/spatial_stats_no_fibroblast"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load, exclude Fibroblast, renormalize, clean
# ============================================================

df_raw <- read.csv(input_file, stringsAsFactors = FALSE)

df_excl <- exclude_and_renormalize(
  df = df_raw,
  exclude = "Fibroblast",
  abundance_col = "rel_abundance",
  roi_col = "ROI_ID",
  celltype_col = "celltype"
)

df <- prepare_spatial_proportion_data(
  df = df_excl,
  abundance_col = "rel_abundance",
  disease_col = "disease_status",
  pathology_col = "pathology",
  region_col = "region",
  scan_col = "Scan_ID",
  celltype_col = "celltype"
)

# ============================================================
# Run contrasts (identical framework to run_spatial_stats.R)
# ============================================================

amyloid_res <- run_amyloid_effect(df = df, abundance_col = "rel_abundance")
disease_res <- run_disease_effect(df = df, abundance_col = "rel_abundance")
overall_res <- run_weighted_overall_effect(df = df, abundance_col = "rel_abundance")
maxpath_res <- run_max_pathology_effect(df = df, abundance_col = "rel_abundance")

write.csv(
  amyloid_res$contrasts,
  file.path(output_dir, "amyloid_effect_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  disease_res$contrasts,
  file.path(output_dir, "disease_effect_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  overall_res$contrasts,
  file.path(output_dir, "overall_effect_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  maxpath_res$contrasts,
  file.path(output_dir, "max_pathology_effect_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  overall_res$weights,
  file.path(output_dir, "overall_effect_region_weights.csv"),
  row.names = FALSE
)

# ============================================================
# Region heterogeneity test
# ============================================================

region_het_res <- run_region_heterogeneity(df = df, abundance_col = "rel_abundance")

write.csv(
  region_het_res$interaction_tests,
  file.path(output_dir, "region_heterogeneity_interaction_tests.csv"),
  row.names = FALSE
)

write.csv(
  region_het_res$posthoc,
  file.path(output_dir, "region_heterogeneity_posthoc.csv"),
  row.names = FALSE
)

# ============================================================
# Combined summary table
# ============================================================

amyloid_df <- format_contrast_summary(amyloid_res$contrasts, contrast_type = "Amyloid")
disease_df <- format_contrast_summary(disease_res$contrasts, contrast_type = "Disease")

overall_df <- format_contrast_summary(
  overall_res$contrasts %>% filter(contrast == "AD_overall_vs_Control"),
  contrast_type = "Overall"
)

maxpath_df <- format_contrast_summary(maxpath_res$contrasts, contrast_type = "MaxPathology")

combined_contrasts <- bind_rows(amyloid_df, disease_df, overall_df, maxpath_df)

write.csv(
  combined_contrasts,
  file.path(output_dir, "combined_spatial_contrast_summary.csv"),
  row.names = FALSE
)

# ============================================================
# Diff against the original (Fibroblast-included) run, if available --
# shows how much excluding Fibroblast moved each contrast.
# ============================================================

original_file <- "results/spatial_stats/combined_spatial_contrast_summary.csv"

if (file.exists(original_file)) {
  original <- read.csv(original_file, stringsAsFactors = FALSE)

  comparison <- combined_contrasts |>
    dplyr::inner_join(
      original,
      by = c("celltype", "region", "contrast", "type"),
      suffix = c("_no_fibroblast", "_original")
    ) |>
    dplyr::mutate(
      log2_OR_delta = .data$log2_OR_no_fibroblast - .data$log2_OR_original
    ) |>
    dplyr::arrange(dplyr::desc(abs(.data$log2_OR_delta)))

  write.csv(
    comparison,
    file.path(output_dir, "comparison_vs_original.csv"),
    row.names = FALSE
  )
} else {
  message(
    "Original run not found at ", original_file,
    " -- skipping comparison_vs_original.csv (run run_spatial_stats.R first if you want it)."
  )
}

message("Fibroblast-excluded sensitivity analysis completed successfully.")
