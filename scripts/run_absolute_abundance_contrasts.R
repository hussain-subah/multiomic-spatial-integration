# ============================================================
# Absolute cell-type abundance contrasts (total-count offset)
#
# Mirrors run_spatial_stats_no_fibroblast.R's shape: same 4-contrast
# framework as run_spatial_stats.R, reused here against log(abs_abundance)
# with a per-ROI total-count offset instead of rel_abundance.
#
# Requires scripts/run_calibration_gate_check.R to have already run --
# every output row is stamped with the gate's gate_status, and
# Disease/Overall/MaxPathology should not be reported as confirmatory
# unless gate_status == "proceed".
#
# Runs twice:
#   - offset_col = "size_factor_mor"  (primary, DESeq2-style)  -> results/spatial_stats_absolute/
#   - offset_col = "total_counts"     (raw-sum sensitivity)    -> results/spatial_stats_absolute_rawlibsize/
# Divergence between the two is itself diagnostic, not discarded.
# ============================================================

library(glmmTMB)
library(emmeans)
library(dplyr)
library(readr)

source("R/contrast_utils.R")
source("R/absolute_abundance_utils.R")

setwd("D:/repo/multiomic-spatial-integration-branchv2")

gate_decision_file <- "results/spatial_stats_absolute/calibration_gate_decision.txt"

if (!file.exists(gate_decision_file)) {
  stop(
    "Calibration gate decision not found at ", gate_decision_file, ". ",
    "Run scripts/run_calibration_gate_check.R first -- every contrast in ",
    "this pipeline is stamped with that gate's decision, and ",
    "Disease/Overall/MaxPathology must not be reported without it.",
    call. = FALSE
  )
}

# ============================================================
# Parse the gate decision file into a per-contrast-type lookup
# ============================================================

decision_lines <- readLines(gate_decision_file)
decision_lines <- decision_lines[grepl("^(Amyloid|Disease|Overall|MaxPathology):", decision_lines)]

gate_status_lookup <- setNames(
  trimws(sub("^[A-Za-z]+:\\s*([A-Za-z_]+).*$", "\\1", decision_lines)),
  trimws(sub("^([A-Za-z]+):.*$", "\\1", decision_lines))
)

message("Calibration gate decisions loaded: ")
print(gate_status_lookup)

# ============================================================
# Load abundance + total-count data (shared across both offset runs)
# ============================================================

abs_long_file <- "results/cell_proportions/roi_celltype_abundance_long_abs.csv"
total_counts_file <- "results/cell_proportions/roi_total_counts.csv"

abs_df_raw <- read.csv(abs_long_file, stringsAsFactors = FALSE)
total_counts_df <- read.csv(total_counts_file, stringsAsFactors = FALSE)

abs_only <- setdiff(abs_df_raw$ROI_ID, total_counts_df$ROI_ID)
counts_only <- setdiff(total_counts_df$ROI_ID, abs_df_raw$ROI_ID)

if (length(abs_only) > 0 || length(counts_only) > 0) {
  stop(
    "Incomplete join between abs_abundance ROIs and total_counts ROIs.\n",
    "In abs_abundance but not total_counts (", length(abs_only), "): ",
    paste(utils::head(abs_only, 10), collapse = ", "), "\n",
    "In total_counts but not abs_abundance (", length(counts_only), "): ",
    paste(utils::head(counts_only, 10), collapse = ", "),
    call. = FALSE
  )
}

merged_df <- dplyr::inner_join(abs_df_raw, total_counts_df, by = "ROI_ID")

# ============================================================
# Helper: run all four contrasts for one offset choice, stamp gate_status,
# write the same file set as spatial_stats_no_fibroblast (per-contrast CSVs
# + region weights + fit status + combined summary).
# ============================================================

run_and_write_absolute_pipeline <- function(offset_col, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  df <- prepare_absolute_abundance_data(
    df = merged_df,
    offset_col = offset_col
  )

  amyloid_res <- run_amyloid_effect_absolute(df)
  disease_res <- run_disease_effect_absolute(df)
  overall_res <- run_weighted_overall_effect_absolute(df)
  maxpath_res <- run_max_pathology_effect_absolute(df)

  stamp_gate <- function(contrast_df, type_label) {
    contrast_df$gate_status <- rep(gate_status_lookup[[type_label]], nrow(contrast_df))
    contrast_df
  }

  amyloid_contrasts <- stamp_gate(amyloid_res$contrasts, "Amyloid")
  disease_contrasts <- stamp_gate(disease_res$contrasts, "Disease")
  overall_contrasts <- stamp_gate(overall_res$contrasts, "Overall")
  maxpath_contrasts <- stamp_gate(maxpath_res$contrasts, "MaxPathology")

  write.csv(
    amyloid_contrasts,
    file.path(output_dir, "amyloid_effect_contrasts.csv"),
    row.names = FALSE
  )

  write.csv(
    disease_contrasts,
    file.path(output_dir, "disease_effect_contrasts.csv"),
    row.names = FALSE
  )

  write.csv(
    overall_contrasts,
    file.path(output_dir, "overall_effect_contrasts.csv"),
    row.names = FALSE
  )

  write.csv(
    maxpath_contrasts,
    file.path(output_dir, "max_pathology_effect_contrasts.csv"),
    row.names = FALSE
  )

  write.csv(
    overall_res$weights,
    file.path(output_dir, "overall_effect_region_weights.csv"),
    row.names = FALSE
  )

  write.csv(
    maxpath_res$fit_status,
    file.path(output_dir, "max_pathology_fit_status.csv"),
    row.names = FALSE
  )

  amyloid_df <- format_absolute_contrast_summary(amyloid_contrasts, contrast_type = "Amyloid")
  disease_df <- format_absolute_contrast_summary(disease_contrasts, contrast_type = "Disease")

  overall_df <- format_absolute_contrast_summary(
    overall_contrasts %>% filter(contrast == "AD_overall_vs_Control"),
    contrast_type = "Overall"
  )

  maxpath_df <- format_absolute_contrast_summary(maxpath_contrasts, contrast_type = "MaxPathology")

  amyloid_df$gate_status <- rep(gate_status_lookup[["Amyloid"]], nrow(amyloid_df))
  disease_df$gate_status <- rep(gate_status_lookup[["Disease"]], nrow(disease_df))
  overall_df$gate_status <- rep(gate_status_lookup[["Overall"]], nrow(overall_df))
  maxpath_df$gate_status <- rep(gate_status_lookup[["MaxPathology"]], nrow(maxpath_df))

  combined_contrasts <- bind_rows(amyloid_df, disease_df, overall_df, maxpath_df)

  write.csv(
    combined_contrasts,
    file.path(output_dir, "combined_spatial_contrast_summary.csv"),
    row.names = FALSE
  )

  message("Absolute-abundance pipeline (offset_col = '", offset_col, "') written to ", output_dir)
}

# ============================================================
# Primary: median-of-ratios size factor
# ============================================================

run_and_write_absolute_pipeline(
  offset_col = "size_factor_mor",
  output_dir = "results/spatial_stats_absolute"
)

# ============================================================
# Sensitivity variant: raw total-count offset (mirrors the Pyro model's
# internal l_r exactly, but is dominated by highly-expressed genes --
# divergence from the median-of-ratios run is itself diagnostic)
# ============================================================

run_and_write_absolute_pipeline(
  offset_col = "total_counts",
  output_dir = "results/spatial_stats_absolute_rawlibsize"
)

message("Absolute-abundance contrasts completed successfully.")
