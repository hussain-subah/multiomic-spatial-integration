# ============================================================
# Figures for the absolute cell-type abundance (total-count offset) pipeline
#
# Covers, for both offset variants (results/spatial_stats_absolute and
# results/spatial_stats_absolute_rawlibsize):
#   - per-pipeline effect-size heatmaps + forest plots, reusing
#     plot_contrast_effsize_heatmap() / plot_contrast_forest()
#     (R/viz_spatial_stats.R) unchanged
#   - the calibration gate diagnostic (results/spatial_stats_absolute only --
#     the gate is run once, upstream of both offset variants)
#   - the flagship absolute-vs-relative log2FC scatter, against
#     results/spatial_stats (the least-confounded proportion-based pipeline)
# ============================================================

library(ggplot2)
library(dplyr)
setwd("D:/repo/multiomic-spatial-integration-branchv2")

source("R/viz_theme.R")
source("R/viz_spatial_stats.R")
source("R/viz_absolute_abundance.R")

fig_root <- "results/figures"

# helper: run a figure call only if its input exists
if_exists <- function(path, expr) {
  if (all(file.exists(path))) {
    force(expr)
  } else {
    message("Skipping (missing input): ", paste(path, collapse = ", "))
  }
}

cell_classes <- c("all", "neuronal", "non_neuronal")

contrast_files <- list(
  Amyloid = "amyloid_effect_contrasts.csv",
  Disease = "disease_effect_contrasts.csv",
  Overall = "overall_effect_contrasts.csv",
  MaxPathology = "max_pathology_effect_contrasts.csv"
)

absolute_dirs <- list(
  "Absolute (offset-adjusted)" = "results/spatial_stats_absolute",
  "Absolute (raw libsize offset)" = "results/spatial_stats_absolute_rawlibsize"
)

# ============================================================
# Per-pipeline heatmaps + forest plots
# ============================================================

for (model_label in names(absolute_dirs)) {
  stats_dir <- absolute_dirs[[model_label]]
  safe_label <- gsub("[^A-Za-z0-9]+", "_", model_label)
  out_dir <- file.path(fig_root, paste0("spatial_stats_", safe_label))

  combined_csv <- file.path(stats_dir, "combined_spatial_contrast_summary.csv")
  if_exists(
    combined_csv,
    for (cc in cell_classes) {
      plot_contrast_effsize_heatmap(combined_csv, out_dir, cell_class = cc)
    }
  )

  for (label in names(contrast_files)) {
    csv <- file.path(stats_dir, contrast_files[[label]])
    for (cc in cell_classes) {
      if_exists(csv, plot_contrast_forest(csv, label, out_dir, cell_class = cc))
    }
  }
}

# ============================================================
# Calibration gate diagnostic
# ============================================================

if_exists(
  c(
    "results/spatial_stats_absolute/calibration_ratio_table.csv",
    "results/spatial_stats_absolute/calibration_gate_summary.csv"
  ),
  plot_calibration_gate_diagnostic(
    "results/spatial_stats_absolute/calibration_ratio_table.csv",
    "results/spatial_stats_absolute/calibration_gate_summary.csv",
    file.path(fig_root, "spatial_stats_absolute")
  )
)

# ============================================================
# Flagship: absolute vs. relative log2FC scatter, against the
# least-confounded proportion-based pipeline (raw beta-regression)
# ============================================================

relative_csv <- "results/spatial_stats/combined_spatial_contrast_summary.csv"

for (model_label in names(absolute_dirs)) {
  stats_dir <- absolute_dirs[[model_label]]
  safe_label <- gsub("[^A-Za-z0-9]+", "_", model_label)
  out_dir <- file.path(fig_root, paste0("spatial_stats_", safe_label))

  if_exists(
    c(file.path(stats_dir, "combined_spatial_contrast_summary.csv"), relative_csv),
    plot_absolute_vs_relative_fc_scatter(
      file.path(stats_dir, "combined_spatial_contrast_summary.csv"),
      relative_csv,
      out_dir
    )
  )
}

message("\nAbsolute-abundance figures written under ", fig_root, ".")
