# ============================================================
# Figures for the Fibroblast-excluded beta-regression sensitivity pipeline,
# plus a cross-model comparison against the original raw-proportion model
#
# scripts/run_all_figures.R already covers the original beta-regression
# pipeline (results/spatial_stats/). This script covers:
#   - per-pipeline effect-size heatmaps + forest plots for the
#     Fibroblast-excluded pipeline -- reusing plot_contrast_effsize_heatmap()
#     / plot_contrast_forest() from R/viz_spatial_stats.R
#   - cross-model comparison figures: does the significant-hit rate shrink
#     once Fibroblast is excluded, and where did the effect sizes actually
#     move, across all four contrasts at once (not just the disease effect
#     that originally motivated the fix)
# ============================================================

library(ggplot2)
library(dplyr)
setwd("D:/repo/multiomic-spatial-integration-branchv2")

source("R/viz_theme.R")
source("R/viz_spatial_stats.R")
source("R/viz_model_comparison.R")

fig_root <- "results/figures"

# helper: run a figure call only if its input exists
if_exists <- function(path, expr) {
  if (all(file.exists(path))) {
    force(expr)
  } else {
    message("Skipping (missing input): ", paste(path, collapse = ", "))
  }
}

# ============================================================
# Config: the two comparable beta pipelines and their contrast files
# ============================================================

model_dirs <- list(
  "Beta (raw proportions)" = "results/spatial_stats",
  "Beta (Fibroblast excluded)" = "results/spatial_stats_no_fibroblast",
  "Absolute (offset-adjusted)" = "results/spatial_stats_absolute"
)

combined_csvs <- lapply(model_dirs, function(d) file.path(d, "combined_spatial_contrast_summary.csv"))

contrast_files <- list(
  Amyloid = "amyloid_effect_contrasts.csv",
  Disease = "disease_effect_contrasts.csv",
  Overall = "overall_effect_contrasts.csv",
  MaxPathology = "max_pathology_effect_contrasts.csv"
)

cell_classes <- c("all", "neuronal", "non_neuronal")

# ============================================================
# Per-pipeline figures for the Fibroblast-excluded sensitivity pipeline
# (Beta raw proportions was already done by run_all_figures.R)
# ============================================================

model_label <- "Beta (Fibroblast excluded)"
stats_dir <- model_dirs[[model_label]]
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

# ============================================================
# Cross-model comparison: does the significant-hit rate shrink once
# Fibroblast is excluded, across all four contrasts (not just Disease_effect)?
# ============================================================

if_exists(
  unlist(combined_csvs),
  plot_significant_hit_summary(unlist(combined_csvs), file.path(fig_root, "model_comparison"))
)

# ============================================================
# Cross-model comparison: where did the effect sizes actually move?
# ============================================================

if_exists(
  c(combined_csvs[["Beta (raw proportions)"]], combined_csvs[["Beta (Fibroblast excluded)"]]),
  plot_effect_size_comparison_scatter(
    combined_csvs[["Beta (raw proportions)"]], combined_csvs[["Beta (Fibroblast excluded)"]],
    "Beta (raw proportions)", "Beta (Fibroblast excluded)",
    file.path(fig_root, "model_comparison")
  )
)

# Beta (raw proportions) vs. Absolute (offset-adjusted): does modeling
# absolute abundance with a total-count offset move effect sizes/hit rates
# the same way as excluding Fibroblast did, or differently?
if_exists(
  c(combined_csvs[["Beta (raw proportions)"]], combined_csvs[["Absolute (offset-adjusted)"]]),
  plot_effect_size_comparison_scatter(
    combined_csvs[["Beta (raw proportions)"]], combined_csvs[["Absolute (offset-adjusted)"]],
    "Beta (raw proportions)", "Absolute (offset-adjusted)",
    file.path(fig_root, "model_comparison")
  )
)

message(
  "\nModel-comparison figures written under ", file.path(fig_root, "model_comparison"),
  " and per-pipeline directories under ", fig_root, "."
)
