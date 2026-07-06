# ============================================================
# Generate all analysis figures from the R-side result CSVs
#
# Reads the outputs of the spatial stats, region heterogeneity,
# co-occurrence, pathway linkage, marker concordance, robustness, and
# deconvolution comparison analyses, and writes PDF+PNG figures under
# results/figures/. Each section is wrapped so a missing input file skips
# that figure with a message rather than aborting the whole run.
#
# (Pseudobulk-validation figures are produced on the Python side by
# notebooks/06_pseudobulk_validation.py, which now also has a per-cell-type
# recovery bar via python/pseudobulk_validation_utils.plot_recovery_bar.)
# ============================================================

library(ggplot2)

source("R/viz_theme.R")
source("R/viz_spatial_stats.R")
source("R/viz_pathway.R")
source("R/viz_marker_concordance.R")
source("R/viz_robustness.R")
source("R/viz_deconvolution.R")

# ============================================================
# Config -- adjust paths for your environment
# ============================================================

stats_dir <- "results/spatial_stats"
pathway_dir <- "results/pathway_proportion_link"
marker_csv <- "results/marker_concordance/marker_concordance.csv"
deconv_dir <- "results/deconvolution_comparison"
baseline_csv <- "results/cell_proportions/spatial_celltype_proportions_for_R.csv"

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
# Spatial stats: contrast effect sizes + forests
# ============================================================

if_exists(
  file.path(stats_dir, "combined_spatial_contrast_summary.csv"),
  plot_contrast_effsize_heatmap(
    file.path(stats_dir, "combined_spatial_contrast_summary.csv"),
    file.path(fig_root, "spatial_stats")
  )
)

contrast_files <- list(
  Amyloid = "amyloid_effect_contrasts.csv",
  Disease = "disease_effect_contrasts.csv",
  Overall = "overall_effect_contrasts.csv",
  MaxPathology = "max_pathology_effect_contrasts.csv"
)

for (label in names(contrast_files)) {
  csv <- file.path(stats_dir, contrast_files[[label]])
  if_exists(csv, plot_contrast_forest(csv, label, file.path(fig_root, "spatial_stats")))
}

# ============================================================
# Region heterogeneity + co-occurrence
# ============================================================

if_exists(
  file.path(stats_dir, "region_heterogeneity_interaction_tests.csv"),
  plot_region_interaction(
    file.path(stats_dir, "region_heterogeneity_interaction_tests.csv"),
    file.path(fig_root, "spatial_stats")
  )
)

if_exists(
  file.path(stats_dir, "celltype_cooccurrence_by_region.csv"),
  {
    plot_cooccurrence_heatmap(
      file.path(stats_dir, "celltype_cooccurrence_by_region.csv"),
      file.path(fig_root, "cooccurrence")
    )
    plot_cooccurrence_network(
      file.path(stats_dir, "celltype_cooccurrence_by_region.csv"),
      file.path(fig_root, "cooccurrence")
    )
  }
)

# ============================================================
# Pathway linkage: one GSEA dotplot + volcano per cell type
# ============================================================

if (dir.exists(pathway_dir)) {
  ranked_files <- list.files(pathway_dir, pattern = "ranked_genes\\.csv$",
                             recursive = TRUE, full.names = TRUE)
  for (rf in ranked_files) {
    celltype <- basename(dirname(rf))
    plot_ranked_gene_volcano(rf, celltype, file.path(fig_root, "pathway", celltype))

    enr <- file.path(dirname(rf), "pathway_enrichment.csv")
    if (file.exists(enr)) {
      plot_gsea_dotplot(enr, celltype, file.path(fig_root, "pathway", celltype))
    }
  }
} else {
  message("Skipping pathway figures (missing dir): ", pathway_dir)
}

# ============================================================
# Marker concordance
# ============================================================

if_exists(marker_csv, plot_marker_concordance(marker_csv, file.path(fig_root, "marker_concordance")))

# ============================================================
# Robustness caterpillars
# ============================================================

robustness_dir <- file.path(stats_dir, "robustness")
if (dir.exists(robustness_dir)) {
  rob_files <- list.files(robustness_dir, pattern = "_robustness\\.csv$", full.names = TRUE)
  for (rf in rob_files) {
    label <- sub("_robustness\\.csv$", "", basename(rf))
    plot_robustness_caterpillar(rf, label, file.path(fig_root, "robustness"))
  }
} else {
  message("Skipping robustness figures (missing dir): ", robustness_dir)
}

# ============================================================
# Deconvolution comparison
# ============================================================

comp_dir <- file.path(deconv_dir, "comparison_summary")

if_exists(
  file.path(comp_dir, "cross_method_correlation.csv"),
  plot_cross_method_heatmap(
    file.path(comp_dir, "cross_method_correlation.csv"),
    file.path(fig_root, "deconvolution")
  )
)

if_exists(
  file.path(comp_dir, "concordance_vs_cell2location.csv"),
  {
    plot_method_agreement_box(
      file.path(comp_dir, "concordance_vs_cell2location.csv"),
      file.path(fig_root, "deconvolution")
    )
    plot_concordance_heatmap(
      file.path(comp_dir, "concordance_vs_cell2location.csv"),
      file.path(fig_root, "deconvolution")
    )
  }
)

if_exists(
  file.path(comp_dir, "mean_abs_diff_vs_cell2location.csv"),
  plot_mean_abs_diff_bar(
    file.path(comp_dir, "mean_abs_diff_vs_cell2location.csv"),
    file.path(fig_root, "deconvolution")
  )
)

if_exists(
  file.path(comp_dir, "cross_method_correlation.csv"),
  plot_method_dendrogram(
    file.path(comp_dir, "cross_method_correlation.csv"),
    file.path(fig_root, "deconvolution")
  )
)

if_exists(
  baseline_csv,
  if (length(list.files(deconv_dir, pattern = "_proportions\\.csv$")) > 0) {
    plot_method_composition_heatmap(deconv_dir, baseline_csv, file.path(fig_root, "deconvolution"))
    plot_celltype_method_spread(deconv_dir, baseline_csv, file.path(fig_root, "deconvolution"))
    plot_celltype_scatter_vs_baseline(deconv_dir, baseline_csv, file.path(fig_root, "deconvolution"))
  } else {
    message("Skipping composition/spread/scatter figures: no method _proportions.csv files yet.")
  }
)

message("\nAll available figures written under ", fig_root)
