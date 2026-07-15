# ============================================================
# Generate all analysis figures from the R-side result CSVs
#
# Reads the outputs of the spatial stats, region heterogeneity,
# co-occurrence, pathway linkage, marker concordance, robustness, and
# deconvolution comparison analyses, and writes PDF+PNG figures under
# results/figures/. Each section is wrapped so a missing input file skips
# that figure with a message rather than aborting the whole run.
#
# (The true-vs-inferred scatter and abundance-bias figures are produced on
# the Python side by notebooks/06_pseudobulk_validation.py via
# python/pseudobulk_validation_utils.py. The R-side figures below add a
# per-cell-type recovery bar, a rarity-vs-accuracy check, cross-condition
# consistency, and composition stacked bars.)
# ============================================================

library(ggplot2)
setwd("D:/repo/multiomic-spatial-integration-branchv2")

source("R/viz_theme.R")
source("R/viz_spatial_stats.R")
source("R/viz_pathway.R")
source("R/viz_marker_concordance.R")
source("R/viz_robustness.R")
source("R/viz_deconvolution.R")
source("R/viz_pseudobulk_validation.R")

# ============================================================
# Config -- adjust paths for your environment
# ============================================================

stats_dir <- "results/spatial_stats"
pathway_dir <- "results/pathway_proportion_link"
marker_dir <- "results/marker_concordance"
marker_csv <- file.path(marker_dir, "marker_concordance.csv")
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
  file.path(stats_dir, "celltype_cooccurrence_by_group.csv"),
  {
    plot_cooccurrence_heatmap(
      file.path(stats_dir, "celltype_cooccurrence_by_group.csv"),
      file.path(fig_root, "cooccurrence")
    )
    plot_cooccurrence_heatmap_by_stratum(
      file.path(stats_dir, "celltype_cooccurrence_by_group.csv"),
      file.path(fig_root, "cooccurrence", "by_stratum")
    )
    plot_cooccurrence_network(
      file.path(stats_dir, "celltype_cooccurrence_by_group.csv"),
      file.path(fig_root, "cooccurrence")
    )
    plot_cooccurrence_network_by_stratum(
      file.path(stats_dir, "celltype_cooccurrence_by_group.csv"),
      file.path(fig_root, "cooccurrence", "by_stratum")
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
    stratum <- basename(dirname(dirname(rf)))
    label <- if (identical(stratum, "all")) celltype else paste0(celltype, "_", stratum)
    fig_dir <- if (identical(stratum, "all")) {
      file.path(fig_root, "pathway", celltype)
    } else {
      file.path(fig_root, "pathway", celltype, stratum)
    }

    plot_ranked_gene_volcano(rf, label, fig_dir)

    enr <- file.path(dirname(rf), "pathway_enrichment.csv")
    if (file.exists(enr)) {
      plot_gsea_dotplot(enr, label, fig_dir)
    }
  }
} else {
  message("Skipping pathway figures (missing dir): ", pathway_dir)
}

# ============================================================
# Marker concordance
# ============================================================

if_exists(marker_csv, plot_marker_concordance(marker_csv, file.path(fig_root, "marker_concordance")))

# Cross-strata replication (pooled vs. disease-specific) -- surfaces cell
# types whose "significant" pooled concordance is actually carried by a
# single disease stratum, or whose sign disagrees between strata.
strata_csvs <- c(
  Combined = marker_csv,
  "AD+CAA" = file.path(marker_dir, "marker_concordance_AD_CAA.csv"),
  Control = file.path(marker_dir, "marker_concordance_Control.csv")
)

if (sum(file.exists(strata_csvs)) >= 2) {
  plot_marker_concordance_strata(strata_csvs, file.path(fig_root, "marker_concordance"))
} else {
  message("Skipping marker concordance strata figure (need >=2 of): ",
          paste(strata_csvs, collapse = ", "))
}

# Factor/marker confound check (scripts/check_marker_factor_mapping.R output)
factor_corr_csv <- file.path(marker_dir, "factor_marker_correlation_matrix.csv")
factor_best_csv <- file.path(marker_dir, "factor_marker_best_matches.csv")

if_exists(
  c(factor_corr_csv, factor_best_csv),
  plot_factor_marker_heatmap(factor_corr_csv, factor_best_csv, file.path(fig_root, "marker_concordance"))
)

if_exists(
  factor_best_csv,
  plot_factor_marker_best_match_bar(factor_best_csv, file.path(fig_root, "marker_concordance"))
)

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

# ============================================================
# Pseudobulk deconvolution validation
# ============================================================

pseudobulk_dir <- "results/pseudobulk_validation"
pseudobulk_conditions <- c("AD+CAA", "Control")

for (cond in pseudobulk_conditions) {
  metrics_csv <- file.path(pseudobulk_dir, cond, "recovery_metrics.csv")
  gt_csv <- file.path(pseudobulk_dir, cond, "ground_truth_proportions.csv")
  inferred_csv <- file.path(pseudobulk_dir, cond, "inferred_proportions.csv")

  if_exists(
    metrics_csv,
    plot_recovery_metric_bar(metrics_csv, cond, file.path(fig_root, "pseudobulk_validation"))
  )

  if_exists(
    c(metrics_csv, gt_csv),
    plot_rarity_vs_accuracy(metrics_csv, gt_csv, cond, file.path(fig_root, "pseudobulk_validation"))
  )

  if_exists(
    c(gt_csv, inferred_csv),
    plot_composition_stacked_bars(gt_csv, inferred_csv, cond, file.path(fig_root, "pseudobulk_validation"))
  )
}

metrics_csvs <- file.path(pseudobulk_dir, pseudobulk_conditions, "recovery_metrics.csv")
if_exists(
  metrics_csvs,
  plot_condition_consistency(
    metrics_csvs[1], metrics_csvs[2],
    pseudobulk_conditions[1], pseudobulk_conditions[2],
    file.path(fig_root, "pseudobulk_validation")
  )
)

message("\nAll available figures written under ", fig_root)
