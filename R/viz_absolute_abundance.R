#' Figures for the absolute cell-type abundance (total-count offset) pipeline
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2 and dplyr.
#'
#' @keywords internal
NULL


#' Calibration gate diagnostic figure
#'
#' The figure a reviewer needs before trusting anything else from the
#' absolute-abundance pipeline: per-ROI `calibration_ratio` (total
#' `abs_abundance` / `total_counts`) by `disease_status`, faceted by region,
#' annotated with the fitted `disease_status` coefficient and the
#' plausibility threshold used to make the gate decision.
#'
#' @param calib_csv Path to `calibration_ratio_table.csv`
#'   (`compute_calibration_ratio_table()` output).
#' @param gate_summary_csv Path to `calibration_gate_summary.csv`
#'   (`run_calibration_gate_check()` output, written by
#'   `scripts/run_calibration_gate_check.R`).
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly), or NULL if inputs are missing.
#' @export
plot_calibration_gate_diagnostic <- function(calib_csv, gate_summary_csv, output_dir) {
  if (!file.exists(calib_csv) || !file.exists(gate_summary_csv)) {
    message("plot_calibration_gate_diagnostic: missing input(s); skipping.")
    return(invisible(NULL))
  }

  calib_df <- utils::read.csv(calib_csv, stringsAsFactors = FALSE)
  gate_summary <- utils::read.csv(gate_summary_csv, stringsAsFactors = FALSE)

  decision <- gate_summary$decision[1]
  coef_log2 <- gate_summary$disease_coefficient_log2[1]
  threshold_log2 <- gate_summary$biological_threshold_log2[1]

  disease_levels <- c("Control", "AD-CAA")
  calib_df$disease_status <- factor(calib_df$disease_status, levels = disease_levels)

  p <- ggplot2::ggplot(calib_df, ggplot2::aes(x = disease_status, y = calibration_ratio, color = disease_status)) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.4) +
    ggplot2::geom_boxplot(width = 0.4, alpha = 0, outlier.shape = NA, color = "grey30") +
    ggplot2::facet_wrap(~ region) +
    ggplot2::scale_color_manual(values = okabe_ito_palette(2), name = NULL) +
    ggplot2::labs(
      title = "Calibration gate: total abs_abundance / total_counts by disease status",
      subtitle = paste0(
        "disease_status coefficient = ", round(coef_log2, 3), " log2 ",
        "(plausibility threshold = ", threshold_log2, ") -- decision: ", decision
      ),
      x = NULL, y = "Calibration ratio (total abs_abundance / total_counts)"
    ) +
    theme_analysis() +
    ggplot2::theme(legend.position = "none")

  save_figure(p, "calibration_gate_diagnostic", output_dir, width = 10, height = 5)
  invisible(p)
}


#' Absolute vs. relative log2 fold-change scatter (flagship figure)
#'
#' Joins this pipeline's `combined_spatial_contrast_summary.csv` against a
#' relative-abundance (proportion-based) pipeline's, on celltype x region x
#' contrast `type`, and plots absolute log2FC (x) against relative log2FC
#' (y), with a y = x reference diagonal:
#' - Near the diagonal: genuine cell-type-specific change, consistent both ways.
#' - Near the x-axis: "everything expanded together" (total RNA/cellularity
#'   shifted; this cell type's *share* didn't really change).
#' - Near the y-axis: the closure-artifact signature this pipeline exists to
#'   catch -- some *other* category's absolute swing mechanically dragging
#'   this one's proportion.
#'
#' @param absolute_csv Path to the absolute-abundance
#'   `combined_spatial_contrast_summary.csv` (has `log2_fold_change`).
#' @param relative_csv Path to a relative-abundance
#'   `combined_spatial_contrast_summary.csv` (has `log2_OR`). Default is
#'   `results/spatial_stats/combined_spatial_contrast_summary.csv`, the
#'   least-confounded proportion-based pipeline (raw beta-regression, no
#'   Fibroblast exclusion) -- the most direct "before" to compare against.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Significance threshold.
#' @return The ggplot object (invisibly), or NULL if inputs are missing or
#'   nothing joins.
#' @export
plot_absolute_vs_relative_fc_scatter <- function(absolute_csv,
                                                 relative_csv = "results/spatial_stats/combined_spatial_contrast_summary.csv",
                                                 output_dir,
                                                 padj_cutoff = 0.05) {
  if (!file.exists(absolute_csv) || !file.exists(relative_csv)) {
    message("plot_absolute_vs_relative_fc_scatter: missing input(s); skipping.")
    return(invisible(NULL))
  }

  abs_df <- utils::read.csv(absolute_csv, stringsAsFactors = FALSE)
  rel_df <- utils::read.csv(relative_csv, stringsAsFactors = FALSE)

  abs_df <- abs_df[, c("celltype", "region", "contrast", "type", "log2_fold_change", "p_adj")]
  colnames(abs_df) <- c("celltype", "region", "contrast", "type", "effect_absolute", "p_adj_absolute")

  rel_df <- rel_df[, c("celltype", "region", "contrast", "type", "log2_OR", "p_adj")]
  colnames(rel_df) <- c("celltype", "region", "contrast", "type", "effect_relative", "p_adj_relative")

  joined <- merge(abs_df, rel_df, by = c("celltype", "region", "contrast", "type"))

  if (nrow(joined) == 0) {
    message("plot_absolute_vs_relative_fc_scatter: no overlapping rows.")
    return(invisible(NULL))
  }

  joined$sig_absolute <- !is.na(joined$p_adj_absolute) & joined$p_adj_absolute < padj_cutoff
  joined$sig_relative <- !is.na(joined$p_adj_relative) & joined$p_adj_relative < padj_cutoff

  cat_levels <- c(
    "Not significant",
    "Significant in absolute only",
    "Significant in relative only",
    "Significant in both"
  )

  joined$category <- ifelse(
    joined$sig_absolute & joined$sig_relative, cat_levels[4],
    ifelse(joined$sig_absolute, cat_levels[2],
    ifelse(joined$sig_relative, cat_levels[3], cat_levels[1]))
  )
  joined$category <- factor(joined$category, levels = cat_levels)

  cat_colors <- stats::setNames(c("grey70", okabe_ito_palette(3)), cat_levels)

  lim <- max(abs(c(joined$effect_absolute, joined$effect_relative)), na.rm = TRUE)

  p <- ggplot2::ggplot(joined, ggplot2::aes(x = effect_absolute, y = effect_relative, color = category)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    ggplot2::geom_point(size = 1.6, alpha = 0.75) +
    ggplot2::facet_wrap(~ type) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::scale_color_manual(values = cat_colors, name = NULL, drop = FALSE) +
    ggplot2::labs(
      title = "Absolute (offset-adjusted) vs. relative (proportion) log2 fold-change",
      subtitle = paste0(
        "Dashed line = y = x. Near diagonal = genuine cell-type-specific change; ",
        "near x-axis = total RNA/cellularity shift; near y-axis = closure artifact. ",
        "Color = which model(s) call it significant (BH p < ", padj_cutoff, ")"
      ),
      x = "Absolute log2 fold-change (offset-adjusted)",
      y = "Relative log2 fold-change (proportion-based)"
    ) +
    theme_analysis()

  save_figure(p, "absolute_vs_relative_fc_scatter", output_dir, width = 11, height = 9)
  invisible(p)
}
