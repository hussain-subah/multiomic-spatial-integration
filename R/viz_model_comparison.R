#' Cross-model comparison figures for the spatial cell-type contrasts
#'
#' Two beta-regression variants were built for the same four contrasts
#' (Amyloid/Disease/Overall/MaxPathology): raw-proportion and
#' Fibroblast-excluded (see `R/contrast_utils.R`,
#' `R/exclude_fibroblast_utils.R`). Both write a
#' `combined_spatial_contrast_summary.csv` in the same long shape (celltype,
#' region, contrast, type, log2_OR, p_adj), so they can be compared directly
#' here.
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2 and dplyr.
#'
#' @keywords internal
NULL


#' Auto-detect a combined-summary file's effect-size column
#'
#' @return A list with `col` (column name) and `label` (axis/legend text).
#' @keywords internal
.effect_size_column_info <- function(df) {
  if ("log2_OR" %in% colnames(df)) {
    list(col = "log2_OR", label = "log2 OR")
  } else {
    stop("No log2_OR column found.", call. = FALSE)
  }
}


#' Significant-hit-rate summary across models x regions, faceted by contrast
#'
#' Reads a named list of `combined_spatial_contrast_summary.csv` paths (one
#' per model; names become the model/legend labels) and plots the percent of
#' cell types called significant (BH-adjusted p < `padj_cutoff`), grouped by
#' region and faceted by contrast `type` (Amyloid/Disease/Overall/
#' MaxPathology) -- so shrinkage in the hit rate as each fix is applied is
#' visible across every contrast at once, not just whichever one motivated
#' the fixes.
#'
#' @param model_csvs Named character vector of
#'   `combined_spatial_contrast_summary.csv` paths; names become model
#'   labels (max 8, per `okabe_ito_palette()`).
#' @param output_dir Figures output directory.
#' @param padj_cutoff Significance threshold.
#' @return The ggplot object (invisibly), or NULL if no input files exist.
#' @export
plot_significant_hit_summary <- function(model_csvs, output_dir, padj_cutoff = 0.05) {
  present <- model_csvs[file.exists(model_csvs)]

  if (length(present) == 0) {
    message("plot_significant_hit_summary: none of the input files exist.")
    return(invisible(NULL))
  }

  # Select only the columns this figure needs before rbinding, in case any
  # input ever picks up an extra/reordered column -- base rbind() errors on
  # data frames whose column sets don't match exactly.
  df <- do.call(rbind, lapply(names(present), function(nm) {
    d <- utils::read.csv(present[[nm]], stringsAsFactors = FALSE)
    d <- d[, c("celltype", "region", "contrast", "type", "p_adj")]
    d$model <- nm
    d
  }))

  summary_df <- df |>
    dplyr::group_by(.data$model, .data$region, .data$type) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_sig = sum(!is.na(.data$p_adj) & .data$p_adj < padj_cutoff),
      pct_sig = 100 * .data$n_sig / .data$n,
      .groups = "drop"
    ) |>
    dplyr::mutate(model = factor(.data$model, levels = names(model_csvs)))

  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = region, y = pct_sig, fill = model)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::facet_wrap(~ type) +
    ggplot2::scale_fill_manual(values = okabe_ito_palette(length(present)), name = NULL) +
    ggplot2::labs(
      title = "Percent of cell types called significant, by model x region",
      subtitle = paste0("BH-adjusted p < ", padj_cutoff, "; faceted by contrast"),
      x = "Region", y = "% cell types significant"
    ) +
    theme_analysis() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  save_figure(p, "model_comparison_significant_hit_rate", output_dir, width = 11, height = 7)
  invisible(p)
}


#' Effect-size comparison scatter between two models, faceted by contrast
#'
#' Reads two `combined_spatial_contrast_summary.csv` files (auto-detecting
#' each one's effect-size column, `log2_OR` or `log2_fold_change`), joins
#' them on celltype x region x contrast, and plots model A's effect size
#' against model B's, colored by which model(s) called that point
#' significant. Faceted by contrast `type` so all four contrasts are visible
#' at once.
#'
#' @param model_a_csv,model_b_csv Paths to the two
#'   `combined_spatial_contrast_summary.csv` files being compared.
#' @param label_a,label_b Axis/legend labels for the two models.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Significance threshold.
#' @return The ggplot object (invisibly), or NULL if inputs are missing or
#'   nothing joins.
#' @export
plot_effect_size_comparison_scatter <- function(model_a_csv, model_b_csv,
                                                label_a, label_b, output_dir,
                                                padj_cutoff = 0.05) {
  if (!file.exists(model_a_csv) || !file.exists(model_b_csv)) {
    message("plot_effect_size_comparison_scatter: missing input for ", label_a, " vs ", label_b, ".")
    return(invisible(NULL))
  }

  a <- utils::read.csv(model_a_csv, stringsAsFactors = FALSE)
  b <- utils::read.csv(model_b_csv, stringsAsFactors = FALSE)

  a_info <- .effect_size_column_info(a)
  b_info <- .effect_size_column_info(b)

  a <- a[, c("celltype", "region", "contrast", "type", a_info$col, "p_adj")]
  colnames(a) <- c("celltype", "region", "contrast", "type", "effect_a", "p_adj_a")

  b <- b[, c("celltype", "region", "contrast", "type", b_info$col, "p_adj")]
  colnames(b) <- c("celltype", "region", "contrast", "type", "effect_b", "p_adj_b")

  joined <- merge(a, b, by = c("celltype", "region", "contrast", "type"))

  if (nrow(joined) == 0) {
    message("plot_effect_size_comparison_scatter: no overlapping rows for ", label_a, " vs ", label_b, ".")
    return(invisible(NULL))
  }

  joined$sig_a <- !is.na(joined$p_adj_a) & joined$p_adj_a < padj_cutoff
  joined$sig_b <- !is.na(joined$p_adj_b) & joined$p_adj_b < padj_cutoff

  cat_levels <- c(
    "Not significant",
    paste("Significant in", label_a, "only"),
    paste("Significant in", label_b, "only"),
    "Significant in both"
  )

  joined$category <- ifelse(
    joined$sig_a & joined$sig_b, cat_levels[4],
    ifelse(joined$sig_a, cat_levels[2],
    ifelse(joined$sig_b, cat_levels[3], cat_levels[1]))
  )
  joined$category <- factor(joined$category, levels = cat_levels)

  cat_colors <- stats::setNames(c("grey70", okabe_ito_palette(3)), cat_levels)

  p <- ggplot2::ggplot(joined, ggplot2::aes(x = effect_a, y = effect_b, color = category)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_point(size = 1.6, alpha = 0.75) +
    ggplot2::facet_wrap(~ type) +
    ggplot2::scale_color_manual(values = cat_colors, name = NULL, drop = FALSE) +
    ggplot2::labs(
      title = paste0(label_a, " vs. ", label_b, ": effect size per celltype x region"),
      subtitle = paste0("Color = which model(s) call it significant (BH p < ", padj_cutoff, ")"),
      x = paste0(label_a, " (", a_info$label, ")"),
      y = paste0(label_b, " (", b_info$label, ")")
    ) +
    theme_analysis()

  safe_a <- gsub("[^A-Za-z0-9]+", "_", label_a)
  safe_b <- gsub("[^A-Za-z0-9]+", "_", label_b)
  save_figure(
    p, paste0("model_comparison_scatter_", safe_a, "_vs_", safe_b),
    output_dir, width = 11, height = 8
  )
  invisible(p)
}
