#' Figure for the marker-gene concordance check
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2.
#'
#' @keywords internal
NULL


#' Marker-concordance bar chart (QC view)
#'
#' Reads `marker_concordance.csv` and draws, per cell type, the correlation
#' between its inferred proportion and its own marker score -- sorted, with
#' cell types that fail (no positive, significant concordance) flagged. A
#' clean go/no-go QC figure for the deconvolution's per-cell-type quality.
#'
#' @param concordance_csv Path to `marker_concordance.csv`.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Significance threshold defining a pass.
#' @return The ggplot object (invisibly).
#' @export
plot_marker_concordance <- function(concordance_csv, output_dir, padj_cutoff = 0.05) {
  df <- utils::read.csv(concordance_csv, stringsAsFactors = FALSE)

  df$pass <- !is.na(df$statistic) & df$statistic > 0 &
    !is.na(df$p_adj) & df$p_adj < padj_cutoff
  df$statistic[is.na(df$statistic)] <- 0

  df <- df[order(df$statistic), ]
  df$celltype <- factor(df$celltype, levels = df$celltype)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = statistic, y = celltype, fill = pass)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, color = "grey50") +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "#D55E00", "TRUE" = "#009E73"),
      labels = c("FALSE" = "flag (no pos. sig. concordance)", "TRUE" = "pass"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Marker-gene concordance by cell type",
      subtitle = "Correlation of inferred proportion with own marker score",
      x = "Concordance correlation", y = NULL
    ) +
    theme_analysis()

  save_figure(p, "marker_concordance", output_dir, width = 8, height = 9)
  invisible(p)
}
