#' Figure for the leave-one-Scan-out robustness check
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2.
#'
#' @keywords internal
NULL


#' Robustness caterpillar plot for one contrast
#'
#' Reads a `<contrast>_robustness.csv` (full estimate plus the min/max of the
#' leave-one-Scan-out refits) and plots, per region x cell-type, the full
#' estimate as a point with the leave-one-out range as a horizontal bar.
#' Points are colored by severity: a result whose *direction* flips when a
#' single scan is dropped (`any_direction_flip`) is a materially worse finding
#' than one that merely crosses the significance threshold
#' (`any_sig_flip`), so the two are shown as distinct tiers rather than
#' collapsed into one flag. Only the most-affected rows are shown to keep the
#' figure legible.
#'
#' Some robustness files (e.g. `overall_effect_robustness.csv`) stack more
#' than one `contrast` (Disease_effect/Amyloid_effect/AD_overall_vs_Control)
#' for the same region x cell-type pair -- the y-axis label includes the
#' contrast whenever the file has more than one, so those rows don't collide.
#'
#' @param robustness_csv Path to a `<contrast>_robustness.csv`.
#' @param contrast_label Label for the title/filename.
#' @param output_dir Figures output directory.
#' @param top_n Number of most-variable rows (by max_abs_delta) to show.
#' @return The ggplot object (invisibly).
#' @export
plot_robustness_caterpillar <- function(robustness_csv, contrast_label, output_dir,
                                        top_n = 30) {
  df <- utils::read.csv(robustness_csv, stringsAsFactors = FALSE)

  if (nrow(df) == 0) {
    message("plot_robustness_caterpillar: empty input for ", contrast_label, ".")
    return(invisible(NULL))
  }

  if (!"any_sig_flip" %in% colnames(df)) df$any_sig_flip <- FALSE
  if (!"any_direction_flip" %in% colnames(df)) df$any_direction_flip <- FALSE

  df$flip_status <- ifelse(
    df$any_direction_flip, "direction flips",
    ifelse(df$any_sig_flip, "significance flips", "stable")
  )
  df$flip_status <- factor(
    df$flip_status,
    levels = c("stable", "significance flips", "direction flips")
  )

  df <- df[order(-df$max_abs_delta), ]
  df <- utils::head(df, top_n)

  df$label <- if ("contrast" %in% colnames(df) && length(unique(df$contrast)) > 1) {
    paste(df$region, df$celltype, df$contrast, sep = " | ")
  } else {
    paste(df$region, df$celltype, sep = " | ")
  }

  if (anyDuplicated(df$label) > 0) {
    df$label <- make.unique(df$label, sep = " #")
  }

  df <- df[order(df$full_estimate), ]
  df$label <- factor(df$label, levels = df$label)

  p <- ggplot2::ggplot(df, ggplot2::aes(y = label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = min_loo_estimate, xmax = max_loo_estimate),
      color = "grey65", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = full_estimate, color = flip_status), size = 2
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "stable" = "grey30",
        "significance flips" = "#E69F00",
        "direction flips" = "#B2182B"
      ),
      name = NULL, drop = FALSE
    ) +
    ggplot2::labs(
      title = paste0(contrast_label, ": leave-one-Scan-out robustness"),
      subtitle = paste0("Full estimate (point) vs. leave-one-out range (bar); top ",
                        top_n, " most variable"),
      x = "Contrast estimate", y = NULL
    ) +
    theme_analysis()

  save_figure(p, paste0("robustness_", contrast_label), output_dir, width = 9, height = 8)
  invisible(p)
}
