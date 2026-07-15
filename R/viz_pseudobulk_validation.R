#' Figures for the pseudobulk deconvolution validation
#'
#' Reads the outputs of `notebooks/06_pseudobulk_validation.py`
#' (`results/pseudobulk_validation/<condition>/`: `recovery_metrics.csv`,
#' `ground_truth_proportions.csv`, `inferred_proportions.csv`) and produces
#' diagnostic figures beyond the scatter/bias plots already written on the
#' Python side (`python/pseudobulk_validation_utils.py`).
#'
#' Requires `R/viz_theme.R` to be sourced.
#'
#' @keywords internal
NULL


#' Sorted bar chart of a per-cell-type recovery metric
#'
#' One bar per cell type, sorted by value. For correlation metrics
#' (`pearson_r`, `spearman_rho`) bars are colored by sign (diverging);
#' otherwise a single sequential color is used. Complements the scatter/bias
#' plots with a single at-a-glance ranking of which cell types the model
#' recovers well vs. poorly.
#'
#' @param metrics_csv Path to `recovery_metrics.csv` (output of
#'   `evaluate_recovery()`); the `overall` row is dropped.
#' @param condition_label Label used in the title and output filename.
#' @param output_dir Figures output directory.
#' @param metric Column to plot: `"pearson_r"`, `"spearman_rho"`, `"rmse"`, or
#'   `"mae"`.
#' @return The ggplot object (invisibly).
#' @export
plot_recovery_metric_bar <- function(metrics_csv, condition_label, output_dir,
                                     metric = "pearson_r") {
  df <- utils::read.csv(metrics_csv, stringsAsFactors = FALSE)
  df <- df[df$celltype != "overall", ]
  df <- df[!is.na(df[[metric]]), ]
  df <- df[order(df[[metric]]), ]
  df$celltype <- factor(df$celltype, levels = df$celltype)

  is_corr <- metric %in% c("pearson_r", "spearman_rho")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[metric]], y = celltype))

  if (is_corr) {
    df$sign <- ifelse(df[[metric]] < 0, "negative", "positive")
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[metric]], y = celltype, fill = sign)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::scale_fill_manual(
        values = c(negative = "#2166AC", positive = "#009E73"),
        guide = "none"
      )
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[metric]], y = celltype, fill = .data[[metric]])) +
      ggplot2::geom_col(width = 0.7) +
      scale_fill_magnitude(name = NULL) +
      ggplot2::guides(fill = "none")
  }

  p <- p +
    ggplot2::labs(
      title = paste0("Pseudobulk recovery by cell type (", condition_label, ")"),
      subtitle = metric,
      x = metric, y = NULL
    ) +
    theme_analysis(base_size = 9)

  save_figure(
    p, paste0("recovery_bar_", metric, "_", condition_label),
    output_dir, width = 7, height = max(4, 0.22 * nrow(df))
  )
  invisible(p)
}


#' Rarity vs. accuracy: mean true abundance vs. recovery correlation
#'
#' Tests whether recovery failures are concentrated in rare, low-variance
#' cell types (points near zero mean abundance) rather than being spread
#' uniformly across the abundance range. Cell types with zero-variance
#' ground truth (`pearson_r` undefined) are dropped with a message.
#'
#' @param metrics_csv Path to `recovery_metrics.csv`.
#' @param ground_truth_csv Path to `ground_truth_proportions.csv` (wide:
#'   `roi_id`, `donor_roi_id`, then one column per cell type).
#' @param condition_label Label used in the title and output filename.
#' @param output_dir Figures output directory.
#' @param metric Correlation column to plot against abundance.
#' @return The merged data frame (invisibly).
#' @export
plot_rarity_vs_accuracy <- function(metrics_csv, ground_truth_csv, condition_label,
                                    output_dir, metric = "pearson_r") {
  metrics <- utils::read.csv(metrics_csv, stringsAsFactors = FALSE)
  metrics <- metrics[metrics$celltype != "overall", ]

  gt <- utils::read.csv(ground_truth_csv, stringsAsFactors = FALSE, check.names = FALSE)
  celltype_cols <- setdiff(colnames(gt), c("roi_id", "donor_roi_id"))
  mean_abund <- data.frame(
    celltype = celltype_cols,
    mean_true = colMeans(gt[, celltype_cols, drop = FALSE]),
    stringsAsFactors = FALSE
  )

  merged <- merge(metrics, mean_abund, by = "celltype")

  n_dropped <- sum(is.na(merged[[metric]]))
  if (n_dropped > 0) {
    message(
      n_dropped, " cell type(s) dropped from rarity-vs-accuracy plot (",
      condition_label, "): zero-variance ground truth, ", metric, " undefined."
    )
  }
  plot_df <- merged[!is.na(merged[[metric]]), ]

  ct <- stats::cor.test(plot_df$mean_true, plot_df[[metric]], method = "spearman")
  subtitle <- sprintf("Spearman rho = %.2f (rarity vs. recovery)", unname(ct$estimate))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = mean_true, y = .data[[metric]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#D55E00", fill = "#D55E00", alpha = 0.15, linewidth = 0.7) +
    ggplot2::geom_point(size = 1.6, alpha = 0.75, color = "#0072B2") +
    ggrepel::geom_text_repel(
      data = rbind(utils::head(plot_df[order(-plot_df[[metric]]), ], 3),
                   utils::head(plot_df[order(plot_df[[metric]]), ], 3)),
      ggplot2::aes(label = celltype), size = 2.6, max.overlaps = 20
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = paste0("Recovery accuracy vs. cell-type abundance (", condition_label, ")"),
      subtitle = subtitle,
      x = "Mean true proportion (log scale)", y = metric
    ) +
    theme_analysis()

  save_figure(p, paste0("rarity_vs_accuracy_", condition_label), output_dir, width = 7, height = 5.5)
  invisible(merged)
}


#' Cross-condition consistency of per-cell-type recovery
#'
#' Scatter of a recovery metric in one condition vs. another, one point per
#' shared cell type. Points concordant in sign (same quadrant on both sides
#' of zero) indicate a cell type that is consistently well- or
#' poorly-recovered regardless of condition (a signature/model property);
#' discordant points suggest condition-specific data issues instead.
#'
#' @param metrics_csv_a,metrics_csv_b Paths to the two conditions'
#'   `recovery_metrics.csv`.
#' @param label_a,label_b Condition labels for axis titles and filename.
#' @param output_dir Figures output directory.
#' @param metric Metric column to compare (default `"pearson_r"`).
#' @return The merged data frame (invisibly).
#' @export
plot_condition_consistency <- function(metrics_csv_a, metrics_csv_b, label_a, label_b,
                                       output_dir, metric = "pearson_r") {
  a <- utils::read.csv(metrics_csv_a, stringsAsFactors = FALSE)
  b <- utils::read.csv(metrics_csv_b, stringsAsFactors = FALSE)
  a <- a[a$celltype != "overall", c("celltype", metric)]
  b <- b[b$celltype != "overall", c("celltype", metric)]
  colnames(a)[2] <- "metric_a"
  colnames(b)[2] <- "metric_b"

  merged <- merge(a, b, by = "celltype")
  merged <- merged[stats::complete.cases(merged), ]

  merged$quadrant <- ifelse(
    merged$metric_a >= 0 & merged$metric_b >= 0, "concordant (+/+)",
    ifelse(merged$metric_a < 0 & merged$metric_b < 0, "concordant (-/-)", "discordant")
  )

  ct <- stats::cor.test(merged$metric_a, merged$metric_b, method = "spearman")
  subtitle <- sprintf(
    "Spearman rho = %.2f | %d/%d cell types concordant in sign",
    unname(ct$estimate), sum(merged$quadrant != "discordant"), nrow(merged)
  )

  axis_lim <- max(abs(c(merged$metric_a, merged$metric_b))) * 1.1

  p <- ggplot2::ggplot(merged, ggplot2::aes(x = metric_a, y = metric_b, color = quadrant)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
    ggplot2::geom_point(size = 1.8, alpha = 0.85) +
    ggrepel::geom_text_repel(
      data = merged[merged$quadrant == "discordant", ],
      ggplot2::aes(label = celltype), size = 2.6, max.overlaps = 20, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("concordant (+/+)" = "#009E73", "concordant (-/-)" = "#2166AC", "discordant" = "#D55E00"),
      name = NULL
    ) +
    ggplot2::coord_equal(xlim = c(-axis_lim, axis_lim), ylim = c(-axis_lim, axis_lim)) +
    ggplot2::labs(
      title = paste0("Recovery consistency: ", label_a, " vs. ", label_b),
      subtitle = subtitle,
      x = paste0(metric, " (", label_a, ")"), y = paste0(metric, " (", label_b, ")")
    ) +
    theme_analysis()

  save_figure(
    p, paste0("condition_consistency_", metric, "_", label_a, "_vs_", label_b),
    output_dir, width = 7.5, height = 7
  )
  invisible(merged)
}


#' True vs. inferred composition, stacked bars for representative ROIs
#'
#' More intuitive than the per-cell-type scatter for showing where
#' compositional mass gets misallocated: for a handful of ROIs, draws paired
#' "True" / "Inferred" stacked bars. Only the `top_n_celltypes` most abundant
#' cell types (by mean true proportion) are individually colored; the rest
#' are pooled into "Other" so the categorical palette stays within the
#' colorblind-safe 8-color limit (see `R/viz_theme.R`).
#'
#' @param ground_truth_csv Path to `ground_truth_proportions.csv`.
#' @param inferred_csv Path to `inferred_proportions.csv` (wide, first column
#'   is the ROI id, unnamed in the header).
#' @param condition_label Label used in the title and output filename.
#' @param output_dir Figures output directory.
#' @param n_rois Number of ROIs to show.
#' @param top_n_celltypes Number of cell types individually colored (<= 7,
#'   leaving one slot for "Other").
#' @param seed Random seed for ROI sampling.
#' @return The long-format plotting data frame (invisibly).
#' @export
plot_composition_stacked_bars <- function(ground_truth_csv, inferred_csv, condition_label,
                                          output_dir, n_rois = 6, top_n_celltypes = 7,
                                          seed = 0) {
  gt <- utils::read.csv(ground_truth_csv, stringsAsFactors = FALSE, check.names = FALSE)
  inf <- utils::read.csv(inferred_csv, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

  celltype_cols <- setdiff(colnames(gt), c("roi_id", "donor_roi_id"))
  rownames(gt) <- gt$roi_id

  top_n_celltypes <- min(top_n_celltypes, 7)
  mean_true <- sort(colMeans(gt[, celltype_cols, drop = FALSE]), decreasing = TRUE)
  top_celltypes <- names(utils::head(mean_true, top_n_celltypes))

  set.seed(seed)
  roi_ids <- sample(rownames(gt), min(n_rois, nrow(gt)))

  .bucket_long <- function(mat, roi_ids, source_label) {
    do.call(rbind, lapply(roi_ids, function(rid) {
      row <- mat[rid, celltype_cols, drop = TRUE]
      top_vals <- as.numeric(row[top_celltypes])
      other_val <- sum(as.numeric(row)) - sum(top_vals)
      data.frame(
        roi_id = rid,
        source = source_label,
        celltype = c(top_celltypes, "Other"),
        proportion = c(top_vals, other_val),
        stringsAsFactors = FALSE
      )
    }))
  }

  long_df <- rbind(
    .bucket_long(gt, roi_ids, "True"),
    .bucket_long(inf, roi_ids, "Inferred")
  )
  long_df$source <- factor(long_df$source, levels = c("True", "Inferred"))
  long_df$celltype <- factor(long_df$celltype, levels = c(top_celltypes, "Other"))

  palette <- c(okabe_ito_palette(top_n_celltypes), "grey80")
  names(palette) <- c(top_celltypes, "Other")

  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = source, y = proportion, fill = celltype)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~roi_id, nrow = 1) +
    ggplot2::scale_fill_manual(values = palette, name = NULL) +
    ggplot2::labs(
      title = paste0("True vs. inferred composition (", condition_label, ")"),
      subtitle = paste0("Top ", top_n_celltypes, " cell types by mean true abundance; rest pooled as \"Other\""),
      x = NULL, y = "Proportion"
    ) +
    theme_analysis(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  save_figure(p, paste0("composition_stacked_bars_", condition_label), output_dir, width = 10, height = 5)
  invisible(long_df)
}
