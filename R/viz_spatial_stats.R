#' Figures for the spatial statistics analyses
#'
#' Covers the contrast effect sizes, the region-heterogeneity test, and the
#' cell-type co-occurrence analysis. Cell types are always placed on an axis
#' (never color-encoded), since there are ~46 of them.
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2.
#'
#' @keywords internal
NULL


#' Effect-size heatmap across region x cell-type, faceted by contrast
#'
#' Reads the combined contrast summary (celltype, region, contrast_type,
#' log2_OR, p_adj) and draws a diverging heatmap of log2 odds ratios, with an
#' asterisk on cells passing an adjusted-p threshold.
#'
#' @param combined_csv Path to `combined_spatial_contrast_summary.csv`.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Threshold for the significance asterisk.
#' @return The ggplot object (invisibly).
#' @export
plot_contrast_effsize_heatmap <- function(combined_csv, output_dir, padj_cutoff = 0.05) {
  df <- utils::read.csv(combined_csv, stringsAsFactors = FALSE)
  df$sig <- ifelse(!is.na(df$p_adj) & df$p_adj < padj_cutoff, "*", "")

  lim <- max(abs(df$log2_OR), na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = region, y = celltype, fill = log2_OR)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = sig), color = "grey15", size = 3, vjust = 0.75) +
    ggplot2::facet_wrap(~ contrast_type, nrow = 1) +
    scale_fill_diverging(name = "log2 OR", limits = c(-lim, lim)) +
    ggplot2::labs(
      title = "Cell-type abundance effect sizes by region",
      subtitle = paste0("* adjusted p < ", padj_cutoff, "; diverging scale centered at 0"),
      x = "Region", y = NULL
    ) +
    theme_analysis() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  save_figure(p, "contrast_effsize_heatmap", output_dir, width = 11, height = 10)
  invisible(p)
}


#' Forest plot of one contrast's estimates with confidence intervals
#'
#' Reads an individual contrast CSV (with `estimate`, `lower.CL`, `upper.CL`
#' from emmeans) and plots estimate +/- CI per cell type, colored by region
#' (<= 3 regions, so categorical color is safe here).
#'
#' @param contrast_csv Path to an individual `*_effect_contrasts.csv`.
#' @param contrast_label Human-readable label for the title/filename.
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly), or NULL if CI columns are absent.
#' @export
plot_contrast_forest <- function(contrast_csv, contrast_label, output_dir) {
  df <- utils::read.csv(contrast_csv, stringsAsFactors = FALSE)

  if (!all(c("estimate", "lower.CL", "upper.CL") %in% colnames(df))) {
    message("plot_contrast_forest: no CI columns in ", contrast_csv, "; skipping.")
    return(invisible(NULL))
  }

  regions <- sort(unique(df$region))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = celltype, color = region)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = lower.CL, xmax = upper.CL),
      height = 0, linewidth = 0.5,
      position = ggplot2::position_dodge(width = 0.6)
    ) +
    ggplot2::geom_point(size = 1.8, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_color_manual(values = okabe_ito_palette(length(regions))) +
    ggplot2::labs(
      title = paste0(contrast_label, ": estimate +/- 95% CI"),
      x = "Contrast estimate (logit scale)", y = NULL, color = "Region"
    ) +
    theme_analysis()

  save_figure(p, paste0("forest_", contrast_label), output_dir, width = 8, height = 10)
  invisible(p)
}


#' Region-heterogeneity interaction significance bar
#'
#' Reads `region_heterogeneity_interaction_tests.csv` and plots each cell
#' type's group:region interaction as -log10(adjusted p), sorted, with a
#' reference line at the significance threshold.
#'
#' @param interaction_csv Path to the interaction-tests CSV.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Significance threshold for the reference line/flag.
#' @return The ggplot object (invisibly).
#' @export
plot_region_interaction <- function(interaction_csv, output_dir, padj_cutoff = 0.05) {
  df <- utils::read.csv(interaction_csv, stringsAsFactors = FALSE)

  padj_col <- if ("p_adj" %in% colnames(df)) "p_adj" else "p.value"
  df$neglog_p <- -log10(df[[padj_col]])
  df$significant <- df[[padj_col]] < padj_cutoff
  df <- df[order(df$neglog_p), ]
  df$celltype <- factor(df$celltype, levels = df$celltype)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = neglog_p, y = celltype, fill = significant)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(
      xintercept = -log10(padj_cutoff),
      linetype = "dashed", color = "grey50"
    ) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "grey75", "TRUE" = "#B2182B"),
      labels = c("FALSE" = "n.s.", "TRUE" = paste0("adj p < ", padj_cutoff)),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Region heterogeneity of disease/amyloid effects",
      subtitle = "group:region interaction per cell type",
      x = expression(-log[10]~"(adjusted p)"), y = NULL
    ) +
    theme_analysis()

  save_figure(p, "region_interaction_significance", output_dir, width = 8, height = 9)
  invisible(p)
}


#' Cell-type co-occurrence correlation heatmap, faceted by region
#'
#' @param cooccurrence_csv Path to `celltype_cooccurrence_by_region.csv`.
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly).
#' @export
plot_cooccurrence_heatmap <- function(cooccurrence_csv, output_dir) {
  df <- utils::read.csv(cooccurrence_csv, stringsAsFactors = FALSE)

  # symmetrize so the heatmap is full, not triangular
  df_sym <- df
  swap <- df
  swap$celltype_1 <- df$celltype_2
  swap$celltype_2 <- df$celltype_1
  df_sym <- rbind(df_sym, swap)

  p <- ggplot2::ggplot(df_sym, ggplot2::aes(x = celltype_1, y = celltype_2, fill = r)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ region) +
    scale_fill_diverging(name = "CLR corr.", limits = c(-1, 1)) +
    ggplot2::labs(
      title = "Cell-type co-occurrence (CLR correlation) by region",
      x = NULL, y = NULL
    ) +
    theme_analysis(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      axis.text.y = ggplot2::element_text(size = 5)
    )

  save_figure(p, "cooccurrence_heatmap", output_dir, width = 14, height = 6)
  invisible(p)
}


#' Cell-type co-occurrence network of significant pairs (per region)
#'
#' Optional: requires `igraph` and `ggraph`. Draws one network per region,
#' edges = significant CLR correlations, edge color by sign.
#'
#' @param cooccurrence_csv Path to `celltype_cooccurrence_by_region.csv`.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Edge significance threshold.
#' @return The ggplot object (invisibly), or NULL if packages are unavailable.
#' @export
plot_cooccurrence_network <- function(cooccurrence_csv, output_dir, padj_cutoff = 0.05) {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    message("plot_cooccurrence_network: igraph/ggraph not installed; skipping network ",
            "(the heatmap covers the same data).")
    return(invisible(NULL))
  }

  df <- utils::read.csv(cooccurrence_csv, stringsAsFactors = FALSE)
  df <- df[!is.na(df$p_adj) & df$p_adj < padj_cutoff, ]

  if (nrow(df) == 0) {
    message("plot_cooccurrence_network: no significant pairs at adj p < ", padj_cutoff, ".")
    return(invisible(NULL))
  }

  plots <- lapply(sort(unique(df$region)), function(rg) {
    sub <- df[df$region == rg, ]
    g <- igraph::graph_from_data_frame(
      sub[, c("celltype_1", "celltype_2", "r")],
      directed = FALSE
    )
    ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(ggplot2::aes(edge_color = r), edge_width = 0.6) +
      ggraph::geom_node_point(size = 2, color = "grey30") +
      ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 2.5) +
      ggraph::scale_edge_color_gradient2(
        low = "#2166AC", mid = "grey80", high = "#B2182B", midpoint = 0, name = "corr."
      ) +
      ggplot2::labs(title = rg) +
      ggplot2::theme_void()
  })

  combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plots, nrow = 1)
  } else {
    plots[[1]]
  }

  save_figure(combined, "cooccurrence_network", output_dir, width = 15, height = 5)
  invisible(combined)
}
