#' Figures for the pathway-proportion linkage analysis
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2.
#'
#' @keywords internal
NULL


#' GSEA dot plot of top pathways for one cell type
#'
#' Reads an fgsea result table (`pathway`, `NES`, `padj`, `size`) and plots
#' the top pathways by adjusted p-value: NES on the x-axis, dot size by set
#' size, dot color by NES (diverging).
#'
#' @param enrichment_csv Path to a `pathway_enrichment.csv`.
#' @param celltype_label Label for the title/filename.
#' @param output_dir Figures output directory.
#' @param top_n Number of pathways (by padj) to show.
#' @return The ggplot object (invisibly), or NULL if the file is empty.
#' @export
plot_gsea_dotplot <- function(enrichment_csv, celltype_label, output_dir, top_n = 15) {
  df <- utils::read.csv(enrichment_csv, stringsAsFactors = FALSE)

  if (nrow(df) == 0 || !"NES" %in% colnames(df)) {
    message("plot_gsea_dotplot: nothing to plot for ", celltype_label, ".")
    return(invisible(NULL))
  }

  df <- df[order(df$padj), ]
  df <- utils::head(df, top_n)
  df <- df[order(df$NES), ]
  df$pathway <- factor(df$pathway, levels = df$pathway)

  lim <- max(abs(df$NES), na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = pathway)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_point(ggplot2::aes(size = size, color = NES)) +
    scale_color_diverging(name = "NES", limits = c(-lim, lim)) +
    ggplot2::scale_size_continuous(name = "Set size") +
    ggplot2::labs(
      title = paste0("Pathways covarying with ", celltype_label, " abundance"),
      subtitle = paste0("Top ", top_n, " by adjusted p-value"),
      x = "Normalized enrichment score", y = NULL
    ) +
    theme_analysis()

  save_figure(p, paste0("gsea_", celltype_label), output_dir, width = 9, height = 6)
  invisible(p)
}


#' Volcano plot of genes ranked by association with a cell type's proportion
#'
#' @param ranked_csv Path to a `ranked_genes.csv` (`gene`, `statistic`,
#'   `p.value`, `p_adj`).
#' @param celltype_label Label for the title/filename.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Significance threshold for coloring.
#' @param label_top Number of top genes (by p) to label at each end.
#' @return The ggplot object (invisibly).
#' @export
plot_ranked_gene_volcano <- function(ranked_csv, celltype_label, output_dir,
                                     padj_cutoff = 0.05, label_top = 10) {
  df <- utils::read.csv(ranked_csv, stringsAsFactors = FALSE)
  df$neglog_p <- -log10(df$p.value)
  df$direction <- ifelse(
    df$p_adj >= padj_cutoff, "n.s.",
    ifelse(df$statistic > 0, "pos. assoc.", "neg. assoc.")
  )

  sig <- df[df$p_adj < padj_cutoff, ]
  sig <- sig[order(-sig$neglog_p), ]
  top_pos <- utils::head(sig[sig$statistic > 0, ], label_top)
  top_neg <- utils::head(sig[sig$statistic < 0, ], label_top)
  to_label <- rbind(top_pos, top_neg)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = statistic, y = neglog_p, color = direction)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::scale_color_manual(
      values = c("neg. assoc." = "#2166AC", "n.s." = "grey75", "pos. assoc." = "#B2182B"),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste0("Genes associated with ", celltype_label, " abundance"),
      x = "Spearman correlation with proportion",
      y = expression(-log[10]~"(p)")
    ) +
    theme_analysis()

  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(to_label) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = to_label,
      ggplot2::aes(label = gene), size = 2.5, max.overlaps = 20, show.legend = FALSE
    )
  }

  save_figure(p, paste0("volcano_", celltype_label), output_dir, width = 8, height = 6)
  invisible(p)
}


#' Dot plot(s) of top cell-type-restricted, non-housekeeping pathways
#'
#' Plots the output of `top_specific_pathways_by_stratum()`: one row per
#' (pathway, cell type) hit, x = NES, colored by NES, point size by
#' -log10(padj). Unlike `plot_gsea_dotplot()` (one cell type's top pathways
#' by significance alone), every hit shown here has already been filtered to
#' a small global frequency across all (celltype, stratum) combinations --
#' i.e. it is not a generic pathway recurring everywhere.
#'
#' @param top_specific_df Output of `top_specific_pathways_by_stratum()`
#'   (`stratum, celltype, pathway, pval, padj, NES, size, n_significant`).
#' @param output_dir Figures output directory.
#' @param label_width Max characters of the pathway name before truncating.
#' @param combined If `TRUE`, also write one combined figure faceted by
#'   stratum in addition to the per-stratum files.
#' @return Named list of ggplot objects, one per stratum (plus `"combined"`
#'   if requested), invisibly. `NULL` if there's nothing to plot.
#' @export
plot_specificity_dotplot <- function(top_specific_df, output_dir, label_width = 55,
                                     combined = FALSE) {
  if (nrow(top_specific_df) == 0) {
    message("plot_specificity_dotplot: nothing to plot (no rows).")
    return(invisible(NULL))
  }

  df <- top_specific_df
  df$pathway_short <- ifelse(
    nchar(df$pathway) > label_width,
    paste0(substr(df$pathway, 1, label_width), "..."),
    df$pathway
  )
  df$label <- paste0(df$pathway_short, "  [", df$celltype, "]")

  lim <- max(abs(df$NES), na.rm = TRUE)

  .dotplot <- function(sub, title, subtitle) {
    sub <- sub[order(sub$NES), ]
    sub$label <- factor(sub$label, levels = unique(sub$label))

    ggplot2::ggplot(sub, ggplot2::aes(x = NES, y = label)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(ggplot2::aes(size = -log10(padj), color = NES)) +
      scale_color_diverging(name = "NES", limits = c(-lim, lim)) +
      ggplot2::scale_size_continuous(name = expression(-log[10]~"(padj)")) +
      ggplot2::labs(title = title, subtitle = subtitle, x = "Normalized enrichment score", y = NULL) +
      theme_analysis(base_size = 9)
  }

  results <- list()

  for (stratum_name in unique(df$stratum)) {
    sub <- df[df$stratum == stratum_name, ]
    p <- .dotplot(
      sub,
      title = paste0("Cell-type-restricted pathway hits (", stratum_name, ")"),
      subtitle = "pathway [driving cell type]"
    )
    save_figure(
      p, paste0("specificity_dotplot_", stratum_name), output_dir,
      width = 9, height = max(4, 0.35 * nrow(sub))
    )
    results[[stratum_name]] <- p
  }

  if (combined) {
    p_all <- .dotplot(df, "Cell-type-restricted, non-housekeeping pathway hits",
                       "pathway [driving cell type]") +
      ggplot2::facet_wrap(~stratum, scales = "free_y", ncol = 1)
    n_strata <- length(unique(df$stratum))
    save_figure(p_all, "specificity_dotplot_combined", output_dir,
                width = 10, height = max(8, 1.8 * n_strata))
    results[["combined"]] <- p_all
  }

  invisible(results)
}
