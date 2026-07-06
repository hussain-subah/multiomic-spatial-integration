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
