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


#' Cross-stratum marker-concordance comparison (pooled vs. disease-specific)
#'
#' Reads the pooled `marker_concordance.csv` alongside its disease-stratified
#' siblings (`marker_concordance_AD_CAA.csv`, `marker_concordance_Control.csv`,
#' as written by `scripts/run_marker_concordance_check.R`) and draws a
#' connected dot ("slope") chart of each cell type's concordance statistic
#' across strata. A pooled-significant result driven by only one disease
#' stratum, or whose sign disagrees between strata, is a bigger red flag than
#' the pooled p-value alone shows -- this figure surfaces both, by outlining
#' cell types whose statistic changes sign between the disease-specific
#' strata.
#'
#' @param strata_csvs Named character vector of paths, e.g.
#'   `c(Combined = "results/marker_concordance/marker_concordance.csv",
#'      "AD+CAA" = "results/marker_concordance/marker_concordance_AD_CAA.csv",
#'      Control = "results/marker_concordance/marker_concordance_Control.csv")`.
#'   Names become the stratum labels; entries whose file doesn't exist are
#'   dropped automatically.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Significance threshold.
#' @return The ggplot object (invisibly), or NULL if fewer than 2 strata exist.
#' @export
plot_marker_concordance_strata <- function(strata_csvs, output_dir, padj_cutoff = 0.05) {
  strata_csvs <- strata_csvs[file.exists(strata_csvs)]

  if (length(strata_csvs) < 2) {
    message("plot_marker_concordance_strata: need >=2 existing stratum files; skipping.")
    return(invisible(NULL))
  }

  dfs <- lapply(names(strata_csvs), function(nm) {
    d <- utils::read.csv(strata_csvs[[nm]], stringsAsFactors = FALSE)
    d$stratum <- nm
    d[, c("celltype", "statistic", "p_adj", "stratum")]
  })
  df <- do.call(rbind, dfs)

  df$sig_pos <- !is.na(df$statistic) & df$statistic > 0 &
    !is.na(df$p_adj) & df$p_adj < padj_cutoff

  # Flag cell types whose statistic changes sign across the *disease-specific*
  # strata (i.e. excluding the pooled "Combined" stratum, which will always
  # sit between them and isn't itself evidence of a flip).
  disease_strata <- setdiff(names(strata_csvs), "Combined")
  flip_lookup <- vapply(
    split(df$statistic[df$stratum %in% disease_strata],
          df$celltype[df$stratum %in% disease_strata]),
    function(vals) {
      vals <- vals[!is.na(vals)]
      length(vals) >= 2 && any(vals > 0) && any(vals < 0)
    },
    logical(1)
  )
  df$flip <- flip_lookup[df$celltype]
  df$flip[is.na(df$flip)] <- FALSE

  ref_stratum <- names(strata_csvs)[1]
  ref_df <- df[df$stratum == ref_stratum, ]
  celltype_order <- ref_df$celltype[order(ref_df$statistic)]
  celltype_order <- union(celltype_order, unique(df$celltype))
  df$celltype <- factor(df$celltype, levels = celltype_order)
  df$stratum <- factor(df$stratum, levels = names(strata_csvs))

  n_shapes <- length(strata_csvs)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = statistic, y = celltype)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey50") +
    ggplot2::geom_line(
      ggplot2::aes(group = celltype, color = flip),
      linewidth = 0.4, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = stratum, fill = sig_pos), size = 2.2, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "grey85", "TRUE" = "#B2182B"),
      labels = c("FALSE" = "consistent sign", "TRUE" = "sign disagrees across disease strata"),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "white", "TRUE" = "grey20"),
      labels = c("FALSE" = "n.s.", "TRUE" = paste0("pos., adj p < ", padj_cutoff)),
      name = NULL
    ) +
    ggplot2::scale_shape_manual(values = c(21, 22, 23, 24)[seq_len(n_shapes)], name = "Stratum") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21))) +
    ggplot2::labs(
      title = "Marker concordance across disease strata",
      subtitle = "Pooled vs. disease-specific replication -- red line = sign disagrees between disease strata",
      x = "Concordance correlation", y = NULL
    ) +
    theme_analysis()

  save_figure(p, "marker_concordance_strata", output_dir, width = 9, height = 10)
  invisible(p)
}


#' Factor-vs-marker confound heatmap (`check_marker_factor_mapping.R` output)
#'
#' Reads `factor_marker_correlation_matrix.csv` (inferred cell type x marker
#' score, Spearman rho) and `factor_marker_best_matches.csv` (which marker
#' score best matches each inferred cell type), and draws a heatmap with the
#' best-matching marker score per row outlined in black. Rows whose outline
#' isn't on the diagonal are a red flag: either the deconvolution mislabels
#' that cell type, or -- if many rows converge on the same column -- that
#' marker score is acting as a generic confound axis rather than
#' identity-specific signal (see `plot_factor_marker_best_match_bar()` for a
#' direct view of the latter).
#'
#' @param correlation_csv Path to `factor_marker_correlation_matrix.csv`.
#' @param best_matches_csv Path to `factor_marker_best_matches.csv`.
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly).
#' @export
plot_factor_marker_heatmap <- function(correlation_csv, best_matches_csv, output_dir) {
  mat <- utils::read.csv(correlation_csv, row.names = 1, check.names = FALSE)
  best <- utils::read.csv(best_matches_csv, stringsAsFactors = FALSE)

  df <- as.data.frame(mat)
  df$inferred_celltype <- rownames(df)

  long <- tidyr::pivot_longer(
    df,
    cols = -inferred_celltype,
    names_to = "marker_celltype",
    values_to = "rho"
  )

  long <- merge(
    long,
    best[, c("inferred_celltype", "best_marker_celltype")],
    by = "inferred_celltype", all.x = TRUE
  )
  long$is_best <- long$marker_celltype == long$best_marker_celltype

  ct_order <- best$inferred_celltype[order(best$best_rho)]
  long$inferred_celltype <- factor(long$inferred_celltype, levels = ct_order)
  long$marker_celltype <- factor(long$marker_celltype, levels = ct_order)

  lim <- max(abs(long$rho), na.rm = TRUE)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = marker_celltype, y = inferred_celltype, fill = rho)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.1) +
    ggplot2::geom_tile(
      data = long[long$is_best, , drop = FALSE],
      color = "black", linewidth = 0.6, fill = NA
    ) +
    scale_fill_diverging(name = "Spearman rho", limits = c(-lim, lim)) +
    ggplot2::labs(
      title = "Inferred cell-type abundance vs. marker-gene score, all pairs",
      subtitle = "Black outline = best-matching marker score per row (should sit on the diagonal)",
      x = "Marker-gene score (cell type)", y = "Inferred cell-type abundance"
    ) +
    theme_analysis(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      axis.text.y = ggplot2::element_text(size = 5)
    )

  save_figure(p, "factor_marker_confound_heatmap", output_dir, width = 12, height = 11)
  invisible(p)
}


#' Best-match marker-score frequency bar (confound-axis detector)
#'
#' Reads `factor_marker_best_matches.csv` and counts how often each marker
#' score "wins" as the best match across all inferred cell types. In a clean,
#' identity-specific deconvolution each marker score would win for roughly
#' one inferred cell type (its own); a marker score winning broadly for many
#' unrelated cell types means it's tracking a generic confound (e.g. total
#' RNA content) rather than that cell type's real biology.
#'
#' @param best_matches_csv Path to `factor_marker_best_matches.csv`.
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly).
#' @export
plot_factor_marker_best_match_bar <- function(best_matches_csv, output_dir) {
  best <- utils::read.csv(best_matches_csv, stringsAsFactors = FALSE)

  counts <- as.data.frame(table(best$best_marker_celltype), stringsAsFactors = FALSE)
  colnames(counts) <- c("marker_celltype", "n_wins")
  counts$is_self <- counts$marker_celltype %in%
    best$inferred_celltype[best$self_is_best]
  counts <- counts[order(counts$n_wins), ]
  counts$marker_celltype <- factor(counts$marker_celltype, levels = counts$marker_celltype)

  p <- ggplot2::ggplot(counts, ggplot2::aes(x = n_wins, y = marker_celltype, fill = is_self)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "#999999", "TRUE" = "#009E73"),
      labels = c("FALSE" = "wins for other cell type(s) only", "TRUE" = "also wins for itself"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "How many inferred cell types does each marker score \"win\" as best match?",
      subtitle = "A marker score winning broadly (bar >> 1) is acting as a generic confound axis, not identity signal",
      x = "Number of inferred cell types for which this is the best-matching marker score",
      y = NULL
    ) +
    theme_analysis(base_size = 8)

  save_figure(p, "factor_marker_best_match_frequency", output_dir, width = 8, height = 10)
  invisible(p)
}
