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
#' Reads the combined contrast summary (celltype, region, type,
#' log2_OR, p_adj) and draws a diverging heatmap of log2 odds ratios, with an
#' asterisk on cells passing an adjusted-p threshold.
#'
#' @param combined_csv Path to `combined_spatial_contrast_summary.csv`.
#' @param output_dir Figures output directory.
#' @param padj_cutoff Threshold for the significance asterisk.
#' @param cell_class One of `"all"`, `"neuronal"`, `"non_neuronal"` (see
#'   `classify_cell_class()`). Restricts the y-axis to that subset and writes
#'   to a suffixed filename, so the label list stays short enough to read.
#' @return The ggplot object (invisibly).
#' @export
plot_contrast_effsize_heatmap <- function(combined_csv, output_dir, padj_cutoff = 0.05,
                                           cell_class = c("all", "neuronal", "non_neuronal")) {
  cell_class <- match.arg(cell_class)
  df <- utils::read.csv(combined_csv, stringsAsFactors = FALSE)
  df$sig <- ifelse(!is.na(df$p_adj) & df$p_adj < padj_cutoff, "*", "")

  if (cell_class != "all") {
    df <- df[classify_cell_class(df$celltype) == cell_class, ]
  }
  n_celltypes <- length(unique(df$celltype))

  lim <- max(abs(df$log2_OR), na.rm = TRUE)

  class_label <- c(all = "all cell types", neuronal = "neuronal cell types",
                    non_neuronal = "non-neuronal cell types")[[cell_class]]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = region, y = celltype, fill = log2_OR)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = sig), color = "black", size = 3, vjust = 0.75) +
    ggplot2::facet_wrap(~ type, nrow = 1) +
    scale_fill_diverging(name = "log2 OR", limits = c(-lim, lim)) +
    ggplot2::labs(
      title = "Cell-type abundance effect sizes by region",
      subtitle = paste0("* adjusted p < ", padj_cutoff, "; diverging scale centered at 0; ", class_label),
      x = "Region", y = NULL
    ) +
    theme_analysis() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  filename <- if (cell_class == "all") {
    "contrast_effsize_heatmap"
  } else {
    paste0("contrast_effsize_heatmap_", cell_class)
  }
  save_figure(p, filename, output_dir, width = 11, height = max(4, 3 + n_celltypes * 0.15))
  invisible(p)
}


#' Forest plot of one contrast's estimates with confidence intervals
#'
#' Reads an individual contrast CSV from emmeans and plots log2 odds ratio
#' +/- CI per cell type, colored by region (<= 3 regions, so categorical
#' color is safe here). Handles both response-scale contrasts (`odds.ratio`,
#' as in the amyloid/disease effect files) and link-scale contrasts
#' (`estimate` in log-odds, as in the overall/max-pathology effect files),
#' and both asymptotic (`asymp.LCL`/`asymp.UCL`, the norm for glmmTMB fits
#' with `df = Inf`) and t-based (`lower.CL`/`upper.CL`) CI columns. If the
#' file contains more than one distinct `contrast` (e.g. the overall-effect
#' file, which stacks Disease_effect/Amyloid_effect/AD_overall_vs_Control),
#' facets by contrast so they don't overplot at identical dodge positions.
#'
#' @param contrast_csv Path to an individual `*_effect_contrasts.csv`.
#' @param contrast_label Human-readable label for the title/filename.
#' @param output_dir Figures output directory.
#' @param cell_class One of `"all"`, `"neuronal"`, `"non_neuronal"` (see
#'   `classify_cell_class()`). Restricts the y-axis to that subset and writes
#'   to a suffixed filename, so the label list stays short enough to read.
#' @return The ggplot object (invisibly), or NULL if value/CI columns are absent.
#' @export
plot_contrast_forest <- function(contrast_csv, contrast_label, output_dir,
                                  cell_class = c("all", "neuronal", "non_neuronal")) {
  cell_class <- match.arg(cell_class)
  df <- utils::read.csv(contrast_csv, stringsAsFactors = FALSE)

  if (cell_class != "all") {
    df <- df[classify_cell_class(df$celltype) == cell_class, ]
  }

  value_col <- intersect(c("odds.ratio", "estimate"), colnames(df))[1]
  ci_cols <- if (all(c("asymp.LCL", "asymp.UCL") %in% colnames(df))) {
    c("asymp.LCL", "asymp.UCL")
  } else if (all(c("lower.CL", "upper.CL") %in% colnames(df))) {
    c("lower.CL", "upper.CL")
  } else {
    NA_character_
  }

  if (is.na(value_col) || anyNA(ci_cols)) {
    message("plot_contrast_forest: no value/CI columns in ", contrast_csv, "; skipping.")
    return(invisible(NULL))
  }

  if (value_col == "odds.ratio") {
    df$log2_OR <- log2(df[[value_col]])
    df$log2_lower <- log2(df[[ci_cols[1]]])
    df$log2_upper <- log2(df[[ci_cols[2]]])
  } else {
    df$log2_OR <- df[[value_col]] / log(2)
    df$log2_lower <- df[[ci_cols[1]]] / log(2)
    df$log2_upper <- df[[ci_cols[2]]] / log(2)
  }

  regions <- sort(unique(df$region))
  n_contrasts <- dplyr::n_distinct(df$contrast)
  n_celltypes <- length(unique(df$celltype))

  class_label <- c(all = "all cell types", neuronal = "neuronal cell types",
                    non_neuronal = "non-neuronal cell types")[[cell_class]]
  title_suffix <- if (cell_class == "all") "" else paste0(" (", class_label, ")")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2_OR, y = celltype, color = region)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = log2_lower, xmax = log2_upper),
      height = 0, linewidth = 0.5,
      position = ggplot2::position_dodge(width = 0.6)
    ) +
    ggplot2::geom_point(size = 1.8, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_color_manual(values = okabe_ito_palette(length(regions))) +
    ggplot2::labs(
      title = paste0(contrast_label, ": log2 OR +/- 95% CI", title_suffix),
      x = "log2 odds ratio", y = NULL, color = "Region"
    ) +
    theme_analysis()

  if (n_contrasts > 1) {
    p <- p + ggplot2::facet_wrap(~ contrast)
  }

  filename <- if (cell_class == "all") {
    paste0("forest_", contrast_label)
  } else {
    paste0("forest_", contrast_label, "_", cell_class)
  }
  save_figure(
    p, filename, output_dir,
    width = 8 + 4 * (n_contrasts - 1), height = max(4, 3 + n_celltypes * 0.15)
  )
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


#' Cell-type co-occurrence correlation heatmap, faceted by stratum
#'
#' Facets by whichever stratify_by column(s) are present in the CSV (e.g.
#' just `region`, or `region` x `disease_status` x `pathology`) -- one panel
#' per unique combination.
#'
#' @param cooccurrence_csv Path to a `run_celltype_cooccurrence()` output
#'   csv, e.g. `celltype_cooccurrence_by_group.csv`.
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly).
#' @export
plot_cooccurrence_heatmap <- function(cooccurrence_csv, output_dir) {
  df <- utils::read.csv(cooccurrence_csv, stringsAsFactors = FALSE)

  strat_cols <- setdiff(colnames(df), c("celltype_1", "celltype_2", "r", "p.value", "p_adj", "n"))
  df$stratum <- do.call(paste, c(df[strat_cols], sep = " / "))

  # symmetrize so the heatmap is full, not triangular
  df_sym <- df
  swap <- df
  swap$celltype_1 <- df$celltype_2
  swap$celltype_2 <- df$celltype_1
  df_sym <- rbind(df_sym, swap)

  n_strata <- length(unique(df_sym$stratum))
  ncol_facets <- ceiling(sqrt(n_strata))
  nrow_facets <- ceiling(n_strata / ncol_facets)

  p <- ggplot2::ggplot(df_sym, ggplot2::aes(x = celltype_1, y = celltype_2, fill = r)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ stratum) +
    scale_fill_diverging(name = "CLR corr.", limits = c(-1, 1)) +
    ggplot2::labs(
      title = paste("Cell-type co-occurrence (CLR correlation) by", paste(strat_cols, collapse = " x ")),
      x = NULL, y = NULL
    ) +
    theme_analysis(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      axis.text.y = ggplot2::element_text(size = 5)
    )

  save_figure(
    p, "cooccurrence_heatmap", output_dir,
    width = max(14, 5 * ncol_facets), height = max(6, 4 * nrow_facets)
  )
  invisible(p)
}


#' Cell-type co-occurrence correlation heatmap, one full-size figure per stratum
#'
#' Companion to `plot_cooccurrence_heatmap()`. That figure facets every
#' stratum into one image, which is legible for a single `stratify_by`
#' column (e.g. `region`, 3 panels) but becomes cramped once
#' `celltype_cooccurrence_by_group.csv` stratifies by `region x
#' disease_status x pathology` -- with ~42 cell types per axis, a facet
#' panel that small isn't readable. This writes one full-size heatmap per
#' unique stratum combination instead, so each one stays legible.
#'
#' @param cooccurrence_csv Path to a `run_celltype_cooccurrence()` output
#'   csv, e.g. `celltype_cooccurrence_by_group.csv`.
#' @param output_dir Figures output directory. One file per stratum is
#'   written here, named `cooccurrence_heatmap_<stratum>`.
#' @return Invisibly, a named list of the ggplot objects (one per stratum).
#' @export
plot_cooccurrence_heatmap_by_stratum <- function(cooccurrence_csv, output_dir) {
  df <- utils::read.csv(cooccurrence_csv, stringsAsFactors = FALSE)

  strat_cols <- setdiff(colnames(df), c("celltype_1", "celltype_2", "r", "p.value", "p_adj", "n"))
  df$stratum <- do.call(paste, c(df[strat_cols], sep = " / "))

  # symmetrize so each heatmap is full, not triangular
  swap <- df
  swap$celltype_1 <- df$celltype_2
  swap$celltype_2 <- df$celltype_1
  df_sym <- rbind(df, swap)

  strata <- sort(unique(df_sym$stratum))

  plots <- lapply(strata, function(st) {
    sub <- df_sym[df_sym$stratum == st, , drop = FALSE]

    p <- ggplot2::ggplot(sub, ggplot2::aes(x = celltype_1, y = celltype_2, fill = r)) +
      ggplot2::geom_tile() +
      scale_fill_diverging(name = "CLR corr.", limits = c(-1, 1)) +
      ggplot2::labs(
        title = "Cell-type co-occurrence (CLR correlation)",
        subtitle = st,
        x = NULL, y = NULL
      ) +
      theme_analysis(base_size = 8) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = ggplot2::element_text(size = 6)
      )

    safe_name <- gsub("[^A-Za-z0-9]+", "_", st)
    safe_name <- gsub("^_+|_+$", "", safe_name)

    save_figure(p, paste0("cooccurrence_heatmap_", safe_name), output_dir, width = 11, height = 9)
    p
  })

  names(plots) <- strata
  invisible(plots)
}


#' Cell-type co-occurrence network of significant pairs (per stratum)
#'
#' Optional: requires `igraph` and `ggraph`. Draws one network per stratum
#' (whichever stratify_by column(s) are present in the CSV -- e.g. just
#' `region`, or `region` x `disease_status` x `pathology`), edges =
#' significant CLR correlations, edge color by sign. Filters on the raw
#' (unadjusted) p-value rather than the BH-adjusted `p_adj` -- per-stratum
#' sample sizes are small enough that adjusting across all pairs in a
#' stratum leaves few or no edges, and this plot is exploratory (the heatmap
#' is the figure that carries the adjusted p-values).
#'
#' @param cooccurrence_csv Path to a `run_celltype_cooccurrence()` output
#'   csv, e.g. `celltype_cooccurrence_by_group.csv`.
#' @param output_dir Figures output directory.
#' @param p_cutoff Edge significance threshold, applied to raw `p.value`.
#' @return The ggplot object (invisibly), or NULL if packages are unavailable.
#' @export
plot_cooccurrence_network <- function(cooccurrence_csv, output_dir, p_cutoff = 0.05) {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    message("plot_cooccurrence_network: igraph/ggraph not installed; skipping network ",
            "(the heatmap covers the same data).")
    return(invisible(NULL))
  }

  df <- utils::read.csv(cooccurrence_csv, stringsAsFactors = FALSE)
  df <- df[!is.na(df$p.value) & df$p.value < p_cutoff, ]

  if (nrow(df) == 0) {
    message("plot_cooccurrence_network: no significant pairs at p < ", p_cutoff, ".")
    return(invisible(NULL))
  }

  strat_cols <- setdiff(colnames(df), c("celltype_1", "celltype_2", "r", "p.value", "p_adj", "n"))
  df$stratum <- do.call(paste, c(df[strat_cols], sep = " / "))

  plots <- lapply(sort(unique(df$stratum)), function(st) {
    sub <- df[df$stratum == st, ]
    g <- igraph::graph_from_data_frame(
      sub[, c("celltype_1", "celltype_2", "r")],
      directed = FALSE
    )
    ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(ggplot2::aes(edge_color = r), edge_width = 0.6) +
      ggraph::geom_node_point(size = 2, color = "grey30") +
      ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 2.5) +
      ggraph::scale_edge_color_gradient2(
        low = "#2166AC", mid = "grey80", high = "#B2182B", midpoint = 0, name = "corr.",
        guide = ggraph::guide_edge_colourbar()
      ) +
      ggplot2::labs(title = st) +
      ggplot2::theme_void()
  })

  n_plots <- length(plots)
  ncol_plots <- min(n_plots, 3)
  nrow_plots <- ceiling(n_plots / ncol_plots)

  combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plots, ncol = ncol_plots)
  } else {
    plots[[1]]
  }

  save_figure(
    combined, "cooccurrence_network", output_dir,
    width = 5 * ncol_plots, height = 5 * nrow_plots
  )
  invisible(combined)
}


#' Cell-type co-occurrence network, one full-size figure per stratum
#'
#' Companion to `plot_cooccurrence_network()`. That figure tiles every
#' stratum's network into a single patchwork image, so each panel shrinks as
#' more strata are added -- once `celltype_cooccurrence_by_group.csv`
#' stratifies by `region x disease_status x pathology`, node labels in a
#' tiled panel are no longer legible. This writes one full-size network per
#' unique stratum combination instead. Same significance filter and
#' `igraph`/`ggraph` dependency as `plot_cooccurrence_network()`.
#'
#' @param cooccurrence_csv Path to a `run_celltype_cooccurrence()` output
#'   csv, e.g. `celltype_cooccurrence_by_group.csv`.
#' @param output_dir Figures output directory. One file per stratum is
#'   written here, named `cooccurrence_network_<stratum>`.
#' @param p_cutoff Edge significance threshold, applied to raw `p.value`.
#' @return Invisibly, a named list of the ggplot objects (one per stratum),
#'   or NULL if packages are unavailable or nothing is significant.
#' @export
plot_cooccurrence_network_by_stratum <- function(cooccurrence_csv, output_dir, p_cutoff = 0.05) {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    message("plot_cooccurrence_network_by_stratum: igraph/ggraph not installed; skipping ",
            "(the heatmap covers the same data).")
    return(invisible(NULL))
  }

  df <- utils::read.csv(cooccurrence_csv, stringsAsFactors = FALSE)
  df <- df[!is.na(df$p.value) & df$p.value < p_cutoff, ]

  if (nrow(df) == 0) {
    message("plot_cooccurrence_network_by_stratum: no significant pairs at p < ", p_cutoff, ".")
    return(invisible(NULL))
  }

  strat_cols <- setdiff(colnames(df), c("celltype_1", "celltype_2", "r", "p.value", "p_adj", "n"))
  df$stratum <- do.call(paste, c(df[strat_cols], sep = " / "))

  strata <- sort(unique(df$stratum))

  plots <- lapply(strata, function(st) {
    sub <- df[df$stratum == st, , drop = FALSE]
    g <- igraph::graph_from_data_frame(
      sub[, c("celltype_1", "celltype_2", "r")],
      directed = FALSE
    )

    p <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(ggplot2::aes(edge_color = r), edge_width = 0.6) +
      ggraph::geom_node_point(size = 2, color = "grey30") +
      ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3) +
      ggraph::scale_edge_color_gradient2(
        low = "#2166AC", mid = "grey80", high = "#B2182B", midpoint = 0, name = "corr.",
        guide = ggraph::guide_edge_colourbar()
      ) +
      ggplot2::labs(title = "Cell-type co-occurrence network", subtitle = st) +
      ggplot2::theme_void()

    safe_name <- gsub("[^A-Za-z0-9]+", "_", st)
    safe_name <- gsub("^_+|_+$", "", safe_name)

    save_figure(p, paste0("cooccurrence_network_", safe_name), output_dir, width = 8, height = 7)
    p
  })

  names(plots) <- strata
  invisible(plots)
}
