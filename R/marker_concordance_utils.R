#' Marker-gene concordance utilities
#'
#' Pseudobulk validation (python/pseudobulk_utils.py) checks whether the
#' deconvolution model recovers a *known synthetic* composition. This is
#' the complementary real-data check: for each cell type, does its
#' inferred ROI proportion actually correlate with expression of its own
#' canonical marker genes in that same real ROI? A cell type whose
#' proportion doesn't track its own markers is a red flag independent of
#' any synthetic simulation.
#'
#' Requires `load_roi_expression()` / `normalize_expression_cpm()` from
#' `R/pathway_proportion_utils.R` to already be sourced.
#'
#' @keywords internal
NULL


#' Load per-cell-type marker gene lists from a Seurat-style marker table
#'
#' @param markers_file Path to a marker table (e.g.
#'   `results/AD_CAA_cluster_markers.txt`) with, at minimum, gene, cluster,
#'   and adjusted p-value columns (as produced by Seurat::FindAllMarkers()).
#' @param gene_col Column with gene names.
#' @param cluster_col Column identifying which cell type/cluster a marker
#'   belongs to (values should match the `celltype` column used elsewhere
#'   in the pipeline).
#' @param padj_col Adjusted p-value column.
#' @param padj_cutoff Significance threshold for including a marker.
#' @param top_n If set, keep only the top N markers per cell type (by
#'   ascending `padj_col`, or by `logfc_col` descending if provided).
#' @param logfc_col Optional log-fold-change column used to break ties /
#'   rank within `top_n`.
#'
#' @return Named list, cell type -> character vector of marker genes.
#' @export
load_marker_genes <- function(markers_file,
                              gene_col = "gene",
                              cluster_col = "cluster",
                              padj_col = "p_val_adj",
                              padj_cutoff = 0.05,
                              top_n = NULL,
                              logfc_col = "avg_log2FC") {

  markers <- utils::read.table(
    markers_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  required_cols <- c(gene_col, cluster_col, padj_col)
  missing_cols <- setdiff(required_cols, colnames(markers))

  if (length(missing_cols) > 0) {
    stop(
      "load_marker_genes: missing column(s) ",
      paste(missing_cols, collapse = ", "),
      " in ", markers_file,
      ". Available columns: ",
      paste(colnames(markers), collapse = ", "),
      call. = FALSE
    )
  }

  markers <- markers[
    !is.na(markers[[padj_col]]) &
      markers[[padj_col]] < padj_cutoff,
  ]

  if (logfc_col %in% colnames(markers)) {
    markers <- markers[
      !is.na(markers[[logfc_col]]) &
        markers[[logfc_col]] > 0,
    ]
  }

  if (!is.null(top_n) && logfc_col %in% colnames(markers)) {
    markers <- markers[
      order(
        markers[[cluster_col]],
        -markers[[logfc_col]],
        markers[[padj_col]]
      ),
    ]
  } else {
    markers <- markers[
      order(
        markers[[cluster_col]],
        markers[[padj_col]]
      ),
    ]
  }

  split_markers <- split(
    markers[[gene_col]],
    markers[[cluster_col]]
  )

  if (!is.null(top_n)) {
    split_markers <- lapply(
      split_markers,
      utils::head,
      n = top_n
    )
  }

  lapply(split_markers, unique)
}


#' Compute a per-ROI marker score
#'
#' Mean normalized expression across a gene set, per ROI. Genes not found
#' in `expr_mat` are silently dropped (a warning is issued if fewer than
#' half the requested markers are found).
#'
#' @param expr_mat Numeric matrix, rows = ROI, columns = gene (typically
#'   already normalized via `normalize_expression_cpm()`).
#' @param genes Character vector of marker genes.
#'
#' @return Named numeric vector, one score per ROI.
#' @export
compute_marker_score <- function(expr_mat, genes, min_markers = 3) {
  found_genes <- intersect(genes, colnames(expr_mat))

  if (length(found_genes) < min_markers) {
    stop(
      "compute_marker_score: only ", length(found_genes),
      " usable marker gene(s); need at least ", min_markers, ".",
      call. = FALSE
    )
  }

  marker_mat <- as.matrix(expr_mat[, found_genes, drop = FALSE])

  # Standardize each marker across ROIs.
  marker_z <- scale(marker_mat)

  # Remove genes with zero variance, which become NA after scaling.
  keep <- colSums(is.finite(marker_z)) == nrow(marker_z)
  marker_z <- marker_z[, keep, drop = FALSE]

  if (ncol(marker_z) < min_markers) {
    stop(
      "compute_marker_score: fewer than ", min_markers,
      " variable marker genes remained after scaling.",
      call. = FALSE
    )
  }

  stats::setNames(
    rowMeans(marker_z),
    rownames(expr_mat)
  )
}


#' Evaluate marker-gene concordance for each cell type
#'
#' For each cell type with both an inferred proportion and a marker gene
#' list, correlates the per-ROI marker score against the per-ROI inferred
#' proportion. A cell type whose proportion doesn't positively correlate
#' with its own markers is a real-data red flag for that cell type's
#' deconvolution quality, independent of the synthetic pseudobulk check.
#'
#' @param expr_mat Numeric matrix, rows = ROI, columns = gene (normalized).
#' @param proportions_df Long-format cell-type proportion table.
#' @param marker_genes_list Named list from `load_marker_genes()`.
#' @param roi_id_col Column in `proportions_df` identifying the ROI.
#' @param proportion_col Column with the relative abundance value.
#' @param celltype_col Column with cell-type labels.
#' @param method Correlation method passed to `stats::cor()`.
#' @param min_shared_rois Minimum overlapping ROIs required.
#'
#' @return Data frame: celltype, n_markers_requested, n_markers_found,
#'   statistic, p.value, p_adj, n.
#' @export
evaluate_marker_concordance <- function(expr_mat,
                                        proportions_df,
                                        marker_genes_list,
                                        roi_id_col = "ROI_ID",
                                        proportion_col = "rel_abundance",
                                        celltype_col = "celltype",
                                        method = "spearman",
                                        min_shared_rois = 5,
                                        min_markers = 3) {
  if (!roi_id_col %in% colnames(proportions_df)) {
    stop(
      "evaluate_marker_concordance: roi_id_col '", roi_id_col,
      "' not found in proportions_df. Available columns: ",
      paste(colnames(proportions_df), collapse = ", "),
      call. = FALSE
    )
  }

  celltypes <- intersect(names(marker_genes_list), unique(proportions_df[[celltype_col]]))

  if (length(celltypes) == 0) {
    stop(
      "evaluate_marker_concordance: no overlap between marker_genes_list names ",
      "and proportions_df[[celltype_col]] values.",
      call. = FALSE
    )
  }

  results <- lapply(celltypes, function(ct) {
    genes <- marker_genes_list[[ct]]

    score <- try(
      compute_marker_score(
        expr_mat,
        genes,
        min_markers = min_markers
      ),
      silent = TRUE
   )

    if (inherits(score, "try-error")) {
  return(data.frame(
    celltype = ct,
    n_markers_requested = length(genes),
    n_markers_found = sum(genes %in% colnames(expr_mat)),
    statistic = NA_real_,
    p.value = NA_real_,
    n = NA_integer_,
    stringsAsFactors = FALSE
  ))
}

    dat_ct <- proportions_df[proportions_df[[celltype_col]] == ct, ]
    rownames(dat_ct) <- dat_ct[[roi_id_col]]

    shared_rois <- intersect(names(score), dat_ct[[roi_id_col]])

    if (length(shared_rois) < min_shared_rois) {
      message(
        "evaluate_marker_concordance: skipping '", ct, "' -- only ",
        length(shared_rois), " shared ROI(s) (need >= ", min_shared_rois, ")."
      )
      return(data.frame(
        celltype = ct,
        n_markers_requested = length(genes),
        n_markers_found = sum(genes %in% colnames(expr_mat)),
        statistic = NA_real_,
        p.value = NA_real_,
        n = length(shared_rois),
        stringsAsFactors = FALSE
      ))
    }

    score_sub <- score[shared_rois]
    proportion_sub <- dat_ct[shared_rois, proportion_col]

    if (method == "spearman") {
      test <- stats::cor.test(
        score_sub,
        proportion_sub,
        method = "spearman",
        exact = FALSE
      )
    } else {
      test <- stats::cor.test(
        score_sub,
        proportion_sub,
        method = method
      )
    }

    data.frame(
      celltype = ct,
      n_markers_requested = length(genes),
      n_markers_found = sum(genes %in% colnames(expr_mat)),
      statistic = unname(test$estimate),
      p.value = test$p.value,
      n = length(shared_rois),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  out$p_adj <- stats::p.adjust(out$p.value, method = "BH")

  out[order(out$statistic, decreasing = TRUE, na.last = TRUE), ]
}
