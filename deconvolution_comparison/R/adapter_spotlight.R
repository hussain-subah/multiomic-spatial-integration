#' SPOTlight deconvolution adapter
#'
#' SPOTlight (Elosua-Bayes et al. 2021) seeds an NMF with a single-cell
#' reference and its marker genes, then NNLS-fits each mixture spot/ROI. Marker
#' genes are computed with scran::findMarkers on the reference (SPOTlight's
#' recommended default) unless supplied.
#'
#' @param ref list from `load_reference_mtx()`.
#' @param mixture gene x ROI matrix from `load_mixture()`.
#' @param n_hvg Number of top marker genes per cell type to pass to SPOTlight.
#'
#' @return ROI x celltype proportion matrix, or NULL if SPOTlight isn't installed.
#' @export
run_spotlight <- function(ref, mixture, n_hvg = 100) {
  if (!require_pkg("SPOTlight", "SPOTlight") ||
      !require_pkg("SingleCellExperiment", "SPOTlight") ||
      !require_pkg("scran", "SPOTlight")) {
    return(NULL)
  }

  aligned <- align_genes(mixture, rownames(ref$counts))
  mixture <- aligned$mixture

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = ref$counts),
    colData = S4Vectors::DataFrame(celltype = ref$cell_type)
  )
  # SPOTlight expects a logcounts assay for marker detection
  sce <- scran::computeSumFactors(sce, clusters = ref$cell_type)
  sce <- scuttle::logNormCounts(sce)

  mgs <- scran::findMarkers(sce, groups = ref$cell_type, direction = "up")

  mgs_df <- do.call(rbind, lapply(names(mgs), function(ct) {
    x <- as.data.frame(mgs[[ct]])
    x <- x[order(x$p.value), ]
    x <- utils::head(x, n_hvg)
    data.frame(
      gene = rownames(x),
      celltype = ct,
      weight = -log10(x$p.value + 1e-300),
      stringsAsFactors = FALSE
    )
  }))

  res <- SPOTlight::SPOTlight(
    x = sce,
    y = as.matrix(mixture),
    groups = as.character(ref$cell_type),
    mgs = mgs_df,
    group_id = "celltype",
    gene_id = "gene",
    weight_id = "weight"
  )

  # res$mat is ROI x celltype
  as.matrix(res$mat)
}
