#' MuSiC deconvolution adapter
#'
#' MuSiC (Wang et al. 2019) is a bulk deconvolution method that weights genes
#' by cross-subject consistency, using an scRNA-seq reference. A natural fit
#' for bulk-like GeoMx ROI data.
#'
#' Reference-container note: MuSiC wants a SingleCellExperiment with a raw
#' `counts` assay plus cell-type and subject/sample columns in colData.
#'
#' @param ref list from `load_reference_mtx()`.
#' @param mixture gene x ROI matrix from `load_mixture()` (raw-ish counts).
#'
#' @return ROI x celltype proportion matrix, or NULL if MuSiC isn't installed.
#' @export
run_music <- function(ref, mixture) {
  if (!require_pkg("MuSiC", "MuSiC") ||
      !require_pkg("SingleCellExperiment", "MuSiC")) {
    return(NULL)
  }

  aligned <- align_genes(mixture, rownames(ref$counts))
  mixture <- aligned$mixture

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = ref$counts),
    colData = S4Vectors::DataFrame(
      celltype = ref$cell_type,
      sample = ref$sample
    )
  )

  est <- MuSiC::music_prop(
    bulk.mtx = mixture,
    sc.sce = sce,
    clusters = "celltype",
    samples = "sample"
  )

  # Est.prop.weighted is ROI x celltype
  as.matrix(est$Est.prop.weighted)
}
