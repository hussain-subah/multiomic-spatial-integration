#' RCTD / spacexr deconvolution adapter
#'
#' RCTD (Cable et al. 2021) is designed for spatial spot data and requires a
#' SpatialRNA object with per-spot coordinates. This GeoMx data has no
#' coordinates, so **dummy coordinates are supplied** -- RCTD's core
#' decomposition is per-ROI, so it still runs, but any spatial benefit is lost
#' and it acts here as a robust reference-based bulk-ish method. `doublet_mode
#' = "full"` is used because ROIs are mixtures of many cell types (not 1-2).
#'
#' @param ref list from `load_reference_mtx()`.
#' @param mixture gene x ROI matrix from `load_mixture()` (integer counts).
#' @param max_cores Cores for `run.RCTD`.
#'
#' @return ROI x celltype proportion matrix, or NULL if spacexr isn't installed.
#' @export
run_rctd <- function(ref, mixture, max_cores = 1) {
  if (!require_pkg("spacexr", "RCTD")) {
    return(NULL)
  }

  aligned <- align_genes(mixture, rownames(ref$counts))
  mixture <- aligned$mixture
  mixture <- round(mixture)  # RCTD expects integer counts
  storage.mode(mixture) <- "integer"

  cell_types <- ref$cell_type
  names(cell_types) <- colnames(ref$counts)
  ref_nUMI <- Matrix::colSums(ref$counts)

  reference <- spacexr::Reference(ref$counts, cell_types, ref_nUMI)

  # Dummy coordinates: RCTD requires them structurally but they carry no
  # spatial information for this platform.
  coords <- data.frame(
    x = seq_len(ncol(mixture)),
    y = rep(0, ncol(mixture)),
    row.names = colnames(mixture)
  )
  roi_nUMI <- colSums(mixture)

  puck <- spacexr::SpatialRNA(coords, mixture, roi_nUMI)

  myRCTD <- spacexr::create.RCTD(puck, reference, max_cores = max_cores)
  myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = "full")

  weights <- as.matrix(myRCTD@results$weights)  # ROI x celltype (unnormalized)
  weights
}
