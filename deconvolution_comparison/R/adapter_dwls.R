#' DWLS deconvolution adapter
#'
#' DWLS (Dampened Weighted Least Squares; Tsoucas et al. 2019) solves each
#' bulk sample against a cell-type signature matrix. Here we reuse the
#' already-exported Cell2location signature matrix as the signature, avoiding
#' DWLS's heavy MAST-based signature build, and solve each ROI independently.
#'
#' @param signature gene x celltype signature matrix from
#'   `load_signature_matrix()`.
#' @param mixture gene x ROI matrix from `load_mixture()`.
#'
#' @return ROI x celltype proportion matrix, or NULL if DWLS isn't installed.
#' @export
run_dwls <- function(signature, mixture) {
  if (!require_pkg("DWLS", "DWLS")) {
    return(NULL)
  }

  shared <- intersect(rownames(signature), rownames(mixture))

  if (length(shared) < 100) {
    warning("run_dwls: only ", length(shared), " shared genes between signature and mixture.",
            call. = FALSE)
  }

  signature <- signature[shared, , drop = FALSE]
  mixture <- mixture[shared, , drop = FALSE]

  roi_ids <- colnames(mixture)
  celltypes <- colnames(signature)

  prop_mat <- matrix(
    NA_real_,
    nrow = length(roi_ids),
    ncol = length(celltypes),
    dimnames = list(roi_ids, celltypes)
  )

  for (i in seq_along(roi_ids)) {
    bulk_sample <- mixture[, i]

    # trimData intersects and orders signature/bulk consistently
    tr <- DWLS::trimData(signature, bulk_sample)

    sol <- try(DWLS::solveDampenedWLS(tr$sig, tr$bulk), silent = TRUE)

    if (inherits(sol, "try-error")) {
      # fall back to OLS if the dampened solver fails to converge
      sol <- try(DWLS::solveOLS(tr$sig, tr$bulk), silent = TRUE)
      if (inherits(sol, "try-error")) next
    }

    prop_mat[i, names(sol)] <- sol
  }

  prop_mat[is.na(prop_mat)] <- 0
  prop_mat
}
