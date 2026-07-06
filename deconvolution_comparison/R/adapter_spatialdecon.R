#' SpatialDecon (NanoString) deconvolution adapter
#'
#' SpatialDecon (Danaher et al. 2022) is NanoString's log-normal regression
#' deconvolution built specifically for GeoMx / DSP ROI data -- arguably the
#' most platform-appropriate reference-based method for this dataset. It takes
#' linear-scale *normalized* expression, a per-ROI background estimate (from
#' the negative control probes), and a cell-profile (signature) matrix.
#'
#' Normalization: Q3 (upper-quartile) scaling, the GeoMx-standard linear
#' normalization, is applied here since the shared mixture loader provides raw
#' counts.
#'
#' Background: pass `background` (a named per-ROI vector, e.g. the mean negative
#' -probe value exported alongside the expression matrix). If absent, a
#' 1st-percentile-per-ROI proxy is used with a warning -- fine for a rough run,
#' but the real negative-probe background is preferred for GeoMx.
#'
#' @param mixture gene x ROI raw-count matrix from `load_mixture()`.
#' @param signature gene x celltype signature (cell-profile) matrix.
#' @param background Optional named numeric vector, per-ROI background
#'   (names = ROI ids matching mixture columns).
#'
#' @return ROI x celltype proportion matrix, or NULL if SpatialDecon isn't
#'   installed.
#' @export
run_spatialdecon <- function(mixture, signature, background = NULL) {
  if (!require_pkg("SpatialDecon", "SpatialDecon")) {
    return(NULL)
  }

  aligned <- align_genes(mixture, rownames(signature))
  mixture <- aligned$mixture
  sig <- signature[aligned$shared_genes, , drop = FALSE]

  # --- Q3 (upper-quartile) normalization to linear-scale normalized expr ---
  q3 <- apply(mixture, 2, function(x) stats::quantile(x, 0.75, names = FALSE))
  q3[q3 == 0] <- 1
  norm <- sweep(mixture, 2, q3, "/") * mean(q3)

  # --- Per-ROI background, broadcast to a genes x ROI matrix ---
  if (is.null(background)) {
    per_roi_bg <- apply(norm, 2, function(x) stats::quantile(x, 0.01, names = FALSE))
    warning("run_spatialdecon: no negative-probe background supplied; using a ",
            "1st-percentile-per-ROI proxy. Supply the exported negative-probe ",
            "background for a proper GeoMx run.", call. = FALSE)
  } else {
    missing_bg <- setdiff(colnames(norm), names(background))
    if (length(missing_bg) > 0) {
      stop("run_spatialdecon: background is missing values for ", length(missing_bg),
           " ROI(s) present in the mixture (e.g. ", missing_bg[1], ").", call. = FALSE)
    }
    per_roi_bg <- background[colnames(norm)]
  }

  bg <- sweep(norm * 0, 2, per_roi_bg, "+")

  res <- SpatialDecon::spatialdecon(
    norm = norm,
    bg = bg,
    X = sig,
    align_genes = TRUE
  )

  # res$beta is celltype x ROI -> transpose to ROI x celltype
  t(as.matrix(res$beta))
}
