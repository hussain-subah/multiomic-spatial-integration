#' BayesPrism deconvolution adapter
#'
#' BayesPrism (Chu et al. 2022) jointly models the reference and mixture in a
#' Bayesian framework. A strong methodological comparison to Cell2location
#' (both Bayesian, scRNA-seq reference), and appropriate for bulk-like mixtures.
#'
#' Orientation note: BayesPrism expects cells x genes for the reference and
#' samples x genes for the mixture -- both transposed from this pipeline's
#' usual genes-in-rows convention.
#'
#' @param ref list from `load_reference_mtx()`.
#' @param mixture gene x ROI matrix from `load_mixture()`.
#' @param n_cores Cores for `run.prism`.
#'
#' @return ROI x celltype proportion matrix, or NULL if BayesPrism isn't installed.
#' @export
run_bayesprism <- function(ref, mixture, n_cores = 1) {
  if (!require_pkg("BayesPrism", "BayesPrism")) {
    return(NULL)
  }

  aligned <- align_genes(mixture, rownames(ref$counts))
  mixture <- aligned$mixture
  shared <- aligned$shared_genes

  # BayesPrism wants cells x genes and samples x genes
  ref_mat <- as.matrix(Matrix::t(ref$counts[shared, , drop = FALSE]))
  mix_mat <- t(mixture)

  cell_type_labels <- as.character(ref$cell_type)

  myPrism <- BayesPrism::new.prism(
    reference = ref_mat,
    mixture = mix_mat,
    input.type = "count.matrix",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_type_labels,
    key = NULL
  )

  bp_res <- BayesPrism::run.prism(prism = myPrism, n.cores = n_cores)

  theta <- BayesPrism::get.fraction(
    bp = bp_res,
    which.theta = "final",
    state.or.type = "type"
  )

  # theta is ROI x celltype
  as.matrix(theta)
}
