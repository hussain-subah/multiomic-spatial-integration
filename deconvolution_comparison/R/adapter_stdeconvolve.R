#' STdeconvolve (reference-free) deconvolution adapter
#'
#' STdeconvolve (Miller et al. 2022) is reference-FREE: it fits an LDA topic
#' model to the mixture alone, with no scRNA-seq reference. Its output is
#' therefore latent *topics*, not named cell types -- so the comparison here is
#' necessarily indirect: each topic's gene distribution (beta) is matched to
#' the reference signature matrix by correlation, topics are assigned to their
#' best-correlated cell type, and topic proportions are aggregated per cell
#' type. This assignment step is approximate; treat STdeconvolve's concordance
#' numbers with more caution than the reference-based methods'.
#'
#' @param mixture gene x ROI matrix from `load_mixture()`.
#' @param signature gene x celltype signature matrix (used only for the
#'   post-hoc topic->celltype assignment, not for the deconvolution itself).
#' @param k Number of topics. Defaults to the number of cell types in the
#'   signature; STdeconvolve normally recommends optimizing this.
#'
#' @return ROI x celltype proportion matrix, or NULL if STdeconvolve isn't
#'   installed.
#' @export
run_stdeconvolve <- function(mixture, signature, k = NULL) {
  if (!require_pkg("STdeconvolve", "STdeconvolve")) {
    return(NULL)
  }

  if (is.null(k)) k <- ncol(signature)

  # STdeconvolve works with a genes x spots counts matrix.
  counts <- STdeconvolve::cleanCounts(as.matrix(mixture), min.lib.size = 1)
  corpus <- STdeconvolve::restrictCorpus(
    counts,
    removeAbove = 1.0,
    removeBelow = 0.05,
    nTopOD = 1000
  )

  ldas <- STdeconvolve::fitLDA(t(as.matrix(corpus)), Ks = k)
  optLDA <- STdeconvolve::optimalModel(models = ldas, opt = k)
  results <- STdeconvolve::getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)

  theta <- results$theta   # ROI x topic (proportions)
  beta <- results$beta     # topic x gene

  # --- Post-hoc topic -> celltype assignment via signature correlation ---
  shared <- intersect(colnames(beta), rownames(signature))

  if (length(shared) < 50) {
    warning("run_stdeconvolve: only ", length(shared),
            " shared genes for topic->celltype assignment; results unreliable.",
            call. = FALSE)
  }

  beta_shared <- beta[, shared, drop = FALSE]
  sig_shared <- signature[shared, , drop = FALSE]

  # correlation of each topic's gene profile with each celltype's signature
  topic_ct_cor <- stats::cor(t(beta_shared), sig_shared)  # topic x celltype
  topic_assignment <- colnames(sig_shared)[apply(topic_ct_cor, 1, which.max)]

  celltypes <- colnames(signature)
  prop_mat <- matrix(
    0,
    nrow = nrow(theta),
    ncol = length(celltypes),
    dimnames = list(rownames(theta), celltypes)
  )

  for (topic_idx in seq_len(ncol(theta))) {
    ct <- topic_assignment[topic_idx]
    prop_mat[, ct] <- prop_mat[, ct] + theta[, topic_idx]
  }

  prop_mat
}
