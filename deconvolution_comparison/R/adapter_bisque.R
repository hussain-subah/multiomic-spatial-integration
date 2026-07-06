#' Bisque deconvolution adapter
#'
#' Bisque (Jew et al. 2020) reference-based decomposition. Both the bulk
#' mixture and the scRNA-seq reference are supplied as Biobase ExpressionSets.
#' `use.overlap = FALSE` because the GeoMx ROIs and the snRNA-seq cells are not
#' matched samples.
#'
#' @param ref list from `load_reference_mtx()`.
#' @param mixture gene x ROI matrix from `load_mixture()`.
#'
#' @return ROI x celltype proportion matrix, or NULL if Bisque isn't installed.
#' @export
run_bisque <- function(ref, mixture) {
  if (!require_pkg("BisqueRNA", "Bisque") ||
      !require_pkg("Biobase", "Bisque")) {
    return(NULL)
  }

  aligned <- align_genes(mixture, rownames(ref$counts))
  mixture <- aligned$mixture

  bulk_eset <- Biobase::ExpressionSet(assayData = as.matrix(mixture))

  sc_pheno <- data.frame(
    cellType = ref$cell_type,
    SubjectName = ref$sample,
    row.names = colnames(ref$counts)
  )
  sc_meta <- data.frame(
    labelDescription = c("cellType", "SubjectName"),
    row.names = c("cellType", "SubjectName")
  )
  sc_eset <- Biobase::ExpressionSet(
    assayData = as.matrix(ref$counts),
    phenoData = methods::new("AnnotatedDataFrame", data = sc_pheno, varMetadata = sc_meta)
  )

  res <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset = bulk_eset,
    sc.eset = sc_eset,
    cell.types = "cellType",
    subject.names = "SubjectName",
    use.overlap = FALSE
  )

  # res$bulk.props is celltype x ROI -> transpose to ROI x celltype
  t(as.matrix(res$bulk.props))
}
