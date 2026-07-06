#' Common utilities for the deconvolution method comparison
#'
#' Shared loaders and the standard output contract used by every tool
#' adapter, so the downstream comparison step (comparison_utils.R) is
#' method-agnostic. See deconvolution_comparison/README.md for the input
#' files and output format.
#'
#' @keywords internal
NULL


#' Load the snRNA-seq reference from the MatrixMarket triple
#'
#' Reads the counts/metadata/genes exported by
#' `scripts/run_sn_reference_export.R` into the pieces every R adapter
#' needs, without depending on Seurat.
#'
#' @param counts_mtx MatrixMarket counts file (genes x cells after transpose;
#'   see `transpose`).
#' @param metadata_csv Cell metadata CSV (one row per cell).
#' @param genes_csv Gene annotation CSV.
#' @param celltype_col Metadata column with cell-type labels.
#' @param sample_col Metadata column identifying the biological sample/subject
#'   (used by MuSiC/Bisque for cross-sample variance). If not present, a single
#'   dummy sample is assigned.
#' @param gene_col Column in genes_csv holding gene ids.
#' @param transpose Whether to transpose the matrix after reading (the export
#'   writes cells x genes, so the default TRUE yields genes x cells).
#'
#' @return list(counts = genes x cells sparse Matrix, cell_type = factor,
#'   sample = factor, genes = character).
#' @export
load_reference_mtx <- function(counts_mtx,
                               metadata_csv,
                               genes_csv,
                               celltype_col = "New_Idents",
                               sample_col = "orig.ident",
                               gene_col = "gene_id",
                               transpose = TRUE) {
  counts <- Matrix::readMM(counts_mtx)
  if (transpose) counts <- Matrix::t(counts)
  counts <- methods::as(counts, "CsparseMatrix")

  meta <- utils::read.csv(metadata_csv, row.names = 1, stringsAsFactors = FALSE)
  genes <- utils::read.csv(genes_csv, stringsAsFactors = FALSE)

  gene_ids <- as.character(genes[[gene_col]])
  rownames(counts) <- gene_ids
  colnames(counts) <- rownames(meta)

  if (!celltype_col %in% colnames(meta)) {
    stop(
      "load_reference_mtx: celltype_col '", celltype_col,
      "' not in metadata. Available: ", paste(colnames(meta), collapse = ", "),
      call. = FALSE
    )
  }

  cell_type <- factor(meta[[celltype_col]])

  sample <- if (sample_col %in% colnames(meta)) {
    factor(meta[[sample_col]])
  } else {
    message("load_reference_mtx: sample_col '", sample_col,
            "' not found; assigning a single dummy sample.")
    factor(rep("sample1", ncol(counts)))
  }

  list(counts = counts, cell_type = cell_type, sample = sample, genes = gene_ids)
}


#' Load the GeoMx mixture (ROI x gene CSV) as a gene x ROI matrix
#'
#' Deconvolution tools expect genes x samples; the exported CSV is ROI x gene,
#' so it is transposed here.
#'
#' @param expression_csv ROI x gene CSV (see
#'   `python.geomx_anndata_utils.export_expression_csv`).
#'
#' @return Numeric matrix, rows = gene, columns = ROI.
#' @export
load_mixture <- function(expression_csv) {
  expr <- utils::read.csv(expression_csv, row.names = 1, check.names = FALSE)
  t(as.matrix(expr))
}


#' Load a gene x celltype signature matrix
#'
#' @param signature_csv Gene x celltype CSV (e.g. an inferred_signatures.csv).
#'
#' @return Numeric matrix, rows = gene, columns = celltype.
#' @export
load_signature_matrix <- function(signature_csv) {
  sig <- utils::read.csv(signature_csv, row.names = 1, check.names = FALSE)
  as.matrix(sig)
}


#' Build an independent per-cell-type signature from raw reference counts
#'
#' Aggregates the raw snRNA-seq reference counts (genes x cells) into a
#' genes x celltype mean-expression profile matrix. Unlike the pipeline's
#' `*_inferred_signatures.csv` (produced by cell2location's NB RegressionModel,
#' the same signatures the Cell2location/SpaceJam baseline is built on), this
#' signature never passes through cell2location -- so signature-based tools
#' (DWLS, SpatialDecon, CIBERSORTx) that use it are not sharing a basis with
#' the baseline they are compared against, removing that circularity.
#'
#' Each cell is normalized to a common library size before averaging, so the
#' profile isn't dominated by high-depth cells.
#'
#' @param ref list from `load_reference_mtx()` (`counts` = genes x cells,
#'   `cell_type` = factor).
#' @param scale_to Per-cell library-size normalization target (counts-per-N).
#' @return Numeric matrix, rows = gene, columns = cell type.
#' @export
build_mean_signature <- function(ref, scale_to = 1e4) {
  counts <- ref$counts

  lib <- Matrix::colSums(counts)
  lib[lib == 0] <- 1
  counts <- counts %*% Matrix::Diagonal(x = scale_to / lib)

  celltypes <- levels(ref$cell_type)

  sig <- vapply(
    celltypes,
    function(ct) {
      cols <- which(ref$cell_type == ct)
      as.numeric(Matrix::rowMeans(counts[, cols, drop = FALSE]))
    },
    numeric(nrow(counts))
  )

  rownames(sig) <- rownames(ref$counts)
  colnames(sig) <- celltypes
  sig
}


#' Restrict two gene-indexed objects to their shared genes
#'
#' @param mixture gene x ROI matrix.
#' @param reference_genes Character vector of reference/signature gene ids.
#'
#' @return list(mixture, shared_genes) with the mixture restricted to shared
#'   genes; a message reports the overlap size.
#' @export
align_genes <- function(mixture, reference_genes) {
  shared <- intersect(rownames(mixture), reference_genes)

  if (length(shared) < 100) {
    warning(
      "align_genes: only ", length(shared), " shared genes between mixture and ",
      "reference -- check that both use the same gene-id convention (symbol vs Ensembl).",
      call. = FALSE
    )
  } else {
    message("align_genes: ", length(shared), " shared genes.")
  }

  list(mixture = mixture[shared, , drop = FALSE], shared_genes = shared)
}


#' Convert a wide ROI x celltype proportion matrix to the standard long table
#'
#' @param prop_mat Matrix or data frame, rows = ROI (row names = ROI id),
#'   columns = celltype.
#' @param method Method name to stamp.
#' @param normalize Row-normalize to sum to 1 (most tools already do; a few
#'   return unnormalized weights).
#'
#' @return Long data frame: method, ROI_ID, celltype, proportion.
#' @export
standardize_proportions <- function(prop_mat, method, normalize = TRUE) {
  prop_mat <- as.matrix(prop_mat)

  if (normalize) {
    row_sums <- rowSums(prop_mat)
    row_sums[row_sums == 0] <- 1
    prop_mat <- prop_mat / row_sums
  }

  long <- data.frame(
    method = method,
    ROI_ID = rep(rownames(prop_mat), times = ncol(prop_mat)),
    celltype = rep(colnames(prop_mat), each = nrow(prop_mat)),
    proportion = as.vector(prop_mat),
    stringsAsFactors = FALSE
  )

  long
}


#' Write a standardized proportion table to the comparison output folder
#'
#' @param long_df Long data frame from `standardize_proportions()`.
#' @param method Method name (used for the filename).
#' @param output_dir Output directory.
#'
#' @return The written file path (invisibly).
#' @export
save_method_output <- function(long_df, method, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(output_dir, paste0(method, "_proportions.csv"))
  utils::write.csv(long_df, out_file, row.names = FALSE)
  message("Wrote ", out_file)
  invisible(out_file)
}


#' Load the Cell2location baseline into the standard long format
#'
#' @param proportions_csv The long cell-proportions CSV
#'   (`spatial_celltype_proportions_for_R.csv`).
#' @param roi_id_col ROI id column.
#' @param celltype_col Cell-type column.
#' @param proportion_col Proportion column.
#'
#' @return Long data frame: method, ROI_ID, celltype, proportion.
#' @export
load_cell2location_baseline <- function(proportions_csv,
                                        roi_id_col = "ROI_ID",
                                        celltype_col = "celltype",
                                        proportion_col = "rel_abundance") {
  df <- utils::read.csv(proportions_csv, stringsAsFactors = FALSE)

  data.frame(
    method = "Cell2location",
    ROI_ID = df[[roi_id_col]],
    celltype = df[[celltype_col]],
    proportion = df[[proportion_col]],
    stringsAsFactors = FALSE
  )
}


#' Load per-ROI negative-probe background as a named vector
#'
#' Used by the SpatialDecon adapter. Expects a CSV with an ROI id column and a
#' background value column (see
#' `python.geomx_anndata_utils.export_negprobes_background_csv`).
#'
#' @param background_csv Path to the background CSV.
#' @param background_col Column holding the per-ROI background value.
#'
#' @return Named numeric vector (names = ROI ids), or NULL if the file is absent.
#' @export
load_roi_background <- function(background_csv, background_col = "background") {
  if (!file.exists(background_csv)) {
    message("load_roi_background: '", background_csv, "' not found; ",
            "SpatialDecon will fall back to a proxy background.")
    return(NULL)
  }

  df <- utils::read.csv(background_csv, row.names = 1, stringsAsFactors = FALSE)
  stats::setNames(df[[background_col]], rownames(df))
}


#' Check a package is installed, with a clear skip message if not
#'
#' @param pkg Package name.
#' @param method Method label for the message.
#'
#' @return TRUE if available, FALSE (with a message) otherwise.
#' @export
require_pkg <- function(pkg, method) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(TRUE)
  }
  message("Skipping ", method, ": package '", pkg, "' is not installed.")
  FALSE
}
