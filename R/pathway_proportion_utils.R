#' Pathway-linkage utilities for cell-type proportions
#'
#' Cell-type proportions are compositional (they sum to ~1 within an ROI),
#' so pathway/GSEA tools can't be run on them directly. These utilities
#' instead correlate a cell type's per-ROI proportion with per-ROI gene
#' expression, rank genes by that association, and run rank-based gene set
#' enrichment (fgsea) on the ranked list -- i.e. "what pathways covary with
#' this cell type's abundance shifts."
#'
#' @keywords internal
NULL


#' Load a wide ROI x gene expression matrix
#'
#' @param expression_csv Path to a CSV with ROI ids as row names and genes
#'   as columns (see `python.geomx_anndata_utils.export_expression_csv()`).
#'
#' @return Numeric matrix, rows = ROI, columns = gene.
#' @export
load_roi_expression <- function(expression_csv) {
  expr_df <- read.csv(expression_csv, row.names = 1, check.names = FALSE)
  as.matrix(expr_df)
}


#' Per-ROI log2-CPM normalization
#'
#' Normalizes each ROI (row) by its own total count before log-transforming,
#' so correlation with cell-type proportion isn't confounded by sequencing
#' depth.
#'
#' @param expr_mat Numeric matrix, rows = ROI, columns = gene (raw counts).
#' @param pseudocount Added before taking log2.
#'
#' @return Numeric matrix, same shape as `expr_mat`.
#' @export
normalize_expression_cpm <- function(expr_mat, pseudocount = 1) {
  lib_sizes <- rowSums(expr_mat)

  if (any(lib_sizes == 0)) {
    stop(
      "normalize_expression_cpm: ", sum(lib_sizes == 0),
      " ROI(s) have zero total counts; remove them before normalizing.",
      call. = FALSE
    )
  }

  cpm <- sweep(expr_mat, 1, lib_sizes, "/") * 1e6

  log2(cpm + pseudocount)
}


#' Rank genes by association with a cell type's proportion
#'
#' Joins `proportions_df` (filtered to one cell type) to `expr_mat` by ROI
#' id, then computes a single vectorized correlation of every gene column
#' against the proportion vector (not a per-gene loop), converting
#' coefficients to p-values via the standard t-approximation
#' (`t = r * sqrt(n-2) / sqrt(1-r^2)`), BH-adjusted across genes.
#'
#' @param expr_mat Numeric matrix, rows = ROI, columns = gene (typically
#'   already normalized via `normalize_expression_cpm()`).
#' @param proportions_df Long-format cell-type proportion table (e.g.
#'   `results/cell_proportions/spatial_celltype_proportions_for_R.csv`).
#' @param celltype Cell type to test (must appear in `celltype_col`).
#' @param roi_id_col Column in `proportions_df` identifying the ROI,
#'   matching `rownames(expr_mat)`.
#' @param proportion_col Column with the relative abundance value.
#' @param celltype_col Column with cell-type labels.
#' @param method Correlation method passed to `stats::cor()`.
#' @param min_shared_rois Minimum overlapping ROIs required between
#'   `expr_mat` and `proportions_df` before proceeding.
#'
#' @return Data frame `gene, statistic, p.value, p_adj`, sorted by
#'   `statistic` descending.
#' @export
rank_genes_by_celltype_association <- function(expr_mat,
                                               proportions_df,
                                               celltype,
                                               roi_id_col = "ROI_ID",
                                               proportion_col = "rel_abundance",
                                               celltype_col = "celltype",
                                               method = "spearman",
                                               min_shared_rois = 5) {
  if (!roi_id_col %in% colnames(proportions_df)) {
    stop(
      "rank_genes_by_celltype_association: roi_id_col '", roi_id_col,
      "' not found in proportions_df. Available columns: ",
      paste(colnames(proportions_df), collapse = ", "),
      call. = FALSE
    )
  }

  dat_ct <- proportions_df[proportions_df[[celltype_col]] == celltype, ]

  if (nrow(dat_ct) == 0) {
    stop(
      "rank_genes_by_celltype_association: no rows found for celltype '",
      celltype, "' in column '", celltype_col, "'.",
      call. = FALSE
    )
  }

  shared_rois <- intersect(rownames(expr_mat), dat_ct[[roi_id_col]])

  if (length(shared_rois) < min_shared_rois) {
    stop(
      "rank_genes_by_celltype_association: only ", length(shared_rois),
      " shared ROI(s) between expr_mat (", nrow(expr_mat),
      " rows) and proportions_df (", nrow(dat_ct), " rows for '", celltype,
      "') via roi_id_col = '", roi_id_col, "' -- check that expr_mat's ",
      "row names and proportions_df[[roi_id_col]] use the same ROI ",
      "identifier convention.",
      call. = FALSE
    )
  }

  rownames(dat_ct) <- dat_ct[[roi_id_col]]
  dat_ct <- dat_ct[shared_rois, ]
  expr_sub <- expr_mat[shared_rois, , drop = FALSE]

  proportion_vec <- dat_ct[[proportion_col]]

  r <- as.numeric(stats::cor(expr_sub, proportion_vec, method = method))
  n <- length(proportion_vec)

  t_stat <- r * sqrt(n - 2) / sqrt(1 - r^2)
  p_value <- 2 * stats::pt(-abs(t_stat), df = n - 2)

  result <- data.frame(
    gene = colnames(expr_sub),
    statistic = r,
    p.value = p_value,
    stringsAsFactors = FALSE
  )

  result$p_adj <- stats::p.adjust(result$p.value, method = "BH")

  result[order(-result$statistic), ]
}


#' Default pattern for generic/housekeeping Reactome pathways
#'
#' Translation, ribosome, spliceosome, NMD, oxidative phosphorylation, and
#' the oversized "Olfactory Receptor" gene family recur as GSEA hits for
#' almost any broad transcriptional-activity signal, regardless of the cell
#' type driving it -- flagged here so they can be excluded when looking for
#' cell-type-specific signal.
#'
#' @keywords internal
.housekeeping_pathway_pattern <- paste(
  "Respiratory", "Electron Transport", "Ribosom", "rRNA", "Translation",
  "Spliceosome", "Processing Of Capped", "mRNA Splicing",
  "Nonsense.Mediated", "SRP.dependent", "Citric Acid", "Olfactory",
  "Selenocysteine", "Selenoamino", "Influenza", "Complex I Biogenesis",
  "Cristae Formation",
  sep = "|"
)


#' Load and combine every `pathway_enrichment.csv` under a pathway-linkage dir
#'
#' Recurses the `<stratum>/<celltype>/pathway_enrichment.csv` directory
#' layout written by `scripts/run_pathway_proportion_link.R` and stacks every
#' file into one long data frame, tagging each row with the stratum and cell
#' type its enrichment came from.
#'
#' @param pathway_dir Root directory (e.g. `results/pathway_proportion_link`).
#' @return Data frame: `stratum, celltype, pathway, pval, padj, NES, size`.
#' @export
load_pathway_enrichment_dir <- function(pathway_dir) {
  files <- list.files(
    pathway_dir, pattern = "pathway_enrichment\\.csv$",
    recursive = TRUE, full.names = TRUE
  )

  if (length(files) == 0) {
    stop(
      "load_pathway_enrichment_dir: no pathway_enrichment.csv files found ",
      "under '", pathway_dir, "'.", call. = FALSE
    )
  }

  read_one <- function(f) {
    rel <- sub(paste0("^", pathway_dir, "/?"), "", f)
    parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
    df <- utils::read.csv(f, stringsAsFactors = FALSE)
    if (nrow(df) == 0) return(NULL)
    df$stratum <- parts[1]
    df$celltype <- parts[2]
    df[, c("stratum", "celltype", "pathway", "pval", "padj", "NES", "size")]
  }

  do.call(rbind, lapply(files, read_one))
}


#' Flag generic/housekeeping pathways and count cross-combination frequency
#'
#' For each pathway, counts how many `(celltype, stratum)` combinations it is
#' significant in. A pathway significant in nearly every combination is very
#' likely tracking a shared technical or compositional confound (sequencing
#' depth, overall transcriptional activity, or the closure constraint that
#' compositional proportions sum to ~1) rather than being independent
#' cell-type-specific biology; a pathway significant in only a handful is a
#' much more credible cell-type-restricted hit.
#'
#' @param combined_df Output of `load_pathway_enrichment_dir()`.
#' @param padj_cutoff Significance threshold.
#' @param housekeeping_pattern Regex (OR-joined) flagging generic pathways.
#' @return List with:
#'   - `frequency`: one row per unique pathway (`pathway`, `n_significant`,
#'     `is_housekeeping`), sorted by `n_significant` descending.
#'   - `significant`: the significant subset of `combined_df`, merged with
#'     `n_significant`/`is_housekeeping`.
#' @export
compute_pathway_specificity <- function(combined_df, padj_cutoff = 0.05,
                                        housekeeping_pattern = .housekeeping_pathway_pattern) {
  sig <- combined_df[!is.na(combined_df$padj) & combined_df$padj < padj_cutoff, ]

  freq <- as.data.frame(table(sig$pathway), stringsAsFactors = FALSE)
  colnames(freq) <- c("pathway", "n_significant")
  freq$is_housekeeping <- grepl(housekeeping_pattern, freq$pathway, ignore.case = TRUE)
  freq <- freq[order(-freq$n_significant), ]

  sig <- merge(sig, freq, by = "pathway")

  list(frequency = freq, significant = sig)
}


#' Top cell-type-restricted, non-housekeeping pathways, per stratum
#'
#' Filters `compute_pathway_specificity()`'s significant hits down to
#' pathways that are (a) not flagged as generic/housekeeping and (b)
#' significant in at most `max_frequency` of the total `(celltype, stratum)`
#' combinations, then keeps the top `top_n` per stratum by adjusted p-value.
#'
#' @param specificity_result Output of `compute_pathway_specificity()`.
#' @param max_frequency Maximum global `n_significant` count to keep (i.e.
#'   the specificity threshold).
#' @param top_n Number of rows to keep per stratum.
#' @return Data frame sorted by `stratum`, then `padj` ascending.
#' @export
top_specific_pathways_by_stratum <- function(specificity_result, max_frequency = 10, top_n = 15) {
  sig <- specificity_result$significant
  candidates <- sig[!sig$is_housekeeping & sig$n_significant <= max_frequency, ]

  do.call(rbind, lapply(split(candidates, candidates$stratum), function(df) {
    df <- df[order(df$padj), ]
    utils::head(df, top_n)
  }))
}


#' Run rank-based gene set enrichment on a ranked gene list
#'
#' @param ranked_genes Data frame from `rank_genes_by_celltype_association()`
#'   (`gene`, `statistic` columns used).
#' @param gmt_file Path to a GMT gene-set file (e.g. MSigDB Hallmark/KEGG/
#'   Reactome/GO). If `NULL`, enrichment is skipped (the ranked gene list is
#'   still useful on its own, e.g. via manual Enrichr upload).
#' @param min_size Minimum pathway size passed to `fgsea::fgsea()`.
#' @param max_size Maximum pathway size passed to `fgsea::fgsea()`.
#'
#' @return Data frame of fgsea results sorted by adjusted p-value, or `NULL`
#'   if `gmt_file` is `NULL`.
#' @export
run_pathway_enrichment <- function(ranked_genes,
                                   gmt_file,
                                   min_size = 15,
                                   max_size = 500) {
  if (is.null(gmt_file)) {
    message("run_pathway_enrichment: gmt_file is NULL, skipping enrichment.")
    return(NULL)
  }

  pathways <- fgsea::gmtPathways(gmt_file)

  ranks <- stats::setNames(ranked_genes$statistic, ranked_genes$gene)
  ranks <- sort(ranks, decreasing = TRUE)

  res <- fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = min_size,
    maxSize = max_size
  )

  res <- as.data.frame(res)
  res[order(res$padj), ]
}
