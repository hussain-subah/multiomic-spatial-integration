source("R/pathway_proportion_utils.R")
source("R/marker_concordance_utils.R")

expression_csv <- "results/geomx_exports/CAA-AD_expression_wide.csv"
proportions_csv <- "results/cell_proportions/roi_celltype_abundance_long.csv"
markers_file <- "results/AD_CAA_cluster_markers.txt"
output_dir <- "results/marker_concordance"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expr_raw <- load_roi_expression(expression_csv)
expr_norm <- normalize_expression_cpm(expr_raw)

prop <- read.csv(proportions_csv, stringsAsFactors = FALSE)

marker_lists <- load_marker_genes(
  markers_file,
  padj_cutoff = 0.05,
  top_n = 20
)

# Build marker-score matrix: ROI x marker-celltype
score_list <- lapply(names(marker_lists), function(ct) {
  genes <- marker_lists[[ct]]

  score <- try(
    compute_marker_score(
      expr_norm,
      genes,
      min_markers = 3
    ),
    silent = TRUE
  )

  if (inherits(score, "try-error")) {
    return(NULL)
  }

  score
})

names(score_list) <- names(marker_lists)
score_list <- score_list[!vapply(score_list, is.null, logical(1))]

score_mat <- do.call(cbind, score_list)
colnames(score_mat) <- names(score_list)


# Build abundance matrix: ROI x inferred cell type
abund_wide <- tidyr::pivot_wider(
  prop[, c("ROI_ID", "celltype", "rel_abundance")],
  id_cols = "ROI_ID",
  names_from = "celltype",
  values_from = "rel_abundance"
)

abund_wide <- as.data.frame(abund_wide)

rownames(abund_wide) <- as.character(abund_wide$ROI_ID)
abund_wide$ROI_ID <- NULL

abund_mat <- as.matrix(abund_wide)
storage.mode(abund_mat) <- "numeric"

cat("Score matrix dimensions:", dim(score_mat), "\n")
cat("Abundance matrix dimensions:", dim(abund_mat), "\n")

cat("First score ROI IDs:\n")
print(head(rownames(score_mat)))

cat("First abundance ROI IDs:\n")
print(head(rownames(abund_mat)))

shared_rois <- intersect(
  rownames(score_mat),
  rownames(abund_mat)
)

cat("Shared ROIs:", length(shared_rois), "\n")

if (length(shared_rois) == 0) {
  stop(
    "No shared ROI identifiers between marker-score and abundance matrices."
  )
}

score_mat <- score_mat[shared_rois, , drop = FALSE]
abund_mat <- abund_mat[shared_rois, , drop = FALSE]

cor_mat <- stats::cor(
  abund_mat,
  score_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

write.csv(
  cor_mat,
  file.path(output_dir, "factor_marker_correlation_matrix.csv")
)

# For each inferred factor, identify the marker score with strongest correlation
best_matches <- data.frame(
  inferred_celltype = rownames(cor_mat),
  best_marker_celltype = apply(
    cor_mat,
    1,
    function(x) colnames(cor_mat)[which.max(x)]
  ),
  best_rho = apply(cor_mat, 1, max, na.rm = TRUE),
  self_rho = vapply(
    rownames(cor_mat),
    function(ct) {
      if (ct %in% colnames(cor_mat)) cor_mat[ct, ct] else NA_real_
    },
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

best_matches$self_is_best <- (
  best_matches$inferred_celltype ==
    best_matches$best_marker_celltype
)

best_matches <- best_matches[
  order(best_matches$best_rho, decreasing = TRUE),
]

write.csv(
  best_matches,
  file.path(output_dir, "factor_marker_best_matches.csv"),
  row.names = FALSE
)

print(best_matches)
