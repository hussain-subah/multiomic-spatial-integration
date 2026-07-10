# ============================================================
# Marker-gene concordance check
#
# Real-data complement to the synthetic pseudobulk validation
# (python/pseudobulk_utils.py): for each cell type, does its inferred
# per-ROI proportion correlate with expression of its own canonical
# marker genes in that same real ROI? A cell type whose proportion
# doesn't track its own markers is a red flag independent of any
# synthetic simulation.
# ============================================================

source("R/pathway_proportion_utils.R")
source("R/marker_concordance_utils.R")

# ============================================================
# Config -- adjust paths for your environment
# ============================================================

expression_csv <- "results/geomx_exports/CAA-AD_expression_wide.csv"
proportions_csv <- "results/cell_proportions/roi_celltype_abundance_long.csv"
markers_file <- "results/AD_CAA_cluster_markers.txt"
output_dir <- "results/marker_concordance"

roi_id_col <- "ROI_ID"
proportion_col <- "rel_abundance"
celltype_col <- "celltype"
method <- "spearman"

# Top N markers per cell type to use for the score (NULL = all significant
# markers at padj < 0.05).
top_n_markers <- 20

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================

expr_raw <- load_roi_expression(expression_csv)
expr_norm <- normalize_expression_cpm(expr_raw)

proportions_df <- read.csv(proportions_csv, stringsAsFactors = FALSE)

marker_genes_list <- load_marker_genes(
  markers_file,
  padj_cutoff = 0.05,
  top_n = top_n_markers
)

message(sprintf(
  "Loaded marker gene lists for %d cell type(s).",
  length(marker_genes_list)
))

# ============================================================
# Evaluate concordance
# ============================================================

# ============================================================
# Helper: run concordance on one subset and save results
# ============================================================

run_concordance_subset <- function(prop_subset,
                                   analysis_name,
                                   output_dir) {

  message(sprintf(
    "Running marker concordance for '%s' using %d ROI-celltype rows.",
    analysis_name,
    nrow(prop_subset)
  ))

  result <- evaluate_marker_concordance(
    expr_mat = expr_norm,
    proportions_df = prop_subset,
    marker_genes_list = marker_genes_list,
    roi_id_col = roi_id_col,
    proportion_col = proportion_col,
    celltype_col = celltype_col,
    method = method,
    min_shared_rois = 5,
    min_markers = 3
  )

  output_file <- file.path(
    output_dir,
    paste0("marker_concordance_", analysis_name, ".csv")
  )

  write.csv(
    result,
    file = output_file,
    row.names = FALSE
  )

  evaluable <- !is.na(result$statistic)

  passing <- (
    evaluable &
      result$statistic > 0 &
      !is.na(result$p_adj) &
      result$p_adj < 0.05
  )

  flagged <- result[
    !evaluable |
      result$statistic <= 0 |
      is.na(result$p_adj) |
      result$p_adj >= 0.05,
  ]

  message(sprintf(
    "%s: %d of %d evaluable cell type(s) show positive, significant concordance.",
    analysis_name,
    sum(passing),
    sum(evaluable)
  ))

  if (nrow(flagged) > 0) {
    message(sprintf(
      "%d cell type(s) were non-significant, non-positive, or unevaluable:",
      nrow(flagged)
    ))

    print(
      flagged[
        ,
        c(
          "celltype",
          "n_markers_found",
          "statistic",
          "p_adj",
          "n"
        )
      ]
    )
  }

  result
}

# ============================================================
# Pooled analysis across all ROIs
# ============================================================

concordance_all <- run_concordance_subset(
  prop_subset = proportions_df,
  analysis_name = "all",
  output_dir = output_dir
)

# Preserve the original standard output filename for compatibility
write.csv(
  concordance_all,
  file = file.path(output_dir, "marker_concordance.csv"),
  row.names = FALSE
)

# ============================================================
# Disease-stratified analyses
# ============================================================

if (!"disease_status" %in% colnames(proportions_df)) {
  stop(
    "Column 'disease_status' was not found in proportions_df. ",
    "Available columns: ",
    paste(colnames(proportions_df), collapse = ", "),
    call. = FALSE
  )
}

disease_levels <- unique(proportions_df$disease_status)
disease_levels <- disease_levels[!is.na(disease_levels)]

concordance_by_disease <- list()

for (disease_group in disease_levels) {

  disease_subset <- proportions_df[
    proportions_df$disease_status == disease_group,
    ,
    drop = FALSE
  ]

  safe_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    disease_group
  )

  concordance_by_disease[[disease_group]] <- run_concordance_subset(
    prop_subset = disease_subset,
    analysis_name = safe_name,
    output_dir = output_dir
  )
}

message("Marker-gene concordance check completed.")
