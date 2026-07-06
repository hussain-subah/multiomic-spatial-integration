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

expression_csv <- "data/CAA-AD_expression_wide.csv"
proportions_csv <- "results/cell_proportions/spatial_celltype_proportions_for_R.csv"
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

concordance <- evaluate_marker_concordance(
  expr_mat = expr_norm,
  proportions_df = proportions_df,
  marker_genes_list = marker_genes_list,
  roi_id_col = roi_id_col,
  proportion_col = proportion_col,
  celltype_col = celltype_col,
  method = method
)

write.csv(
  concordance,
  file = file.path(output_dir, "marker_concordance.csv"),
  row.names = FALSE
)

# ============================================================
# Flag cell types worth a closer look: no positive, significant
# concordance between inferred proportion and its own markers.
# ============================================================

flagged <- concordance[
  is.na(concordance$statistic) | concordance$statistic <= 0 | concordance$p_adj >= 0.05,
]

if (nrow(flagged) > 0) {
  message(sprintf(
    "%d of %d cell type(s) show no positive, significant concordance with their own markers:",
    nrow(flagged), nrow(concordance)
  ))
  print(flagged[, c("celltype", "n_markers_found", "statistic", "p_adj", "n")])
} else {
  message("All cell types show positive, significant concordance with their own markers.")
}

message("Marker-gene concordance check completed.")
