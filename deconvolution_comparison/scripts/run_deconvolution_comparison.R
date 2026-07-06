# ============================================================
# Compare all deconvolution method outputs against the Cell2location baseline
#
# Loads every results/deconvolution_comparison/<method>_proportions.csv plus
# the existing Cell2location result, and writes:
#   - per (method, celltype) concordance vs. baseline
#   - a method x method correlation matrix
#   - a per-method mean absolute proportion difference vs. baseline
#
# There is no ground truth for the real ROIs, so these measure agreement,
# not accuracy.
# ============================================================

source("deconvolution_comparison/R/common_utils.R")
source("deconvolution_comparison/R/comparison_utils.R")

# ============================================================
# Config -- adjust paths for your environment
# ============================================================

method_output_dir <- "results/deconvolution_comparison"
baseline_csv <- "results/cell_proportions/spatial_celltype_proportions_for_R.csv"
summary_dir <- file.path(method_output_dir, "comparison_summary")

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load method outputs + baseline
# ============================================================

baseline_long <- load_cell2location_baseline(baseline_csv)

long_df <- load_all_method_outputs(method_output_dir, baseline_long = baseline_long)

methods_present <- unique(long_df$method)
message("Methods loaded: ", paste(methods_present, collapse = ", "))

if (length(setdiff(methods_present, "Cell2location")) == 0) {
  stop("No non-baseline method outputs found in ", method_output_dir,
       " -- run run_all_deconvolution.R first.", call. = FALSE)
}

# ============================================================
# Concordance vs. baseline (per method x celltype)
# ============================================================

concordance <- concordance_vs_baseline(long_df, baseline_method = "Cell2location")

write.csv(
  concordance,
  file = file.path(summary_dir, "concordance_vs_cell2location.csv"),
  row.names = FALSE
)

# ============================================================
# Cross-method correlation matrix
# ============================================================

cross_cor <- cross_method_correlation(long_df)

write.csv(
  as.data.frame(cross_cor),
  file = file.path(summary_dir, "cross_method_correlation.csv"),
  row.names = TRUE
)

# ============================================================
# Per-method mean absolute difference vs. baseline
# ============================================================

mad_summary <- mean_abs_diff_vs_baseline(long_df, baseline_method = "Cell2location")

write.csv(
  mad_summary,
  file = file.path(summary_dir, "mean_abs_diff_vs_cell2location.csv"),
  row.names = FALSE
)

# ============================================================
# Print a per-method headline (median per-celltype correlation)
# ============================================================

headline <- aggregate(
  pearson_r ~ method,
  data = concordance,
  FUN = function(x) median(x, na.rm = TRUE)
)
colnames(headline)[2] <- "median_celltype_pearson_r_vs_baseline"
headline <- merge(headline, mad_summary[, c("method", "mean_abs_diff")], by = "method")
headline <- headline[order(-headline$median_celltype_pearson_r_vs_baseline), ]

message("\nPer-method agreement with Cell2location (higher r / lower diff = closer):")
print(headline, row.names = FALSE)

message("\nComparison summaries written to ", summary_dir)
