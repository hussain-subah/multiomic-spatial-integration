# ============================================================
# Compare deconvolution outputs against Cell2location baseline
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source(
  "deconvolution_comparison/R/common_utils.R"
)

source(
  "deconvolution_comparison/R/comparison_utils.R"
)

# ============================================================
# Configuration
# ============================================================

output_dir <- "results/deconvolution_comparison"

baseline_file <- paste0(
  "results/cell_proportions/",
  "spatial_celltype_proportions_for_R.csv"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Load Cell2location baseline
# ============================================================

message("Loading Cell2location baseline...")

baseline <- load_cell2location_baseline(
  proportions_csv = baseline_file,
  roi_id_col = "ROI_ID",
  celltype_col = "celltype",
  proportion_col = "rel_abundance"
)

message(
  "Baseline: ",
  dplyr::n_distinct(baseline$ROI_ID),
  " ROIs x ",
  dplyr::n_distinct(baseline$celltype),
  " cell types."
)

# ============================================================
# Find comparison-method outputs
# ============================================================

method_files <- list.files(
  output_dir,
  pattern = "_proportions\\.csv$",
  full.names = TRUE
)

method_files <- method_files[
  !grepl(
    "Cell2location",
    basename(method_files),
    ignore.case = TRUE
  )
]

if (length(method_files) == 0) {
  stop(
    "No comparison-method proportion files were found in ",
    output_dir,
    ". Expected files named <Method>_proportions.csv.",
    call. = FALSE
  )
}

message(
  "Found ",
  length(method_files),
  " comparison method output(s): ",
  paste(
    basename(method_files),
    collapse = ", "
  )
)

# ============================================================
# Storage for combined summaries
# ============================================================

all_celltype_results <- list()
all_lineage_results <- list()
all_usage_results <- list()
all_failed_rois <- list()

# ============================================================
# Run comparisons
# ============================================================

for (method_file in method_files) {

  method_name <- sub(
    "_proportions\\.csv$",
    "",
    basename(method_file)
  )

  message(
    "\n============================================================"
  )
  message("Comparing method: ", method_name)
  message(
    "============================================================"
  )

  method_df <- read_csv(
    method_file,
    show_col_types = FALSE
  )

  required_columns <- c(
    "method",
    "ROI_ID",
    "celltype",
    "proportion"
  )

  missing_columns <- setdiff(
    required_columns,
    names(method_df)
  )

  if (length(missing_columns) > 0) {
    warning(
      "Skipping ",
      method_name,
      ": missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )

    next
  }

  # ----------------------------------------------------------
  # Basic integrity checks
  # ----------------------------------------------------------

  if (anyDuplicated(
    method_df[
      c(
        "ROI_ID",
        "celltype"
      )
    ]
  )) {
    warning(
      "Method ",
      method_name,
      " contains duplicate ROI_ID x celltype rows.",
      call. = FALSE
    )
  }

  unexpected_methods <- setdiff(
    unique(method_df$method),
    method_name
  )

  if (length(unexpected_methods) > 0) {
    warning(
      "Method column in ",
      basename(method_file),
      " contains labels different from the filename-derived method name: ",
      paste(
        unique(method_df$method),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Factor usage and sparsity
  # ----------------------------------------------------------

  usage <- summarize_factor_usage(
    method_df
  ) %>%
    mutate(
      method = method_name,
      .before = 1
    )

  write_csv(
    usage,
    file.path(
      output_dir,
      paste0(
        method_name,
        "_factor_usage.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # Failed ROIs
  # ----------------------------------------------------------

  failed_rois <- identify_failed_rois(
    method_df
  ) %>%
    mutate(
      method = method_name,
      .before = 1
    )

  write_csv(
    failed_rois,
    file.path(
      output_dir,
      paste0(
        method_name,
        "_failed_rois.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # Fine cell-type comparison
  # ----------------------------------------------------------

  by_celltype <- compare_methods_by_celltype(
    method_df = method_df,
    baseline_df = baseline
  ) %>%
    mutate(
      method = method_name,
      .before = 1
    )

  write_csv(
    by_celltype,
    file.path(
      output_dir,
      paste0(
        method_name,
        "_vs_Cell2location_by_celltype.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # Broad-lineage comparison
  # ----------------------------------------------------------

  by_lineage <- compare_methods_by_lineage(
    method_df = method_df,
    baseline_df = baseline
  ) %>%
    mutate(
      method = method_name,
      .before = 1
    )

  write_csv(
    by_lineage,
    file.path(
      output_dir,
      paste0(
        method_name,
        "_vs_Cell2location_by_lineage.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # Console summary
  # ----------------------------------------------------------

  evaluable_roi_ids <- method_df %>%
    group_by(ROI_ID) %>%
    summarise(
      evaluable = any(
        is.finite(proportion)
      ),
      .groups = "drop"
    ) %>%
    filter(evaluable) %>%
    pull(ROI_ID)

  message(
    "Evaluable ROIs: ",
    length(evaluable_roi_ids)
  )

  message(
    "Failed ROIs: ",
    nrow(failed_rois)
  )

  message(
    "Factors never used: ",
    sum(
      usage$n_nonzero == 0,
      na.rm = TRUE
    )
  )

  message(
    "Factors used in >=10% of evaluable ROIs: ",
    sum(
      usage$fraction_nonzero >= 0.10,
      na.rm = TRUE
    )
  )

  evaluable_celltypes <- by_celltype %>%
    filter(
      comparison_status == "evaluable"
    )

  fine_state_median <- if (
    nrow(evaluable_celltypes) == 0 ||
    all(is.na(evaluable_celltypes$spearman_rho))
  ) {
    NA_real_
  } else {
    median(
      evaluable_celltypes$spearman_rho,
      na.rm = TRUE
    )
  }

  message(
    "Median fine-state Spearman rho among evaluable factors: ",
    if (is.na(fine_state_median)) {
      "NA"
    } else {
      round(
        fine_state_median,
        3
      )
    }
  )

  lineage_median <- if (
    nrow(by_lineage) == 0 ||
    all(is.na(by_lineage$spearman_rho))
  ) {
    NA_real_
  } else {
    median(
      by_lineage$spearman_rho,
      na.rm = TRUE
    )
  }

  message(
    "Median lineage Spearman rho across all lineages: ",
    if (is.na(lineage_median)) {
      "NA"
    } else {
      round(
        lineage_median,
        3
      )
    }
  )

  well_used_lineages <- by_lineage %>%
    filter(
      method_fraction_nonzero >= 0.10
    )

  message(
    "Lineages used in >=10% of evaluable ROIs: ",
    nrow(well_used_lineages),
    " of ",
    nrow(by_lineage)
  )

  well_used_lineage_median <- if (
    nrow(well_used_lineages) == 0 ||
    all(is.na(well_used_lineages$spearman_rho))
  ) {
    NA_real_
  } else {
    median(
      well_used_lineages$spearman_rho,
      na.rm = TRUE
    )
  }

  message(
    "Median lineage Spearman rho among lineages used in >=10% of ROIs: ",
    if (is.na(well_used_lineage_median)) {
      "NA"
    } else {
      round(
        well_used_lineage_median,
        3
      )
    }
  )

  message(
    "Lineages with Spearman rho >= 0.5: ",
    sum(
      by_lineage$spearman_rho >= 0.5,
      na.rm = TRUE
    ),
    " of ",
    nrow(by_lineage)
  )

  # ----------------------------------------------------------
  # Store for combined summaries
  # ----------------------------------------------------------

  all_usage_results[[method_name]] <- usage
  all_failed_rois[[method_name]] <- failed_rois
  all_celltype_results[[method_name]] <- by_celltype
  all_lineage_results[[method_name]] <- by_lineage
}

# ============================================================
# Combined cross-method summaries
# ============================================================

if (length(all_usage_results) > 0) {
  write_csv(
    bind_rows(all_usage_results),
    file.path(
      output_dir,
      "all_methods_factor_usage.csv"
    )
  )
}

if (length(all_failed_rois) > 0) {
  write_csv(
    bind_rows(all_failed_rois),
    file.path(
      output_dir,
      "all_methods_failed_rois.csv"
    )
  )
}

if (length(all_celltype_results) > 0) {
  write_csv(
    bind_rows(all_celltype_results),
    file.path(
      output_dir,
      "all_methods_vs_Cell2location_by_celltype.csv"
    )
  )
}

if (length(all_lineage_results) > 0) {
  write_csv(
    bind_rows(all_lineage_results),
    file.path(
      output_dir,
      "all_methods_vs_Cell2location_by_lineage.csv"
    )
  )
}

message(
  "\nDeconvolution comparison completed successfully."
)
