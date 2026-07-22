#' Calibration validation gate for cross-model absolute abundance
#'
#' AD-CAA and Control were deconvolved by two **independently trained**
#' SpaceJam/Cell2location models with no shared calibration anchor. A
#' per-ROI total-count offset (see `R/absolute_abundance_utils.R`) makes
#' abundance comparable *within* a condition, but the two posteriors could
#' still carry an arbitrary, condition-specific scale on top of any real
#' biology -- a "gauge freedom" that a within-condition offset cannot detect,
#' let alone fix.
#'
#' This file provides a diagnostic, not a proof: it bounds how much of a
#' cross-condition contrast's absolute-abundance signal could plausibly be
#' this scale artifact, and forces a documented `proceed` /
#' `proceed_with_correction` / `do_not_proceed` decision rather than silently
#' assuming the two posteriors are on the same scale.
#'
#' The Amyloid-effect contrast (Amyloid vs AmyloidFree, entirely within
#' AD-CAA -- a single posterior) is gauge-safe by construction and is exempt
#' from this whole concern regardless of what these checks find.
#'
#' @keywords internal
NULL


# ============================================================
# Calibration ratio table
# ============================================================

#' Build the per-ROI calibration ratio table
#'
#' Sums `abs_abundance` across all cell types per ROI, joins against
#' per-ROI total counts, and computes `calibration_ratio = total_abs_abundance
#' / total_counts` -- the quantity a shared-scale posterior would keep
#' roughly constant within a condition, and that a scale artifact between
#' the two independently-trained models would shift systematically by
#' `disease_status`.
#'
#' @param abs_df Long-format absolute-abundance dataframe (one row per ROI x
#'   celltype), e.g. `results/cell_proportions/roi_celltype_abundance_long_abs.csv`.
#' @param total_counts_df Per-ROI total-count table, e.g.
#'   `results/cell_proportions/roi_total_counts.csv`.
#' @param roi_col ROI identifier column, present in both inputs.
#' @param abundance_col Absolute-abundance column in `abs_df`.
#' @param meta_cols ROI-level metadata columns to carry through from
#'   `abs_df` (duplicated per cell-type row there, so the first value per
#'   ROI is used).
#'
#' @return Dataframe with one row per ROI: `roi_col`, `meta_cols`,
#'   `total_abs_abundance`, `total_counts`, `size_factor_mor`,
#'   `calibration_ratio`, `log_calibration_ratio`.
#' @export
compute_calibration_ratio_table <- function(abs_df,
                                            total_counts_df,
                                            roi_col = "ROI_ID",
                                            abundance_col = "abs_abundance",
                                            meta_cols = c("disease_status", "pathology", "region", "Scan_ID")) {
  missing_meta <- setdiff(meta_cols, colnames(abs_df))
  if (length(missing_meta) > 0) {
    stop(
      "compute_calibration_ratio_table: missing metadata columns in abs_df: ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }

  if (!roi_col %in% colnames(total_counts_df)) {
    stop(
      "compute_calibration_ratio_table: '", roi_col, "' not found in total_counts_df.",
      call. = FALSE
    )
  }

  roi_summary <- abs_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(roi_col, meta_cols)))) |>
    dplyr::summarise(
      total_abs_abundance = sum(.data[[abundance_col]], na.rm = TRUE),
      .groups = "drop"
    )

  abs_only <- setdiff(roi_summary[[roi_col]], total_counts_df[[roi_col]])
  counts_only <- setdiff(total_counts_df[[roi_col]], roi_summary[[roi_col]])

  if (length(abs_only) > 0 || length(counts_only) > 0) {
    stop(
      "compute_calibration_ratio_table: incomplete join between abs_abundance ",
      "ROIs and total_counts ROIs -- these must match exactly.\n",
      "In abs_abundance but not total_counts (", length(abs_only), "): ",
      paste(utils::head(abs_only, 10), collapse = ", "),
      if (length(abs_only) > 10) ", ..." else "", "\n",
      "In total_counts but not abs_abundance (", length(counts_only), "): ",
      paste(utils::head(counts_only, 10), collapse = ", "),
      if (length(counts_only) > 10) ", ..." else "",
      call. = FALSE
    )
  }

  roi_summary |>
    dplyr::inner_join(total_counts_df, by = roi_col) |>
    dplyr::mutate(
      calibration_ratio = .data$total_abs_abundance / .data$total_counts,
      log_calibration_ratio = log(.data$calibration_ratio)
    )
}


# ============================================================
# Gate decision
# ============================================================

#' Run the calibration gate check for the cross-condition contrasts
#'
#' Fits `log_calibration_ratio ~ disease_status + pathology + region +
#' (1 | Scan_ID)` and decides whether the `disease_status` shift in
#' calibration ratio is small enough to ignore, present but plausibly
#' biological (consistent region pattern across conditions), or large/
#' inconsistent enough that cross-condition absolute-abundance contrasts
#' should not be trusted.
#'
#' The `biological_threshold_log2` is a **user-supplied plausibility
#' bound**, not an asserted "correct" value -- how large a scale shift is
#' still plausibly real biology (vs. an artifact of two independently
#' calibrated posteriors) is a domain judgment call, not a statistical one.
#'
#' @param calib_df Output of `compute_calibration_ratio_table()`.
#' @param biological_threshold_log2 Log2 magnitude below which the
#'   `disease_status` coefficient is considered small relative to
#'   within-condition noise (default 0.5, ~40% shift).
#' @param disease_col,pathology_col,region_col,scan_col Column names.
#'
#' @return List: `decision` (`"proceed"`, `"proceed_with_correction"`, or
#'   `"do_not_proceed"`), `disease_coefficient_log2`, `disease_se_log2`,
#'   `biological_threshold_log2`, `within_condition_cv` (per-condition CV of
#'   `calibration_ratio`), `region_pattern_correlation` (correlation of
#'   region-mean `calibration_ratio` between conditions), and `fit` (the
#'   fitted model, for inspection).
#' @export
run_calibration_gate_check <- function(calib_df,
                                       biological_threshold_log2 = 0.5,
                                       disease_col = "disease_status",
                                       pathology_col = "pathology",
                                       region_col = "region",
                                       scan_col = "Scan_ID") {
  calib_df[[disease_col]] <- factor(calib_df[[disease_col]], levels = c("Control", "AD-CAA"))
  calib_df[[pathology_col]] <- factor(calib_df[[pathology_col]])
  calib_df[[region_col]] <- factor(calib_df[[region_col]])
  calib_df[[scan_col]] <- factor(calib_df[[scan_col]])

  fit <- glmmTMB::glmmTMB(
    stats::as.formula(
      paste0(
        "log_calibration_ratio ~ ", disease_col, " + ", pathology_col,
        " + ", region_col, " + (1 | ", scan_col, ")"
      )
    ),
    data = calib_df,
    family = glmmTMB::gaussian()
  )

  coefs <- summary(fit)$coefficients$cond
  disease_row <- grep(paste0("^", disease_col), rownames(coefs))[1]

  disease_coef_log2 <- coefs[disease_row, "Estimate"] / log(2)
  disease_se_log2 <- coefs[disease_row, "Std. Error"] / log(2)

  within_condition_cv <- calib_df |>
    dplyr::group_by(.data[[disease_col]]) |>
    dplyr::summarise(
      cv = stats::sd(.data$calibration_ratio, na.rm = TRUE) /
        mean(.data$calibration_ratio, na.rm = TRUE),
      .groups = "drop"
    )

  region_means <- calib_df |>
    dplyr::group_by(.data[[disease_col]], .data[[region_col]]) |>
    dplyr::summarise(
      mean_ratio = mean(.data$calibration_ratio, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(disease_col), values_from = mean_ratio)

  region_pattern_correlation <- if (nrow(region_means) >= 2 && ncol(region_means) >= 3) {
    stats::cor(region_means[[2]], region_means[[3]], use = "complete.obs")
  } else {
    NA_real_
  }

  max_within_cv <- max(within_condition_cv$cv, na.rm = TRUE)
  abs_coef <- abs(disease_coef_log2)

  decision <- if (abs_coef < biological_threshold_log2 && abs_coef < 2 * max_within_cv) {
    "proceed"
  } else if (!is.na(region_pattern_correlation) && region_pattern_correlation > 0.5) {
    "proceed_with_correction"
  } else {
    "do_not_proceed"
  }

  list(
    decision = decision,
    disease_coefficient_log2 = disease_coef_log2,
    disease_se_log2 = disease_se_log2,
    biological_threshold_log2 = biological_threshold_log2,
    within_condition_cv = within_condition_cv,
    region_pattern_correlation = region_pattern_correlation,
    fit = fit
  )
}
