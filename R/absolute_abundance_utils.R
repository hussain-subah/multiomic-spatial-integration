#' Absolute cell-type abundance contrasts (total-count offset)
#'
#' `rel_abundance` sums to 1 per ROI: a real shift in one dominant cell type
#' (`Fibroblast`, 45-64% of the panel) mechanically drags every other cell
#' type's proportion in the opposite direction through the closure
#' constraint alone (see `R/contrast_utils.R` / `docs/analysis_methodology.md`
#' Sec. 1). This file models absolute abundance instead, with a per-ROI
#' total-count offset (limma/DESeq2-style) so no single cell type can swing
#' the others.
#'
#' Mirrors `R/contrast_utils.R`'s structure exactly -- same per-region x
#' per-celltype loop, same BH-adjustment scoping, same 4-contrast set --
#' swapping the beta/proportion model for `log(abs_abundance) ~ ... +
#' offset(log_offset)`. Fit with `glmmTMB::gaussian()` rather than `lme4`,
#' to stay consistent with every other modeling file in this repo and add
#' no new dependency; `glmmTMB` supports `offset()` in the formula the same
#' way `lme4` does.
#'
#' Requires `R/contrast_utils.R` to be sourced first (reuses
#' `.add_group_factor()` from there).
#'
#' **Before trusting any output**: confirm `emmeans::contrast()` estimates
#' are unaffected by the offset's assumed reference-grid value. Refit with
#' `log_offset` shifted by an arbitrary constant and check that contrast
#' estimates are unchanged (only the reported marginal means should move) --
#' this is a property of how `emmeans` handles a covariate held at a single
#' value across the whole reference grid, but must be verified against a
#' real fit, not assumed.
#'
#' **This is a diagnostic offset fix, not a fix for cross-model gauge
#' freedom** -- see `R/absolute_calibration_utils.R` for the separate gate
#' that bounds the risk of AD-CAA and Control being deconvolved by two
#' independently-trained, uncalibrated posteriors.
#'
#' @keywords internal
NULL


# ============================================================
# Data preparation
# ============================================================

#' Prepare absolute abundance data with a log total-count offset
#'
#' @param df Long-format absolute-abundance dataframe (one row per ROI x
#'   celltype), joined against a total-count table (e.g.
#'   `roi_total_counts.csv`) so `offset_col` is present.
#' @param offset_col Per-ROI offset column. `"size_factor_mor"` (primary,
#'   DESeq2-style median-of-ratios) or `"total_counts"` (raw-sum sensitivity
#'   variant).
#' @param abundance_col Absolute-abundance column.
#' @param disease_col,pathology_col,region_col,scan_col,celltype_col Column
#'   names, as in `prepare_spatial_proportion_data()`.
#' @param eps Small value to keep `log()` finite at zero abundance/offset.
#'
#' @return Cleaned dataframe with `log_abs_abundance` and `log_offset` added.
#' @export
prepare_absolute_abundance_data <- function(df,
                                            offset_col = "size_factor_mor",
                                            abundance_col = "abs_abundance",
                                            disease_col = "disease_status",
                                            pathology_col = "pathology",
                                            region_col = "region",
                                            scan_col = "Scan_ID",
                                            celltype_col = "celltype",
                                            eps = 1e-8) {
  required_cols <- c(
    abundance_col, offset_col, disease_col, pathology_col,
    region_col, scan_col, celltype_col
  )

  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "prepare_absolute_abundance_data: missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df[[abundance_col]] <- as.numeric(df[[abundance_col]])
  df[[offset_col]] <- as.numeric(df[[offset_col]])

  df$log_abs_abundance <- log(pmax(df[[abundance_col]], eps))
  df$log_offset <- log(pmax(df[[offset_col]], eps))

  df[[disease_col]] <- factor(df[[disease_col]], levels = c("Control", "AD-CAA"))
  df[[pathology_col]] <- factor(df[[pathology_col]], levels = c("AmyloidFree", "Amyloid"))
  df[[region_col]] <- factor(df[[region_col]])
  df[[scan_col]] <- factor(df[[scan_col]])
  df[[celltype_col]] <- factor(df[[celltype_col]])

  return(df)
}


# ============================================================
# Internal helper
# ============================================================

#' Run the offset-adjusted absolute-abundance model safely for one subset
#'
#' @keywords internal
.run_absolute_model <- function(tmp, formula, min_rows = 5, min_scans = 2) {
  if (nrow(tmp) < min_rows) return(NULL)
  if (dplyr::n_distinct(tmp$Scan_ID) < min_scans) return(NULL)

  fit <- try(
    glmmTMB::glmmTMB(
      formula,
      data = tmp,
      family = glmmTMB::gaussian()
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) return(NULL)

  return(fit)
}


# ============================================================
# Amyloid effect
# ============================================================

#' Run amyloid effect model (absolute abundance)
#'
#' Tests Amyloid vs AmyloidFree within AD-CAA only -- entirely within one
#' condition's posterior, so gauge-safe by construction regardless of what
#' the calibration gate finds.
#'
#' Model: `log_abs_abundance ~ pathology + offset(log_offset) + (1 | Scan_ID)`
#'
#' @param df Output of `prepare_absolute_abundance_data()`.
#'
#' @return List with means and contrasts. `fold_ratio = exp(estimate)` on
#'   the contrasts is the offset-adjusted absolute fold-change.
#' @export
run_amyloid_effect_absolute <- function(df) {
  dat_ad <- df |>
    dplyr::filter(.data$disease_status == "AD-CAA")

  mean_results <- list()
  contrast_results <- list()

  for (rg in levels(dat_ad$region)) {
    for (ct in levels(dat_ad$celltype)) {

      tmp <- dat_ad |>
        dplyr::filter(.data$region == rg, .data$celltype == ct) |>
        dplyr::filter(!is.na(.data$log_abs_abundance), is.finite(.data$log_abs_abundance)) |>
        droplevels()

      if (dplyr::n_distinct(tmp$pathology) < 2) next

      fit <- .run_absolute_model(
        tmp = tmp,
        formula = log_abs_abundance ~ pathology + offset(log_offset) + (1 | Scan_ID)
      )

      if (is.null(fit)) next

      emm <- emmeans::emmeans(fit, ~ pathology)
      emm_df <- as.data.frame(emm)
      emm_df$mean_abundance <- exp(emm_df$emmean)

      mean_results[[paste(rg, ct, sep = "__")]] <- emm_df |>
        dplyr::mutate(celltype = ct, region = rg)

      cont <- emmeans::contrast(
        emm,
        method = list(Amyloid_effect = c(-1, 1))
      )

      contrast_results[[paste(rg, ct, sep = "__")]] <-
        as.data.frame(summary(cont, infer = TRUE)) |>
        dplyr::mutate(celltype = ct, region = rg, fold_ratio = exp(.data$estimate))
    }
  }

  means_df <- dplyr::bind_rows(mean_results)

  contrast_df <- dplyr::bind_rows(contrast_results) |>
    dplyr::group_by(.data$region) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup()

  list(
    means = means_df,
    contrasts = contrast_df
  )
}


# ============================================================
# Disease effect
# ============================================================

#' Run disease effect model (absolute abundance)
#'
#' Tests AD-CAA vs Control using AmyloidFree ROIs only. Cross-condition --
#' subject to the calibration gate in `R/absolute_calibration_utils.R`.
#'
#' Model: `log_abs_abundance ~ disease_status + offset(log_offset) + (1 | Scan_ID)`
#'
#' @param df Output of `prepare_absolute_abundance_data()`.
#'
#' @return List with means and contrasts.
#' @export
run_disease_effect_absolute <- function(df) {
  dat_dis <- df |>
    dplyr::filter(.data$pathology == "AmyloidFree")

  mean_results <- list()
  contrast_results <- list()

  for (rg in levels(dat_dis$region)) {
    for (ct in levels(dat_dis$celltype)) {

      tmp <- dat_dis |>
        dplyr::filter(.data$region == rg, .data$celltype == ct) |>
        dplyr::filter(!is.na(.data$log_abs_abundance), is.finite(.data$log_abs_abundance)) |>
        droplevels()

      if (dplyr::n_distinct(tmp$disease_status) < 2) next

      fit <- .run_absolute_model(
        tmp = tmp,
        formula = log_abs_abundance ~ disease_status + offset(log_offset) + (1 | Scan_ID)
      )

      if (is.null(fit)) next

      emm <- emmeans::emmeans(fit, ~ disease_status)
      emm_df <- as.data.frame(emm)
      emm_df$mean_abundance <- exp(emm_df$emmean)

      mean_results[[paste(rg, ct, sep = "__")]] <- emm_df |>
        dplyr::mutate(celltype = ct, region = rg)

      cont <- emmeans::contrast(
        emm,
        method = list(Disease_effect = c(-1, 1))
      )

      contrast_results[[paste(rg, ct, sep = "__")]] <-
        as.data.frame(summary(cont, infer = TRUE)) |>
        dplyr::mutate(celltype = ct, region = rg, fold_ratio = exp(.data$estimate))
    }
  }

  means_df <- dplyr::bind_rows(mean_results)

  contrast_df <- dplyr::bind_rows(contrast_results) |>
    dplyr::group_by(.data$region) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup()

  list(
    means = means_df,
    contrasts = contrast_df
  )
}


# ============================================================
# Weighted overall effect
# ============================================================

#' Run weighted overall AD vs Control effect (absolute abundance)
#'
#' Same three-group design and region-specific amyloid-fraction weighting
#' as `run_weighted_overall_effect()` in `R/contrast_utils.R`. Cross-condition
#' -- subject to the calibration gate.
#'
#' Model: `log_abs_abundance ~ group + offset(log_offset) + (1 | Scan_ID)`
#'
#' @param df Output of `prepare_absolute_abundance_data()`. Requires
#'   `.add_group_factor()` from `R/contrast_utils.R` (source that file first).
#'
#' @return List with means, contrasts, and region weights.
#' @export
run_weighted_overall_effect_absolute <- function(df) {
  overall_dat <- .add_group_factor(df)

  mean_results <- list()
  contrast_results <- list()
  weight_results <- list()

  for (rg in levels(overall_dat$region)) {

    w <- overall_dat |>
      dplyr::filter(.data$region == rg, .data$disease_status == "AD-CAA") |>
      dplyr::summarise(w = mean(.data$pathology == "Amyloid")) |>
      dplyr::pull(.data$w)

    weight_results[[rg]] <- data.frame(region = rg, weight = w)

    for (ct in levels(overall_dat$celltype)) {

      tmp <- overall_dat |>
        dplyr::filter(.data$region == rg, .data$celltype == ct) |>
        dplyr::filter(!is.na(.data$log_abs_abundance), is.finite(.data$log_abs_abundance)) |>
        droplevels()

      if (dplyr::n_distinct(tmp$group) < 2) next

      fit <- .run_absolute_model(
        tmp = tmp,
        formula = log_abs_abundance ~ group + offset(log_offset) + (1 | Scan_ID)
      )

      if (is.null(fit)) next

      emm <- emmeans::emmeans(fit, ~ group)
      emm_df <- as.data.frame(emm)
      emm_df$mean_abundance <- exp(emm_df$emmean)

      mean_results[[paste(rg, ct, sep = "__")]] <- emm_df |>
        dplyr::mutate(celltype = ct, region = rg, weight = w)

      contr <- try(
        emmeans::contrast(
          emm,
          method = list(
            Disease_effect = c(-1, 1, 0),
            Amyloid_effect = c(0, -1, 1),
            AD_overall_vs_Control = c(-1, 1 - w, w)
          )
        ),
        silent = TRUE
      )

      if (inherits(contr, "try-error")) next

      contrast_results[[paste(rg, ct, sep = "__")]] <-
        as.data.frame(summary(contr, infer = TRUE)) |>
        dplyr::mutate(celltype = ct, region = rg, weight = w, fold_ratio = exp(.data$estimate))
    }
  }

  means_df <- dplyr::bind_rows(mean_results)

  contrast_df <- dplyr::bind_rows(contrast_results) |>
    dplyr::group_by(.data$region, .data$contrast) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup()

  weights_df <- dplyr::bind_rows(weight_results)

  list(
    means = means_df,
    contrasts = contrast_df,
    weights = weights_df
  )
}


# ============================================================
# Max pathology effect
# ============================================================

#' Run maximal pathology effect (absolute abundance)
#'
#' Tests AD Amyloid ROIs vs Control AmyloidFree ROIs. Cross-condition --
#' subject to the calibration gate.
#'
#' Model: `log_abs_abundance ~ group2 + offset(log_offset) + (1 | Scan_ID)`
#'
#' @param df Output of `prepare_absolute_abundance_data()`.
#'
#' @return List with means, contrasts, and fit status.
#' @export
run_max_pathology_effect_absolute <- function(df) {
  comparison_dat <- df |>
    dplyr::mutate(
      group2 = dplyr::case_when(
        .data$disease_status == "Control" & .data$pathology == "AmyloidFree" ~ "Control_AmyloidFree",
        .data$disease_status == "AD-CAA" & .data$pathology == "Amyloid" ~ "AD_Amyloid",
        TRUE ~ NA_character_
      ),
      group2 = factor(
        .data$group2,
        levels = c("Control_AmyloidFree", "AD_Amyloid")
      )
    ) |>
    dplyr::filter(!is.na(.data$group2))

  mean_results <- list()
  contrast_results <- list()
  fit_status <- list()

  for (rg in levels(comparison_dat$region)) {
    for (ct in levels(comparison_dat$celltype)) {

      tmp <- comparison_dat |>
        dplyr::filter(.data$region == rg, .data$celltype == ct) |>
        dplyr::filter(!is.na(.data$log_abs_abundance), is.finite(.data$log_abs_abundance)) |>
        droplevels()

      status <- data.frame(
        region = rg,
        celltype = ct,
        n_rows = nrow(tmp),
        n_scans = dplyr::n_distinct(tmp$Scan_ID),
        n_groups = dplyr::n_distinct(tmp$group2),
        fit_ok = FALSE,
        pd_hessian = NA,
        stringsAsFactors = FALSE
      )

      if (nrow(tmp) < 5 || dplyr::n_distinct(tmp$Scan_ID) < 2 ||
          dplyr::n_distinct(tmp$group2) < 2) {
        fit_status[[paste(rg, ct, sep = "__")]] <- status
        next
      }

      fit <- .run_absolute_model(
        tmp = tmp,
        formula = log_abs_abundance ~ group2 + offset(log_offset) + (1 | Scan_ID)
      )

      if (is.null(fit)) {
        fit_status[[paste(rg, ct, sep = "__")]] <- status
        next
      }

      status$fit_ok <- TRUE
      status$pd_hessian <- isTRUE(fit$sdr$pdHess)

      emm <- emmeans::emmeans(fit, ~ group2)
      emm_df <- as.data.frame(emm)
      emm_df$mean_abundance <- exp(emm_df$emmean)

      mean_results[[paste(rg, ct, sep = "__")]] <- emm_df |>
        dplyr::mutate(celltype = ct, region = rg)

      contr <- try(
        emmeans::contrast(
          emm,
          method = list(AD_Amyloid_vs_Control = c(-1, 1))
        ),
        silent = TRUE
      )

      if (!inherits(contr, "try-error")) {
        contrast_results[[paste(rg, ct, sep = "__")]] <-
          as.data.frame(summary(contr, infer = TRUE)) |>
          dplyr::mutate(celltype = ct, region = rg, fold_ratio = exp(.data$estimate))
      }

      fit_status[[paste(rg, ct, sep = "__")]] <- status
    }
  }

  means_df <- dplyr::bind_rows(mean_results)

  contrast_df <- dplyr::bind_rows(contrast_results) |>
    dplyr::group_by(.data$region) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup()

  fit_status_df <- dplyr::bind_rows(fit_status)

  list(
    means = means_df,
    contrasts = contrast_df,
    fit_status = fit_status_df
  )
}


# ============================================================
# Summary formatting
# ============================================================

#' Format absolute-abundance contrast output into unified summary table
#'
#' Same shape as `format_contrast_summary()` in `R/contrast_utils.R`, plus
#' `log2_fold_change` (the honestly-named offset-adjusted effect size).
#' Also emits `log2_OR` as an alias of `log2_fold_change` so
#' `plot_contrast_effsize_heatmap()` and `.effect_size_column_info()`
#' (`R/viz_spatial_stats.R`, `R/viz_model_comparison.R`), which key off a
#' column literally named `log2_OR`, work unmodified against this pipeline's
#' output too.
#'
#' @param contrast_df Contrast dataframe.
#' @param contrast_type Label for contrast class.
#'
#' @return Standardized contrast summary.
#' @export
format_absolute_contrast_summary <- function(contrast_df, contrast_type) {
  contrast_df$log2_fold_change <- contrast_df$estimate / log(2)
  contrast_df$log2_OR <- contrast_df$log2_fold_change

  contrast_df %>%
    dplyr::mutate(type = contrast_type) %>%
    dplyr::select(
      dplyr::any_of(c(
        "celltype",
        "region",
        "contrast",
        "type",
        "log2_fold_change",
        "log2_OR",
        "fold_ratio",
        "estimate",
        "SE",
        "p.value",
        "p_adj"
      ))
    )
}
