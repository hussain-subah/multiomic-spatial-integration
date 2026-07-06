#' Leave-one-Scan-out robustness utilities
#'
#' With as few as 3-4 distinct Scan_IDs per condition backing the random
#' intercept in each beta mixed-effects contrast model, a single Scan_ID
#' can have an outsized influence on an effect estimate. These utilities
#' refit a contrast function once per excluded Scan_ID and summarize how
#' much each region x celltype x contrast estimate moves, and whether its
#' significance call flips, when any single scan is dropped.
#'
#' @keywords internal
NULL


#' Refit a contrast function leaving out one Scan_ID at a time
#'
#' @param df Cleaned dataframe (as returned by `prepare_spatial_proportion_data`).
#' @param model_fn One of the contrast functions in `R/contrast_utils.R`
#'   (e.g. `run_amyloid_effect`), taking `df` and `abundance_col` and
#'   returning a list with a `contrasts` data frame.
#' @param abundance_col Relative abundance column.
#' @param scan_col Column identifying the random-effect cluster to leave
#'   out.
#'
#' @return List with `full` (contrasts from the complete data) and
#'   `leave_one_out` (contrasts from each excluded-scan refit, with an
#'   `excluded_scan` column).
#' @export
run_leave_one_scan_out <- function(df,
                                   model_fn,
                                   abundance_col = "rel_abundance",
                                   scan_col = "Scan_ID") {
  full_res <- model_fn(df, abundance_col = abundance_col)
  full_contrasts <- full_res$contrasts

  scans <- unique(df[[scan_col]])
  loo_results <- list()

  for (s in scans) {
    df_loo <- df[df[[scan_col]] != s, ]

    res_loo <- try(model_fn(df_loo, abundance_col = abundance_col), silent = TRUE)

    if (inherits(res_loo, "try-error") || nrow(res_loo$contrasts) == 0) {
      message(
        "run_leave_one_scan_out: skipping excluded_scan = ", s,
        " (refit failed or produced no contrasts)."
      )
      next
    }

    loo_results[[as.character(s)]] <- res_loo$contrasts |>
      dplyr::mutate(excluded_scan = s)
  }

  list(
    full = full_contrasts,
    leave_one_out = dplyr::bind_rows(loo_results)
  )
}


#' Summarize leave-one-Scan-out robustness
#'
#' For each region x celltype x contrast, reports how much the estimate
#' moves and whether the significance call flips across all
#' single-scan-excluded refits.
#'
#' @param full_contrasts Contrasts data frame from the complete-data fit.
#' @param loo_contrasts Contrasts data frame from `run_leave_one_scan_out()`
#'   (must include an `excluded_scan` column).
#' @param join_cols Columns identifying a unique contrast row.
#' @param sig_cutoff Adjusted p-value threshold defining "significant".
#'
#' @return Data frame: `join_cols`, full_estimate, full_p_adj,
#'   min_loo_estimate, max_loo_estimate, max_abs_delta, n_loo,
#'   any_sig_flip -- sorted by `max_abs_delta` descending (most fragile
#'   results first).
#' @export
summarize_robustness <- function(full_contrasts,
                                 loo_contrasts,
                                 join_cols = c("region", "celltype", "contrast"),
                                 sig_cutoff = 0.05) {
  if (nrow(loo_contrasts) == 0) {
    return(data.frame())
  }

  full_sub <- full_contrasts |>
    dplyr::select(dplyr::all_of(c(join_cols, "estimate", "p_adj"))) |>
    dplyr::rename(full_estimate = estimate, full_p_adj = p_adj)

  loo_sub <- loo_contrasts |>
    dplyr::select(dplyr::all_of(c(join_cols, "estimate", "p_adj", "excluded_scan"))) |>
    dplyr::rename(loo_estimate = estimate, loo_p_adj = p_adj)

  joined <- dplyr::inner_join(loo_sub, full_sub, by = join_cols)

  joined |>
    dplyr::mutate(
      sig_full = .data$full_p_adj < sig_cutoff,
      sig_loo = .data$loo_p_adj < sig_cutoff,
      abs_delta = abs(.data$loo_estimate - .data$full_estimate)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(join_cols))) |>
    dplyr::summarise(
      full_estimate = dplyr::first(.data$full_estimate),
      full_p_adj = dplyr::first(.data$full_p_adj),
      min_loo_estimate = min(.data$loo_estimate, na.rm = TRUE),
      max_loo_estimate = max(.data$loo_estimate, na.rm = TRUE),
      max_abs_delta = max(.data$abs_delta, na.rm = TRUE),
      n_loo = dplyr::n(),
      any_sig_flip = any(.data$sig_loo != dplyr::first(.data$sig_full)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$max_abs_delta))
}
