#' Utilities for re-running contrasts with an unstable cell type excluded
#'
#' `Fibroblast` accounts for ~45-64% of the entire cell-type composition
#' (see `results/cell_proportions` sanity checks) -- an implausible share
#' for a rare perivascular cell type in brain tissue -- and its ~19-point
#' swing between AD-CAA (45%) and Control (64%) is large enough to
#' mechanically inflate the apparent share of nearly every other cell type
#' through the sum-to-one (closure) constraint alone, independent of any
#' real biological change. `exclude_and_renormalize()` drops one or more
#' cell types from the composition and renormalizes the remaining types
#' back to sum to 1 within each ROI, so downstream contrasts are not
#' riding on that swing.
#'
#' @keywords internal
NULL


#' Exclude cell type(s) and renormalize remaining proportions to sum to 1
#'
#' @param df Long-format abundance dataframe (one row per ROI x celltype).
#' @param exclude Character vector of celltype(s) to drop.
#' @param abundance_col Column containing proportions.
#' @param roi_col Column identifying each ROI (rows are unique on
#'   `roi_col` x `celltype_col`).
#' @param celltype_col Cell-type column.
#'
#' @return Long-format dataframe with `exclude` rows removed and
#'   `abundance_col` renormalized to sum to 1 per ROI over the remaining
#'   cell types.
#' @export
exclude_and_renormalize <- function(df,
                                    exclude,
                                    abundance_col = "rel_abundance",
                                    roi_col = "ROI_ID",
                                    celltype_col = "celltype") {
  dropped <- df |>
    dplyr::filter(.data[[celltype_col]] %in% exclude)

  dropped_means <- tapply(
    dropped[[abundance_col]],
    dropped[[celltype_col]],
    mean,
    na.rm = TRUE
  )

  message(
    "Excluding ", paste(exclude, collapse = ", "),
    " (mean share before exclusion: ",
    paste(
      sprintf("%s=%.1f%%", names(dropped_means), 100 * dropped_means),
      collapse = ", "
    ),
    ")"
  )

  df |>
    dplyr::filter(!.data[[celltype_col]] %in% exclude) |>
    dplyr::group_by(.data[[roi_col]]) |>
    dplyr::mutate(
      .roi_total = sum(.data[[abundance_col]], na.rm = TRUE),
      "{abundance_col}" := .data[[abundance_col]] / .roi_total
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.roi_total)
}
