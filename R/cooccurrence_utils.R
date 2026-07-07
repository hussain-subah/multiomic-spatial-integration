#' Cell-type co-occurrence utilities
#'
#' Cell-type proportions are compositional (they sum to ~1 within an ROI),
#' so naive Pearson correlation between cell-type columns is biased toward
#' negative values purely from that sum-to-one constraint ("closure"),
#' independent of any real biological relationship. These utilities apply a
#' centered-log-ratio (CLR) transform before correlating, which is the
#' standard fix for compositional data.
#'
#' Composition also differs structurally by vascular compartment (e.g. SMCs
#' are near-exclusive to Arteries), so correlations are computed separately
#' within each `region` stratum rather than pooled — pooling would mostly
#' recover "which region is this ROI" rather than genuine co-occurrence.
#' Stratifying further by disease/pathology `group` is a natural extension
#' but is left out here given the smaller per-cell sample sizes that would
#' result.
#'
#' @keywords internal
NULL


#' Centered log-ratio transform
#'
#' Row-wise CLR: for each row (ROI), subtracts the mean of the log-proportions
#' across columns (cell types) from each log-proportion. Removes the
#' sum-to-one constraint that induces spurious negative correlations in raw
#' compositional data.
#'
#' @param prop_mat Matrix or data frame of proportions (rows = ROIs,
#'   columns = cell types).
#' @param eps Small value proportions are clipped to before taking logs.
#'
#' @return Data frame of CLR-transformed values, same shape as `prop_mat`.
#' @export
compute_clr <- function(prop_mat, eps = 1e-6) {
  prop_mat <- as.matrix(prop_mat)
  prop_mat <- pmin(pmax(prop_mat, eps), 1 - eps)

  log_mat <- log(prop_mat)
  row_geom_mean_log <- rowMeans(log_mat)

  clr_mat <- log_mat - row_geom_mean_log

  as.data.frame(clr_mat)
}


#' Run CLR-transformed cell-type co-occurrence analysis
#'
#' Pivots the long proportion table to one row per ROI (`Scan_ID`) x
#' `stratify_by` with one column per cell type, applies `compute_clr()`,
#' and computes all pairwise Pearson correlations between cell types within
#' each `stratify_by` stratum.
#'
#' @param df Cleaned dataframe (as returned by `prepare_spatial_proportion_data`).
#' @param abundance_col Relative abundance column.
#' @param stratify_by Column to stratify correlations by. Defaults to
#'   `"region"`.
#' @param min_n Minimum number of complete ROIs required within a stratum
#'   before computing correlations for it.
#'
#' @return Long-format data frame: stratify_by column, celltype_1,
#'   celltype_2, r, p.value, p_adj, n.
#' @export
run_celltype_cooccurrence <- function(df,
                                      abundance_col = "rel_abundance",
                                      stratify_by = "region",
                                      min_n = 3) {
  wide_dat <- df |>
    dplyr::select(dplyr::all_of(c("Scan_ID", stratify_by, "celltype", abundance_col))) |>
    dplyr::filter(
      !is.na(.data[[abundance_col]]),
      !is.na(.data[[stratify_by]]),
      !is.na(.data$Scan_ID),
      !is.na(.data$celltype)
    ) |>
    dplyr::group_by(.data$Scan_ID, .data[[stratify_by]], .data$celltype) |>
    dplyr::summarise(
      abundance = mean(.data[[abundance_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(c("Scan_ID", stratify_by)),
      names_from = "celltype",
      values_from = abundance
    )

  celltype_cols <- setdiff(colnames(wide_dat), c("Scan_ID", stratify_by))

  if (length(celltype_cols) < 2) {
    return(data.frame())
  }

  strata <- unique(wide_dat[[stratify_by]])
  result_list <- list()

  for (grp in strata) {

    sub <- wide_dat[wide_dat[[stratify_by]] == grp, celltype_cols, drop = FALSE]

    sub <- sub[, colSums(!is.na(sub)) > 0, drop = FALSE]

    if (nrow(sub) < min_n) next
    if (ncol(sub) < 2) next

    sub[is.na(sub)] <- 0

    clr_mat <- compute_clr(sub)

    pairs <- utils::combn(colnames(clr_mat), 2, simplify = FALSE)

    pair_results <- lapply(pairs, function(pair) {
      ct1 <- pair[1]
      ct2 <- pair[2]

      test <- try(
        stats::cor.test(clr_mat[[ct1]], clr_mat[[ct2]], method = "pearson"),
        silent = TRUE
      )

      if (inherits(test, "try-error")) return(NULL)

      data.frame(
        celltype_1 = ct1,
        celltype_2 = ct2,
        r = unname(test$estimate),
        p.value = test$p.value,
        n = nrow(clr_mat),
        stringsAsFactors = FALSE
      )
    })

    pair_df <- dplyr::bind_rows(pair_results)

    if (nrow(pair_df) == 0) next

    pair_df[[stratify_by]] <- grp
    result_list[[as.character(grp)]] <- pair_df
  }

  cooc_df <- dplyr::bind_rows(result_list)

  if (nrow(cooc_df) == 0) {
    return(cooc_df)
  }

  cooc_df <- cooc_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(stratify_by))) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(c(stratify_by, "celltype_1", "celltype_2", "r", "p.value", "p_adj", "n")))

  cooc_df
}
