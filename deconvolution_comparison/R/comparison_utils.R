#' Concordance utilities for the deconvolution method comparison
#'
#' Consumes the standard long-format proportion tables written by each
#' adapter (plus the Cell2location baseline) and quantifies agreement. There
#' is no ground truth for the real ROIs, so these measure agreement, not
#' accuracy.
#'
#' @keywords internal
NULL


#' Load all standardized method-output tables plus the baseline
#'
#' @param output_dir Directory holding `<method>_proportions.csv` files.
#' @param baseline_long Optional long baseline table (from
#'   `load_cell2location_baseline()`) to bind in.
#'
#' @return Combined long data frame: method, ROI_ID, celltype, proportion.
#' @export
load_all_method_outputs <- function(output_dir, baseline_long = NULL) {
  files <- list.files(output_dir, pattern = "_proportions\\.csv$", full.names = TRUE)

  method_dfs <- lapply(files, utils::read.csv, stringsAsFactors = FALSE)
  combined <- do.call(rbind, method_dfs)

  if (!is.null(baseline_long)) {
    combined <- rbind(combined, baseline_long)
  }

  combined
}


#' Per-celltype concordance of each method against a baseline
#'
#' For every (method, celltype), correlates that method's per-ROI proportion
#' against the baseline method's, across shared ROIs.
#'
#' @param long_df Combined long table from `load_all_method_outputs()`.
#' @param baseline_method Method name to treat as the reference.
#' @param method_col,roi_col,celltype_col,value_col Column names.
#' @param cor_method Correlation method.
#'
#' @return Data frame: method, celltype, pearson_r, spearman_rho, rmse, n_roi.
#' @export
concordance_vs_baseline <- function(long_df,
                                    baseline_method = "Cell2location",
                                    method_col = "method",
                                    roi_col = "ROI_ID",
                                    celltype_col = "celltype",
                                    value_col = "proportion",
                                    cor_method = "pearson") {
  baseline <- long_df[long_df[[method_col]] == baseline_method, ]

  if (nrow(baseline) == 0) {
    stop("concordance_vs_baseline: no rows for baseline_method '", baseline_method, "'.",
         call. = FALSE)
  }

  other_methods <- setdiff(unique(long_df[[method_col]]), baseline_method)
  results <- list()

  for (m in other_methods) {
    method_df <- long_df[long_df[[method_col]] == m, ]
    shared_ct <- intersect(unique(method_df[[celltype_col]]), unique(baseline[[celltype_col]]))

    for (ct in shared_ct) {
      b <- baseline[baseline[[celltype_col]] == ct, c(roi_col, value_col)]
      x <- method_df[method_df[[celltype_col]] == ct, c(roi_col, value_col)]

      merged <- merge(b, x, by = roi_col, suffixes = c("_base", "_method"))

      if (nrow(merged) < 3) next

      base_vals <- merged[[paste0(value_col, "_base")]]
      method_vals <- merged[[paste0(value_col, "_method")]]

      pear <- if (stats::sd(base_vals) > 0 && stats::sd(method_vals) > 0) {
        stats::cor(base_vals, method_vals, method = cor_method)
      } else NA_real_

      spear <- if (stats::sd(base_vals) > 0 && stats::sd(method_vals) > 0) {
        stats::cor(base_vals, method_vals, method = "spearman")
      } else NA_real_

      results[[paste(m, ct, sep = "__")]] <- data.frame(
        method = m,
        celltype = ct,
        pearson_r = pear,
        spearman_rho = spear,
        rmse = sqrt(mean((base_vals - method_vals)^2)),
        n_roi = nrow(merged),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, results)
  out[order(out$method, -out$pearson_r), ]
}


#' Overall method-vs-method correlation matrix
#'
#' Flattens each method's full ROI x celltype proportion vector (over the
#' shared ROI/celltype grid) and correlates every pair of methods.
#'
#' @param long_df Combined long table.
#' @param method_col,roi_col,celltype_col,value_col Column names.
#' @param cor_method Correlation method.
#'
#' @return Symmetric numeric matrix, method x method.
#' @export
cross_method_correlation <- function(long_df,
                                     method_col = "method",
                                     roi_col = "ROI_ID",
                                     celltype_col = "celltype",
                                     value_col = "proportion",
                                     cor_method = "pearson") {
  long_df$key <- paste(long_df[[roi_col]], long_df[[celltype_col]], sep = "||")

  methods_vec <- unique(long_df[[method_col]])

  wide <- Reduce(
    function(acc, m) {
      sub <- long_df[long_df[[method_col]] == m, c("key", value_col)]
      colnames(sub)[2] <- m
      if (is.null(acc)) sub else merge(acc, sub, by = "key", all = TRUE)
    },
    methods_vec,
    accumulate = FALSE,
    right = FALSE,
    init = NULL
  )

  mat <- as.matrix(wide[, methods_vec, drop = FALSE])

  stats::cor(mat, use = "pairwise.complete.obs", method = cor_method)
}


#' Mean absolute proportion difference between each method and the baseline
#'
#' A per-method single-number summary, complementary to the per-celltype
#' correlations: the mean |method - baseline| proportion over all shared
#' (ROI, celltype) cells.
#'
#' @param long_df Combined long table.
#' @param baseline_method Baseline method name.
#' @param method_col,roi_col,celltype_col,value_col Column names.
#'
#' @return Data frame: method, mean_abs_diff, n_cells.
#' @export
mean_abs_diff_vs_baseline <- function(long_df,
                                      baseline_method = "Cell2location",
                                      method_col = "method",
                                      roi_col = "ROI_ID",
                                      celltype_col = "celltype",
                                      value_col = "proportion") {
  baseline <- long_df[long_df[[method_col]] == baseline_method, ]
  baseline$key <- paste(baseline[[roi_col]], baseline[[celltype_col]], sep = "||")

  other_methods <- setdiff(unique(long_df[[method_col]]), baseline_method)
  results <- list()

  for (m in other_methods) {
    md <- long_df[long_df[[method_col]] == m, ]
    md$key <- paste(md[[roi_col]], md[[celltype_col]], sep = "||")

    merged <- merge(
      baseline[, c("key", value_col)],
      md[, c("key", value_col)],
      by = "key",
      suffixes = c("_base", "_method")
    )

    if (nrow(merged) == 0) next

    diffs <- abs(merged[[paste0(value_col, "_base")]] - merged[[paste0(value_col, "_method")]])

    results[[m]] <- data.frame(
      method = m,
      mean_abs_diff = mean(diffs),
      n_cells = nrow(merged),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, results)
  out[order(out$mean_abs_diff), ]
}
