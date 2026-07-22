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

#' Safe correlation for method-comparison summaries
#'
#' Returns NA when there are too few complete observations or either
#' vector has zero variance.
safe_method_cor <- function(x, y, method = "spearman", min_n = 3) {
  keep <- is.finite(x) & is.finite(y)

  x <- x[keep]
  y <- y[keep]

  if (length(x) < min_n) {
    return(NA_real_)
  }

  if (
    stats::sd(x) == 0 ||
    stats::sd(y) == 0
  ) {
    return(NA_real_)
  }

  stats::cor(
    x,
    y,
    method = method
  )
}


#' Map fine cell states to broad biological lineages
#'
#' @param celltype Character vector of cell-type labels.
#'
#' @return Character vector of broad lineage labels.
map_celltype_to_lineage <- function(celltype) {
  dplyr::case_when(
    grepl("^Astrocytes", celltype) ~ "Astrocytes",
    grepl("^Microglia", celltype) ~ "Microglia",
    grepl("^Oligodendrocytes", celltype) ~ "Oligodendrocytes",
    grepl("^OPC", celltype) ~ "OPC",
    grepl("^ExNeuron", celltype) ~ "Excitatory neurons",
    grepl("^InhNeuron", celltype) ~ "Inhibitory neurons",
    celltype == "Endothelial" ~ "Endothelial",
    celltype == "Pericytes" ~ "Pericytes",
    celltype == "SMC" ~ "SMC",
    grepl("^VLMC", celltype) ~ "VLMC",
    celltype == "Fibroblast" ~ "Fibroblast",
    TRUE ~ celltype
  )
}


#' Summarize how often a deconvolution method uses each factor
#'
#' @param method_df Standard long table with ROI_ID, celltype, proportion.
#'
#' @return One row per cell type.
summarize_factor_usage <- function(method_df) {
  method_df %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_evaluable = sum(
        !is.na(proportion)
      ),
      n_nonzero = sum(
        proportion > 0,
        na.rm = TRUE
      ),
      fraction_nonzero = dplyr::if_else(
        n_evaluable > 0,
        n_nonzero / n_evaluable,
        NA_real_
      ),
      mean_proportion = mean(
        proportion,
        na.rm = TRUE
      ),
      median_proportion = stats::median(
        proportion,
        na.rm = TRUE
      ),
      max_proportion = max(
        proportion,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      usage_status = dplyr::case_when(
        n_evaluable == 0 ~ "unevaluable",
        n_nonzero == 0 ~ "all_zero",
        fraction_nonzero < 0.10 ~ "rare",
        TRUE ~ "used"
      )
    ) %>%
    dplyr::arrange(
      dplyr::desc(fraction_nonzero),
      dplyr::desc(mean_proportion)
    )
}


#' Identify ROIs that could not be deconvolved
#'
#' @param method_df Standard long table with ROI_ID and proportion.
#'
#' @return One row per failed ROI.
identify_failed_rois <- function(method_df) {
  method_df %>%
    dplyr::group_by(ROI_ID) %>%
    dplyr::summarise(
      all_na = all(
        is.na(proportion)
      ),
      n_finite = sum(
        is.finite(proportion)
      ),
      proportion_sum = if (
        all(is.na(proportion))
      ) {
        NA_real_
      } else {
        sum(
          proportion,
          na.rm = TRUE
        )
      },
      .groups = "drop"
    ) %>%
    dplyr::filter(
      all_na |
      !is.finite(proportion_sum) |
      proportion_sum <= 0
    )
}


#' Compare two methods at the fine cell-type level
#'
#' @param method_df Standard long table from the comparison method.
#' @param baseline_df Standard long Cell2location table.
#'
#' @return One row per cell type.
compare_methods_by_celltype <- function(
    method_df,
    baseline_df
) {
  method_clean <- method_df %>%
    dplyr::transmute(
      ROI_ID,
      celltype,
      comparison_method = proportion
    )

  baseline_clean <- baseline_df %>%
    dplyr::transmute(
      ROI_ID,
      celltype,
      baseline = proportion
    )

  joined <- dplyr::inner_join(
    method_clean,
    baseline_clean,
    by = c(
      "ROI_ID",
      "celltype"
    )
  )

  joined %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(
      n_comparable = sum(
        is.finite(comparison_method) &
        is.finite(baseline)
      ),
      n_method_nonzero = sum(
        comparison_method > 0,
        na.rm = TRUE
      ),
      method_fraction_nonzero = dplyr::if_else(
        n_comparable > 0,
        n_method_nonzero / n_comparable,
        NA_real_
      ),
      spearman_rho = safe_method_cor(
        comparison_method,
        baseline,
        method = "spearman"
      ),
      pearson_r = safe_method_cor(
        comparison_method,
        baseline,
        method = "pearson"
      ),
      mean_abs_difference = mean(
        abs(
          comparison_method -
          baseline
        ),
        na.rm = TRUE
      ),
      mean_comparison_method = mean(
        comparison_method,
        na.rm = TRUE
      ),
      mean_baseline = mean(
        baseline,
        na.rm = TRUE
      ),
      comparison_status = dplyr::case_when(
        n_comparable < 3 ~ "insufficient_ROIs",
        n_method_nonzero == 0 ~ "method_all_zero",
        method_fraction_nonzero < 0.10 ~ "method_rare",
        TRUE ~ "evaluable"
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      factor(
        comparison_status,
        levels = c(
          "evaluable",
          "method_rare",
          "method_all_zero",
          "insufficient_ROIs"
        )
      ),
      dplyr::desc(spearman_rho)
    )
}


#' Compare two methods after aggregation to broad lineages
#'
#' @param method_df Standard long table from the comparison method.
#' @param baseline_df Standard long Cell2location table.
#'
#' @return One row per broad lineage.
compare_methods_by_lineage <- function(
    method_df,
    baseline_df
) {
  method_lineage <- method_df %>%
    dplyr::mutate(
      lineage = map_celltype_to_lineage(
        celltype
      )
    ) %>%
    dplyr::group_by(
      ROI_ID,
      lineage
    ) %>%
    dplyr::summarise(
      comparison_method = if (
        all(is.na(proportion))
      ) {
        NA_real_
      } else {
        sum(
          proportion,
          na.rm = TRUE
        )
      },
      .groups = "drop"
    )

  baseline_lineage <- baseline_df %>%
    dplyr::mutate(
      lineage = map_celltype_to_lineage(
        celltype
      )
    ) %>%
    dplyr::group_by(
      ROI_ID,
      lineage
    ) %>%
    dplyr::summarise(
      baseline = if (
        all(is.na(proportion))
      ) {
        NA_real_
      } else {
        sum(
          proportion,
          na.rm = TRUE
        )
      },
      .groups = "drop"
    )

  joined <- dplyr::inner_join(
    method_lineage,
    baseline_lineage,
    by = c(
      "ROI_ID",
      "lineage"
    )
  )

  joined %>%
    dplyr::group_by(lineage) %>%
    dplyr::summarise(
      n_comparable = sum(
        is.finite(comparison_method) &
        is.finite(baseline)
      ),
      method_fraction_nonzero = mean(
        comparison_method > 0,
        na.rm = TRUE
      ),
      spearman_rho = safe_method_cor(
        comparison_method,
        baseline,
        method = "spearman"
      ),
      pearson_r = safe_method_cor(
        comparison_method,
        baseline,
        method = "pearson"
      ),
      mean_abs_difference = mean(
        abs(
          comparison_method -
          baseline
        ),
        na.rm = TRUE
      ),
      mean_comparison_method = mean(
        comparison_method,
        na.rm = TRUE
      ),
      mean_baseline = mean(
        baseline,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(spearman_rho)
    )
}
