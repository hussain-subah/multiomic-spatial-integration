#' DWLS deconvolution adapter
#'
#' Runs Dampened Weighted Least Squares independently for each ROI using an
#' independent cell-type signature derived from the raw snRNA-seq reference.
#' ROI fits can be parallelized on Unix-like systems.
#'
#' @param signature Gene x cell-type signature matrix.
#' @param mixture Gene x ROI expression matrix.
#' @param n_cores Number of CPU cores used across ROIs.
#'
#' @return ROI x cell-type proportion matrix, or NULL if DWLS is unavailable.
#' @export
run_dwls <- function(
    signature,
    mixture,
    n_cores = 1
) {
  if (!require_pkg("DWLS", "DWLS")) {
    return(NULL)
  }

  shared <- intersect(
    rownames(signature),
    rownames(mixture)
  )

  if (length(shared) < 100) {
    warning(
      "run_dwls: only ",
      length(shared),
      " shared genes between signature and mixture.",
      call. = FALSE
    )
  } else {
    message(
      "run_dwls: ",
      length(shared),
      " shared genes."
    )
  }

  signature <- signature[
    shared,
    ,
    drop = FALSE
  ]

  mixture <- mixture[
    shared,
    ,
    drop = FALSE
  ]

  roi_ids <- colnames(mixture)
  celltypes <- colnames(signature)

  fit_one_roi <- function(i) {
    roi_id <- roi_ids[i]
    bulk_sample <- mixture[, i]

    tr <- DWLS::trimData(
      signature,
      bulk_sample
    )

    sol <- NULL

    # DWLS prints each solution internally; capture it to avoid huge logs.
    invisible(
      utils::capture.output({
        sol <- try(
          DWLS::solveDampenedWLS(
            tr$sig,
            tr$bulk
          ),
          silent = TRUE
        )
      })
    )

    used_fallback <- FALSE

    if (inherits(sol, "try-error")) {
      used_fallback <- TRUE

      invisible(
        utils::capture.output({
          sol <- try(
            DWLS::solveOLS(
              tr$sig,
              tr$bulk
            ),
            silent = TRUE
          )
        })
      )
    }

    if (inherits(sol, "try-error")) {
      return(
        list(
          roi_id = roi_id,
          solution = rep(
            NA_real_,
            length(celltypes)
          ),
          fallback = used_fallback,
          failed = TRUE
        )
      )
    }

    original_names <- names(sol)
    sol <- as.numeric(sol)

    if (
      !is.null(original_names) &&
      length(original_names) == length(sol)
    ) {
      names(sol) <- original_names
    } else {
      names(sol) <- celltypes
    }

    complete_sol <- stats::setNames(
      rep(0, length(celltypes)),
      celltypes
    )

    matched <- intersect(
      names(sol),
      celltypes
    )

    complete_sol[matched] <- sol[matched]

    # Remove floating-point negative noise.
    complete_sol[
      complete_sol < 0 &
      complete_sol > -1e-8
    ] <- 0

    if (any(complete_sol < 0)) {
      complete_sol[complete_sol < 0] <- 0
    }

    if (
      all(is.finite(complete_sol)) &&
      sum(complete_sol) > 0
    ) {
      complete_sol <- complete_sol /
        sum(complete_sol)
    } else {
      complete_sol[] <- NA_real_
    }

    list(
      roi_id = roi_id,
      solution = complete_sol,
      fallback = used_fallback,
      failed = all(is.na(complete_sol))
    )
  }

  n_cores <- max(
    1L,
    min(
      as.integer(n_cores),
      length(roi_ids)
    )
  )

  message(
    "run_dwls: fitting ",
    length(roi_ids),
    " ROIs using ",
    n_cores,
    " CPU core(s)."
  )

  if (
    n_cores > 1 &&
    .Platform$OS.type == "unix"
  ) {
    fits <- parallel::mclapply(
      seq_along(roi_ids),
      fit_one_roi,
      mc.cores = n_cores,
      mc.preschedule = TRUE
    )
  } else {
    fits <- lapply(
      seq_along(roi_ids),
      fit_one_roi
    )
  }

  prop_mat <- do.call(
    rbind,
    lapply(
      fits,
      function(x) x$solution
    )
  )

  rownames(prop_mat) <- vapply(
    fits,
    function(x) x$roi_id,
    character(1)
  )

  colnames(prop_mat) <- celltypes

  fallback_rois <- vapply(
    fits,
    function(x) x$fallback,
    logical(1)
  )

  failed_rois <- vapply(
    fits,
    function(x) x$failed,
    logical(1)
  )

  message(
    "run_dwls: OLS fallback used for ",
    sum(fallback_rois),
    " ROI(s)."
  )

  message(
    "run_dwls: failed ROIs: ",
    sum(failed_rois)
  )

  if (any(failed_rois)) {
    warning(
      "run_dwls: failed ROI IDs: ",
      paste(
        head(
          rownames(prop_mat)[failed_rois],
          10
        ),
        collapse = ", "
      ),
      if (sum(failed_rois) > 10) " ..." else "",
      call. = FALSE
    )
  }

  prop_mat
}
