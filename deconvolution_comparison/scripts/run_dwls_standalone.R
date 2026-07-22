suppressPackageStartupMessages({
  library(Matrix)
  library(DWLS)
})

source(
  "deconvolution_comparison/R/common_utils.R"
)

source(
  "deconvolution_comparison/R/adapter_dwls.R"
)

project_dir <- "/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration"

output_dir <- Sys.getenv(
  "DWLS_OUTPUT_DIR",
  unset = file.path(
    project_dir,
    "results",
    "deconvolution_comparison"
  )
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

n_cores <- as.integer(
  Sys.getenv(
    "SLURM_CPUS_PER_TASK",
    unset = "8"
  )
)

message("Loading independent reference signature...")

signature <- load_signature_matrix(
  file.path(
    output_dir,
    "independent_reference_signature.csv"
  )
)

message(
  "Signature: ",
  nrow(signature),
  " genes x ",
  ncol(signature),
  " cell types."
)

message("Loading GeoMx mixture...")

mixture <- load_mixture(
  "results/geomx_exports/CAA-AD_expression_wide.csv"
)

message(
  "Mixture: ",
  nrow(mixture),
  " genes x ",
  ncol(mixture),
  " ROIs."
)

message("Running DWLS...")

start_time <- Sys.time()

dwls_prop <- run_dwls(
  signature = signature,
  mixture = mixture,
  n_cores = n_cores
)

elapsed <- difftime(
  Sys.time(),
  start_time,
  units = "mins"
)

message(
  "DWLS elapsed time: ",
  round(
    as.numeric(elapsed),
    2
  ),
  " minutes."
)

if (is.null(dwls_prop)) {
  stop("DWLS returned NULL.")
}

raw_sums <- rowSums(
  dwls_prop,
  na.rm = TRUE
)

message(
  "DWLS result: ",
  nrow(dwls_prop),
  " ROIs x ",
  ncol(dwls_prop),
  " cell types."
)

message(
  "Raw value range: ",
  paste(
    range(
      dwls_prop,
      na.rm = TRUE
    ),
    collapse = " to "
  )
)

message(
  "Zero or failed ROI rows: ",
  sum(
    !is.finite(raw_sums) |
    raw_sums <= 0
  )
)

message(
  "Negative coefficients: ",
  sum(
    dwls_prop < 0,
    na.rm = TRUE
  )
)

stopifnot(
  nrow(dwls_prop) == 190,
  ncol(dwls_prop) == 46,
  !anyDuplicated(rownames(dwls_prop)),
  !anyDuplicated(colnames(dwls_prop))
)

long <- standardize_proportions(
  prop_mat = dwls_prop,
  method = "DWLS",
  normalize = TRUE
)

save_method_output(
  long_df = long,
  method = "DWLS",
  output_dir = output_dir
)

message(
  "Standardized output: ",
  nrow(long),
  " rows; ",
  length(unique(long$ROI_ID)),
  " ROIs; ",
  length(unique(long$celltype)),
  " cell types."
)

message("DWLS run completed successfully.")
