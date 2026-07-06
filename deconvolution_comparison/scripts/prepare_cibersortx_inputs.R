# ============================================================
# Prepare CIBERSORTx input files
#
# CIBERSORTx cannot be run as a plain function call -- it requires a license
# token and is run via its web portal (cibersortx.stanford.edu) or a licensed
# Docker container. This script only writes the two input files in the
# tab-separated format CIBERSORTx expects; submit them separately, then drop
# the returned proportions into the standard output location (see the note at
# the end).
# ============================================================

source("deconvolution_comparison/R/common_utils.R")

# ============================================================
# Config -- adjust paths for your environment
# ============================================================

mixture_csv <- "data/CAA-AD_expression_wide.csv"
# Default to the independent raw-reference signature written by
# run_all_deconvolution.R (removes circularity with the Cell2location
# baseline). Switch to the cell2location inferred signatures
# ("data/Regression-model/AD+CAA_inferred_signatures.csv") if you want the
# shared-basis version instead.
signature_csv <- "data/independent_reference_signature.csv"
output_dir <- "results/deconvolution_comparison/cibersortx_inputs"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Signature matrix: gene (col 1) x celltype, tab-separated
# ============================================================

signature <- load_signature_matrix(signature_csv)

sig_out <- data.frame(
  GeneSymbol = rownames(signature),
  signature,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  sig_out,
  file = file.path(output_dir, "cibersortx_signature_matrix.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Mixture: gene (col 1) x ROI, tab-separated
# ============================================================

mixture <- load_mixture(mixture_csv)  # gene x ROI

mix_out <- data.frame(
  GeneSymbol = rownames(mixture),
  mixture,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  mix_out,
  file = file.path(output_dir, "cibersortx_mixture.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote CIBERSORTx inputs to ", output_dir)
message(
  "Next steps (manual):\n",
  "  1. Submit these to CIBERSORTx (web portal or licensed Docker) in\n",
  "     'Impute Cell Fractions' mode, using the signature matrix as the\n",
  "     signature and the mixture file as the mixture.\n",
  "  2. Download the resulting fractions table.\n",
  "  3. Reshape it to the standard long format (method, ROI_ID, celltype,\n",
  "     proportion) and save as\n",
  "     results/deconvolution_comparison/CIBERSORTx_proportions.csv\n",
  "     so run_deconvolution_comparison.R picks it up automatically."
)
