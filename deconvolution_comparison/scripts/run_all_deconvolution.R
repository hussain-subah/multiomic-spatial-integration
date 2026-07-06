# ============================================================
# Run all R-based deconvolution adapters on the GeoMx ROI data
#
# Each adapter runs only if its package is installed (a missing package is
# skipped with a message, not an error), so one absent dependency doesn't
# block the rest. Outputs are written in the standard long format to
# results/deconvolution_comparison/<method>_proportions.csv.
#
# CIBERSORTx is not run here (license/web submission required) -- see
# scripts/prepare_cibersortx_inputs.R.
# ============================================================

comparison_root <- "deconvolution_comparison/R"

source(file.path(comparison_root, "common_utils.R"))
source(file.path(comparison_root, "adapter_music.R"))
source(file.path(comparison_root, "adapter_bisque.R"))
source(file.path(comparison_root, "adapter_dwls.R"))
source(file.path(comparison_root, "adapter_spotlight.R"))
source(file.path(comparison_root, "adapter_rctd.R"))
source(file.path(comparison_root, "adapter_bayesprism.R"))
source(file.path(comparison_root, "adapter_stdeconvolve.R"))
source(file.path(comparison_root, "adapter_spatialdecon.R"))

# ============================================================
# Config -- adjust paths for your environment
# ============================================================

reference_dir <- "results/sn_reference/anndata_inputs"
counts_mtx <- file.path(reference_dir, "AD_CAA_counts.mtx")
metadata_csv <- file.path(reference_dir, "AD_CAA_metadata.csv")
genes_csv <- file.path(reference_dir, "AD_CAA_genes.csv")

mixture_csv <- "data/CAA-AD_expression_wide.csv"
background_csv <- "data/CAA-AD_negprobe_background.csv"  # for SpatialDecon

output_dir <- "results/deconvolution_comparison"

celltype_col <- "New_Idents"   # cell-type label column in the reference metadata
sample_col <- "orig.ident"     # biological-sample column (for MuSiC/Bisque)

# Signature source for the signature-based tools (DWLS, SpatialDecon):
#   "independent" -> per-cell-type mean aggregated from the raw reference
#                    counts (build_mean_signature); NOT produced by
#                    cell2location, so it does not share a basis with the
#                    Cell2location/SpaceJam baseline (removes circularity).
#   "inferred"    -> the cell2location NB-regression signatures
#                    (AD+CAA_inferred_signatures.csv); shares a basis with the
#                    baseline, so agreement is partly built-in.
signature_source <- "independent"
inferred_signature_csv <- "data/Regression-model/AD+CAA_inferred_signatures.csv"
independent_signature_out <- "data/independent_reference_signature.csv"

# ============================================================
# Load shared inputs once
# ============================================================

message("Loading reference...")
ref <- load_reference_mtx(
  counts_mtx = counts_mtx,
  metadata_csv = metadata_csv,
  genes_csv = genes_csv,
  celltype_col = celltype_col,
  sample_col = sample_col
)

message("Loading mixture...")
mixture <- load_mixture(mixture_csv)

if (signature_source == "independent") {
  message("Building independent per-cell-type signature from raw reference counts...")
  signature <- build_mean_signature(ref)
  # Save it so prepare_cibersortx_inputs.R (and reruns) can reuse the same
  # independent signature rather than the cell2location one.
  utils::write.csv(signature, independent_signature_out)
  message("Wrote independent signature to ", independent_signature_out)
} else {
  message("Loading cell2location inferred signature matrix...")
  signature <- load_signature_matrix(inferred_signature_csv)
}

message("Loading negative-probe background (for SpatialDecon)...")
background <- load_roi_background(background_csv)

# ============================================================
# Run each adapter, standardize, and save
# ============================================================

run_and_save <- function(method, prop_mat, normalize = TRUE) {
  if (is.null(prop_mat)) return(invisible(NULL))
  long <- standardize_proportions(prop_mat, method = method, normalize = normalize)
  save_method_output(long, method = method, output_dir = output_dir)
}

message("\n== MuSiC ==")
run_and_save("MuSiC", tryCatch(run_music(ref, mixture),
             error = function(e) { message("MuSiC failed: ", conditionMessage(e)); NULL }))

message("\n== Bisque ==")
run_and_save("Bisque", tryCatch(run_bisque(ref, mixture),
             error = function(e) { message("Bisque failed: ", conditionMessage(e)); NULL }))

message("\n== DWLS ==")
run_and_save("DWLS", tryCatch(run_dwls(signature, mixture),
             error = function(e) { message("DWLS failed: ", conditionMessage(e)); NULL }))

message("\n== SPOTlight ==")
run_and_save("SPOTlight", tryCatch(run_spotlight(ref, mixture),
             error = function(e) { message("SPOTlight failed: ", conditionMessage(e)); NULL }))

message("\n== RCTD ==")
run_and_save("RCTD", tryCatch(run_rctd(ref, mixture),
             error = function(e) { message("RCTD failed: ", conditionMessage(e)); NULL }))

message("\n== BayesPrism ==")
run_and_save("BayesPrism", tryCatch(run_bayesprism(ref, mixture),
             error = function(e) { message("BayesPrism failed: ", conditionMessage(e)); NULL }))

message("\n== STdeconvolve ==")
run_and_save("STdeconvolve", tryCatch(run_stdeconvolve(mixture, signature),
             error = function(e) { message("STdeconvolve failed: ", conditionMessage(e)); NULL }))

message("\n== SpatialDecon ==")
run_and_save("SpatialDecon", tryCatch(run_spatialdecon(mixture, signature, background),
             error = function(e) { message("SpatialDecon failed: ", conditionMessage(e)); NULL }))

message("\nAll available deconvolution adapters completed.")
