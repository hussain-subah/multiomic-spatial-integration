# 04b — Extract Per-ROI Total Counts and Size Factors

# Separate from 04_extract_cell_proportions.py: that notebook needs both
# trained param stores; this one only needs the raw h5ad, so it can run
# independently and sooner. Produces the per-ROI total-count offset used by
# the absolute-abundance contrast pipeline (R/absolute_abundance_utils.R)
# and the calibration gate (R/absolute_calibration_utils.R).

import scanpy as sc

from python.abundance_extraction_utils import compute_roi_total_counts

adata = sc.read_h5ad("data/CAA-AD_AnnData.h5ad")

roi_total_counts = compute_roi_total_counts(adata)

print("ROI total counts preview:")
print(roi_total_counts.head())

print(
    "\nn_genes_used_for_size_factor:",
    roi_total_counts["n_genes_used_for_size_factor"].iloc[0],
    "/",
    adata.shape[1],
)

roi_total_counts.to_csv(
    "results/cell_proportions/roi_total_counts.csv",
    index=False,
)

print("\nSaved: results/cell_proportions/roi_total_counts.csv")
