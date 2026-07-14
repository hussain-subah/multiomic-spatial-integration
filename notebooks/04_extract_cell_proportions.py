# =============================================================================
# 04 — Extract and Annotate Spatial Cell-Type Proportions
#
# This notebook converts the validated SpaceJam outputs into analysis-ready
# abundance tables.
#
# Inputs:
#   results/spacejam/*.pt
#   results/regression_model/*_inferred_signatures.csv
#
# Outputs:
#   results/cell_proportions/
# =============================================================================

import scanpy as sc
import pandas as pd
import torch

from python.abundance_extraction_utils import (
    make_abundance_dataframe,
    make_long_abundance_table,
    save_abundance_outputs,
)

# =============================================================================
# Load GeoMx metadata
# =============================================================================

adata_wta = sc.read_h5ad("data/CAA-AD_AnnData.h5ad")

print(adata_wta)

# =============================================================================
# Load validated SpaceJam abundance tensors
# =============================================================================

spacejam_results = (
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/"
    "multiomic-spatial-integration/results/spacejam"
)

ad_spot_abs_t = torch.load(
    f"{spacejam_results}/ADCAA_spot_factors_abs.pt",
    map_location="cpu",
    weights_only=True,
)

ad_spot_rel_t = torch.load(
    f"{spacejam_results}/ADCAA_spot_factors_rel.pt",
    map_location="cpu",
    weights_only=True,
)

ctrl_spot_abs_t = torch.load(
    f"{spacejam_results}/CTRL_spot_factors_abs.pt",
    map_location="cpu",
    weights_only=True,
)

ctrl_spot_rel_t = torch.load(
    f"{spacejam_results}/CTRL_spot_factors_rel.pt",
    map_location="cpu",
    weights_only=True,
)

# Convert to NumPy

ad_spot_abs = ad_spot_abs_t.cpu().numpy()
ad_spot_rel = ad_spot_rel_t.cpu().numpy()

ctrl_spot_abs = ctrl_spot_abs_t.cpu().numpy()
ctrl_spot_rel = ctrl_spot_rel_t.cpu().numpy()

# =============================================================================
# Validate tensors
# =============================================================================

print("\n===== AD+CAA =====")
print("Shape:", ad_spot_abs.shape)
print("Absolute range:", ad_spot_abs.min(), ad_spot_abs.max())
print("Relative range:", ad_spot_rel.min(), ad_spot_rel.max())
print(
    "Row sums:",
    ad_spot_rel.sum(axis=1).min(),
    ad_spot_rel.sum(axis=1).max(),
)

print("\n===== Control =====")
print("Shape:", ctrl_spot_abs.shape)
print("Absolute range:", ctrl_spot_abs.min(), ctrl_spot_abs.max())
print("Relative range:", ctrl_spot_rel.min(), ctrl_spot_rel.max())
print(
    "Row sums:",
    ctrl_spot_rel.sum(axis=1).min(),
    ctrl_spot_rel.sum(axis=1).max(),
)

assert (ad_spot_abs > 0).all()
assert (ctrl_spot_abs > 0).all()

assert torch.allclose(
    ad_spot_rel_t.sum(dim=1),
    torch.ones(ad_spot_rel_t.shape[0]),
    atol=1e-5,
)

assert torch.allclose(
    ctrl_spot_rel_t.sum(dim=1),
    torch.ones(ctrl_spot_rel_t.shape[0]),
    atol=1e-5,
)

# Verify that relative abundances equal normalized absolute abundances

ad_check = ad_spot_abs_t / ad_spot_abs_t.sum(dim=1, keepdim=True)
ctrl_check = ctrl_spot_abs_t / ctrl_spot_abs_t.sum(dim=1, keepdim=True)

assert torch.allclose(ad_check, ad_spot_rel_t, atol=1e-6)
assert torch.allclose(ctrl_check, ctrl_spot_rel_t, atol=1e-6)

print("\nSpaceJam abundance tensors validated.")

# =============================================================================
# Load regression signatures
# =============================================================================

signature_dir = (
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/"
    "multiomic-spatial-integration/results/regression_model"
)

ad_signature = pd.read_csv(
    f"{signature_dir}/AD+CAA_inferred_signatures.csv",
    index_col=0,
)

ctrl_signature = pd.read_csv(
    f"{signature_dir}/Control_inferred_signatures.csv",
    index_col=0,
)

ad_labels = ad_signature.columns.tolist()
ctrl_labels = ctrl_signature.columns.tolist()

assert ad_labels == ctrl_labels

print("\nFactor order validation passed.")
print("Number of factors:", len(ad_labels))
print("First five:", ad_labels[:5])
print("Last five:", ad_labels[-5:])

# =============================================================================
# Validate factor dimensions
# =============================================================================

assert ad_spot_abs.shape[1] == len(ad_labels)
assert ctrl_spot_abs.shape[1] == len(ctrl_labels)

print("\nFactor dimensions validated.")

# =============================================================================
# Split metadata
# =============================================================================

ad_obs = (
    adata_wta.obs[
        adata_wta.obs["disease_status"] == "AD-CAA"
    ]
    .copy()
)

ctrl_obs = (
    adata_wta.obs[
        adata_wta.obs["disease_status"] == "Control"
    ]
    .copy()
)

assert ad_spot_abs.shape[0] == ad_obs.shape[0]
assert ctrl_spot_abs.shape[0] == ctrl_obs.shape[0]

print("ROI dimensions validated.")

# =============================================================================
# Build abundance dataframes
# =============================================================================

ad_abs_df = make_abundance_dataframe(
    spot_factors=ad_spot_abs,
    obs=ad_obs,
    celltype_labels=ad_labels,
)

ad_rel_df = make_abundance_dataframe(
    spot_factors=ad_spot_rel,
    obs=ad_obs,
    celltype_labels=ad_labels,
)

ctrl_abs_df = make_abundance_dataframe(
    spot_factors=ctrl_spot_abs,
    obs=ctrl_obs,
    celltype_labels=ctrl_labels,
)

ctrl_rel_df = make_abundance_dataframe(
    spot_factors=ctrl_spot_rel,
    obs=ctrl_obs,
    celltype_labels=ctrl_labels,
)

# =============================================================================
# Long-format relative abundance
# =============================================================================

ad_long_rel = make_long_abundance_table(
    abundance_df=ad_rel_df,
    celltype_labels=ad_labels,
    abundance_col="rel_abundance",
)

ctrl_long_rel = make_long_abundance_table(
    abundance_df=ctrl_rel_df,
    celltype_labels=ctrl_labels,
    abundance_col="rel_abundance",
)

abundance_long_rel = pd.concat(
    [ad_long_rel, ctrl_long_rel],
    ignore_index=True,
)

# =============================================================================
# Long-format absolute abundance
# =============================================================================

ad_long_abs = make_long_abundance_table(
    abundance_df=ad_abs_df,
    celltype_labels=ad_labels,
    abundance_col="abs_abundance",
)

ctrl_long_abs = make_long_abundance_table(
    abundance_df=ctrl_abs_df,
    celltype_labels=ctrl_labels,
    abundance_col="abs_abundance",
)

abundance_long_abs = pd.concat(
    [ad_long_abs, ctrl_long_abs],
    ignore_index=True,
)

# =============================================================================
# Final validation
# =============================================================================

print("\n===== Final validation =====")

print("Relative rows:", abundance_long_rel.shape)
print("Absolute rows:", abundance_long_abs.shape)

print(
    "Unique ROIs:",
    abundance_long_rel["ROI_ID"].nunique(),
)

print(
    "Unique cell types:",
    abundance_long_rel["celltype"].nunique(),
)

assert (
    abundance_long_rel["ROI_ID"].nunique()
    == adata_wta.n_obs
)

assert (
    abundance_long_rel["celltype"].nunique()
    == len(ad_labels)
)

print("Output tables validated.")

print("\nRelative abundance preview")
print(abundance_long_rel.head())

print("\nAbsolute abundance preview")
print(abundance_long_abs.head())

# ============================================================
# Create combined wide-format table (one row per ROI)
# ============================================================

wide_df = pd.concat(
    [ad_rel_df, ctrl_rel_df],
    axis=0,
    ignore_index=True
)

# =============================================================================
# Save outputs
# =============================================================================

output_dir = "results/cell_proportions"

save_abundance_outputs(
    output_dir=output_dir,
    ad_abs_df=ad_abs_df,
    ad_rel_df=ad_rel_df,
    ctrl_abs_df=ctrl_abs_df,
    ctrl_rel_df=ctrl_rel_df,
    long_df=abundance_long_rel,
)

abundance_long_rel.to_csv(
    f"{output_dir}/roi_celltype_abundance_long.csv",
    index=False,
)

abundance_long_abs.to_csv(
    f"{output_dir}/roi_celltype_abundance_long_abs.csv",
    index=False,
)

wide_df.to_csv(
    "results/cell_proportions/cell2location_abundance_wide_with_meta.csv",
    index=False
)
abundance_long_rel.to_csv(
    "results/cell_proportions/cell2location_abundance_long_rel.csv",
    index=False
)

# This is the file consumed by the downstream R notebooks.

abundance_long_rel.to_csv(
    f"{output_dir}/spatial_celltype_proportions_for_R.csv",
    index=False,
)

print("\nSaved:")
print(f"  {output_dir}/roi_celltype_abundance_long.csv")
print(f"  {output_dir}/roi_celltype_abundance_long_abs.csv")
print(f"  {output_dir}/spatial_celltype_proportions_for_R.csv")

print("\nNotebook 4 completed successfully.")
