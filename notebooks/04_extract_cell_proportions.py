# 04 — Extract and Annotate Spatial Cell-Type Proportions

#This notebook extracts posterior cell-type abundance estimates from trained SpaceJam/Pyro models and converts them into analysis-ready tables.

#Outputs include:

#- absolute cell-type abundance per ROI
#- relative cell-type proportions per ROI
#- long-format ROI × cell-type tables
#- metadata-annotated abundance tables

import scanpy as sc
import pandas as pd
import torch

from python.abundance_extraction_utils import (
    load_variational_means,
    extract_spot_factors,
    normalize_spot_factors,
    make_abundance_dataframe,
    make_long_abundance_table,
    save_abundance_outputs
)
## Load GeoMx AnnData metadata

adata_wta = sc.read_h5ad("data/CAA-AD_AnnData.h5ad")
adata_wta.obs.head()

## Load posterior means from trained spatial models
## Extract spot factors
#`spot_factors` represent inferred absolute cell-type abundance per ROI.
## Normalize to relative proportions
## Load SpaceJam posterior parameters and recover constrained abundances

spacejam_results = (
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/SpaceJam_results"
)

ad_param_store = torch.load(
    f"{spacejam_results}/ADCAA_param_store.pt",
    map_location="cpu",
)

ctrl_param_store = torch.load(
    f"{spacejam_results}/CTRL_param_store.pt",
    map_location="cpu",
)

spot_factor_key = "AutoNormal.locs.spot_factors"

if spot_factor_key not in ad_param_store:
    raise KeyError(
        f"{spot_factor_key} not found in AD+CAA parameter store. "
        f"Available matching keys: "
        f"{[k for k in ad_param_store if 'spot_factors' in k]}"
    )

if spot_factor_key not in ctrl_param_store:
    raise KeyError(
        f"{spot_factor_key} not found in Control parameter store. "
        f"Available matching keys: "
        f"{[k for k in ctrl_param_store if 'spot_factors' in k]}"
    )

# AutoNormal stores positive latent variables in unconstrained log space.
# Transform back to the constrained positive abundance scale.
ad_spot_abs_t = torch.exp(
    ad_param_store[spot_factor_key].detach().cpu()
)

ctrl_spot_abs_t = torch.exp(
    ctrl_param_store[spot_factor_key].detach().cpu()
)

# Convert absolute abundances into within-ROI relative proportions.
ad_spot_rel_t = ad_spot_abs_t / ad_spot_abs_t.sum(
    dim=1,
    keepdim=True,
)

ctrl_spot_rel_t = ctrl_spot_abs_t / ctrl_spot_abs_t.sum(
    dim=1,
    keepdim=True,
)

# Convert to NumPy arrays for the existing helper functions.
ad_spot_abs = ad_spot_abs_t.numpy()
ad_spot_rel = ad_spot_rel_t.numpy()

ctrl_spot_abs = ctrl_spot_abs_t.numpy()
ctrl_spot_rel = ctrl_spot_rel_t.numpy()

## Sanity checks

print(
    "AD+CAA absolute abundance range:",
    ad_spot_abs.min(),
    ad_spot_abs.max(),
)

print(
    "Control absolute abundance range:",
    ctrl_spot_abs.min(),
    ctrl_spot_abs.max(),
)

print(
    "AD+CAA first relative row sums:",
    ad_spot_rel.sum(axis=1)[:5],
)

print(
    "Control first relative row sums:",
    ctrl_spot_rel.sum(axis=1)[:5],
)

assert (ad_spot_abs > 0).all()
assert (ctrl_spot_abs > 0).all()

assert torch.allclose(
    torch.tensor(ad_spot_rel.sum(axis=1)),
    torch.ones(ad_spot_rel.shape[0]),
    atol=1e-5,
)

assert torch.allclose(
    torch.tensor(ctrl_spot_rel.sum(axis=1)),
    torch.ones(ctrl_spot_rel.shape[0]),
    atol=1e-5,
)

## Save corrected constrained-space abundance tensors

torch.save(
    ad_spot_abs_t,
    f"{spacejam_results}/ADCAA_spot_factors_abs.pt",
)

torch.save(
    ad_spot_rel_t,
    f"{spacejam_results}/ADCAA_spot_factors_rel.pt",
)

torch.save(
    ctrl_spot_abs_t,
    f"{spacejam_results}/CTRL_spot_factors_abs.pt",
)

torch.save(
    ctrl_spot_rel_t,
    f"{spacejam_results}/CTRL_spot_factors_rel.pt",
)

## Add cell-type labels

#These labels should match the regression-derived signature columns.
## Load condition-specific cell-type labels

signature_dir = (
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/Regression-model"
)

ad_signature = pd.read_csv(
    f"{signature_dir}/AD+CAA_inferred_signatures.csv",
    index_col=0
)

ctrl_signature = pd.read_csv(
    f"{signature_dir}/Control_inferred_signatures.csv",
    index_col=0
)

ad_labels = ad_signature.columns.tolist()
ctrl_labels = ctrl_signature.columns.tolist()

# Sanity checks: factor counts must match signature columns
if ad_spot_abs.shape[1] != len(ad_labels):
    raise ValueError(
        f"AD+CAA factor count mismatch: "
        f"{ad_spot_abs.shape[1]} spot-factor columns vs "
        f"{len(ad_labels)} signature labels."
    )

if ctrl_spot_abs.shape[1] != len(ctrl_labels):
    raise ValueError(
        f"Control factor count mismatch: "
        f"{ctrl_spot_abs.shape[1]} spot-factor columns vs "
        f"{len(ctrl_labels)} signature labels."
    )

print("AD+CAA label order:")
print(ad_labels)

print("Control label order:")
print(ctrl_labels)

## Build condition-specific metadata tables

ad_obs = adata_wta.obs[
    adata_wta.obs["disease_status"] == "AD-CAA"
].copy()

ctrl_obs = adata_wta.obs[
    adata_wta.obs["disease_status"] == "Control"
].copy()

# Sanity checks: ROI counts must match model rows
if ad_spot_abs.shape[0] != ad_obs.shape[0]:
    raise ValueError(
        f"AD+CAA ROI count mismatch: "
        f"{ad_spot_abs.shape[0]} spot-factor rows vs "
        f"{ad_obs.shape[0]} metadata rows."
    )

if ctrl_spot_abs.shape[0] != ctrl_obs.shape[0]:
    raise ValueError(
        f"Control ROI count mismatch: "
        f"{ctrl_spot_abs.shape[0]} spot-factor rows vs "
        f"{ctrl_obs.shape[0]} metadata rows."
    )

## Build abundance data frames

ad_abs_df = make_abundance_dataframe(
    spot_factors=ad_spot_abs,
    obs=ad_obs,
    celltype_labels=ad_labels
)

ad_rel_df = make_abundance_dataframe(
    spot_factors=ad_spot_rel,
    obs=ad_obs,
    celltype_labels=ad_labels
)

ctrl_abs_df = make_abundance_dataframe(
    spot_factors=ctrl_spot_abs,
    obs=ctrl_obs,
    celltype_labels=ctrl_labels
)

ctrl_rel_df = make_abundance_dataframe(
    spot_factors=ctrl_spot_rel,
    obs=ctrl_obs,
    celltype_labels=ctrl_labels
)

## Create long-format relative-abundance tables

ad_long_rel = make_long_abundance_table(
    abundance_df=ad_rel_df,
    celltype_labels=ad_labels,
    abundance_col="rel_abundance"
)

ctrl_long_rel = make_long_abundance_table(
    abundance_df=ctrl_rel_df,
    celltype_labels=ctrl_labels,
    abundance_col="rel_abundance"
)

abundance_long_rel = pd.concat(
    [ad_long_rel, ctrl_long_rel],
    axis=0,
    ignore_index=True
)

## Create long-format absolute-abundance tables

ad_long_abs = make_long_abundance_table(
    abundance_df=ad_abs_df,
    celltype_labels=ad_labels,
    abundance_col="abs_abundance"
)

ctrl_long_abs = make_long_abundance_table(
    abundance_df=ctrl_abs_df,
    celltype_labels=ctrl_labels,
    abundance_col="abs_abundance"
)

abundance_long_abs = pd.concat(
    [ad_long_abs, ctrl_long_abs],
    axis=0,
    ignore_index=True
)

print("\nRelative abundance preview:")
print(abundance_long_rel.head())

print("\nAbsolute abundance preview:")
print(abundance_long_abs.head())

## Save outputs

save_abundance_outputs(
    output_dir="results/cell_proportions",
    ad_abs_df=ad_abs_df,
    ad_rel_df=ad_rel_df,
    ctrl_abs_df=ctrl_abs_df,
    ctrl_rel_df=ctrl_rel_df,
    long_df=abundance_long_rel
)

# Explicitly save the relative-abundance long table used downstream
abundance_long_rel.to_csv(
    "results/cell_proportions/roi_celltype_abundance_long.csv",
    index=False
)

# Save the absolute-abundance long table separately
abundance_long_abs.to_csv(
    "results/cell_proportions/roi_celltype_abundance_long_abs.csv",
    index=False
)

print("\nSaved:")
print("  results/cell_proportions/roi_celltype_abundance_long.csv")
print("  results/cell_proportions/roi_celltype_abundance_long_abs.csv")
