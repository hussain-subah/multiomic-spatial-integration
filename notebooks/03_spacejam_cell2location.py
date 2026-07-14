# ============================================================
# 03 — Bayesian spatial deconvolution with SpaceJam / Pyro
#
# Fits separate AD+CAA and Control spatial models using the
# condition-specific regression signatures produced by Notebook 02.
#
# Inputs:
#   data/CAA-AD_AnnData.h5ad
#   results/regression_model/AD+CAA_inferred_signatures.csv
#   results/regression_model/Control_inferred_signatures.csv
#
# Outputs:
#   results/spacejam/
#     ADCAA_param_store.pt
#     CTRL_param_store.pt
#     ADCAA_spot_factors_abs.pt
#     ADCAA_spot_factors_rel.pt
#     CTRL_spot_factors_abs.pt
#     CTRL_spot_factors_rel.pt
#     ADCAA_training_loss.csv
#     CTRL_training_loss.csv
#     *_manifest_rois.csv
#     *_manifest_factors.csv
#     *_manifest_experiments.csv
#     *_run_metadata.json
# ============================================================

from pathlib import Path
import gc
import json
import sys

import numpy as np
import pandas as pd
import pyro
import scanpy as sc
import torch

from models.LocationModelWTAMultiExperimentHierarchicalGeneLevel_Modified import (
    LocationModelPyro,
)


# ============================================================
# Configuration
# ============================================================

PROJECT_DIR = Path(
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/multiomic-spatial-integration"
)

DATA_DIR = PROJECT_DIR / "data"
RESULTS_ROOT = PROJECT_DIR / "results"

H5AD_FILE = DATA_DIR / "CAA-AD_AnnData.h5ad"

SIGNATURE_DIR = RESULTS_ROOT / "regression_model"
RESULTS_DIR = RESULTS_ROOT / "spacejam"

AD_SIGNATURE_FILE = (
    SIGNATURE_DIR / "AD+CAA_inferred_signatures.csv"
)

CTRL_SIGNATURE_FILE = (
    SIGNATURE_DIR / "Control_inferred_signatures.csv"
)

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

N_STEPS = 9000
LEARNING_RATE = 0.002

DEVICE = (
    "cuda"
    if torch.cuda.is_available()
    else "cpu"
)

CELL_NUMBER_PRIOR = {
    "cells_per_spot": 8.0,
    "combs_per_spot": 2.5,
    "factors_per_combs": 3.0,
    "cells_mean_var_ratio": 1.0,
    "factors_mean_var_ratio": 1.0,
    "combs_mean_var_ratio": 1.0,
}

print("Python:", sys.executable)
print("PyTorch:", torch.__version__)
print("Pyro:", pyro.__version__)
print("CUDA available:", torch.cuda.is_available())
print("Selected device:", DEVICE)
print("Project directory:", PROJECT_DIR)
print("Signature directory:", SIGNATURE_DIR)
print("SpaceJam output directory:", RESULTS_DIR)

if DEVICE == "cuda":
    print("GPU:", torch.cuda.get_device_name(0))
else:
    raise RuntimeError(
        "Notebook 03 requires a GPU allocation, "
        "but CUDA is unavailable."
    )


# ============================================================
# General helpers
# ============================================================

def dense_numpy(x, dtype=np.float32):
    """
    Convert a dense or sparse matrix-like object to NumPy.
    """
    if hasattr(x, "toarray"):
        x = x.toarray()

    return np.asarray(
        x,
        dtype=dtype,
    )


def make_spot2sample_mat(sample_ids):
    """
    Create an ROI × experiment one-hot matrix.

    Parameters
    ----------
    sample_ids
        One Scan_ID or experiment identifier per ROI.

    Returns
    -------
    matrix
        ROI × experiment one-hot matrix.
    categories
        Experiment labels in matrix-column order.
    """
    sample_cat = pd.Categorical(sample_ids)

    if np.any(sample_cat.codes < 0):
        raise ValueError(
            "Missing Scan_ID values were found while "
            "constructing spot2sample_mat."
        )

    matrix = np.zeros(
        (
            len(sample_cat),
            len(sample_cat.categories),
        ),
        dtype=np.float32,
    )

    matrix[
        np.arange(len(sample_cat)),
        sample_cat.codes,
    ] = 1.0

    categories = (
        sample_cat.categories
        .astype(str)
        .tolist()
    )

    return matrix, categories


# ============================================================
# Gene and matrix alignment
# ============================================================

def align_signature_and_spatial(
    signature_df,
    adata,
    counts_layer="counts",
):
    """
    Align a genes × factors signature matrix with spatial counts.

    The signature gene order is preserved and used to subset the
    spatial AnnData object explicitly.
    """
    if counts_layer not in adata.layers:
        raise KeyError(
            f"AnnData layer '{counts_layer}' was not found."
        )

    signature_df = signature_df.copy()

    signature_df.index = (
        signature_df.index
        .astype(str)
    )

    adata.var_names = (
        adata.var_names
        .astype(str)
    )

    if signature_df.index.duplicated().any():
        duplicated = (
            signature_df.index[
                signature_df.index.duplicated()
            ]
            .unique()
            .tolist()
        )

        raise ValueError(
            "The signature matrix contains duplicated genes: "
            f"{duplicated[:10]}"
        )

    if adata.var_names.duplicated().any():
        duplicated = (
            adata.var_names[
                adata.var_names.duplicated()
            ]
            .unique()
            .tolist()
        )

        raise ValueError(
            "The spatial AnnData contains duplicated genes: "
            f"{duplicated[:10]}"
        )

    shared_genes = signature_df.index[
        signature_df.index.isin(
            adata.var_names
        )
    ]

    if len(shared_genes) == 0:
        raise ValueError(
            "No shared genes were found between the "
            "signature and spatial datasets."
        )

    signature_aligned = (
        signature_df
        .loc[shared_genes]
        .copy()
    )

    adata_aligned = (
        adata[:, shared_genes]
        .copy()
    )

    signature_gene_order = (
        signature_aligned.index
        .astype(str)
        .to_numpy()
    )

    spatial_gene_order = (
        adata_aligned.var_names
        .astype(str)
        .to_numpy()
    )

    if not np.array_equal(
        signature_gene_order,
        spatial_gene_order,
    ):
        raise RuntimeError(
            "Gene alignment failed: signature and spatial "
            "gene orders differ."
        )

    counts = dense_numpy(
        adata_aligned.layers[counts_layer],
        dtype=np.float32,
    )

    signatures = (
        signature_aligned
        .to_numpy(dtype=np.float32)
    )

    if counts.shape[1] != signatures.shape[0]:
        raise RuntimeError(
            "Counts/signature dimension mismatch after alignment: "
            f"counts={counts.shape}, "
            f"signatures={signatures.shape}."
        )

    if not np.isfinite(counts).all():
        raise ValueError(
            "Spatial counts contain non-finite values."
        )

    if np.any(counts < 0):
        raise ValueError(
            "Spatial counts contain negative values."
        )

    if not np.isfinite(signatures).all():
        raise ValueError(
            "Signature matrix contains non-finite values."
        )

    if np.any(signatures < 0):
        raise ValueError(
            "Signature matrix contains negative values."
        )

    print(
        f"Aligned {len(shared_genes):,} genes across "
        f"{adata_aligned.n_obs} ROIs and "
        f"{signature_aligned.shape[1]} factors."
    )

    return (
        signature_aligned,
        adata_aligned,
        counts,
    )


def get_negative_probe_matrix(adata):
    """
    Retrieve the GeoMx negative-probe count matrix.
    """
    if "negProbes" not in adata.obsm:
        raise KeyError(
            "adata.obsm['negProbes'] was not found. "
            "Negative controls are required by this model."
        )

    negative_probes = dense_numpy(
        adata.obsm["negProbes"],
        dtype=np.float32,
    )

    if negative_probes.ndim != 2:
        raise ValueError(
            "Negative-probe input must be a two-dimensional matrix."
        )

    if negative_probes.shape[0] != adata.n_obs:
        raise ValueError(
            "Negative-probe matrix row count does not "
            "match the number of spatial ROIs."
        )

    if not np.isfinite(negative_probes).all():
        raise ValueError(
            "Negative-probe matrix contains non-finite values."
        )

    if np.any(negative_probes < 0):
        raise ValueError(
            "Negative-probe matrix contains negative values."
        )

    return negative_probes


# ============================================================
# Saving and posterior extraction
# ============================================================

def save_param_store(output_file):
    """
    Save the complete Pyro parameter-store state on CPU.
    """
    state = (
        pyro.get_param_store()
        .get_state()
    )

    cpu_state = {
        "params": {
            name: tensor.detach().cpu()
            for name, tensor
            in state["params"].items()
        },
        "constraints": state["constraints"],
    }

    torch.save(
        cpu_state,
        output_file,
    )


def extract_constrained_spot_factors(model):
    """
    Extract constrained posterior medians from AutoNormal.

    AutoNormal location parameters are stored in unconstrained
    latent space. model.guide.median() applies each latent site's
    support transform and therefore returns positive Gamma-supported
    spot-factor values.
    """
    if model.guide is None:
        raise RuntimeError(
            "The model guide is unavailable. "
            "Fit the model before posterior extraction."
        )

    with torch.no_grad():
        posterior_median = (
            model.guide.median()
        )

    if "spot_factors" not in posterior_median:
        raise KeyError(
            "The posterior guide did not return a "
            "'spot_factors' site. Available sites: "
            f"{list(posterior_median)}"
        )

    spot_abs = (
        posterior_median["spot_factors"]
        .detach()
        .cpu()
        .to(torch.float32)
    )

    if spot_abs.ndim != 2:
        raise ValueError(
            "spot_factors must be two-dimensional; "
            f"received shape {tuple(spot_abs.shape)}."
        )

    if not torch.isfinite(spot_abs).all():
        raise ValueError(
            "Non-finite posterior spot-factor values were found."
        )

    if torch.any(spot_abs <= 0):
        raise ValueError(
            "Posterior spot factors contain non-positive values."
        )

    row_totals = (
        spot_abs
        .sum(
            dim=1,
            keepdim=True,
        )
    )

    if torch.any(row_totals <= 0):
        raise ValueError(
            "One or more spot-factor rows have "
            "a non-positive total."
        )

    spot_rel = (
        spot_abs /
        row_totals
    )

    return spot_abs, spot_rel


def validate_factor_outputs(
    spot_abs,
    spot_rel,
    n_rois,
    factor_labels,
    condition_name,
):
    """
    Validate output dimensions, positivity, and relative row sums.
    """
    expected_shape = (
        n_rois,
        len(factor_labels),
    )

    if tuple(spot_abs.shape) != expected_shape:
        raise ValueError(
            f"{condition_name}: absolute factor shape "
            f"{tuple(spot_abs.shape)} does not match "
            f"expected shape {expected_shape}."
        )

    if tuple(spot_rel.shape) != expected_shape:
        raise ValueError(
            f"{condition_name}: relative factor shape "
            f"{tuple(spot_rel.shape)} does not match "
            f"expected shape {expected_shape}."
        )

    if torch.any(spot_rel < 0):
        raise ValueError(
            f"{condition_name}: negative relative factors found."
        )

    if torch.any(spot_rel > 1):
        raise ValueError(
            f"{condition_name}: relative factors greater than 1 found."
        )

    row_sums = spot_rel.sum(dim=1)

    if not torch.allclose(
        row_sums,
        torch.ones_like(row_sums),
        atol=1e-5,
        rtol=1e-5,
    ):
        raise ValueError(
            f"{condition_name}: relative factors do not sum to one."
        )

    print(f"\n{condition_name} output validation")
    print("  Shape:", tuple(spot_abs.shape))

    print(
        "  Absolute range:",
        float(spot_abs.min()),
        "to",
        float(spot_abs.max()),
    )

    print(
        "  Relative range:",
        float(spot_rel.min()),
        "to",
        float(spot_rel.max()),
    )

    print(
        "  Relative row-sum range:",
        float(row_sums.min()),
        "to",
        float(row_sums.max()),
    )


# ============================================================
# Manifests
# ============================================================

def save_factor_manifests(
    file_prefix,
    adata,
    factor_labels,
    experiment_names,
):
    """
    Save exact ROI, factor, and experiment ordering.
    """
    roi_manifest = pd.DataFrame(
        {
            "row_index": np.arange(
                adata.n_obs
            ),
            "ROI_ID": (
                adata.obs_names
                .astype(str)
            ),
        }
    )

    for column in [
        "Scan_ID",
        "disease_status",
        "pathology",
        "region",
    ]:
        if column in adata.obs.columns:
            roi_manifest[column] = (
                adata.obs[column]
                .astype(str)
                .to_numpy()
            )

    roi_manifest.to_csv(
        RESULTS_DIR /
        f"{file_prefix}_manifest_rois.csv",
        index=False,
    )

    factor_manifest = pd.DataFrame(
        {
            "factor_index": np.arange(
                len(factor_labels)
            ),
            "celltype": factor_labels,
        }
    )

    factor_manifest.to_csv(
        RESULTS_DIR /
        f"{file_prefix}_manifest_factors.csv",
        index=False,
    )

    experiment_manifest = pd.DataFrame(
        {
            "experiment_index": np.arange(
                len(experiment_names)
            ),
            "Scan_ID": experiment_names,
        }
    )

    experiment_manifest.to_csv(
        RESULTS_DIR /
        f"{file_prefix}_manifest_experiments.csv",
        index=False,
    )


# ============================================================
# Condition-specific training
# ============================================================

def train_condition_model(
    condition_name,
    file_prefix,
    signature_df,
    adata_condition,
):
    """
    Train one condition-specific SpaceJam model and save outputs.
    """
    print("\n" + "=" * 72)
    print(f"Training SpaceJam model: {condition_name}")
    print("=" * 72)

    (
        signature_aligned,
        adata_aligned,
        counts,
    ) = align_signature_and_spatial(
        signature_df=signature_df,
        adata=adata_condition,
        counts_layer="counts",
    )

    negative_probes = get_negative_probe_matrix(
        adata_aligned
    )

    (
        spot2sample,
        experiment_names,
    ) = make_spot2sample_mat(
        adata_aligned.obs["Scan_ID"]
    )

    if negative_probes.shape[0] != counts.shape[0]:
        raise ValueError(
            f"{condition_name}: counts and negative probes "
            "have different ROI counts."
        )

    if spot2sample.shape[0] != counts.shape[0]:
        raise ValueError(
            f"{condition_name}: counts and spot2sample "
            "have different ROI counts."
        )

    factor_labels = (
        signature_aligned.columns
        .astype(str)
        .tolist()
    )

    print("ROIs:", counts.shape[0])
    print("Genes:", counts.shape[1])
    print(
        "Negative probes:",
        negative_probes.shape[1],
    )
    print(
        "Experiments:",
        spot2sample.shape[1],
    )
    print("Factors:", len(factor_labels))
    print(
        "First five factors:",
        factor_labels[:5],
    )

    pyro.clear_param_store()

    if DEVICE == "cuda":
        torch.cuda.empty_cache()

    model = LocationModelPyro(
        cell_state_mat=(
            signature_aligned
            .to_numpy(dtype=np.float32)
        ),
        X_data=counts,
        Y_data=negative_probes,
        spot2sample_mat=spot2sample,
        device=DEVICE,
        cell_number_prior=CELL_NUMBER_PRIOR,
    )

    loss_history = model.fit(
        n_steps=N_STEPS,
        lr=LEARNING_RATE,
    )

    loss_df = pd.DataFrame(
        {
            "step": np.arange(
                len(loss_history)
            ),
            "loss": loss_history,
        }
    )

    loss_df.to_csv(
        RESULTS_DIR /
        f"{file_prefix}_training_loss.csv",
        index=False,
    )

    save_param_store(
        RESULTS_DIR /
        f"{file_prefix}_param_store.pt"
    )

    spot_abs, spot_rel = (
        extract_constrained_spot_factors(
            model
        )
    )

    validate_factor_outputs(
        spot_abs=spot_abs,
        spot_rel=spot_rel,
        n_rois=adata_aligned.n_obs,
        factor_labels=factor_labels,
        condition_name=condition_name,
    )

    torch.save(
        spot_abs,
        RESULTS_DIR /
        f"{file_prefix}_spot_factors_abs.pt",
    )

    torch.save(
        spot_rel,
        RESULTS_DIR /
        f"{file_prefix}_spot_factors_rel.pt",
    )

    save_factor_manifests(
        file_prefix=file_prefix,
        adata=adata_aligned,
        factor_labels=factor_labels,
        experiment_names=experiment_names,
    )

    metadata = {
        "condition": condition_name,
        "n_rois": int(
            adata_aligned.n_obs
        ),
        "n_genes": int(
            adata_aligned.n_vars
        ),
        "n_negative_probes": int(
            negative_probes.shape[1]
        ),
        "n_experiments": int(
            spot2sample.shape[1]
        ),
        "n_factors": int(
            len(factor_labels)
        ),
        "factor_labels": factor_labels,
        "experiment_names": experiment_names,
        "n_steps": N_STEPS,
        "learning_rate": LEARNING_RATE,
        "posterior_summary": (
            "AutoNormal constrained posterior median"
        ),
        "signature_file": str(
            AD_SIGNATURE_FILE
            if condition_name == "AD+CAA"
            else CTRL_SIGNATURE_FILE
        ),
    }

    with open(
        RESULTS_DIR /
        f"{file_prefix}_run_metadata.json",
        "w",
    ) as handle:
        json.dump(
            metadata,
            handle,
            indent=2,
        )

    print(
        f"{condition_name}: outputs saved under "
        f"{RESULTS_DIR}"
    )

    result = {
        "spot_abs": spot_abs,
        "spot_rel": spot_rel,
        "factor_labels": factor_labels,
        "roi_ids": (
            adata_aligned.obs_names
            .astype(str)
            .tolist()
        ),
    }

    del model
    gc.collect()

    if DEVICE == "cuda":
        torch.cuda.empty_cache()

    pyro.clear_param_store()

    return result


# ============================================================
# Validate required inputs
# ============================================================

required_files = [
    H5AD_FILE,
    AD_SIGNATURE_FILE,
    CTRL_SIGNATURE_FILE,
]

missing_files = [
    path
    for path in required_files
    if not path.exists()
]

if missing_files:
    raise FileNotFoundError(
        "Required input file(s) were not found:\n" +
        "\n".join(
            str(path)
            for path in missing_files
        )
    )


# ============================================================
# Load GeoMx spatial data
# ============================================================

adata_wta = sc.read_h5ad(
    H5AD_FILE
)

if "counts" not in adata_wta.layers:
    raise KeyError(
        "adata_wta.layers['counts'] is required. "
        "Do not silently substitute normalized adata.X."
    )

required_obs_columns = {
    "disease_status",
    "Scan_ID",
}

missing_obs_columns = (
    required_obs_columns -
    set(adata_wta.obs.columns)
)

if missing_obs_columns:
    raise KeyError(
        "Missing required AnnData metadata columns: "
        f"{sorted(missing_obs_columns)}"
    )

print("\nLoaded spatial AnnData:")
print(adata_wta)

print("\nDisease-status counts:")
print(
    adata_wta.obs["disease_status"]
    .value_counts(dropna=False)
)


# ============================================================
# Load corrected regression signatures
# ============================================================

X_ref_ADCAA = pd.read_csv(
    AD_SIGNATURE_FILE,
    index_col=0,
)

X_ref_Control = pd.read_csv(
    CTRL_SIGNATURE_FILE,
    index_col=0,
)

X_ref_ADCAA.index = (
    X_ref_ADCAA.index
    .astype(str)
)

X_ref_Control.index = (
    X_ref_Control.index
    .astype(str)
)

ad_labels = (
    X_ref_ADCAA.columns
    .astype(str)
    .tolist()
)

ctrl_labels = (
    X_ref_Control.columns
    .astype(str)
    .tolist()
)

print("\nAD+CAA factor count:", len(ad_labels))
print("Control factor count:", len(ctrl_labels))
print(
    "Same factor-label set:",
    set(ad_labels) == set(ctrl_labels),
)
print(
    "Same factor-label order:",
    ad_labels == ctrl_labels,
)

if ad_labels != ctrl_labels:
    raise ValueError(
        "AD+CAA and Control signature columns do not have "
        "the same authoritative factor order. Notebook 02 "
        "validation must pass before running SpaceJam.\n\n"
        f"AD order:\n{ad_labels}\n\n"
        f"Control order:\n{ctrl_labels}"
    )

if X_ref_ADCAA.columns.duplicated().any():
    raise ValueError(
        "AD+CAA signatures contain duplicated factor names."
    )

if X_ref_Control.columns.duplicated().any():
    raise ValueError(
        "Control signatures contain duplicated factor names."
    )


# ============================================================
# Split spatial ROIs by disease status
# ============================================================

adata_wta_ADCAA = adata_wta[
    adata_wta.obs["disease_status"] == "AD-CAA"
].copy()

adata_wta_Control = adata_wta[
    adata_wta.obs["disease_status"] == "Control"
].copy()

if adata_wta_ADCAA.n_obs == 0:
    raise ValueError(
        "No AD-CAA ROIs were found."
    )

if adata_wta_Control.n_obs == 0:
    raise ValueError(
        "No Control ROIs were found."
    )

print(
    "\nAD+CAA ROI count:",
    adata_wta_ADCAA.n_obs,
)

print(
    "Control ROI count:",
    adata_wta_Control.n_obs,
)


# ============================================================
# Train the AD+CAA model
# ============================================================

ad_results = train_condition_model(
    condition_name="AD+CAA",
    file_prefix="ADCAA",
    signature_df=X_ref_ADCAA,
    adata_condition=adata_wta_ADCAA,
)


# ============================================================
# Train the Control model
# ============================================================

ctrl_results = train_condition_model(
    condition_name="Control",
    file_prefix="CTRL",
    signature_df=X_ref_Control,
    adata_condition=adata_wta_Control,
)


# ============================================================
# Final cross-condition validation
# ============================================================

if (
    ad_results["factor_labels"] !=
    ctrl_results["factor_labels"]
):
    raise RuntimeError(
        "Saved AD+CAA and Control factors have "
        "different label orders."
    )

if (
    ad_results["spot_abs"].shape[1] !=
    ctrl_results["spot_abs"].shape[1]
):
    raise RuntimeError(
        "Saved AD+CAA and Control outputs have "
        "different factor counts."
    )

print("\n" + "=" * 72)
print(
    "SpaceJam training and validation completed successfully"
)
print("=" * 72)

print(
    "AD+CAA shape:",
    tuple(ad_results["spot_abs"].shape),
)

print(
    "Control shape:",
    tuple(ctrl_results["spot_abs"].shape),
)

print(
    "Factor count:",
    len(ad_results["factor_labels"]),
)

print(
    "Output directory:",
    RESULTS_DIR,
)
