# =============================================================================
# 06 — Pseudobulk Validation of SpaceJam / Cell2Location Deconvolution
#
# This script validates deconvolution recovery using synthetic pseudobulk ROIs
# generated from held-out snRNA-seq cells.
#
# Important design:
#   - Regression signatures are trained using train cells only.
#   - Synthetic ROIs are generated using held-out cells only.
#   - Cell types with insufficient held-out cells are excluded from the
#     recovery experiment, but may remain in the regression training model.
#   - Signature columns always follow the exact factor order stored in the
#     regression model registry.
#   - AD+CAA and Control are validated independently.
#
# Outputs:
#   results/pseudobulk_validation/
#       AD+CAA/
#           signature/
#           recovery_metrics.csv
#           ground_truth_proportions.csv
#           inferred_proportions.csv
#           training_loss.csv
#           recovery_scatter.pdf
#           recovery_bias.pdf
#           validation_factor_manifest.csv
#           excluded_celltypes.csv
#           run_metadata.json
#       Control/
#           ...
#       combined_recovery_metrics.csv
# =============================================================================

from pathlib import Path
import gc
import json
import sys

import numpy as np
import pandas as pd
import pyro
import scanpy as sc
import torch

from python.regression_utils import (
    load_seurat_mtx_as_anndata,
    standardize_regression_obs,
    round_counts_layer,
    minimal_filter_anndata,
    disable_lightning_mpi_detection,
    train_regression_model,
)

from python.pseudobulk_utils import (
    split_reference_cells,
    generate_synthetic_dataset,
)

from python.spacejam_pyro_utils import (
    align_genes,
)

from python.pseudobulk_validation_utils import (
    evaluate_recovery,
    plot_recovery_scatter,
    plot_recovery_bias_by_abundance,
)

from models.LocationModelWTAMultiExperimentHierarchicalGeneLevel_Modified import (
    LocationModelPyro,
)


# =============================================================================
# Configuration
# =============================================================================

PROJECT_DIR = Path(
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/multiomic-spatial-integration"
)

DATA_ROOT = Path(
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025"
)

OUTPUT_ROOT = (
    PROJECT_DIR
    / "results"
    / "pseudobulk_validation"
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

REFERENCE_COUNTS_FILE = DATA_ROOT / "AD_CAA_counts.mtx"
REFERENCE_METADATA_FILE = DATA_ROOT / "AD_CAA_meta.csv"
REFERENCE_GENES_FILE = DATA_ROOT / "AD_CAA_genes.csv"

SPATIAL_H5AD_FILE = (
    PROJECT_DIR
    / "data"
    / "CAA-AD_AnnData.h5ad"
)

CONDITION_MAP = {
    "AD+CAA": "AD-CAA",
    "Control": "Control",
}

CELL_NUMBER_PRIOR = {
    "cells_per_spot": 8.0,
    "factors_per_spot": 7.0,
    "combs_per_spot": 2.5,
    "factors_per_combs": 3.0,
    "cells_mean_var_ratio": 1.0,
    "factors_mean_var_ratio": 1.0,
    "combs_mean_var_ratio": 1.0,
}

# Synthetic validation settings
N_MIXTURES = 100
N_CELLS_PER_MIXTURE = 200
DIRICHLET_CONCENTRATION = 0.3

# Reference regression model settings
REGRESSION_MAX_EPOCHS = 100
REGRESSION_BATCH_SIZE = 1024
REGRESSION_LEARNING_RATE = 0.01
REGRESSION_POSTERIOR_SAMPLES = 1000
REGRESSION_POSTERIOR_BATCH_SIZE = 2500

# SpaceJam settings
SPACEJAM_STEPS = 9000
SPACEJAM_LEARNING_RATE = 0.002
POSTERIOR_PREDICTIVE_SAMPLES = 100

# Train/holdout settings
HOLDOUT_FRACTION = 0.3
MIN_HOLDOUT_CELLS = 10
RANDOM_SEED = 0

DEVICE = (
    "cuda"
    if torch.cuda.is_available()
    else "cpu"
)

print("Python:", sys.executable)
print("PyTorch:", torch.__version__)
print("Pyro:", pyro.__version__)
print("CUDA available:", torch.cuda.is_available())
print("Selected device:", DEVICE)
print("Output directory:", OUTPUT_ROOT)

if DEVICE != "cuda":
    raise RuntimeError(
        "Notebook 06 requires a GPU allocation, "
        "but CUDA is unavailable."
    )

print("GPU:", torch.cuda.get_device_name(0))


# =============================================================================
# Helper functions
# =============================================================================

def require_files(paths):
    """
    Confirm that all required input files exist.
    """
    missing = [
        path
        for path in paths
        if not path.exists()
    ]

    if missing:
        raise FileNotFoundError(
            "Required input file(s) were not found:\n"
            + "\n".join(str(path) for path in missing)
        )


def get_model_celltype_order(model):
    """
    Recover the exact categorical factor order used by the trained
    regression model.

    The model registry order is authoritative for the factor axis of
    per_cluster_mu_fg.
    """
    state = None
    successful_key = None

    for key in ("labels", "celltype"):
        try:
            state = model.adata_manager.get_state_registry(key)
            successful_key = key
            break
        except Exception:
            continue

    if state is None:
        raise ValueError(
            "Could not recover the regression model label registry. "
            "Tried registry keys 'labels' and 'celltype'."
        )

    if hasattr(state, "categorical_mapping"):
        labels = list(state.categorical_mapping)

    elif isinstance(state, dict) and "categorical_mapping" in state:
        labels = list(state["categorical_mapping"])

    else:
        try:
            labels = list(state["categorical_mapping"])
        except Exception as exc:
            raise ValueError(
                "A label registry was found, but categorical_mapping "
                f"could not be recovered. Registry key: {successful_key}"
            ) from exc

    labels = [
        str(label)
        for label in labels
    ]

    if len(labels) == 0:
        raise ValueError(
            "The recovered regression factor order is empty."
        )

    if len(labels) != len(set(labels)):
        raise ValueError(
            "The recovered regression factor order contains duplicates."
        )

    return labels


def export_validation_signatures(
    adata,
    model,
    model_name,
    output_folder,
):
    """
    Export the full train-only signature matrix using the exact regression
    model registry order.

    Returns
    -------
    pandas.DataFrame
        Genes × all trained factors.
    """
    output_folder = Path(output_folder)

    output_folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    try:
        posterior_signature = (
            adata.uns["mod"]
            ["post_sample_means"]
            ["per_cluster_mu_fg"]
            .T
        )
    except KeyError as exc:
        raise KeyError(
            f"{model_name}: could not find "
            "adata.uns['mod']['post_sample_means']"
            "['per_cluster_mu_fg']."
        ) from exc

    posterior_signature = np.asarray(
        posterior_signature
    )

    model_order = get_model_celltype_order(
        model
    )

    expected_shape = (
        adata.n_vars,
        len(model_order),
    )

    if posterior_signature.ndim != 2:
        raise ValueError(
            f"{model_name}: expected a two-dimensional signature "
            f"matrix, received {posterior_signature.shape}."
        )

    if posterior_signature.shape != expected_shape:
        raise ValueError(
            f"{model_name}: posterior signature shape "
            f"{posterior_signature.shape} does not match expected "
            f"genes × factors shape {expected_shape}."
        )

    if not np.isfinite(posterior_signature).all():
        raise ValueError(
            f"{model_name}: signature matrix contains non-finite values."
        )

    if np.any(posterior_signature < 0):
        raise ValueError(
            f"{model_name}: signature matrix contains negative values."
        )

    signature_df = pd.DataFrame(
        posterior_signature,
        index=adata.var_names.astype(str),
        columns=model_order,
    )

    signature_df.index.name = "gene"

    output_file = (
        output_folder
        / f"{model_name}_all_train_factors.csv"
    )

    signature_df.to_csv(
        output_file
    )

    print(
        f"Saved full validation signature matrix: {output_file}"
    )
    print(
        "Full signature shape:",
        signature_df.shape,
    )
    print(
        "First five trained factors:",
        model_order[:5],
    )

    return signature_df


def determine_evaluable_celltypes(
    adata_train,
    adata_holdout,
    model_order,
    condition_name,
):
    """
    Determine which model factors have enough held-out cells for recovery
    validation.

    Factors remain ordered according to model_order.
    """
    train_counts = (
        adata_train.obs["celltype"]
        .astype(str)
        .value_counts()
    )

    holdout_counts = (
        adata_holdout.obs["celltype"]
        .astype(str)
        .value_counts()
    )

    eligible = [
        celltype
        for celltype in model_order
        if (
            celltype in holdout_counts.index
            and holdout_counts[celltype] >= MIN_HOLDOUT_CELLS
        )
    ]

    excluded_rows = []

    for celltype in model_order:

        n_train = int(
            train_counts.get(
                celltype,
                0,
            )
        )

        n_holdout = int(
            holdout_counts.get(
                celltype,
                0,
            )
        )

        if celltype not in eligible:

            if n_holdout == 0:
                reason = "no_holdout_cells"
            elif n_holdout < MIN_HOLDOUT_CELLS:
                reason = "insufficient_holdout_cells"
            else:
                reason = "not_evaluable"

            excluded_rows.append(
                {
                    "celltype": celltype,
                    "n_train_cells": n_train,
                    "n_holdout_cells": n_holdout,
                    "reason": reason,
                }
            )

    if len(eligible) == 0:
        raise RuntimeError(
            f"{condition_name}: no factors have sufficient held-out cells."
        )

    print(
        f"{condition_name}: "
        f"{len(eligible)} of {len(model_order)} trained factors "
        "are evaluable."
    )

    print(
        f"{condition_name}: evaluable factors:"
    )
    print(eligible)

    print(
        f"{condition_name}: excluded factors:"
    )
    print(
        [
            row["celltype"]
            for row in excluded_rows
        ]
    )

    excluded_df = pd.DataFrame(
        excluded_rows,
        columns=[
            "celltype",
            "n_train_cells",
            "n_holdout_cells",
            "reason",
        ],
    )

    factor_manifest = pd.DataFrame(
        {
            "factor_index": np.arange(
                len(eligible)
            ),
            "celltype": eligible,
            "n_train_cells": [
                int(train_counts.get(celltype, 0))
                for celltype in eligible
            ],
            "n_holdout_cells": [
                int(holdout_counts.get(celltype, 0))
                for celltype in eligible
            ],
        }
    )

    return (
        eligible,
        excluded_df,
        factor_manifest,
    )


def validate_synthetic_dataset(
    synthetic,
    factor_labels,
    condition_name,
):
    """
    Validate generated synthetic matrices and ground-truth proportions.
    """
    required_keys = {
        "X_data",
        "Y_data",
        "spot2sample_mat",
        "ground_truth_df",
        "donor_roi_ids",
    }

    missing_keys = (
        required_keys
        - set(synthetic)
    )

    if missing_keys:
        raise KeyError(
            f"{condition_name}: synthetic dataset is missing keys: "
            f"{sorted(missing_keys)}"
        )

    x_data = np.asarray(
        synthetic["X_data"]
    )

    y_data = np.asarray(
        synthetic["Y_data"]
    )

    spot2sample = np.asarray(
        synthetic["spot2sample_mat"]
    )

    ground_truth = synthetic[
        "ground_truth_df"
    ]

    if x_data.ndim != 2:
        raise ValueError(
            f"{condition_name}: X_data must be two-dimensional."
        )

    if y_data.ndim != 2:
        raise ValueError(
            f"{condition_name}: Y_data must be two-dimensional."
        )

    if spot2sample.ndim != 2:
        raise ValueError(
            f"{condition_name}: spot2sample_mat must be two-dimensional."
        )

    if x_data.shape[0] != N_MIXTURES:
        raise ValueError(
            f"{condition_name}: expected {N_MIXTURES} synthetic ROIs, "
            f"found {x_data.shape[0]}."
        )

    if y_data.shape[0] != x_data.shape[0]:
        raise ValueError(
            f"{condition_name}: X_data and Y_data have different "
            "numbers of synthetic ROIs."
        )

    if spot2sample.shape[0] != x_data.shape[0]:
        raise ValueError(
            f"{condition_name}: X_data and spot2sample_mat have "
            "different numbers of synthetic ROIs."
        )

    expected_ground_truth_shape = (
        x_data.shape[0],
        len(factor_labels),
    )

    if ground_truth.shape != expected_ground_truth_shape:
        raise ValueError(
            f"{condition_name}: ground-truth shape "
            f"{ground_truth.shape} does not match expected "
            f"{expected_ground_truth_shape}."
        )

    if ground_truth.columns.tolist() != factor_labels:
        raise RuntimeError(
            f"{condition_name}: ground-truth factor order does not "
            "match the validation signature order."
        )

    if not np.isfinite(x_data).all():
        raise ValueError(
            f"{condition_name}: X_data contains non-finite values."
        )

    if not np.isfinite(y_data).all():
        raise ValueError(
            f"{condition_name}: Y_data contains non-finite values."
        )

    ground_truth_values = ground_truth.to_numpy(
        dtype=float
    )

    if np.any(ground_truth_values < 0):
        raise ValueError(
            f"{condition_name}: ground truth contains negative values."
        )

    ground_truth_sums = (
        ground_truth_values.sum(axis=1)
    )

    if not np.allclose(
        ground_truth_sums,
        1.0,
        atol=1e-6,
    ):
        raise ValueError(
            f"{condition_name}: ground-truth rows do not sum to one."
        )

    print(
        f"{condition_name}: synthetic dataset validation passed."
    )
    print("  X_data shape:", x_data.shape)
    print("  Y_data shape:", y_data.shape)
    print("  spot2sample shape:", spot2sample.shape)
    print("  ground truth shape:", ground_truth.shape)


def extract_posterior_proportions(
    model,
    synthetic_roi_ids,
    factor_labels,
    num_samples,
):
    """
    Extract posterior predictive spot factors and convert them to
    within-ROI relative proportions.
    """
    with torch.no_grad():
        posterior = model.posterior_predictive(
            num_samples=num_samples
        )

    if "spot_factors" not in posterior:
        raise KeyError(
            "posterior_predictive() did not return 'spot_factors'. "
            f"Available keys: {list(posterior)}"
        )

    spot_factors = posterior[
        "spot_factors"
    ]

    if hasattr(spot_factors, "detach"):
        spot_factors = (
            spot_factors
            .detach()
            .cpu()
            .numpy()
        )

    spot_factors = np.asarray(
        spot_factors
    )

    expected_3d_shape = (
        num_samples,
        len(synthetic_roi_ids),
        len(factor_labels),
    )

    expected_2d_shape = (
        len(synthetic_roi_ids),
        len(factor_labels),
    )

    if spot_factors.ndim == 3:

        if spot_factors.shape != expected_3d_shape:
            raise ValueError(
                "Unexpected posterior spot-factor shape: "
                f"{spot_factors.shape}. Expected {expected_3d_shape}."
            )

        abundance = spot_factors.mean(
            axis=0
        )

    elif spot_factors.ndim == 2:

        if spot_factors.shape != expected_2d_shape:
            raise ValueError(
                "Unexpected posterior spot-factor shape: "
                f"{spot_factors.shape}. Expected {expected_2d_shape}."
            )

        abundance = spot_factors

    else:
        raise ValueError(
            "Posterior spot factors must be two- or three-dimensional. "
            f"Received {spot_factors.shape}."
        )

    if not np.isfinite(abundance).all():
        raise ValueError(
            "Posterior abundance matrix contains non-finite values."
        )

    if np.any(abundance < 0):
        raise ValueError(
            "Posterior abundance matrix contains negative values."
        )

    row_totals = abundance.sum(
        axis=1,
        keepdims=True,
    )

    if np.any(row_totals <= 0):
        raise ValueError(
            "One or more posterior abundance rows have "
            "non-positive totals."
        )

    proportions = (
        abundance
        / row_totals
    )

    inferred_df = pd.DataFrame(
        proportions,
        index=synthetic_roi_ids,
        columns=factor_labels,
    )

    inferred_df.index.name = "roi_id"

    if inferred_df.columns.tolist() != factor_labels:
        raise RuntimeError(
            "Inferred factor order does not match the validation "
            "signature order."
        )

    if not np.allclose(
        inferred_df.sum(axis=1),
        1.0,
        atol=1e-6,
    ):
        raise ValueError(
            "Inferred proportion rows do not sum to one."
        )

    return inferred_df


# =============================================================================
# Validate required files
# =============================================================================

require_files(
    [
        REFERENCE_COUNTS_FILE,
        REFERENCE_METADATA_FILE,
        REFERENCE_GENES_FILE,
        SPATIAL_H5AD_FILE,
    ]
)


# =============================================================================
# Load snRNA-seq reference
# =============================================================================

adata_ref = load_seurat_mtx_as_anndata(
    counts_mtx=str(
        REFERENCE_COUNTS_FILE
    ),
    metadata_csv=str(
        REFERENCE_METADATA_FILE
    ),
    genes_csv=str(
        REFERENCE_GENES_FILE
    ),
    gene_col="gene_id",
    transpose_counts=True,
    counts_layer="counts",
)

adata_ref = standardize_regression_obs(
    adata_ref,
    celltype_col="New_Idents",
    new_celltype_col="celltype",
    experiment_col="Experiment",
    experiment_value="Batch1",
)

adata_ref = round_counts_layer(
    adata_ref,
    counts_layer="counts",
    replace_x=True,
)

adata_ref = minimal_filter_anndata(
    adata_ref,
    min_genes=1,
    min_cells=1,
    counts_layer="counts",
)

required_reference_obs = {
    "FDX",
    "celltype",
    "Experiment",
}

missing_reference_obs = (
    required_reference_obs
    - set(adata_ref.obs.columns)
)

if missing_reference_obs:
    raise KeyError(
        "Reference AnnData is missing required columns: "
        f"{sorted(missing_reference_obs)}"
    )

print("\nLoaded reference AnnData:")
print(adata_ref)

print("\nReference condition counts:")
print(
    adata_ref.obs["FDX"]
    .value_counts(dropna=False)
)

disable_lightning_mpi_detection()


# =============================================================================
# Load spatial GeoMx AnnData
# =============================================================================

adata_wta = sc.read_h5ad(
    SPATIAL_H5AD_FILE
)

if "counts" not in adata_wta.layers:
    raise KeyError(
        "Spatial AnnData does not contain layers['counts']. "
        "Do not substitute normalized adata.X."
    )

if "negProbes" not in adata_wta.obsm:
    raise KeyError(
        "Spatial AnnData does not contain obsm['negProbes']."
    )

if "disease_status" not in adata_wta.obs.columns:
    raise KeyError(
        "Spatial AnnData is missing disease_status."
    )

print("\nLoaded spatial AnnData:")
print(adata_wta)

print("\nSpatial disease-status counts:")
print(
    adata_wta.obs["disease_status"]
    .value_counts(dropna=False)
)


# =============================================================================
# Run validation separately by condition
# =============================================================================

results_by_condition = {}

for ref_condition, spatial_condition in CONDITION_MAP.items():

    print("\n" + "=" * 80)
    print(f"Pseudobulk validation: {ref_condition}")
    print("=" * 80)

    condition_output_dir = (
        OUTPUT_ROOT
        / ref_condition
    )

    signature_output_dir = (
        condition_output_dir
        / "signature"
    )

    condition_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    signature_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    pyro.clear_param_store()
    torch.cuda.empty_cache()

    # -------------------------------------------------------------------------
    # 1. Subset reference cells and create train/holdout split
    # -------------------------------------------------------------------------

    adata_ref_cond = adata_ref[
        adata_ref.obs["FDX"] == ref_condition
    ].copy()

    if adata_ref_cond.n_obs == 0:
        raise ValueError(
            f"No reference cells found for {ref_condition}."
        )

    adata_train, adata_holdout = split_reference_cells(
        adata_ref_cond,
        celltype_col="celltype",
        holdout_frac=HOLDOUT_FRACTION,
        min_holdout_cells=MIN_HOLDOUT_CELLS,
        seed=RANDOM_SEED,
    )

    print("Reference cells:", adata_ref_cond.n_obs)
    print("Train cells:", adata_train.n_obs)
    print("Holdout cells:", adata_holdout.n_obs)

    if adata_train.n_obs == 0:
        raise ValueError(
            f"{ref_condition}: training split is empty."
        )

    if adata_holdout.n_obs == 0:
        raise ValueError(
            f"{ref_condition}: holdout split is empty."
        )

    # -------------------------------------------------------------------------
    # 2. Train the reference regression model
    # -------------------------------------------------------------------------

    model_regression = train_regression_model(
        adata_train,
        labels_key="celltype",
        batch_key="Experiment",
        layer="counts",
        max_epochs=REGRESSION_MAX_EPOCHS,
        batch_size=REGRESSION_BATCH_SIZE,
        lr=REGRESSION_LEARNING_RATE,
        accelerator="gpu",
    )

    model_regression.export_posterior(
        adata_train,
        sample_kwargs={
            "num_samples": REGRESSION_POSTERIOR_SAMPLES,
            "batch_size": REGRESSION_POSTERIOR_BATCH_SIZE,
        },
    )

    full_signature_df = export_validation_signatures(
        adata=adata_train,
        model=model_regression,
        model_name=f"{ref_condition}_holdout_validation",
        output_folder=signature_output_dir,
    )

    model_order = get_model_celltype_order(
        model_regression
    )

    if full_signature_df.columns.tolist() != model_order:
        raise RuntimeError(
            f"{ref_condition}: exported full signature order differs "
            "from the model registry."
        )

    # -------------------------------------------------------------------------
    # 3. Select only cell types with enough held-out cells
    # -------------------------------------------------------------------------

    (
        evaluable_celltypes,
        excluded_celltypes_df,
        factor_manifest_df,
    ) = determine_evaluable_celltypes(
        adata_train=adata_train,
        adata_holdout=adata_holdout,
        model_order=model_order,
        condition_name=ref_condition,
    )

    signature_df = (
        full_signature_df
        .loc[:, evaluable_celltypes]
        .copy()
    )

    if signature_df.columns.tolist() != evaluable_celltypes:
        raise RuntimeError(
            f"{ref_condition}: evaluable signature order changed "
            "during subsetting."
        )

    signature_df.to_csv(
        signature_output_dir
        / f"{ref_condition}_holdout_validation_evaluable_factors.csv"
    )

    excluded_celltypes_df.to_csv(
        condition_output_dir
        / "excluded_celltypes.csv",
        index=False,
    )

    factor_manifest_df.to_csv(
        condition_output_dir
        / "validation_factor_manifest.csv",
        index=False,
    )

    # -------------------------------------------------------------------------
    # 4. Align genes across signatures, holdout cells, and spatial data
    # -------------------------------------------------------------------------

    signature_df, adata_holdout_aligned = align_genes(
        signature_df,
        adata_holdout,
    )

    adata_wta_cond = adata_wta[
        adata_wta.obs["disease_status"] == spatial_condition
    ].copy()

    if adata_wta_cond.n_obs == 0:
        raise ValueError(
            f"No spatial ROIs found for {spatial_condition}."
        )

    signature_df, adata_wta_cond_aligned = align_genes(
        signature_df,
        adata_wta_cond,
    )

    shared_genes = (
        adata_wta_cond_aligned
        .var_names
        .astype(str)
    )

    adata_holdout_aligned = (
        adata_holdout_aligned[
            :,
            shared_genes,
        ]
        .copy()
    )

    signature_df = (
        signature_df
        .loc[shared_genes]
        .copy()
    )

    signature_gene_order = (
        signature_df.index
        .astype(str)
        .to_numpy()
    )

    holdout_gene_order = (
        adata_holdout_aligned.var_names
        .astype(str)
        .to_numpy()
    )

    spatial_gene_order = (
        adata_wta_cond_aligned.var_names
        .astype(str)
        .to_numpy()
    )

    if not np.array_equal(
        signature_gene_order,
        holdout_gene_order,
    ):
        raise RuntimeError(
            f"{ref_condition}: signature and holdout gene order differ."
        )

    if not np.array_equal(
        signature_gene_order,
        spatial_gene_order,
    ):
        raise RuntimeError(
            f"{ref_condition}: signature and spatial gene order differ."
        )

    factor_labels = (
        signature_df.columns
        .astype(str)
        .tolist()
    )

    if factor_labels != evaluable_celltypes:
        raise RuntimeError(
            f"{ref_condition}: factor order changed during gene alignment."
        )

    print("Shared genes:", len(shared_genes))
    print("Spatial donor ROIs:", adata_wta_cond_aligned.n_obs)
    print("Evaluable factors:", len(factor_labels))

    # -------------------------------------------------------------------------
    # 5. Generate held-out synthetic pseudobulk ROIs
    # -------------------------------------------------------------------------

    synthetic = generate_synthetic_dataset(
        adata_holdout_aligned,
        celltype_col="celltype",
        adata_real_donor=adata_wta_cond_aligned,
        celltypes=factor_labels,
        n_mixtures=N_MIXTURES,
        n_cells_per_mixture=N_CELLS_PER_MIXTURE,
        concentration=DIRICHLET_CONCENTRATION,
        seed=RANDOM_SEED,
    )

    validate_synthetic_dataset(
        synthetic=synthetic,
        factor_labels=factor_labels,
        condition_name=ref_condition,
    )

    # -------------------------------------------------------------------------
    # 6. Fit SpaceJam on the synthetic ROIs
    # -------------------------------------------------------------------------

    pyro.clear_param_store()
    torch.cuda.empty_cache()

    model_inputs = {
        "cell_state_mat": (
            signature_df
            .to_numpy(dtype=np.float32)
        ),
        "X_data": np.asarray(
            synthetic["X_data"],
            dtype=np.float32,
        ),
        "Y_data": np.asarray(
            synthetic["Y_data"],
            dtype=np.float32,
        ),
        "spot2sample_mat": np.asarray(
            synthetic["spot2sample_mat"],
            dtype=np.float32,
        ),
        "cell_number_prior": CELL_NUMBER_PRIOR,
    }

    model = LocationModelPyro(
        **model_inputs,
        device=DEVICE,
    )

    loss_history = model.fit(
        n_steps=SPACEJAM_STEPS,
        lr=SPACEJAM_LEARNING_RATE,
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
        condition_output_dir
        / "training_loss.csv",
        index=False,
    )

    # -------------------------------------------------------------------------
    # 7. Extract posterior proportions
    # -------------------------------------------------------------------------

    synthetic_roi_ids = [
        f"synthetic_{index:04d}"
        for index in range(
            synthetic["X_data"].shape[0]
        )
    ]

    inferred_df = extract_posterior_proportions(
        model=model,
        synthetic_roi_ids=synthetic_roi_ids,
        factor_labels=factor_labels,
        num_samples=POSTERIOR_PREDICTIVE_SAMPLES,
    )

    ground_truth_df = (
        synthetic["ground_truth_df"]
        .copy()
    )

    if ground_truth_df.columns.tolist() != inferred_df.columns.tolist():
        raise RuntimeError(
            f"{ref_condition}: ground-truth and inferred factor "
            "orders differ."
        )

    if ground_truth_df.shape != inferred_df.shape:
        raise RuntimeError(
            f"{ref_condition}: ground-truth shape "
            f"{ground_truth_df.shape} differs from inferred shape "
            f"{inferred_df.shape}."
        )

    # -------------------------------------------------------------------------
    # 8. Evaluate recovery
    # -------------------------------------------------------------------------

    metrics_df = evaluate_recovery(
        ground_truth_df,
        inferred_df,
    )

    ground_truth_out = (
        ground_truth_df
        .copy()
    )

    ground_truth_out.insert(
        0,
        "roi_id",
        synthetic_roi_ids,
    )

    ground_truth_out.insert(
        1,
        "donor_roi_id",
        synthetic["donor_roi_ids"],
    )

    inferred_out = inferred_df.copy()
    inferred_out.index.name = "roi_id"

    metrics_df.to_csv(
        condition_output_dir
        / "recovery_metrics.csv",
        index=False,
    )

    ground_truth_out.to_csv(
        condition_output_dir
        / "ground_truth_proportions.csv",
        index=False,
    )

    inferred_out.to_csv(
        condition_output_dir
        / "inferred_proportions.csv",
    )

    plot_recovery_scatter(
        ground_truth_df,
        inferred_df,
        save_file=str(
            condition_output_dir
            / "recovery_scatter.pdf"
        ),
    )

    plot_recovery_bias_by_abundance(
        ground_truth_df,
        inferred_df,
        save_file=str(
            condition_output_dir
            / "recovery_bias.pdf"
        ),
    )

    run_metadata = {
        "reference_condition": ref_condition,
        "spatial_condition": spatial_condition,
        "n_reference_cells": int(
            adata_ref_cond.n_obs
        ),
        "n_train_cells": int(
            adata_train.n_obs
        ),
        "n_holdout_cells": int(
            adata_holdout.n_obs
        ),
        "n_trained_factors": int(
            len(model_order)
        ),
        "n_evaluable_factors": int(
            len(factor_labels)
        ),
        "n_excluded_factors": int(
            len(model_order) - len(factor_labels)
        ),
        "trained_factor_order": model_order,
        "evaluable_factor_order": factor_labels,
        "excluded_factors": (
            excluded_celltypes_df["celltype"]
            .tolist()
        ),
        "n_shared_genes": int(
            len(shared_genes)
        ),
        "n_spatial_donor_rois": int(
            adata_wta_cond_aligned.n_obs
        ),
        "n_synthetic_mixtures": N_MIXTURES,
        "n_cells_per_mixture": N_CELLS_PER_MIXTURE,
        "dirichlet_concentration": DIRICHLET_CONCENTRATION,
        "holdout_fraction": HOLDOUT_FRACTION,
        "minimum_holdout_cells": MIN_HOLDOUT_CELLS,
        "random_seed": RANDOM_SEED,
        "regression_max_epochs": REGRESSION_MAX_EPOCHS,
        "regression_posterior_samples": REGRESSION_POSTERIOR_SAMPLES,
        "spacejam_steps": SPACEJAM_STEPS,
        "spacejam_learning_rate": SPACEJAM_LEARNING_RATE,
        "posterior_predictive_samples": POSTERIOR_PREDICTIVE_SAMPLES,
        "device": DEVICE,
    }

    with open(
        condition_output_dir
        / "run_metadata.json",
        "w",
    ) as handle:
        json.dump(
            run_metadata,
            handle,
            indent=2,
        )

    print("\nRecovery metrics:")
    print(metrics_df)

    results_by_condition[
        ref_condition
    ] = {
        "metrics": metrics_df,
        "ground_truth": ground_truth_df,
        "inferred": inferred_df,
        "n_trained_factors": len(model_order),
        "n_evaluable_factors": len(factor_labels),
    }

    # -------------------------------------------------------------------------
    # 9. Clean up GPU and memory before next condition
    # -------------------------------------------------------------------------

    del model
    del model_regression
    del adata_ref_cond
    del adata_train
    del adata_holdout
    del adata_holdout_aligned
    del adata_wta_cond
    del adata_wta_cond_aligned
    del synthetic

    pyro.clear_param_store()

    gc.collect()
    torch.cuda.empty_cache()


# =============================================================================
# Combined final outputs
# =============================================================================

combined_metrics_rows = []
condition_summary_rows = []

for condition, condition_results in results_by_condition.items():

    metrics = (
        condition_results["metrics"]
        .copy()
    )

    metrics.insert(
        0,
        "condition",
        condition,
    )

    combined_metrics_rows.append(
        metrics
    )

    condition_summary_rows.append(
        {
            "condition": condition,
            "n_trained_factors": (
                condition_results["n_trained_factors"]
            ),
            "n_evaluable_factors": (
                condition_results["n_evaluable_factors"]
            ),
        }
    )

combined_metrics = pd.concat(
    combined_metrics_rows,
    axis=0,
    ignore_index=True,
)

combined_metrics.to_csv(
    OUTPUT_ROOT
    / "combined_recovery_metrics.csv",
    index=False,
)

condition_summary_df = pd.DataFrame(
    condition_summary_rows
)

condition_summary_df.to_csv(
    OUTPUT_ROOT
    / "condition_validation_summary.csv",
    index=False,
)

print("\n" + "=" * 80)
print("Pseudobulk validation completed successfully")
print("=" * 80)

print("Conditions:", list(results_by_condition))

print("\nCondition validation summary:")
print(condition_summary_df)

print(
    "\nCombined recovery metrics saved to:",
    OUTPUT_ROOT
    / "combined_recovery_metrics.csv",
)

print(
    "Condition summary saved to:",
    OUTPUT_ROOT
    / "condition_validation_summary.csv",
)
