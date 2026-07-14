import os
import sys
import torch
import scanpy as sc
import pandas as pd
import numpy as np

from python.regression_utils import (
    load_seurat_mtx_as_anndata,
    standardize_regression_obs,
    round_counts_layer,
    minimal_filter_anndata,
    disable_lightning_mpi_detection,
    train_regression_models_by_condition,
    save_regression_models,
    load_and_export_regression_posterior,
)

print("Python executable:", sys.executable)
print("CUDA available:", torch.cuda.is_available())


# ============================================================
# Paths
# ============================================================

data_path = (
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025"
)

project_dir = data_path / "multiomic-spatial-integration"

rresults = project_dir / "results" / "regression_model"

rresults.mkdir(
    parents=True,
    exist_ok=True
)

print("Regression outputs will be saved to:")
print(rresults)

# ============================================================
# Load snRNA-seq reference data
# ============================================================

adata = load_seurat_mtx_as_anndata(
    counts_mtx=str(data_path / "AD_CAA_counts.mtx"),
    metadata_csv=str(data_path / "AD_CAA_meta.csv"),
    genes_csv=str(data_path / "AD_CAA_genes.csv"),
    gene_col="gene_id",
    transpose_counts=True,
    counts_layer="counts",
)

# ============================================================
# Standardize metadata
# ============================================================

adata = standardize_regression_obs(
    adata,
    celltype_col="New_Idents",
    new_celltype_col="celltype",
    experiment_col="Experiment",
    experiment_value="Batch1",
)

adata = round_counts_layer(
    adata,
    counts_layer="counts",
    replace_x=True,
)

adata = minimal_filter_anndata(
    adata,
    min_genes=1,
    min_cells=1,
    counts_layer="counts",
)

print(adata)


# ============================================================
# Train regression models
# ============================================================

disable_lightning_mpi_detection()

trained_models = train_regression_models_by_condition(
    adata,
    condition_col="FDX",
    conditions=["AD+CAA", "Control"],
    labels_key="celltype",
    batch_key="Experiment",
    layer="counts",
    max_epochs=100,
    batch_size=1024,
    lr=0.01,
    accelerator="gpu",
)

save_regression_models(
    trained_models,
    output_dir=str(rresults)
)


# ============================================================
# Create condition-specific AnnData objects
# ============================================================

adata_AD = adata[
    adata.obs["FDX"] == "AD+CAA"
].copy()

adata_CTRL = adata[
    adata.obs["FDX"] == "Control"
].copy()

print("AD+CAA cells:", adata_AD.n_obs)
print("Control cells:", adata_CTRL.n_obs)


# ============================================================
# Load trained models and export posterior means
# ============================================================

adata_AD, mod_AD = load_and_export_regression_posterior(
    adata_AD,
    model_dir=str(rresults / "AD+CAA_regression_model"),
    num_samples=1000,
    batch_size=2500
)

adata_CTRL, mod_CTRL = load_and_export_regression_posterior(
    adata_CTRL,
    model_dir=str(rresults / "Control_regression_model"),
    num_samples=1000,
    batch_size=2500
)


# ============================================================
# Recover authoritative cell-type order from model registry
# ============================================================

def get_model_celltype_order(model):
    """
    Recover the exact categorical label order used internally by
    the trained cell2location/scvi regression model.
    """

    registry_keys_to_try = [
        "labels",
        "celltype",
    ]

    state = None
    successful_key = None

    for key in registry_keys_to_try:
        try:
            state = model.adata_manager.get_state_registry(key)
            successful_key = key
            break
        except Exception:
            continue

    if state is None:
        raise ValueError(
            "Could not recover the label registry from the trained model. "
            f"Tried keys: {registry_keys_to_try}"
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
                "The model state registry was found, but no "
                "'categorical_mapping' was available. "
                f"Registry key: {successful_key}; registry contents: {state}"
            ) from exc

    labels = [str(label) for label in labels]

    if len(labels) == 0:
        raise ValueError(
            "The recovered model cell-type order is empty."
        )

    return labels


ad_model_order = get_model_celltype_order(mod_AD)
ctrl_model_order = get_model_celltype_order(mod_CTRL)

print("\nAD+CAA model factor order:")
print(ad_model_order)

print("\nControl model factor order:")
print(ctrl_model_order)

print(
    "\nSame factor order across conditions:",
    ad_model_order == ctrl_model_order,
)

print("Number of AD+CAA factors:", len(ad_model_order))
print("Number of Control factors:", len(ctrl_model_order))


# This should be true for directly comparable SpaceJam models.
if ad_model_order != ctrl_model_order:
    raise ValueError(
        "AD+CAA and Control regression models have different internal "
        "factor orders. Do not proceed to SpaceJam until this is resolved."
    )


# ============================================================
# Corrected inferred-signature export
# ============================================================

def export_inferred_signatures_corrected(
    adata,
    model,
    covariate_col,
    model_name,
    output_folder,
):
    """
    Export inferred signatures using the exact categorical factor order
    stored in the trained model registry.

    This avoids attaching adata.obs[covariate_col].unique() labels to a
    posterior matrix whose columns follow the model's internal order.
    """

    os.makedirs(output_folder, exist_ok=True)

    required_uns = [
        "mod",
    ]

    for key in required_uns:
        if key not in adata.uns:
            raise KeyError(
                f"{model_name}: adata.uns['{key}'] was not found."
            )

    try:
        mu = (
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

    mu = np.asarray(mu)

    model_celltypes = get_model_celltype_order(model)

    if mu.ndim != 2:
        raise ValueError(
            f"{model_name}: inferred signature array must be 2-dimensional; "
            f"received shape {mu.shape}."
        )

    if mu.shape[0] != adata.n_vars:
        raise ValueError(
            f"{model_name}: inferred signature matrix has {mu.shape[0]} "
            f"genes, but AnnData has {adata.n_vars} genes."
        )

    if mu.shape[1] != len(model_celltypes):
        raise ValueError(
            f"{model_name}: inferred signature matrix contains "
            f"{mu.shape[1]} factors, but the model registry contains "
            f"{len(model_celltypes)} labels."
        )

    inferred_df = pd.DataFrame(
        mu,
        index=adata.var_names.astype(str),
        columns=model_celltypes,
    )

    inferred_df.index.name = "gene"

    inferred_file = os.path.join(
        output_folder,
        f"{model_name}_inferred_signatures.csv",
    )

    inferred_df.to_csv(inferred_file)

    # Analytical mean expression per cell type, placed in the same
    # authoritative factor order.
    average_df = (
        adata.to_df()
        .join(adata.obs[[covariate_col]])
        .groupby(covariate_col, observed=True)
        .mean()
        .T
    )

    missing_average_labels = [
        label
        for label in model_celltypes
        if label not in average_df.columns
    ]

    if missing_average_labels:
        raise ValueError(
            f"{model_name}: analytical mean expression is missing "
            f"cell type(s): {missing_average_labels}"
        )

    average_df = average_df.reindex(
        columns=model_celltypes
    )

    average_df.index.name = "gene"

    average_file = os.path.join(
        output_folder,
        f"{model_name}_mean_cluster_expr.csv",
    )

    average_df.to_csv(average_file)

    # Hard validation that the export order matches the model.
    if inferred_df.columns.tolist() != model_celltypes:
        raise RuntimeError(
            f"{model_name}: exported inferred-signature columns do not "
            "match the model registry order."
        )

    if average_df.columns.tolist() != model_celltypes:
        raise RuntimeError(
            f"{model_name}: exported analytical-average columns do not "
            "match the model registry order."
        )

    print(
        f"\n{model_name}: saved inferred signatures -> "
        f"{inferred_file}"
    )

    print(
        f"{model_name}: saved analytical averages -> "
        f"{average_file}"
    )

    print(
        f"{model_name}: validated {len(model_celltypes)} factors "
        "against the model registry."
    )

    return inferred_df, average_df


# ============================================================
# Export corrected signatures
# ============================================================

inf_AD, avg_AD = export_inferred_signatures_corrected(
    adata=adata_AD,
    model=mod_AD,
    covariate_col="celltype",
    model_name="AD+CAA",
    output_folder=str(rresults)
)

inf_CTRL, avg_CTRL = export_inferred_signatures_corrected(
    adata=adata_CTRL,
    model=mod_CTRL,
    covariate_col="celltype",
    model_name="Control",
    output_folder=str(rresults)
)


# ============================================================
# Extra validation 1:
# exported columns must match model registry exactly
# ============================================================

assert inf_AD.columns.tolist() == ad_model_order
assert inf_CTRL.columns.tolist() == ctrl_model_order

assert avg_AD.columns.tolist() == ad_model_order
assert avg_CTRL.columns.tolist() == ctrl_model_order

print(
    "\nPASS: all exported signature and average-expression columns "
    "match their model registry order."
)


# ============================================================
# Extra validation 2:
# same-position cross-condition signature concordance
# ============================================================

shared_genes = inf_AD.index.intersection(
    inf_CTRL.index
)

ad_test = inf_AD.loc[shared_genes]
ctrl_test = inf_CTRL.loc[shared_genes]

same_position_results = []

for index, celltype in enumerate(ad_model_order):

    ad_label = ad_test.columns[index]
    ctrl_label = ctrl_test.columns[index]

    rho = ad_test.iloc[:, index].corr(
        ctrl_test.iloc[:, index],
        method="spearman",
    )

    same_position_results.append(
        {
            "factor_index": index,
            "AD_label": ad_label,
            "Control_label": ctrl_label,
            "labels_match": ad_label == ctrl_label,
            "same_position_rho": rho,
        }
    )

same_position_df = pd.DataFrame(
    same_position_results
)

same_position_file = os.path.join(
    rresults,
    "signature_same_position_validation.csv",
)

same_position_df.to_csv(
    same_position_file,
    index=False,
)

print("\nSame-position validation summary:")
print(
    same_position_df[
        [
            "factor_index",
            "AD_label",
            "Control_label",
            "labels_match",
            "same_position_rho",
        ]
    ]
)

print(
    "\nMatching labels:",
    same_position_df["labels_match"].sum(),
    "of",
    len(same_position_df),
)

print(
    "Median same-position Spearman correlation:",
    same_position_df["same_position_rho"].median(),
)

print(
    "Minimum same-position Spearman correlation:",
    same_position_df["same_position_rho"].min(),
)


# Every corresponding factor should now have the same biological label.
if not same_position_df["labels_match"].all():
    mismatches = same_position_df.loc[
        ~same_position_df["labels_match"]
    ]

    raise RuntimeError(
        "Corrected exports still contain cross-condition label mismatches:\n"
        f"{mismatches.to_string(index=False)}"
    )


# ============================================================
# Extra validation 3:
# each Control factor should best match the same-named AD factor
# ============================================================

cross_condition_cor = pd.DataFrame(
    index=inf_AD.columns,
    columns=inf_CTRL.columns,
    dtype=float,
)

for ad_celltype in inf_AD.columns:
    for ctrl_celltype in inf_CTRL.columns:
        cross_condition_cor.loc[
            ad_celltype,
            ctrl_celltype,
        ] = inf_AD.loc[
            shared_genes,
            ad_celltype,
        ].corr(
            inf_CTRL.loc[
                shared_genes,
                ctrl_celltype,
            ],
            method="spearman",
        )

best_match_results = []

for ctrl_celltype in inf_CTRL.columns:

    best_ad_celltype = (
        cross_condition_cor[ctrl_celltype]
        .idxmax()
    )

    best_rho = cross_condition_cor.loc[
        best_ad_celltype,
        ctrl_celltype,
    ]

    best_match_results.append(
        {
            "Control_label": ctrl_celltype,
            "best_AD_label": best_ad_celltype,
            "best_rho": best_rho,
            "same_name_best_match": (
                ctrl_celltype == best_ad_celltype
            ),
        }
    )

best_match_df = pd.DataFrame(
    best_match_results
)

cross_condition_cor.to_csv(
    os.path.join(
        rresults,
        "signature_cross_condition_correlation_corrected.csv",
    )
)

best_match_df.to_csv(
    os.path.join(
        rresults,
        "signature_cross_condition_best_matches_corrected.csv",
    ),
    index=False,
)

n_same_name = best_match_df[
    "same_name_best_match"
].sum()

print(
    "\nSame-name best matches:",
    n_same_name,
    "of",
    len(best_match_df),
)

print(
    best_match_df.sort_values(
        "best_rho",
        ascending=False,
    )
)


if n_same_name != len(best_match_df):
    mismatches = best_match_df.loc[
        ~best_match_df["same_name_best_match"]
    ]

    raise RuntimeError(
        "Some corrected Control signatures do not best match their "
        "same-named AD signatures:\n"
        f"{mismatches.to_string(index=False)}"
    )


print(
    "\nAll regression signature-order validation checks passed."
)
