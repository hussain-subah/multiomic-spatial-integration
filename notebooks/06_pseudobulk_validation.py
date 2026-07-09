# 06 — Pseudobulk Validation of SpaceJam / Cell2Location Deconvolution

#This notebook validates the deconvolution model's ability to recover known
#cell-type composition, using synthetic bulk ("pseudobulk") expression
#profiles built from held-out snRNA-seq reference cells.

#Methodology:
#- The reference is split per cell type into train/holdout cells.
#- The cell-type signature matrix is re-derived from train cells only.
#- Synthetic ROIs are built from held-out cells only (never seen by the
#  signature matrix), so recovery isn't circular.
#- Synthetic ROIs borrow negative-probe counts, sequencing depth, and batch
#  identity from real GeoMx ROIs (data/CAA-AD_AnnData.h5ad), so their
#  technical noise is realistic without an invented noise model.
#- Inferred proportions are compared against the known synthetic ground
#  truth.
import os
import torch
import scanpy as sc
import numpy as np
import pandas as pd

from python.regression_utils import (
    load_seurat_mtx_as_anndata,
    standardize_regression_obs,
    round_counts_layer,
    minimal_filter_anndata,
    disable_lightning_mpi_detection,
    train_regression_model,
    export_inferred_signatures,
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

# ============================================================
# Config
# ============================================================

data_path = "data/"
output_root = "results/pseudobulk_validation/"
os.makedirs(output_root, exist_ok=True)

condition_map = {
    "AD+CAA": "AD-CAA",
    "Control": "Control",
}

cell_number_prior = {
    "cells_per_spot": 8.0,
    "factors_per_spot": 7.0,
    "combs_per_spot": 2.5,
    "factors_per_combs": 3.0,
    "cells_mean_var_ratio": 1.0,
    "factors_mean_var_ratio": 1.0,
    "combs_mean_var_ratio": 1.0,
}

# Smoke-test settings
N_MIXTURES = 100
N_CELLS_PER_MIXTURE = 200
DIRICHLET_CONCENTRATION = 0.3
N_STEPS = 9000
LEARNING_RATE = 0.002

# ============================================================
# Load snRNA-seq reference
# ============================================================

adata_ref = load_seurat_mtx_as_anndata(
    counts_mtx=f"{data_path}/AD_CAA_counts.mtx",
    metadata_csv=f"{data_path}/AD_CAA_meta.csv",
    genes_csv=f"{data_path}/AD_CAA_genes.csv",
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

disable_lightning_mpi_detection()

# ============================================================
# Load real GeoMx WTA AnnData
# ============================================================

adata_wta = sc.read_h5ad(f"{data_path}/CAA-AD_AnnData.h5ad")

if "counts" not in adata_wta.layers:
    adata_wta.layers["counts"] = adata_wta.X.copy()

# ============================================================
# Run validation
# ============================================================

results_by_condition = {}

for ref_condition, spatial_condition in condition_map.items():

    print(f"\n==== {ref_condition} ====")

    condition_output_dir = f"{output_root}/{ref_condition}/"
    os.makedirs(condition_output_dir, exist_ok=True)

    # ------------------------------------------------------------
    # 1. Condition subset + train/holdout split
    # ------------------------------------------------------------

    adata_ref_cond = adata_ref[adata_ref.obs["FDX"] == ref_condition].copy()

    adata_train, adata_holdout = split_reference_cells(
        adata_ref_cond,
        celltype_col="celltype",
        holdout_frac=0.3,
        min_holdout_cells=10,
        seed=0,
    )

    print(f"Train cells: {adata_train.n_obs}")
    print(f"Holdout cells: {adata_holdout.n_obs}")

    # ------------------------------------------------------------
    # 2. Re-derive signatures from train cells only
    # ------------------------------------------------------------

    model_regression = train_regression_model(
        adata_train,
        labels_key="celltype",
        batch_key="Experiment",
        layer="counts",
        max_epochs= 100, #100
        batch_size=1024,
        lr=0.01,
        accelerator="gpu", #gpu
    )

    model_regression.export_posterior(
        adata_train,
        sample_kwargs={
            "num_samples": 1000, #100
            "batch_size": 2500, #500
        },
    )

    signature_df, _ = export_inferred_signatures(
        adata_train,
        covariate_col="celltype",
        model_name=f"{ref_condition}_holdout_validation",
        output_folder=f"{condition_output_dir}/signature/",
    )

    # ------------------------------------------------------------
    # 3. Align genes across signature, holdout reference, and WTA
    # ------------------------------------------------------------

    signature_df, adata_holdout_aligned = align_genes(
        signature_df,
        adata_holdout,
    )

    adata_wta_cond = adata_wta[
        adata_wta.obs["disease_status"] == spatial_condition
    ].copy()

    signature_df, adata_wta_cond_aligned = align_genes(
        signature_df,
        adata_wta_cond,
    )

    shared_genes = adata_wta_cond_aligned.var_names

    adata_holdout_aligned = adata_holdout_aligned[:, shared_genes].copy()
    signature_df = signature_df.loc[shared_genes]

    print(f"Shared genes: {len(shared_genes)}")

    # ------------------------------------------------------------
    # 4. Generate synthetic pseudobulk ROIs
    # ------------------------------------------------------------

    synthetic = generate_synthetic_dataset(
        adata_holdout_aligned,
        celltype_col="celltype",
        adata_real_donor=adata_wta_cond_aligned,
        celltypes=list(signature_df.columns),
        n_mixtures=N_MIXTURES,
        n_cells_per_mixture=N_CELLS_PER_MIXTURE,
        concentration=DIRICHLET_CONCENTRATION,
        seed=0,
    )

    # ------------------------------------------------------------
    # 5. Train SpaceJam / Cell2location model on synthetic ROIs
    # ------------------------------------------------------------

    inputs = {
        "cell_state_mat": signature_df.values.astype("float32"),
        "X_data": synthetic["X_data"].astype("float32"),
        "Y_data": synthetic["Y_data"].astype("float32"),
        "spot2sample_mat": synthetic["spot2sample_mat"].astype("float32"),
        "cell_number_prior": cell_number_prior,
    }

    model = LocationModelPyro(**inputs)

    model.fit(
        n_steps=N_STEPS,
        lr=LEARNING_RATE
    )

    if not hasattr(model, "cell_types"):
        model.cell_types = list(signature_df.columns)

    # ------------------------------------------------------------
    # 6. Extract inferred proportions
    # ------------------------------------------------------------


    synthetic_roi_ids = [
        f"synthetic_{i:04d}" for i in range(synthetic["X_data"].shape[0])
    ]

    posterior = model.posterior_predictive(num_samples=100)

    spot_factors = posterior["spot_factors"]

    # Convert torch tensor -> numpy
    if hasattr(spot_factors, "detach"):
        spot_factors = spot_factors.detach().cpu().numpy()

    # Average across posterior samples
    abundance = spot_factors.mean(axis=0)

    # Convert to proportions
    inferred_df = pd.DataFrame(
        abundance,
        index=synthetic_roi_ids,
        columns=list(signature_df.columns),
    )

    inferred_df = inferred_df.div(
        inferred_df.sum(axis=1),
        axis=0
    )
    # ------------------------------------------------------------
    # 7. Evaluate recovery
    # ------------------------------------------------------------

    metrics_df = evaluate_recovery(
        synthetic["ground_truth_df"],
        inferred_df,
    )

    ground_truth_out = synthetic["ground_truth_df"].copy()
    ground_truth_out.insert(0, "roi_id", synthetic_roi_ids)
    ground_truth_out.insert(1, "donor_roi_id", synthetic["donor_roi_ids"])

    metrics_df.to_csv(
        f"{condition_output_dir}/recovery_metrics.csv",
        index=False,
    )

    ground_truth_out.to_csv(
        f"{condition_output_dir}/ground_truth_proportions.csv",
        index=False,
    )

    inferred_df.to_csv(
        f"{condition_output_dir}/inferred_proportions.csv"
    )

    plot_recovery_scatter(
        synthetic["ground_truth_df"],
        inferred_df,
        save_file=f"{condition_output_dir}/recovery_scatter.pdf",
    )

    plot_recovery_bias_by_abundance(
        synthetic["ground_truth_df"],
        inferred_df,
        save_file=f"{condition_output_dir}/recovery_bias.pdf",
    )

    print(metrics_df)

    results_by_condition[ref_condition] = {
        "metrics": metrics_df,
        "ground_truth": synthetic["ground_truth_df"],
        "inferred": inferred_df,
    }

results_by_condition
