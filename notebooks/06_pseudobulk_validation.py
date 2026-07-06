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

import scanpy as sc

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
    train_spacejam_model,
    extract_posterior_abundance,
)
from python.pseudobulk_validation_utils import (
    evaluate_recovery,
    plot_recovery_scatter,
    plot_recovery_bias_by_abundance,
)

## Config -- adjust paths for your environment

data_path = "data/"
output_root = "results/pseudobulk_validation/"

## Load snRNA-seq reference (same inputs as notebooks/02_regression_signatures.py)

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

adata_ref = round_counts_layer(adata_ref, counts_layer="counts", replace_x=True)
adata_ref = minimal_filter_anndata(adata_ref, min_genes=1, min_cells=1, counts_layer="counts")

disable_lightning_mpi_detection()

## Load real GeoMx WTA spatial AnnData (technical donor for synthetic ROIs only)

adata_wta = sc.read_h5ad(f"{data_path}/CAA-AD_AnnData.h5ad")

if "counts" not in adata_wta.layers:
    adata_wta.layers["counts"] = adata_wta.X.copy()

## Reference FDX labels ("AD+CAA"/"Control") and spatial disease_status
## labels ("AD-CAA"/"Control") use different strings for the same conditions.

condition_map = {
    "AD+CAA": "AD-CAA",
    "Control": "Control",
}

# Matches the cell_number_prior actually used for real-data training in
# notebooks/Examples/Cell2Location + SpaceJam.ipynb
cell_number_prior = {
    "cells_per_spot": 8.0,
    "factors_per_spot": 7.0,
    "combs_per_spot": 2.5,
    "factors_per_combs": 3.0,
    "cells_mean_var_ratio": 1.0,
    "factors_mean_var_ratio": 1.0,
    "combs_mean_var_ratio": 1.0,
}

N_MIXTURES = 100
N_CELLS_PER_MIXTURE = 200
DIRICHLET_CONCENTRATION = 0.3
N_STEPS = 9000
LEARNING_RATE = 0.002

results_by_condition = {}

for ref_condition, spatial_condition in condition_map.items():

    print(f"\n==== {ref_condition} ====")

    condition_output_dir = f"{output_root}/{ref_condition}/"

    ## 1. Condition subset + held-out split of the reference

    adata_ref_cond = adata_ref[adata_ref.obs["FDX"] == ref_condition].copy()

    adata_train, adata_holdout = split_reference_cells(
        adata_ref_cond,
        celltype_col="celltype",
        holdout_frac=0.3,
        min_holdout_cells=10,
        seed=0,
    )

    print(f"Train cells: {adata_train.n_obs}, holdout cells: {adata_holdout.n_obs}")

    ## 2. Re-derive the signature matrix from train cells only (never holdout)

    model_regression = train_regression_model(
        adata_train,
        labels_key="celltype",
        batch_key="Experiment",
        layer="counts",
        max_epochs=100,
        batch_size=1024,
        lr=0.01,
        accelerator="gpu",
    )

    model_regression.export_posterior(
        adata_train,
        sample_kwargs={"num_samples": 1000, "batch_size": 2500},
    )

    signature_df, _ = export_inferred_signatures(
        adata_train,
        covariate_col="celltype",
        model_name=f"{ref_condition}_holdout_validation",
        output_folder=f"{condition_output_dir}/signature/",
    )

    ## 3. Align genes: train-derived signature <-> holdout reference <-> real WTA

    signature_df, adata_holdout_aligned = align_genes(signature_df, adata_holdout)

    adata_wta_cond = adata_wta[adata_wta.obs["disease_status"] == spatial_condition].copy()
    signature_df, adata_wta_cond_aligned = align_genes(signature_df, adata_wta_cond)

    # adata_holdout_aligned may still hold genes dropped by the second
    # align_genes() call above -- pin everything to the final shared set.
    shared_genes = adata_wta_cond_aligned.var_names
    adata_holdout_aligned = adata_holdout_aligned[:, shared_genes].copy()
    signature_df = signature_df.loc[shared_genes]

    print(f"Shared genes (signature x holdout x real WTA): {len(shared_genes)}")

    ## 4. Generate the synthetic pseudobulk dataset

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

    ## 5. Train the deconvolution model on the synthetic ROIs

    cell_state_mat = signature_df.values.astype("float32")

    model, loss_hist = train_spacejam_model(
        cell_state_mat,
        synthetic["X_data"],
        synthetic["Y_data"],
        synthetic["spot2sample_mat"],
        n_steps=N_STEPS,
        lr=LEARNING_RATE,
        cell_number_prior=cell_number_prior,
    )

    ## 6. Extract inferred proportions and compare to known ground truth

    synthetic_roi_ids = [f"synthetic_{i:04d}" for i in range(synthetic["X_data"].shape[0])]

    inferred_df = extract_posterior_abundance(
        celltype_labels=synthetic["celltypes"],
        roi_names=synthetic_roi_ids,
        normalize=True,
    )

    metrics_df = evaluate_recovery(synthetic["ground_truth_df"], inferred_df)

    ground_truth_out = synthetic["ground_truth_df"].copy()
    ground_truth_out.insert(0, "roi_id", synthetic_roi_ids)
    ground_truth_out.insert(1, "donor_roi_id", synthetic["donor_roi_ids"])

    metrics_df.to_csv(f"{condition_output_dir}/recovery_metrics.csv", index=False)
    ground_truth_out.to_csv(f"{condition_output_dir}/ground_truth_proportions.csv", index=False)
    inferred_df.to_csv(f"{condition_output_dir}/inferred_proportions.csv")

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
