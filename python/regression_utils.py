"""
Utilities for Cell2Location regression-based cell-type signature inference.

These functions prepare snRNA-seq AnnData objects, train RegressionModel
objects, export posterior means, and save inferred cell-type signatures.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc


def load_seurat_mtx_as_anndata(
    counts_mtx,
    metadata_csv,
    genes_csv,
    gene_col="gene_id",
    metadata_index_col=0,
    transpose_counts=True,
    counts_layer="counts",
):
    """
    Load Seurat-exported counts, metadata, and genes into AnnData.

    Parameters
    ----------
    counts_mtx : str
        Path to Matrix Market count matrix.
    metadata_csv : str
        Path to cell metadata CSV.
    genes_csv : str
        Path to gene annotation CSV.
    gene_col : str
        Column containing gene IDs/names.
    metadata_index_col : int
        Index column for metadata.
    transpose_counts : bool
        Whether to transpose the matrix after reading.
    counts_layer : str
        Layer name for raw counts.

    Returns
    -------
    AnnData
    """
    from scipy.io import mmread
    import anndata

    counts = mmread(counts_mtx)

    if transpose_counts:
        counts = counts.T

    counts = counts.tocsr()

    meta = pd.read_csv(metadata_csv, index_col=metadata_index_col)
    genes = pd.read_csv(genes_csv)
    gene_ids = genes[gene_col].astype(str).values

    adata = anndata.AnnData(X=counts)
    adata.obs = meta
    adata.var = pd.DataFrame(index=gene_ids)
    adata.var.index.name = "gene"
    adata.var_names = gene_ids

    adata.layers[counts_layer] = adata.X.copy()

    return adata


def standardize_regression_obs(
    adata,
    celltype_col="New_Idents",
    new_celltype_col="celltype",
    experiment_col="Experiment",
    experiment_value="Batch1",
):
    """
    Standardize AnnData obs columns for Cell2Location RegressionModel.

    Adds/renames:
    - cell type labels
    - experiment/batch column
    """
    if celltype_col not in adata.obs.columns and new_celltype_col not in adata.obs.columns:
        raise ValueError(
            f"Neither {celltype_col} nor {new_celltype_col} found in adata.obs."
        )

    if celltype_col in adata.obs.columns and new_celltype_col not in adata.obs.columns:
        adata.obs.rename(columns={celltype_col: new_celltype_col}, inplace=True)

    adata.obs[experiment_col] = experiment_value

    return adata


def round_counts_layer(
    adata,
    counts_layer="counts",
    replace_x=True,
    dtype="int64",
):
    """
    Round non-integer counts to integers for Negative Binomial models.

    Useful after SoupX correction, which can produce float values.
    """
    X = adata.layers[counts_layer].copy()

    if sp.issparse(X):
        X.data[X.data < 0] = 0
        X.data = np.rint(X.data)
        X = X.astype(dtype)
    else:
        X = np.array(X, copy=True)
        X[X < 0] = 0
        X = np.rint(X).astype(dtype)

    adata.layers[counts_layer] = X

    if replace_x:
        adata.X = X.copy()

    return adata


def add_gene_detection_metrics(
    adata,
    counts_layer="counts",
):
    """
    Add n_cells and non-zero mean expression metrics to adata.var.
    """
    X = adata.layers[counts_layer]

    if sp.issparse(X):
        n_cells = (X > 0).sum(axis=0).A1
        total_counts = X.sum(axis=0).A1
    else:
        n_cells = (X > 0).sum(axis=0)
        total_counts = X.sum(axis=0)

    adata.var["n_cells"] = n_cells
    adata.var["nonz_mean"] = np.divide(
        total_counts,
        np.where(n_cells == 0, 1, n_cells)
    )

    return adata


def minimal_filter_anndata(
    adata,
    min_genes=1,
    min_cells=1,
    counts_layer="counts",
):
    """
    Apply minimal cell/gene filtering and recompute gene detection metrics.
    """
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    adata = add_gene_detection_metrics(
        adata,
        counts_layer=counts_layer
    )

    return adata


def disable_lightning_mpi_detection():
    """
    Disable Lightning MPI detection.

    Useful on some HPC systems where Lightning incorrectly detects MPI.
    """
    import lightning.fabric.plugins.environments.mpi as lightning_mpi

    def _no_mpi_detect():
        return False

    lightning_mpi.MPIEnvironment.detect = staticmethod(_no_mpi_detect)


def train_regression_model(
    adata,
    labels_key="celltype",
    batch_key="Experiment",
    layer="counts",
    max_epochs=100,
    batch_size=1024,
    lr=0.01,
    accelerator="gpu",
):
    """
    Train a Cell2Location RegressionModel on one AnnData object.

    Returns
    -------
    RegressionModel
    """
    from cell2location.models import RegressionModel

    RegressionModel.setup_anndata(
        adata,
        batch_key=batch_key,
        labels_key=labels_key,
        layer=layer
    )

    model = RegressionModel(adata)

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        lr=lr,
        accelerator=accelerator
    )

    return model


def train_regression_models_by_condition(
    adata,
    condition_col="FDX",
    conditions=None,
    labels_key="celltype",
    batch_key="Experiment",
    layer="counts",
    max_epochs=100,
    batch_size=1024,
    lr=0.01,
    accelerator="gpu",
):
    """
    Train separate regression models by condition.

    Returns
    -------
    dict
        condition -> trained model
    """
    if condition_col not in adata.obs.columns:
        raise ValueError(f"{condition_col} not found in adata.obs.")

    if conditions is None:
        conditions = adata.obs[condition_col].dropna().unique().tolist()

    trained_models = {}

    for condition in conditions:
        print(f"Training RegressionModel for {condition}...")

        adata_sub = adata[adata.obs[condition_col] == condition].copy()

        model = train_regression_model(
            adata_sub,
            labels_key=labels_key,
            batch_key=batch_key,
            layer=layer,
            max_epochs=max_epochs,
            batch_size=batch_size,
            lr=lr,
            accelerator=accelerator
        )

        trained_models[condition] = model

    return trained_models


def save_regression_models(
    trained_models,
    output_dir,
    suffix="_regression_model",
):
    """
    Save trained regression models to disk.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for condition, model in trained_models.items():
        model_dir = output_dir / f"{condition}{suffix}"
        model.save(
          str(model_dir),
          overwrite=True
)


def load_and_export_regression_posterior(
    adata,
    model_dir,
    num_samples=1000,
    batch_size=2500,
):
    """
    Load a saved RegressionModel and export posterior estimates into AnnData.

    Returns
    -------
    AnnData, RegressionModel
    """
    from cell2location.models import RegressionModel

    model = RegressionModel.load(
        model_dir,
        adata=adata
    )

    model.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": num_samples,
            "batch_size": batch_size
        }
    )

    return adata, model


def export_inferred_signatures(
    adata,
    covariate_col="celltype",
    model_name="model",
    output_folder="./",
    key=("mod", "post_sample_means", "per_cluster_mu_fg"),
):
    """
    Export inferred cell-type signatures from Cell2Location RegressionModel posterior.

    Parameters
    ----------
    adata : AnnData
        AnnData after model.export_posterior().
    covariate_col : str
        obs column with cell-type labels.
    model_name : str
        Prefix for output files.
    output_folder : str
        Output directory.
    key : tuple
        Nested key to posterior signatures.

    Returns
    -------
    inferred_df, average_df
    """
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    x = adata.uns
    for k in key:
        x = x[k]

    mu = x.T

    celltypes = adata.obs[covariate_col].astype(str).unique()

    inferred_df = pd.DataFrame(
        mu,
        index=adata.var_names,
        columns=celltypes
    )

    inferred_file = output_folder / f"{model_name}_inferred_signatures.csv"
    inferred_df.to_csv(inferred_file)

    average_df = (
        adata.to_df()
        .join(adata.obs[covariate_col])
        .groupby(covariate_col)
        .mean()
        .T
    )

    average_file = output_folder / f"{model_name}_cluster_average_signatures.csv"
    average_df.to_csv(average_file)

    print(f"Saved inferred signatures -> {inferred_file}")
    print(f"Saved cluster averages -> {average_file}")

    return inferred_df, average_df
