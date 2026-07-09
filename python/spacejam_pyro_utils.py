"""
Utilities for running SpaceJam / Cell2Location-style Pyro models.
"""

import torch
import numpy as np
import pandas as pd

def make_spot2sample_mat(sample_ids):
    """
    sample_ids: array-like of length n_rois, e.g. Scan_ID.
    returns: n_rois x n_exper one-hot matrix and experiment/sample names.
    """
    sample_cat = pd.Categorical(sample_ids)
    n_rois = len(sample_cat)
    n_exper = len(sample_cat.categories)

    mat = np.zeros((n_rois, n_exper), dtype=np.float32)
    mat[np.arange(n_rois), sample_cat.codes] = 1.0

    return mat, sample_cat.categories


def prepare_spacejam_inputs(
    adata,
    signatures,
    neg_probes_key="negProbes",
    counts_layer="counts",
    nuclei_col="Nuclei"
):
    """
    Align AnnData and signature matrix for spatial modeling.
    """

    counts = adata.layers[counts_layer]
    if not isinstance(counts, np.ndarray):
        counts = counts.toarray()

    neg = adata.obsm[neg_probes_key]

    genes_shared = adata.var_names.intersection(signatures.index)

    counts = counts[:, adata.var_names.isin(genes_shared)]
    signatures = signatures.loc[genes_shared]

    return {
        "counts": torch.tensor(counts, dtype=torch.float32),
        "signatures": torch.tensor(signatures.values, dtype=torch.float32),
        "neg_probes": torch.tensor(neg, dtype=torch.float32),
        "nuclei": torch.tensor(adata.obs[nuclei_col].values, dtype=torch.float32)
    }


def train_spacejam_model(
    inputs,
    model_class,
    n_steps=20000,
    learning_rate=0.005,
    device="cpu"
):
    """
    Train Pyro spatial model.
    """

    model = model_class(**inputs)

    model.train_model(
        n_steps=n_steps,
        lr=learning_rate,
        device=device
    )

    return model


def extract_posterior_abundance(
    model,
    adata,
    normalize=True
):
    """
    Extract cell-type abundance from trained model.
    """

    abundance = model.get_cell_abundance()

    df = pd.DataFrame(
        abundance,
        index=adata.obs_names,
        columns=model.cell_types
    )

    if normalize:
        df = df.div(df.sum(axis=1), axis=0)

    return df

def align_genes(signature_df, adata):
    """
    Align a genes x cell-type signature matrix with an AnnData object.

    Parameters
    ----------
    signature_df : pandas.DataFrame
        Rows = genes, columns = cell types/signatures.
    adata : anndata.AnnData
        AnnData object with genes in adata.var_names.

    Returns
    -------
    signature_df_aligned, adata_aligned
    """
    shared_genes = signature_df.index.intersection(adata.var_names)

    if len(shared_genes) == 0:
        raise ValueError("No shared genes between signature_df and adata.var_names.")

    signature_df_aligned = signature_df.loc[shared_genes].copy()
    adata_aligned = adata[:, shared_genes].copy()

    return signature_df_aligned, adata_aligned
