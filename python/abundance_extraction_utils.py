"""
Utilities for extracting and formatting spatial cell-type abundances
from trained Pyro/SpaceJam models.
"""

from pathlib import Path
import torch
import pandas as pd
import numpy as np


def load_variational_means(file):
    """
    Load saved variational posterior means from disk.

    Parameters
    ----------
    file : str
        Path to .pt file.

    Returns
    -------
    dict
        Dictionary of posterior mean tensors.
    """
    return torch.load(file, map_location="cpu")


def extract_spot_factors(
    variational_means,
    key="AutoNormal.locs.spot_factors"
):
    """
    Extract spot_factors from variational posterior means.

    Parameters
    ----------
    variational_means : dict
        Loaded posterior mean dictionary.
    key : str
        Key containing spot factor posterior means.

    Returns
    -------
    numpy.ndarray
        ROI x factor abundance matrix.
    """
    if key not in variational_means:
        available = list(variational_means.keys())
        raise KeyError(
            f"{key} not found in variational_means. "
            f"Available keys include: {available[:10]}"
        )

    spot_factors = variational_means[key]

    if isinstance(spot_factors, torch.Tensor):
        spot_factors = spot_factors.detach().cpu().numpy()

    return np.asarray(spot_factors)


def normalize_spot_factors(spot_factors, eps=1e-12):
    """
    Convert absolute spot factors to relative proportions.

    Parameters
    ----------
    spot_factors : array-like
        ROI x factor matrix.
    eps : float
        Small value to avoid division by zero.

    Returns
    -------
    numpy.ndarray
        Row-normalized abundance matrix.
    """
    spot_factors = np.asarray(spot_factors)
    row_sums = spot_factors.sum(axis=1, keepdims=True)

    return spot_factors / np.maximum(row_sums, eps)


def make_abundance_dataframe(
    spot_factors,
    obs,
    celltype_labels,
    roi_id_col=None
):
    """
    Convert spot factor matrix into metadata-annotated dataframe.

    Parameters
    ----------
    spot_factors : array-like
        ROI x celltype matrix.
    obs : pandas.DataFrame
        ROI metadata.
    celltype_labels : list
        Cell-type names matching columns of spot_factors.
    roi_id_col : str, optional
        Optional ROI ID column to preserve.

    Returns
    -------
    pandas.DataFrame
    """
    spot_factors = np.asarray(spot_factors)

    if spot_factors.shape[1] != len(celltype_labels):
        raise ValueError(
            "Number of spot factor columns does not match number of celltype labels."
        )

    if spot_factors.shape[0] != obs.shape[0]:
        raise ValueError(
            "Number of spot factor rows does not match number of metadata rows."
        )

    abundance_df = pd.DataFrame(
        spot_factors,
        index=obs.index,
        columns=celltype_labels
    )

    meta = obs.copy()

    if roi_id_col is None:
        meta = meta.copy()

        # AnnData obs index is the true ROI identifier.
        # Avoid duplicate ROI_ID columns from GeoMx metadata.
        if "ROI_ID" in meta.columns:
            meta = meta.rename(columns={"ROI_ID": "GeoMx_ROI_ID"})

        meta = meta.reset_index()
        meta = meta.rename(columns={meta.columns[0]: "ROI_ID"})

    elif roi_id_col in meta.columns:
        meta = meta.copy()
        meta["ROI_ID"] = meta[roi_id_col].astype(str)

    else:
       raise ValueError(f"{roi_id_col} not found in obs.")

    abundance_df = abundance_df.reset_index(drop=True)

    return pd.concat([meta.reset_index(drop=True), abundance_df], axis=1)


def make_long_abundance_table(
    abundance_df,
    celltype_labels,
    abundance_col="rel_abundance",
    id_vars=None
):
    """
    Convert wide abundance dataframe to long format.

    Parameters
    ----------
    abundance_df : pandas.DataFrame
        Metadata + cell-type abundance columns.
    celltype_labels : list
        Cell-type abundance columns.
    abundance_col : str
        Name for abundance value column.
    id_vars : list, optional
        Metadata columns to preserve. If None, uses all non-celltype columns.

    Returns
    -------
    pandas.DataFrame
    """
    if id_vars is None:
        id_vars = [c for c in abundance_df.columns if c not in celltype_labels]

    long_df = abundance_df.melt(
        id_vars=id_vars,
        value_vars=celltype_labels,
        var_name="celltype",
        value_name=abundance_col
    )

    return long_df


def save_abundance_outputs(
    output_dir,
    ad_abs_df=None,
    ad_rel_df=None,
    ctrl_abs_df=None,
    ctrl_rel_df=None,
    long_df=None
):
    """
    Save abundance outputs to CSV files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "ADCAA_spot_factors_abs.csv": ad_abs_df,
        "ADCAA_spot_factors_rel.csv": ad_rel_df,
        "CTRL_spot_factors_abs.csv": ctrl_abs_df,
        "CTRL_spot_factors_rel.csv": ctrl_rel_df,
        "roi_celltype_abundance_long.csv": long_df,
    }

    for filename, df in outputs.items():
        if df is not None:
            df.to_csv(output_dir / filename, index=False)
