"""
Evaluation utilities for pseudobulk deconvolution validation: comparing
model-inferred cell-type proportions against known synthetic ground truth
(see python/pseudobulk_utils.py for how the synthetic mixtures are built).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import seaborn as sns


def _recovery_row(label, true_vals, pred_vals):
    n = len(true_vals)

    if n >= 2 and np.std(true_vals) > 0 and np.std(pred_vals) > 0:
        pearson_r, _ = pearsonr(true_vals, pred_vals)
        spearman_rho, _ = spearmanr(true_vals, pred_vals)
    else:
        pearson_r = np.nan
        spearman_rho = np.nan

    rmse = float(np.sqrt(np.mean((true_vals - pred_vals) ** 2)))
    mae = float(np.mean(np.abs(true_vals - pred_vals)))

    return {
        "celltype": label,
        "pearson_r": pearson_r,
        "spearman_rho": spearman_rho,
        "rmse": rmse,
        "mae": mae,
        "n": n,
    }


def evaluate_recovery(ground_truth_df, inferred_df):
    """
    Compare inferred cell-type proportions against known ground truth.

    Parameters
    ----------
    ground_truth_df : pandas.DataFrame
        n_mixtures x celltypes, realized ground-truth proportions.
    inferred_df : pandas.DataFrame
        n_mixtures x celltypes, model-inferred proportions (same row order
        as ground_truth_df).

    Returns
    -------
    pandas.DataFrame
        One row per shared cell type, plus one "overall" row pooling all
        cell-type x mixture pairs: pearson_r, spearman_rho, rmse, mae, n.
    """
    shared_celltypes = [c for c in ground_truth_df.columns if c in inferred_df.columns]

    if len(shared_celltypes) == 0:
        raise ValueError("No shared cell-type columns between ground_truth_df and inferred_df.")

    if ground_truth_df.shape[0] != inferred_df.shape[0]:
        raise ValueError(
            "ground_truth_df and inferred_df must have the same number of "
            "rows (mixtures)."
        )

    rows = [
        _recovery_row(ct, ground_truth_df[ct].to_numpy(), inferred_df[ct].to_numpy())
        for ct in shared_celltypes
    ]

    true_all = ground_truth_df[shared_celltypes].to_numpy().ravel()
    pred_all = inferred_df[shared_celltypes].to_numpy().ravel()
    rows.append(_recovery_row("overall", true_all, pred_all))

    return pd.DataFrame(rows)


def plot_recovery_scatter(ground_truth_df, inferred_df, save_file=None):
    """
    Faceted true-vs-inferred scatter plot, one panel per cell type, with a
    y=x reference line.
    """
    shared_celltypes = [c for c in ground_truth_df.columns if c in inferred_df.columns]

    long_df = pd.concat(
        [
            pd.DataFrame({
                "celltype": ct,
                "true_proportion": ground_truth_df[ct].to_numpy(),
                "inferred_proportion": inferred_df[ct].to_numpy(),
            })
            for ct in shared_celltypes
        ],
        ignore_index=True,
    )

    g = sns.FacetGrid(long_df, col="celltype", col_wrap=6, height=2.5, sharex=True, sharey=True)
    g.map_dataframe(sns.scatterplot, x="true_proportion", y="inferred_proportion", alpha=0.5, s=15)

    axis_max = max(long_df["true_proportion"].max(), long_df["inferred_proportion"].max()) * 1.05

    for ax in g.axes.flat:
        ax.plot([0, axis_max], [0, axis_max], linestyle="--", color="gray", linewidth=1)
        ax.set_xlim(0, axis_max)
        ax.set_ylim(0, axis_max)

    g.set_axis_labels("True proportion", "Inferred proportion")
    g.fig.suptitle("Pseudobulk recovery: true vs. inferred proportion", y=1.02)

    if save_file is not None:
        Path(save_file).parent.mkdir(parents=True, exist_ok=True)
        g.savefig(save_file, dpi=300, bbox_inches="tight")

    return g


def plot_recovery_bias_by_abundance(ground_truth_df, inferred_df, save_file=None):
    """
    Mean(inferred - true) vs. mean(true) per cell type, to surface
    systematic under/over-estimation of rare vs. common cell types.
    """
    shared_celltypes = [c for c in ground_truth_df.columns if c in inferred_df.columns]

    bias_df = pd.DataFrame([
        {
            "celltype": ct,
            "mean_true": np.mean(ground_truth_df[ct].to_numpy()),
            "mean_bias": np.mean(inferred_df[ct].to_numpy() - ground_truth_df[ct].to_numpy()),
        }
        for ct in shared_celltypes
    ])

    plt.figure(figsize=(7, 5))
    sns.scatterplot(data=bias_df, x="mean_true", y="mean_bias")
    plt.axhline(0, linestyle="--", color="gray", linewidth=1)
    plt.xlabel("Mean true proportion")
    plt.ylabel("Mean bias (inferred - true)")
    plt.title("Recovery bias by cell-type abundance")
    plt.tight_layout()

    if save_file is not None:
        Path(save_file).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_file, dpi=300, bbox_inches="tight")

    return bias_df


def plot_recovery_bar(metrics_df, metric="pearson_r", save_file=None):
    """
    Sorted bar chart of a per-cell-type recovery metric.

    Complements the scatter/bias plots with a single at-a-glance ranking of
    which cell types the model recovers well vs. poorly.

    Parameters
    ----------
    metrics_df : pandas.DataFrame
        Output of `evaluate_recovery()` (one row per cell type plus an
        "overall" row, which is dropped here).
    metric : str
        Column to plot ("pearson_r", "spearman_rho", "rmse", or "mae").
    save_file : str, optional
        If given, save the figure there.
    """
    df = metrics_df[metrics_df["celltype"] != "overall"].copy()
    df = df.sort_values(metric, ascending=True)

    # For correlation metrics, color by sign (diverging); for error metrics,
    # a single sequential hue.
    is_corr = metric in ("pearson_r", "spearman_rho")
    if is_corr:
        colors = ["#2166AC" if v < 0 else "#009E73" for v in df[metric]]
    else:
        colors = "#4C72B0"

    plt.figure(figsize=(7, max(4, 0.25 * len(df))))
    plt.barh(df["celltype"], df[metric], color=colors)
    if is_corr:
        plt.axvline(0, linestyle="--", color="gray", linewidth=1)
    plt.xlabel(metric)
    plt.ylabel(None)
    plt.title(f"Pseudobulk recovery by cell type ({metric})")
    plt.tight_layout()

    if save_file is not None:
        Path(save_file).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_file, dpi=300, bbox_inches="tight")

    return df
