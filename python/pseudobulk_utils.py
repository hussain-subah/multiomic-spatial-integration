"""
Pseudobulk synthetic-mixture utilities for validating the SpaceJam/
Cell2Location deconvolution model against known ground-truth cell-type
composition.

Synthetic ROIs are built by summing raw snRNA-seq counts from a held-out
split of the reference (never the cells used to derive the signature
matrix, so the validation isn't circular), then given realistic technical
characteristics -- negative-probe counts, sequencing depth, batch identity
-- borrowed from a real GeoMx ROI. Biological signal always comes from the
held-out reference cells; the donor ROI only contributes noise/depth/batch.
"""

import numpy as np
import pandas as pd

from python.spacejam_pyro_utils import make_spot2sample_mat


def split_reference_cells(
    adata_ref,
    celltype_col="celltype",
    holdout_frac=0.3,
    min_holdout_cells=10,
    seed=0,
):
    """
    Stratified per-cell-type train/holdout split of a reference AnnData.

    Cell types too small to yield at least `min_holdout_cells` holdout
    cells at `holdout_frac` are kept entirely in the train split (they can
    still inform the signature matrix, just can't be represented in
    synthetic mixtures).

    Parameters
    ----------
    adata_ref : AnnData
        Reference snRNA-seq AnnData with a cell-type label column.
    celltype_col : str
        obs column with cell-type labels.
    holdout_frac : float
        Fraction of each cell type's cells to hold out.
    min_holdout_cells : int
        Minimum holdout cells required per cell type to include it in the
        holdout split at all.
    seed : int
        Random seed.

    Returns
    -------
    (AnnData, AnnData)
        (adata_train, adata_holdout)
    """
    rng = np.random.default_rng(seed)

    labels = adata_ref.obs[celltype_col].astype(str).values

    train_idx_parts = []
    holdout_idx_parts = []
    skipped = []

    for ct in pd.unique(labels):
        ct_idx = np.flatnonzero(labels == ct)
        n_holdout = int(round(len(ct_idx) * holdout_frac))

        if n_holdout < min_holdout_cells:
            skipped.append(ct)
            train_idx_parts.append(ct_idx)
            continue

        shuffled = rng.permutation(ct_idx)
        holdout_idx_parts.append(shuffled[:n_holdout])
        train_idx_parts.append(shuffled[n_holdout:])

    train_idx = np.concatenate(train_idx_parts)
    holdout_idx = (
        np.concatenate(holdout_idx_parts) if holdout_idx_parts else np.array([], dtype=int)
    )

    if skipped:
        print(
            f"split_reference_cells: {len(skipped)} cell type(s) had fewer than "
            f"{min_holdout_cells} holdout cells at holdout_frac={holdout_frac} and "
            f"were kept entirely in the train split (excluded from synthetic "
            f"mixtures): {skipped}"
        )

    adata_train = adata_ref[train_idx].copy()
    adata_holdout = adata_ref[holdout_idx].copy()

    return adata_train, adata_holdout


def sample_dirichlet_compositions(celltypes, n_mixtures, concentration=0.3, seed=0):
    """
    Draw ground-truth composition vectors from a Dirichlet distribution.

    A concentration well below 1 yields sparse/skewed mixtures (a few
    dominant cell types), more representative of real ROI composition
    than a uniform (concentration=1) or peaked (concentration>1) mixture.

    Parameters
    ----------
    celltypes : list
        Cell types to draw proportions over.
    n_mixtures : int
        Number of mixtures (rows) to draw.
    concentration : float
        Dirichlet concentration parameter, applied uniformly across
        cell types.
    seed : int
        Random seed.

    Returns
    -------
    pandas.DataFrame
        n_mixtures x len(celltypes), each row summing to 1.
    """
    rng = np.random.default_rng(seed)
    celltypes = list(celltypes)

    alpha = np.full(len(celltypes), concentration)
    draws = rng.dirichlet(alpha, size=n_mixtures)

    return pd.DataFrame(draws, columns=celltypes)


def build_synthetic_mixture(
    adata_holdout,
    celltype_col,
    target_proportions,
    n_cells_per_mixture=200,
    counts_layer="counts",
    seed=None,
):
    """
    Sample cells per a target composition and sum their counts into one
    pseudobulk profile.

    Parameters
    ----------
    adata_holdout : AnnData
        Held-out reference cells (see `split_reference_cells`).
    celltype_col : str
        obs column with cell-type labels.
    target_proportions : pandas.Series
        Target fraction per cell type (index = cell type, values sum to ~1).
    n_cells_per_mixture : int
        Total number of cells to sample into this synthetic ROI.
    counts_layer : str
        Layer holding raw counts.
    seed : int, optional
        Random seed for this mixture's sampling.

    Returns
    -------
    (numpy.ndarray, pandas.Series)
        (pseudobulk_vector over adata_holdout.var_names, realized_proportions
        actually sampled -- the true ground truth -- indexed like
        target_proportions).
    """
    rng = np.random.default_rng(seed)

    labels = adata_holdout.obs[celltype_col].astype(str).values
    counts = adata_holdout.layers[counts_layer]
    sparse = not isinstance(counts, np.ndarray)

    n_genes = adata_holdout.n_vars
    pseudobulk = np.zeros(n_genes, dtype=np.float64)
    realized_counts = pd.Series(0, index=target_proportions.index, dtype=int)

    target_n_cells = (target_proportions * n_cells_per_mixture).round().astype(int)

    for ct, n_needed in target_n_cells.items():
        if n_needed <= 0:
            continue

        ct_idx = np.flatnonzero(labels == ct)

        if len(ct_idx) == 0:
            continue  # cell type absent from holdout pool; can't be represented

        replace = len(ct_idx) < n_needed
        sampled_idx = rng.choice(ct_idx, size=n_needed, replace=replace)

        sampled_counts = counts[sampled_idx]
        if sparse:
            sampled_counts = sampled_counts.toarray()

        pseudobulk += np.asarray(sampled_counts).sum(axis=0)
        realized_counts[ct] = n_needed

    total_cells = int(realized_counts.sum())

    if total_cells == 0:
        raise ValueError(
            "No cells sampled for this mixture -- check that target_proportions' "
            "cell types are present in the holdout pool."
        )

    realized_proportions = realized_counts / total_cells

    return pseudobulk, realized_proportions


def generate_synthetic_dataset(
    adata_holdout,
    celltype_col,
    adata_real_donor,
    celltypes,
    n_mixtures=100,
    n_cells_per_mixture=200,
    concentration=0.3,
    counts_layer="counts",
    real_counts_layer="counts",
    neg_probes_key="negProbes",
    sample_id_col="Scan_ID",
    seed=0,
):
    """
    Build a full synthetic pseudobulk validation dataset.

    Biological signal (ground-truth composition, expression) comes only
    from `adata_holdout`. Technical characteristics -- negative-probe
    counts, total sequencing depth, and batch/experiment identity -- are
    borrowed from randomly chosen real ROIs in `adata_real_donor`, giving
    synthetic ROIs realistic noise without an invented noise model.

    `adata_holdout` and `adata_real_donor` must already share the same
    gene set in the same order (align both with
    `python.spacejam_pyro_utils.align_genes` against the same signature
    matrix before calling this).

    Parameters
    ----------
    adata_holdout : AnnData
        Held-out reference cells (see `split_reference_cells`).
    celltype_col : str
        obs column with cell-type labels in `adata_holdout`.
    adata_real_donor : AnnData
        Real GeoMx WTA AnnData (or a condition subset) donating technical
        characteristics only.
    celltypes : list
        Full cell-type list matching the *signature matrix*'s column
        order (i.e. `cell_state_mat`'s factor order) -- not just the cell
        types present in `adata_holdout`. Cell types absent from the
        signature's train split were never split into holdout at all
        (see `split_reference_cells`), and the model's output has one
        column per signature cell type regardless, so ground truth must
        span the same full list (with 0 proportion for any cell type
        absent from the holdout pool) or column counts won't match the
        model's posterior.
    n_mixtures : int
        Number of synthetic ROIs to generate.
    n_cells_per_mixture : int
        Cells summed per synthetic ROI.
    concentration : float
        Dirichlet concentration for ground-truth compositions.
    seed : int
        Random seed.

    Returns
    -------
    dict
        X_data, Y_data, spot2sample_mat, experiment_names, ground_truth_df,
        donor_roi_ids, celltypes -- ready for
        `python.spacejam_pyro_utils.train_spacejam_model()` alongside a
        gene-aligned `cell_state_mat`.
    """
    if list(adata_holdout.var_names) != list(adata_real_donor.var_names):
        raise ValueError(
            "adata_holdout and adata_real_donor must share the same gene set "
            "in the same order -- align both against the signature matrix "
            "with align_genes() before calling generate_synthetic_dataset()."
        )

    rng = np.random.default_rng(seed)

    celltypes = list(celltypes)

    compositions = sample_dirichlet_compositions(
        celltypes, n_mixtures, concentration=concentration, seed=seed
    )

    n_donors = adata_real_donor.n_obs
    donor_positions = rng.integers(0, n_donors, size=n_mixtures)

    neg_probes = adata_real_donor.obsm[neg_probes_key]
    if not isinstance(neg_probes, np.ndarray):
        neg_probes = np.asarray(neg_probes)

    real_counts = adata_real_donor.layers[real_counts_layer]
    real_sparse = not isinstance(real_counts, np.ndarray)

    X_rows = []
    Y_rows = []
    donor_roi_ids = []
    ground_truth_rows = []

    for i in range(n_mixtures):
        target_proportions = compositions.iloc[i]

        pseudobulk, realized_proportions = build_synthetic_mixture(
            adata_holdout,
            celltype_col=celltype_col,
            target_proportions=target_proportions,
            n_cells_per_mixture=n_cells_per_mixture,
            counts_layer=counts_layer,
            seed=int(rng.integers(0, 2**31 - 1)),
        )

        donor_pos = donor_positions[i]

        donor_row = real_counts[donor_pos]
        if real_sparse:
            donor_row = donor_row.toarray()
        target_depth = int(np.asarray(donor_row).sum())

        pseudobulk_probs = pseudobulk / pseudobulk.sum()
        resampled = rng.multinomial(target_depth, pseudobulk_probs)

        X_rows.append(resampled)
        Y_rows.append(np.asarray(neg_probes[donor_pos]))
        donor_roi_ids.append(adata_real_donor.obs_names[donor_pos])
        ground_truth_rows.append(realized_proportions)

    X_data = np.vstack(X_rows).astype(np.float32)
    Y_data = np.vstack(Y_rows).astype(np.float32)
    ground_truth_df = pd.DataFrame(ground_truth_rows).reset_index(drop=True)

    donor_sample_ids = adata_real_donor.obs[sample_id_col].values[donor_positions]
    spot2sample_mat, experiment_names = make_spot2sample_mat(donor_sample_ids)

    return {
        "X_data": X_data,
        "Y_data": Y_data,
        "spot2sample_mat": spot2sample_mat,
        "experiment_names": experiment_names,
        "ground_truth_df": ground_truth_df,
        "donor_roi_ids": donor_roi_ids,
        "celltypes": celltypes,
    }
