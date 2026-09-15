"""
Principal component analysis over the per-nucleus feature table.

Why the implementation is a bare SVD rather than scikit-learn: a PCA on a
standardised matrix IS the SVD of that matrix, and the repository's
dependency set is deliberately small -- adding scikit-learn for twenty
lines would oblige every collaborator to re-sync their environment. The
results are identical to `sklearn.decomposition.PCA(whiten=False)` up to
component sign, which is fixed here deterministically (see `feature_pca`).

What "all features" can and cannot mean
---------------------------------------
Two groups of columns are excluded by default, and neither is a matter of
taste:

* Position and bookkeeping (`position`, `label`, `mid_z`, the three
  centroid coordinates). Where a nucleus sat in the field is a property of
  the microscope stage, and letting it into the covariance invents a
  component that maps the field layout.
* Per-field background (`*_bg`). These are one value per field, shared by
  every nucleus in it, so they encode acquisition conditions rather than
  nuclei -- and in a multi-condition dataset they alias directly onto
  condition.

The remaining columns are then pruned for redundancy, which matters more
than it sounds. `n_voxels` and `volume_um3` are the same measurement in
different units (r = 1.0000);
each intensity channel contributes median, mean, p90, std and their
background-corrected twins at r > 0.98. Left alone, a feature measured
five times gets five votes, and PC1 becomes "how many ways did we write
down intensity" rather than a property of chromatin. `feature_matrix`
therefore keeps one representative per near-duplicate group, choosing by
the declared priority below rather than by chance.
"""

import numpy as np
import pandas as pd

EXCLUDE_BOOKKEEPING = ("position", "label", "mid_z",
                       "centroid_z_um", "centroid_y_um", "centroid_x_um")
"""Columns describing where the nucleus was, not what it was like."""

PRIORITY = (
    # geometry and chromatin texture first: these are the phenotype
    "volume_um3", "max_area_um2", "solidity", "mid_solidity",
    "dna_cv", "mid_dna_cv_corr", "dna_p90_p10", "dna_p90_p50",
    "dna_top10_fraction", "dna_radial_norm", "mid_dna_radial_norm",
    "dna_radial_um", "mid_dna_radial_um",
    # then intensity, preferring background-corrected and normalised forms
    "DAPI_enrichment", "Cy5_enrichment", "mCherry_enrichment",
    "DAPI_mean_corr", "Cy5_mean_corr", "mCherry_mean_corr",
    "DAPI_integrated_corr", "Cy5_integrated_corr", "mCherry_integrated_corr",
    "DAPI_std", "Cy5_std", "mCherry_std",
)
"""
Order in which competing near-duplicates are resolved.

A feature is kept if it correlates below the threshold with everything
already kept, so earlier entries win. Anything not listed is considered
after these, in table order -- the list biases the outcome towards
interpretable, background-corrected, unit-normalised quantities without
hard-coding the final feature set.
"""


def feature_matrix(table, exclude=(), redundancy=0.98, extra_features=()):
    """
    Numeric feature matrix with bookkeeping and near-duplicates removed.

    Returns (X, dropped) where X is a DataFrame indexed like `table` and
    `dropped` maps each removed column to the reason -- print it; it is the
    audit trail for what the PCA actually saw.

    `redundancy` is the absolute Pearson correlation above which two
    features are treated as the same measurement. Set it to 1.0 to keep
    everything except exact duplicates, or None to disable pruning.
    """
    num = table.select_dtypes(include=[np.number])
    dropped = {}

    for c in num.columns:
        if c in EXCLUDE_BOOKKEEPING and c not in extra_features:
            dropped[c] = "position/bookkeeping"
        elif c.endswith("_bg") and c not in extra_features:
            dropped[c] = "per-field background, not a nucleus property"
        elif c in exclude:
            dropped[c] = "excluded by caller"

    cand = [c for c in num.columns if c not in dropped]
    const = [c for c in cand if num[c].nunique() <= 1]
    for c in const:
        dropped[c] = "constant"
    cand = [c for c in cand if c not in const]

    if redundancy is not None and cand:
        order = ([c for c in PRIORITY if c in cand]
                 + [c for c in cand if c not in PRIORITY])
        corr = num[cand].corr().abs()
        keep = []
        for c in order:
            clash = [k for k in keep if corr.loc[c, k] >= redundancy]
            if clash:
                best = max(clash, key=lambda k: corr.loc[c, k])
                dropped[c] = f"r={corr.loc[c, best]:.4f} with {best}"
            else:
                keep.append(c)
        cand = [c for c in order if c in keep]

    X = num[cand].astype(np.float64)
    if X.isna().any().any():
        n0 = len(X)
        X = X.dropna()
        dropped["__rows__"] = f"{n0 - len(X)} rows with missing values"
    return X, dropped


def feature_pca(table, n_components=6, redundancy=0.98, exclude=(),
                exclude_z_border=True, exclude_xy_border=True):
    """
    Correlation-matrix PCA over the feature table.

    The same two truncation filters as the figure functions apply by
    default, and they are not cosmetic here: a laterally cut nucleus is a
    fragment whose area, volume and integrated intensity are all
    systematically low, so including them adds a truncation axis to the
    decomposition. On the example data the filters also change WHICH
    features survive redundancy pruning -- 797 unfiltered nuclei leave 28,
    the 241 interior ones leave 22 -- so a PCA fit on the unfiltered table
    is not the same model with more points.

    Features are z-scored before decomposition, which is the right choice
    here because the columns differ by orders of magnitude in units: on raw
    values `DAPI_integrated_corr` (up to 2e6) would supply essentially the
    whole variance and the components would be a unit conversion.

    Component signs are fixed so that each component's largest-magnitude
    loading is positive. SVD signs are otherwise arbitrary, and without
    this the same data can produce a mirrored figure on a different machine
    or library version.

    Returns a dict:
        scores      `table` with PC1..PCn appended (ready for plotting)
        loadings    features x components, the orthonormal basis
        explained   per-component variance ratio and its cumulative sum
        features    the columns used, in order
        dropped     column -> why it was not used
        center/scale  the z-scoring applied, for projecting new data
    """
    if exclude_z_border and "mid_at_z_border" in table.columns:
        table = table[~table["mid_at_z_border"]]
    if exclude_xy_border and "touches_xy_border" in table.columns:
        table = table[~table["touches_xy_border"]]

    X, dropped = feature_matrix(table, exclude=exclude, redundancy=redundancy)
    if X.shape[1] < 2:
        raise ValueError(f"need at least 2 features, got {list(X.columns)}")

    center = X.mean()
    scale = X.std(ddof=1).replace(0.0, 1.0)
    Z = ((X - center) / scale).to_numpy()

    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    k = min(n_components, Vt.shape[0])
    U, S, Vt = U[:, :k], S[:k], Vt[:k]

    for i in range(k):                      # deterministic component signs
        j = np.argmax(np.abs(Vt[i]))
        if Vt[i, j] < 0:
            Vt[i], U[:, i] = -Vt[i], -U[:, i]

    var = (S ** 2) / (Z.shape[0] - 1)
    total = ((Z ** 2).sum() / (Z.shape[0] - 1))
    names = [f"PC{i + 1}" for i in range(k)]

    scores = pd.DataFrame(U * S, index=X.index, columns=names)
    out = table.loc[X.index].copy()
    for n in names:
        out[n] = scores[n]

    return dict(
        scores=out,
        loadings=pd.DataFrame(Vt.T, index=X.columns, columns=names),
        explained=pd.DataFrame({"variance_ratio": var / total,
                                "cumulative": np.cumsum(var / total)},
                               index=names),
        features=list(X.columns), dropped=dropped,
        center=center, scale=scale, n_samples=int(Z.shape[0]))
