"""Load the deposited dataset.

Two routes to the same object:

``load_dataset(source="standard")`` (the default) reads the plain formats —
MatrixMarket for the expression matrix, gzipped CSV for the metadata and
embeddings, JSON for the palettes. Nothing here needs ``anndata`` to *read*; the
pieces are assembled into an AnnData only because the plotting code wants one.

``load_dataset(source="h5ad")`` reads ``figure_dataset.h5ad`` instead. Both give
the same values.

The ``logcounts_scaled`` layer that the expression panels are drawn from is not
shipped as a matrix. It is ``logcounts`` multiplied by a per-(dataset, gene)
factor, and reconstructing it that way is bit-identical in float32 while saving
60 MB. :func:`scaled_layer` does the reconstruction.

The 152 genes that were filtered out of the deposited object before it was saved
ride alongside in their own small matrix; :func:`load_extra_genes` reads it.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from anndata import AnnData

from .constants import CLUSTER_COLORS

REPO = Path(__file__).resolve().parent.parent.parent
DATA = REPO / "data"

# obs columns that are categorical, and must not be re-inferred as numbers:
# nr_cluster and cluster are blank for the in vitro cells, so pandas would
# otherwise read them as floats and turn "0" into 0.0.
CATEGORICAL_OBS = (
    "orig_ident",
    "dataset",
    "dataset_type",
    "lc_stage",
    "nr_cluster",
    "cluster",
    "cc_phase",
    "transferred_cc_phase",
    "cell_cycle_group",
)


def read_manifest(datadir: Path = DATA) -> dict:
    return json.loads((datadir / "MANIFEST.json").read_text())


def _apply_categories(obs: pd.DataFrame, manifest: dict) -> pd.DataFrame:
    """Restore the recorded category order, so plots group in the published order."""
    for column, spec in manifest["categoricals"]["obs"].items():
        if column in obs.columns:
            obs[column] = pd.Categorical(
                obs[column], categories=spec["categories"], ordered=spec["ordered"]
            )
    return obs


def _read_matrix(path: Path, n_obs: int, n_vars: int) -> sp.csr_matrix:
    """Read the gzipped MatrixMarket matrix as cells x genes, float32.

    MatrixMarket is decimal text, so scipy hands back float64; the file holds the
    shortest round-trip form of the original float32, and narrowing recovers the
    original bits exactly.
    """
    with gzip.open(path, "rb") as handle:
        matrix = scipy.io.mmread(handle)
    matrix = sp.csr_matrix(matrix)
    if matrix.shape != (n_obs, n_vars):
        raise ValueError(f"{path.name}: expected {(n_obs, n_vars)}, got {matrix.shape}")
    matrix.data = matrix.data.astype(np.float32)
    return matrix


def scaled_layer(adata: AnnData, factors: pd.DataFrame) -> sp.csr_matrix:
    """Rebuild ``logcounts_scaled`` from ``logcounts`` and the per-gene factors.

    The factor depends on the dataset a cell came from as well as the gene — a
    single per-gene factor does not work. Genes with no expression in a dataset
    have no factor; those entries carry no stored value anyway.

    The result is float64, matching the layer in the original analysis object;
    the multiplication reproduces it to within one unit in the last place.
    """
    scaled = adata.layers["logcounts"].tocsr(copy=True).astype(np.float64)
    datasets = pd.Categorical(
        adata.obs["dataset"].astype(str), categories=list(factors.columns)
    ).codes

    # One factor per stored entry: the row's dataset crossed with the entry's gene.
    row_of_entry = np.repeat(np.arange(scaled.shape[0]), np.diff(scaled.indptr))
    table = factors.to_numpy(dtype=np.float64)
    per_entry = table[scaled.indices, datasets[row_of_entry]]
    scaled.data *= np.nan_to_num(per_entry)
    return scaled


def _read_factors(datadir: Path, var_names: pd.Index) -> pd.DataFrame:
    """Long-form factor table -> genes x datasets."""
    long = pd.read_csv(datadir / "logcounts_scaled_factors.csv.gz")
    wide = long.pivot(index="gene_id", columns="dataset", values="factor")
    return wide.reindex(var_names)


def _load_standard(datadir: Path, manifest: dict) -> AnnData:
    n_obs, n_vars = manifest["n_obs"], manifest["n_vars"]

    obs = pd.read_csv(
        datadir / "obs.csv.gz",
        index_col="barcode",
        dtype={column: "str" for column in CATEGORICAL_OBS},
    )
    obs = _apply_categories(obs, manifest)
    var = pd.read_csv(datadir / "var.csv.gz", index_col="gene_id", dtype="str").fillna("")

    adata = AnnData(obs=obs, var=var)
    adata.layers["logcounts"] = _read_matrix(datadir / "logcounts.mtx.gz", n_obs, n_vars)

    for path in sorted(datadir.glob("obsm_*.csv.gz")):
        key = path.name[len("obsm_") : -len(".csv.gz")]
        frame = pd.read_csv(path, index_col=0)
        if not frame.index.equals(adata.obs_names):
            raise ValueError(f"{path.name}: cell order does not match obs.csv.gz")
        dtype = np.float32 if key.startswith("3d_umap") else np.float64
        adata.obsm[key] = frame.to_numpy(dtype=dtype)

    colors = json.loads((datadir / "uns_colors.json").read_text())
    adata.uns.update(colors)

    factors = _read_factors(datadir, adata.var_names)
    adata.varm["logcounts_scaled_factors"] = factors.to_numpy()
    adata.uns["logcounts_scaled_factors_datasets"] = list(factors.columns)
    adata.layers["logcounts_scaled"] = scaled_layer(adata, factors)
    return adata


def _load_h5ad(datadir: Path) -> AnnData:
    import anndata

    adata = anndata.read_h5ad(datadir / "figure_dataset.h5ad")
    factors = pd.DataFrame(
        adata.varm["logcounts_scaled_factors"],
        index=adata.var_names,
        columns=list(adata.uns["logcounts_scaled_factors_datasets"]),
    )
    adata.layers["logcounts_scaled"] = scaled_layer(adata, factors)
    return adata


def _apply_cluster_palette(adata: AnnData) -> AnnData:
    """Put the published cluster palette where scanpy looks for it.

    Neither the deposited object nor the standard files carry
    ``*_cluster_colors``, so ``sc.pl.heatmap`` would colour its groupby bar from
    scanpy's default palette while every other panel passes ``CLUSTER_COLORS``
    explicitly. Setting it here keeps the two agreeing.
    """
    for column in ("nr_cluster", "cluster"):
        values = adata.obs.get(column)
        if values is None or not hasattr(values, "cat"):
            continue
        adata.uns[f"{column}_colors"] = list(CLUSTER_COLORS[: len(values.cat.categories)])
    return adata


def load_dataset(datadir: Path = DATA, source: str = "standard") -> AnnData:
    """Load the deposited dataset, with ``logcounts_scaled`` reconstructed."""
    datadir = Path(datadir)
    manifest = read_manifest(datadir)
    if source == "standard":
        return _apply_cluster_palette(_load_standard(datadir, manifest))
    if source == "h5ad":
        return _apply_cluster_palette(_load_h5ad(datadir))
    raise ValueError(f"source must be 'standard' or 'h5ad', not {source!r}")


def load_extra_genes(datadir: Path = DATA) -> tuple[sp.csr_matrix, pd.DataFrame]:
    """The 152 genes on unplaced contigs that the deposited object drops.

    Same cells in the same order as ``obs.csv.gz``, and normalised the same way
    as ``logcounts``, so the two matrices sit side by side as the 8322-gene
    universe the published differential expression ran on. Only the
    Supplementary 5 volcano panels need them; :mod:`bzfig.de` does the widening.
    """
    datadir = Path(datadir)
    manifest = read_manifest(datadir)
    var = pd.read_csv(datadir / "var_extra.csv.gz", index_col="gene_id", dtype="str").fillna("")
    n_vars = manifest["extra_genes"]["n_vars"]
    if len(var) != n_vars:
        raise ValueError(f"var_extra.csv.gz: expected {n_vars} genes, got {len(var)}")
    matrix = _read_matrix(datadir / "logcounts_extra.mtx.gz", manifest["n_obs"], n_vars)
    return matrix, var


def load_supp1a_markers(datadir: Path = DATA) -> pd.Index:
    """The Supplementary 1A gene list, in the order the panel plots them.

    The shipped file holds a per-cluster score for each of 600 genes and is
    already stored grouped by cluster, so its row order *is* the plotting order.
    Read it off the file rather than re-deriving it from the scores.
    """
    scores = pd.read_csv(Path(datadir) / "supp1a_marker_genes.csv", header=None, index_col=0)
    return scores.index
