#!/usr/bin/env python3
"""Export a lightweight, format-independent data package for the figure repo.

The upstream analysis object
(``integrated_adata_me49_nr_subset_annotated_sparse.h5ad``, 8057 cells x
8170 genes) lives in the *methods* repository and is an ``.h5ad`` file, i.e.
it can only be read comfortably from Python/anndata.  This script distils it
down to the subset of content the figure panels actually need and writes that
subset twice:

1. **Standard, h5ad-independent formats** -- MatrixMarket for the expression
   matrix, gzipped CSV for the cell/gene metadata, the embeddings and the
   scaling factors, and JSON for the colour palettes.  Anything that can read
   MatrixMarket and CSV (R/Matrix, Julia/MatrixMarket.jl, Python/scipy, ...)
   can consume these.
2. **A stripped ``.h5ad``** with exactly the same content, as a convenience
   for scanpy users.

It additionally folds the CCC/MCC cell-cycle assignment (derived from
``nr_pseudotime_BCC_UCC.csv``) into the cell metadata, so downstream code never
has to touch that CSV again.

Every output is listed in ``MANIFEST.json`` with its sha256, size, description
and an (initially empty) ``urls`` array.  Running the script twice on the same
inputs produces byte-identical outputs (gzip streams are written with
``mtime=0``).

Usage
-----
    python scripts/export_dataset.py \
        --h5ad       /path/to/integrated_adata_me49_nr_subset_annotated_sparse.h5ad \
        --pseudotime /path/to/nr_pseudotime_BCC_UCC.csv \
        --outdir     data

    # optional extras, off by default because no panel needs them
    python scripts/export_dataset.py --include-scaled-matrix --include-counts

Package layout
--------------
``logcounts.mtx.gz``
    MatrixMarket coordinate file, **rows = cells, columns = genes** (the
    transpose of the ``.mtx`` CellRanger writes).  Row *i* is row *i* of
    ``obs.csv.gz``, column *j* is row *j* of ``var.csv.gz``.  Values are the
    shortest round-trip decimal form of the source float32.
``logcounts_scaled_factors.csv.gz``
    Long-form ``dataset,gene_id,factor`` table (float64) with
    ``logcounts_scaled[i, j] = logcounts[i, j] * factor[obs.dataset[i], j]``.
``logcounts_extra.mtx.gz`` / ``var_extra.csv.gz``
    The same, for the 152 ToxoDB-65 genes the analysis object was subset away
    from before it was saved -- all of them on unplaced ``KE*`` contigs, and
    among them the apicoplast and mitochondrial transcripts the Supplementary 5
    volcanoes plot.  Same cells in the same order, same normalisation, so the
    two matrices side by side are the 8322-gene universe the published
    differential expression ran on.
``obs.csv.gz`` / ``var.csv.gz``
    Cell / gene metadata; the first column is the index (``barcode`` /
    ``gene_id``).
``obsm_<key>.csv.gz``
    One embedding per file, indexed by barcode, columns ``<key>_1..n``.
``uns_colors.json``
    All ``*_colors`` palettes found in ``adata.uns``, as ``{key: [hex, ...]}``.
``supp1a_marker_genes.csv``
    Verbatim copy of the Supplementary Figure 1A gene list.
``figure_dataset.h5ad``
    The same content as an AnnData object.  ``X`` is deliberately ``None`` --
    the matrix lives in ``.layers`` only, so nothing is stored twice.  The
    scaling factors are in ``.varm['logcounts_scaled_factors']`` (n_vars x
    n_datasets, column order in ``.uns['logcounts_scaled_factors_datasets']``),
    so the file is self-sufficient.
``MANIFEST.json``
    File inventory (filename, sha256, bytes, description, ``urls``) plus the
    categorical level orders and column dtypes, which plain CSV cannot carry.

NOTE ON REDUNDANCY (all numbers measured by this script)
--------------------------------------------------------
* ``logcounts`` is **not** recoverable from ``counts``: the two layers do not
  even share a sparsity pattern (counts has 5,695,710 true non-zeros plus
  115,311 explicitly-stored zeros, logcounts has 4,951,318 non-zeros, and
  115,311 logcounts entries sit exactly on those stored zeros).  ``logcounts``
  was derived upstream from a separate ``raw_counts`` layer (non-integer,
  Seurat-corrected) as ``log1p(1e4 * raw_counts / rowsum(raw_counts))``, and
  ``raw_counts`` is not part of this package.  ``logcounts`` therefore has to
  be shipped as a matrix.
* ``logcounts_scaled`` **is** recoverable: it equals ``logcounts`` times a
  constant that depends only on (``dataset``, gene) -- a gene-wise standard
  deviation reciprocal computed per dataset upstream (for the ``me49``
  dataset the factor matches ``1 / sd(logcounts, ddof=1)`` to float32
  precision; for ``Nonreactivated`` it does not, because the upstream sd was
  computed before that dataset was subset).  Rebuilding it is bit-identical in
  float32 and accurate to 1 ulp (4.4e-16 relative, 1.4e-14 absolute) in
  float64, so the 60 MB matrix is replaced by a 2 x 8170 factor table.  Pass
  ``--include-scaled-matrix`` to ship the matrix anyway.
* ``counts`` and ``raw_counts`` are used by no panel and are not shipped;
  ``--include-counts`` adds ``counts``, and the manifest records how big that
  would be.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread, mmwrite

# --------------------------------------------------------------------------
# Configuration: what goes into the package
# --------------------------------------------------------------------------

#: Expression layers exported by default, with the on-disk dtype to use.
#: ``logcounts_scaled`` drives Fig 1C/1E/1F/1G, Supp 1A and Supp 5D;
#: ``logcounts`` drives Supp 4.
LAYERS: dict[str, str] = {
    "logcounts": "float32",
}

#: Layers that are measured but not shipped, each behind an opt-in flag.
#:
#: ``logcounts_scaled`` (the layer Fig 1C/1E/1F/1G, Supp 1A and Supp 5D use) is
#: reconstructed from ``logcounts`` and ``logcounts_scaled_factors.csv.gz``
#: instead of being shipped -- the reconstruction is bit-exact in float32 and
#: accurate to 1 ulp in float64, and it saves ~60 MB.
#:
#: ``counts`` holds integer UMI counts stored as float64 upstream; writing it
#: as int32 is lossless and roughly halves the MatrixMarket file.  No panel
#: uses it.  ``raw_counts`` is not offered at all -- see the module docstring.
OPTIONAL_LAYERS: dict[str, str] = {
    "logcounts_scaled": "float64",
    "counts": "int32",
}

#: Flag that turns each optional layer on.
OPTIONAL_LAYER_FLAGS = {
    "logcounts_scaled": "--include-scaled-matrix",
    "counts": "--include-counts",
}

#: obs columns to export, in this order.  ``dataset``/``dataset_type``/
#: ``lc_stage`` are redundant with ``orig_ident`` but are kept so the published
#: ``.query(...)`` strings run verbatim.  ``nr_cluster_learned`` is dropped (no
#: panel uses it) and so are the ``2d_projection_*`` / ``3d_umap_harmony_*``
#: obs columns, which duplicate the corresponding obsm entries (the export
#: asserts the duplication before dropping them).
KEEP_OBS = [
    "orig_ident",
    "dataset",
    "dataset_type",
    "lc_stage",
    "nr_cluster",
    "cluster",
    "cc_phase",
    "transferred_cc_phase",
]

#: obs columns that are dropped because an obsm entry carries the same values.
#: Mapping: obs column -> (obsm key, column position).
OBS_EMBEDDING_DUPLICATES = {
    "2d_projection_1": ("2d_projection", 0),
    "2d_projection_2": ("2d_projection", 1),
    "3d_umap_harmony_integration_1": ("3d_umap_harmony_integration", 0),
    "3d_umap_harmony_integration_2": ("3d_umap_harmony_integration", 1),
    "3d_umap_harmony_integration_3": ("3d_umap_harmony_integration", 2),
}

#: var columns to export (on top of the index, the gene id).  ``seqid`` is kept
#: because it documents which contigs are present: the deposited object holds
#: the 14 nuclear chromosomes only, and the 152 genes on the unplaced ``KE*``
#: contigs are shipped beside it -- see ``build_extra_genes``.
KEEP_VAR = ["seqid", "gene_name", "gene_description"]

#: obsm keys to export.  ``3d_umap_harmony_integration`` is required for the
#: bit-exact 3-D renders; ``2d_projection`` is derivable from it but cheap.
KEEP_OBSM = ["3d_umap_harmony_integration", "2d_projection"]

#: Extra small files copied verbatim into the package:
#: destination name -> (path relative to the methods repo, description).
EXTRA_FILES = {
    "supp1a_marker_genes.csv": (
        "notebooks/figures/heatmap_data.csv",
        "Supplementary Figure 1A marker gene list, copied verbatim: 600 rows, "
        "no header; column 1 = gene id, columns 2-7 = per-cluster scores for "
        "nr clusters 0-5. The published code reads it with "
        "pd.read_csv(..., header=None, index_col=0) and takes idxmax(axis=1) "
        "as the cluster assignment",
    ),
}

#: Columns lifted out of nr_pseudotime_BCC_UCC.csv.  ``leiden`` is renamed
#: because a bare "leiden" next to nr_cluster/cluster would be ambiguous; it is
#: a finer (11-cluster) clustering from the pseudotime analysis.
PSEUDOTIME_COLUMNS = {
    "leiden": "nr_leiden",
    "sequence_568349016": "sequence_568349016",
    "pseudotime_568349016": "pseudotime_568349016",
    "sequence_490124": "sequence_490124",
    "pseudotime_490124": "pseudotime_490124",
    "pseudotime_UCC": "pseudotime_UCC",
}

#: Prefix identifying the non-reactivated cells, the only ones the pseudotime
#: analysis covers.
NR_PREFIX = "nr_"

H5AD_NAME = "figure_dataset.h5ad"
MANIFEST_NAME = "MANIFEST.json"
FACTORS_NAME = "logcounts_scaled_factors.csv.gz"
FACTORS_VARM_KEY = "logcounts_scaled_factors"

#: The genes the analysis object was subset away from.  ToxoDB-65 has 8322
#: protein-coding genes; the deposited object keeps the 8170 that sit on the 14
#: nuclear chromosomes.  The other 152 are on unplaced ``KE*`` contigs and carry
#: the apicoplast and mitochondrial transcripts that the Supplementary 5 volcano
#: panels plot at the positive extreme, so they are shipped as a second, narrow
#: matrix rather than folded into ``logcounts.mtx.gz``, which stays exactly the
#: deposited layer.
EXTRA_MATRIX_NAME = "logcounts_extra.mtx.gz"
EXTRA_VAR_NAME = "var_extra.csv.gz"

DEFAULT_METHODS_REPO = Path(
    "/home/agent/workspaces/BradyzoiteHeterogeneity-methods"
)
DEFAULT_H5AD = (
    DEFAULT_METHODS_REPO
    / "data"
    / "integrated_adata_me49_nr_subset_annotated_sparse.h5ad"
)
DEFAULT_PSEUDOTIME = DEFAULT_METHODS_REPO / "data" / "nr_pseudotime_BCC_UCC.csv"
DEFAULT_NR_H5 = DEFAULT_METHODS_REPO / "data" / "non-reactivated_v65.h5"
DEFAULT_ME49_H5AD = DEFAULT_METHODS_REPO / "data" / "011_me49_filtered.h5ad"
DEFAULT_GFF = (
    DEFAULT_METHODS_REPO
    / "data"
    / "input"
    / "genome"
    / "annotations"
    / "ToxoDB-65_TgondiiME49.gff"
)
DEFAULT_OUTDIR = Path(__file__).resolve().parent.parent / "data"

DESCRIPTIONS: dict[str, str] = {
    "counts.mtx.gz": (
        "OPTIONAL (--include-counts): MatrixMarket (integer, coordinate) UMI "
        "count matrix, rows = cells, columns = genes; row order matches "
        "obs.csv.gz, column order matches var.csv.gz. No figure panel uses it"
    ),
    "logcounts.mtx.gz": (
        "MatrixMarket (real, coordinate) log-normalised expression, rows = "
        "cells, columns = genes; float32 values written with shortest "
        "round-trip decimal representation"
    ),
    "logcounts_scaled.mtx.gz": (
        "MatrixMarket (real, coordinate) gene-scaled log expression, rows = "
        "cells, columns = genes; float64 values written with shortest "
        "round-trip decimal representation. Exactly = logcounts * a per-"
        "(dataset, gene) constant -- see logcounts_scaled_factors.csv.gz"
    ),
    EXTRA_MATRIX_NAME: (
        "MatrixMarket (real, coordinate) log-normalised expression of the 152 "
        "ToxoDB-65 genes the deposited object does not carry, rows = cells, "
        "columns = genes; row order matches obs.csv.gz, column order matches "
        "var_extra.csv.gz. Same normalisation as logcounts.mtx.gz, so the two "
        "side by side are the 8322-gene universe the published differential "
        "expression ran on (bzfig.de)"
    ),
    EXTRA_VAR_NAME: (
        "Per-gene metadata for logcounts_extra.mtx.gz, index column 'gene_id'; "
        "same columns as var.csv.gz. Every one of these genes sits on an "
        "unplaced KE* contig"
    ),
    "obs.csv.gz": (
        "Per-cell metadata, index column 'barcode'; includes cell_cycle_group "
        "(CCC/MCC) and the pseudotime columns. Factor level order is in the "
        "'categoricals' section of MANIFEST.json"
    ),
    "var.csv.gz": (
        "Per-gene metadata, index column 'gene_id'; includes gene_name and "
        "gene_description"
    ),
    "uns_colors.json": "All *_colors palettes from adata.uns, as hex strings",
    FACTORS_NAME: (
        "Long-form table (dataset, gene_id, factor), float64. Required: "
        "logcounts_scaled[i, j] = logcounts[i, j] * factor[obs.dataset[i], j]. "
        "The reconstruction is bit-identical to the source layer in float32 "
        "and accurate to 1 ulp (4.4e-16 relative) in float64, which is why "
        "logcounts_scaled.mtx.gz is not shipped"
    ),
    H5AD_NAME: (
        "Stripped AnnData version of this package for scanpy users. X is None "
        "and the matrices live in .layers, so nothing is duplicated. "
        "varm['logcounts_scaled_factors'] (n_vars x n_datasets, column order in "
        "uns['logcounts_scaled_factors_datasets']) lets logcounts_scaled be "
        "rebuilt without any other file"
    ),
}


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    """Return the hex sha256 digest of *path*."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def open_gz(path: Path, level: int) -> gzip.GzipFile:
    """Open *path* for deterministic gzip writing (mtime pinned to 0)."""
    return gzip.GzipFile(filename=str(path), mode="wb", compresslevel=level, mtime=0)


def gzip_csv_kwargs(level: int) -> dict:
    """``compression=`` argument for ``DataFrame.to_csv`` that is reproducible."""
    return {"method": "gzip", "compresslevel": level, "mtime": 0}


def as_object_index(index: pd.Index) -> pd.Index:
    """Force a string index to numpy ``object`` dtype.

    pandas >= 3 defaults to the new ``str`` extension dtype, which anndata
    refuses to write unless ``allow_write_nullable_strings`` is enabled -- and
    the resulting encoding is unreadable by anndata < 0.11 and by most R
    readers.  Coercing back to ``object`` keeps the h5ad maximally portable.
    """
    if index.dtype.kind in "OU" or str(index.dtype) == "str":
        return pd.Index(np.asarray(index, dtype=object), dtype=object, name=index.name)
    return index


def legacy_strings(df: pd.DataFrame) -> pd.DataFrame:
    """Return *df* with every string / categorical column in ``object`` dtype."""
    out = df.copy()
    out.index = as_object_index(df.index)
    for col in out.columns:
        s = out[col]
        if isinstance(s.dtype, pd.CategoricalDtype):
            cats = s.cat.categories
            if cats.dtype.kind in "OU" or str(cats.dtype) == "str":
                cats = pd.Index(np.asarray(cats, dtype=object), dtype=object)
            out[col] = pd.Categorical.from_codes(
                s.cat.codes.to_numpy(), categories=cats, ordered=s.cat.ordered
            )
        elif s.dtype.kind == "O" or str(s.dtype) == "str":
            out[col] = pd.Series(
                np.asarray(s, dtype=object), index=out.index, dtype=object
            )
    return out


#: Factors with more levels than this are recorded by count only -- listing
#: e.g. all 3244 gene_description strings would triple the manifest and the
#: CSV already carries the values verbatim.
MAX_LISTED_LEVELS = 50


def categorical_levels(df: pd.DataFrame) -> dict[str, dict]:
    """Record factor level order / orderedness so CSV consumers can restore it."""
    levels: dict[str, dict] = {}
    for col in df.columns:
        if not isinstance(df[col].dtype, pd.CategoricalDtype):
            continue
        cats = df[col].cat.categories
        entry: dict = {
            "n_categories": int(len(cats)),
            "ordered": bool(df[col].cat.ordered),
        }
        if len(cats) <= MAX_LISTED_LEVELS:
            entry["categories"] = [
                int(c) if isinstance(c, (int, np.integer)) else str(c) for c in cats
            ]
        else:
            entry["categories"] = None
            entry["note"] = (
                "too many levels to list; the level order is alphabetical in "
                "the source object and the CSV carries the values verbatim"
            )
        levels[col] = entry
    return levels


def _string_column_dtypes(df: pd.DataFrame) -> dict[str, str]:
    """``dtype=`` map for ``read_csv`` keeping factor/text columns as text."""
    out = {}
    for col in df.columns:
        dt = df[col].dtype
        if (
            isinstance(dt, pd.CategoricalDtype)
            or dt.kind == "O"
            or str(dt) == "str"
        ):
            out[col] = "str"
    return out


def to_csr(mat) -> sp.csr_matrix:
    """Coerce to CSR with sorted indices, preserving explicitly-stored zeros."""
    out = sp.csr_matrix(mat)
    out.sort_indices()
    return out


def as_sorted_coo_triplets(mat) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (row, col, value) sorted by (row, col); keeps stored zeros."""
    csr = to_csr(mat)
    rows = np.repeat(np.arange(csr.shape[0], dtype=np.int64), np.diff(csr.indptr))
    return rows, csr.indices.astype(np.int64), csr.data


# --------------------------------------------------------------------------
# Building the package contents
# --------------------------------------------------------------------------


def build_obs(adata, pseudotime_csv: Path) -> pd.DataFrame:
    """Assemble the per-cell table: source obs + pseudotime + CCC/MCC group.

    ``cell_cycle_group`` is 'CCC' for the non-reactivated cells that carry a
    non-null ``pseudotime_UCC`` (the unidirectional/committed cell-cycle
    trajectory), 'MCC' for the remaining non-reactivated cells, and missing for
    every cell outside the non-reactivated dataset.
    """
    obs = adata.obs.copy()
    obs.index = pd.Index(np.asarray(adata.obs_names, dtype=object), name="barcode")

    pt = pd.read_csv(pseudotime_csv)
    if "barcode" not in pt.columns:
        raise ValueError(f"{pseudotime_csv} has no 'barcode' column")
    pt = pt.set_index("barcode")
    pt.index = pd.Index(np.asarray(pt.index, dtype=object), name="barcode")

    missing = set(pt.index) - set(obs.index)
    if missing:
        raise ValueError(
            f"{len(missing)} pseudotime barcodes are absent from the h5ad "
            f"(e.g. {sorted(missing)[:3]})"
        )

    is_nr = np.array([str(b).startswith(NR_PREFIX) for b in obs.index])
    nr_barcodes = set(obs.index[is_nr])
    uncovered = nr_barcodes - set(pt.index)
    if uncovered:
        raise ValueError(
            f"{len(uncovered)} '{NR_PREFIX}' cells have no pseudotime row "
            f"(e.g. {sorted(uncovered)[:3]})"
        )

    for src, dst in PSEUDOTIME_COLUMNS.items():
        if src not in pt.columns:
            raise ValueError(f"{pseudotime_csv} has no '{src}' column")
        obs[dst] = pt[src].reindex(obs.index).to_numpy()

    has_ucc = obs["pseudotime_UCC"].notna().to_numpy()
    group = np.full(len(obs), None, dtype=object)
    group[is_nr & has_ucc] = "CCC"
    group[is_nr & ~has_ucc] = "MCC"
    obs["cell_cycle_group"] = pd.Categorical(
        group,
        categories=pd.Index(np.array(["CCC", "MCC"], dtype=object), dtype=object),
        ordered=False,
    )

    n_ccc = int((group == "CCC").sum())
    n_mcc = int((group == "MCC").sum())
    print(
        f"  cell_cycle_group: CCC={n_ccc}  MCC={n_mcc}  "
        f"NA(non-{NR_PREFIX.rstrip('_')})={len(obs) - n_ccc - n_mcc}"
    )

    # Drop the obs columns that merely restate an obsm embedding, after
    # proving that they really are identical.
    for col, (key, pos) in OBS_EMBEDDING_DUPLICATES.items():
        if col not in obs.columns:
            continue
        if key not in adata.obsm:
            raise ValueError(f"cannot drop obs['{col}']: obsm['{key}'] missing")
        a = obs[col].to_numpy(dtype=np.float64)
        b = np.asarray(adata.obsm[key], dtype=np.float64)[:, pos]
        if not np.array_equal(a, b):
            raise ValueError(
                f"obs['{col}'] is not identical to obsm['{key}'][:, {pos}] "
                f"(max abs diff {np.max(np.abs(a - b))}); refusing to drop it"
            )
    dropped = [c for c in OBS_EMBEDDING_DUPLICATES if c in obs.columns]
    print(f"  dropping {len(dropped)} obs columns duplicated by obsm: {dropped}")

    for col in KEEP_OBS:
        if col not in obs.columns:
            raise ValueError(f"required obs column '{col}' missing from source")

    keep = KEEP_OBS + [
        PSEUDOTIME_COLUMNS[c] for c in PSEUDOTIME_COLUMNS
    ] + ["cell_cycle_group"]
    unused = [c for c in obs.columns if c not in keep]
    print(f"  dropping {len(unused)} unused obs columns: {unused}")
    return legacy_strings(obs[keep])


def build_var(adata) -> pd.DataFrame:
    """Assemble the per-gene table (index = gene id)."""
    var = adata.var.copy()
    var.index = pd.Index(np.asarray(adata.var_names, dtype=object), name="gene_id")
    for col in KEEP_VAR:
        if col not in var.columns:
            raise ValueError(f"required var column '{col}' missing from source")
    unused = [c for c in var.columns if c not in KEEP_VAR]
    print(f"  dropping {len(unused)} unused var columns: {unused}")
    return legacy_strings(var[KEEP_VAR])


def read_gff_genes(path: Path) -> pd.DataFrame:
    """seqid / gene_name / gene_description per gene, from a ToxoDB GFF.

    The deposited ``var`` table was built from this file; the export checks that
    this parse reproduces it before trusting it for the 152 extra genes.
    """
    from urllib.parse import unquote

    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "protein_coding_gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in fields[8].split(";") if "=" in kv)
            rows[attrs["ID"]] = (
                fields[0],
                unquote(attrs.get("Name", "")),
                unquote(attrs.get("description", "")),
            )
    frame = pd.DataFrame.from_dict(
        rows, orient="index", columns=["seqid", "gene_name", "gene_description"]
    )
    frame.index = pd.Index(np.asarray(frame.index, dtype=object), name="gene_id")
    return frame


def log_normalise(counts, totals: np.ndarray) -> sp.csr_matrix:
    """``log1p(1e4 * counts / totals)`` in float32, one total per cell."""
    out = to_csr(counts).astype(np.float32)
    out = to_csr(sp.diags((np.float32(1e4) / totals).astype(np.float32)) @ out)
    out.data = np.log1p(out.data).astype(np.float32)
    return out


def max_abs_diff(a, b) -> float:
    """Largest absolute difference between two matrices of the same shape."""
    diff = to_csr(a).astype(np.float64) - to_csr(b).astype(np.float64)
    return float(np.abs(diff.data).max()) if diff.nnz else 0.0


def build_extra_genes(adata, nr_h5: Path, me49_h5ad: Path, gff: Path):
    """Normalised expression of the genes the deposited object does not carry.

    The two objects the deposited one was integrated from still have them: the
    non-reactivated 10x run (8322 genes) and the me49 object (8496).  Upstream
    the normalisation ran per source object and *before* the gene subset -- each
    cell divided by its own total over that object's full gene set, scaled to
    1e4 and log1p'd.  Applied to the 8170 shipped genes that recipe reproduces
    the deposited ``logcounts`` layer to one float32 ulp, which is what makes it
    the right recipe for the 152 that are missing; the returned diagnostics
    record the two denominators that do not work.

    Returns the matrix (cells x 152, rows in ``obs`` order), the gene table and
    the diagnostics.
    """
    import anndata as ad
    import scanpy as sc

    nr = sc.read_10x_h5(nr_h5)
    nr.var_names_make_unique()
    me49 = ad.read_h5ad(me49_h5ad)

    universe = [str(gene) for gene in nr.var_names]
    shipped = [str(gene) for gene in adata.var_names]
    if not set(shipped) <= set(universe):
        raise ValueError(f"{nr_h5.name} is missing genes the deposited object has")
    if not set(universe) <= set(map(str, me49.var_names)):
        raise ValueError(f"{me49_h5ad.name} is missing genes the 10x universe has")
    extra = [gene for gene in universe if gene not in set(shipped)]

    barcodes = np.asarray(list(map(str, adata.obs_names)), dtype=object)
    is_nr = np.array([barcode.startswith(NR_PREFIX) for barcode in barcodes])
    nr_cells = [barcode[len(NR_PREFIX) :] for barcode in barcodes[is_nr]]
    me49_cells = list(barcodes[~is_nr])

    # The two blocks stack nr-first; this puts them back in the deposited order.
    order = np.empty(len(barcodes), dtype=np.int64)
    order[np.flatnonzero(is_nr)] = np.arange(is_nr.sum())
    order[np.flatnonzero(~is_nr)] = is_nr.sum() + np.arange((~is_nr).sum())
    counts = to_csr(
        sp.vstack([to_csr(nr[nr_cells, universe].X), to_csr(me49[me49_cells, universe].X)])
    )[order]

    totals = np.empty(len(barcodes), dtype=np.float32)
    totals[is_nr] = np.asarray(to_csr(nr[nr_cells].X).sum(1)).ravel()
    totals[~is_nr] = np.asarray(to_csr(me49[me49_cells].X).sum(1)).ravel()
    if totals.min() <= 0:
        raise ValueError("a cell has no counts at all; there is nothing to divide by")

    position = {gene: i for i, gene in enumerate(universe)}
    shared_columns = np.array([position[gene] for gene in shipped])
    extra_columns = np.array([position[gene] for gene in extra])
    shared_counts = counts[:, shared_columns]

    # The counts pulled out of the source objects have to be the ones the
    # deposited object kept, or the normalisation below is of something else.
    counts_diff = max_abs_diff(shared_counts, adata.layers["raw_counts"])
    if counts_diff != 0.0:
        raise ValueError(
            f"source counts differ from the deposited raw_counts layer by "
            f"{counts_diff}; the cell or gene alignment is wrong"
        )

    deposited = to_csr(adata.layers["logcounts"])
    denominators = {
        "per source object": totals,
        "all 8322 genes": np.asarray(counts.sum(1)).ravel().astype(np.float32),
        "the 8170 deposited genes": np.asarray(shared_counts.sum(1))
        .ravel()
        .astype(np.float32),
    }
    logcounts_diff = {
        name: max_abs_diff(log_normalise(shared_counts, denominator), deposited)
        for name, denominator in denominators.items()
    }
    if logcounts_diff["per source object"] > 1e-6:
        raise ValueError(
            "the per-source-object normalisation no longer reproduces the "
            f"deposited logcounts layer (max abs diff {logcounts_diff['per source object']})"
        )

    annotation = read_gff_genes(gff)
    missing = [gene for gene in universe if gene not in annotation.index]
    if missing:
        raise ValueError(f"{gff.name} has no gene record for {missing[:3]} ({len(missing)} genes)")
    for column in KEEP_VAR:
        want = np.array(
            [("" if pd.isna(v) else str(v)) for v in adata.var[column]], dtype=object
        )
        got = np.asarray(annotation.loc[shipped, column], dtype=object)
        if not np.array_equal(want, got):
            raise ValueError(f"{gff.name} does not reproduce the deposited var['{column}']")

    matrix = log_normalise(counts[:, extra_columns], totals)
    var = legacy_strings(annotation.loc[extra, KEEP_VAR])
    diagnostics = {
        "n_vars": len(extra),
        "sources": {
            "counts_10x_h5": str(nr_h5),
            "counts_h5ad": str(me49_h5ad),
            "annotation_gff": str(gff),
        },
        "universe": {
            "toxodb_65_genes": len(universe),
            "deposited_genes": len(shipped),
            "extra_genes": len(extra),
        },
        "contigs": {
            "n_contigs": int(var["seqid"].nunique()),
            "genes_per_contig_where_more_than_one": {
                str(seqid): int(n)
                for seqid, n in var["seqid"].astype(str).value_counts().items()
                if n > 1
            },
        },
        "normalisation": (
            "log1p(1e4 * counts / total counts of that cell in its source "
            "object), float32; the same recipe and the same per-cell totals as "
            "the deposited logcounts layer"
        ),
        "row_order": "same as obs.csv.gz",
        "counts_match_deposited_raw_counts": True,
        "deposited_logcounts_max_abs_diff": logcounts_diff,
        "deposited_logcounts_diff_note": (
            "how closely each denominator reproduces the deposited logcounts "
            "layer on the 8170 genes it holds. The per-source-object totals are "
            "within one float32 ulp; a single 8322-gene total is wrong for the "
            "me49 cells, whose source object has 8496 genes; the 8170-gene "
            "total -- all the deposited object can offer -- is wrong for every "
            "cell, which is why the published DE could not be reproduced from it"
        ),
        "stored_entries": int(matrix.nnz),
    }
    return matrix, var, diagnostics


def scaled_factor_table(
    logcounts, logcounts_scaled, dataset_codes: np.ndarray, dataset_names: list[str],
    var_names: np.ndarray,
) -> tuple[pd.DataFrame, dict]:
    """Derive the per-(dataset, gene) constants c with ``scaled = log * c``.

    Returns the long-form table and a dict of diagnostics (how exactly the
    reconstruction reproduces the shipped matrix).
    """
    L = to_csr(logcounts)
    S = to_csr(logcounts_scaled)
    if L.nnz != S.nnz or not np.array_equal(L.indices, S.indices) or not np.array_equal(
        L.indptr, S.indptr
    ):
        return pd.DataFrame(), {"same_pattern": False}

    n_genes = L.shape[1]
    rows = np.repeat(np.arange(L.shape[0], dtype=np.int64), np.diff(L.indptr))
    cols = L.indices.astype(np.int64)
    groups = dataset_codes[rows].astype(np.int64) * n_genes + cols

    order = np.argsort(groups, kind="stable")
    gsorted = groups[order]
    first = np.flatnonzero(np.r_[True, gsorted[1:] != gsorted[:-1]])

    n_datasets = len(dataset_names)
    factors = np.full(n_datasets * n_genes, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = S.data.astype(np.float64) / L.data.astype(np.float64)
    factors[gsorted[first]] = ratio[order][first]

    reconstructed = L.data.astype(np.float64) * factors[groups]
    finite = np.isfinite(reconstructed) & np.isfinite(S.data)
    max_abs = float(np.max(np.abs(reconstructed[finite] - S.data[finite])))
    max_rel = float(
        np.max(np.abs(reconstructed[finite] / S.data[finite] - 1.0))
    )

    table = pd.DataFrame(
        {
            "dataset": np.repeat(
                np.asarray(dataset_names, dtype=object), n_genes
            ),
            "gene_id": np.tile(np.asarray(var_names, dtype=object), n_datasets),
            "factor": factors,
        }
    )
    diagnostics = {
        "same_pattern": True,
        "grouping_obs_column": "dataset",
        "datasets": [str(d) for d in dataset_names],
        "n_groups_determined": int(np.isfinite(factors).sum()),
        "n_groups_total": int(factors.size),
        "n_groups_undetermined": int(factors.size - np.isfinite(factors).sum()),
        "undetermined_note": (
            "a (dataset, gene) pair with no stored logcounts entry has no "
            "observable factor; it is written as an empty field (NaN) and "
            "there is nothing to reconstruct for it"
        ),
        "reconstruction_max_abs_diff": max_abs,
        "reconstruction_max_rel_diff": max_rel,
        "reconstruction_bit_identical_float64": bool(
            np.array_equal(reconstructed, S.data)
        ),
        "reconstruction_bit_identical_float32": bool(
            np.array_equal(
                reconstructed.astype(np.float32), S.data.astype(np.float32)
            )
        ),
    }
    return table, diagnostics


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------


def write_layer_mtx(path: Path, mat, dtype: str, name: str, level: int) -> dict:
    """Write one layer as a gzipped MatrixMarket file; return size statistics."""
    csr = to_csr(mat)
    data = csr.data
    if dtype.startswith("int"):
        if not np.array_equal(data, np.round(data)):
            raise ValueError(f"layer '{name}' is not integral, refusing int export")
        info = np.iinfo(np.dtype(dtype))
        if data.min() < info.min or data.max() > info.max:
            raise ValueError(f"layer '{name}' does not fit into {dtype}")
        out = csr.copy()
        out.data = data.astype(dtype)
    else:
        out = csr.astype(dtype) if csr.dtype != np.dtype(dtype) else csr

    n_stored = int(out.nnz)
    n_true = int(np.count_nonzero(out.data))
    comment = (
        f"\n layer: {name}\n"
        f" orientation: rows = cells (obs), columns = genes (var)\n"
        f" shape: {out.shape[0]} cells x {out.shape[1]} genes\n"
        f" row i <-> row i of obs.csv.gz; column j <-> row j of var.csv.gz\n"
        f" values: {dtype}, written with shortest round-trip decimal form\n"
        f" stored entries: {n_stored} ({n_stored - n_true} of them are "
        f"explicit zeros kept for fidelity with the source object)"
    )

    buf = io.BytesIO()
    mmwrite(buf, out.tocoo(), comment=comment, precision=None)
    raw_bytes = buf.tell()
    with open_gz(path, level) as fh:
        fh.write(buf.getbuffer())

    return {
        "layer": name,
        "dtype_source": str(np.asarray(csr.data).dtype),
        "dtype_written": dtype,
        "stored_entries": n_stored,
        "true_nonzeros": n_true,
        "explicit_zeros": n_stored - n_true,
        "density": n_stored / (out.shape[0] * out.shape[1]),
        "mtx_uncompressed_bytes": int(raw_bytes),
        "mtx_gz_bytes": int(path.stat().st_size),
    }


def write_h5ad(path: Path, obs: pd.DataFrame, var: pd.DataFrame,
               layers: dict, obsm: dict, varm: dict, uns: dict) -> None:
    """Write the stripped AnnData convenience copy (``X`` intentionally None)."""
    import anndata as ad

    adata = ad.AnnData(
        obs=obs,
        var=var,
        layers={k: to_csr(v) for k, v in layers.items()},
        obsm={k: np.asarray(v) for k, v in obsm.items()},
        varm={k: np.asarray(v) for k, v in varm.items()},
        uns=dict(uns),
    )
    adata.write_h5ad(path, compression="gzip")


# --------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------


def verify_supplementary_5(outdir: Path) -> dict:
    """Re-run the Supplementary 5 volcano DE from the package just written.

    The point is that the package is self-sufficient: this loads nothing but
    the files in *outdir* and checks the counts against the ones recorded in
    ``bzfig.constants``.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
    from bzfig import constants, de
    from bzfig.data import load_dataset

    published = {"5C": constants.SUPP5C_DE_COUNTS, "5F": constants.SUPP5F_DE_COUNTS}
    adata = load_dataset(outdir)
    results = {}
    for panel in de.COMPARISONS:
        table = de.volcano_table(adata, panel, outdir)
        counts = de.counts(table)
        expected = constants.SUPP5_DE_COUNTS_REPRODUCED[panel]
        assert counts == expected, f"Supplementary {panel}: {counts} != {expected}"
        results[panel] = {
            **counts,
            "published": published.get(panel),
            "genes_detected_in_both_groups": int(len(table)),
            "log2fc_range": [float(table["log2fc"].min()), float(table["log2fc"].max())],
        }
        print(f"  Supplementary {panel}: up {counts['up']}, down {counts['down']}")
    return results


def verify(outdir: Path, adata, obs: pd.DataFrame, var: pd.DataFrame,
           obsm_keys: list[str], uns_colors: dict,
           export_layers: dict[str, str], extras_root: Path,
           extra_matrix=None, extra_var: pd.DataFrame | None = None,
           check_h5ad: bool = True) -> dict:
    """Re-read every exported file and assert it matches the in-memory source.

    Returns a dict of measured differences.  Raises ``AssertionError`` on any
    mismatch.
    """
    results: dict = {"layers": {}, "obs": {}, "var": {}, "obsm": {}, "h5ad": {}}

    # ---- MatrixMarket layers -------------------------------------------
    for name, dtype in export_layers.items():
        src = to_csr(adata.layers[name])
        with gzip.open(outdir / f"{name}.mtx.gz", "rb") as fh:
            back = to_csr(mmread(fh))
        assert back.shape == src.shape, f"{name}: shape {back.shape} != {src.shape}"
        assert back.nnz == src.nnz, f"{name}: nnz {back.nnz} != {src.nnz}"
        assert np.array_equal(back.indices, src.indices), f"{name}: column indices"
        assert np.array_equal(back.indptr, src.indptr), f"{name}: row pointers"
        # MatrixMarket is decimal text and mmread always hands back float64.
        # The file holds the *shortest round-trip decimal form* of each source
        # value, so casting the parsed value back to the source dtype restores
        # the original bit pattern.  For a float32 layer the intermediate
        # float64 differs from the widened float32 by up to half a float32 ulp
        # (~2.4e-07 here); that is a property of decimal text, not data loss,
        # and it is reported as `max_abs_diff_as_float64`.
        want = src.data
        got_raw = back.data
        got = got_raw.astype(want.dtype)
        max_abs = float(np.max(np.abs(got.astype(np.float64)
                                      - want.astype(np.float64)))) if len(want) else 0.0
        max_abs_f64 = float(np.max(np.abs(got_raw.astype(np.float64)
                                          - want.astype(np.float64)))) if len(want) else 0.0
        exact = bool(np.array_equal(got, want))
        assert exact, f"{name}: values differ, max abs diff {max_abs}"
        results["layers"][name] = {
            "shape": list(src.shape),
            "stored_entries": int(src.nnz),
            "source_dtype": str(want.dtype),
            "max_abs_diff": max_abs,
            "max_abs_diff_as_float64": max_abs_f64,
            "bit_identical": exact,
        }

    # ---- obs -------------------------------------------------------------
    # Categorical columns whose levels look numeric ('0'..'5') would be
    # re-inferred as numbers by read_csv, so read them as text -- exactly what
    # a downstream consumer should do, guided by the manifest's
    # `categoricals` section.
    back_obs = pd.read_csv(
        outdir / "obs.csv.gz",
        index_col=0,
        dtype=_string_column_dtypes(obs),
        float_precision="round_trip",
    )
    assert list(map(str, back_obs.index)) == list(map(str, obs.index)), "obs index"
    for col in obs.columns:
        want_s = obs[col]
        got_s = back_obs[col]
        if isinstance(want_s.dtype, pd.CategoricalDtype) or want_s.dtype.kind == "O":
            w = np.array(
                [None if pd.isna(v) else str(v) for v in want_s], dtype=object
            )
            g = np.array(
                [None if pd.isna(v) else str(v) for v in got_s], dtype=object
            )
            assert np.array_equal(w, g), f"obs['{col}'] categorical values differ"
            results["obs"][col] = "exact (string values)"
        else:
            # CSV holds the shortest round-trip decimal form of each value, so
            # the parsed float64 must be narrowed back to the source dtype
            # (float32 columns) before comparing -- same reasoning as for the
            # MatrixMarket layers above.
            w = want_s.to_numpy(dtype=np.float64)
            g = got_s.to_numpy(dtype=np.float64).astype(want_s.dtype).astype(np.float64)
            assert np.array_equal(np.isnan(w), np.isnan(g)), f"obs['{col}'] NaN mask"
            m = ~np.isnan(w)
            diff = float(np.max(np.abs(g[m] - w[m]))) if m.any() else 0.0
            assert diff == 0.0, f"obs['{col}'] max abs diff {diff}"
            results["obs"][col] = {
                "source_dtype": str(want_s.dtype),
                "max_abs_diff": diff,
            }

    # ---- var -------------------------------------------------------------
    back_var = pd.read_csv(
        outdir / "var.csv.gz",
        index_col=0,
        dtype=_string_column_dtypes(var),
        float_precision="round_trip",
    )
    assert list(map(str, back_var.index)) == list(map(str, var.index)), "var index"
    for col in var.columns:
        w = np.array(
            [None if pd.isna(v) else str(v) for v in var[col]], dtype=object
        )
        g = np.array(
            [None if pd.isna(v) else str(v) for v in back_var[col]], dtype=object
        )
        assert np.array_equal(w, g), f"var['{col}'] values differ"
        results["var"][col] = "exact (string values)"

    # ---- obsm ------------------------------------------------------------
    for key in obsm_keys:
        src_arr = np.asarray(adata.obsm[key])
        want = src_arr.astype(np.float64)
        back = pd.read_csv(
            outdir / f"obsm_{key}.csv.gz", index_col=0,
            float_precision="round_trip",
        )
        assert list(map(str, back.index)) == list(map(str, obs.index)), f"obsm[{key}]"
        got = back.to_numpy(dtype=np.float64).astype(src_arr.dtype).astype(np.float64)
        assert got.shape == want.shape, f"obsm[{key}] shape"
        diff = float(np.max(np.abs(got - want)))
        assert diff == 0.0, f"obsm[{key}] max abs diff {diff}"
        results["obsm"][key] = {
            "shape": list(want.shape),
            "source_dtype": str(src_arr.dtype),
            "max_abs_diff": diff,
        }

    # ---- uns colors ------------------------------------------------------
    with open(outdir / "uns_colors.json") as fh:
        back_colors = json.load(fh)
    assert back_colors == uns_colors, "uns_colors.json differs"
    results["uns_colors"] = "exact"

    # ---- stripped h5ad ---------------------------------------------------
    back_ad = None
    if check_h5ad:
        import anndata as ad

        back_ad = ad.read_h5ad(outdir / H5AD_NAME)
        assert list(map(str, back_ad.obs_names)) == list(map(str, obs.index))
        assert list(map(str, back_ad.var_names)) == list(map(str, var.index))
        for name in export_layers:
            src = to_csr(adata.layers[name])
            got = to_csr(back_ad.layers[name])
            assert np.array_equal(got.indices, src.indices), f"h5ad {name} indices"
            assert np.array_equal(got.indptr, src.indptr), f"h5ad {name} indptr"
            diff = float(np.max(np.abs(got.data.astype(np.float64)
                                       - src.data.astype(np.float64))))
            assert diff == 0.0, f"h5ad layer {name} max abs diff {diff}"
            results["h5ad"][name] = {"max_abs_diff": diff}
        for col in obs.columns:
            w = np.array([None if pd.isna(v) else str(v) for v in obs[col]], dtype=object)
            g = np.array(
                [None if pd.isna(v) else str(v) for v in back_ad.obs[col]], dtype=object
            )
            assert np.array_equal(w, g), f"h5ad obs['{col}']"
            if isinstance(obs[col].dtype, pd.CategoricalDtype):
                assert list(map(str, back_ad.obs[col].cat.categories)) == list(
                    map(str, obs[col].cat.categories)
                ), f"h5ad obs['{col}'] level order"
                assert bool(back_ad.obs[col].cat.ordered) == bool(obs[col].cat.ordered)
        for col in var.columns:
            w = np.array([None if pd.isna(v) else str(v) for v in var[col]], dtype=object)
            g = np.array(
                [None if pd.isna(v) else str(v) for v in back_ad.var[col]], dtype=object
            )
            assert np.array_equal(w, g), f"h5ad var['{col}']"
        for key in obsm_keys:
            diff = float(
                np.max(
                    np.abs(
                        np.asarray(back_ad.obsm[key], dtype=np.float64)
                        - np.asarray(adata.obsm[key], dtype=np.float64)
                    )
                )
            )
            assert diff == 0.0, f"h5ad obsm[{key}] max abs diff {diff}"
        for key, value in uns_colors.items():
            assert list(map(str, back_ad.uns[key])) == list(value), f"h5ad uns['{key}']"
        results["h5ad"]["obs_var_obsm_uns"] = "exact"

    # ---- logcounts_scaled reconstruction --------------------------------
    # The scaled layer is not shipped as a matrix; prove that the shipped
    # logcounts matrix plus the factor table rebuild it.
    factors_path = outdir / FACTORS_NAME
    if "logcounts_scaled" in adata.layers and factors_path.exists():
        with gzip.open(outdir / "logcounts.mtx.gz", "rb") as fh:
            L = to_csr(mmread(fh))
        S = to_csr(adata.layers["logcounts_scaled"])
        assert np.array_equal(L.indices, S.indices), "scaled layer pattern"
        assert np.array_equal(L.indptr, S.indptr), "scaled layer pattern"

        fac = pd.read_csv(
            factors_path, dtype={"dataset": "str", "gene_id": "str"},
            float_precision="round_trip",
        )
        datasets = list(dict.fromkeys(fac["dataset"]))
        n_genes = L.shape[1]
        assert len(fac) == len(datasets) * n_genes, "factor table length"
        assert list(fac["gene_id"][:n_genes]) == list(map(str, var.index))
        wide = fac["factor"].to_numpy(dtype=np.float64).reshape(len(datasets), n_genes)

        ds_pos = {d: i for i, d in enumerate(datasets)}
        codes = np.array([ds_pos[str(d)] for d in back_obs["dataset"]])
        rows = np.repeat(np.arange(L.shape[0], dtype=np.int64), np.diff(L.indptr))
        # L.data came back from decimal text as float64; narrow to the float32
        # the file actually encodes before applying the factor.
        rebuilt = L.data.astype(np.float32).astype(np.float64) * wide[
            codes[rows], L.indices
        ]
        want = S.data.astype(np.float64)
        max_abs = float(np.max(np.abs(rebuilt - want)))
        max_rel = float(np.max(np.abs(rebuilt / want - 1.0)))
        f32_exact = bool(
            np.array_equal(rebuilt.astype(np.float32), want.astype(np.float32))
        )
        assert f32_exact, "logcounts_scaled reconstruction is not float32-exact"
        results["logcounts_scaled_reconstruction"] = {
            "rule": (
                "logcounts_scaled[i, j] = logcounts[i, j] * "
                "factor[obs.dataset[i], j]"
            ),
            "max_abs_diff_float64": max_abs,
            "max_rel_diff_float64": max_rel,
            "bit_identical_float64": bool(np.array_equal(rebuilt, want)),
            "bit_identical_float32": f32_exact,
        }

        # ... and that the .h5ad carries the same factors in varm.
        if check_h5ad:
            varm = np.asarray(back_ad.varm[FACTORS_VARM_KEY], dtype=np.float64)
            assert (
                list(map(str, back_ad.uns[FACTORS_VARM_KEY + "_datasets"]))
                == datasets
            )
            assert np.array_equal(np.isnan(varm.T), np.isnan(wide)), "h5ad varm NaNs"
            m = ~np.isnan(wide)
            assert np.array_equal(varm.T[m], wide[m]), "h5ad varm factors"
            results["logcounts_scaled_reconstruction"]["h5ad_varm"] = "identical"

    # ---- the genes the deposited object dropped -------------------------
    if extra_matrix is not None:
        src = to_csr(extra_matrix)
        with gzip.open(outdir / EXTRA_MATRIX_NAME, "rb") as fh:
            back = to_csr(mmread(fh))
        assert back.shape == src.shape, EXTRA_MATRIX_NAME
        assert np.array_equal(back.indices, src.indices), f"{EXTRA_MATRIX_NAME} columns"
        assert np.array_equal(back.indptr, src.indptr), f"{EXTRA_MATRIX_NAME} rows"
        assert np.array_equal(back.data.astype(np.float32), src.data), EXTRA_MATRIX_NAME
        back_var = pd.read_csv(
            outdir / EXTRA_VAR_NAME, index_col=0, dtype=_string_column_dtypes(extra_var)
        )
        assert list(map(str, back_var.index)) == list(map(str, extra_var.index))
        for col in extra_var.columns:
            # A gene with no name is an empty field in the CSV and comes back as
            # NaN, which is how var.csv.gz already carries it.
            w = np.array(["" if pd.isna(v) else str(v) for v in extra_var[col]], dtype=object)
            g = np.array(["" if pd.isna(v) else str(v) for v in back_var[col]], dtype=object)
            assert np.array_equal(w, g), f"{EXTRA_VAR_NAME}: {col}"
        results["extra_genes"] = {
            "shape": list(src.shape),
            "stored_entries": int(src.nnz),
            "source_dtype": str(src.data.dtype),
            "max_abs_diff": 0.0,
            "bit_identical": True,
        }

    # ---- verbatim extra files -------------------------------------------
    results["extra_files"] = {}
    for dest_name, (rel_src, _) in EXTRA_FILES.items():
        src_bytes = (extras_root / rel_src).read_bytes()
        got_bytes = (outdir / dest_name).read_bytes()
        assert got_bytes == src_bytes, f"{dest_name} differs from {rel_src}"
        results["extra_files"][dest_name] = "byte-identical to source"

    return results


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--h5ad", type=Path, default=DEFAULT_H5AD,
        help="source AnnData object (default: %(default)s)",
    )
    parser.add_argument(
        "--pseudotime", type=Path, default=DEFAULT_PSEUDOTIME,
        help="nr_pseudotime_BCC_UCC.csv (default: %(default)s)",
    )
    parser.add_argument(
        "--nr-h5", type=Path, default=DEFAULT_NR_H5,
        help="10x .h5 of the non-reactivated run, which still has all 8322 "
             "ToxoDB-65 genes (default: %(default)s)",
    )
    parser.add_argument(
        "--me49-h5ad", type=Path, default=DEFAULT_ME49_H5AD,
        help="me49 object the in vitro cells came from (default: %(default)s)",
    )
    parser.add_argument(
        "--gff", type=Path, default=DEFAULT_GFF,
        help="ToxoDB GFF the gene annotation comes from (default: %(default)s)",
    )
    parser.add_argument(
        "--outdir", type=Path, default=DEFAULT_OUTDIR,
        help="directory to write the package into (default: %(default)s)",
    )
    parser.add_argument(
        "--gzip-level", type=int, default=9, choices=range(1, 10), metavar="1-9",
        help="gzip compression level (default: %(default)s)",
    )
    parser.add_argument(
        "--no-h5ad", action="store_true",
        help="skip the convenience figure_dataset.h5ad; it duplicates the "
             "standard-format files and is the largest item after the matrix",
    )
    parser.add_argument(
        "--no-verify", action="store_true",
        help="skip the round-trip verification pass",
    )
    parser.add_argument(
        "--include-counts", action="store_true",
        help="also ship the raw UMI counts layer (no figure panel needs it; "
             "useful only for redoing normalisation/DE from scratch)",
    )
    parser.add_argument(
        "--include-scaled-matrix", action="store_true",
        help="also ship logcounts_scaled.mtx.gz instead of relying on the "
             "(bit-exact in float32) reconstruction from logcounts and "
             "logcounts_scaled_factors.csv.gz",
    )
    parser.add_argument(
        "--no-measure-optional", action="store_true",
        help="skip measuring how big the optional layers would be",
    )
    parser.add_argument(
        "--extras-root", type=Path, default=DEFAULT_METHODS_REPO,
        help="repository the EXTRA_FILES are copied from (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    os.environ.setdefault("MPLCONFIGDIR", os.environ.get("TMPDIR", "/tmp") + "/mpl")
    import anndata as ad  # imported late so --help works without scanpy installed

    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"reading {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes; "
          f"layers={sorted(adata.layers)}")

    export_layers = dict(LAYERS)
    if args.include_scaled_matrix:
        export_layers["logcounts_scaled"] = OPTIONAL_LAYERS["logcounts_scaled"]
    if args.include_counts:
        export_layers["counts"] = OPTIONAL_LAYERS["counts"]
    for name in export_layers:
        if name not in adata.layers:
            raise ValueError(f"source object has no layer '{name}'")

    obs = build_obs(adata, args.pseudotime)
    var = build_var(adata)

    obsm_keys = list(KEEP_OBSM)
    for key in obsm_keys:
        if key not in adata.obsm:
            raise ValueError(f"required obsm key '{key}' missing from source")
    print(f"  dropping {len(adata.obsm) - len(obsm_keys)} unused obsm keys: "
          f"{[k for k in adata.obsm if k not in obsm_keys]}")

    uns_colors = {
        k: [str(x) for x in np.asarray(v).ravel()]
        for k, v in adata.uns.items()
        if k.endswith("_colors")
    }

    written: list[Path] = []
    layer_stats: list[dict] = []

    # ---- matrices --------------------------------------------------------
    for name, dtype in export_layers.items():
        path = outdir / f"{name}.mtx.gz"
        print(f"writing {path.name}")
        stats = write_layer_mtx(path, adata.layers[name], dtype, name, args.gzip_level)
        layer_stats.append(stats)
        written.append(path)
        print(
            f"  stored={stats['stored_entries']:,} "
            f"(true nnz {stats['true_nonzeros']:,}, "
            f"explicit zeros {stats['explicit_zeros']:,})  "
            f"raw={stats['mtx_uncompressed_bytes'] / 1e6:.1f} MB  "
            f"gz={stats['mtx_gz_bytes'] / 1e6:.1f} MB"
        )

    # ---- tables ----------------------------------------------------------
    print("writing obs.csv.gz / var.csv.gz")
    obs.to_csv(outdir / "obs.csv.gz", compression=gzip_csv_kwargs(args.gzip_level))
    var.to_csv(outdir / "var.csv.gz", compression=gzip_csv_kwargs(args.gzip_level))
    written += [outdir / "obs.csv.gz", outdir / "var.csv.gz"]

    for key in obsm_keys:
        arr = np.asarray(adata.obsm[key])
        cols = [f"{key}_{i + 1}" for i in range(arr.shape[1])]
        frame = pd.DataFrame(arr, index=obs.index, columns=cols)
        path = outdir / f"obsm_{key}.csv.gz"
        frame.to_csv(path, compression=gzip_csv_kwargs(args.gzip_level))
        written.append(path)
        DESCRIPTIONS[path.name] = (
            f"Embedding '{key}' ({arr.shape[0]} cells x {arr.shape[1]} dims), "
            f"indexed by barcode; also duplicated as obs columns in the source "
            f"object, which this package drops"
        )
    print(f"writing {len(obsm_keys)} obsm_*.csv.gz files")

    with open(outdir / "uns_colors.json", "w") as fh:
        json.dump(uns_colors, fh, indent=2, sort_keys=True)
        fh.write("\n")
    written.append(outdir / "uns_colors.json")

    # ---- derived factor table -------------------------------------------
    dataset_codes = obs["dataset"].cat.codes.to_numpy()
    dataset_names = [str(c) for c in obs["dataset"].cat.categories]
    factors, factor_diag = scaled_factor_table(
        adata.layers["logcounts"],
        adata.layers["logcounts_scaled"],
        dataset_codes,
        dataset_names,
        np.asarray(var.index, dtype=object),
    )
    h5ad_varm: dict = {}
    h5ad_uns = dict(uns_colors)
    if not factors.empty:
        path = outdir / FACTORS_NAME
        factors.to_csv(
            path, index=False, compression=gzip_csv_kwargs(args.gzip_level)
        )
        written.append(path)
        # Same numbers as a (n_vars x n_datasets) array, so the convenience
        # .h5ad is self-sufficient: the scaled layer can be rebuilt from
        # adata.layers["logcounts"] and adata.varm[FACTORS_VARM_KEY] alone.
        h5ad_varm[FACTORS_VARM_KEY] = (
            factors["factor"]
            .to_numpy(dtype=np.float64)
            .reshape(len(dataset_names), len(var))
            .T.copy()
        )
        h5ad_uns[FACTORS_VARM_KEY + "_datasets"] = [
            str(d) for d in dataset_names
        ]
        print(
            f"writing {FACTORS_NAME}: logcounts_scaled = logcounts * "
            f"c[dataset, gene]; reconstruction max abs diff "
            f"{factor_diag['reconstruction_max_abs_diff']:.3g}, "
            f"bit-exact in float32: "
            f"{factor_diag['reconstruction_bit_identical_float32']}"
        )

    # ---- the genes the deposited object dropped -------------------------
    print(f"writing {EXTRA_MATRIX_NAME} / {EXTRA_VAR_NAME}")
    extra_matrix, extra_var, extra_diag = build_extra_genes(
        adata, args.nr_h5, args.me49_h5ad, args.gff
    )
    extra_stats = write_layer_mtx(
        outdir / EXTRA_MATRIX_NAME,
        extra_matrix,
        "float32",
        "logcounts_extra",
        args.gzip_level,
    )
    extra_var.to_csv(outdir / EXTRA_VAR_NAME, compression=gzip_csv_kwargs(args.gzip_level))
    written += [outdir / EXTRA_MATRIX_NAME, outdir / EXTRA_VAR_NAME]
    print(
        f"  {extra_diag['n_vars']} genes on {len(extra_diag['contigs'])} unplaced "
        f"contigs, {extra_stats['stored_entries']:,} stored entries, "
        f"gz={extra_stats['mtx_gz_bytes'] / 1e6:.2f} MB; deposited logcounts "
        f"reproduced to {extra_diag['deposited_logcounts_max_abs_diff']['per source object']:.3g}"
    )

    # ---- extra verbatim files -------------------------------------------
    for dest_name, (rel_src, description) in EXTRA_FILES.items():
        src_path = args.extras_root / rel_src
        if not src_path.exists():
            raise FileNotFoundError(f"extra file {src_path} not found")
        dest = outdir / dest_name
        dest.write_bytes(src_path.read_bytes())
        DESCRIPTIONS[dest_name] = description + f" (copied from {rel_src})"
        written.append(dest)
        print(f"copying {rel_src} -> {dest_name} ({dest.stat().st_size:,} bytes)")

    # ---- stripped h5ad ---------------------------------------------------
    h5ad_path = outdir / H5AD_NAME
    if args.no_h5ad:
        h5ad_path.unlink(missing_ok=True)
        print(f"skipping {H5AD_NAME} (--no-h5ad)")
    else:
        print(f"writing {H5AD_NAME}")
        write_h5ad(
            h5ad_path,
            obs,
            var,
            {name: adata.layers[name] for name in export_layers},
            {k: adata.obsm[k] for k in obsm_keys},
            h5ad_varm,
            h5ad_uns,
        )
        written.append(h5ad_path)

    # ---- size of the layers we chose not to ship -------------------------
    optional_layers: dict[str, dict] = {}
    if not args.no_measure_optional:
        import tempfile

        for name, dtype in OPTIONAL_LAYERS.items():
            if name in export_layers or name not in adata.layers:
                continue
            with tempfile.TemporaryDirectory() as tmp:
                probe = Path(tmp) / f"{name}.mtx.gz"
                stats = write_layer_mtx(
                    probe, adata.layers[name], dtype, name, args.gzip_level
                )
            optional_layers[name] = {
                **stats,
                "shipped": False,
                "enable_with": OPTIONAL_LAYER_FLAGS[name],
                "note": (
                    "no figure panel uses this layer; adding it would grow the "
                    f"package by {stats['mtx_gz_bytes']:,} bytes of .mtx.gz "
                    "plus roughly the same again inside figure_dataset.h5ad"
                ),
            }
            print(
                f"  (not shipped) {name}: would add "
                f"{stats['mtx_gz_bytes'] / 1e6:.1f} MB as .mtx.gz"
            )

    # ---- verification ----------------------------------------------------
    verification = None
    if not args.no_verify:
        print("verifying round trip")
        verification = verify(
            outdir, adata, obs, var, obsm_keys, uns_colors, export_layers,
            args.extras_root, extra_matrix, extra_var, check_h5ad=not args.no_h5ad,
        )
        print("  all round-trip checks passed")

    # ---- manifest --------------------------------------------------------
    files = []
    for path in written:
        files.append(
            {
                "filename": path.name,
                "sha256": sha256_of(path),
                "bytes": path.stat().st_size,
                "description": DESCRIPTIONS.get(path.name, ""),
                "urls": [],
            }
        )
    files.sort(key=lambda d: d["filename"])

    manifest = {
        "name": "BradyzoiteHeterogeneity figure dataset",
        "description": (
            "Format-independent export of the single-cell object used by the "
            "figure panels. MatrixMarket matrices are cells x genes."
        ),
        "generated_by": "scripts/export_dataset.py",
        "source": {
            "h5ad": str(args.h5ad),
            "h5ad_sha256": sha256_of(args.h5ad),
            "pseudotime_csv": str(args.pseudotime),
            "pseudotime_csv_sha256": sha256_of(args.pseudotime),
        },
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "gzip_level": args.gzip_level,
        "matrix_orientation": "rows = cells (obs), columns = genes (var)",
        "layer_stats": layer_stats,
        "extra_genes": {**extra_diag, **extra_stats},
        "optional_layers": optional_layers,
        "redundancy": {
            "logcounts_from_counts": {
                "recoverable": False,
                "reason": (
                    "counts and logcounts do not share a sparsity pattern; "
                    "logcounts was derived upstream from a separate "
                    "'raw_counts' layer as log1p(1e4 * raw_counts / "
                    "rowsum(raw_counts)), and raw_counts is not shipped"
                ),
            },
            "logcounts_scaled_from_logcounts": {
                "recoverable": True,
                "transform": "logcounts_scaled = logcounts * c[dataset, gene]",
                "factor_file": FACTORS_NAME,
                **factor_diag,
            },
        },
        "categoricals": {
            "obs": categorical_levels(obs),
            "var": categorical_levels(var),
        },
        # CSV is decimal text: float32 columns are written in shortest
        # round-trip float32 form, so a reader that wants bit-identical values
        # must narrow the parsed float64 back to the dtype recorded here.
        "column_dtypes": {
            "obs": {c: str(obs[c].dtype) for c in obs.columns},
            "var": {c: str(var[c].dtype) for c in var.columns},
            "obsm": {k: str(np.asarray(adata.obsm[k]).dtype) for k in obsm_keys},
        },
        "cell_cycle_group": {
            "source": str(args.pseudotime),
            "definition": (
                "CCC = non-reactivated cell with non-null pseudotime_UCC; "
                "MCC = remaining non-reactivated cell; missing otherwise"
            ),
            "counts": {
                "CCC": int((obs["cell_cycle_group"] == "CCC").sum()),
                "MCC": int((obs["cell_cycle_group"] == "MCC").sum()),
                "NA": int(obs["cell_cycle_group"].isna().sum()),
            },
        },
        "verification": verification,
        "files": files,
        "total_bytes": sum(f["bytes"] for f in files),
    }

    def write_manifest() -> None:
        with open(outdir / MANIFEST_NAME, "w") as fh:
            json.dump(manifest, fh, indent=2, sort_keys=False)
            fh.write("\n")

    # The volcano check reads the package back through bzfig, which needs the
    # manifest on disk, so the manifest is written before the check and again
    # with its result.
    write_manifest()
    if verification is not None:
        print("re-running the Supplementary 5 differential expression from the package")
        verification["supplementary_5_de"] = verify_supplementary_5(outdir)
        write_manifest()

    total = manifest["total_bytes"]
    print(f"\n{len(files)} files, {total:,} bytes ({total / 1e6:.1f} MB) "
          f"(+ {MANIFEST_NAME})")
    for f in files:
        print(f"  {f['bytes']:>12,}  {f['sha256'][:16]}...  {f['filename']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
