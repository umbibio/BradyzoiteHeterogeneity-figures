#!/usr/bin/env python
"""Maintain the Figure 3E / 3F integration files in ``data/``.

    python scripts/export_integration.py            # -> figure_3ef_embedding.csv.gz

The integration is a second experiment: 376 in vivo tachyzoites taken from the
peritoneal cavity at 5 dpi (samples S1 and S2) integrated with the 6,505 in vivo
bradyzoites by scVI, and plotted on a UMAP of the 10-dimensional latent space.
Five files carry it, and they follow the same shape as the main package: the
count matrix, its cell and gene tables, the embedding, and the trained
checkpoint.  ``figure_3ef_counts.mtx.gz``, ``figure_3ef_obs.csv.gz``,
``figure_3ef_var.csv.gz`` and ``figure_3ef_scvi_model.pt`` are the analysis
output, assembled in the methods repository by
``notebooks/integrate.S1-S2-S3-2026-06-29.ipynb``; this script does not rebuild
them.

What it does do is regenerate the embedding from the deposited checkpoint and
write the ``figure_3ef`` section of ``MANIFEST.json`` together with fresh
checksums for every file in the package.  Nothing here edits the manifest by
hand; ``export_dataset.refresh_manifest`` owns that file.

Rendering the panels does not need this script -- the embedding it writes is
already in ``data/``.  It needs the ``integration`` extra
(``pip install -e '.[integration]'``), which pulls in torch and scvi-tools.

Loading the checkpoint is deterministic, and so is the UMAP that follows it, so
re-running this reproduces the shipped file.  Training a fresh model does not:
the layout depends on the particular fit, and the published panels were drawn
from a different one.  See docs/reproducibility.md.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

NEIGHBORS = {"n_neighbors": 50, "n_pcs": 10}

NOTE = (
    "A second experiment: in vivo tachyzoites (5 dpi, peritoneal cavity) "
    "integrated with the in vivo bradyzoites by scVI. Sample S3 (55 cells) was "
    "collected but is not part of the published integration."
)


def read_tables(datadir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    obs = pd.read_csv(datadir / "figure_3ef_obs.csv.gz", index_col="barcode", dtype="str")
    var = pd.read_csv(datadir / "figure_3ef_var.csv.gz", index_col="gene_id", dtype="str")
    return obs, var


def build(datadir: Path):
    import anndata

    obs, var = read_tables(datadir)
    with gzip.open(datadir / "figure_3ef_counts.mtx.gz", "rb") as handle:
        counts = sp.csr_matrix(scipy.io.mmread(handle)).astype(np.float32)
    if counts.shape != (len(obs), len(var)):
        raise SystemExit(f"counts shape {counts.shape} does not match {(len(obs), len(var))}")

    adata = anndata.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    return adata


def manifest_section(datadir: Path) -> dict:
    """The ``figure_3ef`` block, read off the deposited tables."""
    obs, var = read_tables(datadir)
    return {
        "n_obs": int(len(obs)),
        "n_vars": int(len(var)),
        "samples": {
            str(sample): int(n)
            for sample, n in sorted(obs["sample"].value_counts().items())
        },
        "note": NOTE,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=REPO / "data")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument(
        "--manifest-only", action="store_true",
        help="refresh MANIFEST.json from the files already in --data and stop; "
             "needs neither torch nor scvi-tools",
    )
    args = parser.parse_args()
    out = args.out or args.data / "figure_3ef_embedding.csv.gz"

    from export_dataset import refresh_manifest

    if not args.manifest_only:
        import scanpy as sc
        import scvi

        adata = build(args.data)
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")
        model = scvi.model.SCVI.load(str(args.data / "figure_3ef_scvi_model.pt"), adata=adata)
        adata.obsm["X_scVI"] = model.get_latent_representation()

        sc.pp.neighbors(adata, use_rep="X_scVI", **NEIGHBORS)
        sc.tl.umap(adata, n_components=2)
        coordinates = adata.obsm["X_umap"][:, ::-1]  # the analysis swaps the two axes

        frame = pd.DataFrame(
            coordinates, index=adata.obs_names, columns=["scvi_umap_1", "scvi_umap_2"]
        )
        frame.index.name = "barcode"

        if out.exists():
            shipped = pd.read_csv(out, index_col="barcode")
            difference = float(np.abs(shipped.to_numpy() - frame.to_numpy()).max())
            print(f"largest difference from the shipped embedding: {difference:.3g}")

        frame.to_csv(out)
        print(f"{out.name}: {len(frame)} cells")

    manifest = refresh_manifest(args.data, {"figure_3ef": manifest_section(args.data)})
    print(
        f"MANIFEST.json: {len(manifest['files'])} files, "
        f"{manifest['total_bytes'] / 1e6:.1f} MB"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
