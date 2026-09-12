#!/usr/bin/env python
"""Regenerate the Figure 3E / 3F embedding from the deposited scVI checkpoint.

    python scripts/integrate_s1_s2.py            # -> data/figure_3ef_embedding.csv.gz

Reads the 6,881 x 8,778 count matrix the model was trained on, runs it through
the deposited checkpoint to get the 10-dimensional latent space, and lays that
out with the same neighbour graph and UMAP call the analysis used.

Rendering the panels does not need this script — the embedding it writes is
already in ``data/``. It needs the ``integration`` extra
(``pip install -e '.[integration]'``), which pulls in torch and scvi-tools.

Loading the checkpoint is deterministic, and so is the UMAP that follows it, so
re-running this reproduces the shipped file. Training a fresh model does not:
the layout depends on the particular fit, and the published panels were drawn
from a different one. See docs/reproducibility.md.
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

NEIGHBORS = {"n_neighbors": 50, "n_pcs": 10}


def build(datadir: Path):
    import anndata
    import scanpy as sc

    obs = pd.read_csv(datadir / "figure_3ef_obs.csv.gz", index_col="barcode", dtype="str")
    var = pd.read_csv(datadir / "figure_3ef_var.csv.gz", index_col="gene_id", dtype="str")
    with gzip.open(datadir / "figure_3ef_counts.mtx.gz", "rb") as handle:
        counts = sp.csr_matrix(scipy.io.mmread(handle)).astype(np.float32)
    if counts.shape != (len(obs), len(var)):
        raise SystemExit(f"counts shape {counts.shape} does not match {(len(obs), len(var))}")

    adata = anndata.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    return adata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=REPO / "data")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    out = args.out or args.data / "figure_3ef_embedding.csv.gz"

    import scanpy as sc
    import scvi

    adata = build(args.data)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")
    model = scvi.model.SCVI.load(str(args.data / "figure_3ef_scvi_model.pt"), adata=adata)
    adata.obsm["X_scVI"] = model.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI", **NEIGHBORS)
    sc.tl.umap(adata, n_components=2)
    coordinates = adata.obsm["X_umap"][:, ::-1]      # the analysis swaps the two axes

    frame = pd.DataFrame(coordinates, index=adata.obs_names,
                         columns=["scvi_umap_1", "scvi_umap_2"])
    frame.index.name = "barcode"

    if out.exists():
        shipped = pd.read_csv(out, index_col="barcode")
        difference = float(np.abs(shipped.to_numpy() - frame.to_numpy()).max())
        print(f"largest difference from the shipped embedding: {difference:.3g}")

    frame.to_csv(out)
    print(f"{out.relative_to(REPO) if out.is_relative_to(REPO) else out}: {len(frame)} cells")
    return 0


if __name__ == "__main__":
    sys.exit(main())
