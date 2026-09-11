#!/usr/bin/env python
"""Render every reproducible figure panel to PNG, SVG and PDF.

    python scripts/make_figures.py                 # everything, into figures/
    python scripts/make_figures.py --panel 1C 1G   # just those
    python scripts/make_figures.py --format png    # one format only

Run ``scripts/fetch_data.py`` first; this reads the deposited dataset from
``data/`` and writes nothing else.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from bzfig import panels  # noqa: E402
from bzfig.data import load_dataset, load_supp1a_markers  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
FORMATS = ("png", "svg", "pdf")
DPI = 300

# Two things make a render differ from the one before it whatever the panel
# draws: matplotlib stamps SVG and PDF with the time of the run, and salts the
# element ids it puts in SVGs with a per-process random value. Dropping the
# stamp and pinning the salt makes all three formats byte-reproducible.
matplotlib.rcParams["svg.hashsalt"] = "bzfig"
METADATA = {"pdf": {"CreationDate": None}, "svg": {"Date": None}}


def save(figure, outdir: Path, name: str, formats) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for extension in formats:
        path = outdir / f"{name}.{extension}"
        figure.savefig(
            path,
            dpi=DPI,
            bbox_inches="tight",
            facecolor=figure.get_facecolor(),
            metadata=METADATA.get(extension),
        )
        written.append(path)
    plt.close(figure)
    return written


def build(adata, markers: pd.Index, datadir: Path) -> dict[str, callable]:
    """Panel name -> zero-argument callable returning a figure.

    The volcano panels need *datadir* as well: they widen the matrix back to the
    8322-gene universe with the extra genes shipped beside it.
    """
    supp4 = panels.supplementary_4_matrices(adata)
    jobs = {
        "Figure_1B_umap": lambda: panels.figure_1b_umap(adata),
        "Figure_1B_cells_per_cluster": lambda: panels.figure_1b_counts(adata),
        "Figure_1C_marker_heatmap": lambda: panels.figure_1c(adata),
        "Figure_1D_unique_markers_per_cluster": lambda: panels.figure_1d(),
        "Figure_1F_cst1_srs44": lambda: panels.figure_1f(adata),
        "Figure_1G_cst1_violin": lambda: panels.figure_1g(adata),
        "Supplementary_1A_all_markers_heatmap": lambda: panels.supplementary_1a(adata, markers),
        "Supplementary_4_CCC": lambda: panels.supplementary_4(
            adata, supp4["CCC"], "Common cell cycle (CCC) expression"
        ),
        "Supplementary_4_MCC": lambda: panels.supplementary_4(
            adata, supp4["MCC"], "Modified cell cycle (MCC) expression"
        ),
        "Supplementary_5A_umap": lambda: panels.supplementary_5a_umap(adata),
        "Supplementary_5A_cells_per_cluster": lambda: panels.supplementary_5a_counts(adata),
        "Supplementary_5B_umap": lambda: panels.supplementary_5b_umap(adata),
        "Supplementary_5B_cells_per_cluster": lambda: panels.supplementary_5b_counts(adata),
        "Supplementary_5C_volcano": lambda: panels.supplementary_5c(adata, datadir),
        "Supplementary_5D_srs22a": lambda: panels.supplementary_5d(adata),
        "Supplementary_5E_volcano": lambda: panels.supplementary_5e(adata, datadir),
        "Supplementary_5F_volcano": lambda: panels.supplementary_5f(adata, datadir),
    }
    for group, label, gene_id in panels.FIG1E_GENES:
        safe = label.replace(" ", "_")
        jobs[f"Figure_1E_{safe}"] = (
            lambda g=gene_id, lbl=label: panels._gene_umap(panels.in_vivo(adata), g, lbl)
        )
    return jobs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=REPO / "data")
    parser.add_argument("--outdir", type=Path, default=REPO / "figures")
    parser.add_argument("--format", nargs="+", choices=FORMATS, default=list(FORMATS))
    parser.add_argument("--panel", nargs="+", default=None, help="substring match on panel names")
    parser.add_argument("--list", action="store_true", help="list panel names and exit")
    args = parser.parse_args()

    adata = load_dataset(args.data)
    markers = load_supp1a_markers(args.data)
    jobs = build(adata, markers, args.data)

    if args.list:
        for name in jobs:
            print(name)
        return 0

    if args.panel:
        wanted = {n for n in jobs for p in args.panel if p.lower() in n.lower()}
        if not wanted:
            print(f"No panel matched {args.panel}; use --list to see the names")
            return 1
        jobs = {n: jobs[n] for n in jobs if n in wanted}

    print(f"Rendering {len(jobs)} panel(s) to {args.outdir}")
    for name, make in jobs.items():
        written = save(make(), args.outdir, name, args.format)
        print(f"  {name}: {', '.join(p.suffix.lstrip('.') for p in written)}")
    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
