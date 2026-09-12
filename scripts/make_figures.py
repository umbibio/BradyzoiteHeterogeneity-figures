#!/usr/bin/env python
"""Render every reproducible figure panel to PNG, SVG and PDF.

    python scripts/make_figures.py                 # everything, into figures/
    python scripts/make_figures.py --panel 1C 1G   # just those
    python scripts/make_figures.py --format png    # one format only

The deposited dataset is fetched automatically if ``data/`` is empty, stale, or
holds git-lfs pointers rather than the files themselves — so a clone made without
git-lfs installed still works. ``--no-fetch`` turns that off.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import fetch_data  # noqa: E402

from bzfig import figure_2h_supplementary_4 as heatmaps  # noqa: E402
from bzfig import figure_3ef  # noqa: E402
from bzfig import panels  # noqa: E402
from bzfig.constants import COHORTS, SUPP1_GENES  # noqa: E402
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


def build(
    adata, markers: pd.Index, datadir: Path, outdir: Path | None = None
) -> dict[str, callable]:
    """Panel name -> zero-argument callable returning a figure.

    The volcano panels need *datadir* as well: they widen the matrix back to the
    8322-gene universe with the extra genes shipped beside it.

    Nothing here touches *adata*, so the registry can be built — and listed —
    without a dataset. The Figure 2H and Supplementary 4 heatmaps additionally
    defer their own (expensive) normalisation until a caller asks for them, and
    write their companion CSV tables under ``outdir/tables`` when *outdir* is
    given.
    """
    # Lazy: listing jobs never reads data, and selecting an existing panel does
    # not compute the Figure 2H / Supplementary 4 heatmaps (or vice versa).

    # The Figure 3E/3F embedding is a small table of its own, read once and
    # shared by the two panels rather than loaded with the main dataset.
    @lru_cache(maxsize=1)
    def integration():
        return figure_3ef.load_integration(datadir)

    # lru_cache so the CCC and MCC panels share one pass over the matrix.
    @lru_cache(maxsize=1)
    def supp4():
        tables = heatmaps.phase_matrices(adata)
        if outdir is not None:
            heatmaps.write_tables(outdir / "tables", phase=tables)
        return tables

    def fig2h():
        table = heatmaps.correlation_matrix(adata)
        if outdir is not None:
            heatmaps.write_tables(outdir / "tables", correlation=table)
        return heatmaps.figure_2h(table)

    jobs = {
        "Figure_1B_umap": lambda: panels.figure_1b_umap(adata),
        "Figure_1B_cells_per_cluster": lambda: panels.figure_1b_counts(adata),
        "Figure_1C_marker_heatmap": lambda: panels.figure_1c(adata),
        "Figure_1D_unique_markers_per_cluster": lambda: panels.figure_1d(),
        "Figure_1F_cst1_srs44": lambda: panels.figure_1f(adata),
        "Figure_1G_cst1_violin": lambda: panels.figure_1g(adata),
        "Figure_2H_correlation_heatmap": fig2h,
        "Figure_3D_in_vitro": lambda: panels.figure_3d(adata),
        "Figure_3E_tachyzoite_highlight": lambda: figure_3ef.figure_3e(integration(), datadir),
        "Figure_3F_tachyzoite_cc_phase": lambda: figure_3ef.figure_3f(integration(), datadir),
        "Supplementary_1A_all_markers_heatmap": lambda: panels.supplementary_1a(adata, markers),
        "Supplementary_4_CCC": lambda: heatmaps.supplementary_4(supp4()["CCC"], "CCC"),
        "Supplementary_4_MCC": lambda: heatmaps.supplementary_4(supp4()["MCC"], "MCC"),
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
    for panel, label, gene_id in SUPP1_GENES:
        jobs[f"Supplementary_{panel}_{label}"] = (
            lambda g=gene_id, lbl=label: panels._gene_umap(panels.in_vivo(adata), g, lbl)
        )
    # Figure 6 redraws Figure 3A-C without the grey layer. Its panels are named
    # for the cohort they hold, not for the column header printed above them,
    # which the published figure does not match to the cells beneath it.
    for letter, name, orig_ident, background in COHORTS:
        jobs[f"Figure_{letter}_umap"] = (
            lambda o=orig_ident, bg=background: panels.figure_3_umap(adata, o, bg)
        )
        jobs[f"Figure_{letter}_cells_per_phase"] = (
            lambda o=orig_ident: panels.figure_3_counts(adata, o)
        )
        jobs[f"Figure_6_{name}"] = lambda o=orig_ident: panels.figure_6_umap(adata, o)
    return jobs


def ensure_data(datadir: Path, no_fetch: bool) -> bool:
    """Make sure ``datadir`` holds the real input files before rendering.

    A clone made without git-lfs leaves pointer stubs in place of the data, which
    would otherwise fail deep inside the loader with something unhelpful. The
    manifest is in plain git precisely so it can be read at this point.
    """
    manifest = fetch_data.read_manifest_file(datadir)
    bad = [e for e in manifest["files"] if not fetch_data.verify(datadir / e["filename"], e)[0]]
    if not bad:
        return True

    names = ", ".join(e["filename"] for e in bad)
    if no_fetch:
        print(f"{len(bad)} input file(s) missing or not intact: {names}")
        print("Run scripts/fetch_data.py, or drop --no-fetch.")
        return False

    print(f"{len(bad)} input file(s) missing or not intact; fetching.")
    bases = fetch_data.base_urls([], manifest)
    missing = [e["filename"] for e in bad if not fetch_data.fetch(e, datadir, bases, force=True)]
    if missing:
        print(f"\nCould not retrieve: {', '.join(missing)}")
        print("Supply a mirror with BZFIG_DATA_URLS, or see data/README.md")
        return False
    print()
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=REPO / "data")
    parser.add_argument("--outdir", type=Path, default=REPO / "figures")
    parser.add_argument("--format", nargs="+", choices=FORMATS, default=list(FORMATS))
    parser.add_argument("--panel", nargs="+", default=None, help="substring match on panel names")
    parser.add_argument("--list", action="store_true", help="list panel names and exit")
    parser.add_argument(
        "--no-fetch", action="store_true", help="fail instead of downloading missing input"
    )
    args = parser.parse_args()

    # Listing is a registry question, not a data question, so it works in a
    # clone that has not fetched the dataset yet.
    if args.list:
        for name in build(None, None, args.data):
            print(name)
        return 0

    if not ensure_data(args.data, args.no_fetch):
        return 1

    adata = load_dataset(args.data)
    markers = load_supp1a_markers(args.data)
    jobs = build(adata, markers, args.data, args.outdir)

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
