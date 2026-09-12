"""Figure 3E and 3F — the in vivo tachyzoites on the scVI integration.

These two panels come from a second experiment: 376 in vivo tachyzoites taken
from the peritoneal cavity at 5 dpi (samples S1 and S2), integrated with the
6,505 in vivo bradyzoites by scVI and plotted on a UMAP of the 10-dimensional
latent space.

The package ships the integration as four files: the 6,881 x 8,778 count matrix
the model was trained on, its cell and gene tables, and the trained scVI
checkpoint. ``data/figure_3ef_embedding.csv.gz`` holds the UMAP coordinates the
panels are drawn from, so rendering needs neither torch nor scvi-tools;
``scripts/integrate_s1_s2.py`` regenerates that file from the checkpoint.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from .constants import BACKGROUND_GRAY, FIG3EF_TACHYZOITE, PHASES
from .data import DATA

TACHYZOITE_SAMPLES = ("S1", "S2")


def phase_palette(datadir: Path = DATA) -> list[str]:
    """The cell-cycle palette that travels with the data, as Figure 3A-C uses."""
    colors = json.loads((Path(datadir) / "uns_colors.json").read_text())
    return list(colors["cc_phase_colors"])


def load_integration(datadir: Path = DATA) -> pd.DataFrame:
    """Cell table and UMAP coordinates for the 6,881 integrated cells."""
    datadir = Path(datadir)
    obs = pd.read_csv(datadir / "figure_3ef_obs.csv.gz", index_col="barcode", dtype="str")
    embedding = pd.read_csv(datadir / "figure_3ef_embedding.csv.gz", index_col="barcode")
    if not obs.index.equals(embedding.index):
        raise ValueError("figure_3ef_embedding.csv.gz: cell order does not match the obs table")
    frame = obs.join(embedding)
    frame["tachyzoite"] = frame["sample"].isin(TACHYZOITE_SAMPLES)
    return frame


def _axes(frame: pd.DataFrame):
    figure, ax = plt.subplots(figsize=(4.6, 3.45), dpi=300)
    ax.set_aspect("equal")
    ax.set_axis_off()
    x, y = frame["scvi_umap_1"], frame["scvi_umap_2"]
    pad = 0.03 * max(x.max() - x.min(), y.max() - y.min())
    ax.set_xlim(x.min() - pad, x.max() + pad)
    ax.set_ylim(y.min() - pad, y.max() + pad)
    return figure, ax


def _scatter(ax, frame, color, size=6.0, zorder=1):
    collection = ax.scatter(
        frame["scvi_umap_1"], frame["scvi_umap_2"],
        s=size, c=color, linewidths=0, zorder=zorder,
    )
    collection.set_rasterized(True)


def figure_3e(frame: pd.DataFrame | None = None, datadir: Path = DATA):
    """The 376 in vivo tachyzoites picked out of the integration."""
    frame = load_integration(datadir) if frame is None else frame
    figure, ax = _axes(frame)
    _scatter(ax, frame[~frame["tachyzoite"]], BACKGROUND_GRAY)
    _scatter(ax, frame[frame["tachyzoite"]], FIG3EF_TACHYZOITE, zorder=2)
    ax.legend(
        handles=[
            plt.Line2D([], [], marker="o", linestyle="", markersize=4,
                       color=FIG3EF_TACHYZOITE, label="in vivo tachyzoite"),
            plt.Line2D([], [], marker="o", linestyle="", markersize=4,
                       color=BACKGROUND_GRAY, label="in vivo bradyzoite"),
        ],
        frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8,
    )
    figure.tight_layout()
    return figure


def figure_3f(frame: pd.DataFrame | None = None, datadir: Path = DATA):
    """The same tachyzoites, coloured by transferred cell-cycle phase."""
    frame = load_integration(datadir) if frame is None else frame
    figure, ax = _axes(frame)
    _scatter(ax, frame[~frame["tachyzoite"]], BACKGROUND_GRAY)
    tachyzoites = frame[frame["tachyzoite"]]
    handles = []
    for phase, color in zip(PHASES, phase_palette(datadir)):
        cells = tachyzoites[tachyzoites["cc_phase"] == phase]
        _scatter(ax, cells, color, zorder=2)
        handles.append(
            plt.Line2D([], [], marker="o", linestyle="", markersize=4, color=color, label=phase)
        )
    ax.legend(handles=handles, frameon=False, loc="center left",
              bbox_to_anchor=(1.0, 0.5), fontsize=8, title="cc_phase",
              title_fontproperties={"size": 8})
    figure.tight_layout()
    return figure


def phase_counts(frame: pd.DataFrame | None = None, datadir: Path = DATA) -> pd.Series:
    """Cells per phase among the in vivo tachyzoites."""
    frame = load_integration(datadir) if frame is None else frame
    tachyzoites = frame[frame["tachyzoite"]]
    return tachyzoites["cc_phase"].value_counts().reindex(PHASES).fillna(0).astype(int)
