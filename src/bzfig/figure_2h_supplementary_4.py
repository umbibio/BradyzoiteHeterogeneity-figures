"""Figure 2H and Supplementary 4 heatmaps, using the shared figure dataset.

The original R scripts used Seurat LogNormalize after the 8,170-gene subset.
The deposited logcounts was normalized before that subset. Closing expm1(logcounts)
to 10,000 over the SAME 8,170 genes recovers the R values to source float32
precision. This panel-specific conversion does not modify adata or other plots.
No additional expression dataset, Seurat installation, or new analysis is used.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from matplotlib.transforms import offset_copy
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .constants import PHASES, RECOVERED_FLAG

METADATA = json.loads(Path(__file__).with_name("figure_2h_supplementary_4_metadata.json").read_text())
PHASE_COLORS = ["#3375a7", "#da8238", "#369539", "#bf3c3e", "#9972b3"]


def normalized_expression(adata):
    """Recover the original R heatmap normalization from shared logcounts.

    Normalize before selecting plotted genes. Zero-total/nonfinite/negative
    inputs are errors, not silently converted to a plausible-looking heatmap.

    The cell total is closed over the analysis object's 8,170 genes, which is
    what the original recipe divided by. The deposited matrix is wider than
    that -- it carries the 152 genes recovered for the volcano panels -- so
    those columns are excluded from the denominator.
    """
    kept = adata.var[RECOVERED_FLAG].to_numpy(dtype=bool)
    if int(kept.sum()) != 8170 or not adata.var_names.is_unique:
        raise ValueError("These panels require the analysis object's 8,170 genes.")
    matrix = sp.csr_matrix(adata.layers["logcounts"], dtype=np.float64, copy=True)
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError("Expected finite nonnegative logcounts.")
    matrix.data = np.expm1(matrix.data)
    totals = np.asarray(matrix[:, kept].sum(axis=1)).ravel()
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise ValueError("Every cell must have a finite positive total.")
    matrix = (sp.diags(10000.0 / totals) @ matrix).tocsr()
    matrix.data = np.log1p(matrix.data)
    return matrix


def selected_cells(adata):
    obs = adata.obs
    selected = ((obs["orig_ident"] == "Nonreactivated")
                & obs["cell_cycle_group"].isin(["CCC", "MCC"])
                & obs["transferred_cc_phase"].isin(PHASES)).to_numpy()
    if int(selected.sum()) != 6505 or not obs.index.is_unique:
        raise ValueError("Expected the 6,505 unique original in-vivo cells with cycle and phase labels.")
    return selected


def positions(adata, genes):
    result = adata.var_names.get_indexer(genes)
    if (result < 0).any():
        raise ValueError(f"Missing required genes: {[g for g, i in zip(genes, result) if i < 0]}")
    return result


def phase_matrices(adata):
    """Same 50 genes, original independent row orders, mean then sample-SD z-score."""
    expression = normalized_expression(adata)
    selected = selected_cells(adata)
    result = {}
    for cycle in ("CCC", "MCC"):
        genes = METADATA["supp4"][cycle]["gene_ids"]
        values = expression[:, positions(adata, genes)].toarray()
        pick = selected & (adata.obs["cell_cycle_group"] == cycle).to_numpy()
        groups = [pick & (adata.obs["transferred_cc_phase"] == p).to_numpy() for p in PHASES]
        counts = [int(g.sum()) for g in groups]
        if counts != METADATA["supp4"][cycle]["phase_counts"]:
            raise ValueError(f"{cycle}: phase cell counts do not match the original panel.")
        mean = np.stack([values[g].mean(axis=0) for g in groups], axis=1)
        sd = mean.std(axis=1, ddof=1, keepdims=True)
        if (sd <= 0).any():
            raise ValueError(f"{cycle}: a required gene has constant phase means.")
        result[cycle] = {
            "mean": pd.DataFrame(mean, index=genes, columns=PHASES),
            "z": pd.DataFrame((mean-mean.mean(axis=1, keepdims=True))/sd, index=genes, columns=PHASES),
            "n_cells": pd.Series(counts, index=PHASES),
        }
    return result


def correlation_matrix(adata):
    """Original 5-by-29 Pearson correlations across cells, not phase means."""
    genes = METADATA["fig2H"]["gene_ids"]
    values = normalized_expression(adata)[selected_cells(adata), :][:, positions(adata, genes)].toarray()
    if not np.isfinite(values).all() or (values.std(axis=0) == 0).any():
        raise ValueError("Correlation input contains nonfinite or constant gene values.")
    matrix = np.corrcoef(values, rowvar=False)[:5, :]
    return pd.DataFrame(matrix, index=genes[:5], columns=genes)


def palette():
    """Fixed LAB-interpolated colors from original circlize::colorRamp2.

    The lookup records a rendering choice, not measured expression. No R
    runtime is required. Breaks are [-2.5,0,2.5] for Supp4 and [-1,0,1] for 2H.
    """
    return ListedColormap(METADATA["color_lut"], name="original_R_blue_white_red")


def figure_2h(matrix):
    labels = METADATA["fig2H"]["labels"]
    fig, ax = plt.subplots(figsize=(12, 4), layout="constrained")
    image = ax.imshow(matrix.to_numpy(), cmap=palette(), vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=9)
    ax.set_yticks(range(5), labels[:5], fontsize=10)
    ax.tick_params(length=0)
    for (row, col), value in np.ndenumerate(matrix.to_numpy()):
        ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("Figure 2H | Pearson correlation across 6,505 cells", fontsize=11)
    fig.colorbar(image, ax=ax, shrink=0.7, label="Pearson r")
    return fig


def supplementary_4(panel, cycle):
    """Draw original-order data and copy the source's explicit label highlights."""
    description = METADATA["supp4"][cycle]
    fig, ax = plt.subplots(figsize=(8.5, 12))
    fig.subplots_adjust(left=0.04, right=0.34, top=0.92, bottom=0.10)
    z = panel["z"].to_numpy()
    image = ax.imshow(z, cmap=palette(), vmin=-2.5, vmax=2.5, aspect="auto")
    labels = [g.replace('_', '-') for g in description["gene_ids"]]
    ax.set_yticks(range(50), labels, fontsize=7)
    ax.yaxis.tick_right()
    ax.set_xticks(range(5), PHASES, fontsize=10)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    # The source highlights descriptions only, leaving gene IDs and text black.
    # Measure the ID column so description strips align without touching the data.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ticks = ax.get_yticklabels()
    id_width = max(t.get_window_extent(renderer).x1 - ax.bbox.x1 for t in ticks)
    offset_points = id_width * 72 / fig.dpi + 3
    transform = offset_copy(ax.get_yaxis_transform(), fig=fig,
                            x=offset_points, units="points")
    font = ticks[0].get_fontproperties()
    text_width = max(renderer.get_text_width_height_descent(d, font, False)[0]
                     for d in description["labels"])
    strip_width = (text_width + 3 * fig.dpi / 72) / ax.bbox.width
    annotation = METADATA["supp4_highlights"]
    highlighted = set(annotation[cycle])
    if not highlighted.issubset(description["gene_ids"]):
        raise ValueError(f"{cycle}: highlighted gene absent from the displayed panel.")
    for row, (gene, label) in enumerate(zip(description["gene_ids"], description["labels"])):
        if gene in highlighted:
            patch = Rectangle((1, row - 0.5), strip_width, 1, transform=transform,
                              facecolor=annotation["background_color"], edgecolor="none",
                              clip_on=False, zorder=0)
            patch.set_gid(f"manual-highlight-{gene}")
            ax.add_patch(patch)
        ax.text(1, row, label, transform=transform, va="center", ha="left",
                fontsize=7, color="black", clip_on=False, zorder=3)
    strip = ax.inset_axes([0, 1.035, 1, 0.012])
    strip.imshow([range(5)], cmap=ListedColormap(PHASE_COLORS), vmin=-0.5, vmax=4.5, aspect="auto")
    strip.set_axis_off()
    strip.set_title(f"{'Common' if cycle == 'CCC' else 'Modified'} cell cycle ({cycle}) expression", fontsize=11, loc="left", pad=12)
    fig.colorbar(image, cax=fig.add_axes([0.04, 0.066, 0.3, 0.012]), orientation="horizontal", label="Within-cycle gene z-score", ticks=[-2.5, 0, 2.5])
    fig.text(0.04, 0.008, "Original R normalization and row order; pink annotations copied from the manuscript.", fontsize=8)
    return fig


def write_tables(outdir, *, correlation=None, phase=None):
    """Companion computed values, regenerated from the shared dataset on every run."""
    outdir.mkdir(parents=True, exist_ok=True)
    if correlation is not None:
        correlation.to_csv(outdir / "Figure_2H_correlations.csv", index_label="gene_id", float_format="%.17g")
    if phase is not None:
        for cycle, panel in phase.items():
            for name in ("mean", "z"):
                panel[name].to_csv(outdir / f"Supplementary_4_{cycle}_{name}.csv", index_label="gene_id", float_format="%.17g")
            panel["n_cells"].to_csv(outdir / f"Supplementary_4_{cycle}_cell_counts.csv", index_label="phase", header=["n_cells"])
