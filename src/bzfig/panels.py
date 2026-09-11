"""One function per published panel.

Each returns a matplotlib figure drawn from the shipped data alone — for the
Supplementary 5 volcanoes that includes the extra genes the deposited object
dropped, which :mod:`bzfig.de` folds back in. Panel assembly — the grid layout,
the row labels beside Figure 1E, the letter labels, the in-plot cluster numbers —
was done in a vector editor and is not reproduced here; these are the individual
panels as the analysis produced them.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from . import de
from .constants import (
    BACKGROUND_GRAY,
    CLUSTER_COLORS,
    EXPRESSION_CMAP,
    EXPRESSION_VMAX,
    EXPRESSION_VMIN,
    FIG1C_GENES,
    FIG1E_GENES,
    FIG1F_GENE,
    FIG3D_IN_VITRO,
    HEATMAP_VMAX,
    PHASES,
    SUPP1_GENES,
    SUPP4_GENES,
    SUPP5D_GENE,
    SUPP5E_CYST_WALL_GENES,
    UNIQUE_MARKERS_PER_CLUSTER,
    VOLCANO_AXES,
    VOLCANO_COLORS,
    VOLCANO_LABEL_RANGE,
)
from .data import DATA
from .scatter3d import AZIMUTH, ELEVATION, plot_3d_preview

UMAP_KEY = "3d_umap_harmony_integration"
EXPRESSION_LAYER = "logcounts_scaled"
CELLCYCLE_LAYER = "logcounts"


# --------------------------------------------------------------------- subsets


def in_vivo(adata):
    """The non-reactivated in vivo bradyzoites — 6505 cells, clusters 0-5."""
    return adata[adata.obs["orig_ident"] == "Nonreactivated"]


def in_vitro_bradyzoites(adata):
    """The in vitro bradyzoites (me49 Day 3) — 950 cells, transferred identities."""
    return adata[adata.obs["orig_ident"] == "me49 Day 3"]


# ------------------------------------------------------------------- utilities


def _categorical_colors(values: pd.Series) -> pd.Series:
    """Map a cluster column onto the published palette, unassigned cells grey."""
    categories = values.cat.categories if hasattr(values, "cat") else sorted(values.dropna().unique())
    lookup = {category: CLUSTER_COLORS[i] for i, category in enumerate(categories)}
    return values.map(lookup).astype("string").fillna("lightgray")


def _cluster_umap(subset, colorby: str, legend: bool = False):
    points = subset.obsm[UMAP_KEY]
    colors = _categorical_colors(subset.obs[colorby])
    ax = plot_3d_preview(
        points, colors, elevation=ELEVATION, angle=AZIMUTH, remove_panels=True
    )
    if legend:
        categories = subset.obs[colorby].cat.categories
        handles = [
            plt.Line2D([], [], marker="o", linestyle="", color=CLUSTER_COLORS[i], label=str(c))
            for i, c in enumerate(categories)
        ]
        ax.legend(handles=handles, loc="upper right", frameon=False, title="Cluster")
    return ax.figure


def _count_bars(counts: pd.Series, xlabel: str = "Count", ylabel: str = "Cluster", palette=None):
    frame = counts.rename_axis(ylabel).rename(xlabel).reset_index()
    frame[ylabel] = frame[ylabel].astype(str)
    figure = plt.figure(figsize=(10, 2))
    sns.barplot(
        frame, x=xlabel, y=ylabel, hue=ylabel, dodge=False,
        palette=palette or CLUSTER_COLORS, legend=False,
    )
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
    return figure


def _gene_umap(subset, gene_id: str, label: str):
    points = subset.obsm[UMAP_KEY]
    values = subset[:, gene_id].layers[EXPRESSION_LAYER]
    values = np.asarray(values.todense()).ravel() if hasattr(values, "todense") else np.asarray(values).ravel()
    ax = plot_3d_preview(
        points,
        values,
        elevation=ELEVATION,
        angle=AZIMUTH,
        cmap=EXPRESSION_CMAP,
        vmin=EXPRESSION_VMIN,
        vmax=EXPRESSION_VMAX,
        remove_panels=True,
        show_colorscale=True,
    )
    ax.set_title(label, loc="center", style="italic", fontsize=16)
    return ax.figure


def _marker_heatmap(subset, gene_ids, labels, swap_axes: bool):
    plot = sc.pl.heatmap(
        subset.copy(),
        list(gene_ids),
        "nr_cluster",
        layer=EXPRESSION_LAYER,
        swap_axes=swap_axes,
        vmax=HEATMAP_VMAX,
        show=False,
    )
    if labels is not None:
        plot["heatmap_ax"].set_yticklabels(labels)
    plot["groupby_ax"].set_xlabel("Cluster") if swap_axes else plot["groupby_ax"].set_ylabel("Cluster")
    return plt.gcf()


# ---------------------------------------------------------------------- Fig 1B


def figure_1b_umap(adata):
    """Cluster UMAP of the in vivo bradyzoites."""
    return _cluster_umap(in_vivo(adata), "nr_cluster")


def figure_1b_counts(adata):
    """Cells per cluster."""
    counts = in_vivo(adata).obs["nr_cluster"].value_counts().sort_index()
    return _count_bars(counts)


# ---------------------------------------------------------------------- Fig 1C


def figure_1c(adata):
    """Selected marker genes per cluster."""
    _, descriptions, gene_ids = map(list, zip(*FIG1C_GENES))
    labels = [f"{d} - {g}" for _, d, g in FIG1C_GENES]
    figure = _marker_heatmap(in_vivo(adata), gene_ids, labels, swap_axes=True)
    plt.subplots_adjust(left=0.4, right=0.95, top=1, bottom=0.1)
    return figure


# ---------------------------------------------------------------------- Fig 1D


def figure_1d(adata=None):
    """Unique markers per cluster — a recorded constant, see constants.py."""
    counts = pd.Series(UNIQUE_MARKERS_PER_CLUSTER).sort_index()
    return _count_bars(counts)


# ------------------------------------------------------------------ Fig 1E, 1F


def figure_1e(adata):
    """Six per-gene expression UMAPs — one figure per gene, keyed by label."""
    subset = in_vivo(adata)
    return {label: _gene_umap(subset, gene_id, label) for _, label, gene_id in FIG1E_GENES}


def figure_1f(adata):
    """cst1 (srs44) expression."""
    label, gene_id = FIG1F_GENE
    return _gene_umap(in_vivo(adata), gene_id, label)


# ---------------------------------------------------------------------- Fig 1G


def figure_1g(adata, seed: int = 0):
    """CST1/SRS44 expression per cluster.

    The overlaid points are jittered, which draws on the global numpy RNG; the
    seed keeps successive renders byte-identical.
    """
    subset = in_vivo(adata)
    ordered = subset.obs.sort_values("nr_cluster", kind="stable").index
    _, gene_id = FIG1F_GENE
    np.random.seed(seed)
    ax = sc.pl.violin(
        subset[ordered].copy(),
        gene_id,
        groupby="nr_cluster",
        layer=EXPRESSION_LAYER,
        show=False,
        dodge=False,
        hue="nr_cluster",
        palette=CLUSTER_COLORS,
    )
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Normalized Expression")
    ax.set_title(f"CST1/SRS44: {gene_id}")
    for collection in ax.collections:
        collection.set_rasterized(True)
    return ax.figure


# ---------------------------------------------------------------- Fig 3, Fig 6


def cohort(adata, orig_ident: str):
    """Boolean mask over all 8057 cells for one orig_ident cohort."""
    return (adata.obs["orig_ident"] == orig_ident).to_numpy()


def _highlight_umap(adata, colors, mask=None):
    """The projection drawn once, coloured per cell; *mask* drops the rest."""
    points = adata.obsm[UMAP_KEY]
    if mask is not None:
        points, colors = points[mask], [c for c, keep in zip(colors, mask) if keep]
    return plot_3d_preview(
        points, colors, elevation=ELEVATION, angle=AZIMUTH, remove_panels=True
    ).figure


def _phase_colors(adata, mask):
    """Phase palette for the masked cells, grey for everything else."""
    palette = dict(zip(PHASES, adata.uns["cc_phase_colors"]))
    phases = adata.obs["transferred_cc_phase"].astype(str)
    return [palette[p] if keep else BACKGROUND_GRAY for p, keep in zip(phases, mask)]


def figure_3_umap(adata, orig_ident: str, background: bool = True):
    """One cohort by transferred cell-cycle phase, the rest of the cells grey.

    Figure 3A is drawn without the grey layer — the published panel is the same
    rendering as Figure 6's first one — so *background* is false for it.
    """
    mask = cohort(adata, orig_ident)
    colors = _phase_colors(adata, mask)
    return _highlight_umap(adata, colors, mask=None if background else mask)


def figure_3_counts(adata, orig_ident: str):
    """Cells per phase in that cohort."""
    mask = cohort(adata, orig_ident)
    counts = adata.obs["transferred_cc_phase"][mask].value_counts().reindex(PHASES)
    return _count_bars(counts, ylabel="CC Phase", palette=adata.uns["cc_phase_colors"])


def figure_3d(adata):
    """The in vitro cells picked out of the projection in flat salmon."""
    in_vitro = (adata.obs["dataset_type"] == "inVitro").to_numpy()
    colors = [FIG3D_IN_VITRO if keep else BACKGROUND_GRAY for keep in in_vitro]
    return _highlight_umap(adata, colors)


def figure_6_umap(adata, orig_ident: str):
    """As Figure 3A-C, with the grey layer dropped rather than drawn."""
    mask = cohort(adata, orig_ident)
    return _highlight_umap(adata, _phase_colors(adata, mask), mask=mask)


# ------------------------------------------------------------- Supp 1A, 1B, 1C


def supplementary_1a(adata, marker_genes: pd.Index):
    """All cluster markers per cluster."""
    subset = in_vivo(adata)
    genes = marker_genes.intersection(subset.var_names)
    return _marker_heatmap(subset, genes, None, swap_axes=False)


def supplementary_1bc(adata):
    """Eleven per-gene expression UMAPs — as Figure 1E, keyed by panel and label."""
    subset = in_vivo(adata)
    return {
        (panel, label): _gene_umap(subset, gene_id, label)
        for panel, label, gene_id in SUPP1_GENES
    }


# -------------------------------------------------------------------- Supp 4


def supplementary_4_matrices(adata) -> dict[str, dict[str, pd.DataFrame]]:
    """Per-phase mean expression and per-gene z-score, for the CCC and MCC panels."""
    panels = {}
    for group in ("CCC", "MCC"):
        subset = adata[adata.obs["cell_cycle_group"] == group, SUPP4_GENES]
        values = subset.layers[CELLCYCLE_LAYER]
        frame = pd.DataFrame(
            np.asarray(values.todense()) if hasattr(values, "todense") else np.asarray(values),
            columns=SUPP4_GENES,
            index=subset.obs_names,
        )
        frame["phase"] = subset.obs["transferred_cc_phase"].astype(str).values
        mean = frame.groupby("phase").mean().reindex(PHASES).T
        z = mean.sub(mean.mean(axis=1), axis=0).div(mean.std(axis=1, ddof=1), axis=0)
        order = (
            pd.DataFrame({"phase": z.idxmax(axis=1).map(PHASES.index), "peak": z.max(axis=1)})
            .sort_values(["phase", "peak"], ascending=[True, False])
            .index
        )
        panels[group] = {
            "mean": mean.loc[order],
            "z": z.loc[order],
            "n_cells": subset.obs["transferred_cc_phase"].value_counts().reindex(PHASES),
        }
    return panels


def supplementary_4(adata, panel: dict, title: str):
    """One of the two cell-cycle heatmaps."""
    z = panel["z"]
    labels = [f"{g} {adata.var.loc[g, 'gene_description']}" for g in z.index]
    figure, ax = plt.subplots(figsize=(7, 12))
    sns.heatmap(
        z, cmap="RdBu_r", center=0, vmin=-2, vmax=2, ax=ax,
        cbar_kws={"label": "Normalized expression (z-score)", "shrink": 0.4},
    )
    ax.set_yticklabels(labels, fontsize=6, rotation=0)
    ax.set_xlabel("Phase")
    ax.set_ylabel(None)
    ax.set_title(title)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    plt.subplots_adjust(left=0.55, right=0.88, top=0.94, bottom=0.04)
    return figure


# -------------------------------------------------------------------- Supp 5


def supplementary_5a_umap(adata):
    """Same panel as Figure 1B, drawn with a legend instead of in-plot numbers."""
    return _cluster_umap(in_vivo(adata), "nr_cluster", legend=True)


def supplementary_5a_counts(adata):
    return figure_1b_counts(adata)


def supplementary_5b_umap(adata):
    """In vitro bradyzoites coloured by transferred cluster identity."""
    return _cluster_umap(in_vitro_bradyzoites(adata), "cluster", legend=True)


def supplementary_5b_counts(adata):
    counts = in_vitro_bradyzoites(adata).obs["cluster"].value_counts().sort_index()
    return _count_bars(counts)


def supplementary_5d(adata):
    """srs22a expression across the in vivo bradyzoites."""
    label, gene_id = SUPP5D_GENE
    return _gene_umap(in_vivo(adata), gene_id, label)


# ------------------------------------------------------- Supp 5 volcano panels


def _gene_label(row) -> str:
    """What the published panels print beside a called-out point."""
    if row["gene_name"]:
        return str(row["gene_name"]).lower()
    text = str(row["gene_description"]).split(",")[0].split("(")[0].strip()
    return text if len(text) <= 22 else text[:21] + "…"


def _volcano(table: pd.DataFrame, panel: str):
    """One volcano: every significant gene, called-out ones coloured and named.

    The published panels colour their points by gene class rather than by
    significance, and only the genes the caption calls out are coloured at all —
    colouring every hypothetical protein would turn half of 5C pink.

    Genes whose adjusted p-value falls off the top of the published axis are not
    plotted, which is what the published panels do: 5C has 100 such genes and
    prints no marker for any of them. They are still in the exported DE table.
    """
    axes = VOLCANO_AXES[panel]
    low, high = VOLCANO_LABEL_RANGE[panel]
    comparison = de.COMPARISONS[panel]

    points = table[table["significant"]]
    with np.errstate(divide="ignore"):
        on_scale = -np.log10(points["pval_adj"].to_numpy()) <= axes["ymax"]
    points = points[on_scale]
    x = points["log2fc"].to_numpy()
    with np.errstate(divide="ignore"):
        y = -np.log10(points["pval_adj"].to_numpy())

    called_out = (x < low) | (x > high)
    ribosomal = points["gene_description"].str.contains("ribosomal protein").to_numpy()
    cyst_wall = (
        points["gene_id"].isin(SUPP5E_CYST_WALL_GENES).to_numpy()
        if panel == "5E"
        else np.zeros(len(points), dtype=bool)
    )
    category = np.where(
        cyst_wall,
        "cyst wall protein",
        np.where(
            called_out & ribosomal,
            "ribosomal protein",
            np.where(called_out, "called out", "other"),
        ),
    )
    legend = {
        "called out": f"log2FC > {high:g} or < {low:g}",
        "ribosomal protein": "ribosomal protein",
        "cyst wall protein": "cyst wall protein",
    }

    figure, ax = plt.subplots(figsize=(11, 3.6))
    for name in ("other", "ribosomal protein", "cyst wall protein", "called out"):
        pick = category == name
        if not pick.any():
            continue
        ax.scatter(
            x[pick], y[pick], s=8, marker="o", linewidths=0,
            color=VOLCANO_COLORS[name], label=legend.get(name),
        )

    # The published panels name the called-out genes that have something to
    # print — not the hypothetical proteins, and not the ribosomal ones, which
    # are what the grey is for. Names are nudged apart from the bottom up, since
    # the called-out genes crowd into a narrow band on the negative side.
    named = [
        (px, py, _gene_label(row))
        for (_, row), px, py in zip(points.iterrows(), x, y)
        if (px < low or px > high)
        and row["gene_description"] != "hypothetical protein"
        and "ribosomal protein" not in row["gene_description"]
    ]
    for on_the_left in (True, False):
        column = sorted(
            (item for item in named if bool(item[0] < 0) == on_the_left), key=lambda i: i[1]
        )
        placed = -np.inf
        for px, py, text in column:
            placed = max(py, placed + 0.045 * axes["ymax"])
            ax.annotate(
                text,
                (px, py),
                xytext=(px - 0.3 if on_the_left else px + 0.3, placed),
                textcoords="data",
                ha="right" if on_the_left else "left",
                va="center",
                fontsize=5,
                style="italic",
                arrowprops={"arrowstyle": "-", "linewidth": 0.3, "color": "0.65",
                            "shrinkA": 0, "shrinkB": 1},
            )

    ticks = list(range(0, axes["ymax"] + 1, axes["ystep"]))
    ax.set_xlim(*axes["xlim"])
    ax.set_xticks(list(axes["xticks"]))
    ax.set_ylim(0, axes["ymax"] * 1.02)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"$10^{{{-tick}}}$" for tick in ticks])
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("pval adj")
    ax.set_title(f"{comparison.label_a} vs {comparison.label_b}", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    if len(ax.collections) > 1:
        ax.legend(loc="upper right", frameon=False, fontsize=6, markerscale=1.5)
    figure.tight_layout()
    return figure


def supplementary_5c(adata, datadir=DATA):
    """In vivo against in vitro bradyzoites — 664 up / 1443 down in the paper."""
    return _volcano(de.volcano_table(adata, "5C", datadir), "5C")


def supplementary_5e(adata, datadir=DATA):
    """In vitro bradyzoite G1 against in vitro tachyzoite G1 — no published counts."""
    return _volcano(de.volcano_table(adata, "5E", datadir), "5E")


def supplementary_5f(adata, datadir=DATA):
    """In vivo against in vitro bradyzoite G1 — 676 up / 1146 down in the paper."""
    return _volcano(de.volcano_table(adata, "5F", datadir), "5F")
