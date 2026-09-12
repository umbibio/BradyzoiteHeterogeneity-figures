#!/usr/bin/env python
"""Build ``preview/index.html`` — one self-contained page showing every panel.

    python scripts/make_preview.py                  # -> preview/index.html
    python scripts/make_preview.py --max-edge 1400  # smaller embedded images

Every rendered panel with the metadata it was drawn from, in one file that makes
**no network request of any kind**. Images are embedded as ``data:`` URIs, the
stylesheet and script are inline, and nothing is fetched at load time.

Everything the page states is read from the repository rather than retyped:

* panel images            ``figures/*.png``
* colour limits, palettes,
  gene lists, constants    ``src/bzfig/constants.py``
* cell counts             ``data/obs.csv.gz`` and ``data/MANIFEST.json``
* panel descriptions      ``docs/reproducibility.md``

Requires ``pillow`` and ``jinja2`` (``pip install -e '.[preview]'``) in addition
to the rendering dependencies.
"""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
import io
import json
import re
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from bzfig import constants as K  # noqa: E402
from bzfig import de  # noqa: E402  (for the published volcano comparisons)

PAGE_TITLE = "Bradyzoite Heterogeneity Figures"

# Image embedding. The renders are 300 dpi; a long edge of 2000 px is more than
# a screen can show and keeps the whole page around 4 MB.
MAX_EDGE = 2000
QUALITY = 95
# Flat, text-heavy panels (heatmaps, bar charts, violins) compress better
# losslessly than lossily, and their small type has to stay crisp; dense scatter
# panels do not. Take the lossless encoding whenever it is not much bigger.
LOSSLESS_RATIO = 1.7


# --------------------------------------------------------------------- helpers


def esc(text: object) -> str:
    return html.escape(str(text), quote=True)


def number(value: int) -> str:
    return f"{value:,}"


# ------------------------------------------------------------------- markdown
#
# Just enough to read the panel table and the section headings out of
# docs/reproducibility.md.


def _split_row(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def _table_rows(lines: list[str]) -> list[list[str]]:
    rows = [_split_row(line) for line in lines]
    return [row for row in rows if not all(set(cell) <= set("-: ") for cell in row)]


def read_tables(markdown: str) -> list[list[dict[str, str]]]:
    """Every pipe table in ``markdown``, each as a list of dicts keyed by its header."""
    blocks: list[list[str]] = []
    block: list[str] = []
    for line in markdown.splitlines():
        if line.lstrip().startswith("|"):
            block.append(line)
        elif block:
            blocks.append(block)
            block = []
    if block:
        blocks.append(block)

    tables = []
    for lines in blocks:
        rows = _table_rows(lines)
        if rows:
            header, *body = rows
            tables.append([dict(zip(header, row)) for row in body])
    return tables


def find_table(markdown: str, *columns: str) -> list[dict[str, str]]:
    """The first table in ``markdown`` that has all of ``columns``, else empty.

    Picked by column rather than by position so that a table added to the docs
    above this one does not silently change what the page reads.
    """
    for table in read_tables(markdown):
        if table and all(column in table[0] for column in columns):
            return table
    return []


def split_sections(markdown: str) -> list[tuple[str, str]]:
    """Every ``## ...`` section as (heading, text including the heading line)."""
    heads = list(re.finditer(r"^##\s+(.*)$", markdown, re.M))
    return [
        (
            head.group(1).strip(),
            markdown[head.start() : heads[i + 1].start() if i + 1 < len(heads) else len(markdown)],
        )
        for i, head in enumerate(heads)
    ]


def doc_section(sections: list[tuple[str, str]], prefix: str) -> str:
    """The section whose heading starts with ``prefix``.

    Matched on a prefix rather than the whole heading so that the docs can be
    reworded without breaking this page.
    """
    for heading, text in sections:
        if heading.startswith(prefix):
            return text
    raise SystemExit(f"make_preview: docs/reproducibility.md has no '## {prefix}…' section")


# ---------------------------------------------------------------------- images


def encode_image(path: Path, max_edge: int, quality: int) -> dict:
    from PIL import Image

    image = Image.open(path).convert("RGB")
    native = image.size
    scale = max_edge / max(native)
    if scale < 1:
        image = image.resize((round(native[0] * scale), round(native[1] * scale)), Image.LANCZOS)

    lossy = io.BytesIO()
    image.save(lossy, format="WEBP", quality=quality, method=6)
    lossless = io.BytesIO()
    image.save(lossless, format="WEBP", lossless=True, method=6)
    best = lossless if lossless.tell() <= LOSSLESS_RATIO * lossy.tell() else lossy

    payload = best.getvalue()
    return {
        "uri": "data:image/webp;base64," + base64.b64encode(payload).decode("ascii"),
        "width": image.width,
        "height": image.height,
        "bytes": len(payload),
        "lossless": best is lossless,
        "native": f"{native[0]}×{native[1]}",
    }


def colour_strip(cmap_name: str, reverse: bool = False) -> str:
    """A 256x1 sample of a matplotlib colormap, as a data URI."""
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np
    from matplotlib import colormaps
    from PIL import Image

    ramp = np.linspace(1, 0, 256) if reverse else np.linspace(0, 1, 256)
    rgba = colormaps[cmap_name](ramp, bytes=True)
    image = Image.fromarray(rgba[np.newaxis, :, :3], "RGB")
    buffer = io.BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


# ------------------------------------------------------------------ dataset


def _rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle))


def dataset_facts(datadir: Path) -> dict:
    """Cell and gene counts, read from the deposited data."""
    manifest = json.loads((datadir / "MANIFEST.json").read_text())
    facts = {
        "n_obs": manifest["n_obs"],
        "n_vars": manifest["n_vars"],
        "cc_counts": manifest["cell_cycle_group"]["counts"],
        "cc_definition": manifest["cell_cycle_group"]["definition"],
    }

    obs = _rows(datadir / "obs.csv.gz")
    in_vivo = [row for row in obs if row["orig_ident"] == "Nonreactivated"]
    in_vitro = [row for row in obs if row["orig_ident"] == "me49 Day 3"]
    facts["in_vivo"] = len(in_vivo)
    facts["in_vivo_clusters"] = Counter(row["nr_cluster"] for row in in_vivo)
    facts["in_vitro"] = len(in_vitro)
    facts["in_vitro_clusters"] = Counter(row["cluster"] for row in in_vitro)
    for group in ("CCC", "MCC"):
        facts[f"{group}_phases"] = Counter(
            row["transferred_cc_phase"] for row in obs if row["cell_cycle_group"] == group
        )

    # The Figure 3 / Figure 6 cohorts, and the phase labels they are coloured by.
    facts["cohorts"] = {}
    for _, _, orig_ident, _ in K.COHORTS:
        cells = [row for row in obs if row["orig_ident"] == orig_ident]
        facts["cohorts"][orig_ident] = {
            "n": len(cells),
            "phases": Counter(row["transferred_cc_phase"] for row in cells),
        }
    facts["in_vitro_all"] = sum(1 for row in obs if row["dataset_type"] == "inVitro")
    facts["phase_colors"] = json.loads((datadir / "uns_colors.json").read_text())[
        "cc_phase_colors"
    ]

    # Figure 3E and 3F come from the separate scVI integration.
    integration = _rows(datadir / "figure_3ef_obs.csv.gz")
    tachyzoites = [row for row in integration if row["sample"] in ("S1", "S2")]
    facts["figure_3ef"] = {
        "n_obs": manifest["figure_3ef"]["n_obs"],
        "n_vars": manifest["figure_3ef"]["n_vars"],
        "samples": Counter(row["sample"] for row in integration),
        "tachyzoite": len(tachyzoites),
        "phases": {
            phase: sum(1 for row in tachyzoites if row["cc_phase"] == phase)
            for phase in K.PHASES
        },
    }

    var = {row["gene_id"]: row for row in _rows(datadir / "var.csv.gz")}
    genes = set(var)
    # Supplementary 1B/1C genes whose printed label the annotation cannot
    # corroborate: no gene_name, and a description that names nothing.
    facts["unnamed_supp1"] = [
        (label, gene_id)
        for _, label, gene_id in K.SUPP1_GENES
        if gene_id in var
        and not var[gene_id]["gene_name"].strip()
        and "hypothetical" in var[gene_id]["gene_description"].lower()
    ]
    with (datadir / "supp1a_marker_genes.csv").open(newline="") as handle:
        markers = [row[0] for row in csv.reader(handle) if row]
    facts["supp1a_listed"] = len(markers)
    facts["supp1a_plotted"] = len(set(markers) & genes)

    # The volcano panels: the two groups of each published comparison, counted
    # the way bzfig.de._groups selects them, and what the test produced.
    facts["universe"] = manifest["recovered_genes"]["universe"]
    facts["de"] = manifest["verification"]["supplementary_5_de"]
    facts["volcano_cells"] = {}
    for panel, comparison in de.COMPARISONS.items():

        def size(sample: str, g1_only: bool = comparison.g1_only) -> int:
            return sum(
                1
                for row in obs
                if row["orig_ident"] == sample
                and (not g1_only or row["cc_phase"] in K.G1_PHASES)
            )

        facts["volcano_cells"][panel] = {
            "a": size(comparison.sample_a),
            "b": size(comparison.sample_b),
        }
    return facts


def phase_counts(counter: Counter) -> str:
    return " · ".join(f"{phase} {number(counter.get(phase, 0))}" for phase in K.PHASES)


def cluster_counts(counter: Counter) -> str:
    keys = sorted(k for k in counter if k)
    return " · ".join(f"{key}: {number(counter[key])}" for key in keys)


# ------------------------------------------------------------------- content

IN_VIVO_SUBSET = "obs.orig_ident == 'Nonreactivated'"
IN_VITRO_SUBSET = "obs.orig_ident == 'me49 Day 3'"
UMAP_ROW = (
    "<code>obsm['3d_umap_harmony_integration']</code> on matplotlib 3-D axes, "
    "elev 60°, azim 0° (<code>bzfig.scatter3d.plot_3d_preview</code>)"
)
MINUS = "−"


def swatches(colors: list[str]) -> str:
    dots = "".join(
        f'<span class="swatch" style="background:{esc(c)}" title="{esc(c)}"></span>' for c in colors
    )
    return f'<span class="swatches">{dots}</span>'


def swatch_key(pairs: list[tuple[str, str]]) -> str:
    """Colours as a small swatch-and-label key: (label, hex) in drawing order."""
    items = "".join(
        f'<span class="key"><span class="swatch" style="background:{esc(colour)}"'
        f' title="{esc(colour)}"></span>{esc(label)}</span>'
        for label, colour in pairs
    )
    return f'<span class="keys">{items}</span>'


def colour_key(names: list[str]) -> str:
    """The volcano point colours, as a small swatch-and-label key."""
    return swatch_key([(name, K.VOLCANO_COLORS[name]) for name in names])


def signed(value: float, places: int = 2, trim: bool = False) -> str:
    """A number for prose, with a typographic minus."""
    text = f"{value:.{places}f}"
    if trim and "." in text:
        text = text.rstrip("0").rstrip(".")
    return text.replace("-", MINUS)


def scale_bar(uri: str, low: str, high: str, label: str) -> str:
    return (
        f'<span class="scale"><span class="scale-end">{esc(low)}</span>'
        f'<img class="scale-bar" src="{uri}" alt="" aria-hidden="true">'
        f'<span class="scale-end">{esc(high)}</span>'
        f'<span class="scale-label">{esc(label)}</span></span>'
    )


def build_panels(facts: dict, verdicts: dict[str, str], figures: Path) -> list[dict]:
    """The panel table. Values come from constants.py, the data and the docs."""
    vmin = (
        f"{MINUS}{abs(K.EXPRESSION_VMIN):g}"
        if K.EXPRESSION_VMIN < 0
        else f"{K.EXPRESSION_VMIN:g}"
    )
    vmax = f"{K.EXPRESSION_VMAX:g}"
    expression_scale = scale_bar(colour_strip(K.EXPRESSION_CMAP), vmin, vmax, K.EXPRESSION_CMAP)
    heatmap_scale = scale_bar(
        colour_strip("viridis"), "0", f"{K.HEATMAP_VMAX:g}", "viridis (scanpy default)"
    )
    supp4_scale = scale_bar(colour_strip("RdBu_r"), f"{MINUS}2", "2", "RdBu_r")
    palette = swatches(K.CLUSTER_COLORS)

    in_vivo_cells = (
        f"{number(facts['in_vivo'])} in vivo bradyzoites — <code>{esc(IN_VIVO_SUBSET)}</code>"
    )
    in_vivo_clusters = f"cells per cluster — {cluster_counts(facts['in_vivo_clusters'])}"
    scaled_layer = (
        "<code>logcounts_scaled</code> — <code>logcounts</code> times a per-(dataset, gene) "
        "factor, rebuilt on load"
    )

    def gene(label: str, gene_id: str) -> str:
        return f"<em>{esc(label)}</em> <code>{esc(gene_id)}</code>"

    # The one-line summary revealed when a panel is hovered: cells, layer, colour.
    dot = " · "
    in_vivo_short = f"{number(facts['in_vivo'])} in vivo cells"
    expression_overlay = dot.join(
        (
            in_vivo_short,
            "logcounts_scaled",
            f"{K.EXPRESSION_CMAP} {MINUS}{abs(K.EXPRESSION_VMIN):g} to {K.EXPRESSION_VMAX:g}, "
            "fixed for every gene panel",
        )
    )
    heatmap_overlay = dot.join(
        (in_vivo_short, "logcounts_scaled", f"viridis, vmax {K.HEATMAP_VMAX}")
    )

    panels: list[dict] = []

    panels.append(
        {
            "id": "fig-1b",
            "overlay": dot.join(
                (in_vivo_short, "coloured by obs.nr_cluster", "cluster palette")
            ),
            "group": "Figure 1",
            "label": "Figure 1B",
            "title": "Cluster UMAP and cells per cluster",
            "lede": verdicts.get("Figure 1B", ""),
            "plates": [
                {"file": "Figure_1B_umap", "caption": "Cluster UMAP of the in vivo bradyzoites"},
                {"file": "Figure_1B_cells_per_cluster", "caption": "Cells per cluster"},
            ],
            "meta": [
                ("Cells", f"{in_vivo_cells}<br>{esc(in_vivo_clusters)}"),
                ("Embedding", UMAP_ROW),
                ("Expression", "none — the panel colours cells by <code>obs.nr_cluster</code>"),
                ("Colour", f"{palette} cluster palette, <code>constants.CLUSTER_COLORS</code>"),
                (
                    "Drawn by",
                    "<code>bzfig.panels.figure_1b_umap</code>, "
                    "<code>bzfig.panels.figure_1b_counts</code>",
                ),
            ],
        }
    )

    panels.append(
        {
            "id": "fig-1c",
            "overlay": f"{heatmap_overlay}{dot}{len(K.FIG1C_GENES)} marker genes",
            "group": "Figure 1",
            "label": "Figure 1C",
            "title": "Selected marker genes per cluster",
            "lede": verdicts.get("Figure 1C", ""),
            "plates": [{"file": "Figure_1C_marker_heatmap", "caption": "Marker gene heatmap"}],
            "meta": [
                ("Cells", f"{in_vivo_cells}, grouped by <code>obs.nr_cluster</code>"),
                ("Data layer", scaled_layer),
                (
                    "Genes",
                    f"{len(K.FIG1C_GENES)} genes in the plotted order, "
                    f"<code>constants.FIG1C_GENES</code> — row groups "
                    + ", ".join(
                        f"{esc(name)} {count}"
                        for name, count in Counter(
                            group for group, _, _ in K.FIG1C_GENES
                        ).items()
                    ),
                ),
                (
                    "Colour",
                    f"{heatmap_scale} vmax {K.HEATMAP_VMAX} "
                    "(<code>constants.HEATMAP_VMAX</code>)",
                ),
                ("Drawn by", "<code>bzfig.panels.figure_1c</code>"),
            ],
        }
    )

    panels.append(
        {
            "id": "fig-1d",
            "overlay": dot.join(
                (
                    "recorded constants, not data",
                    ", ".join(
                        str(value)
                        for _, value in sorted(K.UNIQUE_MARKERS_PER_CLUSTER.items())
                    ),
                    "cluster palette",
                )
            ),
            "group": "Figure 1",
            "label": "Figure 1D",
            "title": "Unique markers per cluster",
            "lede": (
                "Drawn from the six values hard-coded in the analysis notebook, which no "
                "re-run of the marker test recovers."
            ),
            "plates": [
                {
                    "file": "Figure_1D_unique_markers_per_cluster",
                    "caption": "# Unique markers/Cluster, from the recorded constants",
                }
            ],
            "meta": [
                ("Cells", "none — the panel is drawn from recorded values, not from the data"),
                (
                    "Values",
                    "<code>constants.UNIQUE_MARKERS_PER_CLUSTER</code> = "
                    + ", ".join(
                        f"cluster {k} {v}" for k, v in sorted(K.UNIQUE_MARKERS_PER_CLUSTER.items())
                    ),
                ),
                (
                    "Re-derived",
                    "the marker test gives <code>[37, 94, 84, 861, 18, 52]</code> on the "
                    f"deposited {number(facts['n_vars'])} genes and "
                    "<code>[45, 93, 84, 863, 18, 53]</code> on the full "
                    f"{number(facts['universe']['toxodb_65_genes'])}-gene universe; neither is "
                    "the published set",
                ),
                ("Colour", f"{palette} cluster palette, <code>constants.CLUSTER_COLORS</code>"),
                ("Drawn by", "<code>bzfig.panels.figure_1d</code>"),
            ],
        }
    )

    variants = []
    for row_group, label, gene_id in K.FIG1E_GENES:
        variants.append(
            {
                "key": label,
                "label": label,
                "row_group": row_group,
                "file": f"Figure_1E_{label.replace(' ', '_')}",
                "caption": f"{label} — {gene_id} ({row_group.lower()})",
            }
        )
    panels.append(
        {
            "id": "fig-1e",
            "overlay": expression_overlay,
            "group": "Figure 1",
            "label": "Figure 1E",
            "title": "Per-gene expression UMAPs",
            "lede": verdicts.get("Figure 1E", ""),
            "variants": variants,
            "switch_hint": "Six panels, one per gene. Switch between them in place:",
            "meta": [
                ("Cells", in_vivo_cells),
                ("Embedding", UMAP_ROW),
                ("Data layer", scaled_layer),
                (
                    "Genes",
                    "<code>constants.FIG1E_GENES</code><br>"
                    + "<br>".join(
                        f"{esc(row_group)}: {gene(label, gene_id)}"
                        for row_group, label, gene_id in K.FIG1E_GENES
                    ),
                ),
                (
                    "Colour",
                    f"{expression_scale} vmin {vmin}, vmax {vmax} "
                    "(<code>constants.EXPRESSION_VMIN</code> / <code>EXPRESSION_VMAX</code>), "
                    "the same scale for every gene panel",
                ),
                ("Drawn by", "<code>bzfig.panels.figure_1e</code>"),
            ],
        }
    )

    label_1f, gene_1f = K.FIG1F_GENE
    panels.append(
        {
            "id": "fig-1f",
            "overlay": expression_overlay,
            "group": "Figure 1",
            "label": "Figure 1F",
            "title": f"{label_1f} expression",
            "lede": verdicts.get("Figure 1F", ""),
            "plates": [{"file": "Figure_1F_cst1_srs44", "caption": f"{label_1f} — {gene_1f}"}],
            "meta": [
                ("Cells", in_vivo_cells),
                ("Embedding", UMAP_ROW),
                ("Data layer", scaled_layer),
                ("Gene", f"<code>constants.FIG1F_GENE</code> — {gene(label_1f, gene_1f)}"),
                (
                    "Colour",
                    f"{expression_scale} vmin {vmin}, vmax {vmax}, "
                    "shared with Figure 1E and Supplementary 5D",
                ),
                ("Drawn by", "<code>bzfig.panels.figure_1f</code>"),
            ],
        }
    )

    panels.append(
        {
            "id": "fig-1g",
            "overlay": dot.join(
                (in_vivo_short, "logcounts_scaled", f"{gene_1f}", "violins per cluster")
            ),
            "group": "Figure 1",
            "label": "Figure 1G",
            "title": "CST1/SRS44 expression per cluster",
            "lede": verdicts.get("Figure 1G", ""),
            "plates": [{"file": "Figure_1G_cst1_violin", "caption": "Violins per cluster"}],
            "meta": [
                ("Cells", f"{in_vivo_cells}<br>{esc(in_vivo_clusters)}"),
                ("Data layer", scaled_layer),
                ("Gene", gene(label_1f, gene_1f)),
                ("Colour", f"{palette} cluster palette, <code>constants.CLUSTER_COLORS</code>"),
                ("Drawn by", "<code>bzfig.panels.figure_1g</code>"),
                (
                    "Cross-check",
                    "the plotted values were re-derived from the data object and re-plotted from "
                    "the exported table alone; both matched the published panel",
                ),
            ],
        }
    )

    panels.append(
        {
            "id": "fig-2h",
            "overlay": " · ".join(
                (
                    f"{number(facts['in_vivo'])} in vivo cells",
                    "logcounts, re-closed over the 8,170 deposited genes",
                    "Pearson r",
                )
            ),
            "group": "Figure 2",
            "label": "Figure 2H",
            "title": "Cyst wall protein co-expression",
            "lede": (
                "Five cyst wall proteins against the 29 displayed genes, in the published "
                "column order. Contributed by Kourosh Zarringhalam."
            ),
            "plates": [
                {
                    "file": "Figure_2H_correlation_heatmap",
                    "caption": "Pearson correlation across the in vivo bradyzoites",
                }
            ],
            "meta": [
                ("Cells", in_vivo_cells),
                (
                    "Data layer",
                    "<code>logcounts</code>, re-closed over the 8,170 deposited genes as "
                    "<code>log1p(1e4 * expm1(L) / rowSum(expm1(L)))</code> — the "
                    "panel-specific normalisation the original R heatmaps used, recovered",
                ),
                (
                    "Values",
                    "Pearson correlation between each of the five cyst wall proteins and each of "
                    "the 29 displayed genes, in the original column order",
                ),
                (
                    "Genes",
                    "labels and column order from "
                    "<code>figure_2h_supplementary_4_metadata.json</code>",
                ),
                ("Colour", "diverging red/blue, —1 to 1, centred on 0"),
                (
                    "Drawn by",
                    "<code>bzfig.figure_2h_supplementary_4.figure_2h</code>",
                ),
                (
                    "Tables",
                    "the full-precision correlation matrix is in "
                    "<code>figures/tables/Figure_2H_correlations.csv</code>; method in "
                    "<code>docs/figure-2h-supplementary-4.md</code>",
                ),
            ],
        }
    )

    panels.extend(cohort_panels(facts, verdicts, figures))

    panels.append(
        {
            "id": "supp-1a",
            "overlay": f"{heatmap_overlay}{dot}{number(facts['supp1a_plotted'])} marker genes",
            "group": "Supplementary 1",
            "label": "Supplementary 1A",
            "title": "All cluster markers per cluster",
            "lede": verdicts.get("Supplementary 1A", ""),
            "plates": [
                {
                    "file": "Supplementary_1A_all_markers_heatmap",
                    "caption": "All marker genes",
                }
            ],
            "meta": [
                ("Cells", f"{in_vivo_cells}, grouped by <code>obs.nr_cluster</code>"),
                ("Data layer", scaled_layer),
                (
                    "Genes",
                    f"{number(facts['supp1a_listed'])} genes from "
                    "<code>data/supp1a_marker_genes.csv</code>, copied verbatim from the "
                    f"analysis; {number(facts['supp1a_plotted'])} of them are in the deposited "
                    "object and are "
                    "plotted",
                ),
                (
                    "Colour",
                    f"{heatmap_scale} vmax {K.HEATMAP_VMAX} "
                    "(<code>constants.HEATMAP_VMAX</code>)",
                ),
                ("Drawn by", "<code>bzfig.panels.supplementary_1a</code>"),
                (
                    "Rest of the figure",
                    "Supplementary 1B and 1C are the eleven per-gene expression UMAPs below, from "
                    "this same dataset. Panel 1B also carries a BioRender cartoon of the two "
                    "microneme subpopulations, which is figure assembly, not a panel. (The Benke "
                    "et al. caption belongs to Supplementary Figure 2.)",
                ),
            ],
        }
    )

    for supp1_panel in ("1B", "1C"):
        row_group = SUPP1_ROW_GROUPS[supp1_panel]
        members = [row for row in K.SUPP1_GENES if row[0] == supp1_panel]
        unnamed = [
            (label, gene_id)
            for label, gene_id in facts["unnamed_supp1"]
            if any(label == member[1] for member in members)
        ]
        meta = [
            ("Cells", in_vivo_cells),
            ("Embedding", UMAP_ROW),
            ("Data layer", scaled_layer),
            (
                "Genes",
                f"<code>constants.SUPP1_GENES</code>, the {len(members)} of panel "
                f"{supp1_panel}<br>"
                + "<br>".join(gene(label, gene_id) for _, label, gene_id in members),
            ),
            (
                "Colour",
                f"{expression_scale} vmin {vmin}, vmax {vmax}, the same fixed scale as "
                "Figure 1E, 1F and Supplementary 5D",
            ),
            ("Drawn by", "<code>bzfig.panels.supplementary_1bc</code>"),
        ]
        if unnamed:
            meta.append(
                (
                    "Labels not in the annotation",
                    ", ".join(gene(label, gene_id) for label, gene_id in unnamed)
                    + (" is" if len(unnamed) == 1 else " are")
                    + ' "hypothetical protein" in the ToxoDB-65 annotation shipped here, with no '
                    "<code>var.gene_name</code>. Those labels rest on the published caption and "
                    "on the per-gene PNGs the analysis left behind, not on the data.",
                )
            )
        panels.append(
            {
                "id": f"supp-{supp1_panel.lower()}",
                "overlay": expression_overlay,
                "group": "Supplementary 1",
                "label": f"Supplementary {supp1_panel}",
                "title": f"{row_group}, per-gene expression",
                "lede": verdicts.get(f"Supplementary {supp1_panel}", ""),
                "switch_hint": (
                    f"{len(members)} panels, one per gene. Switch between them in place:"
                ),
                "variants": [
                    {
                        "key": label,
                        "label": label,
                        "row_group": row_group,
                        "file": f"Supplementary_{supp1_panel}_{label}",
                        "caption": f"{label} — {gene_id}",
                    }
                    for _, label, gene_id in members
                ],
                "meta": meta,
            }
        )

    panels.append(
        {
            "id": "supp-4",
            "overlay": dot.join(
                (
                    f"CCC {number(facts['cc_counts']['CCC'])} / "
                    f"MCC {number(facts['cc_counts']['MCC'])} cells",
                    "logcounts, re-closed over the 8,170 deposited genes",
                    "per-gene z-score across the five phases",
                    "LAB colour scale",
                )
            ),
            "group": "Supplementary 4",
            "label": "Supplementary 4",
            "title": "Cell-cycle regulators, common and modified cell cycle",
            "lede": (
                "Fifty cell-cycle regulators by phase, with the published gene labels, row "
                "orders, colour scale and row highlights. Contributed by Kourosh Zarringhalam."
            ),
            "plates": [
                {"file": "Supplementary_4_CCC", "caption": "Common cell cycle (CCC)"},
                {"file": "Supplementary_4_MCC", "caption": "Modified cell cycle (MCC)"},
            ],
            "plate_layout": "pair",
            "meta": [
                (
                    "Cells",
                    "<code>obs.cell_cycle_group</code> — CCC "
                    f"{number(facts['cc_counts']['CCC'])} cells "
                    f"({esc(phase_counts(facts['CCC_phases']))}), "
                    f"MCC {number(facts['cc_counts']['MCC'])} cells "
                    f"({esc(phase_counts(facts['MCC_phases']))})",
                ),
                ("Group definition", f"<code>{esc(facts['cc_definition'])}</code>"),
                (
                    "Data layer",
                    "<code>logcounts</code>, re-closed over the 8,170 deposited genes as "
                    "<code>log1p(1e4 * expm1(L) / rowSum(expm1(L)))</code> — the "
                    "panel-specific normalisation the original R heatmaps used, recovered",
                ),
                (
                    "Values",
                    "mean expression per phase, z-scored per gene across the five phases with the "
                    "sample SD (ddof=1), independently within CCC and MCC",
                ),
                (
                    "Genes",
                    f"{len(K.SUPP4_GENES)} regulators — original labels and the two row orders "
                    "in <code>figure_2h_supplementary_4_metadata.json</code>; "
                    "<code>constants.SUPP4_GENES</code> is kept as an independent transcription "
                    "from the printed figure, and the two sets agree exactly",
                ),
                ("Colour", "LAB-interpolated scale from the original figure, breaks at ±2.5"),
                (
                    "Row highlights",
                    "the pale-pink description backgrounds are manual annotations transcribed "
                    "from the printed figure — 19 rows in CCC, 20 in MCC — not a rule "
                    "derived from the data",
                ),
                (
                    "Scale",
                    "the colour bar axis runs −4…4 but the values only span −1.77…+1.79; the "
                    "CCC group is thinly populated, so several genes are detected in one phase "
                    "only and their rows sit at the one-hot extremes of ±1.789 / −0.447",
                ),
                ("Drawn by", "<code>bzfig.figure_2h_supplementary_4.supplementary_4</code>"),
                (
                    "Tables",
                    "phase means, z-scores and cell counts as full-precision CSV in "
                    "<code>figures/tables/</code>; method in "
                    "<code>docs/figure-2h-supplementary-4.md</code>",
                ),
            ],
        }
    )

    panels.append(
        {
            "id": "supp-5a",
            "overlay": dot.join(
                (in_vivo_short, "coloured by obs.nr_cluster", "cluster palette, with a legend")
            ),
            "group": "Supplementary 5",
            "label": "Supplementary 5A",
            "title": "In vivo clusters, with a legend",
            "lede": verdicts.get("Supplementary 5A", ""),
            "plates": [
                {"file": "Supplementary_5A_umap", "caption": "Cluster UMAP with a legend"},
                {"file": "Supplementary_5A_cells_per_cluster", "caption": "Cells per cluster"},
            ],
            "meta": [
                ("Cells", f"{in_vivo_cells}<br>{esc(in_vivo_clusters)}"),
                ("Embedding", UMAP_ROW),
                ("Colour", f"{palette} cluster palette, <code>constants.CLUSTER_COLORS</code>"),
                (
                    "Drawn by",
                    "<code>bzfig.panels.supplementary_5a_umap</code>, "
                    "<code>bzfig.panels.supplementary_5a_counts</code>",
                ),
                (
                    "Relation to Figure 1B",
                    "the same panel, drawn with a legend instead of in-plot numbers",
                ),
            ],
        }
    )

    panels.append(
        {
            "id": "supp-5b",
            "overlay": dot.join(
                (
                    f"{number(facts['in_vitro'])} in vitro cells",
                    "coloured by obs.cluster, the transferred identity",
                    "cluster palette",
                )
            ),
            "group": "Supplementary 5",
            "label": "Supplementary 5B",
            "title": "In vitro bradyzoites by transferred identity",
            "lede": verdicts.get("Supplementary 5B", ""),
            "plates": [
                {
                    "file": "Supplementary_5B_umap",
                    "caption": "In vitro bradyzoites, transferred identities",
                },
                {"file": "Supplementary_5B_cells_per_cluster", "caption": "Cells per cluster"},
            ],
            "meta": [
                (
                    "Cells",
                    f"{number(facts['in_vitro'])} in vitro bradyzoites — "
                    f"<code>{esc(IN_VITRO_SUBSET)}</code><br>"
                    f"transferred clusters {esc(cluster_counts(facts['in_vitro_clusters']))}",
                ),
                ("Embedding", UMAP_ROW),
                (
                    "Expression",
                    "none — cells are coloured by <code>obs.cluster</code>, the transferred "
                    "identity",
                ),
                ("Colour", f"{palette} cluster palette, <code>constants.CLUSTER_COLORS</code>"),
                (
                    "Drawn by",
                    "<code>bzfig.panels.supplementary_5b_umap</code>, "
                    "<code>bzfig.panels.supplementary_5b_counts</code>",
                ),
            ],
        }
    )

    panels.append(volcano_panel("5C", facts))

    label_5d, gene_5d = K.SUPP5D_GENE
    panels.append(
        {
            "id": "supp-5d",
            "overlay": expression_overlay,
            "group": "Supplementary 5",
            "label": "Supplementary 5D",
            "title": f"{label_5d} expression",
            "lede": verdicts.get("Supplementary 5D", ""),
            "plates": [{"file": "Supplementary_5D_srs22a", "caption": f"{label_5d} — {gene_5d}"}],
            "meta": [
                ("Cells", in_vivo_cells),
                ("Embedding", UMAP_ROW),
                ("Data layer", scaled_layer),
                ("Gene", f"<code>constants.SUPP5D_GENE</code> — {gene(label_5d, gene_5d)}"),
                (
                    "Colour",
                    f"{expression_scale} vmin {vmin}, vmax {vmax}, "
                    "shared with Figure 1E and 1F",
                ),
                ("Drawn by", "<code>bzfig.panels.supplementary_5d</code>"),
            ],
        }
    )

    panels.append(volcano_panel("5E", facts))
    panels.append(volcano_panel("5F", facts))

    return panels


# What the two Supplementary 1 gene panels show, as the published figure groups them.
SUPP1_ROW_GROUPS = {"1B": "Microneme transcripts", "1C": "Known cyst wall proteins"}

# How the cohorts of constants.COHORTS are named in prose.
COHORT_TITLES = {
    "Nonreactivated": "In vivo bradyzoites",
    "me49 Day 0": "me49 Day 0",
    "me49 Day 3": "me49 Day 3",
}


def same_bytes(first: Path, second: Path) -> bool:
    import hashlib

    return (
        first.exists()
        and second.exists()
        and hashlib.sha256(first.read_bytes()).digest()
        == hashlib.sha256(second.read_bytes()).digest()
    )


def phase_gaps(panel: str, deposited: Counter) -> dict[str, int]:
    """Published bar minus deposited labels, per phase."""
    published = K.FIG3_PHASE_BARS_PUBLISHED[panel]
    return {phase: published[phase] - deposited.get(phase, 0) for phase in K.PHASES}


def widest_gap(gaps: dict[str, int]) -> tuple[str, int]:
    phase = max(gaps, key=lambda name: abs(gaps[name]))
    return phase, abs(gaps[phase])


def phase_bar_row(panel: str, deposited: Counter, cohort_size: int) -> tuple[str, str]:
    """How the drawn phase bar compares with the one measured off the figure."""
    published = K.FIG3_PHASE_BARS_PUBLISHED[panel]
    drawn = " · ".join(f"{phase} {number(deposited.get(phase, 0))}" for phase in K.PHASES)
    printed = " · ".join(f"{phase} {number(published[phase])}" for phase in K.PHASES)
    gaps = phase_gaps(panel, deposited)
    difference = sum(abs(gap) for gap in gaps.values())
    if difference == 0:
        verdict = "the two agree exactly"
    else:
        phase, size = widest_gap(gaps)
        verdict = (
            f"they differ by at least {number(difference // 2)} of {number(cohort_size)} cells "
            f"({difference / 2 / cohort_size:.1%}), the widest single bar being "
            f"{number(size)} in {phase}"
        )
    return (
        "Phase bar",
        f"drawn from <code>obs.transferred_cc_phase</code> — {esc(drawn)}<br>"
        f"measured off the published bar — {esc(printed)} "
        f"(<code>constants.FIG3_PHASE_BARS_PUBLISHED</code>, scaled so the five bars sum to the "
        f"cohort size); {esc(verdict)}",
    )


def cohort_panels(facts: dict, verdicts: dict[str, str], figures: Path) -> list[dict]:
    """Figure 3A–3D and Figure 6 — one cohort of the projection at a time."""
    phase_key = swatch_key(list(zip(K.PHASES, facts["phase_colors"])))
    grey = f'<span class="swatch" style="background:{K.BACKGROUND_GRAY}"></span>'
    total = facts["n_obs"]
    built: list[dict] = []

    for letter, _, orig_ident, background in K.COHORTS:
        panel = f"Figure {letter}"
        cohort = facts["cohorts"][orig_ident]
        name = COHORT_TITLES[orig_ident]
        layer = (
            f"drawn over the whole {number(total)}-cell projection, the other "
            f"{number(total - cohort['n'])} cells in {grey} <code>{K.BACKGROUND_GRAY}</code> "
            "(<code>constants.BACKGROUND_GRAY</code>)"
            if background
            else "the cohort alone — this panel carries no grey layer"
        )
        meta = [
            (
                "Cells",
                f"{number(cohort['n'])} cells — <code>obs.orig_ident == '{esc(orig_ident)}'</code>"
                f"<br>{layer}",
            ),
            ("Embedding", UMAP_ROW),
            (
                "Colour",
                f"{phase_key}<code>obs.transferred_cc_phase</code>, palette from "
                "<code>uns['cc_phase_colors']</code>, which travels with the data",
            ),
            phase_bar_row(letter, cohort["phases"], cohort["n"]),
            (
                "Drawn by",
                "<code>bzfig.panels.figure_3_umap</code>, "
                "<code>bzfig.panels.figure_3_counts</code>",
            ),
        ]
        if letter == "3A":
            twin = same_bytes(
                figures / "Figure_3A_umap.png", figures / "Figure_6_nonreactivated.png"
            )
            meta.append(
                (
                    "Same as Figure 6",
                    "with no grey layer this is the same rendering as Figure 6's first panel"
                    + (
                        " — the two rendered files are byte-identical."
                        if twin
                        else ", though the two rendered files differ."
                    ),
                )
            )
        if letter == "3C":
            meta.append(
                (
                    "The gap",
                    "the phase vector behind the published panel is in neither repository — not "
                    "<code>transferred_cc_phase</code>, not <code>cc_phase</code>, and nothing in "
                    "the analysis repository. It colours the <em>scatter</em> as well as the bar, "
                    "so this panel puts orange (G1b) where the published one has a blue (G1a) "
                    "group in the upper cluster. The cells and their positions are unaffected.",
                )
            )

        lede = ""
        if letter == "3A":
            phase, size = widest_gap(phase_gaps(letter, cohort["phases"]))
            lede = (
                "The cohort and its positions are the published panel's. Its phase bar is close "
                f"to the published one but not identical — the widest gap is {number(size)} cells "
                f"in {phase} — and the two are compared in the panel metadata."
            )

        caption = f"{name} by cell-cycle phase"
        built.append(
            {
                "id": f"fig-{letter.lower()}",
                "overlay": " · ".join(
                    (
                        f"{number(cohort['n'])} cells",
                        "coloured by obs.transferred_cc_phase",
                        "grey layer" if background else "no grey layer",
                    )
                ),
                "group": "Figure 3",
                "label": panel,
                "title": f"{name} by cell-cycle phase",
                "lede": lede
                or verdicts.get(panel)
                or (
                    "The cohort and its positions are the published panel's; the phase labels are "
                    "not. The vector that coloured the published 3C is in neither repository, so "
                    "the bar and the scatter here both come from the deposited labels."
                ),
                "plates": [
                    {
                        "file": f"Figure_{letter}_umap",
                        "caption": caption
                        + (
                            " — phase labels differ from the published panel"
                            if letter == "3C"
                            else ""
                        ),
                    },
                    {"file": f"Figure_{letter}_cells_per_phase", "caption": "Cells per phase"},
                ],
                "meta": meta,
            }
        )

    built.append(
        {
            "id": "fig-3d",
            "overlay": " · ".join(
                (
                    f"{number(facts['in_vitro_all'])} in vitro cells",
                    "flat highlight, no phase colouring",
                    f"over all {number(total)} cells",
                )
            ),
            "group": "Figure 3",
            "label": "Figure 3D",
            "title": "The in vitro cells in the projection",
            "lede": verdicts.get("Figure 3D", ""),
            "plates": [{"file": "Figure_3D_in_vitro", "caption": "In vitro cells picked out"}],
            "meta": [
                (
                    "Cells",
                    f"{number(facts['in_vitro_all'])} in vitro cells — "
                    f"<code>obs.dataset_type == 'inVitro'</code> — over all {number(total)}",
                ),
                ("Embedding", UMAP_ROW),
                (
                    "Colour",
                    swatch_key(
                        [("in vitro", K.FIG3D_IN_VITRO), ("everything else", K.BACKGROUND_GRAY)]
                    )
                    + "flat, not by phase (<code>constants.FIG3D_IN_VITRO</code> and "
                    "<code>constants.BACKGROUND_GRAY</code>)",
                ),
                ("Drawn by", "<code>bzfig.panels.figure_3d</code>"),
            ],
        }
    )

    tach = facts["figure_3ef"]
    phases = " · ".join(f"{p} {number(n)}" for p, n in tach["phases"].items())
    integration_meta = [
        (
            "Cells",
            f"{number(tach['tachyzoite'])} in vivo tachyzoites (5 dpi, peritoneal cavity) "
            f"— samples S1 {number(tach['samples']['S1'])} and S2 "
            f"{number(tach['samples']['S2'])} — integrated with the "
            f"{number(tach['samples']['NR'])} in vivo bradyzoites",
        ),
        (
            "Embedding",
            "2-D UMAP of the 10-dimensional scVI latent space "
            "(<code>data/figure_3ef_embedding.csv.gz</code>), from the deposited checkpoint "
            "<code>data/figure_3ef_scvi_model.pt</code> over "
            f"{number(tach['n_obs'])} × {number(tach['n_vars'])} counts",
        ),
        (
            "Frame",
            "loading the checkpoint is deterministic and so is the UMAP that follows it, so "
            "this layout reproduces exactly; the published panel was drawn from a different "
            "fit, so its orientation on the page differs",
        ),
    ]
    built.append(
        {
            "id": "fig-3e",
            "overlay": " · ".join(
                (
                    f"{number(tach['tachyzoite'])} in vivo tachyzoites",
                    f"over {number(tach['samples']['NR'])} bradyzoites",
                    "scVI integration",
                )
            ),
            "group": "Figure 3",
            "label": "Figure 3E",
            "title": "In vivo tachyzoites on the integration",
            "lede": (
                f"The {number(tach['tachyzoite'])} in vivo tachyzoites picked out of the "
                "integration. They fall on the common cell cycle: their nearest bradyzoite "
                "neighbours in the latent space are 40.3% CCC, against a 6.9% baseline."
            ),
            "plates": [
                {"file": "Figure_3E_tachyzoite_highlight", "caption": "In vivo tachyzoites in red"}
            ],
            "meta": integration_meta
            + [
                (
                    "Colour",
                    swatch_key(
                        [("in vivo tachyzoite", K.FIG3EF_TACHYZOITE),
                         ("in vivo bradyzoite", K.BACKGROUND_GRAY)]
                    )
                    + "flat (<code>constants.FIG3EF_TACHYZOITE</code>)",
                ),
                ("Drawn by", "<code>bzfig.figure_3ef.figure_3e</code>"),
            ],
        }
    )
    built.append(
        {
            "id": "fig-3f",
            "overlay": " · ".join(
                (
                    f"{number(tach['tachyzoite'])} in vivo tachyzoites",
                    "coloured by obs.cc_phase",
                    "scVI integration",
                )
            ),
            "group": "Figure 3",
            "label": "Figure 3F",
            "title": "The same tachyzoites by cell-cycle phase",
            "lede": (
                "The tachyzoites coloured by their cell-cycle phase, which runs "
                "round the loop in order."
            ),
            "plates": [
                {"file": "Figure_3F_tachyzoite_cc_phase", "caption": "Tachyzoites by cell-cycle phase"}
            ],
            "meta": integration_meta
            + [
                ("Colour", f"{phase_key}<code>obs.cc_phase</code>, the palette Figure 3A\u20133C uses"),
                ("Cells per phase", phases),
                ("Drawn by", "<code>bzfig.figure_3ef.figure_3f</code>"),
            ],
        }
    )

    cohorts = [(COHORT_TITLES[ident], facts["cohorts"][ident]["n"]) for _, _, ident, _ in K.COHORTS]
    built.append(
        {
            "id": "fig-6",
            "overlay": " · ".join(
                (
                    " / ".join(f"{name} {number(size)}" for name, size in cohorts),
                    "coloured by obs.transferred_cc_phase",
                    "no grey layer",
                )
            ),
            "group": "Figure 6",
            "label": "Figure 6",
            "title": "The three cohorts, without the grey layer",
            "lede": (
                "The same three cohorts as Figure 3A–3C, each drawn alone. The first two "
                "reproduce exactly; the third inherits Figure 3C's missing phase vector."
            ),
            "plate_layout": "pair",
            "plates": [
                {
                    "file": f"Figure_6_{name}",
                    "caption": COHORT_TITLES[ident]
                    + {
                        "3A": " — the same rendering as Figure 3A",
                        "3C": " — phase labels differ from the published panel",
                    }.get(letter, ""),
                }
                for letter, name, ident, _ in K.COHORTS
            ],
            "meta": [
                (
                    "Cells",
                    "one cohort of <code>obs.orig_ident</code> per panel — "
                    + " · ".join(f"{esc(name)} {number(size)}" for name, size in cohorts),
                ),
                ("Embedding", UMAP_ROW),
                (
                    "Colour",
                    f"{phase_key}<code>obs.transferred_cc_phase</code>; Figure 6 never draws the "
                    "grey layer",
                ),
                (
                    "Relation to Figure 3",
                    "the same three cohorts as Figure 3A–3C. 3A already carries no grey layer, so "
                    "it and the first panel here are one rendering; 3B and 3C add it.",
                ),
                (
                    "Phase labels",
                    "the third panel is the me49 Day 3 cohort, so it inherits the gap described "
                    "under Figure 3C — positions right, phase colours from the deposited labels.",
                ),
                (
                    "Naming",
                    "the panels are named for the cohort they hold, not for the column header "
                    "printed above them, which the published figure does not match to the cells "
                    "beneath it",
                ),
                ("Drawn by", "<code>bzfig.panels.figure_6_umap</code>"),
            ],
        }
    )
    return built

VOLCANO_IDS = {"5C": "supp-5c", "5E": "supp-5e", "5F": "supp-5f"}

# Notes that belong to one panel only.
VOLCANO_NOTES = {
    "5C": (
        "Points off the top of the axis",
        "100 genes have an adjusted p-value below the published 1e-100 axis top — 20 of them "
        "underflow float64 to zero — and are not plotted, which is what the published panel "
        "does. They stay in the DE table, so the counts are unaffected. 5E and 5F have none.",
    ),
    "5E": (
        "Blue set read off the figure",
        "there is no cyst wall protein list in the data or the paper (two of the blue genes are "
        "a dense granule and a rhoptry protein), so the set was recovered by matching the "
        "published blue points to their fold changes; it is recorded in "
        f"<code>constants.SUPP5E_CYST_WALL_GENES</code>, which holds "
        f"{len(K.SUPP5E_CYST_WALL_GENES)} genes.",
    ),
}


def volcano_panel(panel: str, facts: dict) -> dict:
    """One of the three Supplementary 5 volcano panels."""
    comparison = de.COMPARISONS[panel]
    stats = facts["de"][panel]
    cells = facts["volcano_cells"][panel]
    universe = facts["universe"]
    axes = K.VOLCANO_AXES[panel]
    low, high = K.VOLCANO_LABEL_RANGE[panel]
    plotted = stats["up"] + stats["down"]
    published = stats["published"]
    title = f"{comparison.label_a} vs {comparison.label_b}"

    if published:
        lede = (
            f"{number(stats['up'])} genes up and {number(stats['down'])} down in group A, "
            f"from the full {number(universe['toxodb_65_genes'])}-gene universe; the published "
            f"panel gives {number(published['up_in_vivo'])} / "
            f"{number(published['down_in_vivo'])}."
        )
        counts = (
            f"{number(stats['up'])} up / {number(stats['down'])} down in group A "
            "(<code>constants.SUPP5_DE_COUNTS_REPRODUCED</code>), against a published "
            f"{number(published['up_in_vivo'])} / {number(published['down_in_vivo'])} "
            f"(<code>constants.SUPP5{panel[-1]}_DE_COUNTS</code>)"
        )
    else:
        lede = (
            f"{number(stats['up'])} genes up and {number(stats['down'])} down in group A, "
            f"from the full {number(universe['toxodb_65_genes'])}-gene universe. The paper "
            "quotes no counts for this panel."
        )
        counts = (
            f"{number(stats['up'])} up / {number(stats['down'])} down in group A "
            "(<code>constants.SUPP5_DE_COUNTS_REPRODUCED</code>); the caption quotes none"
        )

    g1 = (
        "<br>G1 only — <code>obs.cc_phase</code> in "
        + ", ".join(f"<code>{phase}</code>" for phase in K.G1_PHASES)
        + " (not <code>transferred_cc_phase</code>, which does not reproduce the counts)"
        if comparison.g1_only
        else ""
    )
    keys = ["other", "called out", "ribosomal protein"]
    if panel == "5E":
        keys.append("cyst wall protein")

    meta = [
        (
            "Cells",
            f"group A — {esc(comparison.label_a)}, <code>obs.orig_ident == "
            f"'{esc(comparison.sample_a)}'</code> ({number(cells['a'])} cells)<br>"
            f"group B (reference) — {esc(comparison.label_b)}, <code>obs.orig_ident == "
            f"'{esc(comparison.sample_b)}'</code> ({number(cells['b'])} cells)" + g1,
        ),
        (
            "Gene universe",
            f"{number(universe['toxodb_65_genes'])} ToxoDB-65 genes, all of them in "
            "<code>logcounts.mtx.gz</code> — the "
            f"{number(universe['deposited_genes'])} of the analysis object and the "
            f"{number(universe['extra_genes'])} recovered after them "
            "(<code>var.in_analysis_object</code>); "
            f"{number(stats['genes_detected_in_both_groups'])} of them are detected in both "
            "groups and tested",
        ),
        (
            "Test",
            "<code>sc.tl.rank_genes_groups(method='wilcoxon', tie_correct=False)</code> against "
            "group B, Benjamini-Hochberg; significant at <code>pvals_adj</code> &lt; "
            f"{K.DE_ALPHA:g} (<code>constants.DE_ALPHA</code>) and |log2FC| &gt; "
            f"{K.DE_LOG2FC_CUTOFF:g} (<code>constants.DE_LOG2FC_CUTOFF</code>)",
        ),
        ("Counts", counts),
        (
            "Plotted",
            f"{number(plotted)} significant genes; log2FC {signed(stats['log2fc_range'][0])} … "
            f"{signed(stats['log2fc_range'][1])}",
        ),
        (
            "Axes",
            f"x {signed(axes['xlim'][0], 1, trim=True)} … "
            f"{signed(axes['xlim'][1], 1, trim=True)}; y adjusted p on a "
            f"reversed log scale to 1e{MINUS}{axes['ymax']} "
            "(<code>constants.VOLCANO_AXES</code>, measured off the published figure)",
        ),
        (
            "Colour",
            f"{colour_key(keys)}pink marks the genes the panel calls out — |log2FC| beyond "
            f"{signed(low, 0)} / {signed(high, 0)} "
            "(<code>constants.VOLCANO_LABEL_RANGE</code>) — not every hypothetical protein, "
            "which is what the published legend says; colouring those would turn half the panel "
            "pink. Ribosomal proteins among them are grey.",
        ),
        VOLCANO_NOTES.get(panel, ("", "")),
        (
            "Drawn by",
            f"<code>bzfig.panels.supplementary_5{panel[-1].lower()}</code> over "
            "<code>bzfig.de.volcano_table</code>",
        ),
        (
            "Method",
            "the full recipe is in <code>docs/reproducibility.md</code>",
        ),
    ]

    return {
        "id": VOLCANO_IDS[panel],
        "overlay": " · ".join(
            (
                f"{number(cells['a'])} vs {number(cells['b'])} cells",
                f"{number(universe['toxodb_65_genes'])}-gene universe",
                f"Wilcoxon, padj < {K.DE_ALPHA:g}, |log2FC| > {K.DE_LOG2FC_CUTOFF:g}",
                f"{number(plotted)} genes plotted",
            )
        ),
        "group": "Supplementary 5",
        "label": f"Supplementary {panel}",
        "title": title[0].upper() + title[1:],
        "lede": lede,
        "plates": [{"file": f"Supplementary_{panel}_volcano", "caption": title}],
        "meta": [row for row in meta if row[0]],
    }


# ------------------------------------------------------------------ rendering

DARK_TOKENS = """
    --bg: #14161a;
    --fg: #e9e7e3;
    --muted: #9aa0a8;
    --faint: #7b818a;
    --line: #2b2f36;
    --line-soft: #23262c;
    --card: #191c21;
    --hover: #20242b;
    --sel: #232830;
    --accent: #8fb6d8;
    --shadow: 0 1px 2px rgba(0,0,0,.5), 0 10px 28px rgba(0,0,0,.42);
    --scrim: rgba(8,9,11,.86);
"""

CSS = """
:root {
  --bg: #f5f4f1;
  --fg: #1a1c1f;
  --muted: #5c6169;
  --faint: #7c828b;
  --line: #dcd8d1;
  --line-soft: #e7e3dc;
  --card: #fbfaf8;
  --hover: #ecebe6;
  --sel: #e5e3dc;
  --accent: #2f5d80;
  --plate: #ffffff;
  --shadow: 0 1px 2px rgba(30,28,24,.07), 0 8px 24px rgba(30,28,24,.06);
  --scrim: rgba(24,23,21,.82);
  --serif: Charter, "Bitstream Charter", "Iowan Old Style", "Palatino Linotype", Palatino,
           Georgia, "Times New Roman", serif;
  --sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, "Liberation Mono", monospace;
  --ok: #3f7d5f;
  --warn: #9a7430;
  --info: #3f6d96;
  --neutral: #7c828b;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) { /*DARK*/ }
}
:root[data-theme="dark"] { /*DARK*/ }

* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--fg);
  font: 16px/1.6 var(--serif);
  overflow-x: hidden;
}
a { color: inherit; }
code { font-family: var(--mono); font-size: .86em; }
kbd {
  font: 11px/1.4 var(--sans); color: var(--fg); background: var(--card);
  border: 1px solid var(--line); border-bottom-width: 2px; border-radius: 3px;
  padding: 1px 5px; white-space: nowrap;
}
code { color: var(--fg); background: var(--sel); padding: .1em .3em; border-radius: 2px; }
em { font-style: italic; }
img { max-width: 100%; }

.layout { display: grid; grid-template-columns: 268px minmax(0, 1fr); align-items: start; }

/* ------------------------------------------------------------------ nav */
.topbar { display: none; }
nav {
  position: sticky; top: 0; align-self: start;
  height: 100vh; height: 100dvh; overflow-y: auto;
  padding: 30px 16px 40px 24px;
  border-right: 1px solid var(--line);
  font: 13px/1.4 var(--sans);
}
.nav-head {
  font-size: 11px; letter-spacing: .1em; text-transform: uppercase;
  color: var(--faint); margin-bottom: 14px;
}
.nav-group { margin-bottom: 18px; }
.nav-group h3 {
  margin: 0 0 6px; padding: 0 8px;
  font: 600 11px/1.4 var(--sans); letter-spacing: .07em; text-transform: uppercase;
  color: var(--faint);
}
nav a {
  display: block; padding: 5px 8px; margin-bottom: 1px; border-radius: 3px;
  color: var(--muted); text-decoration: none;
  box-shadow: inset 2px 0 0 transparent;
}
nav a:hover { background: var(--hover); color: var(--fg); }
nav a.current { background: var(--sel); color: var(--fg); box-shadow: inset 2px 0 0 var(--accent); }
nav a.is-muted { color: var(--faint); }
.nav-label { display: block; font-weight: 600; font-size: 12px; }
.nav-title { display: block; font-size: 12px; color: var(--faint); }
nav a.current .nav-title { color: var(--muted); }
.nav-foot { padding: 12px 8px 0; border-top: 1px solid var(--line-soft); }
.ghost {
  font: 12px/1.4 var(--sans); color: var(--muted);
  background: none; border: 1px solid var(--line); border-radius: 3px;
  padding: 4px 9px; cursor: pointer;
}
.ghost:hover { color: var(--fg); border-color: var(--muted); }

/* ----------------------------------------------------------------- main */
main { min-width: 0; padding: 48px 40px 120px; }
.wrap { max-width: 860px; margin: 0 auto; }
.intro h1 { font-size: 30px; line-height: 1.25; font-weight: 600; margin: 0 0 6px; letter-spacing: -.01em; }
.intro .sub { color: var(--muted); font-size: 15px; margin: 0 0 24px; }
.intro p { margin: 0 0 14px; }
.intro .fine { font: 13px/1.6 var(--sans); color: var(--muted); }
.rule { border: 0; border-top: 1px solid var(--line); margin: 40px 0; }

.panel { scroll-margin-top: 24px; padding-top: 28px; }
.panel + .panel { border-top: 1px solid var(--line-soft); margin-top: 44px; }
.panel-head { display: flex; gap: 16px; align-items: baseline; justify-content: space-between; flex-wrap: wrap; }
.panel h2 { margin: 0; font-size: 20px; font-weight: 600; letter-spacing: -.005em; }
.panel h2 .plabel {
  display: block; font: 600 11px/1.4 var(--sans); letter-spacing: .1em;
  text-transform: uppercase; color: var(--accent); margin-bottom: 3px;
}
.lede { color: var(--muted); font-size: 15px; margin: 8px 0 0; max-width: 62ch; }

/* --------------------------------------------------------------- plates */
.plates { margin: 22px 0 0; display: grid; gap: 22px; }
.plates[data-layout="pair"] { grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); }
.plate {
  position: relative; margin: 0; padding: 16px;
  background: var(--plate); border-radius: 3px; box-shadow: var(--shadow);
  cursor: zoom-in; overflow: hidden;
}
.plate img { display: block; width: 100%; height: auto; }
.plate:focus { outline: none; }
.plate:focus-visible { outline: 2px solid var(--accent); outline-offset: 3px; }
.overlay {
  position: absolute; left: 0; right: 0; bottom: 0;
  background: rgba(255,255,255,.95);
  border-top: 1px solid #e2ded7;
  padding: 9px 14px 10px;
  font: 12px/1.5 var(--sans); color: #33363b;
  display: flex; gap: 14px; align-items: baseline; justify-content: space-between;
  opacity: 0; transform: translateY(8px);
  transition: opacity .16s ease, transform .16s ease;
  pointer-events: none;
}
.overlay .o-main { min-width: 0; }
.overlay .o-caption { font-weight: 600; color: #1a1c1f; }
.overlay .o-facts { color: #5c6169; }
.overlay .o-hint { flex: none; color: #7c828b; white-space: nowrap; }
.plate:hover .overlay, .plate:focus-visible .overlay { opacity: 1; transform: none; }
@media (hover: none) { .overlay { display: none; } }

/* ------------------------------------------------------------ switcher */
.switcher { margin-top: 18px; font: 13px/1.5 var(--sans); }
.switcher .sw-hint { color: var(--muted); margin-bottom: 8px; }
.sw-row { display: flex; flex-wrap: wrap; align-items: center; gap: 8px; margin-bottom: 6px; }
.sw-group {
  font: 600 11px/1.4 var(--sans); letter-spacing: .06em; text-transform: uppercase;
  color: var(--faint); margin-right: 2px;
}
.sw {
  font: 13px/1.4 var(--sans); font-style: italic;
  background: none; color: var(--muted);
  border: 1px solid var(--line); border-radius: 3px; padding: 4px 11px; cursor: pointer;
}
.sw:hover { color: var(--fg); border-color: var(--muted); }
.sw[aria-pressed="true"] { color: var(--fg); background: var(--sel); border-color: var(--muted); font-weight: 600; }

/* ------------------------------------------------------------ metadata */
details.meta { margin-top: 16px; }
details.meta > summary {
  font: 12px/1.5 var(--sans); color: var(--muted);
  cursor: pointer; list-style: none; display: inline-flex; align-items: center; gap: 7px;
  padding: 3px 0;
}
details.meta > summary::-webkit-details-marker { display: none; }
details.meta > summary::before {
  content: ""; width: 0; height: 0;
  border-left: 4px solid currentColor; border-top: 3.5px solid transparent; border-bottom: 3.5px solid transparent;
  transition: transform .15s ease;
}
details.meta[open] > summary::before { transform: rotate(90deg); }
details.meta > summary:hover { color: var(--fg); }
dl.meta-list {
  margin: 10px 0 0; padding: 14px 16px;
  background: var(--card); border: 1px solid var(--line-soft); border-radius: 3px;
  display: grid; grid-template-columns: minmax(120px, 170px) minmax(0, 1fr); gap: 7px 20px;
  font: 13px/1.55 var(--sans);
}
dl.meta-list dt { color: var(--faint); font-weight: 600; }
dl.meta-list dd { margin: 0; color: var(--muted); min-width: 0; overflow-wrap: anywhere; }
dl.meta-list dd code { background: none; padding: 0; color: var(--fg); }
.swatches { display: inline-flex; gap: 3px; vertical-align: -2px; margin-right: 6px; }
.swatch {
  width: 11px; height: 11px; border-radius: 2px; display: inline-block; flex: none;
  box-shadow: inset 0 0 0 1px rgba(128,128,128,.4);
}
.keys { display: flex; flex-wrap: wrap; gap: 3px 14px; margin-bottom: 5px; }
.key { display: inline-flex; align-items: center; gap: 6px; font-size: 12px; }
.scale { display: inline-flex; align-items: center; gap: 6px; margin-right: 8px; }
.scale-bar { width: 96px; height: 9px; border-radius: 2px; display: block; }
.scale-end { font-size: 11px; color: var(--faint); font-variant-numeric: tabular-nums; }
.scale-label { font-size: 12px; color: var(--muted); }

/* Long gene and constant names must break rather than widen the page. */
.lede { overflow-wrap: break-word; }
code { overflow-wrap: break-word; }
footer { margin-top: 60px; padding-top: 18px; border-top: 1px solid var(--line); font: 12px/1.7 var(--sans); color: var(--faint); }

/* --------------------------------------------------------- lightbox */
#lightbox { position: fixed; inset: 0; z-index: 100; }
#lightbox[hidden] { display: none; }
.lb-scrim { position: absolute; inset: 0; background: var(--scrim); }
.lb-body {
  position: absolute; inset: 0; display: grid; grid-template-columns: minmax(0, 1fr) 330px;
  gap: 20px; padding: 24px 24px 24px 28px;
}
.lb-stage {
  min-width: 0; display: flex; align-items: center; justify-content: center;
  overflow: auto; background: #fff; border-radius: 3px; padding: 18px;
}
.lb-stage img { max-width: 100%; max-height: calc(100vh - 100px); width: auto; height: auto; display: block; cursor: zoom-in; }
.lb-stage.zoomed { align-items: flex-start; justify-content: flex-start; }
.lb-stage.zoomed img { max-width: none; max-height: none; width: auto; cursor: zoom-out; }
.lb-side { min-width: 0; overflow-y: auto; color: #fff; font: 13px/1.6 var(--sans); }
.lb-side .lb-plabel { font: 600 11px/1.4 var(--sans); letter-spacing: .1em; text-transform: uppercase; color: #b9c6d3; }
.lb-side h3 { margin: 4px 0 8px; font: 600 17px/1.35 var(--serif); color: #fff; }
.lb-caption { margin: 12px 0 4px; color: #e6e8ea; }
.lb-side dl.meta-list { background: rgba(255,255,255,.06); border-color: rgba(255,255,255,.14); color: #d6dade; grid-template-columns: minmax(0,1fr); gap: 2px 0; }
.lb-side dl.meta-list dt { color: #9fb0c0; margin-top: 9px; }
.lb-side dl.meta-list dd { color: #d6dade; }
.lb-side dl.meta-list dd code, .lb-side code { color: #fff; background: none; }
.lb-controls { position: absolute; top: 24px; right: 24px; display: flex; gap: 8px; }
.lb-controls .ghost, .lb-nav .ghost {
  color: #eef1f4; border-color: rgba(255,255,255,.34); background: rgba(22,24,28,.86);
}
.lb-nav .ghost { font-size: 14px; padding: 3px 13px; }
.lb-controls .ghost:hover, .lb-nav .ghost:hover { color: #fff; border-color: rgba(255,255,255,.7); }
.lb-nav { position: absolute; left: 28px; bottom: 24px; display: flex; gap: 8px; }
.lb-count {
  font: 12px/1.5 var(--sans); color: #cfd5db; background: rgba(22,24,28,.86);
  padding: 5px 9px; border-radius: 3px; font-variant-numeric: tabular-nums;
}

/* --------------------------------------------------------- responsive */
@media (max-width: 1180px) {
  .lb-body { grid-template-columns: minmax(0, 1fr) 280px; }
}
@media (max-width: 900px) {
  .layout { grid-template-columns: minmax(0, 1fr); }
  .topbar {
    display: flex; align-items: center; gap: 12px;
    position: sticky; top: 0; z-index: 50;
    padding: 8px 14px; background: var(--bg); border-bottom: 1px solid var(--line);
  }
  .topbar-title { font: 12px/1.4 var(--sans); color: var(--muted); overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  nav {
    display: none; position: fixed; left: 0; right: 0; top: 41px; bottom: 0;
    height: auto; z-index: 49; background: var(--bg); border-right: 0;
  }
  body[data-nav="open"] nav { display: block; }
  main { padding: 26px 16px 80px; }
  .panel { scroll-margin-top: 52px; }
  .lb-body { grid-template-columns: minmax(0, 1fr); grid-template-rows: minmax(0, 1fr) auto; padding: 56px 14px 14px; gap: 12px; }
  .lb-side { max-height: 34vh; background: #14161a; border-radius: 3px; padding: 12px 14px; }
  .lb-stage img { max-height: 100%; }
  .lb-nav { left: 14px; top: 14px; bottom: auto; }
  .lb-controls { top: 14px; right: 14px; }
  dl.meta-list { grid-template-columns: minmax(0, 1fr); gap: 2px 0; }
  dl.meta-list dt { margin-top: 10px; }
}
@media (prefers-reduced-motion: reduce) {
  * { transition: none !important; scroll-behavior: auto !important; }
}
"""

JS = """
(function () {
  "use strict";

  /* ---------------------------------------------------------- theme */
  var root = document.documentElement;
  var themeButton = document.getElementById("theme-toggle");
  var modes = ["auto", "light", "dark"];
  function readTheme() {
    try { return localStorage.getItem("bzfig-theme"); } catch (e) { return null; }
  }
  function applyTheme(mode) {
    root.setAttribute("data-theme", mode);
    if (themeButton) themeButton.textContent = "Theme: " + mode;
    try { localStorage.setItem("bzfig-theme", mode); } catch (e) { /* private mode */ }
  }
  // Only stamp the root when the reader has actually chosen a theme here. Left
  // alone, whatever stamped it stays -- the host's setting when this page is
  // embedded, and the file's own data-theme="auto" when it is opened directly.
  var stored = readTheme();
  if (stored) {
    applyTheme(stored);
  } else if (themeButton) {
    themeButton.textContent = "Theme: " + (root.getAttribute("data-theme") || "auto");
  }
  if (themeButton) {
    themeButton.addEventListener("click", function () {
      applyTheme(modes[(modes.indexOf(root.getAttribute("data-theme")) + 1) % modes.length]);
    });
  }

  /* ------------------------------------------------------- mobile nav */
  var navToggle = document.getElementById("nav-toggle");
  function closeNav() {
    document.body.removeAttribute("data-nav");
    if (navToggle) navToggle.setAttribute("aria-expanded", "false");
  }
  if (navToggle) {
    navToggle.addEventListener("click", function () {
      var open = document.body.getAttribute("data-nav") === "open";
      if (open) { closeNav(); } else {
        document.body.setAttribute("data-nav", "open");
        navToggle.setAttribute("aria-expanded", "true");
      }
    });
  }

  /* --------------------------------------------------- nav highlighting */
  var links = Array.prototype.slice.call(document.querySelectorAll("nav a[data-target]"));
  var sections = links
    .map(function (a) { return document.getElementById(a.getAttribute("data-target")); })
    .filter(Boolean);
  var linkFor = {};
  links.forEach(function (a) { linkFor[a.getAttribute("data-target")] = a; });
  var current = null;

  links.forEach(function (a) { a.addEventListener("click", closeNav); });

  function setCurrent(id) {
    if (id === current) return;
    current = id;
    links.forEach(function (a) { a.classList.remove("current"); });
    var link = linkFor[id];
    if (link) {
      link.classList.add("current");
      if (link.offsetParent && link.parentNode && link.parentNode.parentNode) {
        var nav = document.getElementById("nav");
        var top = link.offsetTop, bottom = top + link.offsetHeight;
        if (top < nav.scrollTop || bottom > nav.scrollTop + nav.clientHeight) {
          nav.scrollTop = top - nav.clientHeight / 2;
        }
      }
    }
  }

  function refreshCurrent() {
    var best = null, bestTop = -Infinity;
    var line = window.innerHeight * 0.3;
    sections.forEach(function (section) {
      var top = section.getBoundingClientRect().top;
      if (top <= line && top > bestTop) { bestTop = top; best = section; }
    });
    if (!best && sections.length) best = sections[0];
    if (best) setCurrent(best.id);
  }
  var ticking = false;
  window.addEventListener("scroll", function () {
    if (ticking) return;
    ticking = true;
    window.requestAnimationFrame(function () { ticking = false; refreshCurrent(); });
  }, { passive: true });
  window.addEventListener("resize", refreshCurrent);
  refreshCurrent();

  /* ------------------------------------------------------ gene switcher */
  Array.prototype.forEach.call(document.querySelectorAll("[data-switcher]"), function (group) {
    var panel = group.closest(".panel");
    var buttons = Array.prototype.slice.call(group.querySelectorAll(".sw"));
    var plates = Array.prototype.slice.call(panel.querySelectorAll(".plate[data-variant]"));
    function select(key, focus) {
      buttons.forEach(function (button) {
        var on = button.getAttribute("data-variant") === key;
        button.setAttribute("aria-pressed", on ? "true" : "false");
        if (on && focus) button.focus();
      });
      plates.forEach(function (plate) {
        plate.hidden = plate.getAttribute("data-variant") !== key;
      });
    }
    buttons.forEach(function (button, index) {
      button.addEventListener("click", function () { select(button.getAttribute("data-variant"), false); });
      button.addEventListener("keydown", function (event) {
        var step = event.key === "ArrowRight" ? 1 : event.key === "ArrowLeft" ? -1 : 0;
        if (!step) return;
        event.preventDefault();
        event.stopPropagation();
        var next = buttons[(index + step + buttons.length) % buttons.length];
        select(next.getAttribute("data-variant"), true);
      });
    });
  });

  /* ----------------------------------------------------------- lightbox */
  var lightbox = document.getElementById("lightbox");
  var stage = document.getElementById("lb-stage");
  var lbImage = document.getElementById("lb-img");
  var lbLabel = document.getElementById("lb-label");
  var lbTitle = document.getElementById("lb-title");
  var lbCaption = document.getElementById("lb-caption");
  var lbMeta = document.getElementById("lb-meta");
  var lbCount = document.getElementById("lb-count");
  var order = [];
  var position = -1;
  var lastFocus = null;

  function visiblePlates() {
    return Array.prototype.slice.call(document.querySelectorAll(".plate")).filter(function (plate) {
      return !plate.hidden;
    });
  }

  function show(index) {
    position = (index + order.length) % order.length;
    var plate = order[position];
    var panel = plate.closest(".panel");
    var image = plate.querySelector("img");
    stage.classList.remove("zoomed");
    lbImage.src = image.src;
    lbImage.alt = image.alt;
    lbLabel.textContent = panel.getAttribute("data-label") || "";
    lbTitle.textContent = panel.getAttribute("data-title") || "";
    lbCaption.textContent = plate.getAttribute("data-caption") || "";
    var meta = panel.querySelector("dl.meta-list");
    lbMeta.innerHTML = meta ? meta.outerHTML : "";
    lbCount.textContent = (position + 1) + " / " + order.length;
  }

  function open(plate) {
    order = visiblePlates();
    var index = order.indexOf(plate);
    if (index < 0) return;
    lastFocus = document.activeElement;
    lightbox.hidden = false;
    document.body.style.overflow = "hidden";
    show(index);
    document.getElementById("lb-close").focus();
  }

  function close() {
    lightbox.hidden = true;
    lbImage.removeAttribute("src");
    document.body.style.overflow = "";
    if (lastFocus && lastFocus.focus) lastFocus.focus();
  }

  Array.prototype.forEach.call(document.querySelectorAll(".plate"), function (plate) {
    plate.addEventListener("click", function () { open(plate); });
    plate.addEventListener("keydown", function (event) {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        open(plate);
      }
    });
  });

  document.getElementById("lb-close").addEventListener("click", close);
  document.querySelector(".lb-scrim").addEventListener("click", close);
  document.getElementById("lb-prev").addEventListener("click", function () { show(position - 1); });
  document.getElementById("lb-next").addEventListener("click", function () { show(position + 1); });
  stage.addEventListener("click", function (event) {
    if (event.target === lbImage) stage.classList.toggle("zoomed");
  });

  /* ---------------------------------------------------------- keyboard */
  function typing(event) {
    var tag = (event.target.tagName || "").toLowerCase();
    return tag === "input" || tag === "textarea" || tag === "select" || event.target.isContentEditable;
  }

  function goToPanel(step) {
    if (!sections.length) return;
    var index = sections.findIndex(function (section) { return section.id === current; });
    if (index < 0) index = 0;
    var next = sections[Math.min(Math.max(index + step, 0), sections.length - 1)];
    next.scrollIntoView({ behavior: "smooth", block: "start" });
    setCurrent(next.id);
  }

  document.addEventListener("keydown", function (event) {
    if (event.metaKey || event.ctrlKey || event.altKey || typing(event)) return;
    if (!lightbox.hidden) {
      if (event.key === "Escape") { event.preventDefault(); close(); }
      else if (event.key === "ArrowRight") { event.preventDefault(); show(position + 1); }
      else if (event.key === "ArrowLeft") { event.preventDefault(); show(position - 1); }
      else if (event.key === "Tab") {
        var focusable = lightbox.querySelectorAll("button");
        var first = focusable[0], last = focusable[focusable.length - 1];
        if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last.focus(); }
        else if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first.focus(); }
      }
      return;
    }
    if (event.key === "Escape" && document.body.getAttribute("data-nav") === "open") {
      closeNav();
      return;
    }
    if (event.target.closest && event.target.closest(".switcher")) return;
    if (event.key === "ArrowRight") { event.preventDefault(); goToPanel(1); }
    else if (event.key === "ArrowLeft") { event.preventDefault(); goToPanel(-1); }
  });
})();
"""

TEMPLATE = """<!DOCTYPE html>
<html lang="en" data-theme="auto">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="robots" content="noindex">
<title>{{ page_title }}</title>
<style>{{ css }}</style>
</head>
<body>

<div class="topbar">
  <button id="nav-toggle" class="ghost" aria-expanded="false" aria-controls="nav">Panels</button>
  <span class="topbar-title">Reproduced figure panels</span>
</div>

<div class="layout">

<nav id="nav" aria-label="Panels">
  <div class="nav-head">Panels</div>
  {%- for group in nav_groups %}
  <div class="nav-group">
    <h3>{{ group.name }}</h3>
    {%- for item in group.entries %}
    <a href="#{{ item.href }}"{% if item.target %} data-target="{{ item.target }}"{% endif %}
       class="{{ 'is-muted' if item.muted else '' }}">
      <span class="nav-label">{{ item.label }}</span>
      <span class="nav-title">{{ item.title }}</span>
    </a>
    {%- endfor %}
  </div>
  {%- endfor %}
  <div class="nav-foot">
    <button id="theme-toggle" class="ghost" type="button">Theme: auto</button>
  </div>
</nav>

<main>
<div class="wrap">

<header class="intro" id="top">
  <h1>{{ page_title }}</h1>
  <p class="sub">Every figure panel regenerated from the deposited single-cell data.</p>
  <p>The rendered output of <code>scripts/make_figures.py</code> in the
    <code>BradyzoiteHeterogeneity-figures</code> repository, one panel per entry, each with the
    cells, genes, colour scale and code it was drawn from.</p>
  <p>These are the individual panels as the analysis produced them. Panel letters, the Figure 1E
    grid and its row labels, the single shared colour bar across 1E and 1F, the in-plot cluster
    numbers in Figure 1B and the coloured row-group bands beside Figure 1C were all added during
    figure assembly in a vector editor, and are not reproduced here.</p>
  <p class="fine">Hover a panel for its metadata; click it to enlarge, with the full metadata
    beside it. <kbd>←</kbd> and <kbd>→</kbd> move between panels, and between images inside the
    lightbox; <kbd>Esc</kbd> closes it. Figure 1E and Supplementary 1B and 1C put their genes on
    one control each.</p>
  <p class="fine">The embedded images are screen-resolution copies ({{ max_edge }} px on the long
    edge) of the 300 dpi renders. The full-resolution PNG, SVG and PDF for every panel are in
    <code>figures/</code>; the values behind them are in <code>data/</code>, and
    <code>data/MANIFEST.json</code> records the SHA-256 of each input file. This page is a single
    file and fetches nothing.</p>
</header>

<hr class="rule">

{%- for panel in panels %}
<section class="panel" id="{{ panel.id }}"
         data-label="{{ panel.label }}"
         data-title="{{ panel.title }}">
  <div class="panel-head">
    <h2><span class="plabel">{{ panel.label }}</span>{{ panel.title }}</h2>
  </div>
  {%- if panel.lede %}
  <p class="lede">{{ panel.lede }}</p>
  {%- endif %}

  {%- if panel.switcher %}
  <div class="switcher" data-switcher role="group" aria-label="{{ panel.label }} gene">
    <div class="sw-hint">{{ panel.switch_hint }}</div>
    {%- for row in panel.switcher %}
    <div class="sw-row">
      <span class="sw-group">{{ row.group }}</span>
      {%- for item in row.entries %}
      <button class="sw" type="button" data-variant="{{ item.key }}"
              aria-pressed="{{ 'true' if item.active else 'false' }}">{{ item.label }}</button>
      {%- endfor %}
    </div>
    {%- endfor %}
  </div>
  {%- endif %}

  <div class="plates"{% if panel.plate_layout %} data-layout="{{ panel.plate_layout }}"{% endif %}>
    {%- for plate in panel.plates %}
    <figure class="plate" tabindex="0" role="button"
            aria-label="Enlarge: {{ plate.caption }}"
            data-caption="{{ plate.caption }}"
            {%- if plate.variant %} data-variant="{{ plate.variant }}"{% endif %}
            {%- if plate.hidden %} hidden{% endif %}>
      <img src="{{ plate.uri }}" width="{{ plate.width }}" height="{{ plate.height }}"
           alt="{{ panel.label }}: {{ plate.caption }}" decoding="async">
      <figcaption class="overlay">
        <span class="o-main">
          <span class="o-caption">{{ plate.caption }}</span>
          <span class="o-facts"> — {{ plate.facts }}</span>
        </span>
        <span class="o-hint">Click to enlarge</span>
      </figcaption>
    </figure>
    {%- endfor %}
  </div>

  <details class="meta">
    <summary>Panel metadata</summary>
    <dl class="meta-list">
      {%- for term, value in panel.meta %}
      <dt>{{ term }}</dt><dd>{{ value }}</dd>
      {%- endfor %}
      <dt>Rendered files</dt><dd>{{ panel.files }}</dd>
    </dl>
  </details>
</section>
{%- endfor %}

<footer>
  <p>{{ panel_count }} cards covering {{ image_count }} rendered panels, {{ page_size }} in one
    file. Generated by <code>scripts/make_preview.py</code> from <code>figures/</code>,
    <code>src/bzfig/constants.py</code>, <code>docs/reproducibility.md</code> and
    <code>data/</code>. No external resources: every image is embedded, and the page makes no
    network request.</p>
  <p>Repository: <code>BradyzoiteHeterogeneity-figures</code> — run
    <code>python scripts/make_figures.py</code> to re-render the panels and
    <code>python scripts/make_preview.py</code> to rebuild this page.</p>
</footer>

</div>
</main>
</div>

<div id="lightbox" hidden role="dialog" aria-modal="true" aria-label="Enlarged panel">
  <div class="lb-scrim"></div>
  <div class="lb-body">
    <div class="lb-stage" id="lb-stage"><img id="lb-img" alt=""></div>
    <aside class="lb-side">
      <div class="lb-plabel" id="lb-label"></div>
      <h3 id="lb-title"></h3>
      <p class="lb-caption" id="lb-caption"></p>
      <div id="lb-meta"></div>
    </aside>
    <div class="lb-controls">
      <button class="ghost" id="lb-close" type="button">Close (Esc)</button>
    </div>
    <div class="lb-nav">
      <button class="ghost" id="lb-prev" type="button" aria-label="Previous image">←</button>
      <button class="ghost" id="lb-next" type="button" aria-label="Next image">→</button>
      <span class="lb-count" id="lb-count"></span>
    </div>
  </div>
</div>

<script>{{ js }}</script>
</body>
</html>
"""


def human(size: int) -> str:
    return f"{size / 1_000_000:.1f} MB" if size >= 1_000_000 else f"{size / 1_000:.0f} kB"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures", type=Path, default=REPO / "figures")
    parser.add_argument("--data", type=Path, default=REPO / "data")
    parser.add_argument("--docs", type=Path, default=REPO / "docs")
    parser.add_argument("--out", type=Path, default=REPO / "preview" / "index.html")
    parser.add_argument("--max-edge", type=int, default=MAX_EDGE)
    parser.add_argument("--quality", type=int, default=QUALITY)
    args = parser.parse_args()

    import jinja2
    from markupsafe import Markup

    repro = split_sections((args.docs / "reproducibility.md").read_text())
    verdicts = {
        row["Panel"]: row["What it shows"]
        for row in find_table(
            doc_section(repro, "Panels drawn directly"), "Panel", "What it shows"
        )
    }

    facts = dataset_facts(args.data)
    panels = build_panels(facts, verdicts, args.figures)

    total_bytes = 0
    image_count = 0
    for panel in panels:
        overlay = panel["overlay"]

        if "variants" in panel:
            plates = []
            rows: dict[str, list] = {}
            for index, variant in enumerate(panel["variants"]):
                plates.append(
                    {
                        "file": variant["file"],
                        "caption": variant["caption"],
                        "variant": variant["key"],
                        "hidden": index > 0,
                    }
                )
                rows.setdefault(variant["row_group"], []).append(
                    {"key": variant["key"], "label": variant["label"], "active": index == 0}
                )
            panel["plates"] = plates
            panel["switcher"] = [
                {"group": group, "entries": entries} for group, entries in rows.items()
            ]

        for plate in panel["plates"]:
            source = args.figures / f"{plate['file']}.png"
            if not source.exists():
                raise SystemExit(f"make_preview: missing render {source}")
            encoded = encode_image(source, args.max_edge, args.quality)
            plate.update(encoded)
            plate["facts"] = overlay
            total_bytes += encoded["bytes"]
            image_count += 1

        panel["files"] = Markup(
            ", ".join(f"<code>figures/{esc(p['file'])}.png</code>" for p in panel["plates"])
            + " (also <code>.svg</code> and <code>.pdf</code>)"
        )
        panel["meta"] = [(term, Markup(value)) for term, value in panel["meta"]]

    nav_groups: list[dict] = []
    for panel in panels:
        if not nav_groups or nav_groups[-1]["name"] != panel["group"]:
            nav_groups.append({"name": panel["group"], "entries": []})
        nav_groups[-1]["entries"].append(
            {
                "href": panel["id"],
                "target": panel["id"],
                "label": panel["label"],
                "title": panel["title"],
                "muted": False,
            }
        )

    environment = jinja2.Environment(autoescape=True, trim_blocks=False, lstrip_blocks=True)
    template = environment.from_string(TEMPLATE)

    def render(page_size: str) -> str:
        return template.render(
            page_title=PAGE_TITLE,
            css=Markup(CSS.replace("/*DARK*/", DARK_TOKENS)),
            js=Markup(JS),
            panels=panels,
            nav_groups=nav_groups,
            panel_count=len(panels),
            image_count=image_count,
            page_size=page_size,
            max_edge=args.max_edge,
        )

    # The footer quotes the size of the finished file: render once to measure it,
    # then again with the measurement. The second pass changes only that string,
    # so its own size is stable to the tenth of a megabyte the footer prints.
    page = render(human(len(render("0.0 MB").encode("utf-8"))))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(page, encoding="utf-8")
    out = args.out.resolve()
    shown = out.relative_to(REPO) if out.is_relative_to(REPO) else out
    print(f"{shown}: {human(out.stat().st_size)}")
    print(f"  {len(panels)} cards, {image_count} panels ({human(total_bytes)} of WebP)")
    print(f"  images resized to {args.max_edge} px on the long edge, quality {args.quality}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
