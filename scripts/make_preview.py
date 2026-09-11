#!/usr/bin/env python
"""Build ``preview/index.html`` — one self-contained page showing every panel.

    python scripts/make_preview.py                  # -> preview/index.html
    python scripts/make_preview.py --max-edge 1400  # smaller embedded images

The page is a review aid for the authors: every rendered panel, its metadata and
its reproducibility status, in one file that makes **no network request of any
kind**. Images are embedded as ``data:`` URIs, the stylesheet and script are
inline, and nothing is fetched at load time.

Everything the page states is read from the repository rather than retyped:

* panel images            ``figures/*.png``
* colour limits, palettes,
  gene lists, constants    ``src/bzfig/constants.py``
* cell counts             ``data/obs.csv.gz`` and ``data/MANIFEST.json``
* per-panel verdicts and
  the limitations section  ``docs/reproducibility.md``

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

PAGE_TITLE = "Bradyzoite heterogeneity — reproduced figure panels"

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
# Enough of markdown to render the prose sections of docs/reproducibility.md.
# Links are flattened to their text: the page must not offer anything to fetch.


def _inline(text: str) -> str:
    text = esc(text)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<em>\1</em>", text)
    return text


_BULLET = re.compile(r"^\s*[*-] ")
_ORDERED = re.compile(r"^\s*\d+\. ")


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


def md_to_html(markdown: str) -> str:
    lines = markdown.strip("\n").splitlines()
    out: list[str] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        if not line.strip():
            index += 1
            continue
        if line.lstrip().startswith("```"):
            index += 1
            fenced: list[str] = []
            while index < len(lines) and not lines[index].lstrip().startswith("```"):
                fenced.append(lines[index])
                index += 1
            index += 1
            out.append("<pre><code>" + esc("\n".join(fenced)) + "</code></pre>")
        elif line.startswith("#"):
            depth = len(line) - len(line.lstrip("#"))
            level = min(depth + 1, 5)
            out.append(f"<h{level}>{_inline(line.lstrip('#').strip())}</h{level}>")
            index += 1
        elif _BULLET.match(line) or _ORDERED.match(line):
            ordered = bool(_ORDERED.match(line))
            starts = _ORDERED if ordered else _BULLET
            items: list[str] = []
            while index < len(lines) and starts.match(lines[index]):
                item = lines[index][starts.match(lines[index]).end() :].strip()
                index += 1
                while (
                    index < len(lines)
                    and lines[index].startswith("  ")
                    and lines[index].strip()
                    and not _BULLET.match(lines[index])
                    and not _ORDERED.match(lines[index])
                ):
                    item += " " + lines[index].strip()
                    index += 1
                items.append(f"<li>{_inline(item)}</li>")
            tag = "ol" if ordered else "ul"
            out.append(f"<{tag}>" + "".join(items) + f"</{tag}>")
        elif line.lstrip().startswith("|"):
            block: list[str] = []
            while index < len(lines) and lines[index].lstrip().startswith("|"):
                block.append(lines[index])
                index += 1
            rows = _table_rows(block)
            head = "".join(f"<th>{_inline(cell)}</th>" for cell in rows[0])
            body = "".join(
                "<tr>" + "".join(f"<td>{_inline(cell)}</td>" for cell in row) + "</tr>"
                for row in rows[1:]
            )
            out.append(
                f'<div class="scroll-x"><table><thead><tr>{head}</tr></thead>'
                f"<tbody>{body}</tbody></table></div>"
            )
        else:
            paragraph: list[str] = []
            while (
                index < len(lines)
                and lines[index].strip()
                and not lines[index].startswith("#")
                and not lines[index].lstrip().startswith("|")
                and not _BULLET.match(lines[index])
                and not _ORDERED.match(lines[index])
            ):
                paragraph.append(lines[index].strip())
                index += 1
            out.append(f"<p>{_inline(' '.join(paragraph))}</p>")
    return "\n".join(out)


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

    genes = {row["gene_id"] for row in _rows(datadir / "var.csv.gz")}
    with (datadir / "supp1a_marker_genes.csv").open(newline="") as handle:
        markers = [row[0] for row in csv.reader(handle) if row]
    facts["supp1a_listed"] = len(markers)
    facts["supp1a_plotted"] = len(set(markers) & genes)

    # The volcano panels: the two groups of each published comparison, counted
    # the way bzfig.de._groups selects them, and what the test produced.
    facts["universe"] = manifest["extra_genes"]["universe"]
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

STATUS = {
    "exact": (
        "Reproduced exactly",
        "checked against the published panel (docs/reproducibility.md)",
    ),
    "reconstruction": (
        "Reconstruction",
        "the drawing code does not exist; the recipe was derived by matching the published panel",
    ),
    "wider": (
        "Reproduced, wider gene set",
        "the published test ran on the full ToxoDB-65 gene universe, not on the genes the "
        "deposited object kept; see “The volcano panels”",
    ),
    "constant": (
        "Recorded constant",
        "not derivable from the deposited object; the published values are drawn from constants.py",
    ),
}

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


def colour_key(names: list[str]) -> str:
    """The volcano point colours, as a small swatch-and-label key."""
    items = "".join(
        f'<span class="key"><span class="swatch" style="background:{esc(K.VOLCANO_COLORS[n])}"'
        f' title="{esc(K.VOLCANO_COLORS[n])}"></span>{esc(n)}</span>'
        for n in names
    )
    return f'<span class="keys">{items}</span>'


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


def build_panels(facts: dict, verdicts: dict[str, str]) -> list[dict]:
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
            "status": "exact",
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
                (
                    "Assembly",
                    "the in-plot cluster numbers on the published panel were added in a vector "
                    "editor and are not reproduced here",
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
            "status": "exact",
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
                (
                    "Assembly",
                    "the coloured row-group bands beside the published panel were added in a "
                    "vector editor; the band sizes were read off the printed figure",
                ),
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
                "Drawn from the six values hard-coded in the analysis notebook. It is not a "
                "reproduction — see “What does not reproduce”."
            ),
            "status": "constant",
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
                    "Re-running the obvious test",
                    f"gives <code>[37, 94, 84, 861, 18, 52]</code> on the deposited "
                    f"{number(facts['n_vars'])} genes and "
                    "<code>[45, 93, 84, 863, 18, 53]</code> on the full "
                    f"{number(facts['universe']['toxodb_65_genes'])}-gene universe now that it is "
                    "shipped — neither is the published set, so the gene universe is not the "
                    "explanation and whatever fixed these six numbers is still missing",
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
            "status": "exact",
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
                (
                    "Label check",
                    "the two enolase labels are mapped the other way round in an early notebook "
                    "cell; this repository follows the gene annotation and the published panels — "
                    "see the note below",
                ),
                (
                    "Assembly",
                    "the published figure lays the six out as a grid with row labels and a single "
                    "shared colour bar across 1E and 1F, added in a vector editor",
                ),
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
            "status": "exact",
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
            "status": "exact",
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
            "id": "supp-1a",
            "overlay": f"{heatmap_overlay}{dot}{number(facts['supp1a_plotted'])} marker genes",
            "group": "Supplementary 1",
            "label": "Supplementary 1A",
            "title": "All cluster markers per cluster",
            "lede": verdicts.get("Supplementary 1A", ""),
            "status": "exact",
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
                    "Not included",
                    "Supplementary 1B–F is not derived from this dataset; the caption attributes "
                    "it to Benke et al., mined from ToxoDB",
                ),
            ],
        }
    )

    panels.append(
        {
            "id": "supp-4",
            "overlay": dot.join(
                (
                    f"CCC {number(facts['cc_counts']['CCC'])} / "
                    f"MCC {number(facts['cc_counts']['MCC'])} cells",
                    "logcounts",
                    "per-gene z-score across the five phases",
                    f"RdBu_r {MINUS}2 to 2",
                )
            ),
            "group": "Supplementary 4",
            "label": "Supplementary 4",
            "title": "Cell-cycle regulators, common and modified cell cycle",
            "lede": (
                "A reconstruction: the code that drew this panel does not exist, and the CCC cell "
                "selection was made interactively rather than in a script."
            ),
            "status": "reconstruction",
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
                    "<code>logcounts</code> — not the scaled layer the expression panels use",
                ),
                (
                    "Values",
                    "mean expression per <code>obs.transferred_cc_phase</code>, z-scored per gene "
                    "across the five phases; rows sorted by peak phase, then peak height",
                ),
                ("Genes", f"{len(K.SUPP4_GENES)} regulators, <code>constants.SUPP4_GENES</code>"),
                ("Colour", f"{supp4_scale} vmin {MINUS}2, vmax 2, centred on 0"),
                (
                    "Worth knowing",
                    "the colour bar axis runs −4…4 but the values only span −1.77…+1.79, and the "
                    "CCC group is thinly populated, so several genes are detected in one phase "
                    "only and their rows sit at the one-hot extremes of ±1.789 / −0.447",
                ),
                ("Drawn by", "<code>bzfig.panels.supplementary_4</code>"),
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
            "status": "exact",
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
            "status": "exact",
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
            "status": "exact",
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


VOLCANO_IDS = {"5C": "supp-5c", "5E": "supp-5e", "5F": "supp-5f"}

# Caveats that belong to one panel only. Both are written up in
# docs/reproducibility.md, under "Reproduced from a wider gene set".
VOLCANO_NOTES = {
    "5C": (
        "Points off the top of the axis",
        "100 of the plotted points have an adjusted p-value below the published 1e-100 axis top "
        "— 20 of them underflow float64 to zero — and are drawn as triangles on the axis rather "
        "than dropped. 5E and 5F have none.",
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
            f"Reproduces the published {number(published['up_in_vivo'])} up / "
            f"{number(published['down_in_vivo'])} down as {number(stats['up'])} / "
            f"{number(stats['down'])}, from the full "
            f"{number(universe['toxodb_65_genes'])}-gene universe."
        )
        counts = (
            f"{number(stats['up'])} up / {number(stats['down'])} down in group A "
            "(<code>constants.SUPP5_DE_COUNTS_REPRODUCED</code>), against a published "
            f"{number(published['up_in_vivo'])} / {number(published['down_in_vivo'])} "
            f"(<code>constants.SUPP5{panel[-1]}_DE_COUNTS</code>)"
        )
    else:
        lede = (
            f"Reproduces {number(stats['up'])} up / {number(stats['down'])} down. The paper "
            "quotes no counts for this panel; it is validated by its called-out genes instead."
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
            f"{number(universe['toxodb_65_genes'])} ToxoDB-65 genes — the "
            f"{number(universe['deposited_genes'])} of <code>logcounts.mtx.gz</code> widened with "
            f"the {number(universe['extra_genes'])} of <code>logcounts_extra.mtx.gz</code> "
            "(<code>bzfig.de.expanded</code>); "
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
            'the recipe, the evidence for the cutoff and what is not reproduced are in '
            '<a href="#volcanoes">The volcano panels</a>',
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
        "status": "wider",
        "plates": [{"file": f"Supplementary_{panel}_volcano", "caption": title}],
        "meta": [row for row in meta if row[0]],
    }


def missing_panels(facts: dict) -> list[dict]:
    """The panels that are not regenerated, with their recorded values."""
    return [
        {
            "label": "Figure 1D",
            "status": "Recorded constant",
            "what": "# Unique markers/Cluster",
            "why": (
                "Hard-coded in the analysis notebook with no accompanying computation. Re-running "
                "the obvious test on the deposited object gives 37, 94, 84, 861, 18, 52, and on "
                f"the full {number(facts['universe']['toxodb_65_genes'])}-gene universe "
                "45, 93, 84, 863, 18, 53 — neither is the published set, so the gene universe is "
                "not the explanation here. The panel above is drawn from the recorded values, "
                + ", ".join(str(v) for _, v in sorted(K.UNIQUE_MARKERS_PER_CLUSTER.items()))
                + " (<code>constants.UNIQUE_MARKERS_PER_CLUSTER</code>), so it matches the paper "
                "exactly; it is simply not a reproduction."
            ),
        },
    ]


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

.badge {
  flex: none; font: 11px/1.4 var(--sans); letter-spacing: .02em; color: var(--muted);
  border: 1px solid var(--line); border-radius: 999px; padding: 3px 10px 3px 8px;
  display: inline-flex; align-items: center; gap: 6px; white-space: nowrap;
}
.badge::before { content: ""; width: 6px; height: 6px; border-radius: 50%; background: var(--neutral); }
.badge[data-status="exact"]::before { background: var(--ok); }
.badge[data-status="reconstruction"]::before { background: var(--warn); }
.badge[data-status="wider"]::before { background: var(--info); }

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

/* ------------------------------------------------------------- prose */
.doc h2 { font-size: 22px; font-weight: 600; margin: 0 0 10px; }
.doc h3 { font-size: 17px; font-weight: 600; margin: 30px 0 8px; }
.doc h4 { font-size: 15px; font-weight: 600; margin: 22px 0 6px; }
.doc p { margin: 0 0 12px; max-width: 66ch; }
.doc ul, .doc ol { margin: 0 0 14px; padding-left: 22px; max-width: 66ch; }
.doc li { margin-bottom: 6px; }
.doc .quoted { color: var(--muted); }
.scroll-x { overflow-x: auto; margin: 0 0 16px; }
pre {
  margin: 0 0 16px; padding: 12px 14px; overflow-x: auto;
  background: var(--card); border: 1px solid var(--line-soft); border-radius: 3px;
  font: 12px/1.6 var(--mono); color: var(--muted);
}
pre code { background: none; padding: 0; font-size: inherit; white-space: pre; }
table { border-collapse: collapse; font: 13px/1.5 var(--sans); width: 100%; }
th, td { text-align: left; vertical-align: top; padding: 7px 14px 7px 0; border-bottom: 1px solid var(--line-soft); }
th { color: var(--faint); font-weight: 600; white-space: nowrap; }
td { color: var(--muted); }
.missing { margin: 18px 0 30px; }
.missing-row { padding: 14px 0; border-top: 1px solid var(--line-soft); display: grid; grid-template-columns: minmax(150px, 200px) minmax(0, 1fr); gap: 6px 20px; }
.missing-row:last-child { border-bottom: 1px solid var(--line-soft); }
.missing-label { font: 600 13px/1.5 var(--sans); }
.missing-what { font: 13px/1.5 var(--sans); color: var(--faint); margin-top: 2px; }
.missing-why { font: 13px/1.6 var(--sans); color: var(--muted); }
/* Long gene and constant names must break rather than widen the page. */
.doc p, .doc li, .lede, .missing-what, .missing-why, td, th { overflow-wrap: break-word; }
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
.lb-side .badge { color: #cfd6dd; border-color: rgba(255,255,255,.3); }
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
  .missing-row { grid-template-columns: minmax(0, 1fr); }
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
    try { return localStorage.getItem("bzfig-theme") || "auto"; } catch (e) { return "auto"; }
  }
  function applyTheme(mode) {
    root.setAttribute("data-theme", mode);
    if (themeButton) themeButton.textContent = "Theme: " + mode;
    try { localStorage.setItem("bzfig-theme", mode); } catch (e) { /* private mode */ }
  }
  applyTheme(readTheme());
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
  var lbBadge = document.getElementById("lb-badge");
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
    var status = panel.getAttribute("data-status") || "";
    lbBadge.textContent = panel.getAttribute("data-status-label") || "";
    lbBadge.setAttribute("data-status", status);
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
  <div class="nav-group">
    <h3>Reference</h3>
    <a href="#volcanoes" data-target="volcanoes">
      <span class="nav-label">The volcano panels</span>
      <span class="nav-title">Supplementary 5C / 5E / 5F, and the wider gene set</span>
    </a>
    <a href="#limitations" data-target="limitations">
      <span class="nav-label">What does not reproduce</span>
      <span class="nav-title">Figure 1D</span>
    </a>
    <a href="#notes" data-target="notes">
      <span class="nav-label">Notes</span>
      <span class="nav-title">Label correction, scope, assembly</span>
    </a>
  </div>
  <div class="nav-foot">
    <button id="theme-toggle" class="ghost" type="button">Theme: auto</button>
  </div>
</nav>

<main>
<div class="wrap">

<header class="intro" id="top">
  <h1>{{ page_title }}</h1>
  <p class="sub">Every panel that could be regenerated from the deposited single-cell data,
    for review.</p>
  <p>This page is a preview of the figure panels produced by the
    <code>BradyzoiteHeterogeneity-figures</code> repository: the rendered output of
    <code>scripts/make_figures.py</code>, drawn from the deposited dataset, one panel per entry.
    Each carries the metadata it was drawn with and its status — reproduced exactly, reproduced
    from a wider gene set, a reconstruction, or a recorded constant — taken from
    <code>docs/reproducibility.md</code>.</p>
  <p>These are the individual panels as the analysis produced them. Panel letters, the Figure 1E
    grid and its row labels, the single shared colour bar across 1E and 1F, the in-plot cluster
    numbers in Figure 1B and the coloured row-group bands beside Figure 1C were all added during
    figure assembly in a vector editor, and are not reproduced here.</p>
  <p class="fine">Hover a panel for its metadata; click it to enlarge, with the full metadata
    beside it. <kbd>←</kbd> and <kbd>→</kbd> move between panels, and between images inside the
    lightbox; <kbd>Esc</kbd> closes it. Figure 1E has six genes on one control.</p>
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
         data-title="{{ panel.title }}"
         data-status="{{ panel.status }}"
         data-status-label="{{ panel.status_label }}">
  <div class="panel-head">
    <h2><span class="plabel">{{ panel.label }}</span>{{ panel.title }}</h2>
    <span class="badge" data-status="{{ panel.status }}" title="{{ panel.status_note }}">{{ panel.status_label }}</span>
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
      <dt>Status</dt><dd>{{ panel.status_label }} — {{ panel.status_note }}</dd>
    </dl>
  </details>

  {%- if panel.extra %}
  <details class="meta">
    <summary>{{ panel.extra_summary }}</summary>
    <div class="doc quoted" style="margin-top:10px">{{ panel.extra }}</div>
  </details>
  {%- endif %}
</section>
{%- endfor %}

<hr class="rule">

<section class="doc" id="volcanoes">
  <h2>The volcano panels</h2>
  <p class="lede">Supplementary 5C, 5E and 5F were written off as unreproducible when this
    repository was first put together. They are reproducible — the published test simply ran on a
    wider gene set than the deposited object holds. Quoted from
    <code>docs/reproducibility.md</code>, including what still is not reproduced.</p>
  <div class="quoted">{{ volcano_doc }}</div>
</section>

<hr class="rule">

<section class="doc" id="limitations">
  <h2>What does not reproduce</h2>
  <p class="lede">One published panel is not a reproduction. This is the honest list, summarised
    from <code>docs/reproducibility.md</code>; the full text of that section follows.</p>
  <div class="missing">
    {%- for row in missing %}
    <div class="missing-row">
      <div>
        <div class="missing-label">{{ row.label }}</div>
        <div class="missing-what">{{ row.what }}</div>
        <div class="missing-what">{{ row.status }}</div>
      </div>
      <div class="missing-why">{{ row.why }}</div>
    </div>
    {%- endfor %}
  </div>
  <div class="quoted">{{ limitations_doc }}</div>
</section>

<hr class="rule">

<section class="doc" id="notes">
  <h2>Notes</h2>
  <p class="lede">Also quoted from <code>docs/reproducibility.md</code>. Open questions for the
    authors are collected separately in <code>docs/todo-for-authors.md</code>.</p>
  <div class="quoted">{{ notes_doc }}</div>
</section>

<footer>
  <p>{{ panel_count }} cards covering {{ image_count }} rendered panels, {{ page_size }} in one
    file. Generated by
    <code>scripts/make_preview.py</code> from <code>figures/</code>,
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
      <span class="badge" id="lb-badge"></span>
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
        for row in find_table(doc_section(repro, "Reproduced exactly"), "Panel", "What it shows")
    }

    facts = dataset_facts(args.data)
    panels = build_panels(facts, verdicts)

    total_bytes = 0
    image_count = 0
    for panel in panels:
        status_label, status_note = STATUS[panel["status"]]
        panel["status_label"] = status_label
        panel["status_note"] = status_note
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

    panels[[p["id"] for p in panels].index("supp-4")].update(
        extra_summary="Why this is a reconstruction, and how closely it matches",
        extra=Markup(md_to_html(doc_section(repro, "Reproduced as a reconstruction"))),
    )

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

    limitations_doc = md_to_html(doc_section(repro, "Not reproducible"))
    volcano_doc = md_to_html(doc_section(repro, "Reproduced from a wider gene set"))
    # Everything the page has not already quoted, in the order the document has it,
    # so a section added to the docs turns up here rather than being dropped.
    quoted = (
        "Reproduced exactly",
        "Reproduced as a reconstruction",
        "Reproduced from a wider gene set",
        "Not reproducible",
    )
    notes_doc = "\n".join(
        md_to_html(text) for heading, text in repro if not heading.startswith(quoted)
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
            missing=[
                {key: Markup(value) for key, value in row.items()}
                for row in missing_panels(facts)
            ],
            volcano_doc=Markup(volcano_doc),
            limitations_doc=Markup(limitations_doc),
            notes_doc=Markup(notes_doc),
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
