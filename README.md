# Bradyzoite heterogeneity — figure reproduction

Everything needed to regenerate the figure panels derived from the single-cell
RNA-seq data: the deposited dataset, the code that draws each panel, and notes on
how each panel relates to the published figure.

```bash
uv sync                              # or: pip install -e .
python scripts/make_figures.py       # render every panel to PNG, SVG and PDF
python scripts/make_preview.py       # build preview/index.html from what was rendered
```

**You do not need git-lfs.** The input data is stored in LFS, but if your clone
left pointer stubs in `data/` instead of the files, `make_figures.py` notices and
downloads them from the public mirrors before rendering. Run
`python scripts/fetch_data.py` to do that step on its own.

Rendered panels land in `figures/`. `python scripts/make_figures.py --list` shows
the panel names — it needs no data, so it answers in a fresh clone — and
`--panel 1C 1G` renders a subset.

Figure 2H and Supplementary 4 also write full-precision companion tables to
`figures/tables/`, and carry a test suite:

```bash
python scripts/make_figures.py --panel Figure_2H Supplementary_4
python -m unittest discover -s tests -v
```

`preview/index.html` is a single self-contained page showing every rendered panel
with its caption and source — no network access, so it can be opened from disk or
served as a static file. It needs the `preview` extra
(`pip install -e '.[preview]'`); `--max-edge` controls the size of the embedded
images.

## Start here

| | |
| --- | --- |
| [`docs/reproducibility.md`](docs/reproducibility.md) | How each panel is drawn, and where it departs from the published figure |
| [`docs/figure-2h-supplementary-4.md`](docs/figure-2h-supplementary-4.md) | Figure 2H and Supplementary 4: method, outputs and validation |
| [`data/MANIFEST.json`](data/MANIFEST.json) | Every input file with its SHA-256 and download mirrors |

## Panels

Drawn from the deposited data:

* **Figure 1B** — cluster UMAP of the in vivo bradyzoites, and cells per cluster
* **Figure 1C** — selected marker genes per cluster
* **Figure 1E / 1F** — per-gene expression UMAPs, and cst1 (srs44)
* **Figure 1G** — CST1/SRS44 expression violins per cluster
* **Figure 2H** — the 5 × 29 Pearson correlation heatmap across the in vivo cells
* **Figure 3A / 3B / 3C** — each cohort by cell-cycle phase, and its phase bar
* **Figure 3D** — the in vitro cells picked out of the projection
* **Figure 3E / 3F** — the in vivo tachyzoites on the scVI integration, plain and
  by cell-cycle phase
* **Figure 6** — the same three cohorts without the grey background layer
* **Supplementary 1A** — all cluster markers per cluster
* **Supplementary 1B / 1C** — eleven per-gene expression UMAPs, micronemes and cyst wall proteins
* **Supplementary 4** — cell-cycle regulators per phase, common (CCC) and modified (MCC) cell cycles
* **Supplementary 5A / 5B / 5D** — in vivo clusters, in vitro transferred identities, srs22a
* **Supplementary 5C / 5E / 5F** — the three differential-expression volcanoes

Not regenerated, for the reason given in `docs/reproducibility.md`: **Figure 1D**,
which is a recorded constant.

## The data

`scripts/fetch_data.py` reads `data/MANIFEST.json`, tries each file's mirrors in
turn and verifies the SHA-256 of whatever it gets, falling through to the next
mirror on a mismatch. `MANIFEST.json` is deliberately kept out of LFS so it can
always be read, whatever state the rest of `data/` is in. Four mirrors are
listed: the public GitHub copy of this repository, pinned to the commit holding
the current data, and then the two institutional servers.

A file is re-fetched if it is missing, the wrong size, fails its checksum, or is
a git-lfs pointer stub — so a clone made without git-lfs repairs itself. Extra
mirrors can be supplied without editing anything, and take priority over the
published ones:

```bash
python scripts/fetch_data.py --base-url https://example.org/bradyzoite/
BZFIG_DATA_URLS="https://mirror-a/ https://mirror-b/" python scripts/fetch_data.py
python scripts/fetch_data.py --check      # verify what is already on disk
python scripts/make_figures.py --no-fetch # fail rather than download
```

The dataset ships in plain formats that need no special library — MatrixMarket
for the expression matrix, gzipped CSV for the cell and gene metadata and the
embeddings, JSON for the palettes — alongside an `.h5ad` for convenience. Either
route loads through `bzfig.data.load_dataset`.

The matrix is the full 8,322-gene ToxoDB-65 universe. The analysis object had
been subset to the 8,170 genes on the nuclear chromosomes before it was saved;
the other 152 sit on unplaced contigs and carry the apicoplast and mitochondrial
transcripts, and they are recovered from the two objects the analysis integrated
and appended after the 8,170. `var.in_analysis_object` says which is which, and
the export checks that the first 8,170 columns are still bit-identical to the
deposited layer. Only the volcano panels look at the recovered genes.

Figure 3E and 3F come from a second experiment — in vivo tachyzoites
integrated with the in vivo bradyzoites by scVI. The package carries that
integration in full: the count matrix the model was trained on, its cell and gene
tables, the trained checkpoint, and the UMAP coordinates the two panels are drawn
from. Rendering them needs nothing extra; `scripts/export_integration.py` regenerates
the coordinates from the checkpoint and needs the `integration` extra
(`pip install -e '.[integration]'`).

`scripts/export_dataset.py` documents how the deposited dataset was cut down from
the full analysis object, and writes `MANIFEST.json`. It is the only thing that
does: `export_integration.py` adds the Figure 3E/3F section through the same
code, so the checksums and the file list are never edited by hand.

## Layout

```
data/           deposited input (LFS) + MANIFEST.json
scripts/        fetch_data.py, make_figures.py, export_dataset.py, export_integration.py
src/bzfig/      panels.py, figure_2h_supplementary_4.py, figure_3ef.py,
                constants.py, data.py, de.py
tests/          unit tests for the Figure 2H / Supplementary 4 heatmaps
figures/        rendered output
docs/           reproducibility notes
preview/        self-contained HTML preview of every panel
```

`src/bzfig/constants.py` holds every value the figures depend on that is not
derivable from the data — palettes, colour limits, gene lists, and the counts
that were hard-coded when the figures were made. Each is annotated with where it
came from. `src/bzfig/figure_2h_supplementary_4_metadata.json` does the same for
the two heatmaps it is named for: their published gene order, row labels, colour
lookup and the manually highlighted rows, each with its source recorded.

Figure 2H and Supplementary 4 were contributed by Kourosh Zarringhalam.
