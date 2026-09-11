# Bradyzoite heterogeneity — figure reproduction

Everything needed to regenerate the figure panels derived from the single-cell
RNA-seq data: the deposited dataset, the code that draws each panel, and a record
of what reproduces exactly and what does not.

```bash
uv sync                              # or: pip install -e .
python scripts/fetch_data.py         # download + verify the input data
python scripts/make_figures.py       # render every panel to PNG, SVG and PDF
python scripts/make_preview.py       # build preview/index.html from what was rendered
```

Rendered panels land in `figures/`. `python scripts/make_figures.py --list` shows
the panel names; `--panel 1C 1G` renders a subset.

`preview/index.html` is a single self-contained page showing every rendered panel
with its metadata and its reproducibility status — no network access, so it can be
opened from disk or served as a static file. It needs the `preview` extra
(`pip install -e '.[preview]'`); `--max-edge` controls the size of the embedded
images.

## Start here

| | |
| --- | --- |
| [`docs/reproducibility.md`](docs/reproducibility.md) | What reproduces exactly, what is a reconstruction, what is a recorded constant — and why |
| [`docs/todo-for-authors.md`](docs/todo-for-authors.md) | Open questions and missing inputs |
| [`data/MANIFEST.json`](data/MANIFEST.json) | Every input file with its SHA-256 and download mirrors |

## Panels

Drawn from the deposited data:

* **Figure 1B** — cluster UMAP of the in vivo bradyzoites, and cells per cluster
* **Figure 1C** — selected marker genes per cluster
* **Figure 1E / 1F** — per-gene expression UMAPs, and cst1 (srs44)
* **Figure 1G** — CST1/SRS44 expression violins per cluster
* **Supplementary 1A** — all cluster markers per cluster
* **Supplementary 4** — cell-cycle regulators per phase, common (CCC) and modified (MCC) cell cycles
* **Supplementary 5A / 5B / 5D** — in vivo clusters, in vitro transferred identities, srs22a

Not regenerated, for reasons given in `docs/reproducibility.md`: **Figure 1D**
(recorded constant) and the **Supplementary 5C / 5E / 5F** volcano plots
(the underlying gene set is not in the deposited object).

## The data

`scripts/fetch_data.py` reads `data/MANIFEST.json`, tries each file's mirrors in
turn and verifies the SHA-256 of whatever it gets, falling through to the next
mirror on a mismatch. Extra mirrors can be supplied without editing anything:

```bash
python scripts/fetch_data.py --base-url https://example.org/bradyzoite/
BZFIG_DATA_URLS="https://mirror-a/ https://mirror-b/" python scripts/fetch_data.py
python scripts/fetch_data.py --check      # verify what is already on disk
```

The dataset ships in plain formats that need no special library — MatrixMarket
for the expression matrix, gzipped CSV for the cell and gene metadata and the
embeddings, JSON for the palettes — alongside an `.h5ad` for convenience. Either
route loads through `bzfig.data.load_dataset`.

`scripts/export_dataset.py` documents how the deposited dataset was cut down from
the full analysis object.

## Layout

```
data/           deposited input (LFS) + MANIFEST.json
scripts/        fetch_data.py, make_figures.py, export_dataset.py
src/bzfig/      panels.py, constants.py, scatter3d.py, data.py
figures/        rendered output
docs/           reproducibility notes, open questions
preview/        self-contained HTML preview of every panel
```

`src/bzfig/constants.py` holds every value the figures depend on that is not
derivable from the data — palettes, colour limits, gene lists, and the counts
that were hard-coded when the figures were made. Each is annotated with where it
came from.
