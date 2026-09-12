# Bradyzoite heterogeneity — figure reproduction

Data and code for the computational panels listed below, with a record of
numerical matches and outstanding reconstruction/assembly questions. This is
not yet a claim of coverage of every single-cell panel in the manuscript.

```bash
uv sync --extra preview               # or: pip install -e '.[preview]'
uv run python scripts/make_figures.py # render every panel to PNG, SVG and PDF
uv run python scripts/make_preview.py # build preview/index.html from what was rendered
```

**You do not need git-lfs.** The input data is stored in LFS, but if your clone
left pointer stubs in `data/` instead of the files, `make_figures.py` notices and
downloads them from the public mirrors before rendering. Run
`python scripts/fetch_data.py` to do that step on its own.

Rendered panels land in `figures/`. `python scripts/make_figures.py --list` shows
the panel names; `--panel 1C 1G` renders a subset.

If using pip instead of uv, activate that environment and use `python` directly.
To render only Kourosh's contribution, without rerunning any other panels:

```bash
python scripts/make_figures.py --panel Figure_2H Supplementary_4
python -m unittest discover -s tests -v
```

These heatmaps read the same shared dataset, recover the original R normalization
and write full-precision companion tables in `figures/tables/`. No Seurat object
or additional expression download is needed. See `docs/kourosh-contribution.md`.

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
* **Figure 2H** — original 5 × 29 Pearson correlations across the 6,505 cells
* **Supplementary 1A** — all cluster markers per cluster
* **Supplementary 4** — cell-cycle regulators per phase, common (CCC) and modified (MCC) cell cycles
* **Supplementary 5A / 5B / 5D** — in vivo clusters, in vitro transferred identities, srs22a
* **Supplementary 5C / 5E / 5F** — the three differential-expression volcanoes

**Figure 1D** is currently drawn from recorded counts. The manuscript workbook's
final lists independently support these counts; the historical selection recipe
is distinct from replotting them. See `docs/reproducibility.md`.

## The data

`scripts/fetch_data.py` reads `data/MANIFEST.json`, tries each file's mirrors in
turn and verifies the SHA-256 of whatever it gets, falling through to the next
mirror on a mismatch. `MANIFEST.json` is deliberately kept out of LFS so it can
always be read, whatever state the rest of `data/` is in.

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

The manifest tries a commit-pinned public GitHub LFS mirror first, then the
institutional mirrors, avoiding their timeout when the public copy is available.
Update the pin when publishing a new data release; hashes reject stale files.
TLS certificate verification remains enabled; a Python installation lacking its
CA bundle must have that local installation repaired, not disable verification.

The dataset ships in plain formats that need no special library — MatrixMarket
for the expression matrix, gzipped CSV for the cell and gene metadata and the
embeddings, JSON for the palettes — alongside an `.h5ad` for convenience. Either
route loads through `bzfig.data.load_dataset`.

A second, narrow matrix carries the 152 genes that were filtered out of the
deposited object before it was saved — the unplaced contigs, with the apicoplast
and mitochondrial transcripts on them. Only the volcano panels need it, through
`bzfig.data.load_extra_genes`; `bzfig.de` puts the two matrices side by side to
run its current reconstruction. The historical testing universe and small count
differences still need confirmation; see the reproducibility notes.

`scripts/export_dataset.py` documents how the deposited dataset was cut down from
the full analysis object.

## Layout

```
data/           deposited input (LFS) + MANIFEST.json
scripts/        fetch_data.py, make_figures.py, export_dataset.py
src/bzfig/      panels.py, constants.py, scatter3d.py, data.py, de.py
figures/        rendered output
docs/           reproducibility notes, open questions
preview/        self-contained HTML preview of every panel
```

`src/bzfig/constants.py` and `src/bzfig/kourosh_metadata.json` hold values that are not
derivable from the data — palettes, colour limits, gene lists, and the counts
that were hard-coded when the figures were made. Each is annotated with where it
came from.
