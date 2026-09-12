# Figure 2H and Supplementary Figure 4

These panels are rendered by `src/bzfig/figure_2h_supplementary_4.py` through
the repository's standard figure and preview workflow. They use the same
deposited `logcounts`, cell metadata and gene identifiers as the other panels;
no additional expression data, Seurat object or R installation is required.

## Run

```bash
uv sync --extra preview
uv run python scripts/make_figures.py --panel Figure_2H Supplementary_4
uv run python -m unittest discover -s tests -v
uv run python scripts/make_preview.py
```

For pip, install `pip install -e '.[preview]'` in an environment and use
`python` without `uv run`. The preview combines the new renders with the other
existing PNGs in `figures/`; it does not require rerunning other panels.

## Inputs and methods

The original R heatmaps used Seurat LogNormalize after selecting the 8,170
deposited genes. The shared `logcounts` was normalized before this subset.
The panel-specific normalization is recovered as:

```
L_R = log1p(10000 * expm1(L_shared) / rowSum(expm1(L_shared)))
```

The denominator includes **all 8,170 genes**, before selecting plotted genes,
not the expanded 8,322-gene differential-expression matrix. The conversion
does not modify any shared input or expression layer.

- **Figure 2H:** Pearson correlations across 6,505 in-vivo cells; five genes
  against all 29 displayed genes, in the original order.
- **Supplementary 4:** phase means for 50 genes, standardized using the sample
  SD (`ddof=1`) of the five means, independently within CCC and MCC. CCC phase
  cell counts are 11/230/106/56/44; MCC counts are 1517/2789/1219/271/262.

Original gene labels, independent CCC/MCC row orders and the LAB-interpolated
color palette are stored in `src/bzfig/figure_2h_supplementary_4_metadata.json`.
The pale-pink description backgrounds reproduce manual annotations in the
source supplementary figure: 19 CCC rows and 20 MCC rows, with AP2XI-4
highlighted only in MCC. These are explicit rendering metadata, not selections
inferred from expression. Source filename, page and SHA-256 are recorded.

## Outputs

Each panel is saved in PNG, SVG and PDF formats under `figures/`:

- `Figure_2H_correlation_heatmap`
- `Supplementary_4_CCC`
- `Supplementary_4_MCC`

Seven full-precision CSV tables are regenerated under `figures/tables/`:
the Figure 2H correlation matrix and, for each cycle, phase means, z-scores and
phase cell counts. Reference values in `tests/reference_heatmaps.json` are
test fixtures only and are never used as plotting inputs.

## Validation

All **145/145** Figure 2H correlations match the original displayed values at
two-decimal precision; the largest full-precision difference is below **2e-9**.
All **500/500** Supplementary 4 z-scores agree with the original R calculations
within **6e-8**, inside the automated tolerance of **1e-6**.

The tests cover numerical agreement, row order and phase counts, non-mutation
of shared data, rejection of the wrong gene universe, equivalent plain/H5AD
input routes, CSV outputs, data-free panel listing and the original highlight
identities, color and row placement.

Validated with Python 3.12, NumPy 2.2.5, SciPy 1.15.2, pandas 2.2.3,
anndata 0.12.7, Scanpy 1.11.5 and Matplotlib 3.10.3. Numerical reproduction
does not imply pixel-identical manual figure assembly or byte-identical
graphics across library versions. Broader panel coverage and reconstruction
limits are described in [reproducibility.md](reproducibility.md).
