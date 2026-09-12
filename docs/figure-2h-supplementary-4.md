# Figure 2H and Supplementary Figure 4

These panels are rendered by `src/bzfig/figure_2h_supplementary_4.py` through the
repository's standard figure workflow. They use the same deposited `logcounts`,
cell metadata and gene identifiers as every other panel; no additional expression
data, Seurat object or R installation is required.

Contributed by **Kourosh Zarringhalam**. They replace the earlier Supplementary 4
reconstruction, which was fitted to the printed figure rather than derived from
the original calculation.

## Run

```bash
uv sync                                # or: pip install -e .
python scripts/make_figures.py --panel Figure_2H Supplementary_4
python -m unittest discover -s tests -v
```

Under `uv`, prefix the two Python commands with `uv run`. Rendering these three
panels does not require rerunning the other forty-two.

## Inputs and methods

The original R heatmaps ran Seurat `LogNormalize` *after* selecting the 8,170
deposited genes. The shared `logcounts` layer was normalised before that subset.
The panel-specific normalisation is recovered as:

```
L_R = log1p(1e4 * expm1(L_shared) / rowSum(expm1(L_shared)))
```

The denominator includes **all 8,170 genes**, before any plotted gene is
selected — not the expanded 8,322-gene differential-expression matrix, which the
module refuses outright. The conversion copies the matrix; it does not modify any
shared input or expression layer, and no other panel sees it.

- **Figure 2H:** Pearson correlations across the 6,505 in vivo cells; five genes
  against all 29 displayed genes, in the original order.
- **Supplementary 4:** phase means for 50 genes, standardised using the sample SD
  (`ddof=1`) of the five means, independently within CCC and MCC. CCC phase cell
  counts are 11/230/106/56/44; MCC counts are 1517/2789/1219/271/262.

Original gene labels, the independent CCC/MCC row orders and the LAB-interpolated
colour lookup are stored in `src/bzfig/figure_2h_supplementary_4_metadata.json`.
The pale-pink description backgrounds reproduce manual annotations in the source
supplementary figure: 19 CCC rows and 20 MCC rows, with AP2XI-4 highlighted in
MCC only. These are explicit rendering metadata, not selections inferred from
expression. Source filename, page and SHA-256 are recorded with them.

## Outputs

Each panel is saved as PNG, SVG and PDF under `figures/`:

- `Figure_2H_correlation_heatmap`
- `Supplementary_4_CCC`
- `Supplementary_4_MCC`

Seven full-precision CSV tables are regenerated under `figures/tables/`: the
Figure 2H correlation matrix and, for each cycle, the phase means, the z-scores
and the phase cell counts. Reference values in `tests/reference_heatmaps.json`
are test fixtures only and are never used as plotting inputs.

## Validation

Measured against the original values in `tests/reference_heatmaps.json`:

| | |
| --- | --- |
| Figure 2H, displayed two-decimal values matched | **145 / 145** |
| Figure 2H, largest full-precision difference | **1.87e-09** |
| Supplementary 4, z-scores within the 1e-6 test tolerance | **500 / 500** |
| Supplementary 4, largest z-score difference | **5.60e-08** |
| Supplementary 4, largest phase-mean difference | **1.07e-08** |

The residual is float32 rounding in the deposited layer. The denominator is what
decides it: the same correlations computed from `logcounts` as shipped, or
re-closed over the 8,322-gene universe, land 1.0e-03 away and match only 143 of
the 145 printed values.

The tests cover numerical agreement, row order and phase counts, non-mutation of
the shared data, rejection of the wrong gene universe, the equivalence of the
plain and H5AD input routes, the CSV outputs, data-free panel listing, and the
identity, colour and row placement of the manual highlights.

Two cross-checks independent of the reference file: the 50-gene list here is
set-identical to the one transcribed from the printed figure in
`bzfig.constants.SUPP4_GENES`, and the manuscript row orders recorded here are
exactly those the earlier reconstruction derived by sorting on peak phase then
peak height — 0 of 50 rows differ, in either panel.

Validated in this tree with Python 3.12, NumPy 2.1.3, SciPy 1.17.1, pandas 3.0.1,
anndata 0.12.6, Scanpy 1.12 and Matplotlib 3.10.8. Numerical reproduction does
not imply pixel-identical manual figure assembly, nor byte-identical graphics
across library versions; the `%.17g` CSVs in particular track the local BLAS in
their last two or three digits. Broader panel coverage and reconstruction limits
are in [reproducibility.md](reproducibility.md).
