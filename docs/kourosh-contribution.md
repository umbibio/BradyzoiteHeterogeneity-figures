# Kourosh's contribution: Figure 2H and Supplementary 4

Review-branch contribution based on Argenis's commit
`1c25d514bb4adf00b15f5b467ec68571aefd33c3`, 11 September 2026.
Numerical checks pass and Kourosh authorized sharing this contribution for
collaborator review. Final manuscript/assembly decisions remain with the authors
before merging or replacing manuscript figures. This is not a new upstream analysis.

## Scope and workflow

The existing runner now adds `Figure_2H_correlation_heatmap` and routes
`Supplementary_4_CCC` / `Supplementary_4_MCC` through the original R recipe,
implemented in Python. Other panel implementations, DE settings, embeddings,
cluster labels and shared expression layers are unchanged.

```bash
uv sync --extra preview
uv run python scripts/make_figures.py --panel Figure_2H Supplementary_4
uv run python -m unittest discover -s tests -v
uv run python scripts/make_preview.py
```

For pip, install `pip install -e '.[preview]'` in an environment and use `python`
without `uv run`. The preview needs the other existing PNGs present; there is
no reason to recompute their analyses solely to add these three panels.

Outputs are PNG/SVG/PDF plus seven full-precision CSV tables in `figures/tables`.
The tables contain computed correlations, CCC/MCC phase means and z-scores, and
the two phase cell-count vectors. Gene IDs, labels and row orders are in
`src/bzfig/kourosh_metadata.json`. These are rendering metadata, not additional
expression measurements. The original R reference values in
`tests/reference_heatmaps.json` are test fixtures only, never plot inputs.

## Why the normalization conversion is necessary

Argenis's `logcounts` was normalized before the 8,170-gene subset. The original
Kourosh R scripts (`heatmap_phase_mean.R`, `fig2H_heatmap.R`) used Seurat RNA/data
normalized after that subset. The original values can be recovered from the
same uploaded data, without counts, a separate RDS, or another data download:

```
L_R = log1p(10000 * expm1(L_shared) / rowSum(expm1(L_shared)))
```

The denominator includes **all 8,170 deposited genes**, before selecting plotted
genes; it must not include the expanded 8,322-gene DE universe. This conversion
is confined to the new module and does not overwrite any shared input.

Figure 2H uses Pearson correlations across the 6,505 original in-vivo cells,
five rows against 29 columns. Supplementary 4 uses 50 original genes, averages
within the five transferred phases, and uses the sample SD (`ddof=1`) of those
five means, independently for CCC and MCC. The original independent row orders,
manuscript aliases, and LAB-interpolated color palette are retained.

## QC result and limits

- All 11 public input files downloaded and passed the existing manifest hashes.
- All **145/145** Figure 2H values match at the displayed two-decimal precision;
  the largest full-precision correlation difference is below **2e-9**.
- All **500/500** Supplementary 4 z-scores match the audited R calculations
  within **6e-8**, well inside the test tolerance of **1e-6**.
- The eight tests check original numbers and row orders, phase cell counts,
  non-mutation of the shared layer, refusal of the wrong gene universe,
  equivalent plain/H5AD input routes, CSV production, data-free listing, and
  the source's exact per-panel highlight identities, color and row placement.
- This validates the original numerical recipe, not pixel-identical manual
  assembly. Pale-pink description backgrounds have been copied from the final
  supplementary PDF (page 4), with black text and unhighlighted gene IDs.
  The source's 19 CCC / 20 MCC asymmetry is retained: AP2XI-4 (TGME49_315760)
  is highlighted only in MCC. No expression-based selection rule was inferred.
  See the author checklist for this specific source discrepancy and the
  Figure 2H legend wording.

The tests were run with Python 3.12, NumPy 2.2.5, SciPy 1.15.2, pandas 2.2.3,
anndata 0.12.7, Scanpy 1.11.5 and Matplotlib 3.10.3. This records the tested
environment; it does not claim the repository already has a frozen dependency
lock or that byte-identical graphics are guaranteed across library versions.

## Bounded fixes accompanying the contribution

- `--list` works without downloading/loading input data; new panel calculations
  are lazy and do not run when an unrelated panel is selected.
- Added a commit-pinned public GitHub LFS mirror with the existing SHA-256
  verification. Download progress now flushes immediately. TLS verification is
  never disabled.
- Corrected the Supplementary 1A loader's misleading row-sort comment without
  changing the returned marker order.
- Corrected Supplementary 1/2 numbering and added Figure 2H to the preview.
- Updated Supplementary 4 provenance, normalization, scale and status; qualified
  unsupported exact-reproduction claims about the existing volcano reruns.
- Restored the original Supplementary 4 label backgrounds as explicit metadata;
  source filename, page and SHA-256 are recorded. No expression data changed.
- Corrected the preview dependency install command and recorded Figure 1D's
  workbook-list provenance. Neither panel's underlying calculation was changed.

Supplementary 5 historical settings/count reconciliation, final repository
licensing and archival identifiers remain author decisions. This contribution
does not retune statistics, reprocess raw reads, or replace manuscript figures.
