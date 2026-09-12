# What reproduces, and what does not

This document distinguishes panels reproduced from the deposited data,
original-recipe numerical matches, reconstructed analyses and recorded
constants. Figure 2H and Supplementary 4 methods and validation are detailed in
[figure-2h-supplementary-4.md](figure-2h-supplementary-4.md).

## Reproduced exactly

| Panel | What it shows |
| --- | --- |
| Figure 1B | Cluster UMAP of the in vivo bradyzoites, plus cells per cluster |
| Figure 1C | Selected marker genes per cluster |
| Figure 1E | Six per-gene expression UMAPs (tachyzoite and bradyzoite markers) |
| Figure 1F | cst1 (srs44) expression |
| Figure 1G | CST1/SRS44 expression violins per cluster |
| Figure 3A | In vivo bradyzoites by cell-cycle phase — the UMAP; see the bar below |
| Figure 3B | me49 Day 0 by cell-cycle phase, and its phase bar |
| Figure 3D | The in vitro cells picked out of the projection |
| Figure 6 | The first two of the three UMAPs |
| Supplementary 1A | All cluster markers per cluster |
| Supplementary 1B | Five per-gene expression UMAPs (microneme transcripts) |
| Supplementary 1C | Six per-gene expression UMAPs (known cyst wall proteins) |
| Supplementary 5A | As Figure 1B, with a legend |
| Supplementary 5B | In vitro bradyzoites coloured by transferred cluster identity |
| Supplementary 5D | srs22a expression |

For Figure 1C, Figure 1G and Supplementary 1A the exported values were
additionally re-derived from the data object and re-plotted from the exported
tables alone; both matched the published panels.

Supplementary 1B and 1C are the same panel type as Figure 1E, drawn by the same
function over the same 6505 in vivo cells. Every one of the eleven genes was
checked against `var` and against the per-gene PNG the analysis left behind in
`notebooks/figures/`. Nine agree with their annotation; the two that cannot are
`cst4` (`TGME49_261650`) and `cst10` (`TGME49_312330`), both "hypothetical
protein" with no `gene_name` in the ToxoDB-65 annotation shipped here, so their
labels rest on the caption and on those filenames. Panel 1B also holds a
BioRender cartoon of the two microneme subpopulations, which is figure assembly,
not a panel.

## Reproduced from the original recipe — Figure 2H and Supplementary Figure 4

**Figure 2H** (5 × 29 Pearson correlation heatmap) and **Supplementary Figure 4**
(CCC and MCC cell-cycle regulator heatmaps), drawn by
`bzfig.figure_2h_supplementary_4`. This pair was contributed by **Kourosh
Zarringhalam**; the method notes are in
[figure-2h-supplementary-4.md](figure-2h-supplementary-4.md).

These two panels are not a reconstruction fitted to the printed image. They
reproduce the arithmetic of the original R code, from the deposited `logcounts`
alone — no Seurat object, no R runtime, no extra expression download.

**The normalisation is the whole trick.** The original heatmap scripts ran
Seurat `LogNormalize` *after* subsetting to the 8,170 deposited genes, whereas
the deposited `logcounts` layer was normalised before that subset, over the
wider pre-subset universe. Re-closing the deposited layer over the 8,170 genes
recovers the original values:

```
L_R = log1p(1e4 * expm1(L_shared) / rowSum(expm1(L_shared)))
```

The row sum runs over all 8,170 genes, before any plotted gene is selected. It
is local to these two panels and copies the matrix, so no shared layer moves.
The denominator matters and the 8,322-gene differential-expression universe is
the wrong one for it — measured on the Figure 2H matrix against the original
values:

| Denominator | max abs. difference | displayed 2-dp values matched |
| --- | --- | --- |
| none (deposited `logcounts` as shipped) | 1.0e-03 | 143 / 145 |
| the 8,322-gene DE universe | 1.0e-03 | 143 / 145 |
| **the 8,170 deposited genes** | **1.9e-09** | **145 / 145** |

The middle row is a no-op rather than a near miss. For these 6,505 in vivo cells
`expm1(logcounts)` already sums to 10,000 over the 8,322 genes (9999.998 …
10000.002 across the cohort), so re-closing over that set returns the deposited
layer to within 1.6e-07 — float32 rounding. Over the 8,170 deposited genes the
same totals run 7,625 … 10,000, and it is that missing mass, gene by gene and
cell by cell, that the R scripts divided out after subsetting. This is the same
finding as *Why the 8322-gene universe is the right one* below, arrived at from
the other end.

**What was measured here**, re-running the port against the original values in
`tests/reference_heatmaps.json`:

* Figure 2H — all **145 / 145** displayed correlations agree at two decimals;
  largest full-precision difference **1.87e-09**.
* Supplementary 4 — all **500 / 500** z-scores agree to within **5.60e-08**,
  inside the suite's 1e-6 tolerance; phase means to within 1.07e-08.
* Phase cell counts are exactly the panels': CCC 11 / 230 / 106 / 56 / 44,
  MCC 1517 / 2789 / 1219 / 271 / 262, over the 6,505 in vivo cells.

The residuals are float32 rounding in the deposited layer, not a difference of
method. Bit-identical equality is not claimed and could not be.

**What this settles.** The earlier version of these notes carried Supplementary 4
as a reconstruction matched to the printed figure (mean per-gene r = 0.9998 CCC,
0.9997 MCC). Two independent routes now agree on its inputs, which is worth more
than either alone:

* the 50-gene list transcribed from the printed figure
  (`bzfig.constants.SUPP4_GENES`) is **set-identical** to the list carried in the
  source data (`figure_2h_supplementary_4_metadata.json`);
* the independent CCC and MCC row orders transcribed from the manuscript are
  **exactly** the orders the old sort rule derived (peak phase, then peak
  height) — 0 of 50 rows differ, in either panel.

Row labels are now the manuscript's own — mostly short forms of the ToxoDB-65
descriptions (`AP2 IX5` for "AP2 domain transcription factor AP2IX-5"), 22 of the
50 differing in wording. One differs in substance and is still open:
`TGME49_291590` is "MYND-like zinc finger protein" in the figure and
"hypothetical protein" in the annotation shipped here.

The pale-pink description backgrounds are transcribed from page 4 of the
supplementary PDF rather than inferred from expression, which exposes something a
rule could not have: the source highlights **19** rows in CCC and **20** in MCC —
`TGME49_315760` (AP2XI-4) is highlighted in MCC only. Of the 20 highlighted genes,
17 peak in the same phase in both panels and no unhighlighted gene does, so the
caption's "highest within the same cell cycle state" accounts for the set apart
from three genes (`TGME49_207900`, `TGME49_262730`, `TGME49_318470`) that peak one
phase apart. The asymmetry is reproduced as found.

Two things worth knowing when reading the published panel: its colour bar axis
runs −4…4 but no value in either panel exceeds **±1.79** (the render here uses
the original ±2.5 breaks), and the CCC group is thinly populated, so several
genes are detected in one phase only and their rows sit at the one-hot extremes
of ±1.789 / −0.447.

**Companion tables.** Every run rewrites seven full-precision CSVs under
`figures/tables/` — the Figure 2H correlation matrix, and per cycle the phase
means, the z-scores and the phase cell counts. They are outputs, never inputs;
`tests/reference_heatmaps.json` holds the original values and is read only by the
test suite. Because the CSVs are written at `%.17g`, their last two or three
digits track the local BLAS: regenerated here they differ from the committed
copies by at most 4.4e-14 absolute (1.7e-13 relative), with row order and cell
counts identical.

The remaining gap is provenance, not arithmetic: the code that drew these panels
is still not in the analysis repository, and the CCC cell selection was made
interactively in a dashboard rather than in a script — it survives only as a
`pseudotime_UCC` column, carried into this dataset as `obs["cell_cycle_group"]`.

## Figure 3 and Figure 6 — the cohort panels

Figure 3A–3D and Figure 6's three UMAPs all draw the same picture: one cohort of
`obs["orig_ident"]` picked out of the integrated projection, coloured either by
`obs["transferred_cc_phase"]` or, in 3D, flat.

**The embedding is the 3-D one.** These are not 2-D scatters. They are
`bzfig.scatter3d.plot_3d_preview` over `obsm["3d_umap_harmony_integration"]` at
`elev=60, azim=0` — the same rendering path as Figure 1B and Supplementary 5A.
`obsm["2d_projection"]` holds that same view projected by hand, and it is a
tempting substitute, but matplotlib normalises each axis to the unit cube before
projecting: a 3-D render is stretched by the range of the points it is given,
and a scatter of the stored projection is not. Against the published Figure 3D the
3-D render matches at an occupancy IoU of **0.96** and a scatter of
`2d_projection` at **0.69**; the same test on Figure 3B gives 0.81 against 0.69.
The 3-D render is also what the analysis notebook did
(`integrated_adata_me49_nr_subset-2026-06-29.ipynb`, cells 131–137), and its
saved output `notebooks/figures/2026-06-29/umap_dataset.png` is the published
Figure 3D, matching it at IoU 0.98.

**Figure 3A has no grey layer.** 3B and 3C are drawn over the whole projection
with the other cells in `#d3d3d3`; 3A is the cohort alone, which makes it the
same image as Figure 6's first panel — the two published panels are the same
rendering, and the two files here are byte-identical. Figure 6 drops the grey
layer from all three.

**Figure 6's panels are named for their cohort, not for the column header above
them.** The headers and the cells beneath them disagree; that is item 8 of
`docs/todo-for-authors.md`, and the names used here
(`Figure_6_nonreactivated`, `Figure_6_me49_day0`, `Figure_6_me49_day3`) follow
the cells.

### The phase labels — two gaps, and the panel that explains them

Nothing below is tuned to. The phase bars here are drawn from the deposited
labels and the differences are recorded rather than closed.

The published bars were measured off the printed figure at 300 dpi, scaling each
bar's pixel length so the five bars sum to the cohort size.

| Panel | Published bar implies | Deposited `transferred_cc_phase` |
| --- | --- | --- |
| 3A (6505 cells) | 1516 / 3022 / 1350 / 318 / 299 | 1528 / 3019 / 1325 / 327 / 306 |
| 3B (602 cells) | 162 / 161 / 190 / 64 / 25 | **162 / 161 / 190 / 64 / 25** |
| 3C (950 cells) | 186 / 494 / 191 / 50 / 29 | 133 / 620 / 130 / 40 / 27 |

1. **Figure 3C's phase vector is not in either repository.** The published bar
   is about 50 cells off in G1a and 130 in G1b, and nothing deposited matches —
   not `transferred_cc_phase`, not `cc_phase` (which is `G1a 50 / G1b 230 /
   S 66 / M 63 / C 26` with 515 cells unassigned for this cohort), and no object
   in the analysis repository either.
   `notebooks/figures/me49d3_cc_phase_counts.png` reproduces the published
   ratios, so the vector existed; it was never saved. Both 3C's bar **and 3C's
   UMAP colours** are drawn from the deposited labels, so the panel here puts
   orange (G1b) where the published one has a blue (G1a) group in the upper
   cluster. The cells and their positions are unaffected — only which phase each
   is called. Figure 6's third panel inherits the same gap.
2. **Figure 3A's bar is about 1% off**, concentrated in the S bar (about 25
   cells). Almost certainly a slightly different vintage of the label transfer;
   `notebooks/` holds four dated `nr_transferred_cc_phase*.csv` variants and none
   of them match either. A handful of points change colour; the panel is
   otherwise the published one.
3. **Figure 3B matches the deposited labels exactly**, which is what makes the
   other two informative: the recipe is right and the label vector is not.

### Figure 3E and 3F are not from this dataset

They are not reproducible here, and not because something is missing from the
data package: **they come from a different experiment.** Both were drawn by
`notebooks/integrate.S1-S2-S3-2026-06-29.ipynb`, which scVI-integrates the 6505
non-reactivated cells with two further samples, S1 (166 cells) and S2 (210), and
plots the result on `obsm["X_scVI_2Dumap"]`:

* **3E** is that notebook's `NR-S1-S2-highlight` — the 376 S1/S2 cells picked out
  in red against the non-reactivated cells in grey. It matches the published
  panel at an occupancy IoU of **0.96** (0.87 on the red points alone).
* **3F** is its `NR-S1-S2-cc_phaseS1S2` — the same 376 cells coloured by
  `cc_phase`, everything else grey. IoU **0.95** (0.79 on the coloured points).

Neither panel involves `obs["cell_cycle_group"]`, `2d_projection`, or the
8057-cell object this repository ships. Reading 3E as the 447 CCC cells is a
near miss — the CCC selection does trace a similar arc through the in vivo
projection — but the embeddings are different point clouds: the best view of the
3-D UMAP reaches only IoU 0.53 against the published 3E, where a true match
scores 0.84. The S1 and S2 count matrices are not in either repository, so
nothing here can draw these two panels. Figure 3G (flow cytometry) and Figure 3H
are outside this dataset as well.

## Reproduced from a wider gene set — the volcano plots

**Supplementary Figure 5C, 5E and 5F.**

The existing reconstruction uses the full **8322-gene** ToxoDB-65 universe;
the deposited object has **8170** genes. The results below are a reconstruction,
not an exact match to the manuscript counts. Agreement with the historical
gene-level results is not established. The extra 152 sit on unplaced `KE*`
contigs, and among them are the apicoplast and mitochondrial transcripts — ORF F,
two cytochrome b's, cytochrome c oxidase III — that the published panels label at
the positive extreme. They now ship beside the deposited matrix as
`logcounts_extra.mtx.gz` (8057 cells x 152 genes, 86,934 stored values, 391 KB)
and `var_extra.csv.gz`; `bzfig.de` puts the two matrices side by side and runs
the test.

The recipe:

```
universe     all 8322 ToxoDB-65 genes (logcounts.mtx.gz + logcounts_extra.mtx.gz)
test         sc.tl.rank_genes_groups(groups=[A], reference=B, method="wilcoxon",
                                     tie_correct=False)
filter       genes whose mean expm1(logcounts) is > 0 in *both* groups
significant  pvals_adj < 0.05                      (bzfig.constants.DE_ALPHA)
cutoff       |logfoldchanges| > 2                  (bzfig.constants.DE_LOG2FC_CUTOFF)
```

| | Published | Reproduced |
| --- | --- | --- |
| 5C up / down in vivo | 664 / 1443 | **666 / 1441** |
| 5C log2FC extremes | −10.40 … +13.79 (measured) | −10.99 … +13.80 |
| 5E up / down | not quoted | 55 / 69 |
| 5E log2FC extremes | −5.68 … +9.77 (measured) | −7.96 … +9.77 |
| 5F up / down in vivo | 676 / 1146 | **678 / 1145** |
| 5F log2FC extremes | −12.20 … +13.13 (measured) | −12.19 … +13.07 |

The published extremes were measured off the printed figure at 300 dpi against
the axis ticks (60.35 px per log2 unit in 5C, 50.1 in 5F). Individual called-out
points land within 0.06 of the reproduced fold change — sag1 at −10.40 against
−10.42 in 5C, srs22a at +8.19 against +8.19, ORF F at +8.30 against +8.30,
cytochrome c oxidase III at +13.79 against +13.80. The following are the original
reconstruction's explanations for differences in the measured extremes; they
are not independent proof of historical settings:

* **5C**, one gene: a hypothetical protein at −10.99 whose adjusted p-value
  underflows to zero, so it is drawn as a triangle at the axis top — and the
  point detector keeps only components with circularity > 0.65, which excludes
  triangles by construction.
* **5E**, thirteen genes between −7.96 and −5.71 have adjusted p-values between
  0.7 and 1.0. They are in the full test table but are excluded by the plot's
  significance filter, not hidden on its baseline. The full-table and plotted
  fold-change ranges must not be conflated.
* **5F**, none at all — and its negative extreme, the one case with no low-y or
  off-axis points involved, matches to the last digit printed (−12.19 against a
  measured −12.20).

Raster comparisons cannot establish historical settings or gene-level agreement.
The measured plot extremes are approximate image-derived values, not source
data or criteria for adjusting the reconstructed analysis.

**The subsets.** `orig_ident` carries all three cohorts (Nonreactivated 6505,
me49 Day 3 950, me49 Day 0 602).

| Panel | Group A (numerator) | Group B (reference) |
| --- | --- | --- |
| 5C | in vivo bradyzoites, `Nonreactivated` (6505) | in vitro bradyzoites, `me49 Day 3` (950) |
| 5E | in vitro bradyzoite G1, `me49 Day 3` (280) | in vitro tachyzoite G1, `me49 Day 0` (272) |
| 5F | in vivo bradyzoite G1, `Nonreactivated` (1976) | in vitro bradyzoite G1, `me49 Day 3` (280) |

"G1" is `cc_phase` in G1a or G1b, not `transferred_cc_phase`: the transferred
column gives 4547 against 753 cells for 5F and **673 / 1338** genes, 192 off the
published 676 / 1146. This supports the reconstructed subset but does not
independently establish the historical selection.

### Evidence for the upstream normalization denominator

Normalisation ran *before* the gene subset, and per source object: each cell
divided by its own total over that object's full gene set, scaled to 1e4 and
`log1p`'d. Applied to the 8170 genes the deposited object kept, that reproduces
its `logcounts` layer to within one float32 ulp. Nothing else does.

| Denominator | Largest difference from the deposited `logcounts` |
| --- | --- |
| each cell's total in its own source object — 8322 genes for the in vivo 10x run, 8496 for the me49 object | **4.8e-07**, one float32 ulp |
| one 8322-gene total for every cell | 5.3e-03 — wrong for the me49 cells |
| the 8170 deposited genes | 0.27 — wrong for every cell |

This establishes the normalization denominator used upstream of the saved
subset; it does not independently establish the later differential-expression
testing universe.
`scripts/export_dataset.py` measures all three on every run, records them in
`MANIFEST.json` under `extra_genes`, and refuses to write the package if the
first row stops holding.

### Evidence supporting the reconstructed cutoff and filter

All three published panels have a clean empty band around zero — a zero-ink
rectangle spanning the full height, not a thinning. Measured at 600 dpi its edges
are −1.99 / +2.01 (5C), −2.20 / +2.23 (5E) and −2.02 / +2.03 (5F): the cutoff was
applied to the plotted fold change, and it is 2.0. A 1.5-fold cutoff (log2 =
0.585) is excluded by that measurement and by the counts — it gives 2326 up /
2885 down for 5C against a published 664 / 1443. It could still have been an
upstream pre-filter: a gene that fails |log2FC| > 0.585 also fails |log2FC| > 2,
so a 1.5-fold pre-filter would change none of these numbers.

The reconstruction uses a both-groups filter. Scanpy scores a gene with a
zero group mean against a 1e-9 pseudocount, which drops it into a degenerate
bucket: without the filter 5C's `logfoldchanges` run −28.19 … +18.24 and the down
count goes from 1441 to 1879, while the up count does not move at all. Dropping
exactly the genes undetected in one group brings the axis range closer to the
printed panel. This supports the reconstruction but does not uniquely identify
the original filtering recipe.

### What is *not* reproduced

* **The point colouring is only partly rule-derivable.** The published legends
  say pink = hypothetical protein and grey = ribosomal protein (blue = cyst wall
  protein in 5E), but the pink is applied to the *called-out* genes only: 1024 of
  5C's 2107 plotted genes are hypothetical proteins and the published panel has
  about 25 pink points, several of them named genes (sag1, srs2, srs22a, the
  cytochromes). The panels here colour the called-out set — the genes beyond the
  caption's annotation thresholds — pink, except ribosomal proteins, which are
  grey; that rule puts grey on exactly the five points the published 5F has grey.
  5E's blue set is not derivable at all (two of its genes are a dense granule and
  a rhoptry protein), so the seven genes were read back off the published panel
  by matching points to fold changes, and are recorded in
  `bzfig.constants.SUPP5E_CYST_WALL_GENES`.
* **The labels.** The caption's rule — annotate above log2FC 8 or below −7 for
  5C, 6 and −4 for 5E — is what `bzfig.constants.VOLCANO_LABEL_RANGE` uses (5F's
  caption gives no rule; its labelled genes fit 5C's). The published panels also
  label a hand-picked set of transcription factors (bfd1, the ap2s) and the two
  enolases, which no threshold recovers, and shorten long descriptions by hand.
* **Points off the top of the axis.** The published y limits are 1e-100, 1e-60
  and 1e-200, which is what these panels use. 100 of 5C's 2107 points have an
  adjusted p-value below 1e-100 — 20 of them underflow float64 to zero — and are
  not plotted, which is what the published panel does: it prints no marker for
  any of them. They remain in the DE table, so nothing is lost from the numbers;
  only the drawing clips. 5E and 5F have no such points.

  One consequence worth knowing when comparing panels side by side: srs2/p35 is
  among 5C's clipped genes here, while the published panel places it at about
  1e-88 and labels it. The reproduction puts its adjusted p-value slightly lower
  than the original run did; the fold change agrees.
* **ORF F in 5F.** `TGME49_302005` reaches log2FC +8.41 there but only
  pvals_adj = 0.23, so it is not among the plotted genes. The published 5F does
  not label it either; it labels ORF F in 5C, where the reproduction gives
  pvals_adj = 1.7e-04 and puts the point at +8.30, which is where the published
  panel has it. If the published 5F does plot a point for ORF F, that one gene is
  unexplained — an unlabelled point cannot be identified from a printed figure.
* **5E has no published counts** to check against, and its negative extreme
  cannot be checked either, for the reason given above. It is validated instead
  by its called-out genes, which land where the published panel has them: sag1 at −5.67 (published −5.68, both at 1e-60), sag4/srs35a at +9.46,
  bag1 at +9.77, the helicase at −4.01, the CMGC kinase at +6.23, and every one
  of the blue cyst wall points within 0.02 of a published point.

`scripts/export_dataset.py` re-runs all three comparisons from the data package
it has just written and fails if the counts move. They are recorded in the
manifest under `verification.supplementary_5_de` and in
`bzfig.constants.SUPP5_DE_COUNTS_REPRODUCED`; the published counts stay in
`SUPP5C_DE_COUNTS` and `SUPP5F_DE_COUNTS`.

## Not reproducible — recorded as constants

**Figure 1D**, "# Unique markers/Cluster" (38 / 93 / 80 / 881 / 14 / 47).

Hard-coded in the analysis notebook. The underlying quantity is evidently "genes
that are a significant Wilcoxon marker of exactly one cluster", but the original
run used a larger gene universe (8322 genes) than the deposited object (8170).
Re-running that test on the deposited object gives `[37, 94, 84, 861, 18, 52]`:

| Cluster | 0 | 1 | 2 | 3 | 4 | 5 |
| --- | --- | --- | --- | --- | --- | --- |
| Published | 38 | 93 | 80 | 881 | 14 | 47 |
| Re-derived here | 37 | 94 | 84 | 861 | 18 | 52 |

Re-running the same test on the full 8322-gene matrix gives
`[45, 93, 84, 863, 18, 53]`, also distinct from the manuscript counts. The final
source-data workbook's F1D gene lists contain exactly 38 / 93 / 80 / 881 / 14 / 47
genes. Replotting these counts is distinct from reproducing the historical
marker-selection calculation. This repository draws the panel from
`bzfig.constants.UNIQUE_MARKERS_PER_CLUSTER`; the workbook is not a runtime
input and the historical selection calculation is not reproduced here.

## Not included

**Supplementary Figure 2A–F** (cell-cycle-dependent expression of marker genes)
is not derived from this dataset at all — the caption attributes it to Benke et
al., mined from ToxoDB.

An earlier version of this page said that of *Supplementary Figure 1B–F*. That
was wrong twice over: Supplementary Figure 1 has only panels A, B and C, and 1B
and 1C are eleven per-gene expression UMAPs drawn from this dataset — the same
panel type as Figure 1E. They are rendered here.

**Figure 3E and 3F** are drawn from a different object, and **Figure 3G and 3H**
from flow cytometry and imaging; see the Figure 3 section above. The cartoons
above Figure 6's UMAPs and beside Supplementary 1B are BioRender artwork.

## The enolase labels — the figure is right

The repository uses the following gene labels, consistent with the manuscript
and the deposited gene metadata:

* `TGME49_268850` = **enolase 2**, the tachyzoite isoform — detected in 0.6% of
  in vivo cells, the near-blank panel in the tachyzoite-marker row
* `TGME49_268860` = **enolase 1**, the bradyzoite isoform — detected in 68.9%,
  the broadly stained panel in the bradyzoite-marker row

This is confirmed from two independent directions. The ToxoDB annotations used in
this study, releases 65 and 68, in both GFF and GTF, all carry
`TGME49_268850 = enolase 2` and `TGME49_268860 = enolase 1`, as do the `var`
table of the deposited object and NCBI Gene. And the primary literature gives the
same assignment: Dzierszinski et al. (1999) originally characterised ENO2 as the
tachyzoite isoform and ENO1 as bradyzoite-specific; Ferguson et al. (2002),
*Int J Parasitol* 32:1399–1410, localised ENO1 to brain tissue cyst bradyzoites
and ENO2 to replicating stages by immuno-EM in vivo; and Ngô et al. (2015),
*Acta Cryst D*, state the gene IDs explicitly — TgENO1 = TGME49_268860,
TgENO2 = TGME49_268850.

One point for readers coming from mammalian work: the human ENO1/ENO2 numbering
means something entirely different (ENO1 is the ubiquitous isoform there). The
*Toxoplasma* literature uses one convention throughout; the apparent conflict is
only with the human gene names.

## Panel assembly

The published figures were assembled in a vector editor. Panel letters, the
Figure 1E grid and its row labels, the single shared colour bar across 1E/1F, the
in-plot cluster numbers in Figure 1B and the coloured row-group bands beside
Figure 1C are all added at that stage. This repository renders the individual
panels as the analysis produced them; it does not attempt the assembly.
