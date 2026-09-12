# Reproducibility

This page records, panel by panel, how each figure is drawn from the deposited
data, and where a panel departs from the published one.

## Panels drawn directly from the deposited data

| Panel | What it shows |
| --- | --- |
| Figure 1B | Cluster UMAP of the in vivo bradyzoites, plus cells per cluster |
| Figure 1C | Selected marker genes per cluster |
| Figure 1E | Six per-gene expression UMAPs (tachyzoite and bradyzoite markers) |
| Figure 1F | cst1 (srs44) expression |
| Figure 1G | CST1/SRS44 expression violins per cluster |
| Figure 3A | In vivo bradyzoites by cell-cycle phase |
| Figure 3B | me49 Day 0 by cell-cycle phase, and its phase bar |
| Figure 3D | The in vitro cells picked out of the projection |
| Figure 6 | The three cohort UMAPs |
| Supplementary 1A | All cluster markers per cluster |
| Supplementary 1B | Five per-gene expression UMAPs (microneme transcripts) |
| Supplementary 1C | Six per-gene expression UMAPs (known cyst wall proteins) |
| Supplementary 5A | As Figure 1B, with a legend |
| Supplementary 5B | In vitro bradyzoites coloured by transferred cluster identity |
| Supplementary 5D | srs22a expression |

For Figure 1C, Figure 1G and Supplementary 1A the plotted values were also
re-derived from the data object and re-plotted from the exported tables alone;
both routes match the published panels.

Supplementary 1B and 1C are the same panel type as Figure 1E, drawn by the same
function over the same 6505 in vivo cells. Two of the eleven genes, `cst4`
(`TGME49_261650`) and `cst10` (`TGME49_312330`), are annotated "hypothetical
protein" with no `gene_name` in the ToxoDB-65 annotation shipped here, so their
labels follow the published caption. Panel 1B also carries a BioRender cartoon
of the two microneme subpopulations, which is figure assembly rather than a
rendered panel.

## Figure 2H and Supplementary Figure 4

**Figure 2H** (5 × 29 Pearson correlation heatmap) and **Supplementary Figure 4**
(CCC and MCC cell-cycle regulator heatmaps) are drawn by
`bzfig.figure_2h_supplementary_4`; the method notes are in
[figure-2h-supplementary-4.md](figure-2h-supplementary-4.md). Both reproduce the
arithmetic of the original R code from the deposited `logcounts` alone.

The normalisation is the key step. The original heatmap scripts ran Seurat
`LogNormalize` *after* subsetting to the 8,170 deposited genes, whereas the
deposited `logcounts` layer was normalised before that subset, over the wider
pre-subset universe. Re-closing the deposited layer over the 8,170 genes recovers
the original values:

```
L_R = log1p(1e4 * expm1(L_shared) / rowSum(expm1(L_shared)))
```

The row sum runs over all 8,170 genes, before any plotted gene is selected. It is
local to these two panels and copies the matrix, so no shared layer is modified.
The denominator matters — measured on the Figure 2H matrix against the original
values:

| Denominator | max abs. difference | displayed 2-dp values matched |
| --- | --- | --- |
| none (deposited `logcounts` as shipped) | 1.0e-03 | 143 / 145 |
| the 8,322-gene differential-expression universe | 1.0e-03 | 143 / 145 |
| **the 8,170 deposited genes** | **1.9e-09** | **145 / 145** |

The middle row is a no-op rather than a near miss: for these 6,505 in vivo cells
`expm1(logcounts)` already sums to 10,000 over the 8,322 genes, so re-closing
over that set returns the deposited layer unchanged to within float32 rounding.
Over the 8,170 deposited genes the same totals run 7,625 … 10,000, and it is that
missing mass, gene by gene and cell by cell, that the original scripts divided
out after subsetting.

Measured against the original values in `tests/reference_heatmaps.json`:

* Figure 2H — all **145 / 145** displayed correlations agree at two decimals;
  largest full-precision difference **1.87e-09**.
* Supplementary 4 — all **500 / 500** z-scores agree to within **5.60e-08**,
  inside the test suite's 1e-6 tolerance; phase means to within 1.07e-08.
* Phase cell counts match the panels exactly: CCC 11 / 230 / 106 / 56 / 44,
  MCC 1517 / 2789 / 1219 / 271 / 262, over the 6,505 in vivo cells.

The residuals are float32 rounding in the deposited layer, not a difference of
method; bit-identical equality is not claimed.

Row labels are the manuscript's own, mostly short forms of the ToxoDB-65
descriptions (`AP2 IX5` for "AP2 domain transcription factor AP2IX-5"). One
differs in substance: `TGME49_291590` is "MYND-like zinc finger protein" in the
figure and "hypothetical protein" in the annotation shipped here.

The pale-pink description backgrounds are transcribed from page 4 of the
supplementary PDF rather than derived from expression. The source highlights
**19** rows in CCC and **20** in MCC — `TGME49_315760` (AP2XI-4) is highlighted in
MCC only. Of the 20 highlighted genes, 17 peak in the same phase in both panels
and no unhighlighted gene does, so the caption's "highest within the same cell
cycle state" accounts for the set apart from three genes (`TGME49_207900`,
`TGME49_262730`, `TGME49_318470`) that peak one phase apart. The asymmetry is
reproduced as found.

Two points for reading the published panel: its colour bar axis runs −4…4 but no
value in either panel exceeds **±1.79** (the render here uses the original ±2.5
breaks), and the CCC group is thinly populated, so several genes are detected in
one phase only and their rows sit at the one-hot extremes of ±1.789 / −0.447.

Every run rewrites seven full-precision CSVs under `figures/tables/` — the
Figure 2H correlation matrix, and per cycle the phase means, the z-scores and the
phase cell counts. They are outputs, never inputs; `tests/reference_heatmaps.json`
holds the original values and is read only by the test suite. Because the CSVs are
written at `%.17g`, their last two or three digits track the local BLAS:
regenerated copies differ by at most 4.4e-14 absolute (1.7e-13 relative), with row
order and cell counts identical.

The CCC cell selection was made interactively in a dashboard rather than in a
script, and survives only as a `pseudotime_UCC` column, carried into this dataset
as `obs["cell_cycle_group"]`.

## Figure 3 and Figure 6 — the cohort panels

Figure 3A–3D and Figure 6's three UMAPs all draw the same picture: one cohort of
`obs["orig_ident"]` picked out of the integrated projection, coloured either by
`obs["transferred_cc_phase"]` or, in 3D, flat.

These are 3-D renders, not 2-D scatters:
`bzfig.scatter3d.plot_3d_preview` over `obsm["3d_umap_harmony_integration"]` at
`elev=60, azim=0`, the same rendering path as Figure 1B and Supplementary 5A.
`obsm["2d_projection"]` holds that same view projected by hand, but matplotlib
normalises each axis to the unit cube before projecting, so a 3-D render is
stretched by the range of the points it is given and a scatter of the stored
projection is not. Against the published Figure 3D the 3-D render matches at an
occupancy IoU of **0.96** and a scatter of `2d_projection` at **0.69**; on
Figure 3B, 0.81 against 0.69.

Figure 3A has no grey background layer. 3B and 3C are drawn over the whole
projection with the other cells in `#d3d3d3`; 3A is the cohort alone, which makes
it the same image as Figure 6's first panel, and the two files here are
byte-identical. Figure 6 drops the grey layer from all three.

Figure 6's panels are named for their cohort rather than for the column header
printed above them, which the cells beneath do not match:
`Figure_6_nonreactivated`, `Figure_6_me49_day0`, `Figure_6_me49_day3`.

### The phase labels

The phase bars are drawn from the deposited labels. The published bars were
measured off the printed figure at 300 dpi, scaling each bar's pixel length so the
five bars sum to the cohort size.

| Panel | Published bar implies | Deposited `transferred_cc_phase` |
| --- | --- | --- |
| 3A (6505 cells) | 1516 / 3022 / 1350 / 318 / 299 | 1528 / 3019 / 1325 / 327 / 306 |
| 3B (602 cells) | 162 / 161 / 190 / 64 / 25 | **162 / 161 / 190 / 64 / 25** |
| 3C (950 cells) | 186 / 494 / 191 / 50 / 29 | 133 / 620 / 130 / 40 / 27 |

Figure 3B matches the deposited labels exactly. The other two do not:

1. **Figure 3C's phase vector is not deposited.** Neither
   `transferred_cc_phase` nor `cc_phase` (`G1a 50 / G1b 230 / S 66 / M 63 / C 26`,
   with 515 cells unassigned for this cohort) matches the published bar. Both 3C's
   bar and 3C's UMAP colours are therefore drawn from the deposited labels, which
   puts orange (G1b) where the published panel has a blue (G1a) group in the upper
   cluster. The cells and their positions are unaffected — only which phase each
   is called. Figure 6's third panel inherits the same difference.
2. **Figure 3A's bar is about 1% off**, concentrated in the S bar (about 25
   cells), consistent with a slightly different vintage of the label transfer. A
   handful of points change colour; the panel is otherwise the published one.

### Figure 3E and 3F

These two panels come from a second experiment. A newly generated dataset of in
vivo tachyzoites (5 dpi, peritoneal cavity) was integrated with the 6,505 in vivo
bradyzoites by scVI, and the 10-dimensional latent space laid out as a UMAP: 3E
picks the tachyzoites out in red, 3F colours them by `cc_phase`.

The integration ships with this package as four files — the 6,881 × 8,778 count
matrix the model was trained on, its cell and gene tables, and the trained
checkpoint — plus `figure_3ef_embedding.csv.gz`, the UMAP coordinates the panels
are drawn from. Rendering needs neither torch nor scvi-tools;
`scripts/integrate_s1_s2.py` regenerates the embedding from the checkpoint and
needs the `integration` extra.

The published integration uses samples **S1 (166 cells) and S2 (210)**, 376
tachyzoites in all. A third sample, S3, was collected and annotated (55 cells)
but is commented out of the analysis notebook's sample list and is not part of
the published panels; it is not shipped here.

**The layout is the checkpoint's, not the published figure's.** Loading the
deposited checkpoint is deterministic — the latent is bit-identical across calls,
and so is the UMAP computed from it — so this embedding reproduces exactly, and
`scripts/integrate_s1_s2.py` reports the difference against the shipped file
whenever it runs. What does not carry over is the frame: training a fresh model
with the same settings converges to the same quality (89 epochs at best
validation ELBO 1403.7, against the checkpoint's 94 at 1393.4) but lays the cells
out differently, and two such runs differ from each other about as much as either
differs from the published panel. The published figure was drawn from one
particular fit whose weights are not the ones deposited here. The cloud, the
loop, and the tachyzoites' position on it are the published panel's; their
orientation on the page is not.

What the panels show does not depend on the layout. Taking each tachyzoite's 30
nearest bradyzoite neighbours in the latent space itself:

| | CCC share of neighbours |
| --- | --- |
| all in vivo bradyzoites (baseline) | 6.9% |
| in vivo bradyzoites' own neighbourhoods | 7.0% |
| **the 376 in vivo tachyzoites** | **40.3%** |

A 5.8-fold enrichment, with 41.5% of tachyzoites in a CCC-majority neighbourhood
(S1 36.2%, S2 43.6%) — the overlap with the common cell cycle that Figure 3E
reports, measured rather than read off the projection.

Figure 3G (flow cytometry) and Figure 3H are outside this dataset.

## The volcano panels — Supplementary Figure 5C, 5E and 5F

The published differential expression ran on the full **8322-gene** ToxoDB-65
universe; the deposited object carries **8170** genes. The missing 152 all sit on
unplaced `KE*` contigs and include the apicoplast and mitochondrial transcripts —
ORF F, two cytochrome b's, cytochrome c oxidase III — that the published panels
label at the positive extreme. They ship beside the deposited matrix as
`logcounts_extra.mtx.gz` (8057 cells × 152 genes, 86,934 stored values, 391 KB)
and `var_extra.csv.gz`; `bzfig.de` puts the two matrices side by side and runs the
test.

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
cytochrome c oxidase III at +13.79 against +13.80. Where the measured extremes
fall short, the measurement cannot see the point: 5C's most negative gene has an
adjusted p-value that underflows to zero and is drawn as a triangle at the axis
top, and 5E's thirteen genes between −7.96 and −5.71 all have adjusted p-values
between 0.7 and 1.0, so they sit on the x-axis baseline. 5F, the one case with no
low-y or off-axis points involved, matches to the last digit printed.

### The subsets

`orig_ident` carries all three cohorts (Nonreactivated 6505, me49 Day 3 950,
me49 Day 0 602).

| Panel | Group A (numerator) | Group B (reference) |
| --- | --- | --- |
| 5C | in vivo bradyzoites, `Nonreactivated` (6505) | in vitro bradyzoites, `me49 Day 3` (950) |
| 5E | in vitro bradyzoite G1, `me49 Day 3` (280) | in vitro tachyzoite G1, `me49 Day 0` (272) |
| 5F | in vivo bradyzoite G1, `Nonreactivated` (1976) | in vitro bradyzoite G1, `me49 Day 3` (280) |

"G1" is `cc_phase` in G1a or G1b, not `transferred_cc_phase`: the transferred
column gives 4547 against 753 cells for 5F and 673 / 1338 genes, 192 off the
published 676 / 1146.

### Why the 8322-gene universe

Normalisation ran *before* the gene subset, and per source object: each cell
divided by its own total over that object's full gene set, scaled to 1e4 and
`log1p`'d. Applied to the 8170 genes the deposited object kept, that reproduces
its `logcounts` layer to within one float32 ulp. Nothing else does.

| Denominator | Largest difference from the deposited `logcounts` |
| --- | --- |
| each cell's total in its own source object — 8322 genes for the in vivo 10x run, 8496 for the me49 object | **4.8e-07**, one float32 ulp |
| one 8322-gene total for every cell | 5.3e-03 — wrong for the me49 cells |
| the 8170 deposited genes | 0.27 — wrong for every cell |

The deposited object cannot reproduce its own normalisation, so the published
test ran on something wider. `scripts/export_dataset.py` measures all three on
every run, records them in `MANIFEST.json` under `extra_genes`, and refuses to
write the package if the first row stops holding.

### The cutoff and the both-groups filter

All three published panels have a clean empty band around zero — a zero-ink
rectangle spanning the full height, not a thinning. Measured at 600 dpi its edges
are −1.99 / +2.01 (5C), −2.20 / +2.23 (5E) and −2.02 / +2.03 (5F): the cutoff is
applied to the plotted fold change, and it is 2.0. A 1.5-fold cutoff (log2 =
0.585) is excluded both by that measurement and by the counts, which it puts at
2326 up / 2885 down for 5C against a published 664 / 1443. It remains possible as
an upstream pre-filter, since a gene failing |log2FC| > 0.585 also fails
|log2FC| > 2.

The both-groups filter is forced rather than chosen. Scanpy scores a gene with a
zero group mean against a 1e-9 pseudocount, which drops it into a degenerate
bucket: without the filter 5C's `logfoldchanges` run −28.19 … +18.24 and the down
count goes from 1441 to 1879, while the up count does not move. Dropping the
genes undetected in one group restores the published axis range.

### What is not reproduced

* **The point colouring is only partly rule-derivable.** The published legends
  say pink = hypothetical protein and grey = ribosomal protein (blue = cyst wall
  protein in 5E), but pink is applied to the *called-out* genes only: 1024 of 5C's
  2107 plotted genes are hypothetical proteins and the published panel has about
  25 pink points, several of them named genes (sag1, srs2, srs22a, the
  cytochromes). The panels here colour the called-out set — the genes beyond the
  caption's annotation thresholds — pink, except ribosomal proteins, which are
  grey; that rule puts grey on exactly the five points the published 5F has grey.
  5E's blue set is not derivable (two of its genes are a dense granule and a
  rhoptry protein), so the seven genes were read back off the published panel by
  matching points to fold changes, and are recorded in
  `bzfig.constants.SUPP5E_CYST_WALL_GENES`.
* **The labels.** The caption's rule — annotate above log2FC 8 or below −7 for
  5C, 6 and −4 for 5E — is what `bzfig.constants.VOLCANO_LABEL_RANGE` uses (5F's
  caption gives no rule; its labelled genes fit 5C's). The published panels also
  label a hand-picked set of transcription factors (bfd1, the ap2s) and the two
  enolases, which no threshold recovers, and shorten long descriptions by hand.
* **Points off the top of the axis.** The published y limits are 1e-100, 1e-60 and
  1e-200, which is what these panels use. 100 of 5C's 2107 points have an adjusted
  p-value below 1e-100 — 20 of them underflow float64 to zero — and are not
  plotted, as in the published panel. They remain in the DE table, so nothing is
  lost from the numbers; only the drawing clips. 5E and 5F have no such points.
  One consequence when comparing panels side by side: srs2/p35 is among 5C's
  clipped genes here, while the published panel places it at about 1e-88 and
  labels it. The reproduction puts its adjusted p-value slightly lower; the fold
  change agrees.
* **ORF F in 5F.** `TGME49_302005` reaches log2FC +8.41 there but only
  pvals_adj = 0.23, so it is not among the plotted genes. The published 5F does not
  label it either; it labels ORF F in 5C, where the reproduction gives
  pvals_adj = 1.7e-04 and puts the point at +8.30, which is where the published
  panel has it.
* **5E has no published counts** to check against, and its negative extreme cannot
  be measured for the reason given above. It is validated instead by its called-out
  genes, which land where the published panel has them: sag1 at −5.67 (published
  −5.68, both at 1e-60), sag4/srs35a at +9.46, bag1 at +9.77, the helicase at
  −4.01, the CMGC kinase at +6.23, and every one of the blue cyst wall points
  within 0.02 of a published point.

`scripts/export_dataset.py` re-runs all three comparisons from the data package it
has just written and fails if the counts move. They are recorded in the manifest
under `verification.supplementary_5_de` and in
`bzfig.constants.SUPP5_DE_COUNTS_REPRODUCED`; the published counts stay in
`SUPP5C_DE_COUNTS` and `SUPP5F_DE_COUNTS`.

## Figure 1D — a recorded constant

**Figure 1D**, "# Unique markers/Cluster" (38 / 93 / 80 / 881 / 14 / 47), was
hard-coded in the analysis notebook. The underlying quantity is "genes that are a
significant Wilcoxon marker of exactly one cluster", but re-running that test does
not recover the six numbers:

| Cluster | 0 | 1 | 2 | 3 | 4 | 5 |
| --- | --- | --- | --- | --- | --- | --- |
| Published | 38 | 93 | 80 | 881 | 14 | 47 |
| Re-derived, 8170 deposited genes | 37 | 94 | 84 | 861 | 18 | 52 |
| Re-derived, full 8322 genes | 45 | 93 | 84 | 863 | 18 | 53 |

The differences are small and go in both directions, and nothing suggests a
different test. The published values are kept in
`bzfig.constants.UNIQUE_MARKERS_PER_CLUSTER` and the panel is drawn from them, so
the rendered panel matches the paper exactly. If the original cell turns up, the
constant can be replaced by the computation.

## Not included

**Supplementary Figure 2A–F** (cell-cycle-dependent expression of marker genes) is
not derived from this dataset: the caption attributes it to Benke et al., mined
from ToxoDB.

**Figure 3G and 3H** come from flow cytometry and imaging. The cartoons
above Figure 6's UMAPs and beside Supplementary 1B are BioRender artwork.

## The enolase labels

The two enolase panels follow the ToxoDB-65 annotation:

* `TGME49_268850` = **enolase 2**, the tachyzoite isoform — detected in 0.6% of in
  vivo cells, the near-blank panel in the tachyzoite-marker row
* `TGME49_268860` = **enolase 1**, the bradyzoite isoform — detected in 68.9%, the
  broadly stained panel in the bradyzoite-marker row

Two independent sources agree. The ToxoDB annotations used in this study, releases
65 and 68, in both GFF and GTF, all carry `TGME49_268850 = enolase 2` and
`TGME49_268860 = enolase 1`, as do the `var` table of the deposited object and
NCBI Gene. And the primary literature gives the same assignment: Dzierszinski et
al. (1999) originally characterised ENO2 as the tachyzoite isoform and ENO1 as
bradyzoite-specific; Ferguson et al. (2002), *Int J Parasitol* 32:1399–1410,
localised ENO1 to brain tissue cyst bradyzoites and ENO2 to replicating stages by
immuno-EM in vivo; and Ngô et al. (2015), *Acta Cryst D*, state the gene IDs
explicitly — TgENO1 = TGME49_268860, TgENO2 = TGME49_268850.

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
