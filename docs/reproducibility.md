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
| Supplementary 1A | All cluster markers per cluster |
| Supplementary 5A | As Figure 1B, with a legend |
| Supplementary 5B | In vitro bradyzoites coloured by transferred cluster identity |
| Supplementary 5D | srs22a expression |

For Figure 1C, Figure 1G, Supplementary 1A and Supplementary 4 the exported
values were additionally re-derived from the data object and re-plotted from the
exported tables alone; both matched the published panels.

## Reproduced from the original R recipe

**Figure 2H and Supplementary Figure 4** use the original R methods,
implemented in `bzfig.figure_2h_supplementary_4` using the same deposited `logcounts`, cell
metadata and gene identifiers. No Seurat object or additional expression
dataset is required. Original row orders and gene labels are recorded in
`figure_2h_supplementary_4_metadata.json`; measured values are calculated, not hard-coded.

The original Seurat normalization was applied after the 8,170-gene subset.
It is recovered by `log1p(10000 * expm1(logcounts) / rowSum(expm1(logcounts)))`,
with row sums over all 8,170 genes before selecting plotted genes. This
conversion is local to these panels and does not modify any shared layer.

Figure 2H calculates Pearson correlations across the 6,505 cells, the first
five genes against all 29. All 145 displayed two-decimal values match the
original. Supplementary 4 computes phase means then sample-SD z-scores of the
five means, separately for CCC/MCC; its 500 z-scores match the audited R values
within 1e-6. Floating-point equality is not claimed: source logcounts is float32.

The independent publication row orders and original LAB-interpolated color
breaks (-2.5/0/2.5) are retained. The CCC counts remain 11/230/106/56/44; MCC
counts are 1517/2789/1219/271/262. The manual pale-pink description backgrounds
are now copied explicitly from page 4 of the final supplementary PDF, with
black text and unhighlighted gene IDs. The source has 19 highlighted CCC rows
and 20 MCC rows: TGME49_315760 (AP2XI-4) is highlighted only in MCC. This source
asymmetry is preserved, not silently corrected or inferred from expression.
Numerical reproduction is distinct from exact page layout.

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
panel type as Figure 1E. They are a gap in this repository, not out of scope.

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
