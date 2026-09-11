# What reproduces, and what does not

Every panel below was checked against the published figure. This page records
the verdict for each, so nobody has to rediscover the awkward cases.

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

## Reproduced as a reconstruction

**Supplementary Figure 4** (CCC and MCC cell-cycle regulator heatmaps).

The code that drew this panel is not in the analysis repository, and the CCC cell
selection was made interactively in a dashboard rather than in a script — it
survives only as a `pseudotime_UCC` column, which is carried into this dataset as
`obs["cell_cycle_group"]`. The recipe here was derived by matching against the
published panel:

* mean per-gene Pearson r = **0.9998** (CCC) and **0.9997** (MCC)
* the published 50-row order of **both** panels is reproduced exactly by the
  sort rule used here (by peak phase, then peak height)

The 50-gene list was transcribed from the published figure; the gene list itself
lives in a package that is not part of this analysis. `TGME49_291590` is labelled
"MYND-like zinc finger protein" in the figure but "hypothetical protein" in the
annotation shipped here.

Two things worth knowing when reading the panel: the colour bar axis runs −4…4
but the values only span **−1.77…+1.79**, and the CCC group is thinly populated
(G1a 11, G1b 230, S 106, M 56, C 44 cells), so several genes are detected in one
phase only and their rows sit at the one-hot extremes of ±1.789 / −0.447.

## Not reproducible — recorded as constants

**Figure 1D**, "# Unique markers/Cluster" (38 / 93 / 80 / 881 / 14 / 47).

Hard-coded in the analysis notebook. The underlying quantity is evidently "genes
that are a significant Wilcoxon marker of exactly one cluster", but the original
run used a larger gene universe (8322 genes) than the deposited object (8170).
Re-running the obvious test on the deposited object gives `[37, 94, 84, 861, 18,
52]` — close for clusters 0 and 3, not for the others. The published values are
kept in `bzfig.constants.UNIQUE_MARKERS_PER_CLUSTER` and the panel is drawn from
them.

**Supplementary Figure 5C, 5E and 5F** — the three volcano plots.

These are not regenerated. Three independent reasons:

1. **The gene universe is wrong.** Panels C and F label `ORF F`, `cytochrome b`
   and `cytochrome c` at the positive extreme. The deposited object contains only
   the 14 nuclear chromosomes (`var["seqid"]`), so the mitochondrial transcripts
   carrying those points are absent. The published test ran before that
   reduction.
2. **The fold-change metric is not scanpy's.** Scanpy's `logfoldchanges` for
   panel C span −28…+8; the published axis spans about −11…+14. A log2 ratio of
   mean normalised expression with a small pseudocount reproduces the negative
   extreme and the labelled genes well, but cannot reach +14 without the
   mitochondrial genes.
3. **An additional fold-change filter was applied** that is not described in the
   caption — the published volcanoes have an empty band from about −1.7 to +1.9.
   Threshold combinations can be fitted to land near the stated counts (705/1460
   against a target of 664/1443), but no round-number threshold reproduces both
   panels, so any such choice would be fitted rather than derived.

No differential-expression code for these panels exists in the analysis
repository; the volcanoes were drawn in a separate plotting program from an
exported table that was not kept. The published counts are recorded in
`bzfig.constants.SUPP5C_DE_COUNTS` and `SUPP5F_DE_COUNTS`.

The subsets are documented, for anyone who wants to run their own test:

| Panel | Group A | Group B |
| --- | --- | --- |
| 5C | in vivo bradyzoites (6505) | in vitro bradyzoites, me49 Day 3 (950) |
| 5E | in vitro bradyzoite G1 | in vitro tachyzoite G1, me49 Day 0 |
| 5F | in vivo bradyzoite G1 | in vitro bradyzoite G1 |

For 5E and 5F it is undetermined whether G1 was taken from `cc_phase` or
`transferred_cc_phase`.

## Not included

**Supplementary Figure 1B–F** (cell-cycle-dependent expression of marker genes)
is not derived from this dataset at all — the caption attributes it to Benke et
al., mined from ToxoDB.

## A label correction worth knowing

An early cell of the analysis notebook maps the two enolase labels the wrong way
round (`eno-1 → TGME49_268850`, `eno-2 → TGME49_268860`). The gene annotation and
the published panels both use the opposite mapping, which is what this repository
follows:

* `TGME49_268850` = **enolase 2**, detected in 0.6% of in vivo cells — the
  near-blank panel in the tachyzoite-marker row
* `TGME49_268860` = **enolase 1**, detected in 68.9% — the broadly stained panel
  in the bradyzoite-marker row

Some stale PNGs in the analysis repository carry the swapped names. They are not
the naming authority.

## Panel assembly

The published figures were assembled in a vector editor. Panel letters, the
Figure 1E grid and its row labels, the single shared colour bar across 1E/1F, the
in-plot cluster numbers in Figure 1B and the coloured row-group bands beside
Figure 1C are all added at that stage. This repository renders the individual
panels as the analysis produced them; it does not attempt the assembly.
