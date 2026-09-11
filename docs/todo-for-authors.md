# What I need from the authors

Items found while making the figure panels reproducible. Some are now closed;
the ones that remain say what is missing, why it matters, and what would close
it. Nothing here blocks the repository from building the figures.

## Status

| | Item | State |
| --- | --- | --- |
| 1 | Supplementary 5C/5E/5F volcano plots | **Resolved — regenerated; two things to confirm** |
| 2 | Figure 1D unique-marker counts | Accepted as a recorded constant |
| 3 | Enolase gene assignment | **Resolved — no error in the figure** |
| 4 | Supplementary Figure 4 reconstruction | Needs a check |
| 5 | Figure 1C row-group bands | Needs a check |
| 6 | CCC/MCC vs UCC/BCC naming | Needs a consistency pass |
| 7 | Public data URL | Waiting on you, nothing blocked |

---

## 1. The differential expression behind Supplementary Figure 5C, 5E, 5F

**Status:** resolved. The three volcanoes are regenerated from the data package
and rendered by `scripts/make_figures.py`; no exported table is needed any more.

**What was wrong before.** The gene universe. The published test ran on all 8322
ToxoDB-65 genes; the deposited object keeps the 8170 on the nuclear chromosomes,
and the other 152 — on unplaced `KE*` contigs — carry the apicoplast and
mitochondrial transcripts (`ORF F`, the cytochrome b's, cytochrome c oxidase III)
that the panels label at the positive extreme. Without them the fold-change axis
could not reach +14 and the test was being run on a different question. Those 152
genes now ship as `logcounts_extra.mtx.gz`, recovered from the two source objects
the deposited one was integrated from.

The decisive point is that the deposited object cannot reproduce its own
normalisation: dividing each cell by its total over the 8170 genes it holds is
wrong by up to 0.27 in `logcounts`, while dividing by its total in its own source
object is right to one float32 ulp. Normalisation therefore ran before the gene
subset, and the test after it.

| | Published | Reproduced |
| --- | --- | --- |
| 5C up / down in vivo | 664 / 1443 | 666 / 1441 |
| 5F up / down in vivo | 676 / 1146 | 678 / 1145 |
| 5E | no counts quoted | 55 / 69 |

Two genes in each panel, out of two thousand. The full recipe, the evidence for
the fold-change cutoff of 2, and the fold changes of individual labelled genes
measured off the printed figure are in `docs/reproducibility.md`.

**Your question about "G1" is answered.** It is `cc_phase`, not
`transferred_cc_phase`: the transferred column gives 673 / 1338 genes for 5F,
192 off the published 676 / 1146. `cc_phase` puts 1976 in vivo against 280 in
vitro cells in that comparison.

**Two things to confirm, neither of them blocking:**

- **ORF F in panel F.** `TGME49_302005` comes out at log2FC +8.41 but
  pvals_adj = 0.23 in the 5F comparison, so the panel here does not plot it. The
  published 5F does not label it either — ORF F is labelled in 5C, where it is
  significant (pvals_adj = 1.7e-04) and lands at +8.30, exactly where the printed
  panel has it. If your 5F does include a point for ORF F, that one gene is
  unexplained; an unlabelled point cannot be identified from the printed figure.
- **Panel E has no published counts.** The caption quotes counts for 5C and 5F
  but not for 5E, so there is nothing to check 55 / 69 against. The panel agrees
  with the published one everywhere it can be checked — sag1 at −5.67 and
  1e-60, sag4 at +9.46, bag1 at +9.77, the helicase at −4.01, the cyst wall
  points all within 0.02 — but a count from you would settle it.

**One thing worth saying in the methods** if the panels are described again: a
1.5-fold pre-filter is compatible with everything here but is not the plotted
cutoff. The empty band around zero measures −1.99 … +2.01, i.e. |log2FC| > 2; a
1.5-fold cutoff (log2 = 0.585) would put 2326 / 2885 genes on panel C instead of
664 / 1443. A 1.5-fold pre-filter upstream would change none of the counts,
since everything it removes the 2.0 cutoff removes as well.

---

## 2. Where the Figure 1D numbers came from

**Status:** accepted as a recorded constant. Recorded here for completeness — no
action needed unless the provenance turns up.

The "# Unique markers/Cluster" bars — **38 / 93 / 80 / 881 / 14 / 47** — are
hard-coded in the analysis notebook with no accompanying computation.

The quantity is evidently "genes that are a significant Wilcoxon marker of
exactly one cluster". Re-running that on the deposited object gives
`[37, 94, 84, 861, 18, 52]` — close for clusters 0 and 3, off for the others.
The gene universe was the suspect, but it is not the answer: now that the full
8322 genes are shipped (see item 1), the same test on them gives
`[45, 93, 84, 863, 18, 53]`, which is no closer.

**Decision:** the published values stay as a recorded constant in
`bzfig.constants.UNIQUE_MARKERS_PER_CLUSTER`, and the panel is drawn from them.
`docs/reproducibility.md` shows the re-derived numbers alongside so a reader can
see how close they are. If the original cell ever turns up, drop it in and the
panel becomes fully derived.

---

## 3. Enolase gene assignment — RESOLVED, no error in the figure

This was raised as a possible error in the published figure. It is not: the
figure is correct.

An early cell of the analysis notebook (cell 142) maps the enolase labels the
wrong way round. That cell was never used for the published panel — cell 145
was, and it agrees with everything else:

| | `TGME49_268850` | `TGME49_268860` |
| --- | --- | --- |
| Notebook cell 142 (unused) | eno-1 | eno-2 |
| Notebook cell 145 (used for the figure) | enolase 2 | enolase 1 |
| ToxoDB release 65, GFF and GTF | enolase 2 | enolase 1 |
| ToxoDB release 68, GFF and GTF | enolase 2 | enolase 1 |
| NCBI Gene | enolase 2 | enolase 1 |
| Stage specificity in the literature | **tachyzoite** | **bradyzoite** |
| Detected in in vivo bradyzoites | **0.6%** | **68.9%** |
| Printed in the published panel under | Tachyzoite markers | Bradyzoite markers |

Every line agrees. The assignment is also confirmed independently of the
annotation, in the primary literature:

- Dzierszinski et al. (1999) — the original characterisation: ENO2 is expressed
  in tachyzoites, ENO1 solely in bradyzoites.
- Ferguson et al. (2002), *Int J Parasitol* 32:1399–1410 — ENO1 in brain tissue
  cyst bradyzoites, ENO2 in replicating stages, by immuno-EM in vivo.
- Ngô et al. (2015), *Acta Cryst D* — the bradyzoite enolase structure paper,
  which gives the gene IDs explicitly: **TgENO1 = TGME49_268860**,
  **TgENO2 = TGME49_268850**.

No action needed on the manuscript. Two small tidy-ups if you want them:

- Some PNGs in the analysis repository are named from the swapped run
  (`expression_eno-1_TGME49_268850.png` and similar). They are the only
  artefacts carrying the wrong names; worth deleting so nobody is misled later.
- If the readership includes non-parasitologists, a footnote may help: in humans
  the ENO1/ENO2 numbering means something entirely different (ENO1 is the
  ubiquitous isoform there), which is a known trap for readers crossing over
  from mammalian work. It is not an inconsistency within the *Toxoplasma*
  literature, which uses one convention throughout.

---

## 4. Supplementary Figure 4 — confirm the reconstruction

The code that drew this panel does not exist anywhere, and the CCC cell
selection was made interactively in a dashboard rather than in a script. I
reconstructed it and matched it against the published image: mean per-gene
Pearson **r = 0.9998** (CCC) and **0.9997** (MCC), with the published 50-row
order of both panels reproduced exactly. I am confident it is right, but three
details deserve a human eye:

- **The 50-gene list was transcribed from the printed figure**, because the list
  itself lives in a package that is not part of this analysis. Please check it
  against whatever you originally used. It is in
  `src/bzfig/constants.py::SUPP4_GENES`.
- **`TGME49_291590`** is labelled "MYND-like zinc finger protein" in the figure
  but "hypothetical protein" in the annotation we ship. Which is right?
- **The pink highlighting is not fully rule-derivable.** 17 of the 20 pink genes
  are exactly those peaking in the same phase in both panels. The other three —
  `TGME49_207900`, `TGME49_262730`, `TGME49_318470` — peak one phase apart. The
  caption's "e.g. highest within the same cell cycle state" covers this, but it
  was a judgement call; please confirm the intended rule, or supply the list.

**One thing to consider mentioning in the caption:** the colour bar axis runs
−4…4, but no value in either panel exceeds **±1.79**. And the CCC group is thin
— G1a 11, G1b 230, S 106, M 56, C 44 cells — so several genes are detected in a
single phase and their rows sit at the one-hot extremes. A reader cannot tell
that from the figure.

---

## 5. Confirm the Figure 1C row-group bands

The coloured bands beside the heatmap were added during figure assembly, so
they are not in any code. I read them off the published panel by measuring the
pixel runs: **CWPs 5, Micronemes-1 7, Micronemes-2 6, ~Tachyzoites 6,
Rhoptry-1 5, Rhoptry-2 5**.

The boundary worth double-checking is `TGME49_306270` (a hypothetical protein),
which my measurement puts as the last row of **Micronemes-1**. Please confirm,
since a one-row error here would misattribute a gene in the source data table.

Also: the panel prints the first band as **CWPs**; the source data spells it
"Cyst wall proteins". Say if you want the abbreviation used verbatim.

---

## 6. Naming: CCC/MCC vs UCC/BCC

The deposited pseudotime file uses **UCC** and **BCC**; the published figure uses
**CCC** (common cell cycle) and **MCC** (modified cell cycle). I have carried the
published names into the dataset as `obs["cell_cycle_group"]` and kept the
original pseudotime columns under their own names.

Supplementary Figure 5's caption still says "regulators of the BCC". Worth a
consistency pass across the manuscript so one pair of names is used throughout.

---

## 7. Where the data will live

The pipeline fetches its input from a list of mirrors with SHA-256 verification,
so it works with any number of URLs per file. Right now the data is in git LFS
on the internal server, which is enough for us but not for a reader.

**What I need when you have it:** the public URL or DOI — Zenodo, an
institutional server, whatever you choose. Send me the base URL and I will fill
in `data/MANIFEST.json`; readers then need nothing but the repository.

Nothing is blocked on this in the meantime.

---

## 8. For information — no action needed

- **Supplementary Figure 1B–F** is not derived from this dataset at all; the
  caption attributes it to Benke et al., mined from ToxoDB. It is out of scope
  for this repository.
- **Panel assembly** — letters, the Figure 1E grid and its row labels, the shared
  colour bar across 1E/1F, the in-plot cluster numbers — was done in a vector
  editor. This repository renders the individual panels as the analysis produced
  them and does not attempt the assembly.
