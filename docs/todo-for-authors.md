# What I need from the authors

Items found while making the figure panels reproducible. Some are now closed;
the ones that remain say what is missing, why it matters, and what would close
it. Nothing here blocks the repository from building the figures.

## Status

| | Item | State |
| --- | --- | --- |
| 1 | Supplementary 5C/5E/5F volcano plots | Accepted as external — see below |
| 2 | Figure 1D unique-marker counts | Accepted as a recorded constant |
| 3 | Enolase gene assignment | **Resolved — no error in the figure** |
| 4 | Supplementary Figure 4 reconstruction | Needs a check |
| 5 | Figure 1C row-group bands | Needs a check |
| 6 | CCC/MCC vs UCC/BCC naming | Needs a consistency pass |
| 7 | Public data URL | Waiting on you, nothing blocked |

---

## 1. The differential-expression table behind Supplementary Figure 5C, 5E, 5F

**Status:** accepted as external. You have said the genes at the extremes are not
of particular significance and the panels are fine as they stand, so nothing here
is blocking — this is a record of why they are not regenerated, and of the one
thing that would change that.

**Why they cannot be regenerated.** Three independent problems, any one of which
would be enough:

- **The gene universe is gone.** Panels C and F label `ORF F`, `cytochrome b`
  and `cytochrome c` at the positive extreme of the fold-change axis. The
  deposited object contains only the 14 nuclear chromosomes — the mitochondrial
  transcripts carrying those points were filtered out before the object was
  saved. The published test therefore ran on a larger gene set than anything we
  still have.
- **The fold-change metric is unknown.** Scanpy's `logfoldchanges` for panel C
  span −28…+8; the published axis spans about −11…+14. Seurat's default
  `FindMarkers` spans −5.0…+5.1. Neither is it.
- **An undocumented filter was applied.** The published volcanoes have an empty
  band from about −1.7 to +1.9 around zero that the caption does not mention.

Fitted threshold combinations can land near the stated counts (705/1460 against
the published 664/1443), but that is curve-fitting, not reproduction.

**How they are handled.** The three panels are documented as external in
`docs/reproducibility.md`, with the published counts (664/1443 and 676/1146)
recorded in `bzfig.constants.SUPP5C_DE_COUNTS` and `SUPP5F_DE_COUNTS`, and the
cell subsets of each comparison described in prose so anyone can run their own
test.

**The one thing that would change this:** the exported DE table the volcanoes
were plotted from — one row per gene with gene ID, fold change, p-value and
adjusted p-value, per comparison. With it the panels become exactly reproducible,
and it also becomes the source data the publisher wants for them. If it is on a
disk somewhere it is worth ten minutes of looking; if not, no matter.

**Small question, worth answering either way:** for panels 5E and 5F, was "G1"
taken from `cc_phase` or from `transferred_cc_phase`? The two give very different
cell counts (280 vs 753 for the in vitro bradyzoites), and the caption does not
say. This affects how the subsets are described in the methods, independently of
whether the panels are ever re-plotted.

---

## 2. Where the Figure 1D numbers came from

**Status:** accepted as a recorded constant. Recorded here for completeness — no
action needed unless the provenance turns up.

The "# Unique markers/Cluster" bars — **38 / 93 / 80 / 881 / 14 / 47** — are
hard-coded in the analysis notebook with no accompanying computation.

The quantity is evidently "genes that are a significant Wilcoxon marker of
exactly one cluster". Re-running that on the deposited object gives
`[37, 94, 84, 861, 18, 52]` — close for clusters 0 and 3, off for the others.
The gap is almost certainly the gene universe again: the working object had
**8322** genes, the deposited one has **8170**.

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
