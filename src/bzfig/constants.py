"""Values the published figures depend on that are not derivable from the data.

Anything here is a *recorded constant*: it was fixed by the authors when the
figure was made, and reproducing the panel requires knowing it. Each entry says
where it came from.
"""

from __future__ import annotations

# Cluster palette, clusters 0-5. From the analysis notebook
# (integrated_adata_me49_nr_subset-2026-06-29.ipynb, cell 132).
CLUSTER_COLORS = ["#e0007a", "#4f7f6a", "#9c7c38", "#6f5e8d", "#00a6d6", "#7d5a5a"]

# Cell cycle phases, in plotting order.
PHASES = ["G1a", "G1b", "S", "M", "C"]

# The two phases the published "G1" subsets of Supplementary 5E and 5F are made
# of, taken from obs["cc_phase"] (not transferred_cc_phase: that gives 4547 vs
# 753 cells for 5F and 673/1338 genes, against a published 676/1146).
G1_PHASES = ("G1a", "G1b")

# Expression colour scale shared by every per-gene UMAP panel (Fig 1E, 1F, Supp 5D).
# vmin is -vmax * 0.2; vmax is pinned to 5 regardless of the data range.
EXPRESSION_CMAP = "PuRd"
EXPRESSION_VMIN = -1.0
EXPRESSION_VMAX = 5.0

# Heatmap colour ceiling for the marker panels (Fig 1C, Supp 1A).
HEATMAP_VMAX = 6

# Figure 1D, "# Unique markers/Cluster".
#
# NOT DERIVABLE from the published data object. The notebook hard-codes these
# (cell 126). The underlying test is evidently "genes that are a significant
# Wilcoxon marker of exactly one cluster", but the original run used a larger
# gene universe (8322 genes) than the deposited object (8170), so the counts
# cannot be regenerated exactly. Re-running on the deposited object gives
# [37, 94, 84, 861, 18, 52]; see docs/reproducibility.md.
UNIQUE_MARKERS_PER_CLUSTER = {0: 38, 1: 93, 2: 80, 3: 881, 4: 14, 5: 47}

# Supplementary Figure 5C and 5F differential-expression counts, quoted in the
# published caption. bzfig.de reproduces them to within two genes -- see
# docs/reproducibility.md.
SUPP5C_DE_COUNTS = {"up_in_vivo": 664, "down_in_vivo": 1443}
SUPP5F_DE_COUNTS = {"up_in_vivo": 676, "down_in_vivo": 1146}

# Significance and fold-change cutoffs for the Supplementary 5 volcano panels.
#
# DE_ALPHA is the caption's own ("only displayed are those with a p value <
# 0.05 from a Wilcoxon rank-sum test"); the axis is adjusted p, so it is applied
# to pvals_adj. DE_LOG2FC_CUTOFF is not stated anywhere, but every published
# panel has an empty band around zero, measured on the printed 5C axis at
# -1.9..+1.9 (300 dpi, 60.35 px per log2 unit). 2.0 reproduces the published
# counts to within two genes in both panels that quote them; a 1.5-fold cutoff
# (log2 = 0.585) would give 2326 up / 2885 down for 5C against 664 / 1443.
DE_ALPHA = 0.05
DE_LOG2FC_CUTOFF = 2.0

# What bzfig.de produces from the shipped data package, "up" meaning higher in
# group A. Recorded so that a change in the data or in scanpy shows up as a
# failed check in scripts/export_dataset.py; 5E has no published counts.
SUPP5_DE_COUNTS_REPRODUCED = {
    "5C": {"up": 666, "down": 1441},
    "5E": {"up": 55, "down": 69},
    "5F": {"up": 678, "down": 1145},
}

# Volcano panel axes, measured off the published figure at 300 dpi from the tick
# spacing and the ends of the axis lines. The y axis is adjusted p on a reversed
# log scale, labelled in powers of ten.
VOLCANO_AXES = {
    "5C": {"xlim": (-12.1, 14.1), "xticks": (-10, 0, 10), "ymax": 100, "ystep": 20},
    "5E": {"xlim": (-11.0, 15.1), "xticks": (-10, -5, 0, 5, 10, 15), "ymax": 60, "ystep": 20},
    "5F": {"xlim": (-16.3, 15.1), "xticks": (-10, 0, 10), "ymax": 200, "ystep": 50},
}

# Which genes the published panels call out. The caption gives the rule for 5C
# ("all genes with a log2fold change > 8 or <-7 are annotated") and for 5E
# ("log2FC > 6 or < -4"); 5F states none and its labelled genes fit 5C's rule.
# The published panels additionally label a hand-picked set of transcription
# factors (bfd1, the ap2s) and the two enolases, which no rule recovers.
VOLCANO_LABEL_RANGE = {"5C": (-7.0, 8.0), "5E": (-4.0, 6.0), "5F": (-7.0, 8.0)}

# Point colours, from the published legends: the called-out genes are pink,
# except ribosomal proteins, which are grey, and the cyst wall proteins of
# panel 5E, which are blue. Everything else is black.
VOLCANO_COLORS = {
    "other": "#000000",
    "called out": "#fc4ba0",
    "ribosomal protein": "#aaaaaf",
    "cyst wall protein": "#76b8fd",
}

# The genes drawn blue in the published 5E. There is no cyst wall protein list
# in the data or the paper (two of these are a dense granule and a rhoptry
# protein), so the set was read back off the published panel by matching the
# blue points to their fold changes -- the same way the Figure 1C row-group
# bands were recovered.
SUPP5E_CYST_WALL_GENES = [
    "TGME49_209755",
    "TGME49_213067",
    "TGME49_242110",
    "TGME49_251540",
    "TGME49_260520",
    "TGME49_264660",
    "TGME49_312330",
]

# Figure 1E / 1F genes. The label is the one printed on the panel.
#
# Note: an earlier cell of the analysis notebook (cell 142) maps the two enolase
# labels the other way round. Cell 145 -- used here -- agrees with the gene
# annotation and with the published panels: TGME49_268850 ("enolase 2") is
# detected in 0.6% of cells, matching the near-blank tachyzoite-row panel, while
# TGME49_268860 ("enolase 1") is detected in 68.9%, matching the broadly
# labelled bradyzoite-row panel.
FIG1E_GENES = [
    ("Tachyzoite markers", "sag1", "TGME49_233460"),
    ("Tachyzoite markers", "enolase2", "TGME49_268850"),
    ("Tachyzoite markers", "ldh1", "TGME49_232350"),
    ("Bradyzoite markers", "bag1", "TGME49_259020"),
    ("Bradyzoite markers", "enolase1", "TGME49_268860"),
    ("Bradyzoite markers", "ldh2", "TGME49_291040"),
]

FIG1F_GENE = ("cst1 (srs44)", "TGME49_264660")
SUPP5D_GENE = ("srs22a", "TGME49_238440")

# Supplementary Figure 1B and 1C: eleven more per-gene expression UMAPs, the
# same panel type as Figure 1E, in the order the panels print them. Every id was
# checked against `var` and against the per-gene PNGs the analysis left behind
# (notebooks/figures/expression_<NAME>_<id>.png). cst4 and cst10 carry no name
# in the ToxoDB-65 annotation shipped here -- both are "hypothetical protein" --
# so their labels rest on those filenames and on the caption. Panel 1B also
# carries a BioRender cartoon, which is figure assembly, not a panel.
SUPP1_GENES = [
    ("1B", "mic2", "TGME49_201780"),
    ("1B", "m2ap", "TGME49_214940"),
    ("1B", "mic3", "TGME49_319560"),
    ("1B", "ama4", "TGME49_294330"),
    ("1B", "ama2", "TGME49_300130"),
    ("1C", "cst4", "TGME49_261650"),
    ("1C", "cst10", "TGME49_312330"),
    ("1C", "mcp4", "TGME49_208730"),
    ("1C", "bpk1", "TGME49_253330"),
    ("1C", "mag1", "TGME49_270240"),
    ("1C", "mag2", "TGME49_209755"),
]

# Figure 3 and Figure 6 pick one cohort out of the projection. The grey layer
# and the flat in vitro highlight of Figure 3D are from the analysis notebook
# (integrated_adata_me49_nr_subset-2026-06-29.ipynb, cell 172); the phase
# palette is not a constant, it travels with the data as uns["cc_phase_colors"].
BACKGROUND_GRAY = "#d3d3d3"
FIG3D_IN_VITRO = "#FA8072"

# The three cohorts of orig_ident that Figure 3A-C and Figure 6 draw, with the
# Figure 3 panel letter, the name used for the rendered file, and whether that
# Figure 3 panel carries the grey layer. 3A does not: it is the same rendering
# as Figure 6's first panel, the cohort alone. 3B and 3C are drawn over the
# whole projection in grey. Figure 6 never draws the grey layer.
COHORTS = [
    ("3A", "nonreactivated", "Nonreactivated", False),
    ("3B", "me49_day0", "me49 Day 0", True),
    ("3C", "me49_day3", "me49 Day 3", True),
]

# What the printed Figure 3 phase bars say, measured at 300 dpi and scaled so the
# five bars sum to the cohort size. 3B agrees with the deposited
# transferred_cc_phase exactly; 3A is about 1% off and 3C matches nothing in
# either repository (see docs/reproducibility.md). The panels here are drawn
# from the deposited labels, so this is recorded only to say by how much.
FIG3_PHASE_BARS_PUBLISHED = {
    "3A": {"G1a": 1516, "G1b": 3022, "S": 1350, "M": 318, "C": 299},
    "3B": {"G1a": 162, "G1b": 161, "S": 190, "M": 64, "C": 25},
    "3C": {"G1a": 186, "G1b": 494, "S": 191, "M": 50, "C": 29},
}

# Figure 1C: genes in plotting order, with the row-group bands drawn beside the
# heatmap (analysis notebook cell 156; band assignment read off the published
# panel, giving band sizes 5/7/6/6/5/5).
FIG1C_GENES: list[tuple[str, str, str]] = [
    ("Cyst wall proteins", "SAG-related sequence SRS44", "TGME49_264660"),
    ("Cyst wall proteins", "cyst matrix protein MAG2", "TGME49_209755"),
    ("Cyst wall proteins", "bradyzoite pseudokinase BPK1", "TGME49_253330"),
    ("Cyst wall proteins", "Toxoplasma gondii family A protein", "TGME49_278080"),
    ("Cyst wall proteins", "SAG-related sequence SRS49D", "TGME49_207160"),
    ("Micronemes-1", "hypothetical protein", "TGME49_287040"),
    ("Micronemes-1", "hypothetical protein", "TGME49_289370"),
    ("Micronemes-1", "microneme protein MIC12", "TGME49_267680"),
    ("Micronemes-1", "SAG-related sequence SRS22A", "TGME49_238440"),
    ("Micronemes-1", "AMA2", "TGME49_300130"),
    ("Micronemes-1", "AMA4", "TGME49_294330"),
    ("Micronemes-1", "hypothetical protein", "TGME49_306270"),
    ("Micronemes-2", "microneme protein MIC10", "TGME49_250710"),
    ("Micronemes-2", "hypothetical protein", "TGME49_257970"),
    ("Micronemes-2", "dense granule protein GRA64", "TGME49_202620"),
    ("Micronemes-2", "perforin-like protein PLP1", "TGME49_204130"),
    ("Micronemes-2", "MIC2", "TGME49_201780"),
    ("Micronemes-2", "MIC2-associated protein M2AP", "TGME49_214940"),
    ("~Tachyzoites", "Ctr copper transporter family protein", "TGME49_262710"),
    ("~Tachyzoites", "microtubule associated protein SPM2", "TGME49_286590"),
    ("~Tachyzoites", "beta-tubulin, putative", "TGME49_221620"),
    ("~Tachyzoites", "inner membrane complex protein IMC29", "TGME49_243200"),
    ("~Tachyzoites", "LDH1", "TGME49_232350"),
    ("~Tachyzoites", "hypothetical protein", "TGME49_266300"),
    ("Rhoptry-1", "rhoptry protein ROP1", "TGME49_309590"),
    ("Rhoptry-1", "WAVE complex interacting protein WIP", "TGME49_247520"),
    ("Rhoptry-1", "rhoptry protein ROP10", "TGME49_315490"),
    ("Rhoptry-1", "protein phosphatase PP2C-hn", "TGME49_282055"),
    ("Rhoptry-1", "rhoptry protein ROP54", "TGME49_210370"),
    ("Rhoptry-2", "rhoptry kinase family protein", "TGME49_308096"),
    ("Rhoptry-2", "Kinase-like domain-containing protein", "TGME49_321700"),
    ("Rhoptry-2", "rhoptry protein ROP42", "TGME49_209985"),
    ("Rhoptry-2", "rhoptry kinase family protein", "TGME49_308093"),
    ("Rhoptry-2", "rhoptry neck protein RON2L1", "TGME49_294400"),
]

# Supplementary Figure 4: the known cell-cycle regulators shown in both panels,
# transcribed from the published figure (the gene list lives in a package that is
# not part of this analysis).
SUPP4_GENES = [
    "TGME49_" + gene
    for gene in (
        "289710 239010 239910 250800 224230 293280 280470 207900 223120 266910 "
        "290630 313040 311510 218220 260250 239440 223050 243780 291590 267580 "
        "254910 220440 229200 216220 261410 249880 236240 318470 266900 219100 "
        "240460 203710 305330 281450 251740 290020 270330 264690 262730 229020 "
        "256070 237090 254630 249260 309410 247700 315760 224050 271280 292010"
    ).split()
]
