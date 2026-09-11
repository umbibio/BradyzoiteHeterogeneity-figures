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
# published caption. The volcano panels themselves are not reproducible from the
# deposited object -- see docs/reproducibility.md.
SUPP5C_DE_COUNTS = {"up_in_vivo": 664, "down_in_vivo": 1443}
SUPP5F_DE_COUNTS = {"up_in_vivo": 676, "down_in_vivo": 1146}

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
