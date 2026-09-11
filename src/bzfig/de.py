"""The differential expression behind the Supplementary 5 volcano panels.

The published test ran on the full 8322-gene ToxoDB-65 universe, not on the 8170
genes of the deposited object: the missing 152 sit on unplaced contigs and carry
the apicoplast and mitochondrial transcripts the panels label at the positive
extreme. ``logcounts_extra.mtx.gz`` ships them, and widening the matrix back out
is the only thing these panels need that the others do not.

The recipe — Wilcoxon rank-sum on the log-normalised matrix, Benjamini-Hochberg,
genes detected in both groups, ``|log2FC| > 2`` — reproduces the published counts
to within two genes. ``docs/reproducibility.md`` has the evidence, including why
the both-groups filter is forced rather than chosen.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

from .constants import DE_ALPHA, DE_LOG2FC_CUTOFF, G1_PHASES
from .data import DATA, load_extra_genes

GROUP_A, GROUP_B = "A", "B"

COLUMNS = [
    "gene_id",
    "log2fc",
    "pval",
    "pval_adj",
    "mean_a",
    "mean_b",
    "gene_name",
    "gene_description",
    "significant",
]


@dataclass(frozen=True)
class Comparison:
    """One published volcano: ``sample_a`` against ``sample_b`` as reference."""

    label_a: str
    label_b: str
    sample_a: str
    sample_b: str
    g1_only: bool = False


COMPARISONS: dict[str, Comparison] = {
    "5C": Comparison(
        "in vivo bradyzoites", "in vitro bradyzoites", "Nonreactivated", "me49 Day 3"
    ),
    "5E": Comparison(
        "in vitro bradyzoite G1",
        "in vitro tachyzoite G1",
        "me49 Day 3",
        "me49 Day 0",
        g1_only=True,
    ),
    "5F": Comparison(
        "in vivo bradyzoite G1",
        "in vitro bradyzoite G1",
        "Nonreactivated",
        "me49 Day 3",
        g1_only=True,
    ),
}


def expanded(adata: AnnData, datadir: Path = DATA) -> AnnData:
    """``logcounts`` widened back to the 8322 genes the published test used."""
    matrix, var = load_extra_genes(datadir)
    return AnnData(
        X=sp.hstack([adata.layers["logcounts"], matrix], format="csr"),
        obs=adata.obs.copy(),
        var=pd.concat([adata.var[list(var.columns)], var]),
    )


def _groups(obs: pd.DataFrame, comparison: Comparison) -> pd.Series:
    """Group label per cell; the empty string for cells in neither group."""
    sample = obs["orig_ident"].astype(str)
    group = pd.Series("", index=obs.index, dtype=object)
    group[sample == comparison.sample_a] = GROUP_A
    group[sample == comparison.sample_b] = GROUP_B
    if comparison.g1_only:
        group[~obs["cc_phase"].astype(str).isin(G1_PHASES)] = ""
    return group


def volcano_table(adata: AnnData, panel: str, datadir: Path = DATA) -> pd.DataFrame:
    """One row per gene detected in both groups of *panel*'s comparison.

    ``mean_a`` / ``mean_b`` are the mean normalised expression of the gene in
    each group, on the linear scale scanpy's fold change is defined against. A
    gene with a zero mean in one group is dropped: scanpy scores those against a
    1e-9 pseudocount, which throws fold changes out to ±28 and is what put the
    published axis limits out of reach when the panels were first attempted.
    """
    comparison = COMPARISONS[panel]
    wide = expanded(adata, datadir)
    group = _groups(wide.obs, comparison)

    subset = wide[group != ""].copy()
    subset.obs["group"] = pd.Categorical(group[group != ""])
    sc.tl.rank_genes_groups(
        subset,
        "group",
        groups=[GROUP_A],
        reference=GROUP_B,
        method="wilcoxon",
        tie_correct=False,
    )
    table = sc.get.rank_genes_groups_df(subset, GROUP_A).rename(
        columns={
            "names": "gene_id",
            "logfoldchanges": "log2fc",
            "pvals": "pval",
            "pvals_adj": "pval_adj",
        }
    )

    expression = sp.csr_matrix(subset.X).expm1()
    in_a = (subset.obs["group"] == GROUP_A).to_numpy()
    table = table.merge(
        pd.DataFrame(
            {
                "gene_id": list(subset.var_names),
                "mean_a": np.asarray(expression[in_a].mean(0)).ravel(),
                "mean_b": np.asarray(expression[~in_a].mean(0)).ravel(),
            }
        ),
        on="gene_id",
    )
    table = table[(table["mean_a"] > 0) & (table["mean_b"] > 0)]
    table = table.join(wide.var[["gene_name", "gene_description"]], on="gene_id")
    table["significant"] = (table["pval_adj"] < DE_ALPHA) & (
        table["log2fc"].abs() > DE_LOG2FC_CUTOFF
    )
    return table[COLUMNS].reset_index(drop=True)


def counts(table: pd.DataFrame) -> dict[str, int]:
    """Genes up and down in group A — the numbers the caption quotes."""
    significant = table[table["significant"]]
    return {
        "up": int((significant["log2fc"] > 0).sum()),
        "down": int((significant["log2fc"] < 0).sum()),
    }
