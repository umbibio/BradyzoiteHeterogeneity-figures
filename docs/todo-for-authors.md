# Author checks before release

Updated with Kourosh's contribution, 11 September 2026. Rendering a panel
successfully is not the same as confirming all historical analysis settings.
See [contribution and QC](kourosh-contribution.md) for the validation boundary.

| Item | Status / smallest remaining action |
| --- | --- |
| Figure 2H | Numerically reproduced from the shared data; visual approval pending. Confirm that the legend describes Pearson correlations across cells, not across phase means. |
| Supplementary 4 | Original recipe and pink description backgrounds reproduced. Review the source's AP2XI-4 highlight asymmetry noted below. |
| Supplementary 5C/E/F | Existing reconstruction retained; reconcile historical settings and small count differences before calling it exact. |
| Figure 1D | Six recorded counts agree with the final workbook lists; historical marker-selection computation remains unverified. Existing panel unchanged. |
| Figure 1C | Confirm manually assembled row-group boundaries, especially the last Micronemes-1 row. |
| Naming | Check CCC/MCC versus older UCC/BCC terminology in captions. |
| Release | Supply final archive DOI/accessions and choose licensing. A working public download is not a substitute for archival deposition. |

## Supplementary 5: confirmation, not numerical tuning

| Panel | Manuscript up/down | Existing reconstruction up/down |
| --- | --- | --- |
| 5C | 664 / 1443 | 666 / 1441 |
| 5E | not quoted | 55 / 69 |
| 5F | 676 / 1146 | 678 / 1145 |

Argenis: please confirm the original DE table or code and historical gene
universe, G1 column, filtering and plotted cutoff. If the reconstruction is the
intended final analysis, the manuscript/source data need an explicit author
decision to align with it. Do not tune the analysis merely to force the counts.
Evidence for the normalization denominator alone does not prove the DE testing
universe. Raster comparisons do not resolve this provenance question.

The previous notes also mixed a full-table fold-change range with the plotted,
significance-filtered range for 5E, and gave two measured positive extremes for
5F (+13.53 and +13.13). The distinction is corrected in the documentation; the
unverified measurement is not used to change the analysis. ORF F is not plotted
in the current 5F because its adjusted p-value is 0.23.

## Supplementary 4: assembly decisions remain for this contribution

The original R scripts and audited source values do exist in Kourosh's project.
The shared data now reproduce the original 50 genes, independent CCC/MCC row
orders, phase means and sample-SD z-scores. Original labels and the -2.5/0/2.5
color breaks are retained. The CCC phase cell counts are 11/230/106/56/44;
MCC counts are 1517/2789/1219/271/262.

Pink description backgrounds were a manual assembly annotation and are now
copied from page 4 of the final supplementary PDF. It contains 19 highlighted
CCC rows and 20 MCC rows, with 19 shared identities. TGME49_315760 (AP2XI-4) is
highlighted only in MCC. This supersedes the earlier blanket "20 highlighted
genes" wording. The source caption states 40%, but CCC visibly has 19/50.

The code preserves that asymmetry exactly. Kourosh/Robyn: the only remaining
annotation question is whether AP2XI-4 should also be highlighted in CCC;
changing that would be a manuscript correction, not faithful copying. No new
expression-based rule is inferred. Original manuscript aliases, including the
MYND-like label, are preserved for reproduction, not asserted as new annotation.

## Other retained author checks

- **Figure 1C bands:** the earlier reconstruction read these as CWPs 5,
  Micronemes-1 7, Micronemes-2 6, ~Tachyzoites 6, Rhoptry-1 5, Rhoptry-2 5.
  Confirm whether `TGME49_306270` is the last Micronemes-1 row.
- **Enolase:** the existing Figure 1E mapping is retained. The original review
  resolved the suspected swap as an unused notebook-cell issue; no panel fix
  was needed. Historical files have not been deleted or renamed.
- **Coverage:** external Benke/ToxoDB time courses are Supplementary 2A-F, not
  Supplementary 1B-F. Supplementary 1B/C expression maps and computational
  Figure 3 panels are not yet in this repository's registered panel list.
  Assign these explicitly if they are required for the final release.
- **Data/code availability:** all 11 packaged inputs were retrieved through
  the commit-pinned public GitHub LFS mirror and verified by SHA-256 on
  11 September 2026. The two institutional mirrors timed out in this test.
  Final Zenodo DOI, raw-read accession and license still require the authors;
  do not describe planned deposition as complete.

Panel letters, composite layout and other hand-added annotations remain
manuscript-assembly work. These checks do not imply final manuscript approval
or completed archival deposition.
