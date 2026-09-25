# A10 public data search: is an independent cohort available

25 September 2026. Searched because stage 1 left A10's biological unit unresolved,
which caps the screen at within-screen description. An independent cohort would be
the other route to a defensible result. Sources: PubMed and the GEO DataSets index.

## Conclusion

**No comparable public dataset exists.** Nothing found pairs well-level RNA with
measured organoid growth under perturbation. External validation of A10 is therefore
not available from public data today, and the internal preparation structure is the
only validation route. That raises, rather than lowers, the importance of the
unresolved unit question.

## What was searched

| Query | Result |
|---|---|
| Alveolosphere or alveolar organoid with a CRISPR or knockout screen, in lung | 2 articles, neither relevant |
| Alveolar type 2 organoid, fibroblast co-culture, imaging, growth, transcriptome | 0 articles |
| GEO DataSets: alveolar organoid or alveolosphere, mouse | 932 series, top 30 inspected |

The two article hits were a human stem-cell-derived lung organoid stress reporter for
toxicity screening ([10.1016/j.mtbio.2026.102972](https://doi.org/10.1016/j.mtbio.2026.102972))
and a deubiquitinase screen in human alveolar organoids
([10.7150/thno.105994](https://doi.org/10.7150/thno.105994)). Both are human
stem-cell-derived systems without a fibroblast niche partner and without well-level
RNA linked to growth, so neither can validate A10. Article details came from PubMed.

## Nearest GEO candidates, and why each falls short

| Series | Samples | Why it is not a replication cohort |
|---|--:|---|
| GSE271971 | 72 | Alveolar epithelial growth supported by microvascular endothelium. The niche partner is endothelial, not fibroblast, so the dialogue being modelled differs. The nearest candidate of the set. |
| GSE213975 | 14 | A type 2 organoid platform for lung adenocarcinoma subtypes. Too few samples, and the outcome is tumour modelling. |
| GSE162859 | 8 | An alveolar organoid model on arrays. Too few samples, no perturbation series. |
| GSE313671 to GSE313673 | 2 to 6 | Alveolar regeneration mechanism studies, not organoid growth screens. |
| GSE221402, GSE223664, GSE183423 | 1 to 6 | Fibroblast biology in lung; no organoid growth endpoint. |

The three top GEO hits were this screen's own series and its superseries.

## What this means for the plan

1. **The unit question is now the binding constraint.** With no external cohort, the
   only available validation is holding out whole preparations inside this screen, and
   whether replicate indices are independent preparations is exactly what the deposit
   does not say.
2. **A descriptive result is not automatically worthless.** It would still be the
   register's first analysis anchored to a measured growth endpoint, provided it is
   labelled as within-screen association with no independent replication.
3. **One candidate is worth a later look, for a different purpose.** GSE271971 pairs
   alveolar epithelial growth with a supporting niche cell type. It could not validate
   A10, but if A10 finds a programme associated with growth, asking whether the same
   programme appears when the supporting cell is endothelial rather than fibroblast is
   a distinct and interesting question. That is a proposal, not a plan.

## What was not done

No dataset was downloaded from this search. No author was contacted. The source
paper for this screen remains unread by the owner and was not consulted.
