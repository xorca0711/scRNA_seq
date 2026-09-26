# A0: recovery of an executable pilot

26 September 2026. This is an eligibility and source-identity report. The
[fixed specification](../config/pilot_v1.json) precedes discovery. Final scientific
results are reported separately; coverage is not a conservation result.

The initial A0 audit was local and untracked. Commit `4ed38ff` preserves it.
The continuation checked the original plan, A1's regulatory/fate evidence, A5's
completed external-signature test, A5/A11's pairwise component contract and A10's
outcome follow-up. None already supplies A0's three-context test.

## Biological sequence

1. A transition can share expression with another transition because of common
   differentiation work, generic stress, proliferation, or changing identity.
2. Compare an independently annotated intermediate with **both** its earlier and
   later reference states within each biological unit. A change against the
   starting state alone could simply measure acquisition of destination identity.
3. Discover genes meeting the same directional criteria in repair and normal
   development. Distinct lung branches test a broader proposition than a second
   analysis of alveolar injury alone.
4. Freeze the resulting module before looking at its effects in another tissue.
5. Challenge a transferred signal with cycling, stress, depth and matched-gene
   controls. Even success establishes a state association, not causal fate control.

The broad hypothesis is a conserved process that might modulate fate. This pilot
tests only one operational prediction: a shared, intermediate-enriched RNA module.
Spatial position, RNA trajectories and named states do not trace individual cells.

## Candidate decisions

| Candidate | New check | Decision |
|---|---|---|
| Original Negretti developmental atlas | Separate endpoint comparisons permitted by the original plan | Counting pairwise comparisons cannot supply three intermediate-versus-AT1 units; animal mapping remains missing |
| Guo P1, reused as a list by A5 | [Source Figure 1](https://www.nature.com/articles/s41467-018-07770-1) has two mice | Useful list source, below A0's three-unit discovery floor |
| Zhang/Chuang Hippo 2026 | [GSE319370](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319370) has one control deposit at each of two ages, plus mutants | Two normal observations cannot provide three independent normal comparisons; do not count mutants as normal development |
| He 2022 fetal lung | Downloaded author-linked 9,096-cell epithelial object; read metadata only | Late-tip/late-stalk/AT1 coverage fails, including separate endpoint comparisons; no expression effects inspected |
| Quach 2024 fetal lung | [Author resource](https://github.com/The-HQQ/Human_fetal_lung_atlas) identifies a smaller epithelial object | Anonymous Synapse file download returned 403; GEO whole object is 7,981,428,892 bytes. Local free memory was about 1.8 GB at inspection. Eligibility remains unknown; this is not evidence of insufficient replication |
| Sountoulidis 2023 | [Paper-linked UCSC atlas](https://cells.ucsc.edu/?ds=lung-dev) supplies donor IDs, annotations and matched integer counts | Select the declared early airway axis; 11 complete donor triplets before expression-based discovery |
| Miller 2020 | A similarly named UCSC entry was initially reached; its description identifies a different study | Retain the retrieval trail, but do not misattribute it or select it as D2 |
| Original Haber immature-proximal branch | Four mapped mice; two complete triplets | Historical proposal remains ineligible |
| Revised Haber late-progenitor branch | Same four mapped mice; three complete triplets | Freeze this earlier commitment-stage comparison before any programme effects; do not claim the progenitors are uniquely proximal |

## Selected measurement scope

**D1:** Strunz day 10–15 injured mice, excluding source-labeled normal controls.
Compare Krt8+ ADI with AT2 and AT1. Keep activated AT2 as a biologically relevant
alternative comparison. Previous A5 expression work in these mice is disclosed;
the repair dataset is not untouched validation.

**D2:** Sountoulidis author clusters `epi_cl2`, `epi_cl1`, `epi_cl0`, within donor.
The [source](https://doi.org/10.1038/s41556-022-01064-x) combines state annotation
with spatial evidence and inferred developmental relationships. These populations
describe a distal-to-proximal airway axis, not alveolar AT2-to-AT1 maturation.
Use the matching UCSC donor annotations as the biological-unit definition;
captures do not add n. The source's spatial evidence supports an intermediate
position, not direct proof of each cell's ancestry or future fate.

**V1:** Haber source-labeled stem cells, late enterocyte progenitors and mature
proximal enterocytes. The [source](https://doi.org/10.1038/nature24489) places
enterocyte commitment within the crypt-to-villus differentiation system. This
comparison is neither individual-cell lineage tracing nor a claim that every
late progenitor takes the proximal branch. Only explicitly mapped control mice
enter the pilot; unmatched atlas batches stay excluded.

The human developmental role requires a conservative one-to-one mouse–human
mapping. Mapping and assay coverage are fixed before discovery; native counts
remain separate. A small paired pilot does not identify the effect of species,
tissue, platform or age separately. Annotation-gene exclusions reduce direct
reuse of marker panels but cannot make expression-derived labels independent.

## Retrieval and execution limits

- The small Synapse folder metadata were publicly readable; the actual epithelial
  file returned 403. No account restriction was bypassed and no author contacted.
- Europe PMC XML retrieval for Miller and Haber returned HTTP 500. Haber biology
  was checked against its primary PMC text; no inaccessible supplement is assumed
  to lack information.
- An initial guessed GSE215896 download was rejected by its unrelated title.
  GSE215898 is the paper's authoritative SuperSeries (scRNA child GSE215895).
- A 15 MB homology-download bound was slightly too small. The same official
  source was subsequently retrieved with a declared 100 MB bound (15.11 MB file).
- Public downloaded sources remain ignored. Their URLs and hashes are retained
  in the continuation source manifest; the original source manifest is historical.

No result-dependent replacement of a failed module is permitted. In particular,
the original requirement of at least 20 selected genes is retained. If that gate
fails, omit programme transfer/specificity scoring and report a stopped pilot,
not a negative test of every possible universal regulatory process.
