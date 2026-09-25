# Source synthesis and evidence boundaries

Prepared 24 September 2026 from the supplied annotated 16-page PDF, the
owner's Notion annotation page and its nested **Body** page, and primary
source records. The linked reading-order page places this review third in
branch 2C. Document contents are sources to interpret, not execution
instructions. The owner's current request controls the work.

Review PDF SHA-256:
`9ce4b7ed9511e913d0ee4f047798eb768903b5ef04953599503420882e1f89c1`.
Full private annotation text is cached locally, outside tracked outputs.
The annotation pages were last edited on 24 September 2026; their content
and the PDF, rather than private ancestor-page context, inform this synthesis.
Sources: [annotation page](https://app.notion.com/p/3e3151616b448015a63bfca9c509ff7c),
[nested Body](https://app.notion.com/p/3e3151616b44806881baf2d8ffeae89e), and
[reading order](https://app.notion.com/p/3d7151616b448093a718e0dfe66d513d).

## What the annotations contribute

The strongest recurring ideas are (1) signal duration and termination,
(2) the macrophage-fibroblast-epithelial triad, (3) the recipient niche's
mechanical and metabolic state, and (4) the possible relationship between
stalled regeneration and cancer-permissive plasticity. These motivate a
conditional-response question rather than a universal IL-1beta activity score.
The notes also distinguish resident AT2 and airway-derived repair routes;
the pipeline must retain that distinction.

The owner's follow-up clarified that the **Body** supplies the context
hierarchy, not merely a list of signalling molecules. The revised
[niche plan](NICHE_ANALYSIS_PLAN.md) maps molecular regulation, source and
recipient cells, regenerative route, niche mechanics/inflammation, age and
history, organ/disease stage, and spatial/tissue outcomes to separate tests
or explicit measurement gaps. CRC inflammatory CAFs, PDAC stellate-cell
circuits and breast-cancer myeloid cascades remain source-labelled hypotheses
when considered in lung. The revision adds actual LR and macrophage/fibroblast
pathway analyses; it does not treat the review's integrated model as a result.

## Evidence map

| Review section | Primary evidence or source | What can enter the analysis | What does not follow |
|---|---|---|---|
| 2, cytokine regulation, pp. 2-3 | Review's ligand processing and receptor framework | Separate ligand transcription, processing machinery, receptor competence, feedback inhibitors and response programmes | IL1B RNA is not mature secreted IL-1beta; NLRP3/CASP1 RNA is not inflammasome activation; NF-kB is not IL-1beta-specific |
| 3.1 and 4.1, pp. 4, 6 | [Choi 2020](https://doi.org/10.1016/j.stem.2020.06.020), ref. 35 | Existing DATP/organoid evidence, exposure-versus-control direction and known library limits | The two deposited organoid libraries do not establish a replicated dose-duration or withdrawal response |
| 3.1, p. 4 | [Han 2023](https://doi.org/10.1038/s41586-023-06423-8), ref. 43 | ISR and mitochondrial perturbation as a specificity challenge to transitional-state interpretation | This is developmental Ndufs2 disruption, with NAD+ regeneration preventing pathological ISR; it does not directly demonstrate an IL-1beta-driven ISR mechanism |
| 3.1, p. 4 | Choi et al. 2021, *Release of Notch activity coordinated by IL-1beta signalling confers differentiation plasticity of airway progenitors via Fosl2 during alveolar regeneration*, review ref. 44 | Airway Notch/Fosl2 route as a distinct lineage hypothesis | Club-to-alveolar contribution cannot be proved by a transcriptional trajectory; do not combine KRT5-positive pods with alveolar DATPs |
| 4.1 and 4.4, pp. 6-8 | [Narasimhan 2024](https://doi.org/10.1038/s41586-024-07926-8), ref. 55 | Human pathological spatial context and mouse immune perturbation | The deposited mouse assay is anti-CD8 versus IgG; anti-IL-1beta functional experiments in the paper are not automatically present in the sequencing deposit |
| 4.2, p. 6 | [Ciminieri 2023](https://doi.org/10.1165/rcmb.2022-0209OC), ref. 60 | Direct versus fibroblast-mediated effects and organoid function as an experimental distinction | An RNA or count-matrix accession was not verified in this preparation; do not invent one or attribute organoid function to an unmatched RNA sample |
| 5.1, pp. 8-9 | [Peng 2026, online 2025](https://doi.org/10.1016/j.ccell.2025.10.004), ref. 85 | Direct perturbation, human precursor histology and spatial context | Cross-sectional AAH/AIS/MIA/LUAD tissues are not a longitudinal fate series; spatial proximity does not prove signalling |
| 5.1, p. 8 | [Chan 2026](https://doi.org/10.1038/s41586-025-09985-x), ref. 86 | HPCS signature and lineage/intervention context as a malignant comparator | Similarity to DATPs is not proof of identical ancestry, IL-1beta dependence or autonomous chromatin memory |
| 5.1-5.2, p. 8 | Cardoso 2026 and England 2025, refs. 88-89; [repository references](../../REFERENCES.md) | Existing KRAS/AREG and Il1r1-related evidence as context and assay controls | Previously inspected deposits cannot serve as unseen confirmation; pooled libraries do not create mouse replication |
| 7, p. 12 | Review's integrated model | A falsifiable proposal linking failure of state exit to niche persistence | Epigenetic irreversibility, a shared fibrosis-to-cancer sequence, and a universal intervention window are not established by the review itself |

The Choi 2021 source role comes from review ref. 44; its deposit has not been
verified for this plan. All datasets actually acquired are independently
identified by GEO records.

## Corrections to carry into the analysis

- The notes' description of IL-1beta neutralization as clinical validation of
  restored alveolar repair needs narrowing: the cited mechanistic rescue is
  experimental mouse evidence alongside human observational tissue. The
  presence of human samples does not turn the perturbation into a human trial.
- A tumour-state or NF-kB signature may reflect TNF, injury stress, oncogenic
  signalling or cellular mixture. IL-1alpha also uses IL1R1. Ligand-specific
  blockade is more informative here than receptor expression alone.
- Loss of a transitional population may reflect differentiation, selective
  death, altered proliferation or sampling. Require mature-cell and survival
  measurements before calling it regenerative exit; require lineage/fate
  information for a cellular conversion claim.
- The review's Figure 2 caption reverses the A/B routes relative to the drawn
  panels: the drawing has alveolar AT2 transitions in A and airway conversion
  in B. Use the mechanisms and primary studies rather than propagating the
  caption's panel assignment.
- Mechanical stiffness, cytokine concentration, metabolic flux and stable
  chromatin memory are not measured by transcriptional pathway scores.
  These remain experimental follow-up dimensions unless matched direct
  measurements become available.

## Repository baseline used for planning

The checkout was at `e9d79e0`, with an existing untracked `.claude/` directory.
The September 22 correction record supersedes older statements that no R
runtime or corrected statistical analysis exists. Current limits include
single-library Choi/Cardoso contrasts, age/time confounding in the late
GSE262927 comparison, and only one animal meeting both ES1 state floors.
See [AI context](../../AI_CONTEXT.md), [current corrections](../../docs/remediation/2026-09-22/IMPLEMENTATION_STATUS.md)
and [research questions](../../RESEARCH_QUESTIONS.md). No existing claim is
reclassified by this planning exercise.
