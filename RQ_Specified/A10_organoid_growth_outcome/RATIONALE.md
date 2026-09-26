# A10 rationale: connect compartment RNA to measured organoid size

Originally declared 25 September 2026; revised 26 September after the completed
fits and logical review. This is the current biological argument, not a
retrospective preregistration. The [original plan](PLAN.md), frozen specifications
and numerical outputs preserve the execution record. Read the
[stage 4 report](reports/STAGE4_REVISED_REPORT.md) and
[logical review](../../docs/LOGICAL_RATIONALE_REVIEW.md#a10-corrections).

**Execution update:** the diagnostic-first follow-up is now complete. The
[results](reports/FOLLOWUP_RESULTS.md) show conditional information beyond the
two nominated proliferation scores and relative gains in all plate holdouts,
alongside negative absolute R-squared on three of four plates. Target allocation
and plate cannot be separated adequately; independent preparations remain unknown.
The sequence below explains the logic, with the completed follow-up supplying
steps 5/6. Further model expansion was pruned pending design/calibration evidence.

## Biological premise

Alveolar type 2 cells produce surfactant and include progenitors that can replenish
the alveolar epithelium. Type 1 cells form its thin gas-exchange surface. Restoring
these functions requires more than an increase in a transitional RNA signature.

Niche signals help regulate AT2 stemness and differentiation. Nabhan et al. found
fibroblast Wnt niches supporting an AT2 stem-cell subset, alongside injury-induced
epithelial autocrine Wnt recruitment of additional progenitors. Thus the source
and role of Wnt depend on context
([Nabhan 2018](https://pubmed.ncbi.nlm.nih.gov/29420258/)). This motivates measuring
both compartments. It does not guarantee that fibroblast RNA adds information
after epithelial RNA, or that every AT2 organoid culture requires fibroblasts.

## What this design measures

The deposit cultures mouse AT2 cells with human fibroblasts, with perturbations
introduced in the mouse epithelium. Species assignment identifies the origin of
RNA. Both endpoint profiles may nevertheless reflect feedback, cell abundance,
selection and shared culture conditions. Epithelial RNA is not a pure
cell-autonomous readout; fibroblast RNA does not isolate a causal mediator.

Imaging supplies organoid count, mean area and area proportion at days 7 and 14.
The fixed primary endpoint is **day-14 mean area conditional on day-7 mean area**.
It measures average size and morphology, not total tissue yield, mature AT1
function or repair in a living lung. Area can change without proportional changes
in cell number. A10 complements A1's existing lineage and microscopy outcomes;
it is not the only repository question with a non-RNA endpoint.

## Question and logical sequence

The bounded question is whether epithelial scores add information about the
recorded area endpoint after day-7 size and RNA composition proxies, and whether
fibroblast scores add more. These are conditional associations. Day-7 size is
already post-perturbation, while day-14 RNA and imaging are concurrent; controlling
for them does not identify the perturbation's total or mediated effect. Species
read fractions are not direct measurements of cell fractions.

1. **Join the measurements.** Link design, imaging and species-assigned RNA by
   well. The deposit has 886 RNA libraries; 885 enter the completed area analysis.
2. **Identify independent units.** Fifteen deposited plate-replicate groups span
   four plates. Those labels do not establish separate epithelial isolations or
   fibroblast donors/lots. The result remains descriptive while this is unresolved.
3. **Specify the comparison.** Keep the original repair-module model and its
   result. The later growth-module grid changes both programme set and centring.
   It was declared before its own fit, after the original fit had been inspected.
   It is a same-screen adaptation, not independent confirmation.
4. **Read the metric literally.** Within-unit centring uses each held-out group's
   own outcome mean. It evaluates variation around an observed group mean; it
   cannot forecast an unseen group's size. Changing centring changes the
   R-squared denominator. The same 0.02 margin is not automatically equivalent
   in biological importance or more conservative across those scales.
5. **Check heterogeneity before expanding the model.** The revised epithelial
   increment is 0.0426, but positive in only 8/15 groups and concentrated in two
   plates. The fibroblast increment is -0.0154. A negative conditional increment
   does not exclude niche dependence: epithelial RNA may already reflect niche
   action, or the nominated fibroblast scores may miss it.
6. **Separate explanatory description from a stronger test.** Inspect plate and
   target allocation, then predeclare a simpler growth-block comparison and
   plate-level holdouts. A metric that uses held-out outcome means must retain
   its within-group label. Four plates do not establish independent preparations.
7. **Require different evidence for mechanism.** Compartment-specific interventions,
   verified perturbation efficiency and independent biological replication are
   needed to distinguish direct effects, feedback and mediation. More scoring
   cannot supply those design features.

Target-transcript reduction in 58/76 tested genes is a limited consistency check,
not complete assay validation. The comparisons pool other targets and do not
establish editing efficiency in every well or recover a prespecified imaging
positive control. The review records this and the secondary-model implementation
limits; all original numerical outputs remain intact.

## Connections to A5, A11 and niche questions

The original source-list partition overlaps the screen's perturbation targets:

| Original module | Genes | Perturbed members |
|---|--:|---|
| Shared development and injury | 5 | Areg, Itgb6 |
| Development-list-exclusive | 94 | 11, including Egfr, Erbb3, Fzd5, Tfcp2l1 |
| Lesion-list-exclusive | 91 | Tgif1, Tigit |
| Shared injury and lesion | 7 | none |

These are source-list overlaps, not functional validation. A5's revised primary
uses additional external Guo modules; the old Strunz-filtered variants remain
descriptive. A marker may be downstream, redundant or active in another context,
so membership does not imply that its knockout changes organoid area. Conversely,
an area effect need not validate developmental reuse or lesion specificity.

Epithelial Areg or receptor perturbations cannot by themselves resolve A2's
source contributions or A9's fibroblast-specific receptor competence. Those
questions require their own compartment-specific interventions and endpoints.
Likewise A10 does not fill A1's regulation-to-fate gap, A5's lineage/function gap
or A11's missing non-neoplastic injury comparator.

## Source-paper boundary

The follow-up includes a bounded post-fit design check: the main article identifies
TIGIT as its in-plate control, while the detailed methods supplement was not
successfully retrieved. The [source record](tables/followup_v1/source_design_check.json)
states what was and was not inspected. This does not create a study note, update
the owner's reading status or establish preparation independence.

The source paper is roadmap paper 14 (DOI 10.1073/pnas.2606113123). Its owner
reading remains recorded as pending. This review uses the deposit and existing
analysis records, not a new study note or acceptance of that paper's claims
(DEVELOPMENT decision 21). Author requests remain unsent. Preparation identity
is a missing fact, not a scientific choice that can be resolved by relabelling.

## Current conclusion

Under the revised specification, epithelial growth scores add information about
relative mean organoid area inside observed groups, with substantial plate
heterogeneity. Fibroblast growth scores do not improve that model. The descriptive
margin rule is met; independent replication, precise absence of a niche role,
cell fate and an in vivo repair mechanism remain unestablished.
