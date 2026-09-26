# Research questions

> Which epithelial programmes and immune–stromal interactions distinguish
> productive lung repair from persistent remodelling after injury?

This project uses public lung RNA, chromatin and complementary assays to develop
and test its own biological hypotheses. Our observations motivate developmental
programme reuse, context-dependent loss of epithelial identity, and recipient
and fibroblast contributions to injury responses. Those interpretations are
scientific propositions to distinguish from alternatives, not just criticisms
of markers or analysis tools.

**Rewritten 25 September 2026; A0 registered 26 September.** All A1–A14 identifiers remain; A12-S1 remains an
enabling source-identity question. There is no fixed number of “surviving” RQs.
Related questions share evidence but retain separate tests and decisions.

The [claim register](CLAIMS.md) grades the original measurements. A new biological
interpretation does not inherit a historical claim's status. Here **motivation**
describes the connection to observations; **readiness** describes whether a
discriminating test can run. This rewrite establishes no new mechanism and
changes no result, threshold or claim grade. These are source-informed and often
post-analysis hypotheses, not retrospective preregistrations or novelty claims.

The [measurement contracts](docs/RQ_MEASUREMENT_CONTRACTS.md) collect technical
checks that can change a decision. The [figure gallery](analysis/figures/rq/README.md)
preserves measured panels, full captions, diagnostics and labelled designs.
The [migration record](docs/migrations/2026-09-25-rq-reframing/README.md) preserves
the previous register. Paper-specific evidence stays under `Research Article/`;
question-specific execution stays under `RQ_Specified/`.

## Biological hypotheses and execution priorities

The collections do not share an independently measured repair outcome. Analyse
development, infection, fibrosis and neoplasia within their own designs; do not
order them along an assumed repair-to-cancer trajectory. Late sampling is not a
recovery endpoint, and population persistence does not trace the same cells.

| ID | Biological hypothesis / decisive endpoint | Motivation | Readiness and next task |
|---|---|---|---|
| [A0](#a0) | A conserved transition-associated programme may contribute to epithelial fate modulation | Related transitional RNA states across repair and development; shared fate control remains a hypothesis | [Pilot complete](RQ_Specified/A0_conserved_epithelial_transition_program/reports/PILOT_V1_RESULTS.md): 50-gene lung signature fails the mature intestinal endpoint; P4 pruned; universal/causal fate claims unresolved |
| [A1](#a1) | Regulatory features distinguish RNA-similar transitional states and responses | Direct marks, lineage, perturbation and paired CD44 context contrasts | Adaptive batch complete; matched replicated regulatory/fate linkage remains missing |
| [A2](#a2) | AREG changes a fibroblast response through delivery to a competent recipient rather than through source abundance | Seven register rows close the abundance version (C37, C39, C40, C45, C48, C49, C50); the mechanism places the rate-limiting step in the recipient | [Both legs run](RQ_Specified/A2_areg_source_delivery/reports/STAGE5_SYNTHESIS.md): the Areg arm does not move and its null is bounded by an unremoved autocrine source; epithelial Itgb6 moves by -0.938 log2 CPM in three of three readable units; the donor-level leg, once put on a common molecule budget, reads null at rho 0.293 |
| [A3](#a3) | Injury leaves a macrophage programme beyond normal aging | Late population composition | Age-matched controls and comparable sampling needed |
| [A4](#a4) | Wnt maintenance and IL-1 response occur sequentially in an AT2 lineage | Transcript/source observations; sequence untested | Measured activity/history and lineage-linked response needed |
| [A5](#a5) | Adult repair reuses a developmental epithelial component | Neonatal coexpression and label-excluded ADI enrichment; outside developmental list now sourced | Revised external-signature test positive in 24 mice and after identity/control exclusions; lineage/function untested |
| [A6](#a6) | IPF changes shared macrophage states beyond subtype abundance | Cell fractions and RNA contributions differ | Harmonize states and audit donors before a within-state fit |
| [A7](#a7) | Cebpa loss attenuates identity across AT2 states | Reference and transitional contrasts both change | Replicated genotype-by-state design needed |
| [A8](#a8) | A maturation component adds information about mature AT1 contribution | Score dependence motivates separation; limited biological support | Independent mature endpoints needed |
| [A9](#a9) | Fibroblast receptor context modifies AREG response | RNA/resource observations nominate a competence question | RNA screen possible; protein/function data needed |
| [A10](#a10) | Epithelial programmes add information about measured organoid growth | Public RNA and imaging design | Follow-up complete: added information beyond E2F/G2M and positive plate-shift error gains, but negative absolute R-squared on 3/4 plates; independent units unresolved |
| [A11](#a11) | Lesion-associated programmes add to shared plasticity | Reduced HPCS signal across repair, development, IPF and LUAD | Kim lesion association replicates in 8 patients; beyond-shared criterion unresolved (BH q=0.0547) |
| [A12](#a12) | Recipient context explains responses beyond ligand RNA | Cohort/recipient heterogeneity | Conditional component model; activation unmeasured |
| [A13](#a13) | Fibroblast programmes add information beyond macrophage IL1B | Niche heterogeneity motivates joint association | Audit complete triads; current IPF arms below joint-model floor |
| [A14](#a14) | Exposure duration and fibroblast reception separately affect recovery | Mechanistic follow-up to the repair/persistence question | Two decisions; withdrawal and recipient-specific data needed |

A0 completed its bounded pilot: the frozen programme fails the mature intestinal
endpoint, and its conditional specificity work is pruned. The A1 continuation,
A5/A11 revised tests and A10 outcome join/fits are complete.
The A10 plate/design diagnostics and separately specified growth-block test have
also run. Further A10 expansion now needs preparation identities and imaging/
validation design evidence; stronger A1, A5 and A11 conclusions require the
missing evidence named in their cards. Existing joins and scores need not be
rerun. The other questions remain conditional; their design gates still apply.
The [logical review](docs/LOGICAL_RATIONALE_REVIEW.md) checks this sequence against
the actual measurements and implementation, with an adaptive follow-up order.

## Hypothesis cards

<a id="a0"></a>

### A0. Is a conserved programme reused across epithelial transitions, and could it modulate fate?

**Hypothesis.** Different epithelia may reuse part of the cellular work needed to
leave an established identity and acquire another. A shared RNA component is one
possible observable consequence. Whether that component controls maturation,
persistence or reversibility is a separate causal question.

**Biological logic.** First establish the starting, intermediate and destination
populations from each study's biological context. Compare the intermediate with
both endpoints within mice or donors: a difference from the starting population
alone could measure ordinary acquisition of destination identity. Seek one shared
intermediate-enriched module in injury and normal development, then freeze it
before testing another epithelium. Challenge any transferable association with
stress, proliferation and measurement controls. Only matched perturbation and
fate endpoints could subsequently support modulation of fate.

**Current measurement.** The [A0 workspace](RQ_Specified/A0_conserved_epithelial_transition_program/README.md)
contains the original feasibility audit and the separately specified scientific
pilot. D1 is Strunz mouse alveolar repair; D2 is Sountoulidis human early airway
development; V1 is Haber mouse intestinal enterocyte commitment. These are
biologically distinct axes, not a single pooled trajectory. Source spatial and
developmental evidence supports the comparison but does not trace the future of
each sequenced cell. Read the [executed result and stopping decision](RQ_Specified/A0_conserved_epithelial_transition_program/reports/PILOT_V1_RESULTS.md)
before proposing any further score or module. Of 134 qualifying genes, 50 were
frozen. In three intestinal mice, the intermediate-minus-mature median is −0.0914
score points with only one positive mouse; the primary transfer rule fails.
Specificity work is pruned under the rule declared before V1 scores.

**Relation to established work.** A5's external developmental signature supports
partial recruitment in adult repair, not a universal transition machinery.
A5/A11's shared contract contains pairwise overlaps rather than one three-context
core. A1 addresses regulatory and response distinctions; A10 measures organoid
growth. Those results motivate A0 without supplying its cross-tissue transfer or
causal fate test. The [source recovery report](RQ_Specified/A0_conserved_epithelial_transition_program/reports/SOURCE_RECOVERY.md)
records the analysis precedents, cohort decisions and prior exposure.

**Scope of the decision.** The pilot requires 20–50 qualifying one-to-one orthologs
under fixed four-contrast criteria. Failing that operational prediction cannot
exclude a smaller shared component, conserved regulation with different RNA
outputs, or all possible transition processes. A positive transfer would remain
a state association. Neither outcome by itself resolves universality or fate.

<a id="a1-which-regulatory-and-phenotypic-features-distinguish-transitional-epithelial-states-beyond-rna-markers"></a>
<a id="a1"></a>

### A1. Do regulatory programmes distinguish RNA-similar transitional states and their functional responses?

**Hypothesis.** Transitional epithelia share part of an RNA response but differ
in regulatory programmes that help explain maturation, persistence or perturbation
response. DATP, PATS, Krt8 ADI, ABI/aberrant basaloid and HPCS remain source-defined
states; their names predetermine neither equivalence nor separate cell types.

**Our observation.** AT2 RNA loss survives the specified depth controls in two
same-laboratory deposits (C118). Distal-accessibility comparisons do not establish
closure or its timing (C131/C133). The completed
[first batch](RQ_Specified/A1_transitional_epithelial_state_distinction/reports/FIRST_BATCH_REPORT.md)
adds measured PATS endpoints, a separate ten-mouse IRE1α RiboTag contrast and
descriptive ATAC/CD44 profiles. Four of 14,811 genes pass whole-family FDR;
none of the predefined markers or eligible pathways does. The verified
[second batch](RQ_Specified/A1_transitional_epithelial_state_distinction/reports/SECOND_BATCH_REPORT.md)
adds directionally stable but significance-sensitive IRE1α effects, direct
histone profiles at 23 loci in two induced-cell preparations, and one-donor
methylation-domain context. Histone directions can depend on H3/window choice.
Recovered HPCS trace-linked RNA composition covers 5,333 cells / 22 sources;
the latest primary animal table resolves the mouse identities, while current
mScarlet remains unavailable.
The [robustness batch](RQ_Specified/A1_transitional_epithelial_state_distinction/reports/ROBUSTNESS_REPORT.md)
shows promoter-dependent CDKN1A acetylation, source-sensitive HPCS fractions and
chase/library aliasing that prevents the proposed fixed-library-adjusted temporal
contrast. The [closure batch](RQ_Specified/A1_transitional_epithelial_state_distinction/reports/EVIDENCE_CLOSURE_REPORT.md)
recovers CD44's exact eight-mouse crosswalk and fits paired genotype contrasts
plus their direct interaction. Seven transported markers change in both
genotypes; four have detected effect-size interactions, although Sftpc's effect
nearly vanishes when WT2 is omitted. HPCS's author biological
map is now verified, and all stringent K12 changes are confidence abstentions.
The [regulatory/outcome continuation](RQ_Specified/A1_transitional_epithelial_state_distinction/reports/REGULATORY_FATE_REPORT.md)
resolves HPCS harvest timing and Tsutsui perturbation-library identities, and
reanalyzes AP-1 microscopy with fields nested within mice. HOPX responses have
opposite directions by lung region; culture program suppression and
differentiation capacity remain separate experiments. TP53 selected RNA lists
include shared and opposing responses across AT1/AT2 origins. These results
make origin and local context necessary parts of the proposed comparison.
The [reference map](RQ_Specified/A1_transitional_epithelial_state_distinction/reports/ANALYSIS_REFERENCE_MAP.md)
connects published precedents, dependent reanalyses and the next decision.
These within-context measurements do not establish an epigenetic taxonomy.

**Rivals and test.** A common regulatory continuum, distinct branches, and
different routes sharing stress RNA remain alternatives. Compare frozen
loci/programmes in compatible direct histone-mark or accessibility assays using
independent animals/donors and independent state definitions. Where observations
can be linked, test added information beyond RNA about measured protein,
lineage-descendant or perturbation endpoints. Separate studies triangulate;
they do not constitute a paired multiomic fate test. ATAC does not measure
histone marks, methylation or chromosome conformation.

**Decision / readiness.** A precise independent exclusion of the predefined
meaningful regulatory effect or added endpoint information retires that specified
discriminator. Incompatible callers, pooled identities or low replication leave
it inconclusive. The [plan](RQ_Specified/A1_transitional_epithelial_state_distinction/PLAN.md),
[assay map](RQ_Specified/A1_transitional_epithelial_state_distinction/STUDY_MAP.md)
and [lineage audit](RQ_Specified/A1_transitional_epithelial_state_distinction/LINEAGE_AUDIT.md)
govern the input gate. PATS track scaling, TIGIT pool membership and HPCS
source/age reconciliation remain specific holds; CD44 count identities are
resolved. The newer IRE1 cell-resolved cohort has one pooled library per
condition and cannot supply replicated treatment inference. Temporal closure or memory additionally needs actual
time/fate evidence (A14).

**Figures / contracts.** [Existing assay and endpoint panels](RQ_Specified/A1_transitional_epithelial_state_distinction/figures/README.md);
Descriptive direct-mark and descendant-source panels are available; independently
replicated regulatory/endpoint tests remain future work. [Depth panels](analysis/figures/rq/README.md#a1) are
supporting diagnostics. [MC1–MC2](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1).

<a id="a2-which-cells-express-areg-and-how-sensitive-are-candidate-rankings-to-the-resource"></a>
<a id="a2-does-the-functional-contribution-of-areg-sources-depend-on-tissue-context"></a>
<a id="a2"></a>

### A2. Does the fibroblast response to AREG depend on delivery or on abundance?

**Hypothesis.** AREG's contribution to a fibroblast response is set by where the
ligand is released relative to a competent recipient, rather than by how much of it a
compartment transcribes. The mechanism is short-range and recipient-licensed, so the
quantity that matters is the recipient's state and the ligand's point of release.

**Why the abundance question is closed.** This card previously asked which source
dominates. Seven register rows answer that question and none establishes a
depth-independent epithelial hierarchy: the epithelial over myeloid contrast is
sensitive to annotation and molecule matching (C37), an epithelial source in
adenocarcinoma is not established because dendritic cells sit above both tumour states
(C39), those cells are a major source in their own right, a ranking in the opposite
direction (C45), tumour enrichment is refuted and fibrosis enrichment is not
established (C40, C48), and no donor-level correlation between epithelial AREG and
fibroblast EGFR or fibroblast activation was established (C49, C50). C51 records why a
further correlation on those variables would be unreadable: the one significant pair
in that trial tracked sequencing depth, and a frozen rule refused it.

**Mechanism.** Amphiregulin acts on a mesenchymal recipient through that recipient's
own EGFR, and the signal then activates integrin alphaV to release bioactive TGF-beta
from latent complexes, driving myofibroblast differentiation. The recipient in that
work is a PDGFRB-positive pericyte and macrophages are a critical source
([Minutti 2019](https://doi.org/10.1016/j.immuni.2019.01.008)). Silencing
amphiregulin or inhibiting EGFR reduces TGF-beta1-driven fibroblast proliferation,
smooth-muscle actin and collagen, and the amphiregulin silenced there is the
fibroblast's own ([Zhou 2012](https://doi.org/10.1074/jbc.M112.356824)). A
short-range ligand that converts a store the recipient already holds does not require
tissue-level ligand to be rate-limiting, so C49 and C50 do not refute the axis;
neither row is evidence for this framing either, and C49's point estimate is positive
and underpowered. Leukocyte amphiregulin can be non-redundant for lung protection
([Arpaia 2015](https://doi.org/10.1016/j.cell.2015.08.021)), and a review establishes
both epithelial and leukocyte sources
([Zaiss 2015](https://doi.org/10.1016/j.immuni.2015.01.020)), which is why C45
constrains an epithelium-only reading rather than refuting the axis.

**Rivals and test.** [Analysis A2](RQ_Specified/A2_areg_source_delivery/README.md)
holds the plan, the contract and the freeze that governs what may be computed.
Rivals: the recipient supplies the same ligand itself; another EGFR ligand carries the
response; the epithelial state rather than its ligand changes the fibroblast; the
epithelium activates TGF-beta through its own integrin; the fibroblast profile moves
with well composition or read depth. The organoid screen perturbs the mouse epithelium
only and leaves the human fibroblasts unedited, so it can compare an epithelial source
contribution against loss of epithelial reception and against loss of epithelial
integrin-mediated TGF-beta activation. It cannot separate delivery from abundance: a
single well holds one source compartment with no spatial variation, and the audit found
these fibroblasts transcribe AREG at a higher within-compartment level than the
epithelium, so a positive result is equally consistent with the abundance version and
bounds an increment rather than establishing necessity. That contrast needs the spatial
layer.

**Decision / readiness.** Conditional and descriptive; both legs have now run and
neither settles the question. In the organoid screen the frozen fibroblast activation
score does not move on epithelial Areg knockout, a median of +0.036 log2 CPM against
depth-matched controls with two of four units in the predicted direction and an endpoint
standard deviation of 0.482. That null is bounded rather than absent: the recipient
transcribes AREG at 9.846 against the 7.267 removed, the knockout is partial, and the
culture medium is not in the deposit. The arm that moved is epithelial **Itgb6**, at a
median of -0.938 in three of three readable units, with organoid size and fibroblast
content unchanged, which proposes epithelial integrin-mediated TGF-beta activation rather
than the ligand and needs its own design. The donor-level leg was refused by the frozen
C51 depth rule at a composite depth coupling of 0.770, and has since been read: on a common
molecule budget the coupling falls to 0.232 and the correlation falls from 0.433 at nominal
p 0.044 to 0.293 at p 0.186, so the rule was protecting against depth and no coupling
survives it. The first freeze, which
declared a rank test with an exact null, is
[withdrawn](RQ_Specified/A2_areg_source_delivery/reports/STAGE2_WITHDRAWN.md). Read the
[synthesis](RQ_Specified/A2_areg_source_delivery/reports/STAGE5_SYNTHESIS.md). Secreted
ligand, receptor engagement (C36) and proximity remain outside the repository, and the
decisive delivery-versus-abundance contrast needs the spatial layer.
[Current figures](analysis/figures/rq/README.md#a2) are diagnostics.
[MC2 to MC4](docs/RQ_MEASUREMENT_CONTRACTS.md#mc2).

<a id="a3-which-macrophage-programmes-vary-with-phase-and-what-explains-the-differences"></a>
<a id="a3"></a>

### A3. Does prior injury leave a macrophage programme that differs from normal aging?

**Hypothesis.** Prior infection is associated with a late macrophage programme
or state distribution beyond changes expected with age alone.

**Our observation.** Late myeloid states and altered capillary composition
motivate a lasting tissue-response question. G1/W1 retain age and processing
confounding (C158). Reference CAMERA yields no significant W1 sets; its
seven-gene ornithine set is ineligible, not negative. See the
[sample-level analysis](analysis/corrections/statistics/README.md).

**Rivals and test.** Normal aging, recruitment/replacement, subtype mixture and
processing compete with an injury-history effect. Compare age-matched uninjured
and previously injured animals using harmonized states and an estimable
injury-history contrast. Report state fractions and within-state programmes.
Same-cell persistence requires tracing; metabolic flux requires a metabolic
endpoint. Neither follows from late RNA or relative fractions.

**Decision / readiness.** A matched design precisely excluding a predefined
meaningful effect retires the specified injury-associated programme. The current
confounded contrast cannot decide it. Age-matched data are the gate.
[Existing temporal panels](analysis/figures/rq/README.md#a3) describe sampling
after infection; desired primary plots compare animal-level effects with matched
controls. [MC1](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1),
[MC5](docs/RQ_MEASUREMENT_CONTRACTS.md#mc5).

<a id="a4-how-do-current-wnt-activity-and-il-1-responsiveness-overlap-in-at2-cells"></a>
<a id="a4"></a>

### A4. Can Wnt-supported maintenance precede an IL-1-responsive transition in the same AT2 lineage?

**Hypothesis.** Wnt-associated maintenance and IL-1-associated transition can be
sequential states of a lineage rather than fixed opposing subsets.

**Our observation.** Transcript co-detection and Wnt/source profiles provide
feasibility information, not a measured sequence. [Nb1](Research%20Article/gate1_03_nabhan_2018/nb1/README.md)
starts at day 6 and has only one eligible baseline and one day-11 AT2 unit;
it cannot test an acute switch. Current activity, reporter history and sparse
Axin2/Il1r1 transcripts are different measurements.

**Rivals and test.** Stable subsets, concurrent signalling and selection of
different cells remain alternatives. Use a Wnt-history pulse/chase with reporter
washout, present pathway activity, IL-1 challenge and traced descendant/function
endpoints across independent animals. Separate fibroblast ligand sources from
epithelial receiver activity.

**Decision / readiness.** A precise absence of the predicted transition among
verified history-labelled responsive cells weakens the sequential model.
Sparse RNA overlap or an unverified reporter is inconclusive. Suitable
history/activity data are required. [Current figures](analysis/figures/rq/README.md#a4)
are screens; the desired figure is a lineage-linked activity/response time course.
[MC1–MC3](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1).

<a id="a5-which-transitional-signatures-are-specific-to-injury-rather-than-development-or-genotype"></a>
<a id="a5"></a>

### A5. Does adult alveolar repair reuse part of a developmental epithelial programme?

**Hypothesis.** Development and adult repair recruit a shared epithelial
remodelling component, with context-specific additions contributing to different
outcomes.

**Our observation.** Neonatal controls contain Krt8/Cldn4 co-detection at the
common depth budget (3.69%; C119). ES1's label-excluded ADI enrichment is
+1.04/+1.21 detection points in neonatal controls and +7.61/+6.53 in injured
adult controls across technical seeds. These positive within-well observations
motivate reuse without equating neonatal and adult injury states. Only one of
25 external-study animals passes both group floors. See
[ES1](Research%20Article/epithelial_state_specificity/README.md).

The [shared component contract](RQ_Specified/A5_A11_shared_component_contract/README.md)
sourced an outside developmental list, the signature of a mixed type 1 and type 2
population in normal lung at postnatal day 1 (Guo et al. 2019), and froze a
94-gene development-specific module, retained by the owner. At list level the
developmental list shares only five genes with the adult injury list, four of them
type 1 identity genes, and none with the lesion list. List overlap is a
conservative measure, so this does not show the programmes differ.

**Rivals and test.** Generic stress, cycling, age/genotype imbalance and shared
label genes can explain overlap. Freeze source-defined shared and context-specific
components excluding selection genes; test effects beyond generic stress/cycling
in independent developmental and injury contrasts. Use animals as units, retain
genotype and validate outside the component-selection data. Type 1-directed
identity is now the leading rival: 43 of the 94 development-specific genes are
type 1 or type 2 identity genes. A shared transcriptional component would also not
imply a shared route, since neonatal injury is reported to regenerate by type 1 to
type 2 reprogramming (Penkala et al. 2021).

**Decision / readiness.** The revised external Guo test is complete: 24 mice,
transitional versus activated AT2, mean +0.735 detection percentage points
(95% CI 0.579–0.892). Identity-excluded (57 genes) and identity/control-excluded
(53 genes) effects remain positive after Holm correction. This supports partial
transcriptional recruitment, not a shared lineage or repair outcome. The original
94/51-gene modules were filtered with Strunz test-cohort markers and are descriptive
there. Strunz's published poor overall developmental correspondence remains relevant
counterevidence to global equivalence. [Plan](RQ_Specified/A5_developmental_programme_reuse/PLAN.md);
[results and biological interpretation](RQ_Specified/A5_A11_shared_component_contract/reports/REVISED_TEST_RESULTS.md).
[Current figures](analysis/figures/rq/README.md#a5) motivate animal-level
shared-versus-specific effects and held-out evaluation.
[MC1–MC2](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1), [MC5](docs/RQ_MEASUREMENT_CONTRACTS.md#mc5).

<a id="a6-how-much-of-the-ipf-macrophage-proliferation-signal-is-composition-and-what-remains-within-a-shared-noncycling-state"></a>
<a id="a6"></a>

### A6. Does IPF alter shared macrophage states beyond changing their abundance?

**Hypothesis.** IPF has a within-state macrophage programme component in addition
to changes in resident, recruited and proliferating cell proportions.

**Our observation.** Validation-cohort proliferating macrophage fractions average
3.77% in IPF versus 2.11% in controls. Cell and programme-RNA contributions differ,
motivating a biological decomposition. Only three IPF donors and one control
exceed the proliferating-state floor, and cohort labels are not harmonized.
[Existing enrichment](analysis/corrections/statistics/README.md) does not establish
the within-state hypothesis.

**Rivals and test.** Mixture alone, inconsistent annotation and batch effects
compete with a within-state change. Harmonize states independently of the tested
programme, compare donor pseudobulks under a frozen composition standard, and
replicate across cohorts. Do not define “noncycling” solely with the score tested.

**Decision / readiness.** A replicated effect supports the within-state component.
A narrow interval inside a prespecified equivalence margin supports a
composition-only explanation for that programme; nonsignificance does not.
Harmonization and donor coverage are gates. Extend the
[composition panels](analysis/figures/rq/README.md#a6) to donor-level standardized
effects and intervals. [MC1](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1),
[MC3](docs/RQ_MEASUREMENT_CONTRACTS.md#mc3), [MC5](docs/RQ_MEASUREMENT_CONTRACTS.md#mc5).

<a id="a7-does-cebpa-genotype-shift-the-reference-at2-population-and-compress-the-apparent-transitional-contrast"></a>
<a id="a7"></a>

### A7. Does Cebpa loss attenuate AT2 identity across states or preferentially within a transitional state?

**Hypothesis.** Cebpa loss reduces identity across AT2 states, so the smaller
transitional–reference contrast partly reflects a changed reference rather than
selective preservation of intermediates.

**Our observation.** Both reference identity and labelled–reference differences
decrease in mutant wells (C167). The P9 contrast changes from -4.59/-4.79 points
in controls to -0.70/-0.31 in mutants across technical seeds; adult injured wells
show the same qualitative compression. This motivates a genotype-wide effect,
but one well per condition cannot establish an interaction. See
[ES1](Research%20Article/epithelial_state_specificity/README.md).

**Rivals and test.** State-selective effects, genotype-dependent reference
selection and sampling remain alternatives. Estimate genotype effects in both
states and their interaction in replicated age-matched animals. Define states
independently of Sftpc and scored identity genes.

**Decision / readiness.** Comparable shifts with a precisely bounded interaction
support broad attenuation; a replicated meaningful interaction supports a
state-selective component. Both can coexist. Current one-well measurements are
descriptive. Extend [two-population plots](analysis/figures/rq/README.md#a7) to
animal effects and interaction intervals. [MC1–MC3](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1).

<a id="a8-does-a-broad-at1-score-capture-shared-transition-programmes-rather-than-late-maturation"></a>
<a id="a8"></a>

### A8. Does a maturation-specific programme add information about mature AT1 contribution beyond shared transition?

**Hypothesis.** A component beyond shared transition provides additional
information about measured mature AT1 contribution.

**Our observation.** ADI and AT1 holdout lists share 119 genes, and the small
late-AT1 panel is seed-sensitive (C168). This definitional dependence motivates
separation but is weak biological evidence. It does not show that mature AT1
identity is merely an extension of transition.

**Rivals and test.** Shared transition alone, timing, mixture or measurement
quality may explain the endpoint. Freeze disjoint/shared components and test
incremental association or prediction against measured AT1 protein, morphology
or traced descendant yield in independent animals. Do not define both predictor
and “mature” outcome with the same RNA panel.

**Decision / readiness.** Added information that transports supports the nominated
component; a precise absence of a meaningful increment weakens it. Missing
independent endpoints or imprecision is inconclusive. The
[A1 lineage audit](RQ_Specified/A1_transitional_epithelial_state_distinction/LINEAGE_AUDIT.md)
guides sourcing. [Overlap figures](analysis/figures/rq/README.md#a8) are diagnostics;
desired panels compare component effects and held-out endpoint performance.
[MC1](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1), [MC5](docs/RQ_MEASUREMENT_CONTRACTS.md#mc5).

<a id="a9-does-apparent-egfr-ligand-specificity-reflect-receiver-biology-or-receptor-representation-and-coverage"></a>
<a id="a9"></a>

### A9. Does fibroblast receptor context determine the response to AREG?

**Hypothesis.** Receptor abundance and complex competence modify fibroblast
responses to a defined AREG exposure.

**Our observation.** Exact EGFR and EGFR_ERBB2 definitions have different donor
coverage: the complex is scored in 7/22 donors, below the 11-donor retention rule.
This nominates a competence question and identifies a measurement prerequisite;
it does not measure dimer composition or support EGFR-homodimer predominance.

**Rivals and test.** Resource representation, depth, fibroblast mixture and
other pathways compete with receiver biology. Following a within-state RNA
screen, measure receptor protein/activation and matched ligand responses under
receptor-specific perturbation across independent preparations. Complex
composition needs an appropriate direct assay.

**Decision / readiness.** A receptor-specific outcome change with verified
engagement supports context dependence; its precise absence weakens the nominated
mechanism. RNA/resource discordance cannot decide it. The
[ligand report](analysis/corrections/ligand/RESULTS.md) and
[coverage panels](analysis/figures/rq/README.md#a9) support screening; primary
future panels are protein/activation and functional contrasts.
[MC3–MC4](docs/RQ_MEASUREMENT_CONTRACTS.md#mc3).

<a id="a10-do-epithelial-perturbation-responses-predict-organoid-growth-and-fibroblast-responses-across-independent-preparations"></a>
<a id="a10"></a>

### A10. Do epithelial programmes add information about measured organoid growth?

**Hypothesis.** Epithelial perturbation programmes add information about growth
beyond baseline imaging and plate effects; fibroblast response programmes may
provide a separate increment.

**Observation.** GSE307112 provides species-separated epithelial/fibroblast RNA
and well-linked imaging. The [revised analysis](RQ_Specified/A10_organoid_growth_outcome/reports/STAGE4_REVISED_REPORT.md)
uses 885 wells in 15 deposited plate-replicate groups. Epithelial growth programmes
add 0.0426 under the within-unit metric, above the declared 0.02 margin in all
four sensitivity settings. The fibroblast increment is -0.0154. Only 8/15 groups
have a positive epithelial increment, with the pooled result concentrated in
plates 1 and 3. The [dataset gate](docs/NEXT_DATASET_GATE.md) preserves the initial
selection rationale; it is no longer the next unexecuted task.

The [completed follow-up](RQ_Specified/A10_organoid_growth_outcome/reports/FOLLOWUP_RESULTS.md)
finds only four targets shared between plates; 886 GEO sample records still lack
preparation IDs. The six growth scores add beyond E2F/G2M: 4.84% less reference
squared error within groups and 36.96% less under joint plate/target shift.
Full-growth plate-shift improvement is 42.39%, positive in each held-out plate,
but absolute R-squared is negative on 3/4 plates. These percentages use a new
reference-error metric, not the original delta-R-squared. Both fixed outcome
scales and the drop-both sensitivity support the conditional addition.

**Rivals and test.** Initial size, plate/preparation, guide effects or mixture
may explain growth. The completed analysis joined well metadata and fixed day-14
area conditional on day-7 area; number/coverage were nominated as secondary
endpoints, not substituted for the primary. Entire deposited groups were held
out, but independent preparations remain unidentified. Outcome-time RNA supports
concurrent association. The target-transcript check is a limited diagnostic,
not validation of all perturbations or of a positive-control imaging response.

**Decision / readiness.** The revised rule is met descriptively, not as independent
confirmation: it followed the first analysis on the same screen, and centring
uses each held-out group's own outcome mean. The grouping labels do not establish
independent preparations. No confidence bound establishes precise absence of a
fibroblast contribution. Organoid mean area measures size/morphology, not total
tissue production, mature AT1 fate or in vivo repair. The follow-up comparisons
are complete and remain descriptive. Additional model/target extensions were
pruned: resolve preparation mapping, imaging calibration and a matched validation
design before another model-selection cycle.
[Original design schematic](analysis/figures/rq/README.md#a10);
[MC1](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1),
[MC5](docs/RQ_MEASUREMENT_CONTRACTS.md#mc5).

<a id="a11-shared-plasticity-versus-neoplasia-associated-context"></a>
<a id="a11"></a>

### A11. Which lesion-associated programmes add to a shared epithelial plasticity component?

**Hypothesis.** Neoplasia-associated epithelial states contain a shared
remodelling component plus context-associated additions distinguishable from
non-neoplastic repair and fibrosis.

**Our observation.** The overlap-reduced HPCS score rises in 7/7 pooled
repair/development libraries, 3/3 IPF donor pairs and 19/23 LUAD paired comparisons
(mean +0.327 log2 CPM). This positively motivates shared transcriptional biology;
it establishes neither one cell identity nor the absence of lesion-specific
additions. The [evidence review](Research%20Article/gate2_C3_yu_lee_choi_min_2026/EVIDENCE_REVIEW.md)
preserves each context's units and limitations. The 0.327 figure is the discovery
run's broad type 2 compartment; its narrower type 2 label gives 0.323 with the same
19 positive patients.

The [shared component contract](RQ_Specified/A5_A11_shared_component_contract/README.md)
froze a lesion-specific module that is identical, gene for gene, to this
overlap-reduced score, so the 23-patient result is discovery and cannot also
evaluate. At list level the adult injury programme shares seven genes with the
lesion list, mostly stress genes such as the p53 targets Bax and Gdf15, and the
developmental list shares none.

**Rivals and test.** Generic stress, cycling, annotation and composition can
mimic sharing or specificity. Define shared and candidate context-associated
components in discovery data, then evaluate frozen additions in independent
within-study contrasts with paired patients and supported epithelial states.
Require independent identity evidence before calling a cell malignant;
histology groups are not longitudinal progression stages.

**Decision / readiness.** The amended Kim test is complete after reproducing all
46 broad/narrow discovery contrasts. All eight patient differences are positive:
HL +0.681 log2 CPM (exact 95% CI 0.386–0.976; p=0.0078125), with the interval above
the pragmatic 0.10 margin. Stress exclusion remains positive (BH q=0.0234), but
beyond-shared does not meet its declared threshold (HL +0.549; CI −0.034–1.063;
BH q=0.0547). Thus lesion association replicates and the stronger relative-activation
criterion remains unresolved. Changed population definition and lack of a
non-neoplastic injury comparator limit specificity claims. [Plan](RQ_Specified/A11_lesion_programme_addition/PLAN.md);
[results](RQ_Specified/A5_A11_shared_component_contract/reports/REVISED_TEST_RESULTS.md).
[Current PCA/paired panels](analysis/figures/rq/README.md#a11)
motivate disjoint-component heatmaps and held-out patient effects.
[MC1](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1), [MC3](docs/RQ_MEASUREMENT_CONTRACTS.md#mc3),
[MC5](docs/RQ_MEASUREMENT_CONTRACTS.md#mc5).

<a id="a12-recipient-context-and-il-1-specificity"></a>
<a id="a12"></a>

### A12. Does recipient receptor and inhibitor context explain responses beyond IL-1 ligand RNA?

**Hypothesis.** Fibroblast and epithelial receptor/inhibitor context contributes
to recipient programme variation beyond source IL1A/IL1B RNA and subtype mixture.

**Our observation.** IL1B compatibility differs between IPF cohorts and human
target fits vary by recipient. These patterns motivate a recipient model;
different rank denominators cannot establish stronger signalling. Human pathway
results yield 0/279 significant primary estimated-correlation tests versus
148/279 under fixed-0.01 sensitivity. The primary result remains primary.

**Rivals and test.** Ligand amount, mixture, shared inflammation and alternative
ligand families compete with the recipient explanation. Within supported states,
compare a frozen source-only model with a parsimonious receptor/inhibitor
extension using verified patient/donor units. Include IL1R1/IL1RAP and regulatory
components; retain alternative ligands, pathway eligibility and omission
sensitivity. Ligand–receptor inference and pathway enrichment explicitly support
this analysis. IL-1 specificity needs activation or selective perturbation data.

**Decision / readiness.** Replicated added information supports an association;
a precisely excluded meaningful increment weakens it. Ineligible targets or
insufficient matched units leave it unresolved. The
[methods/report](analysis/figures/rq/il1b_context/REPORT.md) and
[component, enrichment and target figures](analysis/figures/rq/README.md#a12)
are complete; the proposed joint model has not run. Future primary panels show
within-state effects and independent recipient responses.
[MC1](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1), [MC3–MC5](docs/RQ_MEASUREMENT_CONTRACTS.md#mc3).

<a id="a12-s1"></a>

**A12-S1, an enabling source-identity question.** Unassigned cells carry median
52–72% of recovered IL1B counts across human histologies. Resolving their identity
can change source attribution; an unknown label does not define a new macrophage
state. Retain this question in the [annotation contract](docs/RQ_MEASUREMENT_CONTRACTS.md#a12-s1)
and [source gallery](analysis/figures/rq/README.md#a12-s1). Cytokine secretion is
a separate measurement.

<a id="a13-fibroblast-context-and-reciprocal-niche-associations"></a>
<a id="a13"></a>

### A13. Do fibroblast programmes add information about epithelial plasticity beyond macrophage IL1B?

**Hypothesis.** Fibroblast inflammatory/recruitment, matrix and trophic programmes
contribute information about epithelial plasticity beyond macrophage IL1B RNA.

**Motivation.** Heterogeneous niche associations motivate a joint model, but
completed analyses establish neither fibroblast mediation nor reciprocal feedback.
See the [human niche report](Research%20Article/gate2_C3_yu_lee_choi_min_2026/trials/u5_human_niche/REPORT.md).

**Rivals and test.** Shared inflammation, histology, capture composition and
macrophage state can explain apparent fibroblast information. Count complete
macrophage–fibroblast–epithelial triads for the exact variables and contrast.
Freeze a small programme/multiplicity family, then compare parsimonious joint
models using within-patient changes where available. Marginal correlations do
not test “beyond IL1B.” Spatial work requires independently defined pathology
regions and a patient-level null; feedback requires intervention.

**Decision / readiness.** Replicated added information supports the association;
precise absence weakens the nominated programme. IPF has six and three complete
triads, below the ten-unit joint-model floor; 23 human cohort patients does not
mean 23 complete paired triads. Coverage is next, not mediation fitting.
[Figure plan/spatial evidence](analysis/figures/rq/README.md#a13): coverage heatmap
first, then patient effects if eligible. [MC1](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1),
[MC3–MC5](docs/RQ_MEASUREMENT_CONTRACTS.md#mc3).

<a id="a14-resolution-versus-persistence-after-signal-withdrawal"></a>
<a id="a14"></a>

### A14. Do exposure duration and fibroblast IL-1 reception separately determine recovery after withdrawal?

**Two hypotheses, separate decisions.** Longer IL-1β exposure reduces mature
epithelial recovery after withdrawal. Separately, fibroblast IL-1 reception
modifies recovery beyond direct epithelial reception. Either can hold without
the other.

**Motivation.** Shared RNA cannot distinguish reversible repair from persistent
dysfunction. This is a mechanistic follow-up, not an existing withdrawal-fate
result. Published tracing informs the design but cannot replace that experiment.

**Rivals and test.** Residual exposure, toxicity, death/replacement and an
epithelial-only response compete with duration and fibroblast dependence.
Compare transient/sustained exposure and verified withdrawal with time-matched
controls and epithelial- versus fibroblast-specific IL1R1 perturbation. Use
independent animals or culture preparations. Measure mature-cell yield, viability,
traced descendants and function alongside exposure and target engagement, with
comparable post-withdrawal intervals. Loss of a transitional score cannot by
itself distinguish maturation, reversion, death or replacement.

**Decision / readiness.** A recovery difference after verified withdrawal
supports the duration hypothesis; an outcome change under fibroblast-specific
intervention with controlled direct epithelial reception supports the second.
Precise absence weakens each separately. Failed engagement or no withdrawal
observation is inconclusive. Appropriate new data are needed; neither result
alone establishes malignant transformation.

**Figures / contracts.** [Existing design schematic](analysis/figures/rq/README.md#a14);
future primary panels: replicate-level recovery, mature-cell yield and
recipient-specific contrasts with exposure/engagement controls. No anticipated
response curve is a result. [MC1](docs/RQ_MEASUREMENT_CONTRACTS.md#mc1),
[MC3–MC4](docs/RQ_MEASUREMENT_CONTRACTS.md#mc3).

## Execution and interpretation rules

Before a new fit, freeze its estimand, biological unit, eligible observations,
primary comparison, multiplicity family, meaningful effect/prediction margin
and validation split in the question-specific plan. Sample floors are eligibility
rules, not power guarantees. Separate confirmation, exploratory discovery and
tests left inconclusive by inadequate data.

Reuse completed measurement checks within their recorded scope; repeat only if
inputs, estimands or a concrete unresolved risk change. Do not broaden a
sensitivity sweep merely to obtain significance. The [methods guide](REPRODUCIBILITY.md),
[measurement contracts](docs/RQ_MEASUREMENT_CONTRACTS.md) and
[ID crosswalk](docs/RQ_MEASUREMENT_CONTRACTS.md#crosswalk) locate requirements
without making diagnostics the biological questions.
