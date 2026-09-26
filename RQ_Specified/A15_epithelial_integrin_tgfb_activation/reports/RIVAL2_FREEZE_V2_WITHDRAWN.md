# Why the second rival-2 freeze is withdrawn, also before anything was scored

27 September 2026. [`config/a15_rival2_freeze_v2.json`](../config/a15_rival2_freeze_v2.json)
is **withdrawn and preserved unedited**.
[`config/a15_rival2_freeze_v3.json`](../config/a15_rival2_freeze_v3.json) is the authority.
**Nothing had been scored under v2.**

A second three-lens red team read v2 against the fatal findings that withdrew v1. All three
returned *needs revision*: v2 fixed the ratio decisively, installed a real engagement
control, and demoted the A0 set correctly, but it relocated rather than removed the construct
problem and introduced three defects of its own.

| Lens | Verdict | Fixed v1 fatal 1 | fatal 2 | fatal 3 |
|---|---|---|---|---|
| Construct validity | needs revision | partially | partially | yes |
| Inference | needs revision | partially | partially | yes |
| Coherence | needs revision | partially | partially | yes |

Two of their quantitative claims were **reproduced by simulation in this repository before
being accepted**, because a reviewer's number is not evidence until it is checked.

## 1. The omnibus statistic confounded location with dispersion

v2's omnibus was the mean between-arm distance minus the mean within-arm distance. Simulated
here over 200 draws with eight mice and 200 genes:

| Scenario | v2 statistic rejects | Centroid distance rejects | Median within-arm distance ratio |
|---|--:|--:|--:|
| Pure dispersion difference, no gene's mean moving | **1.000** | 0.635 | 2.77 |
| True null, equal dispersion | 0.035 | 0.035 | 1.01 |
| True location shift | 1.000 | 1.000 | 1.01 |

A statistic that rejects every time on a pure between-mouse variance change cannot be read as
"the epithelium changed". v3 makes the **centroid distance** the primary, and gates it on a
**dispersion diagnostic** with a declared threshold of 1.5, which the same simulation shows
separates the pathological case cleanly. The v2 statistic is retained and reported with no
decision attached.

## 2. The Hodges-Lehmann interval under-covered

v2 trimmed one pairwise difference from each end and labelled the result 97.14 per cent.
Simulated coverage at four against four, 20,000 draws:

| Trim from each end | Empirical coverage | 1 - 2 P(U <= k) |
|---|--:|--:|
| 0 | 0.9716 | 0.9714 |
| 1 | 0.9417 | 0.9429 |

So the interval v2 reported as 97.14 per cent actually covered **94.18 per cent**, and the
interval that attains the level at these sizes is the **full range** of the sixteen pairwise
differences. v3 selects the trim as the largest k satisfying 1 - 2 P(U <= k) >= 0.95, which
is k = 0 here, and reports the attained coverage.

## 3. Three defects v2 introduced

- **No direction for the contrast actually tested.** v2's composites carried a direction for
  bleomycin against saline, while the contrast tested is 3G9 against Axum8, both bleomycin.
  The reading then demanded separation "in the declared direction" for a contrast with no
  declared direction, and the executor ignored direction entirely. v3 declares both primary
  endpoints **two-sided**, because neither the A15 mechanism nor rival 2 predicts a sign, and
  keeps a direction only for the engagement control, where Horan 2008 supplies one.
- **`Col1a1` sat in both the engagement control and the purity downgrade rule**, so a
  successful engagement result would have pushed the run toward inconclusive through its own
  purity check. v3 carries mesenchymal de-enrichment on `Col3a1` instead.
- **The positive branch waived the engagement control**, declaring that a change establishes
  engagement by itself. An off-target antibody effect would also produce a change. v3 removes
  the waiver.

## 4. The branch names overclaimed

v2 called its positive branch "rival 2 stays live". Both the A15 mechanism and rival 2
predict an epithelial change, so a positive discriminates neither, and the name implied
otherwise. v3 renames it
**`epithelium_differs_but_does_not_discriminate`** and states in the freeze, before
execution, that the likely outcome is the uninformative one and that only the
`weak_bound_on_rival_2` branch discriminates.

## 5. One further defect, found by execution rather than review

The first execution attempt produced `purity breaches: none` when the purity markers **had
never been mapped**: the Ensembl symbol lookup was built during stage 1 and never extended to
the panels v3 declares, so every covariate silently returned "not recovered" and the run
reported vacuous silence as a negative result. Script 01 now requests every symbol any freeze
version declares, script 02 **fails closed** if a declared panel does not map, and an
uncomputed covariate now forces inconclusive rather than passing unnoticed. The engagement
control was refused in that attempt for the same reason; the corrected run is the one
reported.

## What v3 keeps

The design audit, the arms, the batch and sex matching, the join and its agreement with A1's
records, the paired-input purity covariates, the engagement control's gene set and direction,
the label-free filters, the no-gate single-pass execution and the quantified power statement
all carry over unchanged.
