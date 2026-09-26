# A0 continuation: eligibility before expression discovery

26 September 2026. The owner requested completion of the scientific pilot and
a PR. Commit `4ed38ff` preserves the previously untracked audit unchanged.
This continuation does not treat an eligibility failure as a negative biological
result or permission to manufacture replication.

## References checked before launch

- The original A0 plan permits separately eligible paired endpoint contrasts;
  the later reproduction guide's complete-triplet requirement was stronger.
  Reassess the two contrasts at the unchanged 30-cell / three-unit floors.
- A5's revised Guo signature test is complete in Strunz repair mice. It tests
  partial signature recruitment, not A0's discovery-plus-cross-tissue transfer.
  Guo's source paper contains two mice, below the A0 discovery-unit floor.
- A1 supplies regulatory and fate context, but no matched cross-tissue discovery
  dataset. A5/A11 overlap is pairwise, not a universal three-context programme.
- A10's growth scores do not test a conserved transition programme or mature
  fate. Its plate holdouts do not fill A0's developmental or transfer roles.

## Recovery sequence and limits

1. Recalculate both endpoint eligibility counts from original A0 tables.
2. Inspect additional source-supported mouse developmental data. The 2026
   Hippo study is a candidate; establish normal biological units before counts.
3. If mouse candidates fail, inspect human developmental atlases as allowed by
   the original plan's cross-species option. He et al. 2022 is the first candidate,
   with late tip / late stalk / AT1 as the alveolar branch. Inspect donor, library,
   age, dissection and processing metadata, retaining their original labels.
   Quach et al. 2024 is the next candidate if the first fails.
4. Inspect an independent non-lung comparison only after a viable developmental
   role exists. Do not rerun the unchanged old transfer audit unnecessarily.
5. Freeze P1 only if each required contrast has at least three identifiable
   independent units and 30 cells/state, and the branch has biological support
   beyond pseudotime. Use matched processing/region strata where needed. No
   gene expression effects have been inspected for the new cohorts at this stage.
6. If eligible, commit the exact ortholog, gene-filter, label-exclusion,
   ranking, scoring and sensitivity specification before reading expression
   effects. Run original P2-P4 conditionally: a discovery failure may make
   transfer and specificity analyses unnecessary. Preserve every outcome.

Selection may use study design and state coverage, never candidate programme
scores. Changes of branch require source evidence and a dated decision before
scoring; they cannot silently substitute injury for normal development, pooled
cells for biological units, or early commitment for a mature endpoint.

The search is bounded to concrete candidates with accessible evidence. If the
input gate remains unsatisfied, preserve the completed recovery work in a PR
and explicitly state that the scientific pilot is incomplete. The instruction
to finish cannot turn unavailable observations into a scientific result.

## Metadata-driven amendment before discovery

The first two human candidates did not yield an executable cohort: He's alveolar
contrasts failed coverage; Quach's smaller epithelial object returned HTTP 403
to anonymous download, while the public whole-lung object is 7.98 GB. A third
candidate explicitly cited by Quach, Sountoulidis et al. (2023), exposes matched
counts and donor/state metadata through the paper-linked UCSC `lung-dev` entry.
The initial similarly named `fetal-lung` URL actually describes Miller 2020;
its metadata are not attributed to Sountoulidis or used as the selected cohort.
An initially guessed GSE215896 accession was likewise rejected on its title;
the authoritative SuperSeries is GSE215898, with scRNA child GSE215895.

The selected D2 compares author `epi_cl2` (SOX9-high/ETV5-medium distal),
`epi_cl1` (intermediate), and `epi_cl0` (proximal secretory) within each donor.
Eleven donors pass all three 30-cell floors. This is an early airway developmental
axis, replacing the unworkable alveolar branch, not proof of an alveolar lineage
shared with repair. The source supports position by spatial assays and a
developmental interpretation by trajectories; it does not trace each sequenced
intermediate to its endpoint. The original pilot allowed different starting and
destination identities, and this source supplies evidence beyond pseudotime.
The chosen operational test remains a state association, not a transition rate.

For V1, retain the four explicitly mapped Haber control mice and prospectively
change the intermediate to author `Enterocyte.Progenitor.Late`. The start remains
`Stem`, and the destination remains `Enterocyte.Mature.Proximal`. Three mice pass
coverage. This tests an earlier enterocyte commitment stage than the historical
immature-proximal-cell proposal. Late progenitors are not claimed to be exclusively
proximal. The source's differentiation analysis and crypt/villus biology motivate
the comparison; donor linkage of a particular progenitor to its progeny is absent.
This choice uses labels and coverage only, with all gene-programme effects unseen.
The old two-mouse result remains in the historical audit, unchanged.

`config/pilot_v1.json` fixes these choices and all discovery/transfer rules.
The exact programme may be absent: require positive median differences of at
least 0.1 log2(CPM+1), positive directions in at least two-thirds of units, and
intermediate CPM >=1 in at least two-thirds of units, separately in both lung
contexts and against both endpoints. The 0.1 margin is a pragmatic selection
criterion, not a validated biological effect threshold. Do not change it after
seeing how many genes pass. Fewer than 20 qualifying genes stops this specified
single-module pilot; V1 scoring would then have no frozen instrument to test.
