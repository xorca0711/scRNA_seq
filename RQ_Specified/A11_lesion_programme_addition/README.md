# A11: lesion-associated programme beyond shared plasticity

Question-specific work for [A11](../../RESEARCH_QUESTIONS.md#a11). The canonical
hypothesis stays in the register.

**Status: amended test completed.** Read the [plan](PLAN.md),
[biological rationale](../A5_A11_shared_component_contract/BIOLOGICAL_LOGIC.md) and
[results](../A5_A11_shared_component_contract/reports/REVISED_TEST_RESULTS.md).

Lesion-associated expression replicates in all eight paired patients: HL +0.681
log2 CPM (exact 95% CI 0.386–0.976; p=0.0078125). The stress-excluded sensitivity
is positive (BH q=0.0234). The beyond-shared comparison is unresolved (BH q=0.0547),
so the stronger relative-activation criterion is not met. Cancer specificity,
malignant identity and a separate mechanism remain unestablished.

## Why the test moved to a new cohort

The lesion-specific module frozen by the
[shared component contract](../A5_A11_shared_component_contract/README.md) is
identical, gene for gene, to a module the completed human run had already scored in
23 patients. That cohort is therefore the discovery, and a fresh test needs fresh
patients. The Kim 2020 deposit has ten verified tumour-normal pairs and had never
been scored with these modules.

## Layout

| Path | Contents |
|---|---|
| [PLAN.md](PLAN.md) | The pre-registration: populations, instrument, estimands, decision rules, power |
| [config/kim2020_test_contract.json](config/kim2020_test_contract.json) | The same, machine-readable |
| `scripts/00_power_calculation.py` | Power from the discovery's tracked differences only |
| `scripts/01_eligibility_gates.py` | Gene coverage and patient eligibility from names and labels only |
| `tables/` | Power, pairing, cell counts, coverage and run records |

Original gate/power outputs are retained. New results use `tables/test_v2/`.

## Reproduce the completed test

Use a clean output directory and the tracked amended configuration. Run
`Rscript scripts/02_reproduce_discovery.R REPO_ROOT DATA_ROOT`, then the Python
launcher with `scripts/03_prepare_kim.py --data-root DATA_ROOT`, then
`Rscript scripts/04_score_kim.R REPO_ROOT DATA_ROOT`. These scripts refuse to
overwrite results. The Python step verifies all original Kim input hashes and
requires the successful discovery reproduction before reading counts for scoring.
