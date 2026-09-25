# A11: lesion-associated programme beyond shared plasticity

Question-specific work for [A11](../../RESEARCH_QUESTIONS.md#a11). The canonical
hypothesis stays in the register.

**Status: plan amended before scoring; execution authorized.** Read the [plan](PLAN.md)
and [biological rationale](../A5_A11_shared_component_contract/BIOLOGICAL_LOGIC.md).
Eight patient pairs are eligible. Direction, magnitude and biological claim are
reported separately; a positive primary alone means replicated lesion association.

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
