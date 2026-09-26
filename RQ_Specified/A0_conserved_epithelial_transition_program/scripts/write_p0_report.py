"""Generate the P0 report from audited counts plus the reviewed interpretation template."""
import json
from pathlib import Path
from string import Template
import pandas as pd

BASE = Path(__file__).resolve().parents[1]


def main():
    local = json.loads((BASE / "local_audit_record.json").read_text())
    public = json.loads((BASE / "public_audit_record.json").read_text())
    values = {"local_cells": f"{local['retained_cells']:,}",
              "local_libraries": local["samples"], "local_epi": f"{local['epithelial_cells']:,}",
              "local_triplets": local["complete_triplet_units_all_times"]["Alveolar_transitional"],
              "alternative_triplets": local["complete_triplet_units_all_times"]["AT1_AT2"],
              "matrix_cells": f"{local['matrix']['n_cells']:,}",
              "matrix_genes": f"{local['matrix']['n_genes']:,}",
              "matrix_omissions": local["matrix"]["published_state_cells_missing_from_cached_subset"],
              "four_state_cells": f"{local['matrix']['all_four_published_state_cells_in_metadata']:,}"}
    for key, label in [("at2", "AT2"), ("transitional", "Alveolar_transitional"),
                       ("at1_at2", "AT1_AT2"), ("at1", "AT1")]:
        values[key] = f"{local['epithelial_labels'][label]:,}"
    for prefix, source in [("hi", "strunz_high_resolution"), ("whole", "strunz_whole_lung")]:
        for key in ["cells", "samples"]:
            values[f"{prefix}_{key}"] = f"{public[source][key]:,}"
    values["hi_triplets"] = public["strunz_high_resolution"]["triplet_units_all_times"]
    values["whole_at1"] = public["strunz_whole_lung"]["AT1_cells_total"]
    values["gut_cells"] = f"{public['haber_atlas']['header_cells']:,}"
    values["gut_batches"] = public["haber_atlas"]["batches"]
    values["gut_triplets"] = public["haber_atlas"]["triplet_units"]
    gut = pd.read_csv(BASE / "tables/v1_haber_verified_mouse_coverage.csv")
    lines = ["| Verified mouse | Batch | Stem | Immature proximal | Mature proximal | Complete at 30 cells/state |",
             "|---|---|---:|---:|---:|---|"]
    for row in gut.to_dict("records"):
        lines.append(f"| {row['mouse_id']} | {row['batch']} | {row['Stem']} | "
                     f"{row['Enterocyte.Immature.Proximal']} | {row['Enterocyte.Mature.Proximal']} | "
                     f"{'Yes' if row['triplet_eligible'] else 'No'} |")
    values["gut_coverage_table"] = "\n".join(lines)
    template = Template((BASE / "scripts/p0_report.md.in").read_text(encoding="utf-8"))
    output = template.substitute(values)
    (BASE / "reports/P0_ELIGIBILITY_REPORT.md").write_text(output, encoding="utf-8")
    print("Generated reports/P0_ELIGIBILITY_REPORT.md from audit records and reviewed interpretation template")


if __name__ == "__main__":
    main()
