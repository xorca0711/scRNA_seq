"""P0: audit existing source labels and sample coverage, without program scoring."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
FLOOR = 30
MIN_UNITS = 3
STATES = ["AT2", "Alveolar_transitional", "AT1_AT2", "AT1"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    out = BASE / "tables"
    out.mkdir(exist_ok=True)
    metadata = ROOT / "Thesis/gate1_01_niethamer_2025/GSE262927/tables/cell_metadata.csv"
    columns = ["cell_id", "sample_id", "gsm", "orig_ident", "experimental_group",
               "condition", "sacrifice_day", "tamoxifen_start_day", "sex",
               "author_lineage", "author_celltype"]
    cells = pd.read_csv(metadata, usecols=columns, low_memory=False)
    unannotated_condition_cells = int(cells.condition.eq("unannotated").sum())
    cells["condition"] = cells.condition.replace("unannotated", pd.NA)
    if not cells.cell_id.is_unique:
        raise ValueError("Duplicate cell IDs in existing retained-cell metadata")
    fields = ["gsm", "orig_ident", "experimental_group", "condition", "sacrifice_day",
              "tamoxifen_start_day", "sex"]
    if cells.groupby("sample_id")[fields].nunique(dropna=True).gt(1).any().any():
        raise ValueError("Sample metadata not constant within sample_id")
    # Retained cells absent from author metadata have missing annotation fields.
    # Propagate only consistent non-missing sample facts, never cell-state labels.
    samples = cells.groupby("sample_id")[fields].first()
    epi = cells.loc[cells.author_lineage.eq("Epithelium")].copy()
    long = epi.groupby(["sample_id", "author_celltype"], dropna=False).size().rename("n_cells").reset_index()
    long.merge(samples, on="sample_id", validate="many_to_one").to_csv(
        out / "d1_published_state_counts.csv", index=False)
    wide = pd.crosstab(epi.sample_id, epi.author_celltype).reindex(
        index=samples.index, columns=STATES, fill_value=0).fillna(0).astype(int)
    wide = samples.join(wide)
    # Do not merge the two distinct author intermediate labels to improve coverage.
    summaries = []
    for intermediate in ["Alveolar_transitional", "AT1_AT2"]:
        for endpoint in ["AT2", "AT1"]:
            key = f"{intermediate}_vs_{endpoint}_eligible"
            wide[key] = wide[intermediate].ge(FLOOR) & wide[endpoint].ge(FLOOR)
            for (condition, day), part in wide.groupby(["condition", "sacrifice_day"], dropna=False):
                summaries.append({"condition": condition, "sacrifice_day": day,
                    "intermediate": intermediate, "endpoint": endpoint,
                    "total_samples_in_stratum": len(part), "eligible_units": int(part[key].sum()),
                    "minimum_units": MIN_UNITS,
                    "passes_numeric_floor": bool(part[key].sum() >= MIN_UNITS)})
        wide[f"{intermediate}_triplet_eligible"] = (
            wide[[intermediate, "AT2", "AT1"]].ge(FLOOR).all(axis=1))
    wide.reset_index().to_csv(out / "d1_sample_coverage.csv", index=False)
    contrasts = pd.DataFrame(summaries)
    contrasts.to_csv(out / "d1_stratum_eligibility.csv", index=False)

    h5 = ROOT / "Research Article/gate1_01_niethamer_2025/GSE262927/regeneration_focus/alveolar_trajectory.h5ad"
    a = ad.read_h5ad(h5, backed="r")
    try:
        present = set(a.obs_names.astype(str))
        meta_ids = set(cells.cell_id)
        matrix = {"path": str(h5.relative_to(ROOT)), "n_cells": a.n_obs,
                  "n_genes": a.n_vars, "layers": list(a.layers),
                  "raw_shape": None if a.raw is None else list(a.raw.shape),
                  "obs_columns": list(a.obs.columns),
                  "obs_ids_missing_from_metadata": len(present - meta_ids),
                  "all_four_published_state_cells_in_metadata": int(epi.author_celltype.isin(STATES).sum()),
                  "published_state_cells_missing_from_cached_subset": int(
                      (~epi.loc[epi.author_celltype.isin(STATES), "cell_id"].isin(present)).sum()),
                  "expression_values_read": False}
        available = epi.loc[epi.cell_id.isin(present)]
        cached_counts = pd.crosstab(available.sample_id, available.author_celltype).reindex(
            index=samples.index, columns=STATES, fill_value=0).fillna(0).astype(int)
        samples.join(cached_counts).reset_index().to_csv(out / "d1_cached_matrix_state_counts.csv", index=False)
    finally:
        a.file.close()
    if matrix["obs_ids_missing_from_metadata"]:
        raise ValueError("Cached matrix has cell IDs absent from source metadata")
    result = {"stage": "P0", "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)), "metadata_sha256": sha256(metadata),
        "matrix_sha256": sha256(h5), "python": platform.python_version(),
        "packages": {p: importlib.metadata.version(p) for p in ["pandas", "anndata", "h5py"]},
        "cell_floor": FLOOR, "unit_floor": MIN_UNITS, "retained_cells": len(cells),
        "missing_author_celltype": int(cells.author_celltype.isna().sum()),
        "unannotated_condition_cells_treated_as_missing": unannotated_condition_cells,
        "missing_sample_facts_before_sample_aggregation": cells[fields].isna().sum().to_dict(),
        "samples": len(samples), "epithelial_cells": len(epi),
        "epithelial_labels": epi.author_celltype.value_counts().to_dict(),
        "complete_triplet_units_all_times": {s: int(wide[f"{s}_triplet_eligible"].sum())
                                              for s in ["Alveolar_transitional", "AT1_AT2"]},
        "eligible_contrasts": contrasts.loc[contrasts.passes_numeric_floor].to_dict("records"),
        "matrix": matrix, "program_scoring_performed": False,
        "interpretation": "Coverage only; sample independence and biological definitions require source audit."}
    (BASE / "local_audit_record.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ["retained_cells", "samples", "epithelial_cells",
          "epithelial_labels", "complete_triplet_units_all_times", "eligible_contrasts"]}, indent=2))
    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
