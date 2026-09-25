"""Reconstruct Kobayashi ED4 endpoints without treating fields as mice."""
from __future__ import annotations

import hashlib
import json
import re
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def fraction(numerator, denominator):
    if numerator < 0 or denominator < 0 or numerator > denominator:
        raise ValueError("Invalid labelled-cell count")
    return 100 * numerator / denominator if denominator else np.nan


def extract_fields(path):
    rows = []
    with ZipFile(path) as z:
        strings = ["".join(n.itertext()) for n in ET.fromstring(z.read("xl/sharedStrings.xml"))]
        tree = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        for row in tree.findall(".//m:row", NS):
            cells = {}
            for cell in row:
                v = cell.find("m:v", NS)
                if v is not None:
                    cells[re.sub(r"\d", "", cell.get("r"))] = strings[int(v.text)] if cell.get("t") == "s" else v.text
            match = re.fullmatch(r"K19-creER-tdt (ctr|BleoD12)_(AGER|KRT8)_tdt_Mouse([123])", cells.get("A", ""), re.I)
            if not match:
                continue
            condition, marker, mouse = match.groups()
            numerator, denominator = int(cells["C"]), int(cells["B"])
            pct = fraction(numerator, denominator)
            stored = float(cells["D"])
            if denominator and not np.isclose(pct, stored):
                raise ValueError(f"Stored/recomputed mismatch in source row {row.get('r')}")
            rows.append(dict(source_row=int(row.get("r")), source_label=cells["A"],
                             condition="BleoD12" if condition.lower() == "bleod12" else "Control",
                             marker=marker.upper(), mouse=mouse, numerator=numerator,
                             denominator=denominator, source_percent=stored, percent=pct,
                             valid_denominator=bool(denominator)))
    if not rows:
        raise ValueError("No explicitly mouse-labelled source rows found")
    return pd.DataFrame(rows)


def main():
    if (BASE / "reports/lineage_run.json").exists():
        raise SystemExit("Existing lineage run: archive its outputs before a new numerical run.")
    inventory = json.loads((BASE / "reports/lineage_source_inventory.json").read_text())
    source = next(f for f in inventory["files"] if f["study"] == "pats_nature")
    path = ROOT / source["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == source["sha256"]
    fields = extract_fields(path)
    rows = []
    for (condition, marker, mouse), group in fields.groupby(["condition", "marker", "mouse"]):
        rows.append(dict(condition=condition, marker=marker, mouse=mouse, fields=len(group),
                         fields_with_denominator=int(group.valid_denominator.sum()),
                         numerator=int(group.numerator.sum()), denominator=int(group.denominator.sum()),
                         mean_field_percent=group.percent.mean(),
                         pooled_count_percent=fraction(group.numerator.sum(), group.denominator.sum())))
    mice = pd.DataFrame(rows)
    valid = mice[mice.condition == "BleoD12"]
    assert valid.groupby("marker").size().to_dict() == {"AGER": 3, "KRT8": 3}
    table = BASE / "tables/lineage"; table.mkdir(parents=True, exist_ok=True)
    fields.to_csv(table / "pats_source_fields.tsv", sep="\t", index=False, na_rep="NA")
    mice.to_csv(table / "pats_mouse_endpoints.tsv", sep="\t", index=False, na_rep="NA")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.6), sharey=True)
    for ax, marker, color in zip(axes, ["KRT8", "AGER"], ["#1F718B", "#BE6835"]):
        sub = valid[valid.marker == marker].sort_values("mouse")
        x = np.arange(3)
        ax.scatter(x - .07, sub.mean_field_percent, color=color, s=58, label="Mean of field fractions", zorder=3)
        ax.scatter(x + .07, sub.pooled_count_percent, facecolors="white", edgecolors=color, s=58, marker="s", label="Pooled counts within mouse", zorder=3)
        ax.set(xticks=x, xticklabels=[f"Mouse {m}" for m in sub.mouse], ylim=(0, 100), title=f"{marker}+ among alveolar tdTomato+ cells")
        ax.grid(axis="y", alpha=.18)
    axes[0].set_ylabel("Labelled-cell endpoint fraction (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, .065), ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Krt19 lineage endpoint after bleomycin", fontsize=14)
    fig.text(.5, .91, "Kobayashi et al. 2020 · source labels: BleoD12 · 3 mice per marker", ha="center", fontsize=10)
    fig.text(.5, .02, "Control fractions are undefined (0 labelled cells); fields are nested within mice. No between-marker pairing assumed.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .17, 1, .87))
    for ext in ("png", "svg"):
        fig.savefig(BASE / f"figures/a1_pats_lineage_endpoints.{ext}", dpi=220)
    plt.close(fig)
    summary = {"status": "completed", "utc": datetime.now(timezone.utc).isoformat(), "source": source,
               "input_fields": len(fields), "zero_denominator_fields": int((~fields.valid_denominator).sum()),
               "mean_mouse_field_percent": valid.groupby("marker").mean_field_percent.mean().to_dict(),
               "claim": "Reconstruction of published labelled-cell endpoint composition, not a new fate experiment or exit-rate estimate",
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (BASE / "reports/lineage_run.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("source", "script_sha256")}, indent=2))


if __name__ == "__main__":
    main()
