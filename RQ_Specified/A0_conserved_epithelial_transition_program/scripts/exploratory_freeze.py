"""Freeze the exploratory estimators using feature IDs and metadata only."""
import gzip
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from audit_extended_coverage import decode_annotations
from exploratory_common import BASE, ROOT, TABLES, save_json, sha, excluded


def main():
    path = BASE / "exploratory_config.json"
    if path.exists():
        print("Existing frozen configuration retained; no changes made.")
        return
    folder = BASE / "cache/expression"
    repair_genes = gzip.open(folder / "GSE141259_HighResolution_genes.txt.gz", "rt").read().splitlines()
    dev = decode_annotations((folder / "Negretti_var.bin").read_bytes())
    dev.to_csv(TABLES / "developmental_feature_index.csv", index=True, index_label="feature_index")
    with gzip.open(folder / "GSE92332_atlas_UMIcounts.txt.gz", "rt") as f:
        header = next(f)
        intestinal_genes = [line.split("\t", 1)[0].strip('"\n') for line in f]
    labels = "Krt8 Krt18 Krt19 Cldn4 Cdkn1a Lgals3 Sfn Sftpc Sftpa1 Sftpb Sftpd Abca3 Etv5 Ager Pdpn Hopx Aqp5 Cav1 Rtkn2 Clic5 Lmo7 Lcn2 Il33 Mdk Sox9 Scgb1a1 Scgb3a2 Foxj1 Dynlrb2 Ascl1 Lgr5 Olfm4 Ascl2 Smoc2 Lrig1 Alpi Apoa1 Apoa4 Fabp1 Fabp2 Sis Slc5a1 Mki67 Top2a".split()
    config = {"version": "A0-exploratory-2026-09-26", "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "expression_effects_seen_at_freeze": False, "seed": 0, "label_genes": sorted(labels),
        "label_gene_scope": "Conservative known source/endpoint/transfer marker exclusion; not an exhaustive removal of annotation circularity",
        "minimum_adi_cell_detection": .10, "minimum_log2cpm_difference": .25,
        "minimum_program_genes": 20, "maximum_program_genes": 50, "random_gene_sets": 200,
        "cell_floor_for_qualified_groups": 30, "minimum_feature_coverage": .80,
        "score": "mean probability program genes outrank fixed reference genes within each cell; ties=0.5; remove module/reference overlap",
        "E1_gate": {"minimum_positive_heldout_mice_each_endpoint": 6, "minimum_median_standardized_effect_each_endpoint": .25},
        "transfer": {"developmental_ages": ["P0","P3","P5","P7","P14"],
            "developmental_states": ["AT2","Transitional","AT1"],
            "intestinal_states": ["Stem","Enterocyte.Immature.Proximal","Enterocyte.Mature.Proximal"],
            "intestinal_batches_to_mouse": {"B3":"Control-Mouse1","B4":"Control-Mouse2","B7":"Control-Mouse3","B8":"Control-Mouse4"}},
        "generic_control_names": ["HALLMARK_HYPOXIA","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_P53_PATHWAY",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","Immediate_early_AP1_panel"],
        "source_hashes": {str(p.relative_to(ROOT)):sha(p) for p in [
            folder / "GSE141259_HighResolution_rawcounts.mtx.gz", folder / "GSE141259_HighResolution_genes.txt.gz",
            folder / "GSE141259_HighResolution_barcodes.txt.gz", folder / "Negretti_var.bin",
            folder / "GSE92332_atlas_UMIcounts.txt.gz", BASE / "cache/sources/GSE141259_HighResolution_cellinfo.csv.gz",
            ROOT / "Research Article/epithelial_state_specificity/modules.json", ROOT / "raw_data/msigdb/mh.all.v2024.1.Mm.symbols.gmt"]}}
    common = sorted((set(repair_genes) & set(dev.gene) & set(intestinal_genes)) - set(labels))
    common = [g for g in common if not excluded(g,config)]
    assert len(common) >= 1000
    config["reference_genes"] = sorted(np.random.default_rng(0).choice(common, 1000, replace=False).tolist())
    config["reference_universe"] = {"repair_features":len(repair_genes), "developmental_features":len(dev),
        "intestinal_features":len(intestinal_genes), "common_nonexcluded_features":len(common),
        "selection": "Feature identity intersection only; no expression values inspected for reference selection"}
    coverage = pd.read_csv(BASE / "tables/d1_strunz_proposed_window.csv")
    config["repair_units"] = sorted(coverage.loc[coverage.coverage_pass_in_proposed_window,"identifier"].tolist())
    assert len(config["repair_units"]) == 9
    modules = json.loads((ROOT / "Research Article/epithelial_state_specificity/modules.json").read_text())["modules"]
    control = {x["name"]:x["genes"] for x in modules}
    for line in (ROOT / "raw_data/msigdb/mh.all.v2024.1.Mm.symbols.gmt").read_text().splitlines():
        name, url, *genes = line.split("\t")
        if name in config["generic_control_names"]:
            control[name] = sorted(set(genes))
    control["Immediate_early_AP1_panel"] = ["Fos","Fosb","Jun","Junb","Jund","Atf3","Egr1","Dusp1"]
    config["control_modules"] = control
    config["custom_panel_note"] = "Immediate_early_AP1_panel is an analyst-prespecified short marker panel, not a complete pathway or causal activity assay"
    save_json(path, config)
    print(json.dumps({"frozen":True,"repair_mice":9,"common_reference_universe":len(common),"config_sha256":sha(path)}))


if __name__ == "__main__":
    main()
