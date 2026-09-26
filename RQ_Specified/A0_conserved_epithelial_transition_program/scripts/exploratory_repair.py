"""E1: learn a repair-associated candidate and test leave-one-mouse-out stability."""
import gzip
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from exploratory_common import BASE, PROC, TABLES, save_json, sha, rank_scores, paired_effects, pseudobulk, learn_program


def main():
    config = json.loads((BASE / "exploratory_config.json").read_text())
    folder = BASE / "cache/expression"
    genes = gzip.open(folder / "GSE141259_HighResolution_genes.txt.gz", "rt").read().splitlines()
    barcodes = gzip.open(folder / "GSE141259_HighResolution_barcodes.txt.gz", "rt").read().splitlines()
    assert len(genes) == len(set(genes)) and len(barcodes) == len(set(barcodes))
    annotation = pd.read_csv(BASE / "cache/sources/GSE141259_HighResolution_cellinfo.csv.gz", sep="\t").set_index("cell_barcode")
    assert set(barcodes) == set(annotation.index)
    obs = annotation.loc[barcodes].reset_index().rename(columns={"identifier":"unit","cell_type":"state","time_point":"time"})
    with gzip.open(folder / "GSE141259_HighResolution_rawcounts.mtx.gz", "rb") as stream:
        matrix = mmread(stream).tocsr().astype(np.float32)
    assert matrix.shape == (len(obs),len(genes)), matrix.shape
    assert np.isfinite(matrix.data).all() and (matrix.data >= 0).all() and np.equal(matrix.data,np.floor(matrix.data)).all()
    matrix.eliminate_zeros()
    selected = obs.unit.isin(config["repair_units"]).to_numpy()
    matrix, obs = matrix[selected], obs.loc[selected].reset_index(drop=True)
    obs["measured_umi"] = np.asarray(matrix.sum(axis=1)).ravel()
    obs["measured_genes"] = np.diff(matrix.indptr)
    mt = [i for i,g in enumerate(genes) if g.lower().startswith("mt-")]
    obs["measured_mito_fraction"] = np.asarray(matrix[:,mt].sum(axis=1)).ravel() / obs.measured_umi
    assert (obs.measured_umi > 0).all()
    sparse.save_npz(PROC / "repair_counts.npz",matrix)
    obs.to_csv(PROC / "repair_obs.csv",index=False)
    save_json(PROC / "repair_genes.json",genes)
    qc = obs.groupby(["unit","state"]).agg(cells=("state","size"), median_umi=("measured_umi","median"),
        median_genes=("measured_genes","median"), median_mito_fraction=("measured_mito_fraction","median"))
    qc.to_csv(TABLES / "repair_state_qc.csv")
    facts, logcpm, detection = pseudobulk(matrix,obs,genes,["AT2","Krt8+ ADI","AT1"])
    assert len(facts) == 27 and (facts.cells >= 30).all()
    facts.to_csv(TABLES / "repair_pseudobulk_inventory.csv",index=False)
    np.savez_compressed(PROC / "repair_pseudobulk.npz", logcpm=logcpm, detection=detection)
    program, evidence = learn_program(facts,logcpm,detection,genes,config)
    evidence.to_csv(TABLES / "repair_gene_evidence.csv",index=False)
    evidence[evidence.selected].sort_values("minimum_median_effect",ascending=False).to_csv(TABLES / "program_genes.csv",index=False)
    folds, fold_genes = [], []
    if len(program) >= config["minimum_program_genes"]:
        scores, coverage = rank_scores(matrix,genes,program,config["reference_genes"])
        paired_effects(obs,scores).to_csv(TABLES / "repair_apparent_effects.csv",index=False)
        for unit in config["repair_units"]:
            selected_genes, _ = learn_program(facts,logcpm,detection,genes,config,omit=unit)
            if len(selected_genes) < config["minimum_program_genes"]:
                continue
            mask = obs.unit.eq(unit).to_numpy()
            heldout, _ = rank_scores(matrix[mask],genes,selected_genes,config["reference_genes"])
            effect = paired_effects(obs.loc[mask],heldout)
            effect["training_program_genes"] = len(selected_genes)
            effect["jaccard_with_full_program"] = len(set(selected_genes)&set(program))/len(set(selected_genes)|set(program))
            folds.append(effect)
            fold_genes += [{"heldout_unit":unit,"gene":g} for g in selected_genes]
        cv = pd.concat(folds,ignore_index=True) if folds else pd.DataFrame()
        cv.to_csv(TABLES / "repair_heldout_effects.csv",index=False)
        pd.DataFrame(fold_genes).to_csv(TABLES / "repair_fold_programs.csv",index=False)
        summaries = {e:{"positive_mice":int((cv.loc[cv.endpoint==e,"difference"] > 0).sum()),
            "median_difference":float(cv.loc[cv.endpoint==e,"difference"].median()),
            "median_standardized_difference":float(cv.loc[cv.endpoint==e,"standardized_difference"].median())} for e in ["AT2","AT1"]}
        passed = len(cv) == 18 and all(v["positive_mice"] >= 6 and v["median_standardized_difference"] >= .25 for v in summaries.values())
        save_json(BASE / "frozen_repair_program.json", {"genes":program,"reference_genes":config["reference_genes"],
            "config_sha256":sha(BASE/"exploratory_config.json"),"training_context":"repair_only",
            "transfer_outcomes_seen":False,"score":config["score"],"coverage_in_repair":coverage})
    else:
        summaries, passed = {}, False
    result = {"stage":"E1", "input_matrix_shape":[len(barcodes),len(genes)], "analyzed_mice":9,
        "cells_in_selected_mouse_libraries":len(obs), "eligible_genes":int(evidence.eligible.sum()),
        "program_genes":len(program), "heldout_results":summaries,
        "checkpoint_passed":bool(passed),"next_stage_recommended":"E2 specificity" if passed else "stop_candidate_transfer",
        "limitations":["Author cell labels reused; internal validation is not independent state validation",
                      "Repair-associated candidate, not a common two-context program", "No causal or universal inference"]}
    save_json(BASE / "stage_E1_result.json",result)
    print(json.dumps(result,indent=2),flush=True)


if __name__ == "__main__":
    main()
