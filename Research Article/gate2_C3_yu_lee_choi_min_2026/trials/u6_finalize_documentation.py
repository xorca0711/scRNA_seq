"""Update current navigation after the evidence reports exist; preserve original plans."""
from pathlib import Path
import json
import re

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]


def main():
    assert (PAPER / 'EVIDENCE_REVIEW.md').exists()
    h = json.loads((PAPER / 'trials/u5_human_niche/validation.json').read_text())
    t = json.loads((PAPER / 'trials/u5_human_ligand_targets/validation.json').read_text())
    assert h['status'] == t['status'] == 'passed'
    p = PAPER / 'README.md'
    s = p.read_text(encoding='utf-8')
    start = s.index('**Execution status:')
    end = s.index('## Read and review')
    s = s[:start] + '**Execution status: feasible analyses complete; evidence ready for joint review.**\n' + \
        'The completed run includes both IPF cohorts, eligible early mouse niches,\n' + \
        'all 75 human lesion libraries, 56 human spatial sections, nine post-viral\n' + \
        'matrices, and the source-program specificity extensions.\n' + \
        'Read the [evidence review](EVIDENCE_REVIEW.md) and\n' + \
        '[completion register](WORK_PACKAGES.md) for findings, sensitivities and\n' + \
        'questions that remain unidentifiable. No criterion was relaxed after\n' + \
        'viewing results. The [initial report](INITIAL_RUN_REPORT.md) preserves\n' + \
        'the first batch; it is not the current whole-run status.\n\n' + s[end:]
    s = s.replace('| [Analysis pipeline plan]', '| [Completed evidence review](EVIDENCE_REVIEW.md) | Cross-context findings and limits for interpretation review |\n| [Evaluation rules](EVALUATION_RULES.md) | Distinguish null, undercovered, unidentifiable and unmeasured questions |\n| [Reproducibility](REPRODUCIBILITY.md) | Runtime, stage dependencies, provenance and preserved corrections |\n| [Analysis pipeline plan]', 1)
    status_rows = {
        'F01': 'Completed design and biological-unit coverage figure',
        'F02': 'Completed IPF and all-QC human sources; receptor/inhibitor tables retain assay limits',
        'F03': 'Completed IPF, mouse and paired human RNA compatibility with native LR sensitivities',
        'F04': 'Completed IPF, mouse and paired human pathways; eligible IPF/human ligand targets',
        'F05': 'Early niche measurements complete; KAC identity and primary phenotype remain unidentifiable',
        'F06': 'All human libraries and spatial/context matrices processed; independent ROI and post-viral coordinates unavailable',
        'F07': 'Completed specificity panels and evidence matrix; interpretation review ready',
    }
    lines = s.splitlines()
    for i, line in enumerate(lines):
        for code, status in status_rows.items():
            if line.startswith('| ' + code + ' '):
                first = line.split('|')[1].strip()
                lines[i] = f'| {first} | {status} |'
    s = '\n'.join(lines) + '\n'
    marker = '### F02/F03: full-cell LR robustness and IL1B sources'
    s = s.replace(marker, '### F01: context, design and biological-unit coverage\n\n![Context and coverage](figures/context_design_and_coverage.png)\n\nCounts preserve animals, donors and paired patients; companion assays and\nrepeated histologies do not add independent replication.\n[Design table](trials/u6_completion/context_design_map.csv) ·\n[generating script](trials/u6_report_context_coverage.py).\n\n' + marker)
    gallery = '''### F02 extension: human IL1B sources

![Human IL1B sources](figures/human_IL1B_sources.png)

All QC cells remain in the source scan. Unassigned fine labels carry median
IL1B count fractions of 52–72% across histologies. Confident-subtype results
therefore cannot establish the dominant source across all recovered cells.
[Report and tables](trials/u5_human_sources/REPORT.md) ·
[generating script](trials/u5_report_human_sources.py).

### F04/F06: paired human pathways and RNA niches

![Human paired recipient pathways](figures/human_paired_recipient_pathways.png)

![Human paired niche RNA](figures/human_paired_niche_RNA_contrasts.png)

Individual patients remain the unit of comparison. Filled pathway symbols
indicate the declared primary global q<0.05 criterion; RNA compatibility
points are descriptive and do not establish communication.
[Report and all sensitivities](trials/u5_human_niche/REPORT.md) ·
[annotation review](trials/u5_human_full/ANNOTATION_REVIEW.md) ·
[generating script](trials/u5_report_human_niche.py).

### F04 extension: paired human ligand targets

![Human ligand target eligibility and fit](figures/human_ligand_target_eligibility_and_fit.png)

Target eligibility precedes source/receiver expression filtering and ranking.
An unsigned-prior fit to downregulated targets does not establish inhibition.
[Report and omission stability](trials/u5_human_ligand_targets/REPORT.md) ·
[generating script](trials/u5_report_human_ligand_targets.py).

### F07 component: human epithelial programs

![Human epithelial program specificity](figures/human_epithelial_program_specificity.png)

Paired source-program changes in reference-compatible AT2-like populations
remain distinct from source-defined KAC/HPCS identity and malignant status.
[Report and sensitivities](trials/u6_human_specificity/REPORT.md) ·
[generating script](trials/u6_report_human_specificity.py).

### F07: cross-context evidence for review

![Cross-context evidence review](figures/cross_context_evidence_review.png)

[Evidence review](EVIDENCE_REVIEW.md) ·
[complete table](trials/u6_completion/evidence_matrix.csv) ·
[generating script](trials/u6_summarize_evidence.py).

'''
    s = s.replace('## Layout and execution boundary', gallery + '## Layout and execution boundary')
    s = s.replace('U0/U1 are preparation utilities. U2/U3 and the initial IPF pathway arm now\nhave executable stages and run records. Later biological stages are released\nonly when their design, annotation and resource specifications are satisfied;\nlaunch authorization is recorded in the analysis contracts.', 'The [reproduction guide](REPRODUCIBILITY.md) indexes the completed stages.\nThe [release checks](trials/u6_completion/release_validation.json) bind\nvalidated stages, figure reviews and released file hashes. Unidentifiable\nendpoints remain explicit in the analysis contracts and evidence register.')
    p.write_text(s, encoding='utf-8')
    p = PAPER / 'WORK_PACKAGES.md'; s = p.read_text(encoding='utf-8')
    replacements = {
        'All 75 library inputs and 1,214 native LR calls complete; all-QC source profiles complete; paired aggregation/models running': f'Completed; 75 libraries, 1,214 native LR calls, {h["primary_RNA_contrasts"]} primary paired RNA contrasts and {h["primary_pathway_tests"]} primary pathway tests',
        'Completed for eligible IPF receivers; 274 candidate results and full donor-omission checks; paired human refits pending': f'Completed for eligible IPF and human receivers; 274 IPF and {t["candidate_rows"]} human candidate rows with omission checks',
        'HPCS/repair, source ISR and IPF state comparisons complete; paired human extension pending': 'Completed HPCS/repair, source ISR, IPF state and paired human extensions; original Han expression arm remains conditional on missing matrices',
        'Completed package figures released; human figures and final synthesis follow paired models': 'Completed gallery and evidence register; ready for interpretation review',
    }
    for old, new in replacements.items():
        assert old in s, old
        s = s.replace(old, new)
    s += f'''\n## Completed: paired human niches, targets and specificity (U5/U6)

All primary and declared sensitivity stages completed. The primary pathway
family contains {h['primary_pathway_tests']} tests with {h['primary_pathway_q05']} global q<0.05 discoveries.
The {h['primary_RNA_contrasts']} primary RNA contrasts retain individual paired patients.
Conditional target fits yield {t['candidate_rows']} expression-supported candidate rows
across overlapping source scopes and receiver views. These are not independent
replications or causal activation tests. Source-program contrasts remain
separate from the unavailable KAC/NF-κB endpoint.
[Niche report](trials/u5_human_niche/REPORT.md) ·
[target report](trials/u5_human_ligand_targets/REPORT.md) ·
[program report](trials/u6_human_specificity/REPORT.md).

## Completed: gallery and evidence register

The [paper gallery](README.md#figure-gallery) includes context, source,
compatibility, pathway, perturbation, human/spatial and specificity panels.
The root README links to paper galleries without selective figure embeds.
The [evidence review](EVIDENCE_REVIEW.md) separates supported descriptive
findings, tested results without primary discoveries, and questions that
the public data cannot identify. The original criteria and historical claims
are unchanged. [Reproduction and provenance](REPRODUCIBILITY.md).
'''
    p.write_text(s, encoding='utf-8')
    p = PAPER / 'FIGURE_GALLERY_PLAN.md'; s = p.read_text(encoding='utf-8')
    s = s.replace('## Reading order and current readiness', '## Reading order and current readiness\n\nThe completed gallery and evidence register are in [README](README.md#figure-gallery)\nand [EVIDENCE_REVIEW](EVIDENCE_REVIEW.md). The panel specifications below\npreserve the original presentation plan; unavailable endpoints are explicitly\nrepresented by eligibility tables or narrower measured panels.')
    lines = s.splitlines()
    for i, line in enumerate(lines):
        for code, status in status_rows.items():
            if line.startswith('| ' + code + ' |'):
                parts = line.split('|'); parts[3] = ' ' + status + ' '; lines[i] = '|'.join(parts)
    p.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    for name in ['analysis_contract.json', 'niche_analysis_contract.json']:
        p = PAPER / name; x = json.loads(p.read_text()); x['status'] = 'completed_feasible_scope_with_documented_scientific_gates'
        x['completion_report'] = 'EVIDENCE_REVIEW.md'; x['completion_register'] = 'WORK_PACKAGES.md'
        if name == 'analysis_contract.json':
            for stage in x.get('stages', []):
                if stage['id'] == 'U5': stage['status'] = 'IPF_human_and_spatial_measurements_complete_identity_and_region_endpoints_unidentifiable'
                if stage['id'] == 'U6': stage['status'] = 'specificity_and_evidence_synthesis_complete_original_Han_expression_arm_unavailable'
        else:
            x['computed_scope'] = 'IPF_mouse_human_niche_measurements_native_LR_sensitivities_eligible_ligand_targets_spatial_whole_section_context_and_epithelial_specificity'
            x['ligand_receptor']['exact_edges_and_complex_subunits_status'] = 'frozen_resources_and_biological_role_tables_available_for_IPF_mouse_and_human'
            x['pathways']['exact_gene_sets_versions_hashes_status'] = 'all_contexts_frozen_and_coverage_sensitivity_verified'
            x['human_completion'] = dict(report='trials/u5_human_niche/REPORT.md', pathways=h['primary_pathway_tests'], primary_q05=h['primary_pathway_q05'], native_calls=1214, ligand_target_candidates=t['candidate_rows'])
        p.write_text(json.dumps(x, indent=2) + '\n', encoding='utf-8')
    p = PAPER / 'DATASETS.md'; s = p.read_text(encoding='utf-8')
    s = s.replace('## New verified deposits', '## Current execution\n\nThe [completion register](WORK_PACKAGES.md) records the completed public-data\nanalyses. All 75 human RNA libraries, 56 human spatial sections and nine\npost-viral matrices were processed. Both IPF cohorts, eligible early mouse\nniches and specificity extensions have measured outputs. The shortlist below\npreserves source-design questions; endpoint-level resolutions and remaining\nlimits are in the [evidence review](EVIDENCE_REVIEW.md).\n\n## New verified deposits', 1)
    s = s.replace('For GSE300288 and GSE308103, actual eligible cell coverage in each compartment\nis **not yet verified** by the two technical matrix checks. Require same-sample', 'For GSE300288 and GSE308103, actual eligible cell coverage is now measured\nin the [completed context figure](README.md#figure-gallery) and biological-unit\ntables. Comparisons require same-sample')
    p.write_text(s, encoding='utf-8')
    p = ROOT / 'README.md'; s = p.read_text(encoding='utf-8')
    s = s.replace('Initial IPF pathway results and the context-specific gallery plan', 'IPF, mouse and paired human niches; spatial context, specificity and evidence review')
    start = s.index('The [IL-1beta review branch]'); end = s.index('Deposits used in the existing biological analyses', start)
    s = s[:start] + 'The [IL-1beta review branch](Research Article/gate2_C3_yu_lee_choi_min_2026/README.md)\nhas completed its feasible public-data analyses across IPF, early mouse\nblockade, human lesions, spatial context and epithelial program specificity.\nIts [evidence review](Research Article/gate2_C3_yu_lee_choi_min_2026/EVIDENCE_REVIEW.md)\ndistinguishes measured findings from inconclusive and unidentifiable endpoints.\nThe [public-data shortlist](Research Article/gate2_C3_yu_lee_choi_min_2026/DATASETS.md)\nrecords dataset roles; the table below describes the earlier analyses.\n\n' + s[end:]
    p.write_text(s, encoding='utf-8')
    print('Updated current documentation and gallery navigation; original plans and criteria retained.')


if __name__ == '__main__':
    main()
