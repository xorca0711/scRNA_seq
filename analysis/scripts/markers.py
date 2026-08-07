"""
Canonical marker panels for lung / airway scRNA-seq, in both mouse and human
symbol conventions.

Symbols are stored explicitly per species.  Nothing here converts a symbol from
one species to the other by changing its case - orthology is not a case
transformation, and silently "converting" symbols is a well known source of
wrong marker calls.  The pipeline looks up whichever panel matches the species
that was detected from the data, then intersects with the genes that are
actually present in the matrix.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Lineage panels.  Order inside each group is deliberate (most specific first)
# because it drives the order of genes on the dot plots.
# ---------------------------------------------------------------------------

MOUSE_PANEL: dict[str, list[str]] = {
    "AT2": ["Sftpc", "Sftpb", "Abca3", "Lamp3", "Slc34a2", "Cxcl15"],
    "AT1": ["Ager", "Pdpn", "Hopx", "Cav1", "Aqp5", "Rtkn2"],
    "Transitional": ["Krt8", "Krt18", "Krt19", "Cldn4", "Krt17", "Cdkn1a", "Sfn"],
    "Basal": ["Krt5", "Krt14", "Trp63", "Dapl1"],
    "Club": ["Scgb1a1", "Scgb3a2", "Scgb3a1", "Cyp2f2"],
    "Ciliated": ["Foxj1", "Tppp3", "Dynlrb2", "Ccdc153"],
    "Neuroendocrine": ["Ascl1", "Calca", "Chga"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr", "Emcn", "Cldn5"],
    "Endothelial_cap": ["Car4", "Ednrb", "Gpihbp1", "Aplnr"],
    "Lymphatic_EC": ["Mmrn1", "Prox1", "Ccl21a"],
    "Fibroblast": ["Col1a1", "Col1a2", "Dcn", "Pdgfra", "Col3a1"],
    "Fibro_activated": ["Cthrc1", "Tnc", "Postn", "Acta2"],
    "SMC_pericyte": ["Acta2", "Myh11", "Tagln", "Pdgfrb", "Notch3"],
    "Mesothelium": ["Msln", "Upk3b", "Wt1"],
    "Alveolar_Mac": ["Marco", "Pparg", "Fabp4", "Siglecf", "Ear1"],
    "Mono_Mac": ["Lyz2", "Cd68", "Csf1r", "Adgre1", "Itgam", "C1qa"],
    "Neutrophil": ["S100a8", "S100a9", "Retnlg"],
    "DC": ["Cd209a", "Itgax", "H2-Ab1", "Flt3"],
    "T_cell": ["Cd3d", "Cd3e", "Cd3g", "Trbc2"],
    "CD4_T": ["Cd4"],
    "CD8_T": ["Cd8a", "Cd8b1"],
    "Treg": ["Foxp3", "Il2ra", "Ctla4", "Ikzf2"],
    "NK": ["Ncr1", "Klrb1c", "Gzma"],
    "B_cell": ["Cd79a", "Cd79b", "Ms4a1", "Cd37", "Cd74"],
    "Plasma": ["Jchain", "Mzb1", "Xbp1"],
    "Proliferating": ["Mki67", "Top2a", "Ccnb1"],
    "Erythroid": ["Hba-a1", "Hbb-bs", "Alas2"],
    "Platelet": ["Ppbp", "Pf4", "Itga2b"],
    "Pan_epithelial": ["Epcam", "Cdh1", "Krt7"],
    "Pan_immune": ["Ptprc"],
}

HUMAN_PANEL: dict[str, list[str]] = {
    "AT2": ["SFTPC", "SFTPB", "ABCA3", "LAMP3", "SLC34A2", "NAPSA"],
    "AT1": ["AGER", "PDPN", "HOPX", "CAV1", "AQP5", "RTKN2"],
    "Transitional": ["KRT8", "KRT18", "KRT19", "CLDN4", "KRT17", "CDKN1A", "SFN"],
    "Basal": ["KRT5", "KRT14", "TP63", "DAPL1"],
    "Club": ["SCGB1A1", "SCGB3A2", "SCGB3A1", "MGP"],
    "Ciliated": ["FOXJ1", "TPPP3", "DNAH5", "CAPS"],
    "Neuroendocrine": ["ASCL1", "CALCA", "CHGA"],
    "Endothelial": ["PECAM1", "CDH5", "KDR", "EMCN", "CLDN5"],
    "Endothelial_cap": ["CA4", "EDNRB", "GPIHBP1", "APLNR"],
    "Lymphatic_EC": ["MMRN1", "PROX1", "CCL21"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "PDGFRA", "COL3A1"],
    "Fibro_activated": ["CTHRC1", "TNC", "POSTN", "ACTA2"],
    "SMC_pericyte": ["ACTA2", "MYH11", "TAGLN", "PDGFRB", "NOTCH3"],
    "Mesothelium": ["MSLN", "UPK3B", "WT1"],
    "Alveolar_Mac": ["MARCO", "PPARG", "FABP4", "MCEMP1"],
    "Mono_Mac": ["LYZ", "CD68", "CSF1R", "ITGAM", "C1QA"],
    "Neutrophil": ["S100A8", "S100A9", "FCGR3B"],
    "DC": ["CD1C", "ITGAX", "HLA-DRA", "FLT3"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRBC2"],
    "CD4_T": ["CD4"],
    "CD8_T": ["CD8A", "CD8B"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "NK": ["NKG7", "KLRD1", "GNLY"],
    "B_cell": ["CD79A", "CD79B", "MS4A1", "CD37", "CD74"],
    "Plasma": ["JCHAIN", "MZB1", "XBP1"],
    "Proliferating": ["MKI67", "TOP2A", "CCNB1"],
    "Erythroid": ["HBB", "HBA1", "ALAS2"],
    "Platelet": ["PPBP", "PF4", "ITGA2B"],
    "Pan_epithelial": ["EPCAM", "CDH1", "KRT7"],
    "Pan_immune": ["PTPRC"],
}

# Genes requested explicitly for feature-plot UMAPs (task phase 22).
MOUSE_FEATURE_GENES = [
    "Sftpc", "Ager", "Krt5", "Trp63", "Krt8", "Krt19", "Krt17", "Scgb1a1",
    "Foxj1", "Cd3d", "Cd4", "Cd8a", "Foxp3", "Lyz2", "Col1a1", "Cthrc1",
    "Pecam1",
]
HUMAN_FEATURE_GENES = [
    "SFTPC", "AGER", "KRT5", "TP63", "KRT8", "KRT19", "KRT17", "SCGB1A1",
    "FOXJ1", "CD3D", "CD4", "CD8A", "FOXP3", "LYZ", "COL1A1", "CTHRC1",
    "PECAM1",
]

# Compact panel used for the epithelial sub-analysis dot plot.
MOUSE_EPI_PANEL: dict[str, list[str]] = {
    "AT2": ["Sftpc", "Sftpb", "Abca3", "Lamp3", "Slc34a2"],
    "AT2_activated": ["Lcn2", "Il33", "Cxcl15"],
    "Transitional_Krt8": ["Krt8", "Cldn4", "Sfn", "Cdkn1a", "Krt18"],
    "Transitional_Krt19": ["Krt19", "Krt7"],
    "Krt17_high": ["Krt17", "Tnc"],
    "AT1": ["Ager", "Pdpn", "Hopx", "Cav1", "Aqp5"],
    "Basal": ["Krt5", "Krt14", "Trp63"],
    "Club": ["Scgb1a1", "Scgb3a2", "Cyp2f2"],
    "Ciliated": ["Foxj1", "Tppp3", "Dynlrb2"],
    "Proliferating": ["Mki67", "Top2a"],
}
HUMAN_EPI_PANEL: dict[str, list[str]] = {
    "AT2": ["SFTPC", "SFTPB", "ABCA3", "LAMP3", "SLC34A2"],
    "Transitional_KRT8": ["KRT8", "CLDN4", "SFN", "CDKN1A", "KRT18"],
    "Transitional_KRT19": ["KRT19", "KRT7"],
    "KRT17_high": ["KRT17", "TNC"],
    "AT1": ["AGER", "PDPN", "HOPX", "CAV1", "AQP5"],
    "Basal": ["KRT5", "KRT14", "TP63"],
    "Club": ["SCGB1A1", "SCGB3A2", "MGP"],
    "Ciliated": ["FOXJ1", "TPPP3", "CAPS"],
    "Proliferating": ["MKI67", "TOP2A"],
}

# Marker sets used to decide, from coherent evidence rather than one gene,
# whether a cluster belongs to the epithelial compartment.
EPITHELIAL_EVIDENCE = {
    "mouse": {
        "positive": ["Epcam", "Cdh1", "Krt8", "Krt18", "Sftpc", "Scgb1a1",
                     "Foxj1", "Ager", "Krt5"],
        "negative": ["Ptprc", "Pecam1", "Col1a1", "Dcn", "Cdh5"],
    },
    "human": {
        "positive": ["EPCAM", "CDH1", "KRT8", "KRT18", "SFTPC", "SCGB1A1",
                     "FOXJ1", "AGER", "KRT5"],
        "negative": ["PTPRC", "PECAM1", "COL1A1", "DCN", "CDH5"],
    },
}

# Gene-name classes that must not, on their own, drive an annotation.
UNINFORMATIVE_PREFIXES = {
    "mouse": {
        "mitochondrial": ("mt-",),
        "ribosomal": ("Rps", "Rpl"),
        "heat_shock": ("Hsp",),
        "hemoglobin": ("Hba", "Hbb"),
    },
    "human": {
        "mitochondrial": ("MT-",),
        "ribosomal": ("RPS", "RPL"),
        "heat_shock": ("HSP",),
        "hemoglobin": ("HBA", "HBB"),
    },
}


def get_panel(species: str) -> dict[str, list[str]]:
    return MOUSE_PANEL if species == "mouse" else HUMAN_PANEL


def get_epi_panel(species: str) -> dict[str, list[str]]:
    return MOUSE_EPI_PANEL if species == "mouse" else HUMAN_EPI_PANEL


def get_feature_genes(species: str) -> list[str]:
    return MOUSE_FEATURE_GENES if species == "mouse" else HUMAN_FEATURE_GENES


def mito_prefix(species: str) -> str:
    return "mt-" if species == "mouse" else "MT-"


def ribo_prefixes(species: str) -> tuple[str, ...]:
    return ("Rps", "Rpl") if species == "mouse" else ("RPS", "RPL")


def hb_prefixes(species: str) -> tuple[str, ...]:
    return ("Hba", "Hbb") if species == "mouse" else ("HBA", "HBB")


def is_uninformative(gene: str, species: str) -> bool:
    """True for mitochondrial / ribosomal / heat-shock / haemoglobin genes."""
    for prefixes in UNINFORMATIVE_PREFIXES[species].values():
        if gene.startswith(tuple(prefixes)):
            return True
    return False


def filter_panel(panel: dict[str, list[str]], available: set[str]
                 ) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Split a panel into the genes present in the matrix and those absent.

    Returns (present_panel, absent_by_group).  Groups that end up empty are
    dropped from the present panel so that no empty row is ever plotted.
    Duplicate genes are removed keeping first occurrence, because scanpy's
    dotplot raises on repeated var names across groups.
    """
    present: dict[str, list[str]] = {}
    absent: dict[str, list[str]] = {}
    seen: set[str] = set()
    for group, genes in panel.items():
        keep, miss = [], []
        for g in genes:
            if g in available:
                if g not in seen:
                    keep.append(g)
                    seen.add(g)
            else:
                miss.append(g)
        if keep:
            present[group] = keep
        if miss:
            absent[group] = miss
    return present, absent
