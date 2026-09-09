#!/usr/bin/env python
"""Build the thesis-aware portfolio PDF from tracked analysis artefacts."""

from __future__ import annotations

from datetime import date
from pathlib import Path

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "pdf" / "lung_scrna_portfolio_thesis_context.pdf"

TEAL = colors.HexColor("#087F8C")
BLUE = colors.HexColor("#1D4E89")
NAVY = colors.HexColor("#16324F")
GREEN = colors.HexColor("#4C956C")
GOLD = colors.HexColor("#E6A817")
CORAL = colors.HexColor("#D95D4F")
INK = colors.HexColor("#24313B")
MUTED = colors.HexColor("#667580")
PALE = colors.HexColor("#F2F7F8")
LIGHT_BLUE = colors.HexColor("#EDF3F9")
WHITE = colors.white

PAGE_W, PAGE_H = landscape(A4)
MARGIN_X = 15 * mm
MARGIN_TOP = 15 * mm
MARGIN_BOTTOM = 13 * mm


def p(path: str) -> Path:
    result = ROOT / path
    if not result.exists():
        raise FileNotFoundError(result)
    return result


def image(path: str, max_w: float, max_h: float) -> Image:
    source = p(path)
    with PILImage.open(source) as raster:
        width, height = raster.size
    scale = min(max_w / width, max_h / height)
    return Image(str(source), width=width * scale, height=height * scale)


class PortfolioDoc(BaseDocTemplate):
    def __init__(self, filename: Path):
        super().__init__(
            str(filename), pagesize=(PAGE_W, PAGE_H),
            leftMargin=MARGIN_X, rightMargin=MARGIN_X,
            topMargin=MARGIN_TOP, bottomMargin=MARGIN_BOTTOM,
            title="Independent lung scRNA-seq reanalysis in thesis context",
            author="Independent Python/Scanpy reanalysis",
            subject="Human distal lung epithelial states, AT0, and additional lung regeneration results",
        )
        frame = Frame(
            self.leftMargin, self.bottomMargin,
            self.width, self.height, id="body",
        )
        self.addPageTemplates(PageTemplate(id="content", frames=[frame], onPage=draw_page))


def draw_page(canvas, doc) -> None:
    canvas.saveState()
    if doc.page == 1:
        canvas.setFillColor(NAVY)
        canvas.rect(0, 0, PAGE_W, PAGE_H, fill=1, stroke=0)
    else:
        canvas.setStrokeColor(colors.HexColor("#D7E2E5"))
        canvas.line(MARGIN_X, PAGE_H - 10 * mm, PAGE_W - MARGIN_X, PAGE_H - 10 * mm)
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(MARGIN_X, PAGE_H - 7.4 * mm, "Lung scRNA-seq portfolio | thesis-aware synthesis")
        canvas.drawRightString(PAGE_W - MARGIN_X, 7 * mm, f"{doc.page}")
    canvas.restoreState()


styles = getSampleStyleSheet()
styles.add(ParagraphStyle(
    name="CoverTitle", parent=styles["Title"], fontName="Helvetica-Bold",
    fontSize=28, leading=32, textColor=WHITE, alignment=TA_LEFT, spaceAfter=10,
))
styles.add(ParagraphStyle(
    name="CoverSub", parent=styles["BodyText"], fontName="Helvetica",
    fontSize=13, leading=18, textColor=colors.HexColor("#DDEFF1"),
))
styles.add(ParagraphStyle(
    name="Section", parent=styles["Heading1"], fontName="Helvetica-Bold",
    fontSize=20, leading=24, textColor=NAVY, spaceAfter=7, keepWithNext=True,
))
styles.add(ParagraphStyle(
    name="Subsection", parent=styles["Heading2"], fontName="Helvetica-Bold",
    fontSize=12.5, leading=15, textColor=TEAL, spaceBefore=4, spaceAfter=4,
    keepWithNext=True,
))
styles.add(ParagraphStyle(
    name="BodySmall", parent=styles["BodyText"], fontName="Helvetica",
    fontSize=8.7, leading=12, textColor=INK, spaceAfter=5,
))
styles.add(ParagraphStyle(
    name="Body", parent=styles["BodyText"], fontName="Helvetica",
    fontSize=9.6, leading=13.5, textColor=INK, spaceAfter=6,
))
styles.add(ParagraphStyle(
    name="Caption", parent=styles["BodyText"], fontName="Helvetica",
    fontSize=7.4, leading=10, textColor=MUTED, alignment=TA_LEFT, spaceBefore=3,
))
styles.add(ParagraphStyle(
    name="Callout", parent=styles["BodyText"], fontName="Helvetica-Bold",
    fontSize=11, leading=15, textColor=NAVY, alignment=TA_LEFT,
))
styles.add(ParagraphStyle(
    name="Metric", parent=styles["BodyText"], fontName="Helvetica-Bold",
    fontSize=18, leading=20, textColor=TEAL, alignment=TA_CENTER,
))
styles.add(ParagraphStyle(
    name="MetricLabel", parent=styles["BodyText"], fontName="Helvetica",
    fontSize=7.8, leading=10, textColor=MUTED, alignment=TA_CENTER,
))
styles.add(ParagraphStyle(
    name="BulletSmall", parent=styles["BodyText"], fontName="Helvetica",
    fontSize=8.6, leading=11.6, textColor=INK, leftIndent=10, firstLineIndent=-6,
    bulletIndent=0, spaceAfter=3,
))


def para(text: str, style: str = "Body") -> Paragraph:
    return Paragraph(text, styles[style])


def bullet(text: str) -> Paragraph:
    return Paragraph(f"• {text}", styles["BulletSmall"])


def section(title: str, kicker: str | None = None) -> list:
    items = [para(title, "Section")]
    if kicker:
        items.append(para(kicker, "Body"))
    return items


def callout(text: str, background=PALE, accent=TEAL) -> Table:
    table = Table([[para(text, "Callout")]], colWidths=[PAGE_W - 2 * MARGIN_X])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), background),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor("#C9DADC")),
        ("LINEBEFORE", (0, 0), (0, -1), 4, accent),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 10),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
    ]))
    return table


def metrics(items: list[tuple[str, str]], widths: list[float] | None = None) -> Table:
    widths = widths or [(PAGE_W - 2 * MARGIN_X) / len(items)] * len(items)
    cells = [[para(value, "Metric") for value, _ in items],
             [para(label, "MetricLabel") for _, label in items]]
    table = Table(cells, colWidths=widths)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), PALE),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor("#D6E2E5")),
        ("INNERGRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#D6E2E5")),
        ("TOPPADDING", (0, 0), (-1, 0), 7),
        ("BOTTOMPADDING", (0, 1), (-1, 1), 7),
    ]))
    return table


def figure_block(path: str, caption: str, max_w: float, max_h: float) -> KeepTogether:
    return KeepTogether([image(path, max_w, max_h), para(caption, "Caption")])


def figure_items(path: str, caption: str, max_w: float, max_h: float) -> list:
    """Image/caption pair safe for use inside a table cell."""
    return [image(path, max_w, max_h), para(caption, "Caption")]


def two_col(left: list, right: list, left_w: float | None = None) -> Table:
    gap = 8 * mm
    usable = PAGE_W - 2 * MARGIN_X - gap
    left_w = left_w or usable / 2
    right_w = usable - left_w
    table = Table([[left, right]], colWidths=[left_w, right_w], hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (0, -1), gap),
        ("RIGHTPADDING", (1, 0), (1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
    ]))
    return table


def build_story() -> list:
    story: list = []

    # Cover
    cover = Table([
        [para("INDEPENDENT REANALYSIS", "CoverSub")],
        [para("Human distal lung epithelial states<br/>in the context of AT0 biology", "CoverTitle")],
        [para(
            "A portfolio synthesis of GSE178360 with complementary mouse lung-regeneration results from GSE262927",
            "CoverSub",
        )],
        [Spacer(1, 10 * mm)],
        [para(
            "Raw-count reconstruction • candidate identities • reference-oriented epithelial map • "
            "KRT8 / CLDN4 / KRT17 / SFN projections • cross-species regeneration context",
            "CoverSub",
        )],
        [Spacer(1, 14 * mm)],
        [para(f"Portfolio edition | {date.today().isoformat()}", "CoverSub")],
    ], colWidths=[PAGE_W - 2 * MARGIN_X], rowHeights=[None] * 7)
    cover.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), NAVY),
        ("LEFTPADDING", (0, 0), (-1, -1), 18 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 18 * mm),
        ("TOPPADDING", (0, 0), (-1, 0), 16 * mm),
        ("BOTTOMPADDING", (0, -1), (-1, -1), 20 * mm),
    ]))
    story += [cover, PageBreak()]

    # Executive summary
    story += section(
        "Executive summary",
        "The thesis paper asks how human distal-airway epithelial states connect to alveolar lineages, with AT0 proposed as a bipotent progenitor-like population. This report tests how much of that landscape is recoverable from deposited raw counts without importing the authors' processed objects or labels.",
    )
    story.append(metrics([
        ("27,729", "human cells after QC + doublet removal"),
        ("31", "human whole-atlas Leiden clusters"),
        ("6,386", "cells in the focused epithelial display"),
        ("14", "annotated epithelial Leiden regions"),
        ("165.51°", "presentation-only CCW rotation"),
    ]))
    story.append(Spacer(1, 5 * mm))
    left = [
        para("What was recovered", "Subsection"),
        bullet("A distinct <b>SFTPC+SCGB3A2+ AT0 candidate analogue</b> (epithelial Leiden 4), separated from a conventional SFTPC-high AT2 candidate (Leiden 0)."),
        bullet("Basal, differentiating-basal, secretory/TRB-like, ciliated, neuroendocrine, AT2 and AT1 regions supported by coherent marker expression."),
        bullet("A KRT8/CLDN4/KRT17/SFN programme concentrated in basal/transitional regions, with KRT17 and SFN showing the clearest restriction."),
        bullet("The complementary mouse time course independently recovers an AT2-to-Krt8-transitional-to-AT1 progression and its temporal resolution."),
    ]
    right = [
        para("What is not claimed", "Subsection"),
        bullet("The candidate analogues are <b>not transferred author annotations</b> and do not reproduce the paper's 18 clusters one-to-one."),
        bullet("A healthy cross-sectional human dataset does not prove lineage direction. The AT0 interpretation rests on marker convergence and correspondence, not trajectory causality."),
        bullet("The display rotation changes orientation only. It does not alter neighbours, cluster membership, distances or the internal UMAP geometry."),
        bullet("Eight Fig. 1c states are not separately resolved and are listed as unresolved rather than assigned by appearance."),
    ]
    story.append(two_col(left, right))
    story.append(Spacer(1, 4 * mm))
    story.append(callout(
        "Bottom line: the raw-count reanalysis supports an AT0-like human epithelial population and a broader distal-airway-to-alveolar state continuum, while keeping the central distinction between marker-supported correspondence and proven lineage hierarchy.",
        LIGHT_BLUE, BLUE,
    ))
    story.append(PageBreak())

    # Context and design
    story += section("Study context and analytical design")
    context_data = [
        [para("Dataset", "BodySmall"), para("Design", "BodySmall"), para("Role in this report", "BodySmall")],
        [para("GSE178360", "BodySmall"), para("Human healthy distal lung; 3 donors; 10x 3'", "BodySmall"), para("Thesis-focused whole-atlas and epithelial-state reconstruction", "BodySmall")],
        [para("GSE262927", "BodySmall"), para("Mouse H1N1 injury time course; 33 samples", "BodySmall"), para("Dynamic regeneration and vascular-injury context", "BodySmall")],
    ]
    t = Table(context_data, colWidths=[42 * mm, 88 * mm, 128 * mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("TEXTCOLOR", (0, 0), (-1, 0), WHITE),
        ("BACKGROUND", (0, 1), (-1, -1), PALE),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#C9D5D9")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    story += [t, Spacer(1, 5 * mm)]
    design_left = [
        para("Human pipeline", "Subsection"),
        bullet("36,464 barcodes loaded; 29,605 retained after per-sample QC; 27,729 after 1,876 Scrublet calls."),
        bullet("Samples intersected on Ensembl gene ID because donor DD073R used a different GRCh38 annotation build."),
        bullet("log1p(CP10K) expression; 2,500 HVGs; 50 PCs; Leiden resolution 1.0."),
        bullet("Harmony on donor was an explicit override because same-type cells fragmented into donor-private islands; the unintegrated view is retained."),
    ]
    design_right = [
        para("Independence safeguards", "Subsection"),
        bullet("The authors' processed Seurat objects were not used for model fitting or coordinate transfer."),
        bullet("Numeric Leiden labels remain primary; every biological identity is suffixed or described as a candidate."),
        bullet("Reference names are used only when defining marker patterns are observed."),
        bullet("Expression panels use the identical epithelial subset and identical rigid display transform."),
    ]
    story.append(two_col(design_left, design_right))
    story.append(Spacer(1, 4 * mm))
    story.append(callout(
        "Interpretive hierarchy: measured expression and cluster membership > candidate identity > closest thesis analogue. Visual resemblance alone is never used as annotation evidence.",
        PALE, TEAL,
    ))
    story.append(PageBreak())

    # Whole human atlas
    story += section(
        "Human whole-atlas result: integration creates a usable shared donor space",
        "The three human samples are healthy biological replicates. Harmony was selected to reduce donor-private fragmentation while retaining the unintegrated embedding for audit.",
    )
    story.append(figure_block(
        "analysis/GSE178360/figures/umap/UMAP_integration_before_after.png",
        "Before/after integration comparison. The decision is dataset-specific: correcting donor in this healthy replicate design is defensible, unlike correcting sample in the confounded mouse injury time course.",
        PAGE_W - 2 * MARGIN_X, 116 * mm,
    ))
    story.append(Spacer(1, 3 * mm))
    story.append(metrics([
        ("4 / 31", "clusters >75% single-donor after Harmony"),
        ("false", "automated donor-driven clustering flag"),
        ("1.0", "selected whole-atlas Leiden resolution"),
        ("0", "clusters below 20 cells at selected resolution"),
    ]))
    story.append(PageBreak())

    # Epithelial map
    story += section(
        "Thesis-focused epithelial landscape",
        "The epithelial object contains 7,446 cells. The reference-oriented figure retains 6,386 epithelial candidates and excludes 1,060 immune, endothelial, plasma and mesothelial carryover cells.",
    )
    story.append(figure_block(
        "analysis/GSE178360/epithelial_subanalysis/figures/reference_aligned/epithelial_UMAP_proposed_reference_aligned.png",
        "All 14 retained Leiden regions are annotated. Literal Fig. 1c names are used only where their marker combination is supported. Eight source-paper states that are not separately resolved are listed on the figure.",
        PAGE_W - 2 * MARGIN_X, 137 * mm,
    ))
    story.append(para(
        "Orientation note: the underlying Scanpy embedding is mean-centred and rotated 165.51 degrees counter-clockwise for display. There is no reflection, scaling, nonlinear warping, neighbour recalculation or UMAP rerun.",
        "Caption",
    ))
    story.append(PageBreak())

    # AT0 interpretation
    story += section("AT0 as the central thesis-facing result")
    at0_left = [
        para("Why epithelial Leiden 4 is separated", "Subsection"),
        bullet("<b>SFTPC mean 4.41</b>; expressed in 88% of cluster cells."),
        bullet("<b>SCGB3A2 mean 3.49</b>; expressed in 95% of cluster cells."),
        bullet("SFTPB mean 5.13; expressed in 100% of cluster cells."),
        bullet("Its location is between secretory and conventional alveolar regions in the unchanged embedding."),
        para("Together these observations support the label <b>SFTPC+SCGB3A2+ (AT0) candidate analogue</b>.", "Body"),
    ]
    at0_right = [
        para("Why epithelial Leiden 0 remains AT2", "Subsection"),
        bullet("<b>SFTPC mean 6.75</b>; expressed in 100% of cluster cells."),
        bullet("SCGB3A2 is substantially lower: mean 0.71; expressed in 70% of cells."),
        bullet("NAPSA, LAMP3, SFTPA1/2, SFTPD and ABCA3 are among the strongest cluster markers."),
        para("This is the more conventional mature AT2-like state and should not be collapsed with the AT0 candidate.", "Body"),
    ]
    story.append(two_col(at0_left, at0_right))
    story.append(Spacer(1, 5 * mm))
    story.append(metrics([
        ("328", "cells in AT0 candidate (Leiden 4)"),
        ("354", "cells in conventional AT2 candidate (Leiden 0)"),
        ("159", "cells in AT1 candidate (Leiden 18)"),
        ("8", "reference states not separately resolved"),
    ]))
    story.append(Spacer(1, 6 * mm))
    story.append(callout(
        "Evidence grade: strong marker-supported correspondence, but not identity transfer. The healthy human dataset supplies cell-state structure; lineage direction remains a thesis-derived model requiring perturbation, time, or lineage tracing for direct proof.",
        colors.HexColor("#FFF7E3"), GOLD,
    ))
    story.append(Spacer(1, 6 * mm))
    unresolved = [
        "Deuterosomal cells", "Proliferating airway cells", "SFTPB+KRT5- Distal-BC-2",
        "SCGB1A1+MUC5B+", "MUC5AC+MUC5B+", "SFTPB+SCGB3A2+SCGB1A1- TRB-SC",
        "Proliferating AT2", "Immature AT1",
    ]
    story.append(para("States explicitly left unresolved", "Subsection"))
    unresolved_table = Table(
        [[para(unresolved[i], "BodySmall"), para(unresolved[i + 1], "BodySmall")]
         for i in range(0, len(unresolved), 2)],
        colWidths=[(PAGE_W - 2 * MARGIN_X) / 2] * 2,
    )
    unresolved_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), PALE),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#D1DEE1")),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(unresolved_table)
    story.append(PageBreak())

    # Primary feature plots
    story += section(
        "Primary transitional-state marker projections",
        "KRT8, CLDN4, KRT17 and SFN are projected on the same 6,386 cells and the same presentation transform. Colour limits are per-gene 99th percentiles among positive cells, so intensity should be compared spatially within a panel, not numerically between genes.",
    )
    story.append(figure_block(
        "analysis/GSE178360/epithelial_subanalysis/figures/reference_aligned/epithelial_featureplots_KRT8_CLDN4_KRT17_SFN_reference_aligned.png",
        "KRT8 and CLDN4 are broad epithelial/state markers in this object. KRT17 and SFN provide more selective support for basal/transitional regions. AT0, AT2 and AT1 are labelled separately to prevent the thesis-focused alveolar states from being visually merged.",
        PAGE_W - 2 * MARGIN_X, 136 * mm,
    ))
    story.append(PageBreak())

    # Reference markers
    story += section(
        "Marker triangulation across the epithelial landscape",
        "No single marker defines the inferred states. The interpretation relies on complementary airway, secretory, basal and alveolar programs.",
    )
    story.append(figure_block(
        "analysis/GSE178360/epithelial_subanalysis/figures/reference_aligned/epithelial_reference_marker_featureplots_reference_aligned.png",
        "SFTPC distinguishes alveolar regions; SCGB3A2 and SFTPB bridge secretory and AT0-like regions; SCGB1A1 marks secretory cells; KRT5/TP63 mark basal states; FOXJ1 marks ciliated cells; AGER marks AT1; KRT19 is broad in epithelial transition/basal regions.",
        PAGE_W - 2 * MARGIN_X, 136 * mm,
    ))
    story.append(PageBreak())

    # Mouse dynamic result
    story += section(
        "Additional result: dynamic evidence from post-influenza regeneration",
        "The mouse dataset addresses what the cross-sectional human dataset cannot: temporal progression. Author labels were held out of trajectory construction and used only after fitting as an external check.",
    )
    left_fig = figure_items(
        "analysis/GSE262927/regeneration_focus/figures/UMAP_alveolar_pseudotime.png",
        "Diffusion pseudotime on 5,694 alveolar epithelial cells, rooted in AT2.",
        122 * mm, 112 * mm,
    )
    right_fig = figure_items(
        "analysis/GSE262927/regeneration_focus/figures/transitional_timecourse.png",
        "The transitional population peaks after injury and nearly resolves by one year.",
        122 * mm, 112 * mm,
    )
    story.append(two_col(left_fig, right_fig, left_w=126 * mm))
    story.append(Spacer(1, 4 * mm))
    story.append(metrics([
        ("0.013", "median AT2 pseudotime"),
        ("0.179", "median transitional pseudotime"),
        ("0.327", "median AT1 pseudotime"),
        ("27.4%", "median transitional fraction at 11 dpi"),
        ("0.3%", "median transitional fraction at 366 dpi"),
    ]))
    story.append(PageBreak())

    # iCAP and lineage
    story += section(
        "Additional result: epithelial repair resolves, vascular injury persists",
        "The capillary analysis provides a biologically informative contrast: the epithelial transitional programme contracts, while the injury-associated capillary state remains elevated long after infection.",
    )
    left_fig = figure_items(
        "analysis/GSE262927/regeneration_focus/figures/icap_timecourse.png",
        "iCAP is nearly absent in homeostasis, surges after injury and remains elevated at one year.",
        122 * mm, 108 * mm,
    )
    right_fig = figure_items(
        "analysis/GSE262927/lineage_tracing_cohort/figures/icap_origin_by_cre_line.png",
        "Reporter-positive fractions support a CAP1 origin in the informative Kit line; CAP2 lines are under-labelled and therefore uninformative rather than negative.",
        122 * mm, 108 * mm,
    )
    story.append(two_col(left_fig, right_fig, left_w=126 * mm))
    story.append(Spacer(1, 4 * mm))
    story.append(metrics([
        ("2.0%", "median iCAP fraction in homeostasis"),
        ("37.5%", "median iCAP fraction at 25 dpi"),
        ("21.7%", "median iCAP fraction at 366 dpi"),
        ("33±3%", "Kit-line iCAP tracing per animal"),
    ]))
    story.append(PageBreak())

    # Validation and limits
    story += section("Validation, self-audits and limitations")
    qa_left = [
        para("What strengthens confidence", "Subsection"),
        bullet("All results are regenerated from raw deposited counts with deterministic seeds and machine-logged decisions."),
        bullet("All 31 selected human clusters satisfy the minimum size and marker-support criteria at resolution 1.0."),
        bullet("The portfolio validator passes 165 checks: Markdown links, JSON artefacts, headline counts and biological results."),
        bullet("Off-compartment marker projections are retained as negative controls rather than omitted."),
        bullet("The stricter AT0 doublet audit found 3.9% AT0 flagged versus a 6.3% baseline, overturning an earlier misleading co-expression gate."),
    ]
    qa_right = [
        para("What constrains interpretation", "Subsection"),
        bullet("The human embedding is independent of the authors' Seurat UMAP; exact cluster geometry cannot be recovered without their original coordinates and preprocessing state."),
        bullet("Two of three human Scrublet thresholds are expected-rate ranking cuts because the score histograms lacked a trustworthy automatic separation."),
        bullet("No ambient-RNA correction was possible from filtered matrices alone; soup-prone genes such as SFTPC require contextual interpretation."),
        bullet("Some epithelial subclusters are donor-skewed or low-depth and should be validated in an independent cohort."),
        bullet("The mouse and human analyses are complementary, not a pooled cross-species model."),
    ]
    story.append(two_col(qa_left, qa_right))
    story.append(Spacer(1, 6 * mm))
    story.append(callout(
        "The principal portfolio value is not perfect visual replication. It is a reproducible, inspectable reconstruction that distinguishes direct measurements, candidate correspondences, external validation and unresolved biology.",
        LIGHT_BLUE, BLUE,
    ))
    story.append(Spacer(1, 7 * mm))
    story.append(para("Confidence map", "Subsection"))
    confidence = [
        [para("Claim", "BodySmall"), para("Evidence", "BodySmall"), para("Confidence", "BodySmall")],
        [para("Distinct human AT0 candidate analogue", "BodySmall"), para("SFTPC/SCGB3A2/SFTPB co-expression + separate cluster", "BodySmall"), para("Moderate-high", "BodySmall")],
        [para("Human lineage direction", "BodySmall"), para("Cross-sectional marker structure only", "BodySmall"), para("Not established", "BodySmall")],
        [para("Mouse AT2-transitional-AT1 ordering", "BodySmall"), para("Pseudotime + held-out labels + injury time course", "BodySmall"), para("High within dataset", "BodySmall")],
        [para("Persistent mouse iCAP state", "BodySmall"), para("Per-animal longitudinal composition", "BodySmall"), para("High within dataset", "BodySmall")],
    ]
    ct = Table(confidence, colWidths=[72 * mm, 135 * mm, 51 * mm])
    ct.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("TEXTCOLOR", (0, 0), (-1, 0), WHITE),
        ("BACKGROUND", (0, 1), (-1, -1), PALE),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#CAD8DB")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    story.append(ct)
    story.append(PageBreak())

    # Conclusion/references
    story += section("Conclusion and reproducibility")
    story.append(callout(
        "The independent analysis recovers a coherent human distal-lung epithelial landscape with a distinct SFTPC+SCGB3A2+ AT0 candidate analogue, while the complementary injury time course demonstrates a transient Krt8-associated epithelial repair route. Together they support the thesis focus on intermediate distal-lung epithelial states without treating candidate correspondence as lineage proof.",
        colors.HexColor("#EAF6F0"), GREEN,
    ))
    story.append(Spacer(1, 6 * mm))
    repro_left = [
        para("Reproduce the report", "Subsection"),
        para("From the repository root:", "BodySmall"),
        para("<font name='Courier'>python analysis/scripts/09_write_portfolio_pdf.py</font>", "BodySmall"),
        para("Regenerate the reference-oriented human figures:", "BodySmall"),
        para("<font name='Courier'>python analysis/scripts/08_reference_aligned_epithelial_umap.py</font>", "BodySmall"),
        para("Validate tracked artefacts:", "BodySmall"),
        para("<font name='Courier'>python analysis/scripts/validate_portfolio.py</font>", "BodySmall"),
    ]
    repro_right = [
        para("Key artefacts", "Subsection"),
        bullet("analysis/GSE178360/epithelial_subanalysis/figures/reference_aligned/"),
        bullet("analysis/GSE178360/README.md"),
        bullet("analysis/GSE262927/regeneration_focus/"),
        bullet("analysis/GSE262927/lineage_tracing_cohort/"),
        bullet("docs/PIPELINE_AS_RUN.md and docs/ANALYSIS_RATIONALE.md"),
    ]
    story.append(two_col(repro_left, repro_right))
    story.append(Spacer(1, 7 * mm))
    story.append(para("Primary references", "Subsection"))
    refs = [
        "Kadur Lakshminarasimha Murthy P, Sontake V, Tata A, et al. <i>Human distal lung maps and lineage hierarchies reveal a bipotent progenitor.</i> Nature. 2022;604:111-119. doi:10.1038/s41586-022-04541-3. GEO: GSE178360.",
        "Niethamer TK, Planer JD, Morley MP, et al. <i>Longitudinal single-cell profiles of lung regeneration after viral infection reveal persistent injury-associated cell states.</i> Cell Stem Cell. 2025;32:302-321.e6. doi:10.1016/j.stem.2024.12.002. GEO: GSE262927.",
        "Wolock SL, Lopez R, Klein AM. <i>Scrublet: computational identification of cell doublets in single-cell transcriptomic data.</i> Cell Systems. 2019;8:281-291.e9. doi:10.1016/j.cels.2018.11.005.",
    ]
    for ref in refs:
        story.append(bullet(ref))
    story.append(Spacer(1, 6 * mm))
    story.append(para(
        "Report scope: synthesis of tracked repository artefacts; no new statistical model is introduced by document generation. All identity language follows the candidate/analogue guardrails recorded in the source analysis.",
        "Caption",
    ))
    return story


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    doc = PortfolioDoc(OUT)
    doc.build(build_story())
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
