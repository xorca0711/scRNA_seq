"""Render P0 author-label coverage only; no expression or gene-program results."""
from pathlib import Path
import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import pandas as pd

BASE = Path(__file__).resolve().parents[1]


def panel(ax, table, columns, ticks, row_labels, title):
    values = table[columns].to_numpy(int)
    ax.imshow(values >= 30, aspect="auto", interpolation="none",
              cmap=ListedColormap(["#edf1f4", "#267c78"]), vmin=0, vmax=1)
    for row in range(len(values)):
        for col in range(len(columns)):
            ax.text(col, row, str(values[row, col]), ha="center", va="center", fontsize=7,
                    color="white" if values[row, col] >= 30 else "#4b5563")
    ax.set_xticks(range(len(columns)), ticks, fontsize=9)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(table)), row_labels, fontsize=7)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=42)


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.facecolor": "white"})
    local = pd.read_csv(BASE / "tables/d1_sample_coverage.csv")
    local = local.loc[local.condition.notna()].sort_values(["sacrifice_day", "sample_id"])
    high = pd.read_csv(BASE / "tables/d1_strunz_high_resolution_coverage.csv")
    high["day"] = high.time_point.str.extract(r"(\d+)").astype(int)
    high = high.sort_values(["day", "identifier"])
    gut = pd.read_csv(BASE / "tables/v1_haber_verified_mouse_coverage.csv").sort_values("mouse_id")
    fig = plt.figure(figsize=(14, 10.5))
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.12], wspace=.95)
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[0, 1])
    right = grid[0, 2].subgridspec(3, 1, height_ratios=[1, .7, 2.4])
    ax3 = fig.add_subplot(right[0, 0])
    panel(ax1, local, ["AT2", "Alveolar_transitional", "AT1"], ["AT2", "Transitional", "AT1"],
          [f"{s.removeprefix('EEM-scRNA-')} | d{int(d)}" for s, d in zip(local.sample_id, local.sacrifice_day)],
          "A  Influenza cohort\n1/25 annotated animals: complete triplet")
    panel(ax2, high, ["AT2", "Krt8+ ADI", "AT1"], ["AT2", "ADI", "AT1"],
          [f"{s} | d{d}" for s, d in zip(high.identifier, high.day)],
          "B  Bleomycin epithelium\n10/32 libraries: complete triplet")
    panel(ax3, gut, ["Stem", "Enterocyte.Immature.Proximal", "Enterocyte.Mature.Proximal"],
          ["Stem", "Immature\nproximal", "Mature\nproximal"],
          [s.replace("Control-", "") for s in gut.mouse_id],
          "C  Intestinal atlas\n2/4 verified mice: proximal branch")
    notes = fig.add_subplot(right[2, 0]); notes.axis("off")
    notes.text(0, 1, "Coverage audit", fontsize=11, weight="bold", va="top")
    paragraphs = ["Numbers are retained cells in published states. No program scoring was performed.",
        "The planned floor is 30 cells per state and at least 3 independent units per contrast.",
        "A: Eight libraries lack the relevant author-annotated cohort and are omitted here.",
        "B: No time point has three complete triplets. Library-to-animal independence and the time window remain to be checked.",
        "C: Ten batches are not ten mice. Six other batches have not been mapped to mice here.",
        "Coverage does not establish transition identity, conservation or causation."]
    notes.text(0, .9, "\n\n".join(textwrap.fill(p, 37) for p in paragraphs),
               fontsize=8.5, linespacing=1.3, va="top")
    fig.suptitle("A0 | Initial dataset eligibility", x=.065, y=.98, ha="left", fontsize=20, weight="bold")
    fig.legend(handles=[Patch(facecolor="#267c78", label="At least 30 cells"),
                        Patch(facecolor="#edf1f4", label="Below 30 cells")],
               loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5, .025))
    fig.subplots_adjust(top=.83, bottom=.09, left=.08, right=.98)
    folder = BASE / "figures"; folder.mkdir(exist_ok=True)
    fig.savefig(folder / "p0_state_coverage.png", dpi=170)
    plt.close(fig)
    print("Saved figures/p0_state_coverage.png")


if __name__ == "__main__":
    main()
