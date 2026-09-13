"""Figure conventions for the Cardoso trials: one validated palette, recessive chrome.

Colours come from the dataviz reference palette and were validated with its
`scripts/validate_palette.js` on the light surface before use (results in the
C5 run record). Colour follows the entity, never its rank; text uses ink
tokens, never a series colour.
"""

from __future__ import annotations

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
DEEMPH = "#cfcec8"
SLOT = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#eda100"}
FLOX_PLUS, FLOX_FLOX = "#5598e7", "#184f95"
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def apply(plt) -> None:
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
        "font.size": 10, "text.color": INK, "axes.labelcolor": INK_2,
        "axes.edgecolor": AXIS, "axes.linewidth": 1, "xtick.color": MUTED, "ytick.color": INK_2,
        "xtick.labelsize": 9, "ytick.labelsize": 9.5, "axes.titlesize": 11.5,
        "axes.titleweight": "semibold", "axes.titlecolor": INK, "axes.titlelocation": "left",
        "axes.titlepad": 10, "legend.frameon": False, "legend.fontsize": 9,
    })


def recessive(ax, grid_axis: str = "x") -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(AXIS)
    ax.spines["bottom"].set_color(AXIS)
    ax.grid(axis=grid_axis, color=GRID, linewidth=1, linestyle="-")
    ax.set_axisbelow(True)
    ax.tick_params(length=0)


def bare(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def px_scale(ax) -> tuple[float, float]:
    """Data units per display pixel, (x, y). Call after limits and position are final."""
    ax.figure.canvas.draw()
    inv = ax.transData.inverted()
    (x0, y0), (x1, y1) = inv.transform([(0, 0), (1, 1)])
    return abs(x1 - x0), abs(y1 - y0)


def hbar(ax, y: float, x0: float, x1: float, thickness_px: float, colour: str,
         scale: tuple[float, float], radius_px: float = 4.0) -> None:
    """Horizontal bar: square at the baseline, a small rounded data-end."""
    from matplotlib.patches import FancyBboxPatch, Rectangle

    width = x1 - x0
    if width <= 0:
        return
    dx, dy = scale
    height = thickness_px * dy
    radius = min(radius_px * dx, width / 2)
    ax.add_patch(FancyBboxPatch((x0, y - height / 2), width, height,
                                boxstyle=f"round,pad=0,rounding_size={radius}",
                                mutation_aspect=dy / dx, facecolor=colour, edgecolor="none",
                                linewidth=0, zorder=3))
    ax.add_patch(Rectangle((x0, y - height / 2), max(width - radius, 0), height,
                           facecolor=colour, edgecolor="none", linewidth=0, zorder=3))
