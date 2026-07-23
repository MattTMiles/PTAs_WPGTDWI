"""easel5_style — shared style module for the main-paper figure set (EASEL-5).

Publication conventions (binding, from the EASEL-5 brief):
  * matplotlib only; serif fonts (STIX mathtext, no external TeX dependency);
  * colorblind-safe palette: Okabe-Ito, assigned in FIXED order, never cycled;
  * AAS/PRD column widths: single 3.5 in, double 7.25 in;
  * every figure saved as vector PDF + 400-dpi PNG via savefig();
  * sequential data uses one perceptually-uniform ramp (viridis); diverging data
    uses a two-hue map with a neutral midpoint (RdBu_r); never a rainbow.

Every generator imports this module and nothing else styles a figure.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# geometry — AAS/PRD column widths (inches)
# ----------------------------------------------------------------------
SINGLE = 3.5      # \columnwidth (AAS 3.52, PRD 3.375 — 3.5 fits both at scale 1)
DOUBLE = 7.25     # full text width (AAS 7.25, PRD 6.75 — scale to fit)
GOLDEN = 0.618

# ----------------------------------------------------------------------
# palette — Okabe-Ito (colorblind-safe), FIXED assignment order
# ----------------------------------------------------------------------
BLUE = "#0072B2"     # primary series / data
ORANGE = "#E69F00"   # secondary series (e.g. the VLBI tier)
GREEN = "#009E73"    # the honest/adopted/confirmed object
VERMIL = "#D55E00"   # the retired/refuted/danger object
PURPLE = "#CC79A7"   # tertiary series
SKY = "#56B4E9"      # quaternary / light accent
YELLOW = "#F0E442"   # highlight fill only, never a line
INK = "#1a1a1a"
GREY = "#8a8a8a"
LIGHT = "#d9d9d9"

# status fills (P12 caveats table): colorblind-separable and label-carrying
OK_FILL, OK_EDGE = "#d5e8d9", "#009E73"        # TESTED
WARN_FILL, WARN_EDGE = "#fdebc8", "#E69F00"    # PENDING
BAD_FILL, BAD_EDGE = "#f6d6c9", "#D55E00"      # OPEN

SEQ_CMAP = "viridis"     # magnitude (heatmaps)
DIV_CMAP = "RdBu_r"      # signed residuals, neutral midpoint

# ----------------------------------------------------------------------
# rcParams — serif, recessive axes, print sizes
# ----------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.0,
    "axes.labelsize": 8.0,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "legend.fontsize": 6.8,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": False,
    "ytick.right": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.color": "#e3e3e3",
    "grid.linewidth": 0.5,
    "legend.frameon": False,
    "figure.dpi": 140,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,          # embed TrueType — journal requirement
    "ps.fonttype": 42,
})

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def savefig(fig, name):
    """Write vector PDF + 400-dpi PNG into PAPER_figures/ and close."""
    fig.savefig(os.path.join(OUT, f"{name}.pdf"))
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=400)
    plt.close(fig)
    print(f"wrote {name}.pdf/.png")


def panel_label(ax, s, x=0.0, y=1.02, **kw):
    """'(a)'-style bold panel tag, axes coords, above the frame."""
    ax.text(x, y, s, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.5, weight="bold", **kw)


def marginal_hatch(ax, x0, y0, w=1.0, h=1.0, **kw):
    """The floor-provenance MARGINAL convention: hatched overlay on a heatmap cell."""
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((x0, y0), w, h, fill=False, hatch="//",
                           edgecolor=kw.pop("edgecolor", "#c9c9c9"),
                           lw=kw.pop("lw", 0.0), zorder=kw.pop("zorder", 3), **kw))


NOTE_BOX = dict(boxstyle="round,pad=0.30", fc="white", ec=GREY, lw=0.6, alpha=0.92)
