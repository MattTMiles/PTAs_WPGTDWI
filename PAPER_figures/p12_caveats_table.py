"""P12 — THE CAVEATS TABLE as a figure.

Rows = the paper's claims; columns = tested domains; cells = TESTED (naming the
campaign that measured it) / PENDING (naming the queued campaign) / OPEN / not
applicable. Built from STORY.md's scope lines (sections cited per row below);
this is the one figure whose source is the tracked claim skeleton, not an npz.

Run:  python p12_caveats_table.py
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import easel5_style as st

# (label, fill, edge)
T = ("TESTED", st.OK_FILL, st.OK_EDGE)
P = ("PENDING", st.WARN_FILL, st.WARN_EDGE)
O = ("OPEN", st.BAD_FILL, st.BAD_EDGE)
NA = ("—", "#f2f2f2", "#cccccc")

COLS = ["noise\nrealism", "floor\nstatistics", "host basis\n(wrong cpart)",
        "waveform\nsystematics", "prior\nmis-centring", "real\ndata"]

# rows: (claim, STORY anchor, [6 cells]); each cell = (status, campaign note)
ROWS = [
    ("Blind fringe-ID no-go ($f=6.9\\times10^{-7}$)", "S4.1",
     [(NA, "zero-noise,\ninfo-theoretic"), (NA, "no floor term"), (NA, ""),
      (P, "EOB tier"), (P, "ROBUST"), (O, "REAL-ARRAY")]),
    ("Chirp-mass wall (tiers fail; $>20\\times$)", "S4.2",
     [(NA, "zero-noise,\ninfo-theoretic"), (NA, "no floor term"), (NA, ""),
      (P, "EOB tier"), (P, "ROBUST"), (O, "REAL-ARRAY")]),
    ("Certified-count distr. + selection fn", "S3.2–3",
     [(T, "FORGE"), (T, "RECUT"), (O, "D1 open"),
      (NA, "circular only"), (P, "ROBUST"), (O, "REAL-ARRAY")]),
    ("The criterion (floors, purity)", "S5, S6.5",
     [(T, "ANCHOR"), (T, "ANCHOR§4/RECUT"), (O, "D3+D4\nrejected"),
      (T, "CHORUS refits\nper mix"), (P, "ROBUST"), (O, "REAL-ARRAY")]),
    ("The onset surface ($N_{\\rm onset}=59$)", "S6.3",
     [(T, "ANCHOR"), (T, "RECUT (v2.2)"), (O, "D1 open"),
      (NA, "circular pop."), (P, "ROBUST"), (O, "REAL-ARRAY")]),
    ("The switch ($e=0.5$ single / $e=0.3$ pair)", "S7.6",
     [(T, "CHORUS\n(drawn noise)"), (T, "RECUT (v2.2)"), (O, "D1 open"),
      (P, "EOB tier\n(toy comb)"), (P, "ROBUST"), (O, "REAL-ARRAY")]),
    ("Loop safety (holds/inert/motion)", "S8.5",
     [(T, "IGNITE-2 /\nEMBER"), (T, "per-cell refits"), (T, "EMBER scrambled\n(measures, not closes)"),
      (T, "CHORUS loops"), (P, "ROBUST\n(tol grid)"), (O, "REAL-ARRAY")]),
    ("Cascade arithmetic (selectivity)", "S8.6",
     [(T, "SPARK nulls\n(N=1300)"), (T, "v2.2 gate\npassed"), (P, "EPOCH\n(EM-mediated)"),
      (NA, "circular twin"), (P, "ROBUST"), (O, "REAL-ARRAY")]),
    ("Siren payoff ($\\sigma_{D_L}/D_L=6–12\\%$)", "S9",
     [(P, "Cramér–Rao,\nzero-noise"), (NA, "given seeds"), (O, "seeds assume\ntrue cpart"),
      (P, "EOB tier\n(self-siren)"), (T, "SIREN arm C\n(within-fringe)"), (O, "REAL-ARRAY")]),
]

SAMPLER_ROWS = [
    ("Sampler: coverage", "S10.1.5",
     [(T, "SAMPLER S-4\n(N=3 banked)"), (NA, ""), (NA, ""), (NA, ""),
      (NA, ""), (O, "REAL-ARRAY")]),
    ("Sampler: SBC of fringe posterior", "S10.1",
     [(P, "SAMPLER"), (P, "SAMPLER"), (NA, ""), (NA, ""), (NA, ""), (O, "REAL-ARRAY")]),
    ("Sampler: high-dim adaptation", "S10.1.6",
     [(T, "noted: chaotic\nAdam retired"), (NA, ""), (NA, ""), (NA, ""),
      (NA, ""), (NA, "")]),
]


def main():
    nrow = len(ROWS) + len(SAMPLER_ROWS)
    ncol = len(COLS)
    fig, ax = plt.subplots(figsize=(st.DOUBLE, 0.42 * nrow + 1.0))
    ax.set_xlim(0, ncol + 3.1)
    ax.set_ylim(-nrow - 0.8, 1.1)
    ax.axis("off")

    x0 = 3.1  # claim-column width
    for j, c in enumerate(COLS):
        ax.text(x0 + j + 0.5, 0.55, c, ha="center", va="center",
                fontsize=6.6, weight="bold")

    def draw_row(i, claim, anchor, cells, shade=False):
        y = -i - 1
        if shade:
            ax.add_patch(Rectangle((0, y), ncol + x0, 1, fc="#f7f7f7", ec="none",
                                   zorder=0))
        ax.text(0.02, y + 0.5, claim, ha="left", va="center", fontsize=6.2)
        ax.text(x0 - 0.06, y + 0.5, anchor, ha="right", va="center",
                fontsize=5.4, color=st.GREY)
        for j, ((lab, fc, eccol), note) in enumerate(cells):
            ax.add_patch(Rectangle((x0 + j + 0.035, y + 0.06), 0.93, 0.88,
                                   fc=fc, ec=eccol, lw=0.7, zorder=2))
            dy = 0.30 if note else 0.0
            ax.text(x0 + j + 0.5, y + 0.5 + dy, lab, ha="center", va="center",
                    fontsize=5.9, weight="bold", zorder=3,
                    color={"TESTED": st.OK_EDGE, "PENDING": "#a06800",
                           "OPEN": st.BAD_EDGE}.get(lab, st.GREY))
            if note:
                ax.text(x0 + j + 0.5, y + 0.32, note, ha="center", va="center",
                        fontsize=4.6, zorder=3, color=st.INK)

    for i, (claim, anchor, cells) in enumerate(ROWS):
        draw_row(i, claim, anchor, cells, shade=(i % 2 == 1))

    ysep = -len(ROWS) - 0.02
    ax.plot([0, ncol + x0], [ysep, ysep], color=st.INK, lw=0.9)
    ax.text(0.04, ysep - 0.28, "sampler row-group", fontsize=6.2, weight="bold",
            color=st.INK, va="center", style="italic")
    for i, (claim, anchor, cells) in enumerate(SAMPLER_ROWS):
        draw_row(len(ROWS) + i + 0.55, claim, anchor, cells,
                 shade=(i % 2 == 1))

    fig.suptitle("Claim-by-domain caveat matrix (from the tracked claim skeleton)",
                 fontsize=8.5, weight="bold", y=0.90)
    st.savefig(fig, "P12_caveats_table")


if __name__ == "__main__":
    main()
