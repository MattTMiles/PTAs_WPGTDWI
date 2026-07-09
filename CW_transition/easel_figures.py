"""EASEL — figure regeneration from EXISTING npz. No new computation.
Produces:
  trackB_F2_ladder.png       (F2 lever-arm ladder spectrum + the GAP)
  anchor_census_Ktable.png   (K_real vs K_feather per pulsar)
  anchor_census_Ptrue.png    (P(true) vs N_CW for 3 certifiers + J0437)
Prints the key numbers plotted so they can be checked against the npz.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CWT = "/home/mattm/projects/HSYMT/CW_transition"
INK = "#1a1a2e"; MUT = "#6b7280"
C_REAL = "#2563eb"; C_FEATH = "#e07b39"; C_GAP = "#c0392b"
CERT = {"J1713+0747": "#2563eb", "J1909-3744": "#16a34a",
        "J0711-6830": "#9333ea", "J0437-4715": "#e07b39"}


def fig_f2_ladder():
    d = np.load(f"{CWT}/trackB_F2_ladder.npz", allow_pickle=True)
    rows = d["rows"]; names = d["names"]; tol = d["tol_sky"]
    finite = np.isfinite(tol)
    order = np.argsort(tol[finite])[::-1]
    ts = tol[finite][order]                     # loosest -> tightest
    src = rows[finite, 0][order]; psr = rows[finite, 1][order].astype(int)
    rank = np.arange(1, len(ts) + 1)

    FLOOR_LO, FLOOR_HI = 0.05, 0.21             # blind-float ceiling band
    PULLIN = 1e-4                               # L2c Newton pull-in
    loose = ts[0]
    gap_lo = FLOOR_LO / loose; gap_hi = FLOOR_HI / loose

    print("=== F2 LADDER (checked vs npz) ===")
    print(f"  rungs finite = {finite.sum()}/{len(tol)}")
    print(f"  loosest rung = {loose:.3e} scaled  ({names[psr[0]]}, loud {int(src[0])})")
    print(f"  tightest rung = {ts[-1]:.3e} scaled")
    print(f"  rungs >= 0.05 (in float band) = {int((ts >= FLOOR_LO).sum())}  -> THE GAP")
    print(f"  gap span above loosest: {gap_lo:.0f}x .. {gap_hi:.0f}x")
    print(f"  top rungs: " + ", ".join(f"{names[psr[i]]}={ts[i]:.2e}" for i in range(3)))

    fig, ax = plt.subplots(figsize=(10, 6.2))
    ax.scatter(rank, ts, s=14, color=INK, zorder=3, label="per (pulsar, loud-source) rung")

    # reference bands / lines
    ax.axhspan(FLOOR_LO, FLOOR_HI, color=C_FEATH, alpha=0.18, zorder=0)
    ax.axhline(FLOOR_LO, color=C_FEATH, lw=1.2, ls="--")
    ax.axhline(FLOOR_HI, color=C_FEATH, lw=1.2, ls="--")
    ax.axhline(PULLIN, color=C_REAL, lw=1.4, ls="-.")
    ax.text(len(ts), FLOOR_HI * 1.05, "blind-float ceiling band  (0.05–0.21)",
            color=C_FEATH, ha="right", va="bottom", fontsize=10, fontweight="bold")
    ax.text(len(ts), PULLIN * 0.72, "L2c Newton pull-in  (1e-4)",
            color=C_REAL, ha="right", va="top", fontsize=10, fontweight="bold")

    # THE GAP band (loosest rung -> float floor): zero rungs
    ax.axhspan(loose, FLOOR_LO, color=C_GAP, alpha=0.09, zorder=0)
    ax.annotate("THE GAP — zero rungs\n"
                f"{gap_lo:.0f}–{gap_hi:.0f}× above loosest",
                xy=(len(ts) * 0.5, np.sqrt(loose * FLOOR_LO)),
                color=C_GAP, ha="center", va="center", fontsize=12, fontweight="bold")

    # annotate loosest rung + nearest-pulsar (J0437-class) rungs
    ax.annotate(f"loosest rung {loose:.2e}\n({names[psr[0]]})",
                xy=(1, loose), xytext=(len(ts) * 0.14, loose * 2.6),
                color=INK, fontsize=10, ha="left",
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.1))
    dy = [1.9, 1.25, 0.62]      # stagger to avoid collisions
    for i in range(3):
        ax.annotate(names[psr[i]], xy=(rank[i], ts[i]),
                    xytext=(rank[i] + 10, ts[i] * dy[i]), fontsize=8.5, color=MUT,
                    arrowprops=dict(arrowstyle="-", color=MUT, lw=0.6))

    ax.set_yscale("log")
    ax.set_xlabel("rung rank (sorted loosest → tightest)")
    ax.set_ylabel("sky registration tolerance  (scaled source units, 1-rad shift)")
    ax.set_title("F2 lever-arm ladder — the cascade cannot bridge the blind→pull-in gap",
                 fontsize=12.5, fontweight="bold")
    ax.set_ylim(ts.min() * 0.5, FLOOR_HI * 2.5)
    ax.grid(True, which="both", axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(f"{CWT}/trackB_F2_ladder.png", dpi=140)
    plt.close(fig)


def fig_ktable():
    d = np.load(f"{CWT}/anchor_a1_Ktable.npz", allow_pickle=True)
    jn = d["jn"]; Kr = d["K_canon"]; Kf = d["K_feath"]
    order = np.argsort(Kr)                      # best (smallest K) first
    rank = np.arange(1, len(jn) + 1)

    print("\n=== ANCHOR K-TABLE (checked vs npz) ===")
    print(f"  n pulsars = {len(jn)}")
    print(f"  K_canon (real): min={Kr.min():.2f} ({jn[order[0]]})  max={Kr.max():.0f}")
    print(f"  K_feath: min={Kf.min():.2f}  max={Kf.max():.0f}")
    print(f"  K<=1: {int((Kr<=1).sum())}   K<=3: {int((Kr<=3).sum())}   K<=10: {int((Kr<=10).sum())}")
    print("  best 6 real-prior anchor candidates:")
    for i in order[:6]:
        print(f"    {jn[i]:12s} K_real={Kr[i]:8.2f}  K_feath={Kf[i]:8.2f}")

    fig, ax = plt.subplots(figsize=(10, 6.2))
    ax.scatter(rank, Kr[order], s=18, color=C_REAL, zorder=3, label="K real (canonical prior)")
    ax.scatter(rank, Kf[order], s=14, color=C_FEATH, alpha=0.7, zorder=2,
               marker="s", label="K feather (baseline prior)")

    for K, lab in [(1, "K=1 (single fringe)"), (3, "K=3"), (10, "K=10")]:
        ax.axhline(K, color=MUT, lw=1.0, ls="--")
        ax.text(len(jn), K * 1.06, lab, color=MUT, ha="right", va="bottom", fontsize=9)

    ytxt = [8.0, 22.0, 40.0, 200.0, 340.0, 560.0]   # staggered label heights
    for j, i in enumerate(order[:6]):
        r = j + 1
        ax.annotate(jn[i], xy=(r, Kr[i]), xytext=(r + 6, ytxt[j]),
                    fontsize=8.5, color=INK,
                    arrowprops=dict(arrowstyle="-", color=MUT, lw=0.7))

    ax.set_yscale("log")
    ax.set_xlabel("pulsar rank (sorted by K real, best first)")
    ax.set_ylabel("K  =  wrong-fringe count in ±3σ prior window")
    ax.set_title("Anchor census — no pulsar reaches K≤3 under real EM priors",
                 fontsize=12.5, fontweight="bold")
    ax.grid(True, which="both", axis="y", alpha=0.2)
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{CWT}/anchor_census_Ktable.png", dpi=140)
    plt.close(fig)


def fig_ptrue():
    d = np.load(f"{CWT}/anchor_a2_diag.npz", allow_pickle=True)
    diag = d["diag"][0]; popdiag = d["popdiag"][0]
    ncw = [1, 4, 8, 16]
    targets = ["J1713+0747", "J1909-3744", "J0711-6830", "J0437-4715"]
    names = list(diag[1]["names"]); pnames = list(popdiag["names"])

    print("\n=== P(true) vs N_CW (checked vs npz) ===")
    fig, ax = plt.subplots(figsize=(10, 6.2))
    for t in targets:
        i = names.index(t); c = CERT[t]
        med_l = [np.median(diag[n]["lit_P_true"][:, i]) for n in ncw]
        med_s = [np.median(diag[n]["script_P_true"][:, i]) for n in ncw]
        ax.plot(ncw, med_l, "-o", color=c, lw=1.8, ms=6, label=f"{t} (lit prior)")
        ax.plot(ncw, med_s, "--s", color=c, lw=1.4, ms=5, alpha=0.65,
                markerfacecolor="none")
        pi = pnames.index(t)
        pl = popdiag["lit_P_true"][pi]
        ax.scatter([16], [pl], s=150, marker="*", color=c, edgecolor="k",
                   linewidth=0.6, zorder=5)
        print(f"  {t}: lit med {[f'{x:.2f}' for x in med_l]}  "
              f"script med {[f'{x:.2f}' for x in med_s]}  pop*(lit)={pl:.2f}")

    ax.axhline(0.9, color=C_GAP, lw=1.3, ls="-.")
    ax.text(1, 0.905, "certified-Bayes threshold  P(true) = 0.9",
            color=C_GAP, ha="left", va="bottom", fontsize=10, fontweight="bold")

    # legend keys for line-style + pop marker
    from matplotlib.lines import Line2D
    style = [Line2D([], [], color=MUT, ls="-", marker="o", label="lit prior (median)"),
             Line2D([], [], color=MUT, ls="--", marker="s", markerfacecolor="none",
                    label="script prior (median)"),
             Line2D([], [], color=MUT, ls="none", marker="*", ms=12,
                    markerfacecolor=MUT, label="population draw (3 loud+13 faint)")]
    leg1 = ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.add_artist(leg1)
    ax.legend(handles=style, loc="upper left", frameon=False, fontsize=9)

    ax.set_xscale("log"); ax.set_xticks(ncw); ax.set_xticklabels(ncw)
    ax.set_xlabel("N_CW  (number of injected continuous-wave sources)")
    ax.set_ylabel("P(true fringe)")
    ax.set_ylim(-0.03, 1.05)
    ax.set_title("Anchor certification is confusion-driven — P(true) rises only with N_CW",
                 fontsize=12.5, fontweight="bold")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(f"{CWT}/anchor_census_Ptrue.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    fig_f2_ladder()
    fig_ktable()
    fig_ptrue()
    print("\nsaved 3 figures to CW_transition/")
