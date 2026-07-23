"""P4 — THE SWITCH.

(a) certified counts vs eccentricity for single- AND multi-member arms, both
    tiers, corrected floors, whiskers to the count at floor+error (CONFIRMED
    test), thresholds marked (e=0.5 single / e=0.3 pair);
(b) the mechanism — counts vs total active-channel budget n_active with the
    ~30 switch line and the equal-kappa contrast (kappa fixed at 2.65, grade
    flips with channels);
(c) LOTTERY — P(switch-on) over f_ecc x e_char, conservative-rule boundary,
    validated tier-2 points overlaid.

Sources: reports/ch_analysis_recut.npz (counts/floors, criterion-v2.2);
CHORUS_results/ch_analysis.npz (surface_nactive; the n_active convention is
S7.6.4c: TOTAL channel budget, circular member = 1); reports/MAGPIE_audit.md
J1 (kappa=2.65 at the members' f_orb~10^-8.2, from ATLAS_results/
atlas_kappa_forb{1,2}.npz); reports/lottery_tier1.npz (P_chan/P_thr surfaces)
+ reports/lottery_tier2.npz (validated cells).
"""
import json

import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"
KAPPA_E03 = 2.65   # MAGPIE J1: kappa at e=0.3, members' f_orb ~ 10^-8.2 (ATLAS kappa banks)


def main():
    A = np.load(f"{REPO}/reports/ch_analysis_recut.npz", allow_pickle=True)
    tags = [str(t) for t in A["surface_tags"]]
    corr, corr_lo = A["surface_corr"], A["surface_corr_lo"]
    C = np.load(f"{REPO}/CHORUS_results/ch_analysis.npz", allow_pickle=True)
    ctags = [str(t) for t in C["surface_tags"]]
    nact = C["surface_nactive"]

    def g(a, t, tl=tags):
        return float(a[tl.index(t)])

    fig, (ax, axm, axl) = plt.subplots(
        1, 3, figsize=(st.DOUBLE, 2.35), gridspec_kw={"width_ratios": [1.2, 1, 1.15]})
    fig.subplots_adjust(wspace=0.42)

    # ---------- (a) counts vs e, per n_ecc, both tiers ----------
    evals = [0.3, 0.5, 0.7]
    ecols = {1: st.BLUE, 2: st.GREEN, 3: st.PURPLE}
    for ne in (1, 2, 3):
        for tier, ls, mk in (("lit", "-", "o"), ("vlbi", "--", "s")):
            ts = [f"m{ne}e0{int(10*e)}_{tier}" for e in evals]
            y = [g(corr, t) for t in ts]
            ylo = [g(corr_lo, t) for t in ts]
            ax.plot(evals, y, ls, marker=mk, ms=2.6, lw=0.9, color=ecols[ne],
                    label=f"$n_{{\\rm ecc}}={ne}$, {tier}" if True else None)
            ax.vlines(evals, ylo, y, color=ecols[ne], lw=0.6, alpha=0.7)
    circ = g(corr, "m0e00_lit")
    ax.plot([0.0], [circ], "o", ms=3, color=st.INK)
    ax.plot([0.0], [g(corr, "m0e00_vlbi")], "s", ms=3, color=st.INK)
    ax.annotate("all-circular", xy=(0.0, circ), xytext=(0.02, 1.6), fontsize=5.4,
                arrowprops=dict(arrowstyle="->", color=st.INK, lw=0.5))
    ax.axhline(1.0, color=st.VERMIL, lw=1.0, ls="--")
    ax.text(0.695, 1.12, "onset bar", color=st.VERMIL, fontsize=5.6, ha="right",
            weight="bold")
    ax.axvline(0.5, color=ecols[1], lw=0.6, ls=":", alpha=0.8)
    ax.axvline(0.3, color=ecols[2], lw=0.6, ls=":", alpha=0.8)
    ax.text(0.492, 2.6, "switch-on, 1 member", fontsize=5.2, color=ecols[1],
            rotation=90, va="bottom", ha="right")
    ax.text(0.292, 2.6, "switch-on, $\\geq$2 members", fontsize=5.2,
            color=ecols[2], rotation=90, va="bottom", ha="right")
    ax.set_xlim(-0.04, 0.76)
    ax.set_ylim(0, 7.6)
    ax.set_xlabel("eccentricity of the eccentric member(s)")
    ax.set_ylabel("correct certifications / realisation")
    st.panel_label(ax, "(a)")
    ax.legend(loc="upper left", fontsize=4.6, ncol=2, columnspacing=0.8,
              handlelength=1.6, borderaxespad=0.2)
    ax.text(0.98, 0.02, "whisker = count at floor$+$error\n(the CONFIRMED test)",
            transform=ax.transAxes, ha="right", fontsize=5.0, style="italic")

    # ---------- (b) the mechanism: counts vs channel budget ----------
    for tier, mk in (("lit", "o"), ("vlbi", "s")):
        xs, ys = [], []
        for t in tags:
            if t.endswith(tier) and "eU" not in t:
                xs.append(g(nact, t, ctags))
                ys.append(g(corr, t))
        col = st.BLUE if tier == "lit" else st.ORANGE
        axm.plot(xs, ys, mk, ms=2.8, mfc=col, mec="none", alpha=0.8, label=tier)
    axm.axvline(30, color=st.GREEN, lw=1.0)
    axm.text(32, 0.35, "$n_{\\rm active}\\approx 30$", fontsize=5.6,
             color=st.GREEN, weight="bold")
    axm.axhline(1.0, color=st.VERMIL, lw=0.8, ls="--")
    # the equal-kappa contrast: m1e03 / m2e03 / m3e03 (lit) at fixed kappa
    for t, lab in (("m1e03_lit", "$n_{\\rm ecc}$=1"), ("m2e03_lit", "2"),
                   ("m3e03_lit", "3")):
        x0, y0 = g(nact, t, ctags), g(corr, t)
        axm.plot([x0], [y0], "o", ms=5.5, mfc="none", mec=st.VERMIL, mew=0.9,
                 zorder=6)
        axm.annotate(lab, xy=(x0, y0), xytext=(x0 - 1.5, y0 + 0.55), fontsize=5.2,
                     color=st.VERMIL)
    axm.text(0.36, 0.80,
             f"circled: $e=0.3$ arms —\n$\\kappa$ fixed at {KAPPA_E03}, channels\n"
             f"23$\\to$30$\\to$37, grade flips",
             transform=axm.transAxes, va="top", fontsize=5.4, color=st.VERMIL)
    axm.set_xlabel("active-harmonic channel budget $n_{\\rm active}$")
    axm.set_ylabel("correct certifications / realisation")
    axm.set_xlim(10, 115)
    st.panel_label(axm, "(b)")
    axm.legend(loc="lower right", fontsize=5.6)

    # ---------- (c) LOTTERY P(switch-on) ----------
    L = np.load(f"{REPO}/reports/lottery_tier1.npz", allow_pickle=True)
    fg, eg = L["f_grid"], L["e_grid"]
    Pc, Pt = L["P_chan"], L["P_thr"]
    im = axl.pcolormesh(fg, eg, Pc, cmap=st.SEQ_CMAP, vmin=0, vmax=1,
                        shading="nearest")
    # conservative threshold-rule boundary: where P_thr first becomes nonzero
    axl.contour(fg, eg, Pt, levels=[1e-9], colors="white", linewidths=1.6)
    axl.contour(fg, eg, Pt, levels=[1e-9], colors=[st.VERMIL], linewidths=0.8)
    # validated tier-2 cells (rows_json): plotted at (f_ecc = n_ecc/3 of the loud
    # members, e). Overlay uses the e_char axis directly; f_ecc for the m-cells
    # is the eccentric fraction of the 3 loud members.
    T2 = json.loads(str(np.load(f"{REPO}/reports/lottery_tier2.npz",
                                allow_pickle=True)["rows_json"]))
    for row in T2:
        ne = int(row["n_ecc"])
        e0 = float(row["e"])
        grade = str(row["grade"]).upper()
        carried = bool(row.get("carried", False))
        f_ecc = ne / 3.0 + (0.03 if carried else 0.0)   # dgx-carried offset
        mk = {"CONFIRMED": "o", "MARGINAL": "s", "BELOW": "X"}[grade]
        axl.plot([f_ecc], [e0], mk, ms=3.3,
                 mfc=("white" if grade == "CONFIRMED" else
                      st.ORANGE if grade == "MARGINAL" else st.VERMIL),
                 mec=st.INK, mew=0.6, zorder=6)
    axl.set_xlabel("eccentric fraction of loud members, $f_{\\rm ecc}$")
    axl.set_ylabel("characteristic eccentricity $e_{\\rm char}$")
    st.panel_label(axl, "(c)")
    cb = fig.colorbar(im, ax=axl, fraction=0.05, pad=0.02)
    cb.set_label("$P$(switch-on)", fontsize=6)
    cb.ax.tick_params(labelsize=5)
    axl.text(0.03, 0.97, "line: conservative\nthreshold rule",
             transform=axl.transAxes, va="top", fontsize=5.2, color="white",
             weight="bold")
    axl.text(0.04, 0.03,
             "validated cells: $\\circ$ CONFIRMED, sq. MARGINAL,\n"
             "$\\times$ below; sub-rule cells DO certify —\n"
             "the rule is one-sided conservative",
             transform=axl.transAxes, ha="left", fontsize=4.8, color="white")

    st.savefig(fig, "P4_the_switch")
    print(f"P4: m1e03/e05/e07 lit = "
          f"{[g(corr, f'm1e0{i}_lit') for i in (3, 5, 7)]}, "
          f"m2e03 lit {g(corr, 'm2e03_lit'):.2f}; tier2 rows {len(T2)}")


if __name__ == "__main__":
    main()
