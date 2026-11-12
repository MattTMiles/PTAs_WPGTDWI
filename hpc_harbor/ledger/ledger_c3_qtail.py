#!/usr/bin/env python3
"""LEDGER C3 -- Q-TAIL PAYOFF PROPAGATION: what (1-q) wrong-fringe mass does to the siren.

THE DEFECT UNDER AUDIT. SIREN SS4 quotes sigma(D_L)/D_L "conditional on N_seed certified
pulsar terms". The certification criterion the programme actually uses admits a seed at
q_max > 0.9 (glacier_loop.py:104, QBAR). So a "certified" seed carries UP TO 10% posterior
mass on a WRONG fringe -- and the payoff table propagates NONE of it. The banked payoff is
the width of the ALL-SEEDS-RIGHT component of a mixture, quoted as if it were the width of
the posterior.

The correction is not a fudge factor: it is the mixture the criterion already implies.

THE MIXTURE (stated exactly; nothing here is fitted).
Let the N seeds be independently right with probability q. Write x = ln(D_L/D_L_true).
With k of N seeds on a wrong fringe (probability Binom(N, 1-q)):

    p(x) = SUM_k C(N,k) (1-q)^k q^(N-k) * [ component_k(x) ]

Two DECLARED limits bracket component_k, and the answer is quoted BETWEEN them (the house
two-bound discipline, SPARK-3 SS-Fix-3 / faint_fisher_bounds):

  BENIGN-WRONG (optimistic).  A wrong seed is UNINFORMATIVE about the chirp mass: it
    contributes no valid pulsar-term lever, and the analysis is effectively an
    (N-k)-seed analysis. Unbiased.
        component_k = Normal(0, s[N-k])
    where s[.] is the BANKED frac_DL ladder (SIREN_results/siren_summary.npz,
    A_exact__frac_DL, freq-ranked nested sequence N = 0,1,2,3,5). N = 4 is
    log-interpolated between the banked N = 3 and N = 5 and FLAGGED as interpolated.
    This is a LOWER bound on the damage: it assumes the wrong fringe announces itself
    by carrying no information, which a wrong fringe emphatically does not do.

  ADVERSARIAL-WRONG (pessimistic).  A wrong seed is CONFIDENTLY wrong: it contributes
    full Fisher information, centred on a displaced chirp mass. The displacement is the
    one that reabsorbs the fringe. A fringe is BY CONSTRUCTION one cycle of pulsar-term
    phase (trackB_b1_core.FringeTables: fringe spacing dL is the 2*pi mode spacing), and
    at Fisher level the single-seed chirp-mass width is the phase-precision reciprocal,
    so a one-fringe error displaces log10 Mc by

        Delta_1 = 2*pi * s_mc[1]                      (s_mc[1] = banked single-seed
                                                       sigma(log10 Mc) at this cell)

    resisted by the other seeds in Fisher proportion, so with k wrong out of N the joint
    centre moves by (k/N) * Delta_1. Propagated through discovery's own amplitude
    relation log10 D_L = (5/3) log10 Mc + (2/3) log10 f - log10 h:

        mu_k = +- (5/3) * ln10 * (k/N) * Delta_1      (sign of the fringe error is
                                                       symmetric: a two-point +- mixture)
        component_k = 0.5*Normal(+mu_k, s[N]) + 0.5*Normal(-mu_k, s[N])

    This is an UPPER bound: it gives the wrong fringe its full confidence AND a full
    one-cycle displacement, with no partial cancellation between multiple wrong seeds.

WHAT IS REPORTED (three numbers, because a heavy-tailed mixture has no single width):
  core       -- the all-right component width s[N]. THIS IS WHAT SIREN SS4 QUOTES.
  rms        -- sqrt(Var) of the mixture. Dominated by the rare, very wide components.
  hw90       -- half-width of the CENTRAL 90% credible interval of the mixture. The
                operational number: "where is D_L, at 90% credibility".
  P_allright -- q^N. The probability that the quoted core width is the right answer at all.

DECLARED LIMITATIONS (so the number is not over-read):
  * q is treated as INDEPENDENT per seed. Correlated fringe errors (a biased source
    template mis-registering several pulsars the same way -- exactly the S4.17/D1 wander
    mechanism) would be WORSE than this, not better.
  * The banked ladder is a nested FREQ-RANKED sequence; dropping a wrong seed is modelled
    as dropping the LAST seed of that sequence, which understates the loss when the wrong
    seed is the informative one. Also conservative-in-the-optimistic-direction.
  * q_max is a criterion-side confidence, not a calibrated frequentist coverage. If the
    fringe posterior is over-confident (the ANCHOR D2 finding: the Gumbel floor is biased
    PERMISSIVE at high zero-fraction), the true (1-q) is LARGER than the nominal.
  * Asimov/Fisher throughout, inherited from SIREN. These remain Cramer-Rao cores.

REPORT-ONLY. Nothing here arms a protocol step or changes a banked verdict.
Outputs: reports/ledger_c3_qtail.npz, reports/LEDGER_C3_qtail.png

Usage:  python3 hpc_harbor/ledger/ledger_c3_qtail.py [--repo <path>]
"""
import argparse
import os
import sys

import numpy as np
from scipy.special import comb

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito, assigned in FIXED order and never cycled (identity is also carried by
# line style + direct label, so nothing is colour-alone).
INK = "#1a1a1a"
MUTED = "#6b6b6b"
C_CORE = "#0072B2"      # blue   -- the SIREN-quoted core
C_BEN = "#009E73"       # green  -- benign-wrong bound
C_ADV = "#D55E00"       # orange -- adversarial-wrong bound
GRID = "#d9d9d9"

N_LADDER = np.array([0, 1, 2, 3, 5])          # banked freq-ranked nested sequence
LADDER_COL = {0: 0, 1: 1, 2: 2, 3: 3, 5: 4}   # -> column of A_exact__frac_DL
Q_LIST = (0.90, 0.95, 0.99)
Z90 = 1.6448536269514722            # the Gaussian 90% half-width, so core and
                                   # mixture are compared LIKE FOR LIKE, never sigma-vs-interval


def ladder_width(fr, n):
    """frac_DL at n seeds from the banked nested ladder; log-interpolate n = 4."""
    if n in LADDER_COL:
        return float(fr[LADDER_COL[n]]), False
    lo, hi = 3, 5
    w = (n - lo) / (hi - lo)
    return float(np.exp((1 - w) * np.log(fr[LADDER_COL[lo]])
                        + w * np.log(fr[LADDER_COL[hi]]))), True


def mixture(N, q, fr, s_mc1, mode):
    """Return (weights, means, sds) of the Gaussian mixture in x = ln(D_L/D_L_true)."""
    ws, mus, sds = [], [], []
    sN, _ = ladder_width(fr, N)
    d1 = 2.0 * np.pi * s_mc1                       # one-fringe log10 Mc displacement
    for k in range(N + 1):
        pk = comb(N, k) * (1 - q) ** k * q ** (N - k)
        if pk <= 0:
            continue
        if mode == "benign":
            sk, _ = ladder_width(fr, N - k)
            ws.append(pk); mus.append(0.0); sds.append(sk)
        else:                                       # adversarial
            if k == 0:
                ws.append(pk); mus.append(0.0); sds.append(sN)
            else:
                mu = (5.0 / 3.0) * np.log(10.0) * (k / N) * d1
                ws += [0.5 * pk, 0.5 * pk]
                mus += [+mu, -mu]
                sds += [sN, sN]
    return np.array(ws), np.array(mus), np.array(sds)


def mix_stats(ws, mus, sds, nq=200001, span=12.0):
    """RMS and central-90% half-width of a Gaussian mixture, on a dense grid."""
    mean = float(np.sum(ws * mus))
    var = float(np.sum(ws * (sds ** 2 + (mus - mean) ** 2)))
    hi = span * max(float(np.max(sds)), float(np.max(np.abs(mus))) + 1e-12)
    x = np.linspace(-hi, hi, nq)
    from scipy.special import erf
    cdf = np.zeros_like(x)
    for w, m, s in zip(ws, mus, sds):
        cdf += w * 0.5 * (1.0 + erf((x - m) / (s * np.sqrt(2.0))))
    lo = float(np.interp(0.05, cdf, x)); up = float(np.interp(0.95, cdf, x))
    return np.sqrt(var), 0.5 * (up - lo), x, cdf


def mix_pdf(x, ws, mus, sds):
    p = np.zeros_like(x)
    for w, m, s in zip(ws, mus, sds):
        p += w * np.exp(-0.5 * ((x - m) / s) ** 2) / (s * np.sqrt(2 * np.pi))
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/taylor_group/matt_miles/PTAs_WPGTDWI")
    a = ap.parse_args()
    z = np.load(os.path.join(a.repo, "SIREN_results", "siren_summary.npz"),
                allow_pickle=True)
    FR = z["A_exact__frac_DL"]              # (3 f, 3 mc, 7 tags)
    SMC = z["A_exact__sig_log10_mc"]
    gf, gmc = z["grid_f"], z["grid_mc"]
    DL = z["DL_Mpc"]

    # the SIREN headline cells: the ones SS4.2 quotes as "6-12%" at N_seed = 3-5
    CELLS = [(1, 2), (2, 1), (2, 2)]       # (f=-8.0,mc=9.5), (-7.7,9.0), (-7.7,9.5)
    rows = []
    print("[LEDGER-C3] q-tail propagation through the banked SIREN payoff ladder\n")
    print("  cell (log10 f, log10 Mc, D_L)      N   q     P_all   core    "
          "core_hw90  benign_hw90 (x)  adv_hw90 (x)")
    for (i, j) in CELLS:
        fr = FR[i, j][:5]                   # nested freq-ranked columns only
        s_mc1 = float(SMC[i, j][LADDER_COL[1]])
        for N in (3, 5):
            core, interp = ladder_width(fr, N)
            for q in Q_LIST:
                out = {}
                for mode in ("benign", "adversarial"):
                    ws, mus, sds = mixture(N, q, fr, s_mc1, mode)
                    rms, hw, _, _ = mix_stats(ws, mus, sds)
                    out[mode] = (rms, hw)
                pall = q ** N
                chw = Z90 * core
                print(f"  ({gf[i]:+.1f}, {gmc[j]:+.1f}, {DL[i, j]:6.2f} Mpc)   "
                      f"{N}  {q:.2f}  {pall:.3f}  {core:.4f}   {chw:.4f}     "
                      f"{out['benign'][1]:.4f} ({out['benign'][1]/chw:4.2f}x)   "
                      f"{out['adversarial'][1]:.4f} ({out['adversarial'][1]/chw:4.2f}x)")
                rows.append((float(gf[i]), float(gmc[j]), float(DL[i, j]), N, q,
                             pall, core, chw, out["benign"][0], out["benign"][1],
                             out["adversarial"][0], out["adversarial"][1],
                             out["benign"][1] / chw, out["adversarial"][1] / chw,
                             2 * np.pi * s_mc1, interp))
        print()
    R = np.array(rows, float)

    # ------------------------------------------------------------------ figure
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.6, 4.0))
    i, j = 1, 2                                     # the headline cell
    fr = FR[i, j][:5]; s_mc1 = float(SMC[i, j][LADDER_COL[1]]); N = 3
    core, _ = ladder_width(fr, N)

    # LEFT: the mixture itself, benign limit, at the three q values
    x = np.linspace(-0.9, 0.9, 4001)
    for q, ls in zip(Q_LIST, ("-", "--", ":")):
        ws, mus, sds = mixture(N, q, fr, s_mc1, "benign")
        axL.plot(x, mix_pdf(x, ws, mus, sds), ls, color=C_BEN, lw=2.0,
                 label=f"q = {q:.2f}")
    axL.plot(x, mix_pdf(x, np.array([1.0]), np.array([0.0]), np.array([core])),
             "-", color=C_CORE, lw=2.0, label="core (SIREN §4)")
    axL.set_yscale("log"); axL.set_ylim(1e-3, 30)
    axL.set_xlabel(r"$\ln(D_L / D_L^{\rm true})$", color=INK)
    axL.set_ylabel("posterior density", color=INK)
    axL.set_title(f"the mixture the criterion implies\n"
                  fr"$N_{{\rm seed}}={N}$, benign-wrong limit, "
                  fr"$\log_{{10}}f={gf[i]:+.1f}$, $\log_{{10}}M_c={gmc[j]:+.1f}$",
                  fontsize=9.5, color=INK)
    axL.legend(frameon=False, fontsize=8.5, loc="upper right")
    axL.grid(True, color=GRID, lw=0.6, alpha=0.8); axL.set_axisbelow(True)

    # RIGHT: effective width vs q, both bounds, both N -- no second axis
    qq = np.linspace(0.85, 1.0, 61)
    for N, mk in ((3, "o"), (5, "s")):
        core_n, _ = ladder_width(fr, N)
        hb = [mix_stats(*mixture(N, q, fr, s_mc1, "benign"))[1] for q in qq]
        ha = [mix_stats(*mixture(N, q, fr, s_mc1, "adversarial"))[1] for q in qq]
        axR.plot(qq, hb, "-" if N == 3 else "--", color=C_BEN, lw=2.0)
        axR.plot(qq, ha, "-" if N == 3 else "--", color=C_ADV, lw=2.0)
        axR.axhline(Z90 * core_n, color=C_CORE, lw=1.4,
                    ls="-" if N == 3 else "--", alpha=0.9)
        axR.annotate(f"N={N}", (0.858, hb[0]), fontsize=8, color=MUTED)
    for q in Q_LIST:
        axR.axvline(q, color=MUTED, lw=0.8, ls=":")
        axR.annotate(f"{q:.2f}", (q, 0.92), xycoords=("data", "axes fraction"),
                     fontsize=8, color=MUTED, ha="center")
    axR.set_yscale("log")
    axR.set_xlabel("per-seed fringe confidence  q", color=INK)
    axR.set_ylabel(r"90% half-width of $\ln D_L$", color=INK)
    axR.set_title("what the q-tail costs the siren\n"
                  "(blue = SIREN's quoted core; green = benign; orange = adversarial)",
                  fontsize=9.5, color=INK)
    axR.grid(True, color=GRID, lw=0.6, alpha=0.8); axR.set_axisbelow(True)
    for ax in (axL, axR):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(colors=MUTED, labelsize=8.5)
    fig.suptitle("LEDGER C3 — the payoff is quoted on the all-right component of a mixture",
                 fontsize=11, color=INK, y=1.02)
    fig.tight_layout()
    png = os.path.join(a.repo, "reports", "LEDGER_C3_qtail.png")
    fig.savefig(png, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"[LEDGER-C3] figure -> {png}")

    out = os.path.join(a.repo, "reports", "ledger_c3_qtail.npz")
    np.savez(out, table=R,
             columns=np.array(["log10_f", "log10_mc", "DL_Mpc", "N_seed", "q",
                               "P_allright", "core_fracDL", "core_hw90",
                               "benign_rms", "benign_hw90", "adv_rms", "adv_hw90",
                               "inflation_benign", "inflation_adv",
                               "delta1_log10mc", "ladder_interpolated"]),
             note=("LEDGER C3: sigma(D_L)/D_L with the (1-q) wrong-fringe mass convolved "
                   "in, as the Gaussian mixture the q_max>0.9 criterion implies. Two "
                   "declared bounds (benign / adversarial wrong-seed). core_fracDL is the "
                   "number SIREN SS4 quotes. REPORT-ONLY."))
    print(f"[LEDGER-C3] banked -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
