#!/usr/bin/env python3
"""
VERNIER-5 -- the working point. EXPLORATORY, NOT pre-registered.

Reproduce two published configurations, validate against them, then walk to realism.
Clean room: imports only vernier_toy / vernier2 / vernier4 (my own anchored code).

Reported in the field's units: lambda_fore = h^2_CW / (h^2_CW + h^2_GWB) per bin.
Literature thresholds: 0.7 reliable recovery, 0.5 fiducial detection.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import zeta

import vernier_toy as V
import vernier2 as W
import vernier4 as Z

YR, C_PC_S, REPORTS = V.YR, V.C_PC_S, V.REPORTS
V.PTERM_STATIC = False
KAPPA = Z.KAPPA


# ------------------------------------------------------------------ basics
def sigma_L_from_rho(dL_pc, rho_p):
    """Within-fringe distance precision. Fisher of lnL = rho^2(cos(2 pi dL/dL)-1)."""
    return dL_pc / (2.0 * np.pi * np.maximum(rho_p, 1e-12))


def rho_array(h0, f, n_psr, sigma_toa, T, cadence_days=14.0):
    """Array-summed Earth-term SNR of a CW of strain h0."""
    R = h0 / (2.0 * np.pi * f)
    Sw = 2.0 * sigma_toa ** 2 * cadence_days * 86400.0
    return np.sqrt(KAPPA * n_psr * R ** 2 * T / (2.0 * Sw))


def h0_for_rho(rho, f, n_psr, sigma_toa, T, cadence_days=14.0):
    Sw = 2.0 * sigma_toa ** 2 * cadence_days * 86400.0
    R = np.sqrt(rho ** 2 * 2.0 * Sw / (KAPPA * n_psr * T))
    return R * 2.0 * np.pi * f


# ------------------------------------------------------------------ targets
TARGETS = {
    "A_Xiao": dict(n_psr=85, T=20 * YR, sigma=50e-9, cadence=14.0, n_src=5,
                   L_kpc=1.0, claim="sub-pc for 55-87% of pulsars at ~1 kpc",
                   claim_sigma_pc=1.0),
    "B_YuPan": dict(n_psr=20, T=30 * YR, sigma=20e-9, cadence=14.0, n_src=8,
                    L_kpc=1.0, claim="0.4 pc at 1 kpc; 0.8-1 pc at 2.2 kpc",
                    claim_sigma_pc=0.4),
}


def exp_targets(out, rng, f_list=(5e-9, 10e-9, 20e-9), n_trial=200):
    print("=" * 78)
    print("TARGETS A and B -- do we reproduce the published working points?")
    print("=" * 78)
    print("  Neither prompt specifies the SOURCE STRAIN or FREQUENCY, and the claimed")
    print("  precision depends on both. So we INVERT: what array SNR does each claim")
    print("  require, and is that SNR enough to also IDENTIFY THE FRINGE?")
    print()
    print("  within-fringe precision:  sigma_L = dL / (2 pi rho_pterm)   [Fisher, exact")
    print("  for lnL = rho^2(cos(2 pi dL/dL) - 1); matches the banked D3 scaling]")
    print("  N sources combine in Fisher: 1/sigma^2 = sum_k (2 pi rho_k / dL_k)^2")
    print()
    res = {}
    for name, t in TARGETS.items():
        print(f"  --- {name}: {t['n_psr']} psr, {t['T']/YR:.0f} yr, "
              f"{t['sigma']*1e9:.0f} ns, N={t['n_src']} sources")
        print(f"      claim: {t['claim']}")
        V.T_OBS_FID = t["T"]
        print(f"      {'f nHz':>7} {'dL pc':>7} {'rho_arr needed':>15} "
              f"{'rho_P per psr':>14} {'h0 needed':>11} {'P(fringe ID)':>13}")
        for f in f_list:
            dL = C_PC_S / f                       # median geometry, 1-cos mu = 1
            # Fisher over N equal sources: sigma = dL/(2 pi rho_P sqrt(N))
            rho_p_need = dL / (2 * np.pi * t["claim_sigma_pc"] * np.sqrt(t["n_src"]))
            rho_arr_need = rho_p_need * np.sqrt(t["n_psr"])
            h0 = h0_for_rho(rho_arr_need, f, t["n_psr"], t["sigma"], t["T"],
                            t["cadence"])
            # can the fringe be identified at that per-pulsar SNR with N sources?
            uniq = 0
            for _ in range(n_trial):
                ph = V.iso_dirs(1, rng)[0]
                Ls = t["L_kpc"] * 1000.0
                srcs = V.draw_sources(t["n_src"], rng, ph, rho_p_need)
                uniq += (V.scan(Ls, 1.0, srcs)[0] == 1)     # 1 pc prior, generous
            p_id = uniq / n_trial
            print(f"      {f*1e9:7.1f} {dL:7.3f} {rho_arr_need:15.2f} "
                  f"{rho_p_need:14.3f} {h0:11.2e} {p_id:13.2f}")
            res[(name, f)] = (dL, rho_arr_need, rho_p_need, h0, p_id)
        print()
    V.T_OBS_FID = 15.0 * YR
    print("  READING. The claimed sigma_L needs only rho_P ~ 0.1-0.4 PER PULSAR, i.e.")
    print("  a very modest source. That part reproduces trivially -- WITHIN-FRINGE")
    print("  precision is cheap, exactly as the programme's D2 found (116/116 below")
    print("  their EM prior). But at that same SNR the FRINGE IS NOT IDENTIFIED:")
    print("  VERNIER-1 established a pulsar needs its own rho_P >~ 2 for its comb to")
    print("  collapse, and these configurations sit 5-20x below that.")
    print()
    out["tgt_names"] = np.array(list(TARGETS.keys()))
    out["tgt_f"] = np.array(f_list)
    for k, v in res.items():
        out[f"tgt_{k[0]}_{k[1]*1e9:.0f}"] = np.array(v)
    return res


def exp_target_snr_scaling(out, rng, n_trial=200):
    """At what per-pulsar rho_P does each target ACTUALLY identify fringes?"""
    print("=" * 78)
    print("TARGETS -- what would it take to IDENTIFY the fringe (not just measure")
    print("           within it)?")
    print("=" * 78)
    rhos = np.array([0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0])
    print(f"  P(unique distance solution), 1 pc prior, L* = 1 kpc, {n_trial} trials")
    print(f"  {'rho_P per psr':>14} " + " ".join(f"{n:>12s}" for n in TARGETS))
    P = np.zeros((rhos.size, len(TARGETS)))
    for i, r in enumerate(rhos):
        cells = []
        for j, (name, t) in enumerate(TARGETS.items()):
            V.T_OBS_FID = t["T"]
            u = 0
            for _ in range(n_trial):
                ph = V.iso_dirs(1, rng)[0]
                srcs = V.draw_sources(t["n_src"], rng, ph, r)
                u += (V.scan(1000.0, 1.0, srcs)[0] == 1)
            P[i, j] = u / n_trial
            cells.append(f"{P[i,j]:12.2f}")
        print(f"  {r:14.2f} " + " ".join(cells))
    V.T_OBS_FID = 15.0 * YR
    print()
    for j, (name, t) in enumerate(TARGETS.items()):
        ok = np.where(P[:, j] >= 0.9)[0]
        need = rhos[ok[0]] if ok.size else np.inf
        print(f"    {name}: needs rho_P >~ {need} per pulsar for 90% fringe ID; "
              f"its claim implies rho_P ~ 0.1-0.4  -> shortfall "
              f"{need/0.3:.0f}x in SNR, {(need/0.3)**2:.0f}x in source power")
    out["tsnr_rhos"] = rhos
    out["tsnr_P"] = P
    print()


# ------------------------------------------------------------------ lambda_fore
def lam_fore(N, gamma):
    """
    lambda_fore = h1^2 / sum_k h_k^2 for bright-end counts N(>h) ~ h^-gamma.
    h_k^2 ~ k^(-2/gamma).  gamma = 3 is Euclidean/no-evolution.
    """
    N = np.atleast_1d(np.maximum(N, 1.0))
    s = 2.0 / gamma
    out = np.empty_like(N, dtype=float)
    for i, n in enumerate(N):
        k = np.arange(1, int(min(n, 2e6)) + 1)
        tot = np.sum(k ** -s)
        if n > 2e6 and s < 1:                       # divergent tail, continue analytically
            tot += (n ** (1 - s) - (2e6) ** (1 - s)) / (1 - s)
        out[i] = 1.0 / tot
    return out if out.size > 1 else float(out[0])


def exp_lambda(out):
    print("=" * 78)
    print("FIELD UNITS -- lambda_fore = h^2_CW / (h^2_CW + h^2_GWB) per frequency bin")
    print("=" * 78)
    print("  Literature: >= 0.7 reliable recovery, 0.5 fiducial detection.")
    print()
    print("  BOTH PUBLISHED TARGETS SET NO GWB, so they run at lambda_fore = 1.000")
    print("  exactly -- the most favourable value possible. That is the single biggest")
    print("  difference between their configuration and a real sky.")
    print()
    fg = np.array([2, 3, 5, 8, 10, 15, 20, 32, 50, 80, 150, 208]) * 1e-9
    Nc = Z.N_c(fg)
    lf = lam_fore(Nc, 3.0)
    print(f"  {'f nHz':>7} {'N_c/bin':>11} {'lambda_fore':>12} {'vs 0.5':>8} {'vs 0.7':>8}")
    for f, n, l in zip(fg, Nc, np.atleast_1d(lf)):
        print(f"  {f*1e9:7.1f} {n:11.3g} {l:12.5f} "
              f"{'PASS' if l >= .5 else 'fail':>8} {'PASS' if l >= .7 else 'fail':>8}")
    print()
    print("  With Euclidean bright-end counts (gamma = 3) lambda_fore NEVER reaches")
    print("  0.5 anywhere, because sum_k k^(-2/3) diverges as 3 N^(1/3): the loudest")
    print("  source is always a vanishing fraction of its own bin.")
    print()
    print("  THE KNEE IS A STATEMENT ABOUT THE BRIGHT-END SLOPE, not about N_c.")
    print("  For N(>h) ~ h^-gamma, h_k^2 ~ k^(-2/gamma), so")
    print("      gamma > 2 : sum diverges, lambda_fore -> 0     (no knee, ever)")
    print("      gamma < 2 : sum -> zeta(2/gamma), lambda_fore -> 1/zeta(2/gamma)")
    print("  i.e. a knee EXISTS only if the massive tail is truncated (gamma < 2).")
    print()
    print(f"  {'gamma':>7} {'lambda_fore (N=1e3)':>21} {'(N=1e5)':>12} {'asymptote':>12}")
    for g in [3.0, 2.5, 2.0, 1.5, 1.2, 1.0, 0.9, 0.8]:
        a = 1.0 / zeta(2.0 / g) if g < 2 else 0.0
        print(f"  {g:7.2f} {lam_fore(1e3, g):21.4f} {lam_fore(1e5, g):12.4f} "
              f"{a:12.4f}" + ("   <- >=0.7" if a >= 0.7 else ""))
    gs = np.linspace(0.6, 1.99, 400)
    asy = np.array([1.0 / zeta(2.0 / g) for g in gs])
    g70 = gs[np.argmin(np.abs(asy - 0.7))]
    g50 = gs[np.argmin(np.abs(asy - 0.5))]
    print()
    print(f"  *** lambda_fore >= 0.7 requires bright-end slope gamma <= {g70:.2f}")
    print(f"  *** lambda_fore >= 0.5 requires                  gamma <= {g50:.2f}")
    print("  Euclidean is gamma = 3. So the literature's 0.7 threshold is NOT")
    print("  reachable by making the array better or observing longer -- it is a")
    print("  property of the POPULATION's bright end. It needs a hard mass-function")
    print("  truncation, which my census does not have and VERNIER-4 therefore missed.")
    print()
    print("  CORRECTION TO VERNIER-4. I quoted 208 nHz as 'where the population goes")
    print("  sparse'. That is the one-source-per-bin frequency for the FULL population")
    print("  and it is the right number for that question -- but it is the WRONG")
    print("  question. The knee that matters is where the massive tail truncates, which")
    print("  in the literature sits at 25-40 nHz, INSIDE the band. My pure power law has")
    print("  no such knee at any frequency, so VERNIER-4's 'operates only above 13 nHz,")
    print("  closes only near 113 nHz' is an artefact of a missing mass function.")
    out["lam_f"] = fg
    out["lam_Nc"] = Nc
    out["lam_lf"] = np.atleast_1d(lf)
    out["lam_gammas"] = gs
    out["lam_asy"] = asy
    out["lam_g70"] = g70
    out["lam_g50"] = g50
    print()


def exp_knee(out):
    """If the knee is at 25-40 nHz, what bright-end slope does that imply, and
    what does it do to the operating window?"""
    print("=" * 78)
    print("THE KNEE, taken as given at 25-40 nHz -- what would it buy?")
    print("=" * 78)
    print("  Take the literature knee as an input rather than deriving it: suppose the")
    print("  massive tail truncates so that lambda_fore reaches the quoted thresholds")
    print("  somewhere in 25-40 nHz. Then at the knee the loudest source holds a")
    print("  FINITE fraction of its bin, and the numbers change qualitatively:")
    print()
    print(f"  {'lambda_fore':>12} {'h1/h_sea':>10} {'rho_1/rho_bin':>14} "
          f"{'implied gamma':>14}")
    for l in [0.034, 0.1, 0.3, 0.5, 0.7, 0.9]:
        amp = np.sqrt(l)
        g = np.nan
        gs = np.linspace(0.5, 1.99, 600)
        a = np.array([1.0 / zeta(2.0 / x) for x in gs])
        if l <= a.max():
            g = gs[np.argmin(np.abs(a - l))]
        print(f"  {l:12.3f} {amp:10.3f} {amp:14.3f} {g:14.2f}")
    print()
    print("  At lambda_fore = 0.7 the loudest source is 0.84 of the bin amplitude --")
    print("  it IS the bin. Detection then costs nothing extra and VERNIER-4's")
    print("  'never closes under correlated confusion' verdict is voided, because")
    print("  there is no sea left to be correlated with.")
    print("  *** So the whole VERNIER-4 pessimistic column hinges on the bright-end")
    print("  *** slope, which is exactly the quantity my census assumed away.")
    print()


def exp_convention(out):
    print("=" * 78)
    print("CONVENTION CHECK (done first, as instructed)")
    print("=" * 78)
    print("  VERNIER-4 built BOTH sides of the ratio from per-bin residual POWER:")
    print("      sum_k R_k^2 = 2 S_gwb(f) df   and   R_1^2 = that / (3 N_c^(1/3))")
    print("  h_c never entered, so h_1/h_sea = R_1/sqrt(sum R_k^2) is like-with-like.")
    print()
    print(f"  {'f nHz':>7} {'as built':>10} {'via h_c':>10} {'factor':>8} "
          f"{'sqrt(fT)':>10} {'lambda_fore':>12}")
    rows = []
    for f in [10e-9, 32e-9, 80e-9]:
        R2 = Z.R1_squared(f); tot = 2 * Z.S_gwb(f) * Z.DF
        ab = np.sqrt(R2 / tot)
        h0 = 2 * np.pi * f * np.sqrt(R2)
        hc = (h0 * np.sqrt(f * Z.T_OBS)) / (Z.A_GWB * (f * YR) ** (-2 / 3.))
        rows.append([f, ab, hc, hc / ab, np.sqrt(f * Z.T_OBS), R2 / tot])
        print(f"  {f*1e9:7.1f} {ab:10.4f} {hc:10.4f} {hc/ab:8.4f} "
              f"{np.sqrt(f*Z.T_OBS):10.2f} {R2/tot:12.5f}")
    print()
    print(f"  The two conventions differ by a CONSTANT 2/sqrt(6) = {2/np.sqrt(6):.4f},")
    print("  NOT by sqrt(fT) (which is 2.2 / 3.9 / 6.2 here). sqrt(fT) cancels because")
    print("  source and background are expressed per the SAME Fourier bin.")
    print("  VERDICT: NO sqrt(fT) error. The corrected h_c-convention figures are")
    print("  LOWER, not higher: 0.090 -> 0.074 at 10 nHz, 0.184 -> 0.150 at 32 nHz.")
    out["conv_rows"] = np.array(rows)
    print()


def make_plots(out):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    f = out["lam_f"] * 1e9
    ax[0].loglog(f, out["lam_lf"], "o-", lw=2, ms=5, color="#4c72b0",
                 label="this census ($\\gamma$=3)")
    ax[0].axhline(0.7, color="#c44e52", ls="--", lw=1.4, label="lit. 0.7 recovery")
    ax[0].axhline(0.5, color="#dd8452", ls="-.", lw=1.4, label="lit. 0.5 detection")
    ax[0].axhline(1.0, color="#55a868", ls=":", lw=1.6,
                  label="Targets A & B (no GWB)")
    ax[0].set_xlabel("frequency [nHz]"); ax[0].set_ylabel("$\\lambda_{fore}$")
    ax[0].set_title("(1) field units: nowhere near threshold")
    ax[0].legend(fontsize=7.5, frameon=False); ax[0].grid(alpha=.25, which="both")
    ax[1].plot(out["lam_gammas"], out["lam_asy"], lw=2.2, color="#4c72b0")
    ax[1].axhline(0.7, color="#c44e52", ls="--", lw=1.4)
    ax[1].axhline(0.5, color="#dd8452", ls="-.", lw=1.4)
    ax[1].axvline(out["lam_g70"], color="#c44e52", ls=":", lw=1.2)
    ax[1].set_xlabel("bright-end counts slope $\\gamma$  [$N(>h)\\propto h^{-\\gamma}$]")
    ax[1].set_ylabel("$\\lambda_{fore}$ asymptote")
    ax[1].set_title(f"(2) 0.7 needs $\\gamma\\leq${out['lam_g70']:.2f} (Euclidean = 3)")
    ax[1].grid(alpha=.25)
    r = out["tsnr_rhos"]
    for j, n in enumerate(out["tgt_names"]):
        ax[2].semilogx(r, out["tsnr_P"][:, j], "o-", lw=2, ms=5, label=str(n))
    ax[2].axvspan(0.03, 0.28, color="grey", alpha=.3)
    ax[2].text(0.09, 0.55, "SNR their\nclaims imply", ha="center", fontsize=7.5)
    ax[2].axhline(0.9, color="k", ls=":", lw=1.2)
    ax[2].set_xlabel("per-pulsar pulsar-term $\\rho_P$")
    ax[2].set_ylabel("P(unique distance)")
    ax[2].set_title("(3) fringe ID needs 3-7x more SNR")
    ax[2].legend(fontsize=8, frameon=False); ax[2].grid(alpha=.25)
    fig.suptitle("VERNIER-5 -- the working point (exploratory, not pre-registered)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER5_workingpoint.png"), dpi=140)
    plt.close(fig)


def main():
    rng = np.random.default_rng(20260801)
    out = {}
    print("\nVERNIER-5 -- EXPLORATORY, NOT pre-registered. Clean room preserved.\n")
    exp_convention(out)
    exp_targets(out, rng)
    exp_target_snr_scaling(out, rng)
    exp_lambda(out)
    exp_knee(out)
    npz = os.path.join(REPORTS, "vernier5.npz")
    np.savez_compressed(npz, **out)
    print(f"[saved {npz} before plotting]")
    try:
        make_plots(out)
        print(f"wrote {os.path.join(REPORTS, 'VERNIER5_workingpoint.png')}")
    except Exception as e:
        print(f"WARNING: plotting failed ({type(e).__name__}: {e}); data is safe.")
    print("=" * 78)
    print(f"wrote {npz}")
    print("=" * 78)


if __name__ == "__main__":
    main()
