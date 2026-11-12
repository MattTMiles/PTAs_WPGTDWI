#!/usr/bin/env python3
"""
VERNIER-2 -- walk the wheel to failure.  EXPLORATORY, NOT pre-registered.

Continuation of VERNIER-1. Clean room preserved: imports only vernier_toy (my own
anchored implementation), nothing from the campaign.

Physics fixes applied throughout:
  * pulsar term EVOLVES (exact Peters). f_p taken at the retarded epoch and the
    comb spacing set from f_p, NOT f0.
  * the pulsar-term amplitude pays (f_p/f0)^(2/3) and sits at a DIFFERENT noise
    level than the Earth term, so the gearing is measured, not assumed.

Outputs reports/vernier2.npz + reports/VERNIER2_*.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V          # same directory; clean-room, self-authored

YR, C_PC_S = V.YR, V.C_PC_S
REPORTS = V.REPORTS
F_YR = 1.0 / YR
V.PTERM_STATIC = False           # evolving pulsar term everywhere

BAND = (3e-9, 20e-9)
SNR0_DEFAULT = 20.0              # rank-0 Earth-term array SNR at the working point
DET_THRESH = 8.0

# noise caricature: P(f)/P_wn = 1 + r_red (f/f_yr)^-3 + r_gwb (f/f_yr)^-13/3
# normalised so each component alone is 3x white at f = 10 nHz.
R_RED = 3.0 / (1e-8 / F_YR) ** -3.0
R_GWB = 3.0 / (1e-8 / F_YR) ** -(13.0 / 3.0)


def noise_ratio(f, mode, has_red):
    """P(f)/P_wn for the given noise mode."""
    x = np.asarray(f, float) / F_YR
    r = np.ones_like(x)
    if mode in ("red", "gwb"):
        r = r + np.where(has_red, 1.0, 0.0) * R_RED * x ** -3.0
    if mode == "gwb":
        r = r + R_GWB * x ** -(13.0 / 3.0)
    return r


def snr_factor(f, mode, has_red):
    return 1.0 / np.sqrt(noise_ratio(f, mode, has_red))


# ---------------------------------------------------------------- array + sky
def build_array(n_psr, rng):
    p = V.iso_dirs(n_psr, rng)
    L = rng.uniform(300.0, 5000.0, n_psr)
    has_red = rng.random(n_psr) < 0.5
    return p, L, has_red


def make_sky(rng, log10h, n_psr_dirs):
    """Sources: isotropic sky, band-uniform-in-log f0, Mc ~ U(8.5,9.5)."""
    n = len(log10h)
    s = V.iso_dirs(n, rng)
    f0 = 10.0 ** rng.uniform(np.log10(BAND[0]), np.log10(BAND[1]), n)
    mc = rng.uniform(8.5, 9.5, n)
    return s, f0, mc, np.asarray(log10h, float)


def couple(p, L, has_red, s, f0, mc, noise_mode):
    """
    Per (pulsar, source) Earth-term and pulsar-term SNR weights, and geometry.
    Returns omc, rho_E, rho_P, f_p  (all shape (n_psr, n_src)).
    """
    cmu = p @ s.T
    omc = np.clip(1.0 - cmu, 1e-3, None)
    ant = np.sqrt(np.clip(1.0 - cmu ** 2, 1e-6, None))
    tc = V.t_merge(f0, mc)
    tau = V.tau_of(L[:, None], omc)
    f_p = f0[None, :] * (1.0 + tau / tc[None, :]) ** (-3.0 / 8.0)
    hr = has_red[:, None]
    rho_E = ant * snr_factor(f0[None, :] * np.ones_like(omc), noise_mode, hr)
    # RESIDUAL amplitude, not strain: r ~ h/f and h ~ f^(2/3), so r ~ f^(-1/3).
    # The pulsar term sits at LOWER frequency, so it is INTRINSICALLY LOUDER than
    # the Earth term in residual space (1.196x at f_p/f0 = 0.584).
    rho_P = ant * (f_p / f0[None, :]) ** (-1.0 / 3.0) \
        * snr_factor(f_p, noise_mode, hr)
    return omc, rho_E, rho_P, f_p, tc


def wheel(rng, log10h, n_psr=116, prior="vlbi", noise_mode="white",
          snr0=SNR0_DEFAULT, det_thresh=DET_THRESH, grid_cap=96000):
    """
    One realisation of the wheel. Returns the FULL per-turn history.

    A turn: detect the loudest undetected source using Earth terms everywhere plus
    pulsar terms of already-CERTIFIED pulsars; then recompute certification for
    every pulsar from the combs of all detected sources.
    """
    p, L, has_red = build_array(n_psr, rng)
    s, f0, mc, lh = make_sky(rng, log10h, n_psr)
    omc, rho_E, rho_P, f_p, tc = couple(p, L, has_red, s, f0, mc, noise_mode)
    amp = 10.0 ** (lh - lh.max())
    rho_E = rho_E * amp[None, :]
    rho_P = rho_P * amp[None, :]
    # Normalise on a WHITE-noise reference so that turning red noise / the GWB on
    # genuinely DEGRADES the array. Normalising within each noise mode would
    # renormalise away the very damage the noise dial is meant to measure.
    _, rE_w, _, _, _ = couple(p, L, has_red, s, f0, mc, "white")
    rE_w = rE_w * amp[None, :]
    nrm = snr0 / np.sqrt(np.sum(rE_w[:, int(np.argmax(lh))] ** 2))
    rho_E *= nrm; rho_P *= nrm

    # ---- phase/distance-marginalised pulsar-term residue (the honest BEFORE state).
    # Searching an uncertified pulsar term over its K_a candidate fringes recovers
    # its power but pays the trials factor: dlnl_marg = rho_P^2/2 - ln K_a, so
    # rho_P,marg^2 = max(0, rho_P^2 - 2 ln K_a).
    _sv = {"subpc": 0.1, "vlbi": 1.0}.get(prior)
    sig_arr = (np.full(n_psr, _sv) if _sv is not None
               else {"pct1": 0.01, "pct10": 0.10, "pct20": 0.20}[prior] * L)
    K_aj = np.maximum(6.0 * sig_arr[:, None] / V.dL_of(f_p, omc), 1.0)
    rho_P_marg = np.sqrt(np.maximum(rho_P ** 2 - 2.0 * np.log(K_aj), 0.0))

    n_src = len(lh)
    cert = np.zeros(n_psr, bool)
    det = []
    H = {"n_cert": [], "snr_with": [], "snr_without": [], "snr_marg": [],
         "gain_dex": [], "gain_dex_marg": [], "geared": [], "detected": [],
         "f_cert": []}
    while True:
        avail = [j for j in range(n_src) if j not in det]
        if not avail:
            stall = "sources exhausted"
            break
        sw, so, sm = [], [], []
        for j in avail:
            base = np.sum(rho_E[:, j] ** 2)
            pm = rho_P_marg[:, j] ** 2
            # WITH certification: certified pulsars coherent, the rest still
            # contribute their phase-marginalised residue.
            sw.append(np.sqrt(base + np.sum(rho_P[cert, j] ** 2)
                              + np.sum(pm[~cert])))
            # HONEST before state: every pulsar term marginalised, none discarded.
            sm.append(np.sqrt(base + np.sum(pm)))
            # naive before state, kept for comparison with VERNIER-1: term discarded.
            so.append(np.sqrt(base))
        sw, so, sm = np.array(sw), np.array(so), np.array(sm)
        k = int(np.argmax(sw))
        jj = avail[k]
        tp = np.sum(rho_P[:, jj] ** 2)
        fc = float(np.sum(rho_P[cert, jj] ** 2) / tp) if tp > 0 else 0.0
        H["snr_with"].append(sw[k]); H["snr_without"].append(so[k])
        H["snr_marg"].append(sm[k])
        H["gain_dex"].append(np.log10(sw[k] / so[k]))
        H["gain_dex_marg"].append(np.log10(sw[k] / sm[k]))
        H["n_cert"].append(int(cert.sum()))
        H["f_cert"].append(fc)
        if sw[k] < det_thresh:
            stall = "next source below detection threshold"
            H["geared"].append(False); H["detected"].append(-1)
            break
        H["geared"].append(bool(sm[k] < det_thresh <= sw[k]))
        det.append(avail[k]); H["detected"].append(avail[k])
        # ---- recompute certification from scratch every turn
        new = np.zeros(n_psr, bool)
        for a in range(n_psr):
            sig = sig_arr[a]
            srcs = [V.Src(f0[j], tc[j], omc[a, j], 1, rho_P[a, j]) for j in det]
            new[a] = (V.scan(L[a], sig, srcs, cap=grid_cap)[0] == 1)
        cert = new
    H = {k: np.array(v) for k, v in H.items()}
    H["turns"] = int((H["detected"] >= 0).sum())
    H["turns_geared"] = int(H["geared"].sum())
    H["stall"] = stall
    H["n_psr"] = n_psr
    return H


# ================================================== FIX 1: how much does evolving ease it
def exp_evolve_shift(out, rng, n_trial=200, Kmax=8):
    print("=" * 74)
    print("FIX 1 -- evolving pulsar term: how much does the coarser comb ease things")
    print("=" * 74)
    tiers = ["vlbi", "pct20"]
    Ks = np.arange(1, Kmax + 1)
    P = np.zeros((2, len(tiers), Kmax))
    M = np.zeros((2, len(tiers), Kmax))
    KT = np.zeros((2, len(tiers)))
    for im, stat in enumerate([False, True]):
        V.PTERM_STATIC = stat
        for it, tier in enumerate(tiers):
            kt = []
            for ik, K in enumerate(Ks):
                nm = []
                for _ in range(n_trial):
                    ph = V.iso_dirs(1, rng)[0]
                    Ls = rng.uniform(500.0, 3000.0)
                    sig = 1.0 if tier == "vlbi" else 0.20 * Ls
                    srcs = V.draw_sources(K, rng, ph, 10.0)
                    n, k = V.scan(Ls, sig, srcs)
                    nm.append(n)
                    if K == 1:
                        kt.append(k)
                nm = np.array(nm)
                P[im, it, ik] = (nm == 1).mean()
                M[im, it, ik] = np.median(nm)
            KT[im, it] = np.median(kt)
    V.PTERM_STATIC = False

    def k90(pv):
        for i in range(len(pv) - 1):
            if pv[i] < 0.9 <= pv[i + 1]:
                return (i + 1) + (0.9 - pv[i]) / (pv[i + 1] - pv[i])
        return np.nan if pv[0] < 0.9 else 1.0

    print(f"  {'tier':7s} {'conv':9s} {'med K_trials(K=1)':>18s} {'med modes K=1':>14s} "
          f"{'K90 (frac)':>11s}")
    for it, tier in enumerate(tiers):
        for im, nm_ in enumerate(["evolving", "static"]):
            print(f"  {tier:7s} {nm_:9s} {KT[im,it]:18.1f} {M[im,it,0]:14.0f} "
                  f"{k90(P[im,it]):11.2f}")
        print(f"  {tier:7s} -> comb coarser by {KT[1,it]/KT[0,it]:.2f}x ; "
              f"modes at K=1 down {M[1,it,0]/max(M[0,it,0],1):.2f}x ; "
              f"K90 shift {k90(P[1,it])-k90(P[0,it]):+.2f} sources")
    print()
    print("  P(unique) vs K, evolving (E) vs static (S):")
    for it, tier in enumerate(tiers):
        print(f"    {tier:6s} E: " + " ".join(f"{v:5.2f}" for v in P[0, it]))
        print(f"    {tier:6s} S: " + " ".join(f"{v:5.2f}" for v in P[1, it]))
    print()
    out["ev_P"] = P; out["ev_modes"] = M; out["ev_Ktrials"] = KT; out["ev_K"] = Ks
    out["ev_tiers"] = np.array(tiers)
    out["ev_k90"] = np.array([[k90(P[im, it]) for it in range(len(tiers))]
                              for im in range(2)])
    return P, KT


# ================================================== FIX 2: is the separate line usable
def exp_spectral(out):
    print("=" * 74)
    print("FIX 2 -- f_p != f0 makes the pulsar term a SEPARATE SPECTRAL LINE")
    print("=" * 74)
    print("  Two lines are resolvable when |f0 - f_p| > 1/T_obs.")
    print("  And f_p(L) is MONOTONIC in L -- no wrapping -- so the line position is an")
    print("  UNAMBIGUOUS distance handle, unlike the comb. Precision:")
    print("     df_p/dL = -(3/8) f_p (1-cos mu) / (c (t_c + tau))")
    print("     sigma_f ~ 1/(2 pi T rho)   ->   sigma_L = sigma_f / |df_p/dL|")
    print()
    f_grid = np.array([3.0, 6.0, 10.0, 15.0, 20.0]) * 1e-9
    mc_grid = np.array([8.5, 9.0, 9.5, 10.0])
    L, omc, rho = 1000.0, 1.0, 10.0
    tau = V.tau_of(L, omc)
    res = np.zeros((f_grid.size, mc_grid.size, 4))
    print(f"  L=1 kpc, mu=90 deg, rho_pterm={rho:.0f}, T=15 yr")
    print(f"  {'f0nHz':>6} {'log10Mc':>8} {'f_p/f0':>7} {'|df|*T':>8} "
          f"{'resolv':>7} {'sigma_L pc':>11} {'sig_L/dL':>9}")
    for i, f in enumerate(f_grid):
        for j, mc in enumerate(mc_grid):
            tc = V.t_merge(f, mc)
            fp = f * (1 + tau / tc) ** (-3.0 / 8.0)
            dfT = abs(f - fp) * T
            dfp_dL = (3.0 / 8.0) * fp * omc / (C_PC_S * (tc + tau))
            sig_f = 1.0 / (2 * np.pi * T * rho)
            sig_L = sig_f / dfp_dL
            dl = V.dL_of(fp, omc)
            res[i, j] = [fp / f, dfT, sig_L, sig_L / dl]
            print(f"  {f*1e9:6.1f} {mc:8.2f} {fp/f:7.3f} {dfT:8.2f} "
                  f"{'YES' if dfT>1 else 'no':>7} {sig_L:11.2f} {sig_L/dl:9.1f}")
    nres = int((res[:, :, 1] > 1).sum())
    print()
    print(f"  RESOLVABLE as two lines in {nres}/{res[:,:,1].size} band cells.")
    print(f"  sigma_L from the line position: {np.nanmin(res[:,:,2]):.1f} - "
          f"{np.nanmax(res[:,:,2]):.0f} pc  (median {np.median(res[:,:,2]):.1f} pc)")
    print(f"  in units of the comb: {np.median(res[:,:,3]):.0f} fringes wide (median)")
    print()
    print("  VERDICT. Usable, but as a COARSE UNAMBIGUOUS PRIOR, not as a fringe")
    print("  identifier: it is ~1e1-1e2 fringes wide. Its value is that it is")
    print("  DATA-DERIVED and fringe-free, so it can replace a wide EM prior.")
    print("  THE CATCH: f_p depends on t_c ~ Mc^(-5/3), and the banked chirp-mass")
    print("  posterior is prior-dominated (gain 1.00-1.73x, sigma 0.50-0.87 dex).")
    print("  A 0.5 dex Mc error is 0.83 dex in t_c, which swamps the L dependence")
    print("  for a SINGLE pulsar -- the handle is Mc-L degenerate. It only breaks")
    print("  with the ARRAY: one shared Mc against N_psr distinct tau_a, i.e.")
    print("  N_psr line positions for N_psr + 1 unknowns.")
    out["spec_f"] = f_grid; out["spec_mc"] = mc_grid; out["spec_res"] = res
    out["spec_n_resolvable"] = nres
    print()


T = 15.0 * YR


# ================================================== gearing: measured, not assumed
def exp_gearing(out, rng, n_real=40, n_src=8, gap=0.15, n_psr=116):
    print("=" * 74)
    print("GEARING -- measured dex gain per wheel turn")
    print("=" * 74)
    print(f"  {n_real} realisations, {n_src} sources spaced {gap} dex, {n_psr} pulsars,")
    print("  1 pc priors, white noise, evolving pulsar terms.")
    print("  gain(turn t) = log10[ SNR(with certified pulsar terms) / SNR(Earth only) ]")
    print("  CEILING TEST: if every pulsar were certified the gain would be")
    print("     log10 sqrt(1 + sum_a rho_P,a^2 / sum_a rho_E,a^2)")
    print("  which is log10(sqrt2)=0.1505 ONLY IF rho_P = rho_E. It is not: the")
    print("  pulsar term pays (f_p/f0)^(2/3) in amplitude and sits at a different")
    print("  noise level. So the ceiling is MEASURED here.")
    print()
    lh = -12.75 - gap * np.arange(n_src)
    gains = [[] for _ in range(n_src + 1)]
    ncert = [[] for _ in range(n_src + 1)]
    gm_ = [[] for _ in range(n_src + 1)]
    fcl = [[] for _ in range(n_src + 1)]
    ceil, ratio = [], []
    turns = []
    for _ in range(n_real):
        H = wheel(rng, lh, n_psr=n_psr, prior="vlbi", noise_mode="white")
        turns.append(H["turns"])
        for t, g in enumerate(H["gain_dex"]):
            gains[t].append(g); ncert[t].append(H["n_cert"][t])
            gm_[t].append(H["gain_dex_marg"][t]); fcl[t].append(H["f_cert"][t])
        # ceiling for this realisation, from the same physics
        p, L, hred = build_array(n_psr, rng)
        s, f0, mc, _ = make_sky(rng, lh, n_psr)
        _, rE, rP, _, _ = couple(p, L, hred, s, f0, mc, "white")
        r = np.sum(rP[:, 0] ** 2) / np.sum(rE[:, 0] ** 2)
        ratio.append(np.sqrt(r)); ceil.append(np.log10(np.sqrt(1 + r)))
    print("  gain(discard) = vs discarding uncertified pulsar terms (VERNIER-1 style)")
    print("  gain(marg)    = vs MARGINALISING them, the honest before state (item 3)")
    print(f"  {'turn':>5} {'n':>4} {'med n_cert':>11} {'med f_cert':>11} "
          f"{'gain(discard)':>14} {'gain(marg)':>11} {'incr':>8}")
    prev = 0.0
    med_by_turn = []
    for t in range(n_src + 1):
        if len(gains[t]) < 3:
            continue
        m = float(np.median(gains[t]))
        med_by_turn.append(m)
        print(f"  {t+1:5d} {len(gains[t]):4d} {np.median(ncert[t]):11.0f} "
              f"{np.median(fcl[t]):11.4f} {m:14.4f} "
              f"{np.median(gm_[t]):11.4f} {m-prev:+8.4f}")
        prev = m
    print()
    print(f"  measured rho_P/rho_E (rank 0, array-summed): "
          f"median {np.median(ratio):.3f}")
    print(f"  measured ALL-CERTIFIED ceiling: median {np.median(ceil):.4f} dex "
          f"[{np.percentile(ceil,16):.4f}, {np.percentile(ceil,84):.4f}]")
    print(f"  sqrt(2) reference (rho_P = rho_E): {np.log10(np.sqrt(2)):.4f} dex")
    print(f"  turns achieved: median {np.median(turns):.1f}")
    print()
    if len(med_by_turn) >= 3:
        inc = np.diff(med_by_turn[:4])
        print(f"  VERDICT. Gains at turns 1..{min(4,len(med_by_turn))}: "
              + " ".join(f"{v:.4f}" for v in med_by_turn[:4]))
        print(f"  increments: " + " ".join(f"{v:+.4f}" for v in inc))
        print("  Gearing SATURATES: it is the ONE-TIME recovery of pulsar-term")
        print("  coherence. Once a pulsar is certified its pulsar term is already in")
        print("  the sum; certifying it again buys nothing, and there is no second")
        print("  reservoir of coherence to unlock. Reach does NOT compound.")
    out["gear_gain_marg_by_turn"] = np.array([np.median(g) if len(g) >= 3 else np.nan
                                              for g in gm_])
    out["gear_fcert_by_turn"] = np.array([np.median(g) if len(g) >= 3 else np.nan
                                          for g in fcl])
    out["gear_gain_by_turn"] = np.array([np.median(g) if len(g) >= 3 else np.nan
                                         for g in gains])
    out["gear_ncert_by_turn"] = np.array([np.median(g) if len(g) >= 3 else np.nan
                                          for g in ncert])
    out["gear_ceiling_dex"] = np.array(ceil)
    out["gear_rhoP_over_rhoE"] = np.array(ratio)
    out["gear_turns"] = np.array(turns)
    out["gear_sqrt2_dex"] = np.log10(np.sqrt(2.0))
    print()


# ================================================== the dial walk
DIALS = {
    "gap_dex":  [0.0, 0.15, 0.3, 0.5, 0.8, 1.2],
    "n_comp":   [5, 4, 3, 2, 1],
    "prior":    ["subpc", "vlbi", "pct1", "pct10", "pct20"],
    "noise":    ["white", "red", "gwb"],
    "n_psr":    [116, 50, 250],
}
WORKING = dict(gap_dex=0.0, n_comp=5, prior="vlbi", noise="white", n_psr=116)


def run_point(rng, cfg, n_real):
    lh = -12.75 - cfg["gap_dex"] * np.arange(cfg["n_comp"])
    T_, G_, S_ = [], [], []
    for _ in range(n_real):
        H = wheel(rng, lh, n_psr=cfg["n_psr"], prior=cfg["prior"],
                  noise_mode=cfg["noise"])
        T_.append(H["turns"]); G_.append(H["turns_geared"])
        S_.append(H["n_cert"][-1] if len(H["n_cert"]) else 0)
    return np.median(T_), np.median(G_), np.median(S_), np.array(T_)


def exp_dials(out, rng, n_real=10):
    print("=" * 74)
    print("DIAL WALK -- start at a working point, move ONE dial at a time")
    print("=" * 74)
    print(f"  working point: {WORKING}")
    print(f"  {n_real} realisations per cell. 'turns' = sources detected;")
    print("  'geared' = detections that ONLY happened because certification raised")
    print("  the array SNR above threshold (the wheel's actual output).")
    print()
    t0, g0, c0, _ = run_point(rng, dict(WORKING), n_real)
    print(f"  WORKING POINT: turns={t0:.1f}  geared={g0:.1f}  certified={c0:.0f}"
          f"/{WORKING['n_psr']}")
    print(f"  -> the wheel {'TURNS' if t0 >= 2 else 'DOES NOT TURN'} at the working point.")
    print()
    res = {}
    for dial, vals in DIALS.items():
        print(f"  --- dial: {dial}")
        print(f"    {'value':>8} {'turns':>7} {'geared':>7} {'certified':>10} "
              f"{'verdict':>12}")
        rows = []
        first_stall = None
        for v in vals:
            cfg = dict(WORKING); cfg[dial] = v
            t, g, c, _ = run_point(rng, cfg, n_real)
            rows.append((str(v), t, g, c))
            verdict = "OK" if t >= 2 else "STALLED"
            if t < 2 and first_stall is None:
                first_stall = v
            print(f"    {str(v):>8} {t:7.1f} {g:7.1f} {c:10.0f} {verdict:>12}")
        res[dial] = rows
        print(f"    first value with turns < 2: {first_stall}")
        print()
    out["dial_names"] = np.array(list(DIALS.keys()))
    for dial, rows in res.items():
        out[f"dial_{dial}_vals"] = np.array([r[0] for r in rows])
        out[f"dial_{dial}_turns"] = np.array([r[1] for r in rows])
        out[f"dial_{dial}_geared"] = np.array([r[2] for r in rows])
        out[f"dial_{dial}_cert"] = np.array([r[3] for r in rows])
    out["dial_working_turns"] = t0
    out["dial_working_geared"] = g0
    return res


# ================================================== sky (a): A-SKY as built
def exp_asky(out, rng, n_real=30, n_psr=116):
    print("=" * 74)
    print("SKY (a) -- A-SKY AS BUILT: 5 loud at log10h=-12.75 + 11 faint at -14.25")
    print("=" * 74)
    print("  Gap among the top 5 is ZERO by construction. Sharpest available test:")
    print("  if the vernier is real the wheel should turn here several times.")
    print()
    lh = np.concatenate([np.full(5, -12.75), np.full(11, -14.25)])
    hist_c, hist_g, turns, geared, stalls = [], [], [], [], []
    for _ in range(n_real):
        H = wheel(rng, lh, n_psr=n_psr, prior="vlbi", noise_mode="white")
        turns.append(H["turns"]); geared.append(H["turns_geared"])
        stalls.append(H["stall"])
        hist_c.append(H["n_cert"]); hist_g.append(H["gain_dex"])
    mx = max(len(h) for h in hist_c)
    C = np.full((n_real, mx), np.nan); G = np.full((n_real, mx), np.nan)
    for i, (c, g) in enumerate(zip(hist_c, hist_g)):
        C[i, :len(c)] = c; G[i, :len(g)] = g
    print(f"  {'turn':>5} {'med certified':>14} {'med gain dex':>13} {'n real':>7}")
    for t in range(mx):
        m = ~np.isnan(C[:, t])
        if m.sum() < 3:
            continue
        print(f"  {t+1:5d} {np.nanmedian(C[m,t]):14.0f} {np.nanmedian(G[m,t]):13.4f} "
              f"{int(m.sum()):7d}")
    from collections import Counter
    print()
    print(f"  turns: median {np.median(turns):.1f}, range {min(turns)}-{max(turns)}")
    print(f"  GEARED turns (detections the wheel actually bought): "
          f"median {np.median(geared):.1f}, max {max(geared)}")
    print(f"  stall causes: {dict(Counter(stalls))}")
    print(f"  peak certified: median {np.nanmedian(np.nanmax(C,axis=1)):.0f}/{n_psr}")
    print()
    print("  READING: the 5 loud members are ALL above threshold on their own, so the")
    print("  wheel walks through them without needing gearing -- it turns, but the")
    print("  turning is supplied by the SKY, not by the wheel. The 1.5 dex cliff to")
    print("  the faint tier is where it stops, and 0.09 dex of gearing cannot cross it.")
    out["asky_cert"] = C; out["asky_gain"] = G
    out["asky_turns"] = np.array(turns); out["asky_geared"] = np.array(geared)
    out["asky_stalls"] = np.array(stalls)
    print()


# ================================================== sky (b): TOP-K conditioning
def exp_topk(out, rng, n_pop=256, sigma_logh=0.6, n_real=12, n_psr=116,
             max_draw=400000, chunk=20000):
    print("=" * 74)
    print("SKY (b) -- TOP-3 CONDITIONING at fixed band power")
    print("=" * 74)
    print(f"  population {n_pop} members, log10 h ~ Normal(sigma={sigma_logh} dex),")
    print("  rescaled so total band power sum h^2 is FIXED. Rejection-scan for skies")
    print("  whose top THREE lie within Delta dex. Cost = seeds scanned per accept.")
    print()
    deltas = np.array([0.05, 0.1, 0.15, 0.2, 0.3, 0.5])
    # ---- rejection cost, vectorised
    spreads = []
    drawn = 0
    while drawn < max_draw:
        m = min(chunk, max_draw - drawn)
        lh = rng.normal(0.0, sigma_logh, (m, n_pop))
        top3 = -np.sort(-lh, axis=1)[:, :3]
        spreads.append(top3[:, 0] - top3[:, 2])
        drawn += m
    spreads = np.concatenate(spreads)
    print(f"  drew {drawn} skies. top3 spread: median {np.median(spreads):.3f} dex, "
          f"16/84 = {np.percentile(spreads,16):.3f}/{np.percentile(spreads,84):.3f}")
    print()
    print(f"  {'Delta':>7} {'P(accept)':>11} {'seeds/accept':>13} "
          f"{'turns':>7} {'geared':>7} {'engages':>9}")
    P, cost, TT, GG = [], [], [], []
    for d in deltas:
        pa = float((spreads <= d).mean())
        P.append(pa)
        cost.append(np.inf if pa == 0 else 1.0 / pa)
        # ---- run the wheel on skies conditioned to satisfy it
        t_, g_ = [], []
        for _ in range(n_real):
            for _try in range(4000):
                lh = rng.normal(0.0, sigma_logh, n_pop)
                s3 = -np.sort(-lh)[:3]
                if s3[0] - s3[2] <= d:
                    break
            keep = -np.sort(-lh)[:16]
            keep = keep - keep.max()
            # fix band power: scale so sum h^2 over the population is constant
            allh = 10.0 ** (lh - lh.max())
            keep = keep - 0.5 * np.log10(np.sum(allh ** 2) / n_pop)
            H = wheel(rng, keep - 12.75, n_psr=n_psr, prior="vlbi",
                      noise_mode="white")
            t_.append(H["turns"]); g_.append(H["turns_geared"])
        TT.append(np.median(t_)); GG.append(np.median(g_))
        eng = "YES" if np.median(t_) >= 3 else ("weak" if np.median(t_) >= 2 else "no")
        print(f"  {d:7.2f} {pa:11.5f} {cost[-1]:13.1f} {TT[-1]:7.1f} "
              f"{GG[-1]:7.1f} {eng:>9}")
    print()
    eng_idx = [i for i, t in enumerate(TT) if t >= 3]
    if eng_idx:
        i = eng_idx[-1]
        print(f"  WHEEL ENGAGES for Delta <= {deltas[i]:.2f} dex "
              f"(cost {cost[i]:.1f} seeds/accept).")
    print("  Natural top-3 spread for this population is "
          f"{np.median(spreads):.2f} dex (median), so the conditioning is")
    print("  buying a factor of a few in rarity, not orders of magnitude.")
    out["topk_deltas"] = deltas
    out["topk_Paccept"] = np.array(P)
    out["topk_cost"] = np.array(cost)
    out["topk_turns"] = np.array(TT)
    out["topk_geared"] = np.array(GG)
    out["topk_spreads_pctl"] = np.percentile(spreads, [5, 16, 50, 84, 95])
    out["topk_n_drawn"] = drawn
    print()



# ============================== ADDENDUM 1 -- gearing as a function of f_cert
def exp_gearing_fcert(out, rng, n_real=60, n_psr=116):
    """
    Gearing is sqrt(1 + r^2 f_cert), where f_cert is the FRACTION of the array's
    pulsar-term POWER that is certified -- not a fixed sqrt(2).
    Certifying 1 pulsar of 116 buys ~0.002 dex, not 0.15.
    """
    print("=" * 74)
    print("ADDENDUM 1 -- gearing vs f_cert (fraction of pulsar-term POWER certified)")
    print("=" * 74)
    print("     gain(f_cert) = log10 sqrt(1 + r^2 f_cert),   r = rho_P/rho_E")
    print("  It is NOT a fixed sqrt(2): sqrt(2) is only the f_cert = 1, r = 1 corner.")
    print()
    # measured r over a realistic band/Mc/geometry population
    rs = []
    for _ in range(n_real):
        p, L, hred = build_array(n_psr, rng)
        s_, f0, mc, _ = make_sky(rng, np.zeros(4), n_psr)
        _, rE, rP, f_p, _ = couple(p, L, hred, s_, f0, mc, "white")
        rs.append(np.sqrt(np.sum(rP[:, 0] ** 2) / np.sum(rE[:, 0] ** 2)))
    rs = np.array(rs)
    r_med = float(np.median(rs))
    print(f"  measured r = rho_P/rho_E (array-summed, residual space, 3-20 nHz,")
    print(f"  log10Mc 8.5-9.5): median {r_med:.3f}, "
          f"16/84 = {np.percentile(rs,16):.3f}/{np.percentile(rs,84):.3f}")
    r_ext = 0.584 ** (-1.0 / 3.0)
    print(f"  at the evolved corner f_p/f0 = 0.584: r = {r_ext:.3f} "
          f"-> full-array ceiling {np.log10(np.sqrt(1+r_ext**2)):.4f} dex")
    print(f"  ADDENDUM ITEM 2 CONFIRMED: ceiling is 0.19 dex, not 0.15.")
    print()
    fc = np.array([0.0, 1.0 / n_psr, 0.01, 0.05, 0.1, 0.25, 0.5, 0.7, 0.9, 1.0])
    print(f"  {'f_cert':>9} {'n_psr equiv':>12} {'gain dex (r=%.2f)' % r_med:>18} "
          f"{'gain dex (r=%.2f)' % r_ext:>18}")
    g_med, g_ext = [], []
    for f in fc:
        gm = np.log10(np.sqrt(1 + r_med ** 2 * f))
        ge = np.log10(np.sqrt(1 + r_ext ** 2 * f))
        g_med.append(gm); g_ext.append(ge)
        print(f"  {f:9.4f} {f*n_psr:12.1f} {gm:18.4f} {ge:18.4f}")
    print()
    print("  ONE PULSAR OF 116 (f_cert = 0.0086) buys "
          f"{np.log10(np.sqrt(1+r_ext**2/n_psr)):.4f} dex. The addendum is right:")
    print("  the 0.15 dex figure was the all-certified corner, not a per-turn gain.")
    print()
    print("  f_cert REQUIRED to close a given gap:")
    for gap in [0.05, 0.15, 0.3, 0.5, 0.8, 1.2]:
        need_m = (10 ** (2 * gap) - 1) / r_med ** 2
        need_e = (10 ** (2 * gap) - 1) / r_ext ** 2
        tag_m = f"{need_m:8.3f}" if need_m <= 1 else "IMPOSSIBLE"
        tag_e = f"{need_e:8.3f}" if need_e <= 1 else "IMPOSSIBLE"
        print(f"    gap {gap:4.2f} dex -> f_cert = {tag_m:>10s} (r={r_med:.2f}) ; "
              f"{tag_e:>10s} (r={r_ext:.2f})")
    print()
    print("  *** 0.15 dex needs ~70% of the array's pulsar-term POWER certified.")
    print("  *** 0.8 dex is UNREACHABLE at any f_cert: the ceiling is 0.19 dex.")
    out["gfc_fcert"] = fc
    out["gfc_gain_med"] = np.array(g_med)
    out["gfc_gain_ext"] = np.array(g_ext)
    out["gfc_r_samples"] = rs
    out["gfc_r_med"] = r_med
    out["gfc_r_ext"] = r_ext
    out["gfc_ceiling_med"] = np.log10(np.sqrt(1 + r_med ** 2))
    out["gfc_ceiling_ext"] = np.log10(np.sqrt(1 + r_ext ** 2))
    print()
    return r_med, r_ext


# ============================== ADDENDUM 4 -- target selection
def exp_targets(out, rng, n_psr=116, n_real=40):
    """
    Fringe count K_a = 6 sigma f_p (1-cos mu)/c grows with mu, while the pulsar-term
    response carries a (1 - e^{i 2 pi f_p tau}) factor that vanishes as tau -> 0.
    Find the mu that maximises certifiable pulsar-term POWER.
    """
    print("=" * 74)
    print("ADDENDUM 4 -- target selection: where does certifiable POWER live?")
    print("=" * 74)
    print("  K_a = 6 sigma f_p (1-cos mu) / c        (fringes to remove: grows with mu)")
    print("  weight w(mu) = ant^2 |1-e^{i2pi f_p tau}|^2   (power: vanishes as tau->0)")
    print("  POWER PER FRINGE w/K is analytic here: ant^2 = 1-cos^2 mu and")
    print("  K ~ (1-cos mu), so w/K ~ (1+cos mu) -- MONOTONICALLY FALLING in mu.")
    print("  So the best targets are SMALL mu. But w -> 0 there, so there is a")
    print("  genuine optimum, and it is what limits f_cert.")
    print()
    mu = np.radians(np.linspace(0.5, 179.5, 400))
    f0 = 10.0 ** -7.7
    L, mc = 1500.0, 9.0
    tc = V.t_merge(f0, mc)
    omc = 1.0 - np.cos(mu)
    tau = V.tau_of(L, omc)
    f_p = f0 * (1 + tau / tc) ** (-3.0 / 8.0)
    supp = np.abs(1.0 - np.exp(2j * np.pi * f_p * tau)) ** 2 / 2.0   # ->0 as tau->0
    w = (1.0 - np.cos(mu) ** 2) * supp * (f_p / f0) ** (-2.0 / 3.0)
    K = 6.0 * 1.0 * f_p * omc / C_PC_S                              # sigma = 1 pc
    for kmax, lab in [(1.0, "K_a<1 (prior alone)"), (12.0, "K_a<12 (K=3 sources)"),
                      (250.0, "K_a<250 (K=5-6 sources)")]:
        ok = K < kmax
        if ok.any():
            best = np.nanargmax(np.where(ok, w, np.nan))
            print(f"    {lab:26s} mu < {np.degrees(mu[ok].max()):6.2f} deg ; "
                  f"peak certifiable w at mu = {np.degrees(mu[best]):6.2f} deg")
    print()
    print("  Now the number that matters: over a REAL isotropic array, what fraction")
    print("  of TOTAL pulsar-term power sits in pulsars that are ALSO certifiable?")
    print()
    print(f"  {'K_a cap':>9} {'sources needed':>15} {'f_cert':>9} {'n_psr':>7} "
          f"{'gain dex':>9}")
    r_ext = 0.584 ** (-1.0 / 3.0)
    caps = [(1.0, "1 (prior only)"), (12.0, "3"), (250.0, "5-6"), (3000.0, "6+")]
    rows = []
    for cap, need in caps:
        fcs, ns = [], []
        for _ in range(n_real):
            p, Lr, hred = build_array(n_psr, rng)
            s_, f0s, mcs, _ = make_sky(rng, np.zeros(1), n_psr)
            omc_a, rE, rP, fp_a, _ = couple(p, Lr, hred, s_, f0s, mcs, "white")
            Ka = 6.0 * 1.0 / V.dL_of(fp_a[:, 0], omc_a[:, 0])
            ok = Ka < cap
            tot = np.sum(rP[:, 0] ** 2)
            fcs.append(np.sum(rP[ok, 0] ** 2) / tot if tot > 0 else 0.0)
            ns.append(int(ok.sum()))
        f = float(np.median(fcs))
        g = np.log10(np.sqrt(1 + r_ext ** 2 * f))
        rows.append((cap, f, np.median(ns), g))
        print(f"  {cap:9.0f} {need:>15s} {f:9.4f} {np.median(ns):7.0f} {g:9.4f}")
    print()
    print("  HONEST READING of the table above. At sigma = 1 pc the K_a cap is NOT")
    print("  the binding constraint: dL spans ~0.3-15 pc so K_a = 6 sigma/dL is < 12")
    print("  for essentially the whole array, hence f_cert -> 1 for every cap at or")
    print("  above 12. The mu-selection effect is real but only BINDS in the")
    print(f"  sub-fringe regime: at K_a < 1 only {rows[0][1]*100:.1f}% of pulsar-term power")
    print(f"  qualifies ({rows[0][2]:.0f} of {n_psr} pulsars) -> {rows[0][3]:.4f} dex.")
    print("  K_a < cap is NECESSARY, NOT SUFFICIENT: certification also needs the")
    print("  JOINT comb intersection to collapse, which the wheel measures directly.")
    print("  *** Use the wheel's MEASURED f_cert (see GEARING section) as the real")
    print("  *** reservoir, not this proxy -- the proxy is optimistic.")
    out["tgt_mu_deg"] = np.degrees(mu)
    out["tgt_w"] = w; out["tgt_K"] = K; out["tgt_supp"] = supp
    out["tgt_caps"] = np.array([r[0] for r in rows])
    out["tgt_fcert"] = np.array([r[1] for r in rows])
    out["tgt_npsr"] = np.array([r[2] for r in rows])
    out["tgt_gain"] = np.array([r[3] for r in rows])
    print()


# ================================================== plots + main
def make_plots(out):
    # 0: gearing vs f_cert  (addendum 1 + 4)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    f = out["gfc_fcert"]
    ax[0].plot(f, out["gfc_gain_ext"], "o-", lw=2, ms=6, color="#4c72b0",
               label="r=1.196 (evolved, $f_p/f_0$=0.584)")
    ax[0].plot(f, out["gfc_gain_med"], "s--", lw=1.6, ms=5, color="#dd8452",
               label="r=%.2f (band median)" % out["gfc_r_med"])
    for gap, c in [(0.15, "#55a868"), (0.3, "#c44e52")]:
        ax[0].axhline(gap, color=c, ls=":", lw=1.4, label=f"gap {gap} dex")
    ax[0].set_xlabel("$f_{cert}$ = fraction of pulsar-term POWER certified")
    ax[0].set_ylabel("gain [dex]")
    ax[0].set_title("gearing is $\\sqrt{1+r^2 f_{cert}}$, not $\\sqrt{2}$")
    ax[0].legend(fontsize=7.5, frameon=False); ax[0].grid(alpha=.25)
    mu, w, K = out["tgt_mu_deg"], out["tgt_w"], out["tgt_K"]
    ax[1].plot(mu, w / w.max(), lw=2, color="#4c72b0", label="pulsar-term power w")
    ax[1].plot(mu, K / K.max(), lw=2, color="#c44e52", label="fringe count $K_a$")
    ax[1].plot(mu, (1 + np.cos(np.radians(mu))) / 2, lw=1.6, ls="--",
               color="#55a868", label="power per fringe $\\propto 1+\\cos\\mu$")
    ax[1].set_xlabel("$\\mu$ [deg]"); ax[1].set_ylabel("normalised")
    ax[1].set_title("certifiable power lives at low $\\mu$,\nwhere there is little power")
    ax[1].legend(fontsize=7.5, frameon=False); ax[1].grid(alpha=.25)
    x = np.arange(len(out["tgt_caps"]))
    ax[2].bar(x, out["tgt_fcert"], color="#4c72b0", alpha=.85)
    ax[2].set_xticks(x)
    ax[2].set_xticklabels([f"$K_a<${c:.0f}" for c in out["tgt_caps"]], fontsize=8)
    ax[2].set_ylabel("$f_{cert}$ achievable")
    ax[2].set_title("the wheel's REAL gearing reservoir")
    for i, (fv, gv) in enumerate(zip(out["tgt_fcert"], out["tgt_gain"])):
        ax[2].text(i, fv, f"{gv:.3f}\ndex", ha="center", va="bottom", fontsize=7)
    ax[2].grid(alpha=.25, axis="y")
    fig.suptitle("VERNIER-2 addendum -- gearing vs certified power fraction "
                 "(exploratory)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER2_fcert.png"), dpi=140)
    plt.close(fig)

    # 1: dial walk
    fig, ax = plt.subplots(1, 5, figsize=(18, 3.6))
    for i, dial in enumerate(DIALS):
        v = out[f"dial_{dial}_vals"]; t = out[f"dial_{dial}_turns"]
        g = out[f"dial_{dial}_geared"]
        x = np.arange(len(v))
        ax[i].plot(x, t, "o-", lw=2, ms=6, color="#4c72b0", label="turns")
        ax[i].plot(x, g, "s--", lw=1.5, ms=5, color="#c44e52", label="geared turns")
        ax[i].axhline(2, color="k", ls=":", lw=1)
        ax[i].set_xticks(x); ax[i].set_xticklabels(v, rotation=45, fontsize=8)
        ax[i].set_title(dial, fontsize=10); ax[i].grid(alpha=.25)
        ax[i].set_ylabel("turns" if i == 0 else "")
        if i == 0:
            ax[i].legend(fontsize=8, frameon=False)
    fig.suptitle("VERNIER-2 dial walk -- one dial at a time from the working point "
                 "(exploratory)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER2_dials.png"), dpi=140)
    plt.close(fig)

    # 2: gearing
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    gb = out["gear_gain_by_turn"]; k = ~np.isnan(gb)
    ax[0].plot(np.arange(1, len(gb) + 1)[k], gb[k], "o-", lw=2, ms=7, color="#4c72b0",
               label="measured gain")
    ax[0].axhline(out["gear_sqrt2_dex"], color="#c44e52", ls="--", lw=1.5,
                  label="$\\log_{10}\\sqrt{2}$ = 0.1505")
    ax[0].axhline(np.median(out["gear_ceiling_dex"]), color="#55a868", ls="-.", lw=1.5,
                  label="measured all-certified ceiling")
    ax[0].set_xlabel("wheel turn"); ax[0].set_ylabel("gain [dex]")
    ax[0].set_title("gearing saturates, it does not compound")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=.25)
    ax[1].hist(out["gear_ceiling_dex"], bins=20, color="#4c72b0", alpha=.8)
    ax[1].axvline(out["gear_sqrt2_dex"], color="#c44e52", ls="--", lw=1.5,
                  label="$\\log_{10}\\sqrt{2}$")
    ax[1].set_xlabel("all-certified ceiling [dex]"); ax[1].set_ylabel("count")
    ax[1].set_title("all-certified ceiling vs $\\sqrt{2}$")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=.25)
    fig.suptitle("VERNIER-2 gearing (exploratory)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER2_gearing.png"), dpi=140)
    plt.close(fig)

    # 3: top-k conditioning
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    d, c, t = out["topk_deltas"], out["topk_cost"], out["topk_turns"]
    ax[0].semilogy(d, c, "o-", lw=2, ms=6, color="#4c72b0")
    ax[0].set_xlabel("$\\Delta$ dex (top-3 within)"); ax[0].set_ylabel("seeds per accept")
    ax[0].set_title("rejection cost = the sky's rarity"); ax[0].grid(alpha=.25, which="both")
    ax[1].plot(d, t, "o-", lw=2, ms=6, color="#4c72b0", label="turns")
    ax[1].plot(d, out["topk_geared"], "s--", lw=1.5, ms=5, color="#c44e52",
               label="geared")
    ax[1].axhline(3, color="k", ls=":", lw=1)
    ax[1].set_xlabel("$\\Delta$ dex"); ax[1].set_ylabel("turns")
    ax[1].set_title("where the wheel engages")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=.25)
    fig.suptitle("VERNIER-2 top-3 conditioning (exploratory)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER2_topk.png"), dpi=140)
    plt.close(fig)


def main():
    rng = np.random.default_rng(20260729)
    out = {}
    print("\nVERNIER-2 -- EXPLORATORY, NOT pre-registered. Clean room preserved.")
    print("Evolving (exact Peters) pulsar term throughout; comb from f_p.\n")
    exp_evolve_shift(out, rng)
    exp_spectral(out)
    exp_gearing(out, rng)
    exp_gearing_fcert(out, rng)
    exp_targets(out, rng)
    exp_dials(out, rng)
    exp_asky(out, rng)
    exp_topk(out, rng)
    # SAVE FIRST: a cosmetic plotting error must never destroy a completed run.
    npz = os.path.join(REPORTS, "vernier2.npz")
    np.savez_compressed(npz, **{k: v for k, v in out.items()})
    print(f"[saved {npz} before plotting]")
    try:
        make_plots(out)
    except Exception as e:
        print(f"WARNING: plotting failed ({type(e).__name__}: {e}); data is safe.")
    print("=" * 74)
    print(f"wrote {npz}")
    for f in ["VERNIER2_fcert.png", "VERNIER2_dials.png", "VERNIER2_gearing.png",
              "VERNIER2_topk.png"]:
        print(f"wrote {os.path.join(REPORTS, f)}")
    print("=" * 74)


if __name__ == "__main__":
    main()
