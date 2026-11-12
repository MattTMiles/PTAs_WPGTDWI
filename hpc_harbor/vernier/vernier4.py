#!/usr/bin/env python3
"""
VERNIER-4 -- realistic sky, generous array. EXPLORATORY, NOT pre-registered.

Proof of principle, not a forecast. Clean room: imports only vernier_toy/vernier2
(my own anchored code).

THE CENSUS, DERIVED (not chosen)
--------------------------------
Sources per frequency bin (bin width df = 1/T):
    N_c(f) = N_YR * (f/f_yr)^(-11/3),  N_YR ~ 1000 at f_yr = 1/yr
Strain distribution at fixed f: sources uniform in Euclidean volume, h ~ 1/D, so
    N(>h) ~ h^-3   =>   dN/dh ~ h^-4   =>   rank-k strain  h_k^2 ~ k^(-2/3)
    sum_{k=1..N} k^(-2/3) ~ 3 N^(1/3)
The residual power in a bin equals the GWB PSD times the bin width:
    sum_k R_k^2 = 2 S_gwb(f) df,    R_k = h_k / (2 pi f)   (residual, not strain)
so the bright end is fixed by the GWB amplitude with NO free normalisation:
    R_1^2 = 2 S_gwb(f) df / (3 N_c^(1/3))

Only the top N_BRIGHT sources per bin are ever handled individually; everything
fainter is carried as a Gaussian stochastic background holding the residual power.
Cost is independent of N_c, so N_c = 1e7 costs the same as N_c = 1e3.

THE TWO NUMBERS VERNIER-3 HAD TO ASSUME, NOW DERIVED PER BIN
    C(f)     = S_gwb(f) / S_white           confusion-to-white ratio
    f_rem(f) = (N_sub / N_c(f))^(1/3)       fraction of confusion power removable
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V
import vernier2 as W

YR, REPORTS = V.YR, V.REPORTS
F_YR = 1.0 / YR
V.PTERM_STATIC = False

A_GWB = 2.4e-15          # NANOGrav-15yr-like characteristic strain at f_yr
N_YR = 1000.0            # sources per bin at f = 1/yr
T_OBS = 15.0 * YR
DF = 1.0 / T_OBS
N_BRIGHT = 300           # sources simulated individually per bin
KAPPA = 1.0 / 3.0        # sky/polarisation-averaged antenna power factor
DET_SNR = 8.0


def S_gwb(f):
    """GWB residual PSD, s^3.  S = A^2/(12 pi^2 f_yr^3) (f/f_yr)^(-13/3)."""
    return (A_GWB ** 2 / (12 * np.pi ** 2 * F_YR ** 3)) * (np.asarray(f) / F_YR) ** (-13.0 / 3.0)


def S_white(sigma_toa, cadence_days=14.0):
    """White residual PSD, s^3.  S = 2 sigma^2 dt."""
    return 2.0 * sigma_toa ** 2 * (cadence_days * 86400.0)


def N_c(f):
    """Sources per bin of width 1/T."""
    return N_YR * (np.asarray(f) / F_YR) ** (-11.0 / 3.0)


def C_of_f(f, sigma_toa):
    return S_gwb(f) / S_white(sigma_toa)


def f_removable(f, n_sub):
    """Fraction of bin confusion POWER removed by subtracting the top n_sub."""
    n = N_c(f)
    return np.minimum((np.minimum(n_sub, n) / n) ** (1.0 / 3.0), 1.0)


def rank_R2(f, kmax=N_BRIGHT):
    """Residual amplitude^2 of the top kmax sources in the bin, s^2."""
    n = N_c(f)
    tot = 2.0 * S_gwb(f) * DF                      # sum of R^2 over the WHOLE bin
    k = np.arange(1, int(min(kmax, max(n, 1))) + 1)
    return tot * k ** (-2.0 / 3.0) / (3.0 * n ** (1.0 / 3.0)), tot


def snr_of(R2, f, n_psr, sigma_toa, S_conf):
    """Matched-filter SNR^2 of a source of residual amplitude^2 R2."""
    S_tot = S_white(sigma_toa) + S_conf
    return KAPPA * n_psr * R2 * T_OBS / (2.0 * S_tot)


# ============================================================ deliverable 1
def exp_C_and_f(out, sigma_list=(100e-9, 30e-9, 10e-9), n_sub=(10, 300)):
    print("=" * 78)
    print("DELIVERABLE 1 -- C(f) and f_rem(f), DERIVED from the census + array")
    print("=" * 78)
    print(f"  A_gwb = {A_GWB:.1e}, N_c(f_yr) = {N_YR:.0f}/bin, bin width 1/T = "
          f"{DF:.3e} Hz (T = 15 yr)")
    print()
    fg = np.array([2, 3, 5, 8, 10, 15, 20, 32, 50, 80]) * 1e-9
    print(f"  {'f nHz':>7} {'N_c/bin':>11} {'C(100ns)':>10} {'C(30ns)':>9} "
          f"{'C(10ns)':>9} {'f_rem(10)':>10} {'f_rem(300)':>11} {'R_1 ns':>9}")
    rows = []
    for f in fg:
        n = N_c(f)
        R2, _ = rank_R2(f)
        r1 = np.sqrt(R2[0]) * 1e9
        row = [f, n, C_of_f(f, 100e-9), C_of_f(f, 30e-9), C_of_f(f, 10e-9),
               f_removable(f, 10), f_removable(f, 300), r1]
        rows.append(row)
        print(f"  {f*1e9:7.1f} {n:11.3g} {row[2]:10.3g} {row[3]:9.3g} "
              f"{row[4]:9.3g} {row[5]:10.4f} {row[6]:11.4f} {r1:9.3f}")
    rows = np.array(rows)
    print()
    print("  *** THE SCISSORS, stated correctly. AT FIXED TIMING PRECISION, C(f) ~")
    print("  *** f^(-13/3) FALLS with frequency while f_rem(f) ~ N_c^(-1/3) ~ f^(+11/9)")
    print("  *** RISES. Their product goes as f^(-28/9), so M2 has a single operating")
    print("  *** window and it sits at the GWB / WHITE-NOISE CROSSOVER (C ~ 1): below")
    print("  *** it the sea is huge but 1e5-1e7 sources deep and so irreducible; above")
    print("  *** it few sources hold the power but there is no confusion left to remove.")
    print("  *** CRUCIALLY the crossover is NOT fixed -- better timing lowers S_white,")
    print("  *** raising C at every f and pushing the window UP in frequency, to where")
    print("  *** f_rem is large. That is why M2 is array-dependent, not doomed.")
    for s_ in [100e-9, 30e-9, 10e-9]:
        ff = np.logspace(np.log10(2e-9), np.log10(1e-7), 400)
        cc = C_of_f(ff, s_)
        xo = ff[np.argmin(np.abs(cc - 1.0))]
        print(f"      sigma = {s_*1e9:5.0f} ns -> GWB/white crossover at "
              f"{xo*1e9:5.1f} nHz")
    print()
    # M2 gain, exact: 0.5 log10[(1+C)/(1+C(1-f_rem))]
    print(f"  M2 gain (dex) = 0.5 log10[(1+C)/(1+C(1-f_rem))]")
    print(f"  {'f nHz':>7} " + " ".join(f"{'s=%gns,N=%d' % (s*1e9, ns):>14}"
                                        for s in sigma_list for ns in n_sub))
    g = np.zeros((len(fg), len(sigma_list) * len(n_sub)))
    for i, f in enumerate(fg):
        cells = []
        j = 0
        for s in sigma_list:
            for ns in n_sub:
                C = C_of_f(f, s); fr = f_removable(f, ns)
                gv = 0.5 * np.log10((1 + C) / (1 + C * (1 - fr)))
                g[i, j] = gv; j += 1
                cells.append(f"{gv:14.4f}")
        print(f"  {f*1e9:7.1f} " + " ".join(cells))
    print()
    print(f"  MAX M2 gain anywhere in this table: {g.max():.4f} dex "
          f"(at f = {fg[np.unravel_index(g.argmax(), g.shape)[0]]*1e9:.0f} nHz)")
    print("  Peak M2 gain and the frequency it sits at, per timing precision:")
    j = 0
    for s_ in sigma_list:
        for ns in n_sub:
            i = int(np.argmax(g[:, j]))
            print(f"      sigma={s_*1e9:5.0f} ns, subtract {ns:3d}: max "
                  f"{g[:, j].max():.4f} dex at {fg[i]*1e9:5.1f} nHz")
            j += 1
    print()
    print("  VERNIER-3 ASSUMED C = 3, f_rem = 0.83 -> 0.199 dex.")
    print("  Derived: on a 100 ns array the best M2 can do anywhere in the band is")
    print(f"  {g[:, 1].max():.4f} dex -- {0.199/max(g[:, 1].max(),1e-6):.0f}x smaller, i.e. M2 COLLAPSES.")
    print(f"  On a 10 ns array it recovers to {g[:, 5].max():.4f} dex at the top of the band.")
    out["cf_f"] = fg
    out["cf_rows"] = rows
    out["cf_m2gain"] = g
    out["cf_m2gain_max"] = g.max()
    print()
    return fg, rows, g


# ============================================================ the wheel, realistic sky
def wheel_bin(rng, f_bin, n_psr, sigma_toa, arm="both", n_cert_src=12,
              max_turns=40, grid_cap=64000):
    """
    Run the wheel inside ONE frequency bin on the derived census.
    Only the top N_BRIGHT sources exist individually; the rest are Gaussian sea.
    """
    R2, tot_R2 = rank_R2(f_bin)
    nsrc = len(R2)
    p, L, has_red = W.build_array(n_psr, rng)
    s = V.iso_dirs(nsrc, rng)
    mc = rng.uniform(8.5, 9.5, nsrc)
    f0 = np.full(nsrc, f_bin)
    tc = V.t_merge(f0, mc)
    cmu = p @ s.T
    omc = np.clip(1.0 - cmu, 1e-3, None)
    ant = np.sqrt(np.clip(1.0 - cmu ** 2, 1e-6, None))
    f_p = f0[None, :] * (1.0 + V.tau_of(L[:, None], omc) / tc[None, :]) ** (-3.0 / 8.0)
    sig_arr = np.full(n_psr, 1.0)                      # 1 pc priors
    Kaj = np.maximum(6.0 * sig_arr[:, None] / V.dL_of(f_p, omc), 1.0)
    Sw = S_white(sigma_toa)

    def S_conf(det):
        """Gaussian sea power left after subtracting the detected bright sources."""
        if arm in ("M1", "none") or not det:
            return S_gwb(f_bin)
        removed = R2[list(det)].sum() / tot_R2
        return S_gwb(f_bin) * (1.0 - removed)

    w = ant ** 2 / np.sum(ant ** 2, axis=0, keepdims=True)      # per-pulsar split

    def rho2(det):
        """Per-pulsar rho^2; the array total is KAPPA n_psr R^2 T / (2 S_tot)."""
        St = Sw + S_conf(det)
        base = KAPPA * n_psr * R2 * T_OBS / (2.0 * St)
        return base[None, :] * w

    det, cert = [], np.zeros(n_psr, bool)
    H = {"detected": [], "n_cert": [], "f_cert": [], "snr": [], "gain_dex": [],
         "S_conf_frac": []}
    while len(det) < max_turns:
        r2 = rho2(det)                                  # (n_psr, nsrc)
        rE = np.sqrt(np.maximum(r2, 0))
        rP = rE * (f_p / f0[None, :]) ** (-1.0 / 3.0)
        pm = np.maximum(rP ** 2 - 2.0 * np.log(Kaj), 0.0)
        avail = [j for j in range(nsrc) if j not in det]
        if not avail:
            break
        cand = np.array(avail)
        if arm in ("M2", "none"):
            snr = np.sqrt(np.sum(rE[:, cand] ** 2, axis=0)
                          + np.sum(pm[:, cand], axis=0))
        else:
            snr = np.sqrt(np.sum(rE[:, cand] ** 2, axis=0)
                          + np.sum(rP[cert][:, cand] ** 2, axis=0)
                          + np.sum(pm[~cert][:, cand], axis=0))
        k = int(np.argmax(snr))
        best, sb = int(cand[k]), float(snr[k])
        # no-gearing reference: full sea, nothing certified
        r2_0 = KAPPA * n_psr * R2 * T_OBS / (2.0 * (Sw + S_gwb(f_bin)))
        rE0 = np.sqrt(np.maximum(r2_0[None, :] * w, 0))
        rP0 = rE0 * (f_p / f0[None, :]) ** (-1.0 / 3.0)
        pm0 = np.maximum(rP0 ** 2 - 2.0 * np.log(Kaj), 0.0)
        s0 = np.sqrt(np.sum(rE0[:, best] ** 2) + np.sum(pm0[:, best]))
        tp = np.sum(rP[:, best] ** 2)
        H["snr"].append(sb); H["n_cert"].append(int(cert.sum()))
        H["f_cert"].append(float(np.sum(rP[cert, best] ** 2) / tp) if tp > 0 else 0.0)
        H["gain_dex"].append(float(np.log10(sb / s0)) if s0 > 0 else 0.0)
        H["S_conf_frac"].append(float(S_conf(det) / S_gwb(f_bin)))
        if sb < DET_SNR:
            H["detected"].append(-1)
            break
        det.append(best); H["detected"].append(best)
        if arm in ("M1", "both"):
            use = det[:n_cert_src]                      # loudest few carry the combs
            r2c = rho2(det)
            rPc = np.sqrt(np.maximum(r2c, 0)) * (f_p / f0[None, :]) ** (-1.0 / 3.0)
            new = np.zeros(n_psr, bool)
            for a in range(n_psr):
                srcs = [V.Src(f0[j], tc[j], omc[a, j], 1, rPc[a, j]) for j in use]
                new[a] = (V.scan(L[a], sig_arr[a], srcs, cap=grid_cap)[0] == 1)
            cert = new
    H = {k: np.array(v) for k, v in H.items()}
    H["turns"] = int((H["detected"] >= 0).sum())
    return H


def exp_band(out, rng, n_psr=200, sigmas=(100e-9, 30e-9, 10e-9), n_real=4):
    print("=" * 78)
    print("QUESTION 2 -- WHERE in the band does the wheel operate?")
    print("=" * 78)
    print(f"  {n_psr} pulsars, 1 pc priors, both mechanisms on, {n_real} realisations.")
    print("  'turns' = individually detected sources in that bin (SNR > 8).")
    print()
    fg = np.array([2, 3, 5, 8, 10, 15, 20, 32, 50, 80]) * 1e-9
    res = np.zeros((len(sigmas), len(fg)))
    snr1 = np.zeros_like(res)
    for i, sg in enumerate(sigmas):
        for j, f in enumerate(fg):
            T_, S_ = [], []
            for _ in range(n_real):
                H = wheel_bin(rng, f, n_psr, sg, arm="both")
                T_.append(H["turns"]); S_.append(H["snr"][0] if len(H["snr"]) else 0)
            res[i, j] = np.median(T_); snr1[i, j] = np.median(S_)
        print(f"  sigma={sg*1e9:5.0f} ns  turns/bin: "
              + " ".join(f"{v:4.0f}" for v in res[i]))
        print(f"  {'':16s} SNR_1    : " + " ".join(f"{v:4.1f}" for v in snr1[i]))
        ok = np.where(res[i] >= 1)[0]
        print(f"  {'':16s} -> operates at f >= "
              f"{fg[ok[0]]*1e9:.0f} nHz" if ok.size else
              f"  {'':16s} -> NEVER operates anywhere in the band")
        print()
    out["band_f"] = fg
    out["band_sigmas"] = np.array(sigmas)
    out["band_turns"] = res
    out["band_snr1"] = snr1
    print()
    return fg, res, snr1



# ============================================================ Q1: which mechanism
def exp_mechanisms(out, rng, n_psr=5000, sigma_toa=3e-9, n_real=5):
    print("=" * 78)
    print("QUESTION 1 -- does the wheel still turn, and which mechanism wins?")
    print("=" * 78)
    print(f"  Run at a CLOSING array ({n_psr} pulsars, {sigma_toa*1e9:.0f} ns) so that")
    print("  something is detectable at all -- on today's array nothing is, and the")
    print("  M1/M2 question is moot. 1 pc priors, evolving pulsar terms.")
    print()
    fg = np.array([10, 20, 32, 50, 68]) * 1e-9
    arms = ["none", "M1", "M2", "both"]
    T = np.zeros((len(arms), len(fg)))
    G = np.zeros_like(T)
    FC = np.zeros_like(T)
    print(f"  {'f nHz':>7} " + " ".join(f"{a:>22s}" for a in arms))
    for j, f in enumerate(fg):
        cells = []
        for i, a in enumerate(arms):
            t_, g_, fc_ = [], [], []
            for _ in range(n_real):
                H = wheel_bin(rng, f, n_psr, sigma_toa, arm=a)
                t_.append(H["turns"])
                g_.append(np.median(H["gain_dex"]) if len(H["gain_dex"]) else 0.0)
                fc_.append(H["f_cert"][-1] if len(H["f_cert"]) else 0.0)
            T[i, j], G[i, j], FC[i, j] = np.median(t_), np.median(g_), np.median(fc_)
            cells.append(f"{T[i,j]:5.0f} trn {G[i,j]:+7.4f} dex")
        print(f"  {f*1e9:7.1f} " + " ".join(f"{c:>22s}" for c in cells))
    print()
    for j, f in enumerate(fg):
        i1, i2 = arms.index("M1"), arms.index("M2")
        eng = "M1" if G[i1, j] > G[i2, j] else ("M2" if G[i2, j] > G[i1, j] else "tie")
        print(f"    f={f*1e9:5.1f} nHz: M1={G[i1,j]:.4f}  M2={G[i2,j]:.4f}  "
              f"both={G[arms.index('both'),j]:.4f}  -> engine = {eng}  "
              f"(f_cert={FC[arms.index('both'),j]:.3f})")
    out["mech_f"] = fg; out["mech_arms"] = np.array(arms)
    out["mech_turns"] = T; out["mech_gain"] = G; out["mech_fcert"] = FC
    out["mech_npsr"] = n_psr; out["mech_sigma"] = sigma_toa
    print("  *** THE CERTIFICATION / DETECTION SCISSORS (new here).")
    print("  Detection is an ARRAY statistic: rho_array = sqrt(sum_a rho_a^2), so it")
    print("  grows as sqrt(N_psr). Certification is a PER-PULSAR statistic: VERNIER-1")
    print("  showed a pulsar needs its OWN pulsar-term rho_P >~ 2 for its comb to")
    print("  collapse. At fixed source strength rho_a ~ rho_array / sqrt(N_psr), so")
    print("  growing the array to reach detection DRIVES PER-PULSAR SNR DOWN.")
    print()
    print(f"  {'N_psr':>8} {'rho_array for detection':>24} "
          f"{'per-psr rho at detection':>26} {'rho_array needed to CERTIFY':>29}")
    for N in [50, 116, 500, 1000, 5000, 20000]:
        print(f"  {N:8d} {DET_SNR:24.1f} {DET_SNR/np.sqrt(N):26.3f} "
              f"{2.0*np.sqrt(N):29.1f}")
    print()
    print("  So a 5000-pulsar array detects at rho=8 but would need rho ~ 141 -- a")
    print("  source ~18x louder -- before its distances certify. M1 and detection")
    print("  pull in OPPOSITE directions in N_psr. That is why f_cert collapses to")
    print(f"  ~{out['mech_fcert'].max():.3f} here versus 0.2-0.9 in VERNIER-3's small-array toy.")
    print()

# ============================================================ what array closes it
def R1_squared(f):
    """Residual amplitude^2 of the LOUDEST source in each bin (vectorised in f)."""
    n = np.maximum(N_c(f), 1.0)
    tot = 2.0 * S_gwb(f) * DF
    return tot / (3.0 * n ** (1.0 / 3.0))


def snr1_analytic(f, n_psr, sigma_toa, confusion="uncorrelated", n_sub=0):
    """
    SNR of the LOUDEST source in the bin.

    confusion="uncorrelated": the sea averages down like white noise, 1/N_psr.
        OPTIMISTIC -- a real GWB is Hellings-Downs correlated between pulsars and
        therefore does NOT average away with more pulsars.
    confusion="correlated": the sea sets an N-independent SNR ceiling.
        PESSIMISTIC -- assumes none of it can be separated.
    The truth is between these; they bracket the capability statement.
    """
    R2 = R1_squared(f)
    Sw = S_white(sigma_toa)
    Sc = S_gwb(f) * (1.0 - f_removable(f, n_sub) if n_sub else 1.0)
    if confusion == "uncorrelated":
        return np.sqrt(KAPPA * n_psr * R2 * T_OBS / (2.0 * (Sw + Sc)))
    r2_w = KAPPA * n_psr * R2 * T_OBS / (2.0 * Sw)
    r2_c = KAPPA * 1.0 * R2 * T_OBS / (2.0 * Sc)          # no averaging
    return np.sqrt(1.0 / (1.0 / r2_w + 1.0 / r2_c))


def exp_array(out):
    print("=" * 78)
    print("QUESTION 3 -- WHAT ARRAY CLOSES THE LOOP (capability, not probability)")
    print("=" * 78)
    print("  Closure = the LOUDEST source in some bin reaches SNR 8 UNAIDED. That is")
    print("  the gate: until one source is detected there is nothing to subtract and")
    print("  nothing to certify, so the wheel cannot start. M1/M2 only matter after.")
    print()
    fg = np.logspace(np.log10(2e-9), np.log10(2e-7), 300)
    Ns = np.array([100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000])
    sg = np.array([100e-9, 30e-9, 10e-9, 3e-9, 1e-9])
    for model in ["uncorrelated", "correlated"]:
        print(f"  --- confusion model: {model.upper()}"
              + ("  (optimistic: sea averages as 1/N)" if model == "uncorrelated"
                 else "  (pessimistic: sea is an N-independent floor)"))
        print(f"  {'N_psr':>7} " + " ".join(f"{s*1e9:>8.0f}ns" for s in sg))
        best = np.zeros((Ns.size, sg.size))
        bf = np.zeros_like(best)
        for i, N in enumerate(Ns):
            cells = []
            for j, s in enumerate(sg):
                v = snr1_analytic(fg, N, s, confusion=model)
                best[i, j] = v.max(); bf[i, j] = fg[int(np.argmax(v))]
                cells.append(f"{v.max():10.2f}")
            print(f"  {N:7d} " + " ".join(cells))
        ok = best >= DET_SNR
        print(f"    peak-SNR frequency (nHz), same grid:")
        for i, N in enumerate(Ns):
            print(f"  {N:7d} " + " ".join(f"{bf[i,j]*1e9:10.1f}" for j in range(sg.size)))
        if ok.any():
            idx = np.argwhere(ok)
            # cheapest closure = fewest pulsars, then loosest timing
            i, j = idx[np.lexsort((-sg[idx[:, 1]], Ns[idx[:, 0]]))][0]
            print(f"    *** CLOSES at N_psr = {Ns[i]}, sigma = {sg[j]*1e9:.0f} ns "
                  f"(peak SNR {best[i,j]:.1f} at {bf[i,j]*1e9:.0f} nHz)")
        else:
            print(f"    *** NEVER CLOSES anywhere on this grid "
                  f"(max SNR {best.max():.2f} at N_psr={Ns[np.unravel_index(best.argmax(),best.shape)[0]]}, "
                  f"{sg[np.unravel_index(best.argmax(),best.shape)[1]]*1e9:.0f} ns)")
        print()
        out[f"arr_best_{model}"] = best
        out[f"arr_bestf_{model}"] = bf
    out["arr_Ns"] = Ns; out["arr_sigmas"] = sg
    # ---- closure restricted to the NOMINAL PTA band, where the population law is
    # measured rather than extrapolated.
    print("  --- CLOSURE RESTRICTED TO THE NOMINAL BAND (2-32 nHz)")
    print("  Above ~32 nHz the f^-11/3 law is an extrapolation and N_c -> 1, so a")
    print("  closure that only happens at 113 nHz is a statement about extrapolated")
    print("  population, not about the PTA band.")
    fb = fg[fg <= 32e-9]
    for model in ["uncorrelated", "correlated"]:
        bb = np.array([[snr1_analytic(fb, N, s_, confusion=model).max()
                        for s_ in sg] for N in Ns])
        okb = bb >= DET_SNR
        print(f"    {model:14s}: ", end="")
        if okb.any():
            idx = np.argwhere(okb)
            i, j = idx[np.lexsort((-sg[idx[:, 1]], Ns[idx[:, 0]]))][0]
            print(f"closes IN BAND at N_psr = {Ns[i]}, sigma = {sg[j]*1e9:.0f} ns "
                  f"(SNR {bb[i,j]:.1f})")
        else:
            print(f"NEVER closes in band (max SNR {bb.max():.2f} at "
                  f"N_psr={Ns[np.unravel_index(bb.argmax(),bb.shape)[0]]}, "
                  f"{sg[np.unravel_index(bb.argmax(),bb.shape)[1]]*1e9:.0f} ns)")
        out[f"arr_inband_{model}"] = bb
    print()
    print("  WHY THE TWO MODELS DIVERGE SO HARD. In this census the loudest source in")
    print("  a bin holds a fraction 1/(3 N_c^(1/3)) of the bin power, so its amplitude")
    print("  relative to the sea is 1/sqrt(3 N_c^(1/3)):")
    for f in [10e-9, 32e-9, 80e-9, 200e-9]:
        n = N_c(f)
        print(f"      f = {f*1e9:5.0f} nHz : N_c = {n:10.3g}, "
              f"h_1/h_sea = {1/np.sqrt(3*n**(1/3)):.3f}")
    print("  Until N_c ~ 1 the loudest source is BELOW the sea it sits in, so if the")
    print("  sea cannot be averaged down the source cannot be dug out at any N_psr.")
    print(f"  N_c = 1 only at f = {F_YR*(N_YR)**(3/11)*1e9:.0f} nHz, far above the band.")
    print()


def make_plots(out):
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.4))
    f = out["cf_f"] * 1e9
    rows = out["cf_rows"]
    ax[0].loglog(f, rows[:, 2], "o-", lw=2, ms=5, color="#c44e52", label="C, 100 ns")
    ax[0].loglog(f, rows[:, 3], "s-", lw=2, ms=5, color="#dd8452", label="C, 30 ns")
    ax[0].loglog(f, rows[:, 4], "^-", lw=2, ms=5, color="#4c72b0", label="C, 10 ns")
    ax[0].loglog(f, rows[:, 6], "d--", lw=2, ms=5, color="#55a868",
                 label="$f_{rem}$ (top 300)")
    ax[0].axhline(1.0, color="k", ls=":", lw=1.2)
    ax[0].set_xlabel("frequency [nHz]"); ax[0].set_ylabel("C  /  $f_{rem}$")
    ax[0].set_title("(1) derived: C falls, $f_{rem}$ rises")
    ax[0].legend(fontsize=7.5, frameon=False); ax[0].grid(alpha=.25, which="both")
    g = out["cf_m2gain"]
    for j, lab in enumerate(["100ns/10", "100ns/300", "30ns/10", "30ns/300",
                             "10ns/10", "10ns/300"]):
        ax[1].semilogx(f, g[:, j], "o-", ms=4, lw=1.6, label=lab,
                       color=plt.cm.viridis(j / 5))
    ax[1].axhline(0.199, color="#c44e52", ls="--", lw=1.4,
                  label="VERNIER-3 assumed 0.199")
    ax[1].set_xlabel("frequency [nHz]"); ax[1].set_ylabel("M2 gain [dex]")
    ax[1].set_title("(2) M2 has one narrow window")
    ax[1].legend(fontsize=7, frameon=False, ncol=2); ax[1].grid(alpha=.25)
    for m, ls, c in [("uncorrelated", "-", "#4c72b0"), ("correlated", "--", "#c44e52")]:
        b = out[f"arr_best_{m}"]
        for j in [0, 2, 4]:
            ax[2].loglog(out["arr_Ns"], b[:, j], ls, lw=1.7, marker="o", ms=4,
                         color=c, alpha=1.0 - 0.25 * (j / 4),
                         label=f"{m[:5]}, {out['arr_sigmas'][j]*1e9:.0f}ns")
    ax[2].axhline(DET_SNR, color="k", ls=":", lw=1.4, label="detection SNR 8")
    ax[2].set_xlabel("N pulsars"); ax[2].set_ylabel("SNR of loudest source")
    ax[2].set_title("(3) what array closes the loop")
    ax[2].legend(fontsize=6.5, frameon=False, ncol=2); ax[2].grid(alpha=.25, which="both")
    fig.suptitle("VERNIER-4 -- realistic sky, generous array (exploratory)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER4_realistic.png"), dpi=140)
    plt.close(fig)


def main():
    rng = np.random.default_rng(20260731)
    out = {}
    print("\nVERNIER-4 -- EXPLORATORY, NOT pre-registered. Clean room preserved.")
    print("Realistic sky (dN/df ~ f^-11/3, dN/dh ~ h^-4), generous array.")
    print("Proof of principle, not a forecast.\n")
    exp_C_and_f(out)
    exp_array(out)
    exp_band(out, rng)
    exp_mechanisms(out, rng)
    npz = os.path.join(REPORTS, "vernier4.npz")
    np.savez_compressed(npz, **out)
    print(f"[saved {npz} before plotting]")
    try:
        make_plots(out)
        print(f"wrote {os.path.join(REPORTS, 'VERNIER4_realistic.png')}")
    except Exception as e:
        print(f"WARNING: plotting failed ({type(e).__name__}: {e}); data is safe.")
    print("=" * 78)
    print(f"wrote {npz}")
    print("=" * 78)


if __name__ == "__main__":
    main()
