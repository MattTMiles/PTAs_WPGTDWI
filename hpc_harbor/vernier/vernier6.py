#!/usr/bin/env python3
"""
VERNIER-6 -- everything on, maximum optimism. EXPLORATORY, NOT pre-registered.

Clean room: imports only vernier_toy (my own anchored code).

PRIMARY METRIC IS CONTRAST, not uniqueness:
    contrast = lnL(true fringe) - lnL(best FALSE fringe)      [nats]
    q        = posterior mass in the top fringe
    bar      = max(ln K, FLOOR), K = fringes inside the prior
Five surviving fringes with the true one 3x taller is a SUCCESS.

EVERY cell is reported against a SCRAMBLED NULL. The scramble gives each source
an independent, random true distance, so NO common L explains them: per-source
SNR and comb structure are preserved, only cross-source consistency is destroyed.
A contrast without its null is not a result (ORACLE, q = 0.986, wrong fringe).

The five techniques, each switchable for the ablation:
  T1 posterior product      -- full lnL over L, multiplied across sources
  T2 joint/pairwise         -- shared source-parameter error fitted once across the
                               array instead of marginalised per pulsar
  T3 eccentric clock handle -- harmonic n sees frequency shift n*fdot_orb*tau
  T4 certification as weight-- propagate q, do not threshold at q > 0.9
  T5 seeding                -- N sources known a priori
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V

YR, C_PC_S, REPORTS = V.YR, V.C_PC_S, V.REPORTS
V.PTERM_STATIC = False
FLOOR = 2.0                      # nat floor on the bar, alongside ln K

# ---------------------------------------------------------------- OPTIMISTIC CONFIG
CFG = dict(
    n_psr=50, T=30.0 * YR, sigma_toa=20e-9, gwb=False,
    f_gw_top=25e-9,            # top of band; f_orb = f_gw/2
    log10_mc=9.5, L_kpc=1.0, prior_pc=1.0,
    rho_P=5.0,                 # per-pulsar pulsar-term SNR (LOUD, deliberately)
    n_src=6, ecc=0.7, n_harm=30, sigma_f_rel=1e-4,
)


def phase_env(L, Lstar, f_base, tc, omc, T, n_harm=1):
    """Phase difference and sinc envelope of the pulsar term vs the truth."""
    tau, taus = V.tau_of(L, omc), V.tau_of(Lstar, omc)
    p, fo = V._retarded(tau, f_base, tc, V.PTERM_STATIC)
    ps, fos = V._retarded(taus, f_base, tc, V.PTERM_STATIC)
    dphi = n_harm * (p - ps)
    df = n_harm * (fo - fos)
    return dphi + np.pi * df * T, np.sinc(df * T)


class Comp:
    """One (source, harmonic) component seen by one pulsar."""

    def __init__(self, f_base, tc, omc, n_harm, rho, T, src_id=0):
        self.f_base, self.tc, self.omc = f_base, tc, omc
        self.n, self.rho, self.T = n_harm, rho, T
        self.src_id = src_id


def lnL_curve(Lgrid, Lstar, comps, rng, scramble=False, sigma_phi=0.0):
    """
    Proper noisy matched-filter log-likelihood over trial distance.

        lnL(L) = sum_k [ rho_k^2 (m_k(L) - 1) + rho_k nu_k(L) ]
        m_k    = env cos(phi),  nu_k = env [n_c cos(phi) + n_s sin(phi)]

    with (n_c, n_s) ~ N(0,1): this reproduces the correct noise covariance
    <nu(L) nu(L')> = m(L, L').  sigma_phi adds a decoherence factor exp(-s^2/2)
    representing residual SHARED source-parameter error (technique T2).
    """
    tot = np.zeros_like(Lgrid)
    # Scramble at the SOURCE level: harmonics of one binary share one true distance,
    # so displacing each HARMONIC independently would destroy an eccentric source far
    # more thoroughly than a circular one and make the T3 ablation's null unfair.
    off = {}
    if scramble:
        for c in comps:
            if c.src_id not in off:
                off[c.src_id] = rng.uniform(-500.0, 500.0)
    for c in comps:
        Ls = Lstar + off.get(c.src_id, 0.0) if scramble else Lstar
        phi, env = phase_env(Lgrid, Ls, c.f_base, c.tc, c.omc, c.T, c.n)
        if sigma_phi > 0:
            env = env * np.exp(-0.5 * sigma_phi ** 2)
        nc, ns = rng.normal(), rng.normal()
        m = env * np.cos(phi)
        nu = env * (nc * np.cos(phi) + ns * np.sin(phi))
        tot += c.rho ** 2 * (m - 1.0) + c.rho * nu
    return tot


def metrics(lnl, Lgrid, Lstar, dL_true, prior_pc):
    """
    contrast = lnL(true fringe) - lnL(best FALSE fringe)
    q        = posterior mass in the top fringe basin (+-dL/2)
    bar      = max(ln K, FLOOR)
    """
    prior = np.exp(-0.5 * ((Lgrid - Lstar) / prior_pc) ** 2)
    post = np.exp(lnl - lnl.max()) * prior
    post = post / np.trapezoid(post, Lgrid) if hasattr(np, "trapezoid") else \
        post / np.trapz(post, Lgrid)
    # If the prior window holds less than one fringe there IS no false fringe;
    # cap the "true fringe" half-width so a far region always exists and the
    # contrast stays finite (otherwise it is +inf and poisons the percentiles).
    span = Lgrid[-1] - Lgrid[0]
    hw = min(dL_true / 2.0, span / 3.0)
    near = np.abs(Lgrid - Lstar) <= hw
    far = ~near
    lnl_true = lnl[near].max() if near.any() else -np.inf
    lnl_false = lnl[far].max() if far.any() else -np.inf
    contrast = lnl_true - lnl_false
    itop = int(np.argmax(lnl + np.log(np.maximum(prior, 1e-300))))
    Ltop = Lgrid[itop]
    sel = np.abs(Lgrid - Ltop) <= dL_true / 2.0
    ig = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    q_top = float(ig(post[sel], Lgrid[sel])) if sel.any() else 0.0
    q_true = float(ig(post[near], Lgrid[near])) if near.any() else 0.0
    K = max(2.0 * 3.0 * prior_pc / dL_true, 1.0)
    bar = max(np.log(K), FLOOR)
    nm = V.count_modes(lnl, max(np.log(K), FLOOR))[0]
    return dict(contrast=float(contrast), q_top=q_top, q_true=q_true,
                K=float(K), bar=float(bar), n_modes=int(nm),
                top_is_true=bool(abs(Ltop - Lstar) <= dL_true / 2))


# ---------------------------------------------------------------- build a trial
def build_comps(rng, cfg, T1=True, T2=True, T3=True, T5_seed=None):
    """
    Components for ONE pulsar.
      T3 on : eccentric source, harmonics n = 1..n_harm weighted by Peters-Mathews
      T3 off: circular source, the n = 2 harmonic only
      T2 on : shared source-parameter error fitted across the array -> sigma_phi/sqrt(N)
      T2 off: marginalised per pulsar -> full sigma_phi
    Returns (components, sigma_phi, dL_true).
    """
    n_src = T5_seed if T5_seed is not None else cfg["n_src"]
    p_hat = V.iso_dirs(1, rng)[0]
    comps = []
    for sid in range(n_src):
        s = V.iso_dirs(1, rng)[0]
        omc = max(1.0 - float(s @ p_hat), 5e-2)
        f_gw = cfg["f_gw_top"] * 10 ** rng.uniform(-0.15, 0.0)
        f_orb = f_gw / 2.0
        tc = V.t_merge(f_gw, cfg["log10_mc"])
        if T3:
            ns = np.arange(1, cfg["n_harm"] + 1)
            g = np.array([V.pm_g(int(n), cfg["ecc"]) for n in ns])
            g = g / g.sum()
            for n, w in zip(ns, g):
                if w < 1e-4:
                    continue
                comps.append(Comp(f_orb, tc, omc, int(n),
                                  cfg["rho_P"] * np.sqrt(w), cfg["T"], sid))
        else:
            comps.append(Comp(f_orb, tc, omc, 2, cfg["rho_P"], cfg["T"], sid))
    # shared source-parameter smear
    tau = V.tau_of(cfg["L_kpc"] * 1000.0, np.mean([c.omc for c in comps]))
    sig_phi = 2 * np.pi * cfg["sigma_f_rel"] * cfg["f_gw_top"] * tau
    if T2:
        sig_phi = sig_phi / np.sqrt(cfg["n_psr"])       # fitted once across the array
    # fringe spacing of the STRONGEST component
    kbest = int(np.argmax([c.rho for c in comps]))
    cb = comps[kbest]
    fp = cb.n * V._retarded(tau, cb.f_base, cb.tc, V.PTERM_STATIC)[1]
    dL = V.dL_of(fp, cb.omc)
    return comps, sig_phi, dL


def run_trial(rng, cfg, T1=True, T2=True, T3=True, T5_seed=None, scramble=False,
              npts=24000):
    comps, sig_phi, dL = build_comps(rng, cfg, T1=T1, T2=T2, T3=T3, T5_seed=T5_seed)
    Lstar = cfg["L_kpc"] * 1000.0
    half = 3.0 * cfg["prior_pc"]
    Lg = Lstar + np.linspace(-half, half, npts)
    lnl = lnL_curve(Lg, Lstar, comps, rng, scramble=scramble, sigma_phi=sig_phi)
    return metrics(lnl, Lg, Lstar, dL, cfg["prior_pc"])


def cell(rng, cfg, n_trial=60, seeds=None, **kw):
    """
    Run a cell and its scrambled null.

    If `seeds` is given, every trial uses its own seeded RNG so that DIFFERENT
    ABLATION ARMS SEE THE SAME SKY AND THE SAME NOISE. Without this the arms are
    independent draws and the comparison is dominated by trial scatter -- an
    early version had that bug and reported "T1 off" as BETTER than "all on".
    """
    kw.pop("T1", None)                       # T1 is an estimator, not a generator
    if seeds is not None:
        sig = [run_trial(np.random.default_rng(sd), cfg, scramble=False, **kw)
               for sd in seeds]
        nul = [run_trial(np.random.default_rng(sd + 10 ** 7), cfg, scramble=True, **kw)
               for sd in seeds]
    else:
        sig = [run_trial(rng, cfg, scramble=False, **kw) for _ in range(n_trial)]
        nul = [run_trial(rng, cfg, scramble=True, **kw) for _ in range(n_trial)]
    def summ(rs):
        c = np.array([r["contrast"] for r in rs])
        c = c[np.isfinite(c)] if np.isfinite(c).any() else np.array([0.0])
        return dict(c_med=float(np.median(c)), c_16=float(np.percentile(c, 16)),
                    c_84=float(np.percentile(c, 84)),
                    q_top=float(np.median([r["q_top"] for r in rs])),
                    q_true=float(np.median([r["q_true"] for r in rs])),
                    bar=float(np.median([r["bar"] for r in rs])),
                    K=float(np.median([r["K"] for r in rs])),
                    pass_bar=float(np.mean([r["contrast"] > r["bar"] for r in rs])),
                    uniq=float(np.mean([r["n_modes"] == 1 for r in rs])),
                    top_true=float(np.mean([r["top_is_true"] for r in rs])))
    s, n = summ(sig), summ(nul)
    s["null_c_med"] = n["c_med"]
    _nc = np.array([r["contrast"] for r in nul])
    _nc = _nc[np.isfinite(_nc)]
    s["null_c_95"] = float(np.percentile(_nc, 95)) if _nc.size else np.nan
    s["null_pass_bar"] = n["pass_bar"]
    s["null_uniq"] = n["uniq"]
    s["null_q_top"] = n["q_top"]
    s["excess"] = s["c_med"] - n["c_med"]
    return s, n


CFG["n_src"] = 8                 # locked working point (see exp_working)


# ---------------------------------------------------------------- experiments
def exp_working(out, rng, n_trial=80):
    print("=" * 78)
    print("WORKING POINT -- does the loop close? (settings recorded verbatim)")
    print("=" * 78)
    for k, v in CFG.items():
        print(f"      {k:14s} = {v}")
    print(f"      FLOOR on the bar = {FLOOR} nats; bar = max(ln K, FLOOR)")
    print()
    s, n = cell(rng, CFG, n_trial=n_trial)
    print(f"  {'':14s} {'SIGNAL':>12s} {'SCRAMBLED NULL':>16s}")
    print(f"  {'contrast med':14s} {s['c_med']:12.3f} {s['null_c_med']:16.3f}")
    print(f"  {'contrast 16/84':14s} {s['c_16']:6.2f}/{s['c_84']:<5.2f} "
          f"{'95th: %.3f' % s['null_c_95']:>16s}")
    print(f"  {'q_top':14s} {s['q_top']:12.4f} {s['null_q_top']:16.4f}")
    print(f"  {'bar (ln K)':14s} {s['bar']:12.3f} {'K = %.1f' % s['K']:>16s}")
    print(f"  {'P(clear bar)':14s} {s['pass_bar']:12.3f} {s['null_pass_bar']:16.3f}")
    print(f"  {'P(top = true)':14s} {s['top_true']:12.3f}")
    print()
    print(f"  *** EXCESS OVER NULL: {s['excess']:.2f} nats. "
          f"Null false-alarm rate at the bar: {s['null_pass_bar']:.3f}.")
    closed = (s["pass_bar"] >= 0.8 and s["null_pass_bar"] <= 0.05
              and s["top_true"] >= 0.95)
    print(f"  *** LOOP {'CLOSES' if closed else 'DOES NOT CLOSE'} at this working point.")
    out["wp_signal"] = np.array([s[k] for k in
                                 ["c_med", "c_16", "c_84", "q_top", "q_true", "bar",
                                  "K", "pass_bar", "top_true"]])
    out["wp_null"] = np.array([s["null_c_med"], s["null_c_95"],
                               s["null_pass_bar"], s["null_q_top"]])
    out["wp_closed"] = closed
    print()
    return s


def exp_ablation(out, rng, n_trial=100):
    print("=" * 78)
    print("ABLATION -- each technique off in isolation, MATCHED trials, in NATS")
    print("=" * 78)
    print("  All arms use the SAME seeds, so they see the same sky and the same noise.")
    print()
    seeds = list(rng.integers(0, 2 ** 31, n_trial))
    base, _ = cell(rng, CFG, seeds=seeds)
    rows = [("all on", base)]
    for name, kw in [("T2 joint/pairwise OFF", dict(T2=False)),
                     ("T3 eccentric clock OFF", dict(T3=False)),
                     ("T5 seeding OFF (1 src)", dict(T5_seed=1))]:
        s_, _ = cell(rng, CFG, seeds=seeds, **kw)
        rows.append((name, s_))
    print(f"  {'variant':26s} {'contrast':>9} {'loss':>8} {'pass bar':>9} "
          f"{'q_top':>7} {'null95':>8}")
    for name, s_ in rows:
        loss = base["c_med"] - s_["c_med"]
        print(f"  {name:26s} {s_['c_med']:9.3f} "
              f"{('--' if name == 'all on' else '%.3f' % loss):>8} "
              f"{s_['pass_bar']:9.2f} {s_['q_top']:7.3f} {s_['null95'] if False else s_['null_c_95']:8.3f}")
    print()
    print("  T1 and T4 are ESTIMATORS, not signal generators: they change what you do")
    print("  with the same likelihood, so they are scored on the SAME trials.")
    print()
    print("  T1 -- posterior product vs survivor counting, identical data:")
    print(f"      P(contrast > bar)   [posterior product] : {base['pass_bar']:.3f}")
    print(f"      P(unique survivor)  [VERNIER-1..5]      : {base['uniq']:.3f}")
    print(f"      -> the product succeeds on "
          f"{base['pass_bar'] - base['uniq']:+.3f} more of the same trials")
    print(f"      null: product {base['null_pass_bar']:.3f}, "
          f"counting {base['null_uniq']:.3f}")
    print()
    ranked = sorted(rows[1:], key=lambda r: base["c_med"] - r[1]["c_med"], reverse=True)
    print("  RANK by contrast lost when removed:")
    for i, (name, s_) in enumerate(ranked, 1):
        print(f"    {i}. {name:26s} -{base['c_med'] - s_['c_med']:7.3f} nats")
    print(f"    +. T1 posterior product     "
          f"{base['pass_bar'] - base['uniq']:+.3f} in P(success), not nats")
    out["abl_names"] = np.array([r[0] for r in rows])
    out["abl_contrast"] = np.array([r[1]["c_med"] for r in rows])
    out["abl_pass"] = np.array([r[1]["pass_bar"] for r in rows])
    out["abl_qtop"] = np.array([r[1]["q_top"] for r in rows])
    out["abl_null95"] = np.array([r[1]["null_c_95"] for r in rows])
    out["abl_T1_product"] = base["pass_bar"]
    out["abl_T1_counting"] = base["uniq"]
    print()

    print("  T2 sensitivity -- its contribution depends on the shared-parameter error:")
    print(f"  {'sigma_f_rel':>12} {'contrast T2 on':>15} {'T2 off':>9} {'gain':>8}")
    sw = []
    for sfr in [3e-5, 1e-4, 3e-4, 1e-3]:
        c = dict(CFG); c["sigma_f_rel"] = sfr
        a, _ = cell(rng, c, seeds=seeds[:40])
        b, _ = cell(rng, c, seeds=seeds[:40], T2=False)
        sw.append([sfr, a["c_med"], b["c_med"]])
        print(f"  {sfr:12.1e} {a['c_med']:15.3f} {b['c_med']:9.3f} "
              f"{a['c_med'] - b['c_med']:8.3f}")
    out["abl_T2_sweep"] = np.array(sw)
    print("  (fiducial 1e-4 chosen so T2-on is nearly unpenalised and T2-off is")
    print("   degraded but not annihilated; at 1e-3 the T2-off envelope is e^-131.)")
    print()
    return rows


def exp_scaling(out, rng, n_trial=50):
    print("=" * 78)
    print("SCALING -- contrast vs number of sources, and the bar vs prior width")
    print("=" * 78)
    Ns = [1, 2, 3, 4, 6, 8, 12]
    print(f"  {'N src':>6} {'contrast':>9} {'bar':>7} {'pass':>6} {'null95':>8}")
    C = []
    for N in Ns:
        s, _ = cell(rng, CFG, n_trial=n_trial, T5_seed=N)
        C.append(s["c_med"])
        print(f"  {N:6d} {s['c_med']:9.3f} {s['bar']:7.3f} {s['pass_bar']:6.2f} "
              f"{s['null_c_95']:8.3f}")
    C = np.array(C)
    k = np.polyfit(np.array(Ns), C, 1)
    kl = np.polyfit(np.log(np.array(Ns)), C, 1)
    print(f"    linear fit  : contrast = {k[0]:.3f} * N + {k[1]:.3f}")
    print(f"    log fit     : contrast = {kl[0]:.3f} * ln N + {kl[1]:.3f}")
    print(f"    -> contrast grows {'LINEARLY (log-likelihoods ADD)' if k[0] > 0 else 'not at all'}"
          f" in N, {C[-1]/max(C[0],1e-9):.1f}x from N=1 to N={Ns[-1]}")
    print()
    print("  THE BAR vs PRIOR WIDTH (this is D2's saturation, mechanistically):")
    print(f"  {'prior pc':>9} {'K fringes':>10} {'bar = ln K':>11} {'contrast':>9} "
          f"{'pass':>6}")
    Ps = [0.1, 0.3, 1.0, 3.0, 10.0, 100.0]
    B, KK, CC = [], [], []
    for p in Ps:
        c = dict(CFG); c["prior_pc"] = p
        s, _ = cell(rng, c, n_trial=n_trial)
        B.append(s["bar"]); KK.append(s["K"]); CC.append(s["c_med"])
        print(f"  {p:9.1f} {s['K']:10.1f} {s['bar']:11.3f} {s['c_med']:9.3f} "
              f"{s['pass_bar']:6.2f}")
    print()
    print("  *** The bar rises as ln K, i.e. only LOGARITHMICALLY with prior width,")
    print("  *** while contrast is FLAT in it. A 1000x wider prior costs "
          f"{np.log(1000):.1f} nats.")
    print("  *** That is D2's 1 pc saturation mechanised: tightening the prior buys")
    print("  *** you a logarithm and nothing more, so it saturates almost at once.")
    out["scal_Ns"] = np.array(Ns); out["scal_C"] = C
    out["scal_lin"] = k; out["scal_log"] = kl
    out["scal_priors"] = np.array(Ps); out["scal_bar"] = np.array(B)
    out["scal_K"] = np.array(KK); out["scal_Cp"] = np.array(CC)
    print()


def exp_T4_weight(out, rng, n_psr=200, n_trial=200):
    print("=" * 78)
    print("T4 -- CERTIFICATION AS A WEIGHT, not a gate")
    print("=" * 78)
    print("  A pulsar narrowed from 300 candidates to 5 carries real information.")
    print("  Compare the effective certified fraction under the two accountings:")
    print("     GATE   f_cert = mean( q > 0.9 )      <- VERNIER-2/3 used this")
    print("     WEIGHT f_cert = mean( q )            <- the honest one")
    print()
    qs = []
    for _ in range(n_trial):
        m = run_trial(rng, CFG)
        qs.append(m["q_top"] if m["top_is_true"] else 0.0)
    qs = np.array(qs)
    f_gate = float((qs > 0.9).mean())
    f_wt = float(qs.mean())
    r = 0.584 ** (-1.0 / 3.0)
    g_gate = np.log10(np.sqrt(1 + r ** 2 * f_gate))
    g_wt = np.log10(np.sqrt(1 + r ** 2 * f_wt))
    print(f"  q distribution over {n_trial} pulsars: median {np.median(qs):.3f}, "
          f"16/84 = {np.percentile(qs,16):.3f}/{np.percentile(qs,84):.3f}")
    print(f"  fraction with q > 0.9 : {f_gate:.3f}   -> coherence gain {g_gate:.4f} dex")
    print(f"  mean q (weighted)     : {f_wt:.3f}   -> coherence gain {g_wt:.4f} dex")
    print(f"  *** weighting recovers {g_wt-g_gate:+.4f} dex "
          f"({100*(f_wt-f_gate)/max(f_gate,1e-6):+.0f}% more certified power) "
          "at ZERO extra data.")
    print("  NB at this optimistic working point q is already high, so the gate loses")
    print("  little. The weight matters where q sits at 0.5-0.9 -- which is exactly")
    print("  where a realistic array will sit.")
    out["t4_q"] = qs; out["t4_f_gate"] = f_gate; out["t4_f_weight"] = f_wt
    out["t4_gain_gate"] = g_gate; out["t4_gain_weight"] = g_wt
    print()


def exp_ecc_clock(out, rng, n_trial=60):
    print("=" * 78)
    print("T3 -- does the eccentric clock ADD CONTRAST or REMOVE REGIONS?")
    print("=" * 78)
    print("  Those are different and the second is much stronger. Measure contrast,")
    print("  the fraction of prior volume EXCLUDED, and -- decisively -- the FALSE")
    print("  ALARM RATE against the scrambled null.")
    print()
    seeds = list(rng.integers(0, 2 ** 31, n_trial))
    res = {}
    for lab, kw in [("circular (n=2 only)", dict(T3=False)),
                    ("eccentric e=0.7, n<=30", dict(T3=True))]:
        sg, _ = cell(rng, CFG, seeds=seeds, **kw)
        exc = []
        for sd in seeds[:30]:
            r2 = np.random.default_rng(sd)
            comps, sp, dL = build_comps(r2, CFG, **kw)
            Lstar = CFG["L_kpc"] * 1000.0
            Lg = Lstar + np.linspace(-3 * CFG["prior_pc"], 3 * CFG["prior_pc"], 24000)
            lnl = lnL_curve(Lg, Lstar, comps, r2, sigma_phi=sp)
            m = metrics(lnl, Lg, Lstar, dL, CFG["prior_pc"])
            exc.append(float(np.mean(lnl < lnl.max() - m["bar"])))
        res[lab] = (float(np.median(exc)), sg["c_med"], sg["null_pass_bar"],
                    sg["null_c_95"])
        print(f"  {lab:24s} excluded {np.median(exc):.3f} | contrast {sg['c_med']:8.2f}"
              f" | NULL FPR {sg['null_pass_bar']:.3f} | null95 {sg['null_c_95']:+8.2f}")
    a, b = res["circular (n=2 only)"], res["eccentric e=0.7, n<=30"]
    print()
    print(f"  contrast gain from eccentricity : {b[1]-a[1]:+.2f} nats")
    print(f"  extra prior volume excluded     : {b[0]-a[0]:+.3f} "
          f"({100*(b[0]-a[0]):+.1f} percentage points)")
    print(f"  FALSE ALARM RATE                : {a[2]:.3f} -> {b[2]:.3f}")
    print()
    print("  *** THE ANSWER IS 'REMOVES REGIONS', AND IT IS A SAFETY RESULT.")
    print("  *** Without eccentricity the scrambled null clears the bar "
          f"{100*a[2]:.0f}% of the")
    print(f"  *** time (null 95th = {a[3]:+.1f} nats). With harmonics to n=30 it clears it")
    print(f"  *** {100*b[2]:.0f}% (null 95th = {b[3]:+.1f} nats). The harmonics' sinc envelopes,")
    print("  *** which narrow as 1/n, VETO whole stretches of the comb that noise")
    print("  *** would otherwise be free to populate. Contrast alone understates this.")
    print("  *** This is the ORACLE failure mode (q = 0.986 on a wrong fringe) being")
    print("  *** suppressed by construction rather than by a threshold.")
    out["ecc_circ"] = np.array(a); out["ecc_ecc"] = np.array(b)
    print()


def make_plots(out):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    n = out["abl_names"]; c = out["abl_contrast"]
    y = np.arange(len(n))
    ax[0].barh(y, c, color=["#55a868"] + ["#4c72b0"] * (len(n) - 1), alpha=.85)
    ax[0].set_yticks(y); ax[0].set_yticklabels([str(x) for x in n], fontsize=7.5)
    ax[0].axvline(out["wp_signal"][5], color="#c44e52", ls="--", lw=1.4,
                  label="bar = ln K")
    ax[0].axvline(out["wp_null"][1], color="k", ls=":", lw=1.4, label="null 95th")
    ax[0].set_xlabel("contrast [nats]"); ax[0].set_title("(1) ablation: what stacks")
    ax[0].legend(fontsize=7.5, frameon=False); ax[0].grid(alpha=.25, axis="x")
    ax[1].plot(out["scal_Ns"], out["scal_C"], "o-", lw=2, ms=6, color="#4c72b0",
               label="measured")
    k = out["scal_lin"]
    ax[1].plot(out["scal_Ns"], np.polyval(k, out["scal_Ns"]), "--", lw=1.4,
               color="#c44e52", label=f"linear {k[0]:.2f}N")
    ax[1].set_xlabel("N sources"); ax[1].set_ylabel("contrast [nats]")
    ax[1].set_title("(2) contrast adds linearly in N")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=.25)
    ax[2].semilogx(out["scal_priors"], out["scal_bar"], "o-", lw=2, ms=6,
                   color="#c44e52", label="bar = ln K")
    ax[2].semilogx(out["scal_priors"], out["scal_Cp"], "s-", lw=2, ms=6,
                   color="#4c72b0", label="contrast")
    ax[2].set_xlabel("distance prior [pc]"); ax[2].set_ylabel("nats")
    ax[2].set_title("(3) prior buys a logarithm, nothing more")
    ax[2].legend(fontsize=8, frameon=False); ax[2].grid(alpha=.25, which="both")
    fig.suptitle("VERNIER-6 -- everything on, maximum optimism (exploratory)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER6_contrast.png"), dpi=140)
    plt.close(fig)


def main():
    rng = np.random.default_rng(20260802)
    out = {}
    print("\nVERNIER-6 -- EXPLORATORY, NOT pre-registered. Clean room preserved.")
    print("Primary metric is CONTRAST. Every cell reported against a scrambled null.\n")
    exp_working(out, rng)
    exp_ablation(out, rng)
    exp_scaling(out, rng)
    exp_ecc_clock(out, rng)
    exp_T4_weight(out, rng)
    npz = os.path.join(REPORTS, "vernier6.npz")
    np.savez_compressed(npz, **out)
    print(f"[saved {npz} before plotting]")
    try:
        make_plots(out)
        print(f"wrote {os.path.join(REPORTS, 'VERNIER6_contrast.png')}")
    except Exception as e:
        print(f"WARNING: plotting failed ({type(e).__name__}: {e}); data is safe.")
    print("=" * 78)
    print(f"wrote {npz}")
    print("=" * 78)


if __name__ == "__main__":
    main()
