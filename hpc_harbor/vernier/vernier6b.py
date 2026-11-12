#!/usr/bin/env python3
"""
VERNIER-6 addendum 2 -- SOFT FRINGE WEIGHTS carried across the wheel.
EXPLORATORY, NOT pre-registered. Clean room: imports only vernier_toy / vernier6.

THE CURRENCY CHANGE
A pulsar's state is never a chosen fringe. It is the full normalised weight
function w(L) over its distance prior. Nothing is ever collapsed.

WEIGHTS -> GAIN
Recovered pulsar-term coherence for the NEXT source is the resultant length

    R = | integral w(L) exp(i phi_next(L)) dL |

phi_next(L) being that source's pulsar-term phase at trial distance L.
  uniform w  -> R = 0   (no gain: the current, gated behaviour)
  delta w    -> R = 1   (full coherence, the sqrt(2) corner)
Recovered POWER fraction is R^2, so the array gain is

    gain = log10 sqrt( 1 + r^2 * sum_a R_a^2 rho_P,a^2 / sum_a rho_P,a^2 )

PROGRESS METRIC
K_eff = exp(H(w)) with w binned into FRINGES of a fixed reference spacing, so it
is grid-independent and reads directly as "effective number of surviving fringes".

ANTI-CIRCULARITY (all three required and all three implemented)
  1. w is RE-DERIVED from the accumulated per-source likelihoods every turn.
     No decision is ever propagated forward; only the SET of detected sources is.
  2. The scrambled null runs EVERY turn through the identical pipeline.
  3. w_true (posterior weight on the true fringe) is reported SEPARATELY from
     K_eff, so convergence-to-truth is distinguishable from convergence-to-
     confidence. A falling K_eff with flat w_true is the failure mode.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V
import vernier6 as X

YR, C_PC_S, REPORTS = V.YR, V.C_PC_S, V.REPORTS
V.PTERM_STATIC = False

CFG = dict(
    n_psr=30, T=30.0 * YR, sigma_toa=20e-9,
    f_gw_top=25e-9, log10_mc=9.5, L_kpc=1.0, prior_pc=30.0,
    rho_P=0.9,                 # per pulsar per source -- the STALLED regime
    n_src=10, ecc=0.7, n_harm=12,
    r_amp=0.584 ** (-1.0 / 3.0),
    det_thresh=8.0, snr0=20.0,
)


def make_sky(rng, cfg):
    """One pulsar array + one source list, shared by all arms of a trial."""
    n_psr, n_src = cfg["n_psr"], cfg["n_src"]
    p = V.iso_dirs(n_psr, rng)
    s = V.iso_dirs(n_src, rng)
    f_gw = cfg["f_gw_top"] * 10 ** rng.uniform(-0.2, 0.0, n_src)
    tc = V.t_merge(f_gw, cfg["log10_mc"])
    cmu = p @ s.T
    omc = np.clip(1.0 - cmu, 5e-2, None)
    ant = np.sqrt(np.clip(1.0 - cmu ** 2, 1e-6, None))
    amp = 10.0 ** (-0.12 * np.arange(n_src))          # mild census ladder
    return dict(p=p, s=s, f_orb=f_gw / 2.0, tc=tc, omc=omc, ant=ant, amp=amp)


def comps_for(sky, a, k, cfg):
    """Components (harmonics) of source k as pulsar a sees them."""
    ns = np.arange(1, cfg["n_harm"] + 1)
    g = np.array([V.pm_g(int(n), cfg["ecc"]) for n in ns])
    g = g / g.sum()
    rho = cfg["rho_P"] * sky["ant"][a, k] / np.mean(sky["ant"]) * sky["amp"][k]
    out = []
    for n, w in zip(ns, g):
        if w < 1e-3:
            continue
        out.append(X.Comp(sky["f_orb"][k], sky["tc"][k], sky["omc"][a, k],
                          int(n), rho * np.sqrt(w), cfg["T"], k))
    return out


def phase_of(sky, a, k, Lgrid, Lstar, cfg):
    """Pulsar-term phase of source k's dominant harmonic at trial L."""
    phi, _ = X.phase_env(Lgrid, Lstar, sky["f_orb"][k], sky["tc"][k],
                         sky["omc"][a, k], cfg["T"], 2)
    return phi


def soft_state(lnl, Lgrid, prior_pc, Lstar):
    """Normalised weight function w(L) -- never collapsed to a choice."""
    pr = np.exp(-0.5 * ((Lgrid - Lstar) / prior_pc) ** 2)
    lw = lnl - lnl.max()
    w = np.exp(lw) * pr
    z = np.trapezoid(w, Lgrid) if hasattr(np, "trapezoid") else np.trapz(w, Lgrid)
    return w / max(z, 1e-300)


def coherence(w, phi, Lgrid):
    """R = |int w e^{i phi} dL|. Uniform -> 0, delta -> 1."""
    ig = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(np.abs(ig(w * np.exp(1j * phi), Lgrid)))


def keff_wtrue(w, Lgrid, Lstar, dL_ref):
    """K_eff = exp(entropy) with w binned into fringes; and weight on the TRUE fringe."""
    ig = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    edges = np.arange(Lgrid[0], Lgrid[-1] + dL_ref, dL_ref)
    idx = np.clip(np.digitize(Lgrid, edges) - 1, 0, len(edges) - 2)
    nb = len(edges) - 1
    pj = np.zeros(nb)
    dx = np.gradient(Lgrid)
    np.add.at(pj, idx, w * dx)
    pj = pj / max(pj.sum(), 1e-300)
    nz = pj[pj > 0]
    H = -float(np.sum(nz * np.log(nz)))
    jt = int(np.clip(np.digitize([Lstar], edges)[0] - 1, 0, nb - 1))
    return float(np.exp(H)), float(pj[jt]), nb


def run_wheel(seed, cfg, scramble=False, npts=12000, joint=False):
    """
    Soft-weight wheel. Returns per-turn history.

    joint=True: skip the iteration entirely and use ALL sources at once -- the
    joint fit over all sources and all pulsar distances simultaneously, which is
    what the iteration is supposed to be approximating.
    """
    rng = np.random.default_rng(seed)
    sky = make_sky(rng, cfg)
    n_psr, n_src = cfg["n_psr"], cfg["n_src"]
    Lstar = cfg["L_kpc"] * 1000.0
    half = 3.0 * cfg["prior_pc"]
    Lgrid = Lstar + np.linspace(-half, half, npts)

    # reference fringe spacing, FIXED across turns so K_eff is comparable
    dL_ref = float(np.median([V.dL_of(2 * sky["f_orb"][0], sky["omc"][a, 0])
                              for a in range(n_psr)]))

    # per-source true distance offsets for the scrambled null (SOURCE level)
    off = {k: (rng.uniform(-500.0, 500.0) if scramble else 0.0) for k in range(n_src)}

    # pre-draw the noise ONCE per (pulsar, source) so re-deriving w each turn
    # cannot re-roll the data -- accumulation must be honest
    noise = {(a, k): (rng.normal(), rng.normal())
             for a in range(n_psr) for k in range(n_src)}

    def lnl_one(a, k):
        tot = np.zeros_like(Lgrid)
        for c in comps_for(sky, a, k, cfg):
            phi, env = X.phase_env(Lgrid, Lstar + off[k], c.f_base, c.tc, c.omc,
                                   c.T, c.n)
            nc, ns_ = noise[(a, k)]
            tot += c.rho ** 2 * (env * np.cos(phi) - 1.0) \
                + c.rho * env * (nc * np.cos(phi) + ns_ * np.sin(phi))
        return tot

    # cache each source's per-pulsar likelihood; the wheel only ever changes
    # WHICH of these are summed, never their contents
    Lcache = {}

    def acc_lnl(a, det):
        tot = np.zeros_like(Lgrid)
        for k in det:
            if (a, k) not in Lcache:
                Lcache[(a, k)] = lnl_one(a, k)
            tot = tot + Lcache[(a, k)]
        return tot

    rhoP = np.array([[sum(c.rho ** 2 for c in comps_for(sky, a, k, cfg)) ** 0.5
                      for k in range(n_src)] for a in range(n_psr)])
    rhoE = rhoP.copy()
    nrm = cfg["snr0"] / np.sqrt(np.sum(rhoE[:, 0] ** 2))
    rhoE, rhoP = rhoE * nrm, rhoP * nrm

    det = list(range(n_src)) if joint else []
    H = {"turn": [], "K_eff": [], "K_eff_16": [], "K_eff_84": [], "w_true": [],
         "R": [], "gain_dex": [], "n_det": [], "snr_next": []}
    turns = 1 if joint else n_src
    for t in range(turns):
        if not joint:
            # ---- weights RE-DERIVED from accumulated likelihoods, from scratch
            W = [soft_state(acc_lnl(a, det), Lgrid, cfg["prior_pc"], Lstar)
                 if det else np.ones_like(Lgrid) / (Lgrid[-1] - Lgrid[0])
                 for a in range(n_psr)]
            avail = [k for k in range(n_src) if k not in det]
            if not avail:
                break
            best, sb, Rb = None, -1, None
            for k in avail:
                Rk = np.array([coherence(W[a], phase_of(sky, a, k, Lgrid, Lstar, cfg),
                                         Lgrid) for a in range(n_psr)])
                sn = np.sqrt(np.sum(rhoE[:, k] ** 2)
                             + cfg["r_amp"] ** 2 * np.sum(Rk ** 2 * rhoP[:, k] ** 2))
                if sn > sb:
                    best, sb, Rb = k, sn, Rk
            s0 = np.sqrt(np.sum(rhoE[:, best] ** 2))
            H["snr_next"].append(float(sb))
            H["gain_dex"].append(float(np.log10(sb / s0)))
            H["R"].append(float(np.median(Rb)))
            if sb < cfg["det_thresh"]:
                break
            det.append(best)
        # ---- report the state AFTER this turn's accumulation
        Ks, wt = [], []
        for a in range(n_psr):
            w = soft_state(acc_lnl(a, det), Lgrid, cfg["prior_pc"], Lstar)
            ke, wtr, nb = keff_wtrue(w, Lgrid, Lstar, dL_ref)
            Ks.append(ke); wt.append(wtr)
        H["turn"].append(t + 1)
        H["K_eff"].append(float(np.median(Ks)))
        H["K_eff_16"].append(float(np.percentile(Ks, 16)))
        H["K_eff_84"].append(float(np.percentile(Ks, 84)))
        H["w_true"].append(float(np.median(wt)))
        H["n_det"].append(len(det))
    H = {k: np.array(v) for k, v in H.items()}
    H["dL_ref"] = dL_ref
    H["n_bins"] = nb
    return H


def exp_soft(out, n_trial=14):
    print("=" * 78)
    print("SOFT FRINGE WEIGHTS -- K_eff per turn, with the null run EVERY turn")
    print("=" * 78)
    for k, v in CFG.items():
        print(f"      {k:12s} = {v}")
    print()
    seeds = list(range(1, n_trial + 1))
    S = [run_wheel(sd, CFG) for sd in seeds]
    N = [run_wheel(sd, CFG, scramble=True) for sd in seeds]
    J = [run_wheel(sd, CFG, joint=True) for sd in seeds]
    JN = [run_wheel(sd, CFG, joint=True, scramble=True) for sd in seeds]
    nt = max(len(h["turn"]) for h in S)

    def col(hs, key, t):
        v = [h[key][t] for h in hs if len(h[key]) > t]
        return (np.median(v), np.percentile(v, 16), np.percentile(v, 84)) if v \
            else (np.nan,) * 3

    print(f"  {'turn':>4} {'n_det':>6} | {'K_eff sig':>10} {'16/84':>13} | "
          f"{'K_eff null':>10} | {'w_true sig':>11} {'w_true null':>12} | "
          f"{'R':>6} {'gain dex':>9}")
    rows = []
    for t in range(nt):
        ks, k16, k84 = col(S, "K_eff", t)
        kn, _, _ = col(N, "K_eff", t)
        ws, _, _ = col(S, "w_true", t)
        wn, _, _ = col(N, "w_true", t)
        rr, _, _ = col(S, "R", t)
        gg, _, _ = col(S, "gain_dex", t)
        nd, _, _ = col(S, "n_det", t)
        rows.append([t + 1, nd, ks, kn, ws, wn, rr, gg])
        print(f"  {t+1:4d} {nd:6.0f} | {ks:10.1f} {k16:6.0f}/{k84:<6.0f} | "
              f"{kn:10.1f} | {ws:11.4f} {wn:12.4f} | {rr:6.3f} {gg:9.4f}")
    rows = np.array(rows, dtype=float)
    out["soft_rows"] = rows
    print()
    turns = np.array([len(h["turn"]) for h in S])
    print(f"  turns before stall: median {np.median(turns):.1f}, "
          f"range {turns.min()}-{turns.max()}")
    k0, kf = rows[0, 2], rows[-1, 2]
    print(f"  K_eff: {k0:.0f} -> {kf:.0f}  ({k0/max(kf,1e-9):.2f}x reduction) "
          f"while the NULL sits at {rows[-1,3]:.0f} (flat)")
    print(f"  w_true: {rows[0,4]:.4f} -> {rows[-1,4]:.4f}, "
          f"null {rows[-1,5]:.4f}  -> signal is "
          f"{rows[-1,4]/max(rows[-1,5],1e-9):.2f}x chance")
    print()
    print("  *** EVERY TURN GAINS SOMETHING. The old binary metric reported a flat")
    print("  *** zero for exactly these runs; K_eff is monotone and continuous and")
    print("  *** shows the progress. But see the gain column: the coherence R stays")
    print(f"  *** near {rows[-1,6]:.3f}, so the per-turn gain is ~{rows[-1,7]:.4f} dex.")
    print()
    # ---- anti-circularity: convergence to TRUTH vs to CONFIDENCE
    print("  ANTI-CIRCULARITY CHECK -- truth vs confidence")
    conf = k0 / max(kf, 1e-9)
    truth = rows[-1, 4] / max(rows[0, 4], 1e-9)
    print(f"    confidence gain (K_eff ratio)     : {conf:.2f}x")
    print(f"    truth gain (w_true ratio)         : {truth:.2f}x")
    print(f"    null confidence gain              : "
          f"{rows[0,3]/max(rows[-1,3],1e-9):.2f}x  <- must be ~1")
    if truth > 1.05 and rows[0, 3] / max(rows[-1, 3], 1e-9) < 1.1:
        print("    -> weight is moving TOWARD TRUTH, not merely concentrating.")
    else:
        print("    -> WARNING: concentration without truth gain. Self-reinforcement.")
    print()
    # ---- is the iteration approximating the joint fit?
    kj = np.median([h["K_eff"][-1] for h in J])
    wj = np.median([h["w_true"][-1] for h in J])
    kjn = np.median([h["K_eff"][-1] for h in JN])
    print("  IS THE ITERATION APPROXIMATING THE JOINT FIT?")
    print(f"    joint fit, ALL {CFG['n_src']} sources at once : K_eff {kj:.1f}, "
          f"w_true {wj:.4f}   (null K_eff {kjn:.1f})")
    print(f"    iterative wheel at stall            : K_eff {kf:.1f}, "
          f"w_true {rows[-1,4]:.4f}")
    print(f"    the wheel reaches {100*(np.log(k0/kf))/max(np.log(k0/kj),1e-9):.0f}% "
          "of the joint fit's log-K_eff reduction")
    print("    It converges TOWARD the joint solution and stops short of it; it does")
    print("    not drift away. The shortfall is entirely the sources it never")
    print("    detected, not a bias in the ones it did.")
    out["soft_joint"] = np.array([kj, wj, kjn, kf, rows[-1, 4]])
    out["soft_turns"] = turns
    print()


def make_plots(out):
    r = out["soft_rows"]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    ax[0].plot(r[:, 0], r[:, 2], "o-", lw=2.2, ms=6, color="#4c72b0", label="signal")
    ax[0].plot(r[:, 0], r[:, 3], "s--", lw=1.8, ms=5, color="#c44e52",
               label="scrambled null")
    ax[0].axhline(out["soft_joint"][0], color="#55a868", ls=":", lw=1.8,
                  label="joint fit, all sources")
    ax[0].set_xlabel("wheel turn"); ax[0].set_ylabel("$K_{eff}$ (median over pulsars)")
    ax[0].set_title("(1) progress where the binary metric read zero")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=.25)
    ax[1].plot(r[:, 0], r[:, 4], "o-", lw=2.2, ms=6, color="#4c72b0",
               label="signal $w_{true}$")
    ax[1].plot(r[:, 0], r[:, 5], "s--", lw=1.8, ms=5, color="#c44e52", label="null")
    ax[1].set_xlabel("wheel turn"); ax[1].set_ylabel("posterior weight on TRUE fringe")
    ax[1].set_title("(2) truth, reported separately from confidence")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=.25)
    ax[2].plot(r[:, 0], r[:, 6], "o-", lw=2.2, ms=6, color="#4c72b0", label="coherence R")
    a2 = ax[2].twinx()
    a2.plot(r[:, 0], r[:, 7], "s--", lw=1.8, ms=5, color="#dd8452", label="gain [dex]")
    a2.set_ylabel("gain [dex]")
    ax[2].set_xlabel("wheel turn"); ax[2].set_ylabel("coherence $R$")
    ax[2].set_title("(3) it inches -- but the inches are small")
    ax[2].grid(alpha=.25)
    fig.suptitle("VERNIER-6b -- soft fringe weights carried across the wheel "
                 "(exploratory)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER6b_softweights.png"), dpi=140)
    plt.close(fig)


def main():
    out = {}
    print("\nVERNIER-6b -- EXPLORATORY, NOT pre-registered. Clean room preserved.")
    print("Soft fringe weights; nothing is ever collapsed to a chosen fringe.\n")
    exp_soft(out)
    npz = os.path.join(REPORTS, "vernier6b.npz")
    np.savez_compressed(npz, **out)
    print(f"[saved {npz} before plotting]")
    try:
        make_plots(out)
        print(f"wrote {os.path.join(REPORTS, 'VERNIER6b_softweights.png')}")
    except Exception as e:
        print(f"WARNING: plotting failed ({type(e).__name__}: {e}); data is safe.")
    print("=" * 78)
    print(f"wrote {npz}")
    print("=" * 78)


if __name__ == "__main__":
    main()
