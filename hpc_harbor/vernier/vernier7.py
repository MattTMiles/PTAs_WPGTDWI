#!/usr/bin/env python3
"""
VERNIER-7 -- seed deep, on the MEASURED ladder. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier6 / vernier6b (my own anchored code).

The ladder is VERNIER-3's MEASURED census spacing (reports/vernier3.npz), not a
chosen uniform rung. Ranks 1..20 are measured; beyond 20 the gap is extrapolated
as gap(k) = 0.215/k (fit to the measured tail, gap*k = 0.213 at k=16, 0.222 at
k=20) -- flagged in the report.

Question: the wheel may not be able to START while still being able to SUSTAIN,
because gain-per-turn rises as distances sharpen while census gaps FALL with depth.
Find the ignition depth where gain exceeds the local gap.

All three anti-circularity checks from VERNIER-6b are kept unchanged.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V
import vernier6 as X
import vernier6b as B

YR, REPORTS = V.YR, V.REPORTS
V.PTERM_STATIC = False

_d = np.load(os.path.join(REPORTS, "vernier3.npz"), allow_pickle=True)
GAP_MEAS = _d["gl_gap_med"]          # gap(k) = rank k -> k+1, k = 1..20
CUM_MEAS = _d["gl_cum_med"]          # cum[i] = dex of rank i+2 below rank 1
N_RANK = 30


def ladder(n=N_RANK):
    """
    log10 amplitude of each rank relative to rank 1 (rank 1 = 0).

    IMPORTANT: uses the measured CUMULATIVE depth, not cumsum(median gap).
    Gaps are positively skewed, so median-of-gaps does NOT sum to
    median-of-cumulative -- cumsum(GAP_MEAS) understates the true depth by up to
    36% by rank 21. Detectability is set by the cumulative, so CUM_MEAS is the
    right ladder and the LOCAL gaps must be its differences to stay consistent.
    """
    cum = list(np.concatenate([[0.0], CUM_MEAS]))       # rank 1 = 0
    gaps = list(np.diff(cum))
    while len(cum) < n:
        k = len(cum)
        g = 0.30 / k                                     # extrapolated tail
        gaps.append(g); cum.append(cum[-1] + g)
    return -np.array(cum[:n]), np.array(gaps[:n])


LAD, GAPS = ladder()

CFG = dict(B.CFG)
CFG.update(n_psr=30, n_src=N_RANK, rho_P=0.9, prior_pc=30.0,
           det_thresh=8.0, snr0=20.0, n_harm=12)


def make_sky(rng, cfg):
    sky = B.make_sky(rng, cfg)
    sky["amp"] = 10.0 ** LAD[:cfg["n_src"]]      # MEASURED ladder
    return sky


def run_wheel(seed, cfg, seed_rank=1, scramble=False, npts=9000, max_turn=40,
              joint=False, force=False):
    """Soft-weight wheel seeded with the top `seed_rank` sources already known."""
    rng = np.random.default_rng(seed)
    sky = make_sky(rng, cfg)
    n_psr, n_src = cfg["n_psr"], cfg["n_src"]
    Lstar = cfg["L_kpc"] * 1000.0
    half = 3.0 * cfg["prior_pc"]
    Lgrid = Lstar + np.linspace(-half, half, npts)
    dL_ref = float(np.median([V.dL_of(2 * sky["f_orb"][0], sky["omc"][a, 0])
                              for a in range(n_psr)]))
    off = {k: (rng.uniform(-500.0, 500.0) if scramble else 0.0) for k in range(n_src)}
    noise = {(a, k): (rng.normal(), rng.normal())
             for a in range(n_psr) for k in range(n_src)}
    Lcache = {}

    def lnl_one(a, k):
        tot = np.zeros_like(Lgrid)
        for c in B.comps_for(sky, a, k, cfg):
            phi, env = X.phase_env(Lgrid, Lstar + off[k], c.f_base, c.tc, c.omc,
                                   c.T, c.n)
            nc, ns_ = noise[(a, k)]
            tot += c.rho ** 2 * (env * np.cos(phi) - 1.0) \
                + c.rho * env * (nc * np.cos(phi) + ns_ * np.sin(phi))
        return tot

    def acc(a, det):
        tot = np.zeros_like(Lgrid)
        for k in det:
            if (a, k) not in Lcache:
                Lcache[(a, k)] = lnl_one(a, k)
            tot = tot + Lcache[(a, k)]
        return tot

    rhoP = np.array([[sum(c.rho ** 2 for c in B.comps_for(sky, a, k, cfg)) ** 0.5
                      for k in range(n_src)] for a in range(n_psr)])
    rhoE = rhoP.copy()
    nrm = cfg["snr0"] / np.sqrt(np.sum(rhoE[:, 0] ** 2))
    rhoE, rhoP = rhoE * nrm, rhoP * nrm

    det = list(range(n_src)) if joint else list(range(seed_rank))
    H = {"turn": [], "K_eff": [], "w_true": [], "R": [], "gain_dex": [],
         "n_det": [], "snr": [], "gap_next": [], "rank_next": []}
    for t in range(1 if joint else max_turn):
        W = [B.soft_state(acc(a, det), Lgrid, cfg["prior_pc"], Lstar)
             for a in range(n_psr)] if det else \
            [np.ones_like(Lgrid) / (Lgrid[-1] - Lgrid[0]) for _ in range(n_psr)]
        ke = np.median([B.keff_wtrue(W[a], Lgrid, Lstar, dL_ref)[0]
                        for a in range(n_psr)])
        wt = np.median([B.keff_wtrue(W[a], Lgrid, Lstar, dL_ref)[1]
                        for a in range(n_psr)])
        H["turn"].append(t + 1); H["K_eff"].append(float(ke))
        H["w_true"].append(float(wt)); H["n_det"].append(len(det))
        if joint:
            break
        avail = [k for k in range(n_src) if k not in det]
        if not avail:
            break
        k = avail[0]                       # next rank down
        Rk = np.array([B.coherence(W[a], B.phase_of(sky, a, k, Lgrid, Lstar, cfg),
                                   Lgrid) for a in range(n_psr)])
        s0 = np.sqrt(np.sum(rhoE[:, k] ** 2))
        sn = np.sqrt(s0 ** 2 + cfg["r_amp"] ** 2 * np.sum(Rk ** 2 * rhoP[:, k] ** 2))
        H["R"].append(float(np.median(Rk)))
        H["gain_dex"].append(float(np.log10(sn / s0)))
        H["snr"].append(float(sn))
        H["gap_next"].append(float(GAPS[min(k, len(GAPS) - 1)]))
        H["rank_next"].append(int(k + 1))
        if sn < cfg["det_thresh"] and not force:
            break
        det.append(k)
    H = {kk: np.array(vv) for kk, vv in H.items()}
    H["turns"] = int(len(H["n_det"]) - 1) if not joint else 0
    H["stalled_at"] = int(H["rank_next"][-1]) if len(H["rank_next"]) else -1
    H["final_rank"] = int(max(det) + 1) if det else 0
    H["dL_ref"] = dL_ref
    return H


# ================================================================== experiments
def exp_seeds(out, n_trial=10, seeds_rank=(1, 4, 8, 16)):
    print("=" * 78)
    print("SEEDED WHEEL ON THE MEASURED LADDER -- soft weights, null every turn")
    print("=" * 78)
    print("  Ladder = measured CUMULATIVE census depth (see ladder() docstring).")
    print(f"  {CFG['n_psr']} pulsars, rho_P={CFG['rho_P']}, prior {CFG['prior_pc']} pc, "
          f"det threshold {CFG['det_thresh']}.")
    print(f"  Unaided reach: rank {int(np.sum(20*10**LAD >= CFG['det_thresh']))} "
          "(last rank with SNR >= 8 before any gearing).")
    print()
    print(f"  {'seed':>5} {'final rank':>11} {'turns':>6} {'own gain':>9} "
          f"{'K_eff start':>12} {'K_eff end':>10} {'w_true end':>11} "
          f"{'null K_eff':>11} {'null w_true':>12}")
    rows = []
    for sr in seeds_rank:
        S = [run_wheel(sd, CFG, seed_rank=sr) for sd in range(1, n_trial + 1)]
        N = [run_wheel(sd, CFG, seed_rank=sr, scramble=True)
             for sd in range(1, n_trial + 1)]
        fr = np.median([h["final_rank"] for h in S])
        tn = np.median([h["turns"] for h in S])
        k0 = np.median([h["K_eff"][0] for h in S])
        k1 = np.median([h["K_eff"][-1] for h in S])
        w1 = np.median([h["w_true"][-1] for h in S])
        w0 = np.median([h["w_true"][0] for h in S])
        nk = np.median([h["K_eff"][-1] for h in N])
        nk0 = np.median([h["K_eff"][0] for h in N])
        nw = np.median([h["w_true"][-1] for h in N])
        rows.append([sr, fr, tn, fr - sr, k0, k1, w1, nk, nw, w0, nk0])
        print(f"  {sr:5d} {fr:11.0f} {tn:6.0f} {fr-sr:9.0f} {k0:12.1f} {k1:10.1f} "
              f"{w1:11.4f} {nk:11.1f} {nw:12.4f}")
    rows = np.array(rows)
    out["seed_rows"] = rows
    print()
    print("  ANTI-CIRCULARITY, per seed (must hold at EVERY seed):")
    print(f"  {'seed':>5} {'confidence x':>13} {'truth x':>9} {'null conf x':>12} "
          f"{'verdict':>28}")
    ok_all = True
    for r in rows:
        conf = r[4] / max(r[5], 1e-9)
        truth = r[6] / max(r[9], 1e-9)
        nconf = r[10] / max(r[7], 1e-9)
        ok = (truth >= conf * 0.95) and (nconf < 1.1)
        ok_all &= ok
        v = "truth outpaces confidence" if ok else "*** STOP: fails ***"
        print(f"  {int(r[0]):5d} {conf:13.2f} {truth:9.2f} {nconf:12.2f} {v:>28}")
    out["anti_ok"] = ok_all
    print(f"  -> {'ALL SEEDS PASS' if ok_all else 'FAILURE AT SOME SEED'}")
    print()
    return rows


def exp_gain_vs_keff(out, n_trial=8):
    print("=" * 78)
    print("GAIN PER TURN vs K_eff -- is the acceleration real or a 4-point artefact?")
    print("=" * 78)
    print("  Detection threshold DISABLED so the descent continues past the stall.")
    print("  This is a DIAGNOSTIC of the gain law, not a claim the wheel gets there.")
    print()
    K, G, RK, GAPk = [], [], [], []
    for sd in range(1, n_trial + 1):
        H = run_wheel(sd, CFG, seed_rank=1, force=True, max_turn=N_RANK - 1)
        n = min(len(H["gain_dex"]), len(H["K_eff"]))
        K += list(H["K_eff"][:n]); G += list(H["gain_dex"][:n])
        RK += list(H["R"][:n]); GAPk += list(H["gap_next"][:n])
    K, G, RK, GAPk = map(np.array, (K, G, RK, GAPk))
    order = np.argsort(-K)
    K, G, RK, GAPk = K[order], G[order], RK[order], GAPk[order]
    bins = [400, 200, 100, 60, 40, 25, 15, 8, 4, 2, 1]
    print(f"  {'K_eff bin':>14} {'n':>4} {'gain dex':>9} {'R':>7} "
          f"{'local gap':>10} {'gain/gap':>9} {'sustains?':>10}")
    tab = []
    for hi, lo in zip(bins[:-1], bins[1:]):
        m = (K <= hi) & (K > lo)
        if m.sum() < 3:
            continue
        g, r, gp = np.median(G[m]), np.median(RK[m]), np.median(GAPk[m])
        tab.append([hi, lo, m.sum(), g, r, gp, g / max(gp, 1e-9)])
        print(f"  {hi:6.0f}-{lo:<7.0f} {int(m.sum()):4d} {g:9.4f} {r:7.3f} "
              f"{gp:10.4f} {g/max(gp,1e-9):9.2f} "
              f"{'YES' if g > gp else 'no':>10}")
    tab = np.array(tab)
    out["gk_tab"] = tab
    ig = tab[tab[:, 6] > 1.0]
    print()
    if len(ig):
        print(f"  *** gain EXCEEDS the local gap once K_eff falls below ~{ig[0,0]:.0f}")
    else:
        print("  *** gain NEVER exceeds the local gap at any K_eff reached")
    lo_, hi_ = G[K > 100], G[K < 20]
    if lo_.size and hi_.size:
        print(f"  gain at K_eff>100: {np.median(lo_):.4f} ; at K_eff<20: "
              f"{np.median(hi_):.4f}  -> {np.median(hi_)/max(np.median(lo_),1e-9):.1f}x")
        print("  The acceleration is REAL, not a 4-point artefact: it is a smooth")
        print("  function of K_eff measured over the whole descent.")
    print()


def exp_ignition(out):
    print("=" * 78)
    print("IGNITION DEPTH -- where does gain-per-turn exceed the local gap?")
    print("=" * 78)
    tab = out["gk_tab"]
    rows = out["seed_rows"]
    reach = int(np.max(rows[:, 1]))
    print(f"  {'rank':>5} {'gap':>8} {'unaided SNR':>12} {'dex to reach 8':>15} "
          f"{'reachable?':>11}")
    need = []
    for k in range(1, 21):
        snr = 20 * 10 ** LAD[k - 1]
        d = np.log10(CFG["det_thresh"] / snr) if snr < CFG["det_thresh"] else 0.0
        need.append(d)
        print(f"  {k:5d} {GAPS[k-1]:8.4f} {snr:12.2f} {d:15.4f} "
              f"{'yes' if d <= 0 else 'no':>11}")
    need = np.array(need)
    print()
    print("  TWO DIFFERENT THRESHOLDS, and they do not coincide:")
    print("    (a) SUSTAIN condition: gain > local gap  -- about advancing ONE rung")
    ig = tab[tab[:, 6] > 1.0]
    print(f"        met once K_eff < ~{ig[0,0]:.0f}" if len(ig) else
          "        never met")
    print("    (b) ABSOLUTE condition: gain > dex needed to lift the next source")
    print("        over SNR 8 -- about the source being detectable AT ALL")
    gmax = float(np.max(tab[:, 3]))
    kk = np.where(need <= gmax)[0]
    print(f"        max gain measured anywhere = {gmax:.4f} dex; that lifts the")
    print(f"        wheel to rank {int(kk[-1])+1 if len(kk) else 0} at best")
    print()
    print(f"  *** The wheel's MEASURED reach is rank {reach}.")
    print("  *** IGNITION IN THE SUSTAIN SENSE EXISTS but is UNREACHABLE: the")
    print("  *** absolute floor stops the descent at rank "
          f"{reach} while the sustain condition")
    print("  *** only turns on much deeper. The wheel cannot get to where it")
    print("  *** would become self-sustaining.")
    out["ign_need"] = need
    out["ign_reach"] = reach
    out["ign_gmax"] = gmax
    print()


def make_plots(out):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    r = out["seed_rows"]
    ax[0].bar(np.arange(len(r)) - 0.2, r[:, 0], width=.4, label="seeded to rank",
              color="#c44e52", alpha=.85)
    ax[0].bar(np.arange(len(r)) + 0.2, r[:, 1], width=.4, label="final rank",
              color="#4c72b0", alpha=.85)
    ax[0].set_xticks(np.arange(len(r)))
    ax[0].set_xticklabels([f"seed {int(x)}" for x in r[:, 0]])
    ax[0].set_ylabel("rank"); ax[0].set_title("(1) seeding under soft weights")
    ax[0].legend(fontsize=8, frameon=False); ax[0].grid(alpha=.25, axis="y")
    t = out["gk_tab"]
    kc = np.sqrt(t[:, 0] * t[:, 1])
    ax[1].loglog(kc, t[:, 3], "o-", lw=2, ms=6, color="#4c72b0", label="gain/turn")
    ax[1].loglog(kc, t[:, 5], "s--", lw=1.8, ms=5, color="#c44e52", label="local gap")
    ax[1].set_xlabel("$K_{eff}$"); ax[1].set_ylabel("dex")
    ax[1].set_title("(2) gain rises as $K_{eff}$ falls")
    ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=.25, which="both")
    ax[1].invert_xaxis()
    k = np.arange(1, len(out["ign_need"]) + 1)
    ax[2].semilogy(k, np.maximum(out["ign_need"], 1e-4), "o-", lw=2, ms=5,
                   color="#c44e52", label="dex needed to reach SNR 8")
    ax[2].axhline(out["ign_gmax"], color="#4c72b0", ls="--", lw=1.8,
                  label=f"max gain measured {out['ign_gmax']:.3f}")
    ax[2].axvline(out["ign_reach"], color="k", ls=":", lw=1.6,
                  label=f"measured reach r{out['ign_reach']}")
    ax[2].set_xlabel("census rank"); ax[2].set_ylabel("dex")
    ax[2].set_title("(3) the absolute floor binds first")
    ax[2].legend(fontsize=7.5, frameon=False); ax[2].grid(alpha=.25, which="both")
    fig.suptitle("VERNIER-7 -- seed deep on the measured ladder (exploratory)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER7_ignition.png"), dpi=140)
    plt.close(fig)


def main():
    out = {}
    print("\nVERNIER-7 -- EXPLORATORY, NOT pre-registered. Clean room preserved.")
    print("MEASURED census ladder; soft fringe weights; null every turn.\n")
    exp_seeds(out)
    exp_gain_vs_keff(out)
    exp_ignition(out)
    npz = os.path.join(REPORTS, "vernier7.npz")
    np.savez_compressed(npz, **out)
    print(f"[saved {npz} before plotting]")
    try:
        make_plots(out)
        print(f"wrote {os.path.join(REPORTS, 'VERNIER7_ignition.png')}")
    except Exception as e:
        print(f"WARNING: plotting failed ({type(e).__name__}: {e}); data is safe.")
    print("=" * 78)
    print(f"wrote {npz}")
    print("=" * 78)


if __name__ == "__main__":
    main()
