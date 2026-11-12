#!/usr/bin/env python3
"""
VERNIER-3 -- the seeded wheel. EXPLORATORY, NOT pre-registered.

Clean room preserved: imports only vernier_toy / vernier2 (my own anchored code).

Hand the wheel N sources already known and let it run. Two gearing mechanisms are
separated for the first time:

  M1  PULSAR-TERM COHERENCE from certifying distances.
      gain = sqrt(1 + r^2 f_cert), capped at ~0.19 dex as f_cert -> 1.
  M2  SOURCE SUBTRACTION. Removing a known source removes its contribution to the
      CONFUSION noise budget of the next one. Independent of pulsar distances and
      NOT capped by f_cert.

Noise budget (all as ratios to white):
    P(f)/P_wn = 1 + r_red (f/f_yr)^-3 [if pulsar has red noise]
                  + CONF0 * U(S) / U(empty)
where U(S) = sum of h^2 over population members NOT subtracted and not the target.
CONF0 is the confusion-to-white ratio with nothing subtracted. This is the same
caricature as VERNIER-2's GWB term, but now RESPONSIVE to subtraction -- which is
the whole point of M2.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V
import vernier2 as W

YR, C_PC_S, REPORTS = V.YR, V.C_PC_S, V.REPORTS
F_YR = 1.0 / YR
V.PTERM_STATIC = False

N_POP = 256
SIGMA_LOGH = 0.6
CONF0 = 3.0                  # confusion is 3x white before any subtraction
DET_THRESH = 8.0
SNR0 = 20.0                  # rank-0 array SNR with nothing subtracted


def population(rng, n_pop=N_POP, sigma=SIGMA_LOGH):
    """Ranked log10 h, fixed total band power (sum h^2 normalised to n_pop)."""
    lh = np.sort(rng.normal(0.0, sigma, n_pop))[::-1]
    h2 = 10.0 ** (2 * lh)
    lh = lh - 0.5 * np.log10(h2.sum() / n_pop)      # fix band power
    return lh


# ---------------------------------------------------------------- gap ladder
def exp_gap_ladder(out, rng, n_draw=20000, kmax=21):
    print("=" * 74)
    print("GAP LADDER -- measured brightness gap between adjacent census ranks")
    print("=" * 74)
    print(f"  {N_POP} members, log10 h ~ Normal(sigma={SIGMA_LOGH} dex), "
          f"{n_draw} realisations.")
    print("  gap(k) = log10 h_k - log10 h_(k+1).  The wheel turns at rank k when the")
    print("  gearing available exceeds gap(k).")
    print()
    L = np.sort(rng.normal(0.0, SIGMA_LOGH, (n_draw, N_POP)), axis=1)[:, ::-1]
    gaps = -np.diff(L[:, :kmax], axis=1)                 # (n_draw, kmax-1)
    cum = L[:, 0][:, None] - L[:, 1:kmax]                # dex below rank 1
    gmed = np.median(gaps, axis=0)
    cmed = np.median(cum, axis=0)
    # gearing available: M1 ceiling, and M1 at the f_cert the wheel actually reaches
    r_ext = 0.584 ** (-1.0 / 3.0)
    g_m1_max = np.log10(np.sqrt(1 + r_ext ** 2))
    g_m1_real = np.log10(np.sqrt(1 + r_ext ** 2 * 0.45))   # measured f_cert ~0.45
    print(f"  M1 ceiling (f_cert=1)   = {g_m1_max:.4f} dex")
    print(f"  M1 realistic (f_cert=.45)= {g_m1_real:.4f} dex")
    print()
    print(f"  {'rank k':>7} {'gap(k) dex':>11} {'cum dex below r1':>17} "
          f"{'SNR of rank k+1':>16} {'gap<M1max':>10} {'gap<M1real':>11}")
    for i in range(min(20, gaps.shape[1])):
        k = i + 1
        snr = SNR0 * 10 ** (-cmed[i])
        print(f"  {k:7d} {gmed[i]:11.4f} {cmed[i]:17.4f} {snr:16.2f} "
              f"{'YES' if gmed[i]<g_m1_max else 'no':>10} "
              f"{'YES' if gmed[i]<g_m1_real else 'no':>11}")
    xk = np.where(gmed < g_m1_real)[0]
    xk2 = np.where(gmed < g_m1_max)[0]
    print()
    print(f"  *** CROSSING: gap(k) drops below realistic M1 gearing at rank k = "
          f"{xk[0]+1 if xk.size else 'never'}")
    print(f"  *** CROSSING: gap(k) drops below the M1 CEILING at rank k = "
          f"{xk2[0]+1 if xk2.size else 'never'}")
    # absolute floor
    below = np.where(SNR0 * 10 ** (-cmed) < DET_THRESH)[0]
    print(f"  *** ABSOLUTE FLOOR: rank k+1 falls below SNR {DET_THRESH:.0f} at k = "
          f"{below[0]+1 if below.size else '>20'} "
          f"(cum {cmed[below[0]]:.3f} dex below rank 1)" if below.size else "")
    print()
    print("  READING: the two constraints bind at nearly the same rank. The relative")
    print("  gap stops mattering almost exactly where the absolute floor starts.")
    out["gl_gap_med"] = gmed
    out["gl_cum_med"] = cmed
    out["gl_gap_16"] = np.percentile(gaps, 16, axis=0)
    out["gl_gap_84"] = np.percentile(gaps, 84, axis=0)
    out["gl_m1_max"] = g_m1_max
    out["gl_m1_real"] = g_m1_real
    out["gl_cross_real"] = (xk[0] + 1) if xk.size else -1
    out["gl_cross_max"] = (xk2[0] + 1) if xk2.size else -1
    out["gl_floor_rank"] = (below[0] + 1) if below.size else -1
    print()
    return gmed, cmed


# ---------------------------------------------------------------- seeded wheel
def make_scenario(rng, n_psr=116, n_track=40, prior="vlbi"):
    """
    Build ONE sky + array. All arms are then run on the SAME scenario so that the
    M1/M2 comparison is matched -- otherwise each arm sees a different sky and the
    comparison is meaningless (an early version of this had that bug: 'both' came
    out BELOW 'M1' alone, which is impossible).
    """
    sc = {}
    lh_all = population(rng)
    sc["lh_all"] = lh_all
    sc["lh"] = lh_all[:n_track]
    sc["h2_all"] = 10.0 ** (2 * lh_all)
    sc["U_full"] = sc["h2_all"].sum()
    p, L, has_red = W.build_array(n_psr, rng)
    sc["p"], sc["L"], sc["has_red"] = p, L, has_red
    s_ = V.iso_dirs(n_track, rng)
    f0 = 10.0 ** rng.uniform(np.log10(W.BAND[0]), np.log10(W.BAND[1]), n_track)
    mc = rng.uniform(8.5, 9.5, n_track)
    sc["f0"], sc["tc"] = f0, V.t_merge(f0, mc)
    cmu = p @ s_.T
    sc["omc"] = np.clip(1.0 - cmu, 1e-3, None)
    sc["ant"] = np.sqrt(np.clip(1.0 - cmu ** 2, 1e-6, None))
    sc["f_p"] = f0[None, :] * (1.0 + V.tau_of(L[:, None], sc["omc"])
                               / sc["tc"][None, :]) ** (-3.0 / 8.0)
    sc["amp"] = 10.0 ** sc["lh"]
    _sv = {"subpc": 0.1, "vlbi": 1.0}.get(prior)
    sc["sig_arr"] = (np.full(n_psr, _sv) if _sv is not None
                     else {"pct1": 0.01, "pct10": 0.10, "pct20": 0.20}[prior] * L)
    sc["n_psr"], sc["n_track"] = n_psr, n_track
    sc["Kaj"] = np.maximum(6.0 * sc["sig_arr"][:, None]
                           / V.dL_of(sc["f_p"], sc["omc"]), 1.0)
    return sc


def _weights(sc, subtracted, arm):
    """Earth-term and pulsar-term SNR weights given the subtracted set."""
    Usub = sc["h2_all"][list(subtracted)].sum() if subtracted else 0.0
    conf = CONF0 * (sc["U_full"] - Usub) / sc["U_full"]
    if arm in ("M1", "none"):
        conf = CONF0                                  # subtraction disabled
    nr = 1.0 + conf * np.ones_like(sc["omc"])
    nr = nr + sc["has_red"][:, None] * W.R_RED * (sc["f0"][None, :] / F_YR) ** -3.0
    rE = sc["ant"] * sc["amp"][None, :] / np.sqrt(nr)
    rP = (sc["ant"] * sc["amp"][None, :]
          * (sc["f_p"] / sc["f0"][None, :]) ** (-1.0 / 3.0) / np.sqrt(nr))
    return rE, rP


def run_arm(sc, n_seed, arm="both", det_thresh=DET_THRESH, grid_cap=96000):
    """Run the wheel on a prebuilt scenario, seeded with the top n_seed sources."""
    n_psr, n_track = sc["n_psr"], sc["n_track"]
    rE0, _ = _weights(sc, set(), "none")
    nrm = SNR0 / np.sqrt(np.sum(rE0[:, 0] ** 2))

    def certify(detected):
        new = np.zeros(n_psr, bool)
        if arm in ("M2", "none"):
            return new                                # certification disabled
        _, rP = _weights(sc, set(detected), arm)
        rP = rP * nrm
        for a in range(n_psr):
            srcs = [V.Src(sc["f0"][j], sc["tc"][j], sc["omc"][a, j], 1, rP[a, j])
                    for j in detected]
            new[a] = (V.scan(sc["L"][a], sc["sig_arr"][a], srcs, cap=grid_cap)[0] == 1)
        return new

    # no-gearing reference: nothing subtracted, nothing certified, phase marginalised
    rE_ng, rP_ng = _weights(sc, set(), "none")
    rE_ng, rP_ng = rE_ng * nrm, rP_ng * nrm
    pm_ng = np.maximum(rP_ng ** 2 - 2.0 * np.log(sc["Kaj"]), 0.0)
    snr_ng_all = np.sqrt(np.sum(rE_ng ** 2, axis=0) + np.sum(pm_ng, axis=0))

    det = list(range(n_seed))
    cert = certify(det)
    H = {"n_cert": [], "f_cert": [], "gain_dex": [], "snr_next": [],
         "snr_next_nogear": [], "detected": [], "gap_to_next": [], "conf_ratio": []}
    stall = "ran out of tracked members"
    while True:
        avail = [j for j in range(n_track) if j not in det]
        if not avail:
            break
        rE, rP = _weights(sc, set(det), arm)
        rE, rP = rE * nrm, rP * nrm
        pm = np.maximum(rP ** 2 - 2.0 * np.log(sc["Kaj"]), 0.0)
        sw = np.array([np.sqrt(np.sum(rE[:, j] ** 2)
                               + np.sum(rP[cert, j] ** 2)
                               + np.sum(pm[~cert, j])) for j in avail])
        k = int(np.argmax(sw))
        best, snr_b = avail[k], float(sw[k])
        tp = np.sum(rP[:, best] ** 2)
        fc = float(np.sum(rP[cert, best] ** 2) / tp) if tp > 0 else 0.0
        Usub = sc["h2_all"][list(det)].sum() if det else 0.0
        H["n_cert"].append(int(cert.sum())); H["f_cert"].append(fc)
        H["snr_next"].append(snr_b)
        H["snr_next_nogear"].append(float(snr_ng_all[best]))
        H["gain_dex"].append(float(np.log10(snr_b / snr_ng_all[best])))
        H["gap_to_next"].append(float(sc["lh"][det[-1]] - sc["lh"][best])
                                if det else 0.0)
        # record the confusion actually USED: arms with subtraction disabled keep 1.0
        H["conf_ratio"].append(1.0 if arm in ("M1", "none")
                               else float((sc["U_full"] - Usub) / sc["U_full"]))
        if snr_b < det_thresh:
            stall = "next source below detection threshold"
            H["detected"].append(-1)
            break
        det.append(best); H["detected"].append(best)
        cert = certify(det)
    H = {kk: np.array(vv) for kk, vv in H.items()}
    H["turns"] = int((H["detected"] >= 0).sum())
    H["n_seed"], H["arm"], H["stall"] = n_seed, arm, stall
    return H


def exp_seeded(out, rng, seeds=(1, 2, 4, 8, 16), arms=("none", "M1", "M2", "both"),
               n_real=20, n_psr=116):
    print("=" * 74)
    print("SEEDED WHEEL -- hand it N known sources, separate M1 from M2")
    print("=" * 74)
    print(f"  {n_real} MATCHED realisations (all arms on the SAME sky), {n_psr} pulsars,")
    print(f"  1 pc priors, {N_POP}-member population, confusion = {CONF0:.0f}x white")
    print("  before subtraction. 'turns' = sources detected BEYOND the seed. Gain is")
    print("  measured against no subtraction + no certification, phase MARGINALISED.")
    print()
    acc = {(a, N): {"T": [], "FC": [], "G": [], "CR": [], "ST": [], "TOT": []}
           for a in arms for N in seeds}
    for _ in range(n_real):
        for N in seeds:
            sc = make_scenario(rng, n_psr=n_psr)      # one sky, all arms share it
            for a in arms:
                H = run_arm(sc, N, arm=a)
                d = acc[(a, N)]
                d["T"].append(H["turns"])
                d["TOT"].append(H["turns"] + N)          # TOTAL sources known
                d["FC"].append(H["f_cert"][0] if len(H["f_cert"]) else 0.0)
                d["G"].append(np.median(H["gain_dex"]) if len(H["gain_dex"]) else 0.0)
                d["CR"].append(H["conf_ratio"][-1] if len(H["conf_ratio"]) else 1.0)
                d["ST"].append(H["stall"])
    from collections import Counter
    res = {}
    print(f"  {'arm':>6} {'N seed':>7} {'turns':>7} {'TOTAL known':>12} "
          f"{'f_cert':>8} {'gain dex':>9} {'conf left':>10}")
    for a in arms:
        for N in seeds:
            d = acc[(a, N)]
            res[(a, N)] = (np.median(d["T"]), np.median(d["FC"]),
                           np.median(d["G"]), np.median(d["CR"]))
            print(f"  {a:>6} {N:7d} {np.median(d['T']):7.1f} "
                  f"{np.median(d['TOT']):12.1f} {np.median(d['FC']):8.3f} "
                  f"{np.median(d['G']):9.4f} {np.median(d['CR']):10.3f}")
        print()
    out["seed_Ns"] = np.array(seeds)
    out["seed_arms"] = np.array(arms)
    for key, idx in [("turns", 0), ("fcert", 1), ("gain", 2), ("conf", 3)]:
        out[f"seed_{key}"] = np.array([[res[(a, N)][idx] for N in seeds] for a in arms])
    out["seed_total"] = np.array([[np.median(acc[(a, N)]["TOT"]) for N in seeds]
                                  for a in arms])
    print("  --- SATURATION: TOTAL sources known (seed + turns)")
    for a in arms:
        tt = [np.median(acc[(a, N)]["TOT"]) for N in seeds]
        print(f"    {a:>6}: " + " ".join(f"N{n}->{v:.0f}" for n, v in zip(seeds, tt)))
    print("  If TOTAL is flat in N, seeding does NOT extend reach -- the wheel")
    print("  saturates at a depth set by the absolute floor, not by what it is given.")
    print()
    print("  --- M1 vs M2, matched")
    for N in seeds:
        g1, g2, gb = (res[("M1", N)][2], res[("M2", N)][2], res[("both", N)][2])
        t0, t1, t2, tb = (res[("none", N)][0], res[("M1", N)][0],
                          res[("M2", N)][0], res[("both", N)][0])
        print(f"    N={N:2d}: gain M1={g1:.4f} M2={g2:.4f} both={gb:.4f}  |  turns "
              f"none={t0:.0f} M1={t1:.0f} M2={t2:.0f} both={tb:.0f}  |  engine = "
              f"{'M2' if g2 > g1 else 'M1'}")
    print()
    fc = out["seed_fcert"][list(arms).index("both")]
    print("  --- f_cert vs N (both): " + " ".join(f"N{n}={v:.3f}"
                                                  for n, v in zip(seeds, fc)))
    j = np.diff(fc)
    if j.size:
        i = int(np.argmax(j))
        print(f"  *** CRITICAL N: biggest f_cert jump N={seeds[i]} -> N={seeds[i+1]} "
              f"({fc[i]:.3f} -> {fc[i+1]:.3f})")
    print()
    return res


# ---------------------------------------------------------------- absolute floor
def exp_floor(out, rng, seeds=(1, 2, 4, 8, 16), n_real=12, n_psr=116):
    print("=" * 74)
    print("ABSOLUTE FLOOR vs RELATIVE GAP -- which one binds at each N?")
    print("=" * 74)
    print("  At the stall we compare, for the source the wheel could not reach:")
    print("    gap_to_next  : dex below the last DETECTED source")
    print("    gain_achieved: dex the wheel actually delivered")
    print("    snr_next     : with gearing.  Stall means snr_next < 8.")
    print("  Classification: ABSOLUTE if even the M1 ceiling + full subtraction")
    print("  could not lift it over threshold; GAP if the shortfall is the gearing.")
    print()
    r_ext = 0.584 ** (-1.0 / 3.0)
    g_ceil = np.log10(np.sqrt(1 + r_ext ** 2))
    print(f"  {'N seed':>7} {'gap_next':>9} {'gain':>8} {'snr_next':>9} "
          f"{'snr_nogear':>11} {'max reach':>10} {'binds':>10}")
    rows = []
    for N in seeds:
        G, GA, S, SN, B = [], [], [], [], []
        for _ in range(n_real):
            sc = make_scenario(rng, n_psr=n_psr)
            H = run_arm(sc, N, arm="both")
            if not len(H["gap_to_next"]):
                continue
            gp = H["gap_to_next"][-1]; ga = H["gain_dex"][-1]
            sn = H["snr_next"][-1]; sg = H["snr_next_nogear"][-1]
            # could the ceiling have saved it?  conf already fully subtracted here.
            reach = sg * 10 ** g_ceil
            G.append(gp); GA.append(ga); S.append(sn); SN.append(sg)
            B.append("ABSOLUTE" if reach < DET_THRESH else "GAP")
        from collections import Counter
        b = Counter(B).most_common(1)[0][0] if B else "n/a"
        reach_med = np.median(SN) * 10 ** g_ceil if SN else np.nan
        rows.append((N, np.median(G), np.median(GA), np.median(S),
                     np.median(SN), reach_med, b))
        print(f"  {N:7d} {np.median(G):9.4f} {np.median(GA):8.4f} "
              f"{np.median(S):9.2f} {np.median(SN):11.2f} {reach_med:10.2f} {b:>10}")
    print()
    print("  READING: a small gap to an UNDETECTABLE source is still a stall. Where")
    print("  the classification reads ABSOLUTE, closing the gap would not help -- the")
    print("  next source is simply too faint for the array, and no amount of")
    print("  certification or subtraction reaches it.")
    out["floor_N"] = np.array([r[0] for r in rows])
    out["floor_gap"] = np.array([r[1] for r in rows])
    out["floor_gain"] = np.array([r[2] for r in rows])
    out["floor_snr"] = np.array([r[3] for r in rows])
    out["floor_snr_nogear"] = np.array([r[4] for r in rows])
    out["floor_reach"] = np.array([r[5] for r in rows])
    out["floor_binds"] = np.array([r[6] for r in rows])
    print()


def make_plots(out):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    # gap ladder vs gearing
    k = np.arange(1, len(out["gl_gap_med"]) + 1)
    ax[0].plot(k, out["gl_gap_med"], "o-", lw=2, ms=5, color="#4c72b0",
               label="measured gap(k)")
    ax[0].fill_between(k, out["gl_gap_16"], out["gl_gap_84"], color="#4c72b0",
                       alpha=.2)
    ax[0].axhline(out["gl_m1_max"], color="#c44e52", ls="--", lw=1.5,
                  label="M1 ceiling 0.193 dex")
    ax[0].axhline(out["gl_m1_real"], color="#dd8452", ls="-.", lw=1.5,
                  label="M1 at $f_{cert}$=0.45")
    ax[0].set_yscale("log"); ax[0].set_xlabel("census rank k")
    ax[0].set_ylabel("gap to rank k+1 [dex]")
    ax[0].set_title("the wheel turns where gearing > gap")
    ax[0].legend(fontsize=7.5, frameon=False); ax[0].grid(alpha=.25, which="both")
    # turns vs N per arm
    Ns = out["seed_Ns"]
    for i, a in enumerate(out["seed_arms"]):
        ax[1].plot(Ns, out["seed_turns"][i], "o-", lw=2, ms=6, label=str(a),
                   color=plt.cm.viridis(i / max(len(out["seed_arms"]) - 1, 1)))
    ax[1].set_xscale("log", base=2); ax[1].set_xticks(Ns)
    ax[1].set_xticklabels([str(n) for n in Ns])
    ax[1].set_xlabel("N seeded sources"); ax[1].set_ylabel("turns beyond the seed")
    ax[1].set_title("M1 vs M2 vs both"); ax[1].legend(fontsize=8, frameon=False)
    ax[1].grid(alpha=.25)
    # gain per arm
    for i, a in enumerate(out["seed_arms"]):
        ax[2].plot(Ns, out["seed_gain"][i], "s-", lw=2, ms=6, label=str(a),
                   color=plt.cm.viridis(i / max(len(out["seed_arms"]) - 1, 1)))
    ax[2].axhline(out["gl_m1_max"], color="#c44e52", ls="--", lw=1.3,
                  label="M1 ceiling")
    ax[2].set_xscale("log", base=2); ax[2].set_xticks(Ns)
    ax[2].set_xticklabels([str(n) for n in Ns])
    ax[2].set_xlabel("N seeded sources"); ax[2].set_ylabel("gain [dex]")
    ax[2].set_title("which mechanism is the engine?")
    ax[2].legend(fontsize=8, frameon=False); ax[2].grid(alpha=.25)
    fig.suptitle("VERNIER-3 -- the seeded wheel (exploratory, not pre-registered)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER3_seeded.png"), dpi=140)
    plt.close(fig)


def main():
    rng = np.random.default_rng(20260730)
    out = {}
    print("\nVERNIER-3 -- EXPLORATORY, NOT pre-registered. Clean room preserved.")
    print("Evolving pulsar term; M1 (coherence) and M2 (subtraction) separated.\n")
    exp_gap_ladder(out, rng)
    exp_seeded(out, rng)
    exp_floor(out, rng)
    npz = os.path.join(REPORTS, "vernier3.npz")
    np.savez_compressed(npz, **out)
    print(f"[saved {npz} before plotting]")
    try:
        make_plots(out)
        print(f"wrote {os.path.join(REPORTS, 'VERNIER3_seeded.png')}")
    except Exception as e:
        print(f"WARNING: plotting failed ({type(e).__name__}: {e}); data is safe.")
    print("=" * 74)
    print(f"wrote {npz}")
    print("=" * 74)


if __name__ == "__main__":
    main()
