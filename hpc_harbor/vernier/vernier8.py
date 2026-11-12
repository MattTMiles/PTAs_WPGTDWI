#!/usr/bin/env python3
"""
VERNIER-8 -- one dial. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier6 / vernier6b / vernier7.

Everything from VERNIER-7 is held fixed: the CORRECTED measured census ladder,
soft fringe weights, all three anti-circularity checks, per-turn source-level
scrambled null. The ONLY thing swept is the brightest source's per-pulsar SNR.

Note on what the dial does: rank-1 ARRAY SNR is renormalised to 20 in every cell,
so raising rho_P does NOT make sources easier to detect. It makes DISTANCES
easier to pin. That is the intended experiment -- how much does better
certification buy -- and it isolates the coherence channel cleanly.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V
import vernier6b as B
import vernier7 as S

YR, REPORTS = V.YR, V.REPORTS
V.PTERM_STATIC = False

RHOS = [0.9, 1.5, 2.5, 4.0, 5.0, 8.0]
LAD, GAPS = S.LAD, S.GAPS


def step_needed(rank):
    """dex the next source must be lifted to clear SNR 8 (the ABSOLUTE step)."""
    snr = S.CFG["snr0"] * 10 ** LAD[rank - 1]
    return max(np.log10(S.CFG["det_thresh"] / snr), 0.0)


def exp_sweep(out, n_trial=10):
    print("=" * 78)
    print("ONE DIAL -- sweep the brightest source's per-pulsar SNR")
    print("=" * 78)
    print("  Held fixed from VERNIER-7: measured ladder, soft weights, null every")
    print("  turn, all anti-circularity checks. Rank-1 array SNR renormalised to 20")
    print("  in every cell, so the dial changes CERTIFICATION, not detectability.")
    print()
    unaided = int(np.sum(S.CFG["snr0"] * 10 ** LAD >= S.CFG["det_thresh"]))
    print(f"  Unaided reach (no gearing at all): rank {unaided}")
    print()
    print(f"  {'rho_P':>6} {'final rank':>11} {'turns':>6} {'R median':>9} "
          f"{'max gain':>9} {'K_eff end':>10} {'w_true':>8} {'null w':>8} "
          f"{'sustains?':>10}")
    rows, curves = [], {}
    for rp in RHOS:
        cfg = dict(S.CFG); cfg["rho_P"] = rp
        Sr = [S.run_wheel(sd, cfg, seed_rank=1) for sd in range(1, n_trial + 1)]
        Nr = [S.run_wheel(sd, cfg, seed_rank=1, scramble=True)
              for sd in range(1, n_trial + 1)]
        fr = float(np.median([h["final_rank"] for h in Sr]))
        tn = float(np.median([h["turns"] for h in Sr]))
        Rm = float(np.median(np.concatenate([h["R"] for h in Sr if len(h["R"])])))
        gm = float(np.max([np.max(h["gain_dex"]) for h in Sr if len(h["gain_dex"])]))
        ke = float(np.median([h["K_eff"][-1] for h in Sr]))
        k0 = float(np.median([h["K_eff"][0] for h in Sr]))
        wt = float(np.median([h["w_true"][-1] for h in Sr]))
        w0 = float(np.median([h["w_true"][0] for h in Sr]))
        nw = float(np.median([h["w_true"][-1] for h in Nr]))
        nk = float(np.median([h["K_eff"][-1] for h in Nr]))
        nk0 = float(np.median([h["K_eff"][0] for h in Nr]))
        # gain-vs-rank curve, median across trials
        maxr = max(len(h["gain_dex"]) for h in Sr)
        g_by_rank, s_by_rank = [], []
        for i in range(maxr):
            g = [h["gain_dex"][i] for h in Sr if len(h["gain_dex"]) > i]
            r = [h["rank_next"][i] for h in Sr if len(h["rank_next"]) > i]
            if len(g) >= 3:
                g_by_rank.append(np.median(g))
                s_by_rank.append(step_needed(int(np.median(r))))
        curves[rp] = (np.array(g_by_rank), np.array(s_by_rank))
        sustains = bool(np.any(np.array(g_by_rank) > np.array(s_by_rank)))
        rows.append([rp, fr, tn, Rm, gm, ke, wt, nw, float(sustains),
                     k0, w0, nk, nk0])
        print(f"  {rp:6.1f} {fr:11.0f} {tn:6.0f} {Rm:9.3f} {gm:9.4f} {ke:10.1f} "
              f"{wt:8.4f} {nw:8.4f} {'YES' if sustains else 'no':>10}")
    rows = np.array(rows)
    out["sweep"] = rows
    for rp in RHOS:
        out[f"curve_g_{rp}"], out[f"curve_s_{rp}"] = curves[rp]
    print()
    print("  GAIN vs STEP, per rank (gain / step; > 1 means that rung is cleared):")
    print(f"  {'rho_P':>6} " + " ".join(f"{'r%d' % (i + 2):>7}" for i in range(6)))
    for rp in RHOS:
        g, st = curves[rp]
        cells = []
        for i in range(6):
            if i < len(g):
                cells.append(f"{g[i]/max(st[i],1e-4):7.2f}" if st[i] > 1e-4
                             else f"{'free':>7}")
            else:
                cells.append(f"{'--':>7}")
        print(f"  {rp:6.1f} " + " ".join(cells))
    print()
    sus = rows[rows[:, 8] > 0]
    if len(sus):
        print(f"  *** THE NUMBER: the wheel first sustains at rho_P = {sus[0,0]:.1f}")
    else:
        print(f"  *** THE NUMBER: the wheel NEVER sustains, up to rho_P = {RHOS[-1]}.")
        print("  *** Max gain reached anywhere in the sweep: "
              f"{rows[:,4].max():.4f} dex, against a step of "
              f"{step_needed(int(rows[:,1].max())):.4f} dex at the stall rank.")
    print()
    print("  ANTI-CIRCULARITY (only cells with >=2 turns are real tests):")
    print(f"  {'rho_P':>6} {'turns':>6} {'confidence x':>13} {'truth x':>9} "
          f"{'null conf x':>12} {'verdict':>26}")
    for r in rows:
        conf = r[9] / max(r[5], 1e-9)
        truth = r[6] / max(r[10], 1e-9)
        ncf = r[12] / max(r[11], 1e-9)
        if r[2] < 2:
            v = "(too few turns)"
        elif truth >= conf * 0.95 and ncf < 1.1:
            v = "truth outpaces confidence"
        else:
            v = "*** STOP: fails ***"
        print(f"  {r[0]:6.1f} {r[2]:6.0f} {conf:13.2f} {truth:9.2f} {ncf:12.2f} "
              f"{v:>26}")
    print()
    return rows, curves
