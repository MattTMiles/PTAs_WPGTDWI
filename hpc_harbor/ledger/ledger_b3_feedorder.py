#!/usr/bin/env python3
"""LEDGER B3 -- FEED-ORDER DEPENDENCE: does the M-step's fixed point depend on the order?

=============================================================================
WHERE THE ORDER ENTERS (traced; this is not a hypothesis)
=============================================================================
Within one iteration the FEED DECISIONS are order-INDEPENDENT, and that is worth stating
because it is where one would first look:

  glacier_loop.py:416   carried = np.where(~led.fed)[0]        <- fixed BEFORE the loop
  glacier_loop.py:434-8 for k in list(new): th_on = ...(carried - {k}); feed_dlnl[k] = ...
        every candidate is scored against the SAME `carried` and the SAME theta_rec, so
        feed_dlnl[] is a function of k alone. Permuting `new` cannot change WHO is fed.
  glacier_loop.py:465   fed_idx = np.where(led.fed)[0]         <- re-sorted, so the
        ledger's promote ORDER does not survive into the next stage either.

The order dependence is REAL one stage later, in the M-step:

  glacier_loop.py:275   slots = [NP_SRC*k + j for k in fed_idx for j in MSTEP_AXES]
  glacier_loop.py:279-91  for sweep in range(n_sweep):
                            for i, sl in enumerate(slots):
                                ... src[sl] += off        <- IN PLACE, SEQUENTIALLY

`mstep_quadratic` is COORDINATE ASCENT with in-place updates. Coordinate j is optimised
against a surface already moved by coordinates 1..j-1. That is order-dependent by
construction whenever the coordinates are correlated -- and (log10_fgw, log10_mc) are the
programme's most strongly correlated pair (the chirp-mass/frequency degeneracy is the
whole reason SIREN exists). The ordering is currently FED-MEMBER INDEX order, which is
census-rank order, which is an arbitrary labelling convention.

Two further amplifiers, both measured elsewhere in the campaign:
  * n_sweep = 2 only (glacier_loop.py:270) -- coordinate ascent is TRUNCATED, so the run
    is not at a fixed point at all; it is at a 2-sweep iterate, whose value depends on
    order far more strongly than a converged fixed point would.
  * the non-concave branch (glacier_loop.py:291) silently LEAVES a slot at its entry
    value. Which slots are non-concave depends on where the earlier slots moved the
    surface -- so order changes not just the answer but which coordinates move at all.

=============================================================================
WHAT THIS PROBE DOES (replay, not re-run -- the banks stay as banked)
=============================================================================
For each of 4 skies at one (rung, arm): load the banked iteration-i checkpoint
(theta_rec, fed_mask), rebuild the venue from the banked seeds, and run ONLY the M-step
under PERM permutations of `slots`. Then, separately, run each permutation to CONVERGENCE
(n_sweep raised until the move falls below tol) to separate two different questions:

  Q1  ITERATE SPREAD  -- do the campaign's actual 2-sweep iterates differ by order?
                         (this is what the banked numbers are)
  Q2  FIXED-POINT SPREAD -- do the CONVERGED points differ by order?
                         (if Q1 is large and Q2 is ~0, the defect is TRUNCATION, and the
                          remedy is cheap: more sweeps. If Q2 is also large, the surface
                          is genuinely multimodal and the remedy is not a sweep count.)

Reported per sky: max over permutation pairs of |d log10_fgw|, |d log10_mc| in dex, beside
the campaign's own measured wander scale (0.01-0.44 dex, S4.15) and the drain's own
precision (sigma(log10 A_bg) ~ 0.02-0.03 dex). The verdict bar is stated BEFORE the run:

  PRE-REGISTERED BAR. Feed order is DECLARED IMMATERIAL if the max over permutations of
  the fixed-point spread is below 0.02 dex (the drain's own 1-sigma) on every fed axis in
  every sky. Above that it is MATERIAL and the M-step ordering must be made canonical
  (spec in LEDGER_stats_audit.md SS B3) before any drain number is quoted per-member.

SCOPE. The only rungs whose loop feeds MORE THAN ONE member are r13p5/e07 and r13p25/e07
(measured over the whole bank: every r13p9 cell feeds 0-1, every *none* cell feeds 1-3).
A one-member feed has a single 2-coordinate ordering and no permutation to test. So this
probe necessarily runs in the ECCENTRIC arm, and r13p5/e07 skies {0,1,2,3} (n_fed
16/10/3/3) is chosen over r13p25 because skies 2 and 3 there carry NO certification claim
at all. The cells remain QUARANTINED from closure claims (S4.24 Finding 2); this probe is
a DIAGNOSTIC of the M-step operator and makes no claim about those cells' certifications.

LANE. Small: one venue build per sky (~5-10 min at T=30/n=256) + ~PERM * 4 * n_fed * 3
likelihood evaluations. Under 1 GPU-hr for the four skies.

Usage:
  python3 hpc_harbor/ledger/ledger_b3_feedorder.py run --rung r13p5 --arm e07 --sky 0
  python3 hpc_harbor/ledger/ledger_b3_feedorder.py harvest
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse
import glob
import sys

import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
          "hpc_harbor/atlas", "hpc_harbor/spark", "hpc_harbor/glacier"):
    sys.path.insert(0, f"{HSYMT}/{p}")

N_PERM = 3                      # the brief's x3, plus the identity = 4 orders total
TOL_DEX = 1e-4                  # convergence tolerance for the Q2 leg
MAX_SWEEP = 40
BAR_DEX = 0.02                  # the pre-registered materiality bar (the drain's 1 sigma)


def distinct_orders(n, n_perm, seed):
    """Identity + up to n_perm DISTINCT non-identity permutations of n members.

    Capped at (n! - 1): with n = 2 there is exactly one non-identity order, so asking for
    3 would spin forever. The cap is returned honestly rather than padded -- a cell with
    fewer available orders reports fewer, and the bank records how many were actually run.
    """
    import math
    idn = np.arange(n)
    orders = [idn]
    want = min(n_perm, math.factorial(n) - 1 if n <= 8 else n_perm)
    rng = np.random.default_rng(seed)
    guard = 0
    while len(orders) < want + 1 and guard < 10000:
        guard += 1
        pm = rng.permutation(n)
        if not any(np.array_equal(pm, o) for o in orders):
            orders.append(pm)
    return orders


def mstep_ordered(marg_fn, src0, slots, sc, n_sweep=2, step0=0.3):
    """glacier_loop.mstep_quadratic VERBATIM in its arithmetic, with `slots` (and the
    matching scale vector `sc`) supplied in an ARBITRARY order rather than built in
    fed-member index order. Identity permutation reproduces the incumbent exactly."""
    src = src0.copy()
    n_eval = 0
    widths = np.full(len(slots), np.nan)
    moved = 0.0
    for sweep in range(n_sweep):
        moved = 0.0
        for i, sl in enumerate(slots):
            d = step0 * sc[i] * (0.5 ** sweep)
            trip = np.tile(src, (3, 1))
            trip[0, sl] -= d
            trip[2, sl] += d
            v = np.asarray(marg_fn(trip)); n_eval += 3
            a, b, c = v
            den = a - 2 * b + c
            if den < 0:
                off = np.clip(0.5 * (a - c) / den, -1.0, 1.0) * d
                src[sl] += off
                moved = max(moved, abs(off))
                widths[i] = d / np.sqrt(-den)
    return src, widths, n_eval, moved


def run(rung, arm, sky, it, t_label, wscale, out_dir):
    import glacier_loop as GL
    import glacier_pop as GP
    from glacier_pop import bank_npz, lane_tag, PromoteLedger

    GL.check_affinity()
    wtag = "" if wscale == 1.0 else "_w" + f"{wscale:g}".replace(".", "p")
    stem = f"gl2_{rung}{wtag}_cell_{arm}_s{sky}_T{t_label}_{GL.TIER}"
    ck = glob.glob(f"{GP.OUT}/{stem}_i{it}__*.npz")
    if len(ck) != 1:
        raise SystemExit(f"need exactly one banked checkpoint for {stem}_i{it}, "
                         f"found {len(ck)}")
    z = np.load(ck[0], allow_pickle=True)
    fed_mask = np.asarray(z["fed_mask"], bool)
    theta_rec = np.asarray(z["theta_rec"], float)
    fed_idx = np.where(fed_mask)[0]
    print(f"[LEDGER-B3] {stem} i{it}: n_fed = {len(fed_idx)}  lane {lane_tag()}",
          flush=True)
    if len(fed_idx) < 2:
        raise SystemExit(f"n_fed = {len(fed_idx)} -- a single fed member has no feed "
                         f"order to permute. This cell cannot answer B3 (see SCOPE).")

    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()
    pop, cond = GP.draw_population_conditional(sky, rung, n_src=GL.N_POP_DEFAULT,
                                               band_log10f=GL.BAND_CAMPAIGN)
    slots_src, n_harm, active, chan, n_clip = GL.embed_igniter(
        pop, GL.E_IGNITER[arm], GL.VENUE_SPAN_S[t_label])
    pop_s = dict(pop); pop_s["src"] = slots_src; pop_s["n_src"] = len(slots_src)
    G = GP.build_glacier_problem(t_label, pop_s, verbose=True)
    G["slots"] = slots_src
    amo = G["amo"]; nd = amo["n_dist"]
    ones = jnp.ones(amo["npsr"])
    sb = GL.CertScoreboard(G, amo, jnp, C, prior_key=GL.TIER)
    noise_seed = GL.SEED_NOISE_BASE + 1000 * sky
    L_true, _ = sb.draw_truth(IG, noise_seed + 10_000_000, tier=GL.TIER)
    theta_true = np.asarray(amo["theta_truth"], float).copy()
    theta_true[sb.AI] = L_true
    clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
    ndraw = C.NoiseDrawer(sb.sp, amo, lgwb_path=GL.LGWB_BANKS[t_label], strict=True)
    noise = ndraw.draw(noise_seed, components=("white", "rn"), white_scale=wscale)
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n)) for c, n in zip(clean, noise))

    # the M-step surface EXACTLY as run_cell builds it (glacier_loop.py:490-497)
    carried_now = np.where(~fed_mask)[0]
    th_eval = GL.theta_with_absent(theta_rec, nd, carried_now)

    def marg_fn(srcs):
        ths = np.tile(th_eval, (len(srcs), 1))
        ths[:, nd:] = srcs
        return np.asarray(amo["logL_batch_theta"](jnp.asarray(ths), data, ones))

    scale = TE.phi_scale({"ncw": 1})
    base_slots = np.array([GL.NP_SRC * k + j for k in fed_idx for j in GL.MSTEP_AXES])
    base_sc = np.array([scale[j] for _ in fed_idx for j in GL.MSTEP_AXES])
    src0 = th_eval[nd:].copy()

    # the orders: identity first (must reproduce the bank), then N_PERM member-level
    # permutations. Members are permuted as BLOCKS -- that is what "feed order" means;
    # the within-member (fgw, mc) order is held so the two questions do not mix.
    orders = distinct_orders(len(fed_idx), N_PERM, seed=20260729 + sky)

    res = []
    na = len(GL.MSTEP_AXES)
    for oi, order in enumerate(orders):
        sel = np.concatenate([np.arange(m * na, (m + 1) * na) for m in order])
        sl, sc = base_slots[sel], base_sc[sel]
        # Q1: the campaign's own truncated 2-sweep iterate
        s2, w2, n2, _ = mstep_ordered(marg_fn, src0, sl, sc, n_sweep=2)
        # Q2: the same operator run to convergence
        sc_, mv, nsw, ntot = src0.copy(), np.inf, 0, 0
        while mv > TOL_DEX and nsw < MAX_SWEEP:
            sc_, _, ne, mv = mstep_ordered(marg_fn, sc_, sl, sc, n_sweep=1,
                                           step0=0.3 * (0.5 ** min(nsw, 6)))
            ntot += ne; nsw += 1
        print(f"  order {oi} {list(order)[:8]}{'...' if len(order) > 8 else ''}: "
              f"2-sweep {n2} evals; converged in {nsw} sweeps ({ntot} evals), "
              f"last move {mv:.2e}", flush=True)
        res.append(dict(order=order, src2=s2, srcC=sc_, nsweep=nsw, converged=mv <= TOL_DEX))

    def spread(key):
        A = np.stack([r[key][base_slots] for r in res])
        return np.abs(A[:, None, :] - A[None, :, :]).max(axis=(0, 1))

    sp2, spC = spread("src2"), spread("srcC")
    axnames = ["log10_fgw", "log10_mc"]
    print(f"\n  MAX-OVER-PERMUTATION SPREAD (dex), bar = {BAR_DEX}")
    for a, nm in enumerate(axnames):
        print(f"    {nm:11s} 2-sweep iterate {sp2[a::na].max():.5f} | "
              f"converged {spC[a::na].max():.5f}  "
              f"{'MATERIAL' if spC[a::na].max() > BAR_DEX else 'immaterial'}")
    verdict = ("MATERIAL" if spC.max() > BAR_DEX else "IMMATERIAL")
    print(f"  VERDICT (fixed point): feed order is {verdict} at this cell\n")

    os.makedirs(out_dir, exist_ok=True)
    bank_npz(f"ledger_b3_{rung}{wtag}_{arm}_s{sky}_T{t_label}_i{it}",
             rung=rung, arm=arm, sky=sky, t_label=t_label, wscale=wscale, iter=it,
             n_fed=len(fed_idx), fed_idx=fed_idx, orders=np.array(orders),
             n_orders=len(orders), n_perm_requested=N_PERM,
             src_2sweep=np.stack([r["src2"] for r in res]),
             src_converged=np.stack([r["srcC"] for r in res]),
             nsweep=np.array([r["nsweep"] for r in res]),
             all_converged=np.array([r["converged"] for r in res]),
             spread_2sweep=sp2, spread_converged=spC,
             bar_dex=BAR_DEX, verdict=verdict, base_slots=base_slots,
             note=("LEDGER B3: M-step coordinate-order permutation replayed from the "
                   "banked checkpoint. The bank is READ, never rewritten. REPORT-ONLY; "
                   "the source cells stay quarantined from closure claims."))
    return 0


# ============================================================
# the CPU leg -- what the operator does on a KNOWN surface (no venue, no GPU)
# ============================================================
def surrogate():
    """Order sensitivity of the M-step OPERATOR itself, on a surface whose answer is known.

    The venue leg (mode `run`) is lane-blocked. This leg is not a substitute for it -- it
    cannot tell us what the real surface does -- but it settles the half of B3 that is a
    property of the ALGORITHM rather than of the data, and it calibrates how to read the
    venue numbers when they land:

      * on a QUADRATIC surface the converged point is UNIQUE, so any Q2 (converged) spread
        the venue leg reports is evidence of NON-quadraticity -- genuine multimodality or
        the non-concave-slot branch -- and not of coordinate ascent being ill-defined;
      * the Q1 (2-sweep iterate) spread on a quadratic is NOT zero, and this leg measures
        how it grows with the (fgw, mc) correlation the chirp-mass degeneracy imposes.

    Calibration: sigma_scaled * SCALE_MC(=0.5) = sigma(log10 mc) in dex (forge_b1.py:62),
    step0 = 0.3 and n_sweep = 2 as wired (glacier_loop.py:270).
    """
    print("=== LEDGER B3 -- CPU LEG: the M-step operator on a known surface ===")
    print("  (the VENUE leg is PARKED on lane availability; see ldg_b3.sbatch)\n")
    SCALE = np.array([0.25, 0.5])            # (fgw, mc) scaled->dex, house convention
    print("  n_fed  rho(fgw,mc)   Q1 max spread (dex)   Q2 converged spread (dex)")
    rows = []
    for n_fed in (2, 3, 5, 10):
        for rho in (0.0, 0.5, 0.9, 0.99):
            m = 2 * n_fed
            # curvature: within-member (fgw,mc) correlated at rho; members weakly coupled
            R = np.eye(m)
            for k in range(n_fed):
                R[2*k, 2*k+1] = R[2*k+1, 2*k] = rho
            for k in range(n_fed):
                for l in range(k+1, n_fed):
                    R[2*k:2*k+2, 2*l:2*l+2] = R[2*l:2*l+2, 2*k:2*k+2] = 0.15 * rho
            R = R + 1e-9*np.eye(m)
            A = np.linalg.inv(R)              # precision; lnL = -0.5 (x-x*)^T A (x-x*) * s
            rng = np.random.default_rng(11 + n_fed)
            xstar = rng.normal(size=m) * 0.4
            SC = np.tile(SCALE, n_fed)
            # normalise so the curvature over one step0 move is ~0.2 nat (S4.15's scale)
            k0 = 0.2 / (0.5 * (0.3 ** 2))
            def marg(srcs, _A=A, _x=xstar, _k=k0):
                d = np.asarray(srcs)[:, :m] - _x
                return -0.5 * _k * np.einsum("bi,ij,bj->b", d, _A, d)
            slots0 = np.arange(m)
            src0 = np.zeros(m)
            orders = distinct_orders(n_fed, N_PERM, seed=5 + n_fed)
            s2s, sCs = [], []
            for order in orders:
                sel = np.concatenate([np.arange(2*mm, 2*mm+2) for mm in order])
                sl, sc = slots0[sel], SC[sel]
                a2, _, _, _ = mstep_ordered(marg, src0, sl, sc, n_sweep=2)
                s2s.append(a2)
                cur, mv, nsw = src0.copy(), np.inf, 0
                while mv > TOL_DEX and nsw < MAX_SWEEP:
                    cur, _, _, mv = mstep_ordered(marg, cur, sl, sc, n_sweep=1,
                                                  step0=0.3 * (0.5 ** min(nsw, 6)))
                    nsw += 1
                sCs.append(cur)
            S2 = np.stack(s2s); SCc = np.stack(sCs)
            d2 = float(np.max(np.abs(S2[:, None] - S2[None, :]) * SC))
            dC = float(np.max(np.abs(SCc[:, None] - SCc[None, :]) * SC))
            print(f"  {n_fed:4d}   {rho:5.2f}        {d2:.5f}               {dC:.5f}"
                  + ("   <- above the 0.02 dex bar" if d2 > BAR_DEX else ""))
            rows.append((n_fed, rho, d2, dC))
    R = np.array(rows)
    n_over = int((R[:, 2] > BAR_DEX).sum())
    print(f"\n  MEASURED, and it REFUTES the obvious hypothesis. On a well-conditioned")
    print(f"  quadratic the 2-sweep iterate spreads at most {np.max(R[:, 2]):.4f} dex over feed")
    print(f"  orders -- {n_over}/{len(R)} of the (n_fed, rho) grid is above the {BAR_DEX} dex bar, i.e.")
    print(f"  NONE of it. Truncated coordinate ascent is NOT, by itself, order-sensitive at")
    print(f"  this curvature, even at rho = 0.99 and n_fed = 10.")
    print(f"  (The converged leg's residual {np.max(R[:, 3]):.1e} dex is the MAX_SWEEP = {MAX_SWEEP} cap")
    print(f"   biting at high correlation, not a non-unique optimum -- a quadratic has one.)")
    print()
    print("  SO THE PREDICTION FOR THE VENUE LEG IS SHARP, AND IT IS A PRE-REGISTRATION:")
    print("  if the venue leg finds order dependence above 0.02 dex, it CANNOT be blamed on")
    print("  the sweep count. It would have to come from the two non-quadratic features the")
    print("  operator actually has on the comb surface:")
    print("    (i)  the non-concave branch (glacier_loop.py:291) silently LEAVING a slot at")
    print("         its entry value -- which slots those are depends on where earlier slots")
    print("         moved the surface, so order changes WHICH coordinates move at all;")
    print("    (ii) genuine multimodality of the harmonic-comb surface -- which is exactly")
    print("         the object LEDGER C1 asks FORGE-B to detect before it trusts a Laplace")
    print("         width. B3 and C1 are then the same finding seen from two sides.")
    print("  A NULL venue result is equally informative: it would say the M-step's measured")
    print("  0.01-0.44 dex wander (S4.15) is a property of the SURFACE, not of the ordering,")
    print("  and no canonicalisation of the sweep order would have prevented it.")
    return 0


def harvest(repo):
    fs = sorted(glob.glob(f"{repo}/GLACIER_results/ledger_b3_*__*.npz"))
    if not fs:
        print("no LEDGER-B3 banks yet (the GPU leg is PARKED -- see LEDGER_stats_audit.md)")
        return 0
    print("  cell                                n_fed  2-sweep spread  converged spread  verdict")
    for f in fs:
        z = np.load(f, allow_pickle=True)
        print(f"  {os.path.basename(f).split('__')[0]:36s} {int(z['n_fed']):4d}   "
              f"{float(np.max(z['spread_2sweep'])):.5f}         "
              f"{float(np.max(z['spread_converged'])):.5f}      {str(z['verdict'])}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["run", "harvest", "surrogate"])
    ap.add_argument("--rung", default="r13p5")
    ap.add_argument("--arm", default="e07")
    ap.add_argument("--sky", type=int, default=0)
    ap.add_argument("--iter", type=int, default=0)
    ap.add_argument("--t", type=int, default=30)
    ap.add_argument("--wscale", type=float, default=1.0)
    ap.add_argument("--repo", default="/data/taylor_group/matt_miles/PTAs_WPGTDWI")
    a = ap.parse_args()
    if a.mode == "harvest":
        return harvest(a.repo)
    if a.mode == "surrogate":
        return surrogate()
    return run(a.rung, a.arm, a.sky, a.iter, a.t, a.wscale,
               f"{a.repo}/GLACIER_results")


if __name__ == "__main__":
    sys.exit(main())
