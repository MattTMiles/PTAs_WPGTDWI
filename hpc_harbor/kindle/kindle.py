"""
KINDLE driver — g0 license + off-truth pilot. Lane: A100-80 (HPC_SETUP §2).

modes:
  g0     reproduce a banked IGNITE-2 soft-loop trajectory BIT-IDENTICALLY through the
         KINDLE loop body (start_mode='truth'). The license. Chosen target: the ignition
         run ig_sloop_h1275_T30_lit_g0_n6240011 (n_cert 0->1) — the very event the ratio
         gain NaN'd — so g0 doubles as the demonstration that additive dN reads it as +1.
  pilot  the off-truth instrumentation, at a banked T=30 cell (no new problem build needed):
         arm (a) truth control + arm (b-FIX) 1 x sigma(mc) offset + scrambled through both,
         same seeds. Proves dN/margin-flow move and scrambled stays silent BEFORE the
         launched ladder (Pair B @ T=40, the loudness/structure/eccentric variants) is
         committed to GPU.
"""
import os, sys, argparse, time
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
sys.path.insert(0, f"{HSYMT}/hpc_harbor/ignite2")
sys.path.insert(0, f"{HSYMT}/hpc_harbor/kindle")
KOUT = f"{HSYMT}/KINDLE_results"

import ignite2 as IG2
from kindle_loop import run_kindle_realisation, force_recompute_lgwb


def _cols(d, keys):
    return {k: np.asarray(d[k]) for k in keys}


def mode_g0(args):
    jax, jnp, IG, C, B1, TE, F = IG2._import_stack()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    force_recompute_lgwb(C)   # T=30 IGNITE geometry: recompute L_gwb at cpus=8 (§10.16(e))
    os.makedirs(KOUT, exist_ok=True)
    floors = IG2.load_floors()
    P = IG.build_ignite_problem(30)
    skies = {0: F.load_geo_skies([0])[0]}
    E_cache = {}

    # target banked run: the ignition trajectory (n_cert 0->1) at the lit onset cell
    h, T, tier, gid, ns = -12.75, 30, "lit", 0, 6240011
    fl = floors[(h, T, tier)]["fN"]
    print(f"\ng0: reproducing ig_sloop_h1275_T30_lit_g0_n{ns} through the KINDLE loop "
          f"(start_mode=truth) — floor {fl:.4f}", flush=True)
    tag, msg = run_kindle_realisation(jnp, IG, C, B1, TE, F, P, E_cache, h, tier, fl,
                                      gid, skies[gid], ns, ns + 10_000_000,
                                      scrambled=False, start_mode="truth", start_scale=0.0,
                                      prefix="g0")
    print(f"  {tag}: {msg}", flush=True)

    kp = f"{KOUT}/k_{tag}.npz"
    bp = f"{HSYMT}/reports/ig_sloop_h1275_T30_lit_g0_n{ns}.npz"
    if not os.path.exists(bp):
        bp = f"{HSYMT}/IGNITE2_results/ig_sloop_h1275_T30_lit_g0_n{ns}.npz"
    K = np.load(kp, allow_pickle=True); Bk = np.load(bp, allow_pickle=True)

    print("\n=== g0 GATE — per-iteration raw columns, KINDLE vs banked ===", flush=True)
    shared = ["iter_dlnL", "iter_qmax", "iter_mapk", "dlnL_final", "qmax_final",
              "mapk_final", "cert_final", "on_true_final", "lnK", "n_true_grid",
              "traj_n_cert", "traj_wrong", "traj_W", "traj_sig_mc_dex", "traj_lnmarg",
              "traj_step_norm", "traj_src_mc_off", "x_final"]
    worst = 0.0; nfail = 0
    for k in shared:
        a = np.nan_to_num(np.asarray(K[k], float), nan=-1.0, posinf=1e30)
        b = np.nan_to_num(np.asarray(Bk[k], float), nan=-1.0, posinf=1e30)
        if a.shape != b.shape:
            print(f"  {k:20s} SHAPE {a.shape} vs {b.shape}  *** MISMATCH")
            nfail += 1; continue
        dd = float(np.max(np.abs(a - b))) if a.size else 0.0
        worst = max(worst, dd)
        flag = "" if dd == 0.0 else "  *** NONZERO"
        nfail += (dd != 0.0)
        print(f"  {k:20s} max|diff| = {dd:.3e}{flag}")
    print(f"\n  worst = {worst:.3e} over {len(shared)} columns; {nfail} nonzero", flush=True)

    # the demonstration: additive dN reads the ignition that ratio-gain NaN'd
    nc = np.asarray(K["traj_n_cert"], int)
    g = np.asarray(K["traj_gain"], float)
    dN = np.asarray(K["traj_dN"], float)
    print(f"\n=== the ignition, read three ways ===")
    print(f"  n_cert traj : {nc.tolist()}")
    print(f"  RATIO gain  : {np.round(g,3).tolist()}   <- NaN across the 0->1 ignition")
    print(f"  ADDITIVE dN : {dN.astype(int).tolist()}   <- reads the ignition as +1")
    print(f"  C0(start)={int(K['C0_start'])} C_fix={int(K['C_fix'])} dN_final={int(K['dN_final'])}")

    ok = (worst == 0.0 and nfail == 0)
    print(f"\ng0 {'PASSED — the KINDLE loop body IS the banked loop body (0.0e+00).' if ok else 'FAILED'}",
          flush=True)
    return 0 if ok else 1


def mode_pilot(args):
    jax, jnp, IG, C, B1, TE, F = IG2._import_stack()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    force_recompute_lgwb(C)   # T=30 IGNITE geometry: recompute L_gwb at cpus=8 (§10.16(e))
    os.makedirs(KOUT, exist_ok=True)
    floors = IG2.load_floors()
    P = IG.build_ignite_problem(30)
    # -1 = fiducial POP sky (geo_src=None); 0-3 = GEO draws — mirror mode_softloop
    skies = {g: (None if g < 0 else F.load_geo_skies([g])[0]) for g in set(IG.SKIES)}
    E_cache = {}

    # pilot cell: the VLBI onset cell (tight priors, some truth-start certs to repair)
    h, T, tier = -13.25, 30, "vlbi"
    fl = floors[(h, T, tier)]["fN"]
    reps = list(range(args.nreps))
    print(f"\npilot at ({h},{T},{tier}) floor {fl:.4f} — {len(reps)} reps x "
          f"[truth, fix_mc @ {args.scale}sig] + scrambled\n", flush=True)

    # match the banked soft_plan seed convention so the truth arm reproduces a banked
    # trajectory (free extra validation): ns = 6_000_000 + 100_000*T_idx + 10_000*ci + rep
    ci = IG.cells().index((h, tier))
    base = 6_000_000 + 100_000 * IG.T_GRID.index(T) + 10_000 * ci
    base_scr = 6_500_000 + 100_000 * IG.T_GRID.index(T) + 10_000 * ci
    plan = []
    for rep in reps:
        gid = IG.SKIES[rep % len(IG.SKIES)]
        ns = base + rep               # SAME seeds as the banked truth-start sloop reals
        plan.append(dict(gid=gid, ns=ns, scr=False, sm="truth", sc=0.0))
        plan.append(dict(gid=gid, ns=ns, scr=False, sm="fix_mc", sc=args.scale))
    # scrambled arm = SINGLE perturbation (scrambled source at its own position, truth-start).
    # scrambled + fix_mc is a double perturbation and is numerically pathological (all-NaN
    # estep) — a scrambled isotropic source displaced further in mc leaves the covariance
    # degenerate. The meaningful null (off-truth start + no signal → does it cascade?) is the
    # scrambled-truth arm, exactly as IGNITE-2/CHORUS ran it.
    for rep in range(args.nscr):
        gid = IG.SKIES[rep % len(IG.SKIES)]
        ns = base_scr + rep
        plan.append(dict(gid=gid, ns=ns, scr=True, sm="truth", sc=0.0))

    nskip = 0
    for e in plan:
        try:
            tag, msg = run_kindle_realisation(
                jnp, IG, C, B1, TE, F, P, E_cache, h, tier, fl, e["gid"], skies[e["gid"]],
                e["ns"], e["ns"] + 10_000_000, scrambled=e["scr"],
                start_mode=e["sm"], start_scale=e["sc"], prefix="pilot")
            print(f"  {tag:52s} {msg}", flush=True)
        except RuntimeError as ex:
            # blast-radius discipline: one pathological realisation must not kill the shard
            nskip += 1
            print(f"  SKIP (non-finite) g{e['gid']} n{e['ns']} {e['sm']} scr={e['scr']}: "
                  f"{str(ex)[:80]}", flush=True)
    print(f"\nPILOT COMPLETE ({nskip} skipped)", flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["g0", "pilot"])
    ap.add_argument("--nreps", type=int, default=6)
    ap.add_argument("--nscr", type=int, default=3)
    ap.add_argument("--scale", type=float, default=1.0)
    raise SystemExit({"g0": mode_g0, "pilot": mode_pilot}[ap.parse_args().mode](ap.parse_args()))
