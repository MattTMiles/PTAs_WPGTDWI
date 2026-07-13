"""SAMPLER — G1b: the CHAIN-RULE UNIT TEST (small array, Richardson FD).

WHY THIS EXISTS
---------------
On the full 116-pulsar array, autodiff-vs-FD cannot be certified at 1e-8, and the reason is
the target itself, not the gradient:

  * lnL_marg is assembled as  sum_p M_p - (npsr-1)*lnL_ref  -- a difference of ~4.7e7-sized
    quantities yielding ~4e5, so the VALUE carries ~1e-7 nat of cancellation noise (the same
    ~5e-8 that G0 measures against B1Marg). Central FD then carries ~eps/h of gradient noise.
  * the target's condition number is ~1e10-1e11: the registration dims (log10_fgw, log10_mc)
    need h <~ 1e-7 scaled before FD even resolves their curvature, while the flat dims
    (psi/phase0/log10_h/cos_inc) are already roundoff-dominated by h = 1e-6.
  MEASURED (sampler_graddiag): FD converges ONTO AD as h -> 0 in the sharp dims
  (src0_log10_mc: -6.1e3 -> 5.05e4 -> 6.95e4 -> 6.9703e4 vs AD 6.9705e4) and AWAY from it in
  the flat dims. No single step certifies all 18 dims. That is the pre-registered 4e10 hazard
  showing up IN THE GATE.

So the chain rule is gated where it CAN be gated cleanly: the SAME MargJax code path on a
small sub-array (fewer pulsars => far less cancellation in the value), against a
RICHARDSON-extrapolated central difference (error O(h^4) instead of O(h^2)). The assembly
being tested -- per-pulsar VJP + base VJP + the -(npsr-1)*lnref term + the fringe scan/segment
reduction -- is IDENTICAL code, exercised identically; only the array size changes.

PASS BAR: max relative error < 1e-8 over all active dims (the mission's bar), on the
best-converged Richardson estimate.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time, argparse
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")

import sampler_core as S
import trackB_b1_core as C
import trackB_estimator as TE
import stagec_anchor_a2 as A2
from cw_helpers import load_pulsars, build_enterprise_pta, compute_mode_spacing, \
    _enterprise_to_disco_params
from stagec_fisher import make_geometry_injection, N_COMPONENTS, LOG10_EQUAD, CW_PARAM_BASE

RES = S.RES


def build_small(npsr=8, ncw=2, n_loud=1, seed=3000, verbose=True):
    """build_b1_problem, shrunk. SAME code path (B1Split, MaskedDelay, FringeTables)."""
    t0 = time.time()
    ent_psrs, disco_psrs = load_pulsars(npsr)
    names = [p.name for p in ent_psrs]
    priors = A2.load_real_priors(names)
    pta, cwb, _ = build_enterprise_pta(ent_psrs, ncw, components=N_COMPONENTS,
                                       log10_equad=LOG10_EQUAD)
    inj = make_geometry_injection(pta, ent_psrs, ncw, cwb, seed=seed,
                                  population=(n_loud, -13.25, -14.25))
    amo = C.build_b1_amortised(disco_psrs, ncw, inj, cwb)
    temp = _enterprise_to_disco_params(inj, cwb)
    theta_truth = np.array([temp[k] for k in amo["theta_keys"]], dtype=np.float64)
    L0 = np.array([pe.pdist[0] for pe in ent_psrs])
    idx_of = {k: i for i, k in enumerate(amo["theta_keys"])}
    AI = np.array([idx_of[f"{pe.name}_cw_p_dist"] for pe in ent_psrs])
    cw_list = [{k: inj[f"{nm}_{k}"] for k in CW_PARAM_BASE} for nm in cwb]
    dL = np.array([min(compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"],
                                            cw["log10_fgw"], disco_psrs[a].pos)
                       for cw in cw_list) for a in range(npsr)])
    ones = jnp.ones(npsr)
    data_zero = amo["inject_delay"](jnp.asarray(theta_truth), ones)
    _ = amo["logL"](jnp.asarray(theta_truth), data_zero, ones)
    sp = C.B1Split(amo, theta_truth, data_zero, np.ones(npsr))
    P = dict(ent_psrs=ent_psrs, disco_psrs=disco_psrs, names=names, priors=priors,
             pta=pta, cwb=cwb, inj=inj, ncw=ncw, npsr=npsr, amo=amo, sp=sp,
             nd=C.NoiseDrawer(sp, amo), theta_truth=theta_truth, L0=L0, AI=AI,
             dL_truth=dL, n_dist=amo["n_dist"], n_theta=amo["n_theta"],
             mask_one=np.ones(npsr), mask_zero=np.zeros(npsr))
    P["EV_truth"] = C.build_EV(P, dL)
    if verbose:
        print(f"  build_small(npsr={npsr}, ncw={ncw}): {time.time()-t0:.0f}s", flush=True)
    return P


def richardson(f, x, j, h, scale):
    """Central difference + Richardson: D = (4*C(h/2) - C(h))/3, error O(h^4)."""
    def C_(hh):
        xp = x.copy(); xp[j] += hh * scale
        xm = x.copy(); xm[j] -= hh * scale
        return (f(xp) - f(xm)) / (2 * hh * scale)
    c1 = C_(h)
    c2 = C_(h / 2)
    return (4 * c2 - c1) / 3.0, c1, c2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npsr", type=int, default=8)
    ap.add_argument("--ncw", type=int, default=2)
    ap.add_argument("--chunk", type=int, default=64)
    ap.add_argument("--steps", type=float, nargs="+", default=[1e-3, 1e-4, 1e-5, 1e-6])
    a = ap.parse_args()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    print(f"=== G1b: chain-rule unit test, npsr={a.npsr}, ncw={a.ncw} ===", flush=True)

    S.MargJax.F_CHUNK = a.chunk
    P = build_small(a.npsr, a.ncw)
    src_truth = P["theta_truth"][P["n_dist"]:].copy()
    idx = S.active_index([0], S.ACTIVE_AXES, P["ncw"])
    labels = S.active_labels([0], S.ACTIVE_AXES)
    data = P["amo"]["inject_delay"](jnp.asarray(P["theta_truth"]),
                                    jnp.asarray(P["mask_one"]))
    M = S.MargJax(P, data, P["mask_one"], idx, src_truth, tier="lit")

    # perturb off truth so the gradient is not ~0 (a zero gradient cannot be relatively gated)
    rng = np.random.default_rng(3)
    scale = TE.phi_scale(P)[idx]
    x0 = src_truth[idx] + 1e-3 * scale * rng.normal(size=len(idx))

    v, gad = M.value_and_grad(x0)
    print(f"\n  lnL_marg = {v:.6f}; |grad| in [{np.abs(gad).min():.3e}, "
          f"{np.abs(gad).max():.3e}]", flush=True)

    # ---- G1a: AD vs AD. The hand-assembled per-pulsar VJP vs the monolithic jax.grad of the
    # SAME function. No finite differences => no FD noise => an EXACT gate on the chain rule.
    vm, gm = M.monolith_value_and_grad(x0)
    dv = abs(vm - v)
    rel_ad = np.abs(gad - gm) / np.maximum(np.abs(gm), 1e-300)
    ok_a = bool(np.max(rel_ad) < 1e-10) and dv < 1e-6
    print(f"\n== G1a: per-pulsar VJP  vs  monolithic jax.grad (AD vs AD, exact) ==", flush=True)
    print(f"   |value diff|          = {dv:.3e}", flush=True)
    print(f"   max rel grad diff     = {np.max(rel_ad):.3e}  "
          f"(worst dim {labels[int(np.argmax(rel_ad))]})", flush=True)
    print(f"   [G1a] -> {'PASS' if ok_a else 'FAIL'} (bar 1e-10)", flush=True)

    print(f"\n  {'param':>18s} {'AD':>15s} {'Richardson FD':>15s} {'rel err':>10s} {'h*':>8s}",
          flush=True)
    worst = 0.0
    best = np.zeros(len(idx)); besth = np.zeros(len(idx))
    for j, lab in enumerate(labels):
        rows = []
        for h in a.steps:
            R, c1, c2 = richardson(M.value, x0, j, h, scale[j])
            rows.append((abs(R - gad[j]) / max(abs(gad[j]), 1e-300), R, h))
        rows.sort()
        rel, R, h = rows[0]
        best[j] = R; besth[j] = h
        worst = max(worst, rel)
        print(f"  {lab:>18s} {gad[j]:15.8e} {R:15.8e} {rel:10.2e} {h:8.0e}", flush=True)
    ok = worst < 1e-8
    print(f"\n  [G1b] worst relative error = {worst:.3e} -> {'PASS' if ok else 'FAIL'} "
          f"(bar 1e-8)", flush=True)
    np.savez(f"{RES}/grad_unit.npz", gad=gad, fd=best, h=besth, labels=np.array(labels),
             worst=worst, npsr=a.npsr, ncw=a.ncw, ok=ok)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
