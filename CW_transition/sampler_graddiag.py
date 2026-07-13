"""SAMPLER — G1 anatomy: is the AD-vs-FD gap a CHAIN-RULE BUG or the target's own KINKS?

G1 (first run) showed: gradients agree to high relative precision on the large components but
the WORST-relative dim (src1_psi, |grad| ~ 1e-1 against max|grad| ~ 7e4) missed by 8 %, i.e. an
ABSOLUTE gap of ~1e-2 nat/unit. Two hypotheses, and they are separated by the FD step size:

  BUG   -- an error in the hand-assembled per-pulsar VJP. Then the ABSOLUTE gap is
           h-INDEPENDENT (it is a wrong derivative, not a discretisation error).
  KINK  -- lnL_marg's within-fringe peak is a segment MAX over the B_CERT grid, so the target
           is only piecewise-smooth; FD across a within-fringe argmax swap measures a chord,
           not the derivative. Then the gap SHRINKS with h (and the swap count -> 0), while AD
           reports the (correct) subgradient at the incumbent argmax.

Also reported, because it decides whether KINK is even possible: the per-pulsar GROUP SIZE
(grid points per fringe). If every fringe group holds exactly ONE grid point, segment_max is
the identity, the target is smooth, and KINK is REFUTED by construction.
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

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import sampler_core as S
import trackB_b1_core as C
import trackB_estimator as TE

RES = S.RES


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=float, nargs="+", default=[1e-4, 1e-5, 1e-6, 1e-7])
    ap.add_argument("--chunk", type=int, default=64)
    a = ap.parse_args()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)

    S.MargJax.F_CHUNK = a.chunk
    cell = S.build_cell(seed=None, tier="lit", noise=False, arm="A")
    P = cell["P"]
    src_truth = cell["theta_truth"][P["n_dist"]:].copy()
    idx = S.active_index(list(range(C.N_LOUD)), S.ACTIVE_AXES, P["ncw"])
    labels = S.active_labels(list(range(C.N_LOUD)), S.ACTIVE_AXES)
    M = S.MargJax(P, cell["data"], cell["mask"], idx, src_truth, tier="lit")

    # ---- IS SEGMENT_MAX EVEN ACTIVE? (group sizes: grid points per fringe) ----
    gid = np.asarray(M.gid)
    sizes = np.array([np.bincount(gid[p], minlength=M.Kmax)[:M.K[p]] for p in range(M.npsr)],
                     dtype=object)
    mx = np.array([int(s.max()) for s in sizes])
    mn = np.array([int(s.min()) for s in sizes])
    print(f"\n== fringe-group occupancy (B = {M.B} grid points per pulsar) ==", flush=True)
    print(f"   K (fringes in window): min {M.K.min()}, median {int(np.median(M.K))}, "
          f"max {M.K.max()}", flush=True)
    print(f"   grid points per fringe: max over pulsars = {mx.max()}, "
          f"pulsars with >1 point in some fringe: {int(np.sum(mx > 1))}/{M.npsr}", flush=True)
    kink_possible = bool(np.any(mx > 1))
    print(f"   -> segment_max ACTIVE (kinks possible): {kink_possible}", flush=True)

    x0 = src_truth[idx]
    scale = TE.phi_scale(P)[idx]
    t0 = time.time()
    v, gad = M.value_and_grad(x0)
    print(f"\n   value {v:.6f}; value_and_grad warm {time.time()-t0:.2f}s "
          f"(F_CHUNK = {a.chunk})", flush=True)

    print(f"\n== AD vs central FD, step ladder (absolute gap: h-independent => BUG) ==",
          flush=True)
    hdr = "   " + f"{'param':>22s} {'AD':>14s}" + "".join(f"{f'FD h={h:.0e}':>14s}"
                                                          for h in a.steps)
    print(hdr, flush=True)
    GF = {}
    for h in a.steps:
        g = np.zeros(len(x0))
        for j in range(len(x0)):
            xp = x0.copy(); xp[j] += h * scale[j]
            xm = x0.copy(); xm[j] -= h * scale[j]
            g[j] = (M.value(xp) - M.value(xm)) / (2 * h * scale[j])
        GF[h] = g
    for j, lab in enumerate(labels):
        row = f"   {lab:>22s} {gad[j]:14.6e}" + "".join(f"{GF[h][j]:14.6e}" for h in a.steps)
        print(row, flush=True)

    print(f"\n== absolute gap |AD - FD| per step ==", flush=True)
    for h in a.steps:
        d = np.abs(gad - GF[h])
        print(f"   h = {h:.0e}: max {d.max():.3e}, median {np.median(d):.3e}", flush=True)
    print(f"\n== relative gap (vs max(|AD|,|FD|)) per step ==", flush=True)
    for h in a.steps:
        den = np.maximum(np.abs(gad), np.abs(GF[h]))
        r = np.abs(gad - GF[h]) / np.where(den > 0, den, 1.0)
        print(f"   h = {h:.0e}: max {r.max():.3e} (dim {labels[int(np.argmax(r))]}), "
              f"median {np.median(r):.3e}", flush=True)

    np.savez(f"{RES}/grad_diag.npz", gad=gad, steps=np.array(a.steps),
             gfd=np.array([GF[h] for h in a.steps]), labels=np.array(labels),
             K=M.K, maxgroup=mx, kink_possible=kink_possible)
    print(f"\n   saved {RES}/grad_diag.npz", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
