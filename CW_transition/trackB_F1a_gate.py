"""F1a — retroactive GATE on the 7.7x M-step compaction (compact_configs drops w=0 rows).
One M-step from an identical state, FULL 348-config scan vs COMPACT 72-config scan:
report max|Delta| on the M-step objective Q and on every updated source param. Claim of
'mathematically exact' stands only with a number.
Run: python trackB_F1a_gate.py
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
import trackB_p2_pipeline as PP
from trackB_L0_float import compact_configs


def main():
    print("=== F1a compaction gate (full 348 vs compact 72, one M-step, same state) ===", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False); prior = P["priors"]["lit"]
    mstep = PP.make_soft_mstep(P); ndist = P["n_dist"]
    # representative state: F-stat seed source + prior-mean distances, one E-step
    scan = TE.seed_scan(P); seeds = TE.pick_seeds(scan, twoF_thresh=25.0)
    src0 = TE.seeds_to_phi(P, seeds)[:]
    base = P["L0"].copy()
    post, LNL = PP.soft_estep(P, np.concatenate([base, src0]), prior)
    base = post["map_evalL"].copy()
    dcfg_full, w_full, nanch = PP.build_mstep_configs(P, post, base)
    dcfg_cmp, w_cmp, nz = compact_configs(dcfg_full, w_full)
    print(f"  state ready ({time.time()-t0:.0f}s): anchors={nanch}, full C={len(w_full)} "
          f"nonzero={int(np.sum(w_full>0))}, compact C={len(w_cmp)}", flush=True)
    # identical M-step (same src0, nsteps, lr) on each config set
    src_full, q_full = mstep(src0, dcfg_full, w_full, 25, 1.5e-2)
    src_cmp,  q_cmp  = mstep(src0, dcfg_cmp,  w_cmp,  25, 1.5e-2)
    dq = abs(q_full - q_cmp)
    dsrc = np.abs(np.asarray(src_full) - np.asarray(src_cmp))
    print(f"\n  Q_full={q_full:.6f} Q_compact={q_cmp:.6f}  |dQ|={dq:.3e}  rel={dq/(abs(q_full)+1e-9):.3e}", flush=True)
    print(f"  updated-source max|Delta|={dsrc.max():.3e} median={np.median(dsrc):.3e} "
          f"(argmax param {int(dsrc.argmax())})", flush=True)
    scale = TE.phi_scale(P)
    print(f"  updated-source max|Delta|/scale={np.max(dsrc/scale):.3e}", flush=True)
    ok = dsrc.max() < 1e-6 and dq < 1e-3
    print(f"\n  [F1a GATE] {'PASS' if ok else 'CHECK'}: compaction {'exact to' if ok else 'differs by'} "
          f"|dQ|={dq:.2e}, max|dsrc|={dsrc.max():.2e} ({time.time()-t0:.0f}s)", flush=True)
    np.savez("/home/mattm/projects/HSYMT/CW_transition/trackB_F1a_gate.npz",
             q_full=q_full, q_cmp=q_cmp, dq=dq, dsrc=dsrc, nz=nz, nanch=nanch)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}", flush=True)
    sys.exit(main())
