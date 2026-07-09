"""STEP 1b follow-up — the compaction gate FAILED identically after the nan_to_num fix (no-op ->
no NaN configs). Isolate the real second defect. At the F1a state, with IDENTICAL 72 nonzero
configs, test whether the full-vs-compact difference is (i) config-set mismatch, (ii) run-to-run
nondeterminism, or (iii) FP summation-ORDER sensitivity amplified by the ~4e10-conditioned Adam
M-step. Tests: 1-Adam-step (near-gradient) match; determinism (full twice); order sensitivity
(full with rows shuffled); config-set identity (numpy).
Run: python trackB_step1_diag.py
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
    print("=== STEP1b-diag: real cause of full-vs-compact M-step disagreement ===", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False); prior = P["priors"]["lit"]
    mstep = PP.make_soft_mstep(P)
    scan = TE.seed_scan(P); seeds = TE.pick_seeds(scan, twoF_thresh=25.0)
    src0 = TE.seeds_to_phi(P, seeds)[:]
    base = P["L0"].copy()
    post, _ = PP.soft_estep(P, np.concatenate([base, src0]), prior)
    base = post["map_evalL"].copy()
    dcfg_full, w_full, _ = PP.build_mstep_configs(P, post, base)
    dcfg_cmp, w_cmp, nz = compact_configs(dcfg_full, w_full)
    print(f"  state ready {time.time()-t0:.0f}s; full={len(w_full)} nz={int(np.sum(w_full>0))} compact={len(w_cmp)}", flush=True)

    # (0) config-set identity: is compact's 72 exactly the full's 72 nonzero (as sets)?
    nzmask = w_full > 0
    full_nz_w = np.sort(w_full[nzmask]); cmp_w = np.sort(w_cmp[w_cmp > 0])
    same_w = full_nz_w.shape == cmp_w.shape and np.allclose(full_nz_w, cmp_w)
    # match dist rows by weight
    print(f"  [0] config sets identical (weights): {same_w}; "
          f"full_nz sum(w)={w_full[nzmask].sum():.6e} compact sum(w)={w_cmp.sum():.6e} "
          f"|dsum|={abs(w_full[nzmask].sum()-w_cmp.sum()):.2e}", flush=True)

    def run(dc, w, ns): return mstep(src0, dc, w, ns, 1.5e-2)

    # (1) 1-Adam-step (dominated by the single gradient): full vs compact
    sfa, qfa = run(dcfg_full, w_full, 1); sca, qca = run(dcfg_cmp, w_cmp, 1)
    d1 = np.abs(np.asarray(sfa) - np.asarray(sca)).max()
    print(f"  [1] 1-step  |dsrc| full-vs-compact = {d1:.3e}  (|dQ|={abs(qfa-qca):.2e})", flush=True)

    # (2) determinism: full twice (25 steps)
    sf1, qf1 = run(dcfg_full, w_full, 25); sf2, qf2 = run(dcfg_full, w_full, 25)
    d2 = np.abs(np.asarray(sf1) - np.asarray(sf2)).max()
    print(f"  [2] determinism full-vs-full(25) |dsrc| = {d2:.3e} (|dQ|={abs(qf1-qf2):.2e})", flush=True)

    # (3) order sensitivity: full with ROWS SHUFFLED (same set, different scan order), 25 steps
    rng = np.random.default_rng(0); perm = rng.permutation(len(w_full))
    sfs, qfs = run(dcfg_full[perm], w_full[perm], 25)
    d3 = np.abs(np.asarray(sf1) - np.asarray(sfs)).max()
    print(f"  [3] ORDER-sensitivity full-vs-shuffled(25) |dsrc| = {d3:.3e} (|dQ|={abs(qf1-qfs):.2e})", flush=True)

    # (4) full vs compact at 25 (reproduce F1a)
    scmp, qcmp = run(dcfg_cmp, w_cmp, 25)
    d4 = np.abs(np.asarray(sf1) - np.asarray(scmp)).max()
    print(f"  [4] full-vs-compact(25) |dsrc| = {d4:.3e} (|dQ|={abs(qf1-qcmp):.2e})", flush=True)

    print(f"\n  DIAGNOSIS: config sets identical={same_w}; 1-step diff={d1:.1e}; determinism={d2:.1e}; "
          f"order-alone={d3:.1e}; full-v-compact={d4:.1e}", flush=True)
    if d2 < 1e-10 and d3 > 1e-3:
        print("  -> DEFECT = FP summation-ORDER sensitivity amplified by the ill-conditioned Adam "
              "M-step (chaotic). Compaction is not the bug; the M-step is numerically unstable.", flush=True)
    elif d1 > 1e-6:
        print("  -> gradients differ at step 1 -> config/weight handling differs (not pure order).", flush=True)
    np.savez("/home/mattm/projects/HSYMT/CW_transition/trackB_step1_diag.npz",
             same_w=same_w, d1=d1, d2=d2, d3=d3, d4=d4)
    print(f"  saved ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}", flush=True)
    sys.exit(main())
