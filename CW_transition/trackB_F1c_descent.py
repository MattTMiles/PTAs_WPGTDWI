"""F1c — trace the L0 descent leg. Re-run the truth-blind soft-EM (deterministic: same seeds,
LR, compaction as trackB_L0_float) and log, PER ITERATION, each loud source's LABEL-MATCHED
(best 3x3 assignment) registration offset: sky (deg) and scaled sky/freq. Saves the per-iter
source stack. Answers: did the loud sources hold near the F-stat seed or walk off, and WHEN.
Registration-critical params only (sky idx 0,1; freq idx 4) -- amplitude/phase/psi excluded
(extrinsic, they do NOT move the pulsar-term fringe registration).
Run: python trackB_F1c_descent.py --max_iter 20  (jug-gpu, background)
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import argparse, sys, time, itertools
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
import trackB_p2_pipeline as PP
from trackB_L0_float import compact_configs
from trackB_p1_map import ang_deg
CWT = "/home/mattm/projects/HSYMT/CW_transition"


def loud_registration_offsets(src, src_t, scale, n_loud):
    """Best 3x3 assignment by sky angle; return per-assigned-loud (sky_deg, scaled sky, scaled freq)."""
    def U(s, k):
        cg, gp = s[8*k], s[8*k+1]; import numpy as _np
        ss = _np.sqrt(max(1-cg**2, 0)); return _np.array([ss*_np.cos(gp), ss*_np.sin(gp), cg])
    M = np.array([[np.degrees(np.arccos(np.clip(np.dot(U(src,i), U(src_t,j)), -1, 1)))
                   for j in range(n_loud)] for i in range(n_loud)])
    best = min(itertools.permutations(range(n_loud)), key=lambda p: sum(M[i, p[i]] for i in range(n_loud)))
    rows = []
    for i in range(n_loud):
        j = best[i]
        sky_deg = M[i, j]
        sc_sky = max(abs(src[8*i+0]-src_t[8*j+0])/scale[0], abs(src[8*i+1]-src_t[8*j+1])/scale[1])
        sc_frq = abs(src[8*i+4]-src_t[8*j+4])/scale[4]
        rows.append((j, sky_deg, sc_sky, sc_frq))
    return best, rows


def main(max_iter=20):
    print("=== F1c descent trace (label-matched registration offsets per iter) ===", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False); prior = P["priors"]["lit"]
    ndist = P["n_dist"]; scale = TE.phi_scale(P); src_t = P["theta_truth"][ndist:]
    n_loud = TE.CONFIGS["pop"]["population"][0]
    from trackB_split import Split
    split = Split(P["amo_full"], P["data_obs"], P["theta_truth"])
    mstep = PP.make_soft_mstep(P)
    scan = TE.seed_scan(P); seeds = TE.pick_seeds(scan, twoF_thresh=25.0)
    src = TE.seeds_to_phi(P, seeds)[:]
    best0, rows0 = loud_registration_offsets(src, src_t, scale, n_loud)
    print(f"  SEED (it-pre): perm={best0}", flush=True)
    for i, (j, sd, ss, sf) in enumerate(rows0):
        print(f"    fit{i}->truth{j}: sky={sd:.2f}deg scaled_sky={ss:.4f} scaled_freq={sf:.4f}", flush=True)
    base = P["L0"].copy(); LR = lambda it: 1.5e-2 if it < 8 else (6e-3 if it < 18 else 2e-3)
    src_stack = [src.copy()]
    for it in range(max_iter):
        post, LNL = PP.soft_estep(P, np.concatenate([base, src]), prior, split=split)
        base = post["map_evalL"].copy()
        dcfg, w, _ = PP.build_mstep_configs(P, post, base)
        dcfg, w, _ = compact_configs(dcfg, w)
        src, qval = mstep(src, dcfg, w, 25, LR(it))
        src_stack.append(src.copy())
        best, rows = loud_registration_offsets(src, src_t, scale, n_loud)
        summ = " | ".join(f"f{i}->t{j} sky={sd:.1f} sc={ss:.3f}/{sf:.3f}" for i,(j,sd,ss,sf) in enumerate(rows))
        print(f"  it{it:2d} ({time.time()-t0:.0f}s) perm={best}: {summ}", flush=True)
    np.savez(f"{CWT}/trackB_F1c_descent.npz", src_stack=np.array(src_stack),
             src_truth=src_t, scale=scale, n_loud=n_loud)
    print(f"  saved trackB_F1c_descent.npz ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}", flush=True)
    ap = argparse.ArgumentParser(); ap.add_argument("--max_iter", type=int, default=20)
    a = ap.parse_args(); sys.exit(main(a.max_iter))
