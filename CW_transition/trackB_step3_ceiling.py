"""STEP 3 analysis — the repaired-float REGISTRATION CEILING (the wall's height).
Loads trackB_L0_float.npz (repaired-seeder + hygiene float). For each TRUTH loud source, match
to the NEAREST among ALL N_CW fitted sources (not just the first 3 slots -- loud2's seed sits at
slot ~11), and report the registration-critical (sky, freq) scaled offset + the extrinsics
(log10_h, phase0, psi) SEPARATELY. Also correlate the F2 loosest rungs with pulsar distance
(design-theorem lever i: tol ~ 1/L).
Run: python trackB_step3_ceiling.py
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import sys, itertools
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
CWT = "/home/mattm/projects/HSYMT/CW_transition"


def U(s, k):
    cg, gp = s[8*k], s[8*k+1]; ss = np.sqrt(max(1-cg**2, 0))
    return np.array([ss*np.cos(gp), ss*np.sin(gp), cg])


def main():
    d = np.load(f"{CWT}/trackB_L0_float.npz", allow_pickle=True)
    src = d["src"]; src_t = d["src_truth"]; status = str(d["status"]); it = int(d["iter"]) if "iter" in d else -1
    P = TE.build_problem("pop", build_earth=False, verbose=False)
    scale = TE.phi_scale(P); ncw = P["ncw"]; n_loud = TE.CONFIGS["pop"]["population"][0]
    print(f"=== STEP 3 CEILING (float status={status}, iter={it}) ===", flush=True)
    print(f"  match each TRUTH loud -> nearest of ALL {ncw} fitted sources:", flush=True)
    reg_sky = []; reg_frq = []
    for j in range(n_loud):
        ut = U(src_t, j)
        angs = [np.degrees(np.arccos(np.clip(np.dot(U(src, i), ut), -1, 1))) for i in range(ncw)]
        i = int(np.argmin(angs)); a = angs[i]
        sc_sky = max(abs(src[8*i+0]-src_t[8*j+0])/scale[0], abs(src[8*i+1]-src_t[8*j+1])/scale[1])
        sc_frq = abs(src[8*i+4]-src_t[8*j+4])/scale[4]
        ext = {p: abs(src[8*i+q]-src_t[8*j+q])/scale[q] for q, p in [(5,'h'),(6,'phase'),(7,'psi')]}
        reg_sky.append(sc_sky); reg_frq.append(sc_frq)
        print(f"    truth loud{j} <- fit slot{i}: sky={a:.2f}deg scaled_sky={sc_sky:.4f} "
              f"scaled_freq={sc_frq:.4f} | extrinsic scaled h={ext['h']:.2f} phase={ext['phase']:.2f} psi={ext['psi']:.2f}", flush=True)
    reg_max = max(max(reg_sky), max(reg_frq)); reg_med = np.median(reg_sky + reg_frq)
    print(f"\n  REGISTRATION CEILING (repaired float): max sky/freq scaled = {reg_max:.4f}; "
          f"median = {reg_med:.4f}; per-loud sky = {[round(x,4) for x in reg_sky]}", flush=True)
    loosest = 1.85e-3
    gap_dex = np.log10(reg_max / loosest)
    print(f"  F2 loosest rung = {loosest:.2e} scaled -> WALL = float_ceiling/rung = "
          f"{reg_max/loosest:.1f}x = {gap_dex:.2f} dex (ZERO rungs in the gap)", flush=True)

    # design-theorem lever (i): loosest F2 rungs vs pulsar distance
    L0 = np.asarray(P["L0"]); names = P["names"]
    f2 = np.load(f"{CWT}/trackB_F2_ladder.npz", allow_pickle=True)
    rows = f2["rows"]; ts = f2["tol_sky"]; order = np.argsort(-ts)
    print(f"\n  LEVER (i) tol~1/L check -- loosest 6 rungs vs pulsar distance:", flush=True)
    for k in order[:6]:
        pa = int(rows[k, 1])
        print(f"    tol={ts[k]:.3e} psr {names[pa]} L={L0[pa]:.3f} kpc", flush=True)
    print(f"    (median pulsar L = {np.median(L0):.3f} kpc; nearest = {L0.min():.3f}, "
          f"farthest = {L0.max():.3f})", flush=True)
    np.savez(f"{CWT}/trackB_step3_ceiling.npz", reg_sky=np.array(reg_sky), reg_frq=np.array(reg_frq),
             reg_max=reg_max, gap_dex=gap_dex, loosest=loosest)
    print(f"  saved trackB_step3_ceiling.npz", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
