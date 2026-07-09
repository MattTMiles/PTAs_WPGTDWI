"""STEP 2 — loud2 seeding repair. The single-source F-stat map lets loud0/loud1 sidelobes
(2F 285-350) outrank loud2's true peak (2F~218, was rank-13 > the top-16 cut). Repair: NMS with
a SKY-exclusion radius around each accepted peak (kills a strong source's cross-frequency
sidelobes, which sit at ~same sky / different freq) + a deeper candidate list. Data-driven ONLY;
truth used solely to REPORT capture (offsets, rank), never in selection.
Run: python trackB_step2_seeder.py
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
CWT = "/home/mattm/projects/HSYMT/CW_transition"


def pick_seeds_sky(scan, twoF_thresh=25.0, excl_sky_deg=15.0, min_df_bins=1.5):
    """Repaired NMS: reject a peak if it is within excl_sky_deg of a STRONGER accepted peak
    (sky-only -> kills cross-frequency sidelobes), OR within (min_sep small) AND df. Sky-only
    exclusion is the sidelobe killer. Returns seeds sorted by 2F."""
    twoF = scan["twoF"]; lfgws = scan["lfgws"]; cosgt = scan["cosgt"]; gwphi = scan["gwphi"]
    nf, ns = twoF.shape
    flat = [(twoF[i, j], i, j) for i in range(nf) for j in range(ns) if twoF[i, j] > twoF_thresh]
    flat.sort(reverse=True)
    seeds = []
    for val, i, j in flat:
        u = TE._unit(cosgt[j], gwphi[j]); ok = True
        for s in seeds:
            ang = np.degrees(np.arccos(np.clip(np.dot(u, s["u"]), -1, 1)))
            if ang < excl_sky_deg:                       # sky-only exclusion (any freq)
                ok = False; break
        if ok:
            fr = scan["free"][i, j]
            seeds.append(dict(twoF=val, lfgw=lfgws[i], cosgt=cosgt[j], gwphi=gwphi[j],
                              cos_inc=fr[0], log10_h=fr[1], phase0=fr[2], psi=fr[3],
                              u=u, log10_mc=TE.SEED_MC))
    return seeds


def report_capture(P, seeds, tag, cap_deg=25.0):
    louds = TE.loud_truth(P)
    print(f"  [{tag}] {len(seeds)} seeds; loud capture:", flush=True)
    ncap = 0
    for li, lt in enumerate(louds):
        angs = [np.degrees(np.arccos(np.clip(np.dot(lt["u"], s["u"]), -1, 1))) for s in seeds]
        r = int(np.argmin(angs)); a = angs[r]; cap = a < cap_deg
        ncap += cap
        df = abs(seeds[r]["lfgw"] - lt["lfgw"]) / np.median(np.diff(np.sort(np.unique(
             [s["lfgw"] for s in seeds])))) if len(seeds) > 1 else 0.0
        print(f"     loud{li}: nearest seed rank={r} ang={a:.2f}deg 2F={seeds[r]['twoF']:.1f} "
              f"{'CAPTURED' if cap else 'MISS'}", flush=True)
    print(f"     -> {ncap}/3 loud captured (<{cap_deg}deg)", flush=True)
    return ncap


def main():
    print("=== STEP 2 seeder repair (sky-exclusion NMS + depth) ===", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False)
    scan = TE.seed_scan(P)
    print(f"  scan done {time.time()-t0:.0f}s", flush=True)
    # baseline (current)
    seeds_cur = TE.pick_seeds(scan, twoF_thresh=25.0)
    report_capture(P, seeds_cur[:16], "CURRENT top-16", )
    # repaired: sweep exclusion radius, find smallest giving 3/3 in the top-16
    best = None
    for R in [8, 10, 12, 15, 18, 20, 25, 30]:
        seeds = pick_seeds_sky(scan, excl_sky_deg=float(R))
        nc = report_capture(P, seeds[:16], f"sky-excl R={R}deg top-16")
        if nc == 3 and best is None:
            best = (R, seeds)
    if best is None:
        print("  !! no exclusion radius gave 3/3 in top-16 -> STOP (deeper defect)", flush=True)
    else:
        R, seeds = best
        print(f"\n  CHOSEN excl_sky_deg={R} -> 3/3 loud in top-16", flush=True)
        np.savez(f"{CWT}/trackB_step2_seeder.npz", excl_sky_deg=R,
                 seed_2F=np.array([s["twoF"] for s in seeds[:16]]),
                 seed_cosgt=np.array([s["cosgt"] for s in seeds[:16]]),
                 seed_gwphi=np.array([s["gwphi"] for s in seeds[:16]]))
        print(f"  saved trackB_step2_seeder.npz ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}", flush=True)
    sys.exit(main())
