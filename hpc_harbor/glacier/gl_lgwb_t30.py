#!/usr/bin/env python3
"""Cut the CANONICAL T=30 GLACIER L_gwb bank -- the host-freedom fix.

WHY (Matt's point, 2026-07-23): the physics is hardware-invariant; what varies is
np.linalg.eigh's arbitrary eigenvector basis inside Phi_gwb's near-degenerate subspaces
(backend/thread/host dependent). The T=15 campaign bank (b1_L_gwb.npz) already froze that
freedom; T=30 never got a bank, so GLACIER's fan sat on the RECOMPUTED-UNSAFE branch and
inherited host-pinning it doesn't need. This script cuts the T=30 bank ONCE; thereafter
draw(seed) is host-free up to BLAS matmul reduction order (~1e-13 absolute on the
residuals -- measured by the cross-host gate, gl_lgwb_xhost.py, which this script arms
with reference draws).

Phi_gwb depends ONLY on the GP geometry (npsr=116, ngp_gwb(T=30)=23, span), NOT on the
drawn census -- one bank serves every fan cell. Built here from a minimal n=32 census.

MUST run with BLAS threads pinned to 8 BEFORE numpy import (sbatch exports them): the
banked basis is then the cpus=8 convention, continuous with every other bank in the repo.
"""
import os, sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    if os.environ.get(_v) != "8":
        print(f"REFUSED: {_v}={os.environ.get(_v)!r}, need '8' set before python starts "
              "(the banked convention; export in the sbatch).", flush=True)
        sys.exit(2)

import numpy as np
import hashlib, socket

HSYMT = "/home/mattm/projects/HSYMT"
sys.path.insert(0, f"{HSYMT}/hpc_harbor/glacier")
import argparse
_ap = argparse.ArgumentParser(); _ap.add_argument("--t", type=int, default=30)
T_CUT = _ap.parse_args().t
BANK = f"{HSYMT}/CW_transition/b1_L_gwb_T{T_CUT}.npz"

import glacier_pop as GP
from glacier_pop import draw_population, build_glacier_problem
sys.path.insert(0, f"{HSYMT}/CW_transition")
import trackB_b1_core as C
import trackB_estimator  # noqa: F401  (path warm)
import jax.numpy as jnp


def main():
    print("=" * 92)
    print(f"CANONICAL b1_L_gwb_T{T_CUT} -- the frozen GWB square root for the GLACIER venue")
    print(f"  host={socket.gethostname()} numpy {np.__version__} "
          f"OMP={os.environ['OMP_NUM_THREADS']}")
    print("=" * 92, flush=True)

    pop = draw_population(GP.SEED_POP_BASE, n_src=32, band_log10f=(-8.7, -7.5))
    G = build_glacier_problem(T_CUT, pop, verbose=True)
    amo = G["amo"]
    it = amo["internals"]
    ones = jnp.ones(amo["npsr"])
    theta = np.asarray(amo["theta_truth"], float)
    data0 = amo["inject_delay"](jnp.asarray(theta), ones)
    sp = C.B1Split(amo, theta, data0, np.ones(amo["npsr"]))

    # Phi_gwb exactly as NoiseDrawer.__init__ builds it
    Pinv = np.asarray(it["Pinv_gwb"])
    Phi = np.linalg.inv(0.5 * (Pinv + Pinv.T))
    Phi = 0.5 * (Phi + Phi.T)
    print(f"  Phi_gwb: {Phi.shape} (npsr {amo['npsr']} x ngp {sp.ngp_gwb})", flush=True)

    w, V = np.linalg.eigh(Phi)
    w = np.clip(w, 0.0, None)
    L = V * np.sqrt(w)
    fp = C.lgwb_fingerprint(L)
    recon = float(np.max(np.abs(L @ L.T - Phi)))
    scale = float(np.max(np.abs(Phi)))
    print(f"  eigh: L fp={fp}; max|LL^T - Phi| = {recon:.3e} (scale {scale:.3e})", flush=True)

    np.savez(BANK, L_gwb=L, fingerprint=fp, cpus="8",
             host=socket.gethostname(), T_label=T_CUT,
             npsr=amo["npsr"], ngp=sp.ngp_gwb,
             note=f"GLACIER T={T_CUT} canonical GWB square root; geometry-only (census-free). "
                  "Cut at cpus=8 on the hgx03 lane, 2026-07-23. Loaded via "
                  "NoiseDrawer(lgwb_path=..., strict=True); gated by gl_lgwb_xhost.py.")
    print(f"  banked -> {BANK}", flush=True)

    # arm the cross-host gate: reference draws THROUGH the bank on THIS host
    nd = C.NoiseDrawer(sp, amo, lgwb_path=BANK, strict=True)
    print(f"  drawer provenance: {nd.lgwb_prov}", flush=True)
    seeds = [9_000, 9_100, 424_242]
    refs, hashes = {}, {}
    for s in seeds:
        r = nd.draw(s, components=("white", "rn", "gwb"))
        flat = np.concatenate([np.asarray(x, float) for x in r])
        refs[f"ref_{s}"] = flat
        hashes[s] = hashlib.sha256(flat.tobytes()).hexdigest()[:16]
        print(f"  ref draw seed {s}: n={flat.size} sha={hashes[s]}", flush=True)
    np.savez(f"{HSYMT}/CW_transition/b1_L_gwb_T{T_CUT}_refdraws.npz",
             seeds=np.array(seeds), bank_fp=fp,
             host=socket.gethostname(),
             **refs,
             note=f"Reference realisations drawn through b1_L_gwb_T{T_CUT}.npz on the cutting "
                  "host. Cross-host gate: same draws elsewhere must match to <1e-10 abs "
                  "(BLAS matmul reduction order is the only remaining freedom).")
    print("  reference draws banked (b1_L_gwb_T30_refdraws.npz)", flush=True)
    print(f"\n=== T={T_CUT} L_gwb BANK CUT ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
