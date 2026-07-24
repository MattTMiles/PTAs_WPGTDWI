#!/usr/bin/env python3
"""Cross-host identity gate for the T=30 L_gwb bank (run on the CANDIDATE host, e.g. dgx A100).

Loads the banked square root (strict -- no recompute fallback), rebuilds the GLACIER T=30
problem, draws the reference seeds, and compares against the realisations banked by the
cutting host. With L frozen, the only remaining cross-host freedom is BLAS matmul
reduction order: PASS iff max|d residual| < 1e-10 absolute on every reference seed
(measured headroom expected ~1e-13). Also times the draw and prints the lane tag.

PASS here == draw provenance is HOST-FREE for this bank: floors and scores cut on any
host that passes this gate share the same null realisations to fp, and the campaign's
host-pinning reduces to the cpus=8 affinity convention (BLAS threads never touch the
draw path with a banked L, but the affinity refusal stays as declared provenance).
"""
import os, sys
import numpy as np
import socket, time

HSYMT = "/home/mattm/projects/HSYMT"
sys.path.insert(0, f"{HSYMT}/hpc_harbor/glacier")
import argparse
_ap = argparse.ArgumentParser(); _ap.add_argument("--t", type=int, default=30)
T_CUT = _ap.parse_args().t
BANK = f"{HSYMT}/CW_transition/b1_L_gwb_T{T_CUT}.npz"
REFS = f"{HSYMT}/CW_transition/b1_L_gwb_T{T_CUT}_refdraws.npz"

import glacier_pop as GP
from glacier_pop import draw_population, build_glacier_problem
sys.path.insert(0, f"{HSYMT}/CW_transition")
import trackB_b1_core as C
import jax
import jax.numpy as jnp

TOL = 1e-10


def main():
    print("=" * 92)
    print(f"L_gwb T={T_CUT} CROSS-HOST GATE on host={socket.gethostname()}")
    print(f"  devices: {jax.devices()}")
    print("=" * 92, flush=True)
    z = np.load(REFS)
    print(f"  reference host: {z['host']}, bank fp {z['bank_fp']}", flush=True)

    pop = draw_population(GP.SEED_POP_BASE, n_src=32, band_log10f=(-8.7, -7.5))
    G = build_glacier_problem(T_CUT, pop, verbose=True)
    amo = G["amo"]
    ones = jnp.ones(amo["npsr"])
    theta = np.asarray(amo["theta_truth"], float)
    data0 = amo["inject_delay"](jnp.asarray(theta), ones)
    sp = C.B1Split(amo, theta, data0, np.ones(amo["npsr"]))
    nd = C.NoiseDrawer(sp, amo, lgwb_path=BANK, strict=True)
    print(f"  drawer provenance: {nd.lgwb_prov}", flush=True)

    ok = True
    for s in [int(x) for x in z["seeds"]]:
        t0 = time.time()
        r = nd.draw(s, components=("white", "rn", "gwb"))
        flat = np.concatenate([np.asarray(x, float) for x in r])
        ref = np.asarray(z[f"ref_{s}"], float)
        d = float(np.max(np.abs(flat - ref)))
        b = d < TOL
        ok &= b
        print(f"  seed {s}: max|d residual| = {d:.3e} (<{TOL:.0e}) "
              f"[draw {time.time()-t0:.2f}s] -> {'PASS' if b else 'FAIL'}", flush=True)

    print(f"\n=== CROSS-HOST L_gwb GATE: {'PASS -- draws are host-free' if ok else 'FAIL'} "
          f"===", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
