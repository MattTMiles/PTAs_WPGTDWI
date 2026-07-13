"""Generate the CANONICAL b1_L_gwb.npz — the frozen GWB square root (thread-hazard fix).

WHY THIS FILE MUST BE MADE ON ACCRE AT 8 THREADS
------------------------------------------------
Phi_gwb is HD-correlated and NEAR-DEGENERATE. `np.linalg.eigh` fixes the eigenvector basis
inside a degenerate subspace arbitrarily, and LAPACK's choice depends on the BLAS THREAD
COUNT. Every banked noisy realisation in this repo was drawn at `--cpus-per-task=8`, so the
canonical basis is the one eigh returns AT 8 THREADS ON THIS MACHINE. A bank generated
anywhere else reproduces itself forever but NOT the repo's banks.

This script pins the BLAS thread count to 8 BEFORE numpy is imported (the env vars are read
at load time — setting them afterwards is a no-op), builds the real Phi_gwb from the frozen
b1 problem, factors it, and banks L with a content fingerprint.

IT DOES NOT GATE ITSELF. The bank is only canonical if ANCHOR's g1 replay (80 banked
ig_nullN_* T=15 realisations) comes back BIT-IDENTICAL through it. That is `anchor.py gate`,
run with this file in place. Bit-identical or STOP.

Run (inside the harbor container, CPU):
    python CW_transition/make_lgwb_bank.py
"""
import os

# ---- PIN THE BLAS THREAD COUNT BEFORE numpy IS IMPORTED. This is the whole point. ----
CPUS = "8"                                  # the banked convention: --cpus-per-task=8
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = CPUS
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import sys
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
sys.path.insert(0, f"{HSYMT}/CW_transition")
sys.path.insert(0, f"{HSYMT}/CW_lnL_check")
sys.path.insert(0, f"{HSYMT}/hpc_harbor/forge")
sys.path.insert(0, f"{HSYMT}/hpc_harbor/ignite")

import trackB_b1_core as C
import ignite as IG

BANK = os.path.join(HSYMT, "CW_transition", "b1_L_gwb.npz")


def main():
    print("=" * 92)
    print("CANONICAL b1_L_gwb — the frozen GWB square root")
    print("=" * 92)
    print(f"  BLAS threads pinned to {CPUS} (the banked convention) before numpy import")
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        print(f"    {v} = {os.environ.get(v)}")
    print(f"  numpy {np.__version__}")
    try:
        cfg = np.show_config(mode="dicts")
        bl = cfg.get("Build Dependencies", {}).get("blas", {})
        print(f"  BLAS: {bl.get('name')} {bl.get('version', '')}")
    except Exception:
        pass

    if os.path.exists(BANK):
        z = np.load(BANK)
        print(f"\n  bank ALREADY EXISTS at {BANK} (fp={str(z['fingerprint'])}). "
              f"Refusing to overwrite.\n  Delete it deliberately if you mean to regenerate.")
        return 1

    print("\n  building the frozen b1 problem (T=15, 116 pulsars) ...", flush=True)
    P = IG.build_ignite_problem(15, verbose=False)

    # Phi_gwb exactly as NoiseDrawer derives it — same symmetrisation, same inverse.
    # (build_ignite_problem already constructed a NoiseDrawer; take ITS Phi_gwb so the
    #  matrix we factor is bit-for-bit the one the draw path uses.)
    Pinv = np.asarray(P["amo"]["internals"]["Pinv_gwb"])
    Phi = np.linalg.inv(0.5 * (Pinv + Pinv.T))
    Phi = 0.5 * (Phi + Phi.T)
    assert np.array_equal(Phi, P["nd"].Phi_gwb), "Phi_gwb differs from the drawer's own"
    print(f"  Phi_gwb: {Phi.shape[0]}x{Phi.shape[1]}  (identical to NoiseDrawer's)")

    w = np.linalg.eigvalsh(Phi)
    print(f"  spectrum: min {w.min():.3e}  max {w.max():.3e}  cond {w.max()/max(w.min(),1e-300):.3e}")
    ndeg = int(np.sum(np.abs(np.diff(np.sort(w))) < 1e-12 * max(abs(w.max()), 1e-300)))
    print(f"  near-degenerate adjacent eigenvalue pairs: {ndeg}  <- the hazard's cause")

    print("\n  factoring: w, V = eigh(Phi);  L = V * sqrt(clip(w, 0, None))", flush=True)
    w, V = np.linalg.eigh(Phi)
    w = np.clip(w, 0.0, None)
    L = V * np.sqrt(w)

    recon = float(np.max(np.abs(L @ L.T - Phi)))
    scale = float(np.max(np.abs(Phi)))
    print(f"  reconstruction: max|L L^T - Phi| = {recon:.3e}  (max|Phi| = {scale:.3e}, "
          f"rel {recon/scale:.3e})")

    fp = C.lgwb_fingerprint(L)

    # The drawer built during the problem build took the legacy eigh path at these same 8
    # threads. If our L is not bit-identical to its L, the two eigh calls disagreed and the
    # basis is not even stable within this process — bank nothing.
    fp_drawer = C.lgwb_fingerprint(P["nd"].L_gwb)
    print(f"\n  in-process cross-check (this eigh vs the drawer's legacy eigh, both 8 threads):")
    print(f"    ours    fp {fp}")
    print(f"    drawer  fp {fp_drawer}")
    if fp != fp_drawer:
        print("    *** STOP: the two eigh calls disagree AT THE SAME THREAD COUNT. The basis is "
              "not\n        reproducible even within one process; banking it would freeze an "
              "arbitrary rotation.")
        return 1
    print("    identical — the 8-thread basis is stable in-process.")
    np.savez(BANK, L_gwb=L, fingerprint=fp, cpus=int(CPUS))
    print(f"\n  BANKED -> {BANK}")
    print(f"    fingerprint {fp}")
    print(f"    shape {L.shape}   cpus {CPUS}")
    print("\n  NOT YET CANONICAL. Gate it:  anchor.py gate  (80 banked ig_nullN_* T=15)")
    print("  Bit-identical or STOP.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
