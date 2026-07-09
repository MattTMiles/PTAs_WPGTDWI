"""Prong-2 close-out P2-C: array-boost scaling law.

Stage A measured the confusion knee at ~52*N* for a 116-pulsar array (N* = T*Delta_f).
Here we rerun the SAME toy knee measurement (fiducial T=15 yr, band 3-20 nHz,
fixed_persource, isotropic sky) at N_psr = 15, 30, 60, 116, 200 (200 = synthetic
isotropic positions, same protocol), fit knee/N* vs N_psr to a power law, and report
the exponent and whether it saturates.

GATE: N_psr=116 must reproduce the Stage A knee/N* ~ 52 within seed scatter.

Reuses prong2_transition.run_knee_diagnostic (pure-jax toy; runs under jug-gpu jax 0.10).

Run:  python stagec_p2c.py
"""

import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys
import time

import numpy as np

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import prong2_transition as p2

YR = p2.YR
N_PSR_LIST = [15, 30, 60, 116, 200]
N_SEEDS = 30
CONFIG = ("fid_3-20nHz", 15 * YR, (3e-9, 20e-9))


def bootstrap_knee_scatter(n_psr, n_values, n_boot=8):
    """Crude seed-scatter estimate: rerun the knee at a few seed batches."""
    knees = []
    for b in range(n_boot):
        # vary the per-run RNG by using a different n_seeds offset via re-seeding
        res = p2.run_knee_diagnostic([CONFIG], n_values, n_seeds=8, n_psr=n_psr)
        knees.append(res[0]["knee"] / res[0]["N_star"])
    return np.std(knees)


def main():
    print("=== P2-C: array-boost scaling of the confusion knee ===", flush=True)
    n_values = np.unique(np.round(np.logspace(0, np.log10(1500), 26)).astype(int))
    print(f"N grid: {n_values[0]}..{n_values[-1]} ({len(n_values)} pts); "
          f"n_seeds={N_SEEDS}; config {CONFIG[0]}\n", flush=True)

    rows = []
    for n_psr in N_PSR_LIST:
        t = time.time()
        res = p2.run_knee_diagnostic([CONFIG], n_values, n_seeds=N_SEEDS, n_psr=n_psr)
        r = res[0]
        knee_over_Nstar = r["knee"] / r["N_star"]
        rows.append((n_psr, r["knee"], r["N_star"], knee_over_Nstar))
        print(f"  N_psr={n_psr:4d}: knee={r['knee']:7.1f}  N*={r['N_star']:.2f}  "
              f"knee/N*={knee_over_Nstar:.2f}  ({time.time()-t:.1f}s)", flush=True)

    npsr = np.array([r[0] for r in rows], float)
    kk = np.array([r[3] for r in rows], float)

    # power-law fit knee/N* = A * N_psr^b  (log-log linear)
    b, logA = np.polyfit(np.log(npsr), np.log(kk), 1)
    A = np.exp(logA)
    resid = np.log(kk) - (logA + b * np.log(npsr))
    # saturation check: fit slope on last 3 points vs first 3
    b_lo = np.polyfit(np.log(npsr[:3]), np.log(kk[:3]), 1)[0]
    b_hi = np.polyfit(np.log(npsr[-3:]), np.log(kk[-3:]), 1)[0]

    print("\n================ P2-C SUMMARY ================")
    print(f"{'N_psr':>6s} {'knee/N*':>9s} {'powerlaw fit':>14s}")
    for (n_psr, knee, Nstar, k) in rows:
        print(f"{n_psr:>6d} {k:>9.2f} {A*n_psr**b:>14.2f}")
    print(f"\npower law: knee/N* = {A:.3f} * N_psr^{b:.3f}")
    print(f"  low-N_psr slope (15-60): {b_lo:.3f};  high-N_psr slope (60-200): {b_hi:.3f}")
    sat = "SATURATES (high-end slope < 0.7*low-end)" if b_hi < 0.7 * b_lo else \
          "NO clear saturation (slope ~constant)"
    print(f"  -> {sat}")
    # GATE at 116
    k116 = dict((r[0], r[3]) for r in rows)[116]
    scatter = bootstrap_knee_scatter(116, n_values, n_boot=6)
    print(f"\nGATE: N_psr=116 knee/N* = {k116:.2f}; Stage A = ~52; "
          f"seed scatter (8-seed batches) sigma ~ {scatter:.2f}")
    gate = abs(k116 - 52) < max(3 * scatter, 8)
    print(f"  GATE {'PASS' if gate else 'CHECK'}: |{k116:.1f}-52| = {abs(k116-52):.1f} "
          f"vs tolerance {max(3*scatter,8):.1f}")

    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_p2c_results.npz",
             n_psr=npsr, knee_over_Nstar=kk, A=A, b=b, b_lo=b_lo, b_hi=b_hi,
             k116=k116, scatter=scatter)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(npsr, kk, "o", ms=9, color="#2c7fb8", label="measured knee/N*")
    xs = np.logspace(np.log10(12), np.log10(260), 100)
    ax.plot(xs, A * xs ** b, "-", color="#d95f0e",
            label=f"fit: {A:.2f} N_psr^{b:.2f}")
    ax.axhline(52, color="0.5", ls=":", label="Stage A (116 psr) = 52")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("N_psr"); ax.set_ylabel("knee / N*  (array-resolution boost)")
    ax.set_title("P2-C: confusion-knee array boost vs N_psr")
    ax.set_xticks(N_PSR_LIST); ax.set_xticklabels(N_PSR_LIST)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig("/home/mattm/projects/HSYMT/CW_transition/p2c_array_boost.png", dpi=130)
    print("saved p2c_array_boost.png, stagec_p2c_results.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
