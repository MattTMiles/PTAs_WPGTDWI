#!/usr/bin/env python3
"""SIEVE G-V3c -- verification of the E0 patch to `glacier_loop._null_offenders`.

WHY THIS EXISTS. The patch changes a CERTIFICATION BAR, and it had been compile-checked
but never executed: `_null_offenders` only runs inside a GLACIER cell, on GPU. A bar
change that has never run is not a verified change. This exercises the exact body now in
`glacier_loop._null_offenders` against the canonical `recut_surface.offender` on random
draws, and re-checks on synthetic data the two identities that V1 and V3 lean on.

Pure numpy, seconds, no venue. Run it after any edit to that function.
Usage: python3 hpc_harbor/sieve/sv_v3_gate.py [--trials 20000]
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, "/data/taylor_group/matt_miles/PTAs_WPGTDWI/CW_transition")
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
from recut_surface import offender as canonical      # noqa: E402

QBAR = 0.9


def patched(dlnl, lnK, q_of):
    """The exact body now in glacier_loop._null_offenders. Kept literally parallel --
    if that function changes, this must change with it or the gate is theatre."""
    fin = np.isfinite(dlnl)
    m = (dlnl > lnK) & (q_of > QBAR) & fin
    return float(dlnl[m].max()) if m.any() else 0.0


def superseded(dlnl, lnK, q_of):
    """glacier_loop._null_offenders:594 as it stood BEFORE the E0 correction."""
    o = np.where(q_of > QBAR, np.maximum(dlnl - np.maximum(lnK, 0.0), 0.0), 0.0)
    return float(np.max(o))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=20000)
    ap.add_argument("--npsr", type=int, default=116)
    a = ap.parse_args()
    rng = np.random.default_rng(20260729)

    bad = zero_bad = order_bad = n_inf = n_fin = inf_poison = 0
    for t in range(a.trials):
        dl = rng.normal(5, 20, a.npsr)
        kk = np.abs(rng.normal(6, 3, a.npsr))       # lnK = log(max(K,1)) >= 0 always
        q = rng.uniform(0, 1, a.npsr)
        has_inf = (t % 7 == 0)                       # a K = 1 pulsar: infinite peak gap
        if has_inf:
            dl[rng.integers(0, a.npsr)] = np.inf
            n_inf += 1
        p, c, s = patched(dl, kk, q), canonical(dl, kk, q), superseded(dl, kk, q)

        # (1) equality with the canonical function, EXCEPT for the declared non-finite
        #     guard -- which is spark3.offender:1063's, and rides along with the patch.
        if np.isfinite(c):
            bad += int(p != c)
        else:
            bad += int(not np.isfinite(p))

        # (2)+(3) THE TWO IDENTITIES ARE CLAIMS ABOUT THE *SCALE*, so they are tested on
        # FINITE draws, which is the domain on which both statistics are defined. On a
        # draw carrying an inf the superseded form returns inf and the patched form
        # excludes it -- that is not a violated identity, it is a SEPARATE defect of the
        # superseded code, counted below rather than folded in here. (V1 and V3 measure
        # both statistics with the same finite guard, so their numbers are unaffected.)
        if not has_inf:
            n_fin += 1
            zero_bad += int((p == 0.0) != (s == 0.0))   # the zero-atom identity
            order_bad += int(s > p + 1e-12)             # forced by lnK >= 0
        else:
            inf_poison += int(not np.isfinite(s))

    ok = (bad == 0 and zero_bad == 0 and order_bad == 0)
    print(f"trials={a.trials} npsr={a.npsr}; {n_inf} carry a K=1 inf, {n_fin} finite")
    print(f"  patched == recut_surface.offender          : "
          f"{a.trials-bad}/{a.trials}  {'PASS' if bad==0 else 'FAIL'}")
    print(f"  zero atom identical to superseded (finite) : "
          f"{n_fin-zero_bad}/{n_fin}  {'PASS' if zero_bad==0 else 'FAIL'}")
    print(f"  off_glacier <= off_canonical (finite)      : "
          f"{n_fin-order_bad}/{n_fin}  {'PASS' if order_bad==0 else 'FAIL'}")
    print(f"  [separate defect of the SUPERSEDED code] it returned a non-finite offender "
          f"on {inf_poison}/{n_inf} inf-carrying draws -- one such draw poisons the "
          f"Gumbel fit. The correction closes this too.")
    print("G-V3c " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
