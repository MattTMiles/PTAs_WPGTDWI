#!/usr/bin/env python3
"""SIEVE V3 -- THE E0 RE-CUT OF EVERY GLACIER-LINEAGE BANKED VERDICT
(agent SIEVE-A, 2026-07-29 addendum). REPORT-ONLY on the banks: this script writes ONLY
into SIEVE_results/. It does not rewrite a single campaign npz. The code correction
itself lives in glacier_loop._null_offenders and is forward-acting.

WHAT V3 IS FOR.  The defect is stated in full in the patched
`glacier_loop._null_offenders` docstring. Summary: the floor was cut on
(dlnL - lnK)_+ and applied to dlnL. This script re-cuts every banked GLACIER-lineage
scoreboard verdict and prints the before/after table.

WHAT CAN AND CANNOT BE RE-CUT EXACTLY -- THE HONEST LIMIT, STATED BEFORE THE TABLE.
  The correction changes the FLOOR. Recomputing the corrected floor needs the per-pulsar
  columns of each cell's null draws, and `_null_offenders` banks only the resulting
  scalar (`floor`, `floor_err`, `floor_est`, `zero_fraction`). GLACIER's null draws are
  therefore NOT recoverable from the banks, so the bit-exact canonical re-cut of banked
  verdicts REQUIRES RE-RUNNING THE NULL DRAWS ON A GPU. That is pre-registered below,
  not performed here.
  What IS exactly computable from the banks is the SCALE-MATCHED re-cut: apply the
  banked floor on the scale it was actually measured on,
        cert_E0 = (o > max(TRIALS_NAT, floor)) & (q > QBAR),  o = (dlnL - lnK)_+
  which needs no new information at all. SIEVE V1 measured how good a stand-in that is
  for the true canonical cut, on 38 cells where BOTH are computable: the two agree on
  61329 of 61480 pulsar-rows (99.75%), and disagree in BOTH directions -- so it is a
  close stand-in, NOT an identity, and this script reports it as such.
  A Delta-calibrated canonical band is printed beside it, using V1's measured
  Delta = floor_canonical - floor_glacier (median +6.97, range [+1.96, +15.37] nat over
  38 cells), so the reader sees the sensitivity rather than a single point claim.

GATES
  G-V3a  the re-derived criterion-v2.2 mask must reproduce each bank's own `n_cert`
         column EXACTLY. Any failure stops the run: if the incumbent cannot be
         reproduced, nothing said about the replacement means anything.
  G-V3b  off_glacier <= off_canonical must hold on every GLACIER null bank (it is an
         identity, see the patched docstring). Measured here on GLACIER's OWN venue so
         that the Delta imported from V1 is checked against native data, not assumed.

Output bank: SIEVE_results/sieve_v3_recut.npz
Usage: python3 hpc_harbor/sieve/sv_v3_recut.py [--repo <path>]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

QBAR = 0.9
TRIALS_NAT = 0.578

# V1's measured Delta = floor_canonical - floor_glacier, 38 cells, 3 campaigns.
DELTA_MED, DELTA_LO, DELTA_HI = 6.972, 1.958, 15.370

STEM_RE = re.compile(r"^(gl\d[^.]*?)_i(\d+)__(.+)\.npz$")


def cert_v22(dlnl, lnK, q, fl):
    return (dlnl > np.maximum(lnK + TRIALS_NAT, fl)) & (q > QBAR)


def cert_e0(dlnl, lnK, q, fl):
    o = np.maximum(dlnl - np.maximum(lnK, 0.0), 0.0)
    return (o > max(TRIALS_NAT, fl)) & (q > QBAR)


def off_pair(dlnl, lnK, q):
    fin = np.isfinite(dlnl)
    m = (dlnl > lnK) & (q > QBAR) & fin
    can = float(dlnl[m].max()) if m.any() else 0.0
    o = np.where((q > QBAR) & fin, np.maximum(dlnl - np.maximum(lnK, 0.0), 0.0), 0.0)
    return can, (float(np.max(o)) if o.size else 0.0)


def scan(paths, label, rows, events, nulls_pairs, gate_fail):
    for f in sorted(paths):
        m = STEM_RE.match(os.path.basename(f))
        if not m:
            continue
        stem, it, lane = m.group(1), int(m.group(2)), m.group(3)
        d = np.load(f, allow_pickle=True)
        need = ("dlnL_det", "lnK", "qmax", "on_true", "floor")
        if any(k not in d.files for k in need):
            continue
        dl = np.asarray(d["dlnL_det"], float)
        kk = np.asarray(d["lnK"], float)
        qq = np.asarray(d["qmax"], float)
        ot = np.asarray(d["on_true"], bool)
        fl = float(d["floor"])

        est = str(d["floor_est"]) if "floor_est" in d.files else ""
        # DEGENERATE FLOORS TAKE Delta = 0, AND THIS IS FORCED, NOT A CHOICE.
        # The two statistics vanish on exactly the same null draws (lnK >= 0, so the
        # pulsar attaining the glacier max is in the canonical mask). So an all-zero
        # glacier offender sample implies an all-zero CANONICAL sample, hence
        # floor_canonical = floor_glacier = 0 and Delta = 0. Adding the pooled Delta to
        # a degenerate floor would manufacture a bar the canonical statistic never
        # produces -- and it is not a harmless conservatism: it is what killed both
        # banked true certifications on the first pass of this script.
        degenerate = (fl <= 0.0) or ("degenerate" in est)
        dmed, dlo, dhi = ((0.0, 0.0, 0.0) if degenerate
                          else (DELTA_MED, DELTA_LO, DELTA_HI))

        c22 = cert_v22(dl, kk, qq, fl)
        ce0 = cert_e0(dl, kk, qq, fl)
        cmed = cert_v22(dl, kk, qq, max(fl, 0.0) + dmed)
        clo = cert_v22(dl, kk, qq, max(fl, 0.0) + dlo)
        chi = cert_v22(dl, kk, qq, max(fl, 0.0) + dhi)

        # ---- G-V3a --------------------------------------------------------
        if "n_cert" in d.files:
            nb = int(d["n_cert"])
            if int(c22.sum()) != nb:
                gate_fail.append((f, int(c22.sum()), nb))

        is_null = "null" in stem
        rows.append(dict(label=label, stem=stem, it=it, lane=lane, floor=fl,
                         zf=float(d["zero_fraction"]) if "zero_fraction" in d.files
                         else np.nan,
                         est=est, degenerate=degenerate, is_null=is_null,
                         n22=int(c22.sum()), t22=int((c22 & ot).sum()),
                         ne0=int(ce0.sum()), te0=int((ce0 & ot).sum()),
                         nmed=int(cmed.sum()), tmed=int((cmed & ot).sum()),
                         nlo=int(clo.sum()), tlo=int((clo & ot).sum()),
                         nhi=int(chi.sum()), thi=int((chi & ot).sum())))
        for a in np.flatnonzero(c22):
            events.append((label, stem, it, int(a), bool(ot[a]), bool(ce0[a]),
                           bool(cmed[a]), float(dl[a]), float(kk[a]), float(qq[a]), fl))
        if is_null:
            nulls_pairs.append(off_pair(dl, kk, qq))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/taylor_group/matt_miles/PTAs_WPGTDWI")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    R = a.repo
    out = a.out or f"{R}/SIEVE_results"
    os.makedirs(out, exist_ok=True)

    rows, events, nulls_pairs, gate_fail = [], [], [], []
    scan(glob.glob(f"{R}/GLACIER_results/gl*.npz"), "GLACIER",
         rows, events, nulls_pairs, gate_fail)
    scan(glob.glob(f"{R}/FROZEN_results/gl*.npz"), "FROZEN",
         rows, events, nulls_pairs, gate_fail)

    if not rows:
        raise SystemExit("no GLACIER-lineage scoreboard banks found -- STOP")

    # ---- G-V3a --------------------------------------------------------------
    print("=" * 78)
    if gate_fail:
        for f, got, want in gate_fail[:10]:
            print(f"  G-V3a FAIL {os.path.basename(f)}: recomputed {got} vs banked "
                  f"{want}")
        raise SystemExit(f"G-V3a FAIL on {len(gate_fail)}/{len(rows)} banks -- the "
                         f"incumbent criterion is not reproduced, so nothing said about "
                         f"the replacement is meaningful. STOP.")
    print(f"  G-V3a PASS: criterion-v2.2 reproduced on {len(rows)}/{len(rows)} banks")

    # ---- G-V3b --------------------------------------------------------------
    if nulls_pairs:
        P = np.array(nulls_pairs)
        viol = int(np.sum(P[:, 1] > P[:, 0] + 1e-9))
        d_native = P[:, 0] - P[:, 1]
        nz = d_native[P[:, 0] > 0]
        print(f"  G-V3b {'PASS' if viol == 0 else 'FAIL'}: off_gl <= off_can on "
              f"{len(P)-viol}/{len(P)} GLACIER null banks")
        if nz.size:
            print(f"    GLACIER-native offender gap on non-degenerate null banks "
                  f"(n={nz.size}): median {np.median(nz):+.3f} nat, "
                  f"range [{nz.min():+.3f}, {nz.max():+.3f}]")
        else:
            print("    every GLACIER null bank has a ZERO canonical offender -- the "
                  "native gap is not measurable from these banks (see the pre-"
                  "registration below)")
    print("=" * 78)

    sig = [r for r in rows if not r["is_null"]]
    nul = [r for r in rows if r["is_null"]]

    def T(rs, k):
        return int(sum(r[k] for r in rs))

    nfz = len([r for r in rows if r["label"] == "FROZEN"])
    ndg = len([r for r in sig if r["degenerate"]])
    print(f"\nBANKS: {len(sig)} signal cells + {len(nul)} null/scrambled cells; "
          f"{nfz} of the {len(rows)} come from the frozen arm (FROZEN_results, "
          f"read-only). {ndg}/{len(sig)} signal cells carry a DEGENERATE floor "
          f"(Delta = 0 there by the zero-atom identity).\n")
    print("BEFORE / AFTER  (certification events over ALL banked GLACIER-lineage cells)")
    print(f"{'cut':38s} {'certs':>7s} {'true':>6s} {'wrong':>6s}")
    print("-" * 62)
    for name, kn, kt in (("v2.2 AS BANKED (the defect)", "n22", "t22"),
                         ("E0 scale-matched  [exact from banks]", "ne0", "te0"),
                         (f"canonical, Delta = {DELTA_MED:+.2f} [central]",
                          "nmed", "tmed"),
                         (f"canonical, Delta = {DELTA_LO:+.2f} [most permissive]",
                          "nlo", "tlo"),
                         (f"canonical, Delta = {DELTA_HI:+.2f} [most strict]",
                          "nhi", "thi")):
        n, t = T(sig, kn), T(sig, kt)
        print(f"{name:38s} {n:7d} {t:6d} {n-t:6d}")

    print("\nNULL / SCRAMBLED cells (any certification here is a false positive by "
          "construction)")
    print(f"{'cut':38s} {'certs':>7s}")
    print("-" * 46)
    for name, kn in (("v2.2 AS BANKED", "n22"), ("E0 scale-matched", "ne0"),
                     (f"canonical, Delta = {DELTA_MED:+.2f}", "nmed")):
        print(f"{name:38s} {T(nul, kn):7d}")

    # ---- the event table ----------------------------------------------------
    ev = [e for e in events if "null" not in e[1]]
    print(f"\nCERTIFICATION EVENTS ON SIGNAL CELLS: {len(ev)}")
    print(f"{'cell':44s} {'it':>3s} {'psr':>4s} {'true':>5s} {'E0':>4s} {'can':>4s} "
          f"{'dlnL':>9s} {'lnK':>7s} {'floor':>8s}")
    for lab, stem, it, psr, ot, e0, cm, dl, kk, qq, fl in ev:
        print(f"{stem[:44]:44s} {it:3d} {psr:4d} {str(ot):>5s} {str(e0):>4s} "
              f"{str(cm):>4s} {dl:9.3f} {kk:7.3f} {fl:8.3f}")

    n_true = sum(1 for e in ev if e[4])
    print(f"\n  true certifications: {n_true}; surviving E0 scale-matched: "
          f"{sum(1 for e in ev if e[4] and e[5])}; surviving canonical(Delta_med): "
          f"{sum(1 for e in ev if e[4] and e[6])}")
    print(f"  wrong certifications: {len(ev)-n_true}; surviving E0: "
          f"{sum(1 for e in ev if not e[4] and e[5])}; surviving canonical(Delta_med): "
          f"{sum(1 for e in ev if not e[4] and e[6])}")

    print("\n" + "=" * 78)
    print("PRE-REGISTERED FOR SUMMIT, FROM ITERATION ZERO")
    print("=" * 78)
    print("""  1. glacier_loop._null_offenders now computes the CANONICAL offender. Every
     GLACIER-lineage run started after this commit -- Stage-1, the sky and array
     ladders, the frozen arm, and SUMMIT itself -- uses it from iteration 0. No
     mid-run estimator switch is permitted: a cell is scored under one statistic for
     its whole trajectory, or it is re-run.
  2. The banked verdicts above are re-cut IN REPORT ONLY. The campaign npz files are
     unchanged; SIEVE_results/sieve_v3_recut.npz carries both columns so any later
     analysis can choose, and must say which it chose.
  3. A BIT-EXACT canonical re-cut of the banked verdicts needs the null draws re-run
     (the per-pulsar null columns were never banked -- only the scalar floor). Cost:
     N_NULL_FULL draws x sb.columns per (cell, iteration) on the cell's own venue,
     GPU. Until that is paid, the scale-matched column is the defensible one and its
     99.75% row-agreement with the canonical cut (SIEVE V1) is its error bar.
  4. D2 (R1 2F_coh >= 15.132, R2 Delta-2F > 0) is untouched: it never reads the floor.
     Its manufactured-set kills stand exactly as banked.""")

    keys = list(rows[0].keys())
    bank = {k: np.array([r[k] for r in rows]) for k in keys}
    path = f"{out}/sieve_v3_recut.npz"
    np.savez(path, **bank,
             ev_stem=np.array([e[1] for e in ev]),
             ev_it=np.array([e[2] for e in ev]),
             ev_psr=np.array([e[3] for e in ev]),
             ev_on_true=np.array([e[4] for e in ev]),
             ev_keep_e0=np.array([e[5] for e in ev]),
             ev_keep_can=np.array([e[6] for e in ev]),
             ev_dlnL=np.array([e[7] for e in ev]),
             ev_lnK=np.array([e[8] for e in ev]),
             ev_qmax=np.array([e[9] for e in ev]),
             ev_floor=np.array([e[10] for e in ev]),
             delta_med=DELTA_MED, delta_lo=DELTA_LO, delta_hi=DELTA_HI,
             qbar=QBAR, trials_nat=TRIALS_NAT)
    print(f"\n[V3] banked {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
