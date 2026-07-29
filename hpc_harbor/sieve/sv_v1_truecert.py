#!/usr/bin/env python3
"""SIEVE V1 -- THE E0 CANONICAL-STATISTIC RECUT, TESTED WHERE THE LOOP WORKS
(agent SIEVE-A, 2026-07-29 addendum). REPORT-ONLY: no bank is written outside
SIEVE_results/, no campaign verdict is moved, no likelihood is re-evaluated.

WHY.  T7's E0 finding was measured on the GLACIER banks, which hold n = 2 true
certifications. A tightening that "keeps 2/2 true certs" is not evidence that it keeps
true certs -- with n = 2 the 95% upper bound on the kill rate is 78% (rule of three).
This script re-asks the question on the true-cert-rich banks: CHORUS soft loops,
IGNITE-2 soft loops, and the GENERALISE arm A-SKY units.

THE DEFECT, RESTATED PRECISELY.
  canonical (recut_surface.offender:75, and chorus.py:605, ignite2.py:166,
             surface.py:334, spark3.py:1063, anchor.py:322, kindle_d7_fall.py:67,
             bank_surface_offenders.py:51 -- EVERY campaign except one):
      off_can = max{ dlnL[a] : dlnL[a] > lnK[a] and q[a] > 0.9 }, else 0
  glacier   (glacier_loop._null_offenders:594):
      off_gl  = max_a  ( q[a] > 0.9 ? max(dlnL[a] - max(lnK[a],0), 0) : 0 )
  The floor is the (1-alpha) quantile of the offender over no-CW draws. GLACIER cuts its
  floor on off_gl and then compares dlnL against it -- a bar measured on one scale
  applied on another, in the PERMISSIVE direction.

AN EXACT FACT, WORTH STATING BECAUSE IT BOUNDS THE WHOLE ARGUMENT.
  lnK = log(max(K_counted, 1)) >= 0 always, so max(lnK, 0) == lnK. The two statistics
  are therefore ZERO on exactly the same draws -- both vanish iff no pulsar has
  (dlnL > lnK) and (q > 0.9). The zero atom, and hence `zero_fraction`, is IDENTICAL
  under both. They differ only in the VALUE of the non-zero draws, and there by exactly
  lnK of whichever pulsar attains the max. So the gap is a pure lnK offset, not a
  different event set, and Delta = floor_can - floor_gl is measured cell by cell below.

THE THREE CUTS SCORED ON EVERY REALISATION (QBAR = 0.9, TRIALS_NAT = 0.578)
  A  AS-BANKED / CANONICAL      cert = (dlnL > max(lnK + 0.578, floor_can)) & (q > QBAR)
       what CHORUS / IGNITE-2 / GENERALISE actually did, and what patching
       glacier_loop._null_offenders makes GLACIER do.
  B  THE DEFECT                 cert = (dlnL > max(lnK + 0.578, floor_gl )) & (q > QBAR)
       what GLACIER does today, reconstructed on these banks so the damage is measured
       on independent, true-cert-rich substrate rather than extrapolated from n = 2.
  C  E0 SCALE-MATCHED           cert = (o    > max(0.578,       floor_gl )) & (q > QBAR)
       T7's recut: the banked glacier-scale floor re-applied on the scale it was
       measured on, where o[a] = max(dlnL[a] - lnK[a], 0).

  A vs C IS THE LOAD-BEARING COMPARISON FOR V3. If they agree, then T7's scale-matched
  recut -- which needs NO new null draws -- is a valid stand-in for the canonical fix on
  banks whose canonical floor cannot be recomputed without re-running the nulls on a GPU.
  If they disagree, the GLACIER re-cut needs fresh nulls and V3 must say so.

GATES (all bit-exact; a failure means the null set or the estimator is wrong, and every
number downstream would be meaningless)
  G-V1a  the canonical offender array recomputed from each cell's own null banks must
         reproduce the campaign's BANKED offender array element-wise (gen_ledger
         `offenders`, ch_floors `offN_i`, ig2_floors `offN_i`).
  G-V1b  the floor re-derived from that array with the campaign's own estimator must
         reproduce the campaign's banked adopted floor.
  G-V1c  cut A recomputed must reproduce the campaign's banked certification count where
         one is banked.

Output bank: SIEVE_results/sieve_v1_truecert.npz
Usage: python3 hpc_harbor/sieve/sv_v1_truecert.py [--repo <path>]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

QBAR = 0.9
TRIALS_NAT = 0.578
ALPHA = 0.05
Z_ALPHA = -np.log(-np.log(1.0 - ALPHA))
C_SD = 2.80


# ============================================================
# the two statistics
# ============================================================
def off_canonical(dlnL, lnK, qmax):
    """recut_surface.offender:75 verbatim, with spark3's non-finite exclusion (a pulsar
    with a single fringe peak has K = 1 and an infinite gap; an inf would silently
    poison the extreme-value fit)."""
    fin = np.isfinite(dlnL)
    m = (dlnL > lnK) & (qmax > QBAR) & fin
    return (float(dlnL[m].max()) if m.any() else 0.0), int((~fin).sum())


def off_glacier(dlnL, lnK, qmax):
    """glacier_loop._null_offenders:594 verbatim (same non-finite exclusion applied, so
    the two statistics differ ONLY in the scale and not in the hygiene)."""
    fin = np.isfinite(dlnL)
    o = np.where((qmax > QBAR) & fin, np.maximum(dlnL - np.maximum(lnK, 0.0), 0.0), 0.0)
    return float(np.max(o)) if o.size else 0.0


def gumbel_floor(x):
    from scipy.stats import gumbel_r
    x = np.asarray(x, float)
    mu, beta = gumbel_r.fit(x)
    return (float(mu + beta * Z_ALPHA), float(mu), float(beta),
            float(C_SD * beta / np.sqrt(len(x))))


def emp_q95(x):
    return float(np.quantile(np.asarray(x, float), 1.0 - ALPHA))


def wilson(k, n, z=1.959963984540054):
    """Wilson score interval -- Wald collapses to [0,0] at k = 0, and k = 0 is exactly
    the outcome a survival test is built to detect."""
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(p), float(max(0.0, c - h)), float(min(1.0, c + h))


# ============================================================
# the three cuts
# ============================================================
def cut_A(dlnL, lnK, q, floor_can):
    return (dlnL > np.maximum(lnK + TRIALS_NAT, floor_can)) & (q > QBAR)


def cut_B(dlnL, lnK, q, floor_gl):
    return (dlnL > np.maximum(lnK + TRIALS_NAT, floor_gl)) & (q > QBAR)


def cut_C(dlnL, lnK, q, floor_gl):
    o = np.maximum(dlnL - np.maximum(lnK, 0.0), 0.0)
    return (o > max(TRIALS_NAT, floor_gl)) & (q > QBAR)


# ============================================================
# cell assembly
# ============================================================
def cell_key(path, sig_prefix, null_prefix):
    """Strip the campaign prefix and the per-realisation `_g<geo>_n<seed>` suffix. The
    remainder is the cell tag, and the SAME tag identifies that cell's null banks."""
    b = os.path.basename(path)
    if not b.startswith(sig_prefix):
        raise ValueError(f"{b} does not start with {sig_prefix!r}")
    core = b[len(sig_prefix):]
    core = re.sub(r"_g-?\d+_n\d+\.npz$", "", core)
    return core


def load_nulls(pattern):
    out = []
    for f in sorted(glob.glob(pattern)):
        d = np.load(f, allow_pickle=True)
        if "dlnL_det" not in d.files:
            continue
        out.append((f, np.asarray(d["dlnL_det"], float), np.asarray(d["lnK"], float),
                    np.asarray(d["qmax"], float)))
    return out


def score_cell(name, nulls, sigs, banked_off=None, banked_floor=None,
               banked_ncert=None, verbose=True):
    """Returns a dict of per-cell results, or None if the cell has no usable nulls."""
    if len(nulls) < 20:
        if verbose:
            print(f"  {name:34s} SKIP (only {len(nulls)} null draws; a q95 is not "
                  f"readable)")
        return None

    can, gl, n_inf = [], [], 0
    for _, d, k, q in nulls:
        c, ni = off_canonical(d, k, q)
        can.append(c); n_inf += ni
        gl.append(off_glacier(d, k, q))
    can = np.array(can); gl = np.array(gl)

    # ---- G-V1d: the pointwise inequality is an IDENTITY, not an expectation ----
    # lnK = log(max(K,1)) >= 0, so for the pulsar attaining off_gl we have
    #   dlnL - lnK <= dlnL  and that pulsar satisfies dlnL > lnK, hence is in the
    # canonical mask. Therefore off_gl <= off_can on EVERY draw. A violation means the
    # two statistics were evaluated on different samples -- which is a bug in this
    # script, not a finding, and it must stop the cell rather than print a number.
    bad = int(np.sum(gl > can + 1e-9))
    if bad:
        raise RuntimeError(
            f"{name}: G-V1d FAIL -- off_gl > off_can on {bad}/{len(can)} draws. The two "
            f"statistics are not on the same sample. Refusing to report Delta.")

    # ---- G-V1a -----------------------------------------------------------
    g1a = "n/a"
    if banked_off is not None:
        b = np.asarray(banked_off, float)
        if b.shape == can.shape:
            g1a = "PASS" if np.array_equal(np.sort(b), np.sort(can)) else "FAIL"
        else:
            g1a = f"shape {b.shape} vs {can.shape}"

    zf_can = float(np.mean(can == 0.0))
    zf_gl = float(np.mean(gl == 0.0))

    f_can_g, _, _, sd_can = gumbel_floor(can)
    f_gl_g, _, _, sd_gl = gumbel_floor(gl)
    f_can_e, f_gl_e = emp_q95(can), emp_q95(gl)

    # ---- G-V1b -----------------------------------------------------------
    # BOTH floors come from the SAME loaded sample with the SAME estimator. The banked
    # floor is used ONLY as a gate -- mixing a banked canonical floor with a locally
    # computed glacier floor is exactly how a negative Delta (an impossibility, see
    # G-V1d) gets printed as if it were a measurement.
    g1b = "n/a"
    est = "gumbel"
    if banked_floor is not None:
        cand = {"gumbel": f_can_g, "emp": f_can_e}
        est = min(cand, key=lambda k: abs(cand[k] - banked_floor))
        g1b = (f"PASS[{est}]" if abs(cand[est] - banked_floor) < 1e-6
               else f"FAIL(banked {banked_floor:.4f} vs gumbel {f_can_g:.4f} / "
                    f"emp {f_can_e:.4f})")
    floor_can = f_can_g if est == "gumbel" else f_can_e
    floor_gl = f_gl_g if est == "gumbel" else f_gl_e

    # ---- score the signal / sloop realisations ---------------------------
    nA = nB = nC = 0
    tA = tB = tC = 0
    agreeAC = disagreeAC = 0
    ncert_banked_tot = 0
    for f, d, k, q, ot in sigs:
        a = cut_A(d, k, q, floor_can)
        b = cut_B(d, k, q, floor_gl)
        c = cut_C(d, k, q, floor_gl)
        nA += int(a.sum()); nB += int(b.sum()); nC += int(c.sum())
        tA += int((a & ot).sum()); tB += int((b & ot).sum()); tC += int((c & ot).sum())
        agreeAC += int((a == c).sum()); disagreeAC += int((a != c).sum())

    res = dict(name=name, n_null=len(nulls), n_sig=len(sigs), n_inf=n_inf,
               zf_can=zf_can, zf_gl=zf_gl,
               floor_can=floor_can, floor_gl=floor_gl,
               delta=floor_can - floor_gl,
               f_can_g=f_can_g, f_gl_g=f_gl_g, f_can_e=f_can_e, f_gl_e=f_gl_e,
               sd_can=sd_can, sd_gl=sd_gl,
               nA=nA, nB=nB, nC=nC, tA=tA, tB=tB, tC=tC,
               agreeAC=agreeAC, disagreeAC=disagreeAC, g1a=g1a, g1b=g1b)
    if verbose:
        print(f"  {name:34s} nulls={len(nulls):4d} zf={zf_can:.2f} "
              f"F_can={floor_can:8.3f} F_gl={f_gl_g:8.3f} D={floor_can-floor_gl:+8.3f} | "
              f"cert A/B/C {nA:4d}/{nB:4d}/{nC:4d}  true {tA:3d}/{tB:3d}/{tC:3d} | "
              f"G-V1a {g1a}")
    return res


def load_sigs(paths, dl, lk, qq, ot):
    out = []
    for f in sorted(paths):
        d = np.load(f, allow_pickle=True)
        if dl not in d.files:
            continue
        out.append((f, np.asarray(d[dl], float), np.asarray(d[lk], float),
                    np.asarray(d[qq], float), np.asarray(d[ot], bool)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/taylor_group/matt_miles/PTAs_WPGTDWI")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    R = a.repo
    out = a.out or f"{R}/SIEVE_results"
    os.makedirs(out, exist_ok=True)
    rows = []

    # ---------------------------------------------------------------- A-SKY
    print("\n=== GENERALISE arm A-SKY (truth-anchored: theta_base = theta_src, "
          "generalise.py:376) ===")
    for led in sorted(glob.glob(f"{R}/GENERALISE_results/gen_ledger_AS_*.npz")):
        L = np.load(led, allow_pickle=True)
        tag = str(L["tag"])
        nulls = load_nulls(f"{R}/GENERALISE_results/gen_fnullN_{tag}_g*_n*.npz")
        sigs = load_sigs(glob.glob(f"{R}/GENERALISE_results/gen_sig_{tag}_g*_n*.npz"),
                         "dlnL_det", "lnK", "qmax", "on_true")
        r = score_cell(f"ASKY {tag}", nulls, sigs,
                       banked_off=L["offenders"],
                       banked_floor=float(L["floor_adopted"]))
        if r:
            r["family"] = "ASKY"
            r["survivor"] = ("h1275_k5" in tag)
            rows.append(r)

    # ---------------------------------------------------------------- CHORUS
    print("\n=== CHORUS soft loops (fitted templates; the true-cert reservoir) ===")
    ch_sloops = glob.glob(f"{R}/CHORUS_results/ch_sloop_*.npz")
    by_cell = {}
    for f in ch_sloops:
        by_cell.setdefault(cell_key(f, "ch_sloop_", "ch_fnullA_"), []).append(f)
    for cell, fs in sorted(by_cell.items()):
        # THE FLOOR'S OWN NULL SET IS `fnullN`, not `fnullA` (chorus.py:635 --
        # fnullA/fnullL are the EXTRA draws that make the wider fALL floor). Using
        # fnullA here is what produced the impossible negative Delta on the first pass.
        nulls = load_nulls(f"{R}/CHORUS_results/ch_fnullN_{cell}_g*_n*.npz")
        if len(nulls) < 100 and cell == "m0e00_vlbi":
            # (C4) chorus.py:626 reuses the banked IGNITE-2 fresh-floor cell verbatim
            nulls = load_nulls(f"{R}/IGNITE2_results/ig_fnullN_h1325_T30_vlbi_g*_n*.npz")
            print("    [C4] m0e00_vlbi: nulls reused from IGNITE-2 h1325_T30_vlbi "
                  "(chorus.py:626)")
        sigs = load_sigs(fs, "dlnL_final", "lnK", "qmax_final", "on_true_final")
        fl = float(np.load(fs[0], allow_pickle=True)["floor_cell"])
        r = score_cell(f"CHORUS {cell}", nulls, sigs, banked_floor=fl)
        if r:
            r["family"] = "CHORUS"; r["survivor"] = False
            rows.append(r)

    # ---------------------------------------------------------------- IGNITE-2
    print("\n=== IGNITE-2 soft loops (the seeded trajectories) ===")
    ig = glob.glob(f"{R}/IGNITE2_results/ig_sloop_*.npz")
    by_cell = {}
    for f in ig:
        by_cell.setdefault(cell_key(f, "ig_sloop_", "ig_fnullA_"), []).append(f)
    for cell, fs in sorted(by_cell.items()):
        nulls = load_nulls(f"{R}/IGNITE2_results/ig_fnullN_{cell}_g*_n*.npz")
        sigs = load_sigs(fs, "dlnL_final", "lnK", "qmax_final", "on_true_final")
        fl = float(np.load(fs[0], allow_pickle=True)["floor_cell"])
        r = score_cell(f"IGNITE2 {cell}", nulls, sigs, banked_floor=fl)
        if r:
            r["family"] = "IGNITE2"; r["survivor"] = False
            rows.append(r)

    if not rows:
        raise SystemExit("no cells scored -- STOP rather than report an empty table")

    # ---------------------------------------------------------------- totals
    def tot(key, sel=None):
        return int(sum(r[key] for r in rows if sel is None or sel(r)))

    print("\n" + "=" * 78)
    print("V1 TOTALS  (A = canonical/as-banked, B = the GLACIER defect, "
          "C = E0 scale-matched)")
    print("=" * 78)
    for label, sel in (("ALL cells", None),
                       ("A-SKY survivor (8 skies)",
                        lambda r: r["family"] == "ASKY" and r["survivor"]),
                       ("A-SKY all", lambda r: r["family"] == "ASKY"),
                       ("CHORUS sloops", lambda r: r["family"] == "CHORUS"),
                       ("IGNITE-2 sloops", lambda r: r["family"] == "IGNITE2")):
        n = len([r for r in rows if sel is None or sel(r)])
        if not n:
            continue
        print(f"  {label:26s} cells={n:3d}  "
              f"cert A/B/C = {tot('nA', sel):5d}/{tot('nB', sel):5d}/{tot('nC', sel):5d}"
              f"   TRUE cert A/B/C = {tot('tA', sel):4d}/{tot('tB', sel):4d}/"
              f"{tot('tC', sel):4d}")

    TA, TB, TC = tot("tA"), tot("tB"), tot("tC")
    NA, NB, NC = tot("nA"), tot("nB"), tot("nC")
    print()
    for label, sel in (("ALL", None),
                       ("A-SKY survivor",
                        lambda r: r["family"] == "ASKY" and r["survivor"]),
                       ("A-SKY all", lambda r: r["family"] == "ASKY"),
                       ("CHORUS", lambda r: r["family"] == "CHORUS"),
                       ("IGNITE-2", lambda r: r["family"] == "IGNITE2")):
        ta, tc = tot("tA", sel), tot("tC", sel)
        na, nc = tot("nA", sel), tot("nC", sel)
        lost = ta - tc
        p, lo, hi = wilson(max(lost, 0), max(ta, 1))
        print(f"  TRUE-CERT SURVIVAL {label:15s} {tc:5d}/{ta:5d} = {tc/max(ta,1):.4f}"
              f"   net true-cert change {tc-ta:+5d}   kill rate {p:.4f} "
              f"[{lo:.4f}, {hi:.4f}]"
              f"   wrong-cert rate A {(na-ta)/max(na,1):.4f} -> C "
              f"{(nc-tc)/max(nc,1):.4f}")
    print(f"\n  TRUE-CERT SURVIVAL under the canonical/E0 form: "
          f"{TC}/{TA} = {TC/max(TA,1):.4f}")
    print(f"  FALSE certs added by the defect (B - A): "
          f"{(NB-TB)-(NA-TA):+d}  (cert total {NB-NA:+d})")
    print(f"  A-vs-C per-pulsar agreement: {tot('agreeAC')}/"
          f"{tot('agreeAC')+tot('disagreeAC')} "
          f"({tot('disagreeAC')} disagreements)")
    d = np.array([r["delta"] for r in rows])
    print(f"  Delta = floor_can - floor_gl over {len(d)} cells: "
          f"median {np.median(d):+.3f}  mean {d.mean():+.3f}  "
          f"range [{d.min():+.3f}, {d.max():+.3f}]  negative in {int((d<0).sum())}")
    zc = np.array([r["zf_can"] for r in rows]); zg = np.array([r["zf_gl"] for r in rows])
    print(f"  zero-atom identity check (must be exactly equal): "
          f"max|zf_can - zf_gl| = {np.max(np.abs(zc-zg)):.3e}")
    print(f"  G-V1a: {sum(1 for r in rows if r['g1a']=='PASS')} PASS / "
          f"{sum(1 for r in rows if r['g1a']=='FAIL')} FAIL / "
          f"{sum(1 for r in rows if r['g1a'] not in ('PASS','FAIL'))} n/a")
    print(f"  G-V1b: {sum(1 for r in rows if str(r['g1b']).startswith('PASS'))} PASS / "
          f"{sum(1 for r in rows if str(r['g1b']).startswith('FAIL'))} FAIL")

    keys = ["name", "family", "survivor", "n_null", "n_sig", "n_inf", "zf_can", "zf_gl",
            "floor_can", "floor_gl", "delta", "f_can_g", "f_gl_g", "f_can_e", "f_gl_e",
            "sd_can", "sd_gl", "nA", "nB", "nC", "tA", "tB", "tC",
            "agreeAC", "disagreeAC", "g1a", "g1b"]
    bank = {k: np.array([r[k] for r in rows]) for k in keys}
    path = f"{out}/sieve_v1_truecert.npz"
    np.savez(path, **bank, qbar=QBAR, trials_nat=TRIALS_NAT, alpha=ALPHA)
    print(f"\n[V1] banked {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
