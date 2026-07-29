#!/usr/bin/env python3
"""LEDGER B2 -- THE PERSISTENCE RULE (criterion v2.3 companion term) + BANKED RE-CUT.

THE DEFECT UNDER AUDIT: ITERATION MULTIPLICITY. The scoreboard runs at the end of every
loop iteration (glacier_loop.run_cell step (f)) and a cell is called "certified" if ANY
scored iteration certifies. With n_iter = 6 that is SIX independent looks at a moving
joint template, and the trials term the criterion charges (`lnK` over the carried census
+ TRIALS_NAT = 0.578) prices the FRINGE and CENSUS multiplicity only -- it prices no
part of the ITERATION axis. A statistic that is marginal at the bar therefore gets six
chances to cross it, and the campaign's own record shows the crossings FLICKER (S4.24
Finding 1: the two true certs at r13p5/none/s3 appear at i0 and i4).

THE COMPANION TERM (PRE-REGISTERED HERE, FORWARD, 2026-07-29 -- it does NOT retro-change
any banked verdict):

    criterion v2.3 = criterion v2.2  AND  PERSISTENCE

    PERSISTENCE(p, i) := cert_v22(p, i) AND cert_v22(p, i-1)

A pulsar is certified at iteration i under v2.3 only if it was ALSO certified at the
immediately preceding SCORED iteration. Declared consequences, stated before scoring:

  * i = 0 can never certify under v2.3 (no predecessor). This is deliberate: iteration 0
    is the truth-anchored feed (promote is AT DRAWN TRUTH), the single most optimistic
    state the loop ever occupies, and it is exactly the state a persistence rule should
    refuse to trust on its own.
  * v2.3 is STRICTLY MORE CONSERVATIVE than v2.2: cert_v23 subset cert_v22 always. It can
    only remove certifications. It therefore CANNOT create a new wrong cert, and any
    wrong cert it removes is a genuine gain.
  * The rule is a CONJUNCTION, not a re-derivation of the bar. Every quantity it reads
    (dlnL_det, lnK, qmax, floor) is the banked raw column; no floor is re-cut, no
    estimator is switched. `lean-npz discipline` is what makes this possible at all.
  * A cell whose certifications are all at DISJOINT iterations is exactly the flicker
    class the rule exists to kill. If the campaign's only true certs are in that class,
    that is the result, and it is reported as such rather than argued around.

WHY NOT A HARSHER MULTIPLICITY CHARGE (declared, so the choice is auditable). The
alternative is an additive nat penalty ln(n_iter) = 1.79 on the bar. Rejected here for
two reasons: (a) the six looks are NOT independent -- consecutive iterations share the
data and most of the template, so ln(n_iter) over-charges by an unknown factor; (b) an
additive bar change WOULD retro-alter the meaning of every banked floor, which is
forbidden. Persistence is a conjunction on the existing bar: it prices multiplicity by
demanding the excursion be a PROPERTY OF THE STATE rather than of one look, and it costs
nothing to audit.

GATE (run first, must pass before any re-cut number is read): the re-derived v2.2 cert
mask must reproduce the banked `n_cert` column EXACTLY on every bank. If it does not,
the re-cut machinery is wrong and nothing downstream is quotable.

Output bank: reports/ledger_b2_persistence.npz. REPORT-ONLY.

Usage:  python3 hpc_harbor/ledger/ledger_b2_persistence.py [--repo <path>]
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import numpy as np

TRIALS_NAT = 0.578          # SPARK-2 SS2, as wired in glacier_loop.py:106
QBAR = 0.9                  # glacier_loop.py:104

STEM_RE = re.compile(
    r"^(?P<pre>gl1|gl2)_"
    r"(?:(?P<rung>r13p\d+)(?:_w(?P<w>[0-9p]+))?_)?"
    r"(?P<kind>cell|null\d+)_"
    r"(?P<arm>e07|none)_s(?P<sky>\d+)_T(?P<T>\d+)_(?P<tier>\w+?)_i(?P<iter>\d+)$"
)


def cert_v22(z):
    """glacier_loop.py:533 verbatim, from the banked RAW columns."""
    dlnl = np.asarray(z["dlnL_det"], float)
    lnK = np.asarray(z["lnK"], float)
    q = np.asarray(z["qmax"], float)
    fl = float(z["floor"])
    return (dlnl > np.maximum(lnK + TRIALS_NAT, fl)) & (q > QBAR)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/taylor_group/matt_miles/PTAs_WPGTDWI")
    a = ap.parse_args()
    rdir = os.path.join(a.repo, "GLACIER_results")

    # ---------------- load ----------------
    cells = defaultdict(dict)          # (stage,rung,w,T,arm,sky,kind) -> {iter: record}
    for fn in sorted(os.listdir(rdir)):
        if not fn.endswith(".npz") or "__" not in fn:
            continue
        m = STEM_RE.match(fn.split("__")[0])
        if m is None:
            continue
        z = np.load(os.path.join(rdir, fn), allow_pickle=True)
        if "dlnL_det" not in z.files or "n_cert" not in z.files:
            continue
        g = m.groupdict()
        key = (g["pre"], g["rung"] or "stage1",
               float((g["w"] or "1").replace("p", ".")), int(g["T"]),
               g["arm"], int(g["sky"]), g["kind"])
        cells[key][int(g["iter"])] = dict(
            cert=cert_v22(z), banked=int(z["n_cert"]),
            on_true=np.asarray(z["on_true"], bool),
            dlnl=np.asarray(z["dlnL_det"], float),
            q=np.asarray(z["qmax"], float),
            floor=float(z["floor"]), floor_est=str(z["floor_est"]),
            fn=fn)
    n_bank = sum(len(v) for v in cells.values())
    print(f"[LEDGER-B2] {n_bank} per-iteration banks over {len(cells)} cells\n")

    # ---------------- GATE: re-derived v2.2 == banked n_cert ----------------
    bad = [(k, i, int(r["cert"].sum()), r["banked"], r["fn"])
           for k, it in cells.items() for i, r in it.items()
           if int(r["cert"].sum()) != r["banked"]]
    print("-- GATE: re-derived criterion-v2.2 mask vs the banked n_cert column --")
    if bad:
        print(f"   FAIL: {len(bad)}/{n_bank} banks disagree. First 5:")
        for k, i, a_, b_, fn in bad[:5]:
            print(f"     {fn}: re-cut {a_} vs banked {b_}")
        print("   The re-cut machinery does not reproduce the bank. STOPPING -- no "
              "downstream number from this run is quotable.")
        return 1
    print(f"   PASS: {n_bank}/{n_bank} banks reproduce exactly "
          f"(dlnL > max(lnK+{TRIALS_NAT}, floor)) & (qmax > {QBAR})\n")

    # ---------------- the re-cut ----------------
    print("-- THE RE-CUT: every cell carrying at least one v2.2 certification --")
    print("   cell                                         iters  v2.2 certs (psr@iter)"
          "                     -> v2.3")
    recut = []
    tot22 = tot23 = w22 = w23 = 0
    for k in sorted(cells, key=str):
        its = sorted(cells[k])
        ev22, ev23 = [], []
        for i in its:
            c = cells[k][i]["cert"]
            prev = cells[k].get(i - 1)
            pers = c & (prev["cert"] if prev is not None else np.zeros_like(c))
            for p in np.where(c)[0]:
                ev22.append((int(p), i, bool(cells[k][i]["on_true"][p]),
                             float(cells[k][i]["dlnl"][p])))
            for p in np.where(pers)[0]:
                ev23.append((int(p), i, bool(cells[k][i]["on_true"][p]),
                             float(cells[k][i]["dlnl"][p])))
        if not ev22:
            continue
        tot22 += len(ev22); tot23 += len(ev23)
        w22 += sum(1 for e in ev22 if not e[2]); w23 += sum(1 for e in ev23 if not e[2])
        ks = "/".join(str(x) for x in k)
        s22 = ",".join(f"{p}@i{i}{'T' if t else 'F'}" for p, i, t, _ in ev22)
        s23 = ",".join(f"{p}@i{i}{'T' if t else 'F'}" for p, i, t, _ in ev23) or "--none--"
        print(f"   {ks:44s} {len(its):3d}  {s22[:44]:44s} -> {s23}")
        recut.append(dict(cell=ks, n_iter=len(its), ev22=ev22, ev23=ev23))
    print()
    print(f"   TOTALS  v2.2: {tot22} certification events ({w22} wrong / "
          f"{tot22 - w22} on_true)")
    print(f"           v2.3: {tot23} certification events ({w23} wrong / "
          f"{tot23 - w23} on_true)")
    print(f"           killed by persistence: {tot22 - tot23} "
          f"({w22 - w23} wrong, {(tot22 - w22) - (tot23 - w23)} true)\n")

    # ---------------- the flicker anatomy ----------------
    print("-- FLICKER ANATOMY: per (cell, pulsar), which iterations certified --")
    print("   cell / psr                                        iters certified   "
          "consecutive-run lengths   survives v2.3")
    flick = []
    for k in sorted(cells, key=str):
        its = sorted(cells[k])
        allp = set()
        for i in its:
            allp |= set(np.where(cells[k][i]["cert"])[0].tolist())
        for p in sorted(allp):
            hit = [i for i in its if cells[k][i]["cert"][p]]
            runs, cur = [], 1
            for x, y in zip(hit, hit[1:]):
                cur = cur + 1 if y == x + 1 else (runs.append(cur) or 1)
            runs.append(cur)
            surv = max(runs) >= 2
            ot = bool(cells[k][hit[0]]["on_true"][p])
            ks = "/".join(str(x) for x in k)
            print(f"   {ks}/psr{p:<4d}".ljust(53)
                  + f"{str(hit):16s} {str(runs):12s} "
                  f"{'YES' if surv else 'NO ':3s}  on_true={ot}")
            flick.append(dict(cell=ks, psr=int(p), iters=hit, runs=runs,
                              survives=bool(surv), on_true=ot))
    print()

    out = os.path.join(a.repo, "reports", "ledger_b2_persistence.npz")
    np.savez(
        out,
        note=("LEDGER B2: criterion-v2.3 = v2.2 AND persistence (cert holds at two "
              "consecutive scored iterations). PRE-REGISTERED FORWARD 2026-07-29; this "
              "re-cut of the banked histories is REPORT-ONLY and changes no banked "
              "verdict. Gate passed: re-derived v2.2 == banked n_cert on every bank."),
        gate_pass=True, n_banks=n_bank, n_cells=len(cells),
        trials_nat=TRIALS_NAT, qbar=QBAR,
        n_events_v22=tot22, n_events_v23=tot23,
        n_wrong_v22=w22, n_wrong_v23=w23,
        recut=np.array([(r["cell"], r["n_iter"], str(r["ev22"]), str(r["ev23"]))
                        for r in recut], dtype=object),
        recut_keys=np.array(["cell", "n_iter", "ev22", "ev23"]),
        flicker=np.array([(f["cell"], f["psr"], str(f["iters"]), str(f["runs"]),
                           f["survives"], f["on_true"]) for f in flick], dtype=object),
        flicker_keys=np.array(["cell", "psr", "iters", "runs", "survives", "on_true"]),
    )
    print(f"[LEDGER-B2] banked -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
