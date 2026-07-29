#!/usr/bin/env python3
"""SIEVE T7 -- THE E-PROCESS SCOREBOARD (agent SIEVE-A, 2026-07-29). REPORT-ONLY.

QUESTION (brief): on all banked per-iteration scoreboard histories, put a running
likelihood-ratio e-process per pulsar beside the null floor distribution; compare
against the per-iteration alpha = 0.05 bar and against the persistence rule
(LEDGER B2 / criterion-v2.3). Readouts: (i) does any banked TRUE certification die
under e-values; (ii) do the 13 false cascade certs die EARLIER than under the D2
rigidity gate. Verdict: adopt-for-v2.3 / no-gain.

NOTHING HERE ARMS A PROTOCOL STEP, MOVES A BANKED VERDICT, OR ENTERS A CLOSURE CLAIM.
Every quantity is read from the banked RAW columns (dlnL_det, lnK, qmax, on_true,
floor, zero_fraction, full_floors) -- no floor is re-cut, no estimator is switched,
no likelihood is re-evaluated. `lean-npz discipline` is what makes this affordable.

--------------------------------------------------------------------------------
WHAT THE BANKED STATISTICS ACTUALLY ARE (read from the code, not assumed)

  spark3.score_from_LNL:
      dlnL_det[a] = FT.peak[a] top1 - FT.peak[a] top2       <- a FRINGE PEAK GAP
      lnK[a]      = log(max(K_counted[a], 1))               <- counted fringe trials
      qmax[a]     = fringe posterior mass on the MAP fringe
      on_true[a]  = (map_k == n_true_grid)

  glacier_loop._null_offenders (the floor's own statistic, verbatim):
      o[a] = max(dlnL[a] - max(lnK[a], 0), 0)  if qmax[a] > QBAR else 0
      floor = q95 of  max_a o[a]  over n_null no-CW draws
      zero_fraction = fraction of those draws with max_a o[a] == 0

  glacier_loop.run_cell, criterion-v2.2:
      cert[a] = (dlnL[a] > max(lnK[a] + 0.578, floor)) & (qmax[a] > QBAR)

  NOTE THE SCALE the floor is measured on -- see E0 below.

--------------------------------------------------------------------------------
THE THREE CONSTRUCTIONS SCORED HERE

E0  SCALE-MATCHED FLOOR (diagnostic, not an e-value).  The floor is the q95 of the
    OFFENDER statistic o = (dlnL - lnK)_+, but criterion-v2.2 applies it to dlnL.
    E0 re-applies the same banked floor on the scale it was measured on:
        cert_E0[a] = (o[a] > max(TRIALS_NAT, floor)) & (qmax[a] > QBAR)
    This changes no bank and re-cuts no floor; it is a read of how much of the
    v2.2 bar is carried by the mismatch. Reported as a number, not a proposal.

E1  FLOOR-CALIBRATED OFFENDER E-VALUE (the brief's "vs the null floor distribution").
    Two banked points calibrate the null survival of the CELL-MAX offender:
        S(0+) = 1 - zero_fraction        S(floor) = 0.05   (floor is the q95)
    -> exponential upper tail (the Gumbel tail the floor estimator already assumes):
        S(x) = (1 - zf) * exp(-x / beta),   beta = floor / ln((1 - zf) / 0.05)
    p[a] := S(o[a]) is a CONSERVATIVE (super-uniform) per-pulsar p-value, because it
    is calibrated on the max over 116 pulsars. Then the standard p-to-e calibrator
        e = kappa * p**(kappa - 1),  kappa = 1/2  ->  e = 0.5 / sqrt(p)
    which is a valid e-value for any super-uniform p.
    DEGENERATE branch: if zf >= 0.95 or floor == 0 the null sample has no resolvable
    tail; the rule of three caps P(max o > 0) <= 3/n_null, so p = 3/n_null and no
    smaller p is readable from the bank. RESOLUTION FLAG: any p below 1/(n_null+1)
    is EXTRAPOLATED past the null sample and is labelled as such.

E2  FRINGE-TRIALS E-VALUE.  e2[a] = exp(dlnL[a] - lnK[a]) gated on qmax > QBAR.
    Honest statement: dlnL_det is a top1-top2 peak GAP, not a null-vs-signal
    likelihood ratio, so e2 is not a Bayes factor. It is the criterion's OWN trials
    bar re-expressed on the e-scale: v2.2 demands e2 > exp(0.578) = 1.78, whereas
    alpha = 0.05 demands e2 >= 20 (i.e. dlnL > lnK + 3.00). E2 is quoted to price
    the trials term, not to replace the floor.

RUNNING E-PROCESS across the 6 scored iterations, per (cell, pulsar):
    PRODUCT  prod_i e_i  -- valid only under CONDITIONAL INDEPENDENCE across
        iterations, which this loop VIOLATES (the same data are re-scored at an
        evolving template). Reported as the anticonservative upper bound, flagged.
    MEAN     (1/n) sum_i e_i -- a valid e-value at every fixed n by linearity of
        expectation, with NO independence assumption. The honest readout.
    Bar: e >= 1/alpha = 20.

GATE (must pass before any downstream number is read): the re-derived v2.2 cert mask
must reproduce the banked `n_cert` column EXACTLY on every bank -- LEDGER B2's gate,
same form.

Output bank: reports/sieve_t7_eprocess.npz
Usage: python3 hpc_harbor/ledger/sieve_t7_eprocess.py [--repo <path>]
"""
import argparse
import os
import re
from collections import defaultdict

import numpy as np

TRIALS_NAT = 0.578          # glacier_loop.py:106
QBAR = 0.9                  # glacier_loop.py:104
ALPHA = 0.05
E_BAR = 1.0 / ALPHA         # 20
N_NULL_FULL = 100           # glacier_loop.py:107
N_NULL_CARRY = 32           # glacier_loop.py run_cell step (f)

STEM_RE = re.compile(
    r"^(?P<pre>gl1|gl2)_"
    r"(?:(?P<rung>r13p\d+)(?:_w(?P<w>[0-9p]+))?_)?"
    r"(?P<kind>cell|null\d+)_"
    r"(?P<arm>e07|none)_s(?P<sky>\d+)_T(?P<T>\d+)_(?P<tier>\w+?)_i(?P<iter>\d+)$"
)

# The pre-registered D2 population (S4.21.1): the 13 cascade wrong certs live in
# these two cells. Used only to LABEL rows, never to select them.
CASCADE_CELLS = {("gl2", "r13p25", 1.0, 30, "e07", 0, "cell"),
                 ("gl2", "r13p25", 1.0, 30, "e07", 1, "cell")}


def cert_v22(dlnl, lnK, q, fl):
    """glacier_loop.py:533 verbatim, from the banked RAW columns."""
    return (dlnl > np.maximum(lnK + TRIALS_NAT, fl)) & (q > QBAR)


def offender(dlnl, lnK, q):
    """glacier_loop._null_offenders:594 verbatim -- the statistic the FLOOR is cut on."""
    return np.where(q > QBAR, np.maximum(dlnl - np.maximum(lnK, 0.0), 0.0), 0.0)


def e1_from_floor(o, fl, zf, n_null):
    """Floor-calibrated e-value + its p and a resolution flag.

    Returns (e, p, flag) with flag in {'tail', 'degenerate', 'extrapolated'}.
    """
    p_res = 1.0 / (n_null + 1.0)            # smallest p the null SAMPLE resolves
    if o <= 0.0:
        return 0.5, 1.0, "null"             # e = 0.5/sqrt(1) at p = 1
    tail_mass = 1.0 - zf
    if fl <= 0.0 or tail_mass <= ALPHA:
        # no resolvable tail: rule of three on n_null silent draws
        p = min(1.0, 3.0 / n_null)
        return 0.5 / np.sqrt(p), p, "degenerate"
    beta = fl / np.log(tail_mass / ALPHA)
    p = float(min(1.0, tail_mass * np.exp(-o / beta)))
    flag = "extrapolated" if p < p_res else "tail"
    return 0.5 / np.sqrt(max(p, 1e-300)), p, flag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/taylor_group/matt_miles/PTAs_WPGTDWI")
    a = ap.parse_args()
    rdir = os.path.join(a.repo, "GLACIER_results")

    # ---------------- load ----------------
    cells = defaultdict(dict)
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
        dlnl = np.asarray(z["dlnL_det"], float)
        lnK = np.asarray(z["lnK"], float)
        q = np.asarray(z["qmax"], float)
        fl = float(z["floor"])
        ff = bool(z["full_floors"]) if "full_floors" in z.files else True
        zf = float(z["zero_fraction"]) if "zero_fraction" in z.files else 0.0
        cells[key][int(g["iter"])] = dict(
            dlnl=dlnl, lnK=lnK, q=q, floor=fl, zf=zf, full_floors=ff,
            n_null=N_NULL_FULL if ff else N_NULL_CARRY,
            floor_est=str(z["floor_est"]) if "floor_est" in z.files else "?",
            cert=cert_v22(dlnl, lnK, q, fl), o=offender(dlnl, lnK, q),
            banked=int(z["n_cert"]), on_true=np.asarray(z["on_true"], bool), fn=fn)
    n_bank = sum(len(v) for v in cells.values())
    print(f"[SIEVE-T7] {n_bank} per-iteration banks over {len(cells)} cells\n")

    # ---------------- GATE ----------------
    bad = [(r["fn"], int(r["cert"].sum()), r["banked"])
           for it in cells.values() for r in it.values()
           if int(r["cert"].sum()) != r["banked"]]
    print("-- GATE: re-derived criterion-v2.2 mask vs the banked n_cert column --")
    if bad:
        print(f"   FAIL: {len(bad)}/{n_bank} banks disagree. First 5:")
        for fn, x, y in bad[:5]:
            print(f"     {fn}: re-cut {x} vs banked {y}")
        print("   STOPPING -- no downstream number from this run is quotable.")
        return 1
    print(f"   PASS: {n_bank}/{n_bank} banks reproduce exactly\n")

    # ---------------- E0: the floor's own scale ----------------
    print("-- E0: THE FLOOR IS CUT ON (dlnL - lnK)_+ BUT APPLIED TO dlnL --")
    print("   per certification event: the banked bar vs the scale-matched bar")
    print("   cell/psr@iter                                     dlnL    lnK   floor "
          " v2.2bar  matchedbar  v2.2 E0  on_true")
    events = []
    for k in sorted(cells, key=str):
        for i in sorted(cells[k]):
            r = cells[k][i]
            for p in np.where(r["cert"])[0]:
                bar22 = max(r["lnK"][p] + TRIALS_NAT, r["floor"])
                barM = max(TRIALS_NAT, r["floor"]) + max(r["lnK"][p], 0.0)
                e0 = bool(r["o"][p] > max(TRIALS_NAT, r["floor"]))
                e1, p1, fg = e1_from_floor(r["o"][p], r["floor"], r["zf"], r["n_null"])
                e2 = float(np.exp(r["dlnl"][p] - r["lnK"][p]))
                ev = dict(cell="/".join(str(x) for x in k), key=k, psr=int(p), iter=i,
                          dlnl=float(r["dlnl"][p]), lnK=float(r["lnK"][p]),
                          q=float(r["q"][p]), floor=r["floor"], zf=r["zf"],
                          n_null=r["n_null"], o=float(r["o"][p]),
                          bar22=float(bar22), barM=float(barM), e0=e0,
                          e1=float(e1), p1=float(p1), flag=fg, e2=e2,
                          on_true=bool(r["on_true"][p]),
                          cascade=(k in CASCADE_CELLS))
                events.append(ev)
                print(f"   {ev['cell']}/psr{p}@i{i}".ljust(50)
                      + f"{ev['dlnl']:7.2f} {ev['lnK']:6.2f} {ev['floor']:6.2f} "
                      f"{bar22:7.2f}  {barM:9.2f}  YES {'YES' if e0 else 'NO '} "
                      f"  {ev['on_true']}")
    n_ev = len(events)
    n_true = sum(1 for e in events if e["on_true"])
    n_e0 = sum(1 for e in events if e["e0"])
    n_e0_true = sum(1 for e in events if e["e0"] and e["on_true"])
    print(f"\n   v2.2: {n_ev} events ({n_ev - n_true} wrong / {n_true} on_true)")
    print(f"   E0  : {n_e0} events ({n_e0 - n_e0_true} wrong / {n_e0_true} on_true)"
          f"  -- scale-matching the SAME banked floor removes {n_ev - n_e0}\n")

    # ---------------- E1 / E2 per event ----------------
    print("-- E1 (floor-calibrated) and E2 (fringe-trials) per certification event --")
    print("   bar: e >= 1/alpha = 20      [flag] = null-sample resolution of p")
    print("   cell/psr@iter                                       o    p(E1)      E1"
          "        E2   E1>=20 E2>=20  flag          on_true")
    for e in events:
        print(f"   {e['cell']}/psr{e['psr']}@i{e['iter']}".ljust(50)
              + f"{e['o']:6.2f} {e['p1']:9.2e} {e['e1']:7.2f} {e['e2']:9.3g}  "
              f"{'Y' if e['e1'] >= E_BAR else 'n'}      "
              f"{'Y' if e['e2'] >= E_BAR else 'n'}    {e['flag']:13s} {e['on_true']}")
    for nm, key in (("E1", "e1"), ("E2", "e2")):
        s = [e for e in events if e[key] >= E_BAR]
        st = sum(1 for e in s if e["on_true"])
        print(f"\n   {nm} at alpha=0.05: {len(s)} events survive "
              f"({len(s) - st} wrong / {st} on_true) of {n_ev} "
              f"({n_ev - n_true} wrong / {n_true} on_true)")
    print()

    # ---------------- running e-process ----------------
    print("-- RUNNING E-PROCESS per (cell, pulsar) over the scored iterations --")
    print("   PRODUCT is anticonservative (iterations re-score the SAME data);")
    print("   MEAN is a valid e-value at every fixed n with no independence assumption.")
    print("   cell/psr                                    n_it   max_i E1   "
          "prod E1     mean E1   max_i E2      prod E2      mean E2  on_true")
    procs = []
    for k in sorted(cells, key=str):
        its = sorted(cells[k])
        allp = sorted({int(p) for i in its for p in np.where(cells[k][i]["cert"])[0]})
        for p in allp:
            e1s, e2s = [], []
            for i in its:
                r = cells[k][i]
                ev1, _, _ = e1_from_floor(r["o"][p], r["floor"], r["zf"], r["n_null"])
                e1s.append(ev1)
                e2s.append(float(np.exp(r["dlnl"][p] - r["lnK"][p]))
                           if r["q"][p] > QBAR else 0.0)
            e1s, e2s = np.array(e1s), np.array(e2s)
            hit = [i for i in its if cells[k][i]["cert"][p]]
            ot = bool(cells[k][hit[0]]["on_true"][p])
            with np.errstate(over="ignore"):
                pr1, pr2 = float(np.prod(e1s)), float(np.prod(e2s))
            row = dict(cell="/".join(str(x) for x in k), psr=int(p), n_it=len(its),
                       max_e1=float(e1s.max()), prod_e1=pr1, mean_e1=float(e1s.mean()),
                       max_e2=float(e2s.max()), prod_e2=pr2, mean_e2=float(e2s.mean()),
                       on_true=ot, cascade=(k in CASCADE_CELLS), iters=hit)
            procs.append(row)
            print(f"   {row['cell']}/psr{p}".ljust(46)
                  + f"{len(its):3d} {row['max_e1']:10.2f} {pr1:10.3g} "
                  f"{row['mean_e1']:10.2f} {row['max_e2']:10.3g} {pr2:12.3g} "
                  f"{row['mean_e2']:12.3g}  {ot}")
    print()

    # ---------------- readouts ----------------
    print("-- READOUT 1: do banked TRUE certifications die under e-values? --")
    tr = [e for e in events if e["on_true"]]
    for e in tr:
        print(f"   {e['cell']}/psr{e['psr']}@i{e['iter']}: dlnL {e['dlnl']:.2f} "
              f"lnK {e['lnK']:.2f} floor {e['floor']:.3f} zf {e['zf']:.2f} "
              f"n_null {e['n_null']} -> E1 {e['e1']:.2f} ({e['flag']}), "
              f"E2 {e['e2']:.2f}  ==> "
              f"{'SURVIVES' if max(e['e1'], e['e2']) >= E_BAR else 'DIES'} at e>=20")
    n_true_die = sum(1 for e in tr if max(e["e1"], e["e2"]) < E_BAR)
    print(f"   {n_true_die}/{len(tr)} banked true certifications DIE under e >= 20.\n")

    print("-- READOUT 2: do the 13 cascade wrong certs die EARLIER than the D2 gate? --")
    casc = [e for e in events if e["cascade"]]
    surv1 = [e for e in casc if e["e1"] >= E_BAR]
    surv2 = [e for e in casc if e["e2"] >= E_BAR]
    print(f"   cascade population: {len(casc)} events "
          f"({sum(1 for e in casc if not e['on_true'])} wrong)")
    tag = lambda e: "sky%s/psr%d@i%d" % (e["cell"].split("/")[5], e["psr"], e["iter"])
    print(f"   survive E1>=20: {len(surv1)}  {[tag(e) for e in surv1]}")
    print(f"   survive E2>=20: {len(surv2)}  {[tag(e) for e in surv2]}")
    for e in casc:
        print(f"     sky{e['cell'].split('/')[5]}/psr{e['psr']}@i{e['iter']}: "
              f"E1 {e['e1']:8.2f} {'SURV' if e['e1'] >= E_BAR else 'kill'}   "
              f"E2 {e['e2']:11.3g} {'SURV' if e['e2'] >= E_BAR else 'kill'}   "
              f"(D2 rigidity: KILLED on R2)")
    print()

    out = os.path.join(a.repo, "reports", "sieve_t7_eprocess.npz")
    np.savez(
        out,
        note=("SIEVE T7: e-process scoreboard over the banked GLACIER per-iteration "
              "histories. E0 = the banked floor re-applied on the (dlnL-lnK)_+ scale "
              "it was cut on; E1 = floor-calibrated offender e-value (two-point "
              "exponential tail from (zero_fraction, floor=q95), kappa=1/2 p-to-e "
              "calibrator); E2 = exp(dlnL-lnK), the criterion's trials bar on the "
              "e-scale. REPORT-ONLY: no floor re-cut, no likelihood re-evaluated, no "
              "banked verdict moved. Gate passed: re-derived v2.2 == banked n_cert."),
        gate_pass=True, n_banks=n_bank, n_cells=len(cells),
        trials_nat=TRIALS_NAT, qbar=QBAR, alpha=ALPHA, e_bar=E_BAR,
        n_events=n_ev, n_wrong=n_ev - n_true, n_true=n_true,
        n_e0=n_e0, n_e0_true=n_e0_true,
        n_e1_survive=sum(1 for e in events if e["e1"] >= E_BAR),
        n_e2_survive=sum(1 for e in events if e["e2"] >= E_BAR),
        n_true_die=n_true_die,
        events=np.array([[e["cell"], e["psr"], e["iter"], e["dlnl"], e["lnK"],
                          e["q"], e["floor"], e["zf"], e["n_null"], e["o"],
                          e["bar22"], e["barM"], e["e0"], e["p1"], e["e1"],
                          e["flag"], e["e2"], e["on_true"], e["cascade"]]
                         for e in events], dtype=object),
        events_keys=np.array(["cell", "psr", "iter", "dlnL", "lnK", "qmax", "floor",
                              "zero_fraction", "n_null", "offender", "bar_v22",
                              "bar_scalematched", "cert_E0", "p_E1", "E1", "E1_flag",
                              "E2", "on_true", "cascade"]),
        procs=np.array([[r["cell"], r["psr"], r["n_it"], r["max_e1"], r["prod_e1"],
                         r["mean_e1"], r["max_e2"], r["prod_e2"], r["mean_e2"],
                         r["on_true"], r["cascade"], str(r["iters"])]
                        for r in procs], dtype=object),
        procs_keys=np.array(["cell", "psr", "n_iter", "max_E1", "prod_E1", "mean_E1",
                             "max_E2", "prod_E2", "mean_E2", "on_true", "cascade",
                             "cert_iters"]),
    )
    print(f"banked -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
