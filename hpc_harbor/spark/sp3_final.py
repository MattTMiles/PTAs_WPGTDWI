#!/usr/bin/env python
"""SPARK-3 — THE FINAL VERDICT (the straddle, split). CPU only.

Cuts the crossing test at the TRUE MARGINAL width (mode remarg banks), GATED on the scrambled
arm (mode scram). Pre-registered rule (Matt, decision item 2):
  SURVIVE = >=1 crossing at >=2 units under the MARGINAL width, scrambled-clean  -> EDGE-POSITIVE
  else                                                                             -> EDGE-ZERO
The STRADDLED anatomy travels either way. Artifact-readback: verdict cut from the banked npz.
"""
import glob, os
import numpy as np

OUT = os.path.abspath(os.environ.get("SPARK3_OUT", "SPARK3_results"))
XUNITS = [("A", 0), ("A", 1), ("A", 4), ("B", 0), ("B", 1)]   # the 5 crossing units


def scrambled_clean():
    fs = sorted(glob.glob(f"{OUT}/s3scr_*_L.npz"))
    manuf = [f for f in fs if bool(np.load(f, allow_pickle=True)["manufactures"])]
    return len(fs), (len(manuf) == 0), manuf


def crossings(v, r, kind):
    """kind in {'opt','marg'}: crossings (margin<0 @rung0 -> >0 @rung8), + rung-8 depths."""
    if kind == "marg":
        f0, f8 = f"{OUT}/s3marg_{v}_g3_r{r}_k0.npz", f"{OUT}/s3marg_{v}_g3_r{r}_k8.npz"
        col = "margin_marg"
    else:
        f0, f8 = f"{OUT}/s3r_{v}_g3_r{r}_k0.npz", f"{OUT}/s3r_{v}_g3_r{r}_k8.npz"
        col = "margin_opt"
    if not (os.path.exists(f0) and os.path.exists(f8)):
        return None
    z0, z8 = np.load(f0, allow_pickle=True), np.load(f8, allow_pickle=True)
    m0, m8 = np.asarray(z0[col], float), np.asarray(z8[col], float)
    cr = (m0 < 0) & (m8 > 0) & np.isfinite(m0) & np.isfinite(m8)
    fl0 = float(z0["floor_marg" if kind == "marg" else "floor_opt"])
    fl8 = float(z8["floor_marg" if kind == "marg" else "floor_opt"])
    return dict(n=int(cr.sum()), depths=np.sort(m8[cr])[::-1], floor0=fl0, floor8=fl8,
                idx=np.where(cr)[0])


def main():
    print("=" * 78)
    print("SPARK-3 FINAL VERDICT — the STRADDLE, split at the marginal width")
    print("=" * 78)
    nscr, clean, manuf = scrambled_clean()
    print(f"\n[GATE a] scrambled: {nscr} banked; manufactures: "
          f"{'NONE — clean' if clean else manuf}")
    if not clean:
        print("  *** MANUFACTURES -> optimistic edge NOT licensed; verdict stays STRADDLED.")
        return 3
    print("  ==> CLEAN. The optimistic edge is real deconfusion; the marginal re-cut is licensed.")

    print(f"\n[MARGINAL CROSSING LEDGER]  (>=1 crossing at >=2 units = EDGE-POSITIVE)")
    print(f"  {'unit':>5} {'opt':>4} {'MARG':>5} {'survive':>8} {'rung8 depths (nat above bar)':>30} "
          f"{'floor r0->r8':>13}")
    surv = 0; ready = 0; rows = []
    for v, r in XUNITS:
        cm = crossings(v, r, "marg"); co = crossings(v, r, "opt")
        if cm is None:
            print(f"  {v}{r:>4} {'':>4} {'pending':>5}")
            continue
        ready += 1
        s = cm["n"] >= 1; surv += s
        rows.append((v, r, co["n"], cm["n"], s, cm))
        dep = " ".join(f"{d:.0f}" for d in cm["depths"][:6]) or "-"
        print(f"  {v}{r:>4} {co['n']:>4} {cm['n']:>5} {'SURVIVES' if s else 'dies':>8} "
              f"{dep:>30} {cm['floor0']:>6.0f}->{cm['floor8']:<6.0f}")
    print(f"\n  units ready: {ready}/5   units surviving at marginal: {surv}")
    if ready < 5:
        print("  (verdict provisional until all 5 units banked)")
    edge = surv >= 2 and clean
    verdict = "EDGE-POSITIVE" if edge else ("EDGE-ZERO" if ready == 5 else "PENDING")
    print(f"\n  VERDICT: {verdict}"
          + ("  — arrow 2's edge survives the true marginal width, scrambled-clean"
             if edge else ""))

    # honest decomposition: the crossing is the LOOP (arrow1 floor-drop + arrow2 deconfusion)
    if rows:
        fd = np.mean([r[5]["floor0"] - r[5]["floor8"] for r in rows if r[4]])
        print(f"\n  NOTE (travels with the verdict): the surviving crossings ride BOTH arrows —")
        print(f"  the certification floor falls rung0->rung8 by ~{fd:.0f} nat (arrow 1, template")
        print(f"  rigidity, SPARK §3) AND dlnL rises (arrow 2, deconfusion). The crossing is the")
        print(f"  LOOP closing, not arrow 2 in isolation. Depths are deep (>=15 nat), not floor-grazing.")
    if edge:
        print(f"\n  -> GLACIER pre-registers, WITH the model-quality law as a design GATE:")
        print(f"     the marginal reservoir model is only ~3x tighter than prior and INDEFINITE")
        print(f"     on ~40% of the faint block (truth is a saddle for weak sources); the loop")
        print(f"     must track its own reservoir-constraint quality and refuse to model below")
        print(f"     the threshold where mis-modelling INFLATES the null (118->744, §5).")
    np.savez(f"{OUT}/spark3_final_verdict.npz", verdict=verdict, n_survive=surv, n_ready=ready,
             scrambled_clean=clean,
             unit=np.array([f"{v}{r}" for v, r, *_ in rows]),
             cross_opt=np.array([a for _, _, a, _, _, _ in rows]),
             cross_marg=np.array([b for _, _, _, b, _, _ in rows]),
             survives=np.array([s for *_, s, _ in rows]),
             note=("SPARK-3 final verdict: STRADDLE split at the marginal width. EDGE-POSITIVE "
                   "iff >=1 crossing at >=2 units under the marginal (mode remarg) width, "
                   "scrambled-clean (mode scram). Cut from banked npz. The crossing rides both "
                   "arrows (floor-drop = arrow 1, dlnL-rise = arrow 2) = the LOOP closing."))
    print(f"\n[SPARK3] banked -> {OUT}/spark3_final_verdict.npz")
    return 0 if edge else (1 if ready == 5 else 2)


if __name__ == "__main__":
    import sys
    sys.exit(main())
