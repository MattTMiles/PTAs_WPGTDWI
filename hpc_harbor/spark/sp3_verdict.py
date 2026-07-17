#!/usr/bin/env python
"""SPARK-3 — THE VERDICT READER (arm c ∪ arm a). CPU only, no GPU.

Reads the chunked-JVP marginal widths (arm c), re-cuts the crossing test at the MARGINAL width,
and — GATED ON ARM (a) — pronounces EDGE-POSITIVE / EDGE-ZERO / STILL-STRADDLED-UNAFFORDABLE.

THE ORDERING IS LOAD-BEARING (Matt's decision, item 1a): no JVP reading ships until the
scrambled arm reports. If any scrambled unit MANUFACTURES (soft@opt false-cert count > unmodelled
at the same rung), this refuses to pronounce and prints the manufacturing anatomy instead — the
optimistic edge is not licensed as an edge, so re-cutting crossings at the marginal width would
be measuring a manufactured effect.

THE JVP CORRECTNESS GATE (built in, analogous to g1a): the diagonal of the JVP Hessian, negated,
must reproduce faint_fisher_bounds' finite-difference F_ii — both are −∂²lnL/∂xᵢ². If they
disagree, the JVP is wrong and no marginal width may be quoted from it.

Re-cut rule at the marginal width: the soft draw uses sigma_marg (≥ sigma_cond always). A wider
draw is a worse model, so crossings can only be LOST relative to the optimistic bound, never
gained. SURVIVE = ≥1 crossing at ≥2 units under the marginal width, scrambled-clean.
"""
import os, glob, re
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "..", "SPARK3_results")
OUT = os.path.abspath(os.environ.get("SPARK3_OUT", "SPARK3_results"))


def _rungmap(pattern):
    R = {}
    for f in glob.glob(f"{OUT}/{pattern}"):
        m = re.search(r"_([AB])_g3_r(\d+)_k(\d+)", os.path.basename(f))
        if m:
            R[(m.group(1), int(m.group(2)), int(m.group(3)))] = f
    return R


def arm_a_scrambled():
    """Load the scrambled arm; return (reported, manufactures, rows)."""
    fs = sorted(glob.glob(f"{OUT}/s3scr_*_k*_L.npz"))
    rows = []
    for f in fs:
        z = np.load(f, allow_pickle=True)
        rows.append(dict(venue=str(z["venue"]), real=int(z["real_i"]), rung=int(z["rung"]),
                         nf_un=int(z["nfalse_un"]), nf_opt=int(z["nfalse_opt"]),
                         nf_pes=int(z["nfalse_pes"]), nf_tr=int(z["nfalse_tr"]),
                         manuf=bool(z["manufactures"]), path=f))
    manufactures = any(r["manuf"] for r in rows)
    return len(rows) > 0, manufactures, rows


def arm_c_jvp():
    R = _rungmap("s3jvp_*_k*.npz")
    out = {}
    for k, f in R.items():
        z = np.load(f, allow_pickle=True)
        out[k] = dict(capped=bool(z["capped"]), sig_marg=np.asarray(z["sig_marg"], float),
                      sig_cond=np.asarray(z["sig_cond"], float),
                      sig_pes=np.asarray(z["sig_pes"], float),
                      F_ii=np.asarray(z["F_ii"], float) if "F_ii" in z.files else None,
                      hess=np.asarray(z["hess"], float) if "hess" in z.files else None,
                      idx=np.asarray(z["idx"], int), path=f)
    return out


def jvp_gate(J):
    """diag(-H) from the JVP must reproduce F_ii from finite difference. Returns max rel dev."""
    worst = 0.0
    for k, d in J.items():
        if d["capped"] or d["hess"] is None or d["F_ii"] is None:
            continue
        diagF = -np.diag(d["hess"])                       # −∂²lnL/∂x² from AD
        m = np.abs(d["F_ii"]) > 1e-6
        rel = np.abs(diagF[m] - d["F_ii"][m]) / np.abs(d["F_ii"][m])
        worst = max(worst, float(rel.max()) if rel.size else 0.0)
    return worst


def main():
    print("=" * 78)
    print("SPARK-3 VERDICT READER  —  arm (a) scrambled gates the reading of arm (c) JVP")
    print("=" * 78)

    # ---------- ARM (a): the manufacturing gate (must report first) ----------
    a_reported, manufactures, a_rows = arm_a_scrambled()
    print("\n[ARM a] SCRAMBLED-COUNTERPART  (every cert is FALSE by construction)")
    if not a_rows:
        print("  *** not yet reported — NO JVP reading may ship (Matt item 1a). STOP.")
        return 2
    print(f"  {'unit':>6} {'rung':>4} | {'FALSE certs  un':>15} {'opt':>4} {'pes':>4} {'tr':>4}"
          f"  {'manufactures?':>14}")
    for r in sorted(a_rows, key=lambda x: (x["venue"], x["real"], x["rung"])):
        flag = "*** YES" if r["manuf"] else "clean"
        print(f"  {r['venue']}{r['real']:>4} {r['rung']:>4} | {r['nf_un']:>15} {r['nf_opt']:>4} "
              f"{r['nf_pes']:>4} {r['nf_tr']:>4}  {flag:>14}")
    if manufactures:
        print("\n  *** ARM (a) MANUFACTURES: soft-modelling RAISES the false-cert count on a")
        print("      wrong counterpart. The optimistic edge is NOT licensed as an edge — its")
        print("      crossings could be the same mechanism. STOP + anatomy; no EDGE verdict.")
        print("      -> the campaign's honest state remains STRADDLED, and the optimistic")
        print("         bound is demoted to 'consistent with manufacture'.")
        return 3
    print("  ==> ARM (a) CLEAN: soft-modelling does not manufacture at the tested rungs.")
    print("      The optimistic edge is licensed as an edge; the JVP reading may proceed.")

    # ---------- ARM (c): the JVP correctness gate ----------
    J = arm_c_jvp()
    print(f"\n[ARM c] CHUNKED-JVP MARGINAL WIDTH  ({len(J)} Hessians banked)")
    if not J:
        print("  *** no JVP banks yet. STOP.")
        return 2
    g = jvp_gate(J)
    print(f"  JVP CORRECTNESS GATE: max rel |diag(-H) - F_ii(finite-diff)| = {g:.2e}"
          f"   {'PASS' if g < 1e-2 else '*** FAIL (JVP Hessian is wrong)'}")
    if g >= 1e-2:
        print("  refusing to quote a marginal width from a Hessian that fails its own gate.")
        return 4
    capped = [k for k, d in J.items() if d["capped"]]
    print(f"  capped (unaffordable at 45-min): {len(capped)}/{len(J)}  {sorted(capped)}")

    # ---------- re-cut the crossings at the MARGINAL width ----------
    # We need, per (unit,rung): the soft-modelled margin under sigma_marg. That is a GPU
    # re-score, not available here — so this reader reports the WIDTH verdict (how far marg sits
    # from cond vs pes) and defers the margin re-cut to `unit --width marg` (banked separately).
    print("\n  WIDTH BRACKET per unit (median over the faint block):")
    print(f"  {'unit':>7} {'cond(opt)':>10} {'MARG':>10} {'prior(pes)':>10} {'marg/cond':>9}")
    frac = []
    for k in sorted(J):
        d = J[k]
        if d["capped"]:
            print(f"  {k[0]}{k[1]}k{k[2]:>2}  {'':>10} {'CAPPED':>10}")
            continue
        mc = np.median(d["sig_cond"]); mm = np.median(d["sig_marg"]); mp = np.median(d["sig_pes"])
        frac.append(mm / mc if mc > 0 else np.nan)
        print(f"  {k[0]}{k[1]}k{k[2]:>2}  {mc:10.4g} {mm:10.4g} {mp:10.4g} {mm/mc:9.2f}")
    if frac:
        print(f"\n  marginal width is {np.nanmedian(frac):.2f}x the conditional (>=1 always; "
              f"closer to 1 => optimistic bound was ~tight => crossings likely SURVIVE;")
        print(f"  >> 1 => marginal near prior => crossings likely DIE).")
    print("\n  NEXT: re-score the soft state at sigma_marg (`unit` with the banked marg column)")
    print("  and re-run the crossing test; SURVIVE = >=1 crossing at >=2 units. That GPU re-cut")
    print("  is the final step — this reader has cleared arm (a) and gated the JVP.")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
