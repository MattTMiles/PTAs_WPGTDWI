#!/usr/bin/env python
"""PHOENIX -- FROZEN vs LIVE M-step: the readback.

Discharges GLACIER S4.15.1 item 2. Compares the frozen arm (FROZEN_results/, cut by
frozen_mstep.py) against the INHERITED live arm (GLACIER_results/, 24/24 banked), cell
for cell, at r13p9 / circular / T=30 / w=1, skies 0-3.

THE QUESTION, as pre-registered: does refit-once-then-freeze RESCUE THE DRAIN that the
live M-step gives back? The live arm's signature (S4.15 ACTIVE-REGRESSIVE, SUMMIT S0.1)
is: real first bite at i0, then a regression at i1, then parking somewhere between the
bite and baseline. If freezing the template after its first fit stops the give-back, the
regression is the M-step's doing. If the give-back survives freezing, it is the drain
refit's doing and the M-step is exonerated.

WHY THIS COMPARISON IS IMMUNE TO THE LEDGER-B1 DRAIN DEFECT (stated up front, because it
is the reason this arm is worth reading at all): LEDGER measures the feed-free reference
offset at -0.389 +- 0.374 dex and shows `bite = a_bg - a_eff` is the wrong absolute
reference, two-sided by arm. **This script never uses that reference for its verdict.**
Frozen and live are a PAIRED contrast -- identical census seed, identical noise seed,
identical venue, identical driver; the ONLY difference is whether the M-step writes back
after the first fit. Any common reference offset cancels exactly in the difference
`a_bg_frozen(i) - a_bg_live(i)`. Absolute bites are printed for continuity with S4.15 and
carry LEDGER's caveat; THE VERDICT IS READ OFF THE PAIRED DIFFERENCE ONLY.

CPU only.
"""
import glob, re, sys
import numpy as np

LIVE = "GLACIER_results"
FRZ = "FROZEN_results"
STEM = "gl2_r13p9_cell_none_s{s}_T30_lit"
NP_SRC, I_MC, I_FGW = 8, 3, 4
HOLD_DEX = 0.05


def load(d, s):
    fs = sorted(glob.glob(f"{d}/{STEM.format(s=s)}_i*__*.npz"),
                key=lambda f: int(re.search(r"_i(\d+)__", f).group(1)))
    fs = [f for f in fs if "STOPANAT" not in f]
    if not fs:
        return None, None
    its = [int(re.search(r"_i(\d+)__", f).group(1)) for f in fs]
    return its, [np.load(f, allow_pickle=True) for f in fs]


def series(its, Z):
    abg = np.array([float(z["a_bg"]) for z in Z])
    aeff = float(Z[0]["a_eff_drawn"])
    nres = np.array([int(z["n_resolved"]) for z in Z])
    ncert = np.array([int(z["n_cert"]) for z in Z])
    wrong = np.array([int(z["wrong_cert"]) for z in Z])
    nd = [len(np.asarray(z["theta_rec"], float)) - NP_SRC * int(z["n_slot"]) for z in Z]
    fedth = []
    for z, n in zip(Z, nd):
        th = np.asarray(z["theta_rec"], float)
        fedth.append({int(k): (th[n + NP_SRC*k + I_FGW], th[n + NP_SRC*k + I_MC])
                      for k in np.flatnonzero(z["fed_mask"])})
    w, ncom = [], []
    for a, b in zip(fedth[:-1], fedth[1:]):
        c = set(a) & set(b)
        ncom.append(len(c))
        w.append(np.nan if not c else
                 max(max(abs(b[k][0]-a[k][0]), abs(b[k][1]-a[k][1])) for k in c))
    w = np.array(w, float)
    st = np.array([(not np.isnan(x)) and x <= HOLD_DEX and n > 0 for x, n in zip(w, ncom)], bool)
    hold = any(st[j] and st[j+1] for j in range(len(st)-1)) if len(st) > 1 else False
    return dict(its=its, abg=abg, aeff=aeff, nres=nres, ncert=ncert, wrong=wrong,
                wander=w, hold=hold, fedth=fedth)


rows = []
print("PER-SKY TRAJECTORIES -- a_bg by iteration (LIVE vs FROZEN), and the paired difference")
print("=" * 104)
for s in range(4):
    itL, ZL = load(LIVE, s)
    itF, ZF = load(FRZ, s)
    if ZL is None or ZF is None:
        print(f"sky {s}: MISSING ({'live' if ZL is None else ''}{'frozen' if ZF is None else ''}) -- skipped")
        continue
    L, F = series(itL, ZL), series(itF, ZF)
    n = min(len(L["abg"]), len(F["abg"]))
    d = F["abg"][:n] - L["abg"][:n]
    print(f"\nsky {s}   (a_eff = {L['aeff']:.3f}; fed at i0: live {L['nres'][0]}, frozen {F['nres'][0]})")
    print("  iter        " + "".join(f"{i:>9d}" for i in range(n)))
    print("  LIVE  a_bg  " + "".join(f"{v:>9.3f}" for v in L["abg"][:n]))
    print("  FROZ  a_bg  " + "".join(f"{v:>9.3f}" for v in F["abg"][:n]))
    print("  PAIRED diff " + "".join(f"{v:>+9.3f}" for v in d))
    print("  LIVE wander " + "".join(f"{v:>9.3f}" for v in L["wander"][:n-1]))
    print("  FROZ wander " + "".join(f"{v:>9.3f}" for v in F["wander"][:n-1]))
    # give-back = how much of the i0 bite is handed back by the last iteration
    giveL = float(L["abg"][n-1] - L["abg"][0])
    giveF = float(F["abg"][n-1] - F["abg"][0])
    step01 = float(L["abg"][1] - L["abg"][0]) if n > 1 else float("nan")
    rows.append(dict(sky=s, giveL=giveL, giveF=giveF, rescue=giveL - giveF, step01=step01,
                     dend=float(d[n-1]), holdL=L["hold"], holdF=F["hold"],
                     wL=float(np.nanmax(L["wander"])) if not np.all(np.isnan(L["wander"])) else np.nan,
                     wF=float(np.nanmax(F["wander"])) if not np.all(np.isnan(F["wander"])) else np.nan,
                     certL=int(L["ncert"].max()), certF=int(F["ncert"].max()),
                     wrongL=int(L["wrong"].sum()), wrongF=int(F["wrong"].sum()),
                     nresL=int(L["nres"].max()), nresF=int(F["nres"].max())))

if not rows:
    print("\nNo comparable cells yet -- fan has not drained.")
    sys.exit(0)

print("\n" + "=" * 104)
print("THE VERDICT TABLE -- give-back = a_bg(last) - a_bg(i0); positive = drain handed back")
print(f"{'sky':>4}{'giveback LIVE':>15}{'giveback FROZEN':>17}{'RESCUE':>9}"
      f"{'wanderL':>9}{'wanderF':>9}{'holdL':>7}{'holdF':>7}{'certL':>7}{'certF':>7}{'wrongF':>8}")
for r in rows:
    print(f"{r['sky']:>4}{r['giveL']:>+15.3f}{r['giveF']:>+17.3f}{r['rescue']:>+9.3f}"
          f"{r['wL']:>9.3f}{r['wF']:>9.3f}{int(r['holdL']):>7}{int(r['holdF']):>7}"
          f"{r['certL']:>7}{r['certF']:>7}{r['wrongF']:>8}")

gl = np.array([r["giveL"] for r in rows])
gf = np.array([r["giveF"] for r in rows])
res = gl - gf
print("-" * 104)
print(f"mean give-back  LIVE {gl.mean():+.3f}   FROZEN {gf.mean():+.3f}   "
      f"mean RESCUE {res.mean():+.3f} dex  (n = {len(rows)} skies)")
n_better = int((res > 0).sum())
print(f"skies where freezing retains MORE drain: {n_better}/{len(rows)}")

# WHERE does the give-back happen? If it is complete at the FIRST post-feed refit, the
# M-step's later wander cannot be its cause -- the two arms are identical up to and
# including the first fit by construction (freeze-after-FIRST-fit).
fed_rows = [r for r in rows if r["nresL"] > 0]
print("\nWHERE THE GIVE-BACK HAPPENS (fed skies only):")
print(f"{'sky':>4}{'giveback total':>16}{'at i0->i1':>12}{'frac at i0->i1':>16}{'after i1 (LIVE)':>18}")
for r in fed_rows:
    print(f"{r['sky']:>4}{r['giveL']:>+16.3f}{r['step01']:>+12.3f}"
          f"{100*r['step01']/r['giveL'] if r['giveL'] else float('nan'):>15.1f}%"
          f"{r['giveL']-r['step01']:>+18.3f}")

# Verdict on MAGNITUDE, not sign-count: is the rescue a material fraction of the give-back?
mean_give = float(np.mean([abs(r["giveL"]) for r in fed_rows])) if fed_rows else 0.0
mean_absres = float(np.mean([abs(r["rescue"]) for r in fed_rows])) if fed_rows else 0.0
ratio = mean_absres / mean_give if mean_give else float("nan")
print(f"\nmean |rescue| = {mean_absres:.3f} dex vs mean give-back = {mean_give:.3f} dex "
      f"-> {100*ratio:.1f}% of the effect, signed mean {res.mean():+.3f} dex")

print("\nVERDICT:")
if fed_rows and res.mean() > 0.05 and n_better == len(rows):
    print("  FROZEN RESCUES THE DRAIN in every sky. The live arm's give-back is the")
    print("  M-step writing back a wandering template -- freezing it after the first fit")
    print("  keeps the bite. This licenses the third arm (running-belief M-step).")
elif fed_rows and ratio < 0.35 and abs(res.mean()) < 0.05:
    print("  >>> FREEZING DOES NOT RESCUE THE DRAIN.")
    print(f"  The give-back is {100*np.mean([r['step01']/r['giveL'] for r in fed_rows]):.0f}% "
          "complete at the FIRST post-feed refit (i0->i1), where the two arms are")
    print("  IDENTICAL BY CONSTRUCTION (freeze-after-first-fit). Everything the M-step does")
    print("  afterwards is second-order jitter: signed mean rescue "
          f"{res.mean():+.3f} dex, |rescue| {100*ratio:.0f}% of the effect, both signs.")
    print("  => The M-step is EXONERATED as the cause of ACTIVE-REGRESSIVE. The regression")
    print("     lives in the DRAIN REFIT at constant fed set, not in template wander.")
else:
    print(f"  MIXED ({n_better}/{len(rows)} skies rescued, signed mean {res.mean():+.3f} dex,")
    print(f"  |rescue| {100*ratio:.0f}% of the give-back). Report per-sky; do not average.")

print("\n*** HOLD IS DEGENERATE IN THE FROZEN ARM -- do not read holdF as a result. ***")
print("    A frozen template has wander 0.000 by construction, so HOLD is TRUE for every")
print("    fed frozen cell trivially. holdF is printed for completeness only; the")
print("    M-step-trust contour is meaningful in the LIVE arm alone.")
print("\nSafety: frozen-arm wrong certs =", sum(r["wrongF"] for r in rows),
      "| cert counts live/frozen =", sum(r["certL"] for r in rows), "/",
      sum(r["certF"] for r in rows))

np.savez("reports/frozen_vs_live.npz",
         **{k: np.array([r[k] for r in rows]) for k in rows[0]},
         note="PHOENIX frozen-vs-live M-step, GLACIER S4.15.1 item 2. r13p9/none/T30/w1, "
              "skies 0-3, same lane (dgx03 A100-80GB), SG-F1 wrapper inertness PASS "
              "bit-exact. give-back = a_bg(last)-a_bg(i0); RESCUE = giveL - giveF. PAIRED "
              "contrast -- immune to the LEDGER-B1 feed-free reference offset, which "
              "cancels in the difference.")
print("\nbanked reports/frozen_vs_live.npz")
