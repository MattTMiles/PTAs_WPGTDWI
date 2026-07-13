"""THE FLOOR FIX, applied as far as the banked summaries allow — PROVISIONAL, pending re-score.

THE CONVENTION BEING APPLIED (ANCHOR 2026-07-13 §4, adopted verbatim):

    The D2 Gumbel-MLE floor is VALID ONLY where the nullN zero-fraction is <= 20%.
    Above that, the floor is the empirical (1 - alpha) quantile with a bootstrap error.
    The zero-fraction is a REQUIRED column beside every floor, not a caveat.

Rationale: the offender statistic is 0.0 whenever a realisation has no cell passing layer
1 (+) layer 3. At faint h that is MOST realisations, and a Gumbel fitted to a point mass at
zero is pulled DOWN toward it — understating the alpha = 0.05 bar by up to 2.8x, and erring
PERMISSIVE (a floor that is too low lets pure-noise offenders through).

WHAT THIS SCRIPT CAN AND CANNOT DO — the scope of inference, stated up front
---------------------------------------------------------------------------
It re-derives the FLOOR from banked summaries. It CANNOT re-cut the COUNT, because the
per-realisation signal dlnL columns live in SURFACE_results/ and CHORUS_results/ on ACCRE
and are not in this repo. Therefore:

  * every FLOOR below is final;
  * every COUNT below is the one banked under the OLD floor, and is quoted only to BOUND
    the corrected count using the one thing that is certain — the count is MONOTONE
    NON-INCREASING in the floor. Floor up => count can only fall. Floor down => count can
    only rise. Nothing here interpolates a count to a new floor. That would be a fabricated
    number, and it is exactly what the artifact-readback convention forbids.

Banks (read-only): reports/surface_analysis.npz, reports/ch_floors.npz,
reports/ch_analysis.npz, reports/ig2_floors.npz, reports/anchor_ladder.npz.
Writes: reports/floorfix_provisional.npz.

Run: python CW_transition/floorfix_provisional.py     (CPU only, seconds, no GPU, no jax)
"""
import os
import numpy as np

ALPHA = 0.05
ZF_GATE = 0.20                       # the ANCHOR §4 validity gate
BOOT = 4000
SEED = 20260713                      # fixed: the bootstrap is reproducible
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda f: os.path.join(ROOT, "reports", f)
STRUCT = {2: "2+14", 3: "3+13", 5: "5+11"}


def emp_quantile(off, alpha=ALPHA):
    """The empirical (1-alpha) quantile. Valid regardless of the zero point mass."""
    return float(np.quantile(np.asarray(off, dtype=float), 1.0 - alpha))


def boot_sd(off, alpha=ALPHA, B=BOOT, seed=SEED):
    """Bootstrap standard error of the empirical quantile."""
    rng = np.random.default_rng(seed)
    off = np.asarray(off, dtype=float)
    n = len(off)
    q = [np.quantile(rng.choice(off, n, replace=True), 1.0 - alpha) for _ in range(B)]
    return float(np.std(q))


def adopt(gumbel, gumbel_sd, off, zf):
    """Return (floor, err, estimator) under the ANCHOR §4 convention."""
    if zf <= ZF_GATE:
        return float(gumbel), float(gumbel_sd), "gumbel"
    return emp_quantile(off), boot_sd(off), "emp_q95"


def hr(t):
    print()
    print("=" * 100)
    print(t)
    print("=" * 100)


out = {}

# ==================================================================== SURFACE
# SURFACE banked fN_emp and fN_zerofrac per cell but NOT the raw offender vectors, so the
# corrected floor VALUE is recoverable for all 108 cells and its BOOTSTRAP ERROR is not.
hr("SURFACE — 108 cells. Floor correctable; count NOT (no per-realisation signal dlnL in repo).")
S = np.load(R("surface_analysis.npz"), allow_pickle=True)
cols = [str(c) for c in S["cols"]]
T, tiers = S["table"], np.array([str(t) for t in S["tiers"]])
V = np.array([str(v) for v in S["verdicts"]])
g = lambda c: T[:, cols.index(c)]
h, Tb, k = g("h"), g("T"), g("k")
fN, fN_sd, zf, emp = g("fN"), g("fN_sd"), g("fN_zerofrac"), g("fN_emp")
corr, corr_lo, corr_hi = g("corr"), g("corr_lo"), g("corr_hi")

touched = zf > ZF_GATE
adopted = np.where(touched, emp, fN)
onset = V == "ONSET"

print(f"banked verdicts: ONSET {int((V=='ONSET').sum())}  MARGINAL {int((V=='MARGINAL').sum())}"
      f"  below {int((V=='below').sum())}   (total {len(V)})")
print(f"cells with zero-fraction > {ZF_GATE:.0%} -> Gumbel INVALID, floor := empirical q95 : "
      f"{int(touched.sum())} / {len(V)}")
print(f"cells UNTOUCHED (zero-fraction <= {ZF_GATE:.0%}, Gumbel valid, floor and count stand): "
      f"{int((~touched).sum())} / {len(V)}")
print()

at_risk = touched & onset & (adopted > fN)          # floor up   -> count can only FALL
survives = onset & ~touched                          # untouched  -> stands, definitively
may_ignite = touched & ~onset & (adopted < fN)       # floor down -> count can only RISE

print(f"ONSET cells UNTOUCHED (stand, definitively)              : {int(survives.sum())}")
print(f"ONSET cells whose floor RISES  -> AT RISK, indeterminate : {int(at_risk.sum())}")
print(f"non-onset cells whose floor FALLS -> MAY IGNITE          : {int(may_ignite.sum())}")
print(f"=> corrected onset count is BOUNDED: {int(survives.sum())} <= N_onset <= "
      f"{int(survives.sum() + at_risk.sum() + may_ignite.sum())}   (banked, pre-fix: "
      f"{int(onset.sum())})")
print()
print("THE 15 TOUCHED CELLS (the entire blast radius of the floor fix on SURFACE):")
print(f"{'h':>7} {'T':>3} {'tier':>5} {'struct':>6} | {'zero-f':>6} {'Gumbel':>7} {'EMPq95':>7} "
      f"{'ratio':>6} | {'count':>6} {'@fl+sd':>7} | {'dir':>4} | verdict")
for i in np.argsort(-zf):
    if not touched[i]:
        continue
    d = "UP" if adopted[i] > fN[i] else "down"
    print(f"{h[i]:7.2f} {int(Tb[i]):3d} {tiers[i]:>5} {STRUCT[int(k[i])]:>6} | {zf[i]:6.2f} "
          f"{fN[i]:7.2f} {emp[i]:7.2f} {emp[i]/fN[i]:5.2f}x | {corr[i]:6.2f} {corr_lo[i]:7.2f} "
          f"| {d:>4} | {V[i]}")
print()
print("THE TWO ONSET CELLS AT RISK (both faint; both need the exact re-cut):")
for i in np.where(at_risk)[0]:
    bound = corr_lo[i] if adopted[i] > fN[i] + fN_sd[i] else corr[i]
    tag = "RETRACTED" if bound <= 1.0 else "INDETERMINATE"
    print(f"   h={h[i]:.2f} T={int(Tb[i])} {tiers[i]} {STRUCT[int(k[i])]}: floor "
          f"{fN[i]:.2f} -> {emp[i]:.2f} nat ({emp[i]/fN[i]:.2f}x), banked count {corr[i]:.2f}; "
          f"corrected count <= {bound:.2f} -> {tag}")
print()
print("NOTE: SURFACE banked no raw offender vectors, so the empirical floors above carry NO")
print("      bootstrap error. The convention REQUIRES one. ACCRE must bank the nullN offenders.")

out.update(sf_h=h, sf_T=Tb, sf_k=k, sf_tier=tiers, sf_zf=zf, sf_gumbel=fN, sf_gumbel_sd=fN_sd,
           sf_emp=emp, sf_adopted=adopted, sf_touched=touched, sf_verdict=V,
           sf_corr=corr, sf_corr_lo=corr_lo, sf_corr_hi=corr_hi,
           sf_at_risk=at_risk, sf_survives=survives, sf_may_ignite=may_ignite,
           sf_onset_lo=int(survives.sum()),
           sf_onset_hi=int(survives.sum() + at_risk.sum() + may_ignite.sum()),
           sf_onset_banked=int(onset.sum()))

# ==================================================================== CHORUS
# CHORUS DID bank its raw offender vectors -> full corrected floors WITH bootstrap errors.
hr("CHORUS — 26 cells. Raw offenders banked -> floors FULLY corrected, with bootstrap errors.")
F = np.load(R("ch_floors.npz"), allow_pickle=True)
A = np.load(R("ch_analysis.npz"), allow_pickle=True)
cells = [str(c) for c in F["cells"]]
tags = [str(t) for t in A["surface_tags"]]
ccorr, cwrong = A["surface_corr"], A["surface_wrong"]

# ch_floors names a cell "1,e03,lit"; ch_analysis names the same cell "m1e03_lit". Same order,
# different spelling. Align by CONSTRUCTION and assert it, so a future reorder cannot silently
# marry a floor to the wrong count.
def as_tag(c):
    n, e, tier = c.split(",")
    return f"m{n}{e}_{tier}"


assert [as_tag(c) for c in cells] == tags, "ch_floors / ch_analysis cell order disagrees"

rows = []
for i, c in enumerate(cells):
    off = F[f"offN_{i}"]
    zfc = float(F[f"zero_frac_{i}"])
    gu, gsd = float(F[f"fN_{i}"]), float(F[f"fN_sd_{i}"])
    fl, er, est = adopt(gu, gsd, off, zfc)
    rows.append((c, zfc, gu, gsd, fl, er, est, float(ccorr[i]), float(cwrong[i])))

n_inv = sum(1 for r in rows if r[6] != "gumbel")
n_up = sum(1 for r in rows if r[4] > r[2])
print(f"cells with zero-fraction > {ZF_GATE:.0%} -> Gumbel INVALID : {n_inv} / {len(rows)}")
print(f"   of which the floor RISES (count can only fall)         : {n_up}")
print()
print(f"{'cell':>12} | {'zero-f':>6} | {'Gumbel':>7} {'+-fit':>6} | {'ADOPTED':>8} {'+-bs':>6} "
      f"{'est':>8} | {'ratio':>6} | {'count':>6} {'wrong':>6}")
for c, zfc, gu, gsd, fl, er, est, cc, cw in rows:
    print(f"{c:>12} | {zfc:6.2f} | {gu:7.2f} {gsd:6.2f} | {fl:8.2f} {er:6.2f} {est:>8} "
          f"| {fl/gu:5.2f}x | {cc:6.2f} {cw:6.2f}")
print()
print("THE SWITCH-ON CELLS — the number quoted externally:")
for c, zfc, gu, gsd, fl, er, est, cc, cw in rows:
    if c.startswith("1,e03") or c.startswith("1,e05"):
        marg = cc - 1.0
        print(f"   {c:>12}: floor {gu:.2f} -> {fl:.2f} nat ({fl/gu:.2f}x, "
              f"{abs(fl-gu)/gsd:.1f} sigma of its own Gumbel fit error); banked count {cc:.2f}, "
              f"margin over the >1 bar = {marg:+.2f}")

out.update(ch_cell=np.array([r[0] for r in rows]), ch_zf=np.array([r[1] for r in rows]),
           ch_gumbel=np.array([r[2] for r in rows]), ch_gumbel_sd=np.array([r[3] for r in rows]),
           ch_adopted=np.array([r[4] for r in rows]), ch_adopted_sd=np.array([r[5] for r in rows]),
           ch_estimator=np.array([r[6] for r in rows]), ch_corr_oldfloor=np.array([r[7] for r in rows]),
           ch_wrong_oldfloor=np.array([r[8] for r in rows]), ch_n_invalid=n_inv, ch_n_up=n_up)

# ==================================================================== IGNITE-2
hr("IGNITE-2 — the programme's TWO calibrated cells (STORY S6.2.1), under the convention.")
G = np.load(R("ig2_floors.npz"), allow_pickle=True)
ig = []
for i, c in enumerate([str(x) for x in G["cells"]]):
    off = G[f"offN_{i}"]
    zfc = float(G[f"zero_frac_{i}"])
    gu, gsd = float(G[f"fN_{i}"]), float(G[f"fN_sd_{i}"])
    fl, er, est = adopt(gu, gsd, off, zfc)
    ig.append((c, zfc, gu, gsd, fl, er, est))
    verdict = ("GUMBEL VALID — floor and count stand as banked" if est == "gumbel"
               else "GUMBEL INVALID — floor must be RESTATED as the empirical quantile")
    print(f"   {c:>16}  zero-frac {zfc:.2f} | Gumbel {gu:6.2f} +-{gsd:.2f} | "
          f"EMP q95 {emp_quantile(off):6.2f} +-{boot_sd(off):.2f} | {verdict}")
print()
print("   => ANCHOR §4's 'IGNITE-2's two calibrated cells are safe' needs one word changed.")
print("      (-12.75, 30, lit) is safe because its zero-fraction is 0.00 — the estimator holds.")
print("      (-13.25, 30, vlbi) has zero-fraction 0.45: its Gumbel is INVALID under the very")
print("      convention ANCHOR proposes. Its count (0.54) is safe only because the corrected")
print("      floor FALLS (count can only rise) — not because the estimator is sound there.")
print("      The floor VALUE must be restated as the empirical quantile + bootstrap error.")

out.update(ig2_cell=np.array([r[0] for r in ig]), ig2_zf=np.array([r[1] for r in ig]),
           ig2_gumbel=np.array([r[2] for r in ig]), ig2_adopted=np.array([r[4] for r in ig]),
           ig2_adopted_sd=np.array([r[5] for r in ig]),
           ig2_estimator=np.array([r[6] for r in ig]))

# ==================================================================== ANCHOR
hr("ANCHOR — the ladder re-derived from its raw offenders. A BANK DEFECT, found and fixed here.")
L = np.load(R("anchor_ladder.npz"), allow_pickle=True)
offs, zfl, empb = L["offenders"], L["zero_frac"], L["emp_q95"]
rung, hh, tt = L["rung"], L["h"], L["tier"]
NR, NC = 6, 8                       # 6 rungs (R0,R1,R2,R3,R2o,R3o) x 8 cells (4 h x 2 tier)

# THE DEFECT. `offenders` is stored (rung-major, cell-minor) -- offenders[rung*NC + cell] --
# while EVERY sibling column (rung, h, tier, zero_frac, emp_q95, ...) is (cell-major,
# rung-minor) -- row j = cell*NR + rung. The array is the TRANSPOSE of its own metadata, and
# nothing in the npz says so. ANCHOR §9 advertises "every floor is re-cuttable from the bank
# without rerunning a GPU" -- it is, but a naive re-cut marries offenders[j] to the WRONG
# cell's label and the error is SILENT (both are plausible floors). Caught only by re-deriving
# a banked number from its own raw column and finding it did not match.
PERM = np.array([(j % NR) * NC + (j // NR) for j in range(NR * NC)])

rec = np.array([emp_quantile(offs[PERM[j]]) for j in range(NR * NC)])
zrec = np.array([float((np.asarray(offs[PERM[j]]) == 0).mean()) for j in range(NR * NC)])
dev_q, dev_z = np.abs(rec - empb).max(), np.abs(zrec - zfl).max()

naive = np.array([emp_quantile(offs[j]) for j in range(NR * NC)])
print(f"   48 (rung, cell) rows, 150 offenders each.")
print(f"   NAIVE re-cut (offenders[j] labelled by row j):  max |dev| vs banked emp_q95 = "
      f"{np.abs(naive - empb).max():8.3f} nat   <- WRONG, and silently so")
print(f"   TRANSPOSED re-cut (offenders[PERM[j]]):        max |dev| vs banked emp_q95 = "
      f"{dev_q:8.3e} nat   <- EXACT")
print(f"                                                  max |dev| vs banked zero_frac = "
      f"{dev_z:8.3e}")
assert dev_q == 0.0 and dev_z == 0.0, "ANCHOR ladder does not re-derive from its own offenders"
print()
print("   => ANCHOR's PUBLISHED §3 ladder table is CORRECT. The defect is in the npz layout")
print("      only: `offenders` is the transpose of its own metadata, undeclared. FIX ON ACCRE:")
print("      re-bank `offenders` in metadata-row order, or add an explicit index key. Until")
print("      then, any re-cut from this file MUST apply the permutation asserted above.")
print()
print(f"   rows with zero-fraction > {ZF_GATE:.0%} (Gumbel invalid; verdict issued on EMP q95): "
      f"{int((zfl > ZF_GATE).sum())} / 48")
out.update(an_zf=zfl, an_emp_recomputed=rec, an_emp_banked=empb, an_perm=PERM,
           an_max_dev_transposed=float(dev_q), an_max_dev_naive=float(np.abs(naive - empb).max()),
           an_rung=rung, an_h=hh, an_tier=tt)

# ==================================================================== bank
np.savez(R("floorfix_provisional.npz"), alpha=ALPHA, zf_gate=ZF_GATE, boot=BOOT, seed=SEED, **out)
hr("BANKED -> reports/floorfix_provisional.npz")
print("PROVISIONAL. The floors are final; the COUNTS are not re-cut and no count above is")
print("corrected. Onset flips, the corrected onset table, and the CHORUS e=0.3 verdict all")
print("require the per-realisation signal banks on ACCRE. Nothing here may be quoted as a")
print("corrected count.")
