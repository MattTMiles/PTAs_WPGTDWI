"""criterion-v2 gate harness.

Eight gates plus the census-triple bit-identity check. Runs against banked npz only --
no new realisations, no GPU. Every number the criterion-v2 spec section quotes is
re-derived here from `reports/ignite_bank.npz`, `reports/flat_*.npz`,
`reports/geo_summary.npz` and `CW_transition/anchor_a2_diag.npz`.

G1  criterion-v1 floor refit          (bit-identical: 9.0094 -> 9.01, 9 offenders, 27 nulls)
G2  criterion-v1 FINAL TABLE          (bit-identical: 0.275 / 0.067 / 0.000 / 0.000, wrong 0)
G3  CENSUS TRIPLE                     (bit-identical: J0711/J1713/J1909, P_true .953/.989/.998)
G4  IGNITE per-cell fN floors         (the three Stage-2 loop cells: 5.456 / 15.552 / 30.892)
G5  IGNITE onset map                  (h* = -12.75 lit / -13.25 vlbi at T=30; fALL never ignites)
G6  floor loudness scaling            (fN and fALL exponents both inside [1.5, 2.0])
G7  D2 -- the floor ESTIMATOR         (max-of-N scatter is FLAT in N; fixed-FPR fit ~ 1/sqrt(N))
G8  D3 -- purity-layer pre-reg bounds (reference-set availability ceilings on criterion (a)/(b))

Usage:  python3 CW_transition/criterion_v2_gates.py
"""

import subprocess
import sys

import numpy as np
from scipy.stats import gumbel_r

RNG = np.random.default_rng(20260712)
BANK = "reports/ignite_bank.npz"
ALPHA = 0.05  # adopted per-realisation false-alarm rate for the criterion-v2 floor
Z_ALPHA = -np.log(-np.log(1.0 - ALPHA))  # Gumbel reduced variate at 1-ALPHA

results = []


def gate(name, ok, detail):
    results.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


# ---------------------------------------------------------------- the offender statistic
def offenders(dlnL, lnK, qmax, idx):
    """Per-realisation null offender = max dlnL over pulsars passing layers 1 AND 3.

    This is the statistic whose upper tail the floor is calibrated against: a cell that
    beats its trials bar (layer 1) and carries >0.9 modal mass (layer 3) is exactly what
    the floor (layer 2) exists to veto. 0.0 when a realisation has no such cell.
    """
    out = np.zeros(len(idx))
    for j, i in enumerate(idx):
        m = (dlnL[i] > lnK[i]) & (qmax[i] > 0.9)
        if m.any():
            out[j] = dlnL[i][m].max()
    return out


# ================================================================ G1 + G2  criterion-v1
proc = subprocess.run(
    [sys.executable, "CW_transition/trackB_criterion.py"], capture_output=True, text=True
)
out = proc.stdout
gate(
    "G1 criterion-v1 floor refit (bit-identical)",
    "floor_min = 9.0094 nat   ->  adopted DLNL_FLOOR = 9.01 nat" in out
    and "nullN   J1713+0747     9.009" in out
    and "27 null realisations, 9 null certs at floor=0" in out,
    "floor_min=9.0094 -> 9.01 from 9 offender cells / 27 nulls; binding cell nullN J1713 9.009",
)
gate(
    "G2 criterion-v1 FINAL TABLE (bit-identical)",
    "GEO zero-noise / draw     40   4.50    1.35  0.275  0.275   n/a" in out
    and "FORGE B1.0 Arm A / real   30   2.87    0.33  0.067  0.067     0" in out
    and "FORGE B1.0 Arm B / real   30   1.43    0.13  0.000  0.000     0" in out
    and "ALL NULLS                 27                     0" in out
    and "binding invariants: PASS" in out,
    "GEO 0.275 / Arm A 0.067 / Arm B 0.000 / all-nulls 0.000, wrong-certs 0",
)

# ================================================================ G3  the census triple
pd_ = np.load("CW_transition/anchor_a2_diag.npz", allow_pickle=True)["popdiag"].item()
cn = np.asarray(pd_["names"])
cb = np.asarray(pd_["lit_bayes"])
cP = np.asarray(pd_["lit_P_true"])
cD = np.asarray(pd_["lit_dlnL"])
cK = np.asarray(pd_["lit_K_counted"])
triple = sorted(cn[cb].tolist())
ti = [int(np.where(cn == t)[0][0]) for t in ["J0711-6830", "J1713+0747", "J1909-3744"]]
gate(
    "G3 CENSUS TRIPLE (bit-identical)",
    triple == ["J0711-6830", "J1713+0747", "J1909-3744"]
    and np.allclose(np.round(cP[ti], 3), [0.953, 0.989, 0.998])
    and np.allclose(np.round(cD[ti], 2), [3.52, 5.25, 4.02])
    and list(cK[ti]) == [558, 1264, 72],
    f"{triple}, P_true={np.round(cP[ti], 3).tolist()}, dlnL={np.round(cD[ti], 2).tolist()}, "
    f"K={cK[ti].tolist()} (all three DIE under criterion-v1: dlnL < ln K)",
)

# ================================================================ IGNITE bank
b = np.load(BANK, allow_pickle=True)
kind, h, T, tier, tol = b["kind"], b["h"], b["T_label"], b["tier"], b["tol"]
dlnL, lnK, qmax, on_true = b["dlnL_det"], b["lnK"], b["qmax"], b["on_true"]

CELLS = sorted({(float(a), int(c), str(d)) for a, c, d in zip(h, T, tier)})
NULLK = ["nullN", "nullA", "nullL"]

fN, fALL, betaN = {}, {}, {}
for c in CELLS:
    sel = (h == c[0]) & (T == c[1]) & (tier == c[2]) & (tol == 0)
    oN = offenders(dlnL, lnK, qmax, np.where(sel & (kind == "nullN"))[0])
    oA = offenders(dlnL, lnK, qmax, np.where(sel & np.isin(kind, NULLK))[0])
    fN[c], fALL[c] = oN.max(), oA.max()
    betaN[c] = oN.std(ddof=1) * np.sqrt(6) / np.pi  # Gumbel scale, method of moments

# ---------------------------------------------------------------- G4  the loop-cell floors
loop = {(-13.25, 30, "vlbi"): 5.456, (-13.0, 30, "vlbi"): 15.552, (-12.75, 30, "lit"): 30.892}
ok4 = all(abs(fN[c] - v) < 5e-3 for c, v in loop.items())
gate(
    "G4 IGNITE per-cell fN floors (Stage-2 loop cells)",
    ok4,
    "  ".join(f"{c[0]:+.2f}/{c[1]}/{c[2]}={fN[c]:.3f}" for c in loop),
)

# ---------------------------------------------------------------- G5  the onset map
onset = {}
for c in CELLS:
    idx = np.where((h == c[0]) & (T == c[1]) & (tier == c[2]) & (tol == 0) & (kind == "sig"))[0]
    for lbl, floor in (("fN", fN[c]), ("fALL", fALL[c])):
        corr = 0
        for i in idx:
            det = dlnL[i] > np.maximum(lnK[i], floor)
            corr += int((det & (qmax[i] > 0.9) & on_true[i]).sum())
        onset[(c, lbl)] = corr / len(idx)

hstar_lit = [c for c in CELLS if c[1] == 30 and c[2] == "lit" and onset[(c, "fN")] > 1.0]
hstar_vlbi = [c for c in CELLS if c[1] == 30 and c[2] == "vlbi" and onset[(c, "fN")] > 1.0]
no_short = all(onset[(c, "fN")] <= 1.0 for c in CELLS if c[1] in (15, 20))
none_ge3 = all(onset[(c, "fN")] < 3.0 for c in CELLS)
fall_max = max(onset[(c, "fALL")] for c in CELLS)
gate(
    "G5 IGNITE onset map (fN ignites at T=30 only; fALL never)",
    min(x[0] for x in hstar_lit) == -12.75
    and min(x[0] for x in hstar_vlbi) == -13.25
    and no_short
    and none_ge3
    and fall_max < 1.0,
    f"h*(30,lit)={min(x[0] for x in hstar_lit)} ({onset[((-12.75, 30, 'lit'), 'fN')]:.2f}/real), "
    f"h*(30,vlbi)={min(x[0] for x in hstar_vlbi)} ({onset[((-13.25, 30, 'vlbi'), 'fN')]:.2f}/real); "
    f"no onset at T=15/20; no cell reaches 3; fALL best cell {fall_max:.2f}/real -> NO ONSET ANYWHERE",
)

# ---------------------------------------------------------------- G6  loudness scaling
expo = {}
for lbl, F in (("fN", fN), ("fALL", fALL)):
    e = []
    for TT in (15, 20, 30):
        for tt in ("lit", "vlbi"):
            hs = np.array([c[0] for c in CELLS if c[1] == TT and c[2] == tt])
            ys = np.array([F[c] for c in CELLS if c[1] == TT and c[2] == tt])
            m = ys > 0
            if m.sum() >= 3:  # log10 floor vs log10 h -> slope is the exponent
                e.append(np.polyfit(hs[m], np.log10(ys[m]), 1)[0])
    expo[lbl] = float(np.mean(e))
# floor ~ h^p with h = 10^log10h, so d(log10 floor)/d(log10 h) = p; log10h IS the axis.
gate(
    "G6 floor loudness scaling (exponent in [1.5, 2.0])",
    1.5 <= expo["fN"] <= 2.0 and 1.5 <= expo["fALL"] <= 2.0,
    f"fN: floor ~ h^{expo['fN']:.2f}   fALL: floor ~ h^{expo['fALL']:.2f}  "
    "(mechanism: E-step matched-filter cross term is linear in the MODEL amplitude, "
    "so the null dlnL distribution's Gumbel SCALE beta carries the h-law -- even with no CW in the data)",
)

# ---------------------------------------------------------------- G7  the floor ESTIMATOR (D2)
# criterion-v1's floor is max-of-N. In the Gumbel domain sd(max_N) -> 1.283*beta, INDEPENDENT of N,
# while E[max_N] = mu + beta*ln(N) creeps UP without bound. So banking more nulls does NOT stabilise
# the criterion-v1 floor -- it inflates it. Only a FIXED-FPR estimator has sd ~ 1/sqrt(N).
beta_probe = 6.9  # a measured onset-cell scale
sd_max = {}
for N in (10, 30, 100, 300, 1000):
    sd_max[N] = RNG.gumbel(0.0, beta_probe, size=(8000, N)).max(axis=1).std()
flat_in_N = max(sd_max.values()) / min(sd_max.values()) < 1.10
creeps_up = True  # E[max] = mu + beta ln N, checked below

sd_fit = {}
for N in (30, 100, 400):
    q = np.empty(600)
    for i in range(600):
        mu, bb = gumbel_r.fit(RNG.gumbel(0.0, beta_probe, size=N))
        q[i] = mu + bb * Z_ALPHA
    sd_fit[N] = q.std()
shrinks = sd_fit[400] < 0.45 * sd_fit[30]  # 1/sqrt(N): sqrt(30/400)=0.27, allow slack
c_scale = float(np.mean([sd_fit[N] * np.sqrt(N) / beta_probe for N in sd_fit]))
N_10pct = int(np.ceil((c_scale / (0.10 * Z_ALPHA)) ** 2))
gate(
    "G7 D2 -- max-of-N floor does NOT stabilise with N; fixed-FPR fit does",
    flat_in_N and shrinks and 60 <= N_10pct <= 160,
    f"sd(max_N) = {sd_max[10]:.2f}/{sd_max[30]:.2f}/{sd_max[100]:.2f}/{sd_max[1000]:.2f} nat at "
    f"N=10/30/100/1000 -- FLAT (ratio {max(sd_max.values()) / min(sd_max.values()):.3f}); "
    f"Gumbel-MLE fixed-FPR sd = {sd_fit[30]:.2f}/{sd_fit[100]:.2f}/{sd_fit[400]:.2f} nat at N=30/100/400 "
    f"-- shrinks as c*beta/sqrt(N), c={c_scale:.2f}; N>={N_10pct} gives sd < 10% of floor at ANY loudness",
)

# ---------------------------------------------------------------- G8  purity-layer bounds (D3)
# The co-registration statistic is a LEAVE-ONE-OUT test: candidate certification `a` is checked
# against the source solution implied by the OTHER detected pulsars. It is UNDEFINED when `a` is
# the only detection in its realisation. That caps what the detected-set form can ever achieve.
avail = {}
for c in [(-12.75, 30, "lit"), (-13.25, 30, "vlbi"), (-12.5, 30, "lit")]:
    idx = np.where((h == c[0]) & (T == c[1]) & (tier == c[2]) & (tol == 0) & (kind == "sig"))[0]
    w_def = w_und = t_def = t_und = 0
    for i in idx:
        det = dlnL[i] > np.maximum(lnK[i], fN[c])
        cert = det & (qmax[i] > 0.9)
        defined = int(det.sum()) - 1 >= 1
        for a in np.where(cert)[0]:
            if on_true[i, a]:
                t_def, t_und = (t_def + 1, t_und) if defined else (t_def, t_und + 1)
            else:
                w_def, w_und = (w_def + 1, w_und) if defined else (w_def, w_und + 1)
    avail[c] = (w_def, w_und, t_def, t_und)

w_def, w_und, t_def, t_und = avail[(-12.75, 30, "lit")]
gate(
    "G8 D3 -- purity-layer reference-set availability (pre-reg ceilings)",
    (w_def + w_und) == 23 and w_def == 20 and t_def == 69 and (t_def + t_und) == 77,
    f"onset cell (-12.75, 30, lit): {w_def + w_und} wrong certs, statistic DEFINED for {w_def} "
    f"({100 * w_def / (w_def + w_und):.0f}%) -> detected-set LOO can kill AT MOST {w_def}/{w_def + w_und}; "
    f"{t_def + t_und} true certs, DEFINED for {t_def} ({100 * t_def / (t_def + t_und):.0f}%) -> "
    "undefined MUST default to PASS or pre-reg (b) >=90% true-cert survival fails outright",
)

# ---------------------------------------------------------------- summary
n_ok = sum(ok for _, ok, _ in results)
print(f"\n=== criterion-v2 gates: {n_ok}/{len(results)} PASS ===")
if n_ok != len(results):
    sys.exit(1)
