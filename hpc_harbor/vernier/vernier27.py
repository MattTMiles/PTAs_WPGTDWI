#!/usr/bin/env python3
"""
VERNIER-27 -- FIRM THE K* NUMBER. The certification pass.
EXPLORATORY, NOT pre-registered. Toy, cronus, CPU/numpy. Hours, not minutes.

Clean room: imports only vernier27_lib / vernier25 / vernier24 (and through them
vernier_toy / vernier6 / vernier7 / vernier16). None of those files is edited;
the sky construction, confusion budget, F-statistic, certification statistic and
its incoherent-comb floor are all called out of them unmodified.

WHY
  VERNIER-25's headline -- K* ~ 20-25 recovered sources -- is an extrapolation
  off a sampling-starved tail (3-4 realisations per cell) whose last two points
  are, it turns out, THE SAME POINT. The paper's title claim rests on it, so it
  is bought properly here: a dense K grid, 20 realisations per K, several chains
  per realisation, cluster-bootstrap error bars, and the crossing fraction
  MEASURED rather than read off a fitted line.

DEFECT IN VERNIER-25 ARM K, FOUND AND FIXED HERE
  V25's chain() seeded bin b iff b < n_seed, and build_bin prepended exactly ONE
  seed. BAND holds TEN bins. So every K > 10 placed ten seeds:

      K   =  2  3  5  8  12  16
      seeds  2  3  5  8  10  10        <-- K=12 and K=16 are the SAME cell

  The V25 Arm K rows at K = 12 and K = 16 differ only by chain seed, and the
  contrast-vs-combs curve that gave K* ~ 20-25 was extrapolated off a tail that
  had stopped moving for a construction reason. Every V25 result at K <= 10 is
  unaffected. Fixed here -- and the fix is not unique, which is the whole story:

TWO CO-PRIMARY ARMS, BECAUSE K > 10 HAS NO SINGLE HONEST CONSTRUCTION
  With T = 20 yr the band 10-25 nHz holds exactly TEN 1/T bins. K = 32 sources
  at distinct frequencies inside that band therefore MUST share bins.

  ARM P (PACKED)  -- the brief read literally. K seeds spread over the ten bins
      of 10-25 nHz, several per bin, each on its own sub-bin frequency. Sources
      sharing a bin blend in the F-statistic, so this is the crowded, pessimistic
      construction: the cost of finding K sources in one crowded decade.
  ARM W (WIDE)    -- V25 Arm K's own construction continued. One seed per 1/T
      bin, band widened to K bins and centred on the sensitive decade. No two
      seeds share a bin, so this is the optimistic construction: the cost of
      finding K sources spread across the band.

  The truth is between them and they are reported as a BRACKET, not averaged.
  A single number here would be a construction choice wearing a physics hat.

CHANGES vs V25 (all in vernier27_lib.py, all named there)
  1. seed packing with per-source sub-bin frequencies (Arm P)
  2. noise/sky split, so "several chains per realisation" means what it says
  3. COMB_CAP 8 -> 32, per V25's own caveat 1; priced in ARM C
  4. the GENUINE pool made K-invariant (seeds + top-10 sky, not a global top-10)
  5. the detection-floor RNG isolated from the data RNG -- V25's memoised floor
     made a chain's data depend on how many chains had run before it

THE DELIVERABLE
  K*_50, K*_90 with bootstrap error bars, per arm; the crossing-fraction curve;
  and the local slope d(contrast)/d(comb) at each K, which is the
  extrapolation-credibility check. If contrast saturates below the bar the
  SKA-era claim dies and this is where it dies.

SAFETY KIT, every cell
  M = 1 control; pure-noise null through the whole chain; per-turn source-level
  scrambled null (which vetoes the deepening flag by construction); V18
  genuine/split/blend accounting; wrong-certification ledger.

Code hpc_harbor/vernier/vernier27{,_lib}.py . data reports/vernier27.npz
Staged, not committed.
"""
import os
import time
import numpy as np

import vernier27_lib as L
import vernier25 as V25
import vernier24 as W

REPORTS = W.REPORTS

K_LIST    = (12, 16, 20, 24, 28, 32)
N_REAL    = 20
N_CHAIN_P = 8
N_CHAIN_W = 6
N_BOOT    = 4000
RNG_BOOT  = 20270806

ARMS = (("P", "PACKED  -- K seeds in the ten bins of 10-25 nHz, several per bin",
         None, N_CHAIN_P),
        ("W", "WIDE    -- one seed per 1/T bin, band widened to K bins",
         L.wide_band, N_CHAIN_W))


def max_turn_for(K):
    """Room to work K seeds. Checked: chains stop on their own well inside this."""
    return 2 * int(K) + 16


# ==================================================================== plumbing
def run_cell(K, arm="P", n_real=N_REAL, n_chain=N_CHAIN_P, comb_cap=L.COMB_CAP,
             gamma=1.56, tag=""):
    """
    One (arm, K) cell. n_real independent universes (array + sky); n_chain chains
    per universe differing ONLY in the data realisation and the search's own
    randomness. Every array below is (n_real, n_chain).
    """
    band = L.BAND if arm == "P" else L.wide_band(K)
    mt = max_turn_for(K)
    keys = ("S", "BAR", "NC", "GEN", "SPL", "BLD", "CER", "TRU", "WRG", "NUL",
            "TRN", "HIT", "IGN", "DPN")
    A = {k: np.zeros((n_real, n_chain)) for k in keys}
    AEQ = np.zeros(n_real)
    t0 = time.time()
    for r in range(n_real):
        arr = W.make_array(np.random.default_rng(2500 + r))
        sky = 27000 + 977 * r + 13 * K + (0 if arm == "P" else 500000)
        for c in range(n_chain):
            res = L.chain(sky, arr, gamma, band=band, n_seed=K, max_turn=mt,
                          comb_cap=comb_cap,
                          noise_seed=810000 + 7919 * r + 101 * c + K)
            hh = [h for h in res["hist"] if h["which"] >= 0]
            A["S"][r, c] = max([h["S_max"] for h in hh], default=0.0)
            A["BAR"][r, c] = (float(np.median([h["lnK"] for h in hh]))
                              if hh else np.nan)
            A["NC"][r, c] = min(len(hh), comb_cap)
            A["GEN"][r, c] = res["genuine"]; A["SPL"][r, c] = res["split"]
            A["BLD"][r, c] = res["blend"];   A["TRN"][r, c] = res["turns"]
            A["CER"][r, c] = res["n_cert_tot"]; A["TRU"][r, c] = res["n_true_tot"]
            A["WRG"][r, c] = res["n_wrong_tot"]
            A["NUL"][r, c] = sum(h.get("n_cert_null", 0) for h in hh)
            A["HIT"][r, c] = res["seed_hits"]; A["IGN"][r, c] = float(res["ignite"])
            A["DPN"][r, c] = res["deepen"]
            if c == 0:
                AEQ[r] = res["a_eq"]
        el = time.time() - t0
        print(f"    [{tag}{arm} K={K:2d}] real {r+1:2d}/{n_real}  "
              f"S_end med {np.median(A['S'][r]):7.3f}  bar {np.nanmedian(A['BAR'][r]):5.2f}"
              f"  combs {np.median(A['NC'][r]):4.1f}  cert {A['CER'][r].sum():4.0f}"
              f"  ({el:6.1f}s, eta {el/(r+1)*(n_real-r-1):6.1f}s)", flush=True)
    A["AEQ"] = AEQ
    A["bar_med"] = float(np.nanmedian(A["BAR"]))
    A["K"] = K
    A["arm"] = arm
    return A


# ---------------------------------------------------------------- uncertainty
def _clusters(rng, n):
    return rng.integers(0, n, n)


def boot_median(X, n_boot=N_BOOT, seed=RNG_BOOT):
    """
    CLUSTER bootstrap of the median of X over realisations (rows). Chains within
    a realisation are correlated -- they share a universe -- so the resampling
    unit is the row, not the chain.
    """
    rng = np.random.default_rng(seed)
    X = np.asarray(X)
    n = X.shape[0]
    med = float(np.median(X))
    bs = np.array([float(np.median(X[_clusters(rng, n)])) for _ in range(n_boot)])
    return med, float(bs.std(ddof=1))


def boot_frac(flag, n_boot=N_BOOT, seed=RNG_BOOT + 1):
    rng = np.random.default_rng(seed)
    f = np.asarray(flag, float)
    n = f.size
    bs = np.array([f[_clusters(rng, n)].mean() for _ in range(n_boot)])
    return float(f.mean()), float(bs.std(ddof=1))


def kstar(Ks, flags, thresh, n_boot=N_BOOT, seed=RNG_BOOT + 2):
    """
    Smallest K reaching crossing fraction `thresh`, two ways:
      grid   -- smallest GRID K that reaches it (what the brief asks for)
      interp -- linear interpolation between grid points, which is the thing that
                can carry an error bar
    Bootstrap resamples realisations independently within each K.
    Returns (grid_K|None, interp, se, lo, hi, never_frac).
    """
    rng = np.random.default_rng(seed)
    Ks = np.asarray(Ks, float)

    def _first(p):
        hit = np.where(p >= thresh)[0]
        if hit.size == 0:
            return None, np.nan
        i = int(hit[0])
        if i == 0:
            return Ks[0], float(Ks[0])
        p0, p1 = p[i - 1], p[i]
        if p1 <= p0:
            return Ks[i], float(Ks[i])
        w = (thresh - p0) / (p1 - p0)
        return Ks[i], float(Ks[i - 1] + w * (Ks[i] - Ks[i - 1]))

    g_obs, i_obs = _first(np.array([f.mean() for f in flags]))
    vals, never = [], 0
    for _ in range(n_boot):
        p = np.array([f[_clusters(rng, f.size)].mean() for f in flags])
        _, v = _first(p)
        if np.isnan(v):
            never += 1
        else:
            vals.append(v)
    if len(vals) > 2:
        se = float(np.std(vals, ddof=1))
        lo, hi = (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))
    else:
        se = lo = hi = np.nan
    return g_obs, i_obs, se, lo, hi, never / n_boot


# ======================================================================= arms
def run_arm(out, arm, desc, n_chain):
    print()
    print("=" * 120)
    print(f"DELIVERABLE 1{arm} -- ARM {arm}: {desc}")
    print("=" * 120)
    print(f"  {N_REAL} realisations x {n_chain} chains per K, K in {K_LIST}, "
          f"gamma = 1.56, COMB_CAP = {L.COMB_CAP}")
    if arm == "W":
        for K in K_LIST:
            b = L.wide_band(K)
            print(f"    K = {K:2d}: {len(b)} bins, {b[0]*1e9:5.2f}-{b[-1]*1e9:5.2f} nHz")
    print()
    cells = [run_cell(K, arm=arm, n_chain=n_chain) for K in K_LIST]
    out[f"cells_{arm}"] = cells
    _bank(out)
    return cells


def report_grid(out, arm, cells):
    print()
    print("-" * 120)
    print(f"  ARM {arm} PER-K SUMMARY.  S_end = end-state fringe contrast (max over")
    print("  turns of a per-pulsar max); bar = median trials bar ln K over the turns.")
    print("  x_real = fraction of REALISATIONS whose BEST chain reaches the bar.")
    print("-" * 120)
    print(f"  {'K':>3} {'A_eq':>6} {'K_rec':>6} {'combs':>6} {'S_end':>8} {'+-SE':>6} "
          f"{'IQR':>15} {'bar':>5} {'S/bar':>6} {'x_real':>7} {'+-SE':>6} "
          f"{'x_chain':>8} {'CERT/run':>9} {'true':>6} {'WRONG':>6} {'null':>6}")
    rows = []
    for c in cells:
        med, se = boot_median(c["S"])
        q1, q3 = np.percentile(c["S"], [25, 75])
        bar = c["bar_med"]
        # crossing is judged chain by chain against THAT chain's own trials bar,
        # not against the cell median; bar_med is displayed only
        xr, xr_se = boot_frac(_cross(c).max(axis=1))
        xc = float(_cross(c).mean())
        rows.append([c["K"], float(np.median(c["AEQ"])), float(np.median(c["GEN"])),
                     float(np.median(c["NC"])), med, se, q1, q3, bar,
                     med / max(bar, 1e-9), xr, xr_se, xc, float(c["CER"].mean()),
                     float(c["TRU"].mean()), float(c["WRG"].mean()),
                     float(c["NUL"].mean())])
        print(f"  {c['K']:3d} {np.median(c['AEQ']):6.2f} {np.median(c['GEN']):6.1f} "
              f"{np.median(c['NC']):6.1f} {med:8.3f} {se:6.3f} "
              f"{q1:7.2f}-{q3:7.2f} {bar:5.2f} {med/max(bar,1e-9):6.3f} "
              f"{xr:7.2f} {xr_se:6.2f} {xc:8.2f} {c['CER'].mean():9.2f} "
              f"{c['TRU'].mean():6.2f} {c['WRG'].mean():6.2f} {c['NUL'].mean():6.2f}")
    rows = np.array(rows)
    print()
    print("  SE ADEQUACY -- is SE(median S_end) below half the gap to the nearest")
    print("  neighbouring K median?  (the brief's sampling requirement)")
    ok = True
    for i, c in enumerate(cells):
        med, se = rows[i, 4], rows[i, 5]
        gaps = ([abs(med - rows[i - 1, 4])] if i > 0 else []) + \
               ([abs(med - rows[i + 1, 4])] if i < len(cells) - 1 else [])
        g = min(gaps)
        good = se < 0.5 * g
        ok &= bool(good)
        print(f"    K = {c['K']:2d}: SE {se:6.3f}  vs  half-gap {0.5*g:6.3f}   "
              f"{'OK' if good else 'NOT MET'}")
    print(f"  -> criterion {'MET in every cell' if ok else 'NOT met in every cell'}")
    out[f"grid_{arm}"] = rows
    out[f"se_ok_{arm}"] = np.array([float(ok)])
    return rows


def _cross(c):
    """Per-chain: did the end-state contrast reach THAT chain's own trials bar?"""
    return np.where(np.isnan(c["BAR"]), False, c["S"] >= c["BAR"])


CRIT = (("contrast", lambda c: _cross(c).max(axis=1),
         "best chain's end-state contrast reaches the trials bar"),
        ("certification", lambda c: c["CER"].max(axis=1) >= 1.0,
         "best chain actually CERTIFIES at least one pulsar distance"))


def report_kstar(out, arm, cells):
    print()
    print("=" * 120)
    print(f"DELIVERABLE 2{arm} -- ARM {arm}: K*_50 AND K*_90")
    print("=" * 120)
    Ks = [c["K"] for c in cells]
    for name, fn, desc in CRIT:
        flags = [fn(c) for c in cells]
        print()
        print(f"  criterion: {desc}")
        line = "    crossing fraction vs K: "
        for k, f in zip(Ks, flags):
            p, s = boot_frac(f)
            line += f"  K={k}: {p:.2f}+-{s:.2f}"
        print(line)
        for thr, lbl in ((0.5, "K*_50"), (0.9, "K*_90")):
            g, i_, se, lo, hi, nev = kstar(Ks, flags, thr)
            if g is None:
                print(f"    {lbl} ({name}): NOT REACHED at K <= {max(Ks)}  "
                      f"(unreached in {nev*100:.0f} % of bootstrap resamples)")
            else:
                print(f"    {lbl} ({name}): grid K = {int(g)};  interpolated "
                      f"{i_:.1f} +- {se:.1f}   95 % [{lo:.1f}, {hi:.1f}]"
                      + (f";  unreached in {nev*100:.0f} % of resamples"
                         if nev > 0.01 else ""))
            out[f"kstar_{arm}_{name}_{int(thr*100)}"] = np.array(
                [np.nan if g is None else g, i_, se, lo, hi, nev])
        out[f"xfrac_{arm}_{name}"] = np.array([f.mean() for f in flags])
        out[f"xfrac_{arm}_{name}_se"] = np.array([boot_frac(f)[1] for f in flags])
        # the framing flag the brief asks for
        early = [k for k, f in zip(Ks, flags) if f.mean() > 0.0 and k < 16]
        if early:
            print(f"    ** FLAG: nonzero crossing at K = {early} (< 16). "
                  f"That is MeerKAT-era territory, not SKA-era. **")


def report_slope(out, arm, cells, rows):
    print()
    print("=" * 120)
    print(f"DELIVERABLE 3{arm} -- ARM {arm}: SHAPE. IS CONTRAST STILL LINEAR IN COMBS?")
    print("=" * 120)
    print("  V25 read ~2 nat/comb off six points and extrapolated to the bar at 10-11")
    print("  combs. If the curve saturates BELOW the bar, the SKA-era claim is dead.")
    print()
    nc, sm, se = rows[:, 3], rows[:, 4], rows[:, 5]
    print(f"  {'K':>3} {'combs':>6} {'S_end':>8} {'+-SE':>6} {'local dS/dcomb':>15} "
          f"{'within-K dS/dcomb':>19} {'gap to bar':>11}")
    slopes = []
    for i, c in enumerate(cells):
        # d(contrast)/d(comb) across neighbouring grid points. UNDEFINED when the
        # comb count does not move between them -- which is itself the finding in
        # the packed arm, where combs are pinned near 8 at every K, so guard it
        # rather than dividing by ~0 and printing a number with nine digits.
        lo, hi = max(i - 1, 0), min(i + 1, len(cells) - 1)
        dn = nc[hi] - nc[lo]
        d = (sm[hi] - sm[lo]) / dn if abs(dn) > 0.25 else np.nan
        x, y = c["NC"].reshape(-1), c["S"].reshape(-1)
        w = float(np.polyfit(x, y, 1)[0]) if x.std() > 1e-9 else np.nan
        gap = rows[i, 8] - sm[i]
        slopes.append([c["K"], nc[i], sm[i], se[i], d, w, gap])
        print(f"  {c['K']:3d} {nc[i]:6.1f} {sm[i]:8.3f} {se[i]:6.3f} {d:15.3f} "
              f"{w:19.3f} {gap:+11.3f}")
    slopes = np.array(slopes)
    lo_s = float(np.nanmean(slopes[:2, 4])); hi_s = float(np.nanmean(slopes[-2:, 4]))
    bar, s_hi = rows[-1, 8], sm[-1]
    print()
    print(f"  slope at low K ({K_LIST[0]}-{K_LIST[1]}): {lo_s:+.3f} nat/comb; "
          f"at high K ({K_LIST[-2]}-{K_LIST[-1]}): {hi_s:+.3f} nat/comb")
    if not np.isfinite(hi_s):
        v = ("UNDEFINED -- the comb count does not move with K, so there is no "
             "contrast-vs-comb slope to read. Combs are pinned; see dS/dK instead.")
    elif s_hi >= bar:
        v = "ALREADY AT THE BAR by K = %d -- no extrapolation needed." % K_LIST[-1]
    elif hi_s <= 0.0:
        v = ("SATURATED / TURNED OVER BELOW THE BAR. The extrapolation is not "
             "supported and the SKA-era claim does not follow from this curve.")
    elif hi_s < 0.5 * lo_s:
        v = ("SATURATING below the bar. A linear extrapolation OVERSTATES the "
             "reach; K* is further away than V25's read.")
    elif hi_s > 1.5 * lo_s:
        v = "ACCELERATING. K* is NEARER than a linear read would say."
    else:
        v = "~LINEAR. The extrapolation shape is supported."
    print(f"  VERDICT (arm {arm}): {v}")
    out[f"slopes_{arm}"] = slopes
    out[f"slope_verdict_{arm}"] = np.array([lo_s, hi_s, bar, s_hi])
    print(f"  slope_verdict_{arm} = [{lo_s:.4f}, {hi_s:.4f}, bar {bar:.3f}, "
          f"S_end(Kmax) {s_hi:.3f}]")


def armC_capsens(out, Ks=(16, 24, 32), n_real=10, n_chain=4):
    print()
    print("=" * 120)
    print("DELIVERABLE 4 -- ARM C: WHAT V25's COMB_CAP = 8 COST")
    print("=" * 120)
    print("  V25 caveat 1: 'the comb cap is mine and it binds at exactly the point")
    print("  that matters.' Same universes, same noise, cap 8 vs cap 32.")
    print()
    print(f"  {'arm':>4} {'K':>3} {'S_end cap8':>11} {'S_end cap32':>12} {'ratio':>7} "
          f"{'x_real cap8':>12} {'x_real cap32':>13} {'cert cap8':>10} "
          f"{'cert cap32':>11}")
    rows = []
    for arm in ("P", "W"):
        for K in Ks:
            c8 = run_cell(K, arm=arm, n_real=n_real, n_chain=n_chain,
                          comb_cap=8, tag="cap8 ")
            c32 = run_cell(K, arm=arm, n_real=n_real, n_chain=n_chain,
                           comb_cap=32, tag="cap32 ")
            m8, m32 = float(np.median(c8["S"])), float(np.median(c32["S"]))
            x8 = float(_cross(c8).max(axis=1).mean())
            x32 = float(_cross(c32).max(axis=1).mean())
            rows.append([0 if arm == "P" else 1, K, m8, m32, m32 / max(m8, 1e-9),
                         x8, x32, float(c8["CER"].mean()), float(c32["CER"].mean())])
            print(f"  {arm:>4} {K:3d} {m8:11.3f} {m32:12.3f} "
                  f"{m32/max(m8,1e-9):7.2f} {x8:12.2f} {x32:13.2f} "
                  f"{c8['CER'].mean():10.2f} {c32['CER'].mean():11.2f}")
            out["capsens"] = np.array(rows)
            _bank(out)


def controls(out, n_real=12, gamma=1.56):
    print()
    print("=" * 120)
    print("DELIVERABLE 5 -- THE SAFETY KIT")
    print("=" * 120)
    # ---- M = 1 control
    cs, ct, rx, sm = [], [], [], []
    dete = 0
    for s in range(n_real):
        arr = W.make_array(np.random.default_rng(2500 + s))
        r = L.chain(28100 + s, arr, gamma, band=L.BAND[:1], n_seed=1,
                    sea_keep=0.0, max_turn=32, noise_seed=880000 + s)
        dete += int(r["genuine"] > 0)
        for h in r["hist"]:
            if h["which"] >= 0:
                cs.append(h["n_cert"]); ct.append(h["n_true"])
                rx.append(h["rho_P_max"]); sm.append(h["S_max"])
    f_ = lambda v: (float(np.median(v)) if v else 0.0)
    print(f"  M = 1 control (single bin, seeded, sea removed): detected "
          f"{dete}/{n_real}; CERT {f_(cs):.2f} (true {f_(ct):.2f}); "
          f"rho_P max {f_(rx):.2f}; best fringe contrast {f_(sm):.3f} nat")
    print("    one comb cannot break its own degeneracy at any SNR (V24, V25).")

    # ---- pure-noise null through the whole chain, both arms, largest geometry
    nulls = []
    for arm in ("P", "W"):
        band = L.BAND if arm == "P" else L.wide_band(32)
        nf = cc = ig = 0
        for s in range(n_real):
            arr = W.make_array(np.random.default_rng(2500 + s))
            r = L.chain(28200 + s, arr, gamma, band=band, n_seed=32, null=True,
                        max_turn=max_turn_for(32), noise_seed=890000 + s)
            nf += int(r["turns"] > 0); cc += r["n_cert_tot"]; ig += int(r["ignite"])
        exp = 1.0 - (1.0 - 0.05) ** len(band)
        nulls.append([0 if arm == "P" else 1, len(band), nf, n_real, exp, cc, ig])
        print(f"  pure-noise null, arm {arm} (K = 32 geometry, {len(band)} bins, no "
              f"signal): fired {nf}/{n_real} (expected {exp:.2f} at alpha = 0.05 "
              f"per bin), {cc} certifications, {ig} ignitions")
    out["nulls"] = np.array(nulls, float)

    # ---- scrambled-null and wrong-cert ledger, both arms
    print()
    print("  per-turn SOURCE-LEVEL scrambled null, and the WRONG-CERT ledger:")
    print(f"    {'arm':>4} {'K':>3} {'cert/run':>9} {'true/run':>9} {'WRONG/run':>10} "
          f"{'wrong frac':>11} {'scrambled null/run':>19} {'null frac':>10}")
    led = []
    print("    CERT/run counts pulsars NEWLY certified; true/WRONG/null are per-turn")
    print("    sums, so the fractions are taken against true+WRONG, not against CERT.")
    for arm in ("P", "W"):
        for c in out[f"cells_{arm}"]:
            tot = float(c["TRU"].sum() + c["WRG"].sum())
            wf = float(c["WRG"].sum() / max(tot, 1.0))
            nf = float(c["NUL"].sum() / max(tot, 1.0))
            led.append([0 if arm == "P" else 1, c["K"], c["CER"].mean(),
                        c["TRU"].mean(), c["WRG"].mean(), wf, c["NUL"].mean(), nf])
            print(f"    {arm:>4} {c['K']:3d} {c['CER'].mean():9.2f} "
                  f"{c['TRU'].mean():9.2f} {c['WRG'].mean():10.2f} {wf:11.3f} "
                  f"{c['NUL'].mean():19.2f} {nf:10.3f}")
    out["ledger"] = np.array(led)

    # ---- V18 genuine/split/blend accounting
    print()
    print("  V18 accounting -- what the chain's turns actually landed on:")
    print(f"    {'arm':>4} {'K':>3} {'turns':>6} {'genuine':>8} {'split':>6} "
          f"{'blend':>6} {'seeds hit':>10} {'hit/K':>6} {'blend frac':>11} "
          f"{'deepen':>7} {'ignite':>7}")
    acct = []
    for arm in ("P", "W"):
        for c in out[f"cells_{arm}"]:
            t = float(np.median(c["TRN"])); g = float(np.median(c["GEN"]))
            sp = float(np.median(c["SPL"])); bl = float(np.median(c["BLD"]))
            hit = float(np.median(c["HIT"]))
            acct.append([0 if arm == "P" else 1, c["K"], t, g, sp, bl, hit,
                         hit / c["K"], bl / max(t, 1e-9),
                         float(c["DPN"].mean()), float(c["IGN"].mean())])
            print(f"    {arm:>4} {c['K']:3d} {t:6.1f} {g:8.1f} {sp:6.1f} {bl:6.1f} "
                  f"{hit:10.1f} {hit/c['K']:6.2f} {bl/max(t,1e-9):11.2f} "
                  f"{c['DPN'].mean():7.2f} {c['IGN'].mean():7.2f}")
    out["accounting"] = np.array(acct)

    # ---- chance-match floors (K-dependent: the GENUINE pool is now K-invariant)
    print()
    p1 = (1.0 - np.cos(W.TOL_RAD)) / 2.0
    print("  chance-match floor (GENUINE scored against this bin's seeds + top-10 sky):")
    for K in K_LIST:
        sP = int(np.ceil(K / len(L.BAND)))
        for arm, s in (("P", sP), ("W", 1)):
            M = s + W.TOP_TRUE
            print(f"    arm {arm} K = {K:2d}: pool = {s} seed(s) + {W.TOP_TRUE} sky "
                  f"= {M:2d}, floor {1.0 - (1.0 - p1) ** M:.3f}")
    out["ctrl"] = np.array([dete, n_real, f_(cs), f_(ct), f_(rx), f_(sm)])


def report_bracket(out):
    print()
    print("#" * 120)
    print("# THE DELIVERABLE -- K*_50 AND K*_90 AS A BRACKET OVER CONSTRUCTION")
    print("#" * 120)
    for name, _, desc in CRIT:
        print(f"\n  criterion: {desc}")
        for thr, lbl in ((0.5, "K*_50"), (0.9, "K*_90")):
            parts = []
            for arm, _, _, _ in ARMS:
                v = out.get(f"kstar_{arm}_{name}_{int(thr*100)}")
                if v is None or np.isnan(v[0]):
                    parts.append(f"arm {arm}: > {max(K_LIST)}")
                else:
                    parts.append(f"arm {arm}: {v[1]:.1f} +- {v[2]:.1f}")
            print(f"    {lbl}:  " + "   |   ".join(parts))
    print()
    print("  Arm W (one seed per bin) is the optimistic construction and Arm P")
    print("  (packed into 10-25 nHz) the pessimistic one; the truth is between.")


def _bank(out):
    sav = {}
    for k, v in out.items():
        if isinstance(v, np.ndarray):
            sav[k] = v
    for arm, _, _, _ in ARMS:
        for c in out.get(f"cells_{arm}", []):
            for k, v in c.items():
                if isinstance(v, np.ndarray):
                    sav[f"{arm}{c['K']}_{k}"] = v
                elif isinstance(v, float):
                    sav[f"{arm}{c['K']}_{k}"] = np.array([v])
    np.savez(os.path.join(REPORTS, "vernier27.npz"), **sav)


def main():
    t0 = time.time()
    out = {}
    print()
    print("#" * 120)
    print("# VERNIER-27 -- FIRM THE K* NUMBER.  EXPLORATORY, NOT PRE-REGISTERED")
    print("#" * 120)
    print(f"  arm P band {L.BAND[0]*1e9:.2f}-{L.BAND[-1]*1e9:.2f} nHz, "
          f"{len(L.BAND)} bins of 1/T = {W.DF*1e9:.3f} nHz; array {W.N_PSR} pulsars")
    print(f"  seeds at array SNR {L.RHO_SEED:.0f}; Mc 10^{L.MC_LO}-10^{L.MC_HI}; "
          f"gamma = 1.56 (corrected sky, audited counts)")
    print(f"  COMB_CAP {L.COMB_CAP} (V25 used 8); max_turn = 2K + 16 (V25 used 16)")
    print()
    print("  V25 Arm K DEFECT: BAND holds 10 bins and chain() placed one seed in")
    print("  bins b < n_seed, so K = 12 and K = 16 both placed TEN seeds and are the")
    print("  same cell. The V25 extrapolation to K* ~ 20-25 ran off a tail that had")
    print("  stopped moving for that reason. V25 results at K <= 10 are unaffected.")

    for arm, desc, _, n_chain in ARMS:
        cells = run_arm(out, arm, desc, n_chain)
        rows = report_grid(out, arm, cells)
        report_kstar(out, arm, cells)
        report_slope(out, arm, cells, rows)
    report_bracket(out)
    armC_capsens(out)
    controls(out)

    _bank(out)
    print()
    print(f"  banked -> {os.path.join(REPORTS, 'vernier27.npz')}   "
          f"({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
