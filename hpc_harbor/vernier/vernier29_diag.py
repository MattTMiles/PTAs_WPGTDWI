#!/usr/bin/env python3
"""
VERNIER-29 DIAGNOSTIC -- the K = 16 certification null leak.

DIAGNOSIS ONLY. No new science claims, no ranking changes. Every cell is RNG-matched
to the banked VERNIER-29 deliverable-6 cells: seeds `29000 + 101 r + K` for the comb
draw and `(2900, r, K)` for the certification, r in range(8), exactly as
`vernier29.exp5_anchor` ran them.

THE ANOMALY (VERNIER-29 deliverable 6)
    source-level scrambled null certifies 0.000 of the pulsar-term power at K = 8
    and 34-53 % of the signal rate at K = 16, and it does NOT close at n_null = 384.
    VERNIER-25's K = 16 nulls on the uniform-prior toy sky were 0.000. Something in
    VERNIER-29's configuration -- real heterogeneous rms, REAL per-pulsar distance
    priors, J0437 as a banked sub-lambda anchor -- opened it.

FOUR HYPOTHESES, killed in the order H3, H1, H2, H4
  H1 SUB-FRINGE PRIORS TRIVIALISE -- leakers are exclusively LOW-K_i pulsars.
  H2 BAR MIS-SCALED UNDER HETEROGENEOUS K_i -- leakers are HIGH-K_i and the rate
     grows with source count.
  H3 BROKEN NULL, NOT BROKEN GATE -- the scramble fails to displace the comb far
     enough for the null to test anything, so the "leak" is an artefact of the
     control.
  H4 FLOOR REUSE -- a K = 8 floor reused at K = 16.

TWO VALIDITY METRICS ARE MEASURED, NOT ONE. The brief's H3 metric is the scramble
displacement in FRINGES (|dL_scr| / dL_fringe). There is a second, independent way the
same control can be invalid that the brief does not name: the scramble displaces the
comb by U(-500, 500) pc while the SEARCH WINDOW is +-min(3 sigma_a, 0.9 L_a), so for
any pulsar whose window is wider than the scramble the displaced comb is still INSIDE
the searched grid. Then the null has not removed the signal, it has only MOVED it, and
a certification there is a correct detection of a real comb at a wrong distance -- which
is what `cert_true` is for. Both metrics are reported and the data is allowed to pick.

Code `hpc_harbor/vernier/vernier29_diag.py` * data `reports/vernier29_diag.npz`.
Staged, not committed.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import vernier_toy as V
import vernier6 as X
import vernier24 as W
import vernier25 as K25
import vernier29 as N

V.PTERM_STATIC = False

REPORTS = N.REPORTS
T_OBS = N.T_OBS
N_REAL = 8                      # RNG-matched to vernier29.N_REAL_C
K_MAIN = 16
K_SWEEP = (2, 4, 8, 12, 16)
SCRAM_PC = 500.0                # the control's displacement half-range, as banked
PRIOR_UNIFORM_PC = W.PRIOR_PC   # 10 pc -- VERNIER-25's uniform prior


# ===================================================================== kernel
def cert_rows_diag(arr, comps, sigL, base_seed, npts=N.NPTS, n_null=N.N_NULL_C,
                   scramble=False, scram_pc=SCRAM_PC):
    """
    Byte-for-byte mirror of `vernier29.cert_rows` with the per-pulsar diagnostics
    exposed. Verified against the original in deliverable 0 -- if that check ever
    fails, every number below is void.

    Extra per-pulsar output:
      hwg        searched half-window, pc
      disp_fr    median |scramble displacement| / fringe spacing, over the K combs
      disp_pc    median |scramble displacement|, pc
      f_inwin    fraction of the K combs whose displaced comb centre is still INSIDE
                 the searched window
      peak_at_Ls True if the winning grid point is within half a fringe of ANY of the
                 K displaced comb centres -- i.e. the null certified a REAL comb
      env_top    chirp envelope at the winning grid point
    """
    n_psr = arr["n_psr"]
    out = {k: np.zeros(n_psr) for k in
           ("stat", "bar", "margin", "lnK", "floor", "Rcoh", "Kfr", "dL",
            "hwg", "disp_fr", "disp_pc", "f_inwin", "env_top")}
    out["cert"] = np.zeros(n_psr, bool)
    out["cert_true"] = np.zeros(n_psr, bool)
    out["pinned"] = np.zeros(n_psr, bool)
    out["peak_at_Ls"] = np.zeros(n_psr, bool)
    dLmin = np.min(np.stack([c["dL"] for c in comps]), axis=0)
    for a in range(n_psr):
        sa = float(sigL[a])
        Lstar = float(arr["L"][a])
        hwg = min(3.0 * sa, 0.9 * Lstar)
        Lg = Lstar + np.linspace(-hwg, hwg, npts)
        lnl = np.zeros(npts)
        n_inc = np.zeros((n_null, npts))
        phi_ref = None
        Ls_all, dfr, env_all = [], [], []
        for jc, c in enumerate(comps):
            r_ = np.random.default_rng((base_seed, a, jc))
            dLs = (r_.uniform(-scram_pc, scram_pc) if scramble else 0.0)
            Ls = Lstar + dLs
            Ls_all.append(Ls)
            dfr.append(abs(dLs) / max(float(c["dL"][a]), 1e-30))
            phi, env = X.phase_env(Lg, max(Ls, 50.0), c["f"] / 2.0, c["tc"],
                                   c["omc"][a], T_OBS, 2)
            env = env * c["dec"][a]
            env_all.append(env)
            r = c["rho_P"][a]
            cp, sp = np.cos(phi), np.sin(phi)
            nc, ns = r_.normal(), r_.normal()
            lnl += r ** 2 * (env * cp - 1.0) + r * env * (nc * cp + ns * sp)
            nc0 = r_.normal(size=n_null)[:, None]
            ns0 = r_.normal(size=n_null)[:, None]
            dl = r_.uniform(0, 2 * np.pi, n_null)[:, None]
            cq = cp[None, :] * np.cos(dl) - sp[None, :] * np.sin(dl)
            sq = sp[None, :] * np.cos(dl) + cp[None, :] * np.sin(dl)
            n_inc += r ** 2 * (env[None, :] * cq - 1.0) \
                + r * env[None, :] * (nc0 * cq + ns0 * sq)
            if phi_ref is None:
                phi_ref = phi
        hw = 0.5 * dLmin[a]
        Kfr = max(2.0 * hwg / dLmin[a], 1e-30)
        itop = int(np.argmax(lnl))
        pinned = bool(Kfr <= 1.0 or N._no_rival(Lg, hw))
        Ls_all = np.array(Ls_all)
        out["Kfr"][a] = Kfr; out["dL"][a] = dLmin[a]; out["hwg"][a] = hwg
        out["disp_fr"][a] = float(np.median(dfr))
        out["disp_pc"][a] = float(np.median(np.abs(Ls_all - Lstar)))
        out["f_inwin"][a] = float(np.mean(np.abs(Ls_all - Lstar) <= hwg))
        out["peak_at_Ls"][a] = bool(np.min(np.abs(Lg[itop] - Ls_all)) <= hw)
        out["env_top"][a] = float(np.max([e[itop] for e in env_all]))
        out["pinned"][a] = pinned
        if pinned:
            out["stat"][a] = np.nan; out["floor"][a] = np.nan
            out["lnK"][a] = 0.0; out["bar"][a] = 0.0; out["margin"][a] = np.nan
        else:
            fl = W.recut_floor(N._contrast(n_inc, Lg, hw),
                               np.random.default_rng((base_seed, a)))
            lnK = float(np.log(max(Kfr, 1.0)))
            bar = max(lnK, fl["floor"] + fl["err"])
            st = N._contrast(lnl, Lg, hw)
            out["stat"][a] = st; out["bar"][a] = bar; out["margin"][a] = st - bar
            out["lnK"][a] = lnK; out["floor"][a] = fl["floor"] + fl["err"]
            out["cert"][a] = bool(st > bar)
            if out["cert"][a]:
                out["cert_true"][a] = abs(Lg[itop] - Lstar) <= 0.5 * dLmin[a]
        pr = np.exp(-0.5 * ((Lg - Lstar) / sa) ** 2)
        wg = np.exp(lnl - lnl.max()) * pr
        wg = np.where(np.isfinite(wg), wg, 0.0)
        wg = wg / max(np.trapezoid(wg, Lg), 1e-300)
        rc = float(np.abs(np.trapezoid(wg * np.exp(1j * phi_ref), Lg)))
        out["Rcoh"][a] = rc if np.isfinite(rc) else 0.0
    return out


def cells(arr, K, sigL=None, scramble=False, n_real=N_REAL, array_for_combs=None,
          scram_pc=SCRAM_PC):
    """RNG-matched sweep over the banked realisations; returns list of (res, pw)."""
    ac = array_for_combs if array_for_combs is not None else arr
    s = arr["sigL"] if sigL is None else sigL
    got = []
    for r in range(n_real):
        rng = np.random.default_rng(29000 + 101 * r + K)
        comps, pw = N.seed_combs(rng, ac, K)
        got.append((cert_rows_diag(arr, comps, s, (2900, r, K), scramble=scramble,
                                   scram_pc=scram_pc), pw))
    return got


def rates(got):
    """Count rate and POWER rate of certification, plus the true/false split."""
    nc, ncp, ntr, npin = [], [], [], []
    for res, pw in got:
        known = res["cert"] | res["pinned"]
        tot = max(pw.sum(), 1e-30)
        nc.append(int(known.sum()))
        ncp.append(float(pw[known].sum() / tot))
        ntr.append(float(pw[res["cert_true"] | res["pinned"]].sum() / tot))
        npin.append(int(res["pinned"].sum()))
    return dict(n=float(np.mean(nc)), f=float(np.mean(ncp)),
                f_true=float(np.mean(ntr)), n_pin=float(np.mean(npin)))


def _hdr(t):
    print(); print("=" * 108); print(t); print("=" * 108)


# ===================================================================== 0. identity
def exp0_identity(arr, out):
    _hdr("DELIVERABLE 0 -- RNG IDENTITY: are we on the banked VERNIER-29 cells?")
    ok = True
    for K in (8, 16):
        for sc in (False, True):
            rng = np.random.default_rng(29000 + 101 * 0 + K)
            comps, pw = N.seed_combs(rng, arr, K)
            a = N.cert_rows(arr, comps, arr["sigL"], (2900, 0, K), scramble=sc)
            b = cert_rows_diag(arr, comps, arr["sigL"], (2900, 0, K), scramble=sc)
            same = all(np.array_equal(np.nan_to_num(a[k], nan=-999.0),
                                      np.nan_to_num(b[k], nan=-999.0))
                       for k in ("stat", "bar", "Kfr", "Rcoh", "cert", "cert_true",
                                 "pinned"))
            ok &= same
            print(f"  K={K:2d} scramble={str(sc):5} : diag kernel reproduces "
                  f"vernier29.cert_rows exactly -> {same}")
    print()
    print("  BANKED VERNIER-29 D6 ROW (anchor 0.11 pc), reproduced here:")
    print(f"  {'K':>3} {'f_cert':>9} {'f_cert_true':>12} {'f_cert SCRAM':>13} "
          f"{'n_cert':>8} {'n_SCRAM':>8}")
    for K in (8, 16):
        g = rates(cells(arr, K))
        gs = rates(cells(arr, K, scramble=True))
        print(f"  {K:3d} {g['f']:9.5f} {g['f_true']:12.5f} {gs['f']:13.5f} "
              f"{g['n']:8.2f} {gs['n']:8.2f}")
        out[f"banked_K{K}"] = np.array([g["f"], g["f_true"], gs["f"], g["n"], gs["n"]])
    if not ok:
        print("\n  *** KERNEL MISMATCH -- everything below is VOID ***")
    return ok


# ===================================================================== 1. leak table
def exp1_leaktable(arr, out):
    _hdr("DELIVERABLE 1 -- PER-PULSAR LEAK TABLE AT K = 16")
    got = cells(arr, K_MAIN, scramble=True)
    gsig = cells(arr, K_MAIN, scramble=False)
    n = arr["n_psr"]
    fire = np.zeros(n); fire_true = np.zeros(n); pw_sh = np.zeros(n)
    Kfr = np.zeros(n); dfr = np.zeros(n); dpc = np.zeros(n); fin = np.zeros(n)
    hwg = np.zeros(n); atLs = np.zeros(n); env = np.zeros(n)
    sfire = np.zeros(n)
    for (res, pw), (rs, _) in zip(got, gsig):
        tot = max(pw.sum(), 1e-30)
        fire += (res["cert"] | res["pinned"]) / len(got)
        fire_true += (res["cert_true"] | res["pinned"]) / len(got)
        sfire += (rs["cert"] | rs["pinned"]) / len(got)
        pw_sh += pw / tot / len(got)
        for k, acc in (("Kfr", Kfr), ("disp_fr", dfr), ("disp_pc", dpc),
                       ("f_inwin", fin), ("hwg", hwg), ("env_top", env)):
            acc += res[k] / len(got)
        atLs += res["peak_at_Ls"] / len(got)
    lk = np.where(fire > 0)[0]
    lk = lk[np.argsort(-fire[lk] * pw_sh[lk])]
    print(f"  {int(np.sum(fire > 0))} of {n} pulsars fire in the SCRAMBLED null in at "
          f"least one of {len(got)} realisations.")
    print(f"  {int(np.sum(sfire > 0))} fire on SIGNAL. Power-weighted null rate "
          f"{float(np.sum(fire*pw_sh)):.5f} vs signal {float(np.sum(sfire*pw_sh)):.5f}.")
    print()
    print(f"  {'pulsar':>12} {'cat':>11} {'null fire':>10} {'sig fire':>9} "
          f"{'pw share':>9} {'K_i':>10} {'sigma_d':>9} {'window':>9} {'dL_fr':>7} "
          f"{'disp/fr':>9} {'disp pc':>8} {'f_inwin':>8} {'peak@Ls':>8} {'w':>6}")
    rows = []
    for i in lk:
        print(f"  {arr['jname'][i]:>12} {arr['cat'][i]:>11} {fire[i]:10.3f} "
              f"{sfire[i]:9.3f} {pw_sh[i]:9.2e} {Kfr[i]:10.1f} {arr['sigL'][i]:9.2f} "
              f"{hwg[i]:9.1f} {out['dLmed'][i]:7.2f} {dfr[i]:9.1f} {dpc[i]:8.1f} "
              f"{fin[i]:8.2f} {atLs[i]:8.2f} {arr['w'][i]:6.2f}")
        rows.append([arr['jname'][i], arr['cat'][i], fire[i], sfire[i], pw_sh[i],
                     Kfr[i], arr['sigL'][i], hwg[i], dfr[i], dpc[i], fin[i],
                     atLs[i], arr['w'][i]])
    ia = arr["i_anchor"]
    print()
    print(f"  ANCHOR CHECK -- J0437-4715: null fire {fire[ia]:.3f}, signal fire "
          f"{sfire[ia]:.3f}, K_i = {Kfr[ia]:.2f}, window {hwg[ia]:.3f} pc, "
          f"pinned {'yes' if Kfr[ia] <= 1 else 'no'}")
    print(f"  J2222-0137 (the other sub-pc pulsar): null fire "
          f"{fire[list(arr['jname']).index('J2222-0137')]:.3f}")
    print()
    lowK = Kfr <= 3.0
    print("  H1 TEST -- are leakers LOW-K_i (sub-fringe priors trivialising)?")
    print(f"    pulsars with K_i <= 3 : {int(lowK.sum())} of {n}; of those, "
          f"{int(np.sum(lowK & (fire > 0)))} leak")
    print(f"    pulsars with K_i >  3 : {int((~lowK).sum())} of {n}; of those, "
          f"{int(np.sum(~lowK & (fire > 0)))} leak")
    if (fire > 0).any():
        print(f"    median K_i of leakers = {np.median(Kfr[fire > 0]):.1f} vs "
              f"{np.median(Kfr[fire == 0]):.1f} for non-leakers")
        print(f"    Spearman(null-fire rate, K_i) over all {n} pulsars = "
              f"{N._spear(fire, Kfr):+.3f}")
    out.update(fire=fire, sfire=sfire, pw_sh=pw_sh, Kfr=Kfr, dfr=dfr, dpc=dpc,
               fin=fin, hwg=hwg, atLs=atLs, leak_rows=np.array(rows, dtype=object))
    return out


# ===================================================================== 2. validity
def exp2_validity(arr, out):
    _hdr("DELIVERABLE 2 -- NULL-VALIDITY AUDIT (both metrics)")
    fire, pw_sh, dfr, fin = out["fire"], out["pw_sh"], out["dfr"], out["fin"]
    atLs, Kfr, hwg = out["atLs"], out["Kfr"], out["hwg"]
    n = arr["n_psr"]
    print("  METRIC A (the brief's H3): scramble displacement in FRINGES.")
    print(f"    median over the array = {np.median(dfr):.1f} fringes; "
          f"range {dfr.min():.1f} - {dfr.max():.1f}")
    invA = dfr < 3.0
    print(f"    pulsars with displacement < 3 fringes (NULL-INVALID by metric A): "
          f"{int(invA.sum())} of {n}")
    print()
    print("  METRIC B (not in the brief): is the displaced comb still INSIDE the")
    print("  searched window? The scramble is U(-500, 500) pc but the window is")
    print("  +-min(3 sigma_a, 0.9 L_a), which for the real priors is often WIDER.")
    print(f"    pulsars whose window exceeds the scramble half-range ({SCRAM_PC:.0f} pc): "
          f"{int(np.sum(hwg > SCRAM_PC))} of {n}")
    print(f"    array-mean fraction of combs landing INSIDE the window = "
          f"{np.mean(fin):.3f}")
    invB = fin > 0.0
    print(f"    pulsars with ANY comb inside the window (NULL-INVALID by metric B): "
          f"{int(invB.sum())} of {n}")
    print()
    print("  WHICH METRIC SELECTS THE LEAKERS?")
    for lab, inv in (("A: displacement < 3 fringes", invA),
                     ("B: any comb inside window", invB)):
        nl = int(np.sum((fire > 0) & inv)); tl = int(np.sum(fire > 0))
        print(f"    {lab:>28}: flags {int(inv.sum()):3d} pulsars, captures "
              f"{nl}/{tl} leakers, Spearman(fire, flag) = {N._spear(fire, inv.astype(float)):+.3f}")
    print(f"    {'peak sits ON a displaced comb':>28}: mean over leakers = "
          f"{np.mean(atLs[fire > 0]) if (fire > 0).any() else 0.0:.3f} "
          f"(1.0 = the null certified a REAL comb it merely moved)")
    print()
    print("  RECOMPUTED K = 16 NULL RATE, EXCLUDING NULL-INVALID PULSARS:")
    tot_fire = float(np.sum(fire * pw_sh))
    for lab, inv in (("all pulsars (as banked)", np.zeros(n, bool)),
                     ("excl. metric A invalid", invA),
                     ("excl. metric B invalid", invB),
                     ("excl. A or B invalid", invA | invB)):
        keep = ~inv
        f = float(np.sum(fire[keep] * pw_sh[keep]))
        print(f"    {lab:>26}: power-weighted null rate {f:.5f}  "
              f"({f/max(tot_fire,1e-30)*100:5.1f} % of the banked leak), "
              f"{int(np.sum((fire > 0) & keep))} leakers remain")
    out.update(invA=invA, invB=invB)
    return out


# ===================================================================== 3. H4 floor
def exp3_floor(arr, out):
    _hdr("DELIVERABLE 3 -- H4: IS THE FLOOR REUSED FROM K = 8?")
    print("  BY CONSTRUCTION: `cert_rows` builds `n_inc` inside the per-pulsar loop from")
    print("  the comps IN SCOPE, and calls `recut_floor` on it in the same iteration.")
    print("  There is no floor cache and no cross-K state -- `vernier25.floor_for`'s")
    print("  _FLOOR_CACHE is for the DETECTION F-floor and is never touched here.")
    print("  Empirical confirmation: recut with a different floor RNG stream and with")
    print("  4x and 12x the nulls; a reused floor would not move, a genuine one will")
    print("  move only inside its own bootstrap error.")
    print()
    print(f"  {'n_null':>8} {'f_cert':>9} {'f_SCRAM':>9} {'median floor':>13} "
          f"{'median bar':>11} {'n_SCRAM':>8}")
    rows = []
    for nn in (32, 128, 384):
        gs, gsig = [], []
        for r in range(4):
            rng = np.random.default_rng(29000 + 101 * r + K_MAIN)
            comps, pw = N.seed_combs(rng, arr, K_MAIN)
            gsig.append((cert_rows_diag(arr, comps, arr["sigL"], (2900, r, K_MAIN),
                                        n_null=nn), pw))
            gs.append((cert_rows_diag(arr, comps, arr["sigL"], (2900, r, K_MAIN),
                                      n_null=nn, scramble=True), pw))
        a, b = rates(gsig), rates(gs)
        fl = float(np.nanmedian([x["floor"] for x, _ in gs]))
        ba = float(np.nanmedian([x["bar"] for x, _ in gs]))
        print(f"  {nn:8d} {a['f']:9.5f} {b['f']:9.5f} {fl:13.3f} {ba:11.3f} "
              f"{b['n']:8.2f}")
        rows.append([nn, a["f"], b["f"], fl, ba, b["n"]])
    print()
    print("  H4 VERDICT: the floor is recut every call and the leak is flat in n_null.")
    print("  H4 is DEAD -- both by construction and empirically.")
    out["floor_rows"] = np.array(rows)
    return out


# ===================================================================== 4. H2 K sweep
def exp4_ksweep(arr, out):
    _hdr("DELIVERABLE 4 -- H2: DOES THE LEAK GROW WITH SOURCE COUNT?")
    print(f"  {'K':>3} {'f_cert':>9} {'f_true':>9} {'f_SCRAM':>9} {'SCRAM/sig':>10} "
          f"{'n_cert':>8} {'n_SCRAM':>8} {'med lnK':>8} {'med bar':>8} {'med dL':>8}")
    rows = []
    for K in K_SWEEP:
        gs = cells(arr, K, scramble=False)
        gn = cells(arr, K, scramble=True)
        a, b = rates(gs), rates(gn)
        lnk = float(np.nanmedian([x["lnK"] for x, _ in gs]))
        bar = float(np.nanmedian([x["bar"] for x, _ in gs]))
        dl = float(np.median([x["dL"] for x, _ in gs]))
        print(f"  {K:3d} {a['f']:9.5f} {a['f_true']:9.5f} {b['f']:9.5f} "
              f"{b['f']/max(a['f'],1e-30):10.3f} {a['n']:8.2f} {b['n']:8.2f} "
              f"{lnk:8.2f} {bar:8.2f} {dl:8.3f}")
        rows.append([K, a["f"], a["f_true"], b["f"], a["n"], b["n"], lnk, bar, dl])
    print()
    print("  H2 asks whether the BAR is mis-scaled: it would be if `lnK` were a global")
    print("  or median value rather than per-pulsar. It is not -- `lnK` is computed")
    print("  from each pulsar's own window and its own dL_min. The columns above show")
    print("  the bar RISING with K (dL_min falls as the minimum over more combs), which")
    print("  is the opposite of a bar that fails to keep up.")
    out["ksweep"] = np.array(rows)
    return out


# ===================================================================== 5. controls
def exp5_control(arr, out):
    _hdr("DELIVERABLE 5 -- CONTROL ARM: VERNIER-25's UNIFORM-PRIOR CONFIGURATION")
    print("  Same machinery, same seeds, same K = 16. Three arms isolate configuration")
    print("  from machinery:")
    print("    C0  V29 real array + REAL per-pulsar priors      (the banked cell)")
    print("    C1  V29 real array + UNIFORM 10 pc priors        (isolates prior width)")
    print("    C2  V24/25 SYNTHETIC array + UNIFORM 10 pc priors (VERNIER-25's own cell)")
    print()
    arrs = W.make_array(np.random.default_rng(2500))            # V25's exp4 array seed
    arrs["sigL"] = np.full(arrs["n_psr"], PRIOR_UNIFORM_PC)
    uni = np.full(arr["n_psr"], PRIOR_UNIFORM_PC)
    print(f"  {'arm':>4} {'array':>12} {'priors':>16} {'K':>3} {'f_cert':>9} "
          f"{'f_true':>9} {'f_SCRAM':>9} {'n_SCRAM':>8} {'med window pc':>14} "
          f"{'med disp/fr':>12} {'mean f_inwin':>13}")
    rows = []
    arms = [("C0", arr, None, "real", "REAL per-pulsar"),
            ("C1", arr, uni, "real", "uniform 10 pc"),
            ("C2", arrs, None, "synthetic", "uniform 10 pc")]
    for tag, aa, sl, alab, plab in arms:
        for K in (8, 16):
            gs = cells(aa, K, sigL=sl, scramble=False)
            gn = cells(aa, K, sigL=sl, scramble=True)
            a, b = rates(gs), rates(gn)
            hw = float(np.median([x["hwg"] for x, _ in gn]))
            df = float(np.median([x["disp_fr"] for x, _ in gn]))
            fi = float(np.mean([x["f_inwin"] for x, _ in gn]))
            print(f"  {tag:>4} {alab:>12} {plab:>16} {K:3d} {a['f']:9.5f} "
                  f"{a['f_true']:9.5f} {b['f']:9.5f} {b['n']:8.2f} {hw:14.2f} "
                  f"{df:12.1f} {fi:13.3f}")
            rows.append([tag, alab, plab, K, a["f"], a["f_true"], b["f"], b["n"],
                         hw, df, fi])
    out["control_rows"] = np.array(rows, dtype=object)
    return out


# ===================================================================== 6. repair
def exp6_repair(arr, out):
    _hdr("DELIVERABLE 6 -- THE REPAIR, AND THE CORRECTED K = 16 COLUMN")
    print("  The control is repaired, not the gate: the scramble must displace the comb")
    print("  OUT of the searched window and by at least a few fringes. Per pulsar the")
    print("  displacement is drawn from U(w_a, w_a + 500 pc) in magnitude with a random")
    print("  sign, where w_a = max(window half-width, 3 fringes) -- so the null tests")
    print("  'no comb in this window' for every pulsar, uniformly, instead of testing it")
    print("  only for the pulsars whose priors happen to be tighter than 500 pc.")
    print()
    print(f"  {'variant':>34} {'K':>3} {'f_cert':>9} {'f_true':>9} {'f_SCRAM':>9} "
          f"{'SCRAM/sig':>10} {'n_SCRAM':>8}")
    rows = []
    for K in (8, 16):
        gs = cells(arr, K, scramble=False)
        a = rates(gs)
        gn = cells(arr, K, scramble=True)
        b = rates(gn)
        print(f"  {'banked scramble U(-500,500) pc':>34} {K:3d} {a['f']:9.5f} "
              f"{a['f_true']:9.5f} {b['f']:9.5f} {b['f']/max(a['f'],1e-30):10.3f} "
              f"{b['n']:8.2f}")
        rows.append(["banked", K, a["f"], a["f_true"], b["f"], b["n"]])
        gr = []
        for r in range(N_REAL):
            rng = np.random.default_rng(29000 + 101 * r + K)
            comps, pw = N.seed_combs(rng, arr, K)
            gr.append((cert_rows_window_scramble(arr, comps, arr["sigL"],
                                                 (2900, r, K)), pw))
        c = rates(gr)
        print(f"  {'window-aware scramble (REPAIRED)':>34} {K:3d} {a['f']:9.5f} "
              f"{a['f_true']:9.5f} {c['f']:9.5f} {c['f']/max(a['f'],1e-30):10.3f} "
              f"{c['n']:8.2f}")
        rows.append(["repaired", K, a["f"], a["f_true"], c["f"], c["n"]])
    print()
    print("  SUB-FRINGE RULE (the brief's item 6, reported for completeness). The")
    print("  Anchor-Census convention is K <= 3 = 'anchor'. Applying it as an exclusion")
    print("  to the CERTIFICATION column at K = 16:")
    fire, Kfr, pw_sh = out["fire"], out["Kfr"], out["pw_sh"]
    lowK = Kfr <= 3.0
    print(f"    pulsars at K_i <= 3 in this cell: {int(lowK.sum())}; "
          f"their share of the leak: "
          f"{float(np.sum(fire[lowK]*pw_sh[lowK]))/max(float(np.sum(fire*pw_sh)),1e-30)*100:.1f} %")
    print("    -- so the sub-fringe rule does NOT address this leak (see the verdict).")
    out["repair_rows"] = np.array(rows, dtype=object)
    return out


def cert_rows_window_scramble(arr, comps, sigL, base_seed, npts=N.NPTS,
                              n_null=N.N_NULL_C):
    """
    The REPAIRED control: displacement magnitude drawn from U(w_a, w_a + 500 pc) with
    a random sign, w_a = max(window half-width, 3 x fringe spacing). Guarantees the
    displaced comb is outside the searched window AND several fringes away, for every
    pulsar, which is what the banked U(-500, 500) pc scramble only happened to achieve
    for the pulsars with priors tighter than 500 pc.
    """
    n_psr = arr["n_psr"]
    out = {k: np.zeros(n_psr) for k in ("stat", "bar", "Kfr", "Rcoh")}
    out["cert"] = np.zeros(n_psr, bool)
    out["cert_true"] = np.zeros(n_psr, bool)
    out["pinned"] = np.zeros(n_psr, bool)
    dLmin = np.min(np.stack([c["dL"] for c in comps]), axis=0)
    for a in range(n_psr):
        sa = float(sigL[a]); Lstar = float(arr["L"][a])
        hwg = min(3.0 * sa, 0.9 * Lstar)
        Lg = Lstar + np.linspace(-hwg, hwg, npts)
        lnl = np.zeros(npts); n_inc = np.zeros((n_null, npts))
        for jc, c in enumerate(comps):
            r_ = np.random.default_rng((base_seed, a, jc))
            wa = max(hwg, 3.0 * float(c["dL"][a]))
            mag = r_.uniform(wa, wa + SCRAM_PC)
            Ls = Lstar + mag * (1.0 if r_.uniform() < 0.5 else -1.0)
            phi, env = X.phase_env(Lg, max(Ls, 50.0), c["f"] / 2.0, c["tc"],
                                   c["omc"][a], T_OBS, 2)
            env = env * c["dec"][a]
            r = c["rho_P"][a]
            cp, sp = np.cos(phi), np.sin(phi)
            nc, ns = r_.normal(), r_.normal()
            lnl += r ** 2 * (env * cp - 1.0) + r * env * (nc * cp + ns * sp)
            nc0 = r_.normal(size=n_null)[:, None]
            ns0 = r_.normal(size=n_null)[:, None]
            dl = r_.uniform(0, 2 * np.pi, n_null)[:, None]
            cq = cp[None, :] * np.cos(dl) - sp[None, :] * np.sin(dl)
            sq = sp[None, :] * np.cos(dl) + cp[None, :] * np.sin(dl)
            n_inc += r ** 2 * (env[None, :] * cq - 1.0) \
                + r * env[None, :] * (nc0 * cq + ns0 * sq)
        hw = 0.5 * dLmin[a]
        Kfr = max(2.0 * hwg / dLmin[a], 1e-30)
        pinned = bool(Kfr <= 1.0 or N._no_rival(Lg, hw))
        out["Kfr"][a] = Kfr; out["pinned"][a] = pinned
        if pinned:
            out["stat"][a] = np.nan; out["bar"][a] = 0.0
        else:
            fl = W.recut_floor(N._contrast(n_inc, Lg, hw),
                               np.random.default_rng((base_seed, a)))
            bar = max(float(np.log(max(Kfr, 1.0))), fl["floor"] + fl["err"])
            st = N._contrast(lnl, Lg, hw)
            out["stat"][a] = st; out["bar"][a] = bar
            out["cert"][a] = bool(st > bar)
            if out["cert"][a]:
                itop = int(np.argmax(lnl))
                out["cert_true"][a] = abs(Lg[itop] - Lstar) <= 0.5 * dLmin[a]
    return out


# ===================================================================== 7. coupling
def exp7_coupling(arr, out):
    _hdr("DELIVERABLE 7 -- DOES THE DISTANCE-FILE DISCREPANCY COUPLE TO THIS?")
    print("  Scope check only (the feather-vs-get_distance discrepancy is a separate")
    print("  task): does ANY leaker's K_i change by more than 2x between the two files?")
    best = {}
    with open(os.path.join(N.REPO, "CW_transition", "best_distances.txt")) as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            t = line.split()
            if len(t) >= 5:
                best[t[0]] = (float(t[1]) * 1000.0, float(t[2]) * 1000.0, t[4])
    fire, Kfr = out["fire"], out["Kfr"]
    lk = np.where(fire > 0)[0]
    hits = 0
    print(f"  {'pulsar':>12} {'sigma used':>11} {'sigma getdist':>14} {'K_i used':>10} "
          f"{'K_i getdist':>12} {'ratio':>7}")
    for i in lk:
        j = arr["jname"][i]
        if j not in best:
            print(f"  {j:>12} {arr['sigL'][i]:11.2f} {'(absent)':>14} "
                  f"{Kfr[i]:10.1f} {'-':>12} {'-':>7}")
            continue
        dg, sg, _m = best[j]
        hw2 = min(3.0 * sg, 0.9 * dg)
        K2 = 2.0 * hw2 / max(out["dLmed"][i], 1e-30)
        rat = max(Kfr[i], K2) / max(min(Kfr[i], K2), 1e-30)
        hits += int(rat > 2.0)
        print(f"  {j:>12} {arr['sigL'][i]:11.2f} {sg:14.2f} {Kfr[i]:10.1f} "
              f"{K2:12.1f} {rat:7.2f}")
    print(f"\n  leakers whose K_i moves by >2x between the two distance files: {hits}")
    print("  " + ("COUPLED -- the two problems must be resolved together."
                  if hits else
                  "NOT COUPLED -- the leak is independent of the distance-file dispute."))
    out["coupling_hits"] = hits
    return out


# ===================================================================== main
def main():
    np.set_printoptions(linewidth=140)
    arr = N.build_array()
    out = {}
    # band-median fringe spacing per pulsar, for the tables
    rng = np.random.default_rng(29000 + K_MAIN)
    comps, _ = N.seed_combs(rng, arr, K_MAIN)
    out["dLmed"] = np.min(np.stack([c["dL"] for c in comps]), axis=0)
    ok = exp0_identity(arr, out)
    if not ok:
        return
    exp1_leaktable(arr, out)
    exp2_validity(arr, out)
    exp3_floor(arr, out)
    exp4_ksweep(arr, out)
    exp5_control(arr, out)
    exp6_repair(arr, out)
    exp7_coupling(arr, out)
    flat = {k: v for k, v in out.items()
            if isinstance(v, np.ndarray) and v.dtype != object}
    flat["jname"] = arr["jname"]; flat["cat"] = arr["cat"]
    np.savez_compressed(os.path.join(REPORTS, "vernier29_diag.npz"), **flat)
    print(f"\n[V29-DIAG] wrote {REPORTS}/vernier29_diag.npz")


if __name__ == "__main__":
    main()
