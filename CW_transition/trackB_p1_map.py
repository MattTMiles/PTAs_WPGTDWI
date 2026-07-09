"""Track B — P1: the registration map, batched rebuild.

Decisive question (spec / P1 status): is there a UNIQUE joint-lnL / true-fringe
registration NEEDLE at truth, width ~ the 1-rad pulsar-term-phase prediction
(2*pi*L/dL ~ 1.6e4 rad -> ~0.006 deg sky), rising out of the degrees-wide F-stat
blur? Secondary partial-registration side modes? Does the needle need all 116
pulsars or do the certified few carry it?

Rebuild vs attempt-1 (abandoned at 64/169, ~33 s/trial):
  * SINGLE-VMAP E-step. All trials share the SAME data_obs (fixed truth injection),
    so we batch fisher_logL over a FLAT (trial x pulsar x fringe) point stack against
    the one fixed data leaf -- ONE compiled shape (chunk, n_theta). This avoids the
    pathological (D, 4096, 244) DOUBLE-vmap (scan_all's vmap2, in_axes=(0,0)) that
    made attempt-0/1 slow to compile.
  * The trial is the runtime batch axis; the sky grid streams as chunked tensors.
  * QuickCW low-rank split (spec sec 6): DEFERRED unless a stage exceeds ~45 min
    (decision 2026-07-06). Per-fringe cost stays full-array here; we MEASURE first.

GATE (before any production run): the batched E-step reproduces the D=1 reference
scan (trackB_estimator.scan_current) to < 1e-8 on ~8 deterministic grid trials.
The old trackB_p1_partial.npz salvaged NO lnL values (profiled_lnl all NaN,
registration_count all -1 -- in-RAM-only, never persisted); the reference is the
live D=1 eval, not those dead values. The grid is deterministic -> the map is fully
reproduced.

Run (jug-gpu, XLA autotune off, shared 4090):
  python trackB_p1_map.py gate            # correctness gate + coarse wall-clock estimate
  python trackB_p1_map.py run             # staged zoom + freq zoom, checkpoints, report
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
import stagec_anchor_a2 as A2

CWT = "/home/mattm/projects/HSYMT/CW_transition"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]   # the 3 census-ceiling pulsars

# fixed chunk of flat (trial,pulsar,fringe) points -> ONE compiled fisher_logL_batched
# shape. 4096 = scan_all's BIG_CHUNK (known-good peak memory on the 24 GB card).
CHUNK = 4096


# ============================================================
# batched single-vmap E-step: LNL[trial, pulsar, fringe]
# ============================================================
def batched_scan(P, trials, chunk=CHUNK, verbose=False):
    """LNL (nT, npsr, B) at each trial's theta, sweeping each pulsar's own distance
    over its EV fringe grid (other distances held at the trial's theta). Single vmap
    of fisher_logL over a flat point stack vs the fixed data_obs. Fixed chunk shape."""
    fl_b = P["amo_full"]["fisher_logL_batched"]         # jit(vmap(fisher_logL, (0,None)))
    data = P["data_obs"]
    EV, AI = P["EV"], P["AI"]
    npsr, Bn = EV.shape
    nth = P["amo_full"]["n_theta"]
    M = npsr * Bn
    trials_arr = np.stack([np.asarray(t, float) for t in trials])   # (nT, nth)
    nT = trials_arr.shape[0]
    N = nT * M

    # per-within-trial-point (m in 0..M-1) maps: which pulsar, which theta slot, which eval
    a_grid = np.repeat(np.arange(npsr), Bn)             # (M,)
    f_grid = np.tile(np.arange(Bn), npsr)               # (M,)
    slot_grid = AI[a_grid]                              # (M,) theta index of that pulsar's p_dist
    ev_grid = EV[a_grid, f_grid]                        # (M,) the fringe eval distance

    LNL = np.empty(N)
    ar = np.arange(chunk)
    t0 = time.time()
    for c0 in range(0, N, chunk):
        c1 = min(c0 + chunk, N)
        p = np.arange(c0, c1)
        t_idx = p // M
        m_idx = p % M
        th = trials_arr[t_idx].copy()                  # (len,nth)
        th[np.arange(len(p)), slot_grid[m_idx]] = ev_grid[m_idx]
        if len(p) < chunk:                             # pad last chunk -> fixed shape
            pad = np.repeat(th[-1:], chunk - len(p), 0)
            th = np.concatenate([th, pad], 0)
        out = np.asarray(fl_b(jnp.asarray(th), data))
        LNL[c0:c1] = out[:c1 - c0]
        if verbose and (c0 // chunk) % 200 == 0:
            print(f"    ...{c1}/{N} pts ({time.time()-t0:.0f}s)", flush=True)
    return LNL.reshape(nT, npsr, Bn), time.time() - t0


def eval_stage(P, trials, prior_col, chunk=CHUNK, verbose=False):
    """Batched E-step + per-trial fringe posterior + profiled joint lnL.

    Per trial: snap each pulsar to its MAP fringe centre (fringe_posterior), then
    joint lnL = fisher_logL at all-116-snapped distances (fringe-profiled). Register
    a pulsar if map q_max > 0.5; true-register if additionally its MAP fringe is k=0."""
    LNL, t_scan = batched_scan(P, trials, chunk, verbose)
    nT, npsr, Bn = LNL.shape
    ndist = P["n_dist"]
    fl_b = P["amo_full"]["fisher_logL_batched"]
    best_theta = np.stack([np.asarray(t, float) for t in trials])   # (nT, nth)
    qmax = np.zeros((nT, npsr)); mapk = np.zeros((nT, npsr), int)
    tp0 = time.time()
    for t in range(nT):
        post = TE.fringe_posterior(P, LNL[t], None, prior_col, 1.0)
        best_theta[t, :ndist] = post["map_evalL"]
        qmax[t] = post["map_qmax"]; mapk[t] = post["map_k"]
    t_post = time.time() - tp0
    jlnl = np.asarray(fl_b(jnp.asarray(best_theta), P["data_obs"]))  # (nT,) one (nT,nth) compile
    reg = qmax > 0.5
    nreg = reg.sum(1); nreg_true = (reg & (mapk == 0)).sum(1)
    return dict(jlnl=np.asarray(jlnl, float), nreg=nreg, nreg_true=nreg_true,
                qmax=qmax, mapk=mapk, t_scan=t_scan, t_post=t_post)


# ============================================================
# grid construction (square in cos_gwtheta x gwphi about a centre)
# ============================================================
def make_grid(cg0, gp0, half_deg, n):
    dr = np.radians(half_deg)
    cgs = np.clip(cg0 + np.linspace(-dr, dr, n), -1, 1)
    gps = gp0 + np.linspace(-dr, dr, n)
    return cgs, gps, dr


def grid_trials(P, cg0, gp0, half_deg, n):
    """Trials = truth theta with loud0's sky replaced by each grid cell. Row-major
    (outer=cos_gwtheta, inner=gwphi), deterministic -> reproducible."""
    tt = P["theta_truth"].copy(); ndist = P["n_dist"]
    cgs, gps, dr = make_grid(cg0, gp0, half_deg, n)
    trials, coords = [], []
    for cg in cgs:
        for gp in gps:
            th = tt.copy(); th[ndist + 0] = cg; th[ndist + 1] = gp
            trials.append(th); coords.append((cg, gp))
    return trials, np.array(coords), cgs, gps, dr


def _unit(cg, gp):
    return TE._unit(cg, gp)


def ang_deg(cg, gp, cg0, gp0):
    return np.degrees(np.arccos(np.clip(np.dot(_unit(cg, gp), _unit(cg0, gp0)), -1, 1)))


# ============================================================
# finiteness-checked checkpoint (NaN fraction > 0 aborts the stage)
# ============================================================
def finite_check(tag, name, arr, abort=True):
    a = np.asarray(arr, float)
    fin = np.isfinite(a); frac = float(fin.mean())
    mn = np.nanmin(a) if fin.any() else np.nan
    md = np.nanmedian(a) if fin.any() else np.nan
    mx = np.nanmax(a) if fin.any() else np.nan
    print(f"  [{tag}] {name}: finite {frac*100:.2f}%  min/med/max = "
          f"{mn:.4f} / {md:.4f} / {mx:.4f}", flush=True)
    if frac < 1.0 and abort:
        raise RuntimeError(f"[{tag}] {name} NaN fraction {(1-frac)*100:.2f}% > 0 -> ABORT STAGE")
    return frac


def checkpoint(tag, coords, cgs, gps, res, cg0, gp0, extra=None):
    finite_check(tag, "profiled_lnl", res["jlnl"])
    path = f"{CWT}/trackB_p1_{tag}.npz"
    d = dict(tag=tag, trial_sky_positions=coords, grid_cos_gwtheta=cgs, grid_gwphi=gps,
             center_sky=np.array([cg0, gp0]), profiled_lnl=res["jlnl"],
             registration_count=res["nreg"], n_registered_true=res["nreg_true"],
             qmax=res["qmax"], mapk=res["mapk"],
             t_scan=res["t_scan"], t_post=res["t_post"])
    if extra:
        d.update(extra)
    np.savez(path, **d)
    print(f"  [{tag}] checkpoint -> {path}", flush=True)
    return path


# ============================================================
# per-stage needle analysis
# ============================================================
def analyse_stage(tag, res, coords, cgs, gps, dr, cg0, gp0, census_idx, npsr):
    n = len(cgs)
    jl = res["jlnl"].reshape(n, n)
    nrt = res["nreg_true"].reshape(n, n)
    nrg = res["nreg"].reshape(n, n)
    i0 = np.unravel_index(np.argmax(jl), jl.shape)
    best_cg, best_gp = cgs[i0[0]], gps[i0[1]]
    off = ang_deg(best_cg, best_gp, cg0, gp0)
    cc = n // 2
    print(f"\n  [{tag}] {n}x{n} half={np.degrees(dr):.4g} deg | "
          f"joint lnL max={jl.max():.2f} at cell {i0} (ang {off:.5f} deg from centre); "
          f"centre-cell lnL={jl[cc,cc]:.2f}", flush=True)
    print(f"  [{tag}] true-reg: max {int(nrt.max())}/116 (centre {int(nrt[cc,cc])}); "
          f"reg(q>0.5): max {int(nrg.max())} (centre {int(nrg[cc,cc])})", flush=True)

    # FWHM of the joint-lnL peak along both axes through the peak cell
    def fwhm(line, step_deg):
        y = line - line.max(); idx = np.where(y > -0.5)[0]
        if len(idx) < 2:
            return np.nan
        return step_deg * (idx[-1] - idx[0])
    step_deg = np.degrees(2 * dr) / (n - 1)
    w_cg = fwhm(jl[:, i0[1]], step_deg)
    w_gp = fwhm(jl[i0[0], :], step_deg)
    print(f"  [{tag}] joint-lnL FWHM(0.5 nat): {w_cg:.5f} deg (cos_gwtheta), "
          f"{w_gp:.5f} deg (gwphi)", flush=True)

    # 3-census-only vs all-116 registration map (deliverable 4c)
    reg = res["qmax"] > 0.5
    reg_true = reg & (res["mapk"] == 0)
    nrt_census = reg_true[:, census_idx].sum(1).reshape(n, n)
    print(f"  [{tag}] census-3 true-reg: max {int(nrt_census.max())}/3 "
          f"(centre {int(nrt_census.reshape(-1)[cc*n+cc])})", flush=True)
    return dict(best_cell=i0, best_cg=best_cg, best_gp=best_gp, off=off,
                w_cg=w_cg, w_gp=w_gp, jl=jl, nrt=nrt, nrt_census=nrt_census)


# ============================================================
# GATE
# ============================================================
def run_gate(P, n_trials=8):
    print("\n=== GATE: batched E-step vs D=1 scan_current reference ===", flush=True)
    tt = P["theta_truth"]; ndist = P["n_dist"]
    cg0, gp0 = tt[ndist + 0], tt[ndist + 1]
    # deterministic coarse-grid trials (same construction as production coarse stage)
    trials, coords, cgs, gps, dr = grid_trials(P, cg0, gp0, 20.0, 13)
    # sample n_trials spread across the 169-cell grid
    pick = np.linspace(0, len(trials) - 1, n_trials).astype(int)
    sub = [trials[i] for i in pick]

    LNL_b, t_scan = batched_scan(P, sub)            # warm + timed
    LNL_b2, t_scan2 = batched_scan(P, sub)          # warm-only timing
    maxdiff = 0.0
    for k, i in enumerate(pick):
        ref = TE.scan_current(P, trials[i])          # (npsr,B) D=1 reference
        d = float(np.max(np.abs(ref - LNL_b[k])))
        maxdiff = max(maxdiff, d)
    print(f"  max |batched - reference| over {n_trials} trials x 116 x 512 = {maxdiff:.2e}", flush=True)
    ok = maxdiff < 1e-8
    print(f"  GATE {'PASS' if ok else 'FAIL'} (< 1e-8)", flush=True)

    # wall-clock extrapolation to a 169-trial stage (warm scan, batched)
    npsr, Bn = P["EV"].shape
    pts = n_trials * npsr * Bn
    s_per_pt = t_scan2 / pts
    coarse_pts = 169 * npsr * Bn
    est_min = coarse_pts * s_per_pt / 60.0
    print(f"  warm batched scan: {t_scan2:.1f}s for {pts} pts = {s_per_pt*1e3:.4f} ms/pt", flush=True)
    print(f"  => est 169-trial stage E-step ~ {est_min:.1f} min "
          f"({'OK, < 45 min' if est_min < 45 else 'EXCEEDS 45 min -> low-rank split needed'})", flush=True)
    return ok, est_min


# ============================================================
# production: 9x9 2-D coarse+zoom1 (DATA-DRIVEN peak/side-lobe census) ->
# high-res 1-D cuts through TRUTH (TRUTH-PLACED needle width) ->
# uniqueness patch about the zoom1 argmax (DATA-DRIVEN) -> freq zoom.
# Grids sized so every stage < 45 min at the measured 0.544 ms/pt (decision
# 2026-07-06): batch-only, low-rank split DEFERRED-TO-B1. Each element is
# labelled truth-placed vs data-driven (P1 = physics referendum; the truth-
# blind end-to-end zoom is P2's gate).
# ============================================================
N2D = 9                      # 2-D stage side (81 trials ~ 44 min)
NCUT = 41                    # 1-D cut points (0.06 deg span -> 0.0015 deg steps)


def _local_sky_metric(cg0):
    """deg-per-unit for each sky axis about (cg0,*): d(sky_deg)/d(cos_gwtheta) and
    d(sky_deg)/d(gwphi). sin(theta)=sqrt(1-cg0^2)."""
    s = np.sqrt(max(1 - cg0 ** 2, 1e-12))
    return np.degrees(1.0 / s), np.degrees(s)     # (per d cos_gwtheta, per d gwphi)


def run_2d_stage(P, tag, cg0, gp0, half, n, prior_col, census_idx, truth_sky, placement):
    print(f"\n=== STAGE {tag} [{placement}]: {n}x{n} +/-{half} deg about "
          f"(cg={cg0:.5f}, gp={gp0:.5f}) ===", flush=True)
    trials, coords, cgs, gps, dr = grid_trials(P, cg0, gp0, half, n)
    t0 = time.time()
    res = eval_stage(P, trials, prior_col, verbose=True)
    dt = time.time() - t0
    print(f"  [{tag}] E-step {res['t_scan']:.0f}s + posterior {res['t_post']:.1f}s = stage {dt:.0f}s", flush=True)
    checkpoint(tag, coords, cgs, gps, res, cg0, gp0,
               extra=dict(census_idx=census_idx, truth_sky=truth_sky, placement=placement))
    a = analyse_stage(tag, res, coords, cgs, gps, dr, truth_sky[0], truth_sky[1], census_idx, P["npsr"])
    a["dt"] = dt; a["over45"] = dt > 45 * 60
    return a


def run_cut(P, tag, axis, cg_c, gp_c, half_deg, n, prior_col, census_idx, truth_sky):
    """1-D cut: vary one sky axis (0=cos_gwtheta, 1=gwphi) about the centre (TRUTH),
    holding the other. Returns FWHM in sky-deg + peak offset. TRUTH-PLACED."""
    tt = P["theta_truth"].copy(); ndist = P["n_dist"]
    dr = np.radians(half_deg)
    m_cg, m_gp = _local_sky_metric(cg_c)
    if axis == 0:
        vals = np.clip(cg_c + np.linspace(-dr, dr, n), -1, 1)
        x_sky = (vals - cg_c) * m_cg
    else:
        vals = gp_c + np.linspace(-dr, dr, n)
        x_sky = (vals - gp_c) * m_gp
    trials, coords = [], []
    for v in vals:
        th = tt.copy()
        if axis == 0:
            th[ndist + 0] = v; th[ndist + 1] = gp_c
        else:
            th[ndist + 0] = cg_c; th[ndist + 1] = v
        trials.append(th); coords.append((th[ndist + 0], th[ndist + 1]))
    res = eval_stage(P, trials, prior_col)
    finite_check(tag, "profiled_lnl", res["jlnl"])
    jl = res["jlnl"]; i0 = int(np.argmax(jl))
    y = jl - jl.max(); idx = np.where(y > -0.5)[0]
    fwhm = abs(x_sky[idx[-1]] - x_sky[idx[0]]) if len(idx) > 1 else np.nan
    reg_true = (res["qmax"] > 0.5) & (res["mapk"] == 0)
    axname = "cos_gwtheta" if axis == 0 else "gwphi"
    print(f"  [{tag}] {axname} cut +/-{half_deg} deg: peak at x_sky={x_sky[i0]:+.5f} deg "
          f"(offset from truth), FWHM(0.5nat)={fwhm:.5f} deg; true-reg max "
          f"{int(res['nreg_true'].max())}/116 (centre {int(res['nreg_true'][n//2])})", flush=True)
    np.savez(f"{CWT}/trackB_p1_{tag}.npz", axis=axis, x_sky=x_sky, coords=np.array(coords),
             profiled_lnl=jl, registration_count=res["nreg"], n_registered_true=res["nreg_true"],
             qmax=res["qmax"], mapk=res["mapk"], census_idx=census_idx, truth_sky=truth_sky,
             fwhm_deg=fwhm, peak_offset_deg=x_sky[i0], placement="truth-placed")
    return dict(x_sky=x_sky, jl=jl, fwhm=fwhm, peak_off=x_sky[i0], axname=axname)


def run_production(P):
    prior_col = P["priors"]["lit"]
    census_idx = np.array([P["names"].index(nm) for nm in CENSUS])
    tt = P["theta_truth"]; ndist = P["n_dist"]
    cg_t, gp_t = tt[ndist + 0], tt[ndist + 1]      # loud0 TRUTH sky
    truth_sky = np.array([cg_t, gp_t])
    pred = np.degrees(1.0 / (2 * np.pi * _max_pterm_phase(P)))
    print(f"\n  loud0 truth sky: cos_gwtheta={cg_t:.5f} gwphi={gp_t:.5f}", flush=True)
    print(f"  predicted needle width (1-rad pterm phase, max_p 2pi L/dL) ~ {pred:.5f} deg", flush=True)

    # ---- 2-D coarse (truth-centred scan) + zoom1 (recentred on coarse argmax) ----
    a_coarse = run_2d_stage(P, "coarse", cg_t, gp_t, 20.0, N2D, prior_col, census_idx,
                            truth_sky, placement="data-driven (truth-centred wide scan)")
    a_zoom1 = run_2d_stage(P, "zoom1", a_coarse["best_cg"], a_coarse["best_gp"], 2.0, N2D,
                           prior_col, census_idx, truth_sky, placement="data-driven (recentred on coarse argmax)")
    z1_cg, z1_gp = a_zoom1["best_cg"], a_zoom1["best_gp"]
    z1_off_truth = ang_deg(z1_cg, z1_gp, cg_t, gp_t)
    print(f"\n  zoom1 argmax vs truth: {z1_off_truth:.5f} deg apart", flush=True)

    # ---- high-res 1-D cuts through TRUTH: the needle width (truth-placed) ----
    print("\n=== 1-D WIDTH CUTS through TRUTH (truth-placed) ===", flush=True)
    cuts = {}
    for axis in (0, 1):
        an = "cos" if axis == 0 else "gwphi"
        cuts[f"{an}_wide"] = run_cut(P, f"cut_{an}_wide", axis, cg_t, gp_t, 2.0, NCUT,
                                     prior_col, census_idx, truth_sky)
        cuts[f"{an}_fine"] = run_cut(P, f"cut_{an}_fine", axis, cg_t, gp_t, 0.03, NCUT,
                                     prior_col, census_idx, truth_sky)
    fine_w = [cuts["cos_fine"]["fwhm"], cuts["gwphi_fine"]["fwhm"]]
    W = np.nanmax(fine_w) if np.isfinite(np.nanmax(fine_w)) else 0.02   # patch half-width source

    # ---- uniqueness patch about the zoom1 argmax (DATA-DRIVEN centre) ----
    patch_half = float(min(3.0 * W, 0.05)) if np.isfinite(W) else 0.02
    print(f"\n=== UNIQUENESS PATCH [data-driven centre = zoom1 argmax]: "
          f"9x9 +/-{patch_half:.5f} deg (3x fine-cut FWHM {W:.5f}) ===", flush=True)
    a_patch = run_2d_stage(P, "patch", z1_cg, z1_gp, patch_half, N2D, prior_col, census_idx,
                           truth_sky, placement="data-driven (centre=zoom1 argmax; NOT truth-placed)")
    # is the patch max a 2-D interior maximum (not a ridge/saddle edge)?
    pj = a_patch["jl"]; pi = a_patch["best_cell"]
    interior = (0 < pi[0] < N2D - 1) and (0 < pi[1] < N2D - 1)
    patch_max_off_truth = ang_deg(a_patch["best_cg"], a_patch["best_gp"], cg_t, gp_t)
    print(f"  patch argmax interior={interior}; patch argmax vs truth {patch_max_off_truth:.5f} deg; "
          f"zoom1-argmax-vs-truth {z1_off_truth:.5f} deg", flush=True)

    # ---- freq zoom at the (data-driven) best sky ----
    run_freq_zoom(P, prior_col, z1_cg, z1_gp, census_idx)

    # ---- summary ----
    print("\n=== P1 NEEDLE SUMMARY ===", flush=True)
    print(f"  predicted needle width ~ {pred:.5f} deg  (2pi L/dL lever arm)", flush=True)
    for tag, a in [("coarse", a_coarse), ("zoom1", a_zoom1), ("patch", a_patch)]:
        print(f"  {tag:7s} [{'data-driven'}]: argmax {a['off']:.5f} deg from centre | "
              f"true-reg max {int(a['nrt'].max())}/116 | census-3 max {int(a['nrt_census'].max())}/3 "
              f"| {a['dt']:.0f}s{' OVER45' if a['over45'] else ''}", flush=True)
    for k, c in cuts.items():
        print(f"  cut {k:12s} [truth-placed]: FWHM {c['fwhm']:.5f} deg, peak off truth "
              f"{c['peak_off']:+.5f} deg  (pred {pred:.5f})", flush=True)
    print(f"  zoom1 argmax vs truth {z1_off_truth:.5f} deg; patch argmax interior-max={interior}, "
          f"vs truth {patch_max_off_truth:.5f} deg", flush=True)
    return 0


def _max_pterm_phase(P):
    """max_p 2*pi*L_p/dL_p (rad) -- the lever arm setting the needle-width prediction."""
    return float(np.max(2 * np.pi * P["L0"] / P["dL"]))


def run_freq_zoom(P, prior_col, cg, gp, census_idx, half_bins=3.0, n=25):
    """1-D log10_fgw zoom of loud0 at the (data-driven) best sky. Width in freq bins."""
    print(f"\n=== FREQ ZOOM (1-D log10_fgw) at best sky (cg={cg:.5f}, gp={gp:.5f}) ===", flush=True)
    tt = P["theta_truth"].copy(); ndist = P["n_dist"]
    f0 = tt[ndist + 4]                                  # loud0 truth log10_fgw
    T = 6.992e8; df = 1.0 / (3.0 * T)                   # freq-bin scale (1/3T)
    dlf = half_bins * df / (10 ** f0 * np.log(10))      # +/- half_bins in log10 space (approx)
    lfs = f0 + np.linspace(-dlf, dlf, n)
    trials = []
    for lf in lfs:
        th = tt.copy(); th[ndist + 0] = cg; th[ndist + 1] = gp; th[ndist + 4] = lf
        trials.append(th)
    res = eval_stage(P, trials, prior_col)
    finite_check("freqzoom", "profiled_lnl", res["jlnl"])
    jl = res["jlnl"]; i0 = int(np.argmax(jl))
    off_bins = abs(lfs[i0] - f0) * (10 ** f0 * np.log(10)) / df
    y = jl - jl.max(); idx = np.where(y > -0.5)[0]
    fwhm_bins = ((lfs[idx[-1]] - lfs[idx[0]]) * (10 ** f0 * np.log(10)) / df) if len(idx) > 1 else np.nan
    print(f"  freq peak {off_bins:.4f} bins from truth; FWHM(0.5nat) ~ {fwhm_bins:.4f} freq bins; "
          f"true-reg max {int(res['nreg_true'].max())}/116", flush=True)
    np.savez(f"{CWT}/trackB_p1_freqzoom.npz", log10_fgw=lfs, truth_log10_fgw=f0,
             profiled_lnl=jl, registration_count=res["nreg"], n_registered_true=res["nreg_true"],
             qmax=res["qmax"], mapk=res["mapk"], best_sky=np.array([cg, gp]), df=df,
             placement="data-driven-sky (truth-placed freq axis)")
    print(f"  checkpoint -> {CWT}/trackB_p1_freqzoom.npz", flush=True)


# ============================================================
if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "run"])
    ap.add_argument("--config", default="pop", choices=list(TE.CONFIGS))
    a = ap.parse_args()
    print(f"=== P1 registration map ({a.mode}), config={a.config} ===", flush=True)
    t0 = time.time()
    P = TE.build_problem(a.config)
    print(f"  build_problem: {time.time()-t0:.0f}s", flush=True)
    if a.mode == "gate":
        ok, est = run_gate(P)
        sys.exit(0 if ok else 1)
    else:
        ok, est = run_gate(P)                # gate always runs before production
        if not ok:
            print("  GATE FAILED -> refusing production run.", flush=True); sys.exit(1)
        sys.exit(run_production(P))
