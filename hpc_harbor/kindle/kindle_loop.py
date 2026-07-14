"""
KINDLE — the soft-loop GAIN driver, re-instrumented after the g1 diagnosis.

WHY THIS EXISTS (project_progress.md §10.17 open-1; KINDLE g1):
  The brief asked for the "gain>1 contour". g1 (KINDLE_gain_contour.md §1) proved the
  banked gain statistic (ignite2.py:501, a RATIO |C_fix|/|C_0|) is degenerate:
    - every ignition 0->1 divides by |C_0|=0  -> NaN, filtered OUT of the "finite" set
    - every finite value is an n/n tautology  -> 1.000 by arithmetic
    - the count has NEVER gone n->n+1 with n>=1 in 50 honest runs -> gain>1 unreachable
    - 25/30 honest loops accept ZERO M-steps -> the loop is initialised at its own fixed
      point (it starts at the TRUE source, where the gradient is ~0)
  So "gain=1.000" is not marginality; it is a degenerate readout of a loop run at rest.

THE REDESIGN (Matt, this session — recorded as his call per the IGNITE-2 §1.3 precedent):
  * RATIO GAIN RETIRED. Readouts everywhere:
      additive gain  dN = |C_fix| - |C_0|         (well-defined at zero; counts ignitions)
      margin flow    per-pulsar (dlnL - bar) movement toward/away from the bar
  * ARM (a)  start = TRUTH             -> the CONTROL, gated bit-identical vs banked sloop.
  * ARM (b-MAP) start = the criterion's own MAP source (Earth-term-only fit) -> PRIMARY.
      per-realisation start; the realised start-offset is BANKED AS A COLUMN and supplies
      the x-axis (dN & margin-flow vs realised offset).                     [Phase 2 — GPU]
  * ARM (b-FIX) start = truth + 1 x sigma(log10_mc) along the mc eigendirection -> the
      pre-registered CONTROL: "does a KNOWN displacement of KNOWN size get repaired?"
      same seeds as arm (a). ONE magnitude only.
  * SCRAMBLED NULL through every arm (off-truth + no signal is where the hard-lock cascade
      was born; the soft loop's silence there is a deliverable, not an assumption).
  * SANITY: bank whether the start lies inside/outside the certified-seed basin
      (|C_0| at the start vs at truth).

This module reuses ignite2's numerics VERBATIM (make_marg, fd_hessian,
sigma_profile_batch, constants) so the truth arm is bit-identical to the banked loop BY
CONSTRUCTION: the only divergence from ignite2.run_soft_realisation is the start injection,
which is a NO-OP when start_mode=='truth'. g0 verifies this to 0.0e+00.
"""
import os, sys, time
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
sys.path.insert(0, f"{HSYMT}/hpc_harbor/ignite2")
KOUT = f"{HSYMT}/KINDLE_results"

import ignite2 as IG2          # reuse the banked numerics wholesale
from ignite2 import (make_marg, fd_hessian, sigma_profile_batch,
                     DAMP, H_FD, MAX_ITER)


def force_recompute_lgwb(C):
    """Reproduce IGNITE-2's T=30 noise path EXACTLY.

    The canonical b1_L_gwb.npz (RECUT, §10.16(e)) is shaped (3248,3248) for the T=15 ANCHOR
    array. build_ignite_problem(30) needs (5336,5336), so load_or_build_L_gwb hits its
    shape-mismatch raise. But IGNITE-2 ran 2026-07-12, BEFORE that bank existed (RECUT made
    it 2026-07-13) — its T=30 loops used the RECOMPUTE path (np.linalg.eigh of Phi_gwb at
    cpus=8), which §10.14(g)/§10.16(e) state is the convention every banked noisy number was
    drawn under and is bit-identical AT cpus=8. There is NO canonical bank for the T=30
    geometry; recompute-at-8 is the correct and only reproduction path.

    This wrapper points the loader at a non-existent path so it takes the recompute branch
    (the same branch IGNITE-2 took), rather than tripping on the wrong-shaped bank. The g0
    gate verifies the result is bit-identical to the banked trajectory. cpus-per-task MUST be
    8 (the sbatch sets it); the wrapper prints the recompute provenance + thread count so the
    log proves the path taken.
    """
    _orig = C.load_or_build_L_gwb

    def _recompute(Phi_gwb, path, strict=False):
        L, prov = _orig(Phi_gwb, "/nonexistent_kindle_force_recompute.npz", strict=False)
        print(f"[KINDLE] L_gwb via {prov}", flush=True)
        return L, prov

    C.load_or_build_L_gwb = _recompute


def inject_start(theta_cur, x_ref_scaled, sc_sel, sel, mc_dims, sig_mc_scaled,
                 start_mode, start_scale, rng):
    """Return (x_start_scaled, offset_report) for the chosen start convention.

    x_ref_scaled : truth in scaled sel-coords (the arm-a start).
    offset_report: dict of realised per-parameter offset in units of the mc profile width,
                   banked as a column so the readout is dN/margin vs realised offset.
    """
    x = x_ref_scaled.copy()
    if start_mode == "truth":
        return x, dict(mode="truth", off_mc_sig=0.0, off_max_sig=0.0)
    if start_mode == "fix_mc":
        # pre-registered: displace each loud source by start_scale x sigma(log10_mc) along
        # its OWN mc eigendirection (the axis the pilot showed does the registration work).
        for j, d in enumerate(mc_dims):
            x[d] = x_ref_scaled[d] + start_scale * sig_mc_scaled[j]
        off = np.zeros(len(x))
        for j, d in enumerate(mc_dims):
            off[d] = (x[d] - x_ref_scaled[d]) / max(sig_mc_scaled[j], 1e-30)
        return x, dict(mode="fix_mc", off_mc_sig=float(np.max(np.abs(off[mc_dims]))),
                       off_max_sig=float(np.max(np.abs(off))))
    raise ValueError(f"start_mode {start_mode!r} not implemented in kindle_loop "
                     "(map arm is Phase 2 — needs the Earth-term-only optimiser)")


def run_kindle_realisation(jnp, IG, C, B1, TE, F, P, E_cache, h, tier, floor_cell,
                           geo_id, geo_src, noise_seed, dist_seed, scrambled=False,
                           start_mode="truth", start_scale=0.0, max_iter=MAX_ITER,
                           out_dir=KOUT, prefix="ksloop"):
    """One KINDLE soft-loop realisation. Byte-for-byte ignite2.run_soft_realisation EXCEPT:
      (1) the start injection (no-op for start_mode='truth');
      (2) additive-gain + margin-flow + start-offset instrumentation added to the bank.
    Nothing is ever hard-locked; the per-iteration certified set is a SCORE, not a lock."""
    T = P["T_label"]
    tag = (f"{prefix}{'X' if scrambled else ''}_{IG.cell_tag(h, T, tier)}"
           f"_{start_mode}s{start_scale:g}_g{geo_id}_n{noise_seed}")
    path = f"{out_dir}/k_{tag}.npz"
    if os.path.exists(path):
        return tag, "skip"
    t00 = time.time()
    G = IG.cell_geometry(P, geo_src, h, tier)
    dL, EV, FT = G["dL"], G["EV"], G["FT"]
    nd = P["n_dist"]; L0 = P["L0"]; AI = P["AI"]; npsr = P["npsr"]
    one = jnp.ones(npsr); sp = P["sp"]

    L_true, n_true = IG.draw_true_distances_tier(P, dL, EV, seed=dist_seed, tier=tier)
    n_true_grid = C.n_true_on_grid(L_true, L0, dL)
    theta_true = G["theta"].copy(); theta_true[AI] = L_true
    clean = P["amo"]["inject_delay"](jnp.asarray(theta_true), one)
    noise = P["nd"].draw(noise_seed)
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n)) for c, n in zip(clean, noise))

    theta_cur = G["theta"].copy()
    if scrambled:
        theta_cur = F.scramble_source_sky(theta_cur, nd, P["ncw"], seed=noise_seed + 13,
                                          mode="isotropic", loud_only=False)
    theta_truth_src = theta_cur.copy()          # the arm-a start (post-scramble if scrambled)

    sel = np.array([nd + 8 * k + j for k in range(C.N_LOUD) for j in IG.LOOP_PARS])
    sc_sel, Pi = IG.loop_prior(P, sel)
    sig_box = 1.0 / np.sqrt(np.diag(Pi))
    lo_b, hi_b = TE.phi_bounds(P)
    lo_x = lo_b[sel - nd] / sc_sel; hi_x = hi_b[sel - nd] / sc_sel
    D = len(sel)
    mc_dims = [k * len(IG.LOOP_PARS) + 1 for k in range(C.N_LOUD)]

    E = make_marg(B1, P, data, EV, dL, FT, tier, E_cache)
    if E_cache.get("gate_pending", True):
        E.gate_fast_full(theta_cur[nd:][None])
        E_cache["gate_pending"] = False

    src_bg = theta_cur[nd:].copy()

    def lnpost(X_scaled):
        X_scaled = np.atleast_2d(X_scaled)
        srcs = np.tile(src_bg, (X_scaled.shape[0], 1))
        srcs[:, sel - nd] = X_scaled * sc_sel[None, :]
        lm = E.lnmarg(srcs)
        xc = X_scaled - x_ref[None, :]
        return lm - 0.5 * np.sum(xc @ Pi * xc, axis=1)

    x_ref = theta_cur[sel] / sc_sel               # prior centre = the truth (arm-a) start
    x_truth = x_ref.copy()

    # ---- profile width at truth, for the fix_mc offset magnitude ----
    prof0 = sigma_profile_batch(lnpost, x_truth, mc_dims, sig_box)
    sig_mc_scaled = np.array([prof0[d] for d in mc_dims])   # scaled-unit sigma per loud src

    # ---- START INJECTION (the ONLY divergence from ignite2 for start_mode='truth') ----
    rng = np.random.default_rng(noise_seed + 777)
    x, off_report = inject_start(theta_cur, x_ref, sc_sel, sel, mc_dims, sig_mc_scaled,
                                 start_mode, start_scale, rng)

    lnK = np.log(np.maximum(FT.K_counted.astype(float), 1.0))
    bar = np.maximum(lnK, floor_cell)

    def score(theta_src_full):
        tb = theta_src_full.copy(); tb[AI] = L0
        LNL = sp.estep(tb, EV, AI, data, one)
        C._finite(f"kindle estep {tag}", LNL)
        post = FT.posterior(LNL)
        dlnL = np.array([(np.sort(FT.peak[a])[-1] - np.sort(FT.peak[a])[-2])
                         if np.asarray(FT.peak[a]).size > 1 else np.inf
                         for a in range(npsr)])
        det = dlnL > bar
        cert = det & (post["q_max"] > 0.9)
        ptrue, _, _ = FT.p_of_fringe(n_true_grid)
        return dict(dlnL=dlnL, det=det, cert=cert, qmax=post["q_max"],
                    mapk=post["map_k"], ptrue=ptrue)

    # ---- SANITY: |C_0| at the START vs at TRUTH (is the start inside the seed basin?) ----
    th_start = theta_cur.copy(); th_start[sel] = x * sc_sel
    s0 = score(th_start)
    C0_start = int(s0["cert"].sum())
    th_truthstart = theta_cur.copy(); th_truthstart[sel] = x_truth * sc_sel
    C0_truth = int(score(th_truthstart)["cert"].sum())

    hist = dict(n_cert=[], n_cert_true=[], n_detect=[], wrong=[], gain=[],
                W=[], W_grew=[], sig_mc_dex=[], lnmarg=[], step_norm=[],
                src_mc_off=[], stable=[], dN=[], margin_flow=[], n_toward=[], n_away=[])
    cert_sets = []
    id_rows = []
    full = dict(dlnL=[], qmax=[], mapk=[])
    Hpos = None
    prev_ptrue = None
    m_prev = None
    C0 = None

    for it in range(max_iter):
        theta_score = theta_cur.copy()
        theta_score[sel] = x * sc_sel
        s = score(theta_score)
        full["dlnL"].append(np.nan_to_num(s["dlnL"], posinf=1e30))
        full["qmax"].append(s["qmax"]); full["mapk"].append(s["mapk"])
        certs = set(np.where(s["cert"])[0].tolist())
        cert_sets.append(certs)
        wrong = int(np.sum(s["cert"] & (s["mapk"] != n_true_grid)))
        n_true_cert = int(np.sum(s["cert"] & (s["mapk"] == n_true_grid)))
        qf = 1.0 - np.nan_to_num(s["ptrue"], nan=1.0)
        W = float(np.sum(qf))
        W_grew = -1 if prev_ptrue is None else int(np.sum(
            qf > (1.0 - np.nan_to_num(prev_ptrue, nan=1.0)) + 1e-12))
        prev_ptrue = s["ptrue"]
        # RATIO gain kept ONLY for the trail (degenerate — see module docstring)
        gain = (len(cert_sets[-1]) / len(cert_sets[-2])
                if len(cert_sets) > 1 and len(cert_sets[-2]) else np.nan)
        # ADDITIVE gain: |C_it| - |C_0|
        if C0 is None:
            C0 = len(certs)
        dN = len(certs) - C0
        # MARGIN FLOW: movement of (dlnL - bar) since the previous iteration
        m_now = np.nan_to_num(s["dlnL"], posinf=1e30) - bar
        if m_prev is None:
            mflow, n_tw, n_aw = 0.0, 0, 0
        else:
            dm = m_now - m_prev
            mflow = float(np.sum(dm))
            n_tw = int(np.sum(dm > 1e-9)); n_aw = int(np.sum(dm < -1e-9))
        m_prev = m_now
        for a in sorted(certs):
            id_rows.append((it, a, bool(s["mapk"][a] == n_true_grid[a]),
                            float(s["dlnL"][a]), float(s["qmax"][a]),
                            int(s["mapk"][a] - n_true_grid[a])))
        stable = (len(cert_sets) > 1 and certs == cert_sets[-2])

        # ---- soft M-step: damped Newton on lnL_marg + prior (ignite2 verbatim) ----
        g0 = float(lnpost(x[None])[0])
        prof = sigma_profile_batch(lnpost, x, mc_dims, sig_box)
        sig_mc_dex = float(np.median([prof[i] for i in mc_dims]) * F.SCALE_MC)
        if Hpos is None:
            H = fd_hessian(lnpost, x, D, H_FD)
            Hpos = -0.5 * (H + H.T)
            lam, Evec = np.linalg.eigh(Hpos + Pi)
            lam = np.maximum(lam, 1e-6 * np.abs(lam).max())
            A_inv = (Evec / lam) @ Evec.T
        gpts = []
        for i in range(D):
            e = np.zeros(D); e[i] = H_FD
            gpts += [x + e, x - e]
        gv = lnpost(np.array(gpts))
        grad = np.array([(gv[2 * i] - gv[2 * i + 1]) / (2 * H_FD) for i in range(D)])
        step = A_inv @ grad
        cand = np.array([x + DAMP * step, x + 0.5 * DAMP * step, x + 0.25 * DAMP * step])
        gc = lnpost(cand)
        j = int(np.argmax(gc))
        accepted = bool(gc[j] > g0)
        step_norm = float(np.max(np.abs(cand[j] - x))) if accepted else 0.0
        if accepted:
            x = np.clip(cand[j], lo_x, hi_x)

        hist["n_cert"].append(len(certs)); hist["n_cert_true"].append(n_true_cert)
        hist["n_detect"].append(int(s["det"].sum())); hist["wrong"].append(wrong)
        hist["gain"].append(gain); hist["W"].append(W); hist["W_grew"].append(W_grew)
        hist["sig_mc_dex"].append(sig_mc_dex); hist["lnmarg"].append(g0)
        hist["step_norm"].append(step_norm)
        hist["src_mc_off"].append(float(np.max(np.abs(
            (x * sc_sel)[mc_dims] - theta_truth_src[sel][mc_dims]))))
        hist["stable"].append(bool(stable))
        hist["dN"].append(dN); hist["margin_flow"].append(mflow)
        hist["n_toward"].append(n_tw); hist["n_away"].append(n_aw)
        if stable and step_norm < 1e-4:
            break

    theta_fin = theta_cur.copy(); theta_fin[sel] = x * sc_sel
    s_fin = score(theta_fin)
    C_fix = int(s_fin["cert"].sum())
    for k in hist:
        hist[k] = np.array(hist[k], dtype=float)

    os.makedirs(out_dir, exist_ok=True)
    np.savez(path, kind=("ksloopnull" if scrambled else "ksloop"), h=h, T_label=T,
             tier=tier, floor_cell=floor_cell, geo_id=geo_id, noise_seed=noise_seed,
             dist_seed=dist_seed, scrambled=scrambled,
             start_mode=start_mode, start_scale=start_scale,
             off_mc_sig=off_report["off_mc_sig"], off_max_sig=off_report["off_max_sig"],
             sig_mc_scaled=sig_mc_scaled,
             C0_start=C0_start, C0_truth=C0_truth, C_fix=C_fix,
             dN_final=C_fix - C0_start, start_in_basin=bool(C0_start >= 1),
             n_true=n_true, n_true_grid=n_true_grid,
             dlnL_final=s_fin["dlnL"], lnK=lnK, qmax_final=s_fin["qmax"],
             mapk_final=s_fin["mapk"], on_true_final=(s_fin["mapk"] == n_true_grid),
             cert_final=s_fin["cert"], det_final=s_fin["det"],
             ptrue_final=np.nan_to_num(s_fin["ptrue"], nan=-1.0),
             x_final=x * sc_sel, sel=sel, theta0_sel=theta_truth_src[sel],
             id_iter=np.array([r[0] for r in id_rows], int),
             id_psr=np.array([r[1] for r in id_rows], int),
             id_on_true=np.array([r[2] for r in id_rows], bool),
             id_dlnL=np.array([r[3] for r in id_rows]),
             id_qmax=np.array([r[4] for r in id_rows]),
             id_dk=np.array([r[5] for r in id_rows], int),
             n_iter=len(hist["n_cert"]),
             iter_dlnL=np.array(full["dlnL"]), iter_qmax=np.array(full["qmax"]),
             iter_mapk=np.array(full["mapk"], int),
             **{f"traj_{k}": v for k, v in hist.items()},
             names=np.array(P["names"]), census_idx=P["census_idx"],
             anchor_idx=P["anchor_idx"])
    nc = hist["n_cert"]
    return tag, (f"start={start_mode}s{start_scale:g} C0(start/truth)={C0_start}/{C0_truth} "
                 f"n_cert {int(nc[0])}->{int(nc[-1])} dN_fin={C_fix - C0_start:+d} "
                 f"wrong_fin={int(hist['wrong'][-1])} iters={len(nc)} "
                 f"({time.time()-t00:.0f}s)")
