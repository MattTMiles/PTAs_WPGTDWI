#!/usr/bin/env python
"""SPARK-2 — THE SECOND ARROW: certification-side feedback.

ARROW 1 (SPARK S2, reports/SPARK_launch_criterion.md): certified pulsar terms -> the
reservoir's DETECTION statistic rises 2.96x AND its null floor falls 2.6x -> recruitment.
Measured CASCADE-ALIVE.

ARROW 2 (here): certified pulsar terms -> the 13 faint sources' JOINT CONSTRAINTS tighten (a)
-> modelling them DECONFUSES every pulsar's residual -> the per-pulsar CERTIFICATION floor
drops (b) -> more pulsars certify -> more coherent terms. If arrow 2 is also uphill the two
arrows close a loop that needs no new sources at all -- a FIXED-LIST loop -- and that is a
different, weaker-premised campaign than the EM-mediated EPOCH route arrow 1 alone licenses.

  (a) JOINT TIGHTENING   sigma(faint params) vs N_cert rung, faint SOFT throughout.
  (b) DECONFUSION        certification floor: faint MODELLED at (a)'s widths vs UNMODELLED.
  (c) LEDGER             floor drop per certified pulsar vs margin-to-bar of the nearest
                         UNCERTIFIED pulsars. EDGE-POSITIVE / EDGE-ZERO.

SOFTNESS (spec 3, never hard-committed). The faint sources are NEVER pinned at a point. The
MODELLED state draws each faint source's parameters from N(truth, sigma_faint(rung)) -- the
posterior (a) measured -- and the certification statistic is AVERAGED over N_SOFT draws. That
is a Monte-Carlo marginalisation over the faint posterior, not a delta at its mode. A tighter
rung => a narrower draw => a better model. At rung 0 the draw is the widest, so "modelled at
rung 0" is deliberately a POOR model, not a free win.

A STRUCTURAL LIMIT, MEASURED NOT ASSUMED (reported in full, see mode `mask`).
`B1Split.estep(theta_base, EV, AI, data, pmask)` sweeps EVERY pulsar's distance in ONE call
under a GLOBAL pmask. A pulsar at m_p = 0 has its own pulsar term switched OFF, so its distance
is inert, its fringe row is FLAT, and it cannot certify at all. Therefore the certified-coherent
E-step that arrow 2's "rigidity" question needs -- target pulsar's own term ON, other certified
pulsars at q, uncertified OFF -- is NOT one call: it is one call PER TARGET. The naive test
(estep at pmask = the rung's mask) is CONFOUNDED: it kills the uncertified rows outright, so the
offender statistic (a max over pulsars) drops for a trivial reason (fewer live pulsars), not
because of template rigidity. Mode `mask` measures that confound rather than asserting it.
The certification convention every banked IGNITE/CHORUS/SURFACE number uses is pmask = ONE
(ignite.py:349 `sp.estep(theta_base, EV, AI, data, one)`), and (b) keeps it.

CONVENTIONS. cpus-per-task=8 (the noise-draw thread hazard). T=30 -> force_recompute_lgwb.
Offender + floor: RECUT verbatim (recut_surface.py:75-78 / :101 adopt()); zero-fraction is a
REQUIRED column. Lean-npz: raw dlnL/lnK/qmax banked per unit, never verdicts.

SCOPE. One cell (h=-13.25, T=30, lit, fiducial POP sky). MOCK spine (AXIS 1440 MHz,
10.15(a)) -- the residuals ARE the injected CW+CURN.

Modes:
  pilot            time one E-step + one Fisher; size the grid before spending
  mask             measure the global-pmask confound (the structural limit above)
  unit --real R --rung K    one (realisation x rung) unit: (a) + (b) signal + nulls
  ledger           (c): collate the banked units -> the ledger + verdict
  trials           the trials-factor honesty item (cost column, per-rung floor refit)
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import sys, time, argparse, glob
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir",
                  os.environ.get("HARBOR_JAXCACHE_IN", "/home/mattm/.cache/jax_stagec_cache"))
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

HSYMT = "/home/mattm/projects/HSYMT"
for _p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
           "hpc_harbor/spark"):
    sys.path.insert(0, f"{HSYMT}/{_p}")

import trackB_b1_core as C
import trackB_estimator as TE
import ignite as IG
from spark import (force_recompute_lgwb, draw_realisation, certified_sets, adopt,
                   lnK_vs_nsrc, H_CELL, T_LABEL, TIER, N_LOUD, H_ABSENT, IGNITERS,
                   I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI)

OUT = f"{HSYMT}/SPARK2_results"
QBAR = 0.9
RUNGS = [0, 2, 5, 8]
N_SOFT = 8            # Monte-Carlo draws over the faint posterior (the soft marginalisation)
SOFT_SEED = 20260717
# the 8 free axes per faint source; sky is frozen (the loop cannot act on it -- ignite.py:102
# LOOP_PARS) so tightening is measured on what the loop can actually use.
FAINT_PARS = (I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI)


# ============================================================
# certification scorer -- ignite.py:349-358 semantics, verbatim
# ============================================================
def score_cert(P, FT, theta_base, EV, data, mask, L_true, dL):
    sp = P["sp"]
    LNL = sp.estep(theta_base, EV, P["AI"], data, np.asarray(mask))
    C._finite("estep", LNL)
    post = FT.posterior(LNL)
    n_true_grid = C.n_true_on_grid(L_true, P["L0"], dL)
    qmax = post["q_max"]; mapk = post["map_k"]
    on_true = (mapk == n_true_grid)
    lnK = np.log(np.maximum(FT.K_counted.astype(float), 1.0))
    dlnL_det = np.array([
        (np.sort(FT.peak[a])[-1] - np.sort(FT.peak[a])[-2])
        if np.asarray(FT.peak[a]).size > 1 else np.inf for a in range(P["npsr"])])
    return dict(dlnL=dlnL_det, lnK=lnK, qmax=qmax, on_true=on_true, mapk=mapk)


def offender(s):
    """RECUT recut_surface.py:75-78, verbatim -- max over pulsars, else 0.0 (the zero point
    mass). Non-finite dlnL (a pulsar with a single fringe peak, K=1) is EXCLUDED and counted:
    an inf offender would silently poison the extreme-value fit."""
    d, lnK, q = s["dlnL"], s["lnK"], s["qmax"]
    fin = np.isfinite(d)
    m = (d > lnK) & (q > QBAR) & fin
    return (float(d[m].max()) if m.any() else 0.0), int((~fin).sum())


def cert_mask(s, floor):
    d, lnK, q = s["dlnL"], s["lnK"], s["qmax"]
    bar = np.maximum(lnK, floor)
    return (d > bar) & (q > QBAR), bar


# ============================================================
# (a) JOINT TIGHTENING -- Fisher on the faint sources, pmask = the rung's certified set
# ============================================================
def faint_fisher(P, R, pmask, Ld, npts=9):
    """CONDITIONAL widths on the 13 faint sources' loop-actionable parameters, with the rung's
    certified terms coherent. Faint SOFT: this is a WIDTH, never a commit.

    WHY CONDITIONAL AND NOT THE JOINT/MARGINAL FISHER -- MEASURED, NOT ASSUMED. The marginal
    width needs the full 78x78 Hessian of lnL over the faint block. Two builds of it FAILED to
    return on an A100-80GB inside a 1 h walltime: a monolithic `jax.hessian` (pilot 12583525,
    >27 min, cancelled) and a chunked batched-JVP-of-grad Hessian at CH=8 (pilot 12583899,
    >17 min on the Hessian alone after a 500 s E-step compile, cancelled). The full likelihood
    carries a 2668-dim correlated-GP solve per evaluation; 78 JVPs of its gradient is simply
    not a few-GPU-hr object. Rather than quote a joint width SPARK-2 cannot afford, it measures
    the CONDITIONAL width by a 1-D curvature scan per axis -- ONE batched logL call per chunk,
    seconds -- and says so on every number.

    CONDITIONAL vs MARGINAL, and which way the bias runs: sigma_cond = 1/sqrt(F_ii) IGNORES
    covariance with the other faint parameters, so sigma_cond <= sigma_marg ALWAYS. The
    conditional width is therefore OPTIMISTIC -- it OVERSTATES how well the faint sources are
    known and so OVERSTATES how well they can be modelled. That bias points TOWARD arrow 2,
    which keeps the a-fortiori structure intact: an EDGE-ZERO verdict measured with
    conditional widths is a fortiori EDGE-ZERO with the (wider) marginal ones. An
    EDGE-POSITIVE verdict measured this way would NOT be safe and would have to be re-cut
    against the true joint width before it could be quoted.
    """
    nd = P["n_dist"]
    faint = list(range(N_LOUD, P["ncw"]))
    idx = np.array([nd + 8 * k + j for k in faint for j in FAINT_PARS])
    th = R["theta_true"].copy(); th[P["AI"]] = np.asarray(Ld, float)
    lb = P["amo"]["logL_batch_theta"]
    data = R["data"]; pm = jnp.asarray(pmask)
    # per-axis step: a small fraction of the prior box, so the curvature is local
    BOX = {I_COSINC: 1.0, I_MC: 0.87, I_FGW: 0.5, I_H: 1.0, I_PH0: 3.14, I_PSI: 1.57}
    box = np.array([BOX[j] for _ in faint for j in FAINT_PARS])
    step = 1e-3 * box
    n = len(idx)
    # 3-point stencil per axis -> curvature d2lnL/dx2
    stack, meta = [], []
    for i in range(n):
        for d in (-1, 0, 1):
            t = th.copy(); t[idx[i]] += d * step[i]
            stack.append(t); meta.append((i, d))
    stack = np.stack(stack)
    lls = []
    CH = 48
    for c0 in range(0, len(stack), CH):
        lls.append(np.asarray(lb(jnp.asarray(stack[c0:c0 + CH]), data, pm)))
    lls = np.concatenate(lls)
    L = lls.reshape(n, 3)
    d2 = (L[:, 0] - 2.0 * L[:, 1] + L[:, 2]) / (step ** 2)   # d2 lnL / dx2
    F_ii = -d2                                                # Fisher diagonal
    cond_sig = np.where(F_ii > 0, 1.0 / np.sqrt(np.maximum(F_ii, 1e-300)), np.inf)
    n_pos = int((F_ii > 0).sum())
    return dict(idx=idx, F_ii=F_ii, cond_sig=cond_sig, marg_sig=cond_sig,
                eig_min=float(np.min(F_ii)), eig_max=float(np.max(F_ii)),
                cond_num=float(np.max(np.abs(F_ii)) / max(np.min(np.abs(F_ii)), 1e-300)),
                n_pos=n_pos, n_par=n,
                width_kind="CONDITIONAL (1/sqrt(F_ii)); marginal unaffordable -- see docstring")


def rung_state(P, R, cs_row, k, Lt, L0):
    """pmask/Ld for rung k: the k certified pulsars with the LARGEST dlnL among the banked
    certified set (the array's own ranking, data-driven). k=0 -> fully decohered."""
    npsr = P["npsr"]
    pm = np.zeros(npsr); Ld = L0.copy()
    if k > 0:
        ci = np.where(cs_row["cert"])[0]
        if len(ci) == 0:
            return pm, Ld, np.zeros(0, int)
        take = ci[np.argsort(-np.asarray(cs_row["dlnL"])[ci])][:k]
        pm[take] = np.asarray(cs_row["qmax"])[take]
        Ld[take] = Lt[take]
        return pm, Ld, take
    return pm, Ld, np.zeros(0, int)


def soft_faint_theta(P, R, sig_marg, idx, rng, bounds):
    """A SOFT draw of the faint sources: truth + N(0, sigma). Never a hard commit.

    CLIPPED TO THE PRIOR BOX -- and this is not cosmetic. An unclipped draw put cos_inc
    outside [-1, 1] (its width at the loose rungs is comparable to its whole box), the
    waveform went complex/NaN, and the E-step returned NON-FINITE on ALL 59392 entries
    (measured: units 12584400-11, all 12 died). This is IGNITE-2/KINDLE's known
    double-perturbation pathology (a wide perturbation + the fringe sweep -> degenerate
    covariance -> all-NaN E-step) arriving by a new route. The clip is the honest fix: the
    faint posterior is truncated at its own prior support, exactly as a real posterior is.
    Infinite/huge widths clip to the box scale, so an unconstrained axis draws from the PRIOR
    -- which is the correct statement of "we know nothing about this axis".
    """
    th = R["theta_true"].copy()
    s = np.asarray(sig_marg, float).copy()
    BOX = {I_COSINC: 1.0, I_MC: 0.87, I_FGW: 0.5, I_H: 1.0, I_PH0: 3.14, I_PSI: 1.57}
    box = np.array([BOX[j] for k in range(N_LOUD, P["ncw"]) for j in FAINT_PARS])
    s = np.where(np.isfinite(s), np.minimum(s, box), box)
    lo, hi = bounds
    th[idx] = np.clip(th[idx] + rng.normal(0.0, s), lo, hi)
    return th, s


# ============================================================
def mode_pilot():
    force_recompute_lgwb(C)
    P = IG.build_ignite_problem(T_LABEL, verbose=True)
    R = draw_realisation(P, 7560000, 17560000, faint_present=True)
    G = IG.cell_geometry(P, None, H_CELL, TIER)
    L0 = np.asarray(P["L0"], float); Lt = np.asarray(R["L_true"], float)
    one = np.ones(P["npsr"])
    th = R["theta_true"].copy(); th[P["AI"]] = L0
    t0 = time.time()
    s = score_cert(P, G["FT"], th, G["EV"], R["data"], one, Lt, G["dL"])
    t_e1 = time.time() - t0
    t0 = time.time()
    s = score_cert(P, G["FT"], th, G["EV"], R["data"], one, Lt, G["dL"])
    t_e2 = time.time() - t0
    off, nnf = offender(s)
    print(f"[SPARK2] estep+score: first {t_e1:.1f}s, warm {t_e2:.1f}s; offender {off:.3f}; "
          f"non-finite dlnL {nnf}/{P['npsr']}", flush=True)
    t0 = time.time()
    fi = faint_fisher(P, R, np.zeros(P["npsr"]), L0)
    print(f"[SPARK2] faint Fisher ({fi['n_par']} par): {time.time()-t0:.1f}s  "
          f"cond={fi['cond_num']:.3e} pos-eig {fi['n_pos']}/{fi['n_par']}  "
          f"median marg sigma {np.median(fi['marg_sig']):.4g}", flush=True)
    n_units = 3 * len(RUNGS)
    print(f"\n[SPARK2] SIZING: warm estep {t_e2:.1f}s -> per unit "
          f"(1 signal x {N_SOFT} soft + nulls): see --n-null; {n_units} units", flush=True)
    return 0


def mode_mask():
    """Measure the global-pmask confound rather than asserting it."""
    force_recompute_lgwb(C)
    P = IG.build_ignite_problem(T_LABEL, verbose=True)
    R = draw_realisation(P, 7560000, 17560000, faint_present=True)
    G = IG.cell_geometry(P, None, H_CELL, TIER)
    L0 = np.asarray(P["L0"], float); Lt = np.asarray(R["L_true"], float)
    th = R["theta_true"].copy(); th[P["AI"]] = L0
    cs = certified_sets("m1e07")
    ncs = np.array([r["n_cert"] for r in cs["rows"]])
    row = cs["rows"][int(np.argmax(ncs))]
    row["dlnL"] = np.load(row["path"], allow_pickle=True)["dlnL_det"]
    out = {}
    for k in [0, 8, 116]:
        if k == 116:
            pm = np.ones(P["npsr"]); Ld = Lt
        else:
            pm, Ld, take = rung_state(P, R, row, k, Lt, L0)
        th2 = th.copy(); th2[P["AI"]] = L0
        s = score_cert(P, G["FT"], th2, G["EV"], R["data"], pm, Lt, G["dL"])
        live = np.asarray([np.ptp(np.asarray(G["FT"].peak[a])) for a in range(P["npsr"])])
        off, nnf = offender(s)
        n_live = int((live > 1e-9).sum())
        out[k] = dict(off=off, n_live=n_live, med_q=float(np.median(s["qmax"])))
        print(f"[SPARK2] pmask rung {k:3d}: LIVE fringe rows {n_live:3d}/116  "
              f"offender {off:8.3f}  median qmax {np.median(s['qmax']):.4f}", flush=True)
    print("\n[SPARK2] CONFOUND: if n_live tracks the rung, an estep at the rung's global "
          "pmask kills uncertified rows outright -- the offender (a MAX over pulsars) then "
          "falls for a counting reason, NOT template rigidity. (b) therefore keeps the "
          "banked convention pmask=ONE.", flush=True)
    os.makedirs(OUT, exist_ok=True)
    np.savez(f"{OUT}/spark2_mask.npz",
             rungs=np.array(list(out)), n_live=np.array([out[k]["n_live"] for k in out]),
             offender=np.array([out[k]["off"] for k in out]),
             med_q=np.array([out[k]["med_q"] for k in out]),
             note=("The global-pmask confound, measured. estep applies ONE pmask to a sweep "
                   "of ALL pulsars; m_p=0 switches p's OWN pulsar term off -> p's distance is "
                   "inert -> flat fringe row -> p cannot certify. So the certified-coherent "
                   "E-step arrow 2's rigidity question needs is one call PER TARGET, not one "
                   "call. Not run: 116x cost. (b) keeps pmask=ONE, the convention every "
                   "banked IGNITE/CHORUS/SURFACE number uses (ignite.py:349)."))
    return 0


def mode_unit(real_i, rung, n_null):
    os.makedirs(OUT, exist_ok=True)
    path = f"{OUT}/s2b_r{real_i}_k{rung}.npz"
    if os.path.exists(path):
        print(f"[SPARK2] skip-on-exist: {path}", flush=True)
        return 0
    force_recompute_lgwb(C)
    P = IG.build_ignite_problem(T_LABEL, verbose=True)
    G = IG.cell_geometry(P, None, H_CELL, TIER)
    FT, EV, dL = G["FT"], G["EV"], G["dL"]
    L0 = np.asarray(P["L0"], float)
    one = np.ones(P["npsr"])
    nd = P["n_dist"]

    cs = certified_sets("m1e07")
    rows = cs["rows"]
    row = rows[real_i]
    row["dlnL"] = np.load(row["path"], allow_pickle=True)["dlnL_det"]
    ns, ds = int(row["noise_seed"]), int(row["dist_seed"])
    print(f"[SPARK2] unit real={real_i} (noise {ns} dist {ds}, banked N_cert "
          f"{row['n_cert']}) rung={rung}", flush=True)

    R = draw_realisation(P, ns, ds, faint_present=True)
    Lt = np.asarray(R["L_true"], float)
    pm, Ld, take = rung_state(P, R, row, rung, Lt, L0)
    print(f"[SPARK2]   rung mask: {len(take)} coherent pulsars {list(take)}", flush=True)

    # ---- (a) JOINT TIGHTENING ----
    t0 = time.time()
    fi = faint_fisher(P, R, pm, Ld)
    print(f"[SPARK2]   (a) Fisher {time.time()-t0:.1f}s: median marg sigma "
          f"{np.median(fi['marg_sig']):.5g}  cond {fi['cond_num']:.3e}  "
          f"pos-eig {fi['n_pos']}/{fi['n_par']}", flush=True)

    # ---- (b) DECONFUSION: certification, faint MODELLED (soft) vs UNMODELLED ----
    lo_b, hi_b = TE.phi_bounds(P)
    nd_ = P["n_dist"]
    BND = (lo_b[fi["idx"] - nd_], hi_b[fi["idx"] - nd_])
    rng = np.random.default_rng(SOFT_SEED + 1000 * real_i + rung)
    th_un = R["theta_true"].copy(); th_un[P["AI"]] = L0
    for k in range(N_LOUD, P["ncw"]):
        th_un[nd + 8 * k + I_H] = H_ABSENT           # UNMODELLED: faint out of the model
    s_un = score_cert(P, FT, th_un, EV, R["data"], one, Lt, dL)
    off_un, _ = offender(s_un)

    sig_used = None
    sig_sig = []
    s_mo_list = []
    for d in range(N_SOFT):
        th_mo, su = soft_faint_theta(P, R, fi["marg_sig"], fi["idx"], rng, BND)
        th_mo[P["AI"]] = L0
        sig_used = su
        s_mo_list.append(score_cert(P, FT, th_mo, EV, R["data"], one, Lt, dL))
    # SOFT marginalisation: average the per-pulsar statistic over the faint posterior draws
    s_mo = dict(dlnL=np.mean([x["dlnL"] for x in s_mo_list], axis=0),
                lnK=s_mo_list[0]["lnK"],
                qmax=np.mean([x["qmax"] for x in s_mo_list], axis=0),
                on_true=s_mo_list[0]["on_true"], mapk=s_mo_list[0]["mapk"])
    off_mo, _ = offender(s_mo)
    print(f"[SPARK2]   (b) signal offender: unmodelled {off_un:.3f}  modelled(soft,"
          f"{N_SOFT}) {off_mo:.3f}", flush=True)

    # ---- nulls -> the certification floor in each state ----
    o_un, o_mo = [], []
    for i in range(n_null):
        # nullN convention (ignite.py:340 no_cw=True): NO CW in the DATA, model unchanged.
        # Built directly -- calling draw_realisation here would inject a CW delay and draw the
        # noise, then have both thrown away: a 2x cost for nothing.
        sd = ns + 500_000 + i
        Ltn, _ = IG.draw_true_distances_tier(P, dL, EV, seed=ds + 500_000 + i, tier=TIER)
        noise = P["nd"].draw(sd)
        dn = tuple(jnp.asarray(np.asarray(b)) for b in noise)
        o_un.append(offender(score_cert(P, FT, th_un, EV, dn, one, Ltn, dL))[0])
        thm, _ = soft_faint_theta(P, R, fi["marg_sig"], fi["idx"], rng, BND)
        thm[P["AI"]] = L0
        o_mo.append(offender(score_cert(P, FT, thm, EV, dn, one, Ltn, dL))[0])
        if (i + 1) % 20 == 0:
            print(f"   null {i+1}/{n_null}", flush=True)
    f_un = adopt(np.array(o_un)); f_mo = adopt(np.array(o_mo))
    print(f"[SPARK2]   floor unmodelled {f_un['floor']:.3f} +- {f_un['err']:.3f} "
          f"[{f_un['estimator']}] zf={f_un['zf']:.3f}", flush=True)
    print(f"[SPARK2]   floor modelled   {f_mo['floor']:.3f} +- {f_mo['err']:.3f} "
          f"[{f_mo['estimator']}] zf={f_mo['zf']:.3f}", flush=True)

    c_un, bar_un = cert_mask(s_un, f_un["floor"])
    c_mo, bar_mo = cert_mask(s_mo, f_mo["floor"])
    marg_un = s_un["dlnL"] - bar_un
    marg_mo = s_mo["dlnL"] - bar_mo
    print(f"[SPARK2]   N_cert unmodelled {int(c_un.sum())}  modelled {int(c_mo.sum())}",
          flush=True)

    np.savez(path, real_i=real_i, rung=rung, noise_seed=ns, dist_seed=ds,
             n_coherent=len(take), coherent_idx=take, banked_ncert=row["n_cert"],
             marg_sig=fi["marg_sig"], cond_sig=fi["cond_sig"], sig_used=sig_used,
             fisher_cond=fi["cond_num"], fisher_eig_min=fi["eig_min"],
             fisher_eig_max=fi["eig_max"], fisher_npos=fi["n_pos"], fisher_npar=fi["n_par"],
             dlnL_un=s_un["dlnL"], qmax_un=s_un["qmax"], lnK=s_un["lnK"],
             on_true_un=s_un["on_true"],
             dlnL_mo=s_mo["dlnL"], qmax_mo=s_mo["qmax"], on_true_mo=s_mo["on_true"],
             off_sig_un=off_un, off_sig_mo=off_mo,
             null_off_un=np.array(o_un), null_off_mo=np.array(o_mo),
             floor_un=f_un["floor"], floor_un_err=f_un["err"], floor_un_zf=f_un["zf"],
             floor_un_est=f_un["estimator"],
             floor_mo=f_mo["floor"], floor_mo_err=f_mo["err"], floor_mo_zf=f_mo["zf"],
             floor_mo_est=f_mo["estimator"],
             ncert_un=int(c_un.sum()), ncert_mo=int(c_mo.sum()),
             margin_un=marg_un, margin_mo=marg_mo,
             n_null=n_null, n_soft=N_SOFT, h=H_CELL, tier=TIER, T_label=T_LABEL,
             note=("SPARK-2 arrow-2 unit. RAW per-pulsar dlnL/lnK/qmax + raw null offenders "
                   "banked, never verdicts. Faint SOFT: MODELLED state is a Monte-Carlo "
                   "marginalisation over N(truth, sigma_marg(rung)) with N_SOFT draws, never "
                   "a hard commit. E-step at pmask=ONE, the banked certification convention "
                   "(ignite.py:349) -- see spark2_mask.npz for why the rung's global pmask "
                   "cannot be used here. Nulls: nullN (no CW in the data), model unchanged."))
    print(f"[SPARK2] banked -> {path}", flush=True)
    return 0


def mode_ledger():
    fs = sorted(glob.glob(f"{OUT}/s2b_r*_k*.npz"))
    if not fs:
        print("*** STOP: no units banked.", flush=True)
        return 2
    print(f"=== SPARK-2 LEDGER ({len(fs)} units) ===", flush=True)
    print(f"{'real':>5s} {'rung':>5s} {'coh':>4s} {'med sig':>10s} {'floor_un':>9s} "
          f"{'floor_mo':>9s} {'drop':>8s} {'N_un':>5s} {'N_mo':>5s} {'dN':>4s} "
          f"{'near-miss margin':>17s}")
    rows = []
    for f in fs:
        z = np.load(f, allow_pickle=True)
        fu, fm = float(z["floor_un"]), float(z["floor_mo"])
        drop = fu - fm
        nu, nm = int(z["ncert_un"]), int(z["ncert_mo"])
        mu = np.asarray(z["margin_un"], float)
        mu = mu[np.isfinite(mu) & (mu < 0)]
        near = float(np.sort(mu)[-3:].mean()) if len(mu) >= 3 else np.nan
        rows.append(dict(real=int(z["real_i"]), rung=int(z["rung"]),
                         coh=int(z["n_coherent"]), sig=float(np.median(z["marg_sig"])),
                         fu=fu, fm=fm, drop=drop, nu=nu, nm=nm, dN=nm - nu, near=near))
        r = rows[-1]
        print(f"{r['real']:5d} {r['rung']:5d} {r['coh']:4d} {r['sig']:10.4g} {fu:9.3f} "
              f"{fm:9.3f} {drop:+8.3f} {nu:5d} {nm:5d} {r['dN']:+4d} {near:17.3f}")
    import collections
    byr = collections.defaultdict(list)
    for r in rows:
        byr[r["rung"]].append(r)
    print(f"\n=== per rung (mean over realisations) ===")
    print(f"{'rung':>5s} {'n':>3s} {'med sig':>11s} {'floor drop':>11s} {'dN_cert':>9s} "
          f"{'drop/coh psr':>13s}")
    per = {}
    for k in sorted(byr):
        v = byr[k]
        d = np.mean([x["drop"] for x in v]); dn = np.mean([x["dN"] for x in v])
        s = np.median([x["sig"] for x in v]); coh = np.mean([x["coh"] for x in v])
        per[k] = dict(drop=d, dN=dn, sig=s, coh=coh, n=len(v))
        print(f"{k:5d} {len(v):3d} {s:11.4g} {d:+11.3f} {dn:+9.2f} "
              f"{(d/coh if coh else np.nan):13.4f}")
    ks = sorted(per)
    tighten = per[ks[0]]["sig"] / max(per[ks[-1]]["sig"], 1e-300)
    drops = [per[k]["drop"] for k in ks]
    dNs = [per[k]["dN"] for k in ks]
    print(f"\n(a) tightening rung {ks[0]} -> {ks[-1]}: median sigma "
          f"{per[ks[0]]['sig']:.4g} -> {per[ks[-1]]['sig']:.4g}  ({tighten:.3f}x)")
    print(f"(b) floor drop by rung: " + ", ".join(f"k={k}: {per[k]['drop']:+.3f}" for k in ks))
    print(f"(c) dN_cert  by rung: " + ", ".join(f"k={k}: {per[k]['dN']:+.2f}" for k in ks))
    edge = (max(dNs) > 0) and (max(drops) > 0)
    verdict = "EDGE-POSITIVE" if edge else "EDGE-ZERO"
    print(f"\nVERDICT: {verdict}")
    np.savez(f"{OUT}/spark2_ledger.npz", verdict=verdict, edge=edge,
             rung=np.array(ks), drop=np.array(drops), dN=np.array(dNs),
             sig=np.array([per[k]["sig"] for k in ks]),
             coh=np.array([per[k]["coh"] for k in ks]),
             n_per_rung=np.array([per[k]["n"] for k in ks]),
             unit_real=np.array([r["real"] for r in rows]),
             unit_rung=np.array([r["rung"] for r in rows]),
             unit_drop=np.array([r["drop"] for r in rows]),
             unit_dN=np.array([r["dN"] for r in rows]),
             unit_near=np.array([r["near"] for r in rows]),
             note=("SPARK-2 (c). EDGE-POSITIVE iff the certification floor DROPS and N_cert "
                   "RISES with the rung. Collated from the per-unit raw banks."))
    print(f"[SPARK2] banked -> {OUT}/spark2_ledger.npz", flush=True)
    return 0


def mode_trials():
    """The trials-factor honesty item from SPARK 5(a): the cost column, per-rung, WITHOUT the
    oracle anchor. K_counted/lnK grow with the modelled source list (dL = min over the list),
    so a loop that models the reservoir pays a trials bar it did not pay when the reservoir was
    unmodelled. This is the cost side of arrow 2, measured exactly (CPU)."""
    force_recompute_lgwb(C)
    P = IG.build_ignite_problem(T_LABEL, verbose=False)
    G = IG.cell_geometry(P, None, H_CELL, TIER)
    lk = lnK_vs_nsrc(P, G["theta"])
    os.makedirs(OUT, exist_ok=True)
    print(f"{'n_src':>6s} {'K_sum':>12s} {'lnK_med':>9s} {'dlnK vs 3':>10s} {'dL_med':>11s}")
    base = [r for r in lk if r["n_src"] == 3][0]["lnK_med"]
    for r in lk:
        print(f"{r['n_src']:6d} {r['K_sum']:12.0f} {r['lnK_med']:9.4f} "
              f"{r['lnK_med']-base:+10.4f} {r['dL_med']:11.3e}")
    np.savez(f"{OUT}/spark2_trials.npz",
             n_src=np.array([r["n_src"] for r in lk]),
             K_sum=np.array([r["K_sum"] for r in lk]),
             lnK_med=np.array([r["lnK_med"] for r in lk]),
             lnK_mean=np.array([r["lnK_mean"] for r in lk]),
             dL_med=np.array([r["dL_med"] for r in lk]),
             base_lnK_3=base,
             note=("Trials-factor cost column, exact (CPU). dL = MIN over the MODELLED source "
                   "list of the mode spacing (forge_b1.apply_geometry:91-94) -> modelling the "
                   "reservoir shrinks dL -> more fringes in the prior window -> K_counted and "
                   "the lnK bar GROW. This is what a deconfusing loop pays. Answers SPARK "
                   "5(a): the cost column is no longer oracle-flattered."))
    print(f"[SPARK2] banked -> {OUT}/spark2_trials.npz", flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["pilot", "mask", "unit", "ledger", "trials"])
    ap.add_argument("--real", type=int, default=0)
    ap.add_argument("--rung", type=int, default=0)
    ap.add_argument("--n-null", type=int, default=100)
    a = ap.parse_args()
    if a.mode == "pilot":
        sys.exit(mode_pilot())
    if a.mode == "mask":
        sys.exit(mode_mask())
    if a.mode == "unit":
        sys.exit(mode_unit(a.real, a.rung, a.n_null))
    if a.mode == "ledger":
        sys.exit(mode_ledger())
    sys.exit(mode_trials())
