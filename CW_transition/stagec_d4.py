"""Stage C — Deliverable D4: N_CW sweep on the real likelihood.

Sweep N_CW = 1,2,4,8,16 (equal strain, h=-13.75 all sources), 10 Asimov zero-noise
geometry draws each, plus one population draw (3 loud + 13 faint, 10x) at N_CW=16.

Per N_CW track:
  (a) median marginal sigma_L over Fisher-VALID pulsars (marg < dL/4; else the
      within-fringe width is not Gaussian -> "fringe-unresolved", excluded).
  (b) D3 classification counts (direct 512-eval fringe scan; NEVER the prominence
      peak finder). class-(i) [ident AND marg<EM] is the headline.
  (c) marg/cond distribution (does the conditional optimism open up with N_CW?).
  (d) median dlnL_a vs ln K_a (fringe-ID gap: switches on gradually or sharply?).

Reuses the amortised zero-noise Fisher (one compile per N_CW; HVP-assembled full
Hessian so N_CW=16 fits) and fringe eval at fringe centers.

Run in jug-gpu:  python stagec_d4.py --ndraws 10
"""

import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse
import sys
import time

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import (
    load_pulsars, build_enterprise_pta, compute_mode_spacing, _enterprise_to_disco_params,
)
from stagec_fisher import (
    build_fisher_amortised, make_geometry_injection, sigma_L_from_hessian,
    N_COMPONENTS, LOG10_EQUAD, CW_PARAM_BASE, KPC_TO_PC,
)

EQUAL_H = -13.75
NCW_LIST = [1, 2, 4, 8, 16]
N_FRINGE_EVAL = 512
FR_CHUNK = 256
HVP_CHUNK = 8


def assemble_hessian(amo, theta, data):
    """Full (n_theta,n_theta) Hessian via chunked HVPs (no dense OOM)."""
    n = amo["n_theta"]
    eye = np.eye(n)
    cols = np.empty((n, n))
    for s in range(0, n, HVP_CHUNK):
        e = min(s + HVP_CHUNK, n)
        V = eye[s:e]
        if e - s < HVP_CHUNK:
            V = np.vstack([V, np.zeros((HVP_CHUNK - (e - s), n))])
        out = np.array(amo["hvp_batched"](theta, data, jnp.asarray(V)))
        cols[:, s:e] = out[: e - s].T
    return 0.5 * (cols + cols.T)


def fringe_dlnL(amo, theta, data, a_idx, L0, sig_em, dL, n_sigma=3.0):
    """dlnL(true vs best-wrong fringe) for distance index a_idx; direct 512-eval."""
    lo = max(1e-3, L0 - n_sigma * sig_em)
    hi = L0 + n_sigma * sig_em
    kmin = int(np.floor((lo - L0) / dL)); kmax = int(np.ceil((hi - L0) / dL))
    ks = np.arange(kmin, kmax + 1)
    centers = L0 + ks * dL
    keep = (centers >= lo) & (centers <= hi) & (ks != 0)
    wrong, wrong_ks = centers[keep], ks[keep]
    K = len(wrong)
    if K == 0:
        return np.inf, 0
    if K > N_FRINGE_EVAL:
        near = np.abs(wrong_ks) <= N_FRINGE_EVAL // 2
        far = np.where(~near)[0]
        nfar = N_FRINGE_EVAL - int(near.sum())
        far = far[np.linspace(0, len(far) - 1, max(nfar, 0)).astype(int)]
        evalL = wrong[np.concatenate([np.where(near)[0], far])]
    else:
        evalL = np.concatenate([wrong, np.full(N_FRINGE_EVAL - K, wrong[0])])
    evalL = evalL[:N_FRINGE_EVAL]

    ll_truth = float(amo["fisher_logL"](theta, data))
    theta_np = np.array(theta)
    best = -np.inf
    for s in range(0, N_FRINGE_EVAL, FR_CHUNK):
        blk = evalL[s:s + FR_CHUNK]
        tb = np.repeat(theta_np[None, :], FR_CHUNK, axis=0)
        tb[:, a_idx] = blk
        vals = np.array(amo["fisher_logL_batched"](jnp.asarray(tb), data))
        best = max(best, float(vals.max()))
    return ll_truth - best, K


def run_one_draw(amo, ent_psrs, disco_psrs, inj, cwb, em_pc, sig_em_kpc):
    """Return per-pulsar arrays for one draw: cond,marg (pc), dL(kpc), valid, dlnL, K, cls."""
    temp = _enterprise_to_disco_params(inj, cwb)
    theta = jnp.array([temp[k] for k in amo["theta_keys"]], dtype=jnp.float64)
    data = amo["inject_data"](theta)
    n_dist = amo["n_dist"]
    H = assemble_hessian(amo, theta, data)
    negc = int((np.diag(H)[:n_dist] < 0).sum())
    cond, marg, _ = sigma_L_from_hessian(H, n_dist, rcond=1e-10)
    cond_pc, marg_pc = cond * KPC_TO_PC, marg * KPC_TO_PC

    cw_list = [{k: inj[f"{name}_{k}"] for k in CW_PARAM_BASE} for name in cwb]
    npsr = len(ent_psrs)
    dL = np.array([min(compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"],
                                            cw["log10_fgw"], disco_psrs[a].pos)
                       for cw in cw_list) for a in range(npsr)])
    valid = marg_pc < (dL * KPC_TO_PC / 4.0)   # marg_pc < dL_pc/4 -> Gaussian within-fringe
                                                # (was `marg` [kpc] vs dL_pc/4: units bug, 1000x too loose)

    dlnL = np.zeros(npsr); K = np.zeros(npsr); cls = np.zeros(npsr, int)
    idx_of = {k: i for i, k in enumerate(amo["theta_keys"])}
    for a, pe in enumerate(ent_psrs):
        ai = idx_of[f"{pe.name}_cw_p_dist"]
        dd, Ka = fringe_dlnL(amo, theta, data, ai, pe.pdist[0], pe.pdist[1], dL[a])
        Ka = max(1, Ka); dlnL[a] = dd; K[a] = Ka
        ident = dd > np.log(Ka)
        below = marg_pc[a] < em_pc[a]
        cls[a] = 3 if not ident else (1 if below else 2)
    return dict(cond=cond_pc, marg=marg_pc, dL=dL, valid=valid,
                dlnL=dlnL, K=K, cls=cls, negc=negc)


def summarise(res_list, n_dist):
    """Aggregate a list of per-draw dicts."""
    cls = np.array([r["cls"] for r in res_list])
    valid = np.array([r["valid"] for r in res_list])
    marg = np.array([r["marg"] for r in res_list])
    cond = np.array([r["cond"] for r in res_list])
    dlnL = np.array([r["dlnL"] for r in res_list])
    K = np.array([r["K"] for r in res_list])
    # median marginal sigma_L over Fisher-VALID pulsars (per draw, then median)
    med_marg = []
    for r in res_list:
        v = r["valid"]
        med_marg.append(np.median(r["marg"][v]) if v.any() else np.nan)
    ci = (cls == 1).sum(axis=1); cii = (cls == 2).sum(axis=1); ciii = (cls == 3).sum(axis=1)
    nvalid = valid.sum(axis=1)
    ratio = marg / cond
    finite = np.isfinite(dlnL)
    return dict(
        class_i=(int(np.median(ci)), int(ci.min()), int(ci.max())),
        class_ii=(int(np.median(cii)), int(cii.min()), int(cii.max())),
        class_iii=(int(np.median(ciii)), int(ciii.min()), int(ciii.max())),
        n_valid=(int(np.median(nvalid)), int(nvalid.min()), int(nvalid.max())),
        med_marg_valid=float(np.nanmedian(med_marg)),
        ratio_med=float(np.median(ratio)), ratio_84=float(np.percentile(ratio, 84)),
        ratio_max=float(ratio.max()),
        dlnL_med=float(np.median(dlnL[finite])), dlnL_84=float(np.percentile(dlnL[finite], 84)),
        dlnL_max=float(dlnL[finite].max()),
        lnK_med=float(np.median(np.log(K))),
    )


def run(ndraws):
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)
    ent_psrs, disco_psrs = load_pulsars(116)
    em_pc = np.array([p.pdist[1] for p in ent_psrs]) * KPC_TO_PC
    sig_em_kpc = np.array([p.pdist[1] for p in ent_psrs])
    names = [p.name for p in ent_psrs]

    rows = {}
    pop_row = None
    for ncw in NCW_LIST:
        print(f"\n===== N_CW={ncw} (equal strain h={EQUAL_H}) =====", flush=True)
        pta, cwb, _ = build_enterprise_pta(ent_psrs, ncw, components=N_COMPONENTS,
                                           log10_equad=LOG10_EQUAD)
        t0 = time.time()
        inj0 = make_geometry_injection(pta, ent_psrs, ncw, cwb, seed=1000,
                                       h_range=(EQUAL_H, EQUAL_H))
        amo = build_fisher_amortised(disco_psrs, ncw, inj0, cwb)
        # compile HVP + batched logL on the ref draw
        th0 = amo["theta_truth"]; dt0 = amo["inject_data"](th0)
        _ = assemble_hessian(amo, th0, dt0)
        _ = fringe_dlnL(amo, th0, dt0, 0, ent_psrs[0].pdist[0], ent_psrs[0].pdist[1],
                        max(compute_mode_spacing(inj0[f"{cwb[0]}_cos_gwtheta"],
                            inj0[f"{cwb[0]}_gwphi"], inj0[f"{cwb[0]}_log10_fgw"],
                            disco_psrs[0].pos), 1e-6))
        print(f"  build+compile {time.time()-t0:.1f}s", flush=True)

        res_list = []
        for d in range(ndraws):
            td = time.time()
            inj = make_geometry_injection(pta, ent_psrs, ncw, cwb, seed=2000 + d,
                                          h_range=(EQUAL_H, EQUAL_H))
            r = run_one_draw(amo, ent_psrs, disco_psrs, inj, cwb, em_pc, sig_em_kpc)
            res_list.append(r)
            print(f"  draw {d}: {time.time()-td:.1f}s negc={r['negc']}/{amo['n_dist']} "
                  f"class(i)={int((r['cls']==1).sum())} valid={int(r['valid'].sum())} "
                  f"med_marg_valid={np.median(r['marg'][r['valid']]) if r['valid'].any() else np.nan:.4g}pc",
                  flush=True)
        rows[ncw] = summarise(res_list, amo["n_dist"])

        # population contrast at N_CW=16
        if ncw == 16:
            print("  -- population draw (3 loud + 13 faint, 10x) --", flush=True)
            injp = make_geometry_injection(pta, ent_psrs, 16, cwb, seed=3000,
                                           population=(3, -13.25, -14.25))
            rp = run_one_draw(amo, ent_psrs, disco_psrs, injp, cwb, em_pc, sig_em_kpc)
            pop_row = dict(
                class_i=int((rp["cls"] == 1).sum()), class_ii=int((rp["cls"] == 2).sum()),
                class_iii=int((rp["cls"] == 3).sum()), n_valid=int(rp["valid"].sum()),
                med_marg_valid=float(np.median(rp["marg"][rp["valid"]]) if rp["valid"].any() else np.nan),
                ratio_med=float(np.median(rp["marg"] / rp["cond"])),
                dlnL_med=float(np.median(rp["dlnL"][np.isfinite(rp["dlnL"])])),
            )

    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_d4_results.npz",
             ncw_list=np.array(NCW_LIST),
             **{f"row_{n}": np.array(list(rows[n].values()), dtype=object) for n in NCW_LIST})

    # ---- report ----
    print("\n================ D4 SUMMARY (equal strain) ================")
    hdr = f"{'N_CW':>4s} {'class_i':>14s} {'class_iii':>14s} {'n_valid':>12s} " \
          f"{'medMarg_pc':>10s} {'marg/cond med(84,max)':>22s} {'dlnL med(84,max)':>20s} {'lnK_med':>7s}"
    print(hdr)
    for n in NCW_LIST:
        r = rows[n]
        print(f"{n:>4d} {str(r['class_i']):>14s} {str(r['class_iii']):>14s} "
              f"{str(r['n_valid']):>12s} {r['med_marg_valid']:>10.4g} "
              f"{r['ratio_med']:.2f}({r['ratio_84']:.1f},{r['ratio_max']:.0f}){'':>4s} "
              f"{r['dlnL_med']:.3f}({r['dlnL_84']:.2f},{r['dlnL_max']:.1f}){'':>2s} {r['lnK_med']:>7.2f}")
    if pop_row:
        print("\npopulation draw N_CW=16 (3 loud+13 faint, 10x):")
        print(f"  class_i={pop_row['class_i']} class_ii={pop_row['class_ii']} "
              f"class_iii={pop_row['class_iii']} n_valid={pop_row['n_valid']} "
              f"med_marg_valid={pop_row['med_marg_valid']:.4g}pc "
              f"marg/cond_med={pop_row['ratio_med']:.2f} dlnL_med={pop_row['dlnL_med']:.3f}")
    print("saved stagec_d4_results.npz")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ndraws", type=int, default=10)
    sys.exit(run(ap.parse_args().ndraws))
