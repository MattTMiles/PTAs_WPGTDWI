"""ANCHOR CENSUS — Deliverable A2: reclassification with REAL priors (optimized+amended).

Reruns D3 (N_CW=1) + D4 (N_CW=4,8,16 equal strain + 3-loud/13-faint population)
fringe-identification, REUSING the saved geometry-draw seeds (1000 ref, 2000+d draws,
3000 population) -> controlled comparison vs the feather baseline. Truth distance stays
at the injected (feather) value; only the prior WIDTH sigma_EM changes.

Three prior columns from ONE scan/draw (sigma_L and per-fringe lnL are prior-independent):
  feather ; script-optimal (A1 min-sigma) ; literature (script + published composites:
  J0437 Reardon+2016 0.25pc, J1909 Reardon+2021 3pc).

Classifications (all reported separately):
  ANCHOR-CERTIFIED (flat) : dlnL_a > ln K_a  (K_a = COUNTED wrong fringe centres in the
      +/-3 sigma window; round(6 sigma/dL) kept as a cross-check column). Stricter +3.
  PRIOR-CERTIFIED         : zero wrong fringes in the window (prior alone pins the fringe);
      cross-tagged geometry-wide (dL >> median -> antipodal, low-SNR, spurious).
  CERTIFIED-BAYES         : P(true fringe) > 0.9 (strict 0.99), P = softmax[lnL + ln
      Gaussian-prior-weight] over candidate fringes. Honest small-K criterion (prior tails
      penalise wrong fringes ~ offset^2/2sigma^2).
  IMPROVED                : flat-ident AND marginal sigma_L < sigma_EM(real).

OPTIMIZATIONS:
  - vectorized fringe scan across ALL pulsars and ALL draws (one (D x 116 x B) eval
    tensor, chunked); ~15 GPU calls/N_CW vs ~2300 in the loop.
  - JAX persistent compilation cache (4x ~520s compiles amortize across runs).
  - Hessian/sigma_L computed once per draw, shared across the 3 prior columns.

Run in jug-gpu:
  python stagec_anchor_a2.py gate        # vectorized-vs-loop scan gate (1 draw, N_CW=1)
  python stagec_anchor_a2.py run --ndraws 10
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
# STEP 2c: persistent compilation cache (amortize the 4x ~520s compiles)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import load_pulsars, build_enterprise_pta, compute_mode_spacing, _enterprise_to_disco_params
from stagec_fisher import (build_fisher_amortised, make_geometry_injection,
                           sigma_L_from_hessian, N_COMPONENTS, LOG10_EQUAD,
                           CW_PARAM_BASE, KPC_TO_PC)
from stagec_d4 import assemble_hessian

CWT = "/home/mattm/projects/HSYMT/CW_transition"
EQUAL_H = -13.75
NCW_LIST = [int(x) for x in os.environ.get("A2_NCW", "1,4,8,16").split(",")]
B = 512                    # evals per pulsar
BIG_CHUNK = 4096
DENSE_FRINGE_MAX = 64      # windows with <= this many fringes -> dense uniform grid

LIT_INJECT = {
    "J0437-4715": (0.00025, "Reardon+2016 MNRAS455,1751 orbital-Pbdot kinematic 156.79+/-0.25pc"),
    "J1909-3744": (0.00300, "Reardon+2021 PPTA-DR2 Shklovskii kinematic 1152+/-3pc"),
}


def load_real_priors(names):
    a1 = np.load(CWT + "/anchor_a1_Ktable.npz", allow_pickle=True)
    m = {str(a1["fn"][i]): (a1["sig_opt_pc"][i], a1["sig_feath_pc"][i]) for i in range(len(a1["fn"]))}
    sig_script = np.array([m[n][0] for n in names]) / KPC_TO_PC
    sig_feath = np.array([m[n][1] for n in names]) / KPC_TO_PC
    sig_lit = sig_script.copy()
    for i, n in enumerate(names):
        if n in LIT_INJECT:
            sig_lit[i] = LIT_INJECT[n][0]
    return dict(feather=sig_feath, script=sig_script, lit=sig_lit)


# ---------------- eval grid (fixed length B) ----------------
def eval_grid(L0, dL, R):
    """B eval distances over +/-R about L0: dense uniform if <=DENSE_FRINGE_MAX fringes,
    else fringe centres (D3/D4 selection), padded to length B."""
    lo = max(1e-3, L0 - R); hi = L0 + R
    n_fringe = int((hi - lo) / dL) + 1
    if n_fringe <= DENSE_FRINGE_MAX:
        return np.linspace(lo, hi, B)
    kmin = int(np.floor((lo - L0) / dL)); kmax = int(np.ceil((hi - L0) / dL))
    ks = np.arange(kmin, kmax + 1)
    centers = L0 + ks * dL
    keep = (centers >= lo) & (centers <= hi)
    centers, ks = centers[keep], ks[keep]
    if len(centers) > B:
        near = np.abs(ks) <= B // 2
        far = np.where(~near)[0]
        nfar = B - int(near.sum())
        far = far[np.linspace(0, len(far) - 1, max(nfar, 0)).astype(int)]
        centers = centers[np.concatenate([np.where(near)[0], far])][:B]
    if len(centers) < B:
        centers = np.concatenate([centers, np.full(B - len(centers), L0)])  # pad at truth
    return centers


# ---------------- vectorized scan (all pulsars x all draws) ----------------
def make_vmap2(amo):
    return jax.jit(jax.vmap(amo["fisher_logL_batched"], in_axes=(0, 0)))


def scan_all(vmap2, thetas, data_stack, EV, AI):
    """thetas (D,nth); data_stack: pytree with leading (D,); EV (D,npsr,B); AI (npsr,).
    Returns LNL (D,npsr,B)."""
    D, npsr, Bn = EV.shape
    nth = thetas.shape[1]
    M = npsr * Bn
    theta_all = np.repeat(np.asarray(thetas)[:, None, :], M, axis=1)  # (D,M,nth)
    rows = np.arange(npsr).repeat(Bn)
    cols = AI[rows]
    ari = np.arange(M)
    for d in range(D):
        theta_all[d, ari, cols] = EV[d].ravel()
    LNL = np.empty((D, M))
    for c0 in range(0, M, BIG_CHUNK):
        c1 = min(c0 + BIG_CHUNK, M)
        out = np.asarray(vmap2(jnp.asarray(theta_all[:, c0:c1, :]), data_stack))
        LNL[:, c0:c1] = out
    return LNL.reshape(D, npsr, Bn)


# ---------------- classification (per pulsar per draw) ----------------
def classify_full(evalL, lnL, ll_truth, L0, dL, sig, marg_pc, dL_med):
    # true wrong-fringe centres inside the +/-3 sigma window: k in [1..kmax_win] each side
    kmax_win = int(np.floor(3.0 * sig / dL + 1e-9))
    K_counted = 2 * kmax_win
    K_round = max(1, int(round(2.0 * 3.0 * sig / dL)))
    k_all = np.round((evalL - L0) / dL).astype(int)
    sel = (np.abs(k_all) >= 1) & (np.abs(k_all) <= kmax_win)      # sampled pts on in-window wrong fringes
    if kmax_win >= 1 and sel.any():
        kk, lw = k_all[sel], lnL[sel]
        order = np.argsort(kk)
        kk_s, lw_s = kk[order], lw[order]
        uk, idx0 = np.unique(kk_s, return_index=True)
        lnL_u = np.maximum.reduceat(lw_s, idx0)                    # per-fringe peak lnL
        offs_u = uk * dL
        best_wrong = float(lnL_u.max())
        dlnL = ll_truth - best_wrong
        logw = (lnL_u - ll_truth) - offs_u ** 2 / (2.0 * sig ** 2)
        P_true = 1.0 / (1.0 + float(np.exp(logw).sum()))
    else:
        K_counted = 0; dlnL = np.inf; P_true = 1.0                 # prior alone pins one fringe
    lnKc = np.log(max(K_counted, 1))
    ident = dlnL > lnKc
    ident_strict = dlnL > lnKc + 3.0
    ident_round = dlnL > np.log(K_round)
    prior_cert = K_counted == 0
    geom_wide = dL > 5.0 * dL_med          # antipodal geometry -> huge dL, low SNR
    return dict(dlnL=(dlnL if np.isfinite(dlnL) else 50.0), K_counted=K_counted, K_round=K_round,
                ident=ident, ident_strict=ident_strict, ident_round=ident_round,
                prior_cert=prior_cert, geom_wide=geom_wide, P_true=P_true,
                bayes=P_true > 0.9, bayes_strict=P_true > 0.99,
                improved=ident and (marg_pc < sig * KPC_TO_PC))


# ---------------- one config: build amo, batch draws ----------------
def build_amo(ent_psrs, disco_psrs, ncw):
    pta, cwb, _ = build_enterprise_pta(ent_psrs, ncw, components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
    inj0 = make_geometry_injection(pta, ent_psrs, ncw, cwb, seed=1000, h_range=(EQUAL_H, EQUAL_H))
    amo = build_fisher_amortised(disco_psrs, ncw, inj0, cwb)
    return pta, cwb, amo


def prep_draw(amo, ent_psrs, disco_psrs, inj, cwb, priors):
    """Per-draw: theta, data, dL(kpc), marg_pc, EV(npsr,B), AI(npsr,)."""
    temp = _enterprise_to_disco_params(inj, cwb)
    theta = np.array([temp[k] for k in amo["theta_keys"]], dtype=np.float64)
    data = amo["inject_data"](jnp.asarray(theta))
    H = assemble_hessian(amo, jnp.asarray(theta), data)
    n_dist = amo["n_dist"]
    negc = int((np.diag(H)[:n_dist] < 0).sum())
    cond, marg, _ = sigma_L_from_hessian(H, n_dist, rcond=1e-10)
    marg_pc = marg * KPC_TO_PC
    cw_list = [{k: inj[f"{name}_{k}"] for k in CW_PARAM_BASE} for name in cwb]
    npsr = len(ent_psrs)
    dL = np.array([min(compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"], cw["log10_fgw"],
                       disco_psrs[a].pos) for cw in cw_list) for a in range(npsr)])
    idx_of = {k: i for i, k in enumerate(amo["theta_keys"])}
    AI = np.array([idx_of[f"{pe.name}_cw_p_dist"] for pe in ent_psrs])
    EV = np.empty((npsr, B))
    for a, pe in enumerate(ent_psrs):
        R = 3.0 * max(priors["feather"][a], priors["script"][a], priors["lit"][a])
        EV[a] = eval_grid(pe.pdist[0], dL[a], R)
    return dict(theta=theta, data=data, dL=dL, marg_pc=np.array(marg_pc),
                cond_pc=np.array(cond) * KPC_TO_PC, EV=EV, AI=AI, negc=negc)


BOOL_KEYS = ("ident", "ident_strict", "ident_round", "prior_cert", "geom_wide",
             "bayes", "bayes_strict", "improved")


def classify_draw(pre, LNL_d, ent_psrs, priors):
    npsr = len(ent_psrs)
    dL_med = float(np.median(pre["dL"]))          # kpc
    out = {ps: {k: np.zeros(npsr, bool) for k in BOOL_KEYS} for ps in ("feather", "script", "lit")}
    for ps in out:
        out[ps]["P_true"] = np.zeros(npsr)
        out[ps]["dlnL"] = np.zeros(npsr)
        out[ps]["K_counted"] = np.zeros(npsr, int)
    lltruth = pre["ll_truth"]
    for a, pe in enumerate(ent_psrs):
        evalL, lnL = pre["EV"][a], LNL_d[a]
        for ps in ("feather", "script", "lit"):
            c = classify_full(evalL, lnL, lltruth, pe.pdist[0], pre["dL"][a],
                              priors[ps][a], pre["marg_pc"][a], dL_med)
            o = out[ps]
            for k in BOOL_KEYS:
                o[k][a] = c[k]
            o["P_true"][a] = c["P_true"]; o["dlnL"][a] = c["dlnL"]; o["K_counted"][a] = c["K_counted"]
    return out


# ================= loop scan (reference for the STEP-2 gate) =================
def scan_loop_one(amo, theta, data, AI, EV):
    """Per-pulsar loop scan -> LNL (npsr,B). Same evals as scan_all, no cross-pulsar batching."""
    npsr = EV.shape[0]
    theta_np = np.asarray(theta)
    LNL = np.empty((npsr, B))
    for a in range(npsr):
        tb = np.repeat(theta_np[None, :], B, axis=0)
        tb[:, AI[a]] = EV[a]
        LNL[a] = np.asarray(amo["fisher_logL_batched"](jnp.asarray(tb), data))
    return LNL


def stack_data(datas):
    return jax.tree_util.tree_map(lambda *a: jnp.stack(a), *datas)


# ================= drivers =================
def run(ndraws):
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    ent_psrs, disco_psrs = load_pulsars(116)
    names = [p.name for p in ent_psrs]
    priors = load_real_priors(names)
    print(f"[A2] priors med(pc): feather {np.median(priors['feather'])*KPC_TO_PC:.1f} "
          f"script {np.median(priors['script'])*KPC_TO_PC:.1f} lit-injects {list(LIT_INJECT)}", flush=True)

    results, pop, saved = {}, {}, {}
    for ncw in NCW_LIST:
        print(f"\n===== N_CW={ncw} =====", flush=True)
        t0 = time.time()
        pta, cwb, amo = build_amo(ent_psrs, disco_psrs, ncw)
        vmap2 = make_vmap2(amo)
        # warm compile
        th0 = amo["theta_truth"]; dt0 = amo["inject_data"](th0)
        _ = assemble_hessian(amo, th0, dt0)
        _ = vmap2(jnp.zeros((1, B, amo["n_theta"])), stack_data([dt0]))
        print(f"  build+compile {time.time()-t0:.1f}s", flush=True)

        # prep all draws (Hessian/sigma_L once per draw)
        seeds = list(range(2000, 2000 + ndraws))
        pres = []
        tp = time.time()
        for d, sd in enumerate(seeds):
            inj = make_geometry_injection(pta, ent_psrs, ncw, cwb, seed=sd, h_range=(EQUAL_H, EQUAL_H))
            pre = prep_draw(amo, ent_psrs, disco_psrs, inj, cwb, priors)
            pre["ll_truth"] = float(amo["fisher_logL"](jnp.asarray(pre["theta"]), pre["data"]))
            pres.append(pre)
        print(f"  prep {ndraws} draws (Hessian+geom) {time.time()-tp:.1f}s "
              f"negc={[p['negc'] for p in pres][:3]}...", flush=True)

        # ONE batched scan across all draws
        ts = time.time()
        thetas = np.stack([p["theta"] for p in pres])
        data_stack = stack_data([p["data"] for p in pres])
        EV = np.stack([p["EV"] for p in pres])
        AI = pres[0]["AI"]
        LNL = scan_all(vmap2, thetas, data_stack, EV, AI)   # (D,npsr,B)
        print(f"  batched scan {time.time()-ts:.1f}s", flush=True)

        draws_cls = [classify_draw(pres[d], LNL[d], ent_psrs, priors) for d in range(ndraws)]
        results[ncw] = {ps: agg(draws_cls, ps) for ps in ("feather", "script", "lit")}
        results[ncw]["_named"] = {ps: named(draws_cls, ps, names) for ps in ("script", "lit")}
        # per-draw per-pulsar diagnostics (for cause-tagging, onset fractions, valid counts)
        dg = dict(names=np.array(names),
                  dL_pc=np.array([p["dL"] * KPC_TO_PC for p in pres]),         # (D,116)
                  marg_pc=np.array([p["marg_pc"] for p in pres]),
                  sig_feather_pc=priors["feather"] * KPC_TO_PC,
                  sig_script_pc=priors["script"] * KPC_TO_PC,
                  sig_lit_pc=priors["lit"] * KPC_TO_PC)
        for ps in ("feather", "script", "lit"):
            for key in ("ident", "ident_strict", "bayes", "bayes_strict", "prior_cert", "geom_wide"):
                dg[f"{ps}_{key}"] = np.array([d[ps][key] for d in draws_cls])   # (D,116) bool
            dg[f"{ps}_P_true"] = np.array([d[ps]["P_true"] for d in draws_cls])
            dg[f"{ps}_dlnL"] = np.array([d[ps]["dlnL"] for d in draws_cls])
            dg[f"{ps}_K_counted"] = np.array([d[ps]["K_counted"] for d in draws_cls])
        saved[ncw] = dg

        if ncw == 16:
            print("  -- population draw (3 loud+13 faint) --", flush=True)
            injp = make_geometry_injection(pta, ent_psrs, 16, cwb, seed=3000, population=(3, -13.25, -14.25))
            prep = prep_draw(amo, ent_psrs, disco_psrs, injp, cwb, priors)
            prep["ll_truth"] = float(amo["fisher_logL"](jnp.asarray(prep["theta"]), prep["data"]))
            LNLp = scan_all(vmap2, prep["theta"][None], stack_data([prep["data"]]),
                            prep["EV"][None], prep["AI"])[0]
            pcls = classify_draw(prep, LNLp, ent_psrs, priors)
            pop = {ps: {k: int(pcls[ps][k].sum()) for k in ("ident", "ident_strict", "prior_cert", "bayes", "bayes_strict", "improved")}
                   for ps in ("feather", "script", "lit")}
            pop["_named_lit"] = [names[a] for a in range(len(names)) if pcls["lit"]["ident"][a] or pcls["lit"]["bayes"][a]]
            # population per-pulsar diagnostics
            popdiag = dict(names=np.array(names), dL_pc=prep["dL"] * KPC_TO_PC, marg_pc=prep["marg_pc"],
                           sig_lit_pc=priors["lit"] * KPC_TO_PC, sig_script_pc=priors["script"] * KPC_TO_PC,
                           sig_feather_pc=priors["feather"] * KPC_TO_PC)
            for ps in ("feather", "script", "lit"):
                for key in ("ident", "bayes", "bayes_strict", "prior_cert", "geom_wide"):
                    popdiag[f"{ps}_{key}"] = pcls[ps][key]
                popdiag[f"{ps}_P_true"] = pcls[ps]["P_true"]
                popdiag[f"{ps}_dlnL"] = pcls[ps]["dlnL"]
                popdiag[f"{ps}_K_counted"] = pcls[ps]["K_counted"]
            pop["_diag"] = popdiag

    report(results, pop, NCW_LIST)
    np.savez(CWT + "/anchor_a2_results.npz", ncw_list=np.array(NCW_LIST),
             results=np.array([results], dtype=object), pop=np.array([pop], dtype=object))
    np.savez(CWT + "/anchor_a2_diag.npz", diag=np.array([saved], dtype=object),
             popdiag=np.array([pop.get("_diag")], dtype=object))
    print("\n[A2] saved anchor_a2_results.npz + anchor_a2_diag.npz")
    return 0


def agg(draws_cls, ps):
    def stat(key):
        v = np.array([d[ps][key].sum() for d in draws_cls])
        return (int(np.median(v)), int(v.min()), int(v.max()))
    return {k: stat(k) for k in ("ident", "ident_strict", "ident_round", "prior_cert", "bayes", "bayes_strict", "improved")}


def named(draws_cls, ps, names):
    cnt = {}
    for d in draws_cls:
        for a in range(len(names)):
            if d[ps]["ident"][a] or d[ps]["bayes"][a]:
                cnt[names[a]] = cnt.get(names[a], 0) + 1
    return dict(sorted(cnt.items(), key=lambda x: -x[1]))


def report(results, pop, ncw_list):
    print("\n================ A2 CERTIFICATION COUNTS  med(min,max) ================")
    for metric, lab in [("ident", "ANCHOR-CERTIFIED flat (dlnL>lnK_counted)"),
                        ("ident_strict", "flat strict (+3)"),
                        ("bayes", "CERTIFIED-BAYES (P_true>0.9)"),
                        ("bayes_strict", "bayes strict (>0.99)"),
                        ("prior_cert", "PRIOR-CERTIFIED (0 wrong fringes)"),
                        ("improved", "IMPROVED (ident & marg<sigma)")]:
        print(f"\n-- {lab} --")
        print(f"{'N_CW':>4s} | {'feather':>14s} | {'script-real':>14s} | {'literature':>14s}")
        for n in ncw_list:
            r = results[n]
            print(f"{n:>4d} | {str(r['feather'][metric]):>14s} | {str(r['script'][metric]):>14s} | {str(r['lit'][metric]):>14s}")
    print("\n---- named certified/bayes (script | lit), hits/ndraws ----")
    for n in ncw_list:
        print(f"  N_CW={n}: script={results[n]['_named']['script']}  lit={results[n]['_named']['lit']}")
    print("\n---- population N_CW=16 (3 loud+13 faint) ----")
    for ps in ("feather", "script", "lit"):
        print(f"  {ps:>8s}: {pop.get(ps)}")
    print(f"  named (lit): {pop.get('_named_lit')}")


def gate():
    """STEP-2 gate: vectorized scan == loop scan on one draw to 1e-12."""
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    ent_psrs, disco_psrs = load_pulsars(116)
    names = [p.name for p in ent_psrs]
    priors = load_real_priors(names)
    pta, cwb, amo = build_amo(ent_psrs, disco_psrs, 1)
    vmap2 = make_vmap2(amo)
    inj = make_geometry_injection(pta, ent_psrs, 1, cwb, seed=2000, h_range=(EQUAL_H, EQUAL_H))
    pre = prep_draw(amo, ent_psrs, disco_psrs, inj, cwb, priors)
    LNL_loop = scan_loop_one(amo, jnp.asarray(pre["theta"]), pre["data"], pre["AI"], pre["EV"])
    LNL_vec = scan_all(vmap2, pre["theta"][None], stack_data([pre["data"]]), pre["EV"][None], pre["AI"])[0]
    err = np.abs(LNL_loop - LNL_vec).max()
    print(f"[GATE] max|loop - vectorized| over 116x{B} evals = {err:.3e}  "
          f"({'PASS' if err < 1e-12 else 'FAIL'})", flush=True)
    return 0 if err < 1e-12 else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "run"], default="run", nargs="?")
    ap.add_argument("--ndraws", type=int, default=10)
    a = ap.parse_args()
    sys.exit(gate() if a.mode == "gate" else run(a.ndraws))
