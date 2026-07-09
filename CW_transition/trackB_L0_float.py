"""Track B — LAMBDA feasibility probe, DELIVERABLE L0: the FLOAT LANDING ZONE.

Truth-blind soft-EM cold start (F-stat seed -> soft-EM alternation, 30-iter budget,
--use-split ON, stall rule) exactly as trackB_p2_pipeline.mode_gate, but SAVES THE FULL
FLOAT STATE that L1/L2 need:
  - source params (float) + truth, offset diagnostics
  - per-pulsar FULL fringe posterior q_p(n): (uk integers, offsets, lnL, weights), padded
  - truth integer n_true_p = round((L_true_p - L0_p)/dL_p)  (== 0 iff true dist = prior mean)
  - loud-seed provenance (are the 3 loud sources actually IN the top-N_CW seeds used?)
  - census-3 posteriors, carriers, integer support counts, per-iter traces

Saved to trackB_L0_float.npz. Finiteness checks printed; NaN fraction>0 aborts.
Run: python trackB_L0_float.py --use-split   (jug-gpu, background w/ log)
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

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
from trackB_p1_map import ang_deg
import trackB_p2_pipeline as PP
import itertools


def loud_registration_offsets(src, src_t, scale, n_loud):
    """Best 3x3 sky assignment; per assigned loud src return (truth_idx, sky_deg, scaled_sky, scaled_freq)."""
    def U(s, k):
        cg, gp = s[8*k], s[8*k+1]; ss = np.sqrt(max(1-cg**2, 0))
        return np.array([ss*np.cos(gp), ss*np.sin(gp), cg])
    M = np.array([[np.degrees(np.arccos(np.clip(np.dot(U(src,i), U(src_t,j)), -1, 1)))
                   for j in range(n_loud)] for i in range(n_loud)])
    best = min(itertools.permutations(range(n_loud)), key=lambda p: sum(M[i, p[i]] for i in range(n_loud)))
    rows = []
    for i in range(n_loud):
        j = best[i]
        sc_sky = max(abs(src[8*i+0]-src_t[8*j+0])/scale[0], abs(src[8*i+1]-src_t[8*j+1])/scale[1])
        sc_frq = abs(src[8*i+4]-src_t[8*j+4])/scale[4]
        rows.append((j, M[i, j], sc_sky, sc_frq))
    return best, rows

CWT = "/home/mattm/projects/HSYMT/CW_transition"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]


COMPACT_C = 72   # fixed compact config count (<=24 gated anchors x KTOP=3); w=0 rows dropped


def compact_configs(dcfg, w, C=COMPACT_C):
    """EXACT speedup: gradQ = sum_c w_c*grad(lnL_c); w=0 configs contribute nothing.
    Keep the top-C nonzero-weight rows (padded to fixed C with w=0 rows for a stable
    compiled shape). Identical gradient to the full 348-row scan when nonzero rows <= C."""
    w = np.asarray(w); nz = int(np.sum(w > 0))
    order = np.argsort(-w)[:C]
    d_out = dcfg[order].copy(); w_out = w[order].copy()
    if len(order) < C:                                   # pad with zero-weight rows
        pad = C - len(order)
        d_out = np.concatenate([d_out, np.repeat(dcfg[:1], pad, axis=0)], axis=0)
        w_out = np.concatenate([w_out, np.zeros(pad)])
    if nz > C:
        print(f"    !! WARN {nz} nonzero configs > compact cap {C}; dropped smallest "
              f"{nz-C} (approx)", flush=True)
    return d_out, w_out, nz


def finite_report(name, arr):
    a = np.asarray(arr, float).ravel()
    fin = np.isfinite(a)
    frac = float(fin.mean()) if a.size else 1.0
    nanfrac = float(np.isnan(a).mean()) if a.size else 0.0
    if a[fin].size:
        lo, md, hi = float(a[fin].min()), float(np.median(a[fin])), float(a[fin].max())
    else:
        lo = md = hi = float("nan")
    print(f"    [finite] {name:16s} isfinite={frac:.4f} nan={nanfrac:.4f} "
          f"min/med/max={lo:.4g}/{md:.4g}/{hi:.4g}", flush=True)
    if nanfrac > 0.0:
        print(f"    !! NaN fraction > 0 in {name} -> ABORT", flush=True)
        sys.exit(2)


def loud_seed_provenance(P, seeds):
    """Are the 3 loud truths represented among the top-ncw seeds actually used
    (seeds_to_phi keeps seeds[:ncw])? Report per-loud nearest USED-seed angle + rank."""
    ncw = P["ncw"]
    used = seeds[:ncw]
    louds = TE.loud_truth(P)
    rows = []
    for li, lt in enumerate(louds):
        best_ang, best_rank = 1e9, -1
        for r, s in enumerate(used):
            ang = np.degrees(np.arccos(np.clip(np.dot(lt["u"], s["u"]), -1, 1)))
            if ang < best_ang:
                best_ang, best_rank = ang, r
        rows.append((li, best_ang, best_rank, float(used[best_rank]["twoF"]) if best_rank >= 0 else np.nan))
    return rows, len(seeds)


def pack_float_state(P, post, src, base_dist, n_true, census_idx, status, hist, ev_certified=None):
    """Full float state -> dict for np.savez (used both for per-iter checkpoints and final)."""
    npsr = P["npsr"]; ndist = P["n_dist"]
    Kmax = max(len(post["qlist"][p][0]) for p in range(npsr))
    q_uk = np.zeros((npsr, Kmax), int); q_off = np.full((npsr, Kmax), np.nan)
    q_lnL = np.full((npsr, Kmax), np.nan); q_w = np.zeros((npsr, Kmax)); q_K = np.zeros(npsr, int)
    for p in range(npsr):
        uk, offs_u, lnL_u, wp = post["qlist"][p]
        k = len(uk); q_K[p] = k
        q_uk[p, :k] = uk; q_off[p, :k] = offs_u; q_lnL[p, :k] = lnL_u; q_w[p, :k] = wp
    qmax = post["map_qmax"]; map_k = post["map_k"]
    supp = np.array([int(np.sum(q_w[p, :q_K[p]] > 0.01)) for p in range(npsr)])
    carriers = np.where(qmax > 0.9)[0]
    w_on_true = np.zeros(npsr)
    for p in range(npsr):
        m = np.where(q_uk[p, :q_K[p]] == n_true[p])[0]
        w_on_true[p] = float(q_w[p, m[0]]) if m.size else 0.0
    scale = TE.phi_scale(P); src_t = P["theta_truth"][ndist:]
    n_loud = TE.CONFIGS["pop"]["population"][0]; loud_sl = slice(0, 8 * n_loud)
    scaled_off = np.abs(np.asarray(src)[loud_sl] - src_t[loud_sl]) / scale[loud_sl]
    d = dict(src=src, src_truth=src_t, base_dist=base_dist, L0=np.asarray(P["L0"]),
             dL=np.asarray(P["dL"]), L_true=np.asarray(P["theta_truth"][:ndist]), n_true=n_true,
             names=np.array(P["names"]), census_idx=np.array(census_idx),
             q_uk=q_uk, q_off=q_off, q_lnL=q_lnL, q_w=q_w, q_K=q_K,
             qmax=qmax, map_k=map_k, int_support=supp, carriers=carriers,
             w_on_true=w_on_true, loud_scaled_off=scaled_off, status=status,
             hist_Q=np.array([h["Q"] for h in hist]),
             hist_qtsum=np.array([h["qt_sum"] for h in hist]),
             hist_anchors=np.array([h["anchors"] for h in hist]),
             hist_off=np.array([h["off_max"] for h in hist]),
             hist_census=np.array([h["cen"] for h in hist]))
    if ev_certified is not None:
        d["certified"] = np.array(ev_certified)
    return d


def main(max_iter=30, use_split=True, wall_cap=9000.0):
    print(f"=== L0 FLOAT LANDING ZONE (truth-blind soft-EM, use_split={use_split}, "
          f"wall_cap={wall_cap:.0f}s) ===", flush=True)
    t_all = time.time()
    P = TE.build_problem("pop", build_earth=False); prior = P["priors"]["lit"]
    names = list(P["names"]); ndist = P["n_dist"]; npsr = P["npsr"]
    census_idx = [names.index(n) for n in CENSUS]
    mstep = PP.make_soft_mstep(P)
    split = None
    if use_split:
        from trackB_split import Split
        split = Split(P["amo_full"], P["data_obs"], P["theta_truth"])

    # truth integers: n_true_p = round((L_true_p - L0_p)/dL_p)
    L0 = np.asarray(P["L0"]); dL = np.asarray(P["dL"])
    L_true = np.asarray(P["theta_truth"][:ndist])
    n_true = np.round((L_true - L0) / dL).astype(int)
    print(f"  truth integers: n_true all-zero for {int(np.sum(n_true==0))}/{npsr} pulsars "
          f"(nonzero: {np.where(n_true!=0)[0].tolist()[:10]})", flush=True)

    # ---- 1. F-stat seed (truth-blind) ----
    t0 = time.time()
    scan = TE.seed_scan(P)
    seeds = TE.pick_seeds(scan, twoF_thresh=25.0)
    nrec = TE.report_seed_recovery(P, seeds, 25.0)
    prov, nseed = loud_seed_provenance(P, seeds)
    print(f"  loud-seed provenance (in top-{P['ncw']} USED seeds):", flush=True)
    for li, ang, rank, twoF in prov:
        flag = "OK" if rank >= 0 and ang < 20 else "!! DROPPED/FAR"
        print(f"     loud{li}: nearest used seed rank={rank} ang={ang:.2f}deg 2F={twoF:.1f}  {flag}", flush=True)
    src = TE.seeds_to_phi(P, seeds)[:]
    print(f"  F-stat seed: {nseed} cand, ncw-used={P['ncw']}, loud recovered {nrec}/3, {time.time()-t0:.0f}s", flush=True)

    # ---- 2. soft-EM alternation ----
    base_dist = L0.copy()
    LR = lambda it: 1.5e-2 if it < 8 else (6e-3 if it < 18 else 2e-3)
    hist = []
    q_prev = -np.inf; qt_prev = -np.inf; plateau = 0; status = "budget"
    for it in range(max_iter):
        theta = np.concatenate([base_dist, src])
        post, LNL = PP.soft_estep(P, theta, prior, split=split)
        base_dist = post["map_evalL"].copy()
        Pt = PP.qtrue_trace(P, LNL, prior)
        qt_sum = float(np.sum(Pt)); cen = [float(Pt[i]) for i in census_idx]
        off_max, offs = PP.source_offset_from_truth(P, src)
        dcfg, w, anchors = PP.build_mstep_configs(P, post, base_dist)
        dcfg, w, nz_cfg = compact_configs(dcfg, w)       # drops w=0 rows (~6x faster); NOT bit-exact
        src, qval = mstep(src, dcfg, w, 25, LR(it))       # 25 Adam steps (M-step numerically chaotic)
        ncert = int(np.sum(Pt > 0.9))
        cert_census = all(Pt[census_idx[j]] > 0.9 for j in range(3))
        # label-matched REGISTRATION-critical offsets (sky/freq per loud src) -- the wall-height read
        _scale = TE.phi_scale(P); _st = P["theta_truth"][ndist:]
        _nl = TE.CONFIGS["pop"]["population"][0]
        _perm, _rows = loud_registration_offsets(src, _st, _scale, _nl)
        reg = " ".join(f"f{i}->t{j}:sky{sd:.1f}deg/sc{ss:.3f},fr{sf:.3f}" for i,(j,sd,ss,sf) in enumerate(_rows))
        hist.append(dict(it=it, Q=qval, qt_sum=qt_sum, anchors=anchors, off_max=off_max,
                         ncert=ncert, cen=cen))
        print(f"  it{it:2d} Q={qval:.1f} sumQtrue={qt_sum:.2f} anchors={anchors} "
              f"ncert={ncert} census[{cen[0]:.3f},{cen[1]:.3f},{cen[2]:.3f}] lr={LR(it):.0e} "
              f"({time.time()-t_all:.0f}s)\n       REG {reg}", flush=True)
        dq = abs(qval - q_prev) / (abs(q_prev) + 1e-9); dqt = abs(qt_sum - qt_prev) / (abs(qt_prev) + 1e-9)
        plateau = plateau + 1 if (dq < 1e-3 and dqt < 1e-3) else 0
        q_prev = qval; qt_prev = qt_sum
        # per-iter checkpoint (post/LNL of THIS iter's E-step -> recoverable float state)
        ckpt = pack_float_state(P, post, src, base_dist, n_true, census_idx, "checkpoint", hist)
        ckpt["iter"] = it; ckpt["nrec"] = nrec
        np.savez(f"{CWT}/trackB_L0_float.npz", **ckpt)
        if cert_census:
            status = "certified"; print(f"  -> census-3 certified at it{it}", flush=True); break
        if plateau >= 3:
            status = "STALL"; print(f"  -> STALL (Q & sumQtrue plateau 3 iters)", flush=True); break
        if time.time() - t_all > wall_cap:
            status = "WALLCAP"; print(f"  -> WALLCAP hit ({wall_cap:.0f}s); checkpoint holds", flush=True); break

    # ---- 3. final E-step + full float state ----
    theta = np.concatenate([base_dist, src])
    post, LNL = PP.soft_estep(P, theta, prior, split=split)
    ev = TE.eval_gate(P, LNL, prior); Pt_final = ev["P_true"]

    # pack ragged qlist -> padded arrays
    Kmax = max(len(post["qlist"][p][0]) for p in range(npsr))
    q_uk = np.zeros((npsr, Kmax), int)
    q_off = np.full((npsr, Kmax), np.nan)
    q_lnL = np.full((npsr, Kmax), np.nan)
    q_w = np.zeros((npsr, Kmax))
    q_K = np.zeros(npsr, int)
    for p in range(npsr):
        uk, offs_u, lnL_u, wp = post["qlist"][p]
        k = len(uk); q_K[p] = k
        q_uk[p, :k] = uk; q_off[p, :k] = offs_u; q_lnL[p, :k] = lnL_u; q_w[p, :k] = wp

    qmax = post["map_qmax"]; map_k = post["map_k"]
    supp = np.array([int(np.sum(q_w[p, :q_K[p]] > 0.01)) for p in range(npsr)])
    carriers = np.where(qmax > 0.9)[0]
    # posterior mass the float puts on TRUTH's integer n_true_p, per pulsar
    w_on_true = np.zeros(npsr)
    for p in range(npsr):
        m = np.where(q_uk[p, :q_K[p]] == n_true[p])[0]
        w_on_true[p] = float(q_w[p, m[0]]) if m.size else 0.0

    print(f"\n  === FLOAT STATE SUMMARY ===", flush=True)
    off_max, offs = PP.source_offset_from_truth(P, src)
    print(f"    status={status}  final srcOff(loud sky)={off_max:.4f} deg  per-loud={[f'{o:.2f}' for o in offs]}", flush=True)
    scale = TE.phi_scale(P)
    src_t = P["theta_truth"][ndist:]
    n_loud = TE.CONFIGS["pop"]["population"][0]
    loud_sl = slice(0, 8 * n_loud)
    scaled_off = np.abs(src[loud_sl] - src_t[loud_sl]) / scale[loud_sl]
    print(f"    loud src scaled offset: max={scaled_off.max():.4f} median={np.median(scaled_off):.4f} "
          f"(cert tol ~1e-4)", flush=True)
    print(f"    integer support (#fringes q>1%): median {int(np.median(supp))} min {supp.min()} "
          f"max {supp.max()}; locked(==1) {int(np.sum(supp==1))} ==2 {int(np.sum(supp==2))}", flush=True)
    print(f"    concentrated carriers (qmax>0.9): {len(carriers)}", flush=True)
    for j, ci in enumerate(census_idx):
        print(f"    {CENSUS[j]:14s} qmax={qmax[ci]:.3f} map_k={map_k[ci]} n_true={n_true[ci]} "
              f"support={supp[ci]} w_on_true={w_on_true[ci]:.3e} P_true(final)={Pt_final[ci]:.3f}", flush=True)

    print(f"\n  === FINITENESS ===", flush=True)
    finite_report("src", src); finite_report("base_dist", base_dist)
    finite_report("q_w", q_w); finite_report("q_lnL", q_lnL[np.isfinite(q_lnL)])
    finite_report("qmax", qmax); finite_report("P_true_final", Pt_final)

    np.savez(f"{CWT}/trackB_L0_float.npz",
             src=src, src_truth=src_t, base_dist=base_dist, L0=L0, dL=dL,
             L_true=L_true, n_true=n_true,
             names=np.array(names), census_idx=np.array(census_idx),
             q_uk=q_uk, q_off=q_off, q_lnL=q_lnL, q_w=q_w, q_K=q_K,
             qmax=qmax, map_k=map_k, int_support=supp, carriers=carriers,
             w_on_true=w_on_true, P_true_final=Pt_final,
             certified=np.array(ev["certified"]),
             hist_Q=np.array([h["Q"] for h in hist]),
             hist_qtsum=np.array([h["qt_sum"] for h in hist]),
             hist_anchors=np.array([h["anchors"] for h in hist]),
             hist_off=np.array([h["off_max"] for h in hist]),
             hist_census=np.array([h["cen"] for h in hist]),
             loud_scaled_off=scaled_off, status=status, nrec=nrec)
    print(f"\n  saved trackB_L0_float.npz  (total {time.time()-t_all:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("--max_iter", type=int, default=30)
    ap.add_argument("--use-split", action="store_true")
    ap.add_argument("--wall_cap", type=float, default=9000.0)
    a = ap.parse_args()
    sys.exit(main(a.max_iter, use_split=a.use_split, wall_cap=a.wall_cap))
