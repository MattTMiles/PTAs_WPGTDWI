"""Track B — P2: the truth-blind soft-EM cold-start pipeline + B0.2 re-gate.

Probe 2 verdict: the flat P1 plateau is a gauge conspiracy; the SOFT (marginal/posterior)
surface is navigable across the gap with census-class anchors that EMERGE from the E-step
(loud-source-broken, prior-independent). This pipeline realises the soft-EM descent and
runs the make-or-break B0.2 gate on it, end to end, TRUTH-BLIND.

Pipeline (no truth info in seeds/thresholds/anchor-emergence/zoom-centres):
  1. F-STAT SEED (Earth-term matched filter, amplitude/phase/psi/inc profiled) places the
     loud sources ~6-12 deg / <1 freq bin (B0.2 gate 3a). trackB_estimator.seed_scan.
  2. SOFT-EM ALTERNATION (the funnel descent 12 deg->2 deg AND the gap crossing 2 deg->cusp
     are both driven by the soft M-step Q-gradient, which contains the Earth-term funnel at
     large offset and the anchored soft-fringe gradient across the gap; staged learning
     rate = the "coarse->fine" zoom):
       E-step: soft per-pulsar fringe posteriors q_p(n) (A2 softmax + Gaussian EM prior).
         Anchors are NOT assigned -- the pulsars whose q concentrates (max_n q>0.9) are the
         emergent anchor set, logged.
       M-step: ascend the EM Q-function  Q(src) = sum_p sum_n q_p(n) lnL(src, L_p(n))
         over the 8*N_CW source params (q fixed; the standard EM soft M-step). This is the
         SOFT surface probe 2 showed is navigable, NOT the flat HARD-snap surface the B0.2
         pipeline used.
  3. CERTIFY (E-step P_true, A2 criterion) -> census-3 within 0.05 of ceilings, no false
     certs.

INSTRUMENTATION (per EM iter, the ignition trace -- goes in the paper either way):
  QTRUE per pulsar (weight on the true fringe; DIAGNOSTIC, truth only for the log, never
  fed back), emergent anchor set (max_n q>0.9), source offset from truth (DIAGNOSTIC), the
  M-step Q objective, and P_true for the 3 ceiling pulsars.

STALL RULE (pre-registered, honest): if the Q objective AND sum-QTRUE both plateau
(rel change < 1e-3) for 3 consecutive iters and census-3 not yet certified, declare STALL,
checkpoint, stop -- a FAIL datum, not a tuning invitation. Budget 30 iters.

GATE: same 3 pulsars certified (J0711/J1713/J1909), P_true within ~0.05 of A2 ceilings
(0.953/0.989/0.998), zero false certifications, truth-blind. PASS -> build the low-rank
split, then B1. FAIL/STALL -> report WHICH leg failed (descent / ignition / lock).

Run (jug-gpu, background):
  python trackB_p2_pipeline.py test-mstep   # from-truth + perturbation M-step gate
  python trackB_p2_pipeline.py gate         # full truth-blind pipeline + B0.2 gate
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
from trackB_p1_map import batched_scan, ang_deg

CWT = "/home/mattm/projects/HSYMT/CW_transition"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]
CEILING = TE.CEILING
KTOP = 3                              # top-K fringes per pulsar in the soft M-step mixture


# ============================================================
# soft E-step: per-pulsar posteriors, emergent anchors, QTRUE (diagnostic)
# ============================================================
def soft_estep(P, theta, prior, beta=1.0, split=None):
    # split (a trackB_split.Split, from use_split) gives the ~47x faster E-step, bit-exact
    # vs batched_scan (GATE 1). Default None -> batch-only (the correctness path).
    if split is not None:
        LNL = split.estep(np.asarray(theta), P["EV"], P["AI"])   # (npsr,B)
    else:
        LNL, _ = batched_scan(P, [theta]); LNL = LNL[0]          # (npsr,B)
    post = TE.fringe_posterior(P, LNL, None, prior, beta)
    return post, LNL


def qtrue_trace(P, LNL, prior):
    """DIAGNOSTIC only: per-pulsar posterior weight on the TRUE fringe (k=0, since the true
    distance = L0). Uses truth ONLY for the log; never fed to the algorithm."""
    Pt, _, _ = TE.p_true_at(P, LNL, prior, beta=1.0)
    return Pt


# ============================================================
# soft M-step: ascend Q(src) = sum_{p,n} q_p(n) lnL(src, L_p(n)) over source params
# ============================================================
ANCHOR_TAU = 0.9      # M-step uses only pulsars with q_max > tau (emergent anchors)
ANCHOR_MIN = 3        # if fewer than this qualify, lower tau to admit the top-ANCHOR_MIN


def build_mstep_configs(P, post, base_dist):
    """Fixed-shape (npsr*KTOP) mixture, HARD-GATED to the emergent anchors (RTK 'use only
    resolved satellites'): only pulsars with q_max > ANCHOR_TAU drive the source update.
    The fringe-degenerate majority (spread posteriors) is EXCLUDED -- their summed
    wrong-fringe gradient otherwise outvotes the anchors and drags sky/freq off truth
    (making truth a non-fixed-point). At truth the gated pulsars sit on their true fringe
    (grad lnL = 0) -> truth stationary. Returns dist_configs (C,ndist), weights (C,), n_anchor."""
    ndist = P["n_dist"]; npsr = P["npsr"]; AI = P["AI"]
    qmax = post["map_qmax"]
    tau = ANCHOR_TAU
    if int(np.sum(qmax > tau)) < ANCHOR_MIN:                 # bootstrap: admit the top few
        tau = float(np.sort(qmax)[-ANCHOR_MIN]) - 1e-9 if npsr >= ANCHOR_MIN else 0.0
    gate = qmax > tau
    C = npsr * KTOP
    dist_cfg = np.repeat(base_dist[None, :ndist], C, axis=0)
    w = np.zeros(C)
    for p in range(npsr):
        if not gate[p]:
            continue
        uk, offs_u, lnL_u, wp = post["qlist"][p]
        order = np.argsort(-wp)[:KTOP]
        conf = float(wp.max())
        L0 = P["L0"][p]; dLa = P["dL"][p]
        for j, o in enumerate(order):
            c = p * KTOP + j
            dist_cfg[c, AI[p]] = L0 + uk[o] * dLa            # this pulsar at fringe centre
            w[c] = wp[o] * conf
    return dist_cfg, w, int(gate.sum())


def make_soft_mstep(P):
    """Jitted Adam ascent of the EM Q-function on the source params (distances fixed per
    mixture config). The soft surface probe 2 showed is navigable."""
    fl = P["amo_full"]["fisher_logL"]; data = P["data_obs"]; ndist = P["n_dist"]
    step_scale = jnp.asarray(TE.phi_scale(P))
    plo, phi_hi = TE.phi_bounds(P); plo = jnp.asarray(plo); phi_hi = jnp.asarray(phi_hi)

    # grad of a SUM = sum of grads: scan over mixture configs, one array-lnL BACKWARD at a
    # time (O(1) memory in the config count -- the vmap-all form OOMs at 32 GiB).
    def _cfg_grad(src, dist_c, w_c):
        g = jax.grad(lambda s: fl(jnp.concatenate([dist_c, s]), data))(src)
        # per-config nan_to_num BEFORE the weighted sum: a single degenerate-fringe config
        # (incl. w=0 ones where 0*NaN=NaN) must NOT poison the whole gradQ sum (F1a bug).
        return w_c * jnp.nan_to_num(g)

    def gradQ(src, dist_cfg, weights):
        def body(acc, xs):
            d_c, w_c = xs
            return acc + _cfg_grad(src, d_c, w_c), None
        g, _ = jax.lax.scan(body, jnp.zeros_like(src), (dist_cfg, weights))
        return g

    def Qval(src, dist_cfg, weights):
        def body(acc, xs):
            d_c, w_c = xs
            return acc + w_c * fl(jnp.concatenate([d_c, src]), data), None
        q, _ = jax.lax.scan(body, 0.0, (dist_cfg, weights))
        return q

    b1, b2, eps = 0.9, 0.999, 1e-8

    @jax.jit
    def adam(src0, dist_cfg, weights, nsteps, lr):
        def body(t, c):
            src, m, v = c
            g = jnp.nan_to_num(gradQ(src, dist_cfg, weights))
            m = b1 * m + (1 - b1) * g; v = b2 * v + (1 - b2) * g * g
            mh = m / (1 - b1 ** t); vh = v / (1 - b2 ** t)
            step = lr * step_scale * mh / (jnp.sqrt(vh) + eps)
            step = jnp.clip(step, -0.5 * step_scale, 0.5 * step_scale)
            src = jnp.clip(src + step, plo, phi_hi)
            return (src, m, v)
        z = jnp.zeros_like(src0)
        src, _, _ = jax.lax.fori_loop(1, nsteps + 1, body, (jnp.clip(src0, plo, phi_hi), z, z))
        return src, Qval(src, dist_cfg, weights)

    def mstep(src, dist_cfg, weights, nsteps, lr):
        s, q = adam(jnp.asarray(src), jnp.asarray(dist_cfg), jnp.asarray(weights),
                    int(nsteps), float(lr))
        return np.asarray(s), float(q)
    return mstep


# ============================================================
# from-truth + perturbation M-step gate (de-risk before the full run)
# ============================================================
def _drift_groups(src, src_t, scale, n_loud):
    """max scaled drift split by registration-critical vs nuisance params, per loud src.
    sky={0,1}, freq={4} are registration-critical; nuis={2,3,5,6,7} (inc,mc,h,phase,psi)."""
    def mx(idxset):
        vals = []
        for k in range(n_loud):
            for i in idxset:
                j = 8 * k + i
                vals.append(abs(src[j] - src_t[j]) / scale[j])
        return max(vals)
    return mx([0, 1]), mx([4]), mx([2, 3, 5, 6, 7])


def mode_test_mstep():
    print("=== test-mstep: soft M-step -- from truth (fixed pt) + perturb; cert-from-truth ===", flush=True)
    P = TE.build_problem("pop"); prior = P["priors"]["lit"]; names = P["names"]
    mstep = make_soft_mstep(P)
    ndist = P["n_dist"]; theta_t = P["theta_truth"].copy()
    src_t = theta_t[ndist:].copy(); scale = TE.phi_scale(P)
    n_loud = TE.CONFIGS["pop"].get("population", (P["ncw"],))[0]
    loud_sl = slice(0, 8 * n_loud)
    census_idx = [names.index(n) for n in CENSUS]

    def certify(src):
        post, LNL = soft_estep(P, np.concatenate([theta_t[:ndist], src]), prior)
        Pt, _, _ = TE.p_true_at(P, LNL, prior, beta=1.0)
        return [float(Pt[i]) for i in census_idx], int(np.sum(Pt > 0.9))

    cen_t, nc_t = certify(src_t)
    print(f"  census-3 P_true AT TRUTH (baseline): {cen_t} ncert={nc_t}", flush=True)

    # from truth: ascend Q; sky/freq should HOLD (nuisance may wander harmlessly)
    post, _ = soft_estep(P, theta_t, prior)
    dcfg, w, nanch = build_mstep_configs(P, post, theta_t[:ndist])
    src1, q1 = mstep(src_t, dcfg, w, 40, 8e-3)
    sky, frq, nui = _drift_groups(src1, src_t, scale, n_loud)
    cen1, nc1 = certify(src1)
    print(f"  from TRUTH after 40 Adam (n_anchor={nanch}): drift sky={sky:.4f} freq={frq:.4f} "
          f"nuis={nui:.4f} | census-3 {cen1} ncert={nc1}", flush=True)

    rng = np.random.default_rng(7)
    for delta in [0.05, 0.2, 0.5]:
        src_p = src_t.copy(); src_p[loud_sl] += rng.normal(size=8 * n_loud) * scale[loud_sl] * delta
        post_p, _ = soft_estep(P, np.concatenate([theta_t[:ndist], src_p]), prior)
        dcfg_p, w_p, _ = build_mstep_configs(P, post_p, theta_t[:ndist])
        src_out, q_out = mstep(src_p, dcfg_p, w_p, 60, 8e-3)
        s0, f0, _ = _drift_groups(src_p, src_t, scale, n_loud)
        s1, f1, _ = _drift_groups(src_out, src_t, scale, n_loud)
        cen_o, nc_o = certify(src_out)
        print(f"  delta={delta}: sky {s0:.4f}->{s1:.4f}, freq {f0:.4f}->{f1:.4f} "
              f"({'TOWARD' if (s1+f1)<(s0+f0) else 'AWAY'}) | census-3 {cen_o} ncert={nc_o}", flush=True)
    return 0


# ============================================================
# the truth-blind pipeline + B0.2 gate
# ============================================================
def source_offset_from_truth(P, src):
    """DIAGNOSTIC (truth used only for the log): max over loud sources of sky-angle offset
    (deg) between src and truth loud positions."""
    ndist = P["n_dist"]; tt = P["theta_truth"][ndist:]
    n_loud = TE.CONFIGS["pop"].get("population", (P["ncw"],))[0]
    offs = []
    for k in range(n_loud):
        cg, gp = src[8 * k + 0], src[8 * k + 1]
        cgt, gpt = tt[8 * k + 0], tt[8 * k + 1]
        offs.append(ang_deg(cg, gp, cgt, gpt))
    return float(np.max(offs)), offs


def mode_gate(max_iter=30, save=True, use_split=False):
    print(f"=== P2 GATE: truth-blind soft-EM cold-start + B0.2 re-gate (config=pop"
          f"{', use_split' if use_split else ''}) ===", flush=True)
    P = TE.build_problem("pop", build_earth=False); prior = P["priors"]["lit"]
    names = P["names"]; ndist = P["n_dist"]
    census_idx = [names.index(n) for n in CENSUS]
    mstep = make_soft_mstep(P)
    # E-step accelerator (~47x, bit-exact vs batched_scan; GATE 1). Default OFF.
    split = None
    if use_split:
        from trackB_split import Split
        split = Split(P["amo_full"], P["data_obs"], P["theta_truth"])

    # ---- 1. F-stat seed (truth-blind) ----
    t0 = time.time()
    scan = TE.seed_scan(P)
    seeds = TE.pick_seeds(scan, twoF_thresh=25.0)
    nrec = TE.report_seed_recovery(P, seeds, 25.0)
    src = TE.seeds_to_phi(P, seeds)[:]                          # 8*N_CW source vector
    print(f"  F-stat seed: {len(seeds)} sources, loud recovered {nrec}/3, {time.time()-t0:.0f}s", flush=True)

    # ---- 2. soft-EM alternation (staged LR = coarse->fine) ----
    base_dist = P["L0"].copy()                                  # distances start at prior mean
    LR = lambda it: 1.5e-2 if it < 8 else (6e-3 if it < 18 else 2e-3)
    hist = []
    q_prev = -np.inf; qt_prev = -np.inf; plateau = 0
    status = "budget"
    for it in range(max_iter):
        theta = np.concatenate([base_dist, src])
        post, LNL = soft_estep(P, theta, prior, split=split)    # E-step
        base_dist = post["map_evalL"].copy()                    # distances -> MAP (soft mixture centres kept in M-step)
        anchors = int(np.sum(post["map_qmax"] > 0.9))
        Pt = qtrue_trace(P, LNL, prior)                         # DIAGNOSTIC
        qt_sum = float(np.sum(Pt)); cen = [float(Pt[i]) for i in census_idx]
        off_max, _ = source_offset_from_truth(P, src)           # DIAGNOSTIC
        # M-step: ascend Q on sources (hard-gated to emergent anchors)
        dcfg, w, anchors = build_mstep_configs(P, post, base_dist)
        src, qval = mstep(src, dcfg, w, 40, LR(it))
        ncert = int(np.sum(Pt > 0.9))
        cert_census = all(Pt[census_idx[j]] > 0.9 for j in range(3))
        hist.append(dict(it=it, Q=qval, qt_sum=qt_sum, anchors=anchors, off_max=off_max,
                         ncert=ncert, cen=cen))
        print(f"  it{it:2d} Q={qval:.1f} sumQtrue={qt_sum:.2f} anchors={anchors} "
              f"srcOff={off_max:.4f}deg ncert={ncert} census[{cen[0]:.2f},{cen[1]:.2f},{cen[2]:.2f}] "
              f"lr={LR(it):.0e}", flush=True)
        # stall rule
        dq = abs(qval - q_prev) / (abs(q_prev) + 1e-9); dqt = abs(qt_sum - qt_prev) / (abs(qt_prev) + 1e-9)
        plateau = plateau + 1 if (dq < 1e-3 and dqt < 1e-3) else 0
        q_prev = qval; qt_prev = qt_sum
        if cert_census:
            status = "certified"; print(f"  -> census-3 certified at it{it}", flush=True); break
        if plateau >= 3:
            status = "STALL"; print(f"  -> STALL (Q & sumQtrue plateau 3 iters, srcOff={off_max:.4f}deg)", flush=True); break

    # ---- 3. certify (final E-step) ----
    theta = np.concatenate([base_dist, src])
    post, LNL = soft_estep(P, theta, prior, split=split)
    ev = TE.eval_gate(P, LNL, prior); Pt = ev["P_true"]
    print(f"\n  === B0.2 GATE (truth-blind) ===", flush=True)
    ok = True
    for nm, ceil in CEILING.items():
        a = names.index(nm); diff = abs(Pt[a] - ceil); ok &= diff < 0.05
        print(f"    {nm:14s} P_true={Pt[a]:.3f} ceiling={ceil:.3f} |diff|={diff:.3f}", flush=True)
    extra = sorted(set(ev["certified"]) - set(CEILING))
    print(f"    certified(P>0.9): {ev['certified']}", flush=True)
    print(f"    false certs (beyond census-3): {extra}", flush=True)
    off_max, _ = source_offset_from_truth(P, src)

    # ---- FLOAT LANDING ZONE (LAMBDA's search-space spec, STEP 1b deliverable) ----
    # per-pulsar integer support: # fringes carrying >1% posterior weight at convergence.
    supp = np.array([int(np.sum(post["qlist"][p][3] > 0.01)) for p in range(P["npsr"])])
    carriers = np.where(post["map_qmax"] > 0.9)[0]         # concentrated (near-fixed) pulsars
    census_supp = {names[i]: int(supp[i]) for i in census_idx}
    print(f"\n  === FLOAT LANDING ZONE (LAMBDA search-space spec) ===", flush=True)
    print(f"    per-pulsar integer support (#fringes q>1%): median {int(np.median(supp))}, "
          f"min {supp.min()}, max {supp.max()}; #pulsars with support==1 (locked) "
          f"{int(np.sum(supp==1))}, ==2 {int(np.sum(supp==2))}", flush=True)
    print(f"    concentrated carriers (q_max>0.9): {len(carriers)} pulsars; census-3 support {census_supp}", flush=True)
    print(f"    final source offset from truth: {off_max:.4f} deg (the float landing)", flush=True)

    passed = ok and set(ev["certified"]) == set(CEILING) and status != "STALL"
    verdict = "PASS" if passed else ("STALL" if status == "STALL" else "FAIL")
    # failure leg attribution
    leg = "n/a"
    if not passed:
        if nrec < 3:
            leg = "DESCENT (F-stat seed missed loud sources)"
        elif off_max > 0.02:
            leg = "IGNITION (soft-EM did not cross the gap; srcOff %.3f deg > tol)" % off_max
        else:
            leg = "LOCK (sources near truth but certification off)"
    print(f"\n  [P2 GATE] {verdict} | srcOff={off_max:.4f}deg | status={status} | failed-leg: {leg}", flush=True)
    if save:
        np.savez(f"{CWT}/trackB_p2_gate.npz",
                 P_true=Pt, names=np.array(names), certified=np.array(ev["certified"]),
                 src=src, src_truth=P["theta_truth"][ndist:], base_dist=base_dist,
                 hist_Q=np.array([h["Q"] for h in hist]),
                 hist_qtsum=np.array([h["qt_sum"] for h in hist]),
                 hist_anchors=np.array([h["anchors"] for h in hist]),
                 hist_off=np.array([h["off_max"] for h in hist]),
                 hist_ncert=np.array([h["ncert"] for h in hist]),
                 hist_census=np.array([h["cen"] for h in hist]),
                 int_support=supp, carriers=carriers, final_qmax=post["map_qmax"],
                 status=status, verdict=verdict, failed_leg=leg, nrec=nrec)
        print(f"  saved trackB_p2_gate.npz", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["test-mstep", "gate"])
    ap.add_argument("--max_iter", type=int, default=30)
    ap.add_argument("--use-split", action="store_true", help="E-step via trackB_split (~47x)")
    a = ap.parse_args()
    if a.mode == "test-mstep":
        sys.exit(mode_test_mstep())
    else:
        sys.exit(mode_gate(a.max_iter, use_split=a.use_split))
