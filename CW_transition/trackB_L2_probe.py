"""Track B — LAMBDA feasibility probe, DELIVERABLE L2: the GO/NO-GO probe.

Seeded by the float landing (trackB_L0_float.npz), at the FLOAT source params:
  1. STAGED INTEGER SEARCH via Split.score_candidates (census-3 lattice -> fix -> next tiers).
  2. CONDITIONAL SOURCE RE-SOLVE at winning integers (TE.make_mstep, integers FIXED).
  3. ITERATE ONCE (re-search at re-solved source).
  4. RATIO TEST / A2 P_true certification. Pre-registered gate.

Truth's fringe integer is n_true (=0 for all pulsars; true dist == prior mean L0). "matches
truth" = winner integer == n_true (DIAGNOSTIC ONLY, never fed to the search).

Candidate set per searched pulsar = float-neighborhood: integers within radius R of the
float MAP fringe UNION the >1%-mass integers. Reachability of n_true reported, never injected.

Run: python trackB_L2_probe.py --use-split   (jug-gpu, background w/ log)
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import argparse, sys, time, itertools
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE

CWT = "/home/mattm/projects/HSYMT/CW_transition"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]
SEARCH_R = 4          # integer radius around the float MAP fringe (RTK bounded neighborhood)
MASS_THR = 0.005      # also admit any integer with > this float posterior mass
PAD_TO = 1000         # pad each scorer call to >= this many candidates (overhead floor)


def candidate_ints(uk, w, K, map_k, R=SEARCH_R, mass=MASS_THR):
    """Float-neighborhood integer set for one pulsar: [map_k-R, map_k+R] U {>mass mass}."""
    uk = uk[:K]; w = w[:K]
    nb = set(range(int(map_k) - R, int(map_k) + R + 1))
    nb |= set(uk[w > mass].tolist())
    lo, hi = int(uk.min()), int(uk.max())
    return np.array(sorted(i for i in nb if lo <= i <= hi))


def build_search_EV(P, searched, cand, held_int, L0, dL):
    """EV_search (npsr,B): searched pulsars -> their candidate-integer distances (cols
    0..m-1, padded); fixed pulsars -> single held-integer distance at col 0 (padded).
    Returns EV, and col_int (npsr list of int arrays: column->integer)."""
    npsr = P["npsr"]
    B = max(max(len(cand[p]) for p in searched), 1)
    EV = np.empty((npsr, B)); col_int = [None] * npsr
    for p in range(npsr):
        if p in searched:
            ints = cand[p]
        else:
            ints = np.array([int(held_int[p])])
        dists = L0[p] + ints * dL[p]
        row = np.full(B, dists[-1]); row[:len(dists)] = dists
        ci = np.full(B, ints[-1], int); ci[:len(ints)] = ints
        EV[p] = row; col_int[p] = ci
    return EV, col_int


def enumerate_assignments(P, searched, cand, held_col):
    """Cartesian product over searched pulsars' candidate COLUMN indices; fixed pulsars
    pinned to held_col[p]. Returns assignments (K,npsr) of column indices."""
    npsr = P["npsr"]
    base = np.array([held_col[p] for p in range(npsr)], int)
    grids = [range(len(cand[p])) for p in searched]
    combos = list(itertools.product(*grids))
    K = len(combos)
    A = np.repeat(base[None, :], K, axis=0)
    sp = list(searched)
    for r, combo in enumerate(combos):
        for j, p in enumerate(sp):
            A[r, p] = combo[j]
    return A


def pad_assignments(A, pad_to=PAD_TO):
    """Repeat the first row to reach >= pad_to (overhead-floor batching); return padded A
    and the true count K0 (only first K0 are distinct/meaningful)."""
    K0 = A.shape[0]
    if K0 >= pad_to:
        return A, K0
    reps = int(np.ceil(pad_to / K0))
    return np.repeat(A, reps, axis=0)[:max(pad_to, K0)], K0


def score_stage(split, source, EV, AI, A, K0, col_int, n_true, searched):
    """build_tables at `source`, score assignments, return winner/runner/margin + truth match."""
    split.build_tables(np.asarray(source), EV, AI)
    Apad, _ = pad_assignments(A)
    lnL = split.score_candidates(Apad)[:K0]
    order = np.argsort(-lnL)
    w0, w1 = order[0], order[1] if K0 > 1 else order[0]
    win_cols = A[w0]; run_cols = A[w1]
    win_int = {p: int(col_int[p][win_cols[p]]) for p in searched}
    run_int = {p: int(col_int[p][run_cols[p]]) for p in searched}
    margin = float(lnL[w0] - lnL[w1])
    matches_truth = all(win_int[p] == int(n_true[p]) for p in searched)
    return dict(win_int=win_int, run_int=run_int, margin=margin, lnL_win=float(lnL[w0]),
                matches_truth=matches_truth, win_cols=win_cols, K0=K0)


def main(use_split=True):
    print("=== L2 LAMBDA PROBE (staged integer search + re-solve + ratio test) ===", flush=True)
    t0 = time.time()
    fl = np.load(f"{CWT}/trackB_L0_float.npz", allow_pickle=True)
    src_float = fl["src"].copy(); base_dist_float = fl["base_dist"].copy()
    q_uk, q_w, q_K = fl["q_uk"], fl["q_w"], fl["q_K"]
    qmax = fl["qmax"]; map_k = fl["map_k"]; n_true = fl["n_true"]
    census_idx = list(fl["census_idx"]); float_status = str(fl["status"])
    print(f"  loaded float (status={float_status}); loud scaled off max="
          f"{fl['loud_scaled_off'].max():.4f}", flush=True)

    P = TE.build_problem("pop", build_earth=True); names = list(P["names"])
    L0 = np.asarray(P["L0"]); dL = np.asarray(P["dL"]); AI = P["AI"]; npsr = P["npsr"]
    ndist = P["n_dist"]; prior = P["priors"]["lit"]
    from trackB_split import Split
    split = Split(P["amo_full"], P["data_obs"], P["theta_truth"])
    # CONVERGENT conditional re-solve: LBFGS on -lnL over source (distances FIXED). Adam is
    # WRONG here -- its normalised step m/sqrt(v) does not vanish at a flat optimum, so it
    # random-walks ~lr*scale*nsteps away instead of locking. LBFGS with analytic grad -> grad=0.
    import scipy.optimize as sopt
    _fl = P["amo_full"]["fisher_logL"]; _data = P["data_obs"]
    _plo, _phi = TE.phi_bounds(P)
    _vg = jax.jit(lambda baseL, s: jax.value_and_grad(
        lambda ss: -_fl(jnp.concatenate([baseL, ss]), _data))(s))

    # candidate integer sets (float-neighborhood)
    cand = {p: candidate_ints(q_uk[p], q_w[p], q_K[p], map_k[p]) for p in range(npsr)}
    # carrier ordering: qmax desc; census forced into stage 1
    carrier_order = [int(a) for a in np.argsort(-qmax) if qmax[a] > 0.9]
    for ci in census_idx:
        if ci not in carrier_order:
            carrier_order.insert(0, int(ci))
    # reachability report
    print("  census-3 candidate sets (float source):", flush=True)
    for ci in census_idx:
        reach = int(n_true[ci]) in cand[ci].tolist()
        print(f"    {names[ci]:14s} map_k={map_k[ci]} qmax={qmax[ci]:.3f} "
              f"cand={cand[ci].tolist()} n_true={n_true[ci]} reachable={reach}", flush=True)

    # ---- staged search plan: stage1 = census-3; then tiers of ~<=6 carriers padded ----
    rest = [a for a in carrier_order if a not in census_idx]
    stages = [list(census_idx)] + [rest[i:i + 5] for i in range(0, len(rest), 5)]
    stages = [s for s in stages if s]
    print(f"  {len(stages)} stages; sizes {[len(s) for s in stages]}", flush=True)

    def run_search(source, tag):
        held_int = {p: int(map_k[p]) for p in range(npsr)}   # start: everyone at float MAP
        fixed = {}                                           # p -> winner integer (fixed)
        results = []
        for si, stage in enumerate(stages):
            searched = set(stage)
            EV, col_int = build_search_EV(P, searched, cand, held_int, L0, dL)
            held_col = np.zeros(npsr, int)                   # fixed pulsars use col 0 (=held_int)
            A = enumerate_assignments(P, searched, cand, held_col)
            r = score_stage(split, source, EV, AI, A, A.shape[0], col_int, n_true, searched)
            for p in searched:                               # fix winners -> held for next stage
                held_int[p] = r["win_int"][p]; fixed[p] = r["win_int"][p]
            tru = [p for p in searched if r["win_int"][p] == int(n_true[p])]
            print(f"  [{tag}] stage{si} search={sorted(searched)} K={r['K0']} "
                  f"margin={r['margin']:.3e} matchesTruth={r['matches_truth']} "
                  f"({len(tru)}/{len(searched)} at n_true); "
                  f"census_win={[r['win_int'].get(c,'-') for c in census_idx]}", flush=True)
            results.append(r)
        return fixed, results

    from trackB_p1_map import ang_deg
    scale = TE.phi_scale(P); src_t = P["theta_truth"][ndist:]
    n_loud = TE.CONFIGS["pop"]["population"][0]; loud_sl = slice(0, 8 * n_loud)

    def resolve(source0, fixed_map, maxiter=300):
        """Conditional source re-solve: distances FIXED at fixed_map integers (rest at float
        MAP), LBFGS-maximise full lnL over source -> the smooth-quadratic 'fixed' solution.
        Returns src, sky-deg offset, scaled off max."""
        base_L = base_dist_float.copy()
        for p in range(npsr):
            ni = fixed_map.get(p, int(map_k[p])); base_L[p] = L0[p] + ni * dL[p]
        baseLj = jnp.asarray(base_L)
        def fun(s):
            v, g = _vg(baseLj, jnp.asarray(s))
            return float(v), np.asarray(g, float)
        res = sopt.minimize(fun, np.asarray(source0, float), jac=True, method="L-BFGS-B",
                            bounds=list(zip(_plo, _phi)),
                            options=dict(maxiter=maxiter, ftol=1e-12, gtol=1e-9))
        src_o = res.x
        off = np.abs(src_o[loud_sl] - src_t[loud_sl]) / scale[loud_sl]
        skyd = max(ang_deg(src_o[8*k], src_o[8*k+1], src_t[8*k], src_t[8*k+1]) for k in range(n_loud))
        return src_o, skyd, off

    # ---- ORACLE reachability probe (decisive anatomy; DIAGNOSTIC, not part of the blind search) ----
    # Q: is TRUTH (all carriers at n_true=0) the optimum the blind search SHOULD find? Score the
    # all-zeros carrier assignment at the float source, and re-solve from it. If re-solve locks to
    # truth (mas) with higher lnL than the float-winner -> truth IS reachable-in-principle and any
    # FAIL is a search-radius problem. If not -> truth is not the optimum at the float source (deeper).
    print(f"\n  === ORACLE reachability (truth integers, DIAGNOSTIC) === ({time.time()-t0:.0f}s)", flush=True)
    # ALL pulsars' true fringe = n_true (=0); fix the whole array at truth, re-solve source.
    truth_fix = {p: int(n_true[p]) for p in range(npsr)}
    # (a) pull-in from a MODERATE offset (truth loud +0.1 scaled): does the conditional re-solve
    #     LOCK to truth when the integers are right? Tests the P1-needle fixed point.
    rng = np.random.default_rng(11)
    src_mod = src_t.copy()
    src_mod[loud_sl] = src_t[loud_sl] + rng.normal(size=8*n_loud) * scale[loud_sl] * 0.1
    src_a, sky_a, off_a = resolve(src_mod, truth_fix)
    print(f"    (a) from truth+0.1scaled, truth integers -> sky={sky_a:.5f}deg off_max={off_a.max():.4f} "
          f"(fixed-point/pull-in test; cert tol ~0.006deg)", flush=True)
    # (b) pull-in from the ACTUAL float source (the real landing): does truth remain reachable?
    src_ora, sky_ora, off_ora = resolve(src_float, truth_fix)
    print(f"    (b) from FLOAT source, truth integers -> sky={sky_ora:.5f}deg off_max={off_ora.max():.4f} "
          f"(pull-in from the actual float landing)", flush=True)

    print(f"\n  --- SEARCH ROUND 1 (at FLOAT source) --- ({time.time()-t0:.0f}s)", flush=True)
    fixed1, res1 = run_search(src_float, "R1")

    # ---- conditional source re-solve: distances FIXED at winner integers, ascend source ----
    src_re, sky_deg, off_re = resolve(src_float, fixed1)
    print(f"\n  --- RE-SOLVE at fixed integers: loud scaled off max={off_re.max():.4f} "
          f"median={np.median(off_re):.4f}; sky={sky_deg:.4f}deg (float was "
          f"{fl['loud_scaled_off'].max():.4f}) --- ({time.time()-t0:.0f}s)", flush=True)

    print(f"\n  --- SEARCH ROUND 2 (at RE-SOLVED source) --- ({time.time()-t0:.0f}s)", flush=True)
    fixed2, res2 = run_search(src_re, "R2")

    # ---- ratio test / certification at final fixed solution ----
    base_Lf = base_dist_float.copy()
    for p, ni in fixed2.items():
        base_Lf[p] = L0[p] + ni * dL[p]
    theta_f = np.concatenate([base_Lf, src_re])
    LNL = split.estep(theta_f, P["EV"], AI)
    ev = TE.eval_gate(P, LNL, prior); Pt = ev["P_true"]
    print(f"\n  === RATIO TEST / A2 CERTIFICATION (final fixed solution) ===", flush=True)
    for nm, ceil in TE.CEILING.items():
        a = names.index(nm); diff = abs(Pt[a] - ceil)
        wi = fixed2.get(a, int(map_k[a]))
        print(f"    {nm:14s} P_true={Pt[a]:.3f} ceiling={ceil:.3f} |diff|={diff:.3f} "
              f"win_int={wi} (n_true={n_true[a]})", flush=True)
    false_fix = [names[p] for p, ni in fixed2.items() if ni != int(n_true[p]) and res2 and True]
    print(f"    certified(P>0.9): {ev['certified']}", flush=True)

    # per-carrier dlnL(true fringe vs best wrong) from the E-step
    print(f"    per-census dlnL(n_true vs best-wrong) at final:", flush=True)
    for ci in census_idx:
        uk = q_uk[ci][:q_K[ci]]
        # eval LNL columns correspond to P['EV'] grid; recompute integer of each EV col
        evL = P["EV"][ci]; kcol = np.round((evL - L0[ci]) / dL[ci]).astype(int)
        m_true = np.where(kcol == int(n_true[ci]))[0]
        if m_true.size:
            lt = LNL[ci][m_true[0]]
            wrong = LNL[ci][kcol != int(n_true[ci])]
            dd = float(lt - wrong.max()) if wrong.size else np.inf
        else:
            dd = np.nan
        print(f"      {names[ci]:14s} dlnL={dd:.4f}", flush=True)

    # ---- verdict (pre-registered) ----
    cen_ok = all(abs(Pt[names.index(nm)] - c) < 0.05 for nm, c in TE.CEILING.items())
    census_true = all(fixed2.get(names.index(nm), 1) == int(n_true[names.index(nm)]) for nm in CENSUS)
    false_fix_carrier = [names[p] for p in carrier_order
                         if p in fixed2 and fixed2[p] != int(n_true[p]) and res2[0]["margin"] > 0]
    src_locked = sky_deg < 0.006
    if cen_ok and census_true and src_locked and set(ev["certified"]) == set(TE.CEILING):
        verdict = "PASS"
    elif census_true or any(fixed2.get(c) == int(n_true[c]) for c in census_idx):
        verdict = "PARTIAL"
    else:
        verdict = "FAIL"
    fixedpoint_ok = sky_a < 0.006          # truth is a stable conditional fixed point
    oracle_locks = sky_ora < 0.006         # truth reachable from the ACTUAL float
    if fixedpoint_ok and not oracle_locks:
        anatomy = ("FLOAT-TOO-FAR: truth IS a stable conditional fixed point (locks from "
                   "+0.1 scaled), but the actual float landing is outside its pull-in basin "
                   "-> bounded integer search from this float cannot reach truth's integers")
    elif not fixedpoint_ok:
        anatomy = ("NEEDLE-ABSENT/UNSTABLE: even from +0.1 scaled with TRUE integers fixed the "
                   "conditional re-solve does not lock -> the smooth-quadratic claim fails here")
    else:
        anatomy = "truth reachable from the float; any blind failure is search-radius only"
    print(f"\n  [L2 VERDICT] {verdict} | census_true_int={census_true} src_locked={src_locked} "
          f"({sky_deg:.4f}deg) census_P_ok={cen_ok}", flush=True)
    print(f"    ORACLE: fixedpoint(+0.1)->sky={sky_a:.5f}deg locks={fixedpoint_ok}; "
          f"float->sky={sky_ora:.5f}deg locks={oracle_locks}", flush=True)
    print(f"    ANATOMY: {anatomy}", flush=True)

    np.savez(f"{CWT}/trackB_L2_probe.npz",
             fixed_R1=np.array([[p, ni] for p, ni in sorted(fixed1.items())]),
             fixed_R2=np.array([[p, ni] for p, ni in sorted(fixed2.items())]),
             src_float=src_float, src_re=src_re, off_re=off_re, sky_deg=sky_deg,
             P_true=Pt, certified=np.array(ev["certified"]), names=np.array(names),
             census_idx=np.array(census_idx), n_true=n_true, verdict=verdict,
             sky_ora=sky_ora, off_ora=off_ora, oracle_locks=oracle_locks,
             sky_fixedpoint=sky_a, fixedpoint_ok=fixedpoint_ok, anatomy=anatomy,
             R1_margins=np.array([r["margin"] for r in res1]),
             R2_margins=np.array([r["margin"] for r in res2]))
    print(f"\n  saved trackB_L2_probe.npz ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("--use-split", action="store_true")
    a = ap.parse_args()
    sys.exit(main(use_split=a.use_split))
