"""Track B — low-rank QuickCW split (spec sec 6/6.1). The split PRIMITIVE + E-step entry
point + candidate scorer. Batch-only (`trackB_p1_map.batched_scan`) stays the correctness
ORACLE; this module never shares code with it (Amendment 6).

Physics (verified, spec 6.1): fisher_logL is a two-level Woodbury. Distances enter only via
per-pulsar-LOCAL a_p (scalar), b_p (n_gwb vector); the GWB HD-ORF coupling is the FROZEN
outer factor M=cf_cached. Sweeping pulsar p's fringe changes ONLY a_p, b_p:
    lnL(theta) = sum_p a_dep_p  +  a_const  +  0.5[ b.flat^T M^-1 b.flat - ldP - logdet ]
    a_dep_p = -0.5( y^T N^-1 y_p  -  FtNmy_rn_p . cf_p^-1 FtNmy_rn_p )
    b_p     = T^T N^-1 y_gwb,p    -  TtNmF_p @ (cf_p^-1 FtNmy_rn_p)
cf_p (per-pulsar rednoise factor), TtNmF_p, a_const are FROZEN (noise/rednoise only).

Build order (this file grows through the gates): (1) primitive + scalar recon; (2) E-step
values+soft-posterior gates; (3) gradient gate; (4) small-signal; (5) profiling; (6)
candidate scorer. Reuses discovery.matrix primitives directly (dsm.*) so the algebra is
bit-identical to the reference by construction.

Env: jug-gpu, XLA autotune off, x64. Matt commits; gate + stop after each stage.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
import discovery.matrix as dsm

CWT = "/home/mattm/projects/HSYMT/CW_transition"
B_CERT = 512          # certification fringe count -- NEVER less (asserted at the cert entry)
B_SEARCH = 512        # search phases MAY lower this


class Split:
    """Frozen setup + per-pulsar (a_dep, b) evaluator, reconstructed from amo internals by
    mirroring discovery.matrix.make_kernelterms_vary (line ~2631). a_const is fixed by
    difference vs fisher_logL at truth (sidesteps re-deriving the rednoise logdet)."""

    def __init__(self, amo, data_obs, theta_truth):
        it = amo["internals"]
        self.amo = amo; self.data = data_obs
        self.theta_keys = amo["theta_keys"]; self.n_dist = amo["n_dist"]
        self.frozen = it["frozen"]; self.data_keys = it["data_keys"]
        self.ys = it["ys"]
        self.cf_cached = it["cf_cached"]; self.ldP_cached = it["ldP_sum_cached"]
        self.logdet_cached = it["logdet_cf_cached"]
        self.npsr = it["npsr"]; self.ngp_gwb = it["ngp_gwb"]
        vsm = it["vsm"]
        Ns, Fs_rn, Ts = vsm.Ns, vsm.Fs, it["Fs_gwb"]

        # --- frozen per-pulsar rednoise Woodbury (mirror make_kernelterms_vary 2632-2662) ---
        NmFs, ldNs = zip(*[N.solve_2d(F) for N, F in zip(Ns, Fs_rn)])
        if dsm.regularize_FtNmF:
            FtNmFs = [dsm.zero_out_negative_eigenvalues(F.T @ NmF) for F, NmF in zip(Fs_rn, NmFs)]
        else:
            FtNmFs = [F.T @ NmF for F, NmF in zip(Fs_rn, NmFs)]
        TtNmFs = [T.T @ NmF for T, NmF in zip(Ts, NmFs)]
        self.FtNmF = jnp.asarray(np.array(FtNmFs))          # (npsr, ncomp_rn, ncomp_rn)
        self.TtNmF = jnp.asarray(np.array(TtNmFs))          # (npsr, ngp_gwb, ncomp_rn)
        self.ldN = float(sum(ldNs))
        self.solve1d = [N.make_solve_1d() for N in Ns]
        self.Ft = [jnp.asarray(np.asarray(F.T)) for F in Fs_rn]   # (ncomp_rn, n_toa)
        self.Tt = [jnp.asarray(np.asarray(T.T)) for T in Ts]      # (ngp_gwb, n_toa)

        # frozen rednoise factor cf (depends only on rednoise hyperparams, in `frozen`)
        Pinv, ldP = vsm.P_var.make_inv()(self.frozen)
        if Pinv.ndim == 2:
            i1, i2 = jnp.diag_indices(self.FtNmF.shape[1], ndim=2)
            self.cf = dsm.matrix_factor(self.FtNmF.at[:, i1, i2].add(Pinv))
        else:
            self.cf = dsm.matrix_factor(self.FtNmF + Pinv)

        # a_const by difference at truth
        a_dep_t, b_t = self.per_pulsar_ab(jnp.asarray(theta_truth))
        rest = float(jnp.sum(a_dep_t)) + 0.5 * float(
            b_t.reshape(-1) @ dsm.matrix_solve(self.cf_cached, b_t.reshape(-1))
            - self.ldP_cached - self.logdet_cached)
        self.a_const = float(amo["fisher_logL"](jnp.asarray(theta_truth), data_obs)) - rest

    def _params(self, theta_arr):
        params = dict(self.frozen)
        for k, v in zip(self.theta_keys, theta_arr):
            params[k] = v
        for dk, d in zip(self.data_keys, self.data):
            params[dk] = d
        return params

    def per_pulsar_ab(self, theta_arr):
        """(a_dep (npsr,), b (npsr, ngp_gwb)) for a full theta. Differentiable in theta;
        cf/cf_cached are frozen constants (NOT in the differentiated graph)."""
        params = self._params(theta_arr)
        ytNmy = []; FtNmy = []; TtNmy = []
        for p in range(self.npsr):
            yp = self.ys[p](params)
            Nmy, _ = self.solve1d[p](yp)
            ytNmy.append(yp @ Nmy)
            FtNmy.append(self.Ft[p] @ Nmy)
            TtNmy.append(self.Tt[p] @ Nmy)
        ytNmy = jnp.stack(ytNmy)                              # (npsr,)
        FtNmy = jnp.stack(FtNmy)                              # (npsr, ncomp_rn)
        TtNmy = jnp.stack(TtNmy)                              # (npsr, ngp_gwb)
        sol = dsm.matrix_solve(self.cf, FtNmy)               # (npsr, ncomp_rn)
        a_dep = -0.5 * (ytNmy - jnp.sum(FtNmy * sol, axis=1))
        b = TtNmy - jnp.sum(self.TtNmF * sol[:, jnp.newaxis, :], axis=2)
        return a_dep, b

    def lnL(self, theta_arr):
        a_dep, b = self.per_pulsar_ab(theta_arr)
        bflat = b.reshape(-1)
        return (jnp.sum(a_dep) + self.a_const
                + 0.5 * (bflat @ dsm.matrix_solve(self.cf_cached, bflat)
                         - self.ldP_cached - self.logdet_cached))

    # ---- fast E-step: sweep pulsar p over its fringe grid via rank-block update ----
    def _ensure_minv(self):
        if getattr(self, "_Minv", None) is None:
            ng = self.ngp_gwb; M = self.npsr * ng
            self._Minv = np.asarray(dsm.matrix_solve(self.cf_cached, jnp.eye(M)))   # (M,M)
            self._Minv_pp = np.stack([self._Minv[p * ng:(p + 1) * ng, p * ng:(p + 1) * ng]
                                      for p in range(self.npsr)])                   # (npsr,ng,ng)

    def _pulsar_ab_fn(self, p):
        """Cached jitted fn (theta_base, Lvals) -> (a_dep(B,), b(B,ngp)) for pulsar p;
        theta_base is a RUNTIME arg so the compile is reused across E-steps."""
        if getattr(self, "_ab_fns", None) is None:
            self._ab_fns = {}
        if p in self._ab_fns:
            return self._ab_fns[p]
        Ft = self.Ft[p]; Tt = self.Tt[p]; s1d = self.solve1d[p]; ys_p = self.ys[p]
        cf_p = (self.cf[0][p], self.cf[1]); TtNmF_p = self.TtNmF[p]; AI_p = int(self.AI[p])

        @jax.jit
        def run(theta_base, Lvals):
            def one(L):
                params = self._params(theta_base.at[AI_p].set(L))
                yp = ys_p(params); Nmy, _ = s1d(yp)
                return yp @ Nmy, Ft @ Nmy, Tt @ Nmy
            ytNmy, FtNmy, TtNmy = jax.vmap(one)(Lvals)              # (B,),(B,nc),(B,ngp)
            sol = dsm.matrix_solve(cf_p, FtNmy.T).T                 # (B,nc)
            a_dep = -0.5 * (ytNmy - jnp.sum(FtNmy * sol, axis=1))
            b = TtNmy - sol @ TtNmF_p.T                            # (B,ngp)
            return a_dep, b
        self._ab_fns[p] = run
        return run

    def estep(self, theta_base, EV, AI):
        """LNL (npsr,B): each pulsar's distance swept over EV[p], others at theta_base.
        Matches trackB_p1_map.batched_scan via the rank-block update."""
        self.AI = np.asarray(AI); self._ensure_minv()
        tb = jnp.asarray(theta_base)
        a_base, b_base = self.per_pulsar_ab(tb)
        a_base = np.asarray(a_base); b_base = np.asarray(b_base)      # (npsr,), (npsr,ngp)
        bflat = b_base.reshape(-1)
        u = self._Minv @ bflat                                       # (M,)
        Q_base = float(bflat @ u)
        u_p = u.reshape(self.npsr, self.ngp_gwb)                     # (npsr,ngp)
        sum_a = float(a_base.sum())
        const = self.a_const - 0.5 * (self.ldP_cached + self.logdet_cached)
        npsr, B = EV.shape
        LNL = np.empty((npsr, B))
        for p in range(npsr):
            a_pf, b_pf = self._pulsar_ab_fn(p)(tb, jnp.asarray(EV[p]))
            a_pf = np.asarray(a_pf); b_pf = np.asarray(b_pf)         # (B,),(B,ng)
            db = b_pf - b_base[p]                                    # (B,ng)
            Q_pf = Q_base + 2.0 * (db @ u_p[p]) + np.einsum('bi,ij,bj->b', db, self._Minv_pp[p], db)
            LNL[p] = const + (sum_a - a_base[p] + a_pf) + 0.5 * Q_pf
        return LNL

    # ---- candidate scorer (Amendment 4): score integer fringe assignments at FIXED source ----
    def build_tables(self, source_params, EV, AI):
        """Per-(pulsar,fringe) a/b tables at FIXED source: a_table (npsr,B), b_table
        (npsr,B,ngp). One-time single-pulsar work; candidates then need no rebuilds/solves."""
        self.AI = np.asarray(AI); self._ensure_minv()
        ndist = self.n_dist
        tb = np.zeros(len(self.theta_keys)); tb[ndist:] = source_params
        tb[:ndist] = np.asarray(EV)[:, 0]                          # arbitrary base distances
        tb = jnp.asarray(tb)
        npsr, B = EV.shape
        aT = np.empty((npsr, B)); bT = np.empty((npsr, B, self.ngp_gwb))
        for p in range(npsr):
            a_pf, b_pf = self._pulsar_ab_fn(p)(tb, jnp.asarray(EV[p]))
            aT[p] = np.asarray(a_pf); bT[p] = np.asarray(b_pf)
        self._aT = aT; self._bT = bT
        self._const = self.a_const - 0.5 * (self.ldP_cached + self.logdet_cached)
        return aT, bT

    def score_candidates(self, assignments, chunk=20000):
        """assignments (K,npsr) int fringe indices -> lnL (K,). Pure precomputed-table
        arithmetic (gather + one quadratic form): no likelihood rebuilds, no solves in the
        K loop. Requires build_tables() first. Chunked on the GPU."""
        assignments = np.asarray(assignments)
        K, npsr = assignments.shape
        pr = np.arange(npsr)[None, :]
        if getattr(self, "_Minv_j", None) is None:
            self._Minv_j = jnp.asarray(self._Minv)
        out = np.empty(K)
        for c0 in range(0, K, chunk):
            asg = assignments[c0:c0 + chunk]
            a_sel = self._aT[pr, asg]                              # (k,npsr)
            b_sel = self._bT[pr, asg]                              # (k,npsr,ng)
            bflat = jnp.asarray(b_sel.reshape(len(asg), npsr * self.ngp_gwb))
            quad = jnp.sum(bflat * (bflat @ self._Minv_j), axis=1)
            out[c0:c0 + chunk] = self._const + a_sel.sum(1) + 0.5 * np.asarray(quad)
        return out


# ============================================================
# GATE 0: scalar reconstruction  split.lnL == amo.fisher_logL
# ============================================================
def mode_recon():
    print("=== GATE 0: split scalar lnL == fisher_logL ===", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False)
    print(f"  build_problem(build_earth=False) {time.time()-t0:.0f}s", flush=True)
    amo = P["amo_full"]; data = P["data_obs"]; tt = P["theta_truth"]; ndist = P["n_dist"]
    sp = Split(amo, data, tt)
    fl = amo["fisher_logL"]
    print(f"  a_const = {sp.a_const:.6f}", flush=True)

    # positions: truth, two plateau (~0.5,1 deg sky offset on loud0 gwphi), 5 random
    rng = np.random.default_rng(0)
    scale = TE.phi_scale(P)
    thetas = [tt.copy()]
    from trackB_p1_map import _local_sky_metric
    _, m_gp = _local_sky_metric(tt[ndist + 0])
    for off in (0.5, 1.0):
        th = tt.copy(); th[ndist + 1] += np.radians(off) / m_gp; thetas.append(th)
    plo, phi_hi = TE.phi_bounds(P)
    for _ in range(5):
        th = tt.copy()
        src = th[ndist:] + rng.normal(size=len(tt) - ndist) * scale * 0.05
        th[ndist:] = np.clip(src, plo, phi_hi)               # keep waveform valid
        thetas.append(th)
    labels = ["truth", "plateau0.5deg", "plateau1deg"] + [f"rand{i}" for i in range(5)]

    print(f"  {'position':14s} {'fisher_logL':>16s} {'split.lnL':>16s} {'|abs diff|':>12s} {'rel':>10s}", flush=True)
    worst = 0.0
    for lab, th in zip(labels, thetas):
        a = float(fl(jnp.asarray(th), data)); b = float(sp.lnL(jnp.asarray(th)))
        d = abs(a - b); worst = max(worst, d)
        print(f"  {lab:14s} {a:16.6f} {b:16.6f} {d:12.2e} {d/abs(a):10.2e}", flush=True)
    ok = worst < 1e-8
    print(f"\n  [GATE 0] max |diff| = {worst:.2e} -> {'PASS' if ok else 'FAIL'} (<1e-8)", flush=True)
    return 0 if ok else 1


def mode_estep():
    print("=== GATE 1: split.estep == batched_scan (per-pulsar + soft posteriors) ===", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False)
    print(f"  build_problem {time.time()-t0:.0f}s", flush=True)
    amo = P["amo_full"]; data = P["data_obs"]; tt = P["theta_truth"]; ndist = P["n_dist"]
    sp = Split(amo, data, tt)
    from trackB_p1_map import batched_scan, _local_sky_metric
    _, m_gp = _local_sky_metric(tt[ndist + 0])
    positions = {"truth": tt.copy()}
    th = tt.copy(); th[ndist + 1] += np.radians(0.5) / m_gp; positions["plateau0.5deg"] = th

    ok = True
    for lab, tb in positions.items():
        ref, _ = batched_scan(P, [tb]); ref = ref[0]                # (npsr,B)
        _ = sp.estep(tb, P["EV"], P["AI"])                         # warm (per-pulsar compiles)
        t1 = time.time(); split = sp.estep(tb, P["EV"], P["AI"]); dt = time.time() - t1
        perpsr = np.nanmax(np.abs(ref - split), axis=1)
        worst = float(np.nanmax(perpsr)); wp = int(np.nanargmax(perpsr))
        nbad = int(np.sum(perpsr > 1e-8))
        vok = worst < 1e-8; ok &= vok
        print(f"  [{lab}] VALUES: per-pulsar max|dlnL| worst={worst:.2e} @psr{wp}; "
              f"#pulsars>1e-8={nbad}/116 -> {'PASS' if vok else 'FAIL'}; warm estep {dt:.2f}s", flush=True)
        pr = TE.fringe_posterior(P, ref, None, P["priors"]["lit"], 1.0)
        ps = TE.fringe_posterior(P, split, None, P["priors"]["lit"], 1.0)
        dq = max(float(np.max(np.abs(pr["qlist"][p][3] - ps["qlist"][p][3]))) for p in range(P["npsr"]))
        qok = dq < 1e-10; ok &= qok
        print(f"  [{lab}] SOFT POSTERIOR: max|dq_p(n)|={dq:.2e} -> {'PASS' if qok else 'FAIL'} (<1e-10)", flush=True)
    print(f"\n  [GATE 1] {'PASS' if ok else 'FAIL'}", flush=True)
    return 0 if ok else 1


# ============================================================
# GATE 2/3 + profiling + GATE 4 (scorer) -- one build (Amendments 2,3,5,4)
# ============================================================
def mode_post():
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False)
    print(f"  build_problem {time.time()-t0:.0f}s", flush=True)
    amo = P["amo_full"]; data = P["data_obs"]; tt = P["theta_truth"]; ndist = P["n_dist"]
    EV = P["EV"]; AI = P["AI"]; sp = Split(amo, data, tt)
    from trackB_p1_map import batched_scan, _local_sky_metric
    fl = amo["fisher_logL"]; flb = amo["fisher_logL_batched"]
    ok = True

    # --- GATE 2: gradient wrt source params ---
    print("\n=== GATE 2: source-gradient split vs batch-only ===", flush=True)
    g_ref_f = jax.jit(jax.grad(lambda th: fl(th, data)))
    g_spl_f = jax.jit(jax.grad(sp.lnL))
    _, m_gp = _local_sky_metric(tt[ndist + 0])
    thp = tt.copy(); thp[ndist + 1] += np.radians(0.5) / m_gp
    for lab, th in [("truth", tt), ("plateau0.5deg", thp)]:
        gr = np.asarray(g_ref_f(jnp.asarray(th)))[ndist:]
        _ = g_spl_f(jnp.asarray(th))
        t1 = time.time(); gs = np.asarray(g_spl_f(jnp.asarray(th)))[ndist:]; tg = time.time() - t1
        m = np.abs(gr) > 1e-12 * np.abs(gr).max()
        rel = float(np.max(np.abs(gr[m] - gs[m]) / np.abs(gr[m]))) if m.any() else 0.0
        gok = rel < 1e-8; ok &= gok
        print(f"  [{lab}] max rel err (src grad) = {rel:.2e} -> {'PASS' if gok else 'FAIL'}; "
              f"warm grad {tg*1e3:.1f} ms", flush=True)

    # --- GATE 3: small-signal (adjacent-fringe DIFFERENCE) ---
    print("\n=== GATE 3: small-signal adjacent-fringe difference ===", flush=True)
    ref = batched_scan(P, [tt])[0][0]; spl = sp.estep(tt, EV, AI)
    dref = np.abs(np.diff(ref, axis=1)); target = np.unravel_index(
        np.argmin(np.abs(dref - 1e-3)), dref.shape)
    p0, f0 = target
    d_ref = ref[p0, f0 + 1] - ref[p0, f0]; d_spl = spl[p0, f0 + 1] - spl[p0, f0]
    sok = abs(d_ref - d_spl) < 1e-6; ok &= sok
    print(f"  psr{p0} fringe {f0}->{f0+1}: dlnL_ref={d_ref:.6e} dlnL_split={d_spl:.6e} "
          f"|diff|={abs(d_ref-d_spl):.2e} -> {'PASS' if sok else 'FAIL'} (<1e-6)", flush=True)

    # --- Amendment 5: E-step profiling decomposition ---
    print("\n=== PROFILING: E-step split vs batch-only (~32 s ref) ===", flush=True)
    _ = sp.estep(tt, EV, AI)                                       # warm
    tb = jnp.asarray(tt); npsr = P["npsr"]
    typ = time.time()
    tabs = [sp._pulsar_ab_fn(p)(tb, jnp.asarray(EV[p])) for p in range(npsr)]
    _ = [(np.asarray(a), np.asarray(b)) for a, b in tabs]; t_yprod = time.time() - typ
    t_total = time.time()
    _ = sp.estep(tt, EV, AI); t_full = time.time() - t_total
    print(f"  full split E-step = {t_full:.3f} s (per-pulsar yprod+solve ~ {t_yprod:.3f} s; "
          f"rank-update+overhead ~ {t_full-t_yprod:.3f} s)", flush=True)
    print(f"  speedup vs 32 s batch-only reference = {32.0/t_full:.0f}x "
          f"(bounded by single-pulsar yprod recompute; NOT 1000x for the E-step)", flush=True)

    # --- GATE 4: candidate scorer values + cost ---
    print("\n=== GATE 4: candidate scorer values + cost ===", flush=True)
    src = tt[ndist:]; sp.build_tables(src, EV, AI)
    rng = np.random.default_rng(3); B = EV.shape[1]
    asg = rng.integers(0, B, size=(16, npsr))
    # reference: fisher_logL at each assignment (distances = EV[p, asg_p])
    thetas = np.tile(tt, (16, 1))
    for k in range(16):
        thetas[k, :ndist] = EV[np.arange(npsr), asg[k]]
    ref_sc = np.asarray(flb(jnp.asarray(thetas), data))
    spl_sc = sp.score_candidates(asg)
    dmax = float(np.nanmax(np.abs(ref_sc - spl_sc)))
    cok = dmax < 1e-8; ok &= cok
    print(f"  16 random assignments: max|dlnL| = {dmax:.2e} -> {'PASS' if cok else 'FAIL'} (<1e-8)", flush=True)
    for K in (1000, 100000):
        big = rng.integers(0, B, size=(K, npsr))
        _ = sp.score_candidates(big[:100])                       # warm
        t1 = time.time(); _ = sp.score_candidates(big); dt = time.time() - t1
        print(f"  K={K:>6d}: {dt*1e3:.1f} ms total -> {dt/K*1e6:.2f} us/candidate", flush=True)

    print(f"\n  [GATES 2/3/4] {'ALL PASS' if ok else 'FAIL'}", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else "recon"
    if mode == "recon":
        sys.exit(mode_recon())
    elif mode == "estep":
        sys.exit(mode_estep())
    elif mode == "post":
        sys.exit(mode_post())
    print(f"unknown mode {mode}"); sys.exit(2)
