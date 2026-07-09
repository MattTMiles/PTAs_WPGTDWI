"""Track B — deliverable R: the POSTERIOR REFERENDUM.

Track B settled COLD-START POINT ESTIMATION = NO-GO (F2 wall / L2b-c cusp). R asks the
question point-estimation cannot: WHERE DOES THE POSTERIOR MASS SIT? The needle (a unique
joint optimum at truth) is real but mas-narrow; its volume fraction in the source space is
~e^-(tens), while the plateau gains ~ln K_p per pulsar of FRINGE ENTROPY under
marginalisation. The contest is genuinely close and this module measures it.

The object (the honest fringe-marginalised likelihood; count-once STAR-TOPOLOGY marginal):

    lnL_marg(theta) = lnL_ref(theta)  +  sum_p  LSE_n[ (LNL_p(n) - lnL_ref) + lnprior_p(n) ]

  * theta = the 8*n_loud LOUD source params (faint sources FIXED at truth throughout R).
  * base distances are ALWAYS L0 (the EM prior mean = the true fringe k=0, since
    truth L_p == L0_p on this Asimov draw). lnL_ref(theta) = fisher_logL(theta, all L=L0).
  * LNL_p(n) = full-array fisher_logL with pulsar p at fringe n, others at L0 (the split
    E-step; gated 0.00 vs the batched_scan oracle).
  * lnprior_p(n) = -(n dL_p)^2 / (2 sig_p^2), sig from the lit priors (A2), = fringe_posterior's
    Gaussian tail. Summed over the +-3 sigma window fringes.
  * This is log( sum over the fringe lattice of exp[lnL + lnprior] ) UNDER THE ADDITIVE
    (star-topology) approximation lnL(theta,{L(n_p)}) ~ lnL_ref + sum_p (LNL_p(n_p)-lnL_ref).
    The additive-approx error is GATED (2-pulsar sublattice + MAP-config vs exact).
  * NOT the multi-counted `sum_p logsumexp_n[full lnL]` (that counts lnL_ref 116x and is not
    the log of any integral -> not a valid marginal likelihood; the linescan's -8000/0.003
    numbers were that multi-count). Count-once is the only valid evidence object.

Identity used everywhere: exp(-m_p(truth)) == the A2 P_true (census ceilings), and
lnL_marg(truth) = lnL(truth) - sum_p ln P_true_p.  m_p = LSE_n[...]  >= 0.

Env: jug-gpu, XLA autotune off, x64, jax persistent cache. Matt commits; never commit here.
Diagnostic use of TRUTH is licensed for this Asimov physics referendum and is LABELLED.
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
import trackB_split as TS
import discovery.matrix as dsm

CWT = "/home/mattm/projects/HSYMT/CW_transition"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]


# ============================================================
# batched split E-step: LNL[T, pulsar, fringe] at base distances L0
# ============================================================
def make_pulsar_ab_batch(sp, p):
    """(theta_bases (T,nth), Lvals (B,)) -> a_dep (T,B), b (T,B,ngp) for pulsar p.
    vmap over particles of the split's single-theta per-fringe evaluator."""
    Ft = sp.Ft[p]; Tt = sp.Tt[p]; s1d = sp.solve1d[p]; ys_p = sp.ys[p]
    cf_p = (sp.cf[0][p], sp.cf[1]); TtNmF_p = sp.TtNmF[p]; AI_p = int(sp.AI[p])

    @jax.jit
    def run(theta_bases, Lvals):
        def per_particle(tb):
            def one(L):
                params = sp._params(tb.at[AI_p].set(L))
                yp = ys_p(params); Nmy, _ = s1d(yp)
                return yp @ Nmy, Ft @ Nmy, Tt @ Nmy
            ytNmy, FtNmy, TtNmy = jax.vmap(one)(Lvals)          # (B,),(B,nc),(B,ngp)
            sol = dsm.matrix_solve(cf_p, FtNmy.T).T             # (B,nc)
            a_dep = -0.5 * (ytNmy - jnp.sum(FtNmy * sol, axis=1))
            b = TtNmy - sol @ TtNmF_p.T                         # (B,ngp)
            return a_dep, b
        return jax.vmap(per_particle)(theta_bases)              # (T,B),(T,B,ngp)
    return run


@jax.jit
def _rank_update(a_pf, b_pf, b_base_p, u_p, Minv_pp, Q_base, sum_a, a_base_p, const):
    """Per-pulsar LNL block (T,B) via the split rank-update, on GPU.
    LNL = const + (sum_a - a_base_p + a_pf) + 0.5*Q_pf, Q_pf via Delta-b quadratic form."""
    db = b_pf - b_base_p[:, None, :]                              # (T,B,ng)
    Q_pf = (Q_base[:, None]
            + 2.0 * jnp.einsum('tbi,ti->tb', db, u_p)
            + jnp.einsum('tbi,ij,tbj->tb', db, Minv_pp, db))
    return const + (sum_a[:, None] - a_base_p[:, None] + a_pf) + 0.5 * Q_pf


class MargEngine:
    """Holds the built problem + Split; produces LNL[T,p,n] and lnL_marg[T] for a batch of
    LOUD source-param vectors (faint held at truth), base distances L0. Truth-anchored."""

    def __init__(self, config="pop"):
        t0 = time.time()
        self.P = TE.build_problem(config, build_earth=False)
        P = self.P
        self.sp = TS.Split(P["amo_full"], P["data_obs"], P["theta_truth"])
        self.sp.AI = np.asarray(P["AI"]); self.sp._ensure_minv()
        self.ndist = P["n_dist"]; self.nth = P["amo_full"]["n_theta"]
        self.npsr = P["npsr"]; self.ngp = self.sp.ngp_gwb
        self.EV = P["EV"]; self.AI = np.asarray(P["AI"]); self.L0 = P["L0"]; self.dL = P["dL"]
        self.B = self.EV.shape[1]
        self.theta_truth = P["theta_truth"]
        self.src_truth = self.theta_truth[self.ndist:].copy()
        self.priors = P["priors"]
        self.names = P["names"]
        self.census_idx = np.array([self.names.index(nm) for nm in CENSUS])
        # loud/faint split
        pop = TE.CONFIGS[config].get("population", (P["ncw"],))
        self.n_loud = pop[0]; self.ncw = P["ncw"]
        self.n_loud_par = 8 * self.n_loud
        # ACTIVE loud dims for the referendum (others FIXED at truth). env R_ACTIVE:
        #  sky (default) = 2 sky/loud; skyfreq = +log10_fgw; full = all 8/loud.
        #  Rationale: extrinsic/freq Laplace factors CANCEL in f = Zn/(Zn+Zp); F2 shows SKY
        #  binds the registration. Fixing them at truth isolates the sky-registration contest
        #  and makes the plateau SMC tractable (the 24-D broad box has a narrow good-fit ridge
        #  the isotropic RW cannot follow -> acc collapse). Flagged favorable conditioning.
        mode = os.environ.get("R_ACTIVE", "sky")
        per = {"sky": [0, 1], "skyfreq": [0, 1, 4], "full": list(range(8))}[mode]
        self.r_active_mode = mode
        self.active = np.array([8 * k + j for k in range(self.n_loud) for j in per])
        # per-pulsar k index of each EV sample and the k=0 (base) column
        self.k_all = np.stack([np.round((self.EV[a] - self.L0[a]) / self.dL[a]).astype(int)
                               for a in range(self.npsr)])       # (npsr,B)
        self.base_col = np.array([int(np.argmin(np.abs(self.EV[a] - self.L0[a])))
                                  for a in range(self.npsr)])
        # frozen split arrays as np
        self._Minv = self.sp._Minv                               # (M,M)
        self._Minv_pp = self.sp._Minv_pp                         # (npsr,ng,ng)
        self._const = self.sp.a_const - 0.5 * (self.sp.ldP_cached + self.sp.logdet_cached)
        # per-pulsar batched ab fns (compiled lazily) + jitted vmapped base evaluator
        self._abb = {}
        self._ppab_jit = jax.jit(jax.vmap(self.sp.per_pulsar_ab))
        self._phi_lo, self._phi_hi = TE.phi_bounds(self.P)
        # per-pulsar unique-fringe reduction bookkeeping (offsets & window mask), lit prior
        self._fringe_tables()
        print(f"  [MargEngine] built config={config} in {time.time()-t0:.0f}s "
              f"(npsr={self.npsr}, B={self.B}, n_loud={self.n_loud})", flush=True)

    def _fringe_tables(self, prior_key="lit"):
        """Per pulsar: unique fringe k, offsets (kpc), window mask (|offs|<=3sig), and the
        reduceat index to collapse EV samples -> unique-fringe peak lnL. Prior-dependent
        (sigma); rebuildable for a different column."""
        sig = self.priors[prior_key]
        self.prior_key = prior_key
        self.uk = []; self.offs = []; self.redidx = []; self.order = []; self.lnprior = []
        for a in range(self.npsr):
            ks = self.k_all[a]
            order = np.argsort(ks); self.order.append(order)
            uk, idx0 = np.unique(ks[order], return_index=True)
            self.uk.append(uk); self.redidx.append(idx0)
            offs = uk * self.dL[a]; self.offs.append(offs)
            self.lnprior.append(-offs ** 2 / (2.0 * sig[a] ** 2))

    def estep_batch(self, src_T, chunk_psr=None):
        """src_T (T, n_cw_params) LOUD-varied source params (faint slots must be truth).
        Returns (LNL (T, npsr, B), lnref (T,)): LNL = full-array lnL with pulsar p at
        EV[p,n], others at L0; lnref = EXACT all-L0 lnL (= split.lnL(theta,L0)), assembled
        from a_base/b_base (NOT read off the EV grid -- dense-grid pulsars miss exact L0)."""
        src_T = np.atleast_2d(np.asarray(src_T, float))
        T = src_T.shape[0]
        src_T = np.clip(src_T, self._phi_lo, self._phi_hi)      # keep waveforms valid
        tb = np.zeros((T, self.nth)); tb[:, :self.ndist] = self.L0[None]; tb[:, self.ndist:] = src_T
        tb_j = jnp.asarray(tb)
        # base per-pulsar a,b (jitted vmap over particles)
        a_base, b_base = self._ppab_jit(tb_j)                   # (T,npsr),(T,npsr,ngp)
        a_base = np.asarray(a_base); b_base = np.asarray(b_base)
        bflat = b_base.reshape(T, self.npsr * self.ngp)
        u = bflat @ self._Minv                                   # Minv symmetric -> (T,M)
        Q_base = np.sum(bflat * u, axis=1)                       # (T,)
        u_p = u.reshape(T, self.npsr, self.ngp)
        sum_a = a_base.sum(axis=1)                               # (T,)
        lnref = self._const + sum_a + 0.5 * Q_base               # (T,) exact all-L0 lnL
        # GPU rank-update (host einsum over 116 pulsars was the bottleneck)
        b_base_j = jnp.asarray(b_base); u_p_j = jnp.asarray(u_p)
        Minv_pp_j = jnp.asarray(self._Minv_pp)
        Qb_j = jnp.asarray(Q_base); sa_j = jnp.asarray(sum_a); ab_j = jnp.asarray(a_base)
        LNL = np.empty((T, self.npsr, self.B))
        for p in range(self.npsr):
            if p not in self._abb:
                self._abb[p] = make_pulsar_ab_batch(self.sp, p)
            a_pf, b_pf = self._abb[p](tb_j, jnp.asarray(self.EV[p]))   # jnp (T,B),(T,B,ngp)
            LNL[:, p] = np.asarray(_rank_update(a_pf, b_pf, b_base_j[:, p], u_p_j[:, p],
                                                Minv_pp_j[p], Qb_j, sa_j, ab_j[:, p], self._const))
        return LNL, lnref

    # ---- reduce LNL[p,:] (one particle) -> per-fringe peak lnL over unique fringes ----
    def _fringe_peak(self, a, lnl_row):
        s = lnl_row[self.order[a]]
        return np.maximum.reduceat(s, self.redidx[a])            # (n_uniq_a,)

    def marg_from_LNL(self, LNL_T, lnref_T, return_parts=False):
        """(LNL (T,npsr,B), lnref (T,)) -> lnL_marg (T,) = lnref + sum_p m_p,
        m_p = LSE_n[(peak_p(n) - lnref) + lnprior_p(n)]. Vectorised over particles."""
        LNL_T = np.atleast_3d(LNL_T)
        T = LNL_T.shape[0]
        lnref_T = np.atleast_1d(lnref_T).astype(float)
        msum = np.zeros(T)
        m_all = np.empty((T, self.npsr)) if return_parts else None
        for a in range(self.npsr):
            M = LNL_T[:, a, :][:, self.order[a]]                 # (T,B) reordered by k
            peak = np.maximum.reduceat(M, self.redidx[a], axis=1)  # (T,n_uniq)
            logw = (peak - lnref_T[:, None]) + self.lnprior[a][None, :]
            mmax = logw.max(axis=1)
            m_a = mmax + np.log(np.sum(np.exp(logw - mmax[:, None]), axis=1))  # (T,)
            msum += m_a
            if return_parts:
                m_all[:, a] = m_a
        lnmarg = lnref_T + msum
        if return_parts:
            return lnmarg, m_all, np.exp(-m_all), lnref_T
        return lnmarg, lnref_T

    def lnL_marg(self, src_T, return_parts=False):
        LNL, lnref = self.estep_batch(src_T)
        return self.marg_from_LNL(LNL, lnref, return_parts=return_parts)

    # ---- scaled-coordinate front-end (loud dims only; faint fixed at truth) ----
    @property
    def scale_loud(self):
        return TE.phi_scale(self.P)[:self.n_loud_par]

    def src_from_x(self, X):
        """X (T, n_loud_par) scaled offsets of the LOUD params from truth -> src_T (faint=truth)."""
        X = np.atleast_2d(np.asarray(X, float)); T = X.shape[0]
        src = np.tile(self.src_truth, (T, 1))
        src[:, :self.n_loud_par] = self.src_truth[:self.n_loud_par] + X * self.scale_loud
        return src

    FIXED_CHUNK = 256    # every estep_batch call is padded to this T -> ONE compiled shape

    def lnmarg_x(self, X, chunk=None):
        """lnL_marg at truth + X*scale_loud. ALWAYS evaluated in blocks padded to
        FIXED_CHUNK so the per-pulsar jitted evaluators compile exactly one T-shape
        (avoids 116-fn recompiles for every distinct particle count)."""
        X = np.atleast_2d(np.asarray(X, float)); T = X.shape[0]
        C = self.FIXED_CHUNK
        out = np.empty(T)
        for c0 in range(0, T, C):
            blk = X[c0:c0 + C]; nb = blk.shape[0]
            if nb < C:
                blk = np.concatenate([blk, np.zeros((C - nb, blk.shape[1]))], 0)
            LNL, lnref = self.estep_batch(self.src_from_x(blk))
            lm, _ = self.marg_from_LNL(LNL, lnref)
            out[c0:c0 + nb] = lm[:nb]
        return np.nan_to_num(out, nan=-1e300, posinf=-1e300, neginf=-1e300)  # invalid waveform -> reject

    def reg_ext_idx(self):
        """LOUD param indices: registration (sky cos,gwphi + freq) vs extrinsic (inc,mc,h,phase,psi)."""
        reg = [8 * k + j for k in range(self.n_loud) for j in (0, 1, 4)]
        ext = [i for i in range(self.n_loud_par) if i not in reg]
        return np.array(reg), np.array(ext)


def _lse(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


# ============================================================
# GATES
# ============================================================
def mode_gate(E=None):
    print("=== R GATES: lnL_marg foundation ===", flush=True)
    E = E or MargEngine("pop")
    P = E.P; sp = E.sp; tt = E.theta_truth; ndist = E.ndist
    src_truth = E.src_truth

    # --- GATE A: split scalar recon (reuse the banked split gate on a few points) ---
    fl = P["amo_full"]["fisher_logL"]; data = P["data_obs"]
    rng = np.random.default_rng(0); scale = TE.phi_scale(P)
    plo, phi_hi = TE.phi_bounds(P)
    worst = 0.0
    for lab in ["truth"] + [f"rand{i}" for i in range(4)]:
        th = tt.copy()
        if lab != "truth":
            th[ndist:] = np.clip(th[ndist:] + rng.normal(size=len(tt) - ndist) * scale * 0.05, plo, phi_hi)
        a = float(fl(jnp.asarray(th), data)); b = float(sp.lnL(jnp.asarray(th)))
        worst = max(worst, abs(a - b))
    okA = worst < 1e-8
    print(f"  [GATE A] split.lnL == fisher_logL  max|diff|={worst:.2e} -> {'PASS' if okA else 'FAIL'}", flush=True)

    # --- GATE B: batched estep == single split.estep (per-pulsar) ---
    #   also == batched_scan oracle. Build 4 loud-varied particles.
    from trackB_p1_map import batched_scan, _local_sky_metric
    _, m_gp = _local_sky_metric(src_truth[0])
    srcs = [src_truth.copy()]
    for off in (0.2, 0.6, 1.2):
        s = src_truth.copy(); s[1] += np.radians(off) / m_gp; srcs.append(s)
    srcs = np.array(srcs)
    LNL_b, lnref_b = E.estep_batch(srcs)                         # (4,npsr,B),(4,)
    okB = True
    for t in range(len(srcs)):
        tb = np.concatenate([E.L0, srcs[t]])
        ref = sp.estep(tb, E.EV, E.AI)                           # single split
        orc, _ = batched_scan(P, [tb]); orc = orc[0]            # oracle
        d1 = float(np.max(np.abs(ref - LNL_b[t])))
        d2 = float(np.max(np.abs(orc - LNL_b[t])))
        dref = abs(lnref_b[t] - float(sp.lnL(jnp.asarray(tb))))
        okB &= (d1 < 1e-8 and d2 < 1e-8 and dref < 1e-6)
        print(f"  [GATE B] particle {t}: batched-vs-splitEstep {d1:.2e}, batched-vs-oracle {d2:.2e}, "
              f"lnref-vs-split.lnL {dref:.2e}", flush=True)
    print(f"  [GATE B] -> {'PASS' if okB else 'FAIL'} (<1e-8)", flush=True)

    # --- GATE C: exact lnL_ref == fisher_logL(all-L0) ---
    LNL0, lnref0 = E.estep_batch(src_truth[None])               # (1,npsr,B),(1,)
    LNL0 = LNL0[0]; lnref_truth = float(lnref0[0])
    lnL_truth = float(fl(jnp.asarray(tt), data))                 # truth dist == L0 on this draw
    okC = abs(lnref_truth - lnL_truth) < 1e-6
    grid_spread = float(np.max(np.abs(np.array([LNL0[a, E.base_col[a]] for a in range(E.npsr)]) - lnref_truth)))
    print(f"  [GATE C] exact lnref(truth)={lnref_truth:.4f} vs fisher_logL(truth)={lnL_truth:.4f} "
          f"|diff|={abs(lnref_truth-lnL_truth):.2e} -> {'PASS' if okC else 'FAIL'}", flush=True)
    print(f"        (grid base-col spread {grid_spread:.2e} = dense-grid sampling offset, NOT used as ref)", flush=True)

    # --- GATE D: exp(-m_p) == p_true_at for the CENSUS (fringe-centre, exact); array reported ---
    lnmarg, m_all, pt_all, ref_all = E.marg_from_LNL(LNL0[None], lnref0, return_parts=True)
    Pt_ref, _, _ = TE.p_true_at(P, LNL0, E.priors["lit"], beta=1.0)
    dcen = max(abs(pt_all[0, E.names.index(nm)] - Pt_ref[E.names.index(nm)]) for nm in CENSUS)
    darr = float(np.max(np.abs(pt_all[0] - Pt_ref)))
    okD = dcen < 1e-6
    print(f"  [GATE D] census exp(-m_p) vs p_true_at max|diff|={dcen:.2e} -> {'PASS' if okD else 'FAIL'}; "
          f"array max|diff|={darr:.2e} (dense-grid pulsars: global-ref vs per-pulsar-ref, expected)", flush=True)
    for nm in CENSUS:
        a = E.names.index(nm)
        print(f"        {nm}: P_true={pt_all[0,a]:.4f}  m_p={m_all[0,a]:.4f} "
              f"(ceiling {TE.CEILING[nm]})", flush=True)
    print(f"        lnL_marg(truth) = {lnmarg[0]:.4f} = lnL_ref {ref_all[0]:.2f} "
          f"+ sum_p m_p {m_all[0].sum():.2f}  (sum_p m_p = -sum ln P_true)", flush=True)
    print(f"        needle 'height' above the all-degenerate floor lives in these m_p; "
          f"census m_p small (registered), degenerate pulsars m_p ~ ln K_p (entropy).", flush=True)

    # --- GATE E: additive (star-topology) approximation error ---
    #   (E1) 2-pulsar joint sublattice: exact fisher_logL(both moved) vs lnL_ref + d1 + d2.
    #   (E2) MAP config (all pulsars at their truth-eval best fringe): exact profiled joint
    #        vs the additive sum -> the error at the config that dominates Z near truth.
    print("  [GATE E] additive/star-topology approximation error:", flush=True)
    fl_ = P["amo_full"]["fisher_logL"]
    # pick two strong registrars (census J1909, J1713) + test a few fringe offsets
    a1 = E.names.index("J1909-3744"); a2 = E.names.index("J1713+0747")
    LNLt = LNL0
    e1max = 0.0
    for k1 in (-2, -1, 1, 2):
        for k2 in (-2, -1, 1, 2):
            n1 = np.where(E.k_all[a1] == k1)[0]; n2 = np.where(E.k_all[a2] == k2)[0]
            if len(n1) == 0 or len(n2) == 0:
                continue
            c1 = n1[0]; c2 = n2[0]
            th = tt.copy(); th[E.AI[a1]] = E.EV[a1, c1]; th[E.AI[a2]] = E.EV[a2, c2]
            exact = float(fl_(jnp.asarray(th), data))
            add = lnref_truth + (LNLt[a1, c1] - lnref_truth) + (LNLt[a2, c2] - lnref_truth)
            e1max = max(e1max, abs(exact - add))
    print(f"        E1 2-pulsar sublattice max|exact - additive| = {e1max:.2e} nat "
          f"(cross-term b^T M^-1 b off-diagonal)", flush=True)
    # E2: MAP config at truth (should be all k=0 -> trivially exact); test at a plateau src
    s_pl = src_truth.copy(); s_pl[1] += np.radians(1.0) / m_gp
    LNL_pl_T, lnref_pl_T = E.estep_batch(s_pl[None])
    LNL_pl = LNL_pl_T[0]; lnref_pl = float(lnref_pl_T[0])
    post = TE.fringe_posterior(P, LNL_pl, None, E.priors["lit"], 1.0)
    th = np.concatenate([post["map_evalL"], s_pl])
    exact_map = float(fl_(jnp.asarray(th), data))
    add_map = lnref_pl + sum((LNL_pl[a, int(np.argmin(np.abs(E.EV[a] - post["map_evalL"][a])))] - lnref_pl)
                             for a in range(E.npsr))
    print(f"        E2 MAP-config @1deg-plateau src: exact profiled joint {exact_map:.2f}, "
          f"additive {add_map:.2f}, |diff| {abs(exact_map-add_map):.2f} nat "
          f"({E.npsr} pulsars re-registered)", flush=True)

    ok = okA and okB and okC and okD
    print(f"\n  [R GATES A-D] {'ALL PASS' if ok else 'FAIL'}; E is a reported approximation-error bound.", flush=True)
    return 0 if ok else 1


# ============================================================
# R1 -- Z_NEEDLE by local quadrature (truth-placed; TRUTH USE LABELLED)
#   Laplace baseline over the 24 loud dims + per-Hessian-eigenvector 1-D quadrature
#   (doubling-gated), separability-gated. Sharp (registration) directions get the
#   quadrature correction; soft (extrinsic) directions are near-Gaussian (reported).
# ============================================================
def _box_halfwidths_scaled(E):
    """R2 box half-widths in SCALED loud coords: sky +-2 deg, freq +-0.25 dex(full band),
    extrinsics broad physical ranges. Shared with R2 so the needle/plateau measure matches."""
    scale = E.scale_loud
    half = np.zeros(E.n_loud_par)
    tt = E.src_truth
    from trackB_p1_map import _local_sky_metric
    for k in range(E.n_loud):
        cg = tt[8 * k + 0]
        m_cg, m_gp = _local_sky_metric(cg)
        half[8 * k + 0] = (2.0 / m_cg)                 # dcos for 2 deg (scale=1)
        half[8 * k + 1] = (2.0 / m_gp) / scale[8 * k + 1]   # dgwphi for 2 deg
        half[8 * k + 4] = 0.25 / scale[8 * k + 4]      # +-0.25 dex log10_fgw (full band)
        # extrinsics: broad but kept strictly INSIDE phi_bounds (edges give degenerate/nan waveforms)
        half[8 * k + 2] = 0.9 / scale[8 * k + 2]       # cos_inc +-0.9 (avoid |cos_inc|=1)
        half[8 * k + 3] = 0.9 / scale[8 * k + 3]       # log10_mc +-0.9 (truth 9 -> [8.1,9.9] in [7,10])
        half[8 * k + 5] = 1.5 / scale[8 * k + 5]       # log10_h +-1.5 (truth -13.75 -> [-15.25,-12.25])
        half[8 * k + 6] = 0.9 * np.pi / scale[8 * k + 6]   # phase0 ~+-pi (periodic)
        half[8 * k + 7] = 0.45 * np.pi / scale[8 * k + 7]  # psi +-0.45pi (in [0,pi])
    mask = np.zeros(E.n_loud_par, bool); mask[E.active] = True   # inactive dims FIXED at truth
    half[~mask] = 0.0
    return half


def _lnVbox_active(E, boxhalf):
    """ln volume of the box over ACTIVE dims only (inactive are fixed, measure-1)."""
    return float(np.sum(np.log(2 * boxhalf[E.active])))


def _adaptive_step(E, i, g0, target=(-1.5, -0.3), h0=None, hcap=None):
    """Find scaled step h along dim i so g(+h e_i)-g0 in [target]. Bisect on log h.
    hcap bounds soft directions (flat -> would run to the box edge)."""
    reg, _ = E.reg_ext_idx()
    h = h0 if h0 is not None else (1e-6 if i in reg else 0.3)
    hcap = hcap if hcap is not None else 5.0
    for _ in range(40):
        h = min(h, hcap)
        x = np.zeros(E.n_loud_par); x[i] = h
        d = E.lnmarg_x(x[None])[0] - g0
        if target[0] <= d <= target[1]:
            return h, d
        if d > target[1]:            # too shallow -> bigger step
            if h >= hcap:            # capped and still shallow -> soft direction
                return h, d
            h *= 1.8
        else:                        # too deep -> smaller step
            h *= 0.5
    return h, d


def mode_r1(E=None):
    print("=== R1: Z_NEEDLE by local quadrature (TRUTH-PLACED; diagnostic truth use, labelled) ===", flush=True)
    E = E or MargEngine("pop")
    D = E.n_loud_par
    reg, ext = E.reg_ext_idx()
    g0 = float(E.lnmarg_x(np.zeros((1, D)))[0])
    lnL_truth = float(E.P["amo_full"]["fisher_logL"](jnp.asarray(E.theta_truth), E.P["data_obs"]))
    print(f"  lnL_marg(truth) g0 = {g0:.4f}  (= lnL(truth) {lnL_truth:.2f} + sum_p m_p {g0-lnL_truth:.2f})", flush=True)
    boxhalf = _box_halfwidths_scaled(E)

    HESS_CKPT = f"{CWT}/trackB_R_znaddle_hess.npz"
    if os.environ.get("R1_RESUME") and os.path.exists(HESS_CKPT):
        z = np.load(HESS_CKPT)
        hstep = z["hstep"]; Hess = z["Hess"]; Hpos = z["Hpos"]; evals = z["evals"]; evecs = z["evecs"]
        evals_c = z["evals_c"]; lnZ_laplace = float(z["lnZ_laplace"])
        print(f"  [RESUME] loaded Hessian/eigen from {HESS_CKPT}; Laplace ln Z_needle = {lnZ_laplace:.4f}", flush=True)
        return _r1_quadrature(E, D, evecs.shape[1], reg, ext, g0, lnL_truth, boxhalf, hstep, Hess, Hpos,
                              evals, evecs, evals_c, lnZ_laplace)

    A = E.active; nA = len(A)                             # active loud dims (others fixed at truth)
    print(f"  ACTIVE dims = {E.r_active_mode} ({nA} dims; inactive fixed at truth)", flush=True)

    # 1. BATCHED adaptive per-dim step (ACTIVE dims): geometric h-ladder, pick h with dg ~ -1
    t0 = time.time()
    HLAD = np.geomspace(3e-7, 2.0, 16)
    pts_a = []
    for i in A:
        hc = 0.5 * boxhalf[i]
        for h in HLAD:
            x = np.zeros(D); x[i] = min(h, hc); pts_a.append(x)
    dg_a = (E.lnmarg_x(np.array(pts_a)) - g0).reshape(nA, len(HLAD))
    hstep = np.zeros(D)
    for a_i, i in enumerate(A):
        hl = np.minimum(HLAD, 0.5 * boxhalf[i]); d = dg_a[a_i]
        good = np.where((d <= -0.2) & (d >= -3.0))[0]
        hstep[i] = hl[good[np.argmin(np.abs(d[good] + 1.0))]] if len(good) else (hl[-1] if d[-1] > -0.2 else hl[0])
    print(f"  adaptive steps ({time.time()-t0:.0f}s, batched): h~{np.median(hstep[A]):.2e}", flush=True)

    # 2. Hessian by FD over ACTIVE dims (nA x nA submatrix).
    pts = [np.zeros(D)]
    idx_ph = {}; idx_mh = {}
    for i in A:
        idx_ph[i] = len(pts); x = np.zeros(D); x[i] = hstep[i]; pts.append(x.copy())
        idx_mh[i] = len(pts); x = np.zeros(D); x[i] = -hstep[i]; pts.append(x.copy())
    idx_pp = {}
    for ii in range(nA):
        for jj in range(ii + 1, nA):
            i, j = A[ii], A[jj]
            idx_pp[(i, j)] = len(pts); x = np.zeros(D); x[i] = hstep[i]; x[j] = hstep[j]; pts.append(x.copy())
    gv = E.lnmarg_x(np.array(pts), chunk=256); g0b = gv[0]
    HessA = np.zeros((nA, nA))
    for ii, i in enumerate(A):
        HessA[ii, ii] = (gv[idx_ph[i]] + gv[idx_mh[i]] - 2 * g0b) / hstep[i] ** 2
    for ii in range(nA):
        for jj in range(ii + 1, nA):
            i, j = A[ii], A[jj]
            Hij = (gv[idx_pp[(i, j)]] - gv[idx_ph[i]] - gv[idx_ph[j]] + g0b) / (hstep[i] * hstep[j])
            HessA[ii, jj] = HessA[jj, ii] = Hij
    HposA = -0.5 * (HessA + HessA.T)
    evals, evecsA = np.linalg.eigh(HposA)                 # active-space eigen
    evecs = np.zeros((D, nA)); evecs[A, :] = evecsA       # embed in full loud space
    npos = int((evals > 0).sum())
    print(f"  Hessian: {npos}/{nA} positive curvature eigs; range [{evals.min():.3e}, {evals.max():.3e}]", flush=True)
    eps = 1e-6 * evals.max(); evals_c = np.clip(evals, eps, None)
    lnZ_laplace = g0 + 0.5 * (nA * np.log(2 * np.pi) - np.sum(np.log(evals_c)))
    print(f"  Laplace: ln Z_needle = {lnZ_laplace:.4f}  (0.5*sum ln(2pi/lambda_k), {nA} active dims)", flush=True)
    Hess = np.zeros((D, D))                               # embed for the checkpoint
    for ii, i in enumerate(A):
        for jj, j in enumerate(A):
            Hess[i, j] = HessA[ii, jj]
    Hpos = -0.5 * (Hess + Hess.T)
    np.savez(HESS_CKPT, hstep=hstep, Hess=Hess, Hpos=Hpos, evals=evals, evecs=evecs,
             evals_c=evals_c, lnZ_laplace=lnZ_laplace, active=A)
    return _r1_quadrature(E, D, nA, reg, ext, g0, lnL_truth, boxhalf, hstep, Hess, Hpos,
                          evals, evecs, evals_c, lnZ_laplace)


def _r1_quadrature(E, D, nA, reg, ext, g0, lnL_truth, boxhalf, hstep, Hess, Hpos,
                   evals, evecs, evals_c, lnZ_laplace):
    # 3. per-eigenvector 1-D quadrature, BATCHED across eigvecs; doubling gate n=49 vs 97.
    print("  per-eigenvector 1-D quadrature (batched; doubling gate n=49->97):", flush=True)
    reg_frac = np.array([float(np.sum(evecs[:, k][reg] ** 2)) for k in range(nA)])
    ranges = np.zeros(nA)
    for k in range(nA):
        sig = 1.0 / np.sqrt(evals_c[k])
        ranges[k] = min(6.0 * sig, np.min(boxhalf[E.active] / (np.abs(evecs[E.active, k]) + 1e-30)))

    def eigen_lnJ(nq):
        Xall = []; svs = []
        for k in range(nA):
            s = np.linspace(-ranges[k], ranges[k], nq); svs.append(s)
            Xall.append(s[:, None] * evecs[:, k][None, :])
        p = (E.lnmarg_x(np.concatenate(Xall)) - g0).reshape(nA, nq)
        def _tz(y, x):    # trapezoid (np.trapz removed in numpy 2.x)
            return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))
        return np.array([np.log(_tz(np.exp(p[k]), svs[k])) for k in range(nA)])

    lnJ49 = eigen_lnJ(49); lnJ = eigen_lnJ(97)
    dJmax = float(np.max(np.abs(lnJ - lnJ49)))
    print(f"  doubling-gate: max|lnJ(97)-lnJ(49)| over eigs = {dJmax:.3f} nat "
          f"(sum {abs(lnJ.sum()-lnJ49.sum()):.3f}); {'STABLE' if abs(lnJ.sum()-lnJ49.sum())<0.2 else 'REFINE'}", flush=True)
    rgauss = np.exp(lnJ) / (np.sqrt(2 * np.pi) / np.sqrt(evals_c))
    lnZ_quad = g0 + np.sum(lnJ)
    order = np.argsort(-reg_frac)
    print(f"  {'eig':>4s} {'lambda':>10s} {'sigma':>10s} {'reg_frac':>8s} {'r_gauss':>8s} {'range/sig':>9s}", flush=True)
    for k in order[:8]:
        print(f"  {k:4d} {evals_c[k]:10.2e} {1/np.sqrt(evals_c[k]):10.2e} {reg_frac[k]:8.3f} "
              f"{rgauss[k]:8.3f} {ranges[k]*np.sqrt(evals_c[k]):9.2f}", flush=True)
    print(f"  ln Z_needle (quadrature) = {lnZ_quad:.4f}  vs Laplace {lnZ_laplace:.4f} "
          f"(diff {lnZ_quad-lnZ_laplace:+.3f} nat)", flush=True)

    # 4. separability gate: top-2 sharpest eigvecs, 2-D vs product-of-1D
    ks = order[:2]
    v1, v2 = evecs[:, ks[0]], evecs[:, ks[1]]
    s1 = np.linspace(-ranges[ks[0]], ranges[ks[0]], 9)
    s2 = np.linspace(-ranges[ks[1]], ranges[ks[1]], 9)
    S1, S2 = np.meshgrid(s1, s2)
    X = S1.ravel()[:, None] * v1[None, :] + S2.ravel()[:, None] * v2[None, :]
    g2d = E.lnmarg_x(X, chunk=256).reshape(9, 9) - g0
    # 1-D profiles along each
    p1 = E.lnmarg_x(s1[:, None] * v1[None, :], chunk=64) - g0
    p2 = E.lnmarg_x(s2[:, None] * v2[None, :], chunk=64) - g0
    sep_resid = float(np.max(np.abs(g2d - (p1[None, :] + p2[:, None]))))
    print(f"  [separability] top-2 sharp eigvecs: max|g2d - (g1+g2)| = {sep_resid:.3f} nat "
          f"(<~0.1 -> product-of-1D valid)", flush=True)

    # 5. box + reg-tol connection (report sky/freq needle sigma in physical units)
    #    project the sharpest eigvec back to identify the binding physical scale
    kbind = order[0]; sig_bind = 1.0 / np.sqrt(evals_c[kbind])
    print(f"  needle binding sigma (sharpest eig, scaled) = {sig_bind:.3e}; "
          f"B0.2 cert tol ~1e-4 scaled -> box '+-3 reg-tol' ~ +-3e-4 scaled", flush=True)

    np.savez(f"{CWT}/trackB_R_znaddle.npz",
             g0=g0, lnL_truth=lnL_truth, Hess=Hess, Hpos=Hpos, evals=evals, evecs=evecs,
             hstep=hstep, lnJ=lnJ, rgauss=rgauss, reg_frac=reg_frac, ranges=ranges,
             boxhalf=boxhalf, lnZ_laplace=lnZ_laplace, lnZ_quad=lnZ_quad,
             sep_resid=sep_resid, reg_idx=reg, ext_idx=ext)
    print(f"\n  [R1] ln Z_needle = {lnZ_quad:.4f} (quadrature) / {lnZ_laplace:.4f} (Laplace); "
          f"saved trackB_R_znaddle.npz", flush=True)
    return 0


# ============================================================
# R2 -- Z_PLATEAU by tempered SMC over the credible box (vectorised through the split)
#   Box: truth-centred +-2 deg sky per loud, +-0.25 dex freq, broad extrinsics; faint FIXED
#   at truth. The sampler need NOT find the mas needle (Z_needle is R1); if it stumbles into
#   the needle box those particles are excised. Evidence + ESS + multi-seed error reported.
# ============================================================
def _systematic_resample(w, rng):
    N = len(w); pos = (rng.random() + np.arange(N)) / N
    return np.searchsorted(np.cumsum(w), pos)


def smc_plateau(E, boxhalf, N=2000, seed=0, target_ess=0.5, n_mcmc=5, needle_sig=None,
                verbose=True):
    """Tempered SMC. Returns logZ_density = log( (1/V_box) int exp(lnL_marg) dx ),
    plus diagnostics. logZ_plateau_unnorm = logZ_density + ln V_box."""
    rng = np.random.default_rng(seed)
    D = len(boxhalf)
    X = rng.uniform(-boxhalf, boxhalf, size=(N, D))
    L = E.lnmarg_x(X, chunk=256)
    beta = 0.0; logZ = 0.0; s = 0.5; stages = 0
    ess_hist = []; beta_hist = []
    excised = 0
    while beta < 1.0 - 1e-12:
        def ess_at(db):
            lw = db * L; lw = lw - lw.max()
            w = np.exp(lw); return (w.sum() ** 2) / np.sum(w ** 2)
        hi = 1.0 - beta
        if ess_at(hi) >= target_ess * N:
            db = hi
        else:
            lo = 0.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if ess_at(mid) >= target_ess * N: lo = mid
                else: hi = mid
            db = lo
        lw = db * L
        logZ += _lse(lw) - np.log(N)
        w = np.exp(lw - _lse(lw)); ess = 1.0 / np.sum(w ** 2)
        idx = _systematic_resample(w, rng); X = X[idx]; L = L[idx]
        beta += db; stages += 1; ess_hist.append(ess); beta_hist.append(beta)
        # covariance-preconditioned RW on ACTIVE dims (inactive frozen at truth)
        act = E.active
        Cact = np.cov(X[:, act].T)
        Cact = np.atleast_2d(Cact) + 1e-9 * np.eye(len(act)) * np.trace(np.atleast_2d(Cact)) / len(act)
        try:
            Lchol = np.linalg.cholesky(Cact)
        except np.linalg.LinAlgError:
            Lchol = np.diag(np.sqrt(np.maximum(np.diag(Cact), 1e-30)))
        acc = 0.0
        for _ in range(n_mcmc):
            Xp = X.copy()
            z = rng.normal(size=(N, len(act)))
            Xp[:, act] = X[:, act] + s * (z @ Lchol.T)
            inbox = np.all(np.abs(Xp) <= boxhalf + 1e-30, axis=1)
            Lp = np.where(inbox, E.lnmarg_x(Xp, chunk=256), -np.inf)
            a = beta * (Lp - L)
            accept = inbox & (np.log(rng.random(N)) < a)
            X[accept] = Xp[accept]; L[accept] = Lp[accept]
            acc += accept.mean()
        acc /= n_mcmc
        s *= float(np.exp(acc - 0.234)); s = min(max(s, 1e-6), 1.0)
        if verbose:
            print(f"    [smc s{seed}] stage {stages}: beta={beta:.4f} dbeta={db:.4f} "
                  f"ESS={ess:.0f} acc={acc:.2f} s={s:.3f} Lmax={L.max():.1f} logZ={logZ:.3f}", flush=True)
        np.savez(f"{CWT}/trackB_R_zplateau_s{seed}_ckpt.npz",   # recoverable per-stage progress
                 beta=beta, logZ_density=logZ, lnVbox=_lnVbox_active(E, boxhalf),
                 stages=stages, ess_hist=np.array(ess_hist), beta_hist=np.array(beta_hist),
                 X=X, L=L, seed=seed)
    # excise needle-box particles (active dims within few needle-sigma of truth)
    if needle_sig is not None:
        near = np.all(np.abs(X[:, E.active]) < 5.0 * needle_sig, axis=1)
        excised = int(near.sum())
    lnVbox = _lnVbox_active(E, boxhalf)
    return dict(logZ_density=logZ, lnVbox=lnVbox, logZ_plateau=logZ + lnVbox,
                ess_hist=np.array(ess_hist), beta_hist=np.array(beta_hist),
                stages=stages, Lmax=float(L.max()), Lmed=float(np.median(L)),
                X_final=X, L_final=L, excised=excised, N=N, seed=seed)


def mode_r2(E=None):
    print("=== R2: Z_PLATEAU by tempered SMC (truth-centred +-2deg box; needle excised) ===", flush=True)
    E = E or MargEngine("pop")
    boxhalf = _box_halfwidths_scaled(E)
    D = E.n_loud_par
    g0 = float(E.lnmarg_x(np.zeros((1, D)))[0])
    # needle sigma (from R1 if present) for excision
    needle_sig = None
    try:
        z = np.load(f"{CWT}/trackB_R_znaddle.npz")
        needle_sig = float(1.0 / np.sqrt(np.clip(z["evals"], 1e-6 * z["evals"].max(), None).min()))
    except Exception:
        needle_sig = 1e-4

    # --- R2 GATE: batched lnL_marg == single-particle on 8 random box points (1e-10) ---
    rng = np.random.default_rng(7)
    Xg = rng.uniform(-boxhalf, boxhalf, size=(8, D))
    Lbatch = E.lnmarg_x(Xg, chunk=256)
    Lsingle = np.array([E.lnmarg_x(Xg[i][None])[0] for i in range(8)])
    fin = np.isfinite(Lbatch) & np.isfinite(Lsingle) & (Lbatch > -1e299)
    dg = float(np.max(np.abs(Lbatch[fin] - Lsingle[fin]))) if fin.any() else np.inf
    print(f"  [R2 GATE] batched vs single lnL_marg (finite pts {int(fin.sum())}/8) max|diff| = {dg:.2e} -> "
          f"{'PASS' if dg < 1e-9 else 'FAIL'} (<1e-9)", flush=True)
    print(f"  active={E.r_active_mode}; box (scaled) half-widths: sky~{boxhalf[0]:.3f}/{boxhalf[1]:.3f} "
          f"freq~{boxhalf[4]:.3f} ext~{boxhalf[5]:.2f}; ln V_box(active) = {_lnVbox_active(E, boxhalf):.3f}", flush=True)
    print(f"  lnL_marg(truth) = {g0:.2f} (needle peak, for reference)", flush=True)

    # --- multi-seed SMC for the evidence error ---
    N = int(os.environ.get("R2_N", "2000")); nseed = int(os.environ.get("R2_SEEDS", "4"))
    nmc = int(os.environ.get("R2_MCMC", "5"))
    runs = []
    for sd in range(nseed):
        t0 = time.time()
        r = smc_plateau(E, boxhalf, N=N, seed=sd, n_mcmc=nmc, needle_sig=needle_sig, verbose=True)
        r["wall"] = time.time() - t0
        runs.append(r)
        print(f"  [smc s{sd}] DONE: logZ_plateau={r['logZ_plateau']:.3f} "
              f"(density {r['logZ_density']:.3f} + lnV {r['lnVbox']:.3f}); stages={r['stages']} "
              f"Lmax={r['Lmax']:.1f} Lmed={r['Lmed']:.1f} excised={r['excised']} {r['wall']:.0f}s", flush=True)
        np.savez(f"{CWT}/trackB_R_zplateau_s{sd}.npz",
                 logZ_plateau=r["logZ_plateau"], logZ_density=r["logZ_density"], lnVbox=r["lnVbox"],
                 ess_hist=r["ess_hist"], beta_hist=r["beta_hist"], Lmax=r["Lmax"], Lmed=r["Lmed"],
                 X_final=r["X_final"], L_final=r["L_final"], excised=r["excised"], boxhalf=boxhalf,
                 N=N, seed=sd, g0=g0)
    lz = np.array([r["logZ_plateau"] for r in runs])
    print(f"\n  [R2] logZ_plateau = {lz.mean():.3f} +- {lz.std(ddof=1) if len(lz)>1 else 0:.3f} "
          f"(over {nseed} seeds); min ESS across stages "
          f"{min(r['ess_hist'].min() for r in runs):.0f}/{N}", flush=True)
    # best plateau particle (for R3 MAP census P_true)
    ibest = int(np.argmax([r["Lmax"] for r in runs]))
    Xbest = runs[ibest]["X_final"][np.argmax(runs[ibest]["L_final"])]
    np.savez(f"{CWT}/trackB_R_zplateau_summary.npz",
             logZ_plateau_mean=lz.mean(), logZ_plateau_std=lz.std(ddof=1) if len(lz) > 1 else 0.0,
             logZ_all=lz, lnVbox=runs[0]["lnVbox"], boxhalf=boxhalf, g0=g0,
             Xbest=Xbest, Lbest=runs[ibest]["Lmax"], nseed=nseed, N=N)
    print(f"  saved trackB_R_zplateau_summary.npz + per-seed npz", flush=True)
    return 0


# ============================================================
# R3 -- the referendum number f + break-even volume + census P_true at the peak
# ============================================================
def mode_r3(E=None):
    print("=== R3: the referendum number f = Z_needle/(Z_needle+Z_plateau) ===", flush=True)
    zn = np.load(f"{CWT}/trackB_R_znaddle.npz")
    zp = np.load(f"{CWT}/trackB_R_zplateau_summary.npz")
    lnZn_q = float(zn["lnZ_quad"]); lnZn_L = float(zn["lnZ_laplace"])
    lnZp = float(zp["logZ_plateau_mean"]); lnZp_e = float(zp["logZ_plateau_std"])
    E_needle_used = lnZn_q
    # f (logistic of lnZn - lnZp)
    def frac(lnZn):
        d = lnZn - lnZp
        return 1.0 / (1.0 + np.exp(-d))
    f_q = frac(lnZn_q); f_L = frac(lnZn_L)
    print(f"  ln Z_needle  = {lnZn_q:.3f} (quadrature) / {lnZn_L:.3f} (Laplace)", flush=True)
    print(f"  ln Z_plateau = {lnZp:.3f} +- {lnZp_e:.3f}", flush=True)
    print(f"  ln Z_needle - ln Z_plateau = {lnZn_q - lnZp:+.3f} nat", flush=True)
    print(f"  >>> f = {f_q:.4f} (quadrature Z_needle) / {f_L:.4f} (Laplace)", flush=True)
    # error on f from the plateau evidence error
    f_hi = frac(lnZn_q - (lnZp_e)); f_lo = frac(lnZn_q + (lnZp_e))
    print(f"      f 1-sigma band (from plateau logZ error): [{min(f_lo,f_hi):.4f}, {max(f_lo,f_hi):.4f}]", flush=True)

    # break-even SKY radius: Zp ∝ (theta_sky/2deg)^(2 n_loud); set Zp(theta*) = Zn
    n_loud = 3
    ratio = lnZn_q - lnZp                                   # ln(Zn/Zp) at 2 deg
    theta_star = 2.0 * np.exp(ratio / (2 * n_loud))         # deg per loud source
    print(f"  BREAK-EVEN sky prior box: theta* = {theta_star:.4g} deg per loud source "
          f"(Zn=Zp when the sky prior is +-{theta_star:.3g} deg); current box +-2 deg", flush=True)
    # break-even total volume factor
    print(f"  BREAK-EVEN volume factor: V*/V_2deg = exp(ln Zn - ln Zp) = {np.exp(ratio):.3e}", flush=True)

    # census P_true at the posterior PEAK
    print("  census-3 P(true fringe) at the posterior peak:", flush=True)
    if f_q > 0.5:
        print("    peak = NEEDLE (truth): P_true = census ceilings 0.953/0.989/0.998 (by construction).", flush=True)
    else:
        print("    peak = PLATEAU MAP; evaluating p_true_at the best plateau particle...", flush=True)
        E = E or MargEngine("pop")
        Xbest = zp["Xbest"]
        src = E.src_from_x(Xbest[None])
        LNL, lnref = E.estep_batch(src)
        Pt, _, _ = TE.p_true_at(E.P, LNL[0], E.priors["lit"], beta=1.0)
        for nm in CENSUS:
            print(f"      {nm}: P_true(plateau MAP) = {Pt[E.names.index(nm)]:.4f}", flush=True)

    # verdict (pre-registered R4 language)
    print("\n  === R4 VERDICT ===", flush=True)
    f = f_q
    if f >= 0.95:
        print(f"  (i) f = {f:.3f} >= 0.95: THE DATA IDENTIFY TRUTH. Track B's failure is SEARCH-ONLY;", flush=True)
        print(f"      the fringe-marginalised posterior concentrates at the needle -> guided/anchored", flush=True)
        print(f"      inference (prong 3) is the method; the terminal verdict softens, prong-3 unparks.", flush=True)
    elif f <= 0.05:
        print(f"  (ii) f = {f:.3f} <= 0.05: NO-GO DEEPENS TO INFORMATION-THEORETIC. The honestly", flush=True)
        print(f"       marginalised data do NOT prefer truth; only more data helps. The design", flush=True)
        print(f"       theorem + E-track are the whole story.", flush=True)
    else:
        print(f"  (iii) f = {f:.3f} INTERMEDIATE: report f and the break-even sky radius {theta_star:.3g} deg;", flush=True)
        print(f"        the verdict quotes both; B1 (noise) becomes decisive.", flush=True)
    np.savez(f"{CWT}/trackB_R_referendum_result.npz",
             lnZn_quad=lnZn_q, lnZn_laplace=lnZn_L, lnZp=lnZp, lnZp_err=lnZp_e,
             f_quad=f_q, f_laplace=f_L, theta_star_deg=theta_star, ratio=ratio)
    print(f"\n  saved trackB_R_referendum_result.npz", flush=True)
    return 0


def mode_pipe():
    """Build the engine ONCE, then gate -> R1 -> R2 -> R3 (avoids 4 rebuilds/compiles)."""
    print("=== R PIPELINE: gate -> R1 -> R2 -> R3 (one engine build) ===", flush=True)
    E = MargEngine("pop")
    # throughput probe (sizes R2)
    t0 = time.time(); _ = E.lnmarg_x(np.zeros((256, E.n_loud_par))); dt = time.time() - t0
    print(f"  [throughput] estep_batch(256)+marg = {dt:.2f}s ({dt/256*1e3:.2f} ms/particle)", flush=True)
    rc = mode_gate(E)
    if rc != 0:
        print("  GATE FAILED -> abort pipeline.", flush=True); return rc
    mode_r1(E)
    mode_r2(E)
    mode_r3(E)
    return 0


POSTMORTEM_PSRS = ["J0437-4715", "J2222-0137", "J1713+0747", "J1909-3744", "J0711-6830"]


def mode_postmortem(E=None):
    """Bounded CPU-side postmortem on the BANKED R2 plateau samples (no new sampling): per
    pulsar, (a) distinct best-fringe indices the plateau uses vs census K, (b) marginal
    P(true fringe) under the full plateau+needle posterior (needle weighted by f), (c)
    effective wrong-solution count exp(entropy of the fringe posterior). Q: does any single
    pulsar's EM prior+lever arm near-decide its OWN fringe while the joint sky is not?"""
    print("=== R POSTMORTEM: per-pulsar fringe structure on the banked plateau samples ===", flush=True)
    E = E or MargEngine("pop")
    res = np.load(f"{CWT}/trackB_R_referendum_result.npz")
    f = float(res["f_quad"])
    # load plateau particles (both seeds) -> source-sky offsets X (scaled, active=sky)
    Xs = []
    for sd in (0, 1):
        z = np.load(f"{CWT}/trackB_R_zplateau_s{sd}.npz")
        Xs.append(z["X_final"])
    X = np.concatenate(Xs, axis=0); Npart = X.shape[0]
    print(f"  loaded {Npart} plateau particles ({len(Xs)} seeds); f={f:.3e}", flush=True)
    # E-step at the plateau particles + at truth (needle)
    LNL_pl, _ = _estep_chunked(E, X)
    LNL_nd, _ = E.estep_batch(E.src_from_x(np.zeros((1, E.n_loud_par))))   # truth = needle
    idx = {nm: E.names.index(nm) for nm in POSTMORTEM_PSRS}

    print(f"\n  {'pulsar':12s} {'K_cens':>7s} {'dL_pc':>7s} {'sigEM_pc':>8s} "
          f"{'#distinctMAP':>12s} {'MAP=0frac':>9s} {'Ptrue_pl':>9s} {'Ptrue_full':>10s} {'expH':>7s}",
          flush=True)
    table = {}
    for nm in POSTMORTEM_PSRS:
        a = idx[nm]
        ord_a, red_a, uk_a, lp_a = E.order[a], E.redidx[a], E.uk[a], E.lnprior[a]
        sig = E.priors["lit"][a]; dL = E.dL[a]
        Kc = 2 * int(np.floor(3.0 * sig / dL + 1e-9))
        # plateau: per-particle fringe posterior w (beta=1), MAP fringe, q(true)
        M = LNL_pl[:, a, :][:, ord_a]                          # (Npart,B)
        lnL_u = np.maximum.reduceat(M, red_a, axis=1)          # (Npart,nuniq)
        logw = (lnL_u - lnL_u.max(axis=1, keepdims=True)) + lp_a[None, :]
        w = np.exp(logw - logw.max(axis=1, keepdims=True)); w /= w.sum(axis=1, keepdims=True)
        map_k = uk_a[np.argmax(w, axis=1)]                     # (Npart,)
        true_col = np.where(uk_a == 0)[0]
        q_true_pl = w[:, true_col[0]] if len(true_col) else np.zeros(Npart)
        # needle fringe posterior (at truth)
        Mn = LNL_nd[0, a, :][ord_a]; lnLun = np.maximum.reduceat(Mn, red_a)
        logwn = (lnLun - lnLun.max()) + lp_a; wn = np.exp(logwn - logwn.max()); wn /= wn.sum()
        q_true_nd = float(wn[true_col[0]]) if len(true_col) else 0.0
        # aggregates
        distinct_map = np.unique(map_k); n_distinct = len(distinct_map)
        map0_frac = float(np.mean(map_k == 0))
        Ptrue_pl = float(np.mean(q_true_pl))                   # plateau-averaged P(true)
        Ptrue_full = f * q_true_nd + (1 - f) * Ptrue_pl        # full posterior (needle weighted by f)
        qbar = w.mean(axis=0)                                  # marginal (plateau) fringe posterior
        qbar = qbar / qbar.sum()
        H = -np.sum(qbar[qbar > 0] * np.log(qbar[qbar > 0]))
        expH = float(np.exp(H))
        print(f"  {nm:12s} {Kc:7d} {dL*KPC:7.3f} {sig*KPC:8.2f} {n_distinct:12d} "
              f"{map0_frac:9.3f} {Ptrue_pl:9.4f} {Ptrue_full:10.4f} {expH:7.2f}", flush=True)
        # top-3 populated fringes
        top = np.argsort(-qbar)[:4]
        tops = ", ".join(f"k={uk_a[t]}:{qbar[t]:.2f}" for t in top)
        print(f"    {'':12s} top fringes (marginal): {tops}; census K={Kc}, needle P(true)={q_true_nd:.4f}", flush=True)
        table[nm] = dict(Kc=Kc, n_distinct=n_distinct, map0_frac=map0_frac, Ptrue_pl=Ptrue_pl,
                         Ptrue_full=Ptrue_full, expH=expH, q_true_nd=q_true_nd,
                         dL_pc=dL * KPC, sig_pc=sig * KPC)
    np.savez(f"{CWT}/trackB_R_postmortem.npz", f=f, Npart=Npart,
             **{f"{nm}_{k}": table[nm][k] for nm in table for k in table[nm]})
    print("\n  Q: near-decided own fringe on the plateau <=> expH~1 AND MAP concentrated on ONE fringe.", flush=True)
    print("  saved trackB_R_postmortem.npz", flush=True)
    return 0


def _estep_chunked(E, X, C=256):
    """E-step at scaled offsets X (Npart,nloud_par), chunked; returns (LNL (Npart,npsr,B), None)."""
    Npart = X.shape[0]
    LNL = np.empty((Npart, E.npsr, E.B))
    for c0 in range(0, Npart, C):
        blk = X[c0:c0 + C]; nb = blk.shape[0]
        if nb < C:
            blk = np.concatenate([blk, np.zeros((C - nb, blk.shape[1]))], 0)
        l, _ = E.estep_batch(E.src_from_x(blk))
        LNL[c0:c0 + nb] = l[:nb]
    return LNL, None


KPC = 1000.0  # kpc -> pc for the report


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else "gate"
    if mode == "postmortem":
        sys.exit(mode_postmortem())
    if mode == "gate":
        sys.exit(mode_gate())
    elif mode == "r1":
        sys.exit(mode_r1())
    elif mode == "r2":
        sys.exit(mode_r2())
    elif mode == "r3":
        sys.exit(mode_r3())
    elif mode == "pipe":
        sys.exit(mode_pipe())
    print(f"unknown mode {mode}"); sys.exit(2)
