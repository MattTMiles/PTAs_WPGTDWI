"""Track B — B1: the noisy, TARGETED (counterpart-conditioned) pipeline.

Objects
-------
`B1Marg` -- the honest count-once STAR-TOPOLOGY fringe-marginalised likelihood of deliverable
R, rebuilt on `trackB_b1_core.B1Split` so that (data, pulsar-term mask) are RUNTIME args:

    lnL_marg(theta) = lnL_ref(theta) + sum_p LSE_n[ (LNL_p(n) - lnL_ref) + lnprior_p(n) ]

  * lnL_ref(theta) = full-array lnL at the BASE distances (L0 = the EM prior means).
  * LNL_p(n)       = full-array lnL with pulsar p at fringe centre n, others at base.
  * lnprior_p(n)   = -(L_p(n) - L0_p)^2 / (2 sigma_EM,p^2)   (literature priors, A2).
  * This is log(sum over the fringe lattice) under the additive/star approximation, with
    lnL_ref counted ONCE. It is the SAME object R used (gated there to 0.14 nat), so B1's
    numbers are directly comparable to R's f = 6.9e-7.

The distinction that matters for B1: R integrated over the SKY registration dims and held
frequency + chirp mass at truth. The targeted scenario REMOVES the sky (an EM counterpart
supplies it to <1 arcsec) and leaves whatever registration dims remain. Which dims those are
is measured by `trackB_b1_pilot.py`, not assumed.

Modes
-----
  needle      1-D / 2-D profile of lnL_marg about truth in the registration axes: the needle
              width, the plateau ceiling, and the cost model for everything downstream.
  (further modes are added once the pilot has fixed the registration-axis set)

Env: jug-gpu, autotune off, x64. Matt commits; never commit from here.
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
import trackB_b1_core as C
import trackB_estimator as TE
import discovery.matrix as dsm

CWT = "/home/mattm/projects/HSYMT/CW_transition"
I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)
AXIS_NAME = ["cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h", "phase0", "psi"]


# ============================================================
# batched split E-step over a stack of source params (data + mask fixed per call)
# ============================================================
def make_pulsar_ab_batch(sp, p):
    """(thetas (T,nth), Lvals (B,), data, pmask) -> a (T,B), b (T,B,ngp) for pulsar p."""
    Ft = sp.Ft[p]; Tt = sp.Tt[p]; s1d = sp.solve1d[p]; ys_p = sp.ys[p]
    cf_p = (sp.cf[0][p], sp.cf[1]); TtNmF_p = sp.TtNmF[p]; AI_p = int(sp.AI[p])

    @jax.jit
    def run(thetas, Lvals, data_tuple, pmask):
        def per_particle(tb):
            def one(L):
                params = sp._params(tb.at[AI_p].set(L), data_tuple, pmask)
                yp = ys_p(params); Nmy, _ = s1d(yp)
                return yp @ Nmy, Ft @ Nmy, Tt @ Nmy
            ytNmy, FtNmy, TtNmy = jax.vmap(one)(Lvals)
            sol = dsm.matrix_solve(cf_p, FtNmy.T).T
            a = -0.5 * (ytNmy - jnp.sum(FtNmy * sol, axis=1))
            b = TtNmy - sol @ TtNmF_p.T
            return a, b
        return jax.vmap(per_particle)(thetas)
    return run


@jax.jit
def _rank_update(a_pf, b_pf, b_base_p, u_p, Minv_pp, Q_base, sum_a, a_base_p, const):
    db = b_pf - b_base_p[:, None, :]
    Q_pf = (Q_base[:, None]
            + 2.0 * jnp.einsum('tbi,ti->tb', db, u_p)
            + jnp.einsum('tbi,ij,tbj->tb', db, Minv_pp, db))
    return const + (sum_a[:, None] - a_base_p[:, None] + a_pf) + 0.5 * Q_pf


class B1Marg:
    """lnL_marg for a FIXED (data, mask, base distances, EV grid); batched over source params."""

    T_CHUNK = 32          # particles per compiled call (memory: T x B x n_src x n_toa)

    def __init__(self, P, data, mask, base_L=None, EV=None, dL=None, prior_key="lit"):
        self.P = P; self.sp = P["sp"]
        self.data = data
        self.mask = jnp.asarray(mask)
        # mask == 1 everywhere -> use the pterm-only residual (no wasted Earth-term evaluation
        # per fringe). Gated bit-exactly by `gate_fast_full` before any science uses it.
        self.all_on = bool(np.all(np.asarray(mask) == 1.0))
        self.sp.enable_fast_full(self.all_on)
        self.npsr = P["npsr"]; self.ngp = self.sp.ngp_gwb
        self.ndist = P["n_dist"]; self.nth = P["n_theta"]
        self.AI = np.asarray(P["AI"])
        self.base_L = np.asarray(P["L0"] if base_L is None else base_L, float)
        self.dL = np.asarray(P["dL_truth"] if dL is None else dL, float)
        self.EV = np.asarray(P["EV_truth"] if EV is None else EV, float)
        self.B = self.EV.shape[1]
        self.sp.AI = self.AI; self.sp._ensure_minv()
        self._Minv = self.sp._Minv; self._Minv_pp = self.sp._Minv_pp
        self._const = self.sp.a_const - 0.5 * (self.sp.ldP_cached + self.sp.logdet_cached)
        self._abb = {}
        self._ppab_b = jax.jit(jax.vmap(self.sp._per_pulsar_ab_impl, in_axes=(0, None, None)))
        self.FT = C.FringeTables(P, self.EV, self.dL, prior_key)
        # per-pulsar reduce tables as arrays (prior anchored at L0, NOT at base_L)
        self.order = self.FT.order; self.redidx = self.FT.redidx
        self.lnprior = self.FT.lnprior; self.uk = self.FT.uk
        self.n_calls = 0

    def set_smask(self, w):
        """Forward a per-source feed weight to the split (FORGE-G soft source side) and clear the
        batched evaluators, which close over sp.smask through sp._params at trace time. A non-None
        smask forces the masked residual path; the all-on fast-full optimisation is unavailable
        while a smask is live (MultiSourceDelay does not honour SMASK)."""
        if w is not None and self.all_on:
            self.sp.enable_fast_full(False)
            self.all_on = False
        self.sp.set_smask(w)
        self._abb = {}
        self._ppab_b = jax.jit(jax.vmap(self.sp._per_pulsar_ab_impl, in_axes=(0, None, None)))

    # ---- theta packing ----
    def theta_of(self, src_T):
        src_T = np.atleast_2d(np.asarray(src_T, float))
        T = src_T.shape[0]
        th = np.zeros((T, self.nth))
        th[:, :self.ndist] = self.base_L[None, :]
        th[:, self.ndist:] = src_T
        return th

    def estep_batch(self, thetas):
        """thetas (T,nth) -> LNL (T,npsr,B), lnref (T,)."""
        thetas = np.atleast_2d(np.asarray(thetas, float))
        T = thetas.shape[0]
        tb = jnp.asarray(thetas)
        a_base, b_base = self._ppab_b(tb, self.data, self.mask)
        a_base = np.asarray(a_base); b_base = np.asarray(b_base)
        bflat = b_base.reshape(T, self.npsr * self.ngp)
        u = bflat @ self._Minv
        Q_base = np.sum(bflat * u, axis=1)
        u_p = u.reshape(T, self.npsr, self.ngp)
        sum_a = a_base.sum(axis=1)
        lnref = self._const + sum_a + 0.5 * Q_base
        b_base_j = jnp.asarray(b_base); u_p_j = jnp.asarray(u_p)
        Minv_pp_j = jnp.asarray(self._Minv_pp)
        Qb = jnp.asarray(Q_base); sa = jnp.asarray(sum_a); ab = jnp.asarray(a_base)
        LNL = np.empty((T, self.npsr, self.B))
        for p in range(self.npsr):
            if p not in self._abb:
                self._abb[p] = make_pulsar_ab_batch(self.sp, p)
            a_pf, b_pf = self._abb[p](tb, jnp.asarray(self.EV[p]), self.data, self.mask)
            LNL[:, p] = np.asarray(_rank_update(a_pf, b_pf, b_base_j[:, p], u_p_j[:, p],
                                                Minv_pp_j[p], Qb, sa, ab[:, p], self._const))
        self.n_calls += T
        return LNL, np.asarray(lnref)

    def marg_from_LNL(self, LNL, lnref, return_parts=False):
        T = LNL.shape[0]
        msum = np.zeros(T)
        m_all = np.empty((T, self.npsr)) if return_parts else None
        for a in range(self.npsr):
            M = LNL[:, a, :][:, self.order[a]]
            peak = np.maximum.reduceat(M, self.redidx[a], axis=1)
            logw = (peak - lnref[:, None]) + self.lnprior[a][None, :]
            mmax = logw.max(axis=1)
            m_a = mmax + np.log(np.sum(np.exp(logw - mmax[:, None]), axis=1))
            msum += m_a
            if return_parts:
                m_all[:, a] = m_a
        return (lnref + msum, m_all) if return_parts else (lnref + msum, None)

    def lnmarg(self, src_T):
        """src_T (T, 8*ncw) -> lnL_marg (T,). Chunked to ONE compiled particle shape."""
        src_T = np.atleast_2d(np.asarray(src_T, float))
        T = src_T.shape[0]
        Tc = self.T_CHUNK
        out = np.empty(T)
        lo, hi = TE.phi_bounds(self.P)
        for c0 in range(0, T, Tc):
            blk = src_T[c0:c0 + Tc]; nb = blk.shape[0]
            if nb < Tc:
                blk = np.concatenate([blk, np.tile(blk[-1], (Tc - nb, 1))], 0)
            blk = np.clip(blk, lo, hi)
            LNL, lnref = self.estep_batch(self.theta_of(blk))
            lm, _ = self.marg_from_LNL(LNL, lnref)
            out[c0:c0 + nb] = lm[:nb]
        return np.nan_to_num(out, nan=-1e300, posinf=-1e300, neginf=-1e300)

    def gate_fast_full(self, srcs, tol=1e-8):
        """lnL_marg with the pterm-only fast residual == lnL_marg with the masked residual, at
        pmask == 1. Run before any science that relies on the fast path."""
        if not self.all_on:
            return 0.0
        self.sp.enable_fast_full(True)
        a = self.lnmarg(srcs)
        self.sp.enable_fast_full(False)
        self._abb = {}                       # batched evaluators captured the old residual
        b = self.lnmarg(srcs)
        self.sp.enable_fast_full(True)
        self._abb = {}
        d = float(np.max(np.abs(a - b)))
        print(f"  [GATE fast-full] max|lnL_marg(fast) - lnL_marg(masked)| = {d:.3e} -> "
              f"{'PASS' if d < tol else 'FAIL'} (<{tol:.0e})", flush=True)
        if d >= tol:
            raise RuntimeError(f"fast-full residual path disagrees with the masked path by {d:.3e}")
        return d

    def posterior(self, src):
        """Full per-pulsar fringe posterior at ONE source vector (B_CERT grid)."""
        LNL, lnref = self.estep_batch(self.theta_of(np.atleast_2d(src)))
        C._finite("posterior LNL", LNL[0])
        post = self.FT.posterior(LNL[0])
        return post, LNL[0], float(lnref[0])


# ============================================================
# needle profiles
# ============================================================
def _src_offset(src0, k_list, axis, delta_scaled, scale):
    s = src0.copy()
    for k in k_list:
        s[8 * k + axis] += delta_scaled * scale[axis]
    return s


def mode_needle(args):
    print("=== B1 NEEDLE: lnL_marg profiles about truth in the registration axes ===", flush=True)
    P = C.build_b1_problem()
    tt = P["theta_truth"]; nd = P["n_dist"]
    src0 = tt[nd:].copy()
    scale = TE.phi_scale(P)
    loud = list(range(C.N_LOUD))

    noisy = args.noise
    mask1 = P["mask_one"]
    data_clean = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(mask1))
    if noisy:
        nz = P["nd"].draw(seed=args.seed)
        data = tuple(jnp.asarray(np.asarray(d) + np.asarray(n)) for d, n in zip(data_clean, nz))
        print(f"  data: CW(truth) + noise draw seed {args.seed} (white+rn+gwb)", flush=True)
    else:
        data = data_clean
        print(f"  data: zero-noise Asimov", flush=True)

    E = B1Marg(P, data, mask1)
    t0 = time.time(); g0 = float(E.lnmarg(src0[None])[0]); t1 = time.time() - t0
    lnl0 = float(P["amo"]["logL"](jnp.asarray(tt), data, jnp.asarray(mask1)))
    print(f"  lnL(truth) = {lnl0:.4f} ; lnL_marg(truth) = {g0:.4f} "
          f"(fringe entropy {g0-lnl0:+.2f} nat)", flush=True)
    print(f"  first chunk of {E.T_CHUNK} particles: {t1:.1f}s (incl. compiles)", flush=True)
    t0 = time.time(); _ = E.lnmarg(np.tile(src0, (E.T_CHUNK, 1))); tw = time.time() - t0
    print(f"  warm: {tw:.2f}s / {E.T_CHUNK} particles = {tw/E.T_CHUNK*1e3:.0f} ms per lnL_marg", flush=True)

    axes = [(I_FGW, "log10_fgw"), (I_MC, "log10_mc"), (I_COSGT, "cos_gwtheta")]
    res = {}
    for axis, nm in axes:
        # geometric ray in +/- delta, all loud sources moved together (the F2/B0.2 convention)
        ds = np.concatenate([-np.geomspace(args.dmax, args.dmin, args.n),
                             [0.0],
                             np.geomspace(args.dmin, args.dmax, args.n)])
        srcs = np.array([_src_offset(src0, loud, axis, d, scale) for d in ds])
        t0 = time.time(); g = E.lnmarg(srcs); dt = time.time() - t0
        C._finite(f"lnmarg {nm}", g[np.isfinite(g)])
        dg = g - g0
        res[nm] = (ds, dg)
        # half-width at -1 nat, on the positive side
        pos = ds > 0
        hw = np.nan
        idx = np.where(pos & (dg <= -1.0))[0]
        if len(idx):
            hw = ds[idx[0]]
        print(f"\n  -- axis {nm} ({len(ds)} pts, {dt:.0f}s) --", flush=True)
        print(f"     {'delta(scaled)':>14s} {'dlnL_marg':>12s}", flush=True)
        for d, v in zip(ds, dg):
            if d >= 0 and (d == 0 or d in ds[ds > 0][::max(1, args.n // 8)] or d == ds[-1]):
                print(f"     {d:14.2e} {v:12.2f}", flush=True)
        print(f"     needle half-width at -1 nat: {hw:.3e} scaled ; "
              f"dlnL_marg at delta={args.dmax:.0e}: {dg[-1]:.1f} nat", flush=True)
    np.savez(f"{CWT}/b1_needle_{'noisy' if noisy else 'asimov'}.npz",
             g0=g0, lnl0=lnl0, **{f"{k}_d": v[0] for k, v in res.items()},
             **{f"{k}_g": v[1] for k, v in res.items()})
    print(f"\n  saved b1_needle_{'noisy' if noisy else 'asimov'}.npz", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["needle"])
    ap.add_argument("--noise", action="store_true")
    ap.add_argument("--seed", type=int, default=90001)
    ap.add_argument("--dmin", type=float, default=1e-6)
    ap.add_argument("--dmax", type=float, default=1e-1)
    ap.add_argument("--n", type=int, default=24)
    a = ap.parse_args()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    if a.mode == "needle":
        return mode_needle(a)
    return 2


if __name__ == "__main__":
    sys.exit(main())
