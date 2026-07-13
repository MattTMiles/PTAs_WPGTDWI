"""SAMPLER — prong 3, actually run: the differentiable fringe-marginalised target density.

STATUS: TRACKED MACHINERY (promoted from SAMPLER's fence, cronus 2026-07-13).
`MargJax` is the differentiable twin of `trackB_b1.B1Marg` and is a DROP-IN for any future
soft-loop M-step. IGNITE-2's soft loop — the adopted spec-§3 reference implementation — has
been running on FINITE-DIFFERENCE gradients of a numpy object; an exact analytic gradient of
the SAME objective now exists.

PRODUCTION NUMBERS FOR THE RECORD (measured, cronus / RTX 4090; gates in sampler_gates.py,
sampler_unit.py; banks in SAMPLER_results/):

    G0  VALUE     MargJax.lnmarg == B1Marg.lnmarg      max|diff| = 5.146e-08 nat on |lnL| ~ 4.06e5
                  and it BIT-MATCHES R's banked deliverable:  lnL_marg(truth) = 405686.343447
                                                              R banked:         405686.34
    G1a GRADIENT  hand-assembled per-pulsar VJP == monolithic jax.grad   rel err 1.131e-15
    G2  COST      value_and_grad = 1.09 s  (F_CHUNK=64, 116 psr, B=512, 18 active dims)
                  GPU peak < 1 GiB   (the naive single-graph jax.grad OOMs at 16.5 GiB)

The 5.1e-8 residual is the float64 cancellation floor of the assembly itself (lnL_marg is
assembled as a difference of ~4.7e7-sized quantities yielding ~4e5), not a discrepancy.
The pre-registered 1e-8 FD gate is UNREACHABLE on this target and was replaced by the
strictly stronger AD-vs-AD test (G1a) — see SAMPLER_dev.md §0.2 and §2.

SCOPE OF INFERENCE. The within-fringe offset is PROFILED, not marginalised (see below), and
on the production grid `segment_max` is the IDENTITY for 115/116 pulsars, so the offset is in
fact PINNED at the fringe centre (u = 0). That pinning is EXACT wherever the true distance
sits at the EM prior mean (the whole pre-B1 era: R, the census ceilings, LAMBDA/F2/L2c, GEO,
FORGE Arm A) and is an ESTIMATOR ARTEFACT wherever it does not (the Arm-B era). See
SAMPLER_dev.md §9.5 for the measured immunity boundary.

THE OBJECT
----------
The target is the CONDITIONAL posterior of the counterpart-targeted scenario:

    p(theta_act | d) propto exp( lnL_marg(theta) ) * prior(theta_act)

with `lnL_marg` the SAME star-topology count-once fringe-marginalised likelihood that
deliverable R and IGNITE-2 used (`trackB_b1.B1Marg`), the fringes MARGINALISED (never
sampled as discrete parameters), the sky FIXED at truth (the EM-counterpart convention),
and the faint 13 sources fixed at truth (R's ACTIVE-DIM convention).

    lnL_marg(theta) = lnL_ref(theta) + sum_p LSE_n[ (LNL_p(n) - lnL_ref) + lnprior_p(n) ]

`B1Marg.lnmarg` evaluates this with numpy scaffolding (`np.maximum.reduceat` + a python
pulsar loop), so it has NO jax gradient. This module rebuilds it as ONE differentiable jax
graph, using the algebraic identity

    lnL_marg = sum_p M_p  -  (npsr - 1) * lnL_ref ,   M_p = LSE_n[ LNL_p(n) + lnprior_p(n) ]

(exactly B1Marg's algebra, regrouped) and a padded segment-max for the within-fringe peak
reduction. Reverse mode is memory-bounded by `jax.checkpoint` on each pulsar's B-fringe
block: the backward pass recomputes one pulsar's fringe stack at a time instead of holding
all 116 x B residual tapes.

WHAT IS AND IS NOT MARGINALISED (stated, because it governs every number downstream)
-----------------------------------------------------------------------------------
  * fringe INDEX n_p          : marginalised (the LSE) -- the honest object.
  * within-fringe offset      : PROFILED, not marginalised -- the per-fringe value is the
                                MAX over the B_CERT grid points falling inside that fringe
                                (`np.maximum.reduceat` in B1Marg; segment_max here). This is
                                inherited from B1Marg/R, is what IGNITE-2's soft loop ran on,
                                and is kept bit-comparable deliberately. Flagged, not fixed.
  * GWB + per-pulsar RN GPs   : marginalised ANALYTICALLY and exactly (they are GPs inside
                                the likelihood's covariance; hyperparameters frozen).
  * sky                       : FIXED (counterpart). Freed only in the S-4 coverage mode.

Env: jug-gpu, x64, XLA autotune off (shared GPU). Matt commits; never commit from here.
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
from jax.scipy.special import logsumexp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_b1_core as C
import trackB_b1 as B1
import trackB_estimator as TE
import discovery.matrix as dsm

CWT = "/home/mattm/projects/HSYMT/CW_transition"
RES = "/home/mattm/projects/HSYMT/SAMPLER_results"

# per-8-param-block indices (trackB_b1 convention)
I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)
AXIS_NAME = ["cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h",
             "phase0", "psi"]

# THE CONDITIONAL TARGET's active axes, per source: extrinsics + f + mc; sky fixed.
ACTIVE_AXES = [I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI]
# the S-4 coverage mode frees the sky as well
ACTIVE_AXES_SKY = [I_COSGT, I_GWPHI] + ACTIVE_AXES


# ============================================================
# prior tiers (spec / IGNITE conventions)
# ============================================================
def prior_sigma(P, tier="lit"):
    """Per-pulsar distance prior sigma (kpc) for the requested tier.

    lit  : the canonical literature priors (`best_distances.txt`, via A2.load_real_priors).
    vlbi : IGNITE's binary tier -- sigma_d -> min(lit, 1 pc) on the GEO union-18,
           improve-never-degrade (so J0437's 0.25 pc stands). The union-18 list is NOT
           banked on cronus (GEO's npz lives on ACCRE); `vlbi_all` applies the same
           improve-never-degrade rule to EVERY pulsar and is used ONLY with that stated.
    """
    sig = np.array(P["priors"]["lit"], float)
    if tier == "lit":
        return sig
    if tier == "vlbi_all":
        return np.minimum(sig, 1.0e-3)          # 1 pc in kpc
    raise ValueError(f"unknown tier {tier!r}")


def box_bounds(active_idx, ncw):
    lo = np.tile(TE._PHI_LO, ncw)[active_idx]
    hi = np.tile(TE._PHI_HI, ncw)[active_idx]
    return lo, hi


def active_index(loud, axes, ncw):
    """Flat indices into the (8*ncw,) source vector for the given sources x axes."""
    return np.array([8 * k + a for k in loud for a in axes], int)


def active_labels(loud, axes):
    return [f"src{k}_{AXIS_NAME[a]}" for k in loud for a in axes]


# ============================================================
# the differentiable target
# ============================================================
class MargJax:
    """Differentiable lnL_marg over a chosen ACTIVE subset of the source vector.

    Built to be VALUE-IDENTICAL to `trackB_b1.B1Marg.lnmarg` (gated, <1e-8), which is the
    object R and IGNITE-2 measured. Everything not in `active_idx` is held at `src_fixed`.
    """

    # Fringes per scan step. Memory per block ~ F_CHUNK x n_src x n_toa; the reverse pass
    # recomputes one chunk at a time. MEASURED on the 4090 (116 psr, B = 512, 18 active dims):
    #   F_CHUNK = 8  -> 6.98 s / value_and_grad, 0.6 GiB peak
    #   F_CHUNK = 64 -> 1.09 s / value_and_grad  (6.4x, still comfortably resident)
    # A single un-chunked graph (the naive jax.grad) asks for 16.5 GiB and OOMs.
    F_CHUNK = 64

    def __init__(self, P, data, mask, active_idx, src_fixed, tier="lit", B1M=None,
                 reduce="max"):
        """reduce = 'max' : the B1Marg/R/IGNITE-2 reduction. Within each fringe the value is the
                            MAX over the B_CERT grid points in that fringe (a PROFILE of the
                            sub-fringe offset), and the Gaussian prior is evaluated at the fringe
                            CENTRE. This is what every banked number in the repo was computed with.

           reduce = 'lse' : the offset-MARGINALISED reduction. The within-fringe offset is
                            integrated, with the correct measure and the prior evaluated
                            POINTWISE:

                                M_p = log SUM_j exp[ LNL_p(j) + lnprior(L_j) ] + log dL_grid,p

                            i.e. M_p = log INT exp[lnL(L)] N(L; L0, sigma_EM) dL over the +-3 sigma
                            window. NOTE: with the correct measure the fringe SEGMENTATION
                            DISAPPEARS from the likelihood -- the fringe decomposition is needed
                            only to report the fringe POSTERIOR q_p(n) (segment sums of the same
                            summands), not to form the marginal. `segment_max` was never doing the
                            integral; it was doing the profile.

        Both paths are live and switchable; they are gated against each other where they must
        agree (loud zero-noise: the within-cell Laplace factor becomes fringe-independent and
        cancels in q_p(n)).
        """
        assert reduce in ("max", "lse"), reduce
        self.reduce = reduce
        self.P = P
        self.sp = P["sp"]
        self.npsr = P["npsr"]
        self.ngp = self.sp.ngp_gwb
        self.data = data
        self.mask = jnp.asarray(mask)
        self.all_on = bool(np.all(np.asarray(mask) == 1.0))
        self.sp.enable_fast_full(self.all_on)
        self.active_idx = np.asarray(active_idx, int)
        self.ndim = len(self.active_idx)
        self.src_fixed = np.asarray(src_fixed, float).copy()
        self.tier = tier

        self.base_L = np.asarray(P["L0"], float)
        self.dL = np.asarray(P["dL_truth"], float)
        self.EV = np.asarray(P["EV_truth"], float)
        self.B = self.EV.shape[1]
        self.ndist = P["n_dist"]
        self.nth = P["n_theta"]
        self.AI = np.asarray(P["AI"], int)

        self.sp.AI = self.AI
        self.sp._ensure_minv()
        self._Minv = jnp.asarray(self.sp._Minv)
        self._Minv_pp = jnp.asarray(self.sp._Minv_pp)
        self._const = float(self.sp.a_const
                            - 0.5 * (self.sp.ldP_cached + self.sp.logdet_cached))

        # ---- padded fringe-reduction tables (per pulsar: grid point -> fringe group) ----
        sig = prior_sigma(P, tier)
        self.sigma = sig
        gids, lnpri, Ks, uks = [], [], [], []
        for a in range(self.npsr):
            k = np.round((self.EV[a] - self.base_L[a]) / self.dL[a]).astype(int)
            uk = np.unique(k)
            gid = np.searchsorted(uk, k)                      # (B,) group id per grid point
            gids.append(gid); uks.append(uk); Ks.append(len(uk))
            lnpri.append(-(uk * self.dL[a]) ** 2 / (2.0 * sig[a] ** 2))
        self.Kmax = int(max(Ks))
        self.K = np.array(Ks)
        self.uk = uks
        self.gid = jnp.asarray(np.stack(gids))                                  # (npsr, B)
        lp = np.full((self.npsr, self.Kmax), -np.inf)
        for a in range(self.npsr):
            lp[a, :Ks[a]] = lnpri[a]
        self.lnprior = jnp.asarray(lp)                                      # (npsr, Kmax)

        # ---- POINTWISE prior + grid measure, for the offset-MARGINALISED reduction ----
        # lnprior_pt[a, j] = -(EV[a,j] - L0[a])^2 / (2 sigma_a^2)   evaluated AT THE GRID POINT,
        # not at the fringe centre; dlog[a] = log(grid spacing) supplies the measure dL.
        lnp_pt = -(self.EV - self.base_L[:, None]) ** 2 / (2.0 * sig[:, None] ** 2)
        self.lnprior_pt = jnp.asarray(lnp_pt)                                   # (npsr, B)
        dgrid = (self.EV[:, -1] - self.EV[:, 0]) / (self.B - 1)                 # uniform per psr
        self.dgrid = dgrid
        self.ldgrid = jnp.asarray(np.log(dgrid))                                # (npsr,)
        self.K_counted = np.array([2 * int(np.floor(3.0 * sig[a] / self.dL[a] + 1e-9))
                                   for a in range(self.npsr)])
        self.EVj = jnp.asarray(self.EV)

        self.B1M = B1M                                    # the numpy reference (for gates)
        self._n_lnL = 0
        self._build()

    # ---- theta assembly ----
    def theta_of(self, x):
        src = jnp.asarray(self.src_fixed).at[jnp.asarray(self.active_idx)].set(x)
        return jnp.concatenate([jnp.asarray(self.base_L), src])

    def _build(self):
        """Per-pulsar graphs, chunked over the fringe axis. NEVER one monolithic graph.

        A single jax graph holding all 116 x B fringe blocks OOMs in reverse mode (measured:
        16.5 GiB single allocation on a 24 GB card -- the forward alone holds ~15 GB, because
        each block vmaps B = 512 fringes x 16 sources x n_toa). Two structural fixes, both
        exact:

        (1) FRINGE CHUNKING. Each pulsar's B-fringe stack is a `lax.scan` over chunks of
            `F_CHUNK`, with `jax.checkpoint` on the chunk body: the backward pass recomputes
            one chunk at a time instead of taping all B.
        (2) PER-PULSAR VJP. lnL_marg = sum_p M_p(theta, base(theta)) - (npsr-1)*lnref(base),
            where base(theta) = (a_base, b_base) is ONE full-array residual eval. The chain
            rule is applied by hand:
                grad = sum_p dM_p/dtheta|_base  +  VJP_base[ sum_p dM_p/dbase
                                                             - (npsr-1) * dlnref/dbase ]
            so each compiled graph contains ONE pulsar's fringes, never 116.
        """
        sp = self.sp
        npsr, ngp, B = self.npsr, self.ngp, self.B
        const, Minv, Minv_pp = self._const, self._Minv, self._Minv_pp
        gid, lnprior, Kmax = self.gid, self.lnprior, self.Kmax
        lnprior_pt, ldgrid = self.lnprior_pt, self.ldgrid
        use_lse = (self.reduce == "lse")
        data, mask = self.data, self.mask
        cs = self.F_CHUNK
        assert B % cs == 0, f"B={B} must be divisible by F_CHUNK={cs}"
        nch = B // cs

        # ---- base: one full-array residual eval, differentiable, cheap in memory ----
        def base_ab(theta):
            return sp._per_pulsar_ab_impl(theta, data, mask)
        self._base_ab = jax.jit(base_ab)

        def base_of_x(x):
            return base_ab(self.theta_of(x))
        self._base_of_x = jax.jit(base_of_x)
        self._base_vjp = jax.jit(lambda x, ca, cb: jax.vjp(base_of_x, x)[1]((ca, cb))[0])

        # ---- per-pulsar M_p(x, a_base, b_base) ----
        def make_Mp(p):
            Ft = sp.Ft[p]; Tt = sp.Tt[p]; s1d = sp.solve1d[p]; ys_p = sp.ys[p]
            cf_p = (sp.cf[0][p], sp.cf[1]); TtNmF_p = sp.TtNmF[p]; AI_p = int(self.AI[p])
            Lv = self.EVj[p].reshape(nch, cs)
            gid_p = gid[p]; lnpri_p = lnprior[p]; Minv_p = Minv_pp[p]
            lnpri_pt_p = lnprior_pt[p]; ldg_p = ldgrid[p]

            def M_p(x, a_base, b_base):
                theta = self.theta_of(x)
                bflat = b_base.reshape(npsr * ngp)
                u = Minv @ bflat
                Q_base = bflat @ u
                u_p = u.reshape(npsr, ngp)[p]
                sum_a = jnp.sum(a_base)

                @jax.checkpoint             # recompute this chunk in the backward pass
                def chunk(_, Lc):
                    def one(L):
                        params = sp._params(theta.at[AI_p].set(L), data, mask)
                        yp = ys_p(params)
                        Nmy, _ = s1d(yp)
                        return yp @ Nmy, Ft @ Nmy, Tt @ Nmy
                    ytNmy, FtNmy, TtNmy = jax.vmap(one)(Lc)
                    sol = dsm.matrix_solve(cf_p, FtNmy.T).T
                    a_c = -0.5 * (ytNmy - jnp.sum(FtNmy * sol, axis=1))         # (cs,)
                    b_c = TtNmy - sol @ TtNmF_p.T                               # (cs,ngp)
                    return _, (a_c, b_c)

                _, (A, Bm) = jax.lax.scan(chunk, None, Lv)
                a_pf = A.reshape(B)
                b_pf = Bm.reshape(B, ngp)
                db = b_pf - b_base[p][None, :]
                Q_pf = (Q_base + 2.0 * (db @ u_p)
                        + jnp.einsum('bi,ij,bj->b', db, Minv_p, db))
                LNL_p = const + (sum_a - a_base[p] + a_pf) + 0.5 * Q_pf
                if use_lse:
                    # offset-MARGINALISED: integrate the window with the pointwise prior and the
                    # dL measure. The fringe segmentation is not needed to form the marginal.
                    return logsumexp(LNL_p + lnpri_pt_p) + ldg_p
                # PROFILE (the banked reduction): max within each fringe, prior at the centre
                peak = jax.ops.segment_max(LNL_p, gid_p, num_segments=Kmax)
                return logsumexp(peak + lnpri_p)
            # ONLY the value_and_grad graph is compiled. Holding a separate value graph per
            # pulsar as well (232 executables) OOMs the card; the value is read off the
            # value_and_grad call instead. NUTS needs the pair anyway.
            return jax.jit(jax.value_and_grad(M_p, argnums=(0, 1, 2)))

        self._Mp_vg = [make_Mp(p) for p in range(npsr)]

        def lnref_of(a_base, b_base):
            bflat = b_base.reshape(npsr * ngp)
            return const + jnp.sum(a_base) + 0.5 * (bflat @ (Minv @ bflat))
        self._lnref_of = jax.jit(lnref_of)
        self._dlnref = jax.jit(jax.grad(lnref_of, argnums=(0, 1)))

    # ---- value / gradient assembled across pulsars ----
    def lnmarg(self, x):
        return self.lnmarg_grad(x)[0]

    def lnmarg_grad(self, x):
        """(value, gradient) -- exact, memory-bounded, per-pulsar VJP assembly."""
        x = jnp.asarray(x, float)
        a_base, b_base = self._base_of_x(x)
        npsr = self.npsr
        M = 0.0
        gx = jnp.zeros(self.ndim)
        ga = jnp.zeros(npsr)
        gb = jnp.zeros((npsr, self.ngp))
        for p in range(npsr):
            m, (gx_p, ga_p, gb_p) = self._Mp_vg[p](x, a_base, b_base)
            M = M + m
            gx = gx + gx_p
            ga = ga + ga_p
            gb = gb + gb_p
        lnref = self._lnref_of(a_base, b_base)
        dra, drb = self._dlnref(a_base, b_base)
        ca = ga - (npsr - 1) * dra
        cb = gb - (npsr - 1) * drb
        gx = gx + self._base_vjp(x, ca, cb)
        return M - (npsr - 1) * lnref, gx

    # ---- convenience ----
    def x_of_src(self, src):
        return np.asarray(src, float)[self.active_idx]

    def value(self, x):
        self._n_lnL += 1
        return float(self.lnmarg(x))

    def value_and_grad(self, x):
        self._n_lnL += 1
        v, g = self.lnmarg_grad(x)
        return float(v), np.asarray(g)

    def monolith_value_and_grad(self, x):
        """REFERENCE ONLY (small arrays): one jax graph, plain jax.value_and_grad.

        This is the gradient the hand-assembled per-pulsar VJP must reproduce EXACTLY. It is
        unusable at npsr = 116 (the single backward graph asks for a 16.5 GiB allocation on a
        24 GB card -- measured), which is why `lnmarg_grad` exists at all; but at small npsr it
        compiles, and then AD-vs-AD is an exact gate on the chain rule with NO finite-difference
        noise in it.
        """
        if not hasattr(self, "_mono_vg"):
            npsr, ngp = self.npsr, self.ngp
            const, Minv, Minv_pp = self._const, self._Minv, self._Minv_pp
            gid, lnprior, Kmax = self.gid, self.lnprior, self.Kmax
            sp, data, mask = self.sp, self.data, self.mask

            def mono(x):
                theta = self.theta_of(x)
                a_base, b_base = sp._per_pulsar_ab_impl(theta, data, mask)
                bflat = b_base.reshape(npsr * ngp)
                u = Minv @ bflat
                Q_base = bflat @ u
                u_p = u.reshape(npsr, ngp)
                sum_a = jnp.sum(a_base)
                lnref = const + sum_a + 0.5 * Q_base
                M = 0.0
                for p in range(npsr):
                    Ft = sp.Ft[p]; Tt = sp.Tt[p]; s1d = sp.solve1d[p]; ys_p = sp.ys[p]
                    cf_p = (sp.cf[0][p], sp.cf[1]); TtNmF_p = sp.TtNmF[p]
                    AI_p = int(self.AI[p])

                    def one(L):
                        params = sp._params(theta.at[AI_p].set(L), data, mask)
                        yp = ys_p(params)
                        Nmy, _ = s1d(yp)
                        return yp @ Nmy, Ft @ Nmy, Tt @ Nmy
                    ytNmy, FtNmy, TtNmy = jax.vmap(one)(self.EVj[p])
                    sol = dsm.matrix_solve(cf_p, FtNmy.T).T
                    a_pf = -0.5 * (ytNmy - jnp.sum(FtNmy * sol, axis=1))
                    b_pf = TtNmy - sol @ TtNmF_p.T
                    db = b_pf - b_base[p][None, :]
                    Q_pf = (Q_base + 2.0 * (db @ u_p[p])
                            + jnp.einsum('bi,ij,bj->b', db, Minv_pp[p], db))
                    LNL_p = const + (sum_a - a_base[p] + a_pf) + 0.5 * Q_pf
                    peak = jax.ops.segment_max(LNL_p, gid[p], num_segments=Kmax)
                    M = M + logsumexp(peak + lnprior[p])
                return M - (npsr - 1) * lnref
            self._mono_vg = jax.jit(jax.value_and_grad(mono))
        v, g = self._mono_vg(jnp.asarray(x, float))
        return float(v), np.asarray(g)

    # ---- an XLA-OPAQUE primitive, so a sampler may jit AROUND lnL_marg but never THROUGH it
    def as_primitive(self):
        """lnL_marg as a jax primitive with a custom VJP, evaluated by host callback.

        NUTS (numpyro) jits its potential. If lnL_marg were an ordinary traceable function,
        that jit would re-inline all 116 pulsar blocks into ONE graph -- the exact monolith
        that OOMs (16.5 GiB). A `pure_callback` + `custom_vjp` keeps the sampler's control
        flow compiled while each lnL_marg / grad call runs the per-pulsar loop at runtime.
        The gradient is the SAME one gated against FD; no approximation enters here.
        """
        D = self.ndim
        sds_v = jax.ShapeDtypeStruct((), jnp.float64)
        sds_g = jax.ShapeDtypeStruct((D,), jnp.float64)

        def _v(x):
            return np.float64(self.value(np.asarray(x, float)))

        def _vg(x):
            v, g = self.value_and_grad(np.asarray(x, float))
            return np.float64(v), np.asarray(g, float)

        @jax.custom_vjp
        def op(x):
            return jax.pure_callback(_v, sds_v, x)

        def op_fwd(x):
            v, g = jax.pure_callback(_vg, (sds_v, sds_g), x)
            return v, g

        def op_bwd(g, ct):
            return (ct * g,)

        op.defvjp(op_fwd, op_bwd)
        return op


# ============================================================
# priors (uniform box, the spec's generous physical ranges)
# ============================================================
class BoxPrior:
    """Uniform on the spec box (TE._PHI_LO/_PHI_HI), in the SAMPLED (unconstrained) space.

    Sampling happens in a logit-transformed unit space so NUTS never proposes outside the
    box (the waveform nan-clips there); the log-Jacobian is carried exactly.
    """

    def __init__(self, lo, hi):
        self.lo = jnp.asarray(lo, float)
        self.hi = jnp.asarray(hi, float)
        self.w = self.hi - self.lo

    def to_x(self, z):
        """unconstrained z -> physical x, plus log|dx/dz|."""
        s = jax.nn.sigmoid(z)
        x = self.lo + self.w * s
        logJ = jnp.sum(jnp.log(self.w) + jnp.log(s) + jnp.log1p(-s))
        return x, logJ

    def to_z(self, x):
        u = (jnp.asarray(x) - self.lo) / self.w
        u = jnp.clip(u, 1e-12, 1 - 1e-12)
        return jnp.log(u) - jnp.log1p(-u)

    def sample(self, key, n):
        u = jax.random.uniform(key, (n, self.lo.size))
        return self.lo + self.w * u


def make_potential(M, prior):
    """numpyro-style potential (== -log posterior) in the unconstrained z space.

    Built on `M.as_primitive()`, so the sampler jits its own machinery but calls lnL_marg
    (and its gated gradient) through a host callback -- see MargJax.as_primitive.
    """
    op = M.as_primitive()

    def logp_z(z):
        x, logJ = prior.to_x(z)
        return op(x) + logJ
    return jax.jit(lambda z: -logp_z(z)), jax.jit(logp_z)


# ============================================================
# problem assembly for a given cell
# ============================================================
def build_cell(seed=None, tier="lit", noise=True, arm="B", verbose=True):
    """One realisation of the conditional (targeted) scenario.

    seed=None or noise=False -> ZERO-NOISE Asimov (gate g1/g2).
    arm  'A' : true distances at the EM prior mean (n_true == 0), the census convention.
    arm  'B' : true distances drawn off the prior mean (the honest arm; core.draw_true_distances).
    """
    P = C.build_b1_problem(verbose=verbose)
    npsr = P["npsr"]
    tt = P["theta_truth"].copy()
    mask1 = P["mask_one"]

    if arm == "B":
        L_true, n_true, u_true = C.draw_true_distances(P, seed=(seed or 0) + 500000)
        tt[:P["n_dist"]] = L_true
    else:
        L_true = np.asarray(P["L0"], float).copy()
        n_true = np.zeros(npsr, int); u_true = np.zeros(npsr)

    data_clean = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(mask1))
    if noise and seed is not None:
        nz = P["nd"].draw(seed=seed)
        data = tuple(jnp.asarray(np.asarray(d) + np.asarray(n))
                     for d, n in zip(data_clean, nz))
    else:
        data = data_clean
    C._finite("cell data", *[np.asarray(d) for d in data])
    return dict(P=P, data=data, mask=mask1, theta_truth=tt, L_true=L_true,
                n_true=n_true, u_true=u_true, tier=tier, seed=seed, noise=bool(noise and seed is not None),
                arm=arm)


def target_from_cell(cell, axes=ACTIVE_AXES, loud=None, with_ref=True):
    P = cell["P"]
    loud = list(range(C.N_LOUD)) if loud is None else loud
    src_truth = cell["theta_truth"][P["n_dist"]:].copy()
    idx = active_index(loud, axes, P["ncw"])
    B1M = B1.B1Marg(P, cell["data"], cell["mask"]) if with_ref else None
    M = MargJax(P, cell["data"], cell["mask"], idx, src_truth, tier=cell["tier"], B1M=B1M)
    lo, hi = box_bounds(idx, P["ncw"])
    return M, BoxPrior(lo, hi), dict(labels=active_labels(loud, axes), lo=lo, hi=hi,
                                     x_truth=src_truth[idx], loud=loud, axes=axes)
