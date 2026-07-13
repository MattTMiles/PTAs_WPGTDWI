"""SAMPLER — THE SUB-FRINGE OFFSET TEST: is the within-fringe PROFILE cosmetic or material?

THE QUESTION
------------
Every fringe-marginalised number in this repo -- R's `f = 6.9e-7`, IGNITE/IGNITE-2's floors,
counts and q_max, the soft loop's trajectories -- is computed with `B1Marg`'s reduction, which
takes the MAX over the B_CERT grid points inside each fringe (`np.maximum.reduceat`) and puts the
Gaussian prior at the fringe CENTRE. That is a PROFILE of the sub-fringe offset, not a marginal.
The archaeology addendum (SAMPLER_dev.md §8) showed this is the SAME approximation the 2023-era
`snap_to_peak` made (a +-0.6 dL argmax): the modern formulation inherited it unexamined.

The honest reduction integrates the offset with the correct measure and a POINTWISE prior:

    PROFILE (banked):     M_p = LSE_n [ max_{j in n} LNL_p(j)  +  lnprior(centre_n) ]
    MARGINAL (this):      M_p = LSE_j [ LNL_p(j) + lnprior(L_j) ] + log dL_grid,p
                              = log INT exp[lnL(L)] N(L; L0, sigma_EM) dL

Note what happens: WITH THE CORRECT MEASURE THE FRINGE SEGMENTATION DISAPPEARS from the marginal
likelihood. The fringe decomposition survives only in the reported fringe POSTERIOR q_p(n) (segment
sums of the same summands). `segment_max` was never computing the integral; it was profiling.

THE LOAD-BEARING QUESTION (b.i): does the profile CANCEL in R's needle/plateau ratio?
    Delta(theta) = lnL_marg^LSE(theta) - lnL_marg^MAX(theta)
If Delta is the same at the needle and over the plateau, it cancels in ln Z_plateau - ln Z_needle
and f is untouched (COSMETIC). If it differs, f moves by exp(Delta_plateau - Delta_needle)
(MATERIAL). This is answered from BANKED artifacts -- no new realisations:
    * the needle: Delta at truth (the needle is razor-thin, sigma_sky ~ 1e-6 scaled, so Delta is
      flat across it -- CHECKED, not assumed, by evaluating Delta around the needle box);
    * the plateau: R banked its SMC particles (`trackB_R_zplateau_s{0,1}.npz`, X_final (256,24),
      L_final (256,)). Since those particles approximate the OLD plateau posterior,
          ln( Z_plateau^new / Z_plateau^old ) = ln E_old[ exp(Delta) ] ~= ln mean_i exp(Delta_i).

MODES
-----
  gate    the two reductions must agree on q_p(n) in the LOUD ZERO-NOISE limit (the within-cell
          Laplace factor becomes fringe-independent there and cancels in the posterior).
  ratio   (b.i) recompute R's f under the marginal reduction, from R's banked particles.

Env: jug-gpu, x64. Read-only w.r.t. every banked artifact. Matt commits; never commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time, argparse
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import sampler_core as S
import trackB_b1 as B1
import trackB_b1_core as C
import trackB_estimator as TE

CWT = S.CWT
RES = S.RES


# ============================================================
# the two reductions, applied to the SAME LNL (npsr, B). numpy, exact, no GPU.
# ============================================================
class Reductions:
    def __init__(self, P, EV, dL, prior_key="lit"):
        self.npsr = P["npsr"]
        self.L0 = np.asarray(P["L0"], float)
        self.EV = np.asarray(EV, float)
        self.dL = np.asarray(dL, float)
        self.B = self.EV.shape[1]
        sig = np.asarray(P["priors"][prior_key], float)
        self.sig = sig
        self.order, self.redidx, self.uk, self.lnprior_c = [], [], [], []
        for a in range(self.npsr):
            k = np.round((self.EV[a] - self.L0[a]) / self.dL[a]).astype(int)
            o = np.argsort(k, kind="stable")
            uk, i0 = np.unique(k[o], return_index=True)
            self.order.append(o); self.redidx.append(i0); self.uk.append(uk)
            self.lnprior_c.append(-(uk * self.dL[a]) ** 2 / (2.0 * sig[a] ** 2))
        # pointwise prior + measure
        self.lnprior_pt = -(self.EV - self.L0[:, None]) ** 2 / (2.0 * sig[:, None] ** 2)
        self.dgrid = (self.EV[:, -1] - self.EV[:, 0]) / (self.B - 1)
        self.ldgrid = np.log(self.dgrid)

    @staticmethod
    def _lse(x, axis=-1):
        m = np.max(x, axis=axis, keepdims=True)
        return (m + np.log(np.sum(np.exp(x - m), axis=axis, keepdims=True))).squeeze(axis)

    def totals(self, LNL, lnref):
        """(total_max, total_lse) for ONE realisation's LNL (npsr, B)."""
        m_max = np.zeros(self.npsr)
        m_lse = np.zeros(self.npsr)
        for a in range(self.npsr):
            M = LNL[a][self.order[a]]
            peak = np.maximum.reduceat(M, self.redidx[a])
            m_max[a] = self._lse((peak - lnref) + self.lnprior_c[a])
            m_lse[a] = self._lse((LNL[a] - lnref) + self.lnprior_pt[a]) + self.ldgrid[a]
        return lnref + m_max.sum(), lnref + m_lse.sum(), m_max, m_lse

    def qpost(self, LNL):
        """fringe posteriors q_p(n) under BOTH reductions (lists of length npsr)."""
        q_max, q_lse = [], []
        for a in range(self.npsr):
            M = LNL[a][self.order[a]]
            peak = np.maximum.reduceat(M, self.redidx[a])
            lw = (peak - peak.max()) + self.lnprior_c[a]
            w = np.exp(lw - lw.max()); w /= w.sum()
            q_max.append(w)
            # marginal: segment-SUM the same summands (prior pointwise, measure constant per psr)
            lw2 = LNL[a] + self.lnprior_pt[a]
            lw2 = lw2 - lw2.max()
            wj = np.exp(lw2)
            seg = np.add.reduceat(wj[self.order[a]], self.redidx[a])
            seg = seg / seg.sum()
            q_lse.append(seg)
        return q_max, q_lse


def build(args):
    cell = S.build_cell(seed=None, tier="lit", noise=False, arm="A")
    P = cell["P"]
    B1.B1Marg.T_CHUNK = args.chunk
    B1M = B1.B1Marg(P, cell["data"], cell["mask"])
    R = Reductions(P, P["EV_truth"], P["dL_truth"], "lit")
    return cell, P, B1M, R


def loud_srcs(P, src, h_over):
    """copy of the source vector with the 3 loud log10_h overridden (IGNITE's loudness axis)."""
    s = np.asarray(src, float).copy()
    for k in range(C.N_LOUD):
        s[8 * k + S.I_H] = h_over
    return s


def mode_gate(args):
    """The two reductions must agree on q_p(n) in the LOUD zero-noise limit."""
    print("=== OFFSET GATE: profile vs marginal reduction, q_p(n) agreement vs loudness ===",
          flush=True)
    cell, P, B1M, R = build(args)
    src0 = cell["theta_truth"][P["n_dist"]:].copy()
    rows = []
    for h in args.h:
        src = loud_srcs(P, src0, h)
        LNL, lnref = B1M.estep_batch(B1M.theta_of(src[None]))
        C._finite("gate LNL", LNL)
        LNL = LNL[0]; lnref = float(lnref[0])
        t_max, t_lse, _, _ = R.totals(LNL, lnref)
        q_max, q_lse = R.qpost(LNL)
        dq = max(float(np.max(np.abs(a - b))) for a, b in zip(q_max, q_lse))
        qm_med = float(np.median([q.max() for q in q_max]))
        # the census/anchor pulsars specifically
        cen = {nm: (float(q_max[P["names"].index(nm)].max()),
                    float(q_lse[P["names"].index(nm)].max())) for nm in C.CENSUS}
        rows.append((h, t_max, t_lse, t_lse - t_max, dq, qm_med))
        print(f"  h = {h:6.2f}: lnL_marg  profile {t_max:14.4f}  marginal {t_lse:14.4f}  "
              f"Delta {t_lse-t_max:10.4f}", flush=True)
        print(f"             max_p max_n |q_max - q_lse| = {dq:.3e}   "
              f"median q_max(profile) = {qm_med:.3f}", flush=True)
        for nm, (a, b) in cen.items():
            print(f"             {nm:14s} q_max: profile {a:.4f}  marginal {b:.4f}  "
                  f"diff {abs(a-b):.2e}", flush=True)
    A = np.array(rows)
    np.savez(f"{RES}/offset_gate.npz", h=A[:, 0], total_max=A[:, 1], total_lse=A[:, 2],
             delta=A[:, 3], max_dq=A[:, 4], med_qmax=A[:, 5])
    loud = A[np.argmax(A[:, 0])]
    ok = loud[4] < args.tol
    print(f"\n  [GATE] at the loudest h = {loud[0]:.2f}: max|dq| = {loud[4]:.3e} -> "
          f"{'PASS' if ok else 'FAIL'} (<{args.tol:.0e})", flush=True)
    print(f"  (the reductions must converge as q_p -> delta; they need NOT agree in absolute "
          f"lnL_marg -- the marginal carries a per-pulsar within-cell volume factor)", flush=True)
    return 0 if ok else 1


def mode_ratio(args):
    """(b.i) R's f under the marginal reduction, from R's BANKED particles. No new realisations."""
    print("=== OFFSET IMPACT (b.i): does the profile cancel in R's needle/plateau ratio? ===",
          flush=True)
    cell, P, B1M, R = build(args)
    nd = P["n_dist"]
    src_truth = cell["theta_truth"][nd:].copy()

    def delta_of(src_batch):
        """Delta = lnL_marg^LSE - lnL_marg^MAX at each source vector (T, 8*ncw)."""
        out = np.zeros(len(src_batch)); tm = np.zeros(len(src_batch)); tl = np.zeros(len(src_batch))
        Tc = B1M.T_CHUNK
        for c0 in range(0, len(src_batch), Tc):
            blk = src_batch[c0:c0 + Tc]
            nb = len(blk)
            if nb < Tc:
                blk = np.concatenate([blk, np.tile(blk[-1], (Tc - nb, 1))], 0)
            LNL, lnref = B1M.estep_batch(B1M.theta_of(blk))
            C._finite("ratio LNL", LNL)
            for i in range(nb):
                a, b, _, _ = R.totals(LNL[i], float(lnref[i]))
                tm[c0 + i] = a; tl[c0 + i] = b; out[c0 + i] = b - a
        return out, tm, tl

    # ---- the needle: Delta at truth, and its FLATNESS across the needle box ----
    Zn = np.load(f"{CWT}/trackB_R_znaddle.npz")
    print(f"  R banked: lnZ_needle(quad) = {float(Zn['lnZ_quad']):.4f}, "
          f"lnZ_needle(laplace) = {float(Zn['lnZ_laplace']):.4f}", flush=True)
    reg = np.asarray(Zn["reg_idx"], int)        # the 6 active sky dims (into the 24-vector)
    rng = np.random.default_rng(0)
    # probe Delta at truth and at +-3 needle sigmas along each active eigen-direction
    evals = np.asarray(Zn["evals"], float); evecs = np.asarray(Zn["evecs"], float)  # (24,6)
    probes = [src_truth.copy()]
    for j in range(evecs.shape[1]):
        sig_j = 1.0 / np.sqrt(max(evals[j], 1e-300))
        for s in (-3.0, 3.0):
            v = src_truth.copy()
            v[:24] = v[:24] + s * sig_j * evecs[:, j]
            probes.append(v)
    probes = np.array(probes)
    d_needle, tm_n, tl_n = delta_of(probes)
    print(f"  Delta(truth)            = {d_needle[0]:.4f} nat", flush=True)
    print(f"  Delta over the needle box (+-3 sigma, {len(d_needle)-1} probes): "
          f"mean {d_needle[1:].mean():.4f}, sd {d_needle[1:].std():.4f}, "
          f"max|dev from truth| {np.max(np.abs(d_needle[1:]-d_needle[0])):.4f} nat", flush=True)

    # ---- the plateau: Delta at R's banked SMC particles ----
    Ds, Ns = [], []
    for s in (0, 1):
        Zp = np.load(f"{CWT}/trackB_R_zplateau_s{s}.npz")
        X = np.asarray(Zp["X_final"], float)            # (256, 24) loud-source params
        srcs = np.tile(src_truth, (len(X), 1))
        srcs[:, :24] = X
        d, tm, tl = delta_of(srcs)
        Ds.append(d); Ns.append(len(d))
        print(f"  seed {s}: {len(d)} plateau particles; Delta mean {d.mean():.4f}, "
              f"sd {d.std():.4f}, min {d.min():.4f}, max {d.max():.4f}", flush=True)
    D = np.concatenate(Ds)

    # ---- the reweighting identity ----
    # ln(Z_new/Z_old) = ln E_old[exp(Delta)]  (particles ~ the OLD plateau posterior)
    m = D.max()
    dlnZ_plateau = float(m + np.log(np.mean(np.exp(D - m))))
    dlnZ_needle = float(d_needle[0])
    Zres = np.load(f"{CWT}/trackB_R_referendum_result.npz")
    lnZn_old = float(Zres["lnZn_quad"]); lnZp_old = float(Zres["lnZp"])
    f_old = float(Zres["f_quad"])
    lnZn_new = lnZn_old + dlnZ_needle
    lnZp_new = lnZp_old + dlnZ_plateau
    r_old = lnZp_old - lnZn_old
    r_new = lnZp_new - lnZn_new
    f_new = 1.0 / (1.0 + np.exp(r_new))
    print(f"\n  ln Z_plateau - ln Z_needle :  OLD {r_old:+.3f} nat   NEW {r_new:+.3f} nat   "
          f"shift {r_new-r_old:+.3f} nat", flush=True)
    print(f"  f = Z_n/(Z_n+Z_p)          :  OLD {f_old:.3e}          NEW {f_new:.3e}", flush=True)
    print(f"  R's own plateau evidence error (2 seeds) = +-{float(Zres['lnZp_err']):.3f} nat; "
          f"needle quad-vs-Laplace = {abs(float(Zn['lnZ_quad'])-float(Zn['lnZ_laplace'])):.3f} nat",
          flush=True)
    np.savez(f"{RES}/offset_ratio.npz", delta_needle=d_needle, delta_plateau=D,
             dlnZ_needle=dlnZ_needle, dlnZ_plateau=dlnZ_plateau,
             r_old=r_old, r_new=r_new, f_old=f_old, f_new=f_new,
             lnZn_old=lnZn_old, lnZp_old=lnZp_old, lnZp_err=float(Zres["lnZp_err"]))
    print(f"\n  saved {RES}/offset_ratio.npz", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "ratio"])
    ap.add_argument("--h", type=float, nargs="+", default=[-13.25, -12.5, -12.0])
    ap.add_argument("--chunk", type=int, default=16)
    ap.add_argument("--tol", type=float, default=1e-2)
    a = ap.parse_args()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    os.makedirs(RES, exist_ok=True)
    return mode_gate(a) if a.mode == "gate" else mode_ratio(a)


if __name__ == "__main__":
    sys.exit(main())
