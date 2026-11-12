"""fisher_nested_fast -- nested-N joint-array Fisher for the prong2 CW model, fast.

Model (prong2_transition, booked): per pulsar a and source s the residual is a pure
sinusoid at omega_s on a uniform TOA grid t_j = j*dt, j = 0..n-1:
    r_as(t) = Re[ a_r e^{i w t} ],  a_r = zeta w E (1 - P)
    w = B - iA,  A = Fp(1+c^2),  B = -2 Fx c,  E = e^{i phase0},  P = e^{-i w tau_as}
    d r/d phase0  -> i a_r ;   d r/d log10 h -> ln10 a_r
    d r/d L_a     -> a_L = zeta (i w) E P * w K (1 - cos mu)      (K = kpc/c)
Every Fisher entry is therefore a geometric factor times the Dirichlet kernel
    D(x) = sum_j e^{i x t_j} = e^{i x dt (n-1)/2} sin(n x dt/2)/sin(x dt/2)
    sum_t Re[a1 e^{i w1 t}] Re[a2 e^{i w2 t}] = 1/2 Re[a1 a2* D(w1-w2) + a1 a2 D(w1+w2)]
No timing-model fit, no frequency evolution, uniform sampling: the closed form is exact.

ESTIMATOR = booked prong2 'marg': sources' (phase0, log10 h) marginalised, OTHER pulsars'
distances held fixed, so F_ss is shared by every pulsar and
    I_marg,i(N) = I_cond,i(N) - b_i(N)^T F_ss(N)^{-1} b_i(N).
NESTING: sources ordered by index -> F_ss(N) is the leading 2N block of F_ss(N_max) and
its Cholesky factor is the leading block of chol(F_ss(N_max)). b_i(N) is NOT a prefix of
b_i(N_max): dr_i/dL_i sums over all N sources, so every entry of b_i changes with N (a
nested double sum, like I_cond). I_cond is therefore exact at every integer N; I_marg is
evaluated on a user grid (one Cholesky per sky, bucketed triangular solves).
Cholesky breakdown (numerically singular F_ss) is detected with LAPACK dpotrf's `info`
and reported -- never silently replaced by a pseudo-inverse.
"""
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from scipy.linalg import lapack
import sys
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import prong2_transition as P

LN10 = float(np.log(10.0))
KPC_OVER_C = P.KPC_OVER_C


def antenna_all(theta, phi, psi, p_hats):
    """(n_psr, n_src) Fp, Fx, cos_mu -- vectorised copy of P.antenna."""
    n = jnp.stack([jnp.sin(theta)*jnp.cos(phi), jnp.sin(theta)*jnp.sin(phi), jnp.cos(theta)], 1)
    et = jnp.stack([jnp.cos(theta)*jnp.cos(phi), jnp.cos(theta)*jnp.sin(phi), -jnp.sin(theta)], 1)
    ep = jnp.stack([-jnp.sin(phi), jnp.cos(phi), jnp.zeros_like(phi)], 1)
    cmu = p_hats @ n.T; dt = p_hats @ et.T; dp = p_hats @ ep.T
    den = jnp.clip(1.0 - cmu, 1e-6, None)
    Fp0 = 0.5*(dt*dt - dp*dp)/den; Fx0 = dt*dp/den
    c2, s2 = jnp.cos(2*psi)[None, :], jnp.sin(2*psi)[None, :]
    return Fp0*c2 + Fx0*s2, -Fp0*s2 + Fx0*c2, cmu


def amplitudes(srcs, p_hats, L_arr):
    """complex a_r, a_L of shape (n_psr, n_src); omega (n_src,)."""
    th, ph, ps, ci, l10f, l10h, ph0 = [srcs[:, k] for k in range(7)]
    Fp, Fx, cmu = antenna_all(th, ph, ps, p_hats)
    f = 10.0**l10f; om = 2*jnp.pi*f; zeta = 10.0**l10h/om
    w = (-2.0*ci[None, :]*Fx) - 1j*((1.0 + ci[None, :]**2)*Fp)
    tau = L_arr[:, None]*KPC_OVER_C*(1.0 - cmu)
    E = jnp.exp(1j*ph0)[None, :]; Pp = jnp.exp(-1j*om[None, :]*tau)
    a_r = zeta[None, :]*w*E*(1.0 - Pp)
    a_L = zeta[None, :]*(1j*w)*E*Pp*(om[None, :]*KPC_OVER_C*(1.0 - cmu))
    return a_r, a_L, om


def dirichlet(x, dt, n):
    h = 0.5*x*dt
    s = jnp.sin(h)
    small = jnp.abs(s) < 1e-13
    ratio = jnp.where(small, float(n), jnp.sin(n*h)/jnp.where(small, 1.0, s))
    return jnp.exp(1j*h*(n - 1))*ratio


def _assemble(srcs, p_hats, L_arr, dt, n, sigma, psr_chunk):
    """F_ss (2N,2N) interleaved [phase0_s, log10h_s]; cond_i(N) for every N; kernels."""
    a_r, a_L, om = amplitudes(srcs, p_hats, L_arr)
    Dm = dirichlet(om[:, None] - om[None, :], dt, n)
    Dp = dirichlet(om[:, None] + om[None, :], dt, n)
    Mm = a_r.T @ jnp.conj(a_r); Mp = a_r.T @ a_r
    Xm, Xp = Mm*Dm, Mp*Dp
    inv = 0.5/sigma**2
    Fpp = inv*jnp.real(Xm - Xp)
    Fph = inv*jnp.real(1j*LN10*(Xm + Xp))
    Fhh = inv*jnp.real(LN10**2*(Xm + Xp))
    N = srcs.shape[0]
    F = jnp.zeros((2*N, 2*N))
    F = F.at[0::2, 0::2].set(Fpp).at[0::2, 1::2].set(Fph).at[1::2, 0::2].set(Fph.T).at[1::2, 1::2].set(Fhh)
    conds = []
    for a0 in range(0, p_hats.shape[0], psr_chunk):
        aL = a_L[a0:a0+psr_chunk]
        T = inv*jnp.real(aL[:, :, None]*jnp.conj(aL)[:, None, :]*Dm[None] + aL[:, :, None]*aL[:, None, :]*Dp[None])
        R = jnp.cumsum(T, axis=2)
        v = jnp.diagonal(R, axis1=1, axis2=2)                 # sum_{s'<=s} T[s,s']
        dg = jnp.diagonal(T, axis1=1, axis2=2)
        conds.append(jnp.cumsum(2*v - dg, axis=1))             # (chunk, N): cond_i(N) for every N
    return F, jnp.concatenate(conds, 0), a_r, a_L, Dm, Dp


def _b_grid(a_r, a_L, Dm, Dp, grid, sigma, psr_chunk):
    """b_i(N_k) for grid N_k: (n_grid, 2N_max, n_psr)."""
    inv = 0.5/sigma**2
    idx = jnp.array(np.asarray(grid) - 1)
    out = []
    for a0 in range(0, a_r.shape[0], psr_chunk):
        aL = a_L[a0:a0+psr_chunk]; ar = a_r[a0:a0+psr_chunk]
        Qm = jnp.cumsum(Dm[None]*jnp.conj(aL)[:, None, :], axis=2)[:, :, idx]   # (c, s, g)
        Qp = jnp.cumsum(Dp[None]*aL[:, None, :], axis=2)[:, :, idx]
        Z = ar[:, :, None]*(Qm + Qp)
        bp = inv*jnp.real(1j*Z); bh = inv*jnp.real(LN10*Z)           # (c, s, g)
        b = jnp.stack([bp, bh], axis=2).reshape(aL.shape[0], -1, len(grid))   # (c, 2s, g)
        out.append(jnp.transpose(b, (2, 1, 0)))
    return jnp.concatenate(out, axis=2)


class NestedFisher:
    """Compile once per (N_max, n_psr, grid); call .run(srcs) per sky/arm."""

    def __init__(self, p_hats, T, cadence_days, sigma, N_max, grid, L_arr=None, psr_chunk=29,
                 buckets=(32, 128, 512, 1024, 1536, 2000)):
        self.p_hats = jnp.asarray(p_hats); self.n_psr = self.p_hats.shape[0]
        self.n = int(T/(cadence_days*86400.0)); self.dt = T/(self.n - 1)   # == np.linspace(0,T,n)
        self.sigma = sigma; self.N_max = N_max
        self.grid = np.asarray(sorted(set(int(g) for g in grid)))
        self.L_arr = jnp.ones(self.n_psr) if L_arr is None else jnp.asarray(L_arr)
        bk = sorted(set([min(b, 2*N_max) for b in buckets] + [2*N_max]))
        self.buckets = bk
        self._asm = jax.jit(lambda s: _assemble(s, self.p_hats, self.L_arr, self.dt, self.n, self.sigma, psr_chunk))
        self._bg = jax.jit(lambda ar, aL, Dm, Dp: _b_grid(ar, aL, Dm, Dp, self.grid, self.sigma, psr_chunk))
        # bucket plan: grid points whose 2N fits in bucket m (and not the previous one)
        self.plan = []
        prev = 0
        for m in bk:
            ks = [k for k, N in enumerate(self.grid) if prev < 2*N <= m]
            if ks:
                self.plan.append((m, np.array(ks)))
            prev = m
        def _solve(Lf, b, rows_ok):
            q = []
            for m, ks in self.plan:
                bb = b[ks][:, :m, :]                                   # (g, m, psr)
                g_, _, p_ = bb.shape
                rhs = jnp.transpose(bb, (1, 0, 2)).reshape(m, g_*p_)
                y = jax.scipy.linalg.solve_triangular(Lf[:m, :m], rhs, lower=True)
                y = jnp.transpose(y.reshape(m, g_, p_), (1, 0, 2))
                mask = (jnp.arange(m)[None, :] < 2*jnp.asarray(self.grid[ks])[:, None]).astype(y.dtype)
                q.append(jnp.sum(y**2*mask[:, :, None], axis=1))       # (g, psr)
            return jnp.concatenate(q, 0)
        self._solve = jax.jit(_solve, static_argnums=())

    def run(self, srcs):
        """Returns dict: cond (n_psr, N_max) exact at every N; marg (n_grid, n_psr) on the grid
        (NaN where Cholesky broke down); n_chol_ok = number of sources whose 2N block factored."""
        import time as _t
        t0 = _t.time()
        srcs = jnp.asarray(srcs)
        F, cond, a_r, a_L, Dm, Dp = self._asm(srcs)
        b = self._bg(a_r, a_L, Dm, Dp)
        b.block_until_ready(); t1 = _t.time()
        Lf = jnp.linalg.cholesky(F)
        n_ok = self.N_max
        if not bool(jnp.all(jnp.isfinite(jnp.diagonal(Lf)))):
            c, info = lapack.dpotrf(np.asarray(F), lower=1, clean=1)   # info = order of failing minor
            n_ok = (int(info) - 1)//2
            Lnp = np.tril(c); Lnp[2*n_ok:, :] = 0.0; Lnp[:, 2*n_ok:] = 0.0
            Lnp[np.arange(2*n_ok, 2*self.N_max), np.arange(2*n_ok, 2*self.N_max)] = 1.0
            Lf = jnp.asarray(Lnp)
        q = self._solve(Lf, b, None)
        q.block_until_ready(); t2 = _t.time()
        self.t_asm = getattr(self, "t_asm", 0.0) + (t1 - t0)
        self.t_fac = getattr(self, "t_fac", 0.0) + (t2 - t1)
        condg = cond[:, self.grid - 1].T                                # (g, psr)
        marg = jnp.clip(condg - q, 0.0, None)
        marg = np.array(marg)
        marg[self.grid > n_ok, :] = np.nan
        return dict(cond=np.array(cond), marg=marg, n_chol_ok=n_ok, grid=self.grid)
