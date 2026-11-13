"""JOB 54 core: evolving-CW residual, array, sources, closed-form sky likelihood.

Signal per pulsar (monochromatic within the span, evolution between Earth and pulsar term):
  r(t) = zeta_e * W(Phi_e) - zeta_p * W(Phi_p),  W(Phi) = Fp*(1+ci^2)*sin(Phi) + Fx*(-2ci)*cos(Phi)
  Phi_e = 2 pi f t + phi0 ;  Phi_p = 2 pi f_p t + phi0 - dPhi
  tau = L(1-cos mu)/c ; x = (256/5) (G Mc/c^3)^(5/3) (pi f)^(8/3) tau
  f_p = f (1+x)^(-3/8) ; dPhi = 2 pi f tau * [(1+x)^(5/8)-1]/(5x/8) ; zeta_p = zeta_e (f/f_p)^(1/3)
Mc -> 0 gives prong2_transition.single_source_residual exactly.
"""
import os, sys
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")   # shared GPU: do not grab all memory
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import prong2_transition as p2

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "job54_data")
GMSUN_C3 = 4.925490947e-6            # G Msun / c^3 [s]
PC_OVER_C = p2.KPC_OVER_C / 1000.0   # s per pc
YR = p2.YR
N_PSR = 116
SIGMA = 1e-7
T_SPAN = 15 * YR
CADENCE_D = 7.0
NPH_LIST = [0, 1, 2, 3, 4, 6, 10, 20, 40, 80, 116]
SIGL_PHASED = 0.01
DEG2 = (180.0 / np.pi) ** 2


def build_all(seed_dist=54, seed_perm=5401):
    """Array: prong2 build_array(116, seed=777), toas, census (D, sigma) pairs with replacement."""
    p_hats = np.asarray(p2.build_array(N_PSR))            # seed 777 default
    toas = np.asarray(p2.build_toas(T_SPAN, CADENCE_D))
    import csv
    rows = list(csv.DictReader(open(os.path.join(HERE, "job45_t2_parallax_census.csv"))))
    D = np.array([float(r["D_kpc"]) for r in rows]) * 1000.0
    S = np.array([float(r["sigma_d_now_pc"]) for r in rows])
    names = np.array([r["pulsar"] for r in rows])
    rng = np.random.default_rng(seed_dist)
    idx = rng.integers(0, len(D), N_PSR)
    L_pc, sigL_pc = D[idx], S[idx]
    order_near = np.argsort(L_pc, kind="stable")
    order_rand = np.random.default_rng(seed_perm).permutation(N_PSR)
    return dict(p_hats=p_hats, toas=toas, L_pc=L_pc, sigL_pc=sigL_pc, census_idx=idx,
                census_names=names[idx], order_near=order_near, order_rand=order_rand)


def draw_sources(n=200, seed=5400):
    rng = np.random.default_rng(seed)
    theta = np.arccos(rng.uniform(-1, 1, n))
    phi = rng.uniform(0, 2 * np.pi, n)
    psi = rng.uniform(0, np.pi, n)
    cinc = rng.uniform(-1, 1, n)
    phase0 = rng.uniform(0, 2 * np.pi, n)
    return dict(theta=theta, phi=phi, psi=psi, cinc=cinc, phase0=phase0)


def evol(f, mc_s, tau):
    """Return (f_p/f, dPhi) for Earth GW freq f, chirp mass in seconds, retardation tau."""
    w0 = jnp.pi * f
    x = (256.0 / 5.0) * mc_s ** (5.0 / 3.0) * w0 ** (8.0 / 3.0) * tau
    l1p = jnp.log1p(x)
    ok = x > 1e-200
    xs = jnp.where(ok, x, 1.0)
    ratio = jnp.where(ok, jnp.expm1(0.625 * l1p) / (0.625 * xs), 1.0)
    return jnp.exp(-0.375 * l1p), 2.0 * w0 * tau * ratio


def residual(toas, p_hat, L_pc, theta, phi, psi, cinc, f, h, phase0, mc_msun):
    """Evolving CW residual at one pulsar [s]."""
    Fp, Fx, cos_mu = p2.antenna(theta, phi, psi, p_hat)
    tau = L_pc * PC_OVER_C * (1.0 - cos_mu)
    fr, dPhi = evol(f, mc_msun * GMSUN_C3, tau)
    fp = f * fr
    zeta = h / (2.0 * jnp.pi * f)
    zeta_p = zeta * fr ** (-1.0 / 3.0)
    wp, wx = 1.0 + cinc ** 2, -2.0 * cinc
    Pe = 2.0 * jnp.pi * f * toas + phase0
    Pp = 2.0 * jnp.pi * fp * toas + phase0 - dPhi
    return (zeta * (Fp * wp * jnp.sin(Pe) + Fx * wx * jnp.cos(Pe))
            - zeta_p * (Fp * wp * jnp.sin(Pp) + Fx * wx * jnp.cos(Pp)))


_res_arr = jax.jit(jax.vmap(residual, in_axes=(None, 0, 0) + (None,) * 8))


def snr_unit_h(arr, src_i, f, mc, sigma=SIGMA):
    """h giving total array SNR = 1 (Earth + pulsar terms, white noise)."""
    r = _res_arr(jnp.asarray(arr["toas"]), jnp.asarray(arr["p_hats"]), jnp.asarray(arr["L_pc"]),
                 src_i["theta"], src_i["phi"], src_i["psi"], src_i["cinc"], f, 1.0,
                 src_i["phase0"], mc)
    return float(sigma / jnp.sqrt(jnp.sum(r ** 2)))


def src_at(src, i):
    return {k: float(v[i]) for k, v in src.items()}


# ---------------------------------------------------------------- closed-form sky likelihood
def frame(theta0, phi0):
    n0 = np.array([np.sin(theta0) * np.cos(phi0), np.sin(theta0) * np.sin(phi0), np.cos(theta0)])
    ex = np.array([-np.sin(phi0), np.cos(phi0), 0.0])                       # east
    ey = -np.array([np.cos(theta0) * np.cos(phi0), np.cos(theta0) * np.sin(phi0), -np.sin(theta0)])  # north
    return n0, ex, ey


def sky_par(arr, s, f, h, mc, sigma=SIGMA):
    """Pack traced parameters for parts_fn (so one compiled function serves every source)."""
    n0, ex, ey = frame(s["theta"], s["phi"])
    toas = arr["toas"]
    return dict(n0=jnp.asarray(n0), ex=jnp.asarray(ex), ey=jnp.asarray(ey), we=2.0 * np.pi * f,
                zeta=h / (2.0 * np.pi * f), cinc=s["cinc"], psi=s["psi"], ph0=s["phase0"],
                mcs=mc * GMSUN_C3, P=jnp.asarray(arr["p_hats"]), Lsec=jnp.asarray(arr["L_pc"]) * PC_OVER_C,
                s2=sigma ** 2, N=float(len(toas)), dt=float(toas[1] - toas[0]))


def uv2n(u, v, par):
    return (jnp.cos(v) * jnp.cos(u))[:, None] * par["n0"] + (jnp.cos(v) * jnp.sin(u))[:, None] * par["ex"] \
        + jnp.sin(v)[:, None] * par["ey"]


def parts_fn(u, v, par):
    """Closed-form per-pulsar log-likelihood pieces at rotated-frame coords (u, v) [rad] about the truth.

    Returns coh[C,P], inc[C,P], a[C,P], alpha[C,P], sphi1[C,P], base[C,P]
      coh  : coherent ln L at the true distance (exactly 0 at truth)
      inc  : ln L with the pulsar-term phase marginalised (uniform): base + ln I0(a)
      a, alpha : cross term (d-E).P/sigma^2 = a cos(alpha) ; sphi1 = d(phase)/dL [rad/pc]
    """
    N, dt, we, s2, P = par["N"], par["dt"], par["we"], par["s2"], par["P"]
    wp_i, wx_i = 1.0 + par["cinc"] ** 2, -2.0 * par["cinc"]
    psi, ph0 = par["psi"], par["ph0"]

    def Dk(w):
        x = 0.5 * w * dt
        sx = jnp.sin(x)
        small = jnp.abs(sx) < 1e-10
        ratio = jnp.where(small, N, jnp.sin(N * x) / jnp.where(small, 1.0, sx))
        return ratio * jnp.exp(1j * (N - 1) * x)

    def terms(nv):
        nz = jnp.clip(nv[:, 2], -1.0, 1.0)
        th = jnp.arccos(nz); ph = jnp.arctan2(nv[:, 1], nv[:, 0])
        st, ct, sp, cp = jnp.sin(th), jnp.cos(th), jnp.sin(ph), jnp.cos(ph)
        eth = jnp.stack([ct * cp, ct * sp, -st], 1); eph = jnp.stack([-sp, cp, jnp.zeros_like(sp)], 1)
        cmu = nv @ P.T; d1 = eth @ P.T; d2 = eph @ P.T
        den = jnp.clip(1.0 - cmu, 1e-6, None)
        Fp0 = 0.5 * (d1 * d1 - d2 * d2) / den; Fx0 = d1 * d2 / den
        Fp = Fp0 * jnp.cos(2 * psi) + Fx0 * jnp.sin(2 * psi)
        Fx = -Fp0 * jnp.sin(2 * psi) + Fx0 * jnp.cos(2 * psi)
        g = Fx * wx_i - 1j * Fp * wp_i
        tau = par["Lsec"][None, :] * (1.0 - cmu)
        fr, dPhi = evol(we / (2.0 * jnp.pi), par["mcs"], tau)
        wP = we * fr
        cE = par["zeta"] * g * jnp.exp(1j * ph0)
        cP = -par["zeta"] * fr ** (-1.0 / 3.0) * g * jnp.exp(1j * (ph0 - dPhi))
        return cE, cP, wP, cmu

    cE0, cP0, wP0, _ = terms(par["n0"][None, :])
    cE0, cP0, wP0 = cE0[0], cP0[0], wP0[0]
    D2we = Dk(2 * we)
    Dem, Dep = Dk(we - wP0), Dk(we + wP0)
    PP0 = 0.5 * jnp.abs(cP0) ** 2 * N + 0.5 * jnp.real(cP0 ** 2 * Dk(2 * wP0))
    cE, cP, wP, cmu = terms(uv2n(u, v, par))
    dc = cE0 - cE
    qdd = (0.5 * jnp.abs(dc) ** 2 * N + 0.5 * jnp.real(dc ** 2 * D2we)
           + jnp.real(dc * jnp.conj(cP0) * Dem) + jnp.real(dc * cP0 * Dep) + PP0)
    z0 = 0.5 * (dc * Dk(we + wP) + jnp.conj(dc) * Dk(wP - we)) \
        + 0.5 * (cP0 * Dk(wP0 + wP) + jnp.conj(cP0) * Dk(wP - wP0))
    cz = cP * z0
    nP = 0.5 * jnp.abs(cP) ** 2 * N
    n2 = 0.5 * jnp.real(cP ** 2 * Dk(2 * wP))
    coh = -(qdd - 2 * jnp.real(cz) + nP + n2) / (2 * s2)
    a = jnp.abs(cz) / s2
    base = -(qdd + nP) / (2 * s2)
    inc = base + jnp.log(jax.scipy.special.i0e(a)) + a
    return coh, inc, a, jnp.angle(cz), wP * PC_OVER_C * (1.0 - cmu), base


def make_sky(arr, s, f, h, mc, sigma=SIGMA):
    par = sky_par(arr, s, f, h, mc, sigma)
    return dict(parts=lambda u, v: parts_fn(u, v, par), uv2n=lambda u, v: uv2n(u, v, par), par=par)


def partial_lnl(base, a, alpha, sphi, M=40):
    """ln E_{dphi~N(0,sphi^2)} exp(a cos(alpha+dphi)) + base, via Jacobi-Anger series
    sum_m eps_m I_m(a) exp(-m^2 sphi^2/2) cos(m alpha); I_m by Miller backward recurrence."""
    a_s = jnp.clip(a, 1e-2, None)
    Ip1 = jnp.zeros_like(a_s); Ic = jnp.full_like(a_s, 1e-250)
    norm = jnp.zeros_like(a_s); acc = jnp.zeros_like(a_s)
    for m in range(M, 0, -1):          # Ic = I_m
        norm = norm + 2 * Ic
        acc = acc + 2 * Ic * jnp.exp(-0.5 * (m * sphi) ** 2) * jnp.cos(m * alpha)
        Ic, Ip1 = (2.0 * m / a_s) * Ic + Ip1, Ic
    norm = norm + Ic; acc = acc + Ic
    val = jnp.clip(acc / norm, 1e-300, None)       # = E[...] * exp(-a)
    return base + a_s + jnp.log(val)
