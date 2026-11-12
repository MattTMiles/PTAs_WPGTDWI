"""LOTTERY-2 -- the REALISTIC LOUDNESS PRIOR (agent LOTTERY-2, cronus, CPU only).

THE QUESTION. GLACIER Stage-2a banked `P_cond(onset) = 0 EXACTLY`: under the campaign's
fixed-N power-conserving draw the brightest member of an NG15-normalised sky cannot exceed
log10 h = -13.867, because the draw PINS  sum_i h_i^2 == A^2 (3/4)[x(flo)-x(fhi)]  with
A = -14.6 EXACTLY and N = 256 EXACTLY. That zero is a property of the DRAW, not of nature:
three tails are suppressed by construction.

  (T1) N is fixed at 256. A real band population is a Poisson realisation.
  (T2) The strains are a bare lognormal SHAPE, renormalised. Masses/redshifts are drawn and
       then thrown away -- the loudness distribution carries no synthesis information, and
       in particular no  N(>h) ~ h^-3  Euclidean volume tail.
  (T3) A is pinned at the NG15 point estimate. NG15's own posterior has width.
  (T-pin) And the deep one, the mechanism behind the hard ceiling itself: the REALISED sky
       power is constrained to equal the ENSEMBLE mean. One source carrying 100 % of the
       band power is the ceiling only because the total is a constant.

WHAT THIS MODULE BUILDS
-----------------------
(A) FIXED-CLASS FAMILY (M0..M3) -- the incumbent draw with the three tails switched on one
    at a time, so the ablation says WHICH suppression was load-bearing. M0 reproduces
    reports/gen_stage2_pcond_scan.npz (gate g0).
(B) REALISTIC FAMILY (M4) -- a standard phenomenological SMBHB synthesis:
      d^2 n / (dlog10 Mc dz) = ndot0 (Mc/1e7)^(-alpha) exp(-Mc/Mstar) (1+z)^beta exp(-z/z0)
    [VERIFY: Middleton+2016 (MNRAS 455 L72) / Chen, Sesana & Conselice 2019 (MNRAS 488
     401) / the NG15 astro-interpretation Agazie+2023 ApJL 952 L37 use this phenomenological
     form; the fiducial (alpha, Mstar, beta, z0) below are DECLARED, and swept, not fitted.]
    Sources in the observed band, per unit lnf, weighted by the GW-driven residence time:
      dN/(dlog10Mc dz dlnf) ~ Phi(Mc,z) (dV_c/dz) (1+z)^-2 Mz^(-5/3) f^(-8/3)
    strain (sky/inclination/polarisation averaged, the h_c bookkeeping convention):
      h = (8/sqrt(10)) (G Mz)^(5/3) (pi f)^(2/3) / (c^4 d_L)      [VERIFY: Sesana, Vecchio
      & Colacino 2008 MNRAS 390 192 eq. (3)-(6); this is the h whose quadrature sum over
      the band IS int h_c^2 dlnf, i.e. the SAME h the campaign pins.]
    The rate normalisation ndot0 is NOT an input: it is fixed by requiring the ENSEMBLE
    band power to equal the drawn NG15 target,  lambda * <h^2> = P_target(A).  (ndot0 and A
    are degenerate; the field calibrates the synthesis to the measured GWB, and so do we.)

THE POISSON-PROCESS IDENTITY (why the tail is exact, not sampled)
-----------------------------------------------------------------
Band sources are a Poisson process, so the number above any H is Poisson with mean
mu(H) = lambda * P(h >= H), and
        P(brightest >= H) = 1 - exp(-mu(H))          [exact, no MC noise, any depth]
mu(H) is computed by direct quadrature on the (log10 Mc, z) grid with the f-integral in
closed form -- so 1e-8 tails are as accurate as 1e-1 ones. The >= 1e5-draw MC is run beside
it as a cross-check and to carry the joints and the power-conditioning variants.

THREE CONDITIONING VARIANTS (the NG15-consistency question, stated three ways)
-----------------------------------------------------------------------------
  C1 ENSEMBLE-CALIBRATED : <sum h^2> = P_target(A). The realised sky fluctuates. (Default.)
  C2 TOTAL-POWER-PINNED  : realised sum_i h_i^2 must itself sit inside the NG15 posterior.
                           The incumbent draw is the sigma_A -> 0 limit of this.
  C3 BACKGROUND-PINNED   : the UNRESOLVED sum (all but the brightest) sits inside the NG15
                           posterior. This is the defensible one: a single resolvable source
                           is not HD-correlated, so it is not what NG15's cross-correlation
                           amplitude measures.

Seeds: 21_2xx_xxx (a fresh declared block, disjoint from GLACIER 9_0xx/9_5xx and every
LOTTERY range). Matt commits; nothing is committed from here. CPU only, no GPU touched.
"""
import argparse
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
OUT_BANK = f"{HSYMT}/reports/lottery2_realistic_pop.npz"

# ---------------- physical constants (SI) ----------------
G = 6.67430e-11
C_L = 2.99792458e8
MSUN = 1.98892e30
MPC = 3.0856775814913673e22
YR = 365.25 * 86400.0
F_YR = 1.0 / YR

# ---------------- the campaign's declared numbers (carried verbatim) ----------------
A_TARGET_LOG10 = -14.6          # glacier_pop.A_TARGET_LOG10 (repo NG15 convention)
BAND_SCAN = (-8.7, -7.5)        # the band the banked P_cond scan used (gen_stage2_pcond_scan)
BAND_GEN = (-8.0, -7.5)         # glacier_pop.LOG10_FGW_RANGE (the estimator population box)
N_POP_DEFAULT = 256             # glacier_pop.N_POP_DEFAULT
SIGMA_LOGH = 0.6                # glacier_pop.SIGMA_LOGH (lognormal strain-hierarchy width)
SIGMA_NG15_REPO = 0.05          # glacier_pop.SIGMA_NG15_DEX -- the campaign's DECLARED width
SIGMA_NG15_LIT = 0.11           # [VERIFY] NG15 gamma=13/3: A = 2.4 (+0.7/-0.6) e-15
                                # (Agazie+2023 ApJL 951 L8) -> log10 A = -14.62 +0.11/-0.12.
                                # The repo's declared 0.05 is ~2.2x NARROWER than NG15's own
                                # posterior; both are carried.
RUNGS = (-13.87, -13.5, -13.25)   # the GLACIER Stage-2a ladder, exact
RUNG_TAIL_REALISED = -13.9        # rung 1 as actually realised (tail threshold)
ROBUST_H = -12.75                 # A-SKY's surviving line: 5+11, e>=0.3, h >= -12.75

SEED_BASE = 21_200_000

# ---------------- the phenomenological synthesis (DECLARED, swept) ----------------
# [VERIFY] all four fiducials against Middleton+2016 / Chen+2019 / Agazie+2023 (952 L37).
FID = dict(alpha=0.0,          # dn/dlog10 M slope below the cutoff
           log10_Mstar=8.0,    # exponential mass cutoff, Msun -- TOTAL mass by default
           beta=2.0,           # (1+z)^beta
           z0=2.0,             # exp(-z/z0)
           q_min=0.25,         # mass-ratio range [q_min, 1]
           q_slope=0.0)        # p(q) ~ q^q_slope
# [VERIFY -- THE CONVENTION THAT MOVES THE ANSWER MOST] the phenomenological form is normally
# written in TOTAL mass with Mstar ~ 1e8 Msun, but the strain needs the CHIRP mass. Phi is
# therefore defined in total mass and PUSHED FORWARD to chirp mass through a declared p(q)
# (kappa(q) = q^(3/5)/(1+q)^(6/5); kappa(1) = 0.4353, so an equal-mass 1e8 total is a 10^7.64
# chirp). Reading Mstar as a CHIRP mass instead moves every rung probability by orders of
# magnitude -- which is why (alpha, Mstar, chirp-vs-total) is SWEPT and the sweep, not the
# fiducial, is the deliverable.
MC_RANGE = (5.5, 11.0)          # log10 chirp mass integration range, Msun
D_MIN_MPC = (0.0, 1.0, 5.0, 10.0, 20.0, 50.0)   # near-field validity floors (SS(a2), g4)
Z_RANGE = (1e-7, 5.0)           # log-spaced. THE LOWER EDGE IS LOAD-BEARING and was moved here
                                # from 1e-4 after gate g3 caught it: the brightest-member tail is
                                # dominated by the NEAR FIELD (d_L ~ 0.4-2 Mpc; g4 banks the
                                # quantiles), so truncating z from below truncates the answer --
                                # mu(-13.25) is 2.2x LOW at z_lo = 1e-4. It is converged to <0.1 %
                                # by z_lo = 1e-6 and to <0.01 % by 1e-7 (mu contribution ~ z^3).
H0_KMS = 70.0
OMEGA_M = 0.3


# ============================================================
# cosmology (flat LCDM, quadrature -- no external deps)
# ============================================================
def _cosmo(z_grid, h0_kms=H0_KMS, om=OMEGA_M):
    """Return (d_L [m], dV_c/dz [m^3]) on z_grid. Flat LCDM."""
    H0 = h0_kms * 1e3 / MPC                      # s^-1
    D_H = C_L / H0                               # m
    zz = np.linspace(0.0, z_grid.max(), 20001)
    E = np.sqrt(om * (1 + zz) ** 3 + (1 - om))
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (1.0 / E[1:] + 1.0 / E[:-1]) * np.diff(zz))])
    D_C = D_H * np.interp(z_grid, zz, integ)     # comoving distance, m
    d_L = (1.0 + z_grid) * D_C
    E_z = np.sqrt(om * (1 + z_grid) ** 3 + (1 - om))
    dVdz = 4.0 * np.pi * D_C ** 2 * D_H / E_z    # comoving volume element, m^3
    return d_L, dVdz


# ============================================================
# the band power target -- glacier_pop._band_power_target verbatim
# ============================================================
def band_power_target(a_log10=A_TARGET_LOG10, band=BAND_SCAN):
    A2 = 10.0 ** (2 * a_log10)
    flo, fhi = 10.0 ** np.array(band)
    x = lambda f: (f / F_YR) ** (-4.0 / 3.0)
    return A2 * 0.75 * (x(flo) - x(fhi))


def ceiling_log10h(a_log10=A_TARGET_LOG10, band=BAND_SCAN):
    """The fixed-class HARD CEILING: one source carrying 100 % of the pinned band power."""
    return 0.5 * np.log10(band_power_target(a_log10, band))


# ============================================================
# (A) FIXED-CLASS FAMILY -- the incumbent draw, tails switched on one at a time
# ============================================================
def fixed_class_hmax(n_draw, model, band=BAND_SCAN, n_src=N_POP_DEFAULT,
                     sigma_logh=SIGMA_LOGH, sigma_A=SIGMA_NG15_LIT, seed=SEED_BASE,
                     chunk=20000):
    """Brightest-member log10 h under the incumbent power-conserving draw.

    model: 'M0' incumbent (N fixed, A pinned)
           'M1' + N ~ Poisson(256)
           'M2' + A ~ N(-14.6, sigma_A)          (N still fixed)
           'M3' M1 + M2 together (still renormalised -- the PIN is still on)
    The pin is what makes the ceiling; M0..M3 all keep it, by design.
    """
    rng = np.random.default_rng(seed)
    out = np.empty(n_draw)
    done = 0
    while done < n_draw:
        m = min(chunk, n_draw - done)
        pois = model in ("M1", "M3")
        awid = model in ("M2", "M3")
        a = (A_TARGET_LOG10 + rng.normal(0.0, sigma_A, size=m)) if awid \
            else np.full(m, A_TARGET_LOG10)
        target = band_power_target(0.0, band) * 10.0 ** (2 * a)   # scales as A^2
        if pois:
            nn = rng.poisson(n_src, size=m)
            nn = np.maximum(nn, 1)
        else:
            nn = np.full(m, n_src)
        nmax = int(nn.max())
        raw = 10.0 ** rng.normal(0.0, sigma_logh, size=(m, nmax))
        mask = np.arange(nmax)[None, :] < nn[:, None]
        raw = np.where(mask, raw, 0.0)
        tot = np.sum(raw ** 2, axis=1)
        scale = np.sqrt(target / tot)
        out[done:done + m] = np.log10(raw.max(axis=1) * scale)
        done += m
    return out


def power_share_sample(model, n=1_000_000, n_src=N_POP_DEFAULT, sigma_logh=SIGMA_LOGH,
                       seed=SEED_BASE + 3, chunk=50000):
    """s = log10( max_i h_i / sqrt(sum_i h_i^2) ) for the incumbent draw -- the BRIGHTEST
    MEMBER'S POWER SHARE, in dex, and the whole content of the fixed class:

        log10 h_max  =  ceiling(A0) + (A - A0) + s ,     s <= 0 identically.

    s <= 0 IS the hard ceiling (one source can carry at most 100 % of a pinned total). Sampling
    s once lets the fixed-class tail be evaluated EXACTLY at any depth, instead of counting MC
    hits that are all zero: a rung 0.62 dex above the ceiling needs a 5.6-sigma A excursion, and
    no affordable draw-level MC resolves that. Depends on the model only through N."""
    rng = np.random.default_rng(seed)
    out = np.empty(n)
    done = 0
    while done < n:
        m = min(chunk, n - done)
        nn = np.maximum(rng.poisson(n_src, size=m), 1) if model in ("M1", "M3") \
            else np.full(m, n_src)
        nmax = int(nn.max())
        raw = 10.0 ** rng.normal(0.0, sigma_logh, size=(m, nmax))
        raw = np.where(np.arange(nmax)[None, :] < nn[:, None], raw, 0.0)
        out[done:done + m] = np.log10(raw.max(axis=1) / np.sqrt(np.sum(raw ** 2, axis=1)))
        done += m
    return out


def _gauss_upper(x):
    """Q(x) = P(Z >= x), accurate into the deep tail (erfc, not 1 - Phi)."""
    from math import erfc
    xa = np.atleast_1d(np.asarray(x, dtype=float))
    return np.array([0.5 * erfc(v / np.sqrt(2.0)) for v in xa]).reshape(np.shape(x))


def fixed_class_p_above(rung, model, share, sigma_A, band=BAND_SCAN,
                        a_log10=A_TARGET_LOG10):
    """P(brightest >= rung) in the fixed class, EXACT in the A-tail.

    sigma_A = 0 (M0/M1): P = mean(s >= rung - ceiling) -- identically 0 above the ceiling.
    sigma_A > 0 (M2/M3): P = E_s[ Q( (rung - ceiling - s) / sigma_A ) ].
    """
    ceil0 = ceiling_log10h(a_log10, band)
    need = float(rung) - ceil0
    if sigma_A <= 0:
        return float(np.mean(share >= need))
    return float(np.mean(_gauss_upper((need - share) / sigma_A)))


def sigma_above_ceiling(rung, sigma_A, band=BAND_SCAN, a_log10=A_TARGET_LOG10):
    """How many sigma of the NG15 amplitude posterior the rung sits above the fixed-class
    ceiling -- the GLACIER Stage-2a 'tension' axis, recomputed here."""
    return (float(rung) - ceiling_log10h(a_log10, band)) / sigma_A


# ============================================================
# (B) REALISTIC FAMILY -- the synthesis-based loudness prior
# ============================================================
class RealisticPrior:
    """The per-source loudness distribution of an NG15-calibrated band population, and the
    Poisson-process tail it implies. All heavy lifting is quadrature on a (log10Mc, z) grid
    with the f-integral in closed form."""

    def __init__(self, band=BAND_SCAN, alpha=FID["alpha"], log10_Mstar=FID["log10_Mstar"],
                 beta=FID["beta"], z0=FID["z0"], n_m=260, n_z=340,
                 mc_range=MC_RANGE, z_range=Z_RANGE,
                 q_min=FID["q_min"], q_slope=FID["q_slope"], mass_is_chirp=False):
        self.band = tuple(band)
        self.pars = dict(alpha=alpha, log10_Mstar=log10_Mstar, beta=beta, z0=z0,
                         q_min=q_min, q_slope=q_slope, mass_is_chirp=mass_is_chirp)
        self.flo, self.fhi = 10.0 ** np.array(band)

        # LOG-SPACED in z. The deep tail (the loud sources) is dominated by z ~ 0.1-1, and
        # a linear z grid over [0, 5] puts only ~10 nodes there -- g3 caught exactly that
        # (mu moved 50 % under grid refinement on the first, linear-grid, build).
        lm = np.linspace(*mc_range, n_m)                  # log10 Mc / Msun
        z = np.logspace(np.log10(z_range[0]), np.log10(z_range[1]), n_z)
        self.lm, self.z = lm, z
        LM, Z = np.meshgrid(lm, z, indexing="ij")
        Mc = 10.0 ** LM                                    # Msun
        Mz = Mc * (1.0 + Z)                                # redshifted chirp mass

        d_L, dVdz = _cosmo(z)
        d_L = d_L[None, :] + 0.0 * LM
        dVdz = dVdz[None, :] + 0.0 * LM

        # -- phenomenological merger-rate density (shape only; ndot0 absorbed into lambda).
        #    Phi is written in TOTAL mass and PUSHED FORWARD to chirp mass through p(q):
        #    log10 Mc = log10 M + log10 kappa(q), kappa = q^(3/5)/(1+q)^(6/5). At fixed q the
        #    Jacobian is unity in log-space, so the pushforward is a plain 1-D mixture.
        f_M = lambda M: (M / 1e7) ** (-alpha) * np.exp(-M / 10.0 ** log10_Mstar)
        if mass_is_chirp:
            phi_m = f_M(10.0 ** lm)
            self.q_grid, self.q_w = np.array([1.0]), np.array([1.0])
        else:
            qg = np.linspace(q_min, 1.0, 64)
            wq = qg ** q_slope
            wq = wq / wq.sum()
            lgk = np.log10(qg ** 0.6 / (1.0 + qg) ** 1.2)
            phi_m = np.sum(wq[:, None] * f_M(10.0 ** (lm[None, :] - lgk[:, None])), axis=0)
            self.q_grid, self.q_w = qg, wq
        phi = phi_m[:, None] * ((1.0 + Z) ** beta * np.exp(-Z / z0))

        # -- number weight in the band: residence time * volume * time dilation
        #    dt_r/dlnf_r ~ Mz^(-5/3) f^(-8/3) (1+z)^-1 ; observer-frame rate -> (1+z)^-1
        w = phi * dVdz * (1.0 + Z) ** (-2.0) * (Mz * MSUN) ** (-5.0 / 3.0)
        # trapezoid node weights -- the z grid is NON-UNIFORM (log-spaced), so the cell
        # measure cannot be a scalar dz.
        dlm = lm[1] - lm[0]
        wz = np.empty_like(z)
        wz[1:-1] = 0.5 * (z[2:] - z[:-2])
        wz[0] = 0.5 * (z[1] - z[0])
        wz[-1] = 0.5 * (z[-1] - z[-2])
        self.w = w * dlm * wz[None, :]
        self.w_sum = float(self.w.sum())

        # -- strain prefactor: h = K f^(2/3), K = (8/sqrt10)(G Mz)^(5/3) pi^(2/3)/(c^4 d_L)
        self.K = (8.0 / np.sqrt(10.0)) * (G * Mz * MSUN) ** (5.0 / 3.0) \
            * np.pi ** (2.0 / 3.0) / (C_L ** 4 * d_L)
        self.d_mpc = d_L / MPC

        # -- f-marginal moments under p(f) ~ f^(-11/3) on the band
        self.Ef43 = self._f_moment(4.0 / 3.0)
        self.mean_h2 = float(np.sum(self.w * self.K ** 2) * self.Ef43 / self.w_sum)

    # ---- f helpers: p(f) ~ f^(-11/3) on [flo, fhi] (dN/df of a GW-driven inspiral) ----
    def _f_norm(self):
        return (3.0 / 8.0) * (self.flo ** (-8.0 / 3.0) - self.fhi ** (-8.0 / 3.0))

    def _f_moment(self, p):
        """E[f^p] under p(f) ~ f^(-11/3)."""
        q = p - 11.0 / 3.0 + 1.0
        return (self.flo ** q - self.fhi ** q) / (-q) / self._f_norm() if abs(q) > 1e-12 \
            else np.log(self.fhi / self.flo) / self._f_norm()

    def _f_surv(self, f_thr):
        """P(f >= f_thr) under p(f) ~ f^(-11/3), clipped to the band."""
        ft = np.clip(f_thr, self.flo, self.fhi)
        return (3.0 / 8.0) * (ft ** (-8.0 / 3.0) - self.fhi ** (-8.0 / 3.0)) / self._f_norm()

    # ---- the two headline quantities ----
    def lam(self, a_log10=A_TARGET_LOG10):
        """Expected number of band sources: fixes ndot0 from the drawn GWB amplitude."""
        return band_power_target(a_log10, self.band) / self.mean_h2

    def mu(self, log10_h, a_log10=A_TARGET_LOG10, d_min_mpc=0.0, chunk=64):
        """Expected NUMBER of band sources at or above log10_h (Poisson-process mean).

        d_min_mpc: a DECLARED near-field validity floor. The tail is dominated by the Local
        Volume (see tail_provenance / gate g4), where a smooth comoving rate density is an
        extrapolation, not a measurement. mu(H | no source closer than d_min) is the honest
        way to quote it, so the floor is a first-class argument, never a hidden default."""
        H = np.atleast_1d(10.0 ** np.asarray(log10_h, dtype=float))
        wmask = self.w if d_min_mpc <= 0 else np.where(self.d_mpc >= d_min_mpc, self.w, 0.0)
        out = np.empty(H.size)
        for i in range(0, H.size, chunk):                 # chunked: the grid is (n_m, n_z)
            hs = H[i:i + chunk]
            s = self._f_surv((hs[:, None, None] / self.K[None, :, :]) ** 1.5)
            out[i:i + chunk] = (wmask[None, :, :] * s).sum(axis=(1, 2)) / self.w_sum
        out *= self.lam(a_log10)
        return out if np.ndim(log10_h) else float(out[0])

    def p_brightest_above(self, log10_h, a_log10=A_TARGET_LOG10):
        return 1.0 - np.exp(-self.mu(log10_h, a_log10))

    def tail_provenance(self, log10_h, qs=(0.05, 0.5, 0.95)):
        """WHERE the sources above log10_h live: quantiles of d_L [Mpc] and log10 Mc of the
        mu-integrand. Banked because it is the honest scope of the number: if the tail is
        dominated by d_L of a few Mpc, the answer is a homogeneous-population extrapolation
        into the Local Volume, and must be quoted as such."""
        f_thr = (10.0 ** float(log10_h) / self.K) ** 1.5
        c = self.w * self._f_surv(f_thr)
        tot = float(c.sum())
        if tot <= 0:
            return dict(mu_frac=0.0, z=[np.nan] * 3, d_mpc=[np.nan] * 3, log10_mc=[np.nan] * 3)
        cz = np.cumsum(c.sum(axis=0)) / tot
        cm = np.cumsum(c.sum(axis=1)) / tot
        zq = np.array([self.z[min(np.searchsorted(cz, q), len(self.z) - 1)] for q in qs])
        mq = [float(self.lm[min(np.searchsorted(cm, q), len(self.lm) - 1)]) for q in qs]
        d_l, _ = _cosmo(zq)
        return dict(mu_frac=tot / self.w_sum, z=[float(v) for v in zq],
                    d_mpc=[float(v / MPC) for v in d_l], log10_mc=mq)

    def p_n_above(self, log10_h, k, a_log10=A_TARGET_LOG10):
        """P(N sources above log10_h >= k) -- the structure proxy (5+11 <-> k=5).
        Returned as log10 P as well, because in the deep tail P underflows float64."""
        from math import lgamma
        m = float(self.mu(log10_h, a_log10))
        if m < 1e-3:                      # leading term is exact to O(m): P ~ m^k / k!
            log10p = k * np.log10(m) - lgamma(k + 1) / np.log(10.0)
            return float(10.0 ** log10p), float(log10p)
        tot = sum(np.exp(-m + j * np.log(m) - lgamma(j + 1)) for j in range(k))
        p = float(1.0 - tot)
        return p, float(np.log10(max(p, 1e-323)))

    # ---- marginalising the NG15 posterior width over A ----
    def p_brightest_above_marg(self, log10_h, sigma_A=SIGMA_NG15_LIT,
                               a0=A_TARGET_LOG10, n_q=201, n_sig=6.0):
        a = np.linspace(a0 - n_sig * sigma_A, a0 + n_sig * sigma_A, n_q)
        wq = np.exp(-0.5 * ((a - a0) / sigma_A) ** 2)
        wq /= wq.sum()
        mu0 = self.mu(log10_h, a0)
        # mu scales exactly as 10^(2(A-A0)) (lambda ~ P_target ~ A^2)
        vals = np.array([np.sum(wq * (1.0 - np.exp(-mu0 * 10.0 ** (2 * (a - a0)))))
                         for mu0 in np.atleast_1d(mu0)])
        return vals if np.ndim(log10_h) else float(vals[0])

    # ---- the >=1e5-draw MC (cross-check + joints + power conditioning) ----
    def draw_mc(self, n_draw, sigma_A=SIGMA_NG15_LIT, seed=SEED_BASE + 1,
                h_split=None, chunk=20000):
        """Draw-level MC. Per draw: A ~ N(A0, sigma_A); sources above h_split drawn
        explicitly from the Poisson process (the faint remainder can never be brightest and
        its power is added as its conditional mean + Poisson fluctuation of the split part).
        Returns dict with hmax, n_above(rungs), realised total/background power ratios."""
        rng = np.random.default_rng(seed)
        if h_split is None:
            h_split = self._split_for(target_count=200.0)
        h_split = float(h_split)
        # tail sampling grid for h | h >= h_split
        hg = np.linspace(h_split, h_split + 4.0, 4001)
        mu_g = self.mu(hg)                       # decreasing
        mu0 = float(mu_g[0])
        cdf = 1.0 - mu_g / mu0                   # P(h <= x | h >= split)
        cdf = np.maximum.accumulate(cdf)
        mean_h2_below = self._mean_h2_below(h_split)
        lam_below_1 = None

        hmax = np.empty(n_draw)
        nab = {r: np.zeros(n_draw, dtype=np.int64) for r in RUNGS}
        nab[ROBUST_H] = np.zeros(n_draw, dtype=np.int64)
        p_tot = np.empty(n_draw)                 # realised total / P_target(A)
        p_bg = np.empty(n_draw)                  # realised (total - brightest) / P_target(A)
        done = 0
        while done < n_draw:
            m = min(chunk, n_draw - done)
            a = A_TARGET_LOG10 + rng.normal(0.0, sigma_A, size=m)
            scale = 10.0 ** (2 * (a - A_TARGET_LOG10))
            lam_tot = self.lam(A_TARGET_LOG10) * scale
            mu_split = mu0 * scale
            nb = rng.poisson(mu_split)
            nmax = int(nb.max()) if nb.max() > 0 else 1
            u = rng.random((m, nmax))
            hs = np.interp(u, cdf, hg)
            msk = np.arange(nmax)[None, :] < nb[:, None]
            h2 = np.where(msk, 10.0 ** (2 * hs), 0.0)
            hmx = np.where(nb > 0, np.log10(np.sqrt(np.max(h2, axis=1) + 1e-300)), h_split - 9)
            hmax[done:done + m] = hmx
            for r in list(RUNGS) + [ROBUST_H]:
                nab[r][done:done + m] = np.sum(np.where(msk, hs >= r, False), axis=1)
            # faint remainder: lambda_below * <h^2 | below>, Poisson-fluctuated
            lam_below = (lam_tot - mu_split)
            if lam_below_1 is None:
                lam_below_1 = float(self.lam(A_TARGET_LOG10) - mu0)
            n_below = rng.poisson(np.maximum(lam_below, 0.0))
            pow_below = n_below * mean_h2_below
            tot = np.sum(h2, axis=1) + pow_below
            # REFERENCE IS THE MEASUREMENT, not the draw's own ensemble target: the
            # conditioning question is "would NG15 have measured -14.6 in THIS sky", so the
            # realised power is compared with P_target(A0), not P_target(A_drawn).
            tgt0 = band_power_target(A_TARGET_LOG10, self.band)
            p_tot[done:done + m] = tot / tgt0
            p_bg[done:done + m] = (tot - np.max(h2, axis=1)) / tgt0
            done += m
        return dict(hmax=hmax, n_above=nab, p_tot=p_tot, p_bg=p_bg, h_split=h_split,
                    lam=self.lam(A_TARGET_LOG10), mu_split=mu0)

    def _split_for(self, target_count=200.0):
        hg = np.linspace(-19.0, -13.0, 1201)
        m = self.mu(hg)
        i = int(np.argmin(np.abs(m - target_count)))
        return float(min(hg[i], -14.6))

    def _mean_h2_below(self, h_split):
        """<h^2 | h < h_split> over the per-source distribution."""
        lam_all = self.lam(A_TARGET_LOG10)
        tot_pow = band_power_target(A_TARGET_LOG10, self.band)
        hg = np.linspace(h_split, h_split + 5.0, 3001)
        mu_g = self.mu(hg)
        # power carried above the split: integral of h^2 dmu
        h2 = 10.0 ** (2 * hg)
        pow_above = -np.trapezoid(h2, mu_g) if hasattr(np, "trapezoid") \
            else -np.trapz(h2, mu_g)
        n_below = lam_all - mu_g[0]
        return float(max(tot_pow - pow_above, 0.0) / max(n_below, 1e-30))


# ============================================================
# ECCENTRICITY PRIOR (LOTTERY tier-1 parameterisation, verbatim)
# ============================================================
def p_e_above(e_thr, e_peak=0.3, w=0.2, clip=(0.0, 0.9)):
    """P(e >= e_thr) for e ~ LogNormal(mu = ln(e_peak) + w^2, sigma = w), clipped to `clip` --
    the LOTTERY_tier1 SS2 'Taylor/Farr plug-in' family (mode = e_peak).

    EXACT (erfc), not sampled: at (mode 0.3, w 0.1) the e >= 0.5 tail is 5.1 sigma and a
    2e6-draw MC returns a bare 0.0000, which would silently zero a joint."""
    mu = np.log(e_peak) + w ** 2
    if e_thr >= clip[1]:            # the clip piles all mass >= clip[1] AT clip[1]
        return float(_gauss_upper((np.log(clip[1]) - mu) / w))
    return float(_gauss_upper((np.log(max(e_thr, 1e-12)) - mu) / w))


# ============================================================
# ARM-A h-DEPENDENT ECCENTRICITY BOUNDARY (read from the bank, not retyped)
# ============================================================
def armA_threshold_table(path=f"{HSYMT}/reports/gen_armA_surface.npz"):
    d = np.load(path, allow_pickle=True)
    tc, th = d["thr_cols"], d["thr"]              # cols (h, T, k); thr (e*_conf, e*_marg)
    tab = {}
    for row, t in zip(tc, th):
        tab[(float(row[0]), int(row[1]), int(row[2]))] = (float(t[0]), float(t[1]))
    return tab


def e_star_at(h, k=3, T=30, tab=None):
    """The Arm-A boundary threshold e* AT loudness h, for structure k (3+13 / 5+11).

    Measured box: h in {-13.0, -12.75, -12.5} (Arm A) + h = -13.25 census (banked CHORUS
    re-cut: e* = 0.5 at 3+13/T30). BELOW -13.25 the box has no measurement -- e* is returned
    as np.nan with an EXTRAPOLATION flag, never silently filled."""
    tab = tab or armA_threshold_table()
    if abs(h - (-13.25)) < 1e-9:
        return (0.5, 0.5, "banked-census")        # GENERALISE SS2 table, census column
    key = (round(h, 2), T, k)
    if key in tab:
        return (tab[key][0], tab[key][1], "measured")
    if h < -13.25:
        return (np.nan, np.nan, "below-measured-box")
    return (np.nan, np.nan, "off-grid")


# ============================================================
# GATES
# ============================================================
def gate_g0(verbose=True):
    """g0: reproduce the banked P_cond scan (ceiling, quantiles, P(>=-13.9)) with M0."""
    d = np.load(f"{HSYMT}/reports/gen_stage2_pcond_scan.npz", allow_pickle=True)
    band = tuple(d["band_log10f"])
    ceil_bank = float(d["ceiling_log10h"])
    ceil_here = ceiling_log10h(A_TARGET_LOG10, band)
    hb = d["hmax_log10"]
    hm = fixed_class_hmax(100_000, "M0", band=band, seed=SEED_BASE + 11)
    q_bank = np.percentile(hb, [50, 90, 99.9])
    q_here = np.percentile(hm, [50, 90, 99.9])
    p_bank = float((hb >= -13.9).mean())
    p_here = float((hm >= -13.9).mean())
    ok = (abs(ceil_here - ceil_bank) < 1e-9 and np.all(np.abs(q_here - q_bank) < 0.01)
          and abs(p_here - p_bank) < 0.002)
    if verbose:
        print(f"  g0 ceiling  banked {ceil_bank:.6f}  here {ceil_here:.6f}")
        print(f"  g0 p50/p90/p99.9  banked {q_bank}  here {q_here}")
        print(f"  g0 P(>=-13.9)  banked {p_bank:.5f}  here {p_here:.5f}")
        print(f"  g0 -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate_g1(rp, verbose=True):
    """g1: power identity -- lambda * <h^2> == P_target to machine precision, and the MC's
    realised mean total power reproduces it."""
    lhs = rp.lam() * rp.mean_h2
    rhs = band_power_target(A_TARGET_LOG10, rp.band)
    rel = abs(lhs / rhs - 1.0)
    ok = rel < 1e-12
    if verbose:
        print(f"  g1 power identity |lam<h2>/P_target - 1| = {rel:.2e} -> "
              f"{'PASS' if ok else 'FAIL'}")
    return ok


def gate_g2(rp, mc, verbose=True, min_expected=20.0):
    """g2: the >=1e5-draw MC reproduces the analytic Poisson-process tail wherever the MC
    can RESOLVE it (>= min_expected expected hits). Deeper rungs are UNRESOLVABLE by any
    affordable MC -- that is precisely why the analytic Poisson-process tail is the primary
    estimator here; they are reported, never graded."""
    n = len(mc["hmax"])
    ok, rows = True, []
    # auto-chosen RESOLVABLE calibration rungs: where the MC has >= min_expected hits, so the
    # analytic tail is actually graded somewhere and not merely declared unresolvable.
    hg = np.linspace(-16.0, -13.0, 601)
    mug = rp.mu(hg)
    cal = [float(hg[int(np.argmin(np.abs(mug - t)))]) for t in (5e-2, 5e-3, 1e-3)]
    for r in cal + list(RUNGS) + [RUNG_TAIL_REALISED]:
        an = rp.p_brightest_above_marg(r)
        exp_hits = an * n
        mcp = float((mc["hmax"] >= r).mean())
        se = np.sqrt(max(mcp, 1.0 / n) * (1 - mcp) / n)
        if exp_hits < min_expected:
            tag = f"UNRESOLVABLE (expect {exp_hits:.2f} hits in {n})"
        else:
            good = abs(an - mcp) <= max(4 * se, 0.05 * an)
            ok &= good
            tag = "PASS" if good else "FAIL"
        rows.append((r, an, mcp, se))
        if verbose:
            print(f"  g2 rung {r:+.2f}: analytic {an:.3e}  MC {mcp:.3e} "
                  f"(+-{se:.1e}) -> {tag}")
    return ok


def gate_g3(band=BAND_SCAN, verbose=True):
    """g3: mu(rung) is INVARIANT to the integration box's arbitrary edges and to grid
    resolution. The low-mass end sets lambda (how many faint binaries) but carries negligible
    band power, so it must cancel out of the tail; the z edges and the node count must not
    move it either. If they do, the tail is an artefact of a cutoff nobody measured.

    HISTORY (kept, because it is the reason the default moved): the first build used
    z_lo = 1e-4 and this gate FAILED at 1.08 -- the brightest-member tail is NEAR-FIELD
    dominated (g4), so a z floor at 1e-4 truncates 55 % of mu(-13.25). The truncated box is
    retained below as a REGRESSION row, expected to fail, and is not counted in the verdict."""
    base = RealisticPrior(band=band)
    var = [("mc_lo 5.0", RealisticPrior(band=band, mc_range=(5.0, 10.5)), True),
           ("mc_lo 7.0", RealisticPrior(band=band, mc_range=(7.0, 10.5)), True),
           ("grid 2x", RealisticPrior(band=band, n_m=520, n_z=680), True),
           ("grid 4x", RealisticPrior(band=band, n_m=1040, n_z=1360), True),
           ("z_lo 1e-8", RealisticPrior(band=band, z_range=(1e-8, 5.0)), True),
           ("z_max 8", RealisticPrior(band=band, z_range=(1e-7, 8.0)), True),
           ("z_lo 1e-4 [regression]", RealisticPrior(band=band, z_range=(1e-4, 5.0)), False)]
    ok = True
    for name, v, graded in var:
        rel = [abs(v.mu(r) / base.mu(r) - 1.0) for r in RUNGS]
        good = max(rel) < 0.05
        if graded:
            ok &= good
            tag = "PASS" if good else "FAIL"
        else:
            tag = "(expected FAIL, ungraded)" if not good else "(regression no longer bites?)"
        if verbose:
            print(f"  g3 {name:22s}: lam ratio {v.lam()/base.lam():8.3f} | "
                  f"max |dmu/mu| = {max(rel):.2e} -> {tag}")
    return ok


def gate_g4(rp, verbose=True):
    """g4: NOT pass/fail -- a scope disclosure that is banked so the headline cannot be quoted
    without it. Reports where in (d_L, Mc) the sources above each rung actually live. The
    brightest-member tail of ANY homogeneous band population is near-field dominated; if the
    dominant d_L is a few Mpc, the odds are an extrapolation of a cosmological number density
    into the Local Volume, where the real host count is not a continuum."""
    rows = []
    for r in list(RUNGS) + [RUNG_TAIL_REALISED, ROBUST_H]:
        pv = rp.tail_provenance(r)
        rows.append([r, *pv["d_mpc"], *pv["log10_mc"]])
        if verbose:
            print(f"  g4 rung {r:+.3f}: d_L 5/50/95 % = "
                  f"{pv['d_mpc'][0]:.3g}/{pv['d_mpc'][1]:.3g}/{pv['d_mpc'][2]:.3g} Mpc | "
                  f"log10 Mc 5/50/95 % = {pv['log10_mc'][0]:.2f}/{pv['log10_mc'][1]:.2f}/"
                  f"{pv['log10_mc'][2]:.2f}")
    return np.array(rows)


# ============================================================
# THE DELIVERABLES
# ============================================================
def run(n_draw=200_000, band=BAND_SCAN, sigma_A=SIGMA_NG15_LIT, sweep=True, verbose=True):
    p = print if verbose else (lambda *a, **k: None)
    res = {}

    p("\n== GATES ==")
    g0 = gate_g0(verbose)
    rp = RealisticPrior(band=band)
    g1 = gate_g1(rp, verbose)
    g3 = gate_g3(band, verbose)
    g4 = gate_g4(rp, verbose)
    res["g4_provenance"] = g4
    res["g4_cols"] = np.array(["rung", "d_Mpc_p05", "d_Mpc_p50", "d_Mpc_p95",
                               "log10Mc_p05", "log10Mc_p50", "log10Mc_p95"])

    p(f"\n== REALISTIC PRIOR (band {band}, fiducial {FID}) ==")
    p(f"  lambda (expected band sources at A=-14.6) = {rp.lam():.4g}")
    p(f"  <h^2> per source = {rp.mean_h2:.4g}   P_target = "
      f"{band_power_target(A_TARGET_LOG10, band):.4g}")

    mc = rp.draw_mc(n_draw, sigma_A=sigma_A)
    p(f"  MC: {n_draw} draws, h_split = {mc['h_split']:.2f}, "
      f"mu(split) = {mc['mu_split']:.1f}")
    g2 = gate_g2(rp, mc, verbose)

    # ---------------- (a)+(b) P(brightest >= rung): fixed class vs realistic ----------------
    p("\n== (a)+(b) P(brightest member >= rung) -- fixed class vs realistic prior ==")
    ceil = ceiling_log10h(A_TARGET_LOG10, band)
    p(f"  fixed-class HARD CEILING = {ceil:.4f}  (one source carrying 100 % of the pinned "
      f"band power)")
    share = {m: power_share_sample(m, seed=SEED_BASE + 3 + i)
             for i, m in enumerate(("M0", "M1"))}          # M2 shares M0's N, M3 shares M1's
    share["M2"], share["M3"] = share["M0"], share["M1"]
    sig = {"M0": 0.0, "M1": 0.0, "M2": sigma_A, "M3": sigma_A}
    rows = []
    for r in list(RUNGS) + [RUNG_TAIL_REALISED, ROBUST_H]:
        fx = {m: fixed_class_p_above(r, m, share[m], sig[m], band=band) for m in sig}
        an = rp.p_brightest_above_marg(r, sigma_A=sigma_A)
        an_pin = rp.p_brightest_above(r)          # A pinned, realistic loudness
        mcp = float((mc["hmax"] >= r).mean())
        rows.append((r, fx["M0"], fx["M1"], fx["M2"], fx["M3"], an_pin, an, mcp,
                     float(rp.mu(r)), sigma_above_ceiling(r, sigma_A, band)))
        p(f"  rung {r:+.3f} ({sigma_above_ceiling(r, sigma_A, band):+5.2f} sigma_A above the "
          f"ceiling) | M0 {fx['M0']:.2e} M1 {fx['M1']:.2e} M2 {fx['M2']:.2e} "
          f"M3 {fx['M3']:.2e} | M4(A pinned) {an_pin:.3e} | M4(A~NG15) {an:.3e} "
          f"| MC {mcp:.3e} | mu {rp.mu(r):.3e}")
    res["rung_table"] = np.array(rows)
    res["rung_cols"] = np.array(["rung", "M0", "M1", "M2", "M3", "M4_Apinned",
                                 "M4_Amarg", "M4_MC", "mu", "sigma_above_ceiling"])
    res["ceiling_log10h"] = ceil
    res["share_max"] = np.array([float(share[m].max()) for m in ("M0", "M1")])
    # the ablation, stated as a ratio: which suppressed tail was load-bearing?
    p("  ABLATION at the ladder rungs (P / P_M0, with 0/0 shown as the raw P):")
    for r in RUNGS:
        vals = [fixed_class_p_above(r, m, share[m], sig[m], band=band) for m in
                ("M0", "M1", "M2", "M3")]
        p(f"    rung {r:+.3f}: M0 {vals[0]:.2e} -> +Poisson N {vals[1]:.2e} "
          f"-> +NG15 width {vals[2]:.2e} -> both {vals[3]:.2e} "
          f"-> realistic loudness {rp.p_brightest_above_marg(r, sigma_A=sigma_A):.2e}")

    # ---------------- the near-field floor (g4's number, turned into a curve) ----------
    p("\n== (a2) P(brightest >= rung | NO source closer than d_min) ==")
    p("   g4 says the tail is near-field dominated. This is the honest re-quote: how much")
    p("   of each rung's probability survives if the Local Volume is excluded outright.")
    dmin_rows = []
    for r in RUNGS:
        vals = [1.0 - np.exp(-rp.mu(r, d_min_mpc=dm)) for dm in D_MIN_MPC]
        dmin_rows.append((r, *vals))
        p(f"  rung {r:+.3f} | " + " ".join(f"d>{dm:g}Mpc: {v:.2e}"
                                           for dm, v in zip(D_MIN_MPC, vals)))
    res["dmin_table"] = np.array(dmin_rows)
    res["dmin_cols"] = np.array(["rung"] + [f"d_min_{dm:g}Mpc" for dm in D_MIN_MPC])

    # ---------------- conditioning variants ----------------
    p("\n== conditioning variants (NG15 consistency, three readings) ==")
    lo, hi = 10.0 ** (2 * (-2 * sigma_A)), 10.0 ** (2 * (2 * sigma_A))   # +-2 sigma in power
    c2 = (mc["p_tot"] >= lo) & (mc["p_tot"] <= hi)
    c3 = (mc["p_bg"] >= lo) & (mc["p_bg"] <= hi)
    p(f"  power band kept: realised/P_target(-14.6) in [{lo:.3f}, {hi:.3f}] "
      f"(+-2 sigma_A = {2*sigma_A:.2f} dex in A)")
    cond_rows = []
    for r in list(RUNGS) + [ROBUST_H]:
        a1 = float((mc["hmax"] >= r).mean())
        a2 = float((mc["hmax"][c2] >= r).mean()) if c2.sum() else np.nan
        a3 = float((mc["hmax"][c3] >= r).mean()) if c3.sum() else np.nan
        cond_rows.append((r, a1, a2, a3))
        p(f"  rung {r:+.3f} | C1 ensemble {a1:.3e} | C2 total-pinned {a2:.3e} "
          f"| C3 background-pinned {a3:.3e}")
    p(f"  (C2 keeps {c2.mean():.3f} of draws; C3 keeps {c3.mean():.3f})")
    res["cond_table"] = np.array(cond_rows)
    res["cond_cols"] = np.array(["rung", "C1_ensemble", "C2_total_pinned",
                                 "C3_background_pinned"])
    res["cond_keep"] = np.array([c2.mean(), c3.mean()])

    # ---------------- (c) the eccentricity joint ----------------
    p("\n== (c) eccentricity joint ==")
    tab = armA_threshold_table()
    pe = {(pk, w): p_e_above(0.5, e_peak=pk, w=w) for pk in (0.3,) for w in (0.1, 0.2, 0.4)}
    pe3 = {(pk, w): p_e_above(0.3, e_peak=pk, w=w) for pk in (0.3,) for w in (0.1, 0.2, 0.4)}
    for k, v in pe.items():
        p(f"  P(e >= 0.5) under lnN(mode {k[0]}, w {k[1]}) = {v:.4f}")
    for k, v in pe3.items():
        p(f"  P(e >= 0.3) under lnN(mode {k[0]}, w {k[1]}) = {v:.4f}")
    FLAG = {"measured": 0, "banked-census": 1, "below-measured-box": 2, "off-grid": 3}
    joint = []
    for r in RUNGS:
        es, em, flag = e_star_at(r, k=3, T=30, tab=tab)
        # e* is MONOTONE-HARDENING as h falls across the whole measured Arm-A box
        # (0.3 at h >= -12.75 -> 0.5 at h = -13.0/-13.25, 3+13/T30). Below -13.25 there is
        # NO measurement, so the LAST MEASURED value 0.5 is used as a LOWER BOUND on the
        # requirement -- which makes the joint an UPPER BOUND on the odds. Flag banked.
        e_use = es if np.isfinite(es) else 0.5
        for (pk, w), pv_ in pe.items():
            pv = pv_ if abs(e_use - 0.5) < 1e-9 else p_e_above(e_use, e_peak=pk, w=w)
            ph = rp.p_brightest_above_marg(r, sigma_A=sigma_A)
            joint.append((r, 3, e_use, w, ph, pv, ph * pv, FLAG[flag]))
            p(f"  rung {r:+.3f} 3+13 e*={e_use} ({flag}"
              f"{' -> joint is an UPPER BOUND' if flag == 'below-measured-box' else ''}) "
              f"w={w}: P(h) {ph:.3e} x P(e) {pv:.3e} = {ph * pv:.3e}")
    # the robust line: 5+11, e >= 0.3, h >= -12.75, structure proxy N(>= h) >= 5
    p_h_rob = rp.p_brightest_above_marg(ROBUST_H, sigma_A=sigma_A)
    p_k5, log10_pk5 = rp.p_n_above(ROBUST_H, 5)
    p(f"  ROBUST necessary condition  P(brightest >= {ROBUST_H}) = {p_h_rob:.3e} "
      f"(mu = {float(rp.mu(ROBUST_H)):.3e})")
    es_rob, _, flag_rob = e_star_at(ROBUST_H, k=5, T=30, tab=tab)
    p(f"  ROBUST e* read from the bank at (h={ROBUST_H}, k=5, T=30): {es_rob} ({flag_rob})")
    for (pk, w), pv in pe3.items():
        joint.append((ROBUST_H, 5, es_rob, w, p_k5, pv, p_k5 * pv, FLAG[flag_rob]))
        p(f"  ROBUST 5+11 h>={ROBUST_H} e*>={es_rob} w={w}: log10 P(N>=5) = {log10_pk5:.1f} "
          f"x P(e) {pv:.4f} -> log10 P_joint = {log10_pk5 + np.log10(pv):.1f}")
    res["robust"] = np.array([p_h_rob, float(rp.mu(ROBUST_H)), log10_pk5, es_rob])
    res["robust_cols"] = np.array(["P_brightest_ge_m12p75", "mu_m12p75", "log10_P_N5",
                                   "e_star_k5"])
    res["joint_table"] = np.array(joint)
    res["joint_cols"] = np.array(["rung", "k_struct", "e_star", "e_width", "P_h", "P_e",
                                  "P_joint", "e_star_flag"])
    res["joint_flag_key"] = np.array(["0=measured", "1=banked-census",
                                      "2=below-measured-box(UPPER BOUND)", "3=off-grid"])

    # ---------------- sensitivity sweep ----------------
    if sweep:
        p("\n== sensitivity: the tail is controlled by the mass-function cutoff ==")
        sw = []
        for a_ in (-0.5, 0.0, 0.5):
            for ms in (8.0, 8.5, 9.0):
                for chirp in (False, True):
                    r2 = RealisticPrior(band=band, alpha=a_, log10_Mstar=ms,
                                        mass_is_chirp=chirp)
                    vals = [r2.p_brightest_above_marg(r, sigma_A=sigma_A) for r in RUNGS]
                    sw.append((a_, ms, float(chirp), r2.lam(), *vals))
                    p(f"  alpha {a_:+.1f} log10M* {ms:.1f} "
                      f"{'CHIRP' if chirp else 'TOTAL':5s}: lam {r2.lam():.3g} | "
                      + " ".join(f"{r:+.2f}: {v:.2e}" for r, v in zip(RUNGS, vals)))
        res["sweep"] = np.array(sw)
        res["sweep_cols"] = np.array(["alpha", "log10_Mstar", "mass_is_chirp", "lambda",
                                      "P_m13p87", "P_m13p5", "P_m13p25"])
        # amplitude-convention sensitivity (8/sqrt10 vs 2)
        p("\n== sensitivity: strain-amplitude convention ==")
        conv = []
        for fac, name in ((1.0, "8/sqrt10 (h_c bookkeeping)"), (2.0 / (8 / np.sqrt(10)),
                                                                "h0 = 2(GMz)^5/3.. ")):
            r2 = RealisticPrior(band=band)
            r2.K = r2.K * fac
            r2.mean_h2 = r2.mean_h2 * fac ** 2
            vals = [r2.p_brightest_above_marg(r, sigma_A=sigma_A) for r in RUNGS]
            conv.append((fac, r2.lam(), *vals))
            p(f"  x{fac:.3f} {name}: lam {r2.lam():.3g} | "
              + " ".join(f"{r:+.2f}: {v:.2e}" for r, v in zip(RUNGS, vals)))
        res["conv"] = np.array(conv)
        # NG15 width sensitivity
        p("\n== sensitivity: NG15 posterior width ==")
        wid = []
        for s in (0.0, SIGMA_NG15_REPO, SIGMA_NG15_LIT, 0.15):
            vals = [rp.p_brightest_above_marg(r, sigma_A=max(s, 1e-6)) for r in RUNGS]
            wid.append((s, *vals))
            p(f"  sigma_A {s:.3f}: " + " ".join(f"{r:+.2f}: {v:.2e}"
                                                for r, v in zip(RUNGS, vals)))
        res["width"] = np.array(wid)

    res["gates"] = np.array([g0, g1, g2, g3])
    res["gate_names"] = np.array(["g0_reproduce_banked_scan", "g1_power_identity",
                                  "g2_MC_vs_analytic", "g3_cutoff_invariance"])
    # MC raw columns: hmax kept (f32 -- 1e-7 relative, far finer than any use of it) so the
    # cross-check is re-derivable; the power-ratio columns are banked as HISTOGRAMS, since
    # their only use is the C2/C3 acceptance masks, which are already banked as cond_table.
    # Full-f64 raw for all three was 4.6 MB of committed bank for no readback anyone needs.
    res["hmax_mc"] = mc["hmax"].astype(np.float32)
    pe_ = np.linspace(0.0, 3.0, 601)
    res["p_tot_hist_edges"] = pe_
    res["p_tot_hist"] = np.histogram(mc["p_tot"], bins=pe_)[0]
    res["p_bg_hist"] = np.histogram(mc["p_bg"], bins=pe_)[0]
    res["lam"] = rp.lam()
    res["mean_h2"] = rp.mean_h2
    res["band_log10f"] = np.array(band)
    res["sigma_A"] = sigma_A
    res["fid"] = np.array([FID["alpha"], FID["log10_Mstar"], FID["beta"], FID["z0"]])
    res["n_draw"] = n_draw
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200_000)
    ap.add_argument("--band", default="scan", choices=["scan", "gen"])
    ap.add_argument("--out", default=OUT_BANK)
    ap.add_argument("--no-sweep", action="store_true")
    a = ap.parse_args()
    band = BAND_SCAN if a.band == "scan" else BAND_GEN
    r = run(n_draw=a.n, band=band, sweep=not a.no_sweep)
    np.savez_compressed(a.out, **{k: v for k, v in r.items()})
    print(f"\nbanked -> {a.out}")
