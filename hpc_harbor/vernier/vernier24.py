#!/usr/bin/env python3
"""
VERNIER-24 -- the clean-band wheel on the AUDITED sky. EXPLORATORY, NOT
pre-registered. Clean room: imports only vernier_toy / vernier6 / vernier7 /
vernier16 (all my own anchored code).

WHAT THE CENSUS AUDIT CHANGED (CENSUS_CLASH, ACCRE, staged)
  "One source per bin above 8 nHz" was an artefact of the census's fixed
  n_src = 256 with no detectability floor. The corrected accounting, same dN/df
  law with its proper normalisation:
      sources/bin at 10 / 15 / 20 nHz = 6.9e4 / 1.6e4 / 5.4e3
      effective in-bin competitors N_eff ~ 8-180
      lambda_fore ~ 0.31, FLAT in N_bin (Pareto gamma_eff = 1.56; structural for
      gamma < 2, so NO one-source bins exist anywhere)
      loudest in-bin strain roughly right (0.25-0.67 dex) at 10-20 nHz ONLY.

  So the one-source-per-bin sky is a KNOWN ARTEFACT and is not a physical
  bracket. It is kept here only as the zero-competitor control.

THE ARMS
  ARM B  -- THE PHYSICAL SKY. N_bin from the audited counts (VERNIER-4's
            N_c ~ f^(-11/3) reproduces them: 6.9e4 / 1.6e4 / 5.4e3 at
            10 / 15 / 20 nHz). Strains drawn from a Pareto of index
            gamma_eff = 1.56 through a proper h(Mc, f, D) link, so the bin
            total goes as h^2/bin ~ f^(-7/3) -- NOT h independent of f, which
            is the census's gamma_eff = 5.51 flaw and is not reproduced here.
            The top 300 are modelled individually; the rest is a Gaussian sea.
  ARM Am -- ZERO-COMPETITOR CONTROL: one source at Arm B's own lambda_fore
            share of the bin, everything else Gaussian. Isolates competition
            from loudness.
  ARM A  -- the superseded one-source-per-bin sky (lambda = 0.90). KEPT AND
            RELABELLED as a known-artefact control, never as a bracket.
  M1     -- one source holding the whole bin, no remainder at all.

THE DECIDING TEST (deliverable 1) is REALIZED separation, not dof counting.
VERNIER-23 Part 1's M_max = 17-20 is a counting bound and its Part 2 "no lever"
was computed on the artefact census. At 10-20 nHz the lever (1.4-7.2 bins of
pulsar-term line spread) and the confusion (N_eff ~ 8-180) now coexist. So: fit
M sources sharing one bin, WITH pulsar terms and with the pulsar distances
MARGINALISED over their 10 pc priors, at M = 2, 5, 10, 17, 25, 50 against the
Arm B sky, and measure where recovery actually degrades -- by Fisher
conditioning AND by parameter error -- at 10 / 15 / 20 nHz. The answer is a
SEPARABILITY FRONTIER in-band: which bins are winnable and which are not.

THE WHEEL RUNS ACROSS THE BAND, NOT INSIDE ONE BIN. Seven 1/T bins span
10.0-19.5 nHz; each turn searches every bin, takes the loudest detection
anywhere in the band, and certifies against the ACCUMULATED comb of every
source found so far. That accumulation is the vernier -- sources at different
frequencies and sky positions give different fringe spacings
dL = c/(f (1-cos mu)) in the same pulsar, and it is their commensurability, not
any single comb, that can pin a distance. A single source cannot certify at any
SNR, because every fringe centre of one comb gives an identical template.

SAME ARRAY IN EVERY ARM: 116 pulsars, heterogeneous timing (VERNIER-13's
sigma_log = 0.40 dex spread in 1/sigma_rms), log-uniform distances 0.3-5 kpc,
10 pc distance priors. Sources are NOT pre-resolved -- they are searched for on
a 6000-point sky grid. Full-array JOINT inference (VERNIER-12): the shared
source-parameter smear is cut by sqrt(m) with m = n_psr, and every turn re-fits
ALL found amplitudes simultaneously before subtracting.
lambda_eff-CONSISTENT confusion (the VERNIER-21/22 debt): each interloper enters
the per-pulsar noise budget weighted by its MEASURED squared template overlap
with the target, not by its raw power.

CERTIFICATION is decided by CRITERION-v2.2 -- a null-calibrated floor, Gumbel
where the null zero-fraction is <= 20 % and the empirical q95 with a bootstrap
error above the gate, the zero-fraction carried as a required column -- rather
than by VERNIER-1's hard rho_P >= 2 bar, which is reported beside it.

CONTROLS: M = 1 control, pure-noise null through the whole wheel, per-turn
SOURCE-LEVEL scrambled null, and VERNIER-18's genuine/split/blend accounting
with its chance-match floor printed.

SUB-FRINGE NOTE. Every distance grid has an EVEN point count and is centred on
the prior mean, so it never lands on a fringe centre by construction; and the
certification floor is refit from nulls evaluated on the SAME grid, which
absorbs the residual sampling bias. (Both halves of the banked lesson.)
"""
import os
import numpy as np

import vernier_toy as V
import vernier6 as X
import vernier7 as S7                  # measured census ladder LAD (rank <-> dex)
import vernier16 as S

YR, C_PC_S, REPORTS = V.YR, V.C_PC_S, V.REPORTS
V.PTERM_STATIC = False

# ------------------------------------------------------------------ constants
R_AMP     = 1.196          # pulsar-term / Earth-term residual ratio (VERNIER-8)
CERT_BAR  = 2.0            # VERNIER-1's hard per-pulsar bar, reported alongside
A_GWB     = 2.4e-15        # NANOGrav-15yr-like h_c at f_yr           (VERNIER-4)
N_YR      = 1000.0         # analytic sources per bin at f = 1/yr     (VERNIER-4)
F_YR      = 1.0 / YR
KAPPA     = 1.0 / 3.0      # sky/polarisation-averaged antenna power  (VERNIER-4)
T_OBS     = 20.0 * YR
DF        = 1.0 / T_OBS
BAND      = 10e-9 + np.arange(7) * DF                # seven bins, 10.0-19.5 nHz
F_REPORT  = np.array([10.0, 15.0, 20.0]) * 1e-9      # the audit's report points
GAMMA_EFF = 1.56           # AUDITED Pareto index of the in-bin strain function
LAM_AUDIT = 0.31           # AUDITED lambda_fore, flat in N_bin
N_EFF_AUD = (8.0, 180.0)   # AUDITED effective in-bin competitor range
MC_FID    = 9.0            # log10 Mc used for the h(Mc, f, D) link check
N_PSR     = 116
SIG_LOG_W = 0.40           # VERNIER-13 heterogeneous weight spread, dex
SIG_TOA   = 100e-9         # MEDIAN residual rms
CADENCE   = 14.0           # days
PRIOR_PC  = 10.0           # SKA-era timing-parallax prior (Xiao et al.)
LOG10_MC  = 9.0
N_BRIGHT  = 300            # sources modelled individually in Arm B (VERNIER-4)
LAM_A     = 0.90           # Arm A: fraction of bin power in the modelled source
DET_ALPHA = 0.05
N_NULL_F  = 150            # nulls for the detection floor, per bin
N_NULL_C  = 100            # nulls for the per-pulsar certification floor
NPTS      = 1600           # distance-grid points (EVEN -- see sub-fringe note)
NGRID     = 6000           # sky grid for the F-statistic
TOL_RAD   = 0.15           # VERNIER-16's tightened genuine-match tolerance
TOP_TRUE  = 10             # GENUINE is scored against the top-10 of the ladder
RHO_REF   = 20.0           # rank-1 array SNR, as everywhere in the series
DET_SNR   = 8.0
EULER_G   = 0.5772156649015329


# ============================================================ census / loudness
def S_gwb(f):
    """GWB residual PSD, s^3."""
    return (A_GWB ** 2 / (12 * np.pi ** 2 * F_YR ** 3)) \
        * (np.asarray(f) / F_YR) ** (-13.0 / 3.0)


def S_white(sigma_toa, cadence_days=CADENCE):
    return 2.0 * sigma_toa ** 2 * (cadence_days * 86400.0)


def N_c(f):
    """VERNIER-4's analytic occupancy: sources per 1/T bin."""
    return N_YR * (np.asarray(f) / F_YR) ** (-11.0 / 3.0)


def bin_power(f):
    """Total residual amplitude^2 held by ONE 1/T bin, s^2."""
    return 2.0 * S_gwb(f) * DF


def h2_per_bin(f):
    """
    Strain^2 held by ONE 1/T bin, from the h(Mc, f, D) link.

    h_c ~ f^(-2/3) for a GW-driven SMBHB population, so h_c^2 ~ f^(-4/3) and the
    strain power per LINEAR bin is dh^2/df * (1/T) = h_c^2(f)/f * (1/T), i.e.
        h^2 / bin  ~  f^(-7/3)
    which is checked numerically in deliverable 0. Residual power per bin is the
    same divided by (2 pi f)^2, hence the f^(-13/3) of `bin_power`. The census's
    gamma_eff = 5.51 came from drawing h independently of f; that is exactly what
    this link avoids.
    """
    hc2 = A_GWB ** 2 * (np.asarray(f) / F_YR) ** (-4.0 / 3.0)
    return hc2 / np.asarray(f) * DF


def pareto_ladder(rng, f, n_bright=N_BRIGHT, gamma=GAMMA_EFF):
    """
    In-bin residual amplitudes^2 drawn from the AUDITED strain function.

    N_bin sources per bin from the audited counts (VERNIER-4's N_c reproduces
    them), strains from a Pareto of index gamma_eff = 1.56 --  h = u^(-1/gamma).
    Because gamma < 2 the sum is dominated by the loudest few, so lambda_fore is
    O(0.3) and INDEPENDENT of N_bin: that is the structural reason no one-source
    bin exists, and it is what the audit measured.

    The whole draw is normalised so the bin holds `bin_power(f)`; the top
    n_bright are returned individually and the rest become the Gaussian sea.
    """
    n = int(max(round(float(N_c(f))), 1))
    h = rng.uniform(size=n) ** (-1.0 / gamma)
    p = np.sort(h ** 2)[::-1]
    p = p / p.sum()
    tot = float(bin_power(f))
    lam = float(p[0])
    n_eff = float(1.0 / np.sum(p ** 2))
    n_half = int(np.sum(h >= h.max() / 10 ** 0.5))
    nb = min(n_bright, n)
    return (p[:nb] * tot, float(p[nb:].sum() * tot), n,
            dict(lam_fore=lam, n_eff=n_eff, n_halfdex=n_half, n_bin=n))


def bin_sources(rng, f, arm, n_bright=N_BRIGHT, lam_a=LAM_A, gain=1.0):
    """
    Residual amplitude^2 of each individually-modelled source in the bin, s^2,
    plus the PSD of the Gaussian remainder, plus census diagnostics.

    ARM B  : THE PHYSICAL SKY -- the audited Pareto draw above.
    ARM Am : ZERO-COMPETITOR CONTROL -- one source at Arm B's own lambda_fore
             share, everything else Gaussian.
    ARM A  : the superseded one-source-per-bin artefact sky (lambda = 0.90).
    M1     : one source holding the WHOLE bin, no remainder.

    `gain` multiplies the bright end and is taken back out of the remainder, so
    the bin total is conserved: the loudness sweep REDISTRIBUTES the same GWB
    power between resolved and unresolved rather than adding signal.
    """
    tot = float(bin_power(f))
    if arm == "B":
        R2, _sea, _nb, diag = pareto_ladder(rng, f, n_bright)
        R2 = R2 * gain
        return R2, max(tot - float(R2.sum()), 1e-6 * tot) / (2.0 * DF), diag
    if arm == "Am":
        R2 = np.array([LAM_AUDIT * tot])
        diag = dict(lam_fore=LAM_AUDIT, n_eff=1.0, n_halfdex=1, n_bin=1)
    elif arm == "A":
        R2 = np.array([lam_a * tot])
        diag = dict(lam_fore=lam_a, n_eff=1.0, n_halfdex=1, n_bin=1)
    elif arm == "iso":                        # same loudness as Am, NO sea at all
        return (np.array([LAM_AUDIT * tot * gain]), 1e-6 * tot / (2.0 * DF),
                dict(lam_fore=LAM_AUDIT, n_eff=1.0, n_halfdex=1, n_bin=1))
    elif arm == "M1":
        return (np.array([tot * gain]), 1e-6 * tot / (2.0 * DF),
                dict(lam_fore=1.0, n_eff=1.0, n_halfdex=1, n_bin=1))
    else:
        raise ValueError(arm)
    R2 = R2 * gain
    return R2, max(tot - float(R2.sum()), 1e-6 * tot) / (2.0 * DF), diag


# ============================================================ criterion-v2.2
def _gumbel_q(x, alpha, axis=-1):
    """
    Gumbel (1-alpha) quantile by MOMENTS: scale = sd sqrt(6)/pi, loc = mean - g scale.

    Moments rather than MLE because this floor is refit, with a bootstrap, once
    per pulsar per turn; the fit ERROR is the bootstrap sd either way, which is
    what the convention actually requires.
    """
    m = np.mean(x, axis=axis)
    s = np.std(x, axis=axis)
    sc = np.maximum(s * np.sqrt(6.0) / np.pi, 1e-30)
    return (m - EULER_G * sc) - sc * np.log(-np.log(1.0 - alpha))


def recut_floor(x, rng, alpha=DET_ALPHA, gate=0.20, nboot=80):
    """
    CRITERION-v2.2, BINDING (project_progress.md 10.16a).

    Gumbel floor where the null zero-fraction is <= 20 %; above the gate the
    empirical (1-alpha) quantile with a BOOTSTRAP sd. The zero-fraction is a
    REQUIRED column, not a caveat. The test itself is unchanged: pass iff the
    statistic exceeds floor + its own error.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 10:
        return dict(floor=np.inf, err=np.inf, zero_frac=1.0, kind="degenerate")
    zf = float(np.mean(x <= 0.0))
    xb = x[rng.integers(0, x.size, size=(nboot, x.size))]
    if zf <= gate:
        fl = float(_gumbel_q(x, alpha))
        bs = _gumbel_q(xb, alpha, axis=1)
        kind = "gumbel"
    else:
        fl = float(np.quantile(x, 1 - alpha))
        bs = np.quantile(xb, 1 - alpha, axis=1)
        kind = "emp_q95"
    return dict(floor=fl, err=float(np.std(bs)), zero_frac=zf, kind=kind)


# ============================================================ array and bins
OMC_MIN = 5e-2             # series-wide clip on 1 - cos(mu) (VERNIER-toy)


def antenna_reg(P, G, dmin=OMC_MIN):
    """
    F+, Fx with the 1/(1 - cos mu) divergence REGULARISED at the series-wide
    clip, exactly as V.tau_of / V.dL_of are clipped everywhere else.

    VERNIER-16's antenna_grid clips the denominator at 1e-6, which is fine for an
    F-statistic (a ratio) but not for an SNR or a Fisher: one nearly antipodal
    pulsar then carries a response of ~1e5 and dominates every budget in the run.
    With that clip the per-pulsar rho_P `max` column is a single pulsar's
    coordinate artefact, and the marginalised Fisher reports sky errors of ~1e-3
    rad at rho = 3.8. Both are removed here; the effect is reported.
    """
    z = np.array([0.0, 0.0, 1.0])
    u = np.cross(z[None, :], G)
    nn = np.linalg.norm(u, axis=1, keepdims=True)
    u = np.where(nn > 1e-12, u / np.maximum(nn, 1e-30), np.array([1.0, 0.0, 0.0]))
    v = np.cross(G, u)
    pu, pv = P @ u.T, P @ v.T
    d = np.maximum(1.0 + P @ G.T, dmin)
    return (pu ** 2 - pv ** 2) / (2 * d), (pu * pv) / d



def make_array(rng, n_psr=N_PSR, sig_log=SIG_LOG_W, sigma_toa=SIG_TOA):
    """116 pulsars, heterogeneous timing, log-uniform distances 0.3-5 kpc."""
    P = V.iso_dirs(n_psr, rng)
    w = 10.0 ** rng.normal(0.0, sig_log, n_psr)
    w = w / np.sqrt(np.mean(w ** 2))          # unit mean POWER; w = sigma_med/sigma_a
    L = 10 ** rng.uniform(np.log10(300.0), np.log10(5000.0), n_psr)
    G = S.fib_sphere(NGRID)                   # the search grid, shared by all bins
    fpg, fcg = antenna_reg(P, G)
    return dict(P=P, w=w, L=L, sigma_toa=sigma_toa, n_psr=n_psr,
                G=G, fpg=fpg, fcg=fcg)


def make_bin(rng, arr, f, arm, gain=1.0):
    """One frequency bin: sources, per-(pulsar,source) SNR, whitening, overlaps."""
    P, w, n_psr = arr["P"], arr["w"], arr["n_psr"]
    R2, S_sea, diag = bin_sources(rng, f, arm, gain=gain)
    M = R2.size
    G_all = arr["G"]
    Gt = G_all[rng.choice(NGRID, M, replace=False)]
    fp, fc = antenna_reg(P, Gt)
    ant = np.sqrt(fp ** 2 + fc ** 2)
    ant = ant / np.sqrt(np.mean(ant ** 2)) * np.sqrt(KAPPA)     # mean power = kappa

    Swm = S_white(arr["sigma_toa"])
    s_k = np.sqrt(R2 * T_OBS / (2.0 * Swm))       # per-pulsar SNR scale vs WHITE
    rho0 = ant * w[:, None] * s_k[None, :]        # (n_psr, M) against white alone
    C_sea = (S_sea / Swm) * w ** 2                # sea/white, per pulsar
    omega = w / np.sqrt(1.0 + C_sea)              # whitening incl. the Gaussian sea
    rho0 = rho0 / np.sqrt(1.0 + C_sea)[:, None]   # ... and folded into the SNRs

    psi = rng.uniform(0, np.pi, M)
    ci = rng.uniform(-1, 1, M)
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    Fp = fp * c2[None, :] + fc * s2[None, :]
    Fc = -fp * s2[None, :] + fc * c2[None, :]
    Tm = (Fp * ((1 + ci ** 2) / 2)[None, :] + 1j * Fc * ci[None, :]) * omega[:, None]
    Tm = Tm * np.exp(1j * rng.uniform(0, 2 * np.pi, M))[None, :]
    Tn = Tm / np.maximum(np.sqrt(np.sum(np.abs(Tm) ** 2, axis=0)), 1e-30)
    ov2 = np.abs(Tn.conj().T @ Tn) ** 2           # MEASURED squared overlaps

    omc = np.clip(1.0 - P @ Gt.T, 5e-2, None)
    # search templates over the WHOLE sky grid, whitened
    Fpw = arr["fpg"] * omega[:, None]
    Fcw = arr["fcg"] * omega[:, None]
    return dict(f=f, arm=arm, R2=R2, S_sea=S_sea, diag=diag, Gt=Gt, rho0=rho0, C_sea=C_sea,
                omega=omega, Tn=Tn, ov2=ov2, omc=omc, M=M, G=G_all,
                Fpw=Fpw, Fcw=Fcw)


def rho_budget(B, k, sub):
    """
    lambda_eff-CONSISTENT per-pulsar SNR of source k:

        noise_a = 1 + sum_{i != k, unsubtracted} rho0_ai^2 * ov2_ki

    -- each interloper weighted by its measured squared template overlap with the
    target. That suppression is the VERNIER-21 measurement the VERNIER-22 wheel
    left out; raw |R|^2 is the pessimistic limit ov2 -> 1.
    """
    others = [i for i in range(B["M"]) if i != k and i not in sub]
    if others:
        conf = (B["rho0"][:, others] ** 2 * B["ov2"][k, others][None, :]).sum(axis=1)
    else:
        conf = np.zeros(B["rho0"].shape[0])
    return B["rho0"][:, k] / np.sqrt(1.0 + conf)


# ============================================================ detection
def fstat_floor(rng, B, n_null=N_NULL_F):
    """Recut floor on max-F over the sky grid, from PURE-NOISE realisations."""
    Fpw, Fcw = B["Fpw"], B["Fcw"]
    n = Fpw.shape[0]
    Muu = np.sum(Fpw * Fpw, axis=0)
    Muv = np.sum(Fpw * Fcw, axis=0)
    Mvv = np.sum(Fcw * Fcw, axis=0)
    det = np.where(np.abs(Muu * Mvv - Muv ** 2) < 1e-12, 1e-12,
                   Muu * Mvv - Muv ** 2)
    Dm = (rng.normal(size=(n, n_null)) + 1j * rng.normal(size=(n, n_null))) / np.sqrt(2)
    bu, bv = Fpw.T @ Dm, Fcw.T @ Dm                 # (n_grid, n_null)
    Fv = (Mvv[:, None] * np.abs(bu) ** 2
          - 2 * Muv[:, None] * np.real(np.conj(bu) * bv)
          + Muu[:, None] * np.abs(bv) ** 2) / det[:, None]
    return recut_floor(Fv.max(axis=0), rng)


def inject(rng, B, null=False):
    """Whitened data: sum of the modelled sources + unit complex Gaussian noise."""
    n = B["Tn"].shape[0]
    d = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
    if null:
        return d
    rho_arr = np.sqrt((B["rho0"] ** 2).sum(axis=0))          # per source
    return d + (B["Tn"] * rho_arr[None, :]).sum(axis=1)


# ============================================================ certification
def certify(rng, arr, det, scramble=False, npts=NPTS, n_null=N_NULL_C,
            prior_pc=PRIOR_PC, soft=False):
    """
    Per-pulsar certification against the ACCUMULATED comb of every source found
    so far. `det` is a list of (bin, k, subtracted_set) in discovery order; the
    newest entry is the one whose pulsar-term power f_cert weights.

    STATISTIC. Certification means the distance is pinned to ONE fringe, so the
    statistic is the fringe CONTRAST, anchored on the data and not on the truth:

        S_a = lnL(top) - max{ lnL(L) : |L - L_top| > dL_min/2 }

    the winning fringe against the best rival. (An earlier draft used
    peak-minus-median, which only says a COMB is present: it fired at f_cert =
    0.84 on signal and 0.76 on the scrambled null, because the trough between
    fringes sits near -2 rho^2 whether or not the right fringe was picked. That
    quantity is kept as S_comb for comparison and is NOT used to certify.)

    FLOOR: the null has to be the INCOHERENT-COMB null, not the no-comb one. A
    single source cannot break the fringe degeneracy at all -- every fringe
    centre gives an identical template, so the contrast is exactly zero up to the
    chirp envelope. What can produce a contrast is COMMENSURABILITY between the
    combs of several sources with different dL = c/(f(1-cos mu)). Combs at random
    relative offsets align by chance too, and at high SNR they align well, so the
    floor is built from n_null realisations of the same accumulated statistic with
    each comb's phase randomised. (A no-comb floor was tried first and is reported
    beside it: it lets chance alignment through, and the scrambled-null column
    then certifies at ~70 % of the signal rate.)

    Bar = max(ln K, floor + err), K = fringes in the prior window at the finest
    comb.

    The shared source-parameter smear is the FULL-ARRAY joint one (VERNIER-12):
    sigma_phi_a = 2 pi sigma_f tau_a / sqrt(n_psr).
    """
    n_psr = arr["n_psr"]
    comps = []
    for (B, k, sub) in det:
        rho = rho_budget(B, k, sub)
        rho_arr = float(np.sqrt(np.sum(rho ** 2)))
        sig_f = 1.0 / (2 * np.pi * T_OBS * max(rho_arr, 1e-6))
        omc = B["omc"][:, k]
        tau = V.tau_of(arr["L"], omc)
        dec = np.exp(-0.5 * (2 * np.pi * sig_f * tau / np.sqrt(n_psr)) ** 2)
        mc_k = float(B["mc"][k]) if ("mc" in B) else LOG10_MC
        comps.append(dict(f=B["f"], tc=float(V.t_merge(B["f"], mc_k)),
                          omc=omc, rho_P=R_AMP * rho, dec=dec,
                          dL=V.dL_of(B["f"], omc), rho_arr=rho_arr))
    Bt, kt, _ = det[-1]
    pw = Bt["rho0"][:, kt] ** 2
    dL_min = np.min(np.stack([c["dL"] for c in comps]), axis=0)
    lnK = np.log(np.maximum(6.0 * prior_pc / dL_min, 1.0))

    cert = np.zeros(n_psr, bool); ctru = np.zeros(n_psr, bool)
    stat = np.zeros(n_psr); comb = np.zeros(n_psr); Rcoh = np.zeros(n_psr)
    zfr = np.zeros(n_psr); fl_i = np.zeros(n_psr); fl_n = np.zeros(n_psr)
    kinds = []

    def _contrast(curve, Lg, hw):
        """peak minus the best rival outside +-hw of the peak; rows if 2-D."""
        if curve.ndim == 1:
            it = int(np.argmax(curve))
            far = np.abs(Lg - Lg[it]) > hw
            return float(curve[it] - (curve[far].max() if far.any() else curve.min()))
        it = np.argmax(curve, axis=1)
        far = np.abs(Lg[None, :] - Lg[it][:, None]) > hw
        return curve[np.arange(curve.shape[0]), it] \
            - np.where(far, curve, -np.inf).max(axis=1)

    for a in range(n_psr):
        Lstar = arr["L"][a]
        Lg = Lstar + np.linspace(-3 * prior_pc, 3 * prior_pc, npts)   # EVEN npts
        lnl = np.zeros(npts)
        n_inc = np.zeros((n_null, npts))          # incoherent-comb null
        n_noc = np.zeros((n_null, npts))          # no-comb null (reported only)
        phi_ref = None
        for c in comps:
            Ls = Lstar + (rng.uniform(-500.0, 500.0) if scramble else 0.0)
            phi, env = X.phase_env(Lg, max(Ls, 50.0), c["f"] / 2.0, c["tc"],
                                   c["omc"][a], T_OBS, 2)
            env = env * c["dec"][a]
            r = c["rho_P"][a]
            cp, sp = np.cos(phi), np.sin(phi)
            nc, ns = rng.normal(), rng.normal()
            lnl += r ** 2 * (env * cp - 1.0) + r * env * (nc * cp + ns * sp)
            nc0 = rng.normal(size=n_null)[:, None]
            ns0 = rng.normal(size=n_null)[:, None]
            n_noc += r * env[None, :] * (nc0 * cp[None, :] + ns0 * sp[None, :])
            # same comb, RANDOM relative offset -- chance commensurability
            dl = rng.uniform(0, 2 * np.pi, n_null)[:, None]
            cq = cp[None, :] * np.cos(dl) - sp[None, :] * np.sin(dl)
            sq = sp[None, :] * np.cos(dl) + cp[None, :] * np.sin(dl)
            n_inc += r ** 2 * (env[None, :] * cq - 1.0) \
                + r * env[None, :] * (nc0 * cq + ns0 * sq)
            if phi_ref is None:
                phi_ref = phi
        hw = 0.5 * dL_min[a]
        comb[a] = float(lnl.max() - np.median(lnl))
        stat[a] = _contrast(lnl, Lg, hw)
        itop = int(np.argmax(lnl))
        fl = recut_floor(_contrast(n_inc, Lg, hw), rng)
        fn = recut_floor(_contrast(n_noc, Lg, hw), rng)
        zfr[a] = fl["zero_frac"]; kinds.append(fl["kind"])
        fl_i[a] = fl["floor"] + fl["err"]; fl_n[a] = fn["floor"] + fn["err"]
        cert[a] = stat[a] > max(lnK[a], fl_i[a])
        if cert[a]:
            ctru[a] = abs(Lg[itop] - Lstar) <= 0.5 * dL_min[a]
        if soft:
            pr = np.exp(-0.5 * ((Lg - Lstar) / prior_pc) ** 2)
            wg = np.exp(lnl - lnl.max()) * pr
            wg = wg / max(np.trapezoid(wg, Lg), 1e-300)
            Rcoh[a] = float(np.abs(np.trapezoid(wg * np.exp(1j * phi_ref), Lg)))
    rp = comps[-1]["rho_P"]
    return dict(cert=cert, cert_true=ctru, stat=stat, comb=comb, Rcoh=Rcoh,
                n_cert=int(cert.sum()), n_true=int((cert & ctru).sum()),
                f_cert=float(pw[cert].sum() / max(pw.sum(), 1e-30)),
                f_cert_true=float(pw[cert & ctru].sum() / max(pw.sum(), 1e-30)),
                f_soft=float(np.sum(Rcoh ** 2 * pw) / max(pw.sum(), 1e-30)),
                rho_P_med=float(np.median(rp)), rho_P_max=float(np.max(rp)),
                rho_P_p90=float(np.percentile(rp, 90)),
                rho_arr=comps[-1]["rho_arr"],
                hard_bar=float(pw[rp >= CERT_BAR].sum() / max(pw.sum(), 1e-30)),
                lnK=float(np.median(lnK)), dL=float(np.median(dL_min)),
                S_med=float(np.median(stat)), S_max=float(np.max(stat)),
                comb_med=float(np.median(comb)),
                floor_inc=float(np.median(fl_i)), floor_noc=float(np.median(fl_n)),
                n_cert_nocomb=int(np.sum(stat > np.maximum(lnK, fl_n))),
                zero_frac=float(np.median(zfr)),
                kind=("gumbel" if all(x == "gumbel" for x in kinds) else "mixed"))


# ============================================================ the wheel
def rank_of(gain_dex):
    """Census rank reachable at a sensitivity gain of `gain_dex` (VERNIER-7 ladder)."""
    return int(np.sum(RHO_REF * 10.0 ** (S7.LAD + gain_dex) >= DET_SNR))


def wheel(seed, arr, arm, band=BAND, gain=1.0, max_turn=8, null=False, soft=False):
    """
    detect -> certify (criterion-v2.2) -> recover pulsar terms -> subtract -> repeat,
    ACROSS THE BAND.

    Detection: F-statistic over a 6000-point sky grid in every bin, against that
    bin's own null-calibrated floor; the loudest excess anywhere in the band is
    taken. Every turn re-fits ALL amplitudes found in a bin simultaneously before
    subtracting (VERNIER-18's joint CLEAN -- the VERNIER-12 requirement).
    Recovery is scored GENUINE / SPLIT / BLEND against the injected truth, with
    GENUINE reserved for the top-10 of the ladder (chance floors are printed).

    soft=True is VERNIER-23's un-run Part 3: no certification threshold at all --
    the pulsar term enters MARGINALISED over the distance prior, weighted by the
    fringe coherence R_a of the posterior.
    """
    rng = np.random.default_rng(seed)
    bins = [make_bin(rng, arr, f, arm, gain=gain) for f in band]
    floors = [fstat_floor(rng, B) for B in bins]
    data = [inject(rng, B, null=null) for B in bins]
    found = [[] for _ in bins]                    # grid indices found, per bin
    sub = [set() for _ in bins]                   # true source indices removed
    det, hist = [], []
    gen = spl = blend = 0
    for t in range(max_turn):
        best = None
        for b, B in enumerate(bins):
            d = data[b]
            if found[b]:
                A = np.concatenate([np.stack([B["Fpw"][:, j], B["Fcw"][:, j]], 1)
                                    for j in found[b]], axis=1)
                Mm = A.T @ A
                try:
                    coef = np.linalg.solve(
                        Mm + 1e-10 * np.eye(Mm.shape[0]) * np.trace(Mm), A.T @ d)
                    r = d - A @ coef
                except np.linalg.LinAlgError:
                    continue
            else:
                r = d
            Fv = S.fstat(r, B["Fpw"], B["Fcw"])
            j = int(np.argmax(Fv))
            exc = float(Fv[j]) - (floors[b]["floor"] + floors[b]["err"])
            if exc > 0 and (best is None or exc > best[0]):
                best = (exc, b, j, float(Fv[j]))
        if best is None:
            break
        exc, b, j, Fmax = best
        B = bins[b]
        which, kind = -1, "blend"
        if B["Gt"].shape[0]:
            ca = B["Gt"] @ B["G"][j]
            i = int(np.argmax(ca))
            if np.arccos(np.clip(ca[i], -1, 1)) <= TOL_RAD and i < TOP_TRUE:
                if i in sub[b]:
                    kind, spl = "split", spl + 1
                else:
                    kind, which, gen = "genuine", i, gen + 1
            else:
                blend += 1
        found[b].append(j)
        if which < 0:
            # a blend or a split: CLEAN removes it and the wheel keeps turning
            # (there is nothing to certify against, so no certification is run)
            hist.append(dict(turn=t + 1, bin=b, f=B["f"], F=Fmax, kind=kind,
                             which=-1, n_cert=0, n_true=0, n_cert_null=0,
                             f_cert=0.0, f_cert_true=0.0, f_soft=0.0,
                             rho_P_med=0.0, rho_P_max=0.0, rho_arr=0.0,
                             gain_sub=hist[-1]["gain_sub"] if hist else 0.0,
                             gain_cert=hist[-1]["gain_cert"] if hist else 0.0))
            continue
        det.append((B, which, set(sub[b])))
        C = certify(rng, arr, det, soft=soft)
        Cn = certify(rng, arr, det, scramble=True, soft=soft)
        sub[b].add(which)
        # depth: confusion power removed, in amplitude dex, band- and array-median
        gs = []
        for Bb, sb in zip(bins, sub):
            n_all = 1.0 + (Bb["rho0"] ** 2).sum(axis=1)
            left = [i for i in range(Bb["M"]) if i not in sb]
            n_now = (1.0 + (Bb["rho0"][:, left] ** 2).sum(axis=1) if left
                     else np.ones_like(n_all))
            gs.append(np.median(np.log10(n_all / n_now)))
        g_sub = float(0.5 * np.mean(gs))
        fc = C["f_soft"] if soft else C["f_cert"]
        g_cert = float(0.5 * np.log10(1.0 + R_AMP ** 2 * fc))
        hist.append(dict(turn=t + 1, bin=b, f=B["f"], F=Fmax, kind=kind,
                         which=which, n_cert=C["n_cert"], n_true=C["n_true"],
                         n_cert_null=Cn["n_cert"], f_cert=C["f_cert"],
                         f_cert_true=C["f_cert_true"], f_cert_null=Cn["f_cert"],
                         f_soft=C["f_soft"], rho_P_med=C["rho_P_med"],
                         rho_P_p90=C["rho_P_p90"], rho_P_max=C["rho_P_max"],
                         rho_arr=C["rho_arr"], lnK=C["lnK"], dL=C["dL"],
                         S_med=C["S_med"], S_max=C["S_max"],
                         floor_inc=C["floor_inc"], floor_noc=C["floor_noc"],
                         n_cert_nocomb=C["n_cert_nocomb"],
                         comb_med=C["comb_med"], hard_bar=C["hard_bar"],
                         zero_frac=C["zero_frac"], kind_floor=C["kind"],
                         gain_sub=g_sub, gain_cert=g_cert))
    g_tot = (hist[-1]["gain_sub"] + hist[-1]["gain_cert"]) if hist else 0.0
    return dict(hist=hist, genuine=gen, split=spl, blend=blend, turns=len(hist),
                gain=g_tot, rank=rank_of(g_tot), floors=floors, bins=bins)


# ============================================================ Part 3 machinery
def dvec(P, S_, L=None):
    """
    Vectorised VERNIER-15 data vector: Earth block for every pulsar, plus the
    pulsar-term block. Column 6 of S_ is each source's own frequency, which is
    what makes the pulsar-term block independent information.
    """
    th, ph = S_[:, 0], S_[:, 1]
    A, psi, ci, pha = S_[:, 2], S_[:, 3], S_[:, 4], S_[:, 5]
    om = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
    z = np.array([0.0, 0.0, 1.0])
    u = np.cross(z[None, :], om)
    nn = np.linalg.norm(u, axis=1, keepdims=True)
    u = np.where(nn > 1e-12, u / np.maximum(nn, 1e-30), np.array([1.0, 0.0, 0.0]))
    v = np.cross(om, u)
    pu, pv = P @ u.T, P @ v.T
    dd = np.maximum(1.0 + P @ om.T, OMC_MIN)          # regularised, see antenna_reg
    Fp, Fc = (pu ** 2 - pv ** 2) / (2 * dd), (pu * pv) / dd
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    Fp, Fc = Fp * c2 + Fc * s2, -Fp * s2 + Fc * c2
    E = (Fp * ((1 + ci ** 2) / 2) + 1j * Fc * ci) * (A * np.exp(1j * pha))[None, :]
    earth = E.sum(axis=1)
    if L is None:
        return earth
    omc = np.clip(1.0 - P @ om.T, 5e-2, None)
    tau = L[:, None] * omc / C_PC_S
    return np.concatenate([earth,
                           (E * np.exp(-2j * np.pi * S_[:, 6][None, :] * tau)).sum(axis=1)])


def separability(rng, n_psr, M, f, sigma_L_pc=PRIOR_PC, use_pterm=True,
                 marginalise=True, sigma_toa=SIG_TOA, eps=1e-5, tol=1e-10,
                 rho_fixed=None):
    """
    THE DECIDING TEST -- REALIZED separation, not degree-of-freedom counting.

    M sources share ONE bin at frequency f, with the Arm B strain ladder, seen by
    a real array. The Fisher is built over [source params] + [pulsar distances];
    the distances are then MARGINALISED over their 10 pc priors by the Schur
    complement (VERNIER-23's construction). Two things are measured, not counted:

      CONDITIONING -- the smallest eigenvalue of the correlation-normalised
                      marginalised source block, and its condition number.
      PARAMETER ERROR -- sigma(sky) and sigma(A)/A per source from the INVERSE
                      marginalised Fisher, at the source's real SNR.

    A source counts as RESOLVED when sigma_sky is inside half the typical
    nearest-neighbour separation sqrt(4 pi / M) / 2 AND sigma(A)/A < 0.5 -- the
    VERNIER-15 criterion, errors inside tolerance rather than rank alone.

    The realized frontier is the largest M at which the LOUDEST source is still
    resolved (and, separately, at which half the bin is).
    """
    P = V.iso_dirs(n_psr, rng)
    w = 10.0 ** rng.normal(0.0, SIG_LOG_W, n_psr)
    w = w / np.sqrt(np.mean(w ** 2))
    L = 10 ** rng.uniform(np.log10(300.0), np.log10(5000.0), n_psr)
    R2, S_sea, diag = bin_sources(rng, f, "B", n_bright=M)
    amp = np.sqrt(R2 / R2[0])                      # relative amplitudes
    # source frequencies spread across the bin -- what "same bin" allows
    fk = f + rng.uniform(-0.5, 0.5, M) * DF
    S_ = np.stack([np.arccos(rng.uniform(-1, 1, M)), rng.uniform(0, 2 * np.pi, M),
                   amp, rng.uniform(0, np.pi, M), rng.uniform(-1, 1, M),
                   rng.uniform(0, 2 * np.pi, M), fk], axis=1)
    base = dvec(P, S_, L if use_pterm else None)
    npar = 7 if use_pterm else 6
    ncol = M * npar + (n_psr if (use_pterm and marginalise) else 0)
    J = np.zeros((2 * base.size, ncol))
    for i in range(M):
        for k in range(npar):
            sp = S_.copy()
            st = eps * f if k == 6 else eps * (abs(sp[i, k]) + 1.0)
            sp[i, k] += st
            dd = (dvec(P, sp, L if use_pterm else None) - base) / st
            J[:, i * npar + k] = np.concatenate([dd.real, dd.imag])
    if use_pterm and marginalise:
        for a in range(n_psr):
            Lv = L.copy(); st = eps * (abs(Lv[a]) + 1.0); Lv[a] += st
            dd = (dvec(P, S_, Lv) - base) / st
            J[:, M * npar + a] = np.concatenate([dd.real, dd.imag])
    # noise scale: set so the LOUDEST source has the SNR the Arm B budget gives
    S_tot = S_white(sigma_toa) + S_sea
    rho1 = np.sqrt(KAPPA * n_psr * R2[0] * T_OBS / (2.0 * S_tot))
    if rho_fixed is not None:                 # identifiability-limited regime:
        rho1 = float(rho_fixed)               # take loudness out of the question
    e1 = dvec(P, S_[:1], L if use_pterm else None)
    sig2 = float(np.sum(np.abs(e1) ** 2)) / max(rho1 ** 2, 1e-30)
    F = (J.T @ J) / (sig2 / 2.0)          # real and imag each carry sig2/2
    ns = M * npar
    if use_pterm and marginalise:
        Fdd = F[ns:, ns:] + np.eye(n_psr) / sigma_L_pc ** 2
        Fsd = F[ns:, :ns].T
        try:
            Fm = F[:ns, :ns] - Fsd @ np.linalg.solve(Fdd, Fsd.T)
        except np.linalg.LinAlgError:
            return dict(M=M, f=f, rho1=rho1, resolved=0, top_ok=False,
                        emin=0.0, cond=np.inf, sig_sky=np.inf, sig_A=np.inf,
                        diag=diag)
    else:
        Fm = F[:ns, :ns]
    # invert in the CORRELATION metric: the frequency parameter carries a Fisher
    # entry ~1e24 in Hz^-2 while an amplitude carries ~20, so a trace-scaled ridge
    # on the raw matrix swamps every other parameter (it did, and reported
    # sigma_A/A = 1e-6 at rho = 3).
    dg = np.sqrt(np.clip(np.diag(Fm), 1e-300, None))
    Fc = Fm / np.outer(dg, dg)
    ev = np.linalg.eigvalsh(Fc)
    emin = float(max(ev.min(), 0.0) / max(ev.max(), 1e-300))
    try:
        Cc = np.linalg.inv(Fc + 1e-12 * np.eye(ns))
    except np.linalg.LinAlgError:
        return dict(M=M, f=f, rho1=rho1, resolved=0, top_ok=False, emin=emin,
                    cond=np.inf, sig_sky=np.inf, sig_A=np.inf, diag=diag)
    C = Cc / np.outer(dg, dg)
    d = np.clip(np.diag(C), 0.0, None)
    sth = np.sqrt(d[0::npar]); sph = np.sqrt(d[1::npar]) * np.sin(S_[:, 0])
    sky = np.sqrt(sth ** 2 + sph ** 2)
    sA = np.sqrt(d[2::npar]) / np.maximum(S_[:, 2], 1e-30)
    sep = np.sqrt(4 * np.pi / max(M, 2))
    ok = (sky < 0.5 * sep) & (sA < 0.5)
    return dict(M=M, f=f, rho1=float(rho1), resolved=int(ok.sum()),
                top_ok=bool(ok[0]), emin=emin,
                cond=float(ev.max() / max(ev.min(), 1e-300)),
                sig_sky=float(np.median(sky)), sig_A=float(np.median(sA)),
                sig_sky_top=float(sky[0]), sig_A_top=float(sA[0]),
                sep=float(sep), diag=diag)


def chirp_precision(rng, n_psr, f, competitors=True, sigma_toa=SIG_TOA,
                    log10_mc=MC_FID, n_real=24):
    """
    RE-CUT of VERNIER-23's chirp-mass-from-line-width spin-off, WITH the audited
    in-bin competitors present.

    The pulsar term of pulsar a arrives at the retarded frequency
    f_P = f (1 + tau_a/tc)^(-3/8), so for tau << tc the offset is
    delta_a = (3/8) f tau_a / tc, and tc ~ Mc^(-5/3) makes delta ~ Mc^(5/3).
    Each pulsar measures its own f_P to sigma_f = 1/(2 pi T rho_P,a), so

        1/sigma(ln Mc)^2 = sum_a [ (5/3) delta_a 2 pi T rho_P,a ]^2

    and sigma(log10 Mc) = sigma(ln Mc)/ln 10. With competitors, rho_P,a is the
    lambda_eff-consistent confusion budget of the loudest in-bin source; without,
    it is the same source against white + sea alone. Resolvability is reported
    beside it: the line is only resolvable if its spread exceeds one Fourier bin.
    """
    out = {}
    for tag, use in (("coh", "B"), ("gauss", "Am"), ("iso", "iso")):
        vals, wid = [], []
        for _ in range(n_real):
            arr = make_array(rng, n_psr, sigma_toa=sigma_toa)
            B = make_bin(rng, arr, f, use)
            rho = rho_budget(B, 0, set()) * R_AMP
            omc = B["omc"][:, 0]
            tau = V.tau_of(arr["L"], omc)
            tc = float(V.t_merge(f, log10_mc))
            delta = (3.0 / 8.0) * f * tau / tc
            fis = np.sum(((5.0 / 3.0) * delta * 2 * np.pi * T_OBS * rho) ** 2)
            vals.append(1.0 / (np.sqrt(max(fis, 1e-300)) * np.log(10.0)))
            wid.append((delta.max() - delta.min()) * T_OBS)
        out[tag] = (float(np.median(vals)), float(np.median(wid)))
    return out


def rank_marg(rng, n_psr, M, sigma_L_pc, f, tol=1e-10, eps=1e-5):
    """
    Rank of the SOURCE block after MARGINALISING the pulsar distances over their
    priors (VERNIER-23's Schur complement), vectorised so it runs at N = 116.
    Counting prediction: M_max = 3N/7.
    """
    P = V.iso_dirs(n_psr, rng)
    L = 10 ** rng.uniform(np.log10(300.0), np.log10(5000.0), n_psr)
    S_ = np.stack([np.arccos(rng.uniform(-1, 1, M)), rng.uniform(0, 2 * np.pi, M),
                   np.ones(M), rng.uniform(0, np.pi, M), rng.uniform(-1, 1, M),
                   rng.uniform(0, 2 * np.pi, M), np.full(M, f)], axis=1)
    base = dvec(P, S_, L)
    npar = 7
    J = np.zeros((2 * base.size, M * npar + n_psr))
    for i in range(M):
        for k in range(npar):
            sp = S_.copy()
            st = eps * f if k == 6 else eps * (abs(sp[i, k]) + 1.0)
            sp[i, k] += st
            dd = (dvec(P, sp, L) - base) / st
            J[:, i * npar + k] = np.concatenate([dd.real, dd.imag])
    for a in range(n_psr):
        Lv = L.copy(); st = eps * (abs(Lv[a]) + 1.0); Lv[a] += st
        dd = (dvec(P, S_, Lv) - base) / st
        J[:, M * npar + a] = np.concatenate([dd.real, dd.imag])
    F = J.T @ J
    ns = M * npar
    Fdd = F[ns:, ns:] + np.eye(n_psr) / sigma_L_pc ** 2
    Fss, Fsd = F[:ns, :ns], F[ns:, :ns].T
    try:
        Fm = Fss - Fsd @ np.linalg.solve(Fdd, Fsd.T)
    except np.linalg.LinAlgError:
        return 0, ns
    dg = np.sqrt(np.clip(np.diag(Fm), 1e-300, None))
    Fm = Fm / np.outer(dg, dg)
    ev = np.linalg.eigvalsh(Fm); ev = ev[ev > 0]
    return (int(np.sum(ev > tol * ev.max())) if ev.size else 0), ns


# ============================================================ deliverables
ARMS = [("B", "PHYSICAL SKY: audited Pareto census, 300 explicit + faint sea"),
        ("Am", "zero-competitor control at Arm B's own lambda_fore"),
        ("A", "SUPERSEDED one-source-per-bin artefact sky (kept, relabelled)")]


def chance_match(M, tol=TOL_RAD, top=TOP_TRUE):
    p1 = (1.0 - np.cos(tol)) / 2.0
    return 1.0 - (1.0 - p1) ** min(M, top), 1.0 - (1.0 - p1) ** M


def hh_of(res):
    return [h for h in res["hist"] if h["which"] >= 0]


# ------------------------------------------------------------------ 0. census
def exp0_census(out):
    print("=" * 100)
    print("DELIVERABLE 0 -- THE AUDITED SKY, REBUILT")
    print("=" * 100)
    print(f"  T = 20 yr, 1/T = {DF*1e9:.3f} nHz, band {BAND[0]*1e9:.2f}-"
          f"{BAND[-1]*1e9:.2f} nHz in {len(BAND)} bins; array {N_PSR} pulsars, "
          f"sigma_toa(med) {SIG_TOA*1e9:.0f} ns")
    print(f"  strain function: Pareto gamma_eff = {GAMMA_EFF}, h = u^(-1/gamma), "
          f"normalised to the bin total; h(Mc,f,D) link -> h^2/bin ~ f^(-7/3)")
    print()
    # the f^(-7/3) check, measured not asserted
    ff = np.logspace(np.log10(5e-9), np.log10(5e-8), 40)
    sl = np.polyfit(np.log10(ff), np.log10(h2_per_bin(ff)), 1)[0]
    print(f"  h^2 per bin scaling: measured slope {sl:+.4f} against the required "
          f"-7/3 = {-7/3:+.4f}   [{'OK' if abs(sl + 7/3) < 1e-6 else 'MISMATCH'}]")
    print()
    print(f"  {'f nHz':>6} {'N_bin':>9} {'audit':>9} {'lam_fore':>9} {'audit':>7} "
          f"{'N_eff':>7} {'N<0.5dex':>9} {'h1':>10} {'log10 Mc':>9} "
          f"{'rho_arr':>8} {'rhoP med':>9} {'rhoP p90':>9}")
    rng = np.random.default_rng(24)
    AUD = {10: 6.9e4, 15: 1.6e4, 20: 5.4e3}
    rows = []
    for fn in (10.0, 15.0, 20.0):
        f = fn * 1e-9
        lam, ne, nh, h1, ra, rm, rp = [], [], [], [], [], [], []
        for _ in range(200):                       # census diagnostics: cheap
            _r2, _sea, d = bin_sources(rng, f, "B")
            lam.append(d["lam_fore"]); ne.append(d["n_eff"]); nh.append(d["n_halfdex"])
            h1.append(np.sqrt(_r2[0] / bin_power(f) * h2_per_bin(f)))
        for _ in range(10):                        # SNRs: needs the full array
            arr = make_array(rng)
            B = make_bin(rng, arr, f, "B")
            rho = rho_budget(B, 0, set())
            ra.append(np.sqrt(np.sum(rho ** 2)))
            rm.append(np.median(R_AMP * rho)); rp.append(np.percentile(R_AMP * rho, 90))
        h1m = float(np.median(h1))
        # invert h = 2 (T_sun Mc)^(5/3) (pi f)^(2/3) / (D/c) at D = 1 Gpc
        dl = 1e9 * 3.0856775814913673e16 / 299792458.0
        mc = np.log10(((h1m * dl / 2.0) / (np.pi * f) ** (2 / 3.)) ** 0.6
                      / 4.925490947e-6)
        rows.append([fn, float(np.median([d for d in ne])), float(np.median(lam)),
                     h1m, mc, float(np.median(ra)), float(np.median(rm)),
                     float(np.median(rp))])
        print(f"  {fn:6.1f} {N_c(f):9.3g} {AUD[fn]:9.3g} {np.median(lam):9.3f} "
              f"{LAM_AUDIT:7.2f} {np.median(ne):7.1f} {np.median(nh):9.0f} "
              f"{h1m:10.3e} {mc:9.2f} {np.median(ra):8.2f} {np.median(rm):9.3f} "
              f"{np.median(rp):9.3f}")
    print()
    print(f"  audited N_eff range {N_EFF_AUD[0]:.0f}-{N_EFF_AUD[1]:.0f}; the")
    print("  participation ratio 1/sum(p^2) reproduces its LOWER end, and lambda_fore")
    print("  comes out FLAT in N_bin across 6.9e4 -> 5.4e3, as the audit requires --")
    print("  that flatness is structural for gamma < 2 and is why no one-source bin")
    print("  exists at any frequency.")
    out["cen"] = np.array(rows)


# ------------------------------------------------- 1. THE DECIDING TEST
def exp1_frontier(out, Ms=(2, 5, 10, 17, 25, 50), n_real=4, rho_hi=30.0):
    print()
    print("=" * 100)
    print("DELIVERABLE 1 -- THE DECIDING TEST: REALIZED SEPARATION, DISTANCES MARGINALISED")
    print("=" * 100)
    print("  M sources share ONE bin on the Arm B ladder; Fisher over source parameters")
    print("  + pulsar distances, distances MARGINALISED over 10 pc priors (Schur).")
    print("  RESOLVED = sigma_sky < half the typical separation AND sigma(A)/A < 0.5 --")
    print("  errors inside tolerance, the VERNIER-15 criterion, not rank alone.")
    print("  Run TWICE: at the census's own loudness, and at a fixed rho_1 = "
          f"{rho_hi:.0f} which")
    print("  takes loudness out of the question and leaves pure identifiability.")
    print("  CAVEAT stated up front: a Fisher error is LOCAL. sigma_sky ~ 1e-4 rad is a")
    print("  within-fringe width, not a global uncertainty -- the pulsar-term posterior")
    print("  is multi-modal in exactly the way deliverables 2-3 measure.")
    print("  rho_1 here is the loudest source against white + faint sea with the other")
    print("  M-1 JOINTLY FITTED, so it is larger than deliverable 0's rho_arr, which")
    print("  charges the whole bright end to the noise budget.")
    rows = []
    for f in F_REPORT:
        print()
        print(f"  --- f = {f*1e9:.0f} nHz, N_bin = {N_c(f):.2g}, "
              f"N_eff ~ 8-10 " + "-" * 40)
        print(f"  {'M':>4} | {'rho_1':>6} {'res':>7} {'top':>5} {'sA_top':>8} "
              f"{'emin':>9} {'earth':>6} | {'rho_1=%.0f' % rho_hi:>8} {'res':>7} "
              f"{'top':>5} {'sA_top':>8} {'emin':>9} {'earth':>6}")
        for M in Ms:
            cells = []
            for rf in (None, rho_hi):
                A = [separability(np.random.default_rng(400 + 9 * s + int(f * 1e10)),
                                  N_PSR, M, f, rho_fixed=rf) for s in range(n_real)]
                E = [separability(np.random.default_rng(400 + 9 * s + int(f * 1e10)),
                                  N_PSR, M, f, use_pterm=False, rho_fixed=rf)
                     for s in range(n_real)]
                g = lambda L, k: float(np.median([x[k] for x in L]))
                cells.append([g(A, "rho1"), g(A, "resolved"),
                              float(np.mean([x["top_ok"] for x in A])),
                              g(A, "sig_A_top"), g(A, "emin"), g(E, "resolved")])
            rows.append([f * 1e9, M] + cells[0] + cells[1])
            c0, c1 = cells
            print(f"  {M:4d} | {c0[0]:6.2f} {c0[1]:4.1f}/{M:<2d} {c0[2]:5.2f} "
                  f"{c0[3]:8.3f} {c0[4]:9.1e} {c0[5]:6.1f} | {c1[0]:8.1f} "
                  f"{c1[1]:4.1f}/{M:<2d} {c1[2]:5.2f} {c1[3]:8.3f} {c1[4]:9.1e} "
                  f"{c1[5]:6.1f}")
    rows = np.array(rows)
    print()
    print("  REALIZED FRONTIER -- largest M at which the TOP source resolves, and at")
    print("  which at least HALF the bin resolves:")
    print(f"  {'f nHz':>6} | {'census loudness':>28} | "
          f"{'rho_1 = %.0f' % rho_hi:>28}")
    print(f"  {'':>6} | {'top':>6} {'half':>6} {'earth-half':>12} | "
          f"{'top':>6} {'half':>6} {'earth-half':>12}")
    for f in F_REPORT:
        sel = rows[rows[:, 0] == f * 1e9]
        # row = [f, M, (rho1,res,top,sA,emin,earth)_census, (...)_rho_hi]
        pick = lambda cond: (int(max([r[1] for r in sel if cond(r)]))
                             if any(cond(r) for r in sel) else 0)
        a1 = pick(lambda r: r[4] >= 0.5)
        a2 = pick(lambda r: r[3] >= 0.5 * r[1])
        a3 = pick(lambda r: r[7] >= 0.5 * r[1])
        b1 = pick(lambda r: r[10] >= 0.5)
        b2 = pick(lambda r: r[9] >= 0.5 * r[1])
        b3 = pick(lambda r: r[13] >= 0.5 * r[1])
        print(f"  {f*1e9:6.1f} | {a1:6d} {a2:6d} {a3:12d} | "
              f"{b1:6d} {b2:6d} {b3:12d}")
    print()
    print(f"  counting bounds at N = {N_PSR}: Earth 2N/6 = {2*N_PSR/6:.0f}, "
          f"marginalised 3N/7 = {3*N_PSR/7:.0f}, distances known 4N/7 = "
          f"{4*N_PSR/7:.0f}.")
    print("  VERNIER-23 measured M_max = 17-20 marginalised and 8 Earth-only at N = 40")
    print("  (its bounds were 17 and 13), i.e. it sat AT its counting bound.")
    out["front"] = rows


# ------------------------------------------------------------------ 2. wheel
def exp2_wheels(out, n_real=4):
    print()
    print("=" * 100)
    print("DELIVERABLE 2 -- THE WHEEL IN THE CLEAN BAND, ON THE CORRECTED SKY")
    print("=" * 100)
    print("  per turn: search all bins (F-stat, 6000-pt sky grid, per-bin floor =")
    print("  criterion-v2.2 on 150 pure-noise nulls) -> certify against the ACCUMULATED")
    print("  comb of all sources found so far (per-pulsar v2.2 floor, incoherent-comb")
    print("  null, 100 draws each) -> subtract with a joint amplitude refit -> repeat.")
    print()
    print(f"  {'arm':>4} {'turns':>6} {'gen':>5} {'split':>6} {'blend':>6} "
          f"{'rho_arr':>8} {'rhoP med':>9} {'rhoP p90':>9} {'rhoP max':>9} "
          f"{'CERT':>5} {'true':>5} {'null':>5} {'nocomb':>7} {'f_cert':>7} "
          f"{'lnK':>6} {'S max':>8} {'dex':>8} {'rank':>5}")
    rec = {}
    for arm, desc in ARMS:
        acc = []
        for s in range(n_real):
            arr = make_array(np.random.default_rng(1000 + s))
            acc.append(wheel(2000 + 97 * s, arr, arm))
        hh = [h for r in acc for h in hh_of(r)]
        med = lambda key, d=0.0: (np.median([h[key] for h in hh]) if hh else d)
        print(f"  {arm:>4} {np.median([r['turns'] for r in acc]):6.1f} "
              f"{np.mean([r['genuine'] for r in acc]):5.1f} "
              f"{np.mean([r['split'] for r in acc]):6.1f} "
              f"{np.mean([r['blend'] for r in acc]):6.1f} "
              f"{med('rho_arr'):8.2f} {med('rho_P_med'):9.3f} "
              f"{med('rho_P_p90'):9.3f} {med('rho_P_max'):9.2f} "
              f"{med('n_cert'):5.1f} {med('n_true'):5.1f} {med('n_cert_null'):5.1f} "
              f"{med('n_cert_nocomb'):7.1f} {med('f_cert'):7.3f} {med('lnK'):6.2f} "
              f"{med('S_max'):8.3f} "
              f"{np.median([r['gain'] for r in acc]):+8.4f} "
              f"{int(np.median([r['rank'] for r in acc])):5d}")
        rec[arm] = acc
    print()
    print("  per-turn detail, one realisation per arm:")
    for arm, acc in rec.items():
        for h in acc[0]["hist"]:
            print(f"    {arm:>3} t{h['turn']} f={h['f']*1e9:5.2f} nHz F={h['F']:8.2f} "
                  f"{h['kind']:<8} rho_P med {h['rho_P_med']:.3f} max "
                  f"{h['rho_P_max']:7.2f}  CERT {h['n_cert']:>3}/{h['n_true']:>3}"
                  f"/{h['n_cert_null']:>3} (nocomb {h.get('n_cert_nocomb', 0):>3})  "
                  f"S_max {h.get('S_max', 0.0):8.3f} vs bar "
                  f"{max(h.get('lnK', 0.0), h.get('floor_inc', 0.0)):.2f}  "
                  f"dex {h['gain_sub']+h['gain_cert']:+.4f}")
    out["w_arm"] = np.array([a for a, _ in ARMS])
    out["w_turns"] = np.array([np.median([r["turns"] for r in rec[a]]) for a, _ in ARMS])
    out["w_gain"] = np.array([np.median([r["gain"] for r in rec[a]]) for a, _ in ARMS])
    out["w_rank"] = np.array([np.median([r["rank"] for r in rec[a]]) for a, _ in ARMS])
    return rec


def exp3_controls(out, n_real=4, n_null_real=8):
    print()
    print("=" * 96)
    print("DELIVERABLE 3 -- CONTROLS")
    print("=" * 96)
    fired, cs, ct, rx, sm = 0, [], [], [], []
    for s in range(n_real):
        arr = make_array(np.random.default_rng(1000 + s))
        r = wheel(7000 + s, arr, "M1")
        fired += int(r["genuine"] > 0)
        for h in hh_of(r):
            cs.append(h["n_cert"]); ct.append(h["n_true"])
            rx.append(h["rho_P_max"]); sm.append(h["S_max"])
    f_ = lambda v: (np.median(v) if v else 0.0)
    print(f"  M=1 control (one source per bin holding the WHOLE bin, no remainder):")
    print(f"      detected {fired}/{n_real};  CERT {f_(cs):.1f} of {N_PSR} "
          f"(true {f_(ct):.1f});  rho_P max {f_(rx):.2f};  best fringe contrast "
          f"S_max {f_(sm):.2f} nat")
    for arm in ("B", "A"):
        nf, tt, cc = 0, [], 0
        for s in range(n_null_real):
            arr = make_array(np.random.default_rng(1000 + s))
            r = wheel(8000 + s, arr, arm, null=True)
            nf += int(r["turns"] > 0); tt.append(r["turns"])
            cc += sum(h["n_cert"] for h in r["hist"])
        exp = 1.0 - (1.0 - DET_ALPHA) ** len(BAND)
        print(f"  pure-noise null through the whole wheel, arm {arm}: fired "
              f"{nf}/{n_null_real} (expected {exp:.2f} at alpha = {DET_ALPHA} per bin "
              f"x {len(BAND)} bins), mean {np.mean(tt):.1f} turns, "
              f"{cc} certifications")
    print("  per-turn SOURCE-LEVEL scrambled null: the `null` column of every")
    print("  certification table (same geometry, distance displaced U(-500,500) pc).")
    for arm in ("B", "A"):
        M = bin_sources(np.random.default_rng(0), BAND[0], arm)[0].size
        p10, pany = chance_match(M)
        print(f"  chance-match floor, arm {arm} (M = {M:3d}/bin): top-10 {p10:.3f}, "
              f"any {pany:.3f}  -> GENUINE scored against the top-10 only")
    out["ctrl"] = np.array([fired, n_real, f_(cs), f_(ct), f_(rx), f_(sm)])


def exp4_loudness(out, gains=(1.0, 3.0, 10.0, 30.0, 100.0), n_real=3):
    print()
    print("=" * 96)
    print("DELIVERABLE 4 -- LOUDNESS SWEEP: what would the bright end have to weigh?")
    print("=" * 96)
    print("  `gain` multiplies the bright end and is removed from the Gaussian")
    print("  remainder, so the bin total is conserved -- a redistribution of the same")
    print("  GWB power between resolved and unresolved, not extra signal.")
    print()
    print(f"  {'arm':>4} {'gain':>6} {'rho_arr':>9} {'rhoP med':>9} {'rhoP p90':>9} "
          f"{'rhoP max':>9} {'turns':>6} {'gen':>5} {'CERT':>6} {'true':>6} "
          f"{'null':>6} {'nocomb':>7} {'f_cert':>8} {'hardbar':>8} {'dex':>8} "
          f"{'rank':>5}")
    rows = []
    for arm in ("B", "Am"):
        for g in gains:
            acc = []
            for s in range(n_real):
                arr = make_array(np.random.default_rng(1000 + s))
                acc.append(wheel(9000 + s + int(100 * g), arr, arm, gain=g))
            hh = [h for r in acc for h in hh_of(r)]
            m = lambda k: (np.median([h[k] for h in hh]) if hh else 0.0)
            row = [1 if arm == "B" else 2, g, m("rho_arr"), m("rho_P_med"),
                   m("rho_P_p90"), m("rho_P_max"),
                   np.median([r["turns"] for r in acc]),
                   np.mean([r["genuine"] for r in acc]),
                   m("n_cert"), m("n_true"), m("n_cert_null"),
                   m("n_cert_nocomb"), m("f_cert"),
                   m("hard_bar"), np.median([r["gain"] for r in acc]),
                   int(np.median([r["rank"] for r in acc]))]
            rows.append(row)
            print(f"  {arm:>4} {g:6.0f} {row[2]:9.2f} {row[3]:9.3f} {row[4]:9.3f} "
                  f"{row[5]:9.2f} {row[6]:6.1f} {row[7]:5.1f} {row[8]:6.1f} "
                  f"{row[9]:6.1f} {row[10]:6.1f} {row[11]:7.1f} {row[12]:8.3f} "
                  f"{row[13]:8.3f} {row[14]:+8.4f} {int(row[15]):5d}")
    rows = np.array(rows)
    print()
    print("  bright-end power as a fraction of the measured bin total (>1 means the")
    print("  sweep has left the GWB behind and is answering a counterfactual):")
    for arm in ("B", "Am"):
        fr = [float(bin_sources(np.random.default_rng(0), BAND[0], arm,
                                gain=g)[0].sum() / bin_power(BAND[0]))
              for g in gains]
        print(f"      arm {arm}: " + "  ".join(f"x{g:.0f} -> {v:.2f}"
                                               for g, v in zip(gains, fr)))
    print()
    print("  Loudness needed to reach rho_P = 2, at the bin's own confusion:")
    for arm, tag in (("B", 1), ("Am", 2)):
        sel = rows[(rows[:, 0] == tag) & (rows[:, 3] > 0)]
        if sel.shape[0] < 2:
            print(f"    arm {arm}: not measurable (no certifiable turn at any loudness)")
            continue
        for col, nm in ((3, "median"), (4, "p90"), (5, "best")):
            sl = np.polyfit(np.log10(sel[:, 1]), np.log10(sel[:, col]), 1)
            need = 10 ** ((np.log10(2.0) - sl[1]) / sl[0])
            print(f"    arm {arm}, {nm:>6} pulsar: rho_P ~ gain^{sl[0]:.3f} -> needs "
                  f"gain x{need:.3g} of the bright end "
                  f"({0.5*np.log10(max(need,1e-30)):+.2f} dex in amplitude)")
    out["loud"] = rows
    # --- timing precision, at the census's own loudness
    print()
    print("  Arm B (the physical sky) at its own loudness, versus timing precision:")
    print(f"  {'sigma_toa':>10} {'rho_arr':>9} {'rhoP med':>9} {'rhoP p90':>9} "
          f"{'rhoP max':>9} {'turns':>6} {'gen':>5} {'CERT':>6} {'S max':>8} "
          f"{'lnK':>6} {'dex':>8} {'rank':>5}")
    trows = []
    for st in (100e-9, 30e-9, 10e-9):
        acc = []
        for s in range(n_real):
            arr = make_array(np.random.default_rng(1000 + s), sigma_toa=st)
            acc.append(wheel(13000 + s + int(st * 1e10), arr, "B"))
        hh = [h for r in acc for h in hh_of(r)]
        m = lambda k: (np.median([h[k] for h in hh]) if hh else 0.0)
        row = [st * 1e9, m("rho_arr"), m("rho_P_med"), m("rho_P_p90"),
               m("rho_P_max"), np.median([r["turns"] for r in acc]),
               np.mean([r["genuine"] for r in acc]), m("n_cert"),
               max([h["S_max"] for h in hh], default=0.0), m("lnK"),
               np.median([r["gain"] for r in acc]),
               int(np.median([r["rank"] for r in acc]))]
        trows.append(row)
        print(f"  {row[0]:9.0f}n {row[1]:9.2f} {row[2]:9.3f} {row[3]:9.3f} "
              f"{row[4]:9.2f} {row[5]:6.1f} {row[6]:5.1f} {row[7]:6.1f} "
              f"{row[8]:8.3f} {row[9]:6.2f} {row[10]:+8.4f} {int(row[11]):5d}")
    out["timing"] = np.array(trows)


def exp5_part3(out, n_real=3):
    print()
    print("=" * 100)
    print("DELIVERABLE 5 -- ARM B: VERNIER-23's PART 3, PREMISE RESTORED")
    print("=" * 96)
    print("  Separate with MARGINALISED pulsar terms: no certification threshold at")
    print("  all, the pulsar term weighted by the fringe coherence R_a of the distance")
    print("  posterior under its 10 pc prior.")
    print()
    f = float(BAND[len(BAND) // 2])
    print(f"  (a) counting bound: marginalised Fisher rank at the real array size (N = {N_PSR}, "
          f"sigma_L = {PRIOR_PC:.0f} pc, f = {f*1e9:.2f} nHz):")
    rng = np.random.default_rng(24)
    Ms, rk = [], []
    for M in (10, 20, 30, 40, 50, 60):
        r, ns = rank_marg(rng, N_PSR, M, PRIOR_PC, f)
        Ms.append(M); rk.append(r / ns)
        print(f"      M = {M:3d}   source params {ns:4d}   rank/params = {r/ns:6.3f}"
              + ("   FULL" if r == ns else ""))
    print(f"      counting prediction M_max = 3N/7 = {3*N_PSR/7:.1f}; "
          f"Arm B holds {int(N_c(f)):d} sources in this bin")
    print()
    print("  (b) the wheel in Arm B, hard certification vs marginalised:")
    for nm, sf in (("hard cert (v2.2)", False), ("MARGINALISED", True)):
        acc = []
        for s in range(n_real):
            arr = make_array(np.random.default_rng(1000 + s))
            acc.append(wheel(11000 + s, arr, "B", soft=sf))
        hh = [h for r in acc for h in hh_of(r)]
        got = [r for r in acc if r["genuine"] > 0]
        fs = np.median([h["f_soft"] for h in hh]) if hh else 0.0
        fc = np.median([h["f_cert"] for h in hh]) if hh else 0.0
        f_use = fs if sf else fc
        d_cert = 0.5 * np.log10(1.0 + R_AMP ** 2 * f_use)
        g_all = np.median([r["gain"] for r in acc])
        g_hit = np.median([r["gain"] for r in got]) if got else 0.0
        print(f"      {nm:<18} turns {np.median([r['turns'] for r in acc]):4.1f}  "
              f"genuine/run {np.mean([r['genuine'] for r in acc]):4.1f}  "
              f"runs with a genuine find {len(got)}/{len(acc)}")
        print(f"      {'':<18} f_cert {fc:.4f}  f_marg {fs:.4f}  -> pulsar-term depth "
              f"{d_cert:+.4f} dex")
        print(f"      {'':<18} total depth {g_all:+.4f} dex (all runs) / {g_hit:+.4f} "
              f"(runs that found something)  ->  rank "
              f"{rank_of(g_all):2d} / {rank_of(g_hit):2d}   (against 11)")
    out["p3_M"] = np.array(Ms)
    out["p3_rank"] = np.array(rk)


def exp6_chirp(out, n_real=12):
    print()
    print("=" * 100)
    print("DELIVERABLE 6 -- RE-CUT: CHIRP MASS FROM THE LINE WIDTH, WITH COMPETITORS")
    print("=" * 100)
    print("  VERNIER-23's spin-off assumed 0-1 in-bin competitors. Re-cut here against")
    print("  the audited sky (N_eff ~ 8-9 competitors, lambda_eff-consistent budget),")
    print("  beside the zero-competitor control at the same loudness.")
    print()
    print(f"  {'f nHz':>6} {'width x T':>10} {'resolv':>7} | "
          f"{'COHERENT competitors':>21} {'competitors as Gauss':>21} "
          f"{'isolated (white)':>17} | {'coh/iso':>8}")
    rows = []
    rng = np.random.default_rng(24)
    for f in F_REPORT:
        r = chirp_precision(rng, N_PSR, f, n_real=n_real)
        w = r["coh"][1]
        rows.append([f * 1e9, w, r["coh"][0], r["gauss"][0], r["iso"][0]])
        print(f"  {f*1e9:6.1f} {w:10.2f} {'yes' if w > 1 else 'NO':>7} | "
              f"{r['coh'][0]:21.4f} {r['gauss'][0]:21.4f} {r['iso'][0]:17.4f} | "
              f"{r['coh'][0]/max(r['iso'][0],1e-30):8.2f}")
    print()
    print("  Three columns, because 'zero competitor' is ambiguous once the census is")
    print("  corrected: COHERENT = the audited bin with its bright end modelled and")
    print("  entering the budget through MEASURED template overlaps (VERNIER-21);")
    print("  GAUSS = the same power charged as an unresolved Gaussian sea; ISOLATED =")
    print("  the same source against white noise alone. coh < gauss is the lambda_eff")
    print("  suppression doing real work, not an error.")
    print("  VERNIER-23 banked 0.007-0.029 dex assuming 0-1 competitors; superseded.")
    out["chirp"] = np.array(rows)


def main():
    out = {}
    print()
    print("#" * 100)
    print("# VERNIER-24 -- THE CLEAN-BAND WHEEL ON THE AUDITED SKY.")
    print("# EXPLORATORY, NOT PRE-REGISTERED.")
    print("#" * 100)
    exp0_census(out)
    exp1_frontier(out)
    exp2_wheels(out)
    exp3_controls(out)
    exp4_loudness(out)
    exp5_part3(out)
    exp6_chirp(out)
    np.savez(os.path.join(REPORTS, "vernier24.npz"), **out)
    print()
    print(f"  banked -> {os.path.join(REPORTS, 'vernier24.npz')}")


if __name__ == "__main__":
    main()
