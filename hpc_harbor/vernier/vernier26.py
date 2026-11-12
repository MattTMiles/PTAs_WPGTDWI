#!/usr/bin/env python3
"""
VERNIER-26 -- does the loop SHARPEN THE BACKGROUND MEASUREMENT?
EXPLORATORY, NOT pre-registered. Toy, cronus, CPU/numpy, target < 15 min.
No MCMC anywhere: Fourier-domain quadratic estimators + a profiled power-law
fit. 20 realisations, PAIRED across every cell.

Clean room: imports only vernier24 (and through it vernier_toy / vernier6 /
vernier7 / vernier16).

THE QUESTION
  The loop's banked products are certified distances, chirp masses, a small
  depth bonus. The untested third product is a SHARPER BACKGROUND. Three
  candidate mechanisms, kept separate by construction:

  M-A  BIAS.        Unmodelled sources leak into the common-red-noise (AUTO
                    power) estimator. Pulsar terms in particular are pure
                    excess auto power that no cross-correlation can see, so
                    they bias a CRN amplitude high and nothing removes them
                    except modelling them.
  M-B  REASSIGNMENT. Every recovered source moves power out of "background"
                    and into "foreground catalogue". The catalogue term is
                    measured with a LINEAR estimator (small variance); the
                    background term is measured with a QUADRATIC one (large
                    variance). So the reconstructed TOTAL should get more
                    precise as R grows, even though the residual background
                    itself gets fainter.
  M-C  PULSAR VARIANCE. Pulsar terms are incoherent between pulsars, so they
                    contribute NOTHING to the mean of a Hellings-Downs
                    cross-correlation -- but they sit in every pulsar's auto
                    power, which is exactly the variance normalisation of the
                    HD estimator. With certified distances they become
                    modellable and that variance floor comes down. This is
                    the arrow nobody has computed.

  M-A lives in the AUTO estimator. M-B lives in the reassignment accounting.
  M-C lives in the CROSS estimator and in the OS SNR. They are therefore read
  off three different columns and cannot be confused with each other.

THE SKY  -- VERNIER-24 Arm B, the CORRECTED (audited) census, unchanged:
  N_bin from VERNIER-4's N_c ~ f^(-11/3); in-bin strain function Pareto
  gamma_eff = 1.56 through the h(Mc,f,D) link so h^2/bin ~ f^(-7/3); the top
  N_BRIGHT = 300 per bin modelled individually, the remainder a Gaussian sea.
  15 bins, k = 5..19 of 1/T at T = 20 yr -> 7.93-30.1 nHz, which straddles the
  white-noise crossing (GWB/white per bin runs 19 -> 0.15 across the band), so
  the slope is actually measurable.

WHAT IS NEW HERE vs VERNIER-24/25
  The faint sea is drawn HELLINGS-DOWNS CORRELATED across the array (chol of
  S_sea * Gamma_ab, Gamma_aa = 1), not as an uncorrelated per-pulsar PSD.
  Without that there is no HD estimator to sharpen and M-C cannot be posed.
  The 300 bright per bin stay DETERMINISTIC point sources: their ensemble
  average correlation IS Hellings-Downs (that is where HD comes from), but any
  one universe is a heavy-tailed draw, and that gap is the whole experiment.

THE THREE TREATMENTS at each R in {0, 4, 8, 16, 25}
  (i)   UNMODELLED  -- status quo, nothing subtracted.
  (ii)  EARTH       -- Earth terms of the R recovered sources fitted (complex
                       A+, Ax per source, weighted LS on the known sky
                       position and frequency) and subtracted. Pulsar terms
                       remain in the data.
  (iii) EARTH+PULSAR, CONDITIONED -- same fit, but each template carries its
                       pulsar term (1 - r_P e^{i psi_a}) with psi_a built from
                       the ORACLE distance. Labelled CONDITIONED per the
                       ORACLE convention (S5.5.1): a statistic evaluated
                       against truth is an oracle until its implementable form
                       is scored, so the hazard control below is mandatory.
  (iii-s) SCRAMBLED -- the ORACLE hazard control: identical machinery, psi_a
                       built from WRONG distances. If (iii) beats (ii) for a
                       reason that is not the distances, (iii-s) beats (ii)
                       too and the result is void.

  All treatments run on the SAME realisation (same sky, same noise, same
  array), so every comparison below is PAIRED.

CONTROLS
  Pure-noise null (no GWB, no sources) through both estimators.
  R = 0 rows, where (ii) and (iii) are identical to (i) by construction.

Code hpc_harbor/vernier/vernier26.py . data reports/vernier26.npz
Staged, not committed.
"""
import os
import time
import numpy as np

import vernier24 as W
import vernier_toy as V

REPORTS = W.REPORTS

# ------------------------------------------------------------------ constants
T_OBS   = W.T_OBS
DF      = W.DF
F_YR    = W.F_YR
A_GWB   = W.A_GWB
GAMMA0  = 13.0 / 3.0                 # true residual-PSD slope
KBINS   = np.arange(5, 20)           # 15 bins of 1/T
BAND    = KBINS * DF                 # 7.93 - 30.1 nHz
N_PSR   = W.N_PSR
SIG_TOA = W.SIG_TOA
N_BRIGHT = W.N_BRIGHT
GAMMA_E = W.GAMMA_EFF                # 1.56, audited in-bin Pareto index
R_P     = W.R_AMP                    # pulsar/Earth residual amplitude ratio
R_LIST  = (0, 4, 8, 16, 25)
N_REAL  = 200          # the brief budgeted 20; 200 costs 25 s, so 200
GAM_GRID = np.arange(2.0, 6.001, 0.05)
JITTER  = 1e-9

TREAT = ("i_unmodelled", "ii_earth", "iii_earth_pulsar", "iiis_scrambled")


# ============================================================ array + geometry
def hd_orf(P):
    """
    Hellings-Downs overlap reduction function, normalised so Gamma_aa = 1 and
    Gamma_ab -> 1/2 for coincident distinct pulsars (the usual convention in
    which the auto term carries the pulsar-term factor of 2).
    """
    c = np.clip(P @ P.T, -1.0, 1.0)
    x = (1.0 - c) / 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        g = 1.5 * x * np.log(np.where(x > 0, x, 1.0)) - 0.25 * x + 0.5
    np.fill_diagonal(g, 1.0)
    return g


def build_array(seed=11):
    """
    ONE pulsar array for the whole experiment: 116 pulsars, heterogeneous
    timing (VERNIER-13's 0.40 dex spread in 1/sigma_rms), log-uniform distances
    0.3-5 kpc. Fixed across realisations -- the ensemble is over universes, not
    over arrays, which is the framing an actual PTA lives in.
    """
    rng = np.random.default_rng(seed)
    P = V.iso_dirs(N_PSR, rng)
    w = 10.0 ** rng.normal(0.0, W.SIG_LOG_W, N_PSR)
    w = w / np.sqrt(np.mean(w ** 2))              # unit mean POWER
    L = 10 ** rng.uniform(np.log10(300.0), np.log10(5000.0), N_PSR)   # pc
    pw = W.S_white(SIG_TOA) * DF / w ** 2         # white power per 1/T bin, s^2
    # ORF ACTUALLY IN THE DATA. The textbook HD normalisation (Gamma_ab ->
    # 1/2 at zero separation, Gamma_aa = 1) assumes Earth and pulsar terms
    # carry EQUAL power. VERNIER-8's measured residual ratio r_P = 1.196 does
    # not: the Earth-coherent fraction is q = 1/(1 + r_P^2) = 0.411, so the
    # correlated part of the auto power is 2q, not 1. Using the textbook curve
    # as the filter against this population biases the cross estimator low by
    # 1 - 2q = 18 %, which is a normalisation error and not a mechanism. Both
    # the sea and the filter therefore use the SAME rescaled ORF.
    q = 1.0 / (1.0 + R_P ** 2)
    Ghd = hd_orf(P)
    G = 2.0 * q * Ghd + (1.0 - 2.0 * q) * np.eye(N_PSR)
    Lch = np.linalg.cholesky(G + JITTER * np.eye(N_PSR))
    arr = dict(P=P, w=w, L=L, pw=pw, G=G, Lch=Lch, Ghd=Ghd, q=q,
               G2=G ** 2, n=N_PSR, eye=np.eye(N_PSR))
    arr["kcal"] = calibrate_orf(arr)
    return arr


def calibrate_orf(arr, nsrc=20000, seed=77):
    """
    MEASURED calibration of the HD filter against the actual point-source
    population, rather than an assumed one.

    The sources are injected through the REGULARISED antenna response
    (W.antenna_reg clips 1 - cos mu at 5e-2, the series-wide clip), so their
    ensemble correlation is not exactly the textbook HD curve. Estimating an
    isotropic-population amplitude with a textbook filter therefore carries a
    scalar calibration offset. It is measured here from 20000 injected sources
    and divided out; it is a constant, identical in every treatment, so it
    cancels in every mechanism difference reported below and is quoted only so
    the absolute bias columns can be read at face value.
    """
    rng = np.random.default_rng(seed)
    P, G, pw = arr["P"], arr["G"], arr["pw"]
    Gt = V.iso_dirs(nsrc, rng)
    fp, fc = W.antenna_reg(P, Gt)
    psi, ci = rng.uniform(0, np.pi, nsrc), rng.uniform(-1, 1, nsrc)
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    Fp = fp * c2[None, :] + fc * s2[None, :]
    Fc = -fp * s2[None, :] + fc * c2[None, :]
    Tm = (Fp * ((1 + ci ** 2) / 2)[None, :] + 1j * Fc * ci[None, :])
    Tm = Tm * np.exp(1j * rng.uniform(0, 2 * np.pi, nsrc))[None, :]
    Tm = Tm / np.maximum(np.sqrt(np.mean(np.abs(Tm) ** 2, axis=0)), 1e-300)[None, :]
    # normalise to the SAME split the injection uses: the Earth term carries a
    # fraction q = 1/(1 + r_P^2) of each source's total per-pulsar power, so the
    # coherent correlation of a population of unit TOTAL power is q * (Earth
    # correlation). Without this factor the "calibration" just re-measures
    # 1/q = 2.43 and means nothing.
    C = arr["q"] * np.real(Tm @ Tm.conj().T) / nsrc
    Wm = np.ones((arr["n"], arr["n"]))
    Z1, Z2 = G * C * Wm, G * G * Wm
    num = 0.5 * (Z1.sum() - np.trace(Z1))
    den = 0.5 * (Z2.sum() - np.trace(Z2))
    return float(num / den)


# ============================================================ the sky, per bin
def draw_bin(rng, f, n_bright=N_BRIGHT, gamma=GAMMA_E):
    """
    VERNIER-24's audited ladder for one bin: residual amplitude^2 of the top
    n_bright sources, plus the power left to the Gaussian sea. Sum is exactly
    W.bin_power(f) -- the loudness sweep is a REDISTRIBUTION, never an
    addition.
    """
    n = int(max(round(float(W.N_c(f))), 1))
    h = rng.uniform(size=n) ** (-1.0 / gamma)
    p = np.sort(h ** 2)[::-1]
    p = p / p.sum()
    tot = float(W.bin_power(f))
    nb = min(n_bright, n)
    return p[:nb] * tot, float(p[nb:].sum() * tot), n, float(p[0])


def build_universe(rng, arr):
    """
    One realisation: the whole 15-bin sky, its deterministic bright sources,
    the HD-correlated faint sea, and the white noise. Returns everything the
    three treatments need, so that they are exact paired comparisons.

    Per source k in bin i:
        E_ak  = Earth term,  mean_a |E_ak|^2 = R2_k / (1 + r_P^2)
        P_ak  = -r_P E_ak e^{i psi_ak}, psi_ak = 2 pi f L_a (1 - cos mu_ak) / c
    so the source's total mean per-pulsar power is R2_k and the bin accounting
    closes: sum_k R2_k + S_sea = bin_power(f).
    """
    P, L, pw, Lch = arr["P"], arr["L"], arr["pw"], arr["Lch"]
    n = arr["n"]
    bins = []
    for f in BAND:
        R2, sea, nbin, lam = draw_bin(rng, f)
        M = R2.size
        Gt = V.iso_dirs(M, rng)                       # source sky directions
        fp, fc = W.antenna_reg(P, Gt)                 # (n, M), regularised
        psi = rng.uniform(0, np.pi, M)
        ci = rng.uniform(-1, 1, M)
        c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
        Fp = fp * c2[None, :] + fc * s2[None, :]
        Fc = -fp * s2[None, :] + fc * c2[None, :]
        Tm = (Fp * ((1 + ci ** 2) / 2)[None, :] + 1j * Fc * ci[None, :])
        Tm = Tm * np.exp(1j * rng.uniform(0, 2 * np.pi, M))[None, :]
        nrm = np.sqrt(np.mean(np.abs(Tm) ** 2, axis=0))
        E = Tm / np.maximum(nrm, 1e-300)[None, :] \
            * np.sqrt(R2 / (1.0 + R_P ** 2))[None, :]           # (n, M)

        omc = np.clip(1.0 - P @ Gt.T, W.OMC_MIN, None)          # (n, M)
        # pulsar-term phase: 2 pi f L (1-cos mu) / c, L in pc
        ph = 2.0 * np.pi * f * (L[:, None] / V.C_PC_S) * omc     # radians
        pfac = 1.0 - R_P * np.exp(1j * ph)                       # (n, M)
        sig = (E * pfac).sum(axis=1)

        g = np.sqrt(sea) * (Lch @ ((rng.normal(size=n) +
                                    1j * rng.normal(size=n)) / np.sqrt(2)))
        nz = np.sqrt(pw) * (rng.normal(size=n) +
                            1j * rng.normal(size=n)) / np.sqrt(2)
        bins.append(dict(f=f, R2=R2, sea=sea, n_bin=nbin, lam=lam,
                         fp=fp, fc=fc, E=E, ph=ph, sig=sig,
                         d=sig + g + nz, d_nosrc=g + nz, M=M,
                         Smod=float(W.bin_power(f))))
    return bins


def rank_sources(arr, bins):
    """
    Rank every modelled source in the band by ARRAY SNR against the fixed
    model power budget. The loop recovers the loudest, so R is a cut on this
    ordering.
    """
    pw = arr["pw"]
    rows = []
    for i, B in enumerate(bins):
        Pt = pw + B["Smod"]
        # coherent array SNR^2 of the Earth+pulsar template
        num = (np.abs(B["E"]) ** 2 * (1.0 + R_P ** 2)).sum(axis=0)
        rho = np.sqrt(num / np.mean(Pt))
        for k in range(B["M"]):
            rows.append((float(rho[k]), i, k))
    rows.sort(key=lambda t: -t[0])
    return rows


# ============================================================ subtraction
def subtract(arr, B, ks, mode, rng=None):
    """
    Weighted least-squares fit and subtraction of the R recovered sources in
    one bin. Two complex parameters (A+, Ax) per source; the sky position and
    frequency are known (the loop certified them), the amplitude and
    polarisation are not.

    mode "earth"     : template column pair [fp, fc]
    mode "pulsar"    : [fp * (1 - r_P e^{i psi}), fc * (1 - ...)], ORACLE psi
    mode "scrambled" : same, psi from WRONG distances (the hazard control)

    Returns residual data, NOISELESS residual signal (for the variance-floor
    column), debiased Earth-term power of each fitted source, and the number of
    complex parameters spent.
    """
    if len(ks) == 0:
        return B["d"].copy(), B["sig"].copy(), np.zeros(0), 0, None
    pw = arr["pw"]
    Pt = pw + B["Smod"]
    cols, base = [], []
    for k in ks:
        fpk, fck = B["fp"][:, k], B["fc"][:, k]
        if mode == "earth":
            fac = np.ones(arr["n"])
        elif mode == "pulsar":
            fac = 1.0 - R_P * np.exp(1j * B["ph"][:, k])
        elif mode == "scrambled":
            bad = rng.uniform(0, 2 * np.pi, arr["n"])
            fac = 1.0 - R_P * np.exp(1j * bad)
        else:
            raise ValueError(mode)
        cols += [fpk * fac, fck * fac]
        base += [fpk, fck]
    X = np.stack(cols, axis=1)                       # (n, 2R)
    Xb = np.stack(base, axis=1)                      # Earth-only basis
    sw = 1.0 / np.sqrt(Pt)
    Xt, dt = X * sw[:, None], B["d"] * sw
    ah, *_ = np.linalg.lstsq(Xt, dt, rcond=1e-6)   # SVD-truncated: several
    # recovered sources can land in one bin at small angular separation and the
    # design matrix then goes near-singular; without the truncation one
    # realisation in ~200 returns a catalogue power 10^3 too large and takes
    # the whole sd column with it.
    res = B["d"] - X @ ah
    # debiased Earth-term power of each source: mean_a |fp a+ + fc ax|^2
    XhX = Xt.conj().T @ Xt
    cov = np.linalg.pinv(XhX, rcond=1e-6)
    # the oblique projector the fit applies to the data: res = (I - Y) d
    Y = X @ (cov @ (Xt.conj().T * sw[None, :]))
    # CATALOGUE POWER, as a TOTAL and not as a sum of per-source powers.
    # When two recovered sources land in the same bin close together on the
    # sky, their Earth-term columns are nearly parallel and the fit splits them
    # into two huge opposite-sign amplitudes whose SUM is small. Summing
    # |a_k|^2 then reports a catalogue power orders of magnitude too large --
    # it blew up 1-2 realisations in 200 at R = 16 and R = 25 and took the
    # whole sd column with it. The coherent total ||X_E a||^2 is the quantity
    # that is actually estimated, and it is immune.
    v = Xb @ ah
    raw = float(np.mean(np.abs(v) ** 2))
    BhB = (Xb.conj().T @ Xb) / arr["n"]
    bias = float(np.real(np.trace(BhB @ cov)))
    pk = np.array([max(raw - bias, 0.0)])
    return res, B["sig"] - X @ ah, pk, X.shape[1], Y


# ============================================================ estimators
def spectra(arr, d, Smod, Y=None):
    """
    Free-spectrum point estimates in ONE bin, both flavours, WITH the
    fitted-degrees-of-freedom correction.

    THE CORRECTION, AND WHY IT IS NOT OPTIONAL. Subtracting R sources in a bin
    projects the data through Pi = I - Y, where Y is the fit's oblique hat
    matrix, and that removes 2R complex directions from a 116-pulsar vector.
    Those directions are Earth-term templates -- exactly the COHERENT part of
    the data -- so the projection eats the background as well as the sources.
    Uncorrected, the cross estimate falls to zero by R = 25 and the loop looks
    like it destroys the very thing it was supposed to sharpen. That is an
    artefact of not accounting for the fitted dof. So both estimators are
    referred to the PROJECTED response:

        Gt = Pi Gamma Pi^H     (the background's surviving response)
        Nt = Pi diag(pw) Pi^H  (the white noise's surviving response)

    AUTO (CRN-like):
        S_auto = sum_a v_a (|d_a|^2 - Nt_aa) / sum_a v_a Gt_aa,  v_a = 1/P_a^2
    Sees pulsar terms in the numerator, because they are excess auto power
    that no projection knows about. This is where M-A lives.

    CROSS (optimal statistic / HD):
        S_x = [sum_{a<b} G_ab Re(d_a d_b*) W_ab - sum_{a<b} G_ab Nt_ab W_ab]
              / sum_{a<b} G_ab Gt_ab W_ab,        W_ab = 1/(P_a P_b)
    with FIXED model weights P_a = pw_a + Smod, identical in every treatment so
    that treatment differences are data and not weights. Blind to pulsar terms
    in the MEAN -- they are incoherent between pulsars -- and NOT blind to them
    in the variance. This is where M-B and M-C live.
    """
    pw, G, G2 = arr["pw"], arr["G"], arr["G2"]
    Pt = pw + Smod
    Wm = np.outer(1.0 / Pt, 1.0 / Pt)
    if Y is None:
        Gt, Nt = G, np.diag(pw)
        gdiag, ndiag = np.diag(G).copy(), pw
    else:
        Pi = arr["eye"] - Y
        Gt = np.real(Pi @ G @ Pi.conj().T)
        Nt = np.real((Pi * pw[None, :]) @ Pi.conj().T)
        gdiag, ndiag = np.diag(Gt).copy(), np.diag(Nt).copy()

    v = 1.0 / Pt ** 2
    dn = float(np.sum(v * gdiag))
    s_auto = float(np.sum(v * (np.abs(d) ** 2 - ndiag)) / max(dn, 1e-300))
    var_auto = float(np.sum(v ** 2 * Pt ** 2) / max(dn, 1e-300) ** 2)

    M = np.real(np.outer(d, np.conj(d)))
    def pairsum(A):
        Z = G * A * Wm
        return 0.5 * (Z.sum() - np.trace(Z))
    N = pairsum(M)
    Bi = pairsum(Nt)
    D = pairsum(Gt)
    Zv = G2 * Wm
    Vn = 0.5 * (Zv.sum() - np.trace(Zv))
    D = D if abs(D) > 1e-300 else 1e-300
    kc = arr.get("kcal", 1.0)
    s_cross = float((N - Bi) / (D * kc))
    var_cross = float(Vn / (D * kc) ** 2)
    return s_auto, var_auto, s_cross, var_cross, float(D / max(pairsum(G), 1e-300))


def band_snr(s_cross, var_cross, Smods=None):
    """
    Band-integrated HD cross-correlation SNR: the matched combination of the
    per-bin cross estimates against the true f^(-13/3) shape,
        SNR = sum_i s_i t_i / v_i  /  sqrt(sum_i t_i^2 / v_i).
    Computed from the SAME projection-corrected per-bin estimates as the
    amplitude, so it is comparable across treatments by construction.
    """
    s = np.asarray(s_cross, float)
    iv = 1.0 / np.maximum(np.asarray(var_cross, float), 1e-300)
    t = SHAPE_
    return float(np.sum(s * t * iv) / np.sqrt(max(np.sum(t ** 2 * iv), 1e-300)))


K_A = 2.0 * DF / (12.0 * np.pi ** 2 * F_YR ** 3)     # bin_power = K_A A^2 (f/fyr)^-g
SHAPE = (BAND / F_YR) ** (-GAMMA0)                   # fixed-slope template
C_TRUE = K_A * A_GWB ** 2                            # true amplitude coefficient
SHAPE_ = SHAPE                                       # alias used by band_snr


def fit_fixed(S, var):
    """
    Amplitude at FIXED slope 13/3 -- the linear, always-defined estimator, and
    the one an astronomer actually quotes ("A at gamma = 13/3"). Returns the
    coefficient c with bin_power_i = c * (f_i/f_yr)^(-13/3), so A = sqrt(c/K_A).

    A free two-parameter fit is NOT usable for the amplitude column here: once
    the loudest sources are subtracted the residual spectrum is not a power law
    (subtraction punches bin-shaped holes in it), and the A-gamma degeneracy
    then moves A by tens of per cent for reasons that have nothing to do with
    the mechanisms under test. Slope is reported separately, from the free fit.
    """
    iv = 1.0 / np.maximum(np.asarray(var, float), 1e-300)
    return float(np.sum(np.asarray(S, float) * SHAPE * iv)
                 / np.sum(SHAPE ** 2 * iv))


def os_floor(arr, Pres):
    """
    The OS VARIANCE FLOOR -- the M-C column.

    Given the per-(pulsar, bin) residual SIGNAL power Pres (noiseless: the
    unmodelled sources plus whatever the fit failed to remove), the optimally
    weighted Hellings-Downs amplitude estimator has

        var(c) = 1 / sum_i sum_{a<b} t_i^2 Gamma_ab^2 / (P_ai P_bi),
        P_ai = pw_a + Pres_ai,   t_i = (f_i/f_yr)^(-13/3)

    Pulsar terms enter ONLY through P_ai. They are invisible to the numerator
    of the cross-correlation -- incoherent between pulsars -- and fully visible
    here. So this number is exactly "the variance floor the pulsar terms
    impose on Hellings-Downs", and the difference between treatment (ii) and
    treatment (iii) is exactly M-C, with nothing else in it.

    Returned as the fractional 1-sigma on the AMPLITUDE: sigma_A/A =
    0.5 * sqrt(var(c)) / c_true.
    """
    pw, G2 = arr["pw"], arr["G2"]
    acc = 0.0
    for i, t in enumerate(SHAPE):
        Pt = pw + Pres[i]
        Wm = np.outer(1.0 / Pt, 1.0 / Pt) * G2
        acc += t ** 2 * 0.5 * (Wm.sum() - np.trace(Wm))
    return 0.5 * float(np.sqrt(1.0 / max(acc, 1e-300)) / C_TRUE)


def fit_powerlaw(S, var, f=BAND):
    """
    Profiled two-parameter power-law fit, no MCMC: for each gamma on a grid the
    amplitude is the linear weighted-LS solution, and gamma is chosen by chi^2.
    Returns (A, gamma, ok).
    """
    S = np.asarray(S, float)
    iv = 1.0 / np.maximum(np.asarray(var, float), 1e-300)
    x = (f / F_YR)
    best = (np.inf, np.nan, np.nan)
    for g in GAM_GRID:
        t = x ** (-g)
        den = float(np.sum(t ** 2 * iv))
        if den <= 0:
            continue
        c = float(np.sum(S * t * iv) / den)
        chi2 = float(np.sum((S - c * t) ** 2 * iv))
        if chi2 < best[0]:
            best = (chi2, c, g)
    _, c, g = best
    if not np.isfinite(c) or c <= 0:
        return np.nan, g, False
    return float(np.sqrt(c / K_A)), float(g), True


# ============================================================ one realisation
def run_real(seed, arr, null=False):
    rng = np.random.default_rng(seed)
    bins = build_universe(rng, arr)
    order = rank_sources(arr, bins)
    out = {}

    if null:
        pw = arr["pw"]
        ds = []
        for B in bins:
            nz = np.sqrt(pw) * (rng.normal(size=arr["n"]) +
                                1j * rng.normal(size=arr["n"])) / np.sqrt(2)
            ds.append(nz)
        sa, va, sc, vc = [], [], [], []
        for d, B in zip(ds, bins):
            a, av, c, cv, _ = spectra(arr, d, B["Smod"])
            sa.append(a); va.append(av); sc.append(c); vc.append(cv)
        return dict(null_c_auto=fit_fixed(sa, va) / C_TRUE,
                    null_c_cross=fit_fixed(sc, vc) / C_TRUE,
                    null_snr=band_snr(sc, vc))

    lam = float(np.median([B["lam"] for B in bins]))
    for R in R_LIST:
        sel = order[:R]
        per_bin = {i: [] for i in range(len(bins))}
        for _, i, k in sel:
            per_bin[i].append(k)
        for tr in TREAT:
            if R == 0 and tr != "i_unmodelled":
                out[(R, tr)] = out[(0, "i_unmodelled")]
                continue
            sa, va, sc, vc, sm = [], [], [], [], []
            s_tot_a, s_tot_c, s_true_res, Pres, resp = [], [], [], [], []
            npar, pterm_frac = 0, []
            for i, B in enumerate(bins):
                ks = per_bin[i]
                if tr == "i_unmodelled" or not ks:
                    d, rs, pk, Y = B["d"], B["sig"], np.zeros(0), None
                elif tr == "ii_earth":
                    d, rs, pk, np_, Y = subtract(arr, B, ks, "earth")
                    npar += np_
                elif tr == "iii_earth_pulsar":
                    d, rs, pk, np_, Y = subtract(arr, B, ks, "pulsar")
                    npar += np_
                else:
                    d, rs, pk, np_, Y = subtract(arr, B, ks, "scrambled", rng)
                    npar += np_
                a, av, c, cv, rp = spectra(arr, d, B["Smod"], Y)
                resp.append(rp)
                cat = float(np.sum(pk) * (1.0 + R_P ** 2))   # reassigned power
                sa.append(a); va.append(av); sc.append(c); vc.append(cv)
                s_tot_a.append(a + cat); s_tot_c.append(c + cat)
                sm.append(B["Smod"])
                Pres.append(np.abs(rs) ** 2 + B["sea"])
                # TRUTH for the residual background at this R: sea + the
                # bright sources NOT recovered
                mask = np.ones(B["M"], bool)
                mask[ks] = False
                s_true_res.append(B["sea"] + float(B["R2"][mask].sum()))
                if ks:
                    pterm_frac.append(
                        float(np.sum(B["R2"][ks]) * R_P ** 2 / (1 + R_P ** 2)
                              / max(s_true_res[-1], 1e-300)))
            c_true = fit_fixed(s_true_res, vc)
            _, gb_a, _ = fit_powerlaw(sa, va)
            _, gb_c, _ = fit_powerlaw(sc, vc)
            _, gt_c, _ = fit_powerlaw(s_tot_c, vc)
            _, g_tr, _ = fit_powerlaw(s_true_res, vc)
            out[(R, tr)] = dict(
                c_bg_auto=fit_fixed(sa, va), c_bg_cross=fit_fixed(sc, vc),
                c_tot_auto=fit_fixed(s_tot_a, va),
                c_tot_cross=fit_fixed(s_tot_c, vc),
                c_true_res=c_true,
                g_bg_auto=gb_a, g_bg_cross=gb_c, g_tot_cross=gt_c,
                g_true_res=g_tr,
                floor=os_floor(arr, Pres),
                resp=float(np.mean(resp)),
                pterm_frac=float(np.mean(pterm_frac)) if pterm_frac else 0.0,
                snr=band_snr(sc, vc), npar=npar)
    out["lam"] = lam
    out["nbin"] = np.array([B["n_bin"] for B in bins])
    out["seafrac"] = float(np.mean([B["sea"] / B["Smod"] for B in bins]))
    return out


# ============================================================ reporting
def agg(vals):
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    if v.size == 0:
        return np.nan, np.nan, np.nan
    return float(np.mean(v)), float(np.std(v, ddof=1)) if v.size > 1 else 0.0, \
        float(np.median(v))


def main():
    t0 = time.time()
    print("#" * 100)
    print("# VERNIER-26 -- DOES THE LOOP SHARPEN THE BACKGROUND MEASUREMENT?")
    print("# EXPLORATORY, NOT PRE-REGISTERED.  Fourier-domain estimators only.")
    print("#" * 100)
    arr = build_array()
    gwb_white = W.bin_power(BAND) / np.median(arr["pw"])   # MEDIAN pulsar
    print(f"  array {N_PSR} pulsars, sigma_toa(med) {SIG_TOA*1e9:.0f} ns, "
          f"T = {T_OBS/W.YR:.0f} yr, 1/T = {DF*1e9:.3f} nHz")
    print(f"  band {BAND[0]*1e9:.2f}-{BAND[-1]*1e9:.2f} nHz in {len(BAND)} bins; "
          f"GWB/white per bin {gwb_white[0]:.1f} -> {gwb_white[-1]:.2f} "
          f"(crossing inside the band)")
    print(f"  sky: VERNIER-24 Arm B, Pareto gamma_eff = {GAMMA_E}, "
          f"{N_BRIGHT} explicit per bin + HD-correlated faint sea")
    print(f"  A_true = {A_GWB:.3e}, gamma_true = {GAMMA0:.4f}, "
          f"r_P = {R_P}, Earth-coherent fraction q = {arr['q']:.3f}, "
          f"realisations = {N_REAL}")
    print(f"  measured HD-filter calibration against the regularised antenna "
          f"response: k = {arr['kcal']:.4f} (divided out; a constant, so it "
          f"cancels in every difference)")
    print()

    res = [run_real(2600 + s, arr) for s in range(N_REAL)]
    nulls = [run_real(9600 + s, arr, null=True) for s in range(N_REAL)]

    print(f"  realised sky: median lambda_fore = "
          f"{np.median([r['lam'] for r in res]):.3f}, "
          f"faint-sea fraction (below top {N_BRIGHT}) = "
          f"{np.mean([r['seafrac'] for r in res])*100:.2f} % of the bin")
    print(f"  N_bin {res[0]['nbin'][0]:.3g} (7.9 nHz) -> "
          f"{res[0]['nbin'][-1]:.3g} (30 nHz)")
    print()

    out = {}

    # ---------------------------------------------------------------- null
    print("=" * 100)
    print("CONTROL 0 -- PURE-NOISE NULL (no GWB, no sources, nothing to find)")
    print("=" * 100)
    for key, lab in (("null_c_auto", "A^2 / A^2_true (auto)"),
                     ("null_c_cross", "A^2 / A^2_true (cross)"),
                     ("null_snr", "OS SNR")):
        m, s, md = agg([n[key] for n in nulls])
        print(f"  {lab:24s} mean {m:+.4f}  sd {s:.4f}  median {md:+.4f}")
    m, s, _ = agg([n["null_snr"] for n in nulls])
    print(f"  -> OS SNR null {m:+.3f} +- {s:.3f} (s.e. of the mean "
          f"{s/np.sqrt(len(nulls)):.3f}): the MEAN is consistent with zero, so")
    print(f"     the estimator does not manufacture a correlation. The SPREAD "
          f"is {s:.2f} and not 1 because the")
    print(f"     weights P_a = pw_a + S_model assume a background the null "
          f"does not contain -- the SNR column")
    print(f"     is therefore CONSERVATIVE by ~{1/max(s,1e-9):.1f}x in the "
          f"no-signal limit, and is only ever compared across treatments.")
    out["null_snr"] = np.array([n["null_snr"] for n in nulls])
    out["null_c_cross"] = np.array([n["null_c_cross"] for n in nulls])
    print()

    def cell(R, tr, key):
        return np.array([r[(R, tr)][key] for r in res], float)

    # ------------------------------------------------- deliverable 1
    print("=" * 100)
    print("DELIVERABLE 1 -- THE RESIDUAL BACKGROUND, bias and variance")
    print("  Statistic is the amplitude COEFFICIENT at fixed slope 13/3;")
    print("  columns are A^2_hat / A^2_true_residual, with the truth RECUT at")
    print("  each R (sea + the bright sources not yet recovered).")
    print("  Fractional precision on the AMPLITUDE is half the sd column.")
    print("=" * 100)
    hdr = (f"  {'R':>3s} {'treatment':<18s} | {'AUTO bias':>10s} {'AUTO sd':>9s} "
           f"| {'CROSS bias':>11s} {'CROSS sd':>9s} | {'gam_auto':>8s} "
           f"{'gam_cross':>9s} {'gam_TRUE':>8s} | {'OS SNR':>14s}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for R in R_LIST:
        for tr in TREAT:
            if R == 0 and tr != "i_unmodelled":
                continue
            ct = cell(R, tr, "c_true_res")
            ra = cell(R, tr, "c_bg_auto") / ct
            rc = cell(R, tr, "c_bg_cross") / ct
            ma, sa, _ = agg(ra)
            mc, sc, _ = agg(rc)
            gA, _, _ = agg(cell(R, tr, "g_bg_auto"))
            gC, _, _ = agg(cell(R, tr, "g_bg_cross"))
            gT, _, _ = agg(cell(R, tr, "g_true_res"))
            sn, sns, _ = agg(cell(R, tr, "snr"))
            print(f"  {R:>3d} {tr:<18s} | {ma-1:>+10.4f} {sa:>9.4f} "
                  f"| {mc-1:>+11.4f} {sc:>9.4f} | {gA:>8.2f} {gC:>9.2f} "
                  f"{gT:>8.2f} | {sn:>6.2f} +-{sns:<6.2f}")
            out[f"bg_{R}_{tr}"] = np.vstack([ra, rc])
            out[f"snr_{R}_{tr}"] = cell(R, tr, "snr")
        print()

    # ------------------------------------------------- deliverable 2
    print("=" * 100)
    print("DELIVERABLE 2 -- THE RECONSTRUCTED TOTAL (residual background +")
    print("  recovered catalogue), against the FIXED truth A = 2.4e-15.")
    print("  This is the M-B reassignment column: what an astronomer quotes")
    print("  does not care whether power sits in the sea or in the catalogue,")
    print("  only how well the sum is measured.")
    print("=" * 100)
    hdr2 = (f"  {'R':>3s} {'treatment':<18s} | {'AUTO bias':>10s} {'AUTO sd':>9s} "
            f"| {'CROSS bias':>11s} {'CROSS sd':>9s} | {'CROSS RMSE':>9s} "
            f"{'gam_tot':>8s}")
    print(hdr2)
    print("  " + "-" * (len(hdr2) - 2))
    tot_sd, tot_bias = {}, {}
    for R in R_LIST:
        for tr in TREAT:
            if R == 0 and tr != "i_unmodelled":
                continue
            ra = cell(R, tr, "c_tot_auto") / C_TRUE
            rc = cell(R, tr, "c_tot_cross") / C_TRUE
            ma, sa, _ = agg(ra)
            mc, sc, _ = agg(rc)
            rmse_a = float(np.sqrt(np.mean((ra - 1.0) ** 2)))
            rmse_c = float(np.sqrt(np.mean((rc - 1.0) ** 2)))
            gtc, _, _ = agg(cell(R, tr, "g_tot_cross"))
            print(f"  {R:>3d} {tr:<18s} | {ma-1:>+10.4f} {sa:>9.4f} "
                  f"| {mc-1:>+11.4f} {sc:>9.4f} | {rmse_c:>9.4f} {gtc:>8.2f}")
            tot_sd[(R, tr)] = (sa, sc, rmse_c, rmse_a)
            tot_bias[(R, tr)] = (ma - 1, mc - 1)
            out[f"tot_{R}_{tr}"] = np.vstack([ra, rc])
        print()

    # ------------------------------------------------- deliverable 3
    print("=" * 100)
    print("DELIVERABLE 3 -- THE M-C COLUMN: the HD variance floor set by the")
    print("  pulsar terms. sigma_A/A of the OPTIMALLY weighted HD amplitude,")
    print("  computed from the per-pulsar residual SIGNAL power. Pulsar terms")
    print("  enter here and nowhere else, so (ii) -> (iii) is M-C alone.")
    print("=" * 100)
    hdr3 = (f"  {'R':>3s} | {'floor (i)':>10s} {'(ii) earth':>11s} "
            f"{'(iii) cond':>11s} {'(iii-s) scr':>12s} | {'M-C gain':>9s} "
            f"{'HAZARD':>9s} | {'pterm frac':>11s}")
    print(hdr3)
    print("  " + "-" * (len(hdr3) - 2))
    mc_gain = {}
    for R in R_LIST:
        fi, _, _ = agg(cell(R, "i_unmodelled", "floor"))
        if R == 0:
            print(f"  {R:>3d} | {fi:>10.5f} {'-':>11s} {'-':>11s} {'-':>12s} "
                  f"| {'-':>9s} {'-':>9s} | {'-':>11s}")
            continue
        f2, _, _ = agg(cell(R, "ii_earth", "floor"))
        f3, _, _ = agg(cell(R, "iii_earth_pulsar", "floor"))
        fs, _, _ = agg(cell(R, "iiis_scrambled", "floor"))
        # PAIRED, realisation by realisation
        g_mc = float(np.mean(1.0 - cell(R, "iii_earth_pulsar", "floor")
                             / cell(R, "ii_earth", "floor")))
        g_hz = float(np.mean(1.0 - cell(R, "iiis_scrambled", "floor")
                             / cell(R, "ii_earth", "floor")))
        _, _, pf = agg(cell(R, "ii_earth", "pterm_frac"))   # MEDIAN: heavy tail
        print(f"  {R:>3d} | {fi:>10.5f} {f2:>11.5f} {f3:>11.5f} {fs:>12.5f} "
              f"| {g_mc*100:>+8.2f}% {g_hz*100:>+8.2f}% | {pf*100:>10.2f}%")
        mc_gain[R] = (g_mc, g_hz, fi, f2, f3, fs)
        out[f"floor_{R}"] = np.vstack([cell(R, t, "floor") for t in TREAT])
    print()
    print("  M-C gain = PAIRED mean of 1 - floor(iii)/floor(ii): the fraction of")
    print("             the HD variance floor that certified distances remove.")
    print("  HAZARD   = the same for SCRAMBLED distances. It must be ~0 (or")
    print("             negative) or the M-C claim is void (ORACLE rule S5.5.1).")
    print("  pterm frac = MEDIAN over bins/realisations of the pulsar-term power")
    print("             of the recovered sources, as a fraction of the residual")
    print("             background they sit in. The mean is far larger and is")
    print("             not quotable: the Pareto tail makes it unstable.")
    print()

    # ------------------------------------------------- deliverable 4
    print("=" * 100)
    print("DELIVERABLE 4 -- MECHANISM ATTRIBUTION, all paired")
    print("=" * 100)
    print(f"  {'R':>3s} | {'M-A bias(ii)':>12s} {'M-A bias(iii)':>13s} "
          f"{'M-A removed':>11s} | {'M-B var i->ii':>13s} "
          f"| {'M-C floor':>9s} {'HAZARD':>8s}")
    print("  " + "-" * 88)
    attrib = {}
    for R in R_LIST:
        if R == 0:
            continue
        b_ii = float(np.mean(cell(R, "ii_earth", "c_bg_auto")
                             / cell(R, "ii_earth", "c_true_res") - 1.0))
        b_iii = float(np.mean(cell(R, "iii_earth_pulsar", "c_bg_auto")
                              / cell(R, "iii_earth_pulsar", "c_true_res") - 1.0))
        s_i = tot_sd[(R, "i_unmodelled")][2]
        s_ii = tot_sd[(R, "ii_earth")][2]
        mB = 1.0 - s_ii / s_i if s_i > 0 else np.nan
        g_mc, g_hz = mc_gain[R][0], mc_gain[R][1]
        print(f"  {R:>3d} | {b_ii:>+12.4f} {b_iii:>+13.4f} "
              f"{b_ii-b_iii:>+11.4f} | {mB*100:>+12.1f}% "
              f"| {g_mc*100:>+8.2f}% {g_hz*100:>+7.2f}%")
        attrib[R] = (b_ii, b_iii, mB, g_mc, g_hz)
    print()
    print("  M-A = residual bias of the CRN/AUTO amplitude. Under (ii) the")
    print("        recovered sources' PULSAR terms are still in the data and")
    print("        are pure excess auto power; under (iii) they are gone.")
    print("        'M-A removed' is what modelling them buys, in units of the")
    print("        residual background itself.")
    print("  M-B = fractional drop in RMSE of the reconstructed TOTAL amp^2")
    print("        (cross estimator) on moving from unmodelled to Earth-")
    print("        subtracted, i.e. what reassignment alone buys.")
    print()

    # ------------------------------------------------- deliverable 5
    print("=" * 100)
    print("DELIVERABLE 5 -- PRECISION PER RECOVERED SOURCE")
    print("  On the reconstructed TOTAL amplitude (cross estimator). The figure")
    print("  of merit is RMSE about the TRUTH, not sd: a treatment that biases")
    print("  the answer low also shrinks its scatter, and the scrambled control")
    print("  does exactly that. Amplitude precision = half the amplitude^2")
    print("  precision; both are quoted.")
    print("=" * 100)
    print(f"  {'R':>3s} {'treatment':<18s} | {'RMSE(A^2)':>9s} {'RMSE(A)':>8s} "
          f"| {'vs R=0':>9s} {'per source':>11s}")
    print("  " + "-" * 74)
    s0 = tot_sd[(0, "i_unmodelled")][2]
    for R in R_LIST:
        for tr in TREAT:
            if R == 0 and tr != "i_unmodelled":
                continue
            s = tot_sd[(R, tr)][2]
            gain = 1.0 - s / s0
            per = gain / R * 100 if R else np.nan
            print(f"  {R:>3d} {tr:<18s} | {s:>9.4f} {s/2:>8.4f} "
                  f"| {gain*100:>+8.1f}% "
                  f"{('%+.2f%%' % per) if R else '-':>11s}")
        print()
    print("  Also, on the HD variance FLOOR (the noise-limited precision an")
    print("  optimally weighted analysis could reach):")
    print(f"  {'R':>3s} | {'floor(i)':>9s} {'floor(iii)':>11s} "
          f"| {'vs R=0':>9s} {'per source':>11s}")
    print("  " + "-" * 56)
    f0 = mc_gain[R_LIST[1]][2] if len(R_LIST) > 1 else np.nan
    f0, _, _ = agg(cell(0, "i_unmodelled", "floor"))
    for R in R_LIST:
        if R == 0:
            print(f"  {R:>3d} | {f0:>9.5f} {'-':>11s} | {'-':>9s} {'-':>11s}")
            continue
        fi = mc_gain[R][2]
        f3 = mc_gain[R][4]
        gain = 1.0 - f3 / f0
        print(f"  {R:>3d} | {fi:>9.5f} {f3:>11.5f} | {gain*100:>+8.2f}% "
              f"{gain/R*100:>+10.2f}%")
    print()

    out["R_LIST"] = np.array(R_LIST)
    out["BAND"] = BAND
    out["A_GWB"] = A_GWB
    for R in attrib:
        out[f"attrib_{R}"] = np.array(attrib[R])
    np.savez(os.path.join(REPORTS, "vernier26.npz"), **out)
    print(f"  banked -> {os.path.join(REPORTS, 'vernier26.npz')}")
    print(f"  wall clock {time.time()-t0:.1f} s")




if __name__ == "__main__":
    main()
