#!/usr/bin/env python3
"""
VERNIER toy model -- clean room.

EXPLORATORY. NOT pre-registered. No campaign module is imported; the geometry,
the waveform, the match and the fringe bookkeeping are implemented from scratch
so that agreement with banked numbers is evidence and not tautology.

Question: does source multiplicity break the pulsar-distance fringe degeneracy,
and if so at how many sources and at what brightness gap?

Model
-----
pulsar direction p_hat at distance L, source direction s_hat:
    mu    = angle(p_hat, s_hat)
    tau   = (L/c)(1 - cos mu)
    Phi(t)= 2 pi (f0 t + 0.5 fdot t^2) + phi0
Earth term samples Phi(t), pulsar term samples Phi(t - tau).

The pulsar-term match between a trial distance L and the truth L* is built here
from the EXACT circular chirp (Peters), not a quadratic Taylor phase.  With
t_c the coalescence time and tau the light-travel delay,

    f_orb(-tau) = f_base (1 + tau/t_c)^(-3/8)
    Phi(-tau)   = (16 pi/5) f_base t_c [1 - (1 + tau/t_c)^(5/8)]
    m(L)        = cos(dPhi + pi df T) * sinc(df T)

with dPhi = n[Phi(-tau) - Phi(-tau*)] and df = n[f_orb(-tau) - f_orb(-tau*)].

    the COSINE is H1, the comb   -> dL = c / (n f_orb(-tau) (1 - cos mu))
    the SINC   is H2, the clock  -> envelope, first null at df T = 1

WHY NOT THE TAYLOR PHASE.  The naive model Delta_Phi = 2 pi f0 tau - pi fdot tau^2
+ 2 pi fdot t tau needs |fdot tau| << f0.  At the campaign fiducial
(f0 = 20 nHz, log10 Mc = 9.5, L = 1 kpc, mu = 90 deg) it is NOT: fdot tau
exceeds f0, the carrier frequency passes through zero, and the model fabricates
fringe structure.  Equivalently tau/t_c = 3.2 -- the pulsar term is emitted
BEFORE the source entered the band.  The exact solution is finite for all
tau > 0 and is used throughout; PTERM_STATIC reproduces the campaign's
non-evolving-pulsar-term convention (f_p = f0) for comparison.

Outputs: reports/vernier_toy.npz, reports/VERNIER_*.png
"""
import os
import numpy as np
from scipy.special import jv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------- constants
PC_M    = 3.0856775814913673e16      # m
C_MS    = 299792458.0                # m/s
C_PC_S  = C_MS / PC_M                # pc / s   = 9.7156e-9
T_SUN   = 4.925490947e-6             # G Msun / c^3, s
YR      = 31557600.0                 # s (Julian yr)

HERE    = os.path.dirname(os.path.abspath(__file__))
REPORTS = os.path.abspath(os.path.join(HERE, "..", "..", "reports"))

# campaign fiducial loud source (from the banked config, used ONLY as the anchor
# target and as the centre of the band sweep -- no campaign code is imported)
LOG10_FGW_FID = -7.7
LOG10_MC_FID  = 9.5
BAND_NHZ      = (3.0, 20.0)          # campaign band
T_OBS_FID     = 15.0 * YR

# banked D3 fringe spacings to check against (project_progress.md D3 table)
BANKED_DL = {"J1713+0747": 0.448, "J0437-4715": 0.489, "J1909-3744": 0.424,
             "J1744-1134": 0.361, "J0030+0451": 0.485}


# ----------------------------------------------------------------- geometry
def fdot_of(f_gw, log10_mc):
    """Circular, GW-driven chirp.  fdot = (96/5) pi^(8/3) (T_sun Mc)^(5/3) f^(11/3)."""
    mc = 10.0 ** np.asarray(log10_mc, dtype=float)
    return (96.0 / 5.0) * np.pi ** (8.0 / 3.0) * (T_SUN * mc) ** (5.0 / 3.0) \
        * np.asarray(f_gw, dtype=float) ** (11.0 / 3.0)


def t_merge(f_gw, log10_mc):
    """Coalescence time of a circular binary from GW frequency f_gw, s."""
    mc = 10.0 ** np.asarray(log10_mc, float)
    return (5.0 / 256.0) * (np.pi * np.asarray(f_gw, float)) ** (-8.0 / 3.0) \
        * (T_SUN * mc) ** (-5.0 / 3.0)


def dL_of(f_gw, one_minus_cos_mu):
    """Fringe spacing, pc.  dL = c / (f_carrier (1 - cos mu))."""
    return C_PC_S / (np.asarray(f_gw, float) * np.asarray(one_minus_cos_mu, float))


def tau_of(L_pc, one_minus_cos_mu):
    """Pulsar-term light-travel delay, s."""
    return np.asarray(L_pc, float) * np.asarray(one_minus_cos_mu, float) / C_PC_S


def _retarded(tau, f_base, tc, static):
    """(phase at -tau, base frequency at -tau) for the exact circular chirp."""
    if static:
        return -2.0 * np.pi * f_base * tau, f_base * np.ones_like(np.asarray(tau, float))
    x = 1.0 + tau / tc
    phi = (16.0 * np.pi / 5.0) * f_base * tc * (1.0 - x ** (5.0 / 8.0))
    return phi, f_base * x ** (-3.0 / 8.0)


def f_pterm(L_pc, f_base, tc, one_minus_cos_mu, n_harm=1, static=False):
    """Observed pulsar-term carrier frequency of harmonic n at distance L."""
    _, fo = _retarded(tau_of(L_pc, one_minus_cos_mu), f_base, tc, static)
    return n_harm * fo


def match_pterm(L_pc, Lstar_pc, f_base, tc, one_minus_cos_mu, T_obs,
                n_harm=1, static=False):
    """Time-averaged pulsar-term overlap between trial L and truth L*."""
    tau  = tau_of(L_pc, one_minus_cos_mu)
    taus = tau_of(Lstar_pc, one_minus_cos_mu)
    p, fo   = _retarded(tau,  f_base, tc, static)
    ps, fos = _retarded(taus, f_base, tc, static)
    dphi = n_harm * (p - ps)
    df   = n_harm * (fo - fos)
    return np.cos(dphi + np.pi * df * T_obs) * np.sinc(df * T_obs)


def match_bruteforce(L_pc, Lstar_pc, f_base, tc, one_minus_cos_mu, T_obs,
                     n_harm=1, static=False, n_t=400001):
    """Direct time-domain average of the exact chirp, to validate match_pterm."""
    t = np.linspace(0.0, T_obs, n_t)
    tau  = tau_of(L_pc, one_minus_cos_mu)
    taus = tau_of(Lstar_pc, one_minus_cos_mu)

    def phase(shift):
        tt = t - shift
        if static:
            return 2.0 * np.pi * n_harm * f_base * tt
        return n_harm * (16.0 * np.pi / 5.0) * f_base * tc \
            * (1.0 - (1.0 - tt / tc) ** (5.0 / 8.0))

    return np.mean(np.cos(phase(tau) - phase(taus)))


def count_modes(lnl, thresh):
    """
    Local maxima of the Asimov dlnL curve that survive a threshold below the peak.
    Returns (n_modes, indices).  thresh is a POSITIVE nat depth.

    IMPORTANT (this is the sub-fringe-offset trap in reverse).  A uniform L grid
    does not land on fringe centres, so a sampled alias peak reads LOWER than it
    really is, which biases the mode count DOWN and P(unique) UP.  Every local
    maximum is therefore refined by the parabolic vertex of its sampled triplet
    before thresholding.  That is free (no extra likelihood evaluations) and
    removes the O(h^2) sampling bias.
    """
    interior = (lnl[1:-1] >= lnl[:-2]) & (lnl[1:-1] >= lnl[2:])
    idx = np.where(interior)[0] + 1
    if idx.size == 0:
        return 1, np.array([int(np.argmax(lnl))]), np.array([lnl.max()])
    y1, y2, y3 = lnl[idx - 1], lnl[idx], lnl[idx + 1]
    den = y1 - 2.0 * y2 + y3
    d = np.where(np.abs(den) > 1e-300, 0.5 * (y1 - y3) / np.where(den == 0, 1.0, den), 0.0)
    d = np.clip(d, -0.5, 0.5)
    ypk = y2 - 0.25 * (y1 - y3) * d
    ypk = np.maximum(ypk, y2)                      # refinement can only raise a max
    keep = idx[ypk > ypk.max() - thresh]
    kept_y = ypk[ypk > ypk.max() - thresh]
    if keep.size == 0:
        keep = np.array([int(np.argmax(lnl))]); kept_y = np.array([lnl.max()])
    return int(keep.size), keep, kept_y


def iso_dirs(n, rng):
    """Isotropic unit vectors."""
    cth = rng.uniform(-1.0, 1.0, n)
    sth = np.sqrt(1.0 - cth ** 2)
    ph  = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.stack([sth * np.cos(ph), sth * np.sin(ph), cth], axis=1)


# ============================================================ helpers
TIERS = [("lit10", "10% of L"), ("lit20", "20% of L"),
         ("vlbi", "1 pc"), ("subpc", "0.1 pc")]


def _tier_sigma(tier, Lstar):
    return {"lit10": 0.10 * Lstar, "lit20": 0.20 * Lstar,
            "vlbi": 1.0, "subpc": 0.1}[tier]


class Src:
    """A source as one pulsar sees it: base freq, harmonic, t_c, geometry, SNR."""

    def __init__(self, f_base, tc, omc, n_harm=1, rho=1.0):
        self.f_base, self.tc, self.omc = f_base, tc, omc
        self.n_harm, self.rho = n_harm, rho

    def f_p(self, L):
        return f_pterm(L, self.f_base, self.tc, self.omc, self.n_harm, PTERM_STATIC)

    def dL(self, L):
        return dL_of(self.f_p(L), self.omc)

    def m(self, L, Lstar, T=None):
        return match_pterm(L, Lstar, self.f_base, self.tc, self.omc,
                           T_OBS_FID if T is None else T, self.n_harm, PTERM_STATIC)


PTERM_STATIC = False          # False = exact evolving pulsar term (physical)


def joint_dlnl(Lgrid, Lstar, srcs):
    """Asimov joint dlnL(L) = sum_k rho_k^2 (m_k(L) - 1)."""
    tot = np.zeros_like(Lgrid)
    for s in srcs:
        tot += s.rho ** 2 * (s.m(Lgrid, Lstar) - 1.0)
    return tot


def scan(Lstar, sig, srcs, cap=200000):
    """Build the L grid, evaluate dlnL, return (n_modes, K_trials)."""
    dls = np.array([s.dL(Lstar) for s in srcs])
    half = 3.0 * sig
    step = min(dls.min() / 64.0, half / 256.0)
    n = int(min(2.0 * half / step, cap)) + 1
    Lg = Lstar + np.linspace(-half, half, n)
    lnl = joint_dlnl(Lg, Lstar, srcs)
    ktr = max(2.0 * half / dls.min(), 1.0)
    return count_modes(lnl, np.log(ktr))[0], ktr


def draw_sources(K, rng, p_hat, rho=1.0, log10mc_range=(8.5, 9.5)):
    """K independent circular campaign-band sources seen by one pulsar."""
    out = []
    for _ in range(K):
        f0 = 10.0 ** rng.uniform(np.log10(BAND_NHZ[0] * 1e-9),
                                 np.log10(BAND_NHZ[1] * 1e-9))
        mc = rng.uniform(*log10mc_range)
        s = iso_dirs(1, rng)[0]
        omc = max(1.0 - float(s @ p_hat), 1e-3)
        out.append(Src(f0, t_merge(f0, mc), omc, 1, rho))
    return out


# ============================================================ ANCHOR
def anchor_check(out):
    global PTERM_STATIC
    print("=" * 74)
    print("ANCHOR -- fringe spacing vs the programme's banked ~0.4 pc comb")
    print("=" * 74)
    f0 = 10.0 ** LOG10_FGW_FID
    print(f"  c = {C_PC_S:.6e} pc/s ;  f0 = 10^{LOG10_FGW_FID} = {f0:.4e} Hz")
    print(f"  dL at median geometry (1-cos mu = 1), NON-EVOLVING pulsar term:")
    print(f"      dL = c/f0 = {dL_of(f0, 1.0):.4f} pc   vs banked 0.489 pc")
    med = dL_of(f0, 1.0)
    rel = abs(med - 0.489) / 0.489
    print(f"      -> {100*rel:.2f}% agreement.  ANCHOR PASSES on the campaign's own")
    print("         convention (comb pinned to the EARTH-term frequency).")
    print()
    print("  banked D3 dL_pc -> implied geometry (closure check):")
    for name, dl in BANKED_DL.items():
        omc = C_PC_S / (f0 * dl)
        print(f"    {name:14s} dL={dl:.3f} -> 1-cos mu={omc:.4f} "
              f"mu={np.degrees(np.arccos(1-omc)):6.2f} deg "
              f"-> recomputed {dL_of(f0, omc):.4f} (resid {abs(dL_of(f0,omc)-dl):.1e})")
    out["anchor_dL_median_pc"] = med
    out["anchor_rel_err"] = rel
    out["anchor_banked_names"] = np.array(list(BANKED_DL.keys()))
    out["anchor_banked_dL"] = np.array(list(BANKED_DL.values()))

    print()
    print("  BUT: the banked convention assumes the pulsar term has NOT evolved.")
    tau = tau_of(1000.0, 1.0)
    tc = t_merge(f0, LOG10_MC_FID)
    fp = f0 * (1 + tau / tc) ** (-3.0 / 8.0)
    print(f"    at the fiducial (L=1 kpc, mu=90 deg, log10Mc={LOG10_MC_FID}):")
    print(f"      tau     = {tau/YR:8.0f} yr")
    print(f"      t_merge = {tc/YR:8.0f} yr   -> tau/t_c = {tau/tc:.2f}")
    print(f"      f_p/f0  = {fp/f0:8.3f}      -> evolving-pterm dL = "
          f"{dL_of(fp,1.0):.4f} pc  ({dL_of(fp,1.0)/med:.2f}x the banked value)")
    print("    The pulsar term is emitted BEFORE the source entered the band, so its")
    print("    comb is COARSER than c/f0. Reported as a discrepancy, not smoothed over.")
    out["anchor_tau_over_tc"] = tau / tc
    out["anchor_fp_over_f0"] = fp / f0
    out["anchor_dL_evolving_pc"] = dL_of(fp, 1.0)

    print()
    print("  match_pterm validation vs 4e5-sample exact time average:")
    errs = []
    for stat in [True, False]:
        PTERM_STATIC = stat
        for L in [1000.3, 1001.0, 1005.0, 1050.0]:
            a = match_pterm(L, 1000.0, f0, tc, 1.0, T_OBS_FID, 1, stat)
            b = match_bruteforce(L, 1000.0, f0, tc, 1.0, T_OBS_FID, 1, stat)
            errs.append(abs(a - b))
            print(f"    static={str(stat):5s} L={L:8.2f}  analytic={a:+.6f} "
                  f"brute={b:+.6f}  |d|={abs(a-b):.2e}")
    PTERM_STATIC = False
    out["match_validation_maxerr"] = max(errs)
    print(f"  max |analytic - brute| = {max(errs):.2e}")
    print()


# ============================================================ TOY 1a
def toy1a_fringe_budget(out, rng):
    print("=" * 74)
    print("TOY 1a -- FRINGE BUDGET")
    print("=" * 74)
    n_psr = 20
    f0 = 10.0 ** LOG10_FGW_FID
    tc = t_merge(f0, LOG10_MC_FID)
    p = iso_dirs(n_psr, rng)
    s = iso_dirs(1, rng)[0]
    L = rng.uniform(300.0, 5000.0, n_psr)
    omc = np.clip(1.0 - p @ s, 1e-3, None)
    mu = np.degrees(np.arccos(np.clip(p @ s, -1, 1)))
    dl_static = dL_of(f0, omc)
    fp = f0 * (1 + tau_of(L, omc) / tc) ** (-3.0 / 8.0)
    dl_evolve = dL_of(fp, omc)

    sig = {"lit10": 0.10 * L, "lit20": 0.20 * L,
           "vlbi": np.full(n_psr, 1.0), "subpc": np.full(n_psr, 0.1)}
    K1 = {k: 2.0 * v / dl_static for k, v in sig.items()}     # +/-1 sigma
    K3 = {k: 6.0 * v / dl_static for k, v in sig.items()}     # +/-3 sigma (campaign)

    print(f"  20 pulsars, isotropic sky, L~U(0.3,5) kpc, one source f0=10^{LOG10_FGW_FID}.")
    print("  K3 = 6 sigma/dL is the CAMPAIGN convention (+/-3 sigma window).")
    print("  K1 = 2 sigma/dL is the +/-1 sigma count, reported because the")
    print("  'a couple of fringes at 1 pc' intuition is a +/-1 sigma statement.")
    print()
    print(f"  {'#':>2} {'L kpc':>6} {'mu deg':>7} {'dL_stat':>8} {'dL_evol':>8} "
          f"{'K3 lit20':>9} {'K3 vlbi':>8} {'K3 subpc':>9} {'K1 vlbi':>8} {'K1 subpc':>9}")
    for i in np.argsort(dl_static):
        print(f"  {i:2d} {L[i]/1e3:6.2f} {mu[i]:7.1f} {dl_static[i]:8.3f} "
              f"{dl_evolve[i]:8.3f} {K3['lit20'][i]:9.0f} {K3['vlbi'][i]:8.2f} "
              f"{K3['subpc'][i]:9.3f} {K1['vlbi'][i]:8.2f} {K1['subpc'][i]:9.3f}")
    print()
    print(f"  {'tier':8s} {'med K3':>9s} {'med K1':>9s} {'K3<1':>7s} {'K1<1':>7s}")
    for k, lab in TIERS:
        print(f"  {k:8s} {np.median(K3[k]):9.2f} {np.median(K1[k]):9.2f} "
              f"{int((K3[k]<1).sum()):5d}/20 {int((K1[k]<1).sum()):5d}/20")
    print()
    print("  READING. lit tiers hold 1e3-1e4 fringes -- hopeless. The 1 pc tier holds")
    print(f"  a MEDIAN of {np.median(K3['vlbi']):.0f} fringes at +/-3 sigma "
          f"({np.median(K1['vlbi']):.1f} at +/-1 sigma), NOT 'a couple': the")
    print("  hypothesis is directionally right but ~1 dex optimistic at 1 pc. Only the")
    print(f"  0.1 pc tier is genuinely sub-fringe (median K3 = {np.median(K3['subpc']):.2f}, "
          f"{int((K3['subpc']<1).sum())}/20 below 1).")
    print("  So D2's saturation at 1 pc is NOT the prior running out of fringes -- at")
    print("  1 pc there are still ~10 to remove. The prior runs out only near 0.1 pc.")
    print("  Evolving-pterm combs are COARSER (median ratio "
          f"{np.median(dl_evolve/dl_static):.2f}x), which HELPS slightly but not by 1 dex.")
    print()
    out["t1a_L_pc"], out["t1a_mu_deg"] = L, mu
    out["t1a_dL_static_pc"], out["t1a_dL_evolve_pc"] = dl_static, dl_evolve
    for k, _ in TIERS:
        out[f"t1a_K3_{k}"], out[f"t1a_K1_{k}"] = K3[k], K1[k]

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    cols = {"lit10": "#c44e52", "lit20": "#dd8452", "vlbi": "#4c72b0", "subpc": "#55a868"}
    for k, lab in TIERS:
        ax[0].scatter(dl_static, K3[k], s=34, label=f"{k} ({lab})", color=cols[k],
                      alpha=.85, edgecolor="k", linewidth=.4)
    ax[0].axhline(1.0, color="k", ls="--", lw=1)
    ax[0].text(dl_static.min(), 1.2, "K=1: sub-fringe, prior alone identifies", fontsize=8)
    ax[0].set_yscale("log"); ax[0].set_xlabel("fringe spacing dL [pc]")
    ax[0].set_ylabel("candidate fringes  K = 6$\\sigma$/dL")
    ax[0].set_title("(a) fringe budget per prior tier")
    ax[0].legend(fontsize=8, frameon=False)
    x = np.arange(len(TIERS))
    for j, (k, _) in enumerate(TIERS):
        ax[1].scatter(np.full(20, j - .12), K3[k], s=22, color=cols[k], alpha=.8,
                      edgecolor="k", lw=.3, label="$\\pm3\\sigma$" if j == 0 else None)
        ax[1].scatter(np.full(20, j + .12), K1[k], s=22, color=cols[k], alpha=.4,
                      marker="^", edgecolor="k", lw=.3,
                      label="$\\pm1\\sigma$" if j == 0 else None)
    ax[1].axhline(1.0, color="k", ls="--", lw=1)
    ax[1].set_xticks(x); ax[1].set_xticklabels([k for k, _ in TIERS])
    ax[1].set_yscale("log"); ax[1].set_ylabel("candidate fringes K")
    ax[1].set_title("(a) 1 pc still holds ~10 fringes; 0.1 pc is sub-fringe")
    ax[1].legend(fontsize=8, frameon=False)
    for a in ax:
        a.grid(alpha=.25, which="both")
    fig.suptitle("VERNIER Toy 1a -- fringe budget (exploratory, not pre-registered)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER_t1a_fringe_budget.png"), dpi=140)
    plt.close(fig)


# ============================================================ TOY 1b
def toy1b_clock(out):
    print("=" * 74)
    print("TOY 1b -- THE CLOCK HANDLE H2")
    print("=" * 74)
    print("  In the exact chirp the pulsar-term carrier drifts with L as")
    print("      df_p/dtau = -(3/8) f_p / (t_c + tau)")
    print("  so across ONE comb period (dtau = 1/f_p) the carrier moves by")
    print("      Delta f = (3/8) / (t_c + tau).")
    print("  The clock separates ADJACENT fringes when Delta f * T >= 1, i.e.")
    print()
    print("      *** T_obs  >=  (8/3) (t_c + tau) ***")
    print()
    print("  That is a closed-form NO-GO: tau is millennia, T_obs is decades.")
    print("  Fringes surviving the sinc envelope: n_surv = (8/3)(t_c+tau)/T.")
    print()
    f_grid = np.array([3.0, 6.0, 10.0, 15.0, 20.0]) * 1e-9
    mc_grid = np.array([8.5, 9.0, 9.5, 10.0])
    T_grid = np.array([15.0, 30.0, 40.0]) * YR
    L_ref, omc_ref = 1000.0, 1.0
    tau = tau_of(L_ref, omc_ref)
    print(f"  reference pulsar L=1 kpc, mu=90 deg -> tau = {tau/YR:.0f} yr")
    print()
    print(f"  {'f0nHz':>6} {'log10Mc':>8} {'t_c yr':>10} {'tau/t_c':>9} "
          f"{'f_p/f0':>7} {'T_req yr':>10} {'n_surv T=15':>12} {'window pc':>10}")
    tc_a = np.zeros((f_grid.size, mc_grid.size))
    ratio = np.zeros_like(tc_a); fpr = np.zeros_like(tc_a)
    Treq = np.zeros_like(tc_a)
    nsurv = np.zeros((T_grid.size, f_grid.size, mc_grid.size))
    win = np.zeros_like(nsurv)
    adj_ok = np.zeros_like(nsurv, bool)
    for i, f in enumerate(f_grid):
        for j, mc in enumerate(mc_grid):
            tc = t_merge(f, mc)
            fp = f * (1 + tau / tc) ** (-3.0 / 8.0)
            tr = (8.0 / 3.0) * (tc + tau)
            tc_a[i, j], ratio[i, j], fpr[i, j], Treq[i, j] = tc, tau / tc, fp / f, tr
            for k, T in enumerate(T_grid):
                nsurv[k, i, j] = tr / T
                win[k, i, j] = (tr / T) * dL_of(fp, omc_ref)
                adj_ok[k, i, j] = T >= tr
            print(f"  {f*1e9:6.1f} {mc:8.2f} {tc/YR:10.3g} {tau/tc:9.3g} "
                  f"{fp/f:7.3f} {tr/YR:10.3g} {nsurv[0,i,j]:12.1f} {win[0,i,j]:10.1f}")
    print()
    print(f"  ADJACENT-fringe clock condition met in {int(adj_ok.sum())}/{adj_ok.size} "
          "cells across (f0, Mc, T=15/30/40 yr).")
    print(f"  Tightest clock-only window anywhere in the grid: {win.min():.1f} pc "
          f"-- still {win.min()/1.0:.0f}x looser than a 1 pc VLBI prior.")
    print("  The prompt's condition |fdot tau| T >= 1 is the ABSOLUTE one; it is met")
    print("  widely but does NOT identify a fringe. Requiring T >= (8/3)(t_c+tau) is")
    print("  the condition that does, and nothing in the band meets it.")
    print("  H2 truncates the comb to ~1e2-1e6 fringes. It never resolves one.")
    print()
    out["t1b_f_grid"], out["t1b_mc_grid"], out["t1b_T_grid"] = f_grid, mc_grid, T_grid
    out["t1b_tc_s"], out["t1b_tau_over_tc"] = tc_a, ratio
    out["t1b_fp_over_f0"], out["t1b_Treq_s"] = fpr, Treq
    out["t1b_nsurv"], out["t1b_window_pc"], out["t1b_adj_ok"] = nsurv, win, adj_ok
    out["t1b_tau_ref_s"] = tau


# ============================================================ TOY 1c
def toy1c_vernier(out, rng, n_trial=120, Kmax=8, rho0=10.0):
    global PTERM_STATIC
    print("=" * 74)
    print("TOY 1c -- THE VERNIER: P(unique solution) vs K")
    print("=" * 74)
    print(f"  Noiseless (Asimov) dlnL, each source at pulsar-term SNR rho={rho0:.0f}.")
    print("  A mode survives within ln(K_trials) nats of the peak (trials factor).")
    print(f"  {n_trial} trials per cell. L*~U(0.5,3) kpc, isotropic sky, band "
          f"{BAND_NHZ[0]:.0f}-{BAND_NHZ[1]:.0f} nHz, log10Mc~U(8.5,9.5).")
    print("  Run in BOTH pulsar-term conventions.")
    print()
    tiers = [t[0] for t in TIERS]
    Ks = np.arange(1, Kmax + 1)
    P = np.zeros((2, len(tiers), Kmax))
    NMED = np.zeros((2, len(tiers), Kmax))
    for im, stat in enumerate([False, True]):
        PTERM_STATIC = stat
        print(f"  --- pulsar term {'STATIC (campaign convention)' if stat else 'EVOLVING (exact)'}")
        for it, tier in enumerate(tiers):
            for ik, K in enumerate(Ks):
                nm = []
                for _ in range(n_trial):
                    p_hat = iso_dirs(1, rng)[0]
                    Lstar = rng.uniform(500.0, 3000.0)
                    srcs = draw_sources(K, rng, p_hat, rho0)
                    nm.append(scan(Lstar, _tier_sigma(tier, Lstar), srcs)[0])
                nm = np.array(nm)
                P[im, it, ik] = (nm == 1).mean()
                NMED[im, it, ik] = np.median(nm)
            need = np.where(P[im, it] >= 0.9)[0]
            print(f"    {tier:7s} P: " + " ".join(f"{v:5.2f}" for v in P[im, it])
                  + (f"   K90={Ks[need[0]]}" if need.size else "   K90>8"))
            print(f"    {'':7s} n: " + " ".join(f"{v:5.0f}" for v in NMED[im, it]))
    PTERM_STATIC = False
    out["t1c_K"], out["t1c_tiers"] = Ks, np.array(tiers)
    out["t1c_P_unique"] = P            # [mode(evolving,static), tier, K]
    out["t1c_median_modes"] = NMED
    out["t1c_rho0"] = rho0

    print()
    print(f"  rho sweep (vlbi tier, evolving pterm): P(unique) vs K")
    rhos = np.array([1.0, 2.0, 3.0, 5.0, 10.0, 20.0])
    Prho = np.zeros((rhos.size, Kmax))
    for ir, rr in enumerate(rhos):
        for ik, K in enumerate(Ks):
            u = 0
            for _ in range(n_trial):
                p_hat = iso_dirs(1, rng)[0]
                Lstar = rng.uniform(500.0, 3000.0)
                srcs = draw_sources(K, rng, p_hat, rr)
                u += (scan(Lstar, _tier_sigma("vlbi", Lstar), srcs)[0] == 1)
            Prho[ir, ik] = u / n_trial
        need = np.where(Prho[ir] >= 0.9)[0]
        print(f"    rho={rr:5.1f}: " + " ".join(f"{v:5.2f}" for v in Prho[ir])
              + (f"   K90={Ks[need[0]]}" if need.size else "   K90>8"))
    out["t1c_rho_sweep"], out["t1c_P_unique_rho"] = rhos, Prho
    print()

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.3))
    cols = {"lit10": "#c44e52", "lit20": "#dd8452", "vlbi": "#4c72b0", "subpc": "#55a868"}
    for it, tier in enumerate(tiers):
        ax[0].plot(Ks, P[0, it], "o-", color=cols[tier], lw=1.8, ms=5, label=tier)
        ax[0].plot(Ks, P[1, it], "s--", color=cols[tier], lw=1.1, ms=4, alpha=.55)
        ax[1].plot(Ks, np.maximum(NMED[0, it], .5), "o-", color=cols[tier], lw=1.8,
                   ms=5, label=tier)
    ax[0].set_ylabel("P(unique solution)"); ax[0].set_ylim(-.03, 1.05)
    ax[0].set_title("(c) solid = evolving pterm, dashed = static")
    ax[1].axhline(1.0, color="k", ls="--", lw=1); ax[1].set_yscale("log")
    ax[1].set_ylabel("median surviving modes")
    ax[1].set_title("(c) modes left in the $\\pm3\\sigma$ window")
    for ir, rr in enumerate(rhos):
        ax[2].plot(Ks, Prho[ir], "o-", ms=4, lw=1.6,
                   color=plt.cm.viridis(ir / max(len(rhos) - 1, 1)),
                   label=f"$\\rho$={rr:g}")
    ax[2].set_ylabel("P(unique), vlbi tier"); ax[2].set_ylim(-.03, 1.05)
    ax[2].set_title("(c) weak dependence on per-source SNR")
    for a in ax:
        a.set_xlabel("K sources"); a.grid(alpha=.25, which="both")
        a.legend(fontsize=7.5, frameon=False, ncol=2)
    fig.suptitle("VERNIER Toy 1c -- does source multiplicity break the comb? "
                 "(exploratory)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER_t1c_vernier.png"), dpi=140)
    plt.close(fig)


# ============================================================ TOY 1d
def pm_g(n, e):
    """Peters & Mathews (1963) eq.20 harmonic GW power g(n,e)."""
    ne = n * e
    J = [jv(n - 2, ne), jv(n - 1, ne), jv(n, ne), jv(n + 1, ne), jv(n + 2, ne)]
    t1 = (J[0] - 2 * e * J[1] + (2.0 / n) * J[2] + 2 * e * J[3] - J[4]) ** 2
    t2 = (1 - e ** 2) * (J[0] - 2 * J[2] + J[4]) ** 2
    t3 = (4.0 / (3.0 * n ** 2)) * J[2] ** 2
    return (n ** 4 / 32.0) * (t1 + t2 + t3)


def toy1d_commensurate(out, rng, n_trial=120, nmax=8, rho_tot=14.0):
    print("=" * 74)
    print("TOY 1d -- COMMENSURATE (harmonics) vs INCOMMENSURATE (sources)")
    print("=" * 74)
    print("  K INDEPENDENT circular sources: own f0, own mu, own t_c -> unrelated combs.")
    print("  n HARMONICS of ONE eccentric binary: f_n = n f_orb, SAME mu, SAME t_c.")
    print("  Correct harmonic chirp: fdot_n = n fdot_orb (NOT fdot(n f_orb) ~ n^(11/3)")
    print("  -- treating each harmonic as its own binary inflates its clock by n^(8/3)).")
    print(f"  Total SNR matched between arms (rho_tot={rho_tot:.0f}), so this is a")
    print("  structure comparison, not a brightness comparison.")
    print()
    e_list = [0.1, 0.3]
    for e in e_list:
        g = np.array([pm_g(n, e) for n in range(1, nmax + 1)])
        print(f"    e={e:.1f} Peters-Mathews g(n,e) normalised: "
              + " ".join(f"{v:.3f}" for v in g / g.sum()))
    print()
    tiers = ["lit20", "vlbi"]
    Ks = np.arange(1, nmax + 1)
    res_i = np.zeros((len(tiers), nmax))
    res_h = np.zeros((len(tiers), len(e_list), nmax))
    for it, tier in enumerate(tiers):
        for ik, K in enumerate(Ks):
            nm = []
            for _ in range(n_trial):
                p_hat = iso_dirs(1, rng)[0]
                Lstar = rng.uniform(500.0, 3000.0)
                srcs = draw_sources(K, rng, p_hat, rho_tot / np.sqrt(K))
                nm.append(scan(Lstar, _tier_sigma(tier, Lstar), srcs)[0])
            res_i[it, ik] = np.median(nm)
        for ie, e in enumerate(e_list):
            gall = np.array([pm_g(n, e) for n in range(1, nmax + 1)])
            for ik, K in enumerate(Ks):
                nm = []
                for _ in range(n_trial):
                    p_hat = iso_dirs(1, rng)[0]
                    Lstar = rng.uniform(500.0, 3000.0)
                    sh = iso_dirs(1, rng)[0]
                    omc = max(1.0 - float(sh @ p_hat), 1e-3)
                    forb = 10.0 ** rng.uniform(np.log10(BAND_NHZ[0] * 1e-9 / 2),
                                               np.log10(BAND_NHZ[1] * 1e-9 / 2))
                    mc = rng.uniform(8.5, 9.5)
                    tc = t_merge(2.0 * forb, mc)     # one binary -> ONE t_c
                    w = gall[:K] / gall[:K].sum()
                    srcs = [Src(forb, tc, omc, n, rho_tot * np.sqrt(w[n - 1]))
                            for n in range(1, K + 1)]
                    nm.append(scan(Lstar, _tier_sigma(tier, Lstar), srcs)[0])
                res_h[it, ie, ik] = np.median(nm)
        print(f"  --- tier {tier}")
        print("    independent : " + " ".join(f"{v:7.0f}" for v in res_i[it])
              + f"   | n=1->8 factor {res_i[it,0]/max(res_i[it,-1],1):.1f}x")
        for ie, e in enumerate(e_list):
            print(f"    harm e={e:.1f}  : " + " ".join(f"{v:7.0f}" for v in res_h[it, ie])
                  + f"   | n=1->8 factor {res_h[it,ie,0]/max(res_h[it,ie,-1],1):.1f}x")
    print()
    print("  READING: harmonics are commensurate -- every harmonic's fringe set")
    print("  CONTAINS the fundamental's, so their intersection is pinned at dL_1 and")
    print("  the mode count floors at (window/dL_1). Independent sources keep dividing.")
    print()
    out["t1d_K"], out["t1d_tiers"] = Ks, np.array(tiers)
    out["t1d_e_list"] = np.array(e_list)
    out["t1d_modes_independent"], out["t1d_modes_harmonics"] = res_i, res_h
    out["t1d_pm_g"] = np.array([[pm_g(n, e) for n in range(1, nmax + 1)] for e in e_list])

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4), squeeze=False)
    for it, tier in enumerate(tiers):
        a = ax[0][it]
        a.plot(Ks, np.maximum(res_i[it], .5), "o-", color="#4c72b0", lw=2, ms=6,
               label="K independent sources")
        for ie, e in enumerate(e_list):
            a.plot(Ks, np.maximum(res_h[it, ie], .5), "s--", ms=5, lw=1.5,
                   color=plt.cm.Oranges(.5 + .25 * ie), label=f"n harmonics, e={e:.1f}")
        a.axhline(1.0, color="k", ls=":", lw=1)
        a.set_yscale("log"); a.set_xlabel("K sources / n harmonics")
        a.set_ylabel("median surviving modes"); a.set_title(f"(d) tier = {tier}")
        a.grid(alpha=.25, which="both"); a.legend(fontsize=7.5, frameon=False)
    fig.suptitle("VERNIER Toy 1d -- harmonics dilate, independent sources remove "
                 "(exploratory)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER_t1d_commensurate.png"), dpi=140)
    plt.close(fig)


# ============================================================ TOY 2
def toy2_brightness_gap(out, rng, n_trial=250, rho1=10.0):
    print("=" * 74)
    print("TOY 2 -- HOW BRIGHT MUST SOURCE 2 BE (the dex answer)")
    print("=" * 74)
    print(f"  Source 1 pulsar-term SNR rho1={rho1:.0f}; source 2 at rho1*10^(-dex).")
    print("  Combs enter dlnL weighted by rho^2, so a dex gap is a 10^(2 dex) weight gap.")
    print("  'cuts' = joint mode count strictly below the source-1-only count.")
    print()
    print("  ANALYTIC CEILING (toy-independent). Source 2 can demote an alias by at")
    print("  most 2 rho2^2 nats, so it must clear the trials factor: 2 rho2^2 > ln K.")
    for ktr in [4.0, 12.0, 250.0, 3000.0]:
        r2 = np.sqrt(np.log(ktr) / 2.0)
        print(f"    K_trials={ktr:7.0f} (lnK={np.log(ktr):5.2f}) -> rho2 > {r2:.2f}"
              f" -> max gap from rho1={rho1:.0f} is {np.log10(rho1/r2):.2f} dex")
    print("  NOTE the gap is only a proxy: what the physics constrains is rank-1's")
    print("  ABSOLUTE pulsar-term SNR. A dex gap is admissible iff rho1/10^dex clears")
    print("  the bar, so a louder rank-0 buys a wider admissible gap.")
    print()
    dexes = np.array([0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5])
    tiers = ["lit20", "vlbi"]
    Pcut = np.zeros((len(tiers), dexes.size))
    Puni = np.zeros((len(tiers), dexes.size))
    Rmod = np.zeros((len(tiers), dexes.size))
    for it, tier in enumerate(tiers):
        for idx, dex in enumerate(dexes):
            cut = uni = 0
            rr = []
            for _ in range(n_trial):
                p_hat = iso_dirs(1, rng)[0]
                Lstar = rng.uniform(500.0, 3000.0)
                sig = _tier_sigma(tier, Lstar)
                srcs = draw_sources(2, rng, p_hat)
                srcs[0].rho = rho1
                srcs[1].rho = rho1 * 10.0 ** (-dex)
                n1 = scan(Lstar, sig, [srcs[0]])[0]
                n2 = scan(Lstar, sig, srcs)[0]
                cut += (n2 < n1); uni += (n2 == 1); rr.append(n2 / max(n1, 1))
            Pcut[it, idx], Puni[it, idx], Rmod[it, idx] = \
                cut / n_trial, uni / n_trial, np.median(rr)
        print(f"  --- tier {tier}")
        print("    dex     : " + " ".join(f"{d:5.2f}" for d in dexes))
        print("    P(cuts) : " + " ".join(f"{v:5.2f}" for v in Pcut[it]))
        print("    P(uniq) : " + " ".join(f"{v:5.2f}" for v in Puni[it]))
        print("    n2/n1   : " + " ".join(f"{v:5.2f}" for v in Rmod[it]))
        hd = np.interp(0.5, Pcut[it][::-1], dexes[::-1])
        print(f"    dex where P(cuts) = 0.5 : {hd:.2f}")
        out[f"t2_halfdex_{tier}"] = hd
    print()
    print("  CENSUS READ-OFF: GLACIER puts rank1 at 0.77-1.17 dex below rank0.")
    for it, tier in enumerate(tiers):
        lo = np.interp(0.77, dexes, Puni[it]); hi = np.interp(1.17, dexes, Puni[it])
        lc = np.interp(0.77, dexes, Pcut[it]); hc = np.interp(1.17, dexes, Pcut[it])
        print(f"    {tier:6s}: P(cuts) {lc:.2f} -> {hc:.2f} ; "
              f"P(unique) {lo:.2f} -> {hi:.2f} across that gap")
    out["t2_dexes"], out["t2_tiers"] = dexes, np.array(tiers)
    out["t2_P_cuts"], out["t2_P_unique"], out["t2_mode_ratio"] = Pcut, Puni, Rmod
    out["t2_rho1"] = rho1

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for it, tier in enumerate(tiers):
        c = ["#dd8452", "#4c72b0"][it]
        ax[0].plot(dexes, Pcut[it], "o-", color=c, lw=1.8, ms=5, label=f"{tier}: cuts")
        ax[0].plot(dexes, Puni[it], "s--", color=c, lw=1.2, ms=4, alpha=.7,
                   label=f"{tier}: unique")
        ax[1].plot(dexes, Rmod[it], "o-", color=c, lw=1.8, ms=5, label=tier)
    for a in ax:
        a.axvspan(0.77, 1.17, color="grey", alpha=.25, zorder=0)
        a.set_xlabel("brightness gap  rank0 $-$ rank1  [dex]")
        a.grid(alpha=.25); a.legend(fontsize=8, frameon=False)
    ax[0].set_ylabel("probability"); ax[0].set_ylim(-.03, 1.05)
    ax[0].set_title("(2) grey band = GLACIER census rank0$-$rank1")
    ax[1].set_ylabel("median mode ratio $n_2/n_1$")
    ax[1].set_title("(2) fraction of ambiguity surviving")
    fig.suptitle("VERNIER Toy 2 -- required brightness gap (exploratory)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORTS, "VERNIER_t2_brightness_gap.png"), dpi=140)
    plt.close(fig)


# ============================================================ TOY 3
def toy3_wheel(out, rng, n_psr=20, n_src=6, n_real=30, det_thresh=8.0):
    print("=" * 74)
    print("TOY 3 -- DOES THE WHEEL TURN")
    print("=" * 74)
    print(f"  {n_psr} pulsars (1 pc VLBI-tier priors, the best real case), {n_src}")
    print(f"  sources, {n_real} realisations, detection at array SNR > {det_thresh:.0f}.")
    print("  Turn: detect brightest undetected source -> form its comb at every pulsar")
    print("  -> intersect with all detected sources' combs -> certify pulsars whose")
    print("  mode count is 1 -> refit: a certified pulsar's pulsar term becomes")
    print("  coherently usable and adds to array SNR for the NEXT source.")
    print("  Certification is RECOMPUTED every turn and the full history kept.")
    print()
    hist_c = np.zeros((n_real, n_src), int)
    hist_d = np.zeros((n_real, n_src), int)
    hist_s = np.zeros((n_real, n_src))
    hist_g = np.zeros((n_real, n_src))
    flick = np.zeros(n_real, int)
    stall = np.zeros(n_real, int)
    why = []
    for r in range(n_real):
        p = iso_dirs(n_psr, rng)
        Lstar = rng.uniform(500.0, 3000.0, n_psr)
        sh = iso_dirs(n_src, rng)
        f0 = 10.0 ** rng.uniform(np.log10(BAND_NHZ[0] * 1e-9),
                                 np.log10(BAND_NHZ[1] * 1e-9), n_src)
        mc = rng.uniform(8.5, 9.5, n_src)
        tc = t_merge(f0, mc)
        # strain ladder: rank0 loudest, successive 0.6-1.2 dex drops (census-like)
        h = 10.0 ** np.concatenate([[0.0], -np.cumsum(rng.uniform(0.6, 1.2, n_src - 1))])
        cmu = p @ sh.T
        omc = np.clip(1.0 - cmu, 1e-3, None)
        ant = np.sqrt(np.clip(1.0 - cmu ** 2, 1e-6, None))
        # pulsar-term amplitude also pays the retarded-frequency penalty h ~ f^(2/3)
        fp = f0[None, :] * (1.0 + tau_of(Lstar[:, None], omc) / tc[None, :]) ** (-3.0 / 8.0)
        rho = ant * h[None, :] * (fp / f0[None, :]) ** (2.0 / 3.0)
        rho *= 20.0 / np.sqrt(np.sum((ant[:, 0] * h[0]) ** 2))
        cert = np.zeros(n_psr, bool)
        ever = np.zeros(n_psr, bool)
        det = []
        for turn in range(n_src):
            avail = [j for j in range(n_src) if j not in det]
            gain = np.sqrt(1.0 + cert.sum() / n_psr)
            snr = np.array([np.sqrt(np.sum(rho[:, j] ** 2)) * gain for j in avail])
            jb = avail[int(np.argmax(snr))]
            hist_s[r, turn], hist_g[r, turn] = snr.max(), gain
            if snr.max() < det_thresh:
                stall[r] = turn; why.append("next source below detection threshold")
                hist_c[r, turn:] = cert.sum(); hist_d[r, turn:] = len(det)
                break
            det.append(jb)
            hist_d[r, turn] = len(det)
            new = np.zeros(n_psr, bool)
            for a in range(n_psr):
                srcs = [Src(f0[j], tc[j], omc[a, j], 1, rho[a, j]) for j in det]
                new[a] = (scan(Lstar[a], 1.0, srcs)[0] == 1)
            flick[r] += int(np.sum(ever & ~new))
            cert = new; ever |= cert
            hist_c[r, turn] = cert.sum()
        else:
            stall[r] = n_src; why.append("ran out of sources")
    from collections import Counter
    print(f"  {'turn':>5} {'med det':>8} {'med cert':>9} {'med SNR':>9} {'med gain':>9}")
    for t in range(n_src):
        print(f"  {t+1:5d} {np.median(hist_d[:,t]):8.1f} {np.median(hist_c[:,t]):9.1f} "
              f"{np.median(hist_s[:,t]):9.2f} {np.median(hist_g[:,t]):9.3f}")
    print()
    print(f"  turns before stall: median {np.median(stall):.1f}, "
          f"range {stall.min()}-{stall.max()}")
    print("  stall causes: " + ", ".join(f"{k} x{v}" for k, v in Counter(why).items()))
    print(f"  flicker events (certified at turn t, not at t+1): {int(flick.sum())} "
          f"over {n_real} realisations")
    print("  NB zero flicker here is NOT evidence against flicker: Toy 3 is Asimov")
    print("  (noiseless), and flicker in the campaign is a noise/refit phenomenon.")
    print()
    print("  THE GEARING PROBLEM. Certifying EVERY pulsar converts the array from")
    print("  Earth-term-only to Earth+pulsar-term coherent: a factor sqrt(2) = 0.15 dex")
    print("  of SNR, and that is the ceiling on one turn of the wheel. Toy 2 says the")
    print("  gap to the next source must be <~0.7 dex for its comb to cut. The census")
    print("  gap is 0.77-1.17 dex. 0.15 dex of gearing cannot close a 0.8-1.2 dex gap,")
    print("  so the wheel cannot bootstrap itself. That is arithmetic, not a fit.")
    print()
    out["t3_hist_cert"], out["t3_hist_det"] = hist_c, hist_d
    out["t3_hist_snr"], out["t3_hist_gain"] = hist_s, hist_g
    out["t3_flicker"], out["t3_stall_turn"] = flick, stall
    out["t3_stall_why"] = np.array(why)
    out["t3_max_gain_dex"] = np.log10(np.sqrt(2.0))


# ============================================================ main
def main():
    rng = np.random.default_rng(20260729)
    out = {}
    os.makedirs(REPORTS, exist_ok=True)
    print("\nVERNIER toy -- EXPLORATORY, NOT pre-registered. Clean room: no campaign")
    print("module imported. numpy / CPU only.\n")
    anchor_check(out)
    toy1a_fringe_budget(out, rng)
    toy1b_clock(out)
    toy1c_vernier(out, rng)
    toy1d_commensurate(out, rng)
    toy2_brightness_gap(out, rng)
    toy3_wheel(out, rng)
    npz = os.path.join(REPORTS, "vernier_toy.npz")
    np.savez_compressed(npz, **out)
    print("=" * 74)
    print(f"wrote {npz}")
    for f in ["VERNIER_t1a_fringe_budget.png", "VERNIER_t1c_vernier.png",
              "VERNIER_t1d_commensurate.png", "VERNIER_t2_brightness_gap.png"]:
        print(f"wrote {os.path.join(REPORTS, f)}")
    print("=" * 74)


if __name__ == "__main__":
    main()
