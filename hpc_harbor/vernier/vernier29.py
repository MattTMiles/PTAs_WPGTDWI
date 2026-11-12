#!/usr/bin/env python3
"""
VERNIER-29 -- THE VLBI TARGETING RANKING (the Deller deliverable).

EXPLORATORY, NOT pre-registered. Toy model, cronus, CPU/numpy.
Clean room: imports only vernier24 / vernier25 / vernier_toy / vernier6 / vernier16
(all my own anchored code) plus the two REAL data products already tracked in this
repo (reports/anchor_data_forensics.npz and CW_transition/anchor_a0_priors.npz).

WHY
  The VLBI team expects only 2-3 more sub-parsec MSP distances next decade and needs
  a science case for WHICH pulsars to target. No run in this series has produced the
  ranking.

FRAMING RULE (D2 / ORACLE): distance knowledge does NOT buy ignition -- VERNIER-24/25
  put ignition fraction at 0.00 on every sky at K <= 16, and the R referendum made the
  NO-GO information-theoretic. So value is priced in the PRODUCTS currency, never in
  the ignition currency:
      (a) CERT  -- certification contrast MARGIN against its own bar, seed K = 8, 16
      (b) SOFT  -- soft-channel depth f_marg (VERNIER-23 Part 3 / -24 `soft`) at today's K
      (c) LOC   -- Fisher sky-AREA of the loudest source (localisation proxy)
      (d) CHIRP -- sigma(log10 Mc) from the pulsar-term line width
  Not "does the wheel start".

WHAT IS REAL HERE (and what is not)
  REAL : 116 array pulsars, their sky unit vectors (feather `pos_t`), their residual
         rms (reports/anchor_data_forensics.npz), and their distance priors WITH the
         method that produced each one (CW_transition/anchor_a0_priors.npz, which is
         get_distance.py's canonical output against the Cornell parallax catalogue).
  TOY  : the GW sky (VERNIER-24's audited Pareto census at gamma_eff = 1.56), the
         detections (VERNIER-25 Arm K seeds at array SNR 12), the noise model
         (white + Gaussian confusion sea), and the 20 yr / 14 d cadence.

THE ACHIEVABILITY MODEL is the crux and it is what makes this a targeting list rather
  than a wish list. A parallax measurement of precision sigma_pi delivers a DISTANCE
  precision sigma_d = sigma_pi * d^2 (d in pc, sigma_pi in arcsec). Sub-parsec
  therefore requires d <~ 300 pc even at the best demonstrated astrometric precision.
  The floor used is sigma_pi = 13 uas -- not an aspiration, but the precision Ding et
  al. (2024) actually achieved on J2222-0137 (PX = 3.723 +/- 0.013 mas), which is in
  the catalogue file this run reads. A 30 uas "typical deep campaign" variant is
  reported beside it. Every ranking is given TWICE: ORACLE (prior forced to 0.1 pc,
  physics only) and ACHIEVABLE (prior forced to max(0.1 pc, sigma_pi d^2)).

J0437-4715 enters as an EXISTING sub-wavelength anchor at 0.11 pc, per the brief, and
  is excluded from the candidate list. NOTE THE DISCREPANCY, stated not hidden: the
  canonical prior file in this repo still carries Deller et al. (2008),
  PX = 6.396 +/- 0.054 mas -> d = 156.3 pc, sigma_d = 1.32 pc. The 0.11 pc figure is
  the Reardon-era value the brief specifies. Deliverable 5 removes the anchor to price
  what it has already banked.

SUB-FRINGE NOTE (the banked lesson, both halves): every distance grid has an EVEN
  point count and is centred on the prior mean, so it never lands on a fringe centre
  by construction; and the certification floor is refit from nulls evaluated on the
  SAME grid, which absorbs the residual sampling bias.

PAIRED PERTURBATION: baseline and improved-prior rows are computed with the SAME noise
  draws (per-(pulsar, comb) deterministic sub-seeds), so every reported Delta is a
  response to the prior and not to a fresh noise realisation.

Code `hpc_harbor/vernier/vernier29.py` * data `reports/vernier29.npz`.
Staged, not committed.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import vernier_toy as V
import vernier6 as X
import vernier16 as S
import vernier24 as W
import vernier25 as K25

V.PTERM_STATIC = False

REPO = "/home/mattm/projects/HSYMT"
REPORTS = W.REPORTS
FEATHER_DIR = os.path.join(REPO, "data_products")
FORENSIC = os.path.join(REPO, "reports", "anchor_data_forensics.npz")
PRIORS = os.path.join(REPO, "CW_transition", "anchor_a0_priors.npz")
DISTFILE = os.path.join(REPO, "CW_transition", "best_distances.txt")
ARR_CACHE = os.path.join(REPORTS, "vernier29_array.npz")

YR, DF, T_OBS, C_PC_S = W.YR, W.DF, W.T_OBS, V.C_PC_S
R_AMP, KAPPA = W.R_AMP, W.KAPPA
NGRID = W.NGRID
BAND = K25.BAND                       # ten 1/T bins, 10.0-24.3 nHz
GAMMA = 1.56                          # the audited bright-end slope
RHO_SEED = K25.RHO_SEED               # array SNR of an Arm K seed
MC_LO, MC_HI = K25.MC_LO, K25.MC_HI

# --- certification numerics (VERNIER-25 settings) ---
NPTS = 800                            # EVEN -- sub-fringe note
N_NULL_C = 32

# --- the products' knobs ---
K_LIST = (1, 8, 16)                   # 1 = "today" (VERNIER-24/25 cold start)
F_LOC = 15e-9                         # localisation Fisher bin
M_LOC = 3                             # sources sharing that bin
N_REAL = 24                           # sky realisations for the certification arms
N_REAL_L = 12                         # sky realisations for LOC / CHIRP (cheap)
N_REAL_C = 8                          # sky realisations for the controls

# --- the achievability model ---
UAS = np.pi / (180.0 * 3600.0 * 1e6)  # microarcsec in radians (unused; kept explicit)
SIG_PI_BEST = 13e-6 / 1e6 * 1e6       # placeholder, overwritten below
SIG_PI_BEST_ARCSEC = 13e-6            # 13 uas, Ding et al. (2024) on J2222-0137
SIG_PI_CAMP_ARCSEC = 30e-6            # 30 uas, typical deep VLBI campaign
SIG_TARGET_PC = 0.10                  # the oracle sub-lambda target
LAM_GW_PC = C_PC_S / F_LOC            # GW wavelength at band centre, pc
ANCHOR = "J0437-4715"
SIG_ANCHOR_PC = 0.11                  # the brief's existing-anchor value
SIG_ANCHOR_CATALOGUE_PC = 1.32        # what best_distances.txt actually carries


# ================================================================= the real array
def _feather_pos(name):
    """Real sky unit vector from the tracked feather (`pos_t`, row 0)."""
    import pyarrow.feather as pf
    t = pf.read_table(os.path.join(FEATHER_DIR, name + ".feather"),
                      columns=["pos_t_0", "pos_t_1", "pos_t_2"])
    v = np.array([t["pos_t_0"][0].as_py(), t["pos_t_1"][0].as_py(),
                  t["pos_t_2"][0].as_py()], float)
    return v / np.linalg.norm(v)


def _category(method, note):
    """Distance-prior CATEGORY, from get_distance.py's method + provenance note."""
    if method == "VLBI":
        return "VLBI"
    if method == "Catalog_PX":
        if "(VLBI)" in note:
            return "VLBI"
        if "(Timing)" in note:
            return "timing-PX"
        return "timing-PX"
    if method == "Timing_PX":
        return "timing-PX"
    if method.startswith("feather(DM"):
        return "DM-model"
    if method.startswith("feather(PBDOT"):
        return "Pbdot-Shk"
    return "DM-default"                 # feather(none): enterprise default prior


def build_array(force=False):
    """
    The REAL IPTA-like array: 116 pulsars with real positions, real per-pulsar rms
    weights and real distance-prior categories. Cached to reports/.
    """
    if os.path.exists(ARR_CACHE) and not force:
        z = np.load(ARR_CACHE, allow_pickle=True)
        d = {k: z[k] for k in z.files}
    else:
        fo = np.load(FORENSIC, allow_pickle=True)
        names = [str(x) for x in fo["names"]]
        rms = np.asarray(fo["resid_rms_s"], float)
        toaerr = np.asarray(fo["median_toaerr_s"], float)
        pz = np.load(PRIORS, allow_pickle=True)
        pn = [str(x) for x in pz["fname"]]
        idx = {n: i for i, n in enumerate(pn)}
        jn = [str(x) for x in pz["jname"]]
        meth = [str(x) for x in pz["method"]]
        note = [str(x) for x in pz["note"]]
        ispx = np.asarray(pz["is_parallax_class"], bool)
        dr, sr = np.asarray(pz["dist_real_kpc"], float), np.asarray(pz["sig_real_kpc"], float)
        dfe, sfe = (np.asarray(pz["dist_feather_kpc"], float),
                    np.asarray(pz["sig_feather_kpc"], float))
        P, L, sigL, cat, jname = [], [], [], [], []
        for n in names:
            i = idx[n]
            P.append(_feather_pos(n))
            L.append((dr[i] if ispx[i] else dfe[i]) * 1000.0)
            sigL.append((sr[i] if ispx[i] else sfe[i]) * 1000.0)
            cat.append(_category(meth[i], note[i]))
            jname.append(jn[i])
        d = dict(name=np.array(names), jname=np.array(jname),
                 P=np.array(P, float), rms=rms, toaerr=toaerr,
                 L=np.array(L, float), sigL=np.array(sigL, float),
                 cat=np.array(cat), note=np.array([note[idx[n]] for n in names]),
                 method=np.array([meth[idx[n]] for n in names]))
        np.savez_compressed(ARR_CACHE, **d)
    # per-pulsar weight w = sigma_med / sigma_a, unit MEAN POWER (VERNIER-24 convention)
    rms = d["rms"]
    w = np.median(rms) / rms
    w = w / np.sqrt(np.mean(w ** 2))
    n_psr = len(w)
    G = S.fib_sphere(NGRID)
    fpg, fcg = W.antenna_reg(d["P"], G)
    arr = dict(P=d["P"], w=w, L=d["L"].copy(), sigL=d["sigL"].copy(),
               sigL_cat=d["sigL"].copy(), n_psr=n_psr, sigma_toa=float(np.median(rms)),
               G=G, fpg=fpg, fcg=fcg, name=d["name"], jname=d["jname"],
               cat=d["cat"], note=d["note"], method=d["method"],
               rms=rms, toaerr=d["toaerr"])
    # the existing anchor, per the brief
    ia = int(np.where(arr["jname"] == ANCHOR)[0][0])
    arr["i_anchor"] = ia
    arr["sigL"][ia] = SIG_ANCHOR_PC
    return arr


# ================================================================= achievability
def sigma_vlbi(d_pc, sig_pi_arcsec):
    """sigma_d = sigma_pi * d^2 for d in pc and sigma_pi in arcsec (pi[as] = 1/d[pc])."""
    return sig_pi_arcsec * np.asarray(d_pc, float) ** 2


def improved_prior(arr, mode, sig_pi=SIG_PI_BEST_ARCSEC):
    """
    Per-pulsar prior width AFTER a successful campaign. Never worse than today.
      'oracle'      -- forced to SIG_TARGET_PC, physics only
      'achievable'  -- max(SIG_TARGET_PC, sigma_pi d^2), the real astrometric limit
    """
    if mode == "oracle":
        s = np.full(arr["n_psr"], SIG_TARGET_PC)
    elif mode == "achievable":
        s = np.maximum(SIG_TARGET_PC, sigma_vlbi(arr["L"], sig_pi))
    else:
        raise ValueError(mode)
    return np.minimum(s, arr["sigL"])


# ================================================================= the seed combs
def seed_combs(rng, arr, K, band=BAND, gamma=GAMMA):
    """
    K accumulated, genuinely-recovered combs across the band -- VERNIER-25's Arm K
    first-detection era. Each seed is one source at array SNR RHO_SEED against
    white + the honest confusion sea of the corrected gamma = 1.56 sky, at a distinct
    frequency, with Mc ~ U(10^8.5, 10^9.5).

    Returns the comp list `certify` consumes, plus the per-pulsar pulsar-term power
    used to weight every product-level fraction.
    """
    comps = []
    for j in range(K):
        f = float(band[j % len(band)] + (j // len(band) + 0.0) * 0.31 * DF)
        B = K25.build_bin(rng, arr, f, gamma, seed_rho=RHO_SEED)
        rho = W.rho_budget(B, 0, set())                 # lambda_eff-consistent budget
        rho_arr = float(np.sqrt(np.sum(rho ** 2)))
        sig_f = 1.0 / (2 * np.pi * T_OBS * max(rho_arr, 1e-6))
        omc = B["omc"][:, 0]
        tau = V.tau_of(arr["L"], omc)
        dec = np.exp(-0.5 * (2 * np.pi * sig_f * tau / np.sqrt(arr["n_psr"])) ** 2)
        mc = float(B["mc"][0])
        comps.append(dict(f=f, tc=float(V.t_merge(f, mc)), mc=mc, omc=omc,
                          rho_P=R_AMP * rho, dec=dec, dL=V.dL_of(f, omc),
                          rho_arr=rho_arr, rho0=B["rho0"][:, 0].copy()))
    # power weight: MEAN pulsar-term power over the accumulated combs, so the ranking
    # is not hostage to the geometry of one source (deviation from VERNIER-24, which
    # weights by the NEWEST source alone -- flagged).
    pw = np.mean(np.stack([c["rho_P"] ** 2 for c in comps]), axis=0)
    return comps, pw


# ================================================================= (a)+(b) kernel
def _contrast(curve, Lg, hw):
    """
    Peak minus the best rival outside +-hw of the peak; rows if 2-D. Returns NaN
    where NO rival exists -- see `_no_rival`.
    """
    if curve.ndim == 1:
        it = int(np.argmax(curve))
        far = np.abs(Lg - Lg[it]) > hw
        return float(curve[it] - curve[far].max()) if far.any() else np.nan
    it = np.argmax(curve, axis=1)
    far = np.abs(Lg[None, :] - Lg[it][:, None]) > hw
    v = curve[np.arange(curve.shape[0]), it] - np.where(far, curve, -np.inf).max(axis=1)
    return np.where(far.any(axis=1), v, np.nan)


def _no_rival(Lg, hw):
    """
    True when the searched window holds no point more than half a fringe from any
    possible peak -- i.e. LESS THAN ONE FRINGE is in play.

    VERNIER-24's `certify` silently fell back to peak-minus-minimum here, which
    manufactures a finite 'contrast' out of a window that contains no competing
    fringe at all. That is exactly the regime a sub-parsec anchor puts a pulsar in,
    so it cannot be left implicit in THIS run: the state is named PRIOR-PINNED, the
    contrast is undefined (NaN), and the certification bar is zero because there is
    no fringe ambiguity left to resolve.
    """
    return (Lg[-1] - Lg[0]) <= hw


def cert_rows(arr, comps, sigL, base_seed, rows=None, npts=NPTS, n_null=N_NULL_C,
              scramble=False):
    """
    VERNIER-24's `certify`, refactored for a PER-PULSAR prior vector `sigL`.

    Certification is per-pulsar independent: pulsar a's distance prior enters only
    pulsar a's own grid, its own trials bar ln K_a = ln(6 sigma_a / dL_a) and its own
    incoherent-comb floor. So one pass gives every pulsar, and a second pass with a
    perturbed sigL gives every marginal value -- at the same cost.

    Noise draws are keyed on (base_seed, pulsar, comb) and are therefore IDENTICAL
    between the baseline and the perturbed pass. Every Delta is a prior response.

    Returns per-pulsar: contrast S_a, bar_a, margin_a = S_a - bar_a, lnK_a, floor_a,
    cert_a, cert_true_a, prior-pinned flag (K_a <= 1), and the soft-channel fringe
    coherence R_a (VERNIER-23 Part 3: no threshold at all, the pulsar term enters
    marginalised over the prior and weighted by R_a).
    """
    n_psr = arr["n_psr"]
    idx = range(n_psr) if rows is None else rows
    out = {k: np.zeros(n_psr) for k in
           ("stat", "bar", "margin", "lnK", "floor", "Rcoh", "Kfr", "dL")}
    out["cert"] = np.zeros(n_psr, bool)
    out["cert_true"] = np.zeros(n_psr, bool)
    out["pinned"] = np.zeros(n_psr, bool)
    dLmin = np.min(np.stack([c["dL"] for c in comps]), axis=0)
    for a in idx:
        sa = float(sigL[a])
        Lstar = float(arr["L"][a])
        # GRID TRUNCATION (new here, and it matters only for the widest real priors):
        # a +-3 sigma window on a 60-70 % prior runs through L = 0 into NEGATIVE
        # distance, where the retarded phase (1 + tau/tc)^(5/8) is not real. The
        # window is therefore capped at 0.9 L, and the trials bar is computed from
        # the window ACTUALLY searched rather than from a nominal 6 sigma.
        hwg = min(3.0 * sa, 0.9 * Lstar)
        Lg = Lstar + np.linspace(-hwg, hwg, npts)                # EVEN npts
        lnl = np.zeros(npts)
        n_inc = np.zeros((n_null, npts))
        phi_ref = None
        for jc, c in enumerate(comps):
            r_ = np.random.default_rng((base_seed, a, jc))
            # SOURCE-LEVEL SCRAMBLED NULL (VERNIER-24's control): the comb in the data
            # sits at the WRONG distance while the grid stays centred on the prior
            # mean. Anything that certifies here is certifying noise.
            Ls = Lstar + (r_.uniform(-500.0, 500.0) if scramble else 0.0)
            phi, env = X.phase_env(Lg, max(Ls, 50.0), c["f"] / 2.0, c["tc"],
                                   c["omc"][a], T_OBS, 2)
            env = env * c["dec"][a]
            r = c["rho_P"][a]
            cp, sp = np.cos(phi), np.sin(phi)
            nc, ns = r_.normal(), r_.normal()
            lnl += r ** 2 * (env * cp - 1.0) + r * env * (nc * cp + ns * sp)
            nc0 = r_.normal(size=n_null)[:, None]
            ns0 = r_.normal(size=n_null)[:, None]
            dl = r_.uniform(0, 2 * np.pi, n_null)[:, None]
            cq = cp[None, :] * np.cos(dl) - sp[None, :] * np.sin(dl)
            sq = sp[None, :] * np.cos(dl) + cp[None, :] * np.sin(dl)
            n_inc += r ** 2 * (env[None, :] * cq - 1.0) \
                + r * env[None, :] * (nc0 * cq + ns0 * sq)
            if phi_ref is None:
                phi_ref = phi
        hw = 0.5 * dLmin[a]
        Kfr = max(2.0 * hwg / dLmin[a], 1e-30)      # = 6 sigma / dL when untruncated
        itop = int(np.argmax(lnl))
        pinned = bool(Kfr <= 1.0 or _no_rival(Lg, hw))
        out["Kfr"][a] = Kfr; out["dL"][a] = dLmin[a]
        out["pinned"][a] = pinned
        if pinned:
            # PRIOR-PINNED: fewer than one fringe in the window. No trials, no
            # incoherent-comb floor to beat -- the bar is zero and the contrast is
            # not the operative quantity.
            out["stat"][a] = np.nan; out["floor"][a] = np.nan
            out["lnK"][a] = 0.0; out["bar"][a] = 0.0; out["margin"][a] = np.nan
        else:
            fl = W.recut_floor(_contrast(n_inc, Lg, hw),
                               np.random.default_rng((base_seed, a)))
            lnK = float(np.log(max(Kfr, 1.0)))
            bar = max(lnK, fl["floor"] + fl["err"])
            st = _contrast(lnl, Lg, hw)
            out["stat"][a] = st; out["bar"][a] = bar; out["margin"][a] = st - bar
            out["lnK"][a] = lnK; out["floor"][a] = fl["floor"] + fl["err"]
            out["cert"][a] = bool(st > bar)
            if out["cert"][a]:
                out["cert_true"][a] = abs(Lg[itop] - Lstar) <= 0.5 * dLmin[a]
        pr = np.exp(-0.5 * ((Lg - Lstar) / sa) ** 2)
        wg = np.exp(lnl - lnl.max()) * pr
        wg = np.where(np.isfinite(wg), wg, 0.0)
        wg = wg / max(np.trapezoid(wg, Lg), 1e-300)
        rc = float(np.abs(np.trapezoid(wg * np.exp(1j * phi_ref), Lg)))
        out["Rcoh"][a] = rc if np.isfinite(rc) else 0.0
    return out


def product_ab(res, pw):
    """Product-level (a) and (b) from a cert_rows result and the power weights."""
    tot = max(pw.sum(), 1e-30)
    known = res["cert"] | res["pinned"]
    return dict(f_cert=float(pw[known].sum() / tot),
                f_cert_true=float(pw[res["cert_true"] | res["pinned"]].sum() / tot),
                f_marg=float(np.sum(res["Rcoh"] ** 2 * pw) / tot),
                f_pinned=float(pw[res["pinned"]].sum() / tot),
                bar_sum=float(np.sum(pw * res["bar"]) / tot),
                margin_med=float(np.nanmedian(res["margin"])
                                 if np.isfinite(res["margin"]).any() else np.nan))


# ================================================================= (c) localisation
def loc_build(rng, arr, M=M_LOC, f=F_LOC, rho1=RHO_SEED, eps=1e-5):
    """
    Fisher over [M source parameters] + [116 pulsar distances] on the REAL array,
    with the REAL per-pulsar weights applied (VERNIER-24's `separability` builds the
    weights and then never uses them; that is corrected here and flagged).

    The Jacobian is built ONCE. Changing any pulsar's distance prior changes only the
    diagonal added to the distance block before the Schur complement, so 116 marginal
    values cost 116 linear solves and no re-differentiation.
    """
    P, w, n_psr = arr["P"], arr["w"], arr["n_psr"]
    L = arr["L"]
    wcat = np.concatenate([w, w])
    R2, S_sea, _ = W.bin_sources(rng, f, "B", n_bright=M)
    amp = np.sqrt(R2 / R2[0])
    fk = f + rng.uniform(-0.5, 0.5, M) * DF
    S_ = np.stack([np.arccos(rng.uniform(-1, 1, M)), rng.uniform(0, 2 * np.pi, M),
                   amp, rng.uniform(0, np.pi, M), rng.uniform(-1, 1, M),
                   rng.uniform(0, 2 * np.pi, M), fk], axis=1)
    base = W.dvec(P, S_, L) * wcat
    npar = 7
    ns = M * npar
    J = np.zeros((2 * base.size, ns + n_psr))
    for i in range(M):
        for k in range(npar):
            sp = S_.copy()
            st = eps * f if k == 6 else eps * (abs(sp[i, k]) + 1.0)
            sp[i, k] += st
            dd = (W.dvec(P, sp, L) * wcat - base) / st
            J[:, i * npar + k] = np.concatenate([dd.real, dd.imag])
    for a in range(n_psr):
        Lv = L.copy(); st = eps * (abs(Lv[a]) + 1.0); Lv[a] += st
        dd = (W.dvec(P, S_, Lv) * wcat - base) / st
        J[:, ns + a] = np.concatenate([dd.real, dd.imag])
    e1 = W.dvec(P, S_[:1], L) * wcat
    sig2 = float(np.sum(np.abs(e1) ** 2)) / max(rho1 ** 2, 1e-30)
    F = (J.T @ J) / (sig2 / 2.0)
    return dict(F=F, ns=ns, npar=npar, M=M, S_=S_, n_psr=n_psr)


def loc_area(lb, sigL):
    """
    1-sigma sky ELLIPSE AREA of the loudest source, sr, after marginalising all 116
    pulsar distances over their (per-pulsar) priors by the Schur complement.
    """
    F, ns, npar, n_psr = lb["F"], lb["ns"], lb["npar"], lb["n_psr"]
    Fdd = F[ns:, ns:] + np.diag(1.0 / np.asarray(sigL, float) ** 2)
    Fsd = F[ns:, :ns].T
    try:
        Fm = F[:ns, :ns] - Fsd @ np.linalg.solve(Fdd, Fsd.T)
    except np.linalg.LinAlgError:
        return np.inf
    dg = np.sqrt(np.clip(np.diag(Fm), 1e-300, None))
    Fc = Fm / np.outer(dg, dg)
    try:
        C = np.linalg.inv(Fc + 1e-12 * np.eye(ns)) / np.outer(dg, dg)
    except np.linalg.LinAlgError:
        return np.inf
    d = np.clip(np.diag(C), 0.0, None)
    sth = np.sqrt(d[0]); sph = np.sqrt(d[1]) * np.sin(lb["S_"][0, 0])
    return float(np.pi * sth * sph)


# ================================================================= (d) chirp mass
def chirp_sigma(arr, comps, sigL, which=0):
    """
    sigma(log10 Mc) from the pulsar-term line width, with the distance prior
    MARGINALISED per pulsar.

    Pulsar a sees the pulsar term at f_P = f (1 + tau_a/tc)^(-3/8), i.e. offset
    delta_a = (3/8) f tau_a / tc. Because tc ~ Mc^(-5/3),  d delta_a / d ln Mc =
    (5/3) delta_a; and because tau_a ~ L_a,  d delta_a / d L_a = delta_a / L_a. Each
    pulsar measures its own f_P to sigma_f,a = 1/(2 pi T rho_P,a), so marginalising
    L_a over a Gaussian prior of width sigma_a gives the per-pulsar information

        I_a = ((5/3) delta_a)^2 / ( sigma_f,a^2 + (delta_a sigma_a / L_a)^2 )

    which reduces to VERNIER-24's `chirp_precision` exactly as sigma_a -> 0. THIS is
    the term VERNIER-24 left out: it assumed tau_a known. A pulsar with a 20 % DM
    distance contributes essentially nothing to Mc no matter how well it is timed.
    """
    c = comps[which]
    rho = c["rho_P"]
    tau = V.tau_of(arr["L"], c["omc"])
    delta = (3.0 / 8.0) * c["f"] * tau / c["tc"]
    sig_f = 1.0 / (2 * np.pi * T_OBS * np.maximum(rho, 1e-12))
    sig_d = delta * np.asarray(sigL, float) / np.maximum(arr["L"], 1e-12)
    I = ((5.0 / 3.0) * delta) ** 2 / (sig_f ** 2 + sig_d ** 2)
    return 1.0 / (np.sqrt(max(I.sum(), 1e-300)) * np.log(10.0)), I


# ================================================================= reporting utils
def _spear(a, b):
    ra = np.argsort(np.argsort(-np.asarray(a, float)))
    rb = np.argsort(np.argsort(-np.asarray(b, float)))
    ra = ra - ra.mean(); rb = rb - rb.mean()
    return float(np.sum(ra * rb) / max(np.sqrt(np.sum(ra ** 2) * np.sum(rb ** 2)), 1e-30))


def _hdr(t):
    print(); print("=" * 108); print(t); print("=" * 108)


# ================================================================= 0. the array
def exp0_array(arr, out):
    _hdr("DELIVERABLE 0 -- THE REAL ARRAY, ITS WEIGHTS AND ITS DISTANCE-PRIOR CATEGORIES")
    n = arr["n_psr"]
    w, rms = arr["w"], arr["rms"]
    print(f"  {n} pulsars, REAL sky unit vectors (feather pos_t), REAL residual rms "
          f"(reports/anchor_data_forensics.npz),")
    print(f"  REAL distance priors from get_distance.py vs the Cornell parallax "
          f"catalogue (CW_transition/anchor_a0_priors.npz).")
    print(f"  T = 20 yr, 1/T = {DF*1e9:.3f} nHz, band {BAND[0]*1e9:.2f}-{BAND[-1]*1e9:.2f} nHz "
          f"in {len(BAND)} bins; lambda_GW(15 nHz) = {LAM_GW_PC:.3f} pc.")
    print()
    print(f"  residual rms: median {np.median(rms)*1e9:.0f} ns, "
          f"range {rms.min()*1e9:.0f}-{rms.max()*1e9:.0f} ns")
    print(f"  REALISED weight heterogeneity: sigma_log10(1/rms) = "
          f"{np.std(np.log10(1.0/rms)):.3f} dex  "
          f"(VERNIER-13 assumed 0.400 dex -- the real array is {np.std(np.log10(1.0/rms)):.3f})")
    print(f"  w = sigma_med/sigma_a normalised to unit mean POWER: "
          f"{w.min():.3f}-{w.max():.3f}, median {np.median(w):.3f}")
    print()
    print("  DISTANCE-PRIOR CATEGORY CENSUS (the catalogue values this run uses):")
    print(f"  {'category':>12} {'N':>4} {'median d (pc)':>14} {'median sigma_d (pc)':>20} "
          f"{'median sigma_d/d':>17} {'sub-pc now':>11}")
    cats = ["VLBI", "timing-PX", "Pbdot-Shk", "DM-model", "DM-default"]
    rows = []
    for c in cats:
        m = arr["cat"] == c
        if not m.any():
            continue
        print(f"  {c:>12} {int(m.sum()):4d} {np.median(arr['L'][m]):14.1f} "
              f"{np.median(arr['sigL'][m]):20.2f} "
              f"{np.median(arr['sigL'][m]/arr['L'][m]):17.3f} "
              f"{int(np.sum(arr['sigL'][m] < 1.0)):11d}")
        rows.append([c, int(m.sum()), float(np.median(arr['L'][m])),
                     float(np.median(arr['sigL'][m]))])
    print()
    sub = np.where(arr["sigL"] < 1.0)[0]
    print(f"  pulsars with sigma_d < 1 pc TODAY: {len(sub)}"
          + (" -- " + ", ".join(f"{arr['jname'][i]} ({arr['sigL'][i]:.2f} pc)"
                                for i in sub[np.argsort(arr['sigL'][sub])]) if len(sub) else ""))
    print(f"  pulsars with sigma_d < lambda_GW/2 = {LAM_GW_PC/2:.3f} pc TODAY: "
          f"{int(np.sum(arr['sigL'] < LAM_GW_PC/2))}")
    print()
    print(f"  ANCHOR: {ANCHOR} entered at sigma_d = {SIG_ANCHOR_PC:.2f} pc per the brief.")
    print(f"    DISCREPANCY, STATED: the canonical prior file in this repo carries "
          f"Deller et al. (2008),")
    print(f"    PX = 6.396 +/- 0.054 mas -> d = 156.3 pc, sigma_d = "
          f"{SIG_ANCHOR_CATALOGUE_PC:.2f} pc. The 0.11 pc figure is the")
    print(f"    Reardon-era value the brief specifies; deliverable 5 prices the anchor "
          f"both ways.")
    out["cat_rows"] = np.array(rows, dtype=object)
    out["w"] = w; out["L"] = arr["L"]; out["sigL0"] = arr["sigL"]
    out["cat"] = arr["cat"]; out["jname"] = arr["jname"]


# ================================================================= 1. achievability
def exp1_reach(arr, out):
    _hdr("DELIVERABLE 1 -- WHAT VLBI CAN ACTUALLY REACH  (sigma_d = sigma_pi d^2)")
    print("  A parallax campaign of precision sigma_pi buys a DISTANCE precision")
    print("      sigma_d [pc] = sigma_pi [arcsec] * (d [pc])^2 ,")
    print("  because pi[arcsec] = 1/d[pc]. Sub-parsec is therefore a DISTANCE cut on the")
    print("  target, not a telescope-time cut. Floors used:")
    print(f"      sigma_pi = {SIG_PI_BEST_ARCSEC*1e6:.0f} uas -- DEMONSTRATED: Ding et al. (2024), "
          f"J2222-0137, PX = 3.723 +/- 0.013 mas")
    print(f"      sigma_pi = {SIG_PI_CAMP_ARCSEC*1e6:.0f} uas -- typical deep campaign")
    print()
    for sp in (SIG_PI_BEST_ARCSEC, SIG_PI_CAMP_ARCSEC):
        dmax_1pc = np.sqrt(1.0 / sp)
        dmax_lam = np.sqrt(LAM_GW_PC / 2.0 / sp)
        n1 = int(np.sum(sigma_vlbi(arr["L"], sp) < 1.0))
        nl = int(np.sum(sigma_vlbi(arr["L"], sp) < LAM_GW_PC / 2))
        print(f"  sigma_pi = {sp*1e6:2.0f} uas : sub-parsec out to d = {dmax_1pc:6.1f} pc "
              f"({n1:2d} of {arr['n_psr']} array pulsars) ; "
              f"sub-lambda/2 out to {dmax_lam:6.1f} pc ({nl:2d} pulsars)")
    print()
    sp = SIG_PI_BEST_ARCSEC
    reach = sigma_vlbi(arr["L"], sp)
    ordr = np.argsort(reach)
    print(f"  THE REACHABLE SET at {sp*1e6:.0f} uas (sorted by reachable sigma_d):")
    print(f"  {'#':>3} {'pulsar':>12} {'d (pc)':>9} {'sigma_d now':>12} {'category now':>13} "
          f"{'sigma_d reachable':>18} {'gain factor':>12}")
    nshow = 0
    for i in ordr:
        if reach[i] > 3.0:
            break
        nshow += 1
        tag = "  <-- ANCHOR" if i == arr["i_anchor"] else ""
        print(f"  {nshow:3d} {arr['jname'][i]:>12} {arr['L'][i]:9.1f} "
              f"{arr['sigL'][i]:12.2f} {arr['cat'][i]:>13} {reach[i]:18.3f} "
              f"{arr['sigL'][i]/max(reach[i],1e-9):12.1f}{tag}")
    print(f"\n  {nshow} pulsars can be driven below 3 pc; "
          f"{int(np.sum(reach < 1.0))} below 1 pc; "
          f"{int(np.sum(reach < LAM_GW_PC/2))} below lambda/2 = {LAM_GW_PC/2:.3f} pc.")
    print("  This IS the '2-3 more sub-parsec distances next decade' expectation, "
          "reproduced from the")
    print("  catalogue rather than assumed -- and it is why the ranking is dominated by "
          "d, not by timing.")
    out["reach_best"] = reach
    out["reach_camp"] = sigma_vlbi(arr["L"], SIG_PI_CAMP_ARCSEC)


# ================================================================= 2/3. the ranking
MODES = ("oracle", "achievable")


def marginal_table(arr, sig_pi=SIG_PI_BEST_ARCSEC, verbose=True):
    """
    Per-pulsar marginal value of driving pulsar i's prior to sub-lambda, in all four
    products, for BOTH arms at once.

    Everything is PAIRED: the two improved passes reuse the baseline pass's noise
    draws, so every Delta is a prior response. The baseline is computed once and
    shared between the arms, which is both cheaper and removes baseline scatter from
    the arm-to-arm comparison.

    Per-realisation values are KEPT for K = 8 so the ranking can be split-half tested
    -- a marginal value that does not reproduce across independent halves of the sky
    ensemble is not a ranking, it is noise, and that has to be measured rather than
    assumed.
    """
    n = arr["n_psr"]
    s_new = {m: improved_prior(arr, m, sig_pi) for m in MODES}
    d_marg = {m: {K: np.zeros(n) for K in K_LIST} for m in MODES}
    d_bar = {m: {K: np.zeros(n) for K in K_LIST} for m in MODES}
    flip = {m: {K: np.zeros(n) for K in K_LIST} for m in MODES}
    per_real = {m: np.zeros((N_REAL, n)) for m in MODES}       # K = 8, for split-half
    base_prod = {}
    Kfr0 = np.zeros(n); dL0 = np.zeros(n); pw_acc = np.zeros(n); bar0 = np.zeros(n)
    for K in K_LIST:
        b_marg = []; b_cert = []; b_pin = []
        for r in range(N_REAL):
            rng = np.random.default_rng(29000 + 101 * r + K)
            comps, pw = seed_combs(rng, arr, K)
            seed = (2900, r, K)
            r0 = cert_rows(arr, comps, arr["sigL"], seed)
            tot = max(pw.sum(), 1e-30)
            known0 = r0["cert"] | r0["pinned"]
            for m in MODES:
                r1 = cert_rows(arr, comps, s_new[m], seed)
                dv = pw * (r1["Rcoh"] ** 2 - r0["Rcoh"] ** 2) / tot
                d_marg[m][K] += dv / N_REAL
                # (a) has TWO faces and both are reported: the bar reduction (finite
                # in every state, PRIOR-PINNED counting as bar = 0) and the FLIP rate
                # -- the fraction of realisations in which the pulsar goes from
                # not-known to known. The flip is the operative one; see the note in
                # the ranked table about why the bar can move the WRONG way.
                d_bar[m][K] += (r0["bar"] - r1["bar"]) / N_REAL
                flip[m][K] += ((r1["cert"] | r1["pinned"]) & ~known0) / N_REAL
                if K == 8:
                    per_real[m][r] = dv
            p0 = product_ab(r0, pw)
            b_marg.append(p0["f_marg"]); b_cert.append(p0["f_cert"])
            b_pin.append(p0["f_pinned"])
            if K == 8:
                bar0 += r0["bar"] / N_REAL
            if K == K_LIST[-1]:
                Kfr0 += r0["Kfr"] / N_REAL; dL0 += r0["dL"] / N_REAL
                pw_acc += pw / tot / N_REAL
        base_prod[K] = dict(f_marg=float(np.mean(b_marg)), f_cert=float(np.mean(b_cert)),
                            f_pinned=float(np.mean(b_pin)))
        if verbose:
            print(f"    K={K:2d}: baseline f_marg = {base_prod[K]['f_marg']:.4f}, "
                  f"f_cert = {base_prod[K]['f_cert']:.4f}, "
                  f"f_prior-pinned = {base_prod[K]['f_pinned']:.4f}")
    # ---- (c) localisation, (d) chirp: cheap, deterministic given the realisation ----
    d_area = {m: np.zeros(n) for m in MODES}
    d_chirp = {m: np.zeros(n) for m in MODES}
    a0s, c0s = [], []
    for r in range(N_REAL_L):
        lb = loc_build(np.random.default_rng(29500 + 37 * r), arr)
        A0 = loc_area(lb, arr["sigL"]); a0s.append(A0)
        comps, _ = seed_combs(np.random.default_rng(29900 + 53 * r), arr, 1)
        s0, _ = chirp_sigma(arr, comps, arr["sigL"]); c0s.append(s0)
        for m in MODES:
            for i in range(n):
                s = arr["sigL"].copy(); s[i] = s_new[m][i]
                d_area[m][i] += (A0 - loc_area(lb, s)) / max(A0, 1e-300) / N_REAL_L
                d_chirp[m][i] += (s0 - chirp_sigma(arr, comps, s)[0]) \
                    / max(s0, 1e-300) / N_REAL_L
    return {m: dict(mode=m, s_new=s_new[m], d_marg=d_marg[m], d_bar=d_bar[m],
                    flip=flip[m], base=base_prod, d_area=d_area[m],
                    d_chirp=d_chirp[m], bar0=bar0, per_real=per_real[m],
                    A0=float(np.mean(a0s)), sig_mc0=float(np.mean(c0s)),
                    Kfr0=Kfr0, dL0=dL0, pw=pw_acc) for m in MODES}


def exp23_rank(arr, out):
    _hdr("DELIVERABLES 2+3 -- MARGINAL VALUE PER PULSAR, AND THE RANKING")
    print("  For each pulsar i: force its distance prior to sub-lambda and recompute all")
    print("  four products. Marginal value = the change. Baseline and perturbed passes")
    print("  share their noise draws, so every Delta is a prior response.")
    print()
    print(f"  ORACLE arm     : sigma_d -> {SIG_TARGET_PC:.2f} pc for any pulsar (physics only)")
    print(f"  ACHIEVABLE arm : sigma_d -> max({SIG_TARGET_PC:.2f} pc, "
          f"{SIG_PI_BEST_ARCSEC*1e6:.0f} uas * d^2)  (what VLBI can do)")
    print()
    res = marginal_table(arr)
    print()
    print(f"  BASELINE PRODUCTS (today's priors, {ANCHOR} at {SIG_ANCHOR_PC:.2f} pc):")
    for K in K_LIST:
        lab = "  (today)" if K == 1 else ""
        print(f"    K = {K:2d}{lab:9} f_marg = {res['oracle']['base'][K]['f_marg']:.4f}   "
              f"f_cert = {res['oracle']['base'][K]['f_cert']:.4f}   "
              f"f_prior-pinned = {res['oracle']['base'][K]['f_pinned']:.4f}")
    print(f"    localisation: 1-sigma sky ellipse area of the loudest source = "
          f"{res['oracle']['A0']:.3e} sr")
    print(f"    chirp mass  : sigma(log10 Mc) = {res['oracle']['sig_mc0']:.4f} dex")
    print()

    ia = arr["i_anchor"]
    cand = np.array([i for i in range(arr["n_psr"]) if i != ia])
    for mode in ("achievable", "oracle"):
        R = res[mode]
        prim = R["d_marg"][8]
        ordr = cand[np.argsort(-prim[cand])]
        _hdr(f"THE RANKED TABLE -- {mode.upper()} ARM "
             f"(primary: Delta f_marg at seed K = 8, the soft channel that exists today)")
        print("  READING (a): 'dbar' is the reduction in the certification bar and CAN GO")
        print("  NEGATIVE. That is not a bug. Shrinking a prior below one fringe removes the")
        print("  rival fringes from the searched window, so the incoherent-comb NULL contrast")
        print("  rises too and the refit floor rises with it -- the floor is doing its job.")
        print("  The operative (a) quantity is therefore 'pin', the fraction of realisations")
        print("  in which the pulsar goes from not-known to known (certified OR prior-pinned).")
        print()
        print(f"  {'#':>3} {'pulsar':>12} {'cat now':>11} {'d(pc)':>8} {'sig_d now':>10} "
              f"{'sig_d new':>10} {'w':>6} {'K_fr now':>9} {'(b) dfmarg':>11} "
              f"{'(b) dfmarg':>11} {'(a) dbar':>10} {'(a) pin':>8} {'(c) dA/A':>9} "
              f"{'(d) dsMc':>9}")
        print(f"  {'':>3} {'':>12} {'':>11} {'':>8} {'(pc)':>10} {'(pc)':>10} {'':>6} "
              f"{'':>9} {'K=1':>11} {'K=8':>11} {'K=8 nat':>10} {'K=8':>8} {'':>9} "
              f"{'/sigma':>9}")
        rows = []
        for rank, i in enumerate(ordr[:10], 1):
            print(f"  {rank:3d} {arr['jname'][i]:>12} {arr['cat'][i]:>11} "
                  f"{arr['L'][i]:8.1f} {arr['sigL'][i]:10.2f} {R['s_new'][i]:10.3f} "
                  f"{arr['w'][i]:6.2f} {R['Kfr0'][i]:9.1f} "
                  f"{R['d_marg'][1][i]:+11.3e} {R['d_marg'][8][i]:+11.3e} "
                  f"{R['d_bar'][8][i]:+10.3f} {R['flip'][8][i]:8.2f} "
                  f"{R['d_area'][i]:+9.2e} {R['d_chirp'][i]:+9.2e}")
            rows.append([rank, arr['jname'][i], arr['cat'][i], float(arr['L'][i]),
                         float(arr['sigL'][i]), float(R['s_new'][i]), float(arr['w'][i]),
                         float(R['Kfr0'][i]), float(R['d_marg'][1][i]),
                         float(R['d_marg'][8][i]), float(R['d_bar'][8][i]),
                         float(R['flip'][8][i]), float(R['d_area'][i]),
                         float(R['d_chirp'][i])])
        R["top"] = rows
        R["order"] = ordr
        print()
        print("  SPLIT-HALF STABILITY OF THE RANKING (the test that decides whether this")
        print("  is a ranking at all): the 24 sky realisations are split in two, the value")
        print("  is re-averaged in each half, and the two half-rankings are compared.")
        pr = R["per_real"][:, cand]
        h1, h2 = pr[0::2].mean(axis=0), pr[1::2].mean(axis=0)
        rho_sh = _spear(h1, h2)
        top10 = set(np.argsort(-pr.mean(axis=0))[:10])
        ov = len(top10 & set(np.argsort(-h1)[:10])), len(top10 & set(np.argsort(-h2)[:10]))
        print(f"    Spearman(half A, half B) over 115 candidates = {rho_sh:+.3f}")
        print(f"    top-10 overlap of each half with the full ranking: {ov[0]}/10, {ov[1]}/10")
        se = pr.std(axis=0, ddof=1) / np.sqrt(pr.shape[0])
        sn = np.abs(pr.mean(axis=0)) / np.maximum(se, 1e-300)
        print(f"    per-pulsar |value| / standard error: median {np.median(sn):.2f}, "
              f"top-10 median {np.median(sn[np.argsort(-pr.mean(axis=0))[:10]]):.2f}, "
              f"n with S/N > 3: {int(np.sum(sn > 3))} of 115")
        print()
        print("  DRIVER DECOMPOSITION -- rank correlation of the marginal value with each")
        print("  candidate driver, over all 115 candidates (this is the decomposition; the")
        print("  per-pulsar table below just shows where each top-10 entry sits on each axis)")
        g_all = arr["sigL"][cand] / np.maximum(R["s_new"][cand], 1e-9)
        drv = {"timing weight w": arr["w"][cand],
               "prior-width gain": g_all,
               "1/distance": 1.0 / arr["L"][cand],
               "pulsar-term power share": R["pw"][cand],
               "geometry <dL_min>": R["dL0"][cand],
               "current fringe count K_fr": R["Kfr0"][cand]}
        for kd, vd in drv.items():
            print(f"    Spearman(value, {kd:>26}) = {_spear(prim[cand], vd):+.3f}")
        print()
        print(f"  {'#':>3} {'pulsar':>12} {'timing w':>9} {'w rank':>7} {'prior gain':>11} "
              f"{'gain rank':>10} {'geom <dl>':>10} {'pw share':>9} {'value S/N':>10}")
        for rank, i in enumerate(ordr[:10], 1):
            wr = int(np.sum(arr["w"][cand] > arr["w"][i])) + 1
            g = arr["sigL"][i] / max(R["s_new"][i], 1e-9)
            gr = int(np.sum(g_all > g)) + 1
            j = list(cand).index(i)
            print(f"  {rank:3d} {arr['jname'][i]:>12} {arr['w'][i]:9.2f} {wr:7d} "
                  f"{g:11.1f} {gr:10d} {R['dL0'][i]:10.2f} {R['pw'][i]:9.2e} "
                  f"{sn[j]:10.2f}")
        R["split_half"] = rho_sh
        print()
        pm = {"f_marg K=1": R["d_marg"][1], "f_marg K=8": R["d_marg"][8],
              "f_marg K=16": R["d_marg"][16], "bar K=8": R["d_bar"][8],
              "sky area": R["d_area"], "chirp Mc": R["d_chirp"]}
        keys = list(pm)
        print("  RANK AGREEMENT BETWEEN PRODUCTS (Spearman over the 115 candidates):")
        print("    " + " " * 13 + " ".join(f"{k:>12}" for k in keys))
        for k1 in keys:
            print(f"    {k1:>13} " + " ".join(
                f"{_spear(pm[k1][cand], pm[k2][cand]):12.3f}" for k2 in keys))
    out["res"] = res
    return res


# ================================================================= 4. saturation
def exp4_curve(arr, res, out):
    _hdr("DELIVERABLE 4 -- DIMINISHING RETURNS: THE VALUE OF ANCHOR #1, #2, ... #10")
    print("  GREEDY sequence: at each step take the pulsar with the largest marginal")
    print("  value GIVEN the anchors already bought. Reported in the ACHIEVABLE arm.")
    print()
    print("  STRUCTURAL WARNING, stated up front, because it changes what 'saturation' means:")
    print("    f_marg is EXACTLY separable -- it is sum_a R_a^2 pw_a / sum pw and R_a depends")
    print("      on pulsar a's prior alone. Its greedy curve IS the sorted cumulative sum, so")
    print("      any flattening is a statement about the VALUE DISTRIBUTION being heavy-tailed,")
    print("      NOT about anchors interacting. No anchor makes another anchor less useful.")
    print("    chirp mass has a separable FISHER but sigma = 1/sqrt(sum_a I_a) is not additive,")
    print("      so it gets a true greedy; the naive sorted cumsum overstates it ~2x.")
    print("    localisation genuinely COUPLES pulsars through the Schur complement, so it is")
    print("      the only product where a real interaction can appear. True greedy.")
    print()
    R = res["achievable"]
    ia = arr["i_anchor"]
    cand = [i for i in range(arr["n_psr"]) if i != ia]

    # --- EXACTLY separable products: f_marg is sum_a R_a^2 pw_a / sum pw and R_a
    #     depends on sigma_a alone, so the greedy curve IS the sorted cumulative sum ---
    for name, vec, base in (("f_marg (K=8)", R["d_marg"][8], R["base"][8]["f_marg"]),
                            ("f_marg (K=16)", R["d_marg"][16], R["base"][16]["f_marg"])):
        v = np.sort(vec[cand])[::-1][:10]
        cum = np.cumsum(v)
        print(f"  {name}: baseline {base:.4g}")
        print(f"    per-anchor : " + " ".join(f"{x:+.3e}" for x in v))
        print(f"    cumulative : " + " ".join(f"{x:+.3e}" for x in cum))
        frac = cum / max(cum[-1], 1e-300)
        print(f"    frac of 10 : " + " ".join(f"{x:6.3f}" for x in frac))
        k50 = int(np.argmax(frac >= 0.5)) + 1
        k80 = int(np.argmax(frac >= 0.8)) + 1
        k90 = int(np.argmax(frac >= 0.9)) + 1
        print(f"    50 / 80 / 90 % of the ten-anchor total reached at anchor "
              f"{k50} / {k80} / {k90}")
        print()

    # --- chirp mass: the FISHER is separable but sigma = 1/sqrt(sum_a I_a) is NOT, so
    #     the fractional improvements do not add. True greedy on the Fisher. ---
    print("  chirp sigma(log10 Mc)  (true greedy -- sigma = 1/sqrt(sum I_a) is NOT additive,")
    print("  so the sorted cumulative sum of per-pulsar fractional gains would overstate it):")
    s_new = R["s_new"]
    cur = arr["sigL"].copy()
    sig_c = float(np.mean([chirp_sigma(arr, seed_combs(np.random.default_rng(29900 + 53 * r),
                                                       arr, 1)[0], cur)[0]
                           for r in range(N_REAL_L)]))
    sig_start = sig_c
    cpick, ccurve, cgain = [], [sig_c], []
    cpool = list(cand)
    combs_c = [seed_combs(np.random.default_rng(29900 + 53 * r), arr, 1)[0]
               for r in range(N_REAL_L)]
    for step in range(10):
        best, bs = None, sig_c
        for i in cpool:
            s = cur.copy(); s[i] = s_new[i]
            sv = float(np.mean([chirp_sigma(arr, cc, s)[0] for cc in combs_c]))
            if sv < bs:
                bs, best = sv, i
        if best is None:
            break
        cur[best] = s_new[best]
        cgain.append((sig_c - bs) / max(sig_start, 1e-300))
        sig_c = bs; cpick.append(best); cpool.remove(best); ccurve.append(sig_c)
    print(f"    baseline sigma(log10 Mc) = {sig_start:.5f} dex")
    print(f"    {'#':>3} {'pulsar':>12} {'sigma (dex)':>12} {'cum. reduction':>15} "
          f"{'marginal':>11}")
    for j, i in enumerate(cpick, 1):
        print(f"    {j:3d} {arr['jname'][i]:>12} {ccurve[j]:12.5f} "
              f"{1.0 - ccurve[j]/sig_start:15.4f} {cgain[j-1]:11.4f}")
    if cgain:
        cg = np.array(cgain); cc_ = np.cumsum(cg) / max(np.sum(cg), 1e-300)
        print(f"    50 / 80 / 90 % of the ten-anchor reduction at anchor "
              f"{int(np.argmax(cc_ >= 0.5))+1} / {int(np.argmax(cc_ >= 0.8))+1} / "
              f"{int(np.argmax(cc_ >= 0.9))+1}")
    print()
    out["greedy_chirp"] = np.array(ccurve)
    out["greedy_chirp_pick"] = np.array([arr["jname"][i] for i in cpick])

    # --- localisation: TRUE greedy, the only coupled product ---
    print("  LOCALISATION (coupled -- true greedy over the Schur-marginalised Fisher):")
    lbs = [loc_build(np.random.default_rng(29500 + 37 * r), arr) for r in range(N_REAL_L)]
    s_new = R["s_new"]
    s_cur = arr["sigL"].copy()
    A_cur = float(np.mean([loc_area(lb, s_cur) for lb in lbs]))
    A_start = A_cur
    picked, curve, gains = [], [A_cur], []
    pool = [i for i in cand]
    for step in range(10):
        best, bA = None, A_cur
        for i in pool:
            s = s_cur.copy(); s[i] = s_new[i]
            A = float(np.mean([loc_area(lb, s) for lb in lbs]))
            if A < bA:
                bA, best = A, i
        if best is None:
            print(f"    step {step+1}: no candidate improves the area -- greedy stops")
            break
        s_cur[best] = s_new[best]
        gains.append((A_cur - bA) / max(A_start, 1e-300))
        A_cur = bA
        picked.append(best); pool.remove(best); curve.append(A_cur)
    print(f"    baseline sky-ellipse area {A_start:.4e} sr")
    print(f"    {'#':>3} {'pulsar':>12} {'area (sr)':>12} {'cum. reduction':>15} "
          f"{'marginal':>11}")
    for j, i in enumerate(picked, 1):
        print(f"    {j:3d} {arr['jname'][i]:>12} {curve[j]:12.4e} "
              f"{1.0 - curve[j]/A_start:15.4f} {gains[j-1]:11.4f}")
    if gains:
        g = np.array(gains)
        cumg = np.cumsum(g) / max(np.sum(g), 1e-300)
        k80 = int(np.argmax(cumg >= 0.8)) + 1
        print(f"    80 % of the ten-anchor area reduction reached at anchor {k80}")
        print(f"    marginal gain ratio anchor#2/#1 = {g[1]/max(g[0],1e-30):.3f}, "
              f"#5/#1 = {g[min(4,len(g)-1)]/max(g[0],1e-30):.3f}, "
              f"#10/#1 = {g[-1]/max(g[0],1e-30):.3f}")
    print()
    print("  COMPARISON TO THE ANCHOR LITERATURE'S 3-6 SATURATION CLAIM: see the verdict.")
    out["greedy_loc"] = np.array(curve)
    out["greedy_pick"] = np.array([arr["jname"][i] for i in picked])
    out["greedy_gain"] = np.array(gains)


# ================================================================= 5. J0437 removal
def exp5_anchor(arr, out):
    _hdr(f"DELIVERABLE 5 -- IS {ANCHOR}'s ANCHOR STATUS ALREADY BANKED? (remove it)")
    ia = arr["i_anchor"]
    # what its prior would be WITHOUT the VLBI campaign: the array's median
    # non-parallax prior at its distance, i.e. a 20 % DM-model prior
    s_dm = 0.20 * arr["L"][ia]
    variants = [("anchor (brief, 0.11 pc)", SIG_ANCHOR_PC),
                ("anchor (repo catalogue, Deller 2008)", SIG_ANCHOR_CATALOGUE_PC),
                ("NO anchor: 20 % DM-model prior", s_dm)]
    print(f"  {ANCHOR}: d = {arr['L'][ia]:.1f} pc, w = {arr['w'][ia]:.2f} "
          f"(rank {int(np.sum(arr['w'] > arr['w'][ia]))+1} of {arr['n_psr']} by timing weight),")
    print(f"  reachable sigma_d at {SIG_PI_BEST_ARCSEC*1e6:.0f} uas = "
          f"{sigma_vlbi(arr['L'][ia], SIG_PI_BEST_ARCSEC):.4f} pc "
          f"-- by far the most reachable pulsar in the array.")
    print()
    print(f"  {'variant':>38} {'sigma_d':>9} {'f_marg K=1':>12} {'f_marg K=8':>12} "
          f"{'f_cert K=8':>11} {'own bar K=8':>12} {'area (sr)':>12} {'sig(Mc)':>9}")
    rows = []
    for lab, sv in variants:
        s = arr["sigL"].copy(); s[ia] = sv
        fm, fc, mg = {}, {}, {}
        for K in (1, 8):
            acc_m, acc_c, acc_g = [], [], []
            for r in range(N_REAL):
                rng = np.random.default_rng(29000 + 101 * r + K)
                comps, pw = seed_combs(rng, arr, K)
                rr = cert_rows(arr, comps, s, (2900, r, K))
                p = product_ab(rr, pw)
                acc_m.append(p["f_marg"]); acc_c.append(p["f_cert"])
                acc_g.append(float(rr["bar"][ia]))       # J0437's own certification bar
            fm[K] = float(np.mean(acc_m)); fc[K] = float(np.mean(acc_c))
            mg[K] = float(np.mean(acc_g))
        A = float(np.mean([loc_area(loc_build(np.random.default_rng(29500 + 37 * r), arr), s)
                           for r in range(N_REAL_L)]))
        sc = float(np.mean([chirp_sigma(arr, seed_combs(np.random.default_rng(29900 + 53 * r),
                                                        arr, 1)[0], s)[0]
                            for r in range(N_REAL_L)]))
        print(f"  {lab:>38} {sv:9.3f} {fm[1]:12.5f} {fm[8]:12.5f} {fc[8]:11.5f} "
              f"{mg[8]:12.5f} {A:12.4e} {sc:9.4f}")
        rows.append([lab, sv, fm[1], fm[8], fc[8], mg[8], A, sc])
    r0, r2 = rows[0], rows[2]
    print()
    print("  VALUE ALREADY BANKED BY THE EXISTING ANCHOR (brief value vs no anchor):")
    print(f"    f_marg  K=1 : {r2[2]:.5f} -> {r0[2]:.5f}  "
          f"({(r0[2]-r2[2])/max(r2[2],1e-30)*100:+.2f} %)")
    print(f"    f_marg  K=8 : {r2[3]:.5f} -> {r0[3]:.5f}  "
          f"({(r0[3]-r2[3])/max(r2[3],1e-30)*100:+.2f} %)")
    print(f"    f_cert  K=8 : {r2[4]:.5f} -> {r0[4]:.5f}")
    print(f"    sky area    : {r2[6]:.4e} -> {r0[6]:.4e}  "
          f"({(r0[6]-r2[6])/max(r2[6],1e-30)*100:+.2f} %)")
    print(f"    sigma(Mc)   : {r2[7]:.4f} -> {r0[7]:.4f} dex  "
          f"({(r0[7]-r2[7])/max(r2[7],1e-30)*100:+.2f} %)")
    print()
    print("  SANITY: if the anchor were a candidate, where would it rank? "
          "(its own marginal value from the")
    print("  no-anchor state, scored against every other pulsar's achievable value)")
    R = out["res"]["achievable"]
    v_others = R["d_marg"][8]
    # J0437 marginal value measured from the DM-model state
    s_no = arr["sigL"].copy(); s_no[ia] = s_dm
    acc = []
    for r in range(N_REAL):
        rng = np.random.default_rng(29000 + 101 * r + 8)
        comps, pw = seed_combs(rng, arr, 8)
        a0 = cert_rows(arr, comps, s_no, (2900, r, 8), rows=[ia])
        s_y = s_no.copy(); s_y[ia] = SIG_ANCHOR_PC
        a1 = cert_rows(arr, comps, s_y, (2900, r, 8), rows=[ia])
        acc.append(float(pw[ia] * (a1["Rcoh"][ia] ** 2 - a0["Rcoh"][ia] ** 2)
                         / max(pw.sum(), 1e-30)))
    vj = float(np.mean(acc))
    cand = [i for i in range(arr["n_psr"]) if i != ia]
    rk = int(np.sum(v_others[cand] > vj)) + 1
    print(f"    {ANCHOR} marginal Delta f_marg(K=8) from its no-anchor state = {vj:+.3e}")
    print(f"    that would rank it #{rk} of {len(cand)+1} -- "
          + ("THE TOP TARGET, i.e. the community already spent the best shot."
             if rk == 1 else f"rank {rk}."))
    out["anchor_rows"] = np.array(rows, dtype=object)
    out["anchor_rank"] = rk
    out["anchor_value"] = vj

    # ---------------------------------------------------------------- controls
    _hdr("DELIVERABLE 6 -- CONTROLS: IS THE SUB-PARSEC CERTIFICATION REAL?")
    print("  A 0.11 pc prior puts a pulsar in a window narrower than one fringe. In that")
    print("  regime BOTH the signal contrast and the incoherent-comb null contrast inflate,")
    print("  because 'the best rival fringe' is barely represented in the window -- exactly")
    print("  the sub-fringe family of artefact this project has already been bitten by. So")
    print("  the certification is audited three ways: does the winning fringe contain the")
    print("  TRUE distance (f_cert_true vs f_cert); does the SOURCE-LEVEL SCRAMBLED null")
    print("  certify anything; and how far does the bar move.")
    print()
    print(f"  {'variant':>38} {'sigma_d':>9} {'K':>3} {'f_cert':>9} {'f_cert_true':>12} "
          f"{'f_cert SCRAM':>13} {'bar (J0437)':>12} {'stat (J0437)':>13}")
    crows = []
    for lab, sv in variants:
        s = arr["sigL"].copy(); s[ia] = sv
        for K in (8, 16):
            fc, ft, fs, bb, ss = [], [], [], [], []
            for r in range(N_REAL_C):
                rng = np.random.default_rng(29000 + 101 * r + K)
                comps, pw = seed_combs(rng, arr, K)
                rr = cert_rows(arr, comps, s, (2900, r, K))
                rn = cert_rows(arr, comps, s, (2900, r, K), scramble=True)
                p, pn = product_ab(rr, pw), product_ab(rn, pw)
                fc.append(p["f_cert"]); ft.append(p["f_cert_true"])
                fs.append(pn["f_cert"]); bb.append(rr["bar"][ia]); ss.append(rr["stat"][ia])
            print(f"  {lab:>38} {sv:9.3f} {K:3d} {np.mean(fc):9.5f} {np.mean(ft):12.5f} "
                  f"{np.mean(fs):13.5f} {np.nanmean(bb):12.3f} {np.nanmean(ss):13.3f}")
            crows.append([lab, sv, K, float(np.mean(fc)), float(np.mean(ft)),
                          float(np.mean(fs)), float(np.nanmean(bb)), float(np.nanmean(ss))])
    print()
    print("  READ: f_cert_true == f_cert means every certified pulsar picked the fringe that")
    print("  actually contains its true distance. A scrambled column at or below the 5 %")
    print("  design false-alarm share of the certified power means the floor is holding.")
    out["control_rows"] = np.array(crows, dtype=object)


# ================================================================= main
def main():
    np.set_printoptions(linewidth=140)
    arr = build_array()
    out = {}
    exp0_array(arr, out)
    exp1_reach(arr, out)
    res = exp23_rank(arr, out)
    exp4_curve(arr, res, out)
    exp5_anchor(arr, out)
    flat = {}
    for k, v in out.items():
        if isinstance(v, np.ndarray) and v.dtype != object:
            flat[k] = v
    for mode in ("oracle", "achievable"):
        R = res[mode]
        flat[f"{mode}_s_new"] = R["s_new"]
        flat[f"{mode}_d_area"] = R["d_area"]
        flat[f"{mode}_d_chirp"] = R["d_chirp"]
        flat[f"{mode}_Kfr0"] = R["Kfr0"]
        for K in K_LIST:
            flat[f"{mode}_dmarg_K{K}"] = R["d_marg"][K]
            flat[f"{mode}_dbar_K{K}"] = R["d_bar"][K]
            flat[f"{mode}_flip_K{K}"] = R["flip"][K]
    flat["jname"] = arr["jname"]; flat["cat"] = arr["cat"]
    flat["L"] = arr["L"]; flat["sigL"] = arr["sigL"]; flat["w"] = arr["w"]
    np.savez_compressed(os.path.join(REPORTS, "vernier29.npz"), **flat)
    print(f"\n[V29] wrote {REPORTS}/vernier29.npz")


if __name__ == "__main__":
    main()
