#!/usr/bin/env python3
"""
VERNIER-27 chain library -- VERNIER-25's chain with four named changes and
nothing else. vernier24 / vernier25 / vernier_toy / vernier16 are IMPORTED and
never edited; the sky construction, the confusion budget, the F-statistic, the
certification statistic and its incoherent-comb floor are all called out of
those files unmodified.

The four changes, restated:

  1. SEED PACKING.  V25's chain() seeded bin b iff b < n_seed and build_bin
     prepended exactly ONE seed, so with BAND = 10 bins every K > 10 collapsed
     to ten seeds (K = 12 and K = 16 were the same cell). Here K seeds are
     distributed over the bins, several per bin where needed, and each carries
     its OWN frequency on a regular sub-bin grid spanning +-0.35/T of its bin
     centre. Frequencies stay distinct and stay inside 10-25 nHz.

     Why sub-bin frequencies are legitimate here: the comb spacing that the
     vernier works with is dL = c / (f (1 - cos mu)), a function of the source's
     frequency, and a detected source's frequency is measured to
     sigma_f = 1/(2 pi T rho_arr) -- 1.3 % of a bin at rho = 12. The toy's
     F-statistic searches sky position at a fixed bin frequency, so a sub-bin
     offset costs it no detection SNR; this is a limitation of the toy's search
     and is stated as such, not hidden. ARM W in vernier27.py re-runs the same K
     with one seed per bin in a widened band as the cross-check.

  2. NOISE / SKY SPLIT.  `noise_seed` decouples the data realisation and the
     search's randomness from the sky and array draw, so several chains can be
     run on one universe. With noise_seed=None the RNG stream is exactly V25's.

  3. COMB_CAP is a parameter, defaulting to 32 rather than V25's 8. V25's own
     caveat 1 flagged the cap as binding exactly where the mechanism starts to
     work, and this is the run that quotes K*.

  4. TOP_TRUE MADE K-INVARIANT.  A detection counts as GENUINE if it matches one
     of (the bin's seeds + the top-TOP_TRUE sky sources). V25 used the global
     top-TOP_TRUE, which with s seeds prepended leaves only TOP_TRUE - s sky
     slots, biasing the genuine count down as K rises. The chance-match floor
     rises with the enlarged pool and is reported in the controls.

Everything else -- deepening, ignition (two consecutive deepening turns), the
scrambled-null veto, the stall diagnosis, the gain/rank bookkeeping -- is V25's,
copied verbatim so the two instalments' columns mean the same thing.
"""
import numpy as np

import vernier25 as V25
import vernier24 as W
import vernier_toy as V
import vernier16 as S

YR, DF, T_OBS = W.YR, W.DF, W.T_OBS
KAPPA = W.KAPPA
SIG_TOA, NGRID, TOL_RAD, TOP_TRUE = W.SIG_TOA, W.NGRID, W.TOL_RAD, W.TOP_TRUE
R_AMP = W.R_AMP

BAND      = V25.BAND                 # 10.0-24.3 nHz, ten 1/T bins -- V25's band
N_BRIGHT  = V25.N_BRIGHT
MC_LO, MC_HI = V25.MC_LO, V25.MC_HI
RHO_SEED  = V25.RHO_SEED             # 12.0
N_NULL_C  = V25.N_NULL_C
NPTS      = V25.NPTS
N_NULL_F  = V25.N_NULL_F

COMB_CAP  = 32                       # V25 used 8; change 3
MAX_TURN  = 48                       # V25 used 16; K = 32 cannot be worked in 16
SUB_SPAN  = 0.35                     # seed sub-bin offsets span +- this x DF


# ======================================================== seed placement (1)
def seed_plan(K, band):
    """
    K seeds over the bins of `band`: ceil/floor(K/nb) per bin, each on a regular
    sub-bin frequency grid of half-width SUB_SPAN/T about its bin centre.

    Returns a list (one entry per bin) of tuples of seed frequencies.
    """
    nb = len(band)
    base, extra = divmod(int(K), nb)
    out = []
    for b, f in enumerate(band):
        n = base + (1 if b < extra else 0)
        if n == 0:
            out.append(())
        elif n == 1:
            out.append((float(f),))
        else:
            out.append(tuple(f + np.linspace(-SUB_SPAN, SUB_SPAN, n) * DF))
    return out


F_CENTRE = 17e-9                     # centre of the PTA-sensitive decade
F_FLOOR  = 3e-9                      # do not run the band below this


def wide_band(K, centre=F_CENTRE, f_floor=F_FLOOR):
    """
    ARM W: K bins of 1/T, ONE seed each -- V25 Arm K's construction continued past
    K = 10 by widening the band instead of packing the bins. Centred on the
    sensitive decade and clipped at f_floor, so it does not sit entirely in the
    quiet high-frequency bins (checked: an upward-only band from 10 nHz gives the
    same answer, because the discriminating variable is seeds-per-bin, not where
    the band sits).
    """
    K = int(K)
    f0 = max(f_floor, centre - (K - 1) / 2.0 * DF)
    return f0 + np.arange(K) * DF


# ======================================================== bin construction (1)
def build_bin(rng, arr, f, gamma, sea_keep=1.0, n_bright=N_BRIGHT,
              seed_fs=(), seed_rho=RHO_SEED, mc_lo=MC_LO, mc_hi=MC_HI):
    """
    V25.build_bin, generalised from one seed to len(seed_fs) seeds, and carrying
    a per-source frequency array `fk` (bin centre for every non-seed source).
    The sky draw itself is V25.draw_bin, untouched.
    """
    P, w, n_psr = arr["P"], arr["w"], arr["n_psr"]
    R2, S_sea, diag = V25.draw_bin(rng, f, gamma, n_bright, sea_keep)
    Swm = W.S_white(arr["sigma_toa"])
    ns = len(seed_fs)
    if ns:
        # each seed at array SNR seed_rho against white + sea, as in V25
        S_tot = Swm + S_sea
        R2_seed = seed_rho ** 2 * 2.0 * S_tot / (KAPPA * n_psr * T_OBS)
        R2 = np.concatenate([np.full(ns, R2_seed), R2])
    M = R2.size
    fk = np.concatenate([np.asarray(seed_fs, float), np.full(M - ns, float(f))])

    Gt = arr["G"][rng.choice(NGRID, M, replace=False)]
    fp, fc = W.antenna_reg(P, Gt)
    ant = np.sqrt(fp ** 2 + fc ** 2)
    ant = ant / np.sqrt(np.mean(ant ** 2)) * np.sqrt(KAPPA)

    s_k = np.sqrt(R2 * T_OBS / (2.0 * Swm))
    rho0 = ant * w[:, None] * s_k[None, :]
    C_sea = (S_sea / Swm) * w ** 2
    omega = w / np.sqrt(1.0 + C_sea)
    rho0 = rho0 / np.sqrt(1.0 + C_sea)[:, None]

    psi = rng.uniform(0, np.pi, M)
    ci = rng.uniform(-1, 1, M)
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    Fp = fp * c2[None, :] + fc * s2[None, :]
    Fc = -fp * s2[None, :] + fc * c2[None, :]
    Tm = (Fp * ((1 + ci ** 2) / 2)[None, :] + 1j * Fc * ci[None, :]) * omega[:, None]
    Tm = Tm * np.exp(1j * rng.uniform(0, 2 * np.pi, M))[None, :]
    Tn = Tm / np.maximum(np.sqrt(np.sum(np.abs(Tm) ** 2, axis=0)), 1e-30)
    ov2 = np.abs(Tn.conj().T @ Tn) ** 2

    return dict(f=f, fk=fk, R2=R2, S_sea=S_sea, diag=diag, Gt=Gt, rho0=rho0,
                C_sea=C_sea, omega=omega, Tn=Tn, ov2=ov2, M=M, G=arr["G"],
                omc=np.clip(1.0 - P @ Gt.T, 5e-2, None),
                mc=rng.uniform(mc_lo, mc_hi, M), n_seed=ns,
                Fpw=arr["fpg"] * omega[:, None], Fcw=arr["fcg"] * omega[:, None])


def _bview(B, k):
    """
    Shallow bin view whose "f" is source k's own frequency. certify() reads the
    source frequency only through B["f"] (t_merge, phase_env, dL_of), so this is
    all that is needed to give packed seeds distinct comb spacings -- and it
    leaves vernier24 untouched.
    """
    fk = float(B["fk"][k])
    if fk == B["f"]:
        return B
    Bv = dict(B)
    Bv["f"] = fk
    return Bv


# ================================================================== the chain
def chain(seed, arr, gamma, band=BAND, sea_keep=1.0, max_turn=MAX_TURN,
          null=False, soft=True, n_seed=0, seed_rho=RHO_SEED,
          comb_cap=COMB_CAP, noise_seed=None):
    """
    detect -> certify (incoherent-comb v2.2 floor) -> subtract genuine -> repeat,
    across the band, with per-turn DEEPENING bookkeeping. V25.chain verbatim
    except for the four changes documented at the top of this file.
    """
    sky_rng = np.random.default_rng(seed)
    plan = seed_plan(n_seed, band)
    bins = []
    for b, f in enumerate(band):
        bins.append(build_bin(sky_rng, arr, f, gamma, sea_keep=sea_keep,
                              seed_fs=plan[b], seed_rho=seed_rho))
    rng = sky_rng if noise_seed is None else np.random.default_rng(noise_seed)
    # change 5 (a latent V25 hazard, not a design choice): V25 drew the recut
    # detection floors from the SAME stream as the data, and floor_for is memoised
    # on (array, sea level). A cache HIT consumes no randomness and a MISS does, so
    # a V25 chain's data depended on how many chains had run before it. Floors are a
    # calibration quantity; here they get their own stream and the science stream is
    # independent of cache state. Chains are now exactly reproducible from
    # (seed, noise_seed).
    floor_rng = np.random.default_rng((int(seed) * 2654435761) & 0xFFFFFFFF)
    floors = [V25.floor_for(floor_rng, arr, B) for B in bins]
    data = [W.inject(rng, B, null=null) for B in bins]
    found = [[] for _ in bins]
    sub = [set() for _ in bins]
    det, hist = [], []
    gen = spl = blend = 0
    prevF = [None] * len(bins)
    F0 = [None] * len(bins)
    prev_cert = np.zeros(arr["n_psr"], bool)
    for t in range(max_turn):
        best, Fall = None, [None] * len(bins)
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
            Fall[b] = Fv
            j = int(np.argmax(Fv))
            exc = float(Fv[j]) - (floors[b]["floor"] + floors[b]["err"])
            if exc > 0 and (best is None or exc > best[0]):
                best = (exc, b, j, float(Fv[j]))
        if best is None:
            break
        if F0[0] is None:
            F0 = [None if x is None else x.copy() for x in Fall]
        exc, b, j, Fmax = best
        B = bins[b]
        bar_b = floors[b]["floor"] + floors[b]["err"]
        newly = bool(prevF[b] is not None and prevF[b][j] <= bar_b)
        newly0 = bool(F0[b] is not None and F0[b][j] <= bar_b)
        which, kind = -1, "blend"
        ca = B["Gt"] @ B["G"][j]
        i = int(np.argmax(ca))
        # change 4: the genuine pool is (this bin's seeds + the top-TOP_TRUE sky)
        n_true_pool = B["n_seed"] + TOP_TRUE
        if np.arccos(np.clip(ca[i], -1, 1)) <= TOL_RAD and i < n_true_pool:
            if i in sub[b]:
                kind, spl = "split", spl + 1
            else:
                kind, which, gen = "genuine", i, gen + 1
        else:
            blend += 1
        found[b].append(j)
        prevF = Fall
        if which < 0:
            hist.append(dict(turn=t + 1, bin=b, f=B["f"], F=Fmax, kind=kind,
                             which=-1, is_seed=False, newly=False, newly0=newly0,
                             deepen=False, n_combs=len(det),
                             n_cert=int(prev_cert.sum()), n_true=0, n_wrong=0,
                             n_cert_null=0, f_cert=0.0, f_soft=0.0,
                             rho_P_med=0.0, rho_P_p90=0.0, rho_arr=0.0,
                             gain_sub=hist[-1]["gain_sub"] if hist else 0.0,
                             gain_cert=hist[-1]["gain_cert"] if hist else 0.0))
            continue
        det.append((_bview(B, which), which, set(sub[b])))
        others = sorted(det[:-1], key=lambda d: -float(d[0]["R2"][d[1]]))
        use = others[:comb_cap - 1] + [det[-1]]   # target LAST: certify() weights it
        C = W.certify(rng, arr, use, soft=soft, npts=NPTS, n_null=N_NULL_C)
        new_cert = int(np.sum(C["cert"] & ~prev_cert))
        if C["n_cert"] > 0 or newly:
            Cn = W.certify(rng, arr, use, scramble=True, soft=soft, npts=NPTS,
                           n_null=N_NULL_C)
        else:
            Cn = dict(n_cert=0, f_cert=0.0)
        sub[b].add(which)
        prev_cert = prev_cert | C["cert"]
        deepen = bool((newly or new_cert > 0) and Cn["n_cert"] == 0)
        gs = []
        for Bb, sb in zip(bins, sub):
            n_all = 1.0 + (Bb["rho0"] ** 2).sum(axis=1)
            left = [q for q in range(Bb["M"]) if q not in sb]
            n_now = (1.0 + (Bb["rho0"][:, left] ** 2).sum(axis=1) if left
                     else np.ones_like(n_all))
            gs.append(np.median(np.log10(n_all / n_now)))
        g_sub = float(0.5 * np.mean(gs))
        fc = C["f_soft"] if soft else C["f_cert"]
        g_cert = float(0.5 * np.log10(1.0 + R_AMP ** 2 * fc))
        hist.append(dict(turn=t + 1, bin=b, f=B["f"], F=Fmax, kind=kind,
                         which=which, is_seed=bool(which < B["n_seed"]),
                         newly=newly, newly0=newly0, deepen=deepen,
                         new_cert=new_cert, n_combs=len(det),
                         n_cert=C["n_cert"], n_true=C["n_true"],
                         n_wrong=C["n_cert"] - C["n_true"],
                         n_cert_null=Cn["n_cert"], f_cert=C["f_cert"],
                         f_cert_true=C["f_cert_true"], f_soft=C["f_soft"],
                         rho_P_med=C["rho_P_med"], rho_P_p90=C["rho_P_p90"],
                         rho_P_max=C["rho_P_max"], rho_arr=C["rho_arr"],
                         lnK=C["lnK"], S_max=C["S_max"],
                         floor_inc=C["floor_inc"],
                         gain_sub=g_sub, gain_cert=g_cert))
    dp = [h["deepen"] for h in hist]
    ignite = any(dp[i] and dp[i + 1] for i in range(len(dp) - 1))
    conf = []
    for Bb, sb in zip(bins, sub):
        left = [q for q in range(Bb["M"]) if q not in sb]
        c = (Bb["rho0"][:, left] ** 2).sum(axis=1) if left else np.zeros(arr["n_psr"])
        conf.append(np.median(c))
    conf_frac = float(np.mean(conf) / (np.mean(conf) + 1.0))
    g_tot = (hist[-1]["gain_sub"] + hist[-1]["gain_cert"]) if hist else 0.0
    ncert = sum(h.get("new_cert", 0) for h in hist)
    # A_equivalent label: the seed's loudness against this sky's loudest source
    B0 = bins[0]
    if B0["n_seed"]:
        a_eq = float(np.sqrt(B0["R2"][0] / B0["R2"][B0["n_seed"]]))
    else:
        a_eq = 1.0
    return dict(hist=hist, genuine=gen, split=spl, blend=blend, turns=len(hist),
                gain=g_tot, rank=W.rank_of(g_tot), ignite=ignite,
                deepen=int(sum(dp)), n_cert_tot=ncert,
                n_wrong_tot=sum(h.get("n_wrong", 0) for h in hist),
                n_true_tot=sum(h.get("n_true", 0) for h in hist),
                seed_hits=sum(1 for h in hist if h.get("is_seed")),
                exposed=sum(1 for h in hist if h.get("newly0")),
                beyond=sum(1 for h in hist
                           if h["which"] >= 0 and not h.get("is_seed")),
                conf_frac=conf_frac, a_eq=a_eq, bins=bins, floors=floors)


def stall_reason(res):
    """V25's stall diagnosis, unchanged."""
    return V25.stall_reason(res)
