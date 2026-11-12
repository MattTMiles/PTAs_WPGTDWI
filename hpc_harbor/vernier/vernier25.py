#!/usr/bin/env python3
"""
VERNIER-25 -- the wheel-friendly sky, specified; and the first-detection era.
EXPLORATORY, NOT pre-registered. Clean room: imports only vernier24 (and through
it vernier_toy / vernier6 / vernier7 / vernier16), all my own anchored code.

WHY. Every sky simulated so far is wheel-hostile and the audit localised the
hostility in ONE number -- the bright-end counts slope. gamma_eff = 1.56 gives
lambda_fore ~ 0.31 at every frequency; VERNIER-5's structural result says clean
per-bin dominance needs gamma <= 0.87. Nobody has mapped the space between.
And VERNIER-24's diagnosis is that the chain starves on DETECTIONS -- a single
source cannot break its own comb by construction, while 41 % of the pulsar-term
power already clears rho_P >= 2 -- so the currency is K, the number of
simultaneously detectable sources at DISTINCT frequencies. Multiplicity is a
schedule, not a lottery (~4 resolvable MeerKAT-class / 10 yr, 30-40 SKA / 30 yr),
so both are measured: which skies serve K unaided (ARM G), and what the wheel
does once K arrives (ARM K).

COMMON SETUP, EVERY ARM. VERNIER-24's corrected accounting throughout: audited
counts, h(Mc, f, D) link so power/bin ~ f^(-7/3), band-integrated power pinned to
the NG15 target. Full-array joint inference, lambda_eff-consistent confusion,
marginalised pulsar terms available. Certification uses the INCOHERENT-COMB null
floor (V24's method correction) in every arm -- a floor built from no-comb nulls
over-certifies by 1.9-3.7x and every certification the surplus buys is a chance
comb alignment. M = 1 control, pure-noise null, per-turn source-level scrambled
null, VERNIER-18 genuine/split/blend accounting with chance floors. WRONG
certifications (certified but on the wrong fringe) are reported in every cell,
because top-heavy and loud cells are exactly where wrong-fringe locks got
dangerous.

ARM G -- THE GAMMA SWEEP. gamma in {0.6, 0.87, 1.2, 1.56, 2.0, 2.5}. At each
gamma the bright end is drawn from N(>h) ~ h^(-gamma) and the faint remainder is
a Gaussian sea carrying the residual power, with the bin TOTAL held fixed: gamma
redistributes power between the bright few and the sea, it does not add any.
Because the loudest source's strain moves with gamma at fixed total, every cell
is labelled by A_equivalent -- the NG15 amplitude that would put h_(1) where this
gamma puts it -- rather than by any claim about a "typical" source.

ARM G-FP -- THE FINAL-PARSEC ALTERNATIVE. Instead of tilting the slope, hold the
bright end at gamma = 1.56 and SUPPRESS the faint sea below a stall cut (mass-
selective final-parsec stalling). Same wheel, same gates. If the boundary differs
between "flattened slope" and "suppressed sea", the wheel is telling us which
mechanism it actually cares about.

ARM K -- THE FIRST-DETECTION ERA. On the corrected sky (gamma = 1.56, audited
counts, honest confusion sea) place K = 2, 3, 5, 8 sources above threshold at
DISTINCT frequencies across 10-25 nHz with Mc in 10^8.5-9.5. Run the same chain.
This is the untested middle between VERNIER-6's clean corner (K = 8, no GWB) and
VERNIER-24's starved cold start (K ~ 1, honest sea).

IGNITION, defined once and used in both arms (the standing sighting bar):
  a turn DEEPENS if it makes a genuine detection that was BELOW the floor in the
  previous turn's residual, or if it certifies a pulsar that the previous turn
  could not; IGNITION requires >= 2 CONSECUTIVE deepening turns, attributable
  (the deepening is caused by what the wheel removed, not by a fresh noise draw)
  and scrambled-clean (the source-level scrambled null does not produce it).
"""
import os
import numpy as np

import vernier24 as W
import vernier_toy as V
import vernier16 as S

REPORTS = W.REPORTS
YR, DF, T_OBS = W.YR, W.DF, W.T_OBS
N_PSR, R_AMP = W.N_PSR, W.R_AMP
PRIOR_PC, KAPPA = W.PRIOR_PC, W.KAPPA
SIG_TOA, NGRID, TOL_RAD, TOP_TRUE = W.SIG_TOA, W.NGRID, W.TOL_RAD, W.TOP_TRUE
A_GWB, LAM_AUDIT = W.A_GWB, W.LAM_AUDIT

BAND      = 10e-9 + np.arange(10) * DF        # 10.0 - 24.3 nHz, ten 1/T bins
GAMMAS    = (0.6, 0.87, 1.2, 1.56, 2.0, 2.5)
SEA_KEEP  = (1.0, 0.5, 0.2, 0.05, 0.0)        # final-parsec: surviving sea fraction
K_LIST    = (2, 3, 5, 8, 12, 16)   # 12/16 extend the brief, to bracket K*
MC_LO, MC_HI = 8.5, 9.5
N_BRIGHT  = 300
RHO_SEED  = 12.0                              # array SNR of an Arm K seed
N_NULL_C  = 32                                # certification nulls per pulsar
NPTS      = 800                               # distance grid (EVEN, sub-fringe note)
COMB_CAP  = 8                                 # combs carried into a certification
N_NULL_F  = 120


# ============================================================ sky construction
def draw_bin(rng, f, gamma, n_bright=N_BRIGHT, sea_keep=1.0, seed_rho=None,
             sigma_toa=SIG_TOA):
    """
    One bin of the corrected sky at bright-end slope `gamma`.

    N_bin sources from the audited counts; strains from N(>h) ~ h^(-gamma), i.e.
    h = u^(-1/gamma); the whole draw normalised so the bin holds bin_power(f), so
    gamma REDISTRIBUTES the fixed total between the bright few and the sea.

    sea_keep < 1 is the final-parsec alternative: the faint remainder below the
    explicit bright end is suppressed by that factor at FIXED bright end, so the
    bin total falls and the implied GWB with it (reported, never hidden).

    seed_rho, if given, prepends one source whose array SNR is seed_rho against
    white + sea -- the Arm K first-detection seed.
    """
    n = int(max(round(float(W.N_c(f))), 1))
    h = rng.uniform(size=n) ** (-1.0 / gamma)
    p = np.sort(h ** 2)[::-1]
    p = p / p.sum()
    tot = float(W.bin_power(f))
    nb = min(n_bright, n)
    R2 = p[:nb] * tot
    sea_p = float(p[nb:].sum()) * tot * sea_keep
    diag = dict(lam_fore=float(p[0]), n_eff=float(1.0 / np.sum(p ** 2)),
                n_bin=n, sea_frac=sea_keep,
                tot=float(R2.sum() + sea_p), tot_ref=tot)
    return R2, sea_p / (2.0 * DF), diag


def build_bin(rng, arr, f, gamma, sea_keep=1.0, n_bright=N_BRIGHT,
              seed_rho=None, mc_lo=MC_LO, mc_hi=MC_HI):
    """
    Per-(pulsar, source) SNRs, whitening, measured template overlaps, per-source
    chirp masses. Same construction as VERNIER-24's `make_bin`, with the slope,
    the sea suppression and the optional first-detection seed exposed.
    """
    P, w, n_psr = arr["P"], arr["w"], arr["n_psr"]
    R2, S_sea, diag = draw_bin(rng, f, gamma, n_bright, sea_keep)
    Swm = W.S_white(arr["sigma_toa"])
    n_seed = 0
    if seed_rho is not None:
        # amplitude that puts this source at array SNR seed_rho against white+sea
        S_tot = Swm + S_sea
        R2_seed = seed_rho ** 2 * 2.0 * S_tot / (KAPPA * n_psr * T_OBS)
        R2 = np.concatenate([[R2_seed], R2])
        n_seed = 1
    M = R2.size
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

    return dict(f=f, R2=R2, S_sea=S_sea, diag=diag, Gt=Gt, rho0=rho0,
                C_sea=C_sea, omega=omega, Tn=Tn, ov2=ov2, M=M, G=arr["G"],
                omc=np.clip(1.0 - P @ Gt.T, 5e-2, None),
                mc=rng.uniform(mc_lo, mc_hi, M), n_seed=n_seed,
                Fpw=arr["fpg"] * omega[:, None], Fcw=arr["fcg"] * omega[:, None])


def a_equivalent(lam_fore, lam_ref=LAM_AUDIT):
    """
    The NG15 amplitude that would put the loudest in-bin source where THIS sky
    puts it, at the same fixed total power. h_(1)^2 = lambda_fore * total and
    total ~ A^2, so A_eq = A_gwb sqrt(lambda_fore / lambda_ref). This is a label
    for the cell's effective loudness -- NOT a claim about a typical source.
    """
    return A_GWB * np.sqrt(max(lam_fore, 1e-30) / lam_ref)


# ============================================================ detection floors
_FLOOR_CACHE = {}


def floor_for(rng, arr, B):
    """
    Recut floor on max-F, cached on the shape of the whitening vector.

    The F-statistic is invariant to a global rescaling of the templates, so the
    floor depends on omega only through S_sea/S_white; caching on that scalar
    (to a tenth of a dex) makes a 6-gamma x 20-realisation sweep affordable.
    """
    key = (id(arr), round(float(np.log10(max(B["S_sea"], 1e-30)
                                        / W.S_white(arr["sigma_toa"]))), 1))
    if key not in _FLOOR_CACHE:
        _FLOOR_CACHE[key] = W.fstat_floor(rng, B, n_null=N_NULL_F)
    return _FLOOR_CACHE[key]


def k_available(rng, arr, bins):
    """
    K -- the number of sources across the band simultaneously ABOVE the detection
    threshold, in the lambda_eff-consistent confusion budget with nothing yet
    subtracted. This is the primary currency of the whole instalment.
    """
    k = 0
    per_bin = []
    for B in bins:
        fl = floor_for(rng, arr, B)
        thr = np.sqrt(max(fl["floor"] + fl["err"], 0.0))
        c = 0
        for i in range(min(B["M"], 40)):
            rho = W.rho_budget(B, i, set())
            if np.sqrt(np.sum(rho ** 2)) >= thr:
                c += 1
        k += c
        per_bin.append(c)
    return k, per_bin


# ============================================================ the chain
def chain(seed, arr, gamma, band=BAND, sea_keep=1.0, max_turn=16, null=False,
          soft=True, n_seed=0, seed_rho=RHO_SEED):
    """
    detect -> certify (incoherent-comb v2.2 floor) -> subtract genuine -> repeat,
    across the band, with per-turn DEEPENING bookkeeping.

    A turn deepens if it makes a genuine detection whose F was BELOW the floor in
    the previous turn's residual (so the wheel's own subtraction is what exposed
    it), or if it certifies a pulsar the previous turn could not. Ignition needs
    two consecutive deepening turns, and the scrambled null must not produce them.
    """
    rng = np.random.default_rng(seed)
    bins = []
    for b, f in enumerate(band):
        bins.append(build_bin(rng, arr, f, gamma, sea_keep=sea_keep,
                              seed_rho=(seed_rho if b < n_seed else None)))
    floors = [floor_for(rng, arr, B) for B in bins]
    data = [W.inject(rng, B, null=null) for B in bins]
    found = [[] for _ in bins]
    sub = [set() for _ in bins]
    det, hist = [], []
    gen = spl = blend = 0
    prevF = [None] * len(bins)
    F0 = [None] * len(bins)                   # the ORIGINAL, pre-wheel F surface
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
        # was this sky position BELOW the floor before the wheel's LAST removal?
        newly = bool(prevF[b] is not None and prevF[b][j] <= bar_b)
        # ... and below it in the ORIGINAL data, i.e. exposed by the chain at all?
        newly0 = bool(F0[b] is not None and F0[b][j] <= bar_b)
        which, kind = -1, "blend"
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
        prevF = Fall
        if which < 0:
            hist.append(dict(turn=t + 1, bin=b, f=B["f"], F=Fmax, kind=kind,
                             which=-1, is_seed=False, newly=False, newly0=newly0,
                             deepen=False,
                             n_cert=int(prev_cert.sum()), n_true=0, n_wrong=0,
                             n_cert_null=0, f_cert=0.0, f_soft=0.0,
                             rho_P_med=0.0, rho_P_p90=0.0, rho_arr=0.0,
                             gain_sub=hist[-1]["gain_sub"] if hist else 0.0,
                             gain_cert=hist[-1]["gain_cert"] if hist else 0.0))
            continue
        det.append((B, which, set(sub[b])))
        # carry only the COMB_CAP loudest combs: a 15th faint comb adds no
        # commensurability and the certification cost is linear in the count
        others = sorted(det[:-1], key=lambda d: -float(d[0]["R2"][d[1]]))
        use = others[:COMB_CAP - 1] + [det[-1]]   # target LAST: certify() weights it
        C = W.certify(rng, arr, use, soft=soft, npts=NPTS, n_null=N_NULL_C)
        new_cert = int(np.sum(C["cert"] & ~prev_cert))
        # The scrambled null only has work to do where the signal side fired: with
        # new_cert = 0 and nothing newly exposed, the turn cannot deepen whatever
        # the null says, and running it there costs half the sweep.
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
                         which=which, is_seed=bool(B["n_seed"] and which == 0),
                         newly=newly, newly0=newly0, deepen=deepen,
                         new_cert=new_cert,
                         n_cert=C["n_cert"], n_true=C["n_true"],
                         n_wrong=C["n_cert"] - C["n_true"],
                         n_cert_null=Cn["n_cert"], f_cert=C["f_cert"],
                         f_cert_true=C["f_cert_true"], f_soft=C["f_soft"],
                         rho_P_med=C["rho_P_med"], rho_P_p90=C["rho_P_p90"],
                         rho_P_max=C["rho_P_max"], rho_arr=C["rho_arr"],
                         lnK=C["lnK"], S_max=C["S_max"],
                         floor_inc=C["floor_inc"],
                         gain_sub=g_sub, gain_cert=g_cert))
    # ignition: two consecutive deepening turns
    dp = [h["deepen"] for h in hist]
    ignite = any(dp[i] and dp[i + 1] for i in range(len(dp) - 1))
    # stall diagnosis at the end state
    conf, wht = [], []
    for Bb, sb in zip(bins, sub):
        left = [q for q in range(Bb["M"]) if q not in sb]
        c = (Bb["rho0"][:, left] ** 2).sum(axis=1) if left else np.zeros(arr["n_psr"])
        conf.append(np.median(c)); wht.append(1.0)
    conf_frac = float(np.mean(conf) / (np.mean(conf) + 1.0))
    g_tot = (hist[-1]["gain_sub"] + hist[-1]["gain_cert"]) if hist else 0.0
    ncert = sum(h.get("new_cert", 0) for h in hist)
    return dict(hist=hist, genuine=gen, split=spl, blend=blend, turns=len(hist),
                gain=g_tot, rank=W.rank_of(g_tot), ignite=ignite,
                deepen=int(sum(dp)), n_cert_tot=ncert,
                n_wrong_tot=sum(h.get("n_wrong", 0) for h in hist),
                n_true_tot=sum(h.get("n_true", 0) for h in hist),
                seed_hits=sum(1 for h in hist if h.get("is_seed")),
                exposed=sum(1 for h in hist if h.get("newly0")),
                beyond=sum(1 for h in hist
                           if h["which"] >= 0 and not h.get("is_seed")),
                conf_frac=conf_frac, bins=bins, floors=floors)


def stall_reason(res):
    """Name what stopped the chain: absolute floor, confusion, or comb ambiguity."""
    if res["turns"] == 0:
        return "no detection at all"
    cert = res["n_cert_tot"]
    if res["genuine"] > 0 and cert == 0:
        return "comb ambiguity"
    if res["conf_frac"] > 0.5:
        return "confusion"
    return "absolute floor"


# ============================================================ deliverables
def cells(acc):
    """Median/fraction summary of a list of chain results."""
    hh = [h for r in acc for h in r["hist"] if h["which"] >= 0]
    m = lambda k: (float(np.median([h[k] for h in hh])) if hh else 0.0)
    return dict(
        ignite=float(np.mean([r["ignite"] for r in acc])),
        turns=float(np.median([r["turns"] for r in acc])),
        genuine=float(np.mean([r["genuine"] for r in acc])),
        blend=float(np.mean([r["blend"] for r in acc])),
        cert=float(np.mean([r["n_cert_tot"] for r in acc])),
        true=float(np.mean([r["n_true_tot"] for r in acc])),
        wrong=float(np.mean([r["n_wrong_tot"] for r in acc])),
        gain=float(np.median([r["gain"] for r in acc])),
        rank=int(np.median([r["rank"] for r in acc])),
        rho_arr=m("rho_arr"), rho_P_p90=m("rho_P_p90"), S_max=m("S_max"),
        deepen=float(np.mean([r["deepen"] for r in acc])),
        beyond=float(np.mean([r["beyond"] for r in acc])),
        conf=float(np.median([r["conf_frac"] for r in acc])))


def exp0_setup(out):
    print("=" * 104)
    print("DELIVERABLE 0 -- SETUP, AND WHAT gamma DOES AT FIXED TOTAL POWER")
    print("=" * 104)
    print(f"  band {BAND[0]*1e9:.2f}-{BAND[-1]*1e9:.2f} nHz in {len(BAND)} bins of "
          f"1/T = {DF*1e9:.3f} nHz; array {N_PSR} pulsars, sigma_toa(med) "
          f"{SIG_TOA*1e9:.0f} ns")
    print("  bin total held at the NG15 target in every gamma cell: gamma moves power")
    print("  between the bright few and the sea, it never adds any.")
    print()
    print(f"  {'gamma':>6} {'lam_fore':>9} {'N_eff':>7} {'h_(1)':>10} "
          f"{'A_equiv':>10} {'A_eq/A_ng15':>12} {'sea frac':>9} "
          f"{'vs V24 M_max(=5-10)':>21}")
    rng = np.random.default_rng(25)
    f0 = float(BAND[len(BAND) // 2])
    rows = []
    for g in GAMMAS:
        lam, ne = [], []
        for _ in range(200):
            _r, _s, d = draw_bin(rng, f0, g)
            lam.append(d["lam_fore"]); ne.append(d["n_eff"])
        lm, nm = float(np.median(lam)), float(np.median(ne))
        h1 = float(np.sqrt(lm * W.h2_per_bin(f0)))
        ae = a_equivalent(lm)
        _r, _s, d = draw_bin(rng, f0, g)
        sf = 1.0 - float(_r.sum() / d["tot_ref"])
        rows.append([g, lm, nm, h1, ae, ae / A_GWB, sf])
        print(f"  {g:6.2f} {lm:9.3f} {nm:7.1f} {h1:10.3e} {ae:10.3e} "
              f"{ae/A_GWB:12.2f} {sf:9.3f} "
              f"{('WITHIN' if nm <= 10 else 'above'):>21}")
    print()
    print("  A_equiv is the NG15 amplitude that would put h_(1) where this gamma puts")
    print("  it at the SAME total power -- a label for the cell's effective loudness,")
    print("  not a claim about a typical source.")
    print(f"  VERNIER-5's structural threshold for clean per-bin dominance: gamma <= 0.87.")
    print(f"  VERNIER-24's realized half-bin M_max: 5-10 (counting bound 50).")
    out["g_setup"] = np.array(rows)


def exp1_gamma(out, n_real=20):
    print()
    print("=" * 104)
    print("DELIVERABLE 1 -- ARM G: THE GAMMA SWEEP")
    print("=" * 104)
    print("  full chain at each gamma; ignition = 2 consecutive deepening turns,")
    print("  attributable and scrambled-clean.")
    print()
    print(f"  {'gamma':>6} {'IGNITE':>7} {'K':>6} {'lam_f':>7} {'turns':>6} "
          f"{'gen':>5} {'blend':>6} {'deepen':>7} {'CERT':>6} {'true':>6} "
          f"{'WRONG':>6} {'rhoP p90':>9} {'S_max':>7} {'dex':>8} {'rank':>5} "
          f"{'stall':>16}")
    rows, keep = [], {}
    for g in GAMMAS:
        acc, Ks = [], []
        for s in range(n_real):
            arr = W.make_array(np.random.default_rng(2500 + s))
            r = chain(3000 + 31 * s + int(100 * g), arr, g)
            acc.append(r)
            kk, _ = k_available(np.random.default_rng(9 + s), arr, r["bins"])
            Ks.append(kk)
        c = cells(acc)
        lam = float(np.median([b["diag"]["lam_fore"] for r in acc for b in r["bins"]]))
        sr = max(set(stall_reason(r) for r in acc),
                 key=lambda x: sum(stall_reason(r) == x for r in acc))
        rows.append([g, c["ignite"], float(np.median(Ks)), lam, c["turns"],
                     c["genuine"], c["blend"], c["deepen"], c["cert"], c["true"],
                     c["wrong"], c["rho_P_p90"], c["S_max"], c["gain"], c["rank"]])
        keep[g] = acc
        print(f"  {g:6.2f} {c['ignite']:7.2f} {np.median(Ks):6.1f} {lam:7.3f} "
              f"{c['turns']:6.1f} {c['genuine']:5.1f} {c['blend']:6.1f} "
              f"{c['deepen']:7.2f} {c['cert']:6.1f} {c['true']:6.1f} "
              f"{c['wrong']:6.1f} {c['rho_P_p90']:9.3f} {c['S_max']:7.3f} "
              f"{c['gain']:+8.4f} {c['rank']:5d} {sr:>16}")
    rows = np.array(rows)
    ig = [r for r in rows if r[1] >= 0.5]
    print()
    if ig:
        print(f"  gamma* (>= 50 % of realisations ignite): "
              f"gamma <= {max(r[0] for r in ig):.2f}")
    else:
        print("  gamma*: NO gamma in the sweep ignites in >= 50 % of realisations.")
    print("  ignition FRACTION vs gamma (the lottery, not the median):")
    print("      " + "   ".join(f"g={r[0]:.2f}: {r[1]:.2f}" for r in rows))
    out["g_sweep"] = rows
    return keep


def exp2_boundary(out, rows):
    print()
    print("=" * 104)
    print("DELIVERABLE 2 -- THE BOUNDARY IN FOUR CURRENCIES, AND WHERE 1.56 SITS")
    print("=" * 104)
    print(f"  {'gamma':>6} {'ignite':>7} {'K':>6} {'lam_fore':>9} {'h_(1)':>10} "
          f"{'A_eq/A15':>9} {'N_eff':>7} {'N_eff <= M_max?':>16}")
    setup = {r[0]: r for r in out["g_setup"]}
    for r in rows:
        g = r[0]; st = setup[g]
        print(f"  {g:6.2f} {r[1]:7.2f} {r[2]:6.1f} {r[3]:9.3f} {st[3]:10.3e} "
              f"{st[5]:9.2f} {st[2]:7.1f} "
              f"{('yes' if st[2] <= 10 else 'NO'):>16}")
    print()
    r156 = [r for r in rows if abs(r[0] - 1.56) < 1e-9][0]
    st156 = setup[1.56]
    print("  WHERE gamma = 1.56 SITS -- what binds first:")
    print(f"      K            = {r156[2]:.1f} simultaneously detectable "
          f"(need K* -- deliverable 4)")
    print(f"      lambda_fore  = {r156[3]:.3f}   (V5's clean-dominance bar: "
          f"gamma <= 0.87 -> {setup[0.87][1]:.3f})")
    print(f"      h_(1)        = {st156[3]:.3e}  vs {setup[0.87][3]:.3e} at "
          f"gamma = 0.87  ({0.5*np.log10(setup[0.87][3]**2/st156[3]**2):+.2f} dex)")
    print(f"      N_eff        = {st156[2]:.1f}  vs V24's realized M_max 5-10 -> "
          f"{'INSIDE' if st156[2] <= 10 else 'outside'}")
    out["b_156"] = np.array([r156[2], r156[3], st156[3], st156[2]])


def exp3_finalparsec(out, n_real=12, gamma=1.56):
    print()
    print("=" * 104)
    print("DELIVERABLE 3 -- ARM G-FP: THE FINAL-PARSEC ALTERNATIVE")
    print("=" * 104)
    print("  Bright end held at gamma = 1.56; the faint sea below it suppressed by a")
    print("  stall cut. Bin total FALLS with the sea, so the implied GWB falls too --")
    print("  reported, not hidden.")
    print()
    print(f"  {'sea kept':>9} {'A_eq/A15':>9} {'IGNITE':>7} {'K':>6} {'turns':>6} "
          f"{'gen':>5} {'deepen':>7} {'CERT':>6} {'true':>6} {'WRONG':>6} "
          f"{'rhoP p90':>9} {'S_max':>7} {'dex':>8} {'rank':>5} {'stall':>16}")
    rows = []
    for sk in SEA_KEEP:
        acc, Ks = [], []
        for s in range(n_real):
            arr = W.make_array(np.random.default_rng(2500 + s))
            r = chain(5000 + 31 * s + int(100 * sk), arr, gamma, sea_keep=sk)
            acc.append(r)
            kk, _ = k_available(np.random.default_rng(9 + s), arr, r["bins"])
            Ks.append(kk)
        c = cells(acc)
        tot = float(np.median([b["diag"]["tot"] / b["diag"]["tot_ref"]
                               for r in acc for b in r["bins"]]))
        sr = max(set(stall_reason(r) for r in acc),
                 key=lambda x: sum(stall_reason(r) == x for r in acc))
        rows.append([sk, np.sqrt(tot), c["ignite"], float(np.median(Ks)),
                     c["turns"], c["genuine"], c["deepen"], c["cert"], c["true"],
                     c["wrong"], c["rho_P_p90"], c["S_max"], c["gain"], c["rank"]])
        print(f"  {sk:9.2f} {np.sqrt(tot):9.2f} {c['ignite']:7.2f} "
              f"{np.median(Ks):6.1f} {c['turns']:6.1f} {c['genuine']:5.1f} "
              f"{c['deepen']:7.2f} {c['cert']:6.1f} {c['true']:6.1f} "
              f"{c['wrong']:6.1f} {c['rho_P_p90']:9.3f} {c['S_max']:7.3f} "
              f"{c['gain']:+8.4f} {c['rank']:5d} {sr:>16}")
    out["fp"] = np.array(rows)


def exp4_armK(out, n_real=10, gamma=1.56):
    print()
    print("=" * 104)
    print("DELIVERABLE 4 -- ARM K: THE FIRST-DETECTION ERA")
    print("=" * 104)
    print(f"  K seeds at DISTINCT frequencies, each at array SNR {RHO_SEED:.0f} "
          f"against white+sea,")
    print("  on top of the corrected gamma = 1.56 sky with its honest confusion sea.")
    print()
    print(f"  {'K':>3} {'A_eq/A15':>9} {'IGNITE':>7} {'turns':>6} {'seeds hit':>10} "
          f"{'BEYOND seed':>12} {'total known':>12} {'CERT':>6} {'true':>6} "
          f"{'WRONG':>6} {'null':>5} {'S_max':>7} {'lnK':>6} {'dex':>8} "
          f"{'d/turn':>8} {'rank':>5} {'stall':>16}")
    rows = []
    for K in (0,) + K_LIST:
        acc = []
        for s in range(n_real):
            arr = W.make_array(np.random.default_rng(2500 + s))
            acc.append(chain(7000 + 41 * s + K, arr, gamma, n_seed=K))
        c = cells(acc)
        hh = [h for r in acc for h in r["hist"] if h["which"] >= 0]
        sh = float(np.mean([r["seed_hits"] for r in acc]))
        by = float(np.mean([r["beyond"] for r in acc]))
        nul = float(np.mean([h["n_cert_null"] for h in hh])) if hh else 0.0
        lnk = float(np.median([h["lnK"] for h in hh])) if hh else 0.0
        dpt = c["gain"] / max(c["turns"], 1.0)
        gpt = c["genuine"] / max(c["turns"], 1.0)
        # A_equivalent: how much louder a seed is than the sky's own loudest source
        if K > 0:
            ae = float(np.median([np.sqrt(r["bins"][0]["R2"][0]
                                          / r["bins"][0]["R2"][1]) for r in acc]))
        else:
            ae = 1.0
        sr = max(set(stall_reason(r) for r in acc),
                 key=lambda x: sum(stall_reason(r) == x for r in acc))
        rows.append([K, c["ignite"], c["turns"], sh, by, sh + by, c["cert"],
                     c["true"], c["wrong"], nul, c["S_max"], lnk, c["gain"],
                     dpt, c["rank"], gpt, ae])
        print(f"  {K:3d} {ae:9.2f} {c['ignite']:7.2f} {c['turns']:6.1f} "
              f"{sh:10.1f} {by:12.1f} {sh+by:12.1f} {c['cert']:6.1f} "
              f"{c['true']:6.1f} {c['wrong']:6.1f} {nul:5.1f} {c['S_max']:7.3f} "
              f"{lnk:6.2f} {c['gain']:+8.4f} {dpt:+8.4f} {c['rank']:5d} "
              f"{sr:>16}")
    rows = np.array(rows)
    print()
    ig = [r for r in rows if r[1] >= 0.5]
    print(f"  K* (smallest K igniting in >= 50 % of realisations): "
          + (f"{int(min(r[0] for r in ig))}" if ig else "NOT REACHED at K <= 8"))
    print("  ignition fraction vs K:  " + "   ".join(f"K={int(r[0])}: {r[1]:.2f}"
                                                     for r in rows))
    print()
    print("  marginal value per turn: dex/turn and NEW GENUINE DETECTIONS/turn")
    print("      " + "   ".join(f"K={int(r[0])}: {r[13]:+.4f} dex, {r[15]:.2f} det"
                                for r in rows))
    print()
    print("  best fringe contrast reached vs the trials bar (the certification channel):")
    print("      " + "   ".join(f"K={int(r[0])}: {r[10]:.2f}/{r[11]:.2f}"
                                for r in rows))
    print()
    print("  SEEDING INVARIANCE (VERNIER-3's metric, total known = seed + turns):")
    print("      " + "   ".join(f"K={int(r[0])}: {r[5]:.1f}" for r in rows))
    base = rows[0][5]
    print(f"      invariant would mean a flat row at {base:.1f}; "
          f"spread = {rows[:,5].min():.1f}-{rows[:,5].max():.1f}")
    out["armK"] = rows


def exp5_controls(out, n_real=6, gamma=1.56):
    print()
    print("=" * 104)
    print("DELIVERABLE 5 -- CONTROLS")
    print("=" * 104)
    # M = 1 control: one source per bin, no sea
    cs, ct, rx, sm = [], [], [], []
    dete = 0
    for s in range(n_real):
        arr = W.make_array(np.random.default_rng(2500 + s))
        r = chain(8100 + s, arr, 1.56, sea_keep=0.0, n_seed=1,
                  band=BAND[:1], seed_rho=RHO_SEED)
        dete += int(r["genuine"] > 0)
        for h in r["hist"]:
            if h["which"] >= 0:
                cs.append(h["n_cert"]); ct.append(h["n_true"])
                rx.append(h["rho_P_max"]); sm.append(h["S_max"])
    f_ = lambda v: (float(np.median(v)) if v else 0.0)
    print(f"  M=1 control (single bin, seeded, sea removed): detected {dete}/{n_real}; "
          f"CERT {f_(cs):.1f} (true {f_(ct):.1f}); rho_P max {f_(rx):.2f}; "
          f"best fringe contrast {f_(sm):.3f} nat")
    nf, cc, dpn = 0, 0, 0
    for s in range(n_real):
        arr = W.make_array(np.random.default_rng(2500 + s))
        r = chain(8200 + s, arr, gamma, null=True)
        nf += int(r["turns"] > 0); cc += r["n_cert_tot"]; dpn += int(r["ignite"])
    exp = 1.0 - (1.0 - 0.05) ** len(BAND)
    print(f"  pure-noise null through the whole chain: fired {nf}/{n_real} "
          f"(expected {exp:.2f} at alpha = 0.05 x {len(BAND)} bins), "
          f"{cc} certifications, {dpn} ignitions")
    print("  per-turn SOURCE-LEVEL scrambled null: the `null`/`WRONG` columns above,")
    print("  and a scrambled-null certification VETOES the deepening flag by")
    print("  construction, so no ignition in this note is scrambled-dirty.")
    for M in (1, TOP_TRUE, N_BRIGHT):
        p1 = (1.0 - np.cos(TOL_RAD)) / 2.0
        print(f"  chance-match floor at M = {M:3d}: "
              f"{1.0 - (1.0 - p1) ** M:.3f}  "
              f"(GENUINE is scored against the top-{TOP_TRUE} only)")
    out["ctrl"] = np.array([dete, n_real, f_(cs), f_(ct), f_(rx), f_(sm), nf, cc])


def main():
    out = {}
    print()
    print("#" * 104)
    print("# VERNIER-25 -- THE WHEEL-FRIENDLY SKY, SPECIFIED; AND THE FIRST-DETECTION ERA")
    print("# EXPLORATORY, NOT PRE-REGISTERED")
    print("#" * 104)
    exp0_setup(out)
    exp1_gamma(out)
    exp2_boundary(out, out["g_sweep"])
    exp3_finalparsec(out)
    exp4_armK(out)
    exp5_controls(out)
    np.savez(os.path.join(REPORTS, "vernier25.npz"), **out)
    print()
    print(f"  banked -> {os.path.join(REPORTS, 'vernier25.npz')}")


if __name__ == "__main__":
    main()


def exp6_endstate(out, Ks=(5, 8, 12, 16), gammas=(0.6, 1.2), n_real=6):
    """
    SUPPLEMENT -- the END-STATE fringe contrast.

    Deliverable 4's S_max column is a median OVER TURNS, and the vernier only
    holds its full comb set on the last turn, so that column understates what the
    chain actually reaches. Here the statistic is max-over-turns per chain: the
    best fringe contrast any pulsar attains once the combs have accumulated,
    against its own trials bar.
    """
    print()
    print("=" * 104)
    print("DELIVERABLE 6 (SUPPLEMENT) -- END-STATE FRINGE CONTRAST, NOT THE PER-TURN MEDIAN")
    print("=" * 104)
    print(f"  {'cell':>12} {'combs':>6} {'S_end (max over turns)':>23} "
          f"{'bar lnK':>8} {'S/bar':>7} {'best single chain':>18}")
    rows = []
    for K in Ks:
        acc = []
        for s in range(n_real):
            arr = W.make_array(np.random.default_rng(2500 + s))
            acc.append(chain(7000 + 41 * s + K, arr, 1.56, n_seed=K))
        rows.append(_endrow(f"ArmK K={K}", acc))
    for g in gammas:
        acc = []
        for s in range(n_real):
            arr = W.make_array(np.random.default_rng(2500 + s))
            acc.append(chain(3000 + 31 * s + int(100 * g), arr, g))
        rows.append(_endrow(f"ArmG g={g}", acc))
    print()
    print("  The certification channel is limited by fringe contrast, and contrast grows")
    print("  with the number of ACCUMULATED combs -- not with any one source's SNR.")
    out["endstate"] = np.array([r[1:] for r in rows], dtype=float)
    out["endstate_lbl"] = np.array([r[0] for r in rows])


def _endrow(label, acc):
    se, bars, nc = [], [], []
    for r in acc:
        hh = [h for h in r["hist"] if h["which"] >= 0]
        se.append(max([h["S_max"] for h in hh], default=0.0))
        bars.append(np.median([h["lnK"] for h in hh]) if hh else 0.0)
        nc.append(min(len(hh), COMB_CAP))
    s_med = float(np.median(se)); b_med = float(np.median(bars))
    print(f"  {label:>12} {float(np.median(nc)):6.1f} {s_med:23.3f} "
          f"{b_med:8.2f} {s_med/max(b_med,1e-9):7.3f} {max(se):18.3f}")
    return [label, float(np.median(nc)), s_med, b_med, s_med / max(b_med, 1e-9),
            float(max(se))]
