# ================= JOB 50 -- collector for the 12-arm experiment =====================
# Base arm = Job 49 `p2_mf` (two-term physical datum, matched filter).  Factors:
#   retardation  off | on      (on = arm `p2r`: every pulsar term, target included, at
#                               f_P = f_E (1+tau/t_c)^(-3/8) with the retarded amplitude,
#                               tau = L_i (1 - cos mu_ij)/c for THIS pulsar and source)
#   sources      top12 | all   (all = every inventory candidate: rho>=3, arr>=3, 2-30 nHz,
#                               cap 1300 by rho -- inv_rate's own candidate set)
#   template     oracle | machine (rms dkappa*sigma_d, booked) | full (rms dkappa*L_true)
# The booked rng stream (job43c collect36m) is consumed exactly as booked; the Job 49
# `p2` stream is reused for the booked decode list, so (off, top12, machine) on the
# booked 601-point grid IS p2_mf, record for record.  Everything new comes from a
# separate stream rng([50, seed0, sky, pulsar]).
#
# What is modelled how (stated limits):
#  * booked decode sources: explicit phasors (64 loudest contributors, Earth and pulsar
#    term separately) + exact-power complex Gaussian remainder, as in Job 49.
#  * every OTHER candidate: its own datum u = 1 + CN(0, D X + 2/rho^2), X = two-term
#    confusion ratio for that (pulsar, source).  Each candidate has an independent
#    noise draw although ~1100 candidates share ~45 frequency bins; other candidates
#    are counted as confusion for each one (incoherent stack, not a joint fit).
#  * chirp masses: Job 41 recipe (job41_displaced_excess.py:199-205), mass and redshift
#    assigned by amplitude rank from the presentation-sky population.
#  * Earth-deprojection factor D(f): the booked blended factor, evaluated at the
#    target's frequency, applied to the whole confusion (as p2 does).
#  * a source's own Earth term is treated as removed (as booked).
#  * decode grids: "b" = booked 601 points over +-3 sigma_d; "f" = fine grid, step
#    <= min(0.05 pc, finest tooth/6), evaluated by a gridded FFT (phase error <= 0.15 rad).
import os, sys, time
import numpy as np
import scipy.fft as sfft
from scipy.interpolate import CubicSpline
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job49_common as J
NS = J.NS
TSUN = 4.925490947e-6; MSUN = 1.989e30
RET = ("off", "on"); SRC = ("top12", "all"); TPL = ("oracle", "machine", "full")
ARMS12 = [(r, s, t) for r in RET for s in SRC for t in TPL]
ARMNAME = lambda a: f"{a[0]}|{a[1]}|{a[2]}"
NFFT_MAX = 2**24
_CACHE = {}
# chirp-mass convention for the retardation lag:
#   "z"   redshifted Mc(1+z)  (Job 41 recipe; consistent with observed f and observer-time lag)
#   "src" source-frame Mc     (notebook FP_ALL, cell 14 l.1590)
MASS = os.environ.get("JOB50_MASS", "z")


def _mcz_tc(f2, h2):
    if "ordA" not in _CACHE:
        _CACHE["ordA"] = np.argsort(NS["SRC_AMP"]); _CACHE["AS"] = NS["SRC_AMP"][_CACHE["ordA"]]
    amp = h2/np.maximum(f2, 1e-30)**(2.0/3.0)
    k = np.clip(np.searchsorted(_CACHE["AS"], amp), 0, len(_CACHE["AS"])-1); i = _CACHE["ordA"][k]
    ts = TSUN*NS["SRC_MC"][i]*((1.0+NS["SRC_Z"][i]) if MASS == "z" else 1.0)/MSUN
    return (5.0/256.0)*ts**(-5.0/3.0)*(np.pi*np.maximum(f2, 1e-30))**(-8.0/3.0)


def f_ret(f, tc, tau):
    return f*(1.0+tau/np.maximum(tc, 1e-9))**(-3.0/8.0)


NAT_SEED0 = 5000          # natural sky j uses population seed NAT_SEED0 + j (never used by the notebook: its seeds are 3..52)
NAT_TC = {}               # {sky: chirp time of every source from its own drawn (Mc, z)}
NAT_REJ = {}              # {sky: sources rejected by the CW upper limit during the draw}


def natural_pop(j):
    """Raw population draw, nothing rescaled: (seed, f, h, P0 = h^2 f T_OBS, sea = closed-form faint power per bin)."""
    sd = NAT_SEED0+j
    f, h, _mc, _d, S, nrej = NS["draw_population"](seed=sd)
    NAT_REJ[j] = int(nrej)
    # the drawn catalogue entry, recovered exactly (the amplitude lookup of _mcz_tc lands on a
    # neighbouring entry for about 5% of sources through rounding: job53_diag_mc.py)
    _mcz_tc(f[:1], h[:1]); oA = _CACHE["ordA"]; AS = _CACHE["AS"]
    k = np.clip(np.searchsorted(AS, h/f**(2.0/3.0)), 0, len(AS)-1); ti = np.full(len(f), -1)
    for dk in (0, -1, 1, -2, 2):
        ii = oA[np.clip(k+dk, 0, len(AS)-1)]
        m = (ti < 0) & (NS["SRC_MC"][ii] == _mc) & (NS["SRC_DL"][ii] == _d); ti[m] = ii[m]
    assert (ti >= 0).all()
    ts = TSUN*NS["SRC_MC"][ti]*((1.0+NS["SRC_Z"][ti]) if MASS == "z" else 1.0)/MSUN
    NAT_TC[j] = (5.0/256.0)*ts**(-5.0/3.0)*(np.pi*f)**(-8.0/3.0)
    return sd, f, h, h**2*f*NS["T_OBS"], S


def census_priors():
    import csv
    rows = list(csv.DictReader(open(os.path.join(J.BASE, "job45_t2_parallax_census.csv"))))
    return (np.array([float(r["D_kpc"]) for r in rows]), np.array([float(r["sigma_d_now_pc"]) for r in rows]))


def _decode_fine(c, kap, sd, kmax):
    """argmax over a fine uniform grid of S(L) = Re sum_k c_k e^{-i kap_k L} - L^2/(2 sd^2).
    c: (n_arm, n_src) complex.  Returns L_hat per arm and the grid step used."""
    lam_min = 2*np.pi/max(kmax, 1e-30)
    dL = min(0.05, lam_min/6.0)
    dk_need = 0.3/(3.0*sd)       # nearest-bin phase error <= 0.15 rad at the window edge
    N = int(2**np.ceil(np.log2(max((2*np.pi/dL)/dk_need, 1024))))
    if N > NFFT_MAX:
        N = NFFT_MAX; dL = max(dL, (2*np.pi)/(N*dk_need))
    dk = 2*np.pi/(N*dL)
    M = min(int(np.floor(6*sd/dL))+1, N)
    L0 = -3.0*sd
    q = np.rint(kap/dk).astype(np.int64)
    ok = q < N
    cc = (c*np.exp(-1j*kap*L0)[None, :])[:, ok]; q = q[ok]
    Lg = L0+dL*np.arange(M)
    pr = -Lg**2/(2*sd**2)
    out = np.empty(c.shape[0])
    for a in range(c.shape[0]):
        arr = (np.bincount(q, weights=cc[a].real, minlength=N)
               + 1j*np.bincount(q, weights=cc[a].imag, minlength=N)).astype(np.complex64)
        S = sfft.fft(arr, workers=1).real[:M]+pr
        out[a] = Lg[int(np.argmax(S))]
    return out, dL


def collect50(skies, T=45, marg=True, rate=None, seed0=8000, sens_n=1.0, priors="harness",
              fine=True, arms=ARMS12, verbose=False, tau0=False, hook=None, arr_ovr=None,
              fmax_cand=3e-8, cand_cap=1300, fgrid_max=3.3e-8, sky="constructed", mass="exact", noskip=False, cand_bg=0):
    # noskip (Job 54): keep skies whose booked list has fewer than two usable sources (default: skip, as Jobs 49-53)
    # mass (natural skies only): "exact" = each source's own drawn (Mc, z); "lookup" = _mcz_tc as Jobs 50-52
    # Job 53 addition (default = Jobs 49-52): sky = "constructed" (_pop18: three sources per
    # sky rescaled by build_pop15) or "natural" (raw draw_population, one distinct seed per sky)
    # Job 52 additions (defaults = Job 50/51): fmax_cand upper Earth-frequency edge of the
    # candidate pool, cand_cap its size cap, fgrid_max top of the 1/(8T) field grid
    # Job 51 additions (defaults leave Job 50 behaviour bit-identical):
    #   hook     object with .psr(locals) -> dict merged into the record, .sky(locals) -> dict
    #            merged into the sky info
    #   arr_ovr  {sky: array over candidates} replacing the array SNR in dmu, dfq (template term)
    """Returns (records, sky_info).  One record per red-noise-clean pulsar per sky."""
    _pop18 = NS["_pop18"]; rn_draw = NS["rn_draw"]; inv_rate = NS["inv_rate"]
    depro_sigma = NS["depro_sigma"]; C = NS["C_PC_S"]; YR_S = NS["YR_S"]; T_OBS = NS["T_OBS"]
    TOOTH_W = NS["TOOTH_W"]; be = NS["bin_edges"]; bc = NS["bin_cent"]; alias24 = NS["alias_score24"]
    Ts = T*YR_S; half = 100.0/Ts; sc = Ts/T_OBS
    if priors == "census":
        cD, cS = census_priors()
    # marg amplitude stream is sequential over skies 0..: replay it
    rgA = np.random.default_rng(555); rts = []
    for j in range(max(skies)+1):
        A = 6.4e-15*np.exp(rgA.normal(0, 0.32)); rts.append(float(np.clip((A/6.4e-15)**2, 0.2, 3.5)))
    out = []; info = []
    for j in skies:
        rt = rts[j] if marg else (1.0 if rate is None else rate)
        if sky == "natural": sd, f2, h2, P0, sea = natural_pop(j)
        elif sky == "constructed": sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        else: raise ValueError(sky)
        rng = np.random.default_rng(seed0+sd+977*j); Np = 120
        beta = np.arcsin(rng.uniform(-1, 1, Np)); Dk = rng.uniform(0.5, 3.0, Np)
        px = 2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/T); pb = 3.0*(T/30.0)**(-2.5)
        pri = np.maximum(np.where(rng.random(Np) < 0.7, np.minimum(px, pb), px), 0.05)
        Arn, _g = rn_draw(Np); rncl = Arn == 0.0
        r = inv_rate(f2, h2, P0, sea, T, rate=rt)
        if r is None:
            if not noskip: continue
            ks = np.zeros(0, int); sig = fk = rp = arr = np.zeros(0)     # Job 54: no usable booked list; the joint fit does not need one
        else:
            ks, sig, fk, rp, arr = r
        if priors == "census":
            rc = np.random.default_rng([50, 7, seed0, j]); pick = rc.integers(0, len(cD), Np)
            Dk = cD[pick].copy(); pri = np.maximum(cS[pick], 0.05)
        CB = J._sky_contrib(f2, P0, sea, ks, T, rt)
        sconf = np.sqrt(np.maximum(sig**2-1/np.maximum(rp, 1e-9)**2, 0.0))
        Dfac12 = np.array([depro_sigma(sconf[a], 1e15, fk[a], 0, "loo")**2/max(sconf[a]**2, 1e-300)
                           for a in range(len(ks))])
        # ---- candidates (inv_rate l.490-494, with timing noise / sens_n) --------------
        nt = T*365.25/14.0
        rp_all = (h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100.0*1e-9)*sens_n
        arr_all = rp_all*np.sqrt(0.2*120)
        cand = np.where((arr_all >= 3.0) & (f2 >= 2e-9) & (f2 <= fmax_cand) & (rp_all >= 3.0))[0]
        if len(cand) > cand_cap: cand = cand[np.argsort(rp_all[cand])[::-1][:cand_cap]]
        if cand_bg:      # Job 54: also admit the cand_bg in-band sources loudest against timing noise + background, whatever their single-pulsar SNR
            tg = np.where((arr_all >= 3.0) & (f2 >= 2e-9) & (f2 <= fmax_cand))[0]; gEt = np.zeros(len(tg))
            for lo in range(0, len(f2), 4000):
                gEt += (np.sinc((f2[lo:lo+4000, None]-f2[None, tg])*Ts)**2*P0[lo:lo+4000, None]).sum(0)
            bs_ = sea[np.clip(np.digitize(f2[tg], be)-1, 0, len(bc)-1)]*(T_OBS/Ts)
            oth_ = np.maximum(2*(gEt-P0[tg])*sc+2*bs_, 0)*rt; Nt_ = P0[tg]*sc*2.0/np.maximum(rp_all[tg], 1e-9)**2
            cand = np.union1d(cand, tg[np.argsort(arr_all[tg]/np.sqrt(1.0+oth_/Nt_))[::-1][:cand_bg]])
        cand = np.union1d(cand, ks)                      # booked list always inside
        pos12 = np.searchsorted(cand, ks)                # position of booked sources
        nC = len(cand); fc = f2[cand]; Pc = P0[cand]
        tc_all = NAT_TC[j] if (sky == "natural" and mass == "exact") else _mcz_tc(f2, h2)
        # static Earth-term field gE(f) = sum_j sinc^2((f_j-f)T) P_j (T/T_OBS) on a fine grid
        fg8 = np.arange(0.3e-9, fgrid_max, 1.0/(8*Ts)); fg8T = (fg8*Ts).astype(np.float32)
        gE8 = np.zeros(len(fg8)); f32 = (f2*Ts).astype(np.float32); P32 = P0.astype(np.float32)
        for lo in range(0, len(f2), 6000):
            gE8 += (np.sinc(f32[lo:lo+6000, None]-fg8T[None, :])**2*P32[lo:lo+6000, None]).sum(0, dtype=np.float64)
        gEs = CubicSpline(fg8, gE8*sc)
        bsea = lambda f: sea[np.clip(np.digitize(f, be)-1, 0, len(bc)-1)]*(T_OBS/Ts)
        Dc = lambda f: np.array([depro_sigma(1.0, 1e15, x, 0, "loo")**2 for x in np.atleast_1d(f)])
        D_off = Dc(fc)
        X_off = np.maximum(2*(gEs(fc)-Pc*sc)+2*bsea(fc), 0)*rt/(Pc*sc)
        rho_off = rp_all[cand]
        anchors = np.zeros(nC, int); nmov = []; ndrift = [0, 0]
        for i in range(Np):
            if not rncl[i]: continue
            # ---------------- booked stream, exactly as collect36m ----------------------
            sg = np.array([depro_sigma(np.sqrt(max(sig[a]**2-1/max(rp[a], 1e-9)**2, 0.0)),
                                       rp[a], fk[a], i, "loo") for a in range(len(ks))])
            ok = sg < np.pi; ia_ok = np.where(ok)[0]
            booked = False; Lt = None; ia = np.zeros(0, int); mu12 = None; dmu12 = None
            sgu = kap_b = None
            if ok.sum() >= 2:
                sgu = sg[ok]; fku = fk[ok]; arru = arr[ok]
                mu12 = np.arccos(rng.uniform(-0.9, 0.9, len(sgu)))
                kap_b = 2*np.pi*fku*(1-np.cos(mu12))/C
                md = 1.0/np.maximum(arru, 1e-6)
                dmu12 = rng.rayleigh(md*np.sqrt(2/np.pi))
                dfq = 1.0/(T*YR_S*np.maximum(arru, 1e-6))
                dkap = np.sqrt((2*np.pi*(1-np.cos(mu12))/C)**2*dfq**2
                               + (2*np.pi*fku*np.sin(mu12)/C)**2*dmu12**2)
                sgu = np.sqrt(sgu**2+(dkap*pri[i])**2)
                o2 = sgu < np.pi
                if o2.sum() >= 2:
                    sgu = sgu[o2]; kap_b = kap_b[o2]; ia = ia_ok[o2]
                    Lt = float(np.clip(rng.normal(0, pri[i]), -2.9*pri[i], 2.9*pri[i]))
                    rng.normal(0, sgu)                                  # booked eps draw
                    booked = 1.0/np.sqrt(max(float(np.sum(kap_b**2/sgu**2)), 1e-30)) < TOOTH_W
            # ---------------- auxiliary stream ------------------------------------------
            rx = np.random.default_rng([50, seed0, j, i])
            Lt_aux = float(np.clip(rx.normal(0, 1), -2.9, 2.9))*pri[i]
            if Lt is None: Lt = Lt_aux
            mu_all = np.arccos(rx.uniform(-0.9, 0.9, len(f2)))
            un = rx.normal(0, 1, (3, nC)); uray = rx.rayleigh(1.0, nC)
            if mu12 is not None:
                mu_all[ks[ia_ok]] = mu12
            mu = mu_all[cand]
            arrc = arr_all[cand] if arr_ovr is None else arr_ovr[j]
            dmu = uray*np.sqrt(2/np.pi)/np.maximum(arrc, 1e-6)
            if dmu12 is not None:
                dmu[pos12[ia_ok]] = dmu12 if arr_ovr is None else \
                    dmu12*arr[ia_ok]/np.maximum(arrc[pos12[ia_ok]], 1e-6)
            dfq = 1.0/(T*YR_S*np.maximum(arrc, 1e-6))
            dkap = np.sqrt((2*np.pi*(1-np.cos(mu))/C)**2*dfq**2+(2*np.pi*fc*np.sin(mu)/C)**2*dmu**2)
            Lpc = Dk[i]*1e3
            sfz = {"oracle": np.zeros(nC), "machine": dkap*pri[i], "full": dkap*Lpc}
            # ---------------- retardation geometry --------------------------------------
            tau_all = Lpc*(1-np.cos(mu_all))/C*(0.0 if tau0 else 1.0)
            fP_all = f_ret(f2, tc_all, tau_all); pf_all = (fP_all/f2)**(7.0/3.0)
            fPc = fP_all[cand]; pfc = pf_all[cand]
            mov = np.abs(fP_all-f2)*Ts > 0.02; nmov.append(int(mov.sum()))
            drift = (3.0/8.0)*(fPc/np.maximum(tc_all[cand]+tau_all[cand], 1e-9))*Ts*Ts
            # field of all terms at the retarded target frequencies: Earth field (static
            # spline) + this pulsar's pulsar-term field, summed on a 1/(8T) grid (the
            # field is band-limited at 1/(2T)) and cubic-splined to the targets
            gPg = np.zeros(len(fg8)); wP = (P0*pf_all).astype(np.float32)
            fP32 = (fP_all*Ts).astype(np.float32)
            for lo in range(0, len(f2), 6000):
                gPg += (np.sinc(fP32[lo:lo+6000, None]-fg8T[None, :])**2*wP[lo:lo+6000, None]).sum(0, dtype=np.float64)
            gt = gEs(fPc)+CubicSpline(fg8, gPg*sc)(fPc)
            own = sc*Pc*(np.sinc((fc-fPc)*Ts)**2+pfc)
            Pt_on = Pc*pfc*sc
            X_on = np.maximum(gt-own+2*bsea(fPc), 0)*rt/Pt_on
            D_on = Dc(fPc)
            rho_on = rho_off*(fPc/fc)**(-1.0/3.0)
            kap = {"off": 2*np.pi*fc*(1-np.cos(mu))/C, "on": 2*np.pi*fPc*(1-np.cos(mu))/C}
            anchors += (6*pri[i]*kap["off"]/(2*np.pi) < 1.0)
            # ---------------- data: u for every candidate, both retardation states -------
            u = {}; Pn = {}
            g1 = (un[0]+1j*un[1])/np.sqrt(2)
            for rs, X, D, rho in (("off", X_off, D_off, rho_off), ("on", X_on, D_on, rho_on)):
                Pm = 2.0/np.maximum(rho, 1e-9)**2
                Pn[rs] = D*X+Pm
                u[rs] = 1.0+g1*np.sqrt(Pn[rs])
            # booked decode list: explicit phasors, Job 49 p2 stream
            eta_unit = un[2].copy()
            if len(ia):
                rph = np.random.default_rng([49, seed0, j, i, int(round(rt*1e6)), 2, 0])
                for c_, a in enumerate(ia):
                    p = pos12[a]; k = ks[a]
                    top, rem, cf = CB[a]; D = Dfac12[a]
                    th = rph.uniform(0, 2*np.pi, (2, len(top)))
                    g = rph.normal(0, 1, 2)
                    Pm = 2.0/max(rp[a]*sens_n, 1e-9)**2
                    u["off"][p] = 1.0+np.sum(np.sqrt(D*top)[None, :]*np.exp(1j*th)) \
                        + (g[0]+1j*g[1])*np.sqrt((2*D*rem+Pm)/2)
                    Pn["off"][p] = 2*D*cf+Pm
                    # retarded: exact per-source Earth and pulsar contributions at f_P(target)
                    ft = fPc[p]
                    win = (f2 >= ft-half) & (f2 < ft+half); win[k] = False
                    cE = np.sinc((f2[win]-ft)*Ts)**2*P0[win]*sc*rt
                    cP = np.sinc((fP_all[win]-ft)*Ts)**2*P0[win]*pf_all[win]*sc*rt
                    o = np.argsort(cE+cP)[::-1][:len(top)]
                    Ptk = Pt_on[p]; Dn = D_on[p]
                    seat = 2*float(np.atleast_1d(bsea(ft))[0])*rt
                    remn = max((cE.sum()+cP.sum()-cE[o].sum()-cP[o].sum()+seat), 0)/Ptk
                    Pmn = 2.0/max(rho_on[p], 1e-9)**2
                    u["on"][p] = 1.0+np.sum(np.sqrt(Dn*cE[o]/Ptk)*np.exp(1j*th[0, :len(o)])) \
                        + np.sum(np.sqrt(Dn*cP[o]/Ptk)*np.exp(1j*th[1, :len(o)])) \
                        + (g[0]+1j*g[1])*np.sqrt((Dn*remn+Pmn)/2)
                    Pn["on"][p] = Dn*(cE.sum()+cP.sum()+seat)/Ptk+Pmn
                eta_unit[pos12[ia]] = rph.normal(0, 1, len(ia))
            # ---------------- source lists -----------------------------------------------
            m_all = np.ones(nC, bool)
            m_ok = np.zeros(nC, bool); m_ok[pos12[ia_ok]] = True
            m_ia = np.zeros(nC, bool); m_ia[pos12[ia]] = True
            lists = {("top12", "oracle"): m_ok, ("top12", "machine"): m_ia, ("top12", "full"): m_ok,
                     ("all", "oracle"): m_all, ("all", "machine"): m_all, ("all", "full"): m_all}
            ref = ia if len(ia) else (ia_ok if len(ia_ok) else np.arange(len(ks)))
            kref = kap["off"][pos12[ref]]
            rec = dict(sky=j, psr=i, pri=float(pri[i]), Ltrue=Lt, L_pc=float(Lpc), booked=bool(booked),
                       Kc=float(6*pri[i]*kref.min()/(2*np.pi)) if len(kref) else np.nan,
                       Kf=float(6*pri[i]*kref.max()/(2*np.pi)) if len(kref) else np.nan,
                       Kf_all=float(6*pri[i]*kap["off"].max()/(2*np.pi)), ncand=int(nC),
                       n12=int(len(ia)), nmov=int(mov.sum()))
            if len(ia):
                rec["Pn12_off"] = Pn["off"][pos12[ia]].copy(); rec["Pn12_on"] = Pn["on"][pos12[ia]].copy()
                rec["pf12"] = pfc[pos12[ia]].copy(); rec["dfbin12"] = ((fc-fPc)*Ts)[pos12[ia]].copy()
            if booked:
                rec["alias"] = float(alias24(dict(pri=float(pri[i]), Ltrue=Lt, kap=kap_b, sig=sgu)))
                ndrift[0] += int((drift[pos12[ia]] > 1).sum()); ndrift[1] += len(ia)
            Lg = np.linspace(-3*pri[i], 3*pri[i], 601)
            for rs in RET:
                kp = kap[rs]
                Eb = np.exp(-1j*kp[:, None]*Lg[None, :])
                cs = []; keys = []
                for (r_, s_, t_) in arms:
                    if r_ != rs: continue
                    m = lists[(s_, t_)]
                    w = np.where(m, 2.0*np.exp(-sfz[t_]**2/2)/Pn[rs], 0.0)
                    z = np.exp(1j*(kp*Lt+eta_unit*sfz[t_]))*u[rs]
                    cs.append(w*z); keys.append(ARMNAME((r_, s_, t_)))
                    rec["n|"+keys[-1]] = int(m.sum())
                    rec["i|"+keys[-1]] = float(np.sum(np.where(m, np.exp(-sfz[t_]**2)*2.0/Pn[rs], 0.0)))
                cs = np.array(cs)
                Sb = np.real(cs @ Eb)-Lg[None, :]**2/(2*pri[i]**2)
                Lb = Lg[np.argmax(Sb, axis=1)]
                for n_, kname in enumerate(keys):
                    rec["wb|"+kname] = bool(abs(Lb[n_]-Lt) > TOOTH_W)
                if fine:
                    Lf, dL = _decode_fine(cs, kp, pri[i], float(kp.max()))
                    for n_, kname in enumerate(keys):
                        rec["wf|"+kname] = bool(abs(Lf[n_]-Lt) > TOOTH_W)
                    rec["dL_"+rs] = float(dL)
            if hook is not None: rec.update(hook.psr(locals()))
            out.append(rec)
        info.append(dict(sky=j, rate=rt, ncand=int(nC), n12=int(len(ks)), anchors=anchors,
                         is12=np.isin(cand, ks), nmov_med=float(np.median(nmov)), nsrc=int(len(f2)),
                         drift_gt1=ndrift[0], drift_n=ndrift[1],
                         share_off=1/(1+D_off[pos12]*X_off[pos12])))
        if hook is not None: info[-1].update(hook.sky(locals()))
        if verbose:
            print(f"  sky {j}: cand {nC}, movers/psr {np.median(nmov):.0f}", flush=True)
    return out, info
