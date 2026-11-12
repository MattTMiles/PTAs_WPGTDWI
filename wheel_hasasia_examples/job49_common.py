# ================= JOB 49 -- shared machinery: loader + jitter-model switch =======
# Booked code is NOT edited. This module:
#   (1) loads the machine namespace with the Job 47b/c/d sliced fast loader
#       (verbatim cell slices of active_wheel_demo_v4.ipynb + job43c_rebooking.py
#       spans 485-509 inv_rate, 661-689 book, 712-756 collect36m);
#   (2) defines collect36m_sw(): a copy of the booked collect36m with an added
#       `jitter_model` switch.  jitter_model="incumbent" (DEFAULT) is checked
#       record-for-record against the untouched booked collect36m.  In the same
#       pass it builds the PHYSICAL datum from a SEPARATE seeded stream, so the
#       incumbent rng stream (skies, priors, geometry, L_true, eps) is untouched
#       and every arm sees identical skies/pulsars/sources (common random numbers).
#   (3) book_sw(): the booked book() procedure with a selectable outcome column
#       and an optional fixed threshold (for "old gate applied unchanged").
#
# PHYSICAL ARM (per pulsar i, comb k), in units of the source's own amplitude:
#   u_k = 1 + sum_j a_j e^{i theta_j}  [explicit phasors, the M loudest in-window
#                                        contributors to inv_rate's res; ONE phasor
#                                        each (variant p1) or TWO -- Earth + pulsar
#                                        term, independent phases (variant p2)]
#           + CN(0, R_k)                [remaining contributors + unresolved sea,
#                                        exact power, x1 (p1) or x2 (p2)]
#           + CN(0, 2/rho_k^2)          [measurement noise: small-angle var 1/rho^2]
#   a_j^2 = D_k * c_j / Pt_k with c_j = sinc^2 * P_j * (T/T_OBS) * rate (inv_rate's
#   own per-contributor terms) and D_k = depro_sigma's Earth-deprojection factor
#   (1 - f_E E_K) * LOONOISE, so total confusion power == the incumbent's
#   sigma_conf_dep^2 EXACTLY.  theta_j ~ U(0, 2pi) per (pulsar, contributor): the
#   population carries no per-pulsar phases, so they are SUPPLIED here.
#   Fuzzy template term: a pure phase rotation eta ~ N(0, dkap*sigma_d), as in the
#   incumbent (it is a template error, not a power share).
#   datum z_k = e^{i kappa_k L_true} e^{i eta_k} u_k ;  phase error = arg(e^{i eta} u).
# Decoders on the physical datum:
#   D-inc : incumbent wrapped-quadratic decoder, untouched (uses sigma_inc).
#   D-mf  : matched filter  sum_k 2 e^{-s_fz^2/2} Re(z_k e^{-i kappa_k L}) / P_k,
#           P_k = sigma_conf_dep^2 + 2/rho_k^2 (powers the incumbent already assumes).
import os, json, time
import numpy as np

BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
NS = {"__name__": "__main__"}
_LOADED = False
M_EXPLICIT = 64          # explicit phasors per comb; remainder is an exact-power CN


def load(verbose=True):
    global _LOADED
    if _LOADED:
        return NS
    cwd = os.getcwd(); os.chdir(BASE)
    nb = json.load(open("active_wheel_demo_v4.ipynb"))
    cells = {i: "".join(c["source"]) for i, c in enumerate(nb["cells"])
             if c["cell_type"] == "code"}
    # speed-only patch: do not spawn the jug-gpu gate worker (cell 1 comment:
    # "the backend changes speed, never numbers"; no gate surface is used here)
    cells[1] = cells[1].replace(
        '_JUG_PY = "/home/mattm/miniforge3/envs/jug-gpu/bin/python"',
        '_JUG_PY = "/nonexistent/job49_no_worker"')

    def ex(i, lo=None, hi=None, tag=""):
        src = cells[i]
        if lo is not None:
            src = "\n".join(src.splitlines()[lo-1:hi])
        t0 = time.time()
        exec(compile(src, f"<v4-cell-{i}{':'+tag if tag else ''}>", "exec"), NS)
        if verbose:
            print(f"  [cell {i:2d}{' '+tag if tag else ''}] OK {time.time()-t0:6.1f} s",
                  flush=True)
    import contextlib, io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        for i in (1, 3, 5, 7, 9, 11):
            ex(i)
        ex(14, 116, 116, "EXP"); ex(14, 168, 168, "TOOTH_W")
        ex(14, 780, 786, "rn_draw"); ex(14, 1586, 1591, "FP_ALL")
        ex(14, 1072, 1085, "SEED50"); ex(14, 3023, 3029, "_sea50v")
        ex(14, 4043, 4049, "h_of_rhoP"); ex(14, 4154, 4201, "build_pop15")
        ex(15, 13, 17, "s_to_s30"); ex(17, 1, 127, "depro_sigma")
        ex(18, 12, 17, "_pop18"); ex(27, 67, 78, "alias_score24")
    NS["TOOTH_W"] = 0.25
    lines43 = open("job43c_rebooking.py").read().splitlines()
    for lo, hi in ((485, 509), (661, 689), (712, 756)):
        exec(compile("\n".join(lines43[lo-1:hi]), f"<43c:{lo}-{hi}>", "exec"), NS)
    os.chdir(cwd)
    _LOADED = True
    if verbose:
        print("  machine namespace loaded (47b/c/d sliced fast path; gate worker off)")
    return NS


# ---------------------------------------------------------------------------------
def contributors(f2, P0, sea, k, T_yr, rate=1.0):
    """inv_rate's residual for comb k, returned TERM BY TERM.
    Returns (c_others [per in-window source, sinc^2-weighted, x rate],
             sea_term [x rate], Pt, res_check)."""
    T_OBS = NS["T_OBS"]; YR_S = NS["YR_S"]
    bin_edges = NS["bin_edges"]; bin_cent = NS["bin_cent"]
    T = T_yr*YR_S
    o = np.argsort(f2); fs = f2[o]; Pb = P0[o]; half = 100.0/T
    lo = np.searchsorted(fs, f2[k]-half); hi = np.searchsorted(fs, f2[k]+half)
    w = np.sinc((fs[lo:hi]-f2[k])*T)**2
    c = w*Pb[lo:hi]*(T/T_OBS)
    idx = o[lo:hi]
    b = int(np.clip(np.digitize(f2[k], bin_edges)-1, 0, len(bin_cent)-1))
    res = (c.sum()-P0[k]*(T/T_OBS)+sea[b]*(T_OBS/T))*rate        # inv_rate verbatim
    keep = idx != k
    c_oth = c[keep]*rate
    sea_t = sea[b]*(T_OBS/T)*rate
    Pt = P0[k]*(T/T_OBS)
    return c_oth, sea_t, Pt, res


def _sky_contrib(f2, P0, sea, ks, T_yr, rate):
    """Per-comb (explicit amplitudes^2 / Pt, remainder power / Pt) BEFORE depro."""
    out = []
    for k in ks:
        c_oth, sea_t, Pt, res = contributors(f2, P0, sea, k, T_yr, rate)
        o = np.argsort(c_oth)[::-1]
        top = c_oth[o[:M_EXPLICIT]]
        rem = c_oth[o[M_EXPLICIT:]].sum() + sea_t
        out.append((top/Pt, rem/Pt, max(res, 0)/Pt))
    return out


def collect36m_sw(nsky=200, T=45, fuzzy=True, full_dmu=False, marg=False, seed0=8000,
                  rate=None, jitter_model="incumbent", keep_detail=False):
    """Booked collect36m (job43c_rebooking.py:712-756) + jitter_model switch.
    rate=None reproduces collect36m (rt=1 unless marg); a number reproduces the
    T3 `_fixed` harness (job43c:790-829), which is collect36m(full_dmu=True,
    fuzzy=True) at that rate.  Every record carries every outcome:
      wrong_inc          incumbent Gaussian draw, incumbent decoder -- booked
      wrong_<v>_inc      physical datum variant v, incumbent decoder  = D-inc
      wrong_<v>_mf       physical datum variant v, matched filter     = D-mf
    (v in VARIANTS below) and r['wrong'] is set from the SWITCH
    jitter_model in {"incumbent" (default), "p1_inc", "p1_mf", "p2_inc", ...}."""
    _pop18 = NS["_pop18"]; rn_draw = NS["rn_draw"]; inv_rate = NS["inv_rate"]
    depro_sigma = NS["depro_sigma"]; C_PC_S = NS["C_PC_S"]; YR_S = NS["YR_S"]
    TOOTH_W = NS["TOOTH_W"]
    NS['_CAP32'] = 1300
    out = []; rgA = np.random.default_rng(555)
    for j in range(nsky):
        rt = 1.0 if rate is None else rate
        if marg:
            A = 6.4e-15*np.exp(rgA.normal(0, 0.32)); rt = float(np.clip((A/6.4e-15)**2, 0.2, 3.5))
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        rng = np.random.default_rng(seed0+sd+977*j); Np = 120
        beta = np.arcsin(rng.uniform(-1, 1, Np)); Dk = rng.uniform(0.5, 3.0, Np)
        px = 2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/T); pb = 3.0*(T/30.0)**(-2.5)
        pri = np.maximum(np.where(rng.random(Np) < 0.7, np.minimum(px, pb), px), 0.05)
        Arn, _g = rn_draw(Np); rncl = Arn == 0.0
        r = inv_rate(f2, h2, P0, sea, T, rate=rt)
        if r is None: continue
        ks, sig, fk, rp, arr = r
        CB = _sky_contrib(f2, P0, sea, ks, T, rt)            # physical-arm inputs
        sconf = np.sqrt(np.maximum(sig**2-1/np.maximum(rp, 1e-9)**2, 0.0))
        # depro factor D_k: depro_sigma(.., rho->inf)^2 / sigma_conf^2  (loo: i-free)
        Dfac = np.array([depro_sigma(sconf[a], 1e15, fk[a], 0, "loo")**2
                         / max(sconf[a]**2, 1e-300) for a in range(len(ks))])
        for i in range(Np):
            if not rncl[i]: continue
            sg = np.array([depro_sigma(np.sqrt(max(sig[a]**2-1/max(rp[a], 1e-9)**2, 0.0)),
                                       rp[a], fk[a], i, "loo") for a in range(len(ks))])
            ok = sg < np.pi
            if ok.sum() < 2: continue
            sgu = sg[ok]; fku = fk[ok]; arru = arr[ok]; ia = np.where(ok)[0]
            mu = np.arccos(rng.uniform(-0.9, 0.9, len(sgu)))
            kap = 2*np.pi*fku*(1-np.cos(mu))/C_PC_S
            sfz = np.zeros(len(sgu))
            if fuzzy:
                md = 1.0/np.maximum(arru, 1e-6)
                dmu = rng.rayleigh(md*np.sqrt(2/np.pi)) if full_dmu else md
                dfq = 1.0/(T*YR_S*np.maximum(arru, 1e-6))
                dkap = np.sqrt((2*np.pi*(1-np.cos(mu))/C_PC_S)**2*dfq**2
                               + (2*np.pi*fku*np.sin(mu)/C_PC_S)**2*dmu**2)
                sfz = dkap*pri[i]
                sgu = np.sqrt(sgu**2+(dkap*pri[i])**2)
                o2 = sgu < np.pi
                if o2.sum() < 2: continue
                sgu = sgu[o2]; kap = kap[o2]; ia = ia[o2]; sfz = sfz[o2]
            Lt = float(np.clip(rng.normal(0, pri[i]), -2.9*pri[i], 2.9*pri[i]))
            Lg = np.linspace(-3*pri[i], 3*pri[i], 601)
            eps = kap*Lt+rng.normal(0, sgu)                    # INCUMBENT generator
            r_ = (kap[:, None]*Lg[None, :]-eps[:, None]+np.pi) % (2*np.pi)-np.pi
            S = -(r_**2/(2*sgu[:, None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
            if 1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)), 1e-30)) >= TOOTH_W: continue
            Lh = float(Lg[int(np.argmax(S))])
            w_inc = bool(abs(Lh-Lt) > TOOTH_W)
            # ---------------- PHYSICAL ARMS (separate streams; incumbent untouched) --
            rho_k = rp[ia]
            rec = dict(sky=j, wrong=w_inc, sig=sgu.copy(), kap=kap.copy(),
                       pri=float(pri[i]), Ltrue=Lt, wrong_inc=w_inc,
                       K=6*pri[i]/max(2*np.pi/max(kap.min(), 1e-30), 1e-9),
                       rho=rho_k.copy(), sfz=sfz.copy(),
                       Xc=np.array([Dfac[a]*CB[a][2] for a in ia]))
            if keep_detail:
                rec["eps_res"] = eps-kap*Lt
            for tag, mult, gauss in VARIANTS:
                rph = np.random.default_rng([49, seed0, j, i, int(round(rt*1e6)), mult, gauss])
                nC = len(kap)
                u = np.ones(nC, complex); Pn = np.empty(nC)
                for c_, a in enumerate(ia):
                    top, rem, cf = CB[a]
                    D = Dfac[a]
                    Pm = 2.0/max(rho_k[c_], 1e-9)**2
                    if gauss:                 # same total power, complex Gaussian
                        Rg = mult*D*cf
                    else:                     # explicit phasors: mult per contributor
                        th = rph.uniform(0, 2*np.pi, (mult, len(top)))
                        u[c_] += np.sum(np.sqrt(D*top)[None, :]*np.exp(1j*th))
                        Rg = mult*D*rem
                    g = rph.normal(0, 1, 2)*np.sqrt((Rg+Pm)/2)
                    u[c_] += g[0]+1j*g[1]
                    Pn[c_] = mult*D*cf+Pm
                eta = rph.normal(0, 1, nC)*sfz
                dphi = np.angle(u)+eta
                eps_p = kap*Lt+dphi
                r_ = (kap[:, None]*Lg[None, :]-eps_p[:, None]+np.pi) % (2*np.pi)-np.pi
                Sp = -(r_**2/(2*sgu[:, None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
                rec[f"wrong_{tag}_inc"] = bool(abs(float(Lg[int(np.argmax(Sp))])-Lt) > TOOTH_W)
                z = np.exp(1j*(kap*Lt+eta))*u
                wmf = 2.0*np.exp(-sfz**2/2)/Pn
                Sm = (wmf[:, None]*np.real(z[:, None]*np.exp(-1j*kap[:, None]*Lg[None, :]))).sum(0) \
                    - Lg**2/(2*pri[i]**2)
                rec[f"wrong_{tag}_mf"] = bool(abs(float(Lg[int(np.argmax(Sm))])-Lt) > TOOTH_W)
                rec[f"Pn_{tag}"] = Pn
                if keep_detail:
                    rec[f"dphi_{tag}"] = dphi
            if jitter_model != "incumbent":
                rec["wrong"] = rec[f"wrong_{jitter_model}"]
            out.append(rec)
    return out


# physical variants: (tag, terms per other source, complex-Gaussian instead of phasors)
#   p1  : ONE phasor per other source (the incumbent's power ledger) -> isolates the
#         generator change alone (Rice at the incumbent's own share)
#   p2  : TWO phasors per other source (Earth + pulsar term, D0 as specified), sea x2
#   p2g : p2 power, but the confusion as a single complex Gaussian (no explicit phasors)
VARIANTS = (("p1", 1, 0), ("p2", 2, 0), ("p2g", 2, 1))
ARMS = ("inc", "p1_inc", "p1_mf", "p2_inc", "p2_mf", "p2g_inc", "p2g_mf")


def set_arm(rec, arm):
    key = "wrong_inc" if arm in ("inc", "incumbent") else f"wrong_{arm}"
    for r in rec:
        r["wrong"] = r[key]
    return rec


def book_sw(rec, ntr=100, thr_fixed=None, alias=None):
    """book() (job43c:661-689) verbatim procedure; adds: precomputed alias scores
    (they depend on sig/kap/pri/Ltrue only -- identical across arms), the frozen
    threshold value in the output, and an optional FIXED threshold."""
    alias_score24 = NS["alias_score24"]
    sk = np.array([r["sky"] for r in rec]); wr = np.array([r["wrong"] for r in rec], bool)
    al = np.array([alias_score24(r) for r in rec]) if alias is None else alias
    cal = sk < ntr; tst = sk >= ntr; u = np.unique(sk[tst]); nte = len(u)
    qs = np.linspace(0.5, 60, 120); k = np.ones(7)/7.0
    th = [np.quantile(al[cal], q/100) for q in qs]
    pu = []
    for t in th:
        s = cal & (al <= t); pu.append(1-wr[s].mean() if s.sum() >= 10 else np.nan)
    pu = np.array(pu)
    ps = np.convolve(np.nan_to_num(pu, nan=np.nanmean(pu[np.isfinite(pu)])), k, mode="same")
    thr = np.nan; gate = np.nan
    if thr_fixed is None:
        for a in range(len(ps)-1):
            if ps[a] >= 0.90 >= ps[a+1]:
                f = (ps[a]-0.90)/max(ps[a]-ps[a+1], 1e-9)
                thr = th[a]+f*(th[a+1]-th[a]); gate = qs[a]+f*(qs[a+1]-qs[a]); break
        if not np.isfinite(thr):
            return dict(ok=False, gate=np.nan, pur=np.nan, y=np.nan, thr=np.nan,
                        maxpur=float(np.nanmax(ps)))
    else:
        thr = thr_fixed; gate = float(100*np.mean(al[cal] <= thr))
    s = tst & (al <= thr)
    per = {j: (int((s & (sk == j) & ~wr).sum()), int((s & (sk == j)).sum())) for j in u}
    rg = np.random.default_rng(23); by = []; bp = []
    for b in range(400):
        bs = rg.choice(u, size=nte, replace=True)
        r_ = sum(per[j][0] for j in bs); n_ = sum(per[j][1] for j in bs)
        by.append(r_/nte)
        if n_ > 0: bp.append(r_/n_)
    return dict(ok=True, gate=gate, pur=1-wr[s].mean(), y=(s.sum()-wr[s].sum())/nte,
                yb=np.quantile(by, [.1, .9]), pb=np.quantile(bp, [.1, .9]), n=int(s.sum()),
                thr=float(thr), maxpur=float(np.nanmax(ps)))
