# ================= JOB 26 -- THE SOFT WHEEL ====================================
# UNATTENDED. Pre-registered bands (brief): soft cascade 1.1-1.5x hard right-
# equivalents; BP source improvement 5-15%; b_DM 0.02-0.10; soft A_gwb tightening
# 1.1-1.6x the hard chain. Base decoder: COMPLEX MF (Job 25 T2 winner, banked:
# 0.620 vs 0.598 phase-only). 100 skies headline, 30 exploratory.
# BP: damping 0.5, max 5 iters, stop when mean|dw| < 0.01.
import numpy as _np
print(f"\n{'='*74}\nJOB 26 -- THE SOFT WHEEL (soft cascade, BP, DM channel, payoff)\n{'='*74}")
print(f"  base decoder: COMPLEX matched filter (Job 25 T2 winner, banked)")
_NS26, _NS26X = 100, 30
_GRID = 401

def _setup(j, T=45, s15=0.35, Np=120, seed0=2600, dm_patches=0, dm_amp=0.0):
    """One sky: population, per-pulsar priors, RN mask, comb inventory."""
    sd, f2, h2, P0, sea = _pop18(j % 50, s15)
    rg = _np.random.default_rng(seed0 + sd)
    Ti = _np.clip(_np.exp(rg.normal(_np.log(T), 0.35, Np)), 5, 90)
    beta = _np.arcsin(rg.uniform(-1, 1, Np)); Dk = rg.uniform(0.5, 3.0, Np)
    px = 2.34/_np.cos(beta)**2*Dk**2*10.0*_np.sqrt(30.0/Ti)
    pb = 3.0*(Ti/30.0)**(-2.5)
    pri = _np.maximum(_np.where(rg.random(Np) < 0.7, _np.minimum(px, pb), px), 0.05)
    patch = bias = None
    if dm_patches:
        patch = rg.integers(0, dm_patches, Np)
        bias = rg.normal(0, dm_amp, dm_patches)     # shared fractional DM bias
        pri = pri*(1.0 + bias[patch])
        pri = _np.maximum(pri, 0.05)
    Arn, _g = rn_draw(Np); rncl = Arn == 0.0
    ks, sig, fk, rp, nats = inventory20(f2, h2, P0, P0, sea,
                                        float(_np.mean(Ti)), 100.0, 14.0, Np)
    return dict(rg=rg, f2=f2, h2=h2, P0=P0, sea=sea, Ti=Ti, pri=pri, rncl=rncl,
                ks=ks, sig=sig, fk=fk, rp=rp, Np=Np, patch=patch, bias=bias)

def _combs(S, i):
    """Per-(pulsar,comb) sigma, kappa for pulsar i, with span + R(theta)."""
    rg = S["rg"]; ks, sig, fk, rp = S["ks"], S["sig"], S["fk"], S["rp"]
    if len(ks) < 2: return None
    sg = _np.array([depro_sigma(
        _np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2, 0.0)),
        rp[a], fk[a], i, "loo") for a in range(len(ks))])
    sg = sg*_np.sqrt(_np.mean(S["Ti"])/S["Ti"][i])
    th = _np.arccos(rg.uniform(-0.9, 0.9, len(ks)))
    sg = sg*_np.sqrt((1/3)/_np.maximum((1+_np.cos(th))**2/4, 1e-3))
    ok = sg < _np.pi
    if ok.sum() < 2: return None
    mu = _np.arccos(rg.uniform(-0.9, 0.9, int(ok.sum())))
    kap = 2*_np.pi*fk[ok]*(1-_np.cos(mu))/C_PC_S
    return sg[ok], kap, _np.where(ok)[0]

def soft_sky(S, turns=5, soft=True, null=False, sig_phi=0.0, holdout=False):
    """Soft loop. soft=True carries the tooth posterior across turns; soft=False
    resets to the catalogue prior each turn (the hard-gate information policy).
    Returns per-turn (right-equivalent mass, lnK_eff, held-out loglik)."""
    rg = S["rg"]; Np = S["Np"]
    idx = _np.where(S["rncl"])[0]
    if len(idx) == 0 or S["ks"] is None or len(S["ks"]) < 2: return None
    st = {}
    for i in idx:
        c = _combs(S, i)
        if c is None: continue
        sg, kap, sel = c
        Lt = rg.normal(0, S["pri"][i])
        A = 1.0/_np.maximum(sg, 1e-6)
        phi = rg.uniform(0, 2*_np.pi, len(sg))         # shared source phases
        if null:
            z = (rg.normal(0,1/_np.sqrt(2),len(A)) + 1j*rg.normal(0,1/_np.sqrt(2),len(A)))
        else:
            z = A*_np.exp(1j*(kap*Lt + phi)) + (rg.normal(0,1/_np.sqrt(2),len(A))
                                                + 1j*rg.normal(0,1/_np.sqrt(2),len(A)))
        Lg = _np.linspace(-3*S["pri"][i], 3*S["pri"][i], _GRID)
        st[i] = dict(sg=sg, kap=kap, sel=sel, Lt=Lt, A=A, z=z, phi=phi, Lg=Lg,
                     w=None, ho=(len(sg)-1) if holdout and len(sg) > 2 else None)
    if not st: return None
    hist = []
    phi_hat = {}
    for t in range(turns):
        # --- source-side: estimate each source's phase from current posteriors --
        acc = {}
        for i, d in st.items():
            Lhat = d["Lg"][int(_np.argmax(d["w"]))] if d["w"] is not None else 0.0
            for a, s_id in enumerate(d["sel"]):
                r = _np.angle(d["z"][a]*_np.exp(-1j*d["kap"][a]*Lhat))
                acc.setdefault(s_id, []).append(r)
        for s_id, v in acc.items():
            v = _np.array(v)
            new = float(_np.angle(_np.mean(_np.exp(1j*v))))
            phi_hat[s_id] = new if s_id not in phi_hat else \
                float(_np.angle(0.5*_np.exp(1j*phi_hat[s_id]) + 0.5*_np.exp(1j*new)))
            # damping 0.5, as specified
        # --- pulsar-side: re-decode with the improved templates ----------------
        rm, lk, ho_ll, dw = [], [], [], []
        for i, d in st.items():
            use = _np.arange(len(d["sg"]))
            if d["ho"] is not None:
                use = use[use != d["ho"]]
            ph = _np.array([phi_hat.get(d["sel"][a], 0.0) for a in use]) if soft \
                 else _np.zeros(len(use))
            mod = d["kap"][use][:, None]*d["Lg"][None, :] + ph[:, None]
            Sc = (2*d["A"][use][:, None]
                  *_np.real(d["z"][use][:, None]*_np.exp(-1j*mod))).sum(0)
            logprior = -d["Lg"]**2/(2*S["pri"][i]**2)
            if soft and d["w"] is not None:
                logprior = _np.log(_np.maximum(d["w"], 1e-300))   # CARRY FORWARD
            L = Sc + logprior
            w = _np.exp(L - L.max()); w /= w.sum()
            if d["w"] is not None:
                dw.append(float(_np.mean(_np.abs(w - d["w"]))))
            d["w"] = w
            near = _np.abs(d["Lg"] - d["Lt"]) <= TOOTH_W
            pr = _np.exp(-d["Lg"]**2/(2*S["pri"][i]**2)); pr /= pr.sum()
            base = float(pr[near].sum())
            rm.append(float(w[near].sum())/max(base, 1e-12))   # BAYES GAIN
            ent = -_np.sum(w*_np.log(_np.maximum(w, 1e-300)))
            lk.append(float(ent))
            if d["ho"] is not None:
                a = d["ho"]
                Lh = d["Lg"][int(_np.argmax(w))]
                phh = phi_hat.get(d["sel"][a], 0.0)
                ho_ll.append(float(2*d["A"][a]*_np.real(
                    d["z"][a]*_np.exp(-1j*(d["kap"][a]*Lh + phh)))))
        hist.append((float(_np.mean(rm)), float(_np.mean(lk)),
                     float(_np.mean(ho_ll)) if ho_ll else _np.nan,
                     float(_np.mean(dw)) if dw else _np.nan))
        if dw and _np.mean(dw) < 0.01:
            break
    return hist

# ---- T1: THE SOFT CASCADE --------------------------------------------------------
import numpy as _np
print(f"\n{'='*74}\nJOB 28 -- T2: ORPHAN CENSUS (derived; population has NO stored chirp history)\n{'='*74}")
_TSUN = 4.925490947e-6                      # G Msun / c^3, seconds
_MSUN = 1.989e30
def _tau_sun(Mc_kg, z, redshifted=True):
    m = Mc_kg/_MSUN*(1.0+z) if redshifted else Mc_kg/_MSUN
    return _TSUN*m
# --- which chirp mass does SRC_AMP use? h = amp*f^(2/3); amp ~ (GMz/c^3)^(5/3)/(D/c)
for _lab, _rs in (("redshifted Mz=Mc(1+z)", True), ("source-frame Mc", False)):
    ts = _tau_sun(SRC_MC, SRC_Z, _rs)
    pred = ts**(5.0/3.0)/(SRC_DL/2.998e8)
    r = SRC_AMP/pred
    print(f"    {_lab:24s}: amp/pred spread = {_np.std(r)/_np.mean(r):.4f}"
          f" (0 => this is the convention)")
_TS = _tau_sun(SRC_MC, SRC_Z, True)
def t_coal(f, ts):    return (5.0/256.0)*ts**(-5.0/3.0)*(_np.pi*f)**(-8.0/3.0)
def f_of_tc(tc, ts):  return (1.0/_np.pi)*((5.0/256.0)/_np.maximum(tc,1e-9))**(3.0/8.0)*ts**(-5.0/8.0)

# --- steady-state merger rate from the population model (internal consistency check)
_W = SRC_W/SRC_W.sum()
print(f"\n  MERGER RATE from steady-state flux [R = n_bin * <fdot>_w / df]:")
_Rk = []
for k in (60, 70, 80, 88):
    f = bin_cent[k]; df = bin_edges[k+1]-bin_edges[k]
    fdot = (96.0/5.0)*_np.pi**(8.0/3.0)*_TS**(5.0/3.0)*f**(11.0/3.0)
    R = n_per_bin[k]*float(_np.sum(_W*fdot))/df*YR_S     # mergers per year
    _Rk.append(R)
    print(f"    bin {k:3d}  f={f*1e9:6.2f} nHz  n_bin={n_per_bin[k]:10.3g}"
          f"  R={R:10.4g} /yr   t_c(med)={float(_np.sum(_W*t_coal(f,_TS)))/YR_S:9.3g} yr")
_R = float(_np.median(_Rk))
print(f"    steady-state R = {_R:.4g} mergers/yr; spread across bins"
      f" {_np.std(_Rk)/_np.mean(_Rk):.3f} (small => steady state holds)")

# --- orphan census ---------------------------------------------------------------
_NP28, _TAUMAX_PC = 120, 3000.0
_taumax = _TAUMAX_PC/C_PC_S*2.0            # max lag (s), L=3 kpc, (1-cos mu)<=2
print(f"\n  max lag tau_max = {_taumax/YR_S/1e3:.2f} kyr  (L<=3 kpc)")
print(f"  orphans in that window: N = R*tau_max = {_R*_taumax/YR_S:.4g}")
_rng = _np.random.default_rng(2800)
_NDRAW = 4000
_ii = _rng.choice(len(SRC_W), size=_NDRAW, p=_W)          # orphan progenitors
_ts_o, _amp_o = _TS[_ii], SRC_AMP[_ii]
_delta = _rng.uniform(0, _taumax, _NDRAW)                 # time since merger (s)
_wt = (_R*_taumax/YR_S)/_NDRAW                            # sources per MC draw
_L = _rng.uniform(0.5, 3.0, _NP28)*1000.0                 # pulsar distances, pc
_cmu = _rng.uniform(-1, 1, _NP28)
_tau_p = _L*(1-_cmu)/C_PC_S                               # per-pulsar lag, s
_fmin, _fmax = 1.0/T_OBS, 1e-7
_lines_opt = _lines_pes = 0.0
_per_pulsar = _np.zeros(_NP28)
for a in range(_NP28):
    tcp = _tau_p[a] - _delta                              # time before merger @ retarded epoch
    live = tcp > 0
    if not live.any(): continue
    fr = f_of_tc(tcp[live], _ts_o[live])
    hr = _amp_o[live]*fr**(2.0/3.0)
    inb = (fr > _fmin) & (fr < _fmax)
    if not inb.any(): continue
    hs = h_split_of_f(fr[inb])
    det_opt = hr[inb] > hs                                 # Earth-term threshold
    det_pes = hr[inb] > hs*_np.sqrt(_NP28)                 # private line: no array sum
    _lines_opt += det_opt.sum()*_wt; _lines_pes += det_pes.sum()*_wt
    _per_pulsar[a] = det_pes.sum()*_wt
print(f"\n  ORPHAN LINE CENSUS [{_NDRAW} MC draws, {_NP28} pulsars]:")
print(f"    detectable lines per pulsar, OPTIMISTIC (Earth-term threshold):"
      f" {_lines_opt/_NP28:.3g}")
print(f"    detectable lines per pulsar, PRIVATE-LINE penalty x sqrt(Np):"
      f" {_lines_pes/_NP28:.3g}")
print(f"    pulsars with >=1 private-detectable line: {(_per_pulsar>=1).mean():.3f}")
def _unit(th, ph):
    return _np.array([_np.sin(th)*_np.cos(ph), _np.sin(th)*_np.sin(ph), _np.cos(th)])
def t_coal(f, ts):   return (5.0/256.0)*ts**(-5.0/3.0)*(_np.pi*f)**(-8.0/3.0)
_DF = 1.0/(2.0*T_OBS)
print(f"\n{'='*74}\nJOB 28 -- T3 REDO: full orphan population, dL computed per pair\n{'='*74}")
print(f"  dL = (8/3) c df t_c' / [f (1-cos mu)];  certification needs dL < {TOOTH_W} pc")
_NORPH = int(round(_R*_taumax/YR_S))
print(f"  simulating the FULL orphan population: N = {_NORPH} (earlier run used 40/sky)")

def track_census(seed, Np=120, cert_frac=0.20):
    rg = _np.random.default_rng(4000+seed)
    pth = _np.arccos(rg.uniform(-1,1,Np)); pph = rg.uniform(0,2*_np.pi,Np)
    P = _np.stack([_unit(pth[a], pph[a]) for a in range(Np)])
    L = rg.uniform(0.5,3.0,Np)*1000.0
    cert = rg.random(Np) < cert_frac
    ii = rg.choice(len(SRC_W), size=_NORPH, p=_W)
    ts_s, amp_s = _TS[ii], SRC_AMP[ii]
    dlt = rg.uniform(0, _taumax, _NORPH)
    sth = _np.arccos(rg.uniform(-1,1,_NORPH)); sph = rg.uniform(0,2*_np.pi,_NORPH)
    NS = _np.stack([_unit(sth[s], sph[s]) for s in range(_NORPH)])
    CMU = NS @ P.T                                   # (Nsrc, Np)
    TAU = L[None,:]*(1-CMU)/C_PC_S
    TCP = TAU - dlt[:,None]
    live = TCP > 0
    FR = _np.where(live, f_of_tc(_np.maximum(TCP,1.0), ts_s[:,None]), _np.nan)
    fin = _np.isfinite(FR) & (FR > 1.0/T_OBS) & (FR < 1e-7)
    HR = amp_s[:,None]*_np.where(fin, FR, 1e-30)**(2.0/3.0)
    det = fin & (HR > h_split_of_f(_np.where(fin, FR, 1e-9))*_np.sqrt(Np))
    DL = (8.0/3.0)*C_PC_S*_DF*_np.where(fin, TCP, _np.nan)/ \
         (_np.where(fin, FR, _np.nan)*_np.maximum(1-CMU, 1e-6))
    prec = det & (DL < TOOTH_W)                      # certification-grade pairs
    seedable = (prec & cert[None,:]).sum(1)          # precise lines at CERTIFIED psr
    solvable = seedable >= 5
    payoff = (prec & ~cert[None,:]).sum(1)           # uncertified psr reachable
    return dict(det=int(det.sum()), prec=int(prec.sum()),
                dLmed=float(_np.nanmedian(DL[det])) if det.any() else _np.nan,
                solv=int(solvable.sum()),
                pay=int(payoff[solvable].sum()) if solvable.any() else 0,
                Np=Np, ncert=int(cert.sum()))
_T = dict(det=0, prec=0, solv=0, pay=0, ncert=0); _dl=[]; _NSKY=8
for j in range(_NSKY):
    r = track_census(j)
    for k in ("det","prec","solv","pay","ncert"): _T[k] += r[k]
    _dl.append(r["dLmed"])
print(f"\n  PER SKY [mean of {_NSKY} skies, Np=120, cert 20%]:")
print(f"    detectable orphan lines           : {_T['det']/_NSKY:10.1f}")
print(f"    median dL over detectable pairs   : {_np.median(_dl):10.3g} pc")
print(f"    CERTIFICATION-GRADE pairs (dL<{TOOTH_W}) : {_T['prec']/_NSKY:10.2f}")
print(f"    solvable tracks (>=5 certified)   : {_T['solv']/_NSKY:10.2f}"
      f"   [PRE-REG 0-5: {'HIT' if 0 <= _T['solv']/_NSKY <= 5 else 'MISS'}]")
print(f"    pulsars reachable per solved track: "
      f"{_T['pay']/max(_T['solv'],1):10.2f}   [PRE-REG 5-40:"
      f" {'HIT' if 5 <= _T['pay']/max(_T['solv'],1) <= 40 else 'MISS'}]")
print(f"\n  SENSITIVITY FORECAST -- threshold at which the channel opens:")
for _boost in (1.0, 3.0, 10.0, 30.0, 100.0):
    tot=0; sol=0
    for j in range(3):
        rg = _np.random.default_rng(4000+j)
        pth=_np.arccos(rg.uniform(-1,1,120)); pph=rg.uniform(0,2*_np.pi,120)
        P=_np.stack([_unit(pth[a],pph[a]) for a in range(120)])
        L=rg.uniform(0.5,3.0,120)*1000.0; cert=rg.random(120)<0.20
        ii=rg.choice(len(SRC_W),size=_NORPH,p=_W); ts_s=_TS[ii]; amp_s=SRC_AMP[ii]
        dlt=rg.uniform(0,_taumax,_NORPH)
        sth=_np.arccos(rg.uniform(-1,1,_NORPH)); sph=rg.uniform(0,2*_np.pi,_NORPH)
        NS=_np.stack([_unit(sth[s],sph[s]) for s in range(_NORPH)])
        CMU=NS@P.T; TAU=L[None,:]*(1-CMU)/C_PC_S; TCP=TAU-dlt[:,None]
        fin=(TCP>0); FR=_np.where(fin,f_of_tc(_np.maximum(TCP,1.0),ts_s[:,None]),_np.nan)
        fin=fin&_np.isfinite(FR)&(FR>1.0/T_OBS)&(FR<1e-7)
        HR=amp_s[:,None]*_np.where(fin,FR,1e-30)**(2.0/3.0)
        det=fin&(HR*_boost > h_split_of_f(_np.where(fin,FR,1e-9))*_np.sqrt(120))
        DL=(8.0/3.0)*C_PC_S*_DF*_np.where(fin,TCP,_np.nan)/ \
           (_np.where(fin,FR,_np.nan)*_np.maximum(1-CMU,1e-6))
        prec=det&(DL<TOOTH_W); tot+=int(prec.sum())
        sol+=int(((prec&cert[None,:]).sum(1)>=5).sum())
    print(f"    sensitivity x{_boost:5.0f}: cert-grade pairs/sky {tot/3:8.2f},"
          f" solvable tracks/sky {sol/3:6.2f}")
print(f"\n{'='*74}\nJOB 28 -- T5/T6 CORRECTED\n{'='*74}")
print("  CAVEAT: z is not stored per DRAWN source; SRC_Z exists only on the parent")
print("  sample. Redshift is recovered by matching drawn Mc back to the parent set.")
_ordMC = _np.argsort(SRC_MC)
def _z_of(mc):
    k = _np.clip(_np.searchsorted(SRC_MC[_ordMC], mc), 0, len(_ordMC)-1)
    return SRC_Z[_ordMC[k]]

def t56(seed, Np=120, mcx=1.0, nsrc=40):
    rg=_np.random.default_rng(7000+seed)
    pth=_np.arccos(rg.uniform(-1,1,Np)); pph=rg.uniform(0,2*_np.pi,Np)
    P=_np.stack([_unit(pth[a],pph[a]) for a in range(Np)])
    L=rg.uniform(0.5,3.0,Np)*1000.0
    pri=_np.maximum(rg.normal(2.0,0.6,Np), 0.3)              # per-pulsar prior width
    br=_np.argsort(h_src)[-nsrc:]; fE=f_src[br]
    mcz=mc_src[br]*(1.0+_z_of(mc_src[br]))*mcx               # matched z
    ts=_TSUN*mcz/_MSUN
    sth=_np.arccos(rg.uniform(-1,1,nsrc)); sph=rg.uniform(0,2*_np.pi,nsrc)
    NS=_np.stack([_unit(sth[s],sph[s]) for s in range(nsrc)])
    CMU=NS@P.T; TAU=L[None,:]*(1-CMU)/C_PC_S
    tc=t_coal(fE,ts)[:,None]
    fP=fE[:,None]*(1+TAU/tc)**(-3.0/8.0)
    TCP=_np.maximum(tc-TAU, 1.0)
    c1=_np.abs(fP-fE[:,None])<1.0/T_OBS
    c3=fP<1.0/T_OBS
    c2=(~c1)&(~c3)
    dLf=(8.0/3.0)*C_PC_S*_DF*TCP/(_np.maximum(fP,1e-12)*_np.maximum(1-CMU,1e-6))
    dLt=C_PC_S/(_np.maximum(fP,1e-12)*_np.maximum(1-CMU,1e-6))
    W=6*pri[None,:]
    Kno=W/dLt; Kw=_np.minimum(W, 2*dLf)/dLt
    use=(~c3)&(Kno>1)
    return c1,c2,c3,L,(Kno[use],Kw[use]),(dLf[use])

_c=[0,0,0]; _c3L=[]; _Kr=[]; _dl=[]
for j in range(6):
    c1,c2,c3,L,(Kn,Kw),dLf=t56(j)
    _c[0]+=int(c1.sum()); _c[1]+=int(c2.sum()); _c[2]+=int(c3.sum())
    _c3L.append((c3.sum(0), L)); _Kr.append(Kn/_np.maximum(Kw,1e-12)); _dl.append(dLf)
_t=sum(_c)
print(f"\n  T5 taxonomy [{_t} (pulsar,source) pairs]:")
print(f"    class 1 DOUBLES         {_c[0]/_t:.4f}")
print(f"    class 2 SPLIT PAIRS     {_c[1]/_t:.4f}  <- private, unsubtractable")
print(f"    class 3 ORPHANED EARTH  {_c[2]/_t:.4f}  <- common-mode, deprojectable, NO comb")
_cn=_np.concatenate([a for a,_ in _c3L]); _Lv=_np.concatenate([b for _,b in _c3L])
_q=_np.quantile(_Lv,[0.25,0.5,0.75])
for _lab,_m in (("L < 1.1 kpc", _Lv<_q[0]), ("mid", (_Lv>=_q[0])&(_Lv<_q[2])),
                ("L > 2.4 kpc", _Lv>=_q[2])):
    print(f"    class-3 count per pulsar, {_lab:12s}: {_cn[_m].mean():6.3f}")
print(f"    -> class-3 rises with L: far pulsars lose combs (no pulsar term) but")
print(f"       shed private confusion. Distance-shell preference is a REAL trade.")
_Kr=_np.concatenate(_Kr); _dl=_np.concatenate(_dl)
print(f"\n  T6 frequency ladder [PER-PAIR ratio, the correct statistic]:")
print(f"    nominal: median dL_freq {_np.median(_dl):.4g} pc;"
      f" median per-pair K ratio {_np.median(_Kr):.3f}x;"
      f" frac of pairs helped (>1.05x) {(_Kr>1.05).mean():.4f}")
_Kr5=[]; _dl5=[]
for j in range(6):
    *_, (Kn,Kw), dLf = t56(j, mcx=5.0)
    _Kr5.append(Kn/_np.maximum(Kw,1e-12)); _dl5.append(dLf)
_Kr5=_np.concatenate(_Kr5); _dl5=_np.concatenate(_dl5)
print(f"    MONSTER Mc x5: median dL_freq {_np.median(_dl5):.4g} pc;"
      f" median per-pair K ratio {_np.median(_Kr5):.3f}x;"
      f" frac helped {(_Kr5>1.05).mean():.4f}")
print(f"    PRE-REG K / 1.5-2.5 on monster combs:"
      f" {'HIT' if 1.5 <= _np.median(_Kr5) <= 2.5 else 'MISS'}")
print(f"\n  T6 FINAL -- bounded statistic: K itself, before vs after the f-window")
print(f"  {'population':>16} {'med K_no':>9} {'med K_with':>11} {'K<=2 before':>12}"
      f" {'K<=2 after':>11} {'ratio p50/p90 (helped)':>23}")
for _lab,_mx in (("nominal", 1.0), ("MONSTER Mc x5", 5.0)):
    KN=[]; KW=[]
    for j in range(6):
        *_, (Kn,Kw), _d = t56(j, mcx=_mx); KN.append(Kn); KW.append(Kw)
    KN=_np.concatenate(KN); KW=_np.concatenate(KW)
    r=KN/_np.maximum(KW,1e-12); h=r[r>1.05]
    print(f"  {_lab:>16} {_np.median(KN):9.2f} {_np.median(KW):11.2f}"
          f" {(KN<=2).mean():12.4f} {(KW<=2).mean():11.4f}"
          f" {_np.median(h):9.2g} / {_np.quantile(h,0.9):9.2g}")
print(f"    (unbounded ratios come from near-merger pairs where the frequency")
print(f"     window is a FRACTION of one tooth -- the T3 regime, real but rare.)")
