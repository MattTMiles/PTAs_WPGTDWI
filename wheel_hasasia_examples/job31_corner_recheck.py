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
def rank_mass_mu(A, kap, pri, Lt, z, Lg, mu0, rg, M=99):
    near = _np.abs(Lg - Lt) <= TOOTH_W
    logprior = -(Lg-mu0)**2/(2*pri**2)
    E = _np.exp(-1j*(kap[:, None]*Lg[None, :]))
    def mass(zz):
        Sc = (2*A[:, None]*_np.real(zz[:, None]*E)).sum(0)
        L = Sc + logprior
        w = _np.exp(L - L.max()); w /= w.sum()
        return float(w[near].sum())
    real = mass(z)
    surr = _np.array([mass(z*_np.exp(1j*rg.uniform(0, 2*_np.pi, len(z))))
                      for _ in range(M)])
    return (float((surr < real).sum()) + 0.5*float((surr == real).sum()))/M
print(f"\n{'='*74}\nJOB 29 ADDENDUM 2 (corrected) -- the TRUNCATION is the gate\n{'='*74}")
print("  Gate chain measured: arr>=3 (rp>=0.61) -> rp>=rho_cut=3 [1220 srcs/sky]")
print("  -> TOP-40-BY-rp TRUNCATION [binds] -> sigma<pi -> top-nmax=12 by nats.")
print("  Lowering rho_cut cannot act: the top-40-by-rp cut discards the added")
print("  weak sources first. The real lever is the truncation cap.")

def sig_windowed(f2, P0, sea, rp, cand, T, W=100.0):
    """inventory20's sigma, windowed to |df| < W/T (sinc^2 tails beyond are
    negligible). Validated against inventory20 below."""
    o = _np.argsort(f2); fs = f2[o]; Ps = P0[o]
    out = _np.empty(len(cand))
    half = W/T
    for a, k in enumerate(cand):
        lo = _np.searchsorted(fs, f2[k]-half); hi = _np.searchsorted(fs, f2[k]+half)
        fw = fs[lo:hi]; Pw = Ps[lo:hi]
        w = _np.sinc((fw - f2[k])*T)**2
        c = w*Pw*(T/T_OBS)
        b = int(_np.clip(_np.digitize(f2[k], bin_edges)-1, 0, len(bin_cent)-1))
        res = c.sum() - (P0[k]*(T/T_OBS)) + sea[b]*(T_OBS/T)
        out[a] = _np.sqrt(max(res,0)/max(P0[k]*(T/T_OBS),1e-300)
                          + 1/max(rp[k],1e-9)**2)
    return out

# --- validation against inventory20 on the standard configuration
sd, f2, h2, P0, sea = _pop18(0, 0.35)
T = 45.0*YR_S
ks, sg_ref, fk, rpk, nats = inventory20(f2, h2, P0, P0, sea, 45.0, 100.0, 14.0, 120)
nt = 45.0*365.25/14.0
rp_all = (h2/(2*_np.pi*f2))*_np.sqrt(nt/2)/(100.0*1e-9)
sg_win = sig_windowed(f2, P0, sea, rp_all, ks, T)
print(f"\n  VALIDATION of windowed sigma vs inventory20 on its own combs:")
print(f"    max |rel diff| = {_np.max(_np.abs(sg_win-sg_ref)/sg_ref):.3e}"
      f"  -> {'OK' if _np.max(_np.abs(sg_win-sg_ref)/sg_ref) < 0.02 else 'MISMATCH'}")

def cap_sweep(cap, rho_cut=3.0, nsky=15, nmax=12):
    ncomb=[]; nadd=[]; sgm=[]
    for j in range(nsky):
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        nt = 45.0*365.25/14.0
        rp = (h2/(2*_np.pi*f2))*_np.sqrt(nt/2)/(100.0*1e-9)
        arr = rp*_np.sqrt(0.2*120)
        cand = _np.where((arr>=3.0)&(f2>=2e-9)&(f2<=3e-8)&(rp>=rho_cut))[0]
        if len(cand) > cap: cand = cand[_np.argsort(rp[cand])[::-1][:cap]]
        if len(cand) < 2: continue
        sg = sig_windowed(f2, P0, sea, rp, cand, T)
        good = sg < _np.pi
        ncomb.append(int(good.sum())); nadd.append(len(cand))
        if good.any(): sgm.append(float(_np.median(sg[good])))
    return (_np.mean(ncomb), _np.mean(nadd), _np.mean(sgm) if sgm else _np.nan)

print(f"\n  (1) TRUNCATION-CAP SWEEP (rho_cut fixed at 3):")
print(f"  {'cap':>8} {'candidates':>11} {'combs (sigma<pi)':>17} {'median sigma':>13}")
for _cap in (40, 100, 400, 1300):
    nc, na, sm = cap_sweep(_cap)
    print(f"  {_cap:8d} {na:11.1f} {nc:17.2f} {sm:13.3f}")
print(f"\n  (1b) lowering rho_cut WITH the cap raised to 1300:")
print(f"  {'rho_cut':>8} {'candidates':>11} {'combs (sigma<pi)':>17} {'median sigma':>13}")
for _rc in (3.0, 2.5, 2.0, 1.5):
    nc, na, sm = cap_sweep(1300, rho_cut=_rc)
    print(f"  {_rc:8.1f} {na:11.1f} {nc:17.2f} {sm:13.3f}")
def inv_cap(f2, h2, P0, sea, T_yr, cap, rho_cut=3.0, nmax=24):
    T = T_yr*YR_S; nt = T_yr*365.25/14.0
    rp = (h2/(2*_np.pi*f2))*_np.sqrt(nt/2)/(100.0*1e-9)
    arr = rp*_np.sqrt(0.2*120)
    cand = _np.where((arr>=3.0)&(f2>=2e-9)&(f2<=3e-8)&(rp>=rho_cut))[0]
    if len(cand) > cap: cand = cand[_np.argsort(rp[cand])[::-1][:cap]]
    if len(cand) < 2: return None
    sg = sig_windowed(f2, P0, sea, rp, cand, T)
    good = sg < _np.pi
    cand, sg = cand[good], sg[good]
    if len(cand) < 2: return None
    nats = _np.log(_np.pi/sg); o = _np.argsort(nats)[::-1][:nmax]
    return cand[o], sg[o], f2[cand[o]], rp[cand[o]], arr[cand[o]]
def alias_of(w, Lg, ip):
    m = _np.abs(Lg - Lg[ip]) > TOOTH_W
    return float(w[m].max()/max(w[ip], 1e-300)) if m.any() else 0.0
print(f"\n{'='*74}\nJOB 31 -- RE-MEASURE THE HEADLINE AT RAISED CAP\n{'='*74}")
_INVC = {}
def inv_cached(j, cap, T_yr=45.0):
    key=(j,cap,T_yr)
    if key not in _INVC:
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        _INVC[key] = inv_cap(f2, h2, P0, sea, T_yr, cap)
    return _INVC[key]

def decode_corner(cap, nsky=100, nps=40, T_yr=45.0, null=False, seed0=3100):
    """Full pipeline: Gaussian in-window prior, fuzzy templates, honest truth
    drawn from the prior, alias-coincidence score per certification."""
    H=[]; C=[]; AL=[]; SK=[]
    for j in range(nsky):
        r = inv_cached(j, cap, T_yr)
        if r is None: continue
        ks, sg0, fk, rpk, arrk = r; T=T_yr*YR_S
        rg=_np.random.default_rng(seed0+j)
        Ti=_np.clip(_np.exp(rg.normal(_np.log(T_yr),0.35,120)),5,90)
        beta=_np.arcsin(rg.uniform(-1,1,120)); Dk=rg.uniform(0.5,3.0,120)
        px=2.34/_np.cos(beta)**2*Dk**2*10.0*_np.sqrt(30.0/Ti); pb=3.0*(Ti/30.0)**(-2.5)
        pri=_np.maximum(_np.where(rg.random(120)<0.7,_np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(120); rncl=Arn==0.0
        for i in _np.where(rncl)[0][:nps]:
            sg=_np.array([depro_sigma(
                _np.sqrt(max(sg0[a]**2-1/max(rpk[a],1e-9)**2,0.0)),
                rpk[a], fk[a], i, "loo") for a in range(len(ks))])
            sg=sg*_np.sqrt(_np.mean(Ti)/Ti[i])
            mu=_np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            sg=sg*_np.sqrt((1/3)/_np.maximum((1+_np.cos(mu))**2/4,1e-3))
            kap=2*_np.pi*fk*(1-_np.cos(mu))/C_PC_S
            dmu=1.0/_np.maximum(arrk,1e-6); dfq=1.0/(T*_np.maximum(arrk,1e-6))
            dkap=_np.sqrt((2*_np.pi*(1-_np.cos(mu))/C_PC_S)**2*dfq**2
                          +(2*_np.pi*fk*_np.sin(mu)/C_PC_S)**2*dmu**2)
            sg=_np.sqrt(sg**2+(dkap*pri[i])**2)
            ok=sg<_np.pi
            if ok.sum()<2: continue
            sgu,kapu=sg[ok],kap[ok]; A=1.0/_np.maximum(sgu,1e-6)
            Lt=rg.normal(0,pri[i])                       # HONEST truth from prior
            if null:
                z=(rg.normal(0,1/_np.sqrt(2),len(A))+1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            else:
                z=A*_np.exp(1j*(kapu*Lt))+(rg.normal(0,1/_np.sqrt(2),len(A))
                                           +1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            Lg=_np.linspace(-3*pri[i],3*pri[i],_GRID)
            Sc=(2*A[:,None]*_np.real(z[:,None]*_np.exp(
                -1j*(kapu[:,None]*Lg[None,:])))).sum(0)-Lg**2/(2*pri[i]**2)
            w=_np.exp(Sc-Sc.max()); w/=w.sum()
            ip=int(_np.argmax(Sc))
            H.append(abs(float(Lg[ip])-Lt)<=TOOTH_W)
            C.append(1.0/_np.sqrt(max(float(_np.sum(kapu**2/sgu**2)),1e-30))<TOOTH_W)
            AL.append(alias_of(w,Lg,ip)); SK.append(j)
    return (_np.array(H), _np.array(C), _np.array(AL), _np.array(SK))

def corner(H,C,AL,SK,nsky,target=0.90,thr=None):
    """Yield (right psr/sky) at the tightest gate meeting `target` purity."""
    if not C.any(): return dict(pur=_np.nan,y=_np.nan,thr=_np.nan,q=_np.nan)
    al=AL[C]; hh=H[C]
    if thr is not None:
        m=al<=thr
        return dict(pur=hh[m].mean() if m.any() else _np.nan,
                    y=hh[m].sum()/nsky, thr=thr, q=m.mean())
    best=dict(pur=_np.nan,y=0.0,thr=_np.nan,q=_np.nan)
    for q in _np.linspace(0.01,1.0,100):
        t=_np.quantile(al,q); m=al<=t
        if m.sum()<5: continue
        p=hh[m].mean()
        if p>=target and hh[m].sum()/nsky>best["y"]:
            best=dict(pur=p,y=hh[m].sum()/nsky,thr=float(t),q=float(q))
    return best

print(f"\n  T1 -- THE CORNER at each cap [100 skies, identical realisations]")
print(f"  {'cap':>6} {'ungated pur':>12} {'cert/sky':>9} {'>=90% purity':>14}"
      f" {'psr/sky':>8} {'gate q':>7} | {'20%':>12} {'10%':>12} {'7%':>12}")
_D={}
for _cap in (40, 400, 1300):
    H,C,AL,SK = decode_corner(_cap, nsky=100); _D[_cap]=(H,C,AL,SK)
    b90 = corner(H,C,AL,SK,100,0.90)
    row=f"  {_cap:6d} {H.mean():12.4f} {C.sum()/100:9.2f} {b90['pur']:14.4f}" \
        f" {b90['y']:8.2f} {b90['q']:7.3f} |"
    for q in (0.20,0.10,0.07):
        t=_np.quantile(AL[C],q); r=corner(H,C,AL,SK,100,thr=t)
        row+=f" {r['pur']:.3f}/{r['y']:5.2f}"
    print(row)
def decode_v(cap, nsky=100, nps=40, T_yr=45.0, null=False, seed0=3100,
             couple=True, fuzzy=True):
    H=[];C=[];AL=[];SK=[]
    for j in range(nsky):
        r=inv_cached(j,cap,T_yr)
        if r is None: continue
        ks,sg0,fk,rpk,arrk=r; T=T_yr*YR_S
        rg=_np.random.default_rng(seed0+j)
        Ti=_np.clip(_np.exp(rg.normal(_np.log(T_yr),0.35,120)),5,90)
        beta=_np.arcsin(rg.uniform(-1,1,120)); Dk=rg.uniform(0.5,3.0,120)
        px=2.34/_np.cos(beta)**2*Dk**2*10.0*_np.sqrt(30.0/Ti); pb=3.0*(Ti/30.0)**(-2.5)
        pri=_np.maximum(_np.where(rg.random(120)<0.7,_np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(120); rncl=Arn==0.0
        for i in _np.where(rncl)[0][:nps]:
            sg=_np.array([depro_sigma(
                _np.sqrt(max(sg0[a]**2-1/max(rpk[a],1e-9)**2,0.0)),
                rpk[a],fk[a],i,"loo") for a in range(len(ks))])
            sg=sg*_np.sqrt(_np.mean(Ti)/Ti[i])
            mu=_np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            th=mu if couple else _np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            sg=sg*_np.sqrt((1/3)/_np.maximum((1+_np.cos(th))**2/4,1e-3))
            kap=2*_np.pi*fk*(1-_np.cos(mu))/C_PC_S
            if fuzzy:
                dmu=1.0/_np.maximum(arrk,1e-6); dfq=1.0/(T*_np.maximum(arrk,1e-6))
                dkap=_np.sqrt((2*_np.pi*(1-_np.cos(mu))/C_PC_S)**2*dfq**2
                              +(2*_np.pi*fk*_np.sin(mu)/C_PC_S)**2*dmu**2)
                sg=_np.sqrt(sg**2+(dkap*pri[i])**2)
            ok=sg<_np.pi
            if ok.sum()<2: continue
            sgu,kapu=sg[ok],kap[ok]; A=1.0/_np.maximum(sgu,1e-6)
            Lt=rg.normal(0,pri[i])
            z=(rg.normal(0,1/_np.sqrt(2),len(A))+1j*rg.normal(0,1/_np.sqrt(2),len(A))) \
              if null else A*_np.exp(1j*(kapu*Lt))+(rg.normal(0,1/_np.sqrt(2),len(A))
                                                    +1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            Lg=_np.linspace(-3*pri[i],3*pri[i],_GRID)
            Sc=(2*A[:,None]*_np.real(z[:,None]*_np.exp(
                -1j*(kapu[:,None]*Lg[None,:])))).sum(0)-Lg**2/(2*pri[i]**2)
            w=_np.exp(Sc-Sc.max()); w/=w.sum(); ip=int(_np.argmax(Sc))
            H.append(abs(float(Lg[ip])-Lt)<=TOOTH_W)
            C.append(1.0/_np.sqrt(max(float(_np.sum(kapu**2/sgu**2)),1e-30))<TOOTH_W)
            AL.append(alias_of(w,Lg,ip)); SK.append(j)
    return (_np.array(H),_np.array(C),_np.array(AL),_np.array(SK))

print(f"\n  ANCHOR DIAGNOSTIC -- why cap-40 here is not the banked 6.5 psr/sky:")
print(f"  {'harness variant':>34} {'>=90% pur':>10} {'psr/sky':>8}")
for _lab,_cp,_fz in (("banked-like: decoupled, no fuzzy",False,False),
                     ("+ coupled theta=mu (Job 28 fix)",True,False),
                     ("+ fuzzy templates (Job 30)",True,True)):
    H,C,AL,SK=decode_v(40,nsky=100,couple=_cp,fuzzy=_fz)
    b=corner(H,C,AL,SK,100,0.90)
    print(f"  {_lab:>34} {b['pur']:10.4f} {b['y']:8.2f}")
print(f"    NOTE: banked 6.5 psr/sky came from the NOTEBOOK V10/V14 machinery")
print(f"    (collect24d/alias_score24), not this harness. Absolute anchoring to")
print(f"    6.5 therefore requires a notebook-side re-run; the CAP RATIO below is")
print(f"    like-for-like on identical realisations and is the quotable quantity.")

print(f"\n  T1b -- VALIDATION STACK at each cap:")
print(f"  {'cap':>6} {'in-sample':>10} {'OUT-of-sample (frozen thr)':>27}"
      f" {'null psr/sky':>13} {'100->300 |d|':>13}")
for _cap in (40,400,1300):
    H,C,AL,SK=_D[_cap]
    tr=SK<50; te=SK>=50
    b=corner(H[tr],C[tr],AL[tr],SK[tr],50,0.90)
    o=corner(H[te],C[te],AL[te],SK[te],50,thr=b['thr'])
    Hn,Cn,ALn,SKn=decode_v(_cap,nsky=40,null=True)
    nn=corner(Hn,Cn,ALn,SKn,40,thr=b['thr'])
    H3,C3,AL3,SK3=decode_v(_cap,nsky=300)
    b3=corner(H3,C3,AL3,SK3,300,0.90)
    b1=corner(H,C,AL,SK,100,0.90)
    print(f"  {_cap:6d} {b['pur']:10.4f} {o['pur']:14.4f}/{o['y']:6.2f} psr/sky"
          f" {nn['y']:13.3f} {abs(b3['y']-b1['y']):13.3f}")

print(f"\n  T1c -- P_F CONTRACT NUMBERS (P_F = 1 - purity):")
print(f"  {'cap':>6} {'P_F=0.10 psr/sky':>17} {'P_F=0.01 psr/sky':>17}")
for _cap in (40,400,1300):
    H,C,AL,SK=_D[_cap]
    a=corner(H,C,AL,SK,100,0.90); b=corner(H,C,AL,SK,100,0.99)
    print(f"  {_cap:6d} {a['y']:17.2f} {b['y']:17.2f}")

print(f"\n  T2 -- ONE-SHOT turn-1 at cap 1300:")
print(f"  {'config':>22} {'n cert':>7} {'right':>6} {'purity':>7} {'right/RNc':>10}")
for _T in (45.0, 60.0):
    for _cap in (40, 1300):
        H,C,AL,SK=decode_v(_cap,nsky=60,T_yr=_T)
        n=int(C.sum()); rt=int(H[C].sum())
        print(f"  {'T='+str(int(_T))+', cap '+str(_cap):>22} {n:7d} {rt:6d}"
              f" {rt/max(n,1):7.4f} {rt/max(len(H),1):10.4f}")

print(f"\n  T3 -- JOB 28 T1 ALIGNED RUNGS re-checked at cap 1300:")
print(f"  {'variant':>34} {'K min':>7} {'purity':>8} {'yield':>7} {'right/RNc':>10}")
for _cap in (40, 1300):
    for _lab,_hi in (("baseline (-0.9,0.9)",0.9), ("aligned tail (-0.9,1.0)",1.0)):
        H=[];C=[];KM=[]
        for j in range(30):
            r=inv_cached(j,_cap,45.0)
            if r is None: continue
            ks,sg0,fk,rpk,arrk=r
            rg=_np.random.default_rng(3300+j)
            Ti=_np.clip(_np.exp(rg.normal(_np.log(45.0),0.35,120)),5,90)
            beta=_np.arcsin(rg.uniform(-1,1,120)); Dk=rg.uniform(0.5,3.0,120)
            px=2.34/_np.cos(beta)**2*Dk**2*10.0*_np.sqrt(30.0/Ti); pb=3.0*(Ti/30.0)**(-2.5)
            pri=_np.maximum(_np.where(rg.random(120)<0.7,_np.minimum(px,pb),px),0.05)
            Arn,_g=rn_draw(120); rncl=Arn==0.0
            for i in _np.where(rncl)[0][:30]:
                sg=_np.array([depro_sigma(
                    _np.sqrt(max(sg0[a]**2-1/max(rpk[a],1e-9)**2,0.0)),
                    rpk[a],fk[a],i,"loo") for a in range(len(ks))])
                sg=sg*_np.sqrt(_np.mean(Ti)/Ti[i])
                mu=_np.arccos(rg.uniform(-0.9,_hi,len(ks)))
                sg=sg*_np.sqrt((1/3)/_np.maximum((1+_np.cos(mu))**2/4,1e-3))
                kap=2*_np.pi*fk*(1-_np.cos(mu))/C_PC_S
                ok=sg<_np.pi
                if ok.sum()<2: continue
                sgu,kapu=sg[ok],kap[ok]; A=1.0/_np.maximum(sgu,1e-6)
                KM.append(float((6*pri[i]/(2*_np.pi/_np.maximum(kapu,1e-30))).min()))
                Lt=rg.normal(0,pri[i])
                z=A*_np.exp(1j*(kapu*Lt))+(rg.normal(0,1/_np.sqrt(2),len(A))
                                           +1j*rg.normal(0,1/_np.sqrt(2),len(A)))
                Lg=_np.linspace(-3*pri[i],3*pri[i],_GRID)
                Sc=(2*A[:,None]*_np.real(z[:,None]*_np.exp(
                    -1j*(kapu[:,None]*Lg[None,:])))).sum(0)-Lg**2/(2*pri[i]**2)
                H.append(abs(float(Lg[int(_np.argmax(Sc))])-Lt)<=TOOTH_W)
                C.append(1.0/_np.sqrt(max(float(_np.sum(kapu**2/sgu**2)),1e-30))<TOOTH_W)
        H=_np.array(H);C=_np.array(C)
        print(f"  {('cap '+str(_cap)+': '+_lab):>34} {_np.median(KM):7.2f}"
              f" {H[C].mean():8.4f} {C.mean():7.4f} {H[C].sum()/max(len(H),1):10.4f}")
