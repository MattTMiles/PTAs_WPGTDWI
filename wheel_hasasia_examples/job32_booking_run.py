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
import numpy as np
print(f"\n{'='*74}\nJOB 32 -- T1: BOOK THE CORNER (V10/V14 machinery, notebook-side)\n{'='*74}")
_CAP32 = 40
_inv20_orig = inventory20
def inventory20(f2, h2, Psig, Pbkg, sea, T_yr, sig_ns, cad, Np, depth_dex=0.0,
                rho_cut=3.0, nmax=12):
    """Banked inventory20, verbatim, EXCEPT the hard-coded top-40 truncation is
    replaced by the module-level _CAP32, and sigma is computed on a +/-100/T
    window (validated to 7e-7 against the original)."""
    T = T_yr*YR_S
    nt = T_yr*365.25/cad
    rp = (h2/(2*np.pi*f2))*np.sqrt(nt/2)/(sig_ns*1e-9)
    arr = rp*np.sqrt(0.2*Np)
    det = (arr >= 3.0*10**(-depth_dex)) & (f2 >= 2e-9) & (f2 <= 3e-8)
    ks = np.where(det & (rp >= rho_cut))[0]
    if len(ks) == 0:
        return ks, np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)
    if len(ks) > _CAP32:
        ks = ks[np.argsort(rp[ks])[::-1][:_CAP32]]
    o_ = np.argsort(f2); fs = f2[o_]; Pb = Pbkg[o_]
    half = 100.0/T
    sig = np.empty(len(ks))
    for a, k in enumerate(ks):
        lo = np.searchsorted(fs, f2[k]-half); hi = np.searchsorted(fs, f2[k]+half)
        w = np.sinc((fs[lo:hi] - f2[k])*T)**2
        c = w*Pb[lo:hi]*(T/T_OBS)
        b = int(np.clip(np.digitize(f2[k], bin_edges)-1, 0, len(bin_cent)-1))
        res = c.sum() - Pbkg[k]*(T/T_OBS) + sea[b]*(T_OBS/T)
        Pt = Psig[k]*(T/T_OBS)
        sig[a] = np.sqrt(max(res, 0)/max(Pt, 1e-300) + 1/max(rp[k], 1e-9)**2)
    good = sig < np.pi
    ks, sig = ks[good], sig[good]
    if len(ks) == 0:
        return ks, sig, np.zeros(0), np.zeros(0), np.zeros(0)
    nats = np.log(np.pi/sig)
    o = np.argsort(nats)[::-1][:nmax]
    return ks[o], sig[o], f2[ks[o]], rp[ks[o]], nats[o]

def corner32(cap, nsky=60, ntrain=30, prior="gauss", null=False):
    global _CAP32
    _CAP32 = cap
    c = collect24d(nsky=nsky, inject=1.0, prior=prior, null=null)
    if not c: return None
    w = np.array([r["wrong"] for r in c], bool)
    sk = np.array([r["sky"] for r in c])
    asr = np.array([alias_score24(r) for r in c])
    tr, te = sk < ntrain, sk >= ntrain
    nte = nsky-ntrain
    rows={}; C90=None
    for q in (20,10,7,5,3):
        thr = float(np.percentile(asr[tr], q))
        sel = te & (asr <= thr)
        if sel.sum() < 5: continue
        pur = 1-w[sel].mean(); nright = int(sel.sum()-w[sel].sum())
        rows[q]=(pur, nright/nte, int(sel.sum()))
        if pur >= 0.90 and C90 is None: C90=(q,pur,nright/nte)
    return dict(n=len(c), ung=1-w.mean(), rows=rows, C90=C90, asr=asr, w=w, sk=sk)

print(f"\n  PRE-FLIGHT GATE -- reproduce the banked corner at cap 40:")
_r40 = corner32(40)
print(f"    n {_r40['n']} (banked 4423), ungated purity {_r40['ung']:.3f} (banked 0.544)")
print(f"    {'gate':>6} {'purity':>8} {'psr/sky':>9}   [banked: 20%->15.3, 10%->8.6, 7%->6.5]")
for q,(p,y,nn) in _r40["rows"].items():
    print(f"    {q:5d}% {p:8.3f} {y:9.1f}")
_ok = (_r40["C90"] is not None and abs(_r40["C90"][2]-6.5) <= 1.0
       and abs(_r40["ung"]-0.544) <= 0.02)
print(f"    banked corner 6.5 @ 0.903 (7% gate); reproduced"
      f" {_r40['C90'][2]:.1f} @ {_r40['C90'][1]:.3f} (gate {_r40['C90'][0]}%)"
      if _r40["C90"] else "    NO >=90% corner reached")
print(f"    PRE-FLIGHT: {'PASS -- raised-cap rows may be read' if _ok else 'FAIL -- STOP, diagnose harness gap'}")
if _ok:
    print(f"\n  T1 -- BOOKED CORNER vs cap [collect24d, honest truth, train-frozen]:")
    print(f"  {'cap':>6} {'n':>6} {'ungated':>8} {'>=90% gate':>11} {'purity':>8}"
          f" {'psr/sky':>8} | {'20%':>13} {'10%':>13} {'7%':>13}")
    _B={}
    for _cap in (40,400,1300):
        r = corner32(_cap) if _cap!=40 else _r40
        _B[_cap]=r
        s=f"  {_cap:6d} {r['n']:6d} {r['ung']:8.3f}"
        s+= f" {r['C90'][0]:10d}% {r['C90'][1]:8.3f} {r['C90'][2]:8.1f} |" if r["C90"] else \
            f" {'none':>11} {'--':>8} {'--':>8} |"
        for q in (20,10,7):
            if q in r["rows"]:
                p,y,_n = r["rows"][q]; s+=f" {p:.3f}/{y:6.1f}"
            else: s+=f" {'--':>13}"
        print(s)
    print(f"\n  T1b -- P_F CONTRACT (P_F = 1 - purity), booked:")
    print(f"  {'cap':>6} {'P_F=0.10 psr/sky':>17} {'P_F=0.01 psr/sky':>17}")
    for _cap in (40,400,1300):
        r=_B[_cap]; asr,w,sk=r["asr"],r["w"],r["sk"]
        tr,te = sk<30, sk>=30
        best={0.90:0.0, 0.99:0.0}
        for q in np.linspace(0.5,40,80):
            thr=float(np.percentile(asr[tr],q)); sel=te&(asr<=thr)
            if sel.sum()<5: continue
            pur=1-w[sel].mean(); y=(sel.sum()-w[sel].sum())/30
            for tgt in (0.90,0.99):
                if pur>=tgt and y>best[tgt]: best[tgt]=y
        print(f"  {_cap:6d} {best[0.90]:17.1f} {best[0.99]:17.1f}")
    print(f"\n  T1c -- VALIDATION: null PURITY (not yield), and convergence:")
    print(f"  {'cap':>6} {'signal ungated':>15} {'NULL ungated':>13} {'ratio':>7}"
          f" {'60->120 sky |d| psr/sky':>24}")
    for _cap in (40,400,1300):
        rn = corner32(_cap, nsky=30, ntrain=15, null=True)
        r120 = corner32(_cap, nsky=120, ntrain=60)
        d = abs((r120["C90"][2] if r120["C90"] else np.nan) - _B[_cap]["C90"][2])
        print(f"  {_cap:6d} {_B[_cap]['ung']:15.3f} {rn['ung']:13.3f}"
              f" {_B[_cap]['ung']/max(rn['ung'],1e-9):7.2f} {d:24.2f}")
print(f"\n{'='*74}\nJOB 32 -- T2: BOOK THE ONE-SHOT AT RAISED CAP\n{'='*74}")
def oneshot32(T=60, nsky=25, cap=40):
    globals()['_CAP32'] = cap
    ok=n=wrong=rnc=0; sv=[]; Tim=[]
    rg=np.random.default_rng(2700)
    for j in range(nsky):
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        Np=120
        Ti=np.clip(np.exp(rg.normal(np.log(T),0.5,Np)),5,90); Tim.append(np.median(Ti))
        beta=np.arcsin(rg.uniform(-1,1,Np)); Dk=rg.uniform(0.5,3.0,Np)
        px=2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/Ti); pb=3.0*(Ti/30.0)**(-2.5)
        pri=np.maximum(np.where(rg.random(Np)<0.7,np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(Np); rncl=Arn==0.0
        ks,sig,fk,rp,nats=inventory20(f2,h2,P0,P0,sea,float(np.mean(Ti)),100.0,14.0,Np)
        if len(ks)<2: continue
        rnc+=int(rncl.sum())
        for i in np.where(rncl)[0]:
            sgb=np.array([depro_sigma(
                np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2,0.0)),
                rp[a],fk[a],i,"loo") for a in range(len(ks))])
            sgb=sgb*np.sqrt(np.mean(Ti)/Ti[i])
            mu=np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            sgb=sgb*np.sqrt((1/3)/np.maximum((1+np.cos(mu))**2/4,1e-3))
            okm=sgb<np.pi
            if okm.sum()<2: continue
            sgu=sgb[okm]; sv.append(1.0/(1.0+float(np.median(sgu))**2))
            kap=2*np.pi*fk[okm]*(1-np.cos(mu[okm]))/C_PC_S
            Lt=rg.normal(0,pri[i]); Lg=np.linspace(-3*pri[i],3*pri[i],401)
            A=1.0/np.maximum(sgu,1e-6)
            z=A*np.exp(1j*(kap*Lt))+(rg.normal(0,1/np.sqrt(2),len(A))
                                     +1j*rg.normal(0,1/np.sqrt(2),len(A)))
            Sc=(2*A[:,None]*np.real(z[:,None]*np.exp(
                -1j*(kap[:,None]*Lg[None,:])))).sum(0)-Lg**2/(2*pri[i]**2)
            if 1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)),1e-30))<TOOTH_W:
                n+=1
                if abs(float(Lg[int(np.argmax(Sc))])-Lt)<=TOOTH_W: ok+=1
                else: wrong+=1
    return ok,wrong,n,rnc,float(np.median(sv)) if sv else np.nan, float(np.median(Tim))
print(f"  {'config':>20} {'med T_i':>8} {'s_med':>7} {'n':>6} {'right':>6}"
      f" {'purity':>7} {'right/RNc':>10}")
_R32={}
for _T in (45,60):
    for _cap in (40,1300):
        o,w_,n_,rc,sm,tm=oneshot32(T=_T,cap=_cap); _R32[(_T,_cap)]=(o,n_,rc)
        print(f"  {'T='+str(_T)+', cap '+str(_cap):>20} {tm:8.1f} {sm:7.4f} {n_:6d}"
              f" {o:6d} {o/max(n_,1):7.4f} {o/max(rc,1):10.4f}")
_a=_R32[(60,40)]; _b=_R32[(60,1300)]
print(f"    PRE-FLIGHT: median T_i tracks T (45 vs 60) -> span is live")
print(f"    T=60 purity {_a[0]/_a[1]:.3f} -> {_b[0]/_b[1]:.3f}"
      f"  [PRE-REG 0.68-0.74: {'HIT' if 0.68<=_b[0]/_b[1]<=0.74 else 'MISS'}]")
print(f"    T=60 right/RNc {_a[0]/_a[2]:.3f} -> {_b[0]/_b[2]:.3f}"
      f"  [PRE-REG 0.58-0.66: {'HIT' if 0.58<=_b[0]/_b[2]<=0.66 else 'MISS'}]")

print(f"\n{'='*74}\nJOB 32 -- T3: FREQUENCY-LADDER ASSOCIATION TEST\n{'='*74}")
print("  RULE (stated before running): for each (pulsar, monster) pair, predict")
print("  the retarded frequency band f_P(L) for L over the pulsar's prior support")
print("  [0, 3*sigma_d]; among all detectable pulsar-term lines at that pulsar")
print("  falling in that band, select the LOUDEST. No oracle assignment.")
_TSUN=4.925490947e-6; _MSUN=1.989e30
def t_coal_(f,ts): return (5.0/256.0)*ts**(-5.0/3.0)*(np.pi*f)**(-8.0/3.0)
_ordAMP=np.argsort(SRC_AMP); _AMPS=SRC_AMP[_ordAMP]
def mc_z_(amp):
    k=np.clip(np.searchsorted(_AMPS,amp),0,len(_AMPS)-1); i=_ordAMP[k]
    return SRC_MC[i],SRC_Z[i]
def ladder(nsky=25, scramble=False, Np=30):
    mis=0; tot=0; Krat=[]
    for j in range(nsky):
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        T=45.0*YR_S; nt=45.0*365.25/14.0
        rp=(h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100.0*1e-9)
        arr=rp*np.sqrt(0.2*120)
        band=(arr>=3.0)&(f2>=2e-9)&(f2<=3e-8)&(rp>=3.0)
        cand=np.where(band)[0]
        if len(cand)<5: continue
        sgc=sig_windowed(f2,P0,sea,rp,cand,T)
        solv=cand[sgc<np.pi]
        if len(solv)<5: continue
        amp=h2[solv]/np.maximum(f2[solv],1e-30)**(2.0/3.0)
        mc,zz=mc_z_(amp); ts=_TSUN*mc*(1.0+zz)/_MSUN
        mon=solv[mc>=np.quantile(mc,0.90)]                 # monsters: top Mc decile
        tsA=dict(zip(solv.tolist(), ts.tolist()))
        rg=np.random.default_rng(3200+j)
        L=rg.uniform(0.5,3.0,Np)*1000.0
        pri=np.maximum(rg.normal(2.0,0.6,Np),0.3)
        for p in range(Np):
            cmu={k: float(rg.uniform(-0.9,0.9)) for k in solv}
            fP={}; hP={}
            for k in solv:
                tau=L[p]*(1-cmu[k])/C_PC_S
                tc=t_coal_(max(f2[k],1e-30), tsA[k])
                fpk=f2[k]*(1.0+tau/max(tc,1e-9))**(-3.0/8.0)
                if fpk>1.0/T_OBS and fpk<1e-7:
                    fP[k]=fpk; hP[k]=h2[k]*(fpk/max(f2[k],1e-30))**(2.0/3.0)
            if len(fP)<3: continue
            for k in mon:
                if k not in fP: continue
                tot+=1
                tc=t_coal_(max(f2[k],1e-30), tsA[k])
                lo_,hi_=[],[]
                for Lx in (0.05, 3.0*pri[p]+L[p]*0+3.0):
                    tau=Lx*(1-cmu[k])/C_PC_S
                    v=f2[k]*(1.0+tau/max(tc,1e-9))**(-3.0/8.0); lo_.append(v)
                flo,fhi=min(lo_+[fP[k]])*0.999, max(lo_+[fP[k]])*1.001
                inb=[q for q in fP if flo<=fP[q]<=fhi]
                if not inb: mis+=1; continue
                if scramble:
                    pick=inb[int(rg.integers(len(inb)))]
                else:
                    pick=max(inb, key=lambda q: hP[q])
                if pick!=k: mis+=1; Krat.append(1.0); continue
                dLf=(8.0/3.0)*C_PC_S*(1.0/(2*T_OBS))*max(tc-L[p]*(1-cmu[k])/C_PC_S,1.0)/ \
                    (fP[k]*max(1-cmu[k],1e-6))
                dLt=C_PC_S/(fP[k]*max(1-cmu[k],1e-6))
                Kno=6*pri[p]/dLt; Kw=min(6*pri[p],2*dLf)/dLt
                Krat.append(Kno/max(Kw,1e-12))
    return mis/max(tot,1), np.median(Krat) if Krat else np.nan, tot
_m,_k,_t = ladder()
_ms,_ks,_ts_ = ladder(scramble=True)
print(f"\n  mis-identification rate: {_m:.3f}  (n={_t} monster-pulsar pairs)")
print(f"    PRE-REG 10-40%: {'HIT' if 0.10<=_m<=0.40 else 'MISS'}")
print(f"  retained K-shrink (median): /{_k:.2f}  vs oracle /2.50")
print(f"    PRE-REG /1.5-2.2: {'HIT' if 1.5<=_k<=2.2 else 'MISS'}")
print(f"  NULL scrambled assignment: mis-ID {_ms:.3f}, retained shrink /{_ks:.2f}")
print(f"    shrink must vanish under scrambling:"
      f" {'PASS' if _ks < 1.1 else 'FAIL -- shrink survives scrambling'}")
print(f"  VERDICT: {'ORACLE-ONLY, CLOSED (mis-ID > 60%)' if _m>0.60 else 'association is viable'}")
print(f"\n  T3b -- BOUNDED statistic (K itself), mis-ID folded in:")
def ladder_b(nsky=25, scramble=False, Np=30):
    KN=[]; KW=[]; mis=0; tot=0
    for j in range(nsky):
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        T=45.0*YR_S; nt=45.0*365.25/14.0
        rp=(h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100.0*1e-9)
        arr=rp*np.sqrt(0.2*120)
        cand=np.where((arr>=3.0)&(f2>=2e-9)&(f2<=3e-8)&(rp>=3.0))[0]
        if len(cand)<5: continue
        sgc=sig_windowed(f2,P0,sea,rp,cand,T); solv=cand[sgc<np.pi]
        if len(solv)<5: continue
        amp=h2[solv]/np.maximum(f2[solv],1e-30)**(2.0/3.0)
        mc,zz=mc_z_(amp); ts=_TSUN*mc*(1.0+zz)/_MSUN
        mon=solv[mc>=np.quantile(mc,0.90)]; tsA=dict(zip(solv.tolist(),ts.tolist()))
        rg=np.random.default_rng(3200+j)
        L=rg.uniform(0.5,3.0,Np)*1000.0; pri=np.maximum(rg.normal(2.0,0.6,Np),0.3)
        for p in range(Np):
            cmu={k: float(rg.uniform(-0.9,0.9)) for k in solv}
            fP={}; hP={}
            for k in solv:
                tau=L[p]*(1-cmu[k])/C_PC_S; tc=t_coal_(max(f2[k],1e-30),tsA[k])
                v=f2[k]*(1.0+tau/max(tc,1e-9))**(-3.0/8.0)
                if 1.0/T_OBS<v<1e-7: fP[k]=v; hP[k]=h2[k]*(v/max(f2[k],1e-30))**(2.0/3.0)
            if len(fP)<3: continue
            for k in mon:
                if k not in fP: continue
                tot+=1
                tc=t_coal_(max(f2[k],1e-30),tsA[k])
                dLt=C_PC_S/(fP[k]*max(1-cmu[k],1e-6))
                Kno=6*pri[p]/dLt; KN.append(Kno)
                lo=f2[k]*(1.0+(0.05*(1-cmu[k])/C_PC_S)/max(tc,1e-9))**(-3.0/8.0)
                flo,fhi=min(lo,fP[k])*0.999, max(lo,fP[k])*1.001
                inb=[q for q in fP if flo<=fP[q]<=fhi]
                if not inb: mis+=1; KW.append(Kno); continue
                pick=inb[int(rg.integers(len(inb)))] if scramble else max(inb,key=lambda q:hP[q])
                if pick!=k:
                    mis+=1; KW.append(Kno); continue      # mis-ID -> no shrink
                dLf=(8.0/3.0)*C_PC_S*(1.0/(2*T_OBS))*max(tc-L[p]*(1-cmu[k])/C_PC_S,1.0)/ \
                    (fP[k]*max(1-cmu[k],1e-6))
                KW.append(min(6*pri[p],2*dLf)/dLt)
    KN=np.array(KN); KW=np.array(KW)
    return mis/max(tot,1), np.median(KN), np.median(KW), (KN<=2).mean(), (KW<=2).mean(), tot
for _lab,_sc in (("real association",False), ("NULL scrambled",True)):
    m,kn,kw,f2b,f2a,tt = ladder_b(scramble=_sc)
    print(f"    {_lab:>18}: mis-ID {m:.3f} | med K {kn:6.2f} -> {kw:6.2f}"
          f" | K<=2 {f2b:.3f} -> {f2a:.3f} | n={tt}")
