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
import numpy as np, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT="/home/mattm/projects/HSYMT/wheel_hasasia_examples/"
print(f"\n{'='*74}\nJOB 35 T3(a) -- CONFORMAL TOOTH SETS (the open P_F question)\n{'='*74}")
def collect36(nsky=60, T=45, cap=1300, keep_curves=0):
    globals()['_CAP32']=cap
    rows=[]; curves=[]
    for j in range(nsky):
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        rng=np.random.default_rng(8000+sd+977*j); Np=120
        beta_=np.arcsin(rng.uniform(-1,1,Np)); Dk=rng.uniform(0.5,3.0,Np)
        px=2.34/np.cos(beta_)**2*Dk**2*10.0*np.sqrt(30.0/T); pb=3.0*(T/30.0)**(-2.5)
        pri=np.maximum(np.where(rng.random(Np)<0.7,np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(Np); rncl=Arn==0.0
        ks,sig,fk,rp,nats=inventory20(f2,h2,P0,P0,sea,T,100.0,14.0,Np)
        if len(ks)<2: continue
        for i in np.where(rncl)[0]:
            sg=np.array([depro_sigma(
                np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2,0.0)),
                rp[a],fk[a],i,"loo") for a in range(len(ks))])
            ok=sg<np.pi
            if ok.sum()<2: continue
            sgu=sg[ok]; mu=np.arccos(rng.uniform(-0.9,0.9,len(sgu)))
            sgu=sgu*np.sqrt((1/3)/np.maximum((1+np.cos(mu))**2/4,1e-3))
            ok2=sgu<np.pi
            if ok2.sum()<2: continue
            sgu=sgu[ok2]; kap=2*np.pi*fk[ok][ok2]*(1-np.cos(mu[ok2]))/C_PC_S
            Lt=float(np.clip(rng.normal(0,pri[i]),-2.9*pri[i],2.9*pri[i]))
            Lg=np.linspace(-3*pri[i],3*pri[i],601)
            eps=kap*Lt+rng.normal(0,sgu)
            per=(kap[:,None]*Lg[None,:]-eps[:,None]+np.pi)%(2*np.pi)-np.pi
            Sc=-(per**2/(2*sgu[:,None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
            if 1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)),1e-30))>=TOOTH_W: continue
            w=np.exp(Sc-Sc.max()); w/=w.sum()
            # --- teeth: group grid into TOOTH_W-wide cells, mass per cell
            cell=np.floor((Lg+3*pri[i])/(2*TOOTH_W)).astype(int)
            nc=cell.max()+1
            m=np.bincount(cell,weights=w,minlength=nc)
            ctr=np.bincount(cell,weights=Lg*w,minlength=nc)/np.maximum(m,1e-300)
            tcell=int(np.floor((Lt+3*pri[i])/(2*TOOTH_W)))
            order=np.argsort(m)[::-1]
            cum=np.cumsum(m[order])
            pos=int(np.where(order==tcell)[0][0]) if tcell<nc and (order==tcell).any() else len(order)-1
            aps=float(cum[pos])                       # APS score for the TRUE tooth
            ip=int(np.argmax(Sc))
            rows.append(dict(sky=j,aps=aps,m=m,order=order,cum=cum,ctr=ctr,
                             tcell=tcell,wrong=abs(float(Lg[ip])-Lt)>TOOTH_W,
                             al=alias_score24(dict(pri=pri[i],Ltrue=Lt,kap=kap,sig=sgu)),
                             top=float(m.max())))
            if len(curves)<keep_curves:
                curves.append(dict(Lg=Lg,Sc=Sc,Lt=Lt,pri=pri[i],kap=kap,sgu=sgu,eps=eps,
                                   K=6*pri[i]/max(2*np.pi/max(kap.min(),1e-30),1e-9)))
    return rows,curves
_R36,_CUR=collect36(nsky=60, keep_curves=400)
print(f"  n={len(_R36)}; ungated purity {1-np.mean([r['wrong'] for r in _R36]):.3f}")
sk=np.array([r["sky"] for r in _R36]); cal=sk<30; tst=sk>=30
alv=np.array([r["al"] for r in _R36])
dec=np.digitize(alv,np.quantile(alv[cal],[.2,.4,.6,.8]))   # Mondrian: gate quintile
print(f"\n  {'alpha':>6} {'Mondrian':>9} {'coverage':>9} {'nominal':>8} {'med |set|':>10}"
      f" {'singleton frac':>15} {'singleton purity':>17}")
_CONF=[]
for _al in (0.10,0.01):
    for _mond in (False,True):
        cov=[]; sz=[]; sing=[]; singok=[]
        for d in (range(5) if _mond else [None]):
            c=cal&(dec==d) if _mond else cal
            t=tst&(dec==d) if _mond else tst
            if c.sum()<20 or t.sum()<10: continue
            s=np.array([_R36[i]["aps"] for i in np.where(c)[0]])
            n=len(s); q=np.quantile(s,min(np.ceil((n+1)*(1-_al))/n,1.0))
            for i in np.where(t)[0]:
                r=_R36[i]
                k=int(np.searchsorted(r["cum"],q))+1
                k=min(k,len(r["order"]))
                inset=r["tcell"] in set(r["order"][:k].tolist())
                cov.append(inset); sz.append(k); sing.append(k==1)
                if k==1: singok.append(inset)
        cov=np.array(cov); sz=np.array(sz); sing=np.array(sing)
        se=np.sqrt(max(cov.mean()*(1-cov.mean()),1e-9)/max(len(cov),1))
        okc=abs(cov.mean()-(1-_al))<=max(2*se,0.02)
        sp=np.mean(singok) if singok else np.nan
        print(f"  {_al:6.2f} {str(_mond):>9} {cov.mean():9.4f} {1-_al:8.2f}"
              f" {np.median(sz):10.1f} {sing.mean():15.4f} {sp:17.4f}"
              f"  {'OK' if okc else 'MISCOVERAGE'}")
        _CONF.append((_al,_mond,cov.mean(),float(np.median(sz)),sing.mean(),sp,
                      float(sing.sum())/30.0))
print(f"\n  THE QUESTION: does a conformal singleton set deliver P_F<=0.01?")
for (_al,_m,c,s,sf,sp,ysky) in _CONF:
    if _m and _al==0.01:
        print(f"    alpha=0.01 Mondrian: singleton purity {sp:.4f}"
              f" -> P_F = {1-sp:.4f}; yield {ysky:.1f} psr/sky")
        print(f"    {'REACHES P_F<=0.01' if (1-sp)<=0.01 else 'DOES NOT reach P_F<=0.01'}")
with open(OUT+"fig9_conformal.csv","w",newline="") as fh:
    cw=csv.writer(fh); cw.writerow(["alpha","mondrian","coverage","median_set_size",
                                    "singleton_frac","singleton_purity","singleton_psr_per_sky"])
    for r in _CONF: cw.writerow([f"{x:.4f}" if isinstance(x,float) else x for x in r])
print(f"  wrote fig9_conformal.csv")

# ============ FIG 1: the comb, stacking, and the 2-D pair slice ==============
print(f"\n{'='*74}\nJOB 36 -- FIG 1: THE FRINGE COMB\n{'='*74}")
cand=[c for c in _CUR if 8<=c["K"]<=12 and len(c["kap"])>=4]
c1=cand[0] if cand else _CUR[0]
print(f"  representative pulsar: K={c1['K']:.1f}, {len(c1['kap'])} combs,"
      f" sigma_d={c1['pri']:.2f} pc, L_true={c1['Lt']:.3f} pc")
Lg=c1["Lg"]; kap=c1["kap"]; sgu=c1["sgu"]; eps=c1["eps"]; pri=c1["pri"]
def stack(n):
    per=(kap[:n,None]*Lg[None,:]-eps[:n,None]+np.pi)%(2*np.pi)-np.pi
    return -(per**2/(2*sgu[:n,None]**2)).sum(0)-Lg**2/(2*pri**2)
fig,ax=plt.subplots(1,3,figsize=(13,3.9))
S1=stack(1); ax[0].plot(Lg,S1-S1.max(),color="C0",lw=1.1)
gp=np.exp(-Lg**2/(2*pri**2)); ax[0].fill_between(Lg,-60,-60+18*gp,color="0.85",lw=0)
ax[0].axvline(c1["Lt"],color="0.25",lw=.9)
ax[0].set_ylim(-60,3); ax[0].set_xlabel("pulsar distance offset (pc)")
ax[0].set_ylabel("log likelihood"); ax[0].text(.02,.95,"one source",transform=ax[0].transAxes)
ns=[1,2,4,len(kap)]; off=0
for a,n in enumerate(ns):
    S=stack(n); S=S-S.max()
    ax[1].plot(Lg,S+off,color="C0",lw=1.0,alpha=.45+.55*a/3)
    ax[1].text(Lg[0]+.02*(Lg[-1]-Lg[0]),off+2,f"{'all' if n==len(kap) else n} source"+("s" if n>1 else ""),
               fontsize=8,color="0.25")
    off+=45
ax[1].set_ylim(-55,off+10); ax[1].set_xlabel("pulsar distance offset (pc)")
ax[1].set_yticks([]); ax[1].axvline(c1["Lt"],color="0.25",lw=.9)
c2=[c for c in _CUR if c is not c1 and len(c["kap"])>=2][0]
Lq=np.linspace(-3*c2["pri"],3*c2["pri"],161); Lp=np.linspace(Lg[0],Lg[-1],161)
k_sh=min(len(kap),len(c2["kap"]))-1
A=np.zeros((len(Lq),len(Lp)))
for a in range(len(Lq)):
    pp=(kap[k_sh]*Lp-eps[k_sh]+np.pi)%(2*np.pi)-np.pi
    qq=(c2["kap"][k_sh]*Lq[a]-c2["eps"][k_sh]+np.pi)%(2*np.pi)-np.pi
    A[a]=-(pp**2/(2*sgu[k_sh]**2))-(qq**2/(2*c2["sgu"][k_sh]**2))
im=ax[2].imshow(A-A.max(),origin="lower",aspect="auto",cmap="Blues_r",
                extent=[Lp[0],Lp[-1],Lq[0],Lq[-1]],vmin=-25,vmax=0)
ax[2].plot(c1["Lt"],c2["Lt"],'+',color="0.15",ms=11,mew=1.6)
ax[2].set_xlabel("pulsar p offset (pc)"); ax[2].set_ylabel("pulsar q offset (pc)")
for a in ax: a.grid(alpha=.25,lw=.5)
plt.tight_layout(); plt.savefig(OUT+"fig1_comb.png",dpi=150); plt.close()
np.savez(OUT+"fig1_comb.npz",Lg=Lg,S1=stack(1),S2=stack(2),S4=stack(4),
         Sall=stack(len(kap)),Ltrue=c1["Lt"],pri=pri,K=c1["K"],
         joint=A,Lp=Lp,Lq=Lq,Lt_q=c2["Lt"])
print(f"  wrote fig1_comb.png + fig1_comb.npz")

# ============ FIG 5(a,b): gate distributions + ROC ==========================
print(f"\n{'='*74}\nJOB 36 -- FIG 5(a,b): THE GATE\n{'='*74}")
wr=np.array([r["wrong"] for r in _R36],bool)
fig,ax=plt.subplots(1,2,figsize=(9.5,3.9))
bins=np.linspace(0,1,41)
ax[0].hist(alv[~wr],bins=bins,histtype="step",color="C0",lw=1.4,density=True)
ax[0].hist(alv[wr],bins=bins,histtype="step",color="C3",lw=1.4,density=True)
thr20=np.quantile(alv[cal],0.20); ax[0].axvline(thr20,color="0.3",ls="--",lw=1)
ax[0].text(.05,.9,"right",color="C0",transform=ax[0].transAxes)
ax[0].text(.05,.8,"wrong",color="C3",transform=ax[0].transAxes)
ax[0].set_xlabel("alias-coincidence score"); ax[0].set_ylabel("density")
def roc(sc,w):
    o=np.argsort(sc); w2=w[o]
    tp=np.cumsum(~w2)/max((~w2).sum(),1); fp=np.cumsum(w2)/max(w2.sum(),1)
    return fp,tp
for m,lab,c in ((cal,"train","0.55"),(tst,"test","C0")):
    fp,tp=roc(alv[m],wr[m]); ax[1].plot(fp,tp,color=c,lw=1.4)
    ax[1].text(.55,.25 if lab=="train" else .18,lab,color=c,transform=ax[1].transAxes)
for q,mk in ((20,'o'),(10,'s'),(7,'^')):
    t=np.quantile(alv[cal],q/100); m=tst&(alv<=t)
    ax[1].plot((wr[m].sum())/max(wr[tst].sum(),1),((~wr[m]).sum())/max((~wr[tst]).sum(),1),
               mk,color="0.2",ms=6)
ax[1].plot([0,1],[0,1],ls=":",color="0.7",lw=1)
ax[1].set_xlabel("false-positive rate"); ax[1].set_ylabel("true-positive rate")
for a in ax: a.grid(alpha=.25,lw=.5)
plt.tight_layout(); plt.savefig(OUT+"fig5_gate.png",dpi=150); plt.close()
with open(OUT+"fig5_gate.csv","w",newline="") as fh:
    cw=csv.writer(fh); cw.writerow(["alias_score","wrong","split"])
    for i in range(len(alv)):
        cw.writerow([f"{alv[i]:.6f}",int(wr[i]),"cal" if cal[i] else "test"])
print(f"  wrote fig5_gate.png + fig5_gate.csv (n={len(alv)})")
print(f"    test AUC {auc_of(alv[tst],wr[tst]):.3f}; train AUC {auc_of(alv[cal],wr[cal]):.3f}")
import numpy as np, csv
print(f"\n{'='*74}\nJOB 43B -- T1: DIAGNOSING THE COLLAPSE FAILURE AT A = 5.0e-15\n{'='*74}")
T45=45.0*YR_S
def inv_rate(f2,h2,P0,sea,T_yr,cap=1300,rate=1.0,nmax=12):
    """inv_cap with COUNTS-ONLY scaling: res -> rate*res, Pt untouched."""
    T=T_yr*YR_S; nt=T_yr*365.25/14.0
    rp=(h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100.0*1e-9)
    arr=rp*np.sqrt(0.2*120)
    cand=np.where((arr>=3.0)&(f2>=2e-9)&(f2<=3e-8)&(rp>=3.0))[0]
    if len(cand)>cap: cand=cand[np.argsort(rp[cand])[::-1][:cap]]
    if len(cand)<2: return None
    o=np.argsort(f2); fs=f2[o]; Pb=P0[o]; half=100.0/T
    sg=np.empty(len(cand))
    for a,k in enumerate(cand):
        lo=np.searchsorted(fs,f2[k]-half); hi=np.searchsorted(fs,f2[k]+half)
        w=np.sinc((fs[lo:hi]-f2[k])*T)**2
        c=w*Pb[lo:hi]*(T/T_OBS)
        b=int(np.clip(np.digitize(f2[k],bin_edges)-1,0,len(bin_cent)-1))
        res=(c.sum()-P0[k]*(T/T_OBS)+sea[b]*(T_OBS/T))*rate     # COUNTS ONLY
        Pt=P0[k]*(T/T_OBS)
        sg[a]=np.sqrt(max(res,0)/max(Pt,1e-300)+1/max(rp[k],1e-9)**2)
    good=sg<np.pi; cand,sg=cand[good],sg[good]
    if len(cand)<2: return None
    nats=np.log(np.pi/sg); ordr=np.argsort(nats)[::-1][:nmax]
    return cand[ordr],sg[ordr],f2[cand[ordr]],rp[cand[ordr]],arr[cand[ordr]]
def collect_banked(rate=1.0,nsky=200,nps=30,T_yr=45.0,full_dmu=True,marg=False,seed0=9100):
    H=[];AL=[];SK=[];KK=[]
    rgA=np.random.default_rng(555)
    for j in range(nsky):
        rt=rate
        if marg:
            A=6.4e-15*np.exp(rgA.normal(0,0.32)); rt=float(np.clip((A/6.4e-15)**2,0.2,3.5))
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        r=inv_rate(f2,h2,P0,sea,T_yr,rate=rt)
        if r is None: continue
        ks,sg0,fk,rpk,arrk=r
        rg=np.random.default_rng(seed0+j)
        Ti=np.clip(np.exp(rg.normal(np.log(T_yr),0.35,120)),5,90)
        beta=np.arcsin(rg.uniform(-1,1,120)); Dk=rg.uniform(0.5,3.0,120)
        px=2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/Ti); pb=3.0*(Ti/30.0)**(-2.5)
        pri=np.maximum(np.where(rg.random(120)<0.7,np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(120); rncl=Arn==0.0
        for i in np.where(rncl)[0][:nps]:
            sg=np.array([depro_sigma(np.sqrt(max(sg0[a]**2-1/max(rpk[a],1e-9)**2,0.0)),
                         rpk[a],fk[a],i,"loo") for a in range(len(ks))])
            sg=sg*np.sqrt(np.mean(Ti)/Ti[i])
            mu=np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            sg=sg*np.sqrt((1/3)/np.maximum((1+np.cos(mu))**2/4,1e-3))
            kap=2*np.pi*fk*(1-np.cos(mu))/C_PC_S
            md=1.0/np.maximum(arrk,1e-6)
            dmu=rg.rayleigh(md*np.sqrt(2/np.pi)) if full_dmu else md
            dfq=1.0/(T_yr*YR_S*np.maximum(arrk,1e-6))
            dkap=np.sqrt((2*np.pi*(1-np.cos(mu))/C_PC_S)**2*dfq**2
                         +(2*np.pi*fk*np.sin(mu)/C_PC_S)**2*dmu**2)
            sg=np.sqrt(sg**2+(dkap*pri[i])**2)
            ok=sg<np.pi
            if ok.sum()<2: continue
            sgu,kapu=sg[ok],kap[ok]
            Lt=float(np.clip(rg.normal(0,pri[i]),-2.9*pri[i],2.9*pri[i]))
            Lg=np.linspace(-3*pri[i],3*pri[i],601)
            eps=kapu*Lt+rg.normal(0,sgu)
            per=(kapu[:,None]*Lg[None,:]-eps[:,None]+np.pi)%(2*np.pi)-np.pi
            S=-(per**2/(2*sgu[:,None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
            if 1.0/np.sqrt(max(float(np.sum(kapu**2/sgu**2)),1e-30))>=TOOTH_W: continue
            H.append(abs(float(Lg[int(np.argmax(S))])-Lt)<=TOOTH_W)
            AL.append(alias_score24(dict(pri=pri[i],Ltrue=Lt,kap=kapu,sig=sgu)))
            SK.append(j); KK.append(6*pri[i]/max(2*np.pi/max(kapu.min(),1e-30),1e-9))
    return np.array(H),np.array(AL),np.array(SK),np.array(KK)
def r2f(H,KK):
    x=np.log10(np.maximum(KK,1.01))
    if x.std()<1e-6: x=x+np.random.default_rng(0).normal(0,1e-3,len(x))
    xs=[];ys=[]
    for nb in (6,5,4):
        e=np.unique(np.quantile(x,np.linspace(0,1,nb+1)))
        if len(e)<5: continue
        acc=[]
        for a in range(len(e)-1):
            acc.append((x>=e[a])&((x<e[a+1]) if a<len(e)-2 else (x<=e[a+1])))
        merged=[];cur=None
        for m in acc:
            cur=m if cur is None else (cur|m)
            if cur.sum()>=12: merged.append(cur); cur=None
        if cur is not None and merged: merged[-1]=merged[-1]|cur
        xs=[x[m].mean() for m in merged]; ys=[H[m].mean() for m in merged]
        if len(xs)>=4: break
    if len(xs)<3: return np.nan
    xs=np.array(xs);ys=np.array(ys); p=np.polyfit(xs,ys,1)
    den=((ys-ys.mean())**2).sum()
    return float(1-((ys-np.polyval(p,xs))**2).sum()/den) if den>1e-12 else np.nan
def corner_oos(H,AL,SK,ntr=100):
    tr=SK<ntr; te=SK>=ntr; nte=max(len(np.unique(SK[te])),1)
    qs=np.linspace(0.5,60,120); k=np.ones(7)/7.0
    th=[np.quantile(AL[tr],q/100) for q in qs]
    pu=[]
    for t in th:
        s=tr&(AL<=t); pu.append(H[s].mean() if s.sum()>=10 else np.nan)
    pu=np.array(pu)
    ps=np.convolve(np.nan_to_num(pu,nan=np.nanmean(pu[np.isfinite(pu)])),k,mode="same")
    thr=np.nan
    for a in range(len(ps)-1):
        if ps[a]>=0.90>=ps[a+1]:
            f=(ps[a]-0.90)/max(ps[a]-ps[a+1],1e-9); thr=th[a]+f*(th[a+1]-th[a]); break
    if not np.isfinite(thr): thr=np.quantile(AL[tr],0.20)
    s=te&(AL<=thr)
    per={j:(int((s&(SK==j)&H).sum()),int((s&(SK==j)).sum())) for j in np.unique(SK[te])}
    rg=np.random.default_rng(17); bt=[]
    u=np.unique(SK[te])
    for b in range(300):
        bs=rg.choice(u,size=len(u),replace=True)
        bt.append(sum(per[j][0] for j in bs)/len(u))
    return float(H[s].mean()), float(H[s].sum()/nte), np.quantile(bt,[.1,.9]), int(s.sum())
print(f"  (a) does R^2 recover with statistics on the BANKED harness at rate 0.60?")
for lab,ns in (("banked, 60 skies",60),("banked, 200 skies",200)):
    H,AL,SK,KK=collect_banked(rate=0.60,nsky=ns)
    print(f"    {lab:>20}: n={len(H):5d}  law R^2 = {r2f(H,KK):.4f}")
print(f"  (b) banked-harness R^2 across the CI grid (200 skies each):")
print(f"    {'A (1e-15)':>10} {'rate':>6} {'n':>6} {'law R^2':>9}")
for rate,A in ((0.33,3.7),(0.60,5.0),(1.0,6.4),(1.7,8.3),(2.7,10.6)):
    H,AL,SK,KK=collect_banked(rate=rate,nsky=200)
    print(f"    {A:10.1f} {rate:6.2f} {len(H):6d} {r2f(H,KK):9.4f}")
import numpy as np
print(f"\n{'='*74}\nJOB 43C T1 -- ADJUDICATING THE FUZZY-TEMPLATE TERM\n{'='*74}")
print("  WHAT IT MODELS: the source's sky-position and frequency uncertainty.")
print("    dmu ~ 1/arr (rad), df ~ 1/(T*arr) -> dkappa via kappa = 2*pi*f(1-cos mu)/c")
print("    -> extra template phase error dkappa * L.")
print("  WHICH TERM: it enters SIGMA (the decoder's template noise), i.e. the")
print("    sigma_phi term of the criterion. It does NOT touch res/Pt.")
print("  IS IT DOUBLE-COUNTED? inventory20's sigma = confusion(res/Pt) + 1/rho^2.")
print("    The 1/rho^2 term prices the source's AMPLITUDE/phase SNR. Sky position")
print("    and frequency uncertainty appear NOWHERE in inventory20. So the fuzzy")
print("    term adds an effect that is genuinely absent, not a second copy.")
print("  CONTEXT: Jobs 28 T6 / 33 T1 both flagged that comb frequencies are treated")
print("    as ORACLE. The fuzzy term is exactly the honest price of that oracle.")
def collect36f(nsky=200, T=45, cap=1300, fuzzy=False, full_dmu=False, seed0=8000):
    globals()['_CAP32']=cap
    out=[]
    for j in range(nsky):
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        rng=np.random.default_rng(seed0+sd+977*j); Np=120
        beta=np.arcsin(rng.uniform(-1,1,Np)); Dk=rng.uniform(0.5,3.0,Np)
        px=2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/T); pb=3.0*(T/30.0)**(-2.5)
        pri=np.maximum(np.where(rng.random(Np)<0.7,np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(Np); rncl=Arn==0.0
        ks,sig,fk,rp,nats=inventory20(f2,h2,P0,P0,sea,T,100.0,14.0,Np)
        if len(ks)<2: continue
        arr=rp*np.sqrt(0.2*Np)
        for i in range(Np):
            if not rncl[i]: continue
            sg=np.array([depro_sigma(np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2,0.0)),
                         rp[a],fk[a],i,"loo") for a in range(len(ks))])
            ok=sg<np.pi
            if ok.sum()<2: continue
            sgu=sg[ok]; fku=fk[ok]; arru=arr[ok]
            mu=np.arccos(rng.uniform(-0.9,0.9,len(sgu)))
            kap=2*np.pi*fku*(1-np.cos(mu))/C_PC_S
            if fuzzy:
                md=1.0/np.maximum(arru,1e-6)
                dmu=rng.rayleigh(md*np.sqrt(2/np.pi)) if full_dmu else md
                dfq=1.0/(T*YR_S*np.maximum(arru,1e-6))
                dkap=np.sqrt((2*np.pi*(1-np.cos(mu))/C_PC_S)**2*dfq**2
                             +(2*np.pi*fku*np.sin(mu)/C_PC_S)**2*dmu**2)
                sgu=np.sqrt(sgu**2+(dkap*pri[i])**2)
                ok2=sgu<np.pi
                if ok2.sum()<2: continue
                sgu=sgu[ok2]; kap=kap[ok2]
            Ltrue=float(np.clip(rng.normal(0,pri[i]),-2.9*pri[i],2.9*pri[i]))
            Lg=np.linspace(-3*pri[i],3*pri[i],601)
            eps=kap*Ltrue+rng.normal(0,sgu)
            r_=(kap[:,None]*Lg[None,:]-eps[:,None]+np.pi)%(2*np.pi)-np.pi
            S=-(r_**2/(2*sgu[:,None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
            if 1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)),1e-30))>=TOOTH_W: continue
            Lh=float(Lg[int(np.argmax(S))])
            out.append(dict(sky=j,wrong=bool(abs(Lh-Ltrue)>TOOTH_W),
                            sig=sgu.copy(),kap=kap.copy(),pri=float(pri[i]),
                            Ltrue=Ltrue,npsr=int(rncl.sum())))
    return out
def book(rec,ntr=100,label=""):
    sk=np.array([r["sky"] for r in rec]); wr=np.array([r["wrong"] for r in rec],bool)
    al=np.array([alias_score24(r) for r in rec])
    cal=sk<ntr; tst=sk>=ntr; u=np.unique(sk[tst]); nte=len(u)
    qs=np.linspace(0.5,60,120); k=np.ones(7)/7.0
    th=[np.quantile(al[cal],q/100) for q in qs]
    pu=[]
    for t in th:
        s=cal&(al<=t); pu.append(1-wr[s].mean() if s.sum()>=10 else np.nan)
    pu=np.array(pu)
    ps=np.convolve(np.nan_to_num(pu,nan=np.nanmean(pu[np.isfinite(pu)])),k,mode="same")
    thr=np.nan; gate=np.nan
    for a in range(len(ps)-1):
        if ps[a]>=0.90>=ps[a+1]:
            f=(ps[a]-0.90)/max(ps[a]-ps[a+1],1e-9)
            thr=th[a]+f*(th[a+1]-th[a]); gate=qs[a]+f*(qs[a+1]-qs[a]); break
    if not np.isfinite(thr):
        return dict(ok=False,gate=np.nan,pur=np.nan,y=np.nan,
                    maxpur=float(np.nanmax(ps)))
    s=tst&(al<=thr)
    per={j:(int((s&(sk==j)&~wr).sum()),int((s&(sk==j)).sum())) for j in u}
    rg=np.random.default_rng(23); by=[];bp=[]
    for b in range(400):
        bs=rg.choice(u,size=nte,replace=True)
        r_=sum(per[j][0] for j in bs); n_=sum(per[j][1] for j in bs)
        by.append(r_/nte)
        if n_>0: bp.append(r_/n_)
    return dict(ok=True,gate=gate,pur=1-wr[s].mean(),y=(s.sum()-wr[s].sum())/nte,
                yb=np.quantile(by,[.1,.9]),pb=np.quantile(bp,[.1,.9]),n=int(s.sum()))
print(f"\n  TOGGLE TEST (anchor config, ONLY the fuzzy term changes):")
print(f"  {'arm':>28} {'gate %':>7} {'purity':>8} {'psr/sky':>9} {'10/90 yield':>16}")
A=book(collect36f(fuzzy=False),label="no fuzzy")
B=book(collect36f(fuzzy=True,full_dmu=False),label="fuzzy, mean dmu")
for lab,r in (("ANCHOR: no fuzzy",A),("+ fuzzy (mean dmu)",B)):
    if r["ok"]:
        print(f"  {lab:>28} {r['gate']:7.1f} {r['pur']:8.4f} {r['y']:9.2f}"
              f" [{r['yb'][0]:6.2f},{r['yb'][1]:6.2f}]")
    else:
        print(f"  {lab:>28} {'--':>7} {'no 0.90 crossing':>8}"
              f"  max train purity {r['maxpur']:.4f}")
if A["ok"] and B["ok"]:
    print(f"\n    old->new: purity {A['pur']:.4f} -> {B['pur']:.4f} ({B['pur']-A['pur']:+.4f});"
          f" yield {A['y']:.2f} -> {B['y']:.2f} (x{B['y']/A['y']:.3f})")
    print(f"    does fuzzy explain 43B's 0.883-0.893 purity floor?"
          f" {'YES -- contract still clears here' if B['pur']>=0.90 else 'PARTLY'}")
elif A["ok"] and not B["ok"]:
    print(f"\n    WITH FUZZY THERE IS NO 0.90 CROSSING (max train purity"
          f" {B['maxpur']:.4f}) -- this IS the 43B floor, now reproduced on the")
    print(f"    anchor configuration with everything else frozen.")
import numpy as np, csv
print(f"\n{'='*74}\nJOB 43C T2 -- STEPWISE RE-BOOKING (fuzzy ON = adjudicated baseline)\n{'='*74}")
def collect36m(nsky=200,T=45,fuzzy=True,full_dmu=False,marg=False,seed0=8000):
    globals()['_CAP32']=1300
    out=[]; rgA=np.random.default_rng(555)
    for j in range(nsky):
        rt=1.0
        if marg:
            A=6.4e-15*np.exp(rgA.normal(0,0.32)); rt=float(np.clip((A/6.4e-15)**2,0.2,3.5))
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        rng=np.random.default_rng(seed0+sd+977*j); Np=120
        beta=np.arcsin(rng.uniform(-1,1,Np)); Dk=rng.uniform(0.5,3.0,Np)
        px=2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/T); pb=3.0*(T/30.0)**(-2.5)
        pri=np.maximum(np.where(rng.random(Np)<0.7,np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(Np); rncl=Arn==0.0
        r=inv_rate(f2,h2,P0,sea,T,rate=rt)
        if r is None: continue
        ks,sig,fk,rp,arr=r
        for i in range(Np):
            if not rncl[i]: continue
            sg=np.array([depro_sigma(np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2,0.0)),
                         rp[a],fk[a],i,"loo") for a in range(len(ks))])
            ok=sg<np.pi
            if ok.sum()<2: continue
            sgu=sg[ok]; fku=fk[ok]; arru=arr[ok]
            mu=np.arccos(rng.uniform(-0.9,0.9,len(sgu)))
            kap=2*np.pi*fku*(1-np.cos(mu))/C_PC_S
            if fuzzy:
                md=1.0/np.maximum(arru,1e-6)
                dmu=rng.rayleigh(md*np.sqrt(2/np.pi)) if full_dmu else md
                dfq=1.0/(T*YR_S*np.maximum(arru,1e-6))
                dkap=np.sqrt((2*np.pi*(1-np.cos(mu))/C_PC_S)**2*dfq**2
                             +(2*np.pi*fku*np.sin(mu)/C_PC_S)**2*dmu**2)
                sgu=np.sqrt(sgu**2+(dkap*pri[i])**2)
                o2=sgu<np.pi
                if o2.sum()<2: continue
                sgu=sgu[o2]; kap=kap[o2]
            Lt=float(np.clip(rng.normal(0,pri[i]),-2.9*pri[i],2.9*pri[i]))
            Lg=np.linspace(-3*pri[i],3*pri[i],601)
            eps=kap*Lt+rng.normal(0,sgu)
            r_=(kap[:,None]*Lg[None,:]-eps[:,None]+np.pi)%(2*np.pi)-np.pi
            S=-(r_**2/(2*sgu[:,None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
            if 1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)),1e-30))>=TOOTH_W: continue
            Lh=float(Lg[int(np.argmax(S))])
            out.append(dict(sky=j,wrong=bool(abs(Lh-Lt)>TOOTH_W),sig=sgu.copy(),
                            kap=kap.copy(),pri=float(pri[i]),Ltrue=Lt))
    return out
print(f"  {'step':>40} {'gate %':>7} {'purity':>8} {'psr/sky':>9} {'10/90':>16} {'contract':>9}")
STEPS=[("1. BASELINE (fuzzy, mean dmu)",dict(fuzzy=True,full_dmu=False)),
       ("2. + full dmu (F8)",dict(fuzzy=True,full_dmu=True)),
       ("3. + marginalised A (F7) = FINAL",dict(fuzzy=True,full_dmu=True,marg=True))]
R={}
for lab,kw in STEPS:
    r=book(collect36m(**kw)); R[lab]=r
    if r["ok"]:
        print(f"  {lab:>40} {r['gate']:7.1f} {r['pur']:8.4f} {r['y']:9.2f}"
              f" [{r['yb'][0]:6.2f},{r['yb'][1]:6.2f}]"
              f" {'HOLDS' if r['pur']>=0.90 else 'below':>9}")
    else:
        print(f"  {lab:>40} {'--':>7} no 0.90 crossing (max train {r['maxpur']:.4f})")
b=R[STEPS[0][0]]; d=R[STEPS[1][0]]; f=R[STEPS[2][0]]
if b["ok"] and d["ok"]:
    print(f"\n    step 2 (full dmu): purity {b['pur']:.4f}->{d['pur']:.4f}"
          f" ({d['pur']-b['pur']:+.4f}), yield x{d['y']/b['y']:.3f}"
          f"  [PRE-REG purity +0.02..0.05: {'HIT' if 0.02<=d['pur']-b['pur']<=0.05 else 'MISS'};"
          f" yield +2..5%: {'HIT' if 1.02<=d['y']/b['y']<=1.05 else 'MISS'}]")
if d["ok"] and f["ok"]:
    print(f"    step 3 (marginalised): yield x{f['y']/d['y']:.3f},"
          f" purity {f['pur']-d['pur']:+.4f}")
    print(f"    FINAL: {f['y']:.2f} [{f['yb'][0]:.2f},{f['yb'][1]:.2f}] psr/sky"
          f" @ purity {f['pur']:.4f} [{f['pb'][0]:.3f},{f['pb'][1]:.3f}], gate {f['gate']:.1f}%")
    print(f"    contract >=0.90 out-of-sample:"
          f" {'HOLDS -- this row books as FINAL' if f['pur']>=0.90 else 'BELOW -- corner lives at gate %.1f%% with purity %.4f'%(f['gate'],f['pur'])}")
print(f"\n  T3 -- CORNER(A) ON THE ADJUDICATED CONFIGURATION")
print(f"  {'A (1e-15)':>10} {'rate':>6} {'gate %':>7} {'purity':>8} {'psr/sky':>9}")
rows=[]
for rate,A in ((0.33,3.7),(0.60,5.0),(1.0,6.4),(1.7,8.3),(2.7,10.6)):
    rec=[]
    globals()['_CAP32']=1300
    # reuse collect36m with a fixed rate by temporarily overriding marg draw
    def _fixed(nsky=200,rate=rate):
        out=[];
        for j in range(nsky):
            sd,f2,h2,P0,sea=_pop18(j%50,0.35)
            rng=np.random.default_rng(8000+sd+977*j); Np=120
            beta=np.arcsin(rng.uniform(-1,1,Np)); Dk=rng.uniform(0.5,3.0,Np)
            px=2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/45); pb=3.0*(45/30.0)**(-2.5)
            pri=np.maximum(np.where(rng.random(Np)<0.7,np.minimum(px,pb),px),0.05)
            Arn,_g=rn_draw(Np); rncl=Arn==0.0
            r=inv_rate(f2,h2,P0,sea,45.0,rate=rate)
            if r is None: continue
            ks,sig,fk,rp,arr=r
            for i in range(Np):
                if not rncl[i]: continue
                sg=np.array([depro_sigma(np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2,0.0)),
                             rp[a],fk[a],i,"loo") for a in range(len(ks))])
                ok=sg<np.pi
                if ok.sum()<2: continue
                sgu=sg[ok]; fku=fk[ok]; arru=arr[ok]
                mu=np.arccos(rng.uniform(-0.9,0.9,len(sgu)))
                kap=2*np.pi*fku*(1-np.cos(mu))/C_PC_S
                md=1.0/np.maximum(arru,1e-6)
                dmu=rng.rayleigh(md*np.sqrt(2/np.pi))
                dfq=1.0/(45*YR_S*np.maximum(arru,1e-6))
                dkap=np.sqrt((2*np.pi*(1-np.cos(mu))/C_PC_S)**2*dfq**2
                             +(2*np.pi*fku*np.sin(mu)/C_PC_S)**2*dmu**2)
                sgu=np.sqrt(sgu**2+(dkap*pri[i])**2)
                o2=sgu<np.pi
                if o2.sum()<2: continue
                sgu=sgu[o2]; kap=kap[o2]
                Lt=float(np.clip(rng.normal(0,pri[i]),-2.9*pri[i],2.9*pri[i]))
                Lg=np.linspace(-3*pri[i],3*pri[i],601)
                eps=kap*Lt+rng.normal(0,sgu)
                rr=(kap[:,None]*Lg[None,:]-eps[:,None]+np.pi)%(2*np.pi)-np.pi
                S=-(rr**2/(2*sgu[:,None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
                if 1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)),1e-30))>=TOOTH_W: continue
                Lh=float(Lg[int(np.argmax(S))])
                out.append(dict(sky=j,wrong=bool(abs(Lh-Lt)>TOOTH_W),sig=sgu.copy(),
                                kap=kap.copy(),pri=float(pri[i]),Ltrue=Lt))
        return out
    rr=book(_fixed())
    if rr["ok"]:
        print(f"  {A:10.1f} {rate:6.2f} {rr['gate']:7.1f} {rr['pur']:8.4f} {rr['y']:9.2f}")
        rows.append((A,rate,rr['gate'],rr['pur'],rr['y'],rr['yb'][0],rr['yb'][1]))
    else:
        print(f"  {A:10.1f} {rate:6.2f} {'--':>7} no 0.90 crossing")
        rows.append((A,rate,np.nan,np.nan,np.nan,np.nan,np.nan))
ys=[r[4] for r in rows if np.isfinite(r[4])]
if ys: print(f"    yield span {min(ys):.2f}-{max(ys):.2f} (factor {max(ys)/min(ys):.2f})")
with open("/home/mattm/projects/HSYMT/wheel_hasasia_examples/fig43_corner_vs_A.csv","w",newline="") as fh:
    cw=csv.writer(fh); cw.writerow(["harness","A_1e-15","count_scale","gate_pct","purity","psr_per_sky","q10","q90"])
    for r in rows: cw.writerow(["banked_adjudicated",r[0],r[1],f"{r[2]:.2f}",f"{r[3]:.4f}",f"{r[4]:.3f}",f"{r[5]:.3f}",f"{r[6]:.3f}"])
    if f["ok"]: cw.writerow(["banked_adjudicated","MARGINAL","",f"{f['gate']:.2f}",f"{f['pur']:.4f}",f"{f['y']:.3f}",f"{f['yb'][0]:.3f}",f"{f['yb'][1]:.3f}"])
print(f"  wrote fig43_corner_vs_A.csv")
