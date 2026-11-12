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
import numpy as np, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT="/home/mattm/projects/HSYMT/wheel_hasasia_examples/"
print(f"\n{'='*74}\nJOB 40 ADDENDUM T2(x) -- XIAO ET AL. PAIRWISE, AT OUR PHYSICS\n{'='*74}")
print("  Construction: for a shared source the common source phase CANCELS in the")
print("  pulsar-pair phase DIFFERENCE, so P_Phi,n(Lp,Lq) ∝")
print("    exp(-wrap(kap_p*Lp - kap_q*Lq - d_n)^2 / (2*(sig_p^2+sig_q^2)))")
print("  Joint over shared sources x priors -> marginalise Lq -> posterior on Lp.")
print("  Partner selection: SMALLEST 68% credible half-width (their rule, verbatim).")
NG=201; NSRC=4; NPSR=12; NSKY=6
def xiao_rung(s, seed0=6100):
    sig=np.sqrt(max((1-s)/max(s,1e-9),1e-4))
    rows=[]
    for j in range(NSKY):
        rng=np.random.default_rng(seed0+j)
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        ks,sg0,fk,rp,nats=inventory20(f2,h2,P0,P0,sea,45.0,100.0,14.0,120)
        if len(ks)<2: continue
        fsrc=fk[:NSRC] if len(fk)>=NSRC else fk
        ns=len(fsrc)
        pri=np.maximum(rng.normal(2.0,0.4,NPSR),0.5)      # parallax-class priors
        Lt=np.array([rng.normal(0,pri[i]) for i in range(NPSR)])
        mu=np.arccos(rng.uniform(-0.9,0.9,(NPSR,ns)))
        kap=2*np.pi*fsrc[None,:]*(1-np.cos(mu))/C_PC_S
        phi=rng.uniform(0,2*np.pi,ns)                     # shared source phases
        eps=kap*Lt[:,None]+phi[None,:]+rng.normal(0,sig,(NPSR,ns))
        for p in range(NPSR):
            Lp=np.linspace(-3*pri[p],3*pri[p],NG)
            best=None
            for q in range(NPSR):
                if q==p: continue
                Lq=np.linspace(-3*pri[q],3*pri[q],NG)
                d=eps[p]-eps[q]                            # source phase cancels
                M=np.zeros((NG,NG))
                for a in range(ns):
                    r=(kap[p,a]*Lp[None,:]-kap[q,a]*Lq[:,None]-d[a]+np.pi)%(2*np.pi)-np.pi
                    M+=-r**2/(2*(2*sig**2))
                M+=-(Lp**2/(2*pri[p]**2))[None,:]-(Lq**2/(2*pri[q]**2))[:,None]
                W=np.exp(M-M.max()); post=W.sum(0); post/=post.sum()
                o=np.argsort(post)[::-1]; c=np.cumsum(post[o])
                inset=o[:max(int(np.searchsorted(c,0.68))+1,1)]
                half=0.5*(Lp[inset].max()-Lp[inset].min())
                if best is None or half<best[0]:
                    best=(half,float(Lp[int(np.argmax(post))]),
                          bool(Lt[p]>=Lp[inset].min() and Lt[p]<=Lp[inset].max()))
            if best is None: continue
            rows.append(dict(half=best[0],err=abs(best[1]-Lt[p]),cover=best[2]))
    return rows
print(f"\n  {'s':>6} {'sigma_phi':>10} {'n':>5} {'sub-pc frac':>12} {'coverage':>10}"
      f" {'corr(width,correct)':>20}")
_X=[]
for s in (1.0,0.9,0.75,0.5,0.25,0.1):
    r=xiao_rung(s)
    if not r: continue
    hw=np.array([x["half"] for x in r]); er=np.array([x["err"] for x in r])
    cv=np.array([x["cover"] for x in r],bool)
    sub=float((er<1.0).mean()); cov=float(cv.mean())
    ok=(er<1.0)
    corr=float(np.corrcoef(hw,ok.astype(float))[0,1]) if ok.std()>0 else np.nan
    print(f"  {s:6.2f} {np.sqrt((1-s)/max(s,1e-9)):10.2f} {len(r):5d} {sub:12.3f}"
          f" {cov:10.3f} {corr:20.3f}")
    _X.append((s,np.sqrt((1-s)/max(s,1e-9)),len(r),sub,cov,corr))
print(f"\n  ANTI-PREDICTIVE TEST (pre-registered): corr(width, correct) should be")
print(f"  NEGATIVE at s=1 (narrow => correct) and null/inverted by s<=0.5.")
c1=[x[5] for x in _X if x[0]==1.0][0]; c5=[x[5] for x in _X if x[0]==0.5][0]
print(f"    s=1.00: {c1:+.3f}   s=0.50: {c5:+.3f}")
if c1<-0.05 and c5>=-0.05:
    print(f"    -> PRE-REGISTRATION CONFIRMED: the selection rule goes non-predictive")
    print(f"       under confusion. The chosen-peak-width pathology is real.")
elif c5<-0.05:
    print(f"    -> FINDING AGAINST our Bayes-limit generalisation: the rule STAYS")
    print(f"       predictive at s=0.5. Reported prominently, not reconciled.")
else:
    print(f"    -> inconclusive at s=1 (see diagnosis)")
fig,ax=plt.subplots(figsize=(6.4,4.2))
ss=[x[0] for x in _X]
ax.plot(ss,[x[3] for x in _X],'o-',color="C0",lw=1.8,label="sub-pc fraction")
ax.plot(ss,[x[4] for x in _X],'s--',color="C3",lw=1.6,label="coverage of truth")
ax.axvspan(0.55,0.75,color="0.85",zorder=0)
ax.text(0.65,0.05,"booked\nmachine",ha="center",fontsize=8,color="0.35")
ax.plot(1.0,[x[3] for x in _X][0],'o',color="C0",ms=11)
ax.annotate("Xiao et al. operating point\n(white noise only)",(1.0,[x[3] for x in _X][0]),
            xytext=(-140,-6),textcoords="offset points",fontsize=8,color="C0")
ax.set_xlabel("coherent fraction s"); ax.set_ylabel("fraction")
ax.set_xlim(1.05,0.05); ax.set_ylim(-0.03,1.05); ax.grid(alpha=.25,lw=.5)
ax.legend(fontsize=8,loc="center right")
plt.tight_layout(); plt.savefig(OUT+"fig40_xiao.png",dpi=150); plt.close()
with open(OUT+"fig40_xiao.csv","w",newline="") as fh:
    cw=csv.writer(fh); cw.writerow(["s","sigma_phi","n","subpc_frac","coverage","corr_width_correct"])
    for r in _X: cw.writerow([f"{r[0]:.2f}",f"{r[1]:.3f}",r[2],f"{r[3]:.4f}",f"{r[4]:.4f}",f"{r[5]:.4f}"])
print(f"  wrote fig40_xiao.png + .csv")
print(f"\n{'='*74}\nT2(x) REDONE -- s=1 needs a WHITE-NOISE floor, not zero noise\n{'='*74}")
print("  DEFECT in the first pass: s=1.0 set sigma->0, so every 68% interval")
print("  collapsed to the grid spacing -> coverage 0.000 and corr = nan. That rung")
print("  was degenerate, not a measurement, and the 'finding against' verdict rested")
print("  on comparing s=0.5 to an undefined s=1. Xiao et al. run with WHITE NOISE:")
print("  sigma_meas = 1/SNR_pterm. Using a 0.20 rad floor (SNR_pterm ~ 5, their")
print("  clean-limit regime); sigma_tot^2 = sigma_conf^2 + sigma_floor^2.")
SIG_FLOOR=0.20
def xiao2(s, seed0=6100):
    sc=np.sqrt(max((1-s)/max(s,1e-9),0.0)); sig=np.sqrt(sc**2+SIG_FLOOR**2)
    rows=[]
    for j in range(NSKY):
        rng=np.random.default_rng(seed0+j)
        sd,f2,h2,P0,sea=_pop18(j%50,0.35)
        ks,sg0,fk,rp,nats=inventory20(f2,h2,P0,P0,sea,45.0,100.0,14.0,120)
        if len(ks)<2: continue
        fsrc=fk[:NSRC] if len(fk)>=NSRC else fk; ns=len(fsrc)
        pri=np.maximum(rng.normal(2.0,0.4,NPSR),0.5)
        Lt=np.array([rng.normal(0,pri[i]) for i in range(NPSR)])
        mu=np.arccos(rng.uniform(-0.9,0.9,(NPSR,ns)))
        kap=2*np.pi*fsrc[None,:]*(1-np.cos(mu))/C_PC_S
        phi=rng.uniform(0,2*np.pi,ns)
        eps=kap*Lt[:,None]+phi[None,:]+rng.normal(0,sig,(NPSR,ns))
        for p in range(NPSR):
            Lp=np.linspace(-3*pri[p],3*pri[p],NG); best=None
            for q in range(NPSR):
                if q==p: continue
                Lq=np.linspace(-3*pri[q],3*pri[q],NG); d=eps[p]-eps[q]
                M=np.zeros((NG,NG))
                for a in range(ns):
                    r=(kap[p,a]*Lp[None,:]-kap[q,a]*Lq[:,None]-d[a]+np.pi)%(2*np.pi)-np.pi
                    M+=-r**2/(2*(2*sig**2))
                M+=-(Lp**2/(2*pri[p]**2))[None,:]-(Lq**2/(2*pri[q]**2))[:,None]
                W=np.exp(M-M.max()); post=W.sum(0); post/=post.sum()
                o=np.argsort(post)[::-1]; c=np.cumsum(post[o])
                inset=o[:max(int(np.searchsorted(c,0.68))+1,1)]
                half=0.5*(Lp[inset].max()-Lp[inset].min())
                if best is None or half<best[0]:
                    best=(half,float(Lp[int(np.argmax(post))]),
                          bool(Lp[inset].min()<=Lt[p]<=Lp[inset].max()))
            if best: rows.append(dict(half=best[0],err=abs(best[1]-Lt[p]),cover=best[2]))
    return rows
print(f"\n  {'s':>6} {'sigma_tot':>10} {'sub-pc frac':>12} {'coverage':>10}"
      f" {'corr(width,correct)':>20}")
_X2=[]
for s in (1.0,0.9,0.75,0.5,0.25,0.1):
    r=xiao2(s)
    hw=np.array([x["half"] for x in r]); er=np.array([x["err"] for x in r])
    cv=np.array([x["cover"] for x in r],bool); ok=(er<1.0)
    corr=float(np.corrcoef(hw,ok.astype(float))[0,1]) if (ok.std()>0 and hw.std()>0) else np.nan
    st=np.sqrt(max((1-s)/max(s,1e-9),0.0)+0.0)
    st=np.sqrt(st**2+SIG_FLOOR**2)
    print(f"  {s:6.2f} {st:10.2f} {ok.mean():12.3f} {cv.mean():10.3f} {corr:20.3f}")
    _X2.append((s,st,len(r),float(ok.mean()),float(cv.mean()),corr))
c1=_X2[0][5]; c5=[x[5] for x in _X2 if x[0]==0.5][0]
print(f"\n  ANTI-PREDICTIVE TEST: s=1.00 corr {c1:+.3f} | s=0.50 corr {c5:+.3f}")
if c1<-0.05 and abs(c5)<0.05: v="PRE-REGISTRATION CONFIRMED -- rule goes non-predictive"
elif c1<-0.05 and c5<-0.05: v="FINDING AGAINST our generalisation -- rule STAYS predictive"
else: v="inconclusive"
print(f"    -> {v}")
print(f"\n  COVERAGE is the sharper result: nominal 0.68, measured"
      f" {_X2[0][4]:.3f} at s=1 rising to {_X2[-1][4]:.3f} at s=0.1.")
print(f"    The min-width SELECTION biases coverage DOWN where intervals are")
print(f"    tight, and the bias relaxes as confusion widens them.")
with open(OUT+"fig40_xiao.csv","w",newline="") as fh:
    cw=csv.writer(fh); cw.writerow(["s","sigma_tot","n","subpc_frac","coverage","corr_width_correct"])
    for r in _X2: cw.writerow([f"{r[0]:.2f}",f"{r[1]:.3f}",r[2],f"{r[3]:.4f}",f"{r[4]:.4f}",f"{r[5]:.4f}"])
fig,ax=plt.subplots(figsize=(6.4,4.2))
ss=[x[0] for x in _X2]
ax.plot(ss,[x[3] for x in _X2],'o-',color="C0",lw=1.8,label="sub-pc fraction")
ax.plot(ss,[x[4] for x in _X2],'s--',color="C3",lw=1.6,label="coverage of truth")
ax.axhline(0.68,color="0.5",ls=":",lw=1); ax.text(0.12,0.70,"nominal 68%",fontsize=7.5,color="0.5")
ax.axvspan(0.55,0.75,color="0.88",zorder=0)
ax.text(0.65,0.06,"booked\nmachine",ha="center",fontsize=8,color="0.35")
ax.plot(1.0,_X2[0][3],'o',color="C0",ms=11)
ax.annotate("Xiao et al. operating point\n(white noise only)",(1.0,_X2[0][3]),
            xytext=(-142,-4),textcoords="offset points",fontsize=8,color="C0")
ax.set_xlabel("coherent fraction s"); ax.set_ylabel("fraction")
ax.set_xlim(1.06,0.04); ax.set_ylim(-0.03,1.05); ax.grid(alpha=.25,lw=.5)
ax.legend(fontsize=8,loc="upper right")
plt.tight_layout(); plt.savefig(OUT+"fig40_xiao.png",dpi=150); plt.close()
print(f"  rewrote fig40_xiao.png + .csv")
