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
print(f"\n{'='*74}\nJOB 30 -- THE MIDDLE TIER (Task 4 redo)\n{'='*74}")
print("  GATE-CHAIN CORRECTION: binding gate is the top-40-by-rp TRUNCATION,")
print("  not array SNR>3 (arr>=3 <=> rp>=0.612 at Np=120; non-binding).")
print("  Your rungs 2.5/2.0/1.5 are mapped onto the truncation cap instead.")

def midtier(j, T_yr=45.0, Np=40):
    """combs vs phase-solvable-but-excluded, with coupled-angle R(theta)."""
    sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
    T = T_yr*YR_S; nt = T_yr*365.25/14.0
    rp = (h2/(2*_np.pi*f2))*_np.sqrt(nt/2)/(100.0*1e-9)
    arr = rp*_np.sqrt(0.2*120)
    cand = _np.where((arr>=3.0)&(f2>=2e-9)&(f2<=3e-8)&(rp>=3.0))[0]
    if len(cand) < 2: return None
    sg = sig_windowed(f2, P0, sea, rp, cand, T)
    solv = cand[sg < _np.pi]                       # phase-solvable census
    if len(solv) < 2: return None
    keep = cand[_np.argsort(rp[cand])[::-1][:40]]  # the banked truncation
    sgk = sig_windowed(f2, P0, sea, rp, keep, T)
    kk = keep[sgk < _np.pi]
    nats = _np.log(_np.pi/sgk[sgk < _np.pi])
    combs = kk[_np.argsort(nats)[::-1][:12]]
    mid = _np.setdiff1d(solv, combs)
    rg = _np.random.default_rng(3000+j)
    F = []
    for p in range(Np):                            # coupled angle: R from mu
        muC = _np.arccos(rg.uniform(-0.9,0.9,len(combs)))
        muM = _np.arccos(rg.uniform(-0.9,0.9,len(mid)))
        iC = float((rp[combs]**2*(1+_np.cos(muC))**2/4).sum())
        iM = float((rp[mid]**2*(1+_np.cos(muM))**2/4).sum())
        F.append(iM/max(iC,1e-30))
    return F, len(combs), len(mid)
_F=[]; _nc=[]; _nm=[]
for j in range(100):
    r = midtier(j)
    if r is None: continue
    _F += r[0]; _nc.append(r[1]); _nm.append(r[2])
_F=_np.array(_F)
print(f"\n  T1 -- middle-tier information fraction F = sum_mid rho^2 R / sum_combs rho^2 R")
print(f"    combs/sky {_np.mean(_nc):.2f}   middle-tier/sky {_np.mean(_nm):.2f}"
      f"   [{len(_F)} pulsar cases, 100 skies, T=45]")
print(f"    F: mean {_F.mean():.4f}  median {_np.median(_F):.4f}"
      f"  p10 {_np.quantile(_F,0.1):.4f}  p90 {_np.quantile(_F,0.9):.4f}")
print(f"    PRE-REG 0.15-0.50: {'HIT' if 0.15 <= _F.mean() <= 0.50 else 'MISS'}"
      f"   GATE F>=0.10: {'PASS -> build T2' if _F.mean() >= 0.10 else 'FAIL -> STOP'}")

# ---------------- T2: alias-gate percentiles + scrambled-position null ----------
def alias_of(w, Lg, ip):
    m = _np.abs(Lg - Lg[ip]) > TOOTH_W
    return float(w[m].max()/max(w[ip], 1e-300)) if m.any() else 0.0
def dec30(cap, scramble_mid=False, nsky=30, nps=40, T_yr=45.0, M=99):
    rows=[]; ranks=[]; nrn=0
    for j in range(nsky):
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        r = inv_cap(f2, h2, P0, sea, T_yr, cap)
        if r is None: continue
        ks, sg0, fk, rpk, arrk = r; T = T_yr*YR_S
        base = inv_cap(f2, h2, P0, sea, T_yr, 40)
        bset = set(base[0].tolist()) if base is not None else set()
        ismid = _np.array([k not in bset for k in ks])
        rg = _np.random.default_rng(2950+j)
        Ti=_np.clip(_np.exp(rg.normal(_np.log(T_yr),0.35,120)),5,90)
        beta=_np.arcsin(rg.uniform(-1,1,120)); Dk=rg.uniform(0.5,3.0,120)
        px=2.34/_np.cos(beta)**2*Dk**2*10.0*_np.sqrt(30.0/Ti); pb=3.0*(Ti/30.0)**(-2.5)
        pri=_np.maximum(_np.where(rg.random(120)<0.7,_np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(120); rncl=Arn==0.0
        idx=_np.where(rncl)[0][:nps]; nrn+=len(idx)
        for i in idx:
            sg=_np.array([depro_sigma(
                _np.sqrt(max(sg0[a]**2-1/max(rpk[a],1e-9)**2,0.0)),
                rpk[a], fk[a], i, "loo") for a in range(len(ks))])
            sg=sg*_np.sqrt(_np.mean(Ti)/Ti[i])
            mu=_np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            sg=sg*_np.sqrt((1/3)/_np.maximum((1+_np.cos(mu))**2/4,1e-3))
            mu_used = mu.copy()
            if scramble_mid and ismid.any():        # NULL: scramble middle-tier sky
                mu_used[ismid] = rg.permutation(mu[ismid])
                mu_used[ismid] = _np.arccos(rg.uniform(-0.9,0.9,int(ismid.sum())))
            kap=2*_np.pi*fk*(1-_np.cos(mu_used))/C_PC_S
            dmu=1.0/_np.maximum(arrk,1e-6); dfq=1.0/(T*_np.maximum(arrk,1e-6))
            dkap=_np.sqrt((2*_np.pi*(1-_np.cos(mu))/C_PC_S)**2*dfq**2
                          +(2*_np.pi*fk*_np.sin(mu)/C_PC_S)**2*dmu**2)
            sg=_np.sqrt(sg**2+(dkap*pri[i])**2)
            ok=sg<_np.pi
            if ok.sum()<2: continue
            sgu,kapu=sg[ok],kap[ok]; A=1.0/_np.maximum(sgu,1e-6)
            w_=1.0/sgu**2; neff=float(w_.sum()**2/_np.sum(w_**2))
            Lt=rg.normal(0,pri[i])
            z=A*_np.exp(1j*(kapu*Lt))+(rg.normal(0,1/_np.sqrt(2),len(A))
                                       +1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            Lg=_np.linspace(-3*pri[i],3*pri[i],_GRID)
            Sc=(2*A[:,None]*_np.real(z[:,None]*_np.exp(
                -1j*(kapu[:,None]*Lg[None,:])))).sum(0)-Lg**2/(2*pri[i]**2)
            wpost=_np.exp(Sc-Sc.max()); wpost/=wpost.sum()
            ip=int(_np.argmax(Sc)); hit=abs(float(Lg[ip])-Lt)<=TOOTH_W
            cert = 1.0/_np.sqrt(max(float(_np.sum(kapu**2/sgu**2)),1e-30)) < TOOTH_W
            rows.append((hit, cert, alias_of(wpost,Lg,ip), neff))
            ranks.append(rank_mass_mu(A,kapu,pri[i],Lt,z,Lg,0.0,rg,M=M))
    H=_np.array([r[0] for r in rows]); C=_np.array([r[1] for r in rows])
    AL=_np.array([r[2] for r in rows]); NE=_np.array([r[3] for r in rows])
    out=dict(neff=NE.mean(), pur=H.mean(), rank=float(_np.mean(ranks)),
             y=C.sum()/max(nrn,1), gpur=H[C].mean() if C.any() else _np.nan)
    for q in (0.20,0.10,0.07):
        if C.any():
            thr=_np.quantile(AL[C], q); m=C&(AL<=thr)
            out[f"p{int(q*100)}"]=(H[m].mean() if m.any() else _np.nan, m.sum()/max(nrn,1))
    return out
print(f"\n  T2 -- CONFIGURATION: T=45, heterogeneous spans, truth drawn from prior,")
print(f"        fuzzy templates ON. Baseline ungated purity must be 0.45-0.60.")
print(f"  {'arm':>22} {'N_eff':>6} {'rank':>7} {'ungated':>8} {'gated':>7} {'yield':>7}"
      f" {'20%':>13} {'10%':>13} {'7%':>13}")
_res={}
for _lab,_cap,_sc in (("combs-only (cap 40)",40,False), ("rung 2.5 (cap 100)",100,False),
                      ("rung 2.0 (cap 400)",400,False), ("rung 1.5 (cap 1300)",1300,False),
                      ("NULL scrambled sky",1300,True)):
    o=dec30(_cap, scramble_mid=_sc); _res[_lab]=o
    f=lambda k: f"{o[k][0]:.3f}/{o[k][1]:.3f}" if k in o else "  --  "
    print(f"  {_lab:>22} {o['neff']:6.2f} {o['rank']:7.4f} {o['pur']:8.4f}"
          f" {o['gpur']:7.4f} {o['y']:7.4f} {f('p20'):>13} {f('p10'):>13} {f('p7'):>13}")
_b=_res["combs-only (cap 40)"]; _r2=_res["rung 2.0 (cap 400)"]; _n=_res["NULL scrambled sky"]
print(f"    baseline ungated purity {_b['pur']:.4f}"
      f" -> {'IN 0.45-0.60 (headroom OK)' if 0.45<=_b['pur']<=0.60 else 'OUT OF RANGE'}")
print(f"    N_eff x{_r2['neff']/_b['neff']:.3f} at rung 2.0"
      f"  [PRE-REG 1.3-2.0: {'HIT' if 1.3<=_r2['neff']/_b['neff']<=2.0 else 'MISS'}]")
print(f"    ungated purity x{_r2['pur']/_b['pur']:.3f}"
      f"  [PRE-REG 1.05-1.25: {'HIT' if 1.05<=_r2['pur']/_b['pur']<=1.25 else 'MISS'}]")
print(f"    NULL scrambled-sky purity x{_n['pur']/_b['pur']:.3f}"
      f"  [PRE-REG <1.03: {'PASS' if _n['pur']/_b['pur']<1.03 else 'FAIL -- weight bookkeeping'}]")
print(f"\n  NULL FIX: previous scrambled-sky null permuted mu BEFORE generating z,")
print("  so data and template shared the scrambled geometry -- it tested nothing.")
print("  Correct null: GENERATE with true kappa, DECODE with middle-tier kappa")
print("  scrambled. Middle-tier templates are then wrong and must contribute noise.")
def dec30b(cap, scramble_mid=False, nsky=30, nps=40, T_yr=45.0, M=99):
    H=[]; C=[]; ranks=[]; nrn=0; NE=[]
    for j in range(nsky):
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        r = inv_cap(f2, h2, P0, sea, T_yr, cap)
        if r is None: continue
        ks, sg0, fk, rpk, arrk = r; T = T_yr*YR_S
        base = inv_cap(f2, h2, P0, sea, T_yr, 40)
        bset = set(base[0].tolist()) if base is not None else set()
        ismid = _np.array([k not in bset for k in ks])
        rg = _np.random.default_rng(2950+j)
        Ti=_np.clip(_np.exp(rg.normal(_np.log(T_yr),0.35,120)),5,90)
        beta=_np.arcsin(rg.uniform(-1,1,120)); Dk=rg.uniform(0.5,3.0,120)
        px=2.34/_np.cos(beta)**2*Dk**2*10.0*_np.sqrt(30.0/Ti); pb=3.0*(Ti/30.0)**(-2.5)
        pri=_np.maximum(_np.where(rg.random(120)<0.7,_np.minimum(px,pb),px),0.05)
        Arn,_g=rn_draw(120); rncl=Arn==0.0
        idx=_np.where(rncl)[0][:nps]; nrn+=len(idx)
        for i in idx:
            sg=_np.array([depro_sigma(
                _np.sqrt(max(sg0[a]**2-1/max(rpk[a],1e-9)**2,0.0)),
                rpk[a], fk[a], i, "loo") for a in range(len(ks))])
            sg=sg*_np.sqrt(_np.mean(Ti)/Ti[i])
            mu=_np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            sg=sg*_np.sqrt((1/3)/_np.maximum((1+_np.cos(mu))**2/4,1e-3))
            kap_true=2*_np.pi*fk*(1-_np.cos(mu))/C_PC_S      # GENERATES the data
            kap_tmpl=kap_true.copy()
            if scramble_mid and ismid.any():                  # DECODER sees wrong geometry
                mus=_np.arccos(rg.uniform(-0.9,0.9,int(ismid.sum())))
                kap_tmpl[ismid]=2*_np.pi*fk[ismid]*(1-_np.cos(mus))/C_PC_S
            dmu=1.0/_np.maximum(arrk,1e-6); dfq=1.0/(T*_np.maximum(arrk,1e-6))
            dkap=_np.sqrt((2*_np.pi*(1-_np.cos(mu))/C_PC_S)**2*dfq**2
                          +(2*_np.pi*fk*_np.sin(mu)/C_PC_S)**2*dmu**2)
            sg=_np.sqrt(sg**2+(dkap*pri[i])**2)
            ok=sg<_np.pi
            if ok.sum()<2: continue
            sgu=sg[ok]; kg=kap_true[ok]; kt=kap_tmpl[ok]
            A=1.0/_np.maximum(sgu,1e-6)
            w_=1.0/sgu**2; NE.append(float(w_.sum()**2/_np.sum(w_**2)))
            Lt=rg.normal(0,pri[i])
            z=A*_np.exp(1j*(kg*Lt))+(rg.normal(0,1/_np.sqrt(2),len(A))
                                     +1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            Lg=_np.linspace(-3*pri[i],3*pri[i],_GRID)
            Sc=(2*A[:,None]*_np.real(z[:,None]*_np.exp(
                -1j*(kt[:,None]*Lg[None,:])))).sum(0)-Lg**2/(2*pri[i]**2)
            ip=int(_np.argmax(Sc)); H.append(abs(float(Lg[ip])-Lt)<=TOOTH_W)
            C.append(1.0/_np.sqrt(max(float(_np.sum(kt**2/sgu**2)),1e-30))<TOOTH_W)
            ranks.append(rank_mass_mu(A,kt,pri[i],Lt,z,Lg,0.0,rg,M=M))
    H=_np.array(H); C=_np.array(C)
    return dict(pur=H.mean(), gpur=H[C].mean() if C.any() else _np.nan,
                y=C.sum()/max(nrn,1), rank=float(_np.mean(ranks)), neff=_np.mean(NE))
print(f"\n  {'arm':>26} {'N_eff':>6} {'rank':>7} {'ungated':>8} {'gated':>7} {'yield':>7} {'right/RNc':>10}")
_o={}
for _lab,_cap,_sc in (("combs-only (cap 40)",40,False),
                      ("cap 1300, true geometry",1300,False),
                      ("cap 1300, SCRAMBLED mid",1300,True)):
    o=dec30b(_cap,_sc); _o[_lab]=o
    print(f"  {_lab:>26} {o['neff']:6.2f} {o['rank']:7.4f} {o['pur']:8.4f}"
          f" {o['gpur']:7.4f} {o['y']:7.4f} {o['gpur']*o['y']:10.4f}")
_b=_o["combs-only (cap 40)"]; _t=_o["cap 1300, true geometry"]; _s=_o["cap 1300, SCRAMBLED mid"]
_gain=_t['pur']-_b['pur']; _sg=_s['pur']-_b['pur']
print(f"\n    true gain in ungated purity  : {_gain:+.4f}")
print(f"    gain surviving scrambling    : {_sg:+.4f}"
      f"  ({_sg/max(_gain,1e-9):.1%} of it)")
print(f"    scrambled/baseline ratio {_s['pur']/_b['pur']:.3f}"
      f"  [PRE-REG <1.03: {'PASS -- gain is INFORMATION' if _s['pur']/_b['pur']<1.03 else 'FAIL -- weight bookkeeping'}]")
print(f"    rank: baseline {_b['rank']:.4f} | true {_t['rank']:.4f} | scrambled {_s['rank']:.4f}")
