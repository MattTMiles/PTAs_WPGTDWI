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
print(f"\n{'='*74}\nJOB 29 -- T1: WHAT OWNS THE COMB ELEMENTS (gates everything)\n{'='*74}")
_ordAMP = _np.argsort(SRC_AMP); _AMPS = SRC_AMP[_ordAMP]
def _mc_z_of_amp(amp):
    k = _np.clip(_np.searchsorted(_AMPS, amp), 0, len(_AMPS)-1)
    i = _ordAMP[k]; return SRC_MC[i], SRC_Z[i]
_FE = 0.485                       # measured Earth-term fraction in band (Job 17)
def _unit(th, ph):
    return _np.array([_np.sin(th)*_np.cos(ph), _np.sin(th)*_np.sin(ph), _np.cos(th)])
def t_coal(f, ts): return (5.0/256.0)*ts**(-5.0/3.0)*(_np.pi*f)**(-8.0/3.0)
_TSUN = 4.925490947e-6; _MSUN = 1.989e30

def t1_attrib(j, T_yr=45.0, Np=20, s15=0.35):
    sd, f2, h2, P0, sea = _pop18(j % 50, s15)
    T = T_yr*YR_S
    ks, sig, fk, rp, nats = inventory20(f2, h2, P0, P0, sea, T_yr, 100.0, 14.0, 120)
    if len(ks) < 2: return None
    amp = h2/_np.maximum(f2, 1e-30)**(2.0/3.0)
    mc, zz = _mc_z_of_amp(amp)
    ts = _TSUN*mc*(1.0+zz)/_MSUN
    tc = t_coal(_np.maximum(f2, 1e-30), ts)
    rg = _np.random.default_rng(2900+j)
    pth = _np.arccos(rg.uniform(-1,1,Np)); pph = rg.uniform(0,2*_np.pi,Np)
    P = _np.stack([_unit(pth[a], pph[a]) for a in range(Np)])
    L = rg.uniform(0.5,3.0,Np)*1000.0
    sth = _np.arccos(rg.uniform(-1,1,len(f2))); sph = rg.uniform(0,2*_np.pi,len(f2))
    NS = _np.stack([_unit(sth[s], sph[s]) for s in range(len(f2))])
    CMU = NS @ P.T                                       # (Nsrc, Np)
    TAU = L[None,:]*(1-CMU)/C_PC_S
    FP = f2[:,None]*(1.0 + TAU/_np.maximum(tc[:,None],1e-9))**(-3.0/8.0)
    isres = _np.zeros(len(f2), bool); isres[ks] = True
    out = []
    for a, k in enumerate(ks):
        b = int(_np.clip(_np.digitize(f2[k], bin_edges)-1, 0, len(bin_cent)-1))
        wE = _np.sinc((f2 - f2[k])*T)**2                 # Earth-term weight
        wP = _np.sinc((FP - f2[k])*T)**2                 # pulsar-term weight
        PE = _FE*P0; PP = (1.0-_FE)*P0
        other = _np.ones(len(f2), bool); other[k] = False
        cE = float((wE[other]*PE[other]).sum())*(T/T_OBS)
        aP = ((wP*PP[:,None])[other & isres, :]).sum(0)*(T/T_OBS)
        bP = ((wP*PP[:,None])[other & ~isres, :]).sum(0)*(T/T_OBS)
        seab = sea[b]*(T_OBS/T)
        tot = aP + bP + cE + seab
        out.append((aP/_np.maximum(tot,1e-300), (bP+seab)/_np.maximum(tot,1e-300),
                    cE/_np.maximum(tot,1e-300)))
    return out

_A=[]; _B=[]; _C=[]
for j in range(40):
    r = t1_attrib(j)
    if r is None: continue
    for (a,b,c) in r:
        _A.append(a); _B.append(b); _C.append(_np.full_like(a, c))
_A=_np.concatenate(_A); _B=_np.concatenate(_B); _C=_np.concatenate(_C)
print(f"\n  Attribution of confusion power at comb elements"
      f" [{len(_A)} (element,pulsar) cases, 40 skies, T=45yr]:")
print(f"    (a) RESOLVED-source pulsar terms  : mean {_A.mean():.4f}"
      f"  median {_np.median(_A):.4f}  p90 {_np.quantile(_A,0.9):.4f}")
print(f"    (b) unresolved pulsar terms + sea : mean {_B.mean():.4f}"
      f"  median {_np.median(_B):.4f}")
print(f"    (c) Earth terms (all classes)     : mean {_C.mean():.4f}"
      f"  median {_np.median(_C):.4f}")
print(f"    frac of elements with (a) > 0.15  : {(_A>0.15).mean():.4f}")
print(f"    frac of elements with (a) > 0.30  : {(_A>0.30).mean():.4f}")
_gate = _A.mean()
print(f"\n  PRE-REGISTERED (a) = 0.20-0.45; measured {_gate:.4f} ->"
      f" {'HIT' if 0.20 <= _gate <= 0.45 else 'MISS'}")
print(f"  GATE (a) >= 0.15 to proceed to Task 2:"
      f" {'PASS -- build the joint decoder' if _gate >= 0.15 else 'FAIL -- STOP, sea owns the bins'}")
print(f"\n{'='*74}\nJOB 29 ADDENDUM -- attribution vs SENSITIVITY\n{'='*74}")
print("  CORRECTION: baseline used isres[ks] (top-nmax COMBS) as 'resolved'.")
print("  The detected set is det & rp>=rho_cut BEFORE nmax truncation. Fixed here;")
print("  both are reported so the baseline shift is visible.")
def _detmask(f2, h2, T_yr, sig_ns, cad, Np, depth_dex=0.0, rho_cut=3.0):
    nt = T_yr*365.25/cad
    rp = (h2/(2*_np.pi*f2))*_np.sqrt(nt/2)/(sig_ns*1e-9)
    arr = rp*_np.sqrt(0.2*Np)
    det = (arr >= 3.0*10**(-depth_dex)) & (f2 >= 2e-9) & (f2 <= 3e-8)
    return det & (rp >= rho_cut)

def attrib_sens(j, boost=1.0, T_yr=45.0, Np=20, s15=0.35, use_full=True, nmax=40):
    sig_ns = 100.0/boost
    sd, f2, h2, P0, sea = _pop18(j % 50, s15)
    T = T_yr*YR_S
    ks, sig, fk, rp, nats = inventory20(f2, h2, P0, P0, sea, T_yr, sig_ns, 14.0,
                                        120, nmax=nmax)
    if len(ks) < 2: return None
    full = _detmask(f2, h2, T_yr, sig_ns, 14.0, 120)
    isres = full if use_full else _np.isin(_np.arange(len(f2)), ks)
    amp = h2/_np.maximum(f2, 1e-30)**(2.0/3.0)
    mc, zz = _mc_z_of_amp(amp); ts = _TSUN*mc*(1.0+zz)/_MSUN
    tc = t_coal(_np.maximum(f2, 1e-30), ts)
    rg = _np.random.default_rng(2900+j)
    pth = _np.arccos(rg.uniform(-1,1,Np)); pph = rg.uniform(0,2*_np.pi,Np)
    Pv = _np.stack([_unit(pth[a], pph[a]) for a in range(Np)])
    L = rg.uniform(0.5,3.0,Np)*1000.0
    sth = _np.arccos(rg.uniform(-1,1,len(f2))); sph = rg.uniform(0,2*_np.pi,len(f2))
    NS = _np.stack([_unit(sth[s], sph[s]) for s in range(len(f2))])
    CMU = NS @ Pv.T
    TAU = L[None,:]*(1-CMU)/C_PC_S
    FP = f2[:,None]*(1.0 + TAU/_np.maximum(tc[:,None],1e-9))**(-3.0/8.0)
    out = []
    for a, k in enumerate(ks):
        b = int(_np.clip(_np.digitize(f2[k], bin_edges)-1, 0, len(bin_cent)-1))
        wE = _np.sinc((f2 - f2[k])*T)**2
        wP = _np.sinc((FP - f2[k])*T)**2
        PE = _FE*P0; PP = (1.0-_FE)*P0
        other = _np.ones(len(f2), bool); other[k] = False
        cE = float((wE[other]*PE[other]).sum())*(T/T_OBS)
        aP = ((wP*PP[:,None])[other & isres, :]).sum(0)*(T/T_OBS)
        bP = ((wP*PP[:,None])[other & ~isres, :]).sum(0)*(T/T_OBS)
        tot = aP + bP + cE + sea[b]*(T_OBS/T)
        out.append(aP/_np.maximum(tot, 1e-300))
    # co-location among RESOLVED sources
    bb = _np.digitize(f2[full], bin_edges)
    _, cnt = _np.unique(bb, return_counts=True)
    return out, int(full.sum()), float((cnt > 1).sum()/max(len(cnt),1))

print(f"\n  {'sensitivity':>12} {'resolved/sky':>13} {'co-loc rate':>12}"
      f" {'mean (a)':>9} {'median':>8} {'p90':>7} {'frac>0.15':>10} {'gate':>6}")
_cross = None
for _B in (1.0, 3.0, 10.0, 30.0, 100.0):
    A=[]; NR=[]; CL=[]
    for j in range(40):
        r = attrib_sens(j, boost=_B)
        if r is None: continue
        o, nr, cl = r; A.append(_np.concatenate([x[None,:] if x.ndim else x for x in o]) if False else _np.concatenate(o)); NR.append(nr); CL.append(cl)
    A=_np.concatenate([_np.atleast_1d(x) for x in A])
    m=A.mean(); g = "PASS" if m >= 0.15 else "fail"
    if m >= 0.15 and _cross is None: _cross = _B
    print(f"  {'x'+str(int(_B)):>12} {_np.mean(NR):13.1f} {_np.mean(CL):12.4f}"
          f" {m:9.4f} {_np.median(A):8.4f} {_np.quantile(A,0.9):7.4f}"
          f" {(A>0.15).mean():10.4f} {g:>6}")
print(f"\n  baseline recheck at x1 with COMBS-only 'resolved' (as first reported):")
A=[]
for j in range(40):
    r = attrib_sens(j, boost=1.0, use_full=False, nmax=12)
    if r is None: continue
    A.append(_np.concatenate(r[0]))
A=_np.concatenate([_np.atleast_1d(x) for x in A])
print(f"    mean (a) = {A.mean():.4f}  (first reported 0.0577)")
print(f"\n  VERDICT: {'ERA-DEPENDENT -- crosses 0.15 at sensitivity x'+str(int(_cross)) if _cross else 'CLOSED at all rungs tested (no crossing to x100)'}")
print(f"\n  THREE DEFINITIONS OF 'RESOLVED' -- the gate depends on which:")
print("    (1) COMBS      : top-nmax by nats, sigma<pi          [decoder's set]")
print("    (2) PHASE-SOLVABLE: all detected with sigma<pi       [what a JOINT")
print("                        decoder could actually model]")
print("    (3) DETECTED   : det & rp>=rho_cut, no sigma cut     [formal threshold]")
def _sig_all(f2, h2, P0, sea, T_yr, sig_ns, cad, Np, idx):
    T = T_yr*YR_S; nt = T_yr*365.25/cad
    rp = (h2/(2*_np.pi*f2))*_np.sqrt(nt/2)/(sig_ns*1e-9)
    out = _np.full(len(idx), _np.inf)
    for c0 in range(0, len(idx), 256):                      # chunked
        sub = idx[c0:c0+256]
        W = _np.sinc((f2[None,:] - f2[sub][:,None])*T)**2
        cc = W*P0[None,:]*(T/T_OBS)
        res = cc.sum(1) - cc[_np.arange(len(sub)), sub]
        bb = _np.clip(_np.digitize(f2[sub], bin_edges)-1, 0, len(bin_cent)-1)
        res = res + sea[bb]*(T_OBS/T)
        Pt = P0[sub]*(T/T_OBS)
        out[c0:c0+256] = _np.sqrt(_np.maximum(res,0)/_np.maximum(Pt,1e-300)
                                  + 1/_np.maximum(rp[sub],1e-9)**2)
    return out

def attrib3(j, boost=1.0, T_yr=45.0, Np=20, s15=0.35, nmax=40):
    sig_ns = 100.0/boost
    sd, f2, h2, P0, sea = _pop18(j % 50, s15)
    T = T_yr*YR_S
    ks, sg, fk, rp, nats = inventory20(f2, h2, P0, P0, sea, T_yr, sig_ns, 14.0,
                                       120, nmax=nmax)
    if len(ks) < 2: return None
    full = _np.where(_detmask(f2, h2, T_yr, sig_ns, 14.0, 120))[0]
    if not len(full): return None
    sig_f = _sig_all(f2, h2, P0, sea, T_yr, sig_ns, 14.0, 120, full)
    solv = full[sig_f < _np.pi]
    masks = {}
    masks[1] = _np.isin(_np.arange(len(f2)), ks)
    masks[2] = _np.isin(_np.arange(len(f2)), solv)
    masks[3] = _np.isin(_np.arange(len(f2)), full)
    amp = h2/_np.maximum(f2,1e-30)**(2.0/3.0)
    mc, zz = _mc_z_of_amp(amp); ts = _TSUN*mc*(1.0+zz)/_MSUN
    tc = t_coal(_np.maximum(f2,1e-30), ts)
    rg = _np.random.default_rng(2900+j)
    pth=_np.arccos(rg.uniform(-1,1,Np)); pph=rg.uniform(0,2*_np.pi,Np)
    Pv=_np.stack([_unit(pth[a],pph[a]) for a in range(Np)])
    L=rg.uniform(0.5,3.0,Np)*1000.0
    sth=_np.arccos(rg.uniform(-1,1,len(f2))); sph=rg.uniform(0,2*_np.pi,len(f2))
    NS=_np.stack([_unit(sth[s],sph[s]) for s in range(len(f2))])
    CMU=NS@Pv.T; TAU=L[None,:]*(1-CMU)/C_PC_S
    FP=f2[:,None]*(1.0+TAU/_np.maximum(tc[:,None],1e-9))**(-3.0/8.0)
    res = {}
    for d, m in masks.items():
        acc=[]
        for a, k in enumerate(ks):
            b=int(_np.clip(_np.digitize(f2[k],bin_edges)-1,0,len(bin_cent)-1))
            wE=_np.sinc((f2-f2[k])*T)**2; wP=_np.sinc((FP-f2[k])*T)**2
            PE=_FE*P0; PP=(1.0-_FE)*P0
            oth=_np.ones(len(f2),bool); oth[k]=False
            cE=float((wE[oth]*PE[oth]).sum())*(T/T_OBS)
            aP=((wP*PP[:,None])[oth&m,:]).sum(0)*(T/T_OBS)
            bP=((wP*PP[:,None])[oth&~m,:]).sum(0)*(T/T_OBS)
            tot=aP+bP+cE+sea[b]*(T_OBS/T)
            acc.append(aP/_np.maximum(tot,1e-300))
        res[d]=(_np.concatenate(acc), int(m.sum()))
    return res

print(f"\n  {'sens':>6} {'def':>4} {'N resolved/sky':>15} {'mean (a)':>9}"
      f" {'median':>8} {'frac>0.15':>10} {'gate':>6}")
_CR={1:None,2:None,3:None}
for _B in (1.0, 3.0, 10.0):
    ACC={1:[],2:[],3:[]}; NN={1:[],2:[],3:[]}
    for j in range(25):
        r=attrib3(j, boost=_B)
        if r is None: continue
        for d in (1,2,3): ACC[d].append(r[d][0]); NN[d].append(r[d][1])
    for d in (1,2,3):
        A=_np.concatenate(ACC[d]); m=A.mean()
        if m>=0.15 and _CR[d] is None: _CR[d]=_B
        print(f"  {'x'+str(int(_B)):>6} {d:>4} {_np.mean(NN[d]):15.1f} {m:9.4f}"
              f" {_np.median(A):8.4f} {(A>0.15).mean():10.4f}"
              f" {'PASS' if m>=0.15 else 'fail':>6}")
print(f"\n  crossing sensitivity by definition:"
      f" combs {_CR[1] or '>x10'}, phase-solvable {_CR[2] or '>x10'},"
      f" detected {_CR[3] or '>x10'}")
