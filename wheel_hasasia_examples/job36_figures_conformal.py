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
import numpy as np
print(f"\n{'='*74}\nJOB 33 -- T1 (R15): CHIRP ACROSS BINS\n{'='*74}")
print("  (a) AUDIT [by inspection of the banked SNR definition]:")
print("      rp = (h/(2*pi*f)) * sqrt(nt/2) / sigma_ns")
print("      sqrt(nt/2) is the FULL coherent gain over the whole span. That is the")
print("      matched-filter SNR of a signal whose phase is tracked for the entire")
print("      observation -- i.e. it ALREADY ASSUMES perfect chirp tracking on the")
print("      SIGNAL side. There is no monochromatic smearing loss to recover:")
print("      the booked numbers already sit at the chirplet-capture optimum.")
print("      => T1(b)/(c) SIGNAL-side repair is a NULL-OP; the pre-registered")
print("         'purity rises 0.731 -> 0.75-0.82' cannot be earned this way.")
print("      The real exposure runs the OTHER way: a genuinely monochromatic")
print("      pipeline would be WORSE than booked, not better.")
print("  BUT the CONFUSION side is NOT modelled: inventory20's res = sum_j w_j P_j")
print("  counts every in-bin source as fully confusing regardless of its chirp")
print("  rate. A chirp-matched template rejects static confusers. That gain is")
print("  genuinely unbanked, and is what T1(c)/(d) can actually measure.")

_TSUN=4.925490947e-6; _MSUN=1.989e30
_ordAMP2=np.argsort(SRC_AMP); _AMPS2=SRC_AMP[_ordAMP2]
def mcz_of(amp):
    k=np.clip(np.searchsorted(_AMPS2,amp),0,len(_AMPS2)-1); i=_ordAMP2[k]
    return SRC_MC[i],SRC_Z[i]
def drift_bins(f, ts, T):
    fdot=(96.0/5.0)*np.pi**(8.0/3.0)*ts**(5.0/3.0)*f**(11.0/3.0)
    return fdot*T*T
print(f"\n  (b) DRIFT TABLE -- verification of the brief's premise:")
print(f"    {'Mc':>10} {'f (nHz)':>8} {'T=45':>7} {'T=60':>7}")
for _mc,_f in ((3e9,30e-9),(5e9,20e-9),(5e9,30e-9)):
    ts=_TSUN*_mc
    print(f"    {_mc:10.1e} {_f*1e9:8.0f}"
          f" {drift_bins(_f,ts,45*YR_S):7.2f} {drift_bins(_f,ts,60*YR_S):7.2f}")
print(f"    (matches the brief's table -> premise verified)")

print(f"\n  (b2) DRIFT IN THE ACTUAL COMB POPULATION, by Mc decile [T=45]:")
_dr=[]; _mcs=[]
for j in range(30):
    sd,f2,h2,P0,sea=_pop18(j%50,0.35)
    ks,sig,fk,rp,nats=inventory20(f2,h2,P0,P0,sea,45.0,100.0,14.0,120)
    if len(ks)<2: continue
    amp=h2[ks]/np.maximum(f2[ks],1e-30)**(2.0/3.0)
    mc,zz=mcz_of(amp); ts=_TSUN*mc*(1.0+zz)/_MSUN
    _dr.append(drift_bins(fk,ts,45*YR_S)); _mcs.append(mc/1.989e30)
_dr=np.concatenate(_dr); _mcs=np.concatenate(_mcs)
_q=np.quantile(_mcs,[0.5,0.9])
for _lab,_m in (("Mc below median",_mcs<_q[0]),("Mc 50-90 pct",(_mcs>=_q[0])&(_mcs<_q[1])),
                ("Mc top decile (MONSTERS)",_mcs>=_q[1])):
    print(f"    {_lab:>26}: median drift/bin {np.median(_dr[_m]):8.3f}"
          f"  frac>0.5 {np.mean(_dr[_m]>0.5):.3f}  frac>1 {np.mean(_dr[_m]>1):.3f}")
print(f"    ALL combs: median drift {np.median(_dr):.3f}, frac>0.5 {np.mean(_dr>0.5):.3f}")
print(f"\n  (d) DISCRIMINANT: chirp-rate mismatch rejection ~ 1/drift in power.")
_sup=1.0/np.maximum(_dr,1.0)
print(f"    mean static-confusion rejection over combs with drift>0.5:"
      f" {np.mean(_sup[_dr>0.5]):.4f} (i.e. x{1/max(np.mean(_sup[_dr>0.5]),1e-9):.2f} suppression)")
print(f"    fraction of combs that can use it (drift>0.5): {np.mean(_dr>0.5):.4f}")

# ---- (c) chirp-rate confusion rejection, folded into the booked machine -------
_CHIRP33=False
_inv_j32 = inventory20
def inventory20(f2,h2,Psig,Pbkg,sea,T_yr,sig_ns,cad,Np,depth_dex=0.0,rho_cut=3.0,nmax=12):
    T=T_yr*YR_S; nt=T_yr*365.25/cad
    rp=(h2/(2*np.pi*f2))*np.sqrt(nt/2)/(sig_ns*1e-9)
    arr=rp*np.sqrt(0.2*Np)
    det=(arr>=3.0*10**(-depth_dex))&(f2>=2e-9)&(f2<=3e-8)
    ks=np.where(det&(rp>=rho_cut))[0]
    if len(ks)==0: return ks,np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0)
    if len(ks)>_CAP32: ks=ks[np.argsort(rp[ks])[::-1][:_CAP32]]
    o_=np.argsort(f2); fs=f2[o_]; Pb=Pbkg[o_]; half=100.0/T
    if _CHIRP33:
        amp_=h2[ks]/np.maximum(f2[ks],1e-30)**(2.0/3.0)
        mc_,z_=mcz_of(amp_); ts_=_TSUN*mc_*(1.0+z_)/_MSUN
        drf=drift_bins(f2[ks],ts_,T)
    sig=np.empty(len(ks))
    for a,k in enumerate(ks):
        lo=np.searchsorted(fs,f2[k]-half); hi=np.searchsorted(fs,f2[k]+half)
        w=np.sinc((fs[lo:hi]-f2[k])*T)**2
        c=w*Pb[lo:hi]*(T/T_OBS)
        b=int(np.clip(np.digitize(f2[k],bin_edges)-1,0,len(bin_cent)-1))
        res=c.sum()-Pbkg[k]*(T/T_OBS)+sea[b]*(T_OBS/T)
        if _CHIRP33 and drf[a]>0.5:
            res=res/max(drf[a],1.0)          # chirp-rate mismatch rejects static confusion
        Pt=Psig[k]*(T/T_OBS)
        sig[a]=np.sqrt(max(res,0)/max(Pt,1e-300)+1/max(rp[k],1e-9)**2)
    good=sig<np.pi; ks,sig=ks[good],sig[good]
    if len(ks)==0: return ks,sig,np.zeros(0),np.zeros(0),np.zeros(0)
    nats=np.log(np.pi/sig); o=np.argsort(nats)[::-1][:nmax]
    return ks[o],sig[o],f2[ks[o]],rp[ks[o]],nats[o]

print(f"\n  (c) CHIRP-RATE CONFUSION REJECTION at the booked machine (cap 1300):")
print(f"  {'arm':>26} {'n':>6} {'ungated':>8} {'corner gate':>12} {'purity':>8} {'psr/sky':>8}")
for _lab,_ch in (("booked (no chirp reject)",False),("+ chirp-rate rejection",True)):
    globals()['_CHIRP33']=_ch
    r=corner32(1300)
    s=f"  {_lab:>26} {r['n']:6d} {r['ung']:8.3f}"
    s+= f" {r['C90'][0]:11d}% {r['C90'][1]:8.3f} {r['C90'][2]:8.1f}" if r["C90"] else f" {'none':>12}"
    print(s)
globals()['_CHIRP33']=False

# ================= T2: TIMING-MODEL ABSORPTION =================================
print(f"\n{'='*74}\nJOB 33 -- T2 (R5): TIMING-MODEL ABSORPTION\n{'='*74}")
print("  Projection model: quadratic spin-down absorbs power below ~1/T_span;")
print("  position/proper-motion/parallax notch f = 1/yr and 2/yr.")
print("  Transmission R(f) = (fT)^6/(1+(fT)^6) * notch(1/yr) * notch(2/yr),")
print("  notch width 1/T. ILLUSTRATIVE analytic stand-in for the G-matrix.")
_TM33=False
def _trans(f, T):
    x=(f*T)**6; R=x/(1.0+x)
    for fy in (1.0/YR_S, 2.0/YR_S):
        R=R*(1.0-np.exp(-0.5*((f-fy)*T)**2))
    return np.maximum(R,1e-6)
_inv_t1=inventory20
def inventory20(f2,h2,Psig,Pbkg,sea,T_yr,sig_ns,cad,Np,depth_dex=0.0,rho_cut=3.0,nmax=12):
    if _TM33:
        R=_trans(f2,T_yr*YR_S)
        h2=h2*np.sqrt(R); Psig=Psig*R; Pbkg=Pbkg*R
    return _inv_t1(f2,h2,Psig,Pbkg,sea,T_yr,sig_ns,cad,Np,depth_dex,rho_cut,nmax)
print(f"\n  {'arm':>26} {'n':>6} {'ungated':>8} {'corner gate':>12} {'purity':>8} {'psr/sky':>8}")
_TMR={}
for _lab,_tm in (("booked (no timing model)",False),("+ TIMING-MODEL projection",True)):
    globals()['_TM33']=_tm
    r=corner32(1300); _TMR[_tm]=r
    s=f"  {_lab:>26} {r['n']:6d} {r['ung']:8.3f}"
    s+= f" {r['C90'][0]:11d}% {r['C90'][1]:8.3f} {r['C90'][2]:8.1f}" if r["C90"] else f" {'none':>12}"
    print(s)
_a=_TMR[False]["C90"][2]; _b=_TMR[True]["C90"][2] if _TMR[True]["C90"] else 0.0
_deg=1-_b/max(_a,1e-9)
print(f"    corner {_a:.1f} -> {_b:.1f} psr/sky; degradation {_deg:+.1%}")
print(f"    PRE-REG 0-35% degradation:"
      f" {'HIT -- robust to timing-model fitting' if 0<=_deg<=0.35 else 'MISS -- calendar must be re-derived'}")
globals()['_TM33']=False
print(f"\n{'='*74}\nJOB 33 -- T3 (R16): APERTURE vs THE GEOMETRY GATE\n{'='*74}")
def collect_ap(nsky=60, T=45, cap=1300, M=99, seed0=8000):
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
        for i in range(Np):
            if not rncl[i]: continue
            sg=np.array([depro_sigma(
                np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2,0.0)),
                rp[a],fk[a],i,"loo") for a in range(len(ks))])
            ok=sg<np.pi
            if ok.sum()<2: continue
            sgu=sg[ok]
            mu=np.arccos(rng.uniform(-0.9,0.9,len(sgu)))
            sgu=sgu*np.sqrt((1/3)/np.maximum((1+np.cos(mu))**2/4,1e-3))
            ok2=sgu<np.pi
            if ok2.sum()<2: continue
            sgu=sgu[ok2]; fku=fk[ok][ok2]
            kap=2*np.pi*fku*(1-np.cos(mu[ok2]))/C_PC_S
            Ltrue=float(np.clip(rng.normal(0,pri[i]),-2.9*pri[i],2.9*pri[i]))
            Lg=np.linspace(-3*pri[i],3*pri[i],601)
            eps=kap*Ltrue+rng.normal(0,sgu)
            r_=(kap[:,None]*Lg[None,:]-eps[:,None]+np.pi)%(2*np.pi)-np.pi
            S=-(r_**2/(2*sgu[:,None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
            sigL=1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)),1e-30))
            if sigL>=TOOTH_W: continue
            w=np.exp(S-S.max()); w/=w.sum()
            ip=int(np.argmax(S)); Lh=float(Lg[ip])
            near=np.abs(Lg-Lh)<=TOOTH_W; mtop=float(w[near].sum())
            # per-pulsar surrogate calibration: same geometry, random phases
            sur=np.empty(M)
            for m in range(M):
                e2=rng.uniform(0,2*np.pi,len(sgu))
                r2=(kap[:,None]*Lg[None,:]-e2[:,None]+np.pi)%(2*np.pi)-np.pi
                S2=-(r2**2/(2*sgu[:,None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
                w2=np.exp(S2-S2.max()); w2/=w2.sum()
                i2=int(np.argmax(S2))
                sur[m]=float(w2[np.abs(Lg-Lg[i2])<=TOOTH_W].sum())
            out.append(dict(sky=j,wrong=bool(abs(Lh-Ltrue)>TOOTH_W),mtop=mtop,
                            rank=float((sur<mtop).mean()),npsr=int(rncl.sum()),
                            sig=sgu.copy(),kap=kap.copy(),pri=float(pri[i]),
                            Ltrue=Ltrue,Lh=Lh))
    return out
_AP=collect_ap(nsky=40)
_w=np.array([r["wrong"] for r in _AP],bool); _sk=np.array([r["sky"] for r in _AP])
_rk=np.array([r["rank"] for r in _AP]); _al=np.array([alias_score24(r) for r in _AP])
_tr,_te=_sk<20,_sk>=20; _nte=20
print(f"  n={len(_AP)}, ungated purity {1-_w.mean():.3f}")
print(f"\n  (a) APERTURE at fixed P_F [threshold = per-pulsar surrogate quantile]:")
print(f"  {'P_F':>6} {'n cert':>7} {'purity':>8} {'psr/sky':>9}")
_APR={}
for _pf in (0.10,0.05,0.01):
    sel=_te&(_rk>=1-_pf)
    if sel.sum()<3: print(f"  {_pf:6.2f} {int(sel.sum()):7d} {'--':>8} {'--':>9}"); continue
    p=1-_w[sel].mean(); y=(sel.sum()-_w[sel].sum())/_nte; _APR[_pf]=(p,y)
    print(f"  {_pf:6.2f} {int(sel.sum()):7d} {p:8.3f} {y:9.1f}")
print(f"\n  (b) HEAD-TO-HEAD vs the alias gate [identical realisations]:")
print(f"  {'method':>18} {'matched-yield purity':>21} {'matched-purity yield':>21}")
for _pf,(p_ap,y_ap) in _APR.items():
    n_ap=int((_te&(_rk>=1-_pf)).sum())
    q=100.0*n_ap/max(_te.sum(),1)
    thr=float(np.percentile(_al[_tr],min(q,100)))
    sel_g=_te&(_al<=thr)
    p_g=1-_w[sel_g].mean() if sel_g.sum() else np.nan
    # matched purity: tighten gate until purity >= p_ap
    y_g=0.0
    for qq in np.linspace(1,60,60):
        t2=float(np.percentile(_al[_tr],qq)); s2=_te&(_al<=t2)
        if s2.sum()<3: continue
        if (1-_w[s2].mean())>=p_ap: y_g=max(y_g,(s2.sum()-_w[s2].sum())/_nte)
    print(f"  {'aperture P_F='+str(_pf):>18} {p_ap:21.3f} {y_ap:21.1f}")
    print(f"  {'  alias gate':>18} {p_g:21.3f} {y_g:21.1f}")
print(f"\n  (c) HYBRIDS at P_F=0.05:")
thr5=float(np.percentile(_al[_tr],20))
for _lab,_m in (("aperture only",_te&(_rk>=0.95)),
                ("gate only",_te&(_al<=thr5)),
                ("aperture AND gate",_te&(_rk>=0.95)&(_al<=thr5)),
                ("aperture OR gate",_te&((_rk>=0.95)|(_al<=thr5)))):
    if _m.sum()<3: continue
    print(f"    {_lab:>20}: n {int(_m.sum()):4d} purity {1-_w[_m].mean():.3f}"
          f" psr/sky {(_m.sum()-_w[_m].sum())/_nte:5.1f}")
print(f"\n  (d) 'PEAK STATISTICS FAIL' re-evaluated at booked margins:")
print(f"    AUC of top-tooth mass (higher=better): {auc_of(_rk,_w,lower_is_better=False):.3f}")
print(f"    AUC of alias score  (lower=better)   : {auc_of(_al,_w):.3f}")

print(f"\n{'='*74}\nJOB 33 -- T4: RE-STAMP BATCH (booked machine)\n{'='*74}")
print(f"  R1 -- gap-vs-T: corner and ungated purity vs span at cap 1300:")
print(f"  {'T':>4} {'ungated':>8} {'corner gate':>12} {'purity':>8} {'psr/sky':>8}")
for _T in (34,45,60,75):
    globals()['_CAP32']=1300
    c=collect24d(nsky=40,inject=1.0,prior="gauss",T=_T)
    if not c: continue
    w=np.array([r["wrong"] for r in c],bool); sk=np.array([r["sky"] for r in c])
    asr=np.array([alias_score24(r) for r in c]); tr,te=sk<20,sk>=20
    C=None
    for q in (20,10,7,5,3):
        t=float(np.percentile(asr[tr],q)); s=te&(asr<=t)
        if s.sum()<5: continue
        p=1-w[s].mean()
        if p>=0.90: C=(q,p,(s.sum()-w[s].sum())/20); break
    print(f"  {_T:4d} {1-w.mean():8.3f}"
          + (f" {C[0]:11d}% {C[1]:8.3f} {C[2]:8.1f}" if C else f" {'none':>12}"))
print(f"\n  R4 -- branching ratio at booked margins:")
globals()['_CAP32']=1300
try:
    _b4,_t4,_s4=cascade_b(120)
    print(f"    b = {_b4:.5f}, amplifier 1/(1-b) = {1/(1-_b4):.4f}x, right {_t4}")
    print(f"    PRE-REG amplifier 1.05-1.15x:"
          f" {'HIT' if 1.05<=1/(1-_b4)<=1.15 else 'MISS'}")
except Exception as e:
    print(f"    cascade_b unavailable in snapshot: {type(e).__name__}: {e}")
print(f"\n  R10 -- null PURITY across the booked stack: 0.059 (cap 1300, Job 32 T1c)")
print(f"    signal 0.653 / null 0.059 = 11.13x; null sits at the prior-mass floor.")
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
