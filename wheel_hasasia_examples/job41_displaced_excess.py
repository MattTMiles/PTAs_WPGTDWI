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
print(f"\n{'='*74}\nJOB 41 -- T2: THE DISPLACED AUTO-EXCESS\n{'='*74}")
print("  autos = Earth-term + pulsar-term power; crosses carry Earth terms ONLY")
print("  (pulsar terms are private). auto - cross_implied = the pulsar-term")
print("  population, placed at RETARDED frequencies f_P.")
_TSUN=4.925490947e-6; _MSUN=1.989e30; _FE=0.485
_ordA2=np.argsort(SRC_AMP); _AS2=SRC_AMP[_ordA2]
def mcz(a):
    k=np.clip(np.searchsorted(_AS2,a),0,len(_AS2)-1); i=_ordA2[k]
    return SRC_MC[i],SRC_Z[i]
def tcoal(f,ts): return (5.0/256.0)*ts**(-5.0/3.0)*(np.pi*f)**(-8.0/3.0)
def excess_spectrum(j,Np=60,null=False,seed=0):
    sd,f2,h2,P0,sea=_pop18(j%50,0.35)
    amp=h2/np.maximum(f2,1e-30)**(2.0/3.0)
    mc,z=mcz(amp); ts=_TSUN*mc*(1.0+z)/_MSUN
    tc=tcoal(np.maximum(f2,1e-30),ts)
    rg=np.random.default_rng(4800+j*17+seed)
    L=rg.uniform(0.5,3.0,Np)*1000.0
    cmu=rg.uniform(-1,1,(Np,len(f2)))
    tau=L[:,None]*(1-cmu)/C_PC_S
    fP=f2[None,:] if null else f2[None,:]*(1.0+tau/np.maximum(tc[None,:],1e-9))**(-3.0/8.0)
    hP=amp[None,:]*np.maximum(fP,1e-30)**(2.0/3.0)
    Ppt=(1.0-_FE)*hP**2*np.maximum(fP,1e-30)*T_OBS          # pulsar-term power
    ok=(fP>2e-10)&(fP<1e-7)
    idx=np.clip(np.digitize(fP,bin_edges)-1,0,len(bin_cent)-1)
    spec=np.zeros(len(bin_cent))
    np.add.at(spec,idx[ok],Ppt[ok])
    return spec/Np
IN=(bin_cent>=15e-9)&(bin_cent<=40e-9); AB=(bin_cent>50e-9)
def ratio_of(sp): return sp[IN].sum()/max(sp[AB].sum(),1e-300)
NSKY=20
sig=np.mean([excess_spectrum(j) for j in range(NSKY)],axis=0)
nul=np.mean([excess_spectrum(j,null=True) for j in range(NSKY)],axis=0)
Rs, Rn = ratio_of(sig), ratio_of(nul)
print(f"\n  in-band(15-40 nHz) / above-50-nHz auto-excess power:")
print(f"    SIGNAL (retarded pulsar terms) : {Rs:10.3f}")
print(f"    NULL   (tau=0, terms at f_E)   : {Rn:10.3f}")
print(f"    displacement factor            : {Rs/max(Rn,1e-30):10.3f}x")
print(f"    PRE-REG excess concentrated in 15-40 regardless of f_E, vs null"
      f" tracking f_E (~2x): {'HIT' if Rs/max(Rn,1e-30)>2 else 'MISS'}")
sur=[ratio_of(excess_spectrum(j%NSKY,null=True,seed=100+m)) for m in range(99) for j in [m%NSKY]]
sur=np.array(sur)
rank=float((sur<Rs).mean())
print(f"\n  DETECTABILITY vs M=99 scrambled-retardation surrogates:")
print(f"    surrogate ratio: median {np.median(sur):.3f}, 10/90"
      f" [{np.quantile(sur,.1):.3f},{np.quantile(sur,.9):.3f}]")
print(f"    signal rank among surrogates: {rank:.4f}"
      f"  -> {'DETECTABLE (p < 0.01)' if rank>=0.99 else 'not separable'}")
zs=(Rs-sur.mean())/max(sur.std(),1e-30)
print(f"    z = {zs:+.1f} sigma against the surrogate ensemble")
with open(OUT+"fig41_displaced_excess.csv","w",newline="") as fh:
    cw=csv.writer(fh); cw.writerow(["freq_Hz","excess_signal","excess_null_tau0"])
    for i in range(len(bin_cent)):
        cw.writerow([f"{bin_cent[i]:.4e}",f"{sig[i]:.6e}",f"{nul[i]:.6e}"])
fig,ax=plt.subplots(figsize=(6.8,4.2))
ax.axvspan(15,40,color="0.88",zorder=0)
ax.text(27,ax.get_ylim()[1],"",ha="center")
ax.loglog(bin_cent*1e9,np.maximum(sig,1e-40),color="C0",lw=1.8)
ax.loglog(bin_cent*1e9,np.maximum(nul,1e-40),color="C3",lw=1.4,ls="--")
ax.set_xlim(1,100)
m=np.maximum(sig,1e-40).max()
ax.set_ylim(m*1e-6,m*3)
ax.text(22,m*1.2,"15-40 nHz",fontsize=8,color="0.35",ha="center")
ax.text(60,m*1e-2,"null: terms at $f_E$",fontsize=8,color="C3")
ax.text(6,m*0.3,"signal: retarded",fontsize=8,color="C0")
ax.set_xlabel("frequency (nHz)"); ax.set_ylabel("auto-minus-cross excess power")
ax.grid(alpha=.25,lw=.5,which="both")
plt.tight_layout(); plt.savefig(OUT+"fig41_displaced_excess.png",dpi=150); plt.close()
print(f"  wrote fig41_displaced_excess.png + .csv")
