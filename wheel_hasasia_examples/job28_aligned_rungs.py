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
# ================================================================================
# T1 REDO -- the DM channel PROPERLY SPECIFIED.
# The Job 26 construction put the shared DM term in the prior WIDTH and drew truth
# from the widened prior, so "learning" it narrowed the prior below the true spread
# -- a mis-specification, which the rank metric correctly penalises (b = -0.043).
# A DM-calibration channel is a shared OFFSET in L: patch-mates share a distance
# bias that one certified pulsar can measure and hand to its neighbours. Prior
# WIDTH is untouched; only the prior CENTRE moves. Grid coverage is held identical
# before/after so the comparison cannot be a windowing artefact.
# ================================================================================
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

def dm_offset(npatch=20, off_pc=0.5, shuffle=False, oracle=True, nsky=150, nps=60):
    r1s, r2s = [], []
    for j in range(nsky):
        S = _setup(j)                      # NO width inflation
        if S["ks"] is None or len(S["ks"]) < 2: continue
        rg = S["rg"]
        patch = rg.integers(0, npatch, S["Np"])
        off = rg.normal(0, off_pc, npatch)          # shared per-patch distance bias
        patch_used = rg.permutation(patch) if shuffle else patch
        idx = _np.where(S["rncl"])[0][:nps]
        if not len(idx): continue
        p0 = int(idx[_np.argmin(S["pri"][idx])])
        oh = off[patch[p0]] if oracle else off[patch[p0]] + rg.normal(0, 0.15)
        mates = set(_np.where(patch_used == patch_used[p0])[0].tolist())
        for i in idx:
            if i == p0 or i not in mates: continue
            c = _combs(S, i)
            if c is None: continue
            sg, kap, sel = c; A = 1.0/_np.maximum(sg, 1e-6)
            pri = S["pri"][i]
            Lt = rg.normal(0, pri) + off[patch[i]]   # truth carries the offset
            z = A*_np.exp(1j*(kap*Lt)) + (rg.normal(0,1/_np.sqrt(2),len(A))
                                          + 1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            # IDENTICAL grid both arms -> no windowing artefact
            span = 3*pri + 3*off_pc
            Lg = _np.linspace(-span, span, _GRID)
            r1s.append(rank_mass_mu(A, kap, pri, Lt, z, Lg, 0.0, rg))   # uncalibrated
            r2s.append(rank_mass_mu(A, kap, pri, Lt, z, Lg, oh,  rg))   # calibrated
    r1s, r2s = _np.array(r1s), _np.array(r2s)
    if not len(r1s): return (_np.nan,)*5 + (0,)
    d = r2s-r1s; se = float(d.std()/_np.sqrt(len(d)))
    b = float(d.mean())/max(float(r1s.mean())-0.5, 1e-6)
    bse = se/max(float(r1s.mean())-0.5, 1e-6)
    return b, bse, float(r1s.mean()), float(r2s.mean()), float(d.mean()), len(d)

print(f"\n  T1 REDO -- DM channel as a shared OFFSET (properly specified):")
print(f"    {'off (pc)':>9} {'bias':>10} {'rank pre':>9} {'rank post':>10}"
      f" {'d(rank)':>9} {'b_DM':>9} {'+/-':>7} {'n':>5}")
_OFF = {}
for _op, _sh, _orc, _lab in ((0.50, False, True,  "ORACLE"),
                             (0.50, False, False, "estimated"),
                             (1.00, False, True,  "ORACLE"),
                             (0.25, False, True,  "ORACLE"),
                             (0.50, True,  True,  "SHUFFLE")):
    b, bse, r1, r2, d, n = dm_offset(off_pc=_op, shuffle=_sh, oracle=_orc)
    _OFF[(_op, _lab)] = (b, bse, n)
    print(f"    {_op:9.2f} {_lab:>10} {r1:9.4f} {r2:10.4f} {d:+9.4f}"
          f" {b:+9.4f} {bse:7.4f} {n:5d}")
_bh, _bhse, _ = _OFF[(0.50, "ORACLE")]; _bs2, _, _ = _OFF[(0.50, "SHUFFLE")]
print(f"    headline b_DM {_bh:+.4f} +/- {_bhse:.4f} -> PRE-REG 0.05-0.15:"
      f" {'HIT' if 0.05 <= _bh <= 0.15 else 'MISS'}")
print(f"    patch-shuffle b_DM {_bs2:+.4f} -> PRE-REG < 0.01:"
      f" {'PASS -- spatial correlation does the work' if abs(_bs2) < 0.01 else 'FAIL -- bookkeeping artefact'}")

# ---- T2 column fix: frac of RN-clean needs the RN-clean DENOMINATOR -------------
def oneshot27b(T=60, nsky=25, mcx=1.0):
    ok = n = wrong = rnc = 0; svals = []
    rg = _np.random.default_rng(2700)
    for j in range(nsky):
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        if mcx != 1.0: h2 = h2*mcx**(5.0/3.0)
        Np = 120
        Ti = _np.clip(_np.exp(rg.normal(_np.log(T), 0.5, Np)), 5, 90)
        beta = _np.arcsin(rg.uniform(-1,1,Np)); Dk = rg.uniform(0.5,3.0,Np)
        px = 2.34/_np.cos(beta)**2*Dk**2*10.0*_np.sqrt(30.0/Ti)
        pb = 3.0*(Ti/30.0)**(-2.5)
        pri = _np.maximum(_np.where(rg.random(Np)<0.7, _np.minimum(px,pb), px), 0.05)
        Arn,_g = rn_draw(Np); rncl = Arn == 0.0
        ks, sig, fk, rp, nats = inventory20(f2, h2, P0, P0, sea,
                                            float(_np.mean(Ti)), 100.0, 14.0, Np)
        if len(ks) < 2: continue
        rnc += int(rncl.sum())                       # <-- the DENOMINATOR
        for i in _np.where(rncl)[0]:
            sgb = _np.array([depro_sigma(
                _np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2,0.0)),
                rp[a], fk[a], i, "loo") for a in range(len(ks))])
            sgb = sgb*_np.sqrt(_np.mean(Ti)/Ti[i])
            th = _np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            sgb = sgb*_np.sqrt((1/3)/_np.maximum((1+_np.cos(th))**2/4, 1e-3))
            okm = sgb < _np.pi
            if okm.sum() < 2: continue
            sgu, fku = sgb[okm], fk[okm]
            svals.append(1.0/(1.0 + float(_np.median(sgu))**2))
            mu = _np.arccos(rg.uniform(-0.9,0.9,len(sgu)))
            kap = 2*_np.pi*fku*(1-_np.cos(mu))/C_PC_S
            Lt = rg.normal(0, pri[i]); Lg = _np.linspace(-3*pri[i], 3*pri[i], 401)
            A = 1.0/_np.maximum(sgu,1e-6)
            z = A*_np.exp(1j*(kap*Lt)) + (rg.normal(0,1/_np.sqrt(2),len(A))
                                          + 1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            Sc = (2*A[:,None]*_np.real(z[:,None]*_np.exp(
                -1j*(kap[:,None]*Lg[None,:])))).sum(0) - Lg**2/(2*pri[i]**2)
            if 1.0/_np.sqrt(max(float(_np.sum(kap**2/sgu**2)),1e-30)) < TOOTH_W:
                n += 1
                if abs(float(Lg[int(_np.argmax(Sc))])-Lt) <= TOOTH_W: ok += 1
                else: wrong += 1
    return ok, wrong, n, rnc, float(_np.median(svals)) if svals else _np.nan
print(f"\n  T2 CORRECTED TABLE (frac of RN-clean now uses the RN-clean denominator):")
print(f"    {'config':>30} {'s_med':>7} {'n cert':>7} {'right':>6} {'purity':>7}"
      f" {'RN-clean':>9} {'right/RNc':>10}")
for _lab, _T, _mcx in (("T=45, all levers", 45, 1.0), ("T=60, all levers", 60, 1.0),
                       ("T=60, + concentration Mc x2", 60, 2.0)):
    ok, wr, n, rnc, sm = oneshot27b(T=_T, mcx=_mcx)
    print(f"    {_lab:>30} {sm:7.4f} {n:7d} {ok:6d} {ok/max(n,1):7.3f}"
          f" {rnc:9d} {ok/max(rnc,1):10.3f}")
def combs28b(S, i, lo=-0.9, hi=0.9, couple=False):
    rg = S["rg"]; ks, sig, fk, rp = S["ks"], S["sig"], S["fk"], S["rp"]
    if len(ks) < 2: return None
    sg = _np.array([depro_sigma(
        _np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2, 0.0)),
        rp[a], fk[a], i, "loo") for a in range(len(ks))])
    sg = sg*_np.sqrt(_np.mean(S["Ti"])/S["Ti"][i])
    cm = rg.uniform(lo, hi, len(ks)); mu = _np.arccos(cm)
    th = mu if couple else _np.arccos(rg.uniform(-0.9, 0.9, len(ks)))
    sg = sg*_np.sqrt((1/3)/_np.maximum((1+_np.cos(th))**2/4, 1e-3))
    ok = sg < _np.pi
    if ok.sum() < 2: return None
    return sg[ok], 2*_np.pi*fk[ok]*(1-_np.cos(mu[ok]))/C_PC_S, _np.where(ok)[0]

def t1b(lo, hi, couple, nsky=40, nps=60, M=99):
    ranks, pur, Kmin, ncert, nrn, coarse, nK = [], [], [], 0, 0, 0, 0
    for j in range(nsky):
        S = _setup(j)
        if S["ks"] is None or len(S["ks"]) < 2: continue
        rg = S["rg"]; idx = _np.where(S["rncl"])[0][:nps]; nrn += len(idx)
        for i in idx:
            c = combs28b(S, i, lo, hi, couple)
            if c is None: continue
            sg, kap, sel = c; A = 1.0/_np.maximum(sg, 1e-6); pri = S["pri"][i]
            Kc = 6*pri/(2*_np.pi/_np.maximum(kap, 1e-30))
            Kmin.append(float(Kc.min())); nK += 1
            if Kc.min() <= 2.0: coarse += 1
            Lt = rg.normal(0, pri)
            z = A*_np.exp(1j*(kap*Lt)) + (rg.normal(0,1/_np.sqrt(2),len(A))
                                          + 1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            Lg = _np.linspace(-3*pri, 3*pri, _GRID)
            Sc = (2*A[:,None]*_np.real(z[:,None]*_np.exp(
                -1j*(kap[:,None]*Lg[None,:])))).sum(0) - Lg**2/(2*pri**2)
            hit = abs(float(Lg[int(_np.argmax(Sc))])-Lt) <= TOOTH_W
            if 1.0/_np.sqrt(max(float(_np.sum(kap**2/sg**2)),1e-30)) < TOOTH_W:
                ncert += 1; pur.append(1.0 if hit else 0.0)
            ranks.append(rank_mass_mu(A, kap, pri, Lt, z, Lg, 0.0, rg, M=M))
    p = float(_np.mean(pur)) if pur else _np.nan
    return dict(rank=float(_np.mean(ranks)), pur=p, y=ncert/max(nrn,1),
                e2e=ncert*p/max(nrn,1), Kmin=float(_np.median(Kmin)),
                cf=coarse/max(nK,1))
print(f"\n  T1b -- ALIGNED TAIL ISOLATED (anti-aligned cut held at -0.9):")
print(f"  {'variant':>36} {'K min':>7} {'coarse':>7} {'rank':>7} {'purity':>7}"
      f" {'yield':>7} {'right/RNc':>10}")
_R = {}
for _lab, _lo, _hi, _cp in (
        ("BANKED (-0.9,0.9) decoupled",  -0.9, 0.9, False),
        ("matched baseline, coupled",    -0.9, 0.9, True),
        ("ALIGNED TAIL IN (-0.9,1.0)",   -0.9, 1.0, True),
        ("aligned tail in, decoupled",   -0.9, 1.0, False)):
    r = t1b(_lo, _hi, _cp); _R[_lab] = r
    print(f"  {_lab:>36} {r['Kmin']:7.2f} {r['cf']:7.3f} {r['rank']:7.4f}"
          f" {r['pur']:7.3f} {r['y']:7.3f} {r['e2e']:10.3f}")
_bb, _aa = _R["matched baseline, coupled"], _R["ALIGNED TAIL IN (-0.9,1.0)"]
_kc = _bb['Kmin']/max(_aa['Kmin'],1e-9); _dy = _aa['y']/max(_bb['y'],1e-9)-1
print(f"\n    K cut (matched) {_bb['Kmin']:.2f} -> {_aa['Kmin']:.2f} = {_kc:.2f}x"
      f"  [PRE-REG 2-5x: {'HIT' if 2 <= _kc <= 5 else 'MISS'}]")
print(f"    gated yield {_dy:+.1%}  [PRE-REG +20-60%: {'HIT' if 0.20 <= _dy <= 0.60 else 'MISS'}]")
print(f"    end-to-end right/RN-clean {_bb['e2e']:.3f} -> {_aa['e2e']:.3f}"
      f" = {_aa['e2e']/max(_bb['e2e'],1e-9)-1:+.1%}")
print(f"    NOTE: banked K_min is ALREADY {_R['BANKED (-0.9,0.9) decoupled']['Kmin']:.2f}"
      f" -- the cut boundary (cos=0.9) is itself near-unambiguous, so the")
print(f"    headroom the brief predicts is largely BANKED ALREADY (the falsifier).")
