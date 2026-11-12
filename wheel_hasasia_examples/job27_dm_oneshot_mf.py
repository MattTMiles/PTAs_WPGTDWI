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
# ================= JOB 27 =======================================================
# Metric rebuilt on the campaign's PERMUTATION FLOOR: the statistic is the
# percentile rank of the real posterior mass within its OWN phase-scrambled
# surrogate distribution. Under the null this is Uniform(0,1) BY CONSTRUCTION
# (mean 0.500), so no normalisation bias is possible -- the defect that made
# Job 26's Bayes-gain metric unquotable (hard 0.912, null 0.833, both < 1.0).
print(f"\n{'='*74}\nJOB 27 -- DM CHANNEL / ONE-SHOT REDO / MF SHORTFALL\n{'='*74}")
_M_SURR = 25

def rank_mass(A, kap, pri, Lt, z, Lg, ph, rg, M=_M_SURR):
    near = _np.abs(Lg - Lt) <= TOOTH_W
    logprior = -Lg**2/(2*pri**2)
    mod = kap[:, None]*Lg[None, :]
    if ph is not None: mod = mod + ph[:, None]
    E = _np.exp(-1j*mod)
    def mass(zz):
        Sc = (2*A[:, None]*_np.real(zz[:, None]*E)).sum(0)
        L = Sc + logprior
        w = _np.exp(L - L.max()); w /= w.sum()
        return float(w[near].sum())
    real = mass(z)
    surr = _np.array([mass(z*_np.exp(1j*rg.uniform(0, 2*_np.pi, len(z))))
                      for _ in range(M)])
    return (float((surr < real).sum()) + 0.5*float((surr == real).sum()))/M

# sanity: the metric's own null must sit at 0.500
_nullranks = []
for j in range(8):
    S = _setup(j)
    if S["ks"] is None or len(S["ks"]) < 2: continue
    for i in _np.where(S["rncl"])[0][:40]:
        c = _combs(S, i)
        if c is None: continue
        sg, kap, sel = c; A = 1.0/_np.maximum(sg, 1e-6); rg = S["rg"]
        Lt = rg.normal(0, S["pri"][i])
        z = (rg.normal(0, 1/_np.sqrt(2), len(A))
             + 1j*rg.normal(0, 1/_np.sqrt(2), len(A)))          # PURE NOISE
        Lg = _np.linspace(-3*S["pri"][i], 3*S["pri"][i], _GRID)
        _nullranks.append(rank_mass(A, kap, S["pri"][i], Lt, z, Lg, None, rg))
print(f"  METRIC CHECK -- rank on pure noise: mean {_np.mean(_nullranks):.4f}"
      f" +/- {_np.std(_nullranks)/_np.sqrt(len(_nullranks)):.4f}"
      f" (n={len(_nullranks)}); EXPECT 0.500 by construction")

# ---- T1: THE DM-CALIBRATION CHANNEL, NULL-CALIBRATED ---------------------------
def dm_run(npatch, amp, shuffle=False, oracle=True, nsky=25, nps=60):
    r1s, r2s = [], []
    for j in range(nsky):
        S = _setup(j, dm_patches=npatch, dm_amp=amp)
        if S["ks"] is None or len(S["ks"]) < 2 or S["patch"] is None: continue
        rg = S["rg"]
        # patch-shuffle null: permute the ASSUMED patch labels only. pri is left
        # exactly as built, so every marginal prior width is preserved.
        patch_used = rg.permutation(S["patch"]) if shuffle else S["patch"]
        idx = _np.where(S["rncl"])[0][:nps]
        p0 = int(idx[_np.argmin(S["pri"][idx])]) if len(idx) else None
        if p0 is None: continue
        bt = S["bias"][S["patch"][p0]]
        bhat = bt if oracle else bt + rg.normal(0, 0.10)
        mates = set(_np.where(patch_used == patch_used[p0])[0].tolist())
        for i in idx:
            if i == p0: continue
            c = _combs(S, i)
            if c is None: continue
            sg, kap, sel = c; A = 1.0/_np.maximum(sg, 1e-6)
            Lt = rg.normal(0, S["pri"][i])
            z = A*_np.exp(1j*(kap*Lt)) + (rg.normal(0, 1/_np.sqrt(2), len(A))
                                          + 1j*rg.normal(0, 1/_np.sqrt(2), len(A)))
            Lg = _np.linspace(-3*S["pri"][i], 3*S["pri"][i], _GRID)
            r1s.append(rank_mass(A, kap, S["pri"][i], Lt, z, Lg, None, rg))
            # after: patch-mates inherit the learned calibration -> tighter prior
            pri2 = S["pri"][i]/(1.0 + bhat) if i in mates else S["pri"][i]
            pri2 = max(pri2, 0.05)
            Lg2 = _np.linspace(-3*pri2, 3*pri2, _GRID)
            r2s.append(rank_mass(A, kap, pri2, Lt, z, Lg2, None, rg))
    if not r1s: return _np.nan, _np.nan, _np.nan, 0
    r1, r2 = float(_np.mean(r1s)), float(_np.mean(r2s))
    se = float(_np.std(_np.array(r2s)-_np.array(r1s))/_np.sqrt(len(r1s)))
    b = (r2-r1)/max(r1-0.5, 1e-6)
    return b, r1, r2, se

print(f"\n  T1 -- DM-calibration channel, RANK-CALIBRATED [b_DM pre-reg 0.05-0.15;"
      f"\n        patch-shuffle b_DM < 0.01]")
print(f"    {'patches':>8} {'DM amp':>7} {'corr frac':>10} {'bias':>10}"
      f" {'rank pre':>9} {'rank post':>10} {'d(rank)':>9} {'b_DM':>8}")
_T1 = {}
for _npatch, _amp, _lab, _shuf, _orc in (
        (20, 0.30, "ORACLE",    False, True),
        (20, 0.30, "estimated", False, False),
        (10, 0.30, "ORACLE",    False, True),
        (20, 0.15, "ORACLE",    False, True),
        (20, 0.0577, "ORACLE",  False, True),   # corr frac 0.50 (defensible)
        (20, 0.0315, "ORACLE",  False, True),   # corr frac 0.30 (defensible)
        (20, 0.30, "SHUFFLE",   True,  True)):
    b, r1, r2, se = dm_run(_npatch, _amp, shuffle=_shuf, oracle=_orc)
    corr = _amp/_np.sqrt(_amp**2 + 0.10**2)
    _T1[(_npatch, _amp, _lab)] = (b, r1, r2, se)
    print(f"    {_npatch:8d} {_amp:7.4f} {corr:10.3f} {_lab:>10}"
          f" {r1:9.4f} {r2:10.4f} {r2-r1:+9.4f} {b:8.4f}")
_bo = _T1[(20, 0.30, "ORACLE")][0]; _bs = _T1[(20, 0.30, "SHUFFLE")][0]
print(f"    headline b_DM {_bo:.4f} -> PRE-REG 0.05-0.15:"
      f" {'HIT' if 0.05 <= _bo <= 0.15 else 'MISS'}")
print(f"    patch-shuffle b_DM {_bs:.4f} -> PRE-REG < 0.01:"
      f" {'PASS -- spatial correlation is doing the work' if abs(_bs) < 0.01 else 'FAIL -- channel is a BOOKKEEPING ARTEFACT'}")
print(f"    corr-frac note: 0.949 at amp 0.30 is the TOY mapping"
      f" amp/sqrt(amp^2+0.1^2);\n      defensible rungs (corr 0.30/0.50) reported above."
      f" NO measured\n      NE2001-vs-YMW16 residual set is in my record -- not fabricating one.")
# CEILING AUDIT: if ranks pile at 1.0, the metric cannot register improvement and
# a null b_DM would be an artefact of saturation, not evidence against the channel.
def dm_run2(npatch, amp, shuffle=False, oracle=True, nsky=25, nps=60, M=99):
    r1s, r2s = [], []
    for j in range(nsky):
        S = _setup(j, dm_patches=npatch, dm_amp=amp)
        if S["ks"] is None or len(S["ks"]) < 2 or S["patch"] is None: continue
        rg = S["rg"]
        patch_used = rg.permutation(S["patch"]) if shuffle else S["patch"]
        idx = _np.where(S["rncl"])[0][:nps]
        if not len(idx): continue
        p0 = int(idx[_np.argmin(S["pri"][idx])])
        bt = S["bias"][S["patch"][p0]]
        bhat = bt if oracle else bt + rg.normal(0, 0.10)
        mates = set(_np.where(patch_used == patch_used[p0])[0].tolist())
        for i in idx:
            if i == p0 or i not in mates: continue      # MATES ONLY: where the
            c = _combs(S, i)                            # channel can act at all
            if c is None: continue
            sg, kap, sel = c; A = 1.0/_np.maximum(sg, 1e-6)
            Lt = rg.normal(0, S["pri"][i])
            z = A*_np.exp(1j*(kap*Lt)) + (rg.normal(0,1/_np.sqrt(2),len(A))
                                          + 1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            Lg = _np.linspace(-3*S["pri"][i], 3*S["pri"][i], _GRID)
            r1s.append(rank_mass(A, kap, S["pri"][i], Lt, z, Lg, None, rg, M=M))
            pri2 = max(S["pri"][i]/(1.0+bhat), 0.05)
            Lg2 = _np.linspace(-3*pri2, 3*pri2, _GRID)
            r2s.append(rank_mass(A, kap, pri2, Lt, z, Lg2, None, rg, M=M))
    return _np.array(r1s), _np.array(r2s)

_r1, _r2 = dm_run2(20, 0.30)
_sat = float((_r1 >= 0.999).mean())
print(f"\n  CEILING AUDIT [M=99 surrogates, patch-mates only, n={len(_r1)}]:")
print(f"    fraction of pulsars at rank 1.000 (saturated): {_sat:.3f}")
print(f"    mean rank pre {_r1.mean():.4f} post {_r2.mean():.4f}"
      f" d {_r2.mean()-_r1.mean():+.4f}")
_un = _r1 < 0.999
if _un.sum() > 20:
    d = _r2[_un]-_r1[_un]
    print(f"    UNSATURATED subset only (n={int(_un.sum())}): rank"
          f" {_r1[_un].mean():.4f} -> {_r2[_un].mean():.4f}"
          f" d {d.mean():+.4f} +/- {d.std()/_np.sqrt(len(d)):.4f}"
          f" -> b_DM {(_r2[_un].mean()-_r1[_un].mean())/max(_r1[_un].mean()-0.5,1e-6):+.4f}")
print(f"    prior-width change actually applied to mates: median"
      f" {_np.median([1/(1+0.30)]):.3f}x (oracle bias 0.30 config)")
# ---- T1 DECISIVE: enough n to EXCLUDE the pre-registered band, not just miss it
_r1, _r2 = dm_run2(20, 0.30, nsky=150, M=99)
_d = _r2 - _r1; _se = float(_d.std()/_np.sqrt(len(_d)))
_bse = _se/max(_r1.mean()-0.5, 1e-6)
_bhat = float(_d.mean())/max(_r1.mean()-0.5, 1e-6)
print(f"\n  T1 DECISIVE [n={len(_d)} patch-mates, M=99]:")
print(f"    rank {_r1.mean():.4f} -> {_r2.mean():.4f}; d(rank) {_d.mean():+.4f}"
      f" +/- {_se:.4f}")
print(f"    b_DM = {_bhat:+.4f} +/- {_bse:.4f};  2-sigma upper limit"
      f" {_bhat+2*_bse:+.4f}")
print(f"    PRE-REG 0.05-0.15:"
      f" {'EXCLUDED (2-sigma upper below band)' if _bhat+2*_bse < 0.05 else 'not excluded -- n still insufficient'}")

# ---- T2: ONE-SHOT REDO, SPAN BUG FIXED -----------------------------------------
# Job 25 bug: T was accepted but never used -- span hardcoded at log(25.0), so the
# T=45 and T=60 rows were byte-identical. Fixed here, WITH a pre-flight invariant.
def oneshot27(T=60, nsky=25, mcx=1.0, complexch=True, rweight=True):
    ok = n = wrong = 0; svals = []
    rg = _np.random.default_rng(2700)
    for j in range(nsky):
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        if mcx != 1.0: h2 = h2*mcx**(5.0/3.0)
        Np = 120
        Ti = _np.clip(_np.exp(rg.normal(_np.log(T), 0.5, Np)), 5, 90)   # <-- T USED
        beta = _np.arcsin(rg.uniform(-1,1,Np)); Dk = rg.uniform(0.5,3.0,Np)
        px = 2.34/_np.cos(beta)**2*Dk**2*10.0*_np.sqrt(30.0/Ti)
        pb = 3.0*(Ti/30.0)**(-2.5)
        pri = _np.maximum(_np.where(rg.random(Np)<0.7, _np.minimum(px,pb), px), 0.05)
        Arn,_g = rn_draw(Np); rncl = Arn == 0.0
        ks, sig, fk, rp, nats = inventory20(f2, h2, P0, P0, sea,
                                            float(_np.mean(Ti)), 100.0, 14.0, Np)
        if len(ks) < 2: continue
        for i in _np.where(rncl)[0]:
            sgb = _np.array([depro_sigma(
                _np.sqrt(max(sig[a]**2-1/max(rp[a],1e-9)**2,0.0)),
                rp[a], fk[a], i, "loo") for a in range(len(ks))])
            sgb = sgb*_np.sqrt(_np.mean(Ti)/Ti[i])
            th = _np.arccos(rg.uniform(-0.9,0.9,len(ks)))
            if rweight:
                sgb = sgb*_np.sqrt((1/3)/_np.maximum((1+_np.cos(th))**2/4, 1e-3))
            okm = sgb < _np.pi
            if okm.sum() < 2: continue
            sgu, fku = sgb[okm], fk[okm]
            svals.append(1.0/(1.0 + float(_np.median(sgu))**2))
            mu = _np.arccos(rg.uniform(-0.9,0.9,len(sgu)))
            kap = 2*_np.pi*fku*(1-_np.cos(mu))/C_PC_S
            Lt = rg.normal(0, pri[i])
            Lg = _np.linspace(-3*pri[i], 3*pri[i], 401)
            A = 1.0/_np.maximum(sgu,1e-6)
            z = A*_np.exp(1j*(kap*Lt)) + (rg.normal(0,1/_np.sqrt(2),len(A))
                                          + 1j*rg.normal(0,1/_np.sqrt(2),len(A)))
            mod = kap[:,None]*Lg[None,:]
            if complexch:
                Sc = (2*A[:,None]*_np.real(z[:,None]*_np.exp(-1j*mod))).sum(0)
            else:
                eps = _np.angle(z)
                rr = (mod-eps[:,None]+_np.pi)%(2*_np.pi)-_np.pi
                Sc = -(rr**2/(2*sgu[:,None]**2)).sum(0)
            Sc = Sc - Lg**2/(2*pri[i]**2)
            if 1.0/_np.sqrt(max(float(_np.sum(kap**2/sgu**2)),1e-30)) < TOOTH_W:
                n += 1
                if abs(float(Lg[int(_np.argmax(Sc))])-Lt) <= TOOTH_W: ok += 1
                else: wrong += 1
    return ok, wrong, n, float(_np.median(svals)) if svals else _np.nan

print(f"\n  T2 -- ONE-SHOT REDO [Job 25 Task 5 was VOID: span hardcoded]")
_rows = []
for _lab, _T, _mcx in (("T=45, all levers", 45, 1.0),
                       ("T=60, all levers", 60, 1.0),
                       ("T=60, + concentration Mc x2", 60, 2.0)):
    _rows.append((_lab, _T, _mcx) + oneshot27(T=_T, mcx=_mcx))
# PRE-FLIGHT INVARIANT: s_med MUST differ across T before any result is read
_s45 = [r[6] for r in _rows if r[1] == 45][0]
_s60 = [r[6] for r in _rows if r[1] == 60 and r[2] == 1.0][0]
_gate = abs(_s45-_s60)/max(_s45, 1e-9) > 0.01
print(f"    PRE-FLIGHT GATE: s_med(T=45) {_s45:.4f} vs s_med(T=60) {_s60:.4f}"
      f" -> rel diff {abs(_s45-_s60)/max(_s45,1e-9):.4f}"
      f" -> {'PASS, span is live' if _gate else 'FAIL -- span STILL not propagating; results NOT read'}")
if _gate:
    print(f"    {'config':>30} {'s_med':>7} {'n':>6} {'right':>7} {'purity':>8}"
          f" {'frac RN-clean':>14}")
    for _lab, _T, _mcx, ok, wrong, n, sm in _rows:
        print(f"    {_lab:>30} {sm:7.4f} {n:6d} {ok:7d} {ok/max(n,1):8.3f}"
              f" {ok/max(n,1)*n/max(n,1):14.3f}")
    _p60 = [r for r in _rows if r[1]==60 and r[2]==1.0][0]
    _pf = _p60[3]/max(_p60[5],1)
    print(f"    turn-1 right fraction (T=60, no concentration) {_pf:.3f}"
          f" -> PRE-REG 0.20-0.40: {'HIT' if 0.20 <= _pf <= 0.40 else 'MISS'}")

# ---- T3: COMPLEX-MF SHORTFALL ---------------------------------------------------
# Pre-check by INSPECTION: Job 25 Task 2 set A = 1/max(sigma) as BOTH the generating
# amplitude and the template. The machine was ALREADY in the oracle-A condition, so
# amplitude-template error cannot explain the 1.037x. Confirmed by measurement below,
# then the real candidate is tested: sigma HETEROGENEITY (phase-only already weights
# by 1/sigma^2, so it captures most of the amplitude information the MF adds).
print(f"\n  T3 -- COMPLEX-MF SHORTFALL (1.037x measured vs 1.17x toy)")
def mf_cmp(nsky=25, T=45, sA=0.0, uniform_sigma=False):
    ph_ok = cx_ok = n = 0
    rg = _np.random.default_rng(2500)
    for r in collect24d(nsky=nsky, T=T, s15=0.35, inject=1.0, prior="gauss"):
        sgu, kap, pri, Lt = r["sig"], r["kap"], r["pri"], r["Ltrue"]
        if uniform_sigma: sgu = _np.full_like(sgu, float(_np.mean(sgu)))
        A = 1.0/_np.maximum(sgu, 1e-6)
        z = A*_np.exp(1j*(kap*Lt)) + (rg.normal(0,1/_np.sqrt(2),len(A))
                                      + 1j*rg.normal(0,1/_np.sqrt(2),len(A)))
        At = A*_np.exp(rg.normal(0, sA, len(A))) if sA > 0 else A   # template error
        Lg = _np.linspace(-3*pri, 3*pri, 401)
        mod = kap[:,None]*Lg[None,:]
        eps = _np.angle(z)
        rr = (mod-eps[:,None]+_np.pi)%(2*_np.pi)-_np.pi
        Sp = -(rr**2/(2*sgu[:,None]**2)).sum(0) - Lg**2/(2*pri**2)
        Sx = (2*At[:,None]*_np.real(z[:,None]*_np.exp(-1j*mod))).sum(0) - Lg**2/(2*pri**2)
        n += 1
        if abs(float(Lg[int(_np.argmax(Sp))])-Lt) <= TOOTH_W: ph_ok += 1
        if abs(float(Lg[int(_np.argmax(Sx))])-Lt) <= TOOTH_W: cx_ok += 1
    return ph_ok/max(n,1), cx_ok/max(n,1), n
print(f"    {'condition':>28} {'phase-only':>11} {'complex MF':>11} {'gain':>7}")
for _lab, _sA, _uni in (("oracle A (as banked)", 0.0, False),
                        ("A template err 20%", 0.2, False),
                        ("A template err 40%", 0.4, False),
                        ("A template err 80%", 0.8, False),
                        ("UNIFORM sigma, oracle A", 0.0, True)):
    p, c, n = mf_cmp(sA=_sA, uniform_sigma=_uni)
    print(f"    {_lab:>28} {p:11.4f} {c:11.4f} {c/max(p,1e-9):7.3f}x   (n={n})")
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
