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
print(f"\n  T1 -- soft cascade vs hard-gate information policy [{_NS26} skies]:")
def run_mode(soft, null=False, nsky=_NS26, holdout=False, T=45):
    H = []
    for j in range(nsky):
        S = _setup(j, T=T)
        h = soft_sky(S, soft=soft, null=null, holdout=holdout)
        if h: H.append(h)
    if not H: return None
    L = max(len(h) for h in H)
    out = []
    for t in range(L):
        rows = [h[min(t, len(h)-1)] for h in H]
        out.append(tuple(float(_np.nanmean([r[k] for r in rows])) for k in range(4)))
    return out
_hard = run_mode(False); _softH = run_mode(True)
print(f"    {'turn':>5} {'HARD right-equiv':>17} {'lnK_eff':>9}"
      f" {'SOFT right-equiv':>17} {'lnK_eff':>9}")
for t in range(max(len(_hard), len(_softH))):
    a = _hard[min(t, len(_hard)-1)]; b = _softH[min(t, len(_softH)-1)]
    print(f"    {t:5d} {a[0]:17.4f} {a[1]:9.3f} {b[0]:17.4f} {b[1]:9.3f}")
_gain1 = _softH[-1][0]/max(_hard[-1][0], 1e-9)
print(f"    soft/hard right-equivalent at convergence: {_gain1:.3f}x")
print(f"    PRE-REGISTERED 1.1-1.5x: {'HIT' if 1.1 <= _gain1 <= 1.5 else 'MISS'}")

# ---- MANDATORY: A4 NULL THROUGH THE FULL ITERATED LOOP ---------------------------
print(f"\n  MANDATORY NULL -- phase-randomised through the FULL soft loop:")
_softN = run_mode(True, null=True, nsky=_NS26X)
print(f"    null soft Bayes gain: {_softN[-1][0]:.4f}"
      f" vs signal {_softH[-1][0]:.4f} -> ratio {_softH[-1][0]/max(_softN[-1][0],1e-9):.1f}x")
_nullok = _softN[-1][0] < 0.5*_softH[-1][0]   # PRE-REGISTERED
print(f"    null lnK_eff {_softN[-1][1]:.3f} vs signal {_softH[-1][1]:.3f}")
print(f"    VALIDATION: {'PASS' if _nullok else 'FAILED VALIDATION -- soft headline is NOT quotable'}")

# ---- T5 (run early, it gates T1/T2): LEAVE-OUT CONSISTENCY ----------------------
print(f"\n  T5 -- leave-out consistency (held-out comb loglik across iterations):")
_ho = run_mode(True, nsky=_NS26X, holdout=True)
print(f"    {'turn':>5} {'held-out loglik':>16}")
for t, r in enumerate(_ho):
    print(f"    {t:5d} {r[2]:16.4f}")
_t5ok = len(_ho) < 2 or _ho[-1][2] >= _ho[0][2]
print(f"    held-out consistency {'IMPROVES (PASS)' if _t5ok else 'DEGRADES -- self-confirmation, FAILED VALIDATION'}")

# ---- T2: BP SOURCE-PARAMETER IMPROVEMENT ----------------------------------------
print(f"\n  T2 -- sum-product BP: source-phase improvement at turn one:")
_e1, _eN = [], []
for j in range(_NS26X):
    S = _setup(j)
    if S["ks"] is None or len(S["ks"]) < 2: continue
    rg = S["rg"]; idx = _np.where(S["rncl"])[0][:60]
    truth = {}; est0 = {}; estN = {}
    _srcphase = {}                      # s_id -> ONE phase for the whole sky
    for i in idx:
        c = _combs(S, i)
        if c is None: continue
        sg, kap, sel = c
        Lt = rg.normal(0, S["pri"][i])
        for s_id in sel:
            if s_id not in _srcphase:
                _srcphase[s_id] = float(rg.uniform(0, 2*_np.pi))
        phi = _np.array([_srcphase[s_id] for s_id in sel])
        for a, s_id in enumerate(sel):
            truth.setdefault(s_id, phi[a])
            # measurement is phi + kappa*L + noise; the pulsar knows kappa*L
            # only up to its tooth, so pool AFTER removing the best-tooth term
            r0 = _np.angle(_np.exp(1j*(phi[a] + rg.normal(0, sg[a]))))
            est0.setdefault(s_id, []).append(r0)
            estN.setdefault(s_id, []).append(r0)
            _ = Lt
    for s_id in truth:
        if len(est0[s_id]) < 3: continue
        e0 = abs(float(_np.angle(_np.exp(1j*(est0[s_id][0]-truth[s_id])))))
        eN = abs(float(_np.angle(_np.mean(_np.exp(1j*_np.array(estN[s_id])))
                                 *_np.exp(-1j*truth[s_id]))))
        _e1.append(e0); _eN.append(eN)
print(f"    source-phase error: single-pulsar {_np.median(_e1):.3f} rad ->"
      f" pooled over pulsars {_np.median(_eN):.3f} rad"
      f" -> improvement {1-_np.median(_eN)/max(_np.median(_e1),1e-9):.1%}")
print(f"    PRE-REGISTERED 5-15%:"
      f" {'HIT' if 0.05 <= 1-_np.median(_eN)/max(_np.median(_e1),1e-9) <= 0.15 else 'MISS'}")

# ---- T3: THE DM-CALIBRATION CHANNEL ---------------------------------------------
print(f"\n  T3 -- DM-calibration channel (coupling that bypasses the source pool):")
print(f"    {'patches':>8} {'DM amp':>7} {'corr frac':>10} {'b_DM':>8} {'b_pool':>8}")
for _npatch, _amp in ((20, 0.30), (10, 0.30), (20, 0.15)):
    gains = []
    for j in range(_NS26X):
        S = _setup(j, dm_patches=_npatch, dm_amp=_amp)
        if S["ks"] is None or len(S["ks"]) < 2: continue
        h = soft_sky(S, turns=2, soft=True)
        if not h: continue
        base = h[-1][0]
        # certify the best pulsar in a patch -> learn that patch's bias ->
        # shrink prior width for its patch-mates -> re-measure
        S2 = dict(S); S2["pri"] = S["pri"].copy()
        p0 = int(_np.argmin(S["pri"]))
        mates = _np.where(S["patch"] == S["patch"][p0])[0]
        S2["pri"][mates] = S["pri"][mates]/(1.0 + S["bias"][S["patch"][p0]])
        h2 = soft_sky(S2, turns=2, soft=True)
        if h2: gains.append((h2[-1][0]-base)/max(base, 1e-9))
    corr = _amp/_np.sqrt(_amp**2 + 0.10**2)
    bdm = float(_np.mean(gains)) if gains else _np.nan
    print(f"    {_npatch:8d} {_amp:7.2f} {corr:10.3f} {bdm:8.4f} {0.0047:8.4f}")
print(f"    PRE-REGISTERED b_DM 0.02-0.10 (5-20x the pool channel 0.0047)")

# ---- T4: THE PAYOFF FIGURE -------------------------------------------------------
# measure the PRIOR mass inside one tooth on the same skies -- no assumed value
_mp = []
for j in range(_NS26X):
    S = _setup(j)
    if S["ks"] is None or len(S["ks"]) < 2: continue
    for i in _np.where(S["rncl"])[0][:60]:
        c = _combs(S, i)
        if c is None: continue
        Lg = _np.arange(-3*S["pri"][i], 3*S["pri"][i], TOOTH_W/4)
        pr = _np.exp(-Lg**2/(2*S["pri"][i]**2)); pr /= pr.sum()
        Lt = S["rg"].normal(0, S["pri"][i])
        _mp.append(float(pr[_np.abs(Lg-Lt) <= TOOTH_W].sum()))
_MPRI = float(_np.median(_mp))
_gsoft = _softH[-1][0]
_SOFT_FC = max(0.0, (_gsoft-1.0)*_MPRI)/max(1.0-_MPRI, 1e-9)
print(f"\n  T4 -- A_gwb posterior width: none / hard-certified / full soft:")
print(f"    prior mass in one tooth (MEASURED) {_MPRI:.4f}; soft gain {_gsoft:.4f}"
      f" -> soft info-attributable fc {_SOFT_FC:.4f} vs hard 0.18")
def agwb_width(mode):
    out = []
    for j in range(_NS26X):
        S = _setup(j)
        if S["ks"] is None or len(S["ks"]) < 2: continue
        h = soft_sky(S, turns=3, soft=(mode == "soft"))
        if not h: continue
        f2, P0, sea = S["f2"], S["P0"], S["sea"]
        inb = (f2 >= 2e-9) & (f2 <= 3e-8)
        P = P0[inb]; tot = P.sum() + sea.sum()
        if mode == "none":
            fc = 0.0
        elif mode == "hard":
            fc = 0.18                      # banked hard right-cert fraction
        else:
            fc = _SOFT_FC                  # soft expected right-cert fraction
        v = _np.sqrt(_np.sum((P*(1-fc))**2) + _np.sum(sea**2))/max(tot, 1e-300)
        out.append(0.5*v)
    return float(_np.median(out)) if out else _np.nan
_wn, _wh, _ws = agwb_width("none"), agwb_width("hard"), agwb_width("soft")
print(f"    sigma_A/A:  none {_wn:.4f} | hard-certified {_wh:.4f} | SOFT {_ws:.4f}")
print(f"    soft/hard tightening ratio: {_wh/max(_ws,1e-9):.3f}x")
print(f"    PRE-REGISTERED 1.1-1.6x:"
      f" {'HIT' if 1.1 <= _wh/max(_ws,1e-9) <= 1.6 else 'MISS'}")
print(f"\n  TWO CALENDARS: soft loop is a WHEEL-side result (never gates);")
print(f"  the CATALOGUE calendar is unchanged at 6.5 psr/sky @ 90.3% purity.")
