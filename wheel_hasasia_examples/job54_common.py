# ================= JOB 54 -- ladder, frequency/amplitude split, talk exports ====================
# Joint54 = job53_common.Joint53 (same bins, same draws, same statistic) with
#   * more variants (table VAR): quiet-sky ladder steps and the two split arms
#       fo  known pulsar terms at the retarded FREQUENCY with the unretarded amplitude
#       ao  known pulsar terms at the Earth frequency with the retarded AMPLITUDE
#     (in both, every other source is retarded as in the reference)
#   * the recovered distance saved for every arm
#   * more per-pair diagnostics for N = 12 (phase error, tooth constants, composition on and off)
#   * build() / decode() exposed for the talk exports (job54_export.py)
# With variants ref / A / off / kn / ot it returns Joint53's flags exactly (job54_test.py).
import numpy as np
import scipy.fft as sfft
TW = 0.25
NCAP = 2**22
TL = {"t1": 1.0, "t10": 10.0, "tinf": 1e6}
#        known state, other state, Earth terms removed, sea factor, other discrete sources factor
VAR = {"ref": ("on", "on", False, 1.0, 1.0), "A": ("on", "on", True, 1.0, 1.0), "off": ("off", "off", False, 1.0, 1.0),
       "kn": ("on", "off", False, 1.0, 1.0), "ot": ("off", "on", False, 1.0, 1.0),
       "fo": ("fo", "on", False, 1.0, 1.0), "ao": ("ao", "on", False, 1.0, 1.0),
       "q1": ("off", "off", False, 0.0, 0.0),      # known binaries + timing noise only
       "q2": ("off", "off", False, 1.0, 0.0),      # + sea
       "q2s": ("off", "off", False, 0.0, 1.0)}     # + other binaries (no sea): the swapped order
LADDER = [(v, n, "t1") for n in (7, 12) for v in ("q1", "q2", "q2s", "off", "ref")] + [("ref", n, t) for n in (7, 12) for t in ("t10", "tinf")]
SPLIT = [(v, n, "t1") for n in (7, 12, "all") for v in ("fo", "ao", "ot", "off", "ref")]
def arms_main(): return list(dict.fromkeys(LADDER+SPLIT))
def arms_span(): return [("ref", n, "t1") for n in (7, 12, "all")]


def build(loc, rank, need):
    """every bin of one pulsar: all sources, both terms, sea and timing draws (rng [51, seed0, j, i], as Jobs 51-53)"""
    j, i, seed0 = loc["j"], loc["i"], loc["seed0"]
    f2, P0, cand, fc, Pc = loc["f2"], loc["P0"], loc["cand"], loc["fc"], loc["Pc"]
    sc, Ts, rt = loc["sc"], loc["Ts"], loc["rt"]; rho_off = loc["rho_off"]; bsea = loc["bsea"]
    nC = len(cand); ns = len(f2); nb = int(3.3e-8*Ts)+1; fb = np.arange(1, nb+1)/Ts
    Nb = (Pc*sc*2.0/np.maximum(rho_off, 1e-300)**2/fc**3)[0]*fb**3
    seab = 2.0*np.atleast_1d(bsea(fb))*rt
    reg = np.zeros(ns, bool); reg[cand] = True; non = ~reg
    rg = np.random.default_rng([51, seed0, j, i])
    thE = rg.uniform(0, 2*np.pi, int(non.sum())); thP = rg.uniform(0, 2*np.pi, int(non.sum()))
    thC = rg.uniform(0, 2*np.pi, int(reg.sum())); g = rg.normal(0, 1, (4, nb))
    rk_cand, rk = rank; assert np.array_equal(rk_cand, cand)
    order = np.argsort(rk)[::-1]
    aE = np.sqrt(P0[non]*sc*rt); sE = np.sinc((f2[non][:, None]-fb[None, :])*Ts)
    B = dict(nC=nC, nb=nb, fb=fb, Nb=Nb, seab=seab, order=order, non=non, thE=thE, thP=thP, aE=aE,
             PEn=(aE**2) @ sE**2, dEn=(aE*np.exp(1j*thE)) @ sE, Ts=Ts, sc=sc, rt=rt)
    del sE
    cE = (np.sqrt(Pc*sc)[:, None]*np.sinc((fc[:, None]-fb[None, :])*Ts))[order]
    eE = np.exp(1j*thC[order])[:, None]*cE
    B.update(cE=cE, eE=eE, cumE=np.vstack([np.zeros(nb), np.cumsum(eE, 0)]), cumPE=np.vstack([np.zeros(nb), np.cumsum(cE**2, 0)]),
             d_sea=(g[0]+1j*g[1])*np.sqrt(seab/2), d_tim1=(g[2]+1j*g[3])*np.sqrt(Nb/2), fE=fc[order], thC=thC[order])
    ST = {}
    on = (loc["fP_all"], loc["pf_all"], loc["fPc"], loc["pfc"], loc["kap"]["on"]); off = (f2, np.ones(ns), fc, np.ones(nC), loc["kap"]["off"])
    for st in sorted(need):
        fP_all, pf_all, fPc, pfc, kap = {"on": on, "off": off, "fo": (on[0], off[1], on[2], off[3], on[4]), "ao": (off[0], on[1], off[2], on[3], off[4])}[st]
        D = dict(cP=(np.sqrt(Pc*pfc*sc)[:, None]*np.sinc((fPc[:, None]-fb[None, :])*Ts))[order], kap=kap[order], fP=fPc[order], pf=pfc[order],
                 far=((fc-fPc)*Ts > 1.0)[order])
        D["cumPP"] = np.vstack([np.zeros(nb), np.cumsum(D["cP"]**2, 0)])
        if st in ("on", "off"):
            aP = np.sqrt(P0[non]*pf_all[non]*sc*rt); sP = np.sinc((fP_all[non][:, None]-fb[None, :])*Ts)
            far = (f2[non]-fP_all[non])*Ts > 1.0
            D.update(PPn=(aP**2) @ sP**2, dPn=(aP*np.exp(1j*thP)) @ sP, PPn_far=((aP**2)*far) @ sP**2, aP=aP, fPn=fP_all[non],
                     cumPPfar=np.vstack([np.zeros(nb), np.cumsum(D["cP"]**2*D["far"][:, None], 0)]))
            del sP
        ST[st] = D
    B["ST"] = ST
    return B


def grid(kap, sd):
    kmax = float(kap.max()); dL = min(0.05, (2*np.pi/max(kmax, 1e-30))/6.0)
    N = int(2**np.ceil(np.log2(max(12.0*sd/dL, 2**15)))); capped = False
    if N > NCAP:
        N = NCAP
        if N*dL < 6.0*sd*1.001: dL = 6.0*sd*1.001/N; capped = True
    dk = 2*np.pi/(N*dL); M = min(int(np.floor(6*sd/dL))+1, N); L0 = -3.0*sd
    return dict(dL=dL, N=N, dk=dk, M=M, L0=L0, Lg=L0+dL*np.arange(M), capped=capped)


def phasors(B, Gd, Lt):
    """pulsar-term data phasors of every candidate on this grid, per state.  B["phi"] (optional): a random phase per
    candidate, carried by its pulsar term in the data and, for a known source, by its template (Job 54: without it the
    unknown candidates' pulsar terms all have phase kappa*Lt and line up with the truth when Lt is small)."""
    G = {}; phi = B.get("phi")
    for st, D in B["ST"].items():
        q = np.minimum(np.rint(D["kap"]/Gd["dk"]).astype(np.int64), Gd["N"]//2-1); kq = q*Gd["dk"]
        eP = np.exp(1j*(kq*Lt+(0.0 if phi is None else phi)))[:, None]*D["cP"]
        G[st] = dict(q=q, kq=kq, eP=eP, cumeP=np.vstack([np.zeros(B["nb"]), np.cumsum(eP, 0)]))
    return G


def setup(B, G, var, n, tl):
    """data and unmodelled variance per bin for one arm; n = number known (rank order)"""
    ks_, os_, a_mode, sea_f, dis = VAR[var]; s = TL[tl]; ST = B["ST"]; Do = ST[os_]; Go = G[os_]
    PEu = (B["cumPE"][-1]-B["cumPE"][n]+B["PEn"])*dis; PPu = (Do["cumPP"][-1]-Do["cumPP"][n]+Do["PPn"])*dis
    d = G[ks_]["cumeP"][n]+dis*((Go["cumeP"][-1]-Go["cumeP"][n])+Do["dPn"])+B["d_tim1"]/s
    if a_mode:
        Pun = PPu+B["seab"]/2+B["Nb"]/s**2; d = d+B["d_sea"]/np.sqrt(2)
    else:
        Pun = PEu+PPu+sea_f*B["seab"]+B["Nb"]/s**2; d = d+dis*(B["cumE"][-1]-B["cumE"][n]+B["dEn"])+sea_f*B["d_sea"]
    return d, Pun, PEu, PPu


def stat_fft(B, Gd, G, var, n, d, Pun):
    kst = VAR[var][0]; ck = B["ST"][kst]["cP"][:n]; qo = G[kst]["q"][:n]; N = Gd["N"]; dk = Gd["dk"]; L0 = Gd["L0"]
    if B.get("phi") is not None: ck = ck*np.exp(1j*B["phi"][:n])[:, None]
    W = 2.0*(ck @ (np.conj(d)/Pun))*np.exp(1j*qo*dk*L0)
    arr = np.zeros(N, complex)
    arr += np.bincount(qo, weights=W.real, minlength=N)+1j*np.bincount(qo, weights=W.imag, minlength=N)
    cw = ck/np.sqrt(Pun)[None, :]; Gm = cw @ np.conj(cw).T
    dq = (qo[:, None]-qo[None, :]); ph = np.exp(1j*dq*(dk*L0)); dq = dq % N
    Gp = (Gm*ph).ravel(); dq = dq.ravel()
    arr -= np.bincount(dq, weights=Gp.real, minlength=N)+1j*np.bincount(dq, weights=Gp.imag, minlength=N)
    return sfft.ifft(arr, norm="forward", workers=1).real[:Gd["M"]]


class Joint54:
    def __init__(self, rank, arms, diag=False, sig_arms=(), sigs=(), randph=False):
        self.rank = rank; self.arms = list(arms); self.diag = diag; self.sig_arms = list(sig_arms); self.sigs = tuple(sigs); self.randph = randph

    def sky(self, loc):
        return {}

    def psr(self, loc):
        i = loc["i"]; Lt = loc["Lt"]; sd = float(loc["pri"][i])
        need = {VAR[v][k] for (v, _, _) in self.arms+self.sig_arms for k in (0, 1)} | ({"on", "off"} if self.diag else set())
        B = build(loc, self.rank[loc["j"]], need); nC = B["nC"]; ST = B["ST"]; out = {"nC": int(nC)}
        if self.randph: B["phi"] = np.random.default_rng([54, 7, loc["seed0"], loc["j"], i]).uniform(0, 2*np.pi, nC)
        for sg in self.sigs:           # prior-width slice: every pulsar given the same prior width
            sds = float(sg); Lts = Lt/sd*sds; tag = f"|s{sg:g}"
            for kst in sorted({VAR[a[0]][0] for a in self.sig_arms}, reverse=True):
                Gd = grid(ST[kst]["kap"], sds); G = phasors(B, Gd, Lts); pr = -Gd["Lg"]**2/(2*sds**2); dL = Gd["dL"]; M = Gd["M"]; Lg = Gd["Lg"]
                for (var, n_, tl) in self.sig_arms:
                    if VAR[var][0] != kst: continue
                    n = nC if n_ == "all" else min(n_, nC); name = f"{var}|{n_}|{tl}{tag}"
                    d, Pun, _, _ = setup(B, G, var, n, tl); S = stat_fft(B, Gd, G, var, n, d, Pun)+pr
                    i0 = int(np.argmax(S)); Lh = Lg[i0]; lo = max(i0-int(np.ceil(TW/dL)), 0); hi = min(i0+int(np.ceil(TW/dL))+1, M)
                    sec = max(S[:lo].max() if lo > 0 else -np.inf, S[hi:].max() if hi < M else -np.inf)
                    out["w|"+name] = bool(abs(Lh-Lts) > TW); out["m|"+name] = float(S[i0]-sec)
                if kst == "on":
                    for n_ in sorted({a[1] for a in self.sig_arms}, key=str):
                        n = nC if n_ == "all" else min(n_, nC); out[f"Kco|{n_}{tag}"] = float(6*sds*G["on"]["kq"][:n].min()/(2*np.pi))
        keep = {}
        for kst in sorted({VAR[a[0]][0] for a in self.arms} | ({"on", "off"} if self.diag else set()), reverse=True):
            Gd = grid(ST[kst]["kap"], sd); G = phasors(B, Gd, Lt); keep[kst] = (Gd, G)
            pr = -Gd["Lg"]**2/(2*sd**2); dL = Gd["dL"]; M = Gd["M"]; Lg = Gd["Lg"]
            if kst == "on": out["j_dL"] = float(dL); out["j_cap"] = bool(Gd["capped"])
            for (var, n_, tl) in self.arms:
                if VAR[var][0] != kst: continue
                n = nC if n_ == "all" else min(n_, nC); name = f"{var}|{n_}|{tl}"
                d, Pun, _, _ = setup(B, G, var, n, tl)
                S = stat_fft(B, Gd, G, var, n, d, Pun)+pr
                i0 = int(np.argmax(S)); Lh = Lg[i0]
                lo = max(i0-int(np.ceil(TW/dL)), 0); hi = min(i0+int(np.ceil(TW/dL))+1, M)
                sec = max(S[:lo].max() if lo > 0 else -np.inf, S[hi:].max() if hi < M else -np.inf)
                out["w|"+name] = bool(abs(Lh-Lt) > TW); out["m|"+name] = float(S[i0]-sec); out["Lh|"+name] = float(Lh)
            if kst == "on":
                for n_ in sorted({a[1] for a in self.arms}, key=str):
                    n = nC if n_ == "all" else min(n_, nC)
                    out[f"Kfin|{n_}"] = float(6*sd*G["on"]["kq"][:n].max()/(2*np.pi)); out[f"Kco|{n_}"] = float(6*sd*G["on"]["kq"][:n].min()/(2*np.pi))
        if self.diag:
            n = min(12, nC); nb = B["nb"]; f32 = np.float32
            for st, var in (("on", "ref"), ("off", "off")):
                Gd, G = keep[st]; D = ST[st]; d, Pun, PEu, PPu = setup(B, G, var, n, "t1")
                cp2 = D["cP"][:n]**2; own = cp2.sum(1); w_ = cp2/own[:, None]; unm = w_ @ Pun
                eP = G[st]["eP"][:n]; r = d[None, :]-(eP.sum(0)[None, :]-eP)
                out[f"d_dphi_{st}"] = np.angle((r*np.conj(eP)).sum(1)).astype(f32)
                out[f"d_inf_{st}"] = (cp2/Pun[None, :]).sum(1).astype(f32)      # sum over bins of own power / unmodelled power (inverse-variance combination)
                out[f"d_own_{st}"] = own.astype(f32); out[f"d_unm_{st}"] = unm.astype(f32); out[f"d_bun_{st}"] = (own/(own+unm)).astype(f32)
                PPfar = D["cumPPfar"][-1]-D["cumPPfar"][n]+D["PPn_far"]
                out[f"d_cmp_{st}"] = np.stack([w_ @ PEu, w_ @ PPfar, w_ @ (PPu-PPfar), w_ @ B["seab"], w_ @ B["Nb"]], 1).astype(f32)
                out[f"d_kap_{st}"] = D["kap"][:n].astype(np.float64); out[f"d_fP_{st}"] = D["fP"][:n].astype(np.float64)
                out[f"b_PP_{st}"] = (D["cumPP"][-1]+D["PPn"]).astype(f32)
                out[f"b_cnt_{st}"] = np.bincount(np.clip(np.rint(D["fP"]*B["Ts"]).astype(int)-1, 0, nb-1), minlength=nb).astype(np.int16)
            out["d_mov"] = ((B["fE"][:n]-ST["on"]["fP"][:n])*B["Ts"]).astype(f32); out["d_pf"] = ST["on"]["pf"][:n].astype(f32)
            out["b_PE"] = (B["cumPE"][-1]+B["PEn"]).astype(f32); out["b_sea"] = B["seab"].astype(f32)
        return out
