# ================= JOB 53 -- natural skies, the pile-up split, prior-width slice ===========
# Joint53 is Job 52's joint fit (job52_common.Joint52, linear-phase arms) rewritten so that the
# retardation state of the KNOWN set and of every OTHER source can be chosen separately, and so
# that one pulsar can be decoded under several prior widths.  Same random draws as Jobs 51/52
# (rng [51, seed0, j, i]); with variant "ref"/"A" it returns Job 52's flags exactly (job53_test.py).
#   variant   known set's pulsar terms   other sources' pulsar terms   unknown Earth terms
#   ref       retarded                   retarded                      in the data
#   A         retarded                   retarded                      removed (and Earth half of the sea)
#   off       Earth frequency            Earth frequency               in the data
#   kn        retarded                   Earth frequency               in the data      ("known only")
#   ot        Earth frequency            retarded                      in the data      ("others only")
import numpy as np
import scipy.fft as sfft
TW = 0.25
NCAP = 2**22
NLIST = (1, 2, 3, 5, 7, 12, 20, 30, 60, 100, 300, "all")
N4 = (7, 12, 30, "all")
N3 = (7, 12, "all")
TL = {"t1": 1.0, "t10": 10.0, "tinf": 1e6}
STATE = {"ref": ("on", "on", False), "A": ("on", "on", True), "off": ("off", "off", False),
         "kn": ("on", "off", False), "ot": ("off", "on", False), "Aoff": ("off", "off", True)}
SIGS = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0)


def arms_main():           # Part 3a + Part 1
    a = [(v, n, "t1") for v in ("ref", "A") for n in NLIST]
    a += [("ref", n, t) for n in N4 for t in ("t10", "tinf")]
    a += [(v, n, "t1") for v in ("off", "kn", "ot") for n in N3]
    return a


def arms_sig():            # Part 3b
    return [(v, n, "t1") for v in ("ref", "A") for n in N4]


def arms_span():           # Part 4
    return [("ref", n, "t1") for n in N3]


class Joint53:
    def __init__(self, rank, arms, sig_arms=(), sigs=(), diag=False):
        self.rank = rank; self.arms = list(arms); self.sig_arms = list(sig_arms); self.sigs = tuple(sigs)
        self.diag = diag

    def sky(self, loc):
        return {}

    def psr(self, loc):
        j, i, seed0 = loc["j"], loc["i"], loc["seed0"]
        f2, P0, cand, fc, Pc = loc["f2"], loc["P0"], loc["cand"], loc["fc"], loc["Pc"]
        sc, Ts, rt = loc["sc"], loc["Ts"], loc["rt"]
        Lt0 = loc["Lt"]; sd0 = float(loc["pri"][i]); rho_off = loc["rho_off"]; bsea = loc["bsea"]
        nC = len(cand); ns = len(f2)
        nb = int(3.3e-8*Ts)+1
        fb = np.arange(1, nb+1)/Ts
        Ncal = Pc*sc*2.0/np.maximum(rho_off, 1e-300)**2/fc**3
        Nb = Ncal[0]*fb**3
        seab = 2.0*np.atleast_1d(bsea(fb))*rt
        need = {STATE[v][k] for (v, _, _) in self.arms+self.sig_arms for k in (0, 1)} | ({"on", "off"} if self.diag else set())
        # ---- random draws, exactly as Jobs 51/52 -------------------------------------------------
        reg = np.zeros(ns, bool); reg[cand] = True
        rg = np.random.default_rng([51, seed0, j, i])
        thE = rg.uniform(0, 2*np.pi, int((~reg).sum())); thP = rg.uniform(0, 2*np.pi, int((~reg).sum()))
        thC = rg.uniform(0, 2*np.pi, int(reg.sum())); g = rg.normal(0, 1, (4, nb))
        non = ~reg
        rk_cand, rk = self.rank[j]; assert np.array_equal(rk_cand, cand)
        order = np.argsort(rk)[::-1]                              # loudest first (arr_bg_opt)
        # ---- Earth terms -----------------------------------------------------------------------
        aE = np.sqrt(P0[non]*sc*rt)
        sE = np.sinc((f2[non][:, None]-fb[None, :])*Ts)
        PEn = (aE**2) @ sE**2; dEn = (aE*np.exp(1j*thE)) @ sE
        del sE
        cE = (np.sqrt(Pc*sc)[:, None]*np.sinc((fc[:, None]-fb[None, :])*Ts))[order]
        eE = np.exp(1j*thC[order])[:, None]*cE
        cumE = np.vstack([np.zeros(nb), np.cumsum(eE, 0)]); cumPE = np.vstack([np.zeros(nb), np.cumsum(cE**2, 0)])
        d_sea = (g[0]+1j*g[1])*np.sqrt(seab/2); d_tim1 = (g[2]+1j*g[3])*np.sqrt(Nb/2)
        # ---- pulsar terms, per retardation state ------------------------------------------------
        ST = {}
        for st in sorted(need):
            if st == "on":
                fP_all, pf_all, fPc, pfc, kap = loc["fP_all"], loc["pf_all"], loc["fPc"], loc["pfc"], loc["kap"]["on"]
            else:
                fP_all, pf_all, fPc, pfc, kap = f2, np.ones(ns), fc, np.ones(nC), loc["kap"]["off"]
            aP = np.sqrt(P0[non]*pf_all[non]*sc*rt)
            sP = np.sinc((fP_all[non][:, None]-fb[None, :])*Ts)
            far = (f2[non]-fP_all[non])*Ts > 1.0
            D = dict(PPn=(aP**2) @ sP**2, dPn=(aP*np.exp(1j*thP)) @ sP, PPn_far=((aP**2)*far) @ sP**2,
                     cP=(np.sqrt(Pc*pfc*sc)[:, None]*np.sinc((fPc[:, None]-fb[None, :])*Ts))[order],
                     kap=kap[order], fP=fPc[order], far=((fc-fPc)*Ts > 1.0)[order])
            D["cumPP"] = np.vstack([np.zeros(nb), np.cumsum(D["cP"]**2, 0)])
            D["cumPPfar"] = np.vstack([np.zeros(nb), np.cumsum(D["cP"]**2*D["far"][:, None], 0)])
            ST[st] = D
            del sP
        out = {"nC": int(nC)}
        band = (fb >= 2e-9) & (fb <= 3e-8)
        if "on" in ST:
            gw = cumPE[-1]+ST["on"]["cumPP"][-1]+PEn+ST["on"]["PPn"]+seab
            for n_ in NLIST:
                n = nC if n_ == "all" else min(n_, nC)
                out[f"sh|{n_}"] = float((cumPE[n]+ST["on"]["cumPP"][n])[band].sum()/gw[band].sum())

        # ---- one prior width: grid, data, arms ---------------------------------------------------
        def case(tag, sd, Lt, arms):
            setup = None
            for kst in sorted({STATE[a[0]][0] for a in arms}, reverse=True):      # grid follows the known set's state
                kmax = float(ST[kst]["kap"].max()); lam_min = 2*np.pi/max(kmax, 1e-30)
                dL = min(0.05, lam_min/6.0)
                N = int(2**np.ceil(np.log2(max(12.0*sd/dL, 2**15)))); capped = False
                if N > NCAP:
                    N = NCAP
                    if N*dL < 6.0*sd*1.001: dL = 6.0*sd*1.001/N; capped = True
                dk = 2*np.pi/(N*dL)
                M = min(int(np.floor(6*sd/dL))+1, N); L0 = -3.0*sd; Lg = L0+dL*np.arange(M)
                pr = -Lg**2/(2*sd**2)
                if kst == "on": out["j_dL"+tag] = float(dL); out["j_cap"+tag] = bool(capped)
                G = {}
                for st, D in ST.items():
                    q = np.minimum(np.rint(D["kap"]/dk).astype(np.int64), N//2-1); kq = q*dk
                    eP = np.exp(1j*kq*Lt)[:, None]*D["cP"]
                    G[st] = dict(q=q, kq=kq, cumeP=np.vstack([np.zeros(nb), np.cumsum(eP, 0)]))

                def setup(var, n, tl, G=G):
                    ks_, os_, a_mode = STATE[var]; s = TL[tl]
                    Do = ST[os_]; Go = G[os_]
                    PEu = cumPE[-1]-cumPE[n]+PEn; PPu = Do["cumPP"][-1]-Do["cumPP"][n]+Do["PPn"]
                    d = G[ks_]["cumeP"][n]+(Go["cumeP"][-1]-Go["cumeP"][n])+Do["dPn"]+d_tim1/s
                    if a_mode:
                        Pun = PPu+seab/2+Nb/s**2; d = d+d_sea/np.sqrt(2)
                    else:
                        Pun = PEu+PPu+seab+Nb/s**2; d = d+(cumE[-1]-cumE[n]+dEn)+d_sea
                    return d, Pun, PEu, PPu

                for (var, n_, tl) in arms:
                    if STATE[var][0] != kst: continue
                    n = nC if n_ == "all" else min(n_, nC); name = f"{var}|{n_}|{tl}{tag}"
                    d, Pun, _, _ = setup(var, n, tl)
                    ck = ST[kst]["cP"][:n]; qo = G[kst]["q"][:n]
                    W = 2.0*(ck @ (np.conj(d)/Pun))*np.exp(1j*qo*dk*L0)
                    arr = np.zeros(N, complex)
                    arr += np.bincount(qo, weights=W.real, minlength=N)+1j*np.bincount(qo, weights=W.imag, minlength=N)
                    cw = ck/np.sqrt(Pun)[None, :]; Gm = cw @ cw.T
                    dq = (qo[:, None]-qo[None, :]); ph = np.exp(1j*dq*(dk*L0)); dq = dq % N
                    Gp = (Gm*ph).ravel(); dq = dq.ravel()
                    arr -= np.bincount(dq, weights=Gp.real, minlength=N)+1j*np.bincount(dq, weights=Gp.imag, minlength=N)
                    S = sfft.ifft(arr, norm="forward", workers=1).real[:M]+pr
                    i0 = int(np.argmax(S)); Lh = Lg[i0]
                    lo = max(i0-int(np.ceil(TW/dL)), 0); hi = min(i0+int(np.ceil(TW/dL))+1, M)
                    sec = max(S[:lo].max() if lo > 0 else -np.inf, S[hi:].max() if hi < M else -np.inf)
                    out["w|"+name] = bool(abs(Lh-Lt) > TW); out["m|"+name] = float(S[i0]-sec)
                    out["C|"+name] = float(np.sum(np.log1p((ck**2).sum(0)/Pun)))
                if kst == "on":
                    for n_ in sorted({a[1] for a in arms}, key=str):
                        n = nC if n_ == "all" else min(n_, nC)
                        out[f"Kfin|{n_}{tag}"] = float(6*sd*G["on"]["kq"][:n].max()/(2*np.pi))
                        out[f"Kco|{n_}{tag}"] = float(6*sd*G["on"]["kq"][:n].min()/(2*np.pi))
            return setup

        setup = case("", sd0, Lt0, self.arms)
        # ---- diagnostics (N = 12 pairs) -----------------------------------------------------------
        if self.diag:
            n = min(12, nC); res = {}
            for st, var in (("on", "ref"), ("off", "off")):
                D = ST[st]; d, Pun, PEu, PPu = setup(var, n, "t1")
                cp2 = D["cP"][:n]**2; own = cp2.sum(1); w_ = cp2/own[:, None]
                unm = w_ @ Pun
                res[st] = (own, unm)
                out[f"d_own_{st}"] = own.astype(np.float32); out[f"d_unm_{st}"] = unm.astype(np.float32)
                out[f"d_bun_{st}"] = (own/(own+unm)).astype(np.float32)
                if st == "on":
                    _, Pkn, _, _ = setup("kn", n, "t1")              # own term retarded, every other pulsar term not
                    out["d_unm_kn"] = (w_ @ Pkn).astype(np.float32)
                    PPfar = D["cumPPfar"][-1]-D["cumPPfar"][n]+D["PPn_far"]
                    out["d_cmp"] = np.stack([w_ @ PEu, w_ @ PPfar, w_ @ (PPu-PPfar), w_ @ seab, w_ @ (Nb*np.ones(nb))], 1).astype(np.float32)
                    out["d_fP"] = D["fP"][:n].astype(np.float64); out["d_fE"] = fc[order][:n].astype(np.float64)
                    out["d_mov"] = ((fc[order][:n]-D["fP"][:n])*Ts).astype(np.float32)
                    out["d_bbin12"] = (own/(w_ @ gw)).astype(np.float32)
                    cpa = D["cP"]**2; owa = cpa.sum(1)
                    bb = owa/np.maximum((cpa @ gw)/np.maximum(owa, 1e-300), 1e-300)
                    out["d_bbin_max"] = float(bb.max()); out["d_bbin_max_rank"] = int(np.argmax(bb))
        # ---- prior-width slice ------------------------------------------------------------------------
        for sg in self.sigs:
            case(f"|s{sg:g}", float(sg), Lt0/sd0*float(sg), self.sig_arms)
        return out
