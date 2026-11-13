# ================= JOB 52 -- best case: known-set size N, arms A-D ======================
# Reference = Job 51 joint fit (job51_common.JointHook): bins built once per pulsar, known
# sources' Earth terms subtracted, template = known sources' pulsar-term phasors at trial L,
# unmodelled = everything else + timing noise.  Same draws as Job 51 (rng [51,seed0,j,i]).
#   N     known set = the N candidates ranked highest by arr_bg_opt (Job 51 Part 3)
#   A     Earth terms of every source outside the known set deleted; Earth half of the sea deleted
#   B     candidate pool includes Earth frequencies above 30 nHz; bins to 100 nHz  (ext config)
#   C     each known source's retarded frequency, amplitude, leakage and phase follow trial L
#   D     A + B + C
# Linear-phase arms: S(L) = Re sum_j W_j e^{i k_j L} - sum_jk G_jk e^{i(k_j-k_k)L}  -> one FFT.
import numpy as np
import scipy.fft as sfft
from scipy.interpolate import CubicSpline
import job51_common as H51
TW = 0.25
NCAP = 2**22
NLIST = (1, 2, 3, 5, 7, 12, 20, 30, 60, 100, 300, "all")
N4 = (7, 12, 30, "all")
TL = {"t1": 1.0, "t10": 10.0, "tinf": 1e6}
SEGCAP = 4000


def arms_for(ext, mode="main"):
    if mode == "main":
        if ext: return [(v, n, t) for v in ("B", "D") for n in N4 for t in ("t1", "tinf")]
        a = [("ref", n, "t1") for n in NLIST]
        a += [("ref", n, t) for n in N4 for t in ("t10", "tinf")]
        a += [(v, n, t) for v in ("A", "C") for n in N4 for t in ("t1", "tinf")]
        return a
    if mode == "p5":      # Part 5: one change at a time (N = 12, all), x1
        return [(v, n, "t1") for v in ("ref", "nosea", "nodisc", "earthin") for n in (12, "all")] \
            + [("ref", n, t) for n in (12, "all") for t in ("tinf", "a37", "a106")]
    raise ValueError(mode)
TLX = dict(TL, a37=3.7/6.4, a106=10.6/6.4)     # amplitude x a == timing noise / a (all GW scales together)


class Stat52(H51.StatHook):
    def psr(self, loc):
        if not self.pri: self.fP0 = loc["fPc"].copy()
        return super().psr(loc)

    def sky(self, loc):
        out = super().sky(loc); out.update(fc=loc["fc"].copy(), fP0=self.fP0, cand=loc["cand"].copy())
        return out


class Joint52:
    def __init__(self, rank, ext=False, mode="main", ret=True, part4=False, binmax=None, arms=None):
        self.rank = rank; self.ext = ext; self.mode = mode; self.ret = ret; self.part4 = part4
        self.binmax = binmax                    # diagnostic: ignore bins above this frequency
        self.arms = arms_for(ext, mode) if arms is None else arms

    def sky(self, loc):
        return {}

    # ------------------------------------------------------------------------------------
    def psr(self, loc):
        j, i, seed0 = loc["j"], loc["i"], loc["seed0"]
        f2, P0, cand, fc, Pc = loc["f2"], loc["P0"], loc["cand"], loc["fc"], loc["Pc"]
        sc, Ts, rt = loc["sc"], loc["Ts"], loc["rt"]
        Lt = loc["Lt"]; sd = float(loc["pri"][i]); rho_off = loc["rho_off"]; bsea = loc["bsea"]
        if self.ret:
            fP_all, pf_all, fPc, pfc = loc["fP_all"], loc["pf_all"], loc["fPc"], loc["pfc"]
            kap = loc["kap"]["on"]
        else:
            fP_all, pf_all, fPc, pfc = f2, np.ones(len(f2)), fc, np.ones(len(fc)); kap = loc["kap"]["off"]
        nC = len(cand); ns = len(f2)
        nb0 = int(3.3e-8*Ts)+1
        nb = int(1.0e-7*Ts)+1 if self.ext else nb0
        fb = np.arange(1, nb+1)/Ts
        Ncal = Pc*sc*2.0/np.maximum(rho_off, 1e-300)**2/fc**3
        Nb = Ncal[0]*fb**3
        seab = 2.0*np.atleast_1d(bsea(fb))*rt
        # ---- regular pool (Job 51 list) -> identical random draws in both configs ------------
        if self.ext:
            rp, ar = loc["rp_all"], loc["arr_all"]
            rc = np.where((ar >= 3.0) & (f2 >= 2e-9) & (f2 <= 3e-8) & (rp >= 3.0))[0]
            if len(rc) > 1300: rc = rc[np.argsort(rp[rc])[::-1][:1300]]
            rc = np.union1d(rc, loc["ks"])
        else:
            rc = cand
        reg = np.zeros(ns, bool); reg[rc] = True
        rg = np.random.default_rng([51, seed0, j, i])
        thE = rg.uniform(0, 2*np.pi, int((~reg).sum())); thP = rg.uniform(0, 2*np.pi, int((~reg).sum()))
        thC = rg.uniform(0, 2*np.pi, int(reg.sum())); g = rg.normal(0, 1, (4, nb0))
        if nb > nb0: g = np.hstack([g, rg.normal(0, 1, (4, nb-nb0))])
        phE = np.empty(ns); phE[~reg] = thE; phE[reg] = thC
        phP = np.zeros(ns); phP[~reg] = thP
        pool = np.zeros(ns, bool); pool[cand] = True
        # ---- grid (as Job 51) ------------------------------------------------------------------
        kmax = float(kap.max()); lam_min = 2*np.pi/max(kmax, 1e-30)
        dL = min(0.05, lam_min/6.0)
        N = int(2**np.ceil(np.log2(max(12.0*sd/dL, 2**15)))); capped = False
        if N > NCAP:
            N = NCAP
            if N*dL < 6.0*sd*1.001: dL = 6.0*sd*1.001/N; capped = True
        dk = 2*np.pi/(N*dL)
        q = np.rint(kap/dk).astype(np.int64); q = np.minimum(q, N//2-1); kq = q*dk
        M = min(int(np.floor(6*sd/dL))+1, N); L0 = -3.0*sd; Lg = L0+dL*np.arange(M)
        pr = -Lg**2/(2*sd**2)
        # ---- per-source bin amplitudes -----------------------------------------------------------
        non = ~pool
        aE = np.sqrt(P0[non]*sc*rt); aP = np.sqrt(P0[non]*pf_all[non]*sc*rt)
        sE = np.sinc((f2[non][:, None]-fb[None, :])*Ts); sP = np.sinc((fP_all[non][:, None]-fb[None, :])*Ts)
        PEn = (aE**2) @ sE**2; PPn = (aP**2) @ sP**2
        dEn = (aE*np.exp(1j*phE[non])) @ sE; dPn = (aP*np.exp(1j*phP[non])) @ sP
        del sE, sP
        cE = np.sqrt(Pc*sc)[:, None]*np.sinc((fc[:, None]-fb[None, :])*Ts)
        cP = np.sqrt(Pc*pfc*sc)[:, None]*np.sinc((fPc[:, None]-fb[None, :])*Ts)
        eE = np.exp(1j*phE[cand])[:, None]*cE; eP = np.exp(1j*kq*Lt)[:, None]*cP
        rk_cand, rk = self.rank[j]; assert np.array_equal(rk_cand, cand)
        order = np.argsort(rk)[::-1]                              # loudest first (arr_bg_opt)
        cE, cP, eE, eP, qo, kqo = cE[order], cP[order], eE[order], eP[order], q[order], kq[order]
        cumE = np.vstack([np.zeros(nb), np.cumsum(eE, 0)]); cumPE = np.vstack([np.zeros(nb), np.cumsum(cE**2, 0)])
        cumPP = np.vstack([np.zeros(nb), np.cumsum(cP**2, 0)])
        dP_pool = eP.sum(0)
        d_sea = (g[0]+1j*g[1])*np.sqrt(seab/2); d_tim1 = (g[2]+1j*g[3])*np.sqrt(Nb/2)
        gw = cumPE[-1]+cumPP[-1]+PEn+PPn+seab
        band = (fb >= 2e-9) & (fb <= 3e-8); blo = fb <= 1e-8; bhi = (fb >= 1.5e-8) & (fb <= 4e-8)
        out = {"j_dL": float(dL), "j_cap": bool(capped), "j_M": int(M), "nC": int(nC)}
        # diagnostic: correlation of the unmodelled data between neighbouring bins (every candidate known)
        u_ = (cumE[-1]*0+dEn+dPn+d_sea)/np.sqrt(PEn+PPn+seab)
        hb = fb > 4e-8; lb = (fb > 2e-9) & (fb < 1e-8)
        out["corr_hi"] = float(np.abs(np.mean((u_[1:]*np.conj(u_[:-1]))[hb[1:]]))) if hb.sum() > 3 else np.nan
        out["corr_lo"] = float(np.abs(np.mean((u_[1:]*np.conj(u_[:-1]))[lb[1:]])))

        def setup(var, n, tl):
            """data and unmodelled variance for one arm; n = number known (rank order)"""
            s = TLX[tl]; a_mode = var in ("A", "D")
            PEu = cumPE[-1]-cumPE[n]+PEn; PPu = cumPP[-1]-cumPP[n]+PPn
            dEu = cumE[-1]-cumE[n]+dEn
            sea_f = 1.0; dis = 1.0; ein = 0.0
            if var == "nosea": sea_f = 0.0
            if var == "nodisc": dis = 0.0
            if var == "earthin": ein = 1.0
            d = dP_pool*1.0+d_tim1/s
            if dis == 0.0:
                d = eP[:n].sum(0)+d_tim1/s; PPu = PPu*0; PEu = PEu*0; dEu = dEu*0
            else:
                d = d+dPn
            if a_mode:
                Pun = PPu+seab/2+Nb/s**2; d = d+d_sea/np.sqrt(2)
            else:
                Pun = PEu+PPu+sea_f*seab+Nb/s**2+ein*cumPE[n]; d = d+dEu+sea_f*d_sea+ein*cumE[n]
            if self.binmax is not None: Pun = np.where(fb > self.binmax, Pun*1e30, Pun)
            return d, Pun, PEu

        def finish(name, S, Pk, Pun):
            S = S+pr; i0 = int(np.argmax(S)); Lh = Lg[i0]
            lo = max(i0-int(np.ceil(TW/dL)), 0); hi = min(i0+int(np.ceil(TW/dL))+1, M)
            sec = max(S[:lo].max() if lo > 0 else -np.inf, S[hi:].max() if hi < M else -np.inf)
            out["w|"+name] = bool(abs(Lh-Lt) > TW); out["m|"+name] = float(S[i0]-sec)
            out["C|"+name] = float(np.sum(np.log1p(Pk/Pun)))

        # ---- per-N diagnostics --------------------------------------------------------------------
        ns_used = sorted({nC if n == "all" else min(n, nC) for (_, n, _) in self.arms})
        for n in ns_used:
            tag = "all" if n == nC else str(n)
            pk = cumPE[n]+cumPP[n]
            out["sh|"+tag] = float(pk[band].sum()/gw[band].sum())
            out["shhi|"+tag] = float(np.median((pk/gw)[bhi])) if bhi.any() else np.nan
            out["Kfin|"+tag] = float(6*sd*kqo[:n].max()/(2*np.pi)); out["Kco|"+tag] = float(6*sd*kqo[:n].min()/(2*np.pi))
            PEu = cumPE[-1]-cumPE[n]+PEn; un = PEu+(cumPP[-1]-cumPP[n]+PPn)+seab
            fe = (PEu+seab/2)/un
            out["fe_lo|"+tag] = float(np.median(fe[blo & band])); out["fe_hi|"+tag] = float(np.median(fe[bhi]))
            out["fe_all|"+tag] = float(np.median(fe[band]))
        # ---- linear-phase arms: one FFT each ---------------------------------------------------------
        ph0 = np.exp(1j*qo*dk*L0)
        for (var, n_, tl) in self.arms:
            if var in ("C", "D"): continue
            n = nC if n_ == "all" else min(n_, nC); name = f"{var}|{n_}|{tl}"
            d, Pun, _ = setup(var, n, tl)
            ck = cP[:n]
            W = 2.0*(ck @ (np.conj(d)/Pun))*ph0[:n]
            arr = np.zeros(N, complex)
            arr += np.bincount(qo[:n], weights=W.real, minlength=N)+1j*np.bincount(qo[:n], weights=W.imag, minlength=N)
            cw = ck/np.sqrt(Pun)[None, :]; G = cw @ cw.T
            dq = (qo[:n, None]-qo[None, :n]); ph = np.exp(1j*dq*(dk*L0)); dq = dq % N
            Gp = (G*ph).ravel(); dq = dq.ravel()
            arr -= np.bincount(dq, weights=Gp.real, minlength=N)+1j*np.bincount(dq, weights=Gp.imag, minlength=N)
            S = sfft.ifft(arr, norm="forward", workers=1).real[:M]
            finish(name, S, (ck**2).sum(0), Pun)
        # ---- arms C / D: frequency, amplitude, leakage and phase follow the trial distance ------------
        cd = [a for a in self.arms if a[0] in ("C", "D")]
        if cd:
            C_ = loc["C"]; Lpc = loc["Lpc"]
            mu = loc["mu"][order]; tc = np.maximum(loc["tc_all"][cand][order], 1e-9); f0 = fc[order]; Po = Pc[order]
            gfac = (1-np.cos(mu))/C_
            tau0 = Lpc*gfac
            Labs = np.maximum(Lpc+(Lg-Lt), 0.0)
            fPl = lambda La: f0[:, None]*(1.0+(La[None, :]*gfac[:, None])/tc[:, None])**(-3.0/8.0)
            ends = fPl(np.array([Labs[0], Labs[-1]]))
            dfb = np.abs(ends[:, 1]-ends[:, 0])*Ts
            out["c_mov01_all"] = float(np.mean(dfb > 0.1)); out["c_mov01_12"] = float(np.mean(dfb[:12] > 0.1))
            out["c_mov01_30"] = float(np.mean(dfb[:30] > 0.1)); out["c_dfmax"] = float(dfb.max())
            nseg = int(np.clip(np.ceil(dfb.max()/0.01), 1, SEGCAP)); out["c_nseg"] = nseg
            out["c_segcap"] = bool(dfb.max()/0.01 > SEGCAP)
            levels = sorted({nC if n == "all" else min(n, nC) for (_, n, _) in cd})
            setups = {}
            for (var, n_, tl) in cd:
                n = nC if n_ == "all" else min(n_, nC)
                d, Pun, _ = setup(var, n, tl); setups[(var, n_, tl)] = (n, np.conj(d)/Pun, 1.0/Pun, Pun)
            Sall = {k: np.empty(M) for k in setups}
            base = (8.0/5.0)*2*np.pi*f0*tc; ref0 = (1.0+tau0/tc)**(5.0/8.0); phT = kqo*Lt
            edges = np.linspace(0, M, nseg+1).astype(int)
            CH = max(int(4e6//max(nC, 1)), 256)
            for sgi in range(nseg):
                a0, a1 = edges[sgi], edges[sgi+1]
                if a1 <= a0: continue
                Lm = np.array([Labs[(a0+a1)//2]])
                fPm = fPl(Lm)[:, 0]
                cs = (np.sqrt(Po*sc*(fPm/f0)**(7.0/3.0))[:, None]*np.sinc((fPm[:, None]-fb[None, :])*Ts)).astype(np.complex64)
                for c0 in range(a0, a1, CH):
                    c1 = min(c0+CH, a1); La = Labs[c0:c1]
                    phi = phT[:, None]+base[:, None]*((1.0+(La[None, :]*gfac[:, None])/tc[:, None])**(5.0/8.0)-ref0[:, None])
                    E = np.exp(1j*phi).astype(np.complex64)
                    tcum = {}; t = np.zeros((nb, c1-c0), np.complex64); prev = 0
                    for n in levels:
                        if n > prev: t = t+cs[prev:n].T @ E[prev:n]
                        tcum[n] = t; prev = n
                    for k, (n, v, ip, Pun) in setups.items():
                        tt = tcum[n]
                        Sall[k][c0:c1] = 2.0*np.real(v.astype(np.complex64) @ tt)-(tt.real**2+tt.imag**2).T @ ip
            for k, (n, v, ip, Pun) in setups.items():
                finish(f"{k[0]}|{k[1]}|{k[2]}", Sall[k], (cP[:n]**2).sum(0), Pun)
        # ---- Part 4: phase jitter measured in the bins (N = 12, reference) -----------------------------
        if self.part4:
            n = min(12, nC)
            for tl in ("t1", "tinf"):
                d, Pun, _ = setup("ref", n, tl)
                r = d[None, :]-(eP[:n].sum(0)[None, :]-eP[:n])        # minus the other known pulsar terms
                z = (r*np.conj(eP[:n])).sum(1)                       # amplitude-weighted sum over bins
                own = (cP[:n]**2).sum(1)
                unm = (cP[:n]**2 @ Pun)/own                           # unmodelled power under the pair's leakage profile
                gwp = (cP[:n]**2 @ gw)/own
                out["p4_dphi_"+tl] = np.angle(z).astype(np.float32)
                out["p4_bun_"+tl] = (own/(own+unm)).astype(np.float32)
                out["p4_bbin_"+tl] = (own/gwp).astype(np.float32)
                if getattr(self, "p4_debug", False):                  # Job 52 follow-up: what the unmodelled power is made of
                    w_ = cP[:n]**2/own[:, None]; s_ = TLX[tl]
                    PEu_ = cumPE[-1]-cumPE[n]+PEn
                    out["p4_cmp_"+tl] = np.stack([w_ @ PEu_, w_ @ (Pun-PEu_-seab-Nb/s_**2), w_ @ (seab/2), w_ @ (seab/2), w_ @ (Nb/s_**2),
                                                  w_ @ (cumPE[-1]-cumPE[n]), w_ @ PEn], 1).astype(np.float32)   # unknown Earth, unknown pulsar, sea E, sea P, timing; Earth split cand / non-cand
                    out["p4_own"] = own.astype(np.float32)
                    out["p4_f"] = fc[order][:n].astype(np.float64); out["p4_fP"] = fPc[order][:n].astype(np.float64)
        return out
