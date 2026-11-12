# ================= JOB 51 -- hooks on the Job 50 collector ============================
# Part 2  JointHook : every frequency bin of a pulsar built ONCE (all sources, both terms,
#                     sea, one timing-noise draw); oracle joint fit over all candidates.
# Parts 3/4 StatHook: array SNR of every candidate's Earth term against (a) timing noise,
#                     (b) timing + all other GW power, (c) timing + the part that combining
#                     pulsars cannot remove; pulsar-term share counts.
# Geometry, source lists, distances and sky angles are those of job50_common.collect50
# (same streams), so every arm is on the same skies, pulsars and seeds as Job 50.
import numpy as np
import scipy.fft as sfft
from scipy.interpolate import CubicSpline
TW = 0.25
NCAP = 2**22


class JointHook:
    def __init__(self, noiseless=False):
        self.noiseless = noiseless

    def sky(self, loc):
        return {}

    def psr(self, loc):
        j, i, seed0 = loc["j"], loc["i"], loc["seed0"]
        f2, P0, cand, fc, Pc = loc["f2"], loc["P0"], loc["cand"], loc["fc"], loc["Pc"]
        sc, Ts, rt = loc["sc"], loc["Ts"], loc["rt"]
        fP_all, pf_all, fPc, pfc = loc["fP_all"], loc["pf_all"], loc["fPc"], loc["pfc"]
        kap = loc["kap"]["on"]; Lt = loc["Lt"]; sd = float(loc["pri"][i]); rho_off = loc["rho_off"]
        bsea = loc["bsea"]
        nb = int(3.3e-8*Ts)+1
        fb = np.arange(1, nb+1)/Ts                                 # Fourier bins k/T
        # timing-noise power per bin, same units as P0*sc (calibrated on the candidates:
        # Pm = 2/rho^2 in units of the source's own bin power, job50_common l.211)
        Ncal = Pc*sc*2.0/np.maximum(rho_off, 1e-300)**2/fc**3
        assert np.allclose(Ncal, Ncal[0], rtol=1e-6)
        Nb = Ncal[0]*fb**3
        seab = 2.0*np.atleast_1d(bsea(fb))*rt                      # Earth + pulsar sea
        nonc = np.ones(len(f2), bool); nonc[cand] = False
        aE = np.sqrt(P0[nonc]*sc*rt); aP = np.sqrt(P0[nonc]*pf_all[nonc]*sc*rt)
        sE = np.sinc((f2[nonc][:, None]-fb[None, :])*Ts)           # (n_non, nb)
        sP = np.sinc((fP_all[nonc][:, None]-fb[None, :])*Ts)
        Pnon = (aE**2) @ sE**2+(aP**2) @ sP**2
        Pun = Pnon+seab+Nb                                         # unmodelled variance per bin
        cE = np.sqrt(Pc*sc)[:, None]*np.sinc((fc[:, None]-fb[None, :])*Ts)      # candidates, Earth
        cP = np.sqrt(Pc*pfc*sc)[:, None]*np.sinc((fPc[:, None]-fb[None, :])*Ts)  # candidates, pulsar
        PE = (cE**2).sum(0); PP = (cP**2).sum(0)
        # ---- one realisation of the bins -------------------------------------------------
        rg = np.random.default_rng([51, seed0, j, i])
        thE = rg.uniform(0, 2*np.pi, int(nonc.sum())); thP = rg.uniform(0, 2*np.pi, int(nonc.sum()))
        thC = rg.uniform(0, 2*np.pi, len(cand)); g = rg.normal(0, 1, (4, nb))
        d_non = (aE*np.exp(1j*thE)) @ sE+(aP*np.exp(1j*thP)) @ sP
        d_sea = (g[0]+1j*g[1])*np.sqrt(seab/2); d_tim = (g[2]+1j*g[3])*np.sqrt(Nb/2)
        d_candE = np.exp(1j*thC) @ cE
        # ---- grid ------------------------------------------------------------------------
        kmax = float(kap.max()); lam_min = 2*np.pi/max(kmax, 1e-30)
        dL = min(0.05, lam_min/6.0)
        N = int(2**np.ceil(np.log2(max(12.0*sd/dL, 2**15))))
        capped = False
        if N > NCAP:
            N = NCAP
            if N*dL < 6.0*sd*1.001:
                dL = 6.0*sd*1.001/N; capped = True
        dk = 2*np.pi/(N*dL)
        q = np.rint(kap/dk).astype(np.int64); okq = q < N//2
        kq = q*dk                                                  # the sky's kappa, on the grid
        M = min(int(np.floor(6*sd/dL))+1, N); L0 = -3.0*sd
        Lg = L0+dL*np.arange(M)
        sig_b = (cP*np.exp(1j*kq*Lt)[:, None]).sum(0)              # candidates' pulsar terms, true L
        unm = 0.0 if self.noiseless else (d_non+d_sea+d_tim)
        d_joint = sig_b+unm                                        # candidate Earth terms subtracted
        d_stack = d_joint+(0.0 if self.noiseless else d_candE)     # ... not subtracted
        Ptot = Pun+PE+PP
        ph0 = np.exp(1j*q*dk*L0); qq = q[okq]
        Sj = np.zeros(M); Ss = np.zeros(M)
        for b in range(nb):
            w = (cP[:, b]*ph0)[okq]
            arr = (np.bincount(qq, weights=w.real, minlength=N)
                   + 1j*np.bincount(qq, weights=w.imag, minlength=N)).astype(np.complex64)
            t = sfft.ifft(arr, norm="forward", workers=1)[:M]      # t_b(L) = sum_j c_jb e^{i kq L}
            x = t.real*d_joint[b].real+t.imag*d_joint[b].imag      # Re(d conj t)
            Sj += (2.0*x-(t.real**2+t.imag**2))/Pun[b]
            Ss += 2.0*(t.real*d_stack[b].real+t.imag*d_stack[b].imag)/Ptot[b]
        pr = -Lg**2/(2*sd**2)
        Lj = Lg[int(np.argmax(Sj+pr))]; Lsh = Lg[int(np.argmax(Ss+pr))]
        gw = PE+PP+Pnon+seab; band = (fb >= 2e-9) & (fb <= 3e-8)
        return {"wj": bool(abs(Lj-Lt) > TW), "wsh": bool(abs(Lsh-Lt) > TW),
                "j_dL": float(dL), "j_cap": bool(capped), "j_N": int(N), "j_nq": int(okq.sum()),
                "j_share": float(np.median(((PE+PP)/gw)[band])),
                "j_shareP": float(np.median((PP/gw)[band])),
                "j_share_sum": float((PE+PP)[band].sum()/gw[band].sum()),   # power-weighted over bins
                "j_share_lo": float(np.median(((PE+PP)/gw)[band & (fb <= 1e-8)])),
                "j_unm_gw": float(np.median(((Pnon+seab)/gw)[band])),
                "j_tim_over_un": float(np.median((Nb/Pun)[band])),
                "j_info": float(np.sum(2.0*PP/Pun))}


class StatHook:
    """Per sky: array SNR of each candidate's Earth term three ways; pulsar-term share counts."""
    THR = (0.5, 0.2, 0.1)

    def __init__(self, store=False):
        self.store = store; self._new()

    def _new(self):
        self.oth_bg = []; self.oth_opt = []; self.cp = []; self.sh = []; self.sht = []; self.K = []
        self.pri = []

    def psr(self, loc):
        fc, Pc, sc, rt = loc["fc"], loc["Pc"], loc["sc"], loc["rt"]
        gP = CubicSpline(loc["fg8"], loc["gPg"]*sc)(fc)            # this pulsar's pulsar-term field
        bs = np.atleast_1d(loc["bsea"](fc))
        self.oth_bg.append(np.maximum(loc["gEs"](fc)-Pc*sc+gP+2*bs, 0)*rt)
        self.oth_opt.append(np.maximum(gP+bs, 0)*rt)
        self.cp.append(loc["dkap"]*loc["Lpc"]*loc["arrc"])         # (dkappa * L) * arr
        X = loc["X_on"]
        self.sh.append(1.0/(1.0+X)); self.sht.append(1.0/(1.0+X+2.0/np.maximum(loc["rho_on"], 1e-300)**2))
        self.K.append(6*loc["pri"][loc["i"]]*loc["kap"]["on"]/(2*np.pi))
        self.pri.append(float(loc["pri"][loc["i"]]))
        return {}

    def sky(self, loc):
        Pc, sc, rho_off = loc["Pc"], loc["sc"], loc["rho_off"]
        ob = np.array(self.oth_bg); oo = np.array(self.oth_opt)
        Psig = Pc*sc; Nt = Psig*2.0/np.maximum(rho_off, 1e-300)**2
        arr = lambda oth: np.sqrt(0.2*120*np.mean(2.0*Psig[None, :]/(Nt[None, :]+oth), axis=0))
        a_noise = np.sqrt(0.2*120*2.0*Psig/Nt)
        assert np.allclose(a_noise, loc["arr_all"][loc["cand"]], rtol=1e-6)     # same normalisation as `arr`
        sh = np.array(self.sh); sht = np.array(self.sht); K = np.array(self.K)
        cnt = {}
        for nm, S in (("gw", sh), ("tot", sht)):
            for kn, km in (("any", np.ones_like(K, bool)), ("K1", K < 1), ("K10", K < 10)):
                cnt[f"n_{nm}_{kn}"] = np.array([((S > t) & km).sum(0) for t in self.THR]).astype(np.int16)
        out = dict(arr_noise=a_noise, arr_bg=arr(ob), arr_opt=arr(oo), nclean=len(self.pri), **cnt)
        if self.store:
            out.update(oth_bg=ob.astype(np.float32), oth_opt=oo.astype(np.float32),
                       cpair=np.array(self.cp).astype(np.float32), Psig=Psig, Nt=Nt)
        self._new()
        return out
