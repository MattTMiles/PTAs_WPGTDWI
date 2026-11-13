# JOB 53 Part 5: the four scoreboard cases of one_pulsar_six_worlds.ipynb with the joint fit.
# Geometry (kappa of the 8 sources, prior 2.0 pc, truth +1.42 pc) and per-source sigma_phi are the notebook's
# (printed from its re-executed kernel).  Data here: complex sum  z_k = e^{i kappa_k L}(1 + C_k),
# E|C_k|^2 = sigma_phi_k^2 = (1-s_k)/s_k, complex Gaussian.  Decoder: joint fit = matched filter with weights
# 2/P_k (one unit-amplitude source per bin: see job53_heatmap.py header), fine grid, prior N(0, 2 pc).
import json, re, numpy as np
TW = 0.25; SD = 2.0; LT = 1.42
nbtxt = open("one_pulsar_six_worlds.ipynb").read()
print("check of the stored notebook text: '0.965' present:", "0.965" in nbtxt, "| '0.515':", "0.515" in nbtxt, "| '0.985':", "0.985" in nbtxt, "| '0.540':", "0.540" in nbtxt, "| 'K (coarsest comb)  : 8.2':", "K (coarsest comb)  : 8.2" in nbtxt)
rep = open("ACTIVE_WHEEL_REPORT.md").read().splitlines()
hits = [i+1 for i, l in enumerate(rep) if "0.985" in l and "0.540" in l or ("0.985" in l and "red" in l.lower())]
print("ACTIVE_WHEEL_REPORT.md lines mentioning 0.985:", [i+1 for i, l in enumerate(rep) if "0.985" in l][:8], "; 0.540:", [i+1 for i, l in enumerate(rep) if "0.540" in l][:8])
kap = np.array([28.740931, 28.205128, 15.761829, 25.835907, 26.778100, 12.124330, 4.285693, 16.260948])
print(f"K = 6 sigma_d kappa/2pi: coarsest tooth {6*SD*kap.min()/(2*np.pi):.2f}, finest {6*SD*kap.max()/(2*np.pi):.1f}")
W = {"1 clean": np.full(8, 0.30), "2 + red noise": np.full(8, 0.60), "3 + GWB (s = 0.5)": np.full(8, 1.00),
     "4 + realistic loudness": np.array([0.0645576, 17.62738, 52.95773, 74.14730, 111.2266, 111.2266, 111.2266, 111.2266]),
     "5 real world (4 + red noise)": np.array([0.603463, 17.637591, 52.961132, 74.149728, 111.228189, 111.228189, 111.228189, 111.228189])}
NB = {"1 clean": 1.000, "2 + red noise": 0.965, "3 + GWB (s = 0.5)": 0.515, "4 + realistic loudness": 0.000, "5 real world (4 + red noise)": 0.000}
dL = min(0.05, 2*np.pi/kap.max()/6); Lg = np.arange(-3*SD, 3*SD+dL/2, dL); E = np.exp(-1j*kap[:, None]*Lg[None, :])
def run(sig, nreal, truth, seed, gauss_phase=False, cut=False):
    rg = np.random.default_rng(seed); ok = 0; P = sig**2
    use = sig < np.pi if cut else np.ones(8, bool)
    for t in range(nreal):
        Lt = truth if truth is not None else float(np.clip(rg.normal(0, SD), -2.9*SD, 2.9*SD))
        if gauss_phase:      # the notebook's generator and decoder (real Gaussian phase error, wrapped quadratic)
            eps = kap*Lt+rg.normal(0, sig)
            r = (kap[:, None]*Lg[None, :]-eps[:, None]+np.pi) % (2*np.pi)-np.pi
            S = -((r**2)/(2*P[:, None]))[use].sum(0)-Lg**2/(2*SD**2)
        else:
            C = (rg.normal(0, 1, 8)+1j*rg.normal(0, 1, 8))*np.sqrt(P/2); z = np.exp(1j*kap*Lt)*(1+C)
            S = ((2.0/P)[:, None]*np.real(z[:, None]*E)).sum(0)-Lg**2/(2*SD**2)
        ok += abs(Lg[int(np.argmax(S))]-Lt) <= TW
    return ok/nreal
print(f"\n{'world':<30s} notebook | notebook method rebuilt here (truth 1.42) | joint fit, truth 1.42 | joint fit, truth drawn from the prior | notebook method, truth drawn")
for k, sig in W.items():
    a = run(sig, 4000, LT, 1, gauss_phase=True, cut=True); b = run(sig, 4000, LT, 2); c = run(sig, 4000, None, 3); d = run(sig, 4000, None, 4, gauss_phase=True, cut=True)
    print(f"{k:<30s} {NB[k]:.3f}    | {a:.3f}                                     | {b:.3f}                 | {c:.3f}                                 | {d:.3f}")
print("prior only, truth drawn from the prior: %.3f; truth fixed at 1.42 pc: 0.000" % np.mean(np.abs(np.clip(np.random.default_rng(9).normal(0, SD, 200000), -2.9*SD, 2.9*SD)) <= TW))
