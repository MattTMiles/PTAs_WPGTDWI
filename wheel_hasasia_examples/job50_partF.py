# ================= JOB 50 PART F ======================================================
# F1: share-versus-frequency figure with BOTH terms of every source, pulsar terms at
#     retarded frequencies with retarded amplitudes, true beta on the vertical axis,
#     no readability line.  New file; the old figure is untouched.
# F2: real-likelihood sqrt(2).  stagec_p2b.py:57 takes A_CW as the RMS of the pulsar-term
#     residual and :86 divides the per-quadrature amplitude error by it.  Corrected lines
#     (applied here to the SAVED per-pulsar arrays A_CW, sig_res of stagec_p2b_results.npz,
#     so nothing booked is rerun or overwritten):
#         amp       = np.sqrt(2.0) * A_CW          # sinusoid amplitude, not RMS
#         rho       = amp / sig_res                # standard pulsar-term SNR
#         sigma_phi = 1.0 / rho                    # small-angle phase error
import os, sys, json
import numpy as np
from scipy.special import i0e, i1e
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C
NS = C.J.load(verbose=False); os.chdir(C.J.BASE)
OUT = {}
# ---------------- F1 -------------------------------------------------------------------
f, h, T = NS["f_src"], NS["h_src"], NS["T_OBS"]; edges = NS["bin_edges"]; sea = NS["S_FAINT"]; nb = len(edges)-1
rngA = np.random.default_rng(404)
L = rngA.uniform(500.0, 3000.0, 120); mu = np.arccos(rngA.uniform(-0.9, 0.9, (120, len(f))))
tau = L[41]*(1-np.cos(mu[41]))/NS["C_PC_S"]
PE = h**2*f*T
def beta_true(tc):
    fP = C.f_ret(f, tc, tau); PP = PE*(fP/f)**(7/3)
    bE = np.digitize(f, edges)-1; bP = np.digitize(fP, edges)-1
    tot = 2*sea.copy()
    okE = (bE >= 0) & (bE < nb); okP = (bP >= 0) & (bP < nb)
    np.add.at(tot, bE[okE], PE[okE]); np.add.at(tot, bP[okP], PP[okP])
    b = np.full(len(f), np.nan); b[okP] = PP[okP]/tot[bP[okP]]
    return fP, b
ts_nb = NS["G_SI"]*NS["mc_src"]/NS["C_SI"]**3
tc_nb = (5/256)*ts_nb**(-5/3)*(np.pi*f)**(-8/3)
tc_41 = C._mcz_tc(f, h)
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6), sharey=True)
for a, tc, ttl in ((ax[0], tc_nb, "notebook masses (source frame; 232/44 convention)"),
                   (ax[1], tc_41, "Job 41 masses (redshifted; used by p2r)")):
    fP, b = beta_true(tc); ok = np.isfinite(b)
    inb = ok & (fP >= 15e-9) & (fP < 40e-9)
    a.scatter(fP[ok & ~inb]*1e9, b[ok & ~inb], s=3, color="0.7", lw=0)
    a.scatter(fP[inb]*1e9, b[inb], s=4, color="C0", lw=0)
    a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(1e-5, 1.2); a.set_xlim(1, 100)
    a.set_xlabel("pulsar-term frequency f_P [nHz]"); a.set_title(ttl, fontsize=10); a.grid(alpha=.25, lw=.5)
    top = np.sort(b[inb])[::-1][:20]
    OUT[ttl] = dict(n_inband=int(inb.sum()), med=float(np.median(b[inb])), top20_med=float(np.median(top)),
                    n_gt_0p1=int((b[inb] > 0.1).sum()), n_gt_0p01=int((b[inb] > 0.01).sum()))
    print(f"F1 [{ttl}]: 15-40 nHz pulsar terms {inb.sum()}; true share median {np.median(b[inb]):.5f};"
          f" median of top 20 {np.median(top):.3f}; with share > 0.1: {(b[inb] > 0.1).sum()}, > 0.01: {(b[inb] > 0.01).sum()}")
ax[0].set_ylabel("true share: pulsar-term power / all GW power in its 1/T bin\n(Earth + pulsar terms of every source + sea, this pulsar)")
plt.tight_layout()
plt.savefig("/home/mattm/projects/HSYMT/PAPER_figures/job50_share_vs_frequency_two_term.png", dpi=150); plt.close()
print("wrote PAPER_figures/job50_share_vs_frequency_two_term.png (frozen sky, pulsar 41, L = %.0f pc)" % L[41])
# ---------------- F2 -------------------------------------------------------------------
d = np.load("/home/mattm/projects/HSYMT/CW_transition/stagec_p2b_results.npz")
A_CW, sig_res = d["A_CW"], d["sig_res"]
rho = np.sqrt(2.0)*A_CW/sig_res; sphi = 1.0/rho
g = rho**2/2                                   # A^2 / total noise power
coh = np.sqrt(np.pi)/2*np.sqrt(g)*(i0e(g/2)+i1e(g/2))
print(f"F2 (116 pulsars): booked medians SNR_pterm {np.median(d['snr_pterm']):.3f}, sigma_phi {np.median(d['sigma_phi']):.3f} rad")
print(f"   corrected: pulsar-term SNR rho median {np.median(rho):.3f}; small-angle phase error 1/rho median {np.median(sphi):.3f} rad;"
      f" exact coherence <cos dphi> median {np.median(coh):.3f} (booked-formula coherence e^(-sigma^2/2) median {np.median(np.exp(-d['sigma_phi']**2/2)):.3f});"
      f" pulsars with rho > 1: {(rho > 1).sum()} (booked SNR_pterm > 1: {(d['snr_pterm'] > 1).sum()})")
OUT["F2"] = dict(snr_booked=float(np.median(d["snr_pterm"])), sphi_booked=float(np.median(d["sigma_phi"])),
                 rho=float(np.median(rho)), sphi=float(np.median(sphi)), coh=float(np.median(coh)),
                 n_rho_gt1=int((rho > 1).sum()), n_booked_gt1=int((d["snr_pterm"] > 1).sum()))
np.savez("job50_data/F2_stagec_p2b_corrected.npz", names=d["names"], rho=rho, sigma_phi_small=sphi, coherence=coh)
json.dump(OUT, open("job50_data/partF.json", "w"), indent=1)
