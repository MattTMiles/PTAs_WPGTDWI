# ============ JOB 50 follow-up -- frames of the population, and what "array SNR" means ===========
# Numerical tests behind JOB50_MAIN.md section "Frames and mass convention".
import os, sys, json
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C
J = C.J; os.chdir(J.BASE); J.load(verbose=False); NS = J.NS
out = {}
SRC_MC, SRC_Z, SRC_AMP, SRC_DL = NS["SRC_MC"], NS["SRC_Z"], NS["SRC_AMP"], NS["SRC_DL"]
KH = NS["KH"]
# (1) strain normalisation uses the REDSHIFTED chirp mass and d_L (cell 3 l.100)
a_z = KH*(SRC_MC*(1+SRC_Z))**(5/3)/SRC_DL
a_s = KH*SRC_MC**(5/3)/SRC_DL
out["amp_matches_redshifted"] = bool(np.allclose(SRC_AMP, a_z))
out["amp_matches_source_frame"] = bool(np.allclose(SRC_AMP, a_s))
print("(1) SRC_AMP == KH (Mc(1+z))^(5/3)/d_L :", out["amp_matches_redshifted"],
      "| == KH Mc^(5/3)/d_L :", out["amp_matches_source_frame"])
# (2) presentation sky: mc_src is the INTRINSIC mass (cell 3 l.227), f_src is observed
f_src, h_src, mc_src, d_src = NS["f_src"], NS["h_src"], NS["mc_src"], NS["d_src"]
o = np.argsort(SRC_MC); k = np.clip(np.searchsorted(SRC_MC[o], mc_src), 0, len(o)-1)
out["mc_src_is_SRC_MC_member"] = float(np.mean(SRC_MC[o][k] == mc_src))
# recover z of each presentation source from (mc_src, d_src) and test the strain
key = {(m, d): z for m, d, z in zip(SRC_MC, SRC_DL, SRC_Z)}
z_src = np.array([key[(m, d)] for m, d in zip(mc_src, d_src)])
out["h_src_uses_redshifted_mass"] = bool(np.allclose(h_src, KH*(mc_src*(1+z_src))**(5/3)*f_src**(2/3)/d_src))
print("(2) mc_src drawn from SRC_MC (intrinsic):", out["mc_src_is_SRC_MC_member"],
      "| h_src == KH (mc_src(1+z))^(5/3) f_obs^(2/3)/d_L :", out["h_src_uses_redshifted_mass"])
# (3) chirp rate used for FP_ALL vs the observer-frame chirp rate
fac = (1+z_src)**(5/3)
w = h_src**2
out["z_src_median"] = float(np.median(z_src)); out["z_src_p16_p84"] = [float(x) for x in np.percentile(z_src, [16, 84])]
out["chirp_rate_underestimate_median"] = float(np.median(fac))
top = np.argsort(h_src)[::-1][:200]
out["z_top200_median"] = float(np.median(z_src[top])); out["chirp_rate_underestimate_top200"] = float(np.median(fac[top]))
print(f"(3) presentation sky N={len(f_src)}: z median {np.median(z_src):.2f} "
      f"[16,84%: {out['z_src_p16_p84'][0]:.2f},{out['z_src_p16_p84'][1]:.2f}]; "
      f"(1+z)^(5/3) median {np.median(fac):.2f}; loudest 200: z {np.median(z_src[top]):.2f}, factor {np.median(fac[top]):.2f}")
# (4) harness skies: is the amplitude -> (Mc, z) lookup (Job 41 recipe) exact?
oa = np.argsort(SRC_AMP); AS = SRC_AMP[oa]
mism = []; zc = []; rat = []; arr_w = []; arr_t = []
for j in range(0, 20):
    sd, f2, h2, P0, sea = NS["_pop18"](j % 50, 0.35)
    amp = h2/np.maximum(f2, 1e-30)**(2/3)
    kk = np.clip(np.searchsorted(AS, amp), 0, len(AS)-1)
    km = np.clip(kk-1, 0, len(AS)-1)
    best = np.where(np.abs(AS[km]-amp) < np.abs(AS[kk]-amp), km, kk)
    mism.append(np.abs(AS[kk]/amp-1)); mism.append(-np.abs(AS[best]/amp-1))   # sign-tag nearest
    r = NS["inv_rate"](f2, h2, P0, sea, 45, rate=1.0)
    if r is None: continue
    ks, sig, fk, rp, arr = r
    zc.append(SRC_Z[oa[kk[ks]]])
    # (5) array SNR behind dkappa: white timing noise only (job43c l.491-492) vs with the
    #     per-pulsar confusion ledger (other sources + sea in the bin) added in quadrature
    arr_w.append(arr); arr_t.append(np.sqrt(0.2*120)/sig); rat.append(arr/(np.sqrt(0.2*120)/sig))
m_used = np.concatenate(mism[0::2]); m_near = -np.concatenate(mism[1::2])
zc = np.concatenate(zc); arr_w = np.concatenate(arr_w); arr_t = np.concatenate(arr_t); rat = np.concatenate(rat)
out["lookup_rel_mismatch_used_median"] = float(np.median(m_used)); out["lookup_rel_mismatch_used_p99"] = float(np.percentile(m_used, 99))
out["lookup_exact_fraction_1e-9"] = float(np.mean(m_near < 1e-9))
out["z_booked12_median"] = float(np.median(zc)); out["chirp_factor_booked12_median"] = float(np.median((1+zc)**(5/3)))
print(f"(4) harness skies 0-19: amp->(Mc,z) lookup relative mismatch (as used, searchsorted): median {np.median(m_used):.2e}, "
      f"99% {np.percentile(m_used,99):.2e}; nearest-entry exact (<1e-9) fraction {np.mean(m_near<1e-9):.4f}")
print(f"    booked 12 sources/sky: z median {np.median(zc):.2f}, (1+z)^(5/3) median {np.median((1+zc)**(5/3)):.2f}")
out["arr_white_median"] = float(np.median(arr_w)); out["arr_with_confusion_median"] = float(np.median(arr_t))
out["arr_ratio_median"] = float(np.median(rat)); out["arr_ratio_p16_p84"] = [float(x) for x in np.percentile(rat, [16, 84])]
print(f"(5) booked 12 sources/sky, 20 skies, rate 1: array SNR white-only median {np.median(arr_w):.0f}; "
      f"with per-pulsar confusion ledger {np.median(arr_t):.1f}; ratio median {np.median(rat):.0f} "
      f"[16,84%: {out['arr_ratio_p16_p84'][0]:.0f},{out['arr_ratio_p16_p84'][1]:.0f}]")
json.dump(out, open("job50_data/frames.json", "w"), indent=1)
