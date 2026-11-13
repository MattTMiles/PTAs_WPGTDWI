# JOB 54: sources detectable against the background (version B, Earth terms only) that the candidate cap leaves out
import pickle, numpy as np, os
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import importlib.util, io, contextlib
spec = importlib.util.spec_from_file_location("sens", "job54_sens.py")
CAP = pickle.load(open("job54_data/sens_capture_sky100.pkl", "rb")); f2, h2 = CAP["f2"], CAP["h2"]
nt = 45*365.25/14; Rbg = 1.0+(CAP["gE"]+CAP["gP_mean"]+2*CAP["bsea"])*CAP["rt"]/(CAP["Ncal"]*CAP["fg"]**3)
bgp = (np.interp(f2, CAP["fg"], Rbg)-1)*CAP["Ncal"]*f2**3; own = CAP["P0"]*CAP["sc"]*CAP["rt"]
R = 1+np.maximum(bgp-own, 0)/(CAP["Ncal"]*f2**3)
snr = h2/(2*np.pi*f2)*np.sqrt(nt/2)/1e-7*np.sqrt(0.2*120)/np.sqrt(R)
inb = (f2 >= 2e-9) & (f2 <= 3e-8); isc = np.zeros(len(f2), bool); isc[CAP["cand"]] = True
d = np.where((snr >= 3) & inb)[0]; d = d[np.argsort(-snr[d])]
print(f"sky 100, 120 pulsars, Earth terms, noise = timing + all other GW power (Job 51 arr_bg definition): {len(d)} sources at SNR >= 3 in 2-30 nHz; in the 1,300-candidate list: {int(isc[d].sum())}; left out by the cap: {int((~isc[d]).sum())}")
print("  f [nHz] | SNR against background | timing-only single-pulsar SNR rp | rank by rp among in-band sources | in candidate list")
rp = CAP["rp_all"]; rk = np.argsort(np.argsort(-np.where(inb, rp, -1)))
for q in d: print(f"  {f2[q]*1e9:6.2f} | {snr[q]:5.2f} | {rp[q]:7.1f} | {rk[q]+1:5d} | {bool(isc[q])}")
print(f"candidate list: {isc.sum()} sources, frequencies 10/50/90% {np.percentile(f2[isc]*1e9, [10, 50, 90]).round(2)} nHz, max {f2[isc].max()*1e9:.2f} nHz; smallest rp kept {rp[isc].min():.1f}")
