# JOB 53 tests: (1) constructed switch + Joint53 reproduce Job 52 flag by flag on census sky 100;
# (2) natural switch: strains equal the raw draw, bin powers equal the raw sums; (3) what build_pop15 changed.
import os, sys, pickle, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np
import job50_common as C, job52_common as H52, job53_common as H
import job53_run as R
t0 = time.time(); C.J.load(verbose=False); NS = C.NS
j = 100
# ---- (1) -----------------------------------------------------------------------------------------
src = pickle.load(open(f"job52_data/parts/stat_c_{j:03d}.pkl", "rb"))[1]
rank = {d["sky"]: (d["cand"], d["arr_opt"]) for d in src}
arms = [("ref", n, "t1") for n in H.NLIST]+[("A", n, "t1") for n in H.N4]+[("ref", n, t) for n in H.N4 for t in ("t10", "tinf")]+[("off", n, "t1") for n in (12, "all")]
rec, info = C.collect50([j], hook=H.Joint53(rank, arms), arms=R.ORA, **R.BASE)
old = [r for r in pickle.load(open("job52_data/run_main_c.pkl", "rb"))["rec"] if r["sky"] == j]
oldoff = [r for r in pickle.load(open("job52_data/run_p5off_c.pkl", "rb"))["rec"] if r["sky"] == j]
assert len(rec) == len(old) == 79
for a in arms:
    k = f"w|{a[0]}|{a[1]}|{a[2]}"; ref = oldoff if a[0] == "off" else old
    kk = k.replace("w|off|", "w|ref|") if a[0] == "off" else k
    same = sum(r[k] == o[kk] for r, o in zip(rec, ref)); dm = max(abs(r["m"+k[1:]]-o["m"+kk[1:]]) for r, o in zip(rec, ref))
    print(f"(1) constructed sky {j}, arm {a}: flags equal Job 52 for {same}/79 pulsars; max |margin diff| {dm:.2e}")
    assert same == 79
print(f"(1) PASS  ({time.time()-t0:.0f} s)")
# ---- (2) -----------------------------------------------------------------------------------------
class Cap:
    def psr(self, loc): return {}
    def sky(self, loc):
        self.loc = {k: loc[k] for k in ("f2", "h2", "P0", "sea", "sd")}; return {}
cap = Cap(); C.collect50([j], hook=cap, arms=R.ORA, sky="natural", **R.BASE)
f, h, mc, dl, S, nrej = NS["draw_population"](seed=C.NAT_SEED0+j)
T_OBS = NS["T_OBS"]; be = NS["bin_edges"]; nbin = len(NS["bin_cent"])
L = cap.loc
print(f"(2) natural sky {j}: seed {L['sd']}, {len(f)} sources, {nrej} rejected by the upper limit")
print(f"    max |f - raw| {np.max(np.abs(L['f2']-f)):.1e}; max |h/raw - 1| {np.max(np.abs(L['h2']/h-1)):.1e}; max |P0/(h^2 f T) - 1| {np.max(np.abs(L['P0']/(h**2*f*T_OBS)-1)):.1e}; max |sea - raw sea| {np.max(np.abs(L['sea']-S)):.1e}")
b = np.digitize(f, be)-1; ok = (b >= 0) & (b < nbin)
raw = np.zeros(nbin); np.add.at(raw, b[ok], (h**2*f*T_OBS)[ok]); raw += S
got = np.zeros(nbin); np.add.at(got, b[ok], L["P0"][ok]); got += L["sea"]
print(f"    per-bin total power, decoder input / raw sum: max |ratio - 1| {np.max(np.abs(got[raw>0]/raw[raw>0]-1)):.1e} over {int((raw>0).sum())} bins")
assert np.array_equal(L["f2"], f) and np.array_equal(L["h2"], h) and np.array_equal(L["sea"], S)
# chirp-mass lookup is exact on natural skies
i_ = C._CACHE["ordA"][np.clip(np.searchsorted(C._CACHE["AS"], h/f**(2/3)), 0, len(C._CACHE["AS"])-1)]
print(f"    chirp-mass lookup (job50_common._mcz_tc) on natural sky: max |Mc_lookup/Mc_true - 1| {np.max(np.abs(NS['SRC_MC'][i_]/mc-1)):.1e}")
print("(2) PASS")
# ---- (3) what the construction changed, sky 100 (population 0) -------------------------------------------
cap2 = Cap(); C.collect50([j], hook=cap2, arms=R.ORA, sky="constructed", **R.BASE)
sd, f0, h0 = NS["SEED50"][0][0], NS["SEED50"][0][1], NS["SEED50"][0][2]
fr, hr, mcr, dlr, Sr, _ = NS["draw_population"](seed=sd)
assert np.array_equal(fr, f0) and np.array_equal(hr, h0)
Lc = cap2.loc; ch = np.where(np.abs(Lc["h2"]/h0-1) > 1e-9)[0]
inb = (f0 >= 2e-9) & (f0 <= 1e-8)
print(f"(3) constructed sky {j} = population seed {sd}: {len(ch)} of {len(f0)} sources have a strain different from the raw draw; raised: {int((Lc['h2'][ch] > h0[ch]).sum())}, lowered: {int((Lc['h2'][ch] < h0[ch]).sum())}")
up = ch[Lc["h2"][ch] > h0[ch]]
print(f"    raised sources: f (nHz) {np.round(f0[up]*1e9, 2)}, strain x {np.round(Lc['h2'][up]/h0[up], 2)}")
dn = ch[Lc["h2"][ch] < h0[ch]]
print(f"    lowered sources: {len(dn)}, strain x {np.round(np.percentile(Lc['h2'][dn]/h0[dn], [0, 50, 100]), 3)} (min/median/max)")
print(f"    sea bins changed: {int((np.abs(Lc['sea']/np.maximum(Sr, 1e-300)-1) > 1e-9).sum())} of {nbin}; sea total x {Lc['sea'].sum()/Sr.sum():.3f}")
i_ = C._CACHE["ordA"][np.clip(np.searchsorted(C._CACHE["AS"], Lc["h2"]/f0**(2/3)), 0, len(C._CACHE["AS"])-1)]
mz_true = mcr*(1+NS["SRC_Z"][C._CACHE["ordA"][np.clip(np.searchsorted(C._CACHE["AS"], hr/fr**(2/3)), 0, len(C._CACHE["AS"])-1)]])
mz_look = NS["SRC_MC"][i_]*(1+NS["SRC_Z"][i_])
print(f"    chirp mass (redshifted) used by the machine / the drawn source's: raised sources {np.round(mz_look[up]/mz_true[up], 2)}; chirp time x {np.round((mz_look[up]/mz_true[up])**(-5/3), 2)}")
print(f"TOTAL {time.time()-t0:.0f} s")
