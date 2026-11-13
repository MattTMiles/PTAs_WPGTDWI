# JOB 53: how often does the amplitude lookup in job50_common._mcz_tc return the wrong catalogue entry?
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np
import job50_common as C
C.J.load(verbose=False); NS = C.NS
C._mcz_tc(np.array([1e-8]), np.array([1e-15]))
oA = C._CACHE["ordA"]; AS = C._CACHE["AS"]; MC = NS["SRC_MC"]; Z = NS["SRC_Z"]; DL = NS["SRC_DL"]
print(f"catalogue: {len(AS)} entries; exact duplicate amplitudes: {int((np.diff(AS) == 0).sum())}")
for label, seeds in (("raw draws behind constructed skies (seeds 3-22)", range(3, 23)), ("natural skies 100-119 (seeds 5100-5119)", range(5100, 5120))):
    bad = tot = 0; rat = []; badc = totc = 0; ratc = []
    for sd in seeds:
        f, h, mc, dl, S, _ = NS["draw_population"](seed=sd)
        k = np.clip(np.searchsorted(AS, h/f**(2.0/3.0)), 0, len(AS)-1); i = oA[k]
        wrong = (MC[i] != mc) | (DL[i] != dl)
        # recover the true entry: neighbour whose (Mc, d_L) match exactly
        true_i = np.full(len(f), -1)
        for dk in (0, -1, 1, -2, 2):
            kk = np.clip(k+dk, 0, len(AS)-1); ii = oA[kk]; m = (true_i < 0) & (MC[ii] == mc) & (DL[ii] == dl); true_i[m] = ii[m]
        assert (true_i >= 0).all()
        mzr = (MC[i]*(1+Z[i]))/(MC[true_i]*(1+Z[true_i]))
        bad += int(wrong.sum()); tot += len(f); rat += list(mzr[wrong])
        nt = 45*365.25/14.0; rp = (h/(2*np.pi*f))*np.sqrt(nt/2)/1e-7; cand = (rp*np.sqrt(24) >= 3) & (rp >= 3) & (f >= 2e-9) & (f <= 3e-8)
        badc += int((wrong & cand).sum()); totc += int(cand.sum()); ratc += list(mzr[wrong & cand])
    rat = np.array(rat); ratc = np.array(ratc)
    print(f"{label}: wrong entry for {bad} of {tot} sources ({bad/tot:.3f}); among candidates {badc} of {totc} ({badc/totc:.3f})")
    print(f"   redshifted chirp mass used / true, wrong ones: 10/50/90% {np.percentile(rat, [10, 50, 90]).round(2)}; chirp time used / true: {np.percentile(rat**(-5/3), [10, 50, 90]).round(2)}; candidates: {np.percentile(ratc**(-5/3), [10, 50, 90]).round(2) if len(ratc) else '-'}")
