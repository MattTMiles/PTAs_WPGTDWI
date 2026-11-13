# JOB 53 Part 4: why are skies missing at short span?  (collect50: `if r is None: continue` after inv_rate)
import os, pickle
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np
import job50_common as C
C.J.load(verbose=False); NS = C.NS
for T in (15, 25, 45):
    have = {d["sky"] for d in pickle.load(open(f"job53_data/run_{'main' if T == 45 else 'span%d' % T}_nat.pkl", "rb"))["info"]}
    miss = [j for j in range(100, 300) if j not in have]; none = 0; ncand = []; ngood = []
    for j in miss:
        sd, f2, h2, P0, sea = C.natural_pop(j)
        r = NS["inv_rate"](f2, h2, P0, sea, T, rate=1.0); none += r is None
        nt = T*365.25/14.0; rp = (h2/(2*np.pi*f2))*np.sqrt(nt/2)/1e-7
        ncand.append(int(((rp*np.sqrt(24) >= 3) & (rp >= 3) & (f2 >= 2e-9) & (f2 <= 3e-8)).sum()))
    print(f"T = {T} yr: {len(miss)} skies missing; inv_rate returns None for {none} of them; candidates passing the timing-noise cut on those skies: min {min(ncand) if ncand else '-'}, median {np.median(ncand) if ncand else '-'}")
