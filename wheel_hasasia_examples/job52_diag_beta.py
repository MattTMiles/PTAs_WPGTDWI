# JOB 52 follow-up: what limits beta_unmod? Component breakdown of the unmodelled power (N = 12, reference).
import os, sys, pickle, time
for _v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np
import job50_common as C, job52_common as H, job52_run as R
sky = int(sys.argv[1]); t0 = time.time(); C.J.load(verbose=False)
out = {}
for nm in ("p5_c", "p5off_c"):
    _, kw, hk, rk = R.RUNS[nm]; kw = dict(kw)
    src = pickle.load(open(f"job52_data/parts/{rk}_{sky:03d}.pkl", "rb"))[1]
    hook = H.Joint52({d["sky"]: (d["cand"], d["arr_opt"]) for d in src}, ext=hk[1], mode=hk[2], ret=hk[3], part4=True)
    hook.p4_debug = True; kw["hook"] = hook
    out[nm] = C.collect50([sky], **kw)
    print(nm, sky, f"{time.time()-t0:.0f} s", flush=True)
pickle.dump(out, open(f"job52_data/diag_beta_{sky}.pkl", "wb"))
