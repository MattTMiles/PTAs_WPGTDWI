import os, sys, time, pickle
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
import job53_run as R
R.C.J.load(verbose=False)
for n in sys.argv[2:]:
    t0 = time.time(); R.work((n, int(sys.argv[1]))); print(n, f"{time.time()-t0:.0f} s", flush=True)
