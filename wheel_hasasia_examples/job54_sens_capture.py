# JOB 54: capture the first natural sky's (sky 100) background fields for the detection-reach curves (version B)
import os, sys, pickle
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job52_common as H52, job54_run as R
C.J.load(verbose=False)
class Cap(H52.Stat52):
    gP = []
    def psr(self, loc):
        self.gP.append(loc["gPg"]*loc["sc"]); return super().psr(loc)
    def sky(self, loc):
        out = super().sky(loc); fg = loc["fg8"]
        self.cap = dict(fg=fg, gE=loc["gEs"](fg), gP_mean=np.mean(self.gP, 0), gP_med=np.median(self.gP, 0), bsea=np.atleast_1d(loc["bsea"](fg)), rt=loc["rt"], sc=loc["sc"], Ts=loc["Ts"],
                        Ncal=float((loc["Pc"]*loc["sc"]*2.0/loc["rho_off"]**2/loc["fc"]**3)[0]), f2=loc["f2"], h2=loc["h2"], P0=loc["P0"], tc=loc["tc_all"], cand=loc["cand"],
                        arr_all=loc["arr_all"], rp_all=loc["rp_all"], arr_noise=out["arr_noise"], arr_bg=out["arr_bg"], arr_opt=out["arr_opt"], nclean=out["nclean"])
        return out
kw = dict(R.RUNS["stat15"][0]); kw.pop("T"); kw.pop("noskip"); h = Cap(store=False)
C.collect50([100], hook=h, fgrid_max=1.05e-7, **kw)
old = pickle.load(open("job53_data/parts/stat_nat_100.pkl", "rb"))[1][0]
print("arr_bg equal to Job 53's ranking run (max rel diff):", float(np.max(np.abs(h.cap["arr_bg"]/old["arr_bg"]-1))), "; arr_opt:", float(np.max(np.abs(h.cap["arr_opt"]/old["arr_opt"]-1))))
print("detectable on sky 100: arr_bg >= 3:", int((h.cap["arr_bg"] >= 3).sum()), "; arr_bg_opt >= 3:", int((h.cap["arr_opt"] >= 3).sum()), "; timing only >= 3:", int((h.cap["arr_noise"] >= 3).sum()), "; all sources timing-only arr >= 3:", int((h.cap["arr_all"] >= 3).sum()))
pickle.dump(h.cap, open("job54_data/sens_capture_sky100.pkl", "wb"))
