# JOB 52 follow-up: is the beta_unmod edge the CW upper limit or the constructed dominance (build_pop15)?
import os
for _v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np, pickle
import job50_common as C
C.J.load(verbose=False); NS = C.J.NS if hasattr(C.J, "NS") else C.NS
g = NS; SEED50 = g["SEED50"]; T_OBS = g["T_OBS"]; YR = g["YR_S"]
s30 = g["s_to_s30"](0.35); print(f"s_to_s30(0.35) = {s30:.4f}; T_OBS = {T_OBS/YR:.1f} yr; predicted single-term share at 45 yr = {1/(1+(1-s30)/s30/1.5):.3f}; predicted beta_unmod (own/(own+2*rest)) = {1/(1+2*(1-s30)/s30/1.5):.3f}")
nrej = []; rows = []
for j in range(100, 120):
    sd, f0, h0 = SEED50[j % 50][0], SEED50[j % 50][1], SEED50[j % 50][2]
    nrej.append(g["draw_population"](seed=sd)[5])
    Pn, sea_n, br, ok, ratio, h_n = g["build_pop15"](f0, h0, g["_sea50v"][j % 50], s30, 3, seed=sd)
    ul = g["h95_of_f"](f0); inb = (f0 >= 2e-9) & (f0 <= 3e-8)
    top = np.argsort(np.where(inb, h_n, 0))[::-1][:12]
    rows.append((j, sd, float((h0/ul)[inb].max()), float((h_n/ul)[inb].max()), int(((h_n/ul)[inb] > 1).sum()), [round(float(x), 2) for x in (h_n/h0)[br]], [round(float(x)*1e9, 2) for x in f0[br]], [round(float(x), 2) for x in (h_n/ul)[br]]))
print("upper-limit rejections in the raw draw, per sky:", nrej, "total", sum(nrej))
print("sky seed | max h/h95 raw draw | max h/h95 after build_pop15 | n sources above h95 after | boost h_new/h_old of the 3 constructed | their f (nHz) | their h/h95")
for r in rows: print("  ", r)
a = np.array([r[2] for r in rows]); b = np.array([r[3] for r in rows])
print(f"raw draw: max h/h95 median {np.median(a):.2f}, max {a.max():.2f}; skies within 20% of the limit: {(a>0.8).sum()}/20")
print(f"after construction: median {np.median(b):.2f}, max {b.max():.2f}; skies with a source above the limit: {(b>1).sum()}/20")
