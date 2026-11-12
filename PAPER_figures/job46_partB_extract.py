"""JOB 46 Part B.0 -- extract the frozen machine's presentation sky (draw_population(seed=3),
active_wheel_demo_v4.ipynb cells 1+3) and the audit's retarded pulsar-term frequencies
(cell 14 lines 1585-1590: rng 404, L ~ U(0.5,3) kpc, cos mu ~ U(-0.9,0.9), pulsar 41)."""
import os, json, sys, time, numpy as np
os.environ["MPLBACKEND"] = "Agg"
BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
os.chdir(BASE)
t0 = time.time()
nb = json.load(open("active_wheel_demo_v4.ipynb"))
NS = {"__name__": "__main__"}
c1 = "".join(nb["cells"][1]["source"]).replace(
    '_JUG_PY = "/home/mattm/miniforge3/envs/jug-gpu/bin/python"', '_JUG_PY = "/nonexistent"')
assert '_JUG_PY = "/nonexistent"' in c1
exec(compile(c1, "<cell1>", "exec"), NS)
import matplotlib.pyplot as plt
plt.close("all")
exec(compile("".join(nb["cells"][3]["source"]), "<cell3>", "exec"), NS)
plt.close("all")
f_src, h_src, mc_src = NS["f_src"], NS["h_src"], NS["mc_src"]
N_DISC = len(f_src)
G_SI, C_SI, C_PC_S, T_OBS, YR_S = NS["G_SI"], NS["C_SI"], NS["C_PC_S"], NS["T_OBS"], NS["YR_S"]
EXP = (120, 100)                                             # cell 14 line 115
_rngA = np.random.default_rng(404)
L_AUD = _rngA.uniform(500.0, 3000.0, EXP[0])
_mu_all = np.arccos(_rngA.uniform(-0.9, 0.9, (EXP[0], N_DISC)))
_dt_all = L_AUD[:, None]*(1 - np.cos(_mu_all))/C_PC_S
_KCH_all = (96/5)*np.pi**(8/3)*(G_SI*mc_src/C_SI**3)**(5/3)
FP_ALL = (f_src[None, :]**(-8/3) + (8/3)*_KCH_all[None, :]*_dt_all)**(-3/8)
i = 41
print(f"presentation sky: draw_population(seed=3): N_DISC={N_DISC}  T_OBS={T_OBS/YR_S:.0f} yr  "
      f"F_RANGE={NS['F_RANGE']}  1/T={1e9/T_OBS:.3f} nHz  n_bins={len(NS['bin_cent'])}")
print(f"pulsar 41: L={L_AUD[i]:.0f} pc; cos mu in (-0.9,0.9); look-back dt median {np.median(_dt_all[i])/YR_S:.0f} yr")
print("(0d) table reproduction: Earth vs pulsar terms per band (pulsar 41)")
for lo, hi in ((2e-9, 5e-9), (5e-9, 1e-8), (1e-8, 2e-8), (2e-8, 5e-8), (5e-8, 1e-7)):
    ne = int(((f_src >= lo) & (f_src < hi)).sum()); npu = int(((FP_ALL[i] >= lo) & (FP_ALL[i] < hi)).sum())
    print(f"    {lo*1e9:5.0f}-{hi*1e9:5.0f} nHz: Earth {ne:6d}   pulsar {npu:6d}")
np.savez("/tmp/claude-1005/-home-mattm-soft/3b8c1241-329a-43ca-acaf-8dccedc6cfd9/scratchpad/partB_sky.npz",
         f_src=f_src, h_src=h_src, mc_src=mc_src, fp41=FP_ALL[i], L41=L_AUD[i], mu41=_mu_all[i],
         bin_edges=NS["bin_edges"], S_faint=NS["S_FAINT"], T_OBS=T_OBS, F_RANGE=np.array(NS["F_RANGE"]), fp_all=FP_ALL, L_all=L_AUD)
print(f"[timing] extraction {time.time()-t0:.1f}s")
