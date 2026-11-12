"""JOB 46 Part A.0 -- reproduce the booked prong-2 numbers at N = 1, 100, 398(~410), 1000.
Same seeds, same array, same machinery (prong2_transition.run_sweep), 30 sky realisations."""
import sys, time, numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import prong2_transition as P
import jax.numpy as jnp

d = np.load("/home/mattm/projects/HSYMT/CW_transition/prong2_results.npz", allow_pickle=True)
print("=== BOOKED CONFIG (prong2_results.npz) ===")
print(f"N_psr={int(d['ps_n_psr'])}  T={float(d['ps_T'])/P.YR:.1f} yr  cadence={float(d['ps_cadence_days'])} d  "
      f"sigma_white={float(d['ps_sigma'])*1e9:.0f} ns  band={d['ps_f_band']*1e9} nHz  "
      f"n_seeds(sky realisations)={int(d['ps_n_seeds'])}  n_max={int(d['ps_n_max'])}  rcond={float(d['prod_rcond']):.0e}")
print(f"N grid: {d['ps_n'].tolist()}")
print(f"T*Df = {float(d['N_star']):.3f}   knee(first 0.5 crossing) = {float(d['knee_N']):.1f}")
print("pulsar sky: build_array(116, seed=777) -> isotropic SYNTHETIC directions, all L = 1 kpc (not a catalogue)")
n_toa = int(15*P.YR/(7*86400)); print(f"n_toa per pulsar = {n_toa}")

n_vals = [1, 100, 398, 1000]
p_hats = P.build_array(116)
t0 = time.time()
res_ps = P.run_sweep(n_vals, mode="fixed_persource", n_max=1000, p_hats=p_hats, log10_h=-14.3)
res_tot = P.run_sweep(n_vals, mode="fixed_total", n_max=1000, p_hats=p_hats, h_total=5e-14)
print(f"[timing] 2 modes x 4 N x 30 skies: {time.time()-t0:.1f}s")
print("=== REPRODUCTION vs BOOKED ===")
for lab, r, pfx in (("fixed_persource", res_ps, "ps"), ("fixed_total", res_tot, "tot")):
    for i, n in enumerate(n_vals):
        j = list(d[f"{pfx}_n"]).index(n)
        rat_new = r["marg"][i]/r["cond"][i]; rat_old = d[f"{pfx}_marg"][j]/d[f"{pfx}_cond"][j]
        print(f"  [{lab:15s}] N={n:4d}  cond new/booked = {r['cond'][i]:.4e}/{d[pfx+'_cond'][j]:.4e}  "
              f"marg new/booked = {r['marg'][i]:.4e}/{d[pfx+'_marg'][j]:.4e}  ratio new/booked = {rat_new:.3f}/{rat_old:.3f}")
