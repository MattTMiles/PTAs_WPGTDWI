import numpy as np, glob
for f in sorted(glob.glob("/home/mattm/projects/HSYMT/SPARK2_results/s2b_r*_k*.npz")):
    z = np.load(f, allow_pickle=True)
    print(f"{f.split('/')[-1]}: rung={int(z['rung'])} coh={int(z['n_coherent'])} "
          f"med_cond_sig={np.median(z['marg_sig']):.4g} "
          f"off_sig un={float(z['off_sig_un']):.3f} mo={float(z['off_sig_mo']):.3f} "
          f"floor_un={float(z['floor_un']):.3f}(zf={float(z['floor_un_zf']):.2f}) "
          f"floor_mo={float(z['floor_mo']):.3f}(zf={float(z['floor_mo_zf']):.2f}) "
          f"N_un={int(z['ncert_un'])} N_mo={int(z['ncert_mo'])}")
for f in ["spark2_mask.npz", "spark2_trials.npz"]:
    try:
        z = np.load(f"/home/mattm/projects/HSYMT/SPARK2_results/{f}", allow_pickle=True)
        print(f"\n{f}: " + ", ".join(f"{k}={np.asarray(z[k]).ravel()[:5]}" for k in z.files if k != "note"))
    except Exception as e: print(f, e)
