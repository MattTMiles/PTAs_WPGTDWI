"""Stage C — Deliverable D5: frozen-GWB sensitivity (caveat 7.2).

At N_CW=8 equal strain (the pre-onset point, most sensitive), recompute the
D2/D3 pipeline with the frozen GWB log-amplitude at truth and truth +/- 0.5 dex.
The GWB GP is frozen at truth in the fast-scan / amortised Fisher; this bounds how
optimistic that freezing is. The Asimov data (pure CW injection) is GWB-independent;
only the likelihood's frozen GWB prior amplitude moves, changing how much low-
frequency power the GP absorbs -> the noise weighting -> sigma_L and dlnL.

Report movement in: median marginal sigma_L (Fisher-valid), median dlnL_a, class-i.

Run in jug-gpu:  python stagec_d5.py --ndraws 10
"""

import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse
import sys
import time

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import load_pulsars, build_enterprise_pta, compute_mode_spacing
from stagec_fisher import (
    build_fisher_amortised, make_geometry_injection,
    N_COMPONENTS, LOG10_EQUAD, GWB_LOG10_A, KPC_TO_PC,
)
from stagec_d4 import run_one_draw, EQUAL_H

NCW = 8
SHIFTS = [0.0, +0.5, -0.5]


def run(ndraws):
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)
    ent_psrs, disco_psrs = load_pulsars(116)
    em_pc = np.array([p.pdist[1] for p in ent_psrs]) * KPC_TO_PC
    sig_em_kpc = np.array([p.pdist[1] for p in ent_psrs])
    pta, cwb, _ = build_enterprise_pta(ent_psrs, NCW, components=N_COMPONENTS,
                                       log10_equad=LOG10_EQUAD)

    out = {}
    for shift in SHIFTS:
        gwb = GWB_LOG10_A + shift
        print(f"\n===== frozen GWB log10_A = {gwb:.2f} (truth{shift:+.1f} dex) =====", flush=True)
        inj0 = make_geometry_injection(pta, ent_psrs, NCW, cwb, seed=1000,
                                       h_range=(EQUAL_H, EQUAL_H))
        inj0["gwb_log10_A"] = gwb
        t0 = time.time()
        amo = build_fisher_amortised(disco_psrs, NCW, inj0, cwb)
        th0 = amo["theta_truth"]; dt0 = amo["inject_data"](th0)
        from stagec_d4 import assemble_hessian, fringe_dlnL
        _ = assemble_hessian(amo, th0, dt0)
        _ = fringe_dlnL(amo, th0, dt0, 0, ent_psrs[0].pdist[0], ent_psrs[0].pdist[1],
                        max(compute_mode_spacing(inj0[f"{cwb[0]}_cos_gwtheta"],
                            inj0[f"{cwb[0]}_gwphi"], inj0[f"{cwb[0]}_log10_fgw"],
                            disco_psrs[0].pos), 1e-6))
        print(f"  build+compile {time.time()-t0:.1f}s", flush=True)

        med_marg, med_dlnL, ci = [], [], []
        for d in range(ndraws):
            inj = make_geometry_injection(pta, ent_psrs, NCW, cwb, seed=2000 + d,
                                          h_range=(EQUAL_H, EQUAL_H))
            r = run_one_draw(amo, ent_psrs, disco_psrs, inj, cwb, em_pc, sig_em_kpc)
            v = r["valid"]
            med_marg.append(np.median(r["marg"][v]) if v.any() else np.nan)
            fin = np.isfinite(r["dlnL"])
            med_dlnL.append(np.median(r["dlnL"][fin]))
            ci.append(int((r["cls"] == 1).sum()))
            print(f"  draw {d}: med_marg={med_marg[-1]:.4g}pc med_dlnL={med_dlnL[-1]:.3f} "
                  f"class_i={ci[-1]}", flush=True)
        out[shift] = dict(
            med_marg=float(np.nanmedian(med_marg)),
            med_dlnL=float(np.median(med_dlnL)),
            ci_med=int(np.median(ci)), ci_min=int(min(ci)), ci_max=int(max(ci)),
        )

    print("\n================ D5 SUMMARY (N_CW=8, GWB freeze sensitivity) ================")
    print(f"{'GWB dex':>8s} {'med marg sigma_L pc':>20s} {'med dlnL_a':>12s} {'class_i med(min,max)':>22s}")
    base = out[0.0]
    for shift in SHIFTS:
        o = out[shift]
        print(f"{shift:+8.1f} {o['med_marg']:>20.4g} {o['med_dlnL']:>12.3f} "
              f"{str((o['ci_med'], o['ci_min'], o['ci_max'])):>22s}")
    # movement relative to truth
    dm_p = out[+0.5]['med_marg'] / base['med_marg']
    dm_m = out[-0.5]['med_marg'] / base['med_marg']
    dd_p = out[+0.5]['med_dlnL'] / max(base['med_dlnL'], 1e-9)
    dd_m = out[-0.5]['med_dlnL'] / max(base['med_dlnL'], 1e-9)
    print(f"\nmovement vs truth: marg sigma_L x{dm_m:.2f} (-0.5dex) .. x{dm_p:.2f} (+0.5dex)")
    print(f"                   med dlnL_a x{dd_m:.2f} (-0.5dex) .. x{dd_p:.2f} (+0.5dex)")
    print(f"cross-ref D3: timing+GP absorption already costs ~1.6x in pulsar-term SNR "
          f"(-> ~2.6x in info/dlnL, ~1.6x in sigma_L).")
    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_d5_results.npz",
             shifts=np.array(SHIFTS),
             med_marg=np.array([out[s]['med_marg'] for s in SHIFTS]),
             med_dlnL=np.array([out[s]['med_dlnL'] for s in SHIFTS]),
             ci_med=np.array([out[s]['ci_med'] for s in SHIFTS]))
    print("saved stagec_d5_results.npz")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ndraws", type=int, default=10)
    sys.exit(run(ap.parse_args().ndraws))
