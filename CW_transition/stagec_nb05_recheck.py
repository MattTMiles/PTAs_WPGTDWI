"""Stage C — D5 add-on: blast-radius closure for notebook 05 (the key science nb).

Rerun nb-05's joint-DeltalnL assessment with the CORRECTED direct-evaluation fringe
scanner (D3's fringe-center eval) instead of the prominence-based peak finder that
over-reports identification when the fringe modulation < the 0.5 prominence floor.

Same configuration as nb-05: 116 psr, seed 1234, realistic per-source draws
(h in [-14,-13.5], f in [-8,-7.5], mc [8.5,9.5]), GWB at NG15 (-14.6, 13/3), 1 us
white noise (equad -6), 14 GP components, NOISY simulate. Sweep N_CW = 1..10.

original joint dlnL = lnL(truth) - lnL(all psr at per-pulsar best-wrong mode), best
  wrong via find_best_wrong_mode_in_prior (prominence floor 0.5 -> nan when modulation
  small -> pulsar LEFT AT TRUTH -> joint drop understated -> separability overstated).
corrected joint dlnL = same, but best-wrong distance from direct lnL eval at fringe
  centers L0+k*dL (no prominence floor).

Run in jug-gpu:  python stagec_nb05_recheck.py
"""

import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys
import time

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import (
    load_pulsars, build_enterprise_pta, generate_injection_params, simulate,
    build_fast_scan_likelihood, compute_mode_spacing,
    compute_joint_best_wrong_in_prior, scan_pulsar_distance,
)

# nb-05 configuration (verbatim)
SEED = 1234
NCW_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
LOG10_H_RANGE = (-14.0, -13.5)
LOG10_FGW_RANGE = (-8.0, -7.5)
LOG10_MC_RANGE = (8.5, 9.5)
GWB_LOG10_A = -14.6
GWB_GAMMA = 13.0 / 3.0
LOG10_EQUAD = -6.0
N_COMPONENTS = 14

def best_wrong_dist_direct(logl_fn, bv, pk, dist_key, L0, sig_em, dL, n_sigma=3.0):
    """Best-wrong DISTANCE (others at truth) via a DIRECT FINE GRID scan, max beyond
    0.5 dL from truth. No find_peaks / no prominence floor (the bug). Uses a fine grid
    (not just fringe centers) because nb-05 is NOISY -> the best-wrong peak is shifted
    off the fringe centre; centres alone (D3's Asimov method) under-find it."""
    lo = max(1e-3, L0 - n_sigma * sig_em); hi = L0 + n_sigma * sig_em
    K = max(1, int((hi - lo) / dL))
    n_points = int(np.clip(8 * K, 2000, 8000))
    sd, sl = scan_pulsar_distance(logl_fn, bv, pk, dist_key, lo, hi, n_points,
                                  n_components=N_COMPONENTS)
    mask = np.abs(sd - L0) > 0.5 * dL
    if not mask.any():
        return np.nan
    return float(sd[mask][np.argmax(sl[mask])])


def corrected_joint(logl_fn, bv, pk, ent_psrs, disco_psrs, inj, cwb):
    """Corrected joint dlnL: direct fine-grid best-wrong per pulsar, set all, eval."""
    cw_src = [{"cos_gwtheta": inj[f"{n}_cos_gwtheta"], "gwphi": inj[f"{n}_gwphi"],
               "log10_fgw": inj[f"{n}_log10_fgw"]} for n in cwb]
    ll_truth = float(logl_fn(bv))
    wrong = np.array(bv)
    for pe, pd in zip(ent_psrs, disco_psrs):
        key = f"{pe.name}_cw_p_dist"
        if key not in pk:
            continue
        dL = min(compute_mode_spacing(c["cos_gwtheta"], c["gwphi"], c["log10_fgw"], pd.pos)
                 for c in cw_src)
        bw = best_wrong_dist_direct(logl_fn, bv, pk, key, pe.pdist[0], pe.pdist[1], dL)
        if np.isfinite(bw):
            wrong[pk.index(key)] = bw
    ll_wrong = float(logl_fn(jnp.asarray(wrong)))
    return ll_truth - ll_wrong


def run():
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)
    ent_psrs, disco_psrs = load_pulsars(116)
    rows = []
    for ncw in NCW_VALUES:
        t0 = time.time()
        rng = np.random.default_rng(SEED); np.random.seed(SEED)
        pta, cwb, _ = build_enterprise_pta(ent_psrs, ncw, components=N_COMPONENTS,
                                           log10_equad=LOG10_EQUAD)
        inj = generate_injection_params(
            pta, ent_psrs, ncw, cwb, log10_h=None, scenario="realistic", rng=rng,
            log10_h_range=LOG10_H_RANGE, log10_fgw_range=LOG10_FGW_RANGE,
            log10_mc_range=LOG10_MC_RANGE, gwb_log10_A=GWB_LOG10_A, gwb_gamma=GWB_GAMMA)
        sim = simulate(pta, inj)
        resid_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim)}
        logl_fn, pk, bv = build_fast_scan_likelihood(
            disco_psrs, resid_map, ncw, inj, cwb,
            components=N_COMPONENTS, log10_equad=LOG10_EQUAD)

        orig = compute_joint_best_wrong_in_prior(
            logl_fn, bv, pk, ent_psrs, disco_psrs, inj,
            cw_block_names=cwb, n_sigma=3.0, n_components=N_COMPONENTS)
        orig_joint = orig["delta_lnL"]
        corr_joint = corrected_joint(logl_fn, bv, pk, ent_psrs, disco_psrs, inj, cwb)
        rows.append((ncw, orig_joint, corr_joint))
        print(f"  N_CW={ncw}: original joint dlnL={orig_joint:.2f}  "
              f"corrected joint dlnL={corr_joint:.2f}  ({time.time()-t0:.1f}s)", flush=True)

    print("\n========== nb-05 joint dlnL: original vs corrected ==========")
    print(f"{'N_CW':>4s} {'original':>12s} {'corrected':>12s} {'ratio corr/orig':>16s}")
    for ncw, o, c in rows:
        r = c / o if o != 0 else np.nan
        print(f"{ncw:>4d} {o:>12.2f} {c:>12.2f} {r:>16.3f}")
    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_nb05_recheck.npz",
             ncw=np.array([r[0] for r in rows]),
             orig=np.array([r[1] for r in rows]),
             corr=np.array([r[2] for r in rows]))
    print("saved stagec_nb05_recheck.npz")
    return 0


if __name__ == "__main__":
    sys.exit(run())
