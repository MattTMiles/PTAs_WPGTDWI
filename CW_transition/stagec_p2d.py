"""Prong-2 close-out P2-D: loose ends (two bounded items).

Item 1 (strain reconciliation): per-source OPTIMAL SNR (full Earth+pulsar CW residual,
  marginalised noise) for the D4 draws (h=-13.75, equal strain) vs the old prong-1
  optimizer / nb-03 loud injections. The CW_node_sampling optimizer sandboxes (MK5/MK7)
  inject log10_h = -12.0; nb-03 uses -13.0; nb-05/D4 use -13.75..-13.5. Report the ratio
  -> closes the D6 "were the old runs louder" question.
  SNR_s = sqrt(2[logL(theta_off,0) - logL(theta_off, r_full_s)]) (reuses the amortised
  likelihood; theta_off kills the CW model so the residual IS the data).

Item 2 (nb-01 rerun): nb 01_single (single CW, scenario="single", log10_fgw=-8, mc=9,
  NOISY simulate, seed 1234) at h in {-14,-13.5,-13}. Original joint dlnL
  (compute_joint_best_wrong_in_prior, prominence peak finder) vs corrected joint dlnL
  (D3 direct fine-grid scanner). 00/02 stay flagged-not-rerun per the D6 decision.

Run in jug-gpu:  python stagec_p2d.py
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
    build_fast_scan_likelihood, compute_joint_best_wrong_in_prior,
)
from discovery.deterministic import make_phase_connected_binary
from stagec_fisher import (
    build_fisher_amortised, make_geometry_injection,
    N_COMPONENTS, LOG10_EQUAD, CW_PARAM_BASE,
)
from stagec_nb05_recheck import corrected_joint

EQUAL_H = -13.75
_FULL = jax.vmap(make_phase_connected_binary(pulsarterm=True),
                 in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0, None))
_KEYS = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
         "log10_fgw", "log10_h", "phase0", "psi")


def full_cw_data(disco_psrs, src):
    st = {k: jnp.array([src[k]]) for k in _KEYS}
    out = []
    for p in disco_psrs:
        d = _FULL(p.toas, p.pos, st["cos_gwtheta"], st["gwphi"], st["cos_inc"],
                  st["log10_mc"], st["log10_fgw"], st["log10_h"], st["phase0"],
                  st["psi"], p.pdist[0])
        out.append(jnp.sum(d, axis=0))
    return tuple(out)


def item1_strain():
    print("\n########## ITEM 1: strain reconciliation ##########", flush=True)
    ent, disco = load_pulsars(116)
    pta1, cwb1, _ = build_enterprise_pta(ent, 1, components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
    inj1 = make_geometry_injection(pta1, ent, 1, cwb1, seed=1000, h_range=(EQUAL_H, EQUAL_H))
    t0 = time.time()
    amo = build_fisher_amortised(disco, 1, inj1, cwb1)
    tk = amo["theta_keys"]
    off = np.array(amo["theta_truth"]); off[tk.index("cw_log10_h")] = -30.0
    off = jnp.asarray(off)
    zero = tuple(jnp.zeros(len(p.toas)) for p in disco)
    ll0 = float(amo["fisher_logL"](off, zero))
    print(f"  build {time.time()-t0:.1f}s", flush=True)

    def snr(src):
        ll = float(amo["fisher_logL"](off, full_cw_data(disco, src)))
        return float(np.sqrt(max(2 * (ll0 - ll), 0.0)))

    # collect the D4 equal-strain sources (10 draws x N_CW), take per-source SNR at h=-13.75
    def snrs_for_h(h, seeds=range(2000, 2010), ncw=8):
        pta, cwb, _ = build_enterprise_pta(ent, ncw, components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
        vals = []
        for s in seeds:
            inj = make_geometry_injection(pta, ent, ncw, cwb, seed=s, h_range=(h, h))
            for name in cwb:
                src = {k: inj[f"{name}_{k}"] for k in CW_PARAM_BASE}
                vals.append(snr(src))
        return np.array(vals)

    print("  computing per-source optimal SNR at each injection strain ...", flush=True)
    snr_d4 = snrs_for_h(-13.75)     # D4 / nb05-ish
    snr_nb03 = snrs_for_h(-13.0)    # nb-03 loud
    snr_opt = snrs_for_h(-12.0)     # CW_node_sampling optimizer (MK5/MK7)
    med = lambda a: float(np.median(a))
    print("\n  per-source OPTIMAL SNR (full Earth+pulsar, marginalised noise), median over 80 sources:")
    print(f"    h=-13.75 (D4 equal strain)       : {med(snr_d4):8.1f}")
    print(f"    h=-13.0  (nb-03 loud)            : {med(snr_nb03):8.1f}")
    print(f"    h=-12.0  (CW_node_sampling opt)  : {med(snr_opt):8.1f}")
    print(f"\n  RATIO old/new per-source SNR:")
    print(f"    optimizer(-12.0) / D4(-13.75) = {med(snr_opt)/med(snr_d4):6.1f}x  "
          f"(strain ratio 10^1.75 = {10**1.75:.1f})")
    print(f"    nb-03(-13.0)     / D4(-13.75) = {med(snr_nb03)/med(snr_d4):6.1f}x  "
          f"(strain ratio 10^0.75 = {10**0.75:.1f})")
    print("  => the old optimizer runs (h=-12, often noise-free) were ~50x louder per")
    print("     source than the realistic D4 injections; SNR scales linearly with strain,")
    print("     so their 20/20 distance recovery was a far easier regime (SNR ~ hundreds).")
    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_p2d_item1.npz",
             snr_d4=snr_d4, snr_nb03=snr_nb03, snr_opt=snr_opt)
    return med(snr_d4), med(snr_nb03), med(snr_opt)


def item2_nb01():
    print("\n########## ITEM 2: nb-01 joint-dlnL recheck (corrected scanner) ##########", flush=True)
    ent, disco = load_pulsars(116)
    rows = []
    for h in (-14.0, -13.5, -13.0):
        t0 = time.time()
        rng = np.random.default_rng(1234); np.random.seed(1234)
        pta, cwb, _ = build_enterprise_pta(ent, 1, components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
        inj = generate_injection_params(pta, ent, 1, cwb, log10_h=h, scenario="single", rng=rng)
        sim = simulate(pta, inj)
        resid = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim)}
        logl, pk, bv = build_fast_scan_likelihood(disco, resid, 1, inj, cwb,
                                                  components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
        orig = compute_joint_best_wrong_in_prior(logl, bv, pk, ent, disco, inj,
                                                 cw_block_names=cwb, n_sigma=3.0,
                                                 n_components=N_COMPONENTS)["delta_lnL"]
        corr = corrected_joint(logl, bv, pk, ent, disco, inj, cwb)
        rows.append((h, orig, corr))
        print(f"  nb-01 single CW h={h}: original joint dlnL={orig:.2f}  "
              f"corrected={corr:.2f}  ({time.time()-t0:.1f}s)", flush=True)
    print("\n  nb-01 original vs corrected joint dlnL:")
    print(f"  {'log10_h':>8s} {'original':>12s} {'corrected':>12s}")
    for h, o, c in rows:
        print(f"  {h:>8.1f} {o:>12.2f} {c:>12.2f}")
    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_p2d_item2.npz",
             h=np.array([r[0] for r in rows]),
             orig=np.array([r[1] for r in rows]), corr=np.array([r[2] for r in rows]))
    return rows


def main():
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)
    item1_strain()
    item2_nb01()
    print("\nsaved stagec_p2d_item1.npz, stagec_p2d_item2.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
