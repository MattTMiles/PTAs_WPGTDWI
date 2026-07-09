"""Prong-2 close-out P2-A: DETECTABILITY vs RANGEABILITY (the headline figure).

Two curves vs N_CW on the SAME D4 draws (116 psr, equal strain h=-13.75, 10 Asimov
geometry draws; + the 3-loud+13-faint population draw):

  DETECTABLE = # sources individually detectable via EARTH-TERM-ONLY information.
    Per source s, the matched-filter SNR of its Earth term with ALL pulsar terms
    omitted from the model, against the real marginalised noise (white+timing+frozen
    RN+GWB). SNR_s^2 = (s|s) = r_earth_s^T Sigma^-1 r_earth_s, obtained from the
    amortised zero-noise likelihood: with the CW model amplitude killed (log10_h=-30)
    the residual IS the data, so
        SNR_s^2 = 2*[ logL(theta_off, data=0) - logL(theta_off, data=r_earth_s) ].
    Count SNR_s > 5 (also 3 and 8 as sensitivity).

  RANGEABLE = the D4 class-i count (fringe-identified AND marginal sigma_L < EM prior),
    loaded from stagec_d4_results.npz.

GATE: for one loud source, the Earth-term (s|s) SNR must reproduce the standard
white-noise matched-filter CW SNR sqrt(sum (r_earth/sigma)^2) within ~20%.

Run in jug-gpu:  python stagec_p2a.py
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
from cw_helpers import load_pulsars, build_enterprise_pta
from discovery.deterministic import make_phase_connected_binary
from stagec_fisher import (
    build_fisher_amortised, make_geometry_injection,
    N_COMPONENTS, LOG10_EQUAD, CW_PARAM_BASE,
)

EQUAL_H = -13.75
NCW_LIST = [1, 2, 4, 8, 16]
D4_NPZ = "/home/mattm/projects/HSYMT/CW_transition/stagec_d4_results.npz"

# Earth-term-only CW delay (pulsarterm=False): 10 positional args, no p_dist
_EARTH = jax.vmap(make_phase_connected_binary(pulsarterm=False),
                  in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0))
_KEYS = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
         "log10_fgw", "log10_h", "phase0", "psi")


def earth_term_data(disco_psrs, src_params):
    """Per-pulsar Earth-term-only delay for ONE source (tuple over disco_psrs)."""
    st = {k: jnp.array([src_params[k]]) for k in _KEYS}
    out = []
    for p in disco_psrs:
        d = _EARTH(p.toas, p.pos, st["cos_gwtheta"], st["gwphi"], st["cos_inc"],
                   st["log10_mc"], st["log10_fgw"], st["log10_h"], st["phase0"], st["psi"])
        out.append(jnp.sum(d, axis=0))
    return tuple(out)


def src_list_from_inj(inj, cwb):
    return [{k: inj[f"{name}_{k}"] for k in CW_PARAM_BASE} for name in cwb]


def main():
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)
    ent_psrs, disco_psrs = load_pulsars(116)

    # ONE amortised likelihood (num_cw=1) gives the frozen marginalised noise; we only
    # use logL(theta_off, data) for arbitrary Earth-term data.
    pta1, cwb1, _ = build_enterprise_pta(ent_psrs, 1, components=N_COMPONENTS,
                                         log10_equad=LOG10_EQUAD)
    inj1 = make_geometry_injection(pta1, ent_psrs, 1, cwb1, seed=1000,
                                   h_range=(EQUAL_H, EQUAL_H))
    t0 = time.time()
    amo = build_fisher_amortised(disco_psrs, 1, inj1, cwb1)
    theta_keys = amo["theta_keys"]
    theta_off = np.array(amo["theta_truth"])
    theta_off[theta_keys.index("cw_log10_h")] = -30.0  # kill CW model amplitude
    theta_off = jnp.asarray(theta_off)
    zero_data = tuple(jnp.zeros(len(p.toas)) for p in disco_psrs)
    ll0 = float(amo["fisher_logL"](theta_off, zero_data))  # const (data=0)
    print(f"build+compile {time.time()-t0:.1f}s; logL(theta_off,0)={ll0:.3f}", flush=True)

    def snr_of_source(src):
        data = earth_term_data(disco_psrs, src)
        ll = float(amo["fisher_logL"](theta_off, data))
        s2 = 2.0 * (ll0 - ll)
        return float(np.sqrt(max(s2, 0.0)))

    # ---------- GATE: one loud source, my SNR vs white matched-filter ----------
    print("\n--- GATE: Earth-term (s|s) SNR vs white matched-filter, one loud source ---")
    rng = np.random.default_rng(7)
    gate_src = {"cos_gwtheta": 0.2, "gwphi": 1.5, "cos_inc": 0.3, "log10_mc": 9.5,
                "log10_fgw": -7.7, "log10_h": -13.0, "phase0": 1.0, "psi": 0.5}
    snr_full = snr_of_source(gate_src)
    # white matched filter: sqrt(sum (r_earth/sigma)^2)
    data = earth_term_data(disco_psrs, gate_src)
    s2_white = 0.0
    for p, r in zip(disco_psrs, data):
        s2_white += float(np.sum((np.array(r) / np.asarray(p.toaerrs)) ** 2))
    snr_white = np.sqrt(s2_white)
    ratio = snr_full / snr_white
    print(f"  loud source (log10_h=-13, log10_fgw=-7.7, mc=9.5):")
    print(f"  SNR_full (marginalised white+timing+RN+GWB) = {snr_full:.2f}")
    print(f"  SNR_white (matched filter sum (r/sigma)^2)   = {snr_white:.2f}")
    print(f"  ratio full/white = {ratio:.3f}  (GWB+timing degradation)  "
          f"GATE {'PASS' if 0.8 <= ratio <= 1.25 else 'see note: GWB/timing degradation >20%'}")

    # ---------- sweep ----------
    d4 = np.load(D4_NPZ, allow_pickle=True)
    rangeable = {}
    for n in NCW_LIST:
        ci = d4[f"row_{n}"][0]  # (med,min,max) class-i
        rangeable[n] = ci

    thresholds = [3, 5, 8]
    det = {n: {t: [] for t in thresholds} for n in NCW_LIST}
    print("\n--- DETECTABLE sweep (Earth-term SNR) ---", flush=True)
    for n in NCW_LIST:
        pta, cwb, _ = build_enterprise_pta(ent_psrs, n, components=N_COMPONENTS,
                                           log10_equad=LOG10_EQUAD)
        for d in range(10):
            inj = make_geometry_injection(pta, ent_psrs, n, cwb, seed=2000 + d,
                                          h_range=(EQUAL_H, EQUAL_H))
            snrs = np.array([snr_of_source(s) for s in src_list_from_inj(inj, cwb)])
            for t in thresholds:
                det[n][t].append(int((snrs > t).sum()))
        med5 = int(np.median(det[n][5]))
        print(f"  N_CW={n}: detectable(SNR>5) median {med5}/{n}  "
              f"(SNR>3 {int(np.median(det[n][3]))}, SNR>8 {int(np.median(det[n][8]))}); "
              f"rangeable(class-i) median {rangeable[n][0]}", flush=True)

    # population draw (3 loud + 13 faint, 10x) at N_CW=16
    pta, cwb, _ = build_enterprise_pta(ent_psrs, 16, components=N_COMPONENTS,
                                       log10_equad=LOG10_EQUAD)
    injp = make_geometry_injection(pta, ent_psrs, 16, cwb, seed=3000,
                                   population=(3, -13.25, -14.25))
    snrs_p = np.array([snr_of_source(s) for s in src_list_from_inj(injp, cwb)])
    pop_det = {t: int((snrs_p > t).sum()) for t in thresholds}
    pop_range = 0  # from D4: population class-i = 0
    print(f"\n  POPULATION (3 loud+13 faint) N_CW=16: detectable SNR>5 = {pop_det[5]} "
          f"(>3 {pop_det[3]}, >8 {pop_det[8]}); rangeable(class-i) = {pop_range}")
    print(f"    per-source SNRs: {np.array2string(np.sort(snrs_p)[::-1], precision=1, max_line_width=200)}")

    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_p2a_results.npz",
             ncw=np.array(NCW_LIST),
             det5=np.array([[np.median(det[n][5]), min(det[n][5]), max(det[n][5])] for n in NCW_LIST]),
             det3=np.array([[np.median(det[n][3]), min(det[n][3]), max(det[n][3])] for n in NCW_LIST]),
             det8=np.array([[np.median(det[n][8]), min(det[n][8]), max(det[n][8])] for n in NCW_LIST]),
             rangeable=np.array([list(rangeable[n]) for n in NCW_LIST]),
             pop_det=np.array([pop_det[3], pop_det[5], pop_det[8]]), pop_range=pop_range,
             gate_snr_full=snr_full, gate_snr_white=snr_white)

    # ---------- figure ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.array(NCW_LIST)
    d5 = np.array([[np.median(det[n][5]), min(det[n][5]), max(det[n][5])] for n in NCW_LIST])
    rg = np.array([list(rangeable[n]) for n in NCW_LIST])
    ax.plot(x, d5[:, 0], "o-", color="#2c7fb8", label="detectable (Earth-term SNR>5)")
    ax.fill_between(x, d5[:, 1], d5[:, 2], color="#2c7fb8", alpha=0.2)
    ax.plot(x, rg[:, 0], "s-", color="#d95f0e", label="rangeable (class-i)")
    ax.fill_between(x, rg[:, 1], rg[:, 2], color="#d95f0e", alpha=0.2)
    ax.plot([16], [pop_det[5]], "o", ms=11, mfc="none", color="#2c7fb8", label="detectable (population)")
    ax.plot([16], [pop_range], "s", ms=11, mfc="none", color="#d95f0e", label="rangeable (population)")
    ax.plot(x, x, ":", color="0.6", label="all sources")
    ax.set_xlabel("N_CW"); ax.set_ylabel("number of sources")
    ax.set_title("Detectability vs Rangeability (116 psr, equal strain, Asimov)")
    ax.legend(fontsize=8); ax.set_xscale("log", base=2); ax.set_xticks(x); ax.set_xticklabels(x)
    fig.tight_layout()
    fig.savefig("/home/mattm/projects/HSYMT/CW_transition/p2a_detect_vs_range.png", dpi=130)
    print("saved p2a_detect_vs_range.png, stagec_p2a_results.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
