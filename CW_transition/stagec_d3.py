"""Stage C — Deliverable D3: the fringe split (the honest headline).

D2's marginal sigma_L is WITHIN-FRINGE precision only. A within-fringe sigma below
the EM prior is NOT a distance measurement unless the correct distance FRINGE is also
identified. This deliverable does that classification, per pulsar, on the SAME 10
zero-noise (Asimov) geometry draws as D2.

Per pulsar a (Asimov = signal-only, zero noise -> expected dlnL; optimistic upper bound,
like D2's zero-noise Fisher):
  dL_a   : fringe spacing (compute_mode_spacing; min over CWs)
  K_a    : number of candidate fringes in the +/-3 sigma_EM window = 2*3*sigma_EM/dL_a
  dlnL_a : joint dlnL between true fringe and best wrong fringe in the window
           (find_best_wrong_mode_in_prior; others + CW at truth)
  identifiable_a := dlnL_a > ln(K_a)              (beats the trials factor)
  identifiable_strict := dlnL_a > ln(K_a) + 3      (robustness margin)

Classification (uses D2 marginal sigma_L, same seeds):
  (i)   identifiable AND marginal sigma_L < EM prior  [distance genuinely improved]
  (ii)  identifiable AND marginal sigma_L >= EM prior
  (iii) NOT identifiable                              [sigma_L meaningless]

Run in jug-gpu:  python stagec_d3.py --ndraws 10
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
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import (
    load_pulsars, build_enterprise_pta, build_fast_scan_likelihood,
    inject_noisefree_cw, compute_mode_spacing,
    _enterprise_to_disco_params, MultiSourceDelay,
)

_FRINGE_EVAL_CACHE = {}


N_FRINGE_EVAL = 512  # FIXED eval budget per pulsar (fixed shape -> compile once)
_FRINGE_CHUNK = 256


def fringe_deltalnL(logl_fn, base_values, param_keys, dist_key, L0, sig_em, dL,
                    n_sigma=3.0):
    """Joint dlnL between the true fringe and the best WRONG fringe in +/-n_sigma_EM.

    Evaluates lnL at the fringe CENTERS L0 + k*dL (each fringe's ~peak, since the
    likelihood is periodic in the pulsar-term phase = 2*pi*L/dL) rather than a
    find_peaks scan, whose prominence floor fails when the whole fringe modulation
    is < the floor (the near-degenerate single-CW case). A FIXED N_FRINGE_EVAL
    centers are evaluated (dense near truth + uniform across the window) so the
    batch shape is constant and XLA compiles once. Returns (dlnL, K) with K = the
    TRUE number of wrong fringes in the window (honest trials factor).
    """
    lo = max(1e-3, L0 - n_sigma * sig_em)
    hi = L0 + n_sigma * sig_em
    kmin = int(np.floor((lo - L0) / dL))
    kmax = int(np.ceil((hi - L0) / dL))
    ks = np.arange(kmin, kmax + 1)
    centers = L0 + ks * dL
    keep = (centers >= lo) & (centers <= hi) & (ks != 0)
    wrong, wrong_ks = centers[keep], ks[keep]
    K = len(wrong)
    if K == 0:
        return np.inf, 0

    # select a FIXED-size set: dense nearest fringes + uniform across the rest,
    # then pad by repetition to exactly N_FRINGE_EVAL (dups don't change the max)
    if K > N_FRINGE_EVAL:
        near = np.abs(wrong_ks) <= N_FRINGE_EVAL // 2
        far_idx = np.where(~near)[0]
        n_far = N_FRINGE_EVAL - int(near.sum())
        far_idx = far_idx[np.linspace(0, len(far_idx) - 1, max(n_far, 0)).astype(int)]
        sel = np.concatenate([np.where(near)[0], far_idx])
        evalL = wrong[sel]
    else:
        evalL = np.concatenate([wrong, np.full(N_FRINGE_EVAL - K, wrong[0])])
    evalL = evalL[:N_FRINGE_EVAL]

    idx_key = param_keys.index(dist_key)
    cache = _FRINGE_EVAL_CACHE.get(id(logl_fn))
    if cache is None:
        cache = jax.jit(jax.vmap(logl_fn))
        _FRINGE_EVAL_CACHE[id(logl_fn)] = cache

    ll_truth = float(logl_fn(jnp.asarray(base_values)))
    base_j = jnp.asarray(base_values)
    best = -np.inf
    for s in range(0, N_FRINGE_EVAL, _FRINGE_CHUNK):
        blk = evalL[s:s + _FRINGE_CHUNK]  # always length _FRINGE_CHUNK
        t = jnp.repeat(base_j[None, :], _FRINGE_CHUNK, axis=0)
        x = t.at[:, idx_key].set(jnp.asarray(blk))
        best = max(best, float(np.array(cache(x)).max()))
    return ll_truth - best, K
from discovery.deterministic import make_phase_connected_binary
from stagec_fisher import (
    make_geometry_injection, N_COMPONENTS, LOG10_EQUAD, CW_PARAM_BASE, KPC_TO_PC,
)

D2_NPZ = "/home/mattm/projects/HSYMT/CW_transition/stagec_d2_results.npz"


def cw_list_from_inj(inj, cwb):
    suff = ["" if i == 0 else f"_{i+1}" for i in range(len(cwb))]
    out = []
    for name, s in zip(cwb, suff):
        out.append({k: inj[f"{name}_{k}"] for k in CW_PARAM_BASE})
    return out


def pulsar_term_snr(disco_psrs, inj, cwb):
    """Per-pulsar pulsar-term SNR: sqrt( sum (r_pterm/sigma_toa)^2 ) at truth.

    r_pterm = (full CW delay incl pulsar term) - (Earth-term-only delay).
    Uses the heterogeneous white-noise sigma (toaerrs). Order-unity check only.
    """
    cw_list = cw_list_from_inj(inj, cwb)
    keys = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
            "log10_fgw", "log10_h", "phase0", "psi")
    stacked = {k: jnp.array([cw[k] for cw in cw_list]) for k in keys}
    full = jax.vmap(make_phase_connected_binary(pulsarterm=True),
                    in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0, None))
    # Earth-only term is distance-independent: p_dist is keyword-only, 10 positional args.
    earth = jax.vmap(make_phase_connected_binary(pulsarterm=False),
                     in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0))
    snr = {}
    for psr in disco_psrs:
        base_args = (psr.toas, psr.pos, stacked["cos_gwtheta"], stacked["gwphi"],
                     stacked["cos_inc"], stacked["log10_mc"], stacked["log10_fgw"],
                     stacked["log10_h"], stacked["phase0"], stacked["psi"])
        r_full = np.array(jnp.sum(full(*base_args, psr.pdist[0]), axis=0))
        r_earth = np.array(jnp.sum(earth(*base_args), axis=0))
        r_pterm = r_full - r_earth
        sig = np.asarray(psr.toaerrs)
        snr[psr.name] = float(np.sqrt(np.sum((r_pterm / sig) ** 2)))
    return snr


def run(ndraws):
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)
    print(f"=== D3: fringe split, 116 psr, 1 CW, {ndraws} Asimov draws ===", flush=True)

    ent_psrs, disco_psrs = load_pulsars(116)
    pta, cwb, _ = build_enterprise_pta(ent_psrs, 1, components=N_COMPONENTS,
                                       log10_equad=LOG10_EQUAD)
    names = [p.name for p in ent_psrs]
    em_pc = np.array([p.pdist[1] for p in ent_psrs]) * KPC_TO_PC
    sig_em_kpc = np.array([p.pdist[1] for p in ent_psrs])

    d2 = np.load(D2_NPZ, allow_pickle=True)
    d2_names = list(d2["psr_names"])
    assert d2_names == names, "D2/D3 pulsar order mismatch"
    marg_draws_pc = d2["marg_draws"]  # (ndraws, 116) pc, seeds 2000+d

    npsr = len(names)
    # per (draw, pulsar): K, dlnL, dL, ident, ident_strict, class
    K = np.zeros((ndraws, npsr)); dlnL = np.zeros((ndraws, npsr))
    dL = np.zeros((ndraws, npsr))
    ident = np.zeros((ndraws, npsr), bool); ident_s = np.zeros((ndraws, npsr), bool)
    cls = np.zeros((ndraws, npsr), int)
    snr_draws = np.zeros((ndraws, npsr))

    for d in range(ndraws):
        t0 = time.time()
        inj = make_geometry_injection(pta, ent_psrs, 1, cwb, seed=2000 + d)
        cw_list = cw_list_from_inj(inj, cwb)
        resid_map, _ = inject_noisefree_cw(disco_psrs, cw_list)  # zero-noise Asimov
        logl_fn, pk, bv = build_fast_scan_likelihood(
            disco_psrs, resid_map, 1, inj, cwb,
            components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
        snr = pulsar_term_snr(disco_psrs, inj, cwb)
        _FRINGE_EVAL_CACHE.clear()

        for a, (pe, pd) in enumerate(zip(ent_psrs, disco_psrs)):
            dLa = min(compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"],
                                           cw["log10_fgw"], pd.pos) for cw in cw_list)
            dL[d, a] = dLa
            dd, Ka = fringe_deltalnL(logl_fn, bv, pk, f"{pe.name}_cw_p_dist",
                                     pe.pdist[0], pe.pdist[1], dLa, n_sigma=3.0)
            Ka = max(1, Ka)
            K[d, a] = Ka
            dlnL[d, a] = dd
            snr_draws[d, a] = snr[pe.name]
            idn = dd > np.log(Ka)
            ident[d, a] = idn
            ident_s[d, a] = dd > np.log(Ka) + 3.0
            marg_below = marg_draws_pc[d, a] < em_pc[a]
            if not idn:
                cls[d, a] = 3
            elif marg_below:
                cls[d, a] = 1
            else:
                cls[d, a] = 2
        print(f"  draw {d}: {time.time()-t0:.1f}s  "
              f"class(i)={int((cls[d]==1).sum())} (ii)={int((cls[d]==2).sum())} "
              f"(iii)={int((cls[d]==3).sum())}  ident={int(ident[d].sum())} "
              f"ident_strict={int(ident_s[d].sum())}", flush=True)

    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_d3_results.npz",
             psr_names=np.array(names), em_pc=em_pc, K=K, dlnL=dlnL, dL=dL,
             ident=ident, ident_strict=ident_s, cls=cls, snr=snr_draws,
             marg_draws_pc=marg_draws_pc)

    # -------- report --------
    print("\n================ D3 SUMMARY ================")
    for c, lab in [(1, "(i)   ident & marg<EM  [genuinely improved]"),
                   (2, "(ii)  ident & marg>=EM"),
                   (3, "(iii) NOT identifiable [sigma_L meaningless]")]:
        cnt = (cls == c).sum(axis=1)
        print(f"class {lab}: median {int(np.median(cnt))}/{npsr} "
              f"(range {cnt.min()}-{cnt.max()} across draws)")
    idc = ident.sum(axis=1); idsc = ident_s.sum(axis=1)
    print(f"identifiable (dlnL>lnK):      median {int(np.median(idc))}/{npsr} "
          f"(range {idc.min()}-{idc.max()})")
    print(f"identifiable (dlnL>lnK+3):    median {int(np.median(idsc))}/{npsr} "
          f"(range {idsc.min()}-{idsc.max()})")

    # class-(i) pulsar list: those in class (i) in a MAJORITY of draws
    frac_i = (cls == 1).mean(axis=0)
    order = np.argsort(-frac_i)
    print(f"\nclass-(i) pulsars (in class i for >= half of draws): "
          f"{int((frac_i >= 0.5).sum())}")
    print(f"{'pulsar':13s} {'frac_i':>6s} {'K_a':>8s} {'dlnL_a':>10s} {'lnK_a':>7s} "
          f"{'dL_pc':>8s} {'marg_pc':>8s} {'EM_pc':>8s}")
    Kmed = np.median(K, axis=0); dmed = np.median(dlnL, axis=0)
    dLmed_pc = np.median(dL, axis=0) * KPC_TO_PC
    margmed = np.median(marg_draws_pc, axis=0)
    for a in order:
        if frac_i[a] < 0.5:
            break
        dstr = "inf" if not np.isfinite(dmed[a]) else f"{dmed[a]:10.2f}"
        print(f"{names[a]:13s} {frac_i[a]:6.1f} {Kmed[a]:8.1f} {dstr:>10s} "
              f"{np.log(Kmed[a]):7.2f} {dLmed_pc[a]:8.4g} {margmed[a]:8.4g} {em_pc[a]:8.4g}")

    # J0437 bookkeeping
    if "J0437-4715" in names:
        a = names.index("J0437-4715")
        print(f"\nJ0437-4715: IN array. K={Kmed[a]:.1f} dlnL={dmed[a]:.2f} "
              f"lnK={np.log(Kmed[a]):.2f} ident={dmed[a]>np.log(Kmed[a])} "
              f"marg={margmed[a]:.4g} pc EM={em_pc[a]:.4g} pc class(median)="
              f"{int(np.median(cls[:,a]))}")

    # SNR cross-check on best class-(i) candidate
    cand = [a for a in order if frac_i[a] >= 0.5]
    if cand:
        a = cand[0]
        snr_a = np.median(snr_draws[:, a]); dL_a_kpc = np.median(dL[:, a])
        cond_pc = np.median(d2["cond_draws"][:, a])
        pred_pc = dL_a_kpc / (2 * np.pi * snr_a) * KPC_TO_PC
        print(f"\nSNR cross-check on best class-(i) pulsar {names[a]}:")
        print(f"  SNR_pterm={snr_a:.3g}  dL={dL_a_kpc*KPC_TO_PC:.4g} pc  "
              f"K={Kmed[a]:.1f} = 2*3*sigma_EM/dL = "
              f"2*3*{sig_em_kpc[a]*KPC_TO_PC:.3g}/{dL_a_kpc*KPC_TO_PC:.4g}")
        print(f"  dL/(2pi*SNR) = {pred_pc:.4g} pc   vs  conditional sigma_L(Fisher) = "
              f"{cond_pc:.4g} pc   ratio {pred_pc/cond_pc:.3f}")
    print("saved stagec_d3_results.npz")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ndraws", type=int, default=10)
    sys.exit(run(ap.parse_args().ndraws))
