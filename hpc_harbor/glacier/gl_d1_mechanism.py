#!/usr/bin/env python3
"""D1 -- MECHANISM CUT of the r13p25/e07 wrong-cert cascade (authorized 2026-07-25).

The SPARK selectivity question inverted: are the cascade's certified fringes RIGID
(true two-component Earth+pulsar structure) or manufactured composite power? Three
re-scorings of every scoring iteration of the two cascade cells, differing in ONE
ingredient each, so the 13 wrong certs decompose cleanly:

  V0  (theta_rec banked, real data)   -- reproduces the banked columns (control).
  V1  (theta TRUTH-anchored, real data) -- THE WANDER REMOVED: fed slots at drawn
      truth. If the false fringes vanish here, the M-step wander BUILT them
      (mismatch-manufactured); if they persist, the comb structure alone suffices.
  V2  (theta_rec banked, CLEAN truth data) -- THE NOISE REMOVED: zero-noise injection
      at the cell's truth. If the false fringes persist here, noise is irrelevant --
      the template-mismatch field is the whole mechanism.

Per (sky, iter): full per-pulsar columns banked for all three variants; the paper's
safety figure reads (dlnL, q, on_true) triplets per variant. Rigidity verdict per the
authorization: wrong certs that die under V1 and survive under V2 are
rigidity-DEFICIENT wander products -- SPARK's two-component gate (D2) is then the
in-record defense.

Budget: ~6 iters x 3 variants x 2 skies x ~160 s ~ 1.6 GPU-hr (<= the 8 authorized).
"""
import os, sys
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import glob
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
          "hpc_harbor/atlas", "hpc_harbor/spark", "hpc_harbor/glacier"):
    sys.path.insert(0, f"{HSYMT}/{p}")

import glacier_loop as GL
import glacier_pop as GP
from glacier_pop import bank_npz

RUNG, ARM = "r13p25", "e07"


def main():
    GL.check_affinity()
    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()
    out = {}
    for sky in (0, 1):
        pop, cond = GP.draw_population_conditional(sky, RUNG, n_src=256,
                                                   band_log10f=GL.BAND_CAMPAIGN)
        slots, n_harm, active, chan, n_clip = GL.embed_igniter(
            pop, GL.E_IGNITER[ARM], GL.VENUE_SPAN_S[30])
        pop_s = dict(pop); pop_s["src"] = slots; pop_s["n_src"] = len(slots)
        G = GP.build_glacier_problem(30, pop_s, verbose=True)
        G["slots"] = slots
        amo = G["amo"]; nd = amo["n_dist"]
        ones = jnp.ones(amo["npsr"])
        sb = GL.CertScoreboard(G, amo, jnp, C, prior_key=GL.TIER)
        noise_seed = GL.SEED_NOISE_BASE + 1000*sky
        L_true, n_true = sb.draw_truth(IG, noise_seed + 10_000_000, tier=GL.TIER)
        theta_true = np.asarray(amo["theta_truth"], float).copy()
        theta_true[sb.AI] = L_true
        clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
        ndraw = C.NoiseDrawer(sb.sp, amo, lgwb_path=GL.LGWB_BANKS[30], strict=True)
        noise = ndraw.draw(noise_seed, components=("white", "rn"))
        data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                     for c, n in zip(clean, noise))
        cleanT = tuple(jnp.asarray(np.asarray(c)) for c in clean)

        R = GP.OUT
        stem = f"gl2_{RUNG}_cell_{ARM}_s{sky}_T30_lit"
        prev_ci, prev_q = np.zeros(0, int), np.zeros(amo["npsr"])
        for it in range(6):
            fs = sorted(glob.glob(f"{R}/{stem}_i{it}__*.npz"))
            if not fs:
                break
            z = np.load(fs[0], allow_pickle=True)
            th_rec = np.asarray(z["theta_rec"], float)
            fed = np.where(np.asarray(z["fed_mask"]))[0]
            led = GP.PromoteLedger(slots)
            for k in fed:
                led.promote(int(k), slots[int(k)], iteration=it)
            # truth-anchored variant of the recovery state: fed slots at drawn truth,
            # distances at the banked recovery values (wander removal isolates SOURCE
            # params; theta_rec of this campaign never moves distances)
            th_tru = th_rec.copy()
            for k in fed:
                th_tru[nd + GP.NP_SRC*k : nd + GP.NP_SRC*(k+1)] = slots[int(k)]
            for tag, th, dat in (("V0", th_rec, data), ("V1", th_tru, data),
                                 ("V2", th_rec, cleanT)):
                dl, lnK, q, ot = sb.columns(th, led, dat, ones, prev_ci, prev_q)
                out[f"s{sky}_i{it}_{tag}_dlnL"] = dl
                out[f"s{sky}_i{it}_{tag}_q"] = q
                out[f"s{sky}_i{it}_{tag}_on_true"] = np.asarray(ot).astype(int)
                print(f"  [s{sky} i{it} {tag}] max dlnL {float(np.max(dl)):.2f} "
                      f"n(q>.9) {int((np.asarray(q) > 0.9).sum())}", flush=True)
            prev_ci = np.asarray(z["cert_idx"]).ravel().astype(int)
            prev_q = np.asarray(z["qmax"], float)
        out[f"s{sky}_fed_final"] = fed
    path = bank_npz("gl2_d1_mechanism",
                    note=("D1 mechanism cut (S4.18): three-variant re-scoring of the "
                          "cascade cells. V0 banked-state control; V1 wander removed "
                          "(fed at truth); V2 noise removed (clean data). Wrong certs "
                          "dying under V1 + surviving under V2 = wander-manufactured, "
                          "rigidity-deficient -> D2 tests SPARK's two-component gate."),
                    **out)
    print(f"banked -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
