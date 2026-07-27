#!/usr/bin/env python3
"""FRONTIER-v2 VALIDATION CASE (S4.20b, authorized 2026-07-25) -- gl1_null0_none_s1.

The reference artifact: the Stage-1 control cell whose frontier (v1: curvature ratio
alone) fed the igniter's WRONG-SKY template copy at iter 2 (safety STOP #2, S4.19).
Frontier-v2 (pre-registered S4.20a) adds the DATA-SUPPORT term: feed requires BOTH
ratio < 0.5 AND dlnL_feed > 0, where

  dlnL_feed(k) = logL(theta, k PRESENT at its template, other carried absent)
               - logL(theta, k ABSENT,               other carried absent)

(the own-term-live present-vs-absent contrast; the <d,dh>-facing term the curvature
self-term lacks). This script REPLAYS the cell's three frontier decisions from the
banked checkpoints (data rebuilt from seeds; nothing re-run in place -- the artifact
stays as banked) and scores every v1 candidate under v2.

PRE-REGISTERED EXPECTATION: the wrong-sky copy (slot 0) is REFUSED at the data-support
term (dlnL_feed <= 0); the cell's REAL feeds (the control's first-bite members) carry
positive support (v2 refuses none of them -- no false refusal).

Minutes-class: 3 Fisher sweeps + 2 likelihood evals per candidate.
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
from glacier_loop import I_COSGT, I_GWPHI

ARM, SKY, REAL = "none", 1, 0
STEM = f"gl1_null{REAL}_{ARM}_s{SKY}_T30_{GL.TIER}"


def main():
    GL.check_affinity()
    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()

    # ---- rebuild the cell exactly as run_cell (Stage-1 path, scrambled) ----
    pop = GL.draw_population(GL.SEED_POP_BASE + SKY, n_src=GL.N_POP_DEFAULT,
                             band_log10f=GL.BAND_CAMPAIGN)
    slots, n_harm, active, chan, n_clip = GL.embed_igniter(
        pop, GL.E_IGNITER[ARM], GL.VENUE_SPAN_S[GL.T_LABEL])
    pop_s = dict(pop); pop_s["src"] = slots; pop_s["n_src"] = len(slots)
    G = GP.build_glacier_problem(GL.T_LABEL, pop_s, verbose=True)
    G["slots"] = slots
    amo = G["amo"]; nd = amo["n_dist"]
    ones = jnp.ones(amo["npsr"])
    sb = GL.CertScoreboard(G, amo, jnp, C, prior_key=GL.TIER)
    noise_seed = GL.SEED_NOISE_BASE + 1000*SKY + REAL
    L_true, _ = sb.draw_truth(IG, noise_seed + 10_000_000, tier=GL.TIER)
    theta_true = np.asarray(amo["theta_truth"], float).copy()
    theta_true[sb.AI] = L_true
    clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
    ndraw = C.NoiseDrawer(sb.sp, amo, lgwb_path=GL.LGWB_BANKS[GL.T_LABEL], strict=True)
    noise = ndraw.draw(noise_seed, components=("white", "rn"))
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                 for c, n in zip(clean, noise))

    # the scramble (recovery template only, never the data) -- run_cell verbatim
    theta_rec0 = theta_true.copy()
    rng_s = np.random.default_rng(GL.SEED_SCRAM_BASE + 1000*SKY + REAL)
    theta_rec0[nd + I_COSGT] = rng_s.uniform(-1, 1)
    theta_rec0[nd + I_GWPHI] = rng_s.uniform(0, 2*np.pi)

    box_sigma = np.array([1.0, np.pi, 1.0, 0.5, 0.25, 0.25, np.pi, 0.5*np.pi]) \
        / np.sqrt(3.0)
    lb = amo["logL_batch_theta"]
    R = GP.OUT
    out = dict(arm=ARM, sky=SKY, real=REAL, cell_stem=STEM)

    # iter -> input state: iter 0 from the initial state, iter i>0 from bank i-1
    def state_in(it):
        if it == 0:
            return theta_rec0.copy(), np.zeros(len(slots), bool)
        fs = sorted(glob.glob(f"{R}/{STEM}_i{it-1}__*.npz"))
        if not fs:
            return None, None
        z = np.load(fs[0], allow_pickle=True)
        return (np.asarray(z["theta_rec"], float).copy(),
                np.asarray(z["fed_mask"], bool).copy())

    for it in range(3):                     # banked i0, i1 exist; the STOP was iter 2
        th, fed = state_in(it)
        if th is None:
            print(f"  iter {it}: no input bank -- stop", flush=True); break
        carried = np.where(~fed)[0]
        sig_opt, sig_pes, _ = GL.fisher_conditional(amo, jnp, th, carried, nd,
                                                    box_sigma, data, ones)
        ratio = np.max(sig_opt / sig_pes, axis=1)
        v1_new = np.where((ratio < GL.GATE_RATIO) & ~fed)[0]
        # score v2's data-support term for every v1 candidate, + slot 0 for the record
        cands = sorted(set(v1_new.tolist()) | ({0} if not fed[0] else set()))
        rows = []
        for k in cands:
            th_on = GL.theta_with_absent(th, nd, np.setdiff1d(carried, [k]))
            th_off = GL.theta_with_absent(th, nd, carried)
            ll = np.asarray(lb(jnp.asarray(np.stack([th_on, th_off])), data, ones))
            dsup = float(ll[0] - ll[1])
            v1 = bool(k in v1_new)
            v2 = bool(v1 and dsup > 0.0)
            rows.append((k, float(ratio[k]), dsup, v1, v2))
            tag = ("SCRAMBLED igniter" if k == 0 else "real member")
            print(f"  [iter {it}] slot {k:3d} ({tag}): ratio {ratio[k]:.4f} "
                  f"dlnL_feed {dsup:+.3f} v1_feed={v1} v2_feed={v2}", flush=True)
        rows = np.array(rows, float)
        out[f"i{it}_rows"] = rows            # cols: slot, ratio, dlnL_feed, v1, v2
        # advance nothing -- the replay reads each iter's input from the BANKED chain

    # verdict: the pre-registered expectation, read off iter 2's row for slot 0
    r2 = out.get("i2_rows")
    if r2 is not None and len(r2):
        s0 = r2[r2[:, 0] == 0]
        if len(s0):
            fed_v1, fed_v2 = bool(s0[0, 3]), bool(s0[0, 4])
            refused = fed_v1 and not fed_v2
            out["v1_fed_scrambled"] = fed_v1
            out["v2_refused_scrambled"] = refused
            # false-refusal count over ALL real-member v1 feeds in the replay
            fr = 0
            for it in range(3):
                rr = out.get(f"i{it}_rows")
                if rr is None: continue
                real_feeds = rr[(rr[:, 0] != 0) & (rr[:, 3] == 1)]
                fr += int((real_feeds[:, 4] == 0).sum())
            out["v2_false_refusals"] = fr
            print(f"\nVERDICT: v1 fed scrambled={fed_v1}; v2 refused={refused}; "
                  f"v2 false refusals of real feeds={fr}", flush=True)

    path = bank_npz("gl2_v2_frontier_validation",
                    note=("FRONTIER-v2 first validation case (S4.20b): replay of "
                          "gl1_null0_none_s1's frontier decisions from banked "
                          "checkpoints; v2 = ratio<0.5 AND dlnL_feed>0 (own-term-live "
                          "present-vs-absent data support). Pre-registered "
                          "expectation: the wrong-sky copy is refused at the "
                          "data-support term; no false refusal of real feeds."),
                    **out)
    print(f"banked -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
