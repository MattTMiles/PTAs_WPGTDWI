#!/usr/bin/env python3
"""D2 -- THE RIGIDITY GATE, KILL CHECK (S4.21, gate pre-stated BEFORE this ran).

THE GATE (v2.3-candidate; thresholds from SPARK_results/spark_s2c.npz, NOTHING tuned
on these cells):
  R1 (coherent detection): 2F_coh >= 15.132  (= floor_sC_m1e05_gumbel - 2*sd,
      the LOWEST banked coherent-state floor minus 2 sigma -- kill against the
      lowest bar is the robust kill)
  R2 (lever-arm rigidity):  Delta2F = 2F_coh(k*) - 2F_s0(k*) > 0  at the argmax
      member k* (coherence at the certified lever arms must ADD matched power)
  RETAIN iff R1 AND R2; applied per (cell, scoring iteration) to that iteration's
  certified SET; 2F_coh = MAX over fed members (most generous attribution).

STATISTIC (SPARK's oracle F-stat, spark.py::mode_s2c verbatim conventions):
  2F_k(pmask, Ld) = 2*[ max over (cos_inc, log10_h, phase0, psi) of
                        logL(theta, k free, data, pmask)
                        - logL(theta, k REMOVED, data, pmask) ]
anchored on THE CLAIM AS MADE: the banked (wandered) theta_rec of the cascade cells,
carried census absent, coherent state = (pmask = banked qmax on the certified set,
Ld = FT.posterior map_L recomputed at the banked state); s0 = (pmask = 0, Ld = L0).

KILL CHECK passes iff all 13 cascade wrong certs (r13p25/e07 skies 0-1) are killed.
The FALSE-NEGATIVE check (banked IGNITE-2/CHORUS/A-SKY true certs) is the companion
run -- BOTH must pass for v2.3-candidate status; either failing is reported as-is.
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
from spark import _adam_profile, I_COSINC, I_H, I_PH0, I_PSI, H_ABSENT

RUNG, ARM = "r13p25", "e07"
BAR_R1 = 15.132          # floor_sC_m1e05_gumbel - 2*sd = 15.748 - 0.617 (S4.21)


def main():
    GL.check_affinity()
    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()
    out = dict(bar_r1=BAR_R1)
    n_sets = n_killed = n_certs = n_certs_killed = 0

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
        L_true, _ = sb.draw_truth(IG, noise_seed + 10_000_000, tier=GL.TIER)
        theta_true = np.asarray(amo["theta_truth"], float).copy()
        theta_true[sb.AI] = L_true
        clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
        ndraw = C.NoiseDrawer(sb.sp, amo, lgwb_path=GL.LGWB_BANKS[30], strict=True)
        noise = ndraw.draw(noise_seed, components=("white", "rn"))
        data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                     for c, n in zip(clean, noise))
        lb = amo["logL_batch_theta"]
        fl_lnl = amo["logL"]

        # the oracle profile on THIS problem (spark.py::make_profile16 pattern; the
        # member enters as a traced base offset so fed members vmap in one call)
        def lnl_at(th_base, kb, dat, pmask, free):
            th = th_base.at[kb + I_COSINC].set(free[0])
            th = th.at[kb + I_H].set(free[1])
            th = th.at[kb + I_PH0].set(free[2])
            th = th.at[kb + I_PSI].set(free[3])
            return fl_lnl(th, dat, pmask)
        prof = jax.jit(jax.vmap(_adam_profile(lnl_at),
                                in_axes=(None, 0, None, None)))

        CH = 4    # profile chunk: 20-wide vmap OOMed at 53 GiB (job 12778961); 4-wide fits

        def twoF(th_state, fed_idx, pmask, Ld):
            """2F per fed member at (pmask, Ld); th_state already carries the
            certified claim (wandered sources); distances overwritten by Ld.
            Chunked CH-wide (padded, one compile) -- the full-fed vmap OOMs."""
            th = np.asarray(th_state, float).copy()
            th[sb.AI] = Ld
            kb_all = [nd + GP.NP_SRC*int(k) for k in fed_idx]
            thj, pmj = jnp.asarray(th), jnp.asarray(pmask)
            mx = []
            for c0 in range(0, len(kb_all), CH):
                kb = kb_all[c0:c0+CH]
                pad = CH - len(kb)
                KB = jnp.asarray(kb + [kb[0]]*pad)
                m, _ = prof(thj, KB, data, pmj)
                mx.append(np.asarray(m)[:len(kb)])
            mx = np.concatenate(mx)
            ll_off = []
            for c0 in range(0, len(fed_idx), CH):
                off = []
                for k in fed_idx[c0:c0+CH]:
                    t = th.copy(); t[nd + GP.NP_SRC*int(k) + I_H] = H_ABSENT
                    off.append(t)
                ll_off.append(np.asarray(lb(jnp.asarray(np.stack(off)), data, pmj)))
            ll_off = np.concatenate(ll_off)
            return 2.0 * (mx - ll_off)

        R = GP.OUT
        stem = f"gl2_{RUNG}_cell_{ARM}_s{sky}_T30_{GL.TIER}"
        prev_ci, prev_q = np.zeros(0, int), np.zeros(amo["npsr"])
        for it in range(6):
            fs = sorted(glob.glob(f"{R}/{stem}_i{it}__*.npz"))
            if not fs:
                break
            z = np.load(fs[0], allow_pickle=True)
            ci = np.asarray(z["cert_idx"]).ravel().astype(int)
            th_rec = np.asarray(z["theta_rec"], float)
            fed = np.where(np.asarray(z["fed_mask"]))[0]
            qpsr = np.asarray(z["q_of_psr"], float)
            if len(ci):
                ck = f"{R}/gl2_d2_kill_ck_s{sky}_i{it}.npz"
                if os.path.exists(ck):
                    zc = np.load(ck, allow_pickle=True)
                    for kk in zc.files:
                        out[kk] = zc[kk]
                    v = np.asarray(zc[f"s{sky}_i{it}_verdict"], float)
                    killed = bool(v[4])
                    n_sets += 1; n_killed += int(killed)
                    n_certs += len(ci); n_certs_killed += len(ci)*int(killed)
                    print(f"  [s{sky} i{it}] checkpoint exists -- resumed "
                          f"({'KILLED' if killed else 'RETAINED'})", flush=True)
                    prev_ci, prev_q = ci, qpsr
                    continue
                led = GP.PromoteLedger(slots)
                for k in fed:
                    led.promote(int(k), slots[int(k)], iteration=it)
                # E-step at the banked state -> the claimed (MAP) distances
                sb.columns(th_rec, led, data, ones, prev_ci, prev_q)
                post = sb.FT.posterior(sb._last_LNL)
                map_L = np.asarray(post["map_L"], float)
                carried = np.setdiff1d(np.arange(len(slots)), fed)
                th_state = GL.theta_with_absent(th_rec, nd, carried)
                # coherent state: THE CLAIM AS MADE
                pm_coh = np.zeros(amo["npsr"]); pm_coh[ci] = qpsr[ci]
                Ld_coh = sb.L0.copy(); Ld_coh[ci] = map_L[ci]
                tf_coh = twoF(th_state, fed, pm_coh, Ld_coh)
                # s0 state: earth-only reference
                tf_s0 = twoF(th_state, fed, np.zeros(amo["npsr"]), sb.L0)
                ks = int(np.argmax(tf_coh))
                c2F, d2F = float(tf_coh[ks]), float(tf_coh[ks] - tf_s0[ks])
                r1, r2 = (c2F >= BAR_R1), (d2F > 0.0)
                killed = not (r1 and r2)
                n_sets += 1; n_killed += int(killed)
                n_certs += len(ci); n_certs_killed += len(ci)*int(killed)
                out[f"s{sky}_i{it}_cert_idx"] = ci
                out[f"s{sky}_i{it}_q"] = qpsr[ci]
                out[f"s{sky}_i{it}_mapL"] = map_L[ci]
                out[f"s{sky}_i{it}_twoF_coh"] = tf_coh
                out[f"s{sky}_i{it}_twoF_s0"] = tf_s0
                out[f"s{sky}_i{it}_fed"] = fed
                out[f"s{sky}_i{it}_verdict"] = np.array(
                    [c2F, d2F, float(r1), float(r2), float(killed)])
                print(f"  [s{sky} i{it}] set n={len(ci)}: 2F_coh {c2F:.2f} "
                      f"(bar {BAR_R1}) Delta2F {d2F:+.2f} -> "
                      f"{'KILLED' if killed else '** RETAINED **'}", flush=True)
                np.savez(ck, **{kk: out[kk] for kk in out
                                if kk.startswith(f"s{sky}_i{it}_")})
            prev_ci, prev_q = ci, qpsr

    out["n_sets"] = n_sets; out["n_sets_killed"] = n_killed
    out["n_certs"] = n_certs; out["n_certs_killed"] = n_certs_killed
    out["kill_check_pass"] = bool(n_certs_killed == n_certs and n_certs > 0)
    print(f"\nKILL CHECK: {n_certs_killed}/{n_certs} certs killed "
          f"({n_killed}/{n_sets} sets) -> "
          f"{'PASS' if out['kill_check_pass'] else 'FAIL'}", flush=True)
    path = bank_npz("gl2_d2_rigidity_kill",
                    note=("D2 KILL CHECK (S4.21): SPARK oracle-2F rigidity gate "
                          "(R1: 2F_coh>=15.132 = lowest banked coherent floor - 2sd; "
                          "R2: lever-arm gain > 0 at the argmax member) applied to "
                          "the cascade cells' certified sets at their banked states. "
                          "Thresholds from spark_s2c.npz, stated before scoring, "
                          "never tuned on these cells."),
                    **out)
    print(f"banked -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
