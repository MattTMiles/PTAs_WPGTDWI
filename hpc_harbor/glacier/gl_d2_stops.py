#!/usr/bin/env python3
"""S4.23 -- THE STOP-MATERIAL KILL PASS (Matt 2026-07-26; gate frozen in S4.21).

Applies the IDENTICAL R1/R2 protocol (bar 15.132; lever-arm gain > 0 at the argmax fed
member; set-level; claim-as-made anchoring) to the STOPANAT banks of the 7 STOP cells:
6x r13p25/e07 scrambled nulls (STOP #3) + the capable-corner promote (STOP #4).
REPORT-ONLY / NON-PRE-REGISTERED (D2's pass/fail stays on the 13 signal-arm certs).

DECISION RULE (Matt, verbatim, stated in advance): all killed on R2 -> the STOPs
resolve as "manufacture occurs but is executable", loop protocol gains the kill step;
ANY survivor -> escalate immediately (purity breach inside the switch-on regime).

Promote-STOP cells stopped BEFORE the certification readout: the readout runs here at
the STOP state (fed + the promoted member); "cert material" is gated on the FLOOR-FREE
criterion parts (q > QBAR and dlnL > lnK + trials) since no floor was cut at a
promote-STOP; none expected -> banked as "no rigidity-scorable claim" (feed-layer
defect, frontier-v2's territory). Zero new draws anywhere: frozen seeds, banked states.
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

BAR_R1 = 15.132
# (cell stem, rung, arm, sky, real, T, wscale) -- grouped so builds are shared
GROUPS = {
    ("r13p25", "e07", 0, 30, 1.0): ["gl2_r13p25_null%d_e07_s0_T30_lit" % r
                                    for r in (0, 1, 2)],
    ("r13p25", "e07", 1, 30, 1.0): ["gl2_r13p25_null%d_e07_s1_T30_lit" % r
                                    for r in (0, 1, 2)],
    ("r13p9", "none", 1, 40, 0.25): ["gl2_r13p9_w0p25_null2_none_s1_T40_lit"],
}


def main():
    GL.check_affinity()
    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()
    R = GP.OUT
    out = dict(bar_r1=BAR_R1, report_only=True, pre_registered=False)
    n_sets = n_killed = 0
    survivors = []

    for (rung, arm, sky, T, wscale), stems in GROUPS.items():
        pop, cond = GP.draw_population_conditional(sky, rung, n_src=256,
                                                   band_log10f=GL.BAND_CAMPAIGN)
        slots, n_harm, active, chan, n_clip = GL.embed_igniter(
            pop, GL.E_IGNITER[arm], GL.VENUE_SPAN_S[T])
        pop_s = dict(pop); pop_s["src"] = slots; pop_s["n_src"] = len(slots)
        G = GP.build_glacier_problem(T, pop_s, verbose=True)
        G["slots"] = slots
        amo = G["amo"]; nd = amo["n_dist"]
        ones = jnp.ones(amo["npsr"])
        sb = GL.CertScoreboard(G, amo, jnp, C, prior_key=GL.TIER)
        lb = amo["logL_batch_theta"]; fl_lnl = amo["logL"]

        def lnl_at(th_base, kb, dat, pmask, free):
            th = th_base.at[kb + I_COSINC].set(free[0])
            th = th.at[kb + I_H].set(free[1])
            th = th.at[kb + I_PH0].set(free[2])
            th = th.at[kb + I_PSI].set(free[3])
            return fl_lnl(th, dat, pmask)
        prof = jax.jit(jax.vmap(_adam_profile(lnl_at),
                                in_axes=(None, 0, None, None)))
        CH = 4

        for stem in stems:
            fa = sorted(glob.glob(f"{R}/{stem}_STOPANAT_i*__*.npz"))
            if not fa:
                print(f"  {stem}: NO STOPANAT bank -- skipped", flush=True); continue
            za = np.load(fa[0], allow_pickle=True)
            it = int(za["iter"]); real = int(za["real"])
            kind = str(za["stop_kind"])
            th_rec = np.asarray(za["theta_rec"], float)
            fed_mask = np.asarray(za["fed_mask"], bool).copy()
            # data at THIS cell's frozen seeds (real enters the noise seed only)
            noise_seed = GL.SEED_NOISE_BASE + 1000*sky + real
            L_true, _ = sb.draw_truth(IG, noise_seed + 10_000_000, tier=GL.TIER)
            theta_true = np.asarray(amo["theta_truth"], float).copy()
            theta_true[sb.AI] = L_true
            clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
            ndraw = C.NoiseDrawer(sb.sp, amo, lgwb_path=GL.LGWB_BANKS[T], strict=True)
            noise = ndraw.draw(noise_seed, components=("white", "rn"),
                               white_scale=wscale)
            data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                         for c, n in zip(clean, noise))

            def twoF(fed_idx, pmask, Ld):
                th = GL.theta_with_absent(
                    th_rec, nd, np.setdiff1d(np.arange(len(slots)), fed_idx))
                th[sb.AI] = Ld
                kb_all = [nd + GP.NP_SRC*int(k) for k in fed_idx]
                thj, pmj = jnp.asarray(th), jnp.asarray(pmask)
                mx = []
                for c0 in range(0, len(kb_all), CH):
                    kb = kb_all[c0:c0+CH]
                    KB = jnp.asarray(kb + [kb[0]]*(CH - len(kb)))
                    m, _ = prof(thj, KB, data, pmj)
                    mx.append(np.asarray(m)[:len(kb)])
                mx = np.concatenate(mx)
                offs = []
                for k in fed_idx:
                    t = th.copy(); t[nd + GP.NP_SRC*int(k) + I_H] = H_ABSENT
                    offs.append(t)
                llo = np.asarray(lb(jnp.asarray(np.stack(offs)), data, pmj))
                return 2.0 * (mx - llo)

            if kind == "scrambled_promote":
                fed_mask[int(za["promoted_k"])] = True     # the claim IS the promote
            fed = np.where(fed_mask)[0]
            led = GP.PromoteLedger(slots)
            for k in fed:
                led.promote(int(k), slots[int(k)], iteration=it)
            # previous-iter cohering rows from the regular banked chain
            fp = sorted(glob.glob(f"{R}/{stem}_i{it-1}__*.npz")) if it > 0 else []
            if fp:
                zp = np.load(fp[0], allow_pickle=True)
                prev_ci = np.asarray(zp["cert_idx"], int)
                prev_q = np.asarray(zp["q_of_psr"], float)
            else:
                prev_ci, prev_q = np.zeros(0, int), np.zeros(amo["npsr"])
            dl, lnK, q, ot = sb.columns(th_rec, led, data, ones, prev_ci, prev_q)
            post = sb.FT.posterior(sb._last_LNL)
            map_L = np.asarray(post["map_L"], float)

            if kind == "wrong_cert":
                ci = np.asarray(za["cert_idx"], int)
            else:
                ci = np.where((np.asarray(q) > GL.QBAR) &
                              (np.asarray(dl) > np.asarray(lnK) + GL.TRIALS_NAT))[0]
                if not len(ci):
                    out[f"{stem}_verdict_txt"] = "no rigidity-scorable claim"
                    print(f"  {stem} [{kind} i{it}]: NO cert material at the STOP "
                          f"state -> feed-layer defect (frontier-v2 territory)",
                          flush=True)
                    continue
            pm = np.zeros(amo["npsr"]); pm[ci] = np.asarray(q)[ci]
            Ld = sb.L0.copy(); Ld[ci] = map_L[ci]
            tfc = twoF(fed, pm, Ld)
            tf0 = twoF(fed, np.zeros(amo["npsr"]), sb.L0)
            ks = int(np.argmax(tfc))
            c2F, d2F = float(tfc[ks]), float(tfc[ks] - tf0[ks])
            r1, r2 = c2F >= BAR_R1, d2F > 0.0
            killed = not (r1 and r2)
            n_sets += 1; n_killed += int(killed)
            if not killed:
                survivors.append(stem)
            out[f"{stem}_cert_idx"] = ci
            out[f"{stem}_twoF_coh"] = tfc
            out[f"{stem}_twoF_s0"] = tf0
            out[f"{stem}_verdict"] = np.array(
                [c2F, d2F, float(r1), float(r2), float(killed)])
            print(f"  {stem} [{kind} i{it}] set n={len(ci)}: 2F_coh {c2F:.2f} "
                  f"Delta2F {d2F:+.2f} -> {'KILLED' if killed else '** SURVIVOR **'}",
                  flush=True)

    out["n_sets"] = n_sets; out["n_killed"] = n_killed
    out["all_killed"] = bool(n_killed == n_sets)
    if survivors:
        print(f"\n*** ESCALATE TO MATT: rigidity SURVIVOR(S) in the STOP material: "
              f"{survivors} -- purity breach inside the switch-on regime; closure "
              f"language OFF per the S4.23 decision rule.", flush=True)
    else:
        print(f"\nSTOP-PASS: {n_killed}/{n_sets} scored sets killed; promote cells "
              f"without material banked as feed-layer defects. Per the S4.23 rule the "
              f"STOPs resolve as 'manufacture occurs but is executable'.", flush=True)
    path = bank_npz("gl2_d2_stops",
                    note=("S4.23 STOP-material kill pass -- REPORT-ONLY / "
                          "NON-PRE-REGISTERED (D2 pass/fail stays on the 13 "
                          "signal-arm certs). Frozen S4.21 gate over the STOPANAT "
                          "banks of the 6 e07 nulls + the capable-corner promote; "
                          "zero new draws (frozen seeds, banked states). Decision "
                          "rule per Matt 2026-07-26, recorded in S4.23."),
                    **out)
    print(f"banked -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
