#!/usr/bin/env python3
"""S4.24.1 -- (A) Finding-2 report-only scoring + (B) r13p5/none/s3 floor re-cut.

(A) The three r13p5/e07/s0 wrong certs (i0 psr62; i2 psr62+43) scored under the
FROZEN S4.21 R1/R2 gate as a SECOND report-only population -- direct scoring from the
existing banks, no replay. PREDICTION LOGGED IN THE CAPSTONE BEFORE THIS RAN: R2 kill
for all three; any survivor = immediate escalation, closure language off.

(B) The r13p5/none/s3 cell's floors re-cut at FULL nulls at its two cert iterations
(i0, i4): the banked floors were degenerate (0.0, emp_q95_degenerate on an all-zero
offender sample at 32 nulls). The provisional true certs (psr8, dlnL 2.35/1.60) are
re-evaluated against the re-cut floors; if the floor retires them, that is the result.
Zero new-draw semantics do not apply to (B): null draws for floors are the floor
machinery's OWN draws (the same convention every banked floor used).
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
R = GP.OUT


def build(sky, arm, rung="r13p5", T=30):
    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()
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
    noise_seed = GL.SEED_NOISE_BASE + 1000*sky
    L_true, _ = sb.draw_truth(IG, noise_seed + 10_000_000, tier=GL.TIER)
    theta_true = np.asarray(amo["theta_truth"], float).copy()
    theta_true[sb.AI] = L_true
    clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
    ndraw = C.NoiseDrawer(sb.sp, amo, lgwb_path=GL.LGWB_BANKS[T], strict=True)
    noise = ndraw.draw(noise_seed, components=("white", "rn"))
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                 for c, n in zip(clean, noise))
    return jax, jnp, FL, slots, G, amo, nd, ones, sb, ndraw, data


def banked_iter(stem, it):
    fs = [f for f in sorted(glob.glob(f"{R}/{stem}_i{it}__*.npz"))
          if "STOPANAT" not in f]
    return np.load(fs[0], allow_pickle=True) if fs else None


def main():
    GL.check_affinity()
    out = dict(bar_r1=BAR_R1, report_only=True, pre_registered=False,
               prediction="R2 kill for all three (logged S4.24.1 before scoring)")

    # ---------- (A) Finding-2: r13p5/e07/s0, sets at i0 (n=1) and i2 (n=2) ----------
    jax, jnp, FL, slots, G, amo, nd, ones, sb, ndraw, data = build(0, "e07")
    lb = amo["logL_batch_theta"]; fl_lnl = amo["logL"]

    def lnl_at(th_base, kb, dat, pmask, free):
        th = th_base.at[kb + I_COSINC].set(free[0])
        th = th.at[kb + I_H].set(free[1])
        th = th.at[kb + I_PH0].set(free[2])
        th = th.at[kb + I_PSI].set(free[3])
        return fl_lnl(th, dat, pmask)
    prof = jax.jit(jax.vmap(_adam_profile(lnl_at), in_axes=(None, 0, None, None)))
    CH = 4

    def twoF(th_rec, fed_idx, pmask, Ld):
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

    stem = "gl2_r13p5_cell_e07_s0_T30_lit"
    survivors = []
    for it in (0, 2):
        z = banked_iter(stem, it)
        ci = np.asarray(z["cert_idx"]).ravel().astype(int)
        th_rec = np.asarray(z["theta_rec"], float)
        fed = np.where(np.asarray(z["fed_mask"]))[0]
        qpsr = np.asarray(z["q_of_psr"], float)
        zp = banked_iter(stem, it-1) if it > 0 else None
        prev_ci = np.asarray(zp["cert_idx"], int) if zp is not None else np.zeros(0, int)
        prev_q = (np.asarray(zp["q_of_psr"], float) if zp is not None
                  else np.zeros(amo["npsr"]))
        led = GP.PromoteLedger(slots)
        for k in fed:
            led.promote(int(k), slots[int(k)], iteration=it)
        sb.columns(th_rec, led, data, ones, prev_ci, prev_q)
        map_L = np.asarray(sb.FT.posterior(sb._last_LNL)["map_L"], float)
        pm = np.zeros(amo["npsr"]); pm[ci] = qpsr[ci]
        Ld = sb.L0.copy(); Ld[ci] = map_L[ci]
        tfc = twoF(th_rec, fed, pm, Ld)
        tf0 = twoF(th_rec, fed, np.zeros(amo["npsr"]), sb.L0)
        ks = int(np.argmax(tfc))
        c2F, d2F = float(tfc[ks]), float(tfc[ks] - tf0[ks])
        killed = not (c2F >= BAR_R1 and d2F > 0.0)
        if not killed:
            survivors.append(f"{stem}_i{it}")
        out[f"ext_i{it}_cert_idx"] = ci
        out[f"ext_i{it}_verdict"] = np.array([c2F, d2F, float(killed)])
        print(f"  [ext r13p5/e07/s0 i{it}] set n={len(ci)}: 2F_coh {c2F:.2f} "
              f"Delta2F {d2F:+.2f} -> {'KILLED' if killed else '** SURVIVOR **'}",
              flush=True)
    out["ext_survivors"] = np.array(survivors)
    if survivors:
        print("\n*** ESCALATE TO MATT: survivor(s) in the r13p5 extension "
              "population -- closure language OFF per the standing rule.", flush=True)

    # ---------- (B) floor re-cut: r13p5/none/s3 at full nulls, iters 0 and 4 --------
    jax, jnp, FL, slots, G, amo, nd, ones, sb, ndraw, data = build(3, "none")
    stem = "gl2_r13p5_cell_none_s3_T30_lit"
    for it in (0, 4):
        z = banked_iter(stem, it)
        th_rec = np.asarray(z["theta_rec"], float)
        fed = np.where(np.asarray(z["fed_mask"]))[0]
        zp = banked_iter(stem, it-1) if it > 0 else None
        prev_ci = np.asarray(zp["cert_idx"], int) if zp is not None else np.zeros(0, int)
        prev_q = (np.asarray(zp["q_of_psr"], float) if zp is not None
                  else np.zeros(amo["npsr"]))
        led = GP.PromoteLedger(slots)
        for k in fed:
            led.promote(int(k), slots[int(k)], iteration=it)
        off = GL._null_offenders(sb, ndraw, jnp, th_rec, led, 3, it,
                                 prev_ci, prev_q, n_null=GL.N_NULL_FULL, wscale=1.0)
        zf = float((off == 0.0).mean())
        if float(np.ptp(off)) == 0.0:
            fl, err, est = FL["emp_quantile"](off), 0.0, "emp_q95_degenerate"
        else:
            gu, mu, beta, gsd, nn = FL["gumbel_floor"](off)
            fl, err, est = FL["adopt"](gu, gsd, off, zf)
        dl = float(np.asarray(z["dlnL_det"])[8])
        lnK8 = float(np.asarray(z["lnK"])[8])
        q8 = float(np.asarray(z["qmax"])[8])
        cert = (dl > max(lnK8 + GL.TRIALS_NAT, fl)) and (q8 > GL.QBAR)
        out[f"recut_i{it}"] = np.array([fl, err, zf, dl, float(cert)])
        out[f"recut_i{it}_est"] = est
        print(f"  [recut r13p5/none/s3 i{it}] full-null floor {fl:.3f}+-{err:.3f} "
              f"({est}, zf {zf:.3f}, n={GL.N_NULL_FULL}); psr8 dlnL {dl:.2f} -> "
              f"cert {'SURVIVES' if cert else 'RETIRED'}", flush=True)

    path = bank_npz("gl2_d2_ext",
                    note=("S4.24.1: (A) Finding-2 second report-only population "
                          "(r13p5/e07/s0 wrong certs) under the frozen S4.21 gate, "
                          "prediction logged before scoring; (B) r13p5/none/s3 floor "
                          "re-cut at full nulls, provisional true certs re-evaluated. "
                          "Report-only / non-pre-registered; 13-cert pass/fail "
                          "unchanged."),
                    **out)
    print(f"banked -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
