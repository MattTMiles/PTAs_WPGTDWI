#!/usr/bin/env python3
"""SIEVE V2 -- FALSE-ALARM RATE OF THE T2 TEMPLATE-HEALTH MONITOR ON TRUE-ANCHORED
MEMBERS (agent SIEVE-A, 2026-07-29 addendum).

WHY THIS EXISTS.  T2 killed the even/odd cross-validation split and promoted a cheaper
quantity in its place: the PRESENT-vs-ABSENT data-support contrast on an already-fed
member, which is frontier-v2's own support term re-asked. T2 measured it on THREE
GLACIER cells holding TWO true certifications. A monitor that fires on a healthy
template is worse than no monitor, and n = 2 cannot bound a false-alarm rate. This
script measures the false-alarm rate where the loop WORKS -- on truth-anchored members
of the GENERALISE ARM A-SKY survivor cell, all 8 skies.

WHY A-SKY IS THE RIGHT SUBSTRATE, AND CHORUS/IGNITE-2 SLOOPS ARE NOT.
  generalise.run_gen_realisation builds its search template as
      theta_base = theta_src.copy();  theta_base[AI] = L0
  i.e. EVERY SOURCE PARAMETER IS AT TRUTH and only the distances sit at the fiducial
  L0 -- the pure fringe problem with no source-parameter wander at all. That is exactly
  "true-anchored", by construction rather than by inference.
  The CHORUS `ch_sloop_*` and IGNITE-2 `ig_sloop_*` banks are SOFT-LOOP TRAJECTORIES:
  the template is FITTED, so a member there is anchored only to the extent the fit
  converged. Scoring them would mix genuine wander into the false-alarm count and
  inflate it. They carry the true certifications (V1's substrate) but they cannot
  bound a FALSE-ALARM rate. Declared, not assumed.

WHICH CELL.  gen_armAS_sky.npz verdict column: cell 1 = 'e0.3 h-12.75 5+11' is the ONLY
SURVIVES (n_conf = 6 of 8 skies); cells 0/2/3 FAIL. Default --cell 1.

THE TWO CONTRASTS (both computed; they answer different questions)

  INCREMENTAL (leave-one-out).  The direct analogue of the GLACIER loop, where the
  other members are already in the template:
      dlnL_incr[i] = logL(all 16 present) - logL(all present except i)
  This is what frontier-v2 asks about a candidate against a populated carried set.

  ISOLATED (member alone).  The analogue of the FIRST fed member, when nothing else is
  in the template yet:
      dlnL_iso[i] = logL(only i present) - logL(nothing present)
  T2's promoted quantity sits between these two; reporting both brackets it instead of
  picking one and hoping.

  NB THE T2 BUG THIS INHERITS THE FIX FROM.  In `run_cell`'s frontier-v2 the "present"
  leg may be written setdiff1d(carried, [k]) because there k is a CANDIDATE and is
  still IN carried. For a slot that is ALREADY FED, k is NOT in carried and that form
  collapses to th_on == th_off and returns identically 0.000 -- which in T2 printed a
  clean, wrong, pre-registration-confirming result. Both contrasts here are built by
  explicit index sets and G-V2d asserts the two templates differ.

ABSENCE IS STATED AT H_ABSENT = -30.0, not spark3's -18.0. BASELINE measured -18 to be
NOT absent at T = 30 (the residual member still moves the likelihood). G-V2c re-measures
the plateau on this venue rather than trusting that across campaigns.

FALSE ALARM := dlnL <= 0 on a member that IS in the data. Reported stratified by
injected loudness, because a member 1.5 dex below the loud set carries no support to
find and the monitor is CORRECT to say so -- calling that a false alarm would be
scoring the monitor on a member the loop would never have fed:
    LOUD   the k_loud = 5 members injected at h              <- the promotion number
    FAINT  the 11 members injected at H_FAINT = -14.25
    ALL    all 16
Wilson 95% intervals (Wald is wrong at p near 0, which is the whole point here), plus
the per-sky spread, because the 15 realisations within a sky share one geometry and are
NOT independent.

READS  GENERALISE_results/ (read-only) and the shared venue banks.
WRITES SIEVE_results/ only.  Gate G-V2a refuses to run if the shared GENERALISE L_gwb
bank is absent, because generalise._gen_lgwb would then WRITE into another campaign's
bank directory.

Output bank: SIEVE_results/sieve_v2_monitor_<cell>__<lane>.npz
Usage: python3 hpc_harbor/sieve/sv_v2_monitor.py [--cell 1] [--nreal 15]
"""
import os, sys, time
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import argparse
import glob
import socket
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
          "hpc_harbor/atlas", "hpc_harbor/chorus", "hpc_harbor/generalise"):
    sys.path.insert(0, f"{HSYMT}/{p}")

import chorus as CH
import generalise as GEN

OUT = os.environ.get("GLACIER_OUT", f"{HSYMT}/SIEVE_results")
H_ABSENT = -30.0                      # BASELINE: -18.0 is NOT absent at T = 30
I_H = GEN.I_H
NP_SRC = 8
N_POP = 16


# ---------------------------------------------------------------- helpers
def wilson(k, n, z=1.959963984540054):
    """Wilson score interval. Wald gives [0,0] at k=0, which would be a lie dressed as
    a measurement -- and k=0 is the outcome this test is built to detect."""
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(p), float(max(0.0, c - h)), float(min(1.0, c + h))


def theta_absent(theta, nd, idx):
    """theta with every member in `idx` set to H_ABSENT. Explicit index set -- never
    derived by set arithmetic on a carried list (the T2 degeneracy)."""
    t = np.asarray(theta, float).copy()
    for i in np.atleast_1d(np.asarray(idx, int)):
        t[nd + NP_SRC * int(i) + I_H] = H_ABSENT
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", type=int, default=1,
                    help="index into generalise.AS_CELLS; 1 = the SURVIVES cell")
    ap.add_argument("--nreal", type=int, default=GEN.N_SIG_C)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    lane = f"{os.environ.get('SLURM_JOB_ID', 'local')}"
    host = socket.gethostname()
    ncpu = len(os.sched_getaffinity(0))
    print(f"[V2] host={host} cpus={ncpu} job={lane} start={time.strftime('%FT%T')}",
          flush=True)

    # ---- G-V2a: the shared GENERALISE L_gwb bank must already exist ----------
    banks = sorted(glob.glob(f"{GEN.OUT}/gen_L_gwb_n*.npz"))
    if not banks:
        raise SystemExit("G-V2a FAIL: no gen_L_gwb_n*.npz in GENERALISE_results -- "
                         "running would WRITE into another campaign's bank. STOP.")
    print(f"[G-V2a] PASS shared L_gwb bank present: "
          f"{[os.path.basename(b) for b in banks]}", flush=True)

    jax, jnp, C, B1, TE, F, IG = GEN._import_stack()
    print(f"[V2] jax {jax.__version__} {jax.devices()}", flush=True)

    # keep the GLOBAL unit index -- `unit_entries` seeds off it, so re-deriving it by
    # searching a fresh as_units() list would be one refactor away from silently
    # scoring the wrong seeds.
    units = [(ui, u) for ui, u in enumerate(GEN.as_units()) if u["ci"] == a.cell]
    if not units:
        raise SystemExit(f"no A-SKY units with ci={a.cell}")
    h_cell, e_cell, k_cell = GEN.AS_CELLS[a.cell]
    print(f"[V2] cell {a.cell}: h={h_cell} e={e_cell} k_loud={k_cell}; "
          f"{len(units)} skies", flush=True)

    t0 = time.time()
    P = CH.build_chorus_problem(1, T_label=30)
    nd = P["n_dist"]; AI = P["AI"]; one = jnp.ones(P["npsr"])
    lb = P["amo"]["logL_batch_theta"]
    print(f"[V2] venue built in {time.time()-t0:.0f}s; lgwb: {P['nd'].lgwb_prov}",
          flush=True)
    lgwb_prov = str(P["nd"].lgwb_prov)

    rows = []          # one row per (sky, realisation, member)
    gate_c = None
    gate_d_ok = True
    prov_match = []

    for ui, u in units:
        sky = u["sky"]
        geo = F.load_geo_skies([sky])[0]
        G = GEN.gen_geometry(P, C, geo, u)
        theta_src = G["theta"]; dL = G["dL"]; EV = G["EV"]
        loud = np.array([i < u["k"] for i in range(N_POP)])
        tied = np.zeros(N_POP, bool); tied[list(u["placement"])] = True

        ents = [e for e in GEN.unit_entries("AS", ui, u)
                if e["kind"] == "sig"][:a.nreal]
        tsky = time.time()
        for ri, e in enumerate(ents):
            # ---- reproduce the banked realisation's data exactly -------------
            L_true, _ = CH.draw_true_distances_tier(P, dL, EV, seed=e["dist_seed"],
                                                    tier=u["tier"])
            theta_true = theta_src.copy(); theta_true[AI] = L_true
            clean = P["amo"]["inject_delay"](jnp.asarray(theta_true), one)
            noise = P["nd"].draw(e["noise_seed"])
            data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                         for c, n in zip(clean, noise))

            bank = GEN.real_path(u, e)
            if os.path.exists(bank):
                prov_match.append(str(np.load(bank, allow_pickle=True)["lgwb_prov"])
                                  == lgwb_prov)

            # ---- G-V2c (once): is H_ABSENT = -30 actually absent on this venue? --
            if gate_c is None:
                th30 = theta_absent(theta_true, nd, [0])
                th45 = theta_absent(theta_true, nd, [0]).copy()
                th45[nd + NP_SRC * 0 + I_H] = -45.0
                ll = np.asarray(lb(jnp.asarray(np.stack([th30, th45])), data, one))
                gate_c = float(ll[0] - ll[1])
                verdict = ("PASS" if abs(gate_c) < 1e-3 else
                           "FAIL -- -30 is NOT the absence plateau on this venue")
                print(f"[G-V2c] logL(-30) - logL(-45) = {gate_c:+.3e} nat  ({verdict})",
                      flush=True)

            # ---- the two contrasts, one batched call -------------------------
            all_idx = np.arange(N_POP)
            stack = [theta_true]                                   # 0: all present
            stack += [theta_absent(theta_true, nd, [i]) for i in all_idx]      # 1..16
            stack += [theta_absent(theta_true, nd, np.setdiff1d(all_idx, [i]))
                      for i in all_idx]                            # 17..32: only i
            stack += [theta_absent(theta_true, nd, all_idx)]       # 33: none present
            ths = np.stack(stack)
            # G-V2d: no contrast may be built from two identical templates
            if np.array_equal(ths[1], ths[0]) or np.array_equal(ths[33], ths[17]):
                gate_d_ok = False
            ll = np.asarray(lb(jnp.asarray(ths), data, one))
            d_incr = ll[0] - ll[1:1 + N_POP]
            d_iso = ll[1 + N_POP:1 + 2 * N_POP] - ll[33]

            for i in range(N_POP):
                rows.append((sky, int(e["noise_seed"]), i,
                             float(d_incr[i]), float(d_iso[i]),
                             float(h_cell if loud[i] else GEN.H_FAINT),
                             bool(loud[i]), bool(tied[i])))
            if ri == 0:
                print(f"  [s{sky} r0] incr loud={np.round(d_incr[:u['k']],2)} "
                      f"iso loud={np.round(d_iso[:u['k']],2)}", flush=True)
        print(f"  sky {sky}: {len(ents)} realisations ({time.time()-tsky:.0f}s)",
              flush=True)

    if not gate_d_ok:
        raise SystemExit("G-V2d FAIL: a present/absent pair is identical -- the "
                         "contrast is degenerate (the T2 failure mode). STOP.")
    print("[G-V2d] PASS: every present/absent template pair differs", flush=True)

    sky = np.array([r[0] for r in rows], int)
    seed = np.array([r[1] for r in rows], int)
    mem = np.array([r[2] for r in rows], int)
    d_incr = np.array([r[3] for r in rows], float)
    d_iso = np.array([r[4] for r in rows], float)
    h_inj = np.array([r[5] for r in rows], float)
    loud = np.array([r[6] for r in rows], bool)
    tied = np.array([r[7] for r in rows], bool)

    print(f"\n[V2] lgwb provenance match vs banked realisations: "
          f"{int(np.sum(prov_match))}/{len(prov_match)}", flush=True)

    print("\n=== FALSE-ALARM RATE (false alarm := dlnL <= 0 on a member present in "
          "the data) ===")
    res = {}
    for cname, contrast in (("incr", d_incr), ("iso", d_iso)):
        for sname, m in (("LOUD", loud), ("FAINT", ~loud), ("ALL",
                                                            np.ones_like(loud))):
            n = int(m.sum()); k = int((contrast[m] <= 0).sum())
            p, lo, hi = wilson(k, n)
            per = [float(np.mean(contrast[m & (sky == s)] <= 0))
                   for s in np.unique(sky)]
            res[f"{cname}_{sname}"] = (k, n, p, lo, hi)
            print(f"  {cname:5s} {sname:6s}  FAR = {k:4d}/{n:4d} = {p:.4f}  "
                  f"[{lo:.4f}, {hi:.4f}]  per-sky {min(per):.2f}-{max(per):.2f}  "
                  f"median dlnL = {np.median(contrast[m]):+.2f}")

    stem = f"sieve_v2_monitor_c{a.cell}"
    path = f"{a.out}/{stem}__{host}_{lane}.npz"
    np.savez(path, sky=sky, seed=seed, member=mem, dlnL_incr=d_incr, dlnL_iso=d_iso,
             h_inj=h_inj, loud=loud, tied=tied,
             cell=a.cell, h_cell=h_cell, e_cell=e_cell, k_loud=k_cell,
             h_absent=H_ABSENT, gate_c=(np.nan if gate_c is None else gate_c),
             gate_d_ok=gate_d_ok, lgwb_prov=lgwb_prov,
             prov_match=np.array(prov_match, bool),
             far_keys=np.array(sorted(res)),
             far_vals=np.array([res[k] for k in sorted(res)], float),
             host=host, job=lane, cpus=ncpu, time=time.strftime("%FT%T"))
    print(f"\n[V2] banked {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
