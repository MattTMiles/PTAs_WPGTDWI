#!/usr/bin/env python3
"""SIEVE T2 -- LISA-STYLE EVEN/ODD CROSS-VALIDATION OF THE FED TEMPLATE.

REPORT-ONLY. Nothing here arms a protocol step, moves a banked verdict, or enters a
closure claim. No campaign bank is written or overwritten; output goes to SIEVE_results/.

--------------------------------------------------------------------------------
THE QUESTION (brief): on banked GLACIER fed-cells, split each cell's TOAs into even and
odd epochs. Does the fed template fit the HELD-OUT half, at consistent parameters
(delta-lnL on the held-out half > 0)? Score truth-anchored feeds AND the banked WANDERED
states. PRE-REGISTERED: anchored PASS, wandered FAIL. If confirmed -> a free
wander-detector for the loop.

--------------------------------------------------------------------------------
THE HALF-VENUE, AND WHY IT IS THE SAME VENUE

A half is built by subsetting the DISCOVERY pulsars' per-TOA arrays (toas, stoas,
toaerrs, residuals, freqs, backend_flags, flags, pos_t, sunssb, planetssb, Mmat) to one
epoch parity AFTER the T-extension, then running the stock build. Three properties make
the halves comparable to the full venue and to each other, and all three are GATED, not
assumed:

  (1) IDENTICAL INJECTION. `make_geometry_injection` and the census override run off the
      ENTERPRISE pulsars, which are never subset, so `theta_keys`, `theta_truth`,
      `n_dist` and every source parameter are identical across full/even/odd.
  (2) IDENTICAL GP STRUCTURE. Both parities span the full baseline, so `ds.getspan` and
      hence the span-scaled `rn_comp`/`gwb_comp` counts are unchanged. The halves are not
      shorter datasets; they are sparser ones.
  (3) IDENTICAL SIGNAL. The CW delay is pointwise in TOA, so the half-venue's
      `inject_delay(theta)` must equal the full venue's restricted to that parity.
      GATE G-T2a checks this to 0 (bit-exact) on every pulsar.

DATA. The realisation is generated ONCE on the FULL venue exactly as `run_cell` does
(clean at the drawn distance truth + `NoiseDrawer.draw(noise_seed)`, white+rn, the
banked L_gwb, strict) so it reproduces the banked cell. Each half then takes that same
data vector restricted to its parity. Nothing is re-drawn per half -- re-drawing would
change the realisation and the test would compare two noise draws instead of two halves.

--------------------------------------------------------------------------------
THE STATISTIC

For a fed slot k, on a half H, at a template state theta:

    dlnL_H(k; theta) = logL_H(theta, k PRESENT, other carried ABSENT)
                     - logL_H(theta, k ABSENT,  all   carried ABSENT)

the own-term-live present-vs-absent contrast -- the SAME quantity frontier-v2 uses as
its data-support term (glacier_loop.run_cell step (b)), evaluated on a half instead of
on everything.

Two readings, both reported:

  (A) HELD-OUT WITH REFIT -- the genuine cross-validation. Refit the M-step on half A
      (stock `mstep_quadratic`, stock axes, stock scale) starting from the state under
      test, then score dlnL on half B at that fit. Reported both ways round (fit even /
      test odd, and fit odd / test even). PASS iff dlnL_heldout > 0.
      Parameter consistency: |theta_hat_A - theta_hat_B| on the M-step axes, in units of
      the M-step's own returned width.

  (B) NO-REFIT -- the banked template scored directly on each half. This is the CHEAP
      reading, and it is the one that could become a free wander-detector, because it
      costs two likelihood calls per fed member and no fit at all. It is NOT a held-out
      test for the wandered state (that state was fit on both halves), and is labelled
      accordingly wherever it is quoted.

STATES SCORED
  anchored: theta_true -- the drawn census truth (what `led.promote` feeds AT).
  wandered: the banked `theta_rec` at the requested iteration (post-M-step).

Output bank: SIEVE_results/sieve_t2_xval_<stem>_i<iter>__<lane>.npz
Usage:
  python3 hpc_harbor/sieve/sv_t2_xval.py --rung r13p25 --arm e07 --sky 0 --iter 5
  python3 hpc_harbor/sieve/sv_t2_xval.py --rung r13p5 --arm none --sky 3 --iter 0
"""
import os, sys, time
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import argparse
import glob
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
          "hpc_harbor/atlas", "hpc_harbor/spark", "hpc_harbor/glacier"):
    sys.path.insert(0, f"{HSYMT}/{p}")

import glacier_loop as GL
import glacier_pop as GP

PER_TOA = ("toas", "stoas", "toaerrs", "residuals", "freqs", "backend_flags",
           "pos_t", "sunssb", "planetssb", "Mmat")


def subset_pulsar(p, mask):
    """Subset every per-TOA array on a discovery pulsar in place. The attribute list is
    `ignite.extend_pulsar`'s own list -- if that function grows an array, this must too,
    so the two are kept literally parallel and a missing attribute is an ERROR."""
    n = len(p.toas)
    for attr in PER_TOA:
        v = getattr(p, attr, None)
        if v is None:
            raise RuntimeError(f"{p.name}: per-TOA attribute {attr!r} absent -- the "
                               f"half-venue subset would silently skip it.")
        v = np.asarray(v)
        if v.shape[0] != n:
            raise RuntimeError(f"{p.name}.{attr} has leading dim {v.shape[0]} != "
                               f"{n} TOAs -- refusing to subset blind.")
        setattr(p, attr, v[mask])
    p.flags = {k: (np.asarray(v)[mask] if np.asarray(v).shape[:1] == (n,) else v)
               for k, v in p.flags.items()}
    return p


def build_venue(T, pop_slots, parity=None, verbose=False):
    """`glacier_pop.build_glacier_problem`, with the discovery pulsars subset to one TOA
    parity after the T-extension. parity=None reproduces the stock build exactly."""
    import jax
    jax.config.update("jax_enable_x64", True)
    jax.config.update("jax_compilation_cache_dir",
                      os.environ.get("HARBOR_JAXCACHE_IN",
                                     "/home/mattm/.cache/jax_stagec_cache"))
    jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
    import trackB_b1_core as C
    import ignite as IG
    from cw_helpers import load_pulsars, build_enterprise_pta
    from stagec_fisher import make_geometry_injection, N_COMPONENTS, LOG10_EQUAD
    import discovery as ds

    n_src = pop_slots["n_src"]
    ent_psrs, disco_psrs = load_pulsars(116)
    pta, cwb, _ = build_enterprise_pta(ent_psrs, n_src, components=N_COMPONENTS,
                                       log10_equad=LOG10_EQUAD)
    inj = make_geometry_injection(pta, ent_psrs, n_src, cwb, seed=pop_slots["seed"])
    for i, name in enumerate(cwb):
        s = pop_slots["src"][i]
        inj[f"{name}_cos_gwtheta"] = s[GP.I_COSGT]
        inj[f"{name}_gwphi"] = s[GP.I_GWPHI]
        inj[f"{name}_cos_inc"] = s[GP.I_COSINC]
        inj[f"{name}_log10_mc"] = s[GP.I_MC]
        inj[f"{name}_log10_fgw"] = s[GP.I_FGW]
        inj[f"{name}_log10_h"] = s[GP.I_H]
        inj[f"{name}_phase0"] = s[GP.I_PH0]
        inj[f"{name}_psi"] = s[GP.I_PSI]
    inj["gwb_log10_A"] = GP.A_TARGET_LOG10
    inj["gwb_gamma"] = GP.GAMMA_BG

    span0 = ds.getspan(disco_psrs)
    dT = float(T - 15)
    if dT > 0:
        for p in disco_psrs:
            IG.extend_pulsar(p, dT)
        span1 = ds.getspan(disco_psrs)
        rn_comp = int(round(30 * span1 / span0))
        gwb_comp = int(round(N_COMPONENTS * span1 / span0))
    else:
        rn_comp, gwb_comp = 30, N_COMPONENTS

    masks = None
    if parity is not None:
        masks = []
        for p in disco_psrs:
            m = np.zeros(len(p.toas), bool)
            m[int(parity)::2] = True
            masks.append(m)
            subset_pulsar(p, m)
        span2 = ds.getspan(disco_psrs)
        if abs(span2 / span1 - 1.0) > 0.02:
            raise RuntimeError(f"parity {parity}: span moved {span1:.4e} -> {span2:.4e}; "
                               f"the halves must span the same baseline as the full "
                               f"venue or the GP bases are not comparable.")

    amo = C.build_b1_amortised(disco_psrs, n_src, inj, cwb,
                               components=gwb_comp, rn_components=rn_comp)
    return dict(amo=amo, disco_psrs=disco_psrs, ent_psrs=ent_psrs, cwb=cwb, inj=inj,
                T_label=T, rn_comp=rn_comp, gwb_comp=gwb_comp, masks=masks,
                span_s=float(ds.getspan(disco_psrs)))


def dlnl_member(amo, jnp, theta, nd, carried, k, data, ones):
    """The own-term-live present-vs-absent contrast for slot k.

    NB the absent set must be built as `carried UNION {k}` and then k removed for the
    "on" leg. `run_cell`'s frontier-v2 can write `setdiff1d(carried, [k])` because there
    k is a CANDIDATE and so is still in `carried`; T2 scores slots that are ALREADY FED,
    for which k is NOT in `carried`, and the naive form collapses to th_on == th_off and
    returns identically 0.000. That is exactly what the first scored run produced
    (job 12834899), which is why the degeneracy is asserted against below rather than
    trusted.
    """
    k = int(k)
    absent_k = np.union1d(np.asarray(carried, int), np.array([k]))
    th_on = GL.theta_with_absent(theta, nd, np.setdiff1d(absent_k, [k]))
    th_off = GL.theta_with_absent(theta, nd, absent_k)
    if np.array_equal(th_on, th_off):
        raise RuntimeError(f"slot {k}: present and absent templates are identical -- "
                           f"the contrast is degenerate and would return 0.")
    ll = np.asarray(amo["logL_batch_theta"](
        jnp.asarray(np.stack([th_on, th_off])), data, ones))
    return float(ll[0] - ll[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rung", default="r13p25")
    ap.add_argument("--arm", default="e07")
    ap.add_argument("--sky", type=int, default=0)
    ap.add_argument("--t", type=int, default=30)
    ap.add_argument("--iter", type=int, default=5)
    ap.add_argument("--max-fed", type=int, default=8,
                    help="score at most this many fed slots (brightest-first slot order)")
    a = ap.parse_args()
    GL.check_affinity()

    stem = f"gl2_{a.rung}_cell_{a.arm}_s{a.sky}_T{a.t}_{GL.TIER}"
    # READ the campaign bank from GLACIER_results explicitly. GP.OUT is repointed to
    # SIEVE_results by the launcher so that bank_npz below can never write into the
    # campaign directory; the read path must therefore be stated, not inherited.
    fs = sorted(glob.glob(f"{HSYMT}/GLACIER_results/{stem}_i{a.iter}__*.npz"))
    if not fs:
        print(f"no banked cell {stem}_i{a.iter}", flush=True)
        return 2
    z = np.load(fs[0], allow_pickle=True)
    fed_mask = np.asarray(z["fed_mask"], bool)
    theta_wander = np.asarray(z["theta_rec"], float)
    print(f"[T2] cell {stem} i{a.iter}: {int(fed_mask.sum())} fed of {len(fed_mask)} "
          f"slots  (bank {os.path.basename(fs[0])})", flush=True)

    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()
    pop, cond = GP.draw_population_conditional(a.sky, a.rung, n_src=256,
                                               band_log10f=GL.BAND_CAMPAIGN)
    slots, n_harm, active, chan, n_clip = GL.embed_igniter(
        pop, GL.E_IGNITER[a.arm], GL.VENUE_SPAN_S[a.t])
    pop_s = dict(pop); pop_s["src"] = slots; pop_s["n_src"] = len(slots)
    if len(fed_mask) != len(slots):
        raise RuntimeError(f"rebuilt {len(slots)} slots but the bank has "
                           f"{len(fed_mask)} -- census reproduction broken.")

    t0 = time.time()
    Vf = build_venue(a.t, pop_s, parity=None)
    print(f"[T2] full venue {time.time()-t0:.0f}s", flush=True)
    amo = Vf["amo"]; nd = amo["n_dist"]
    ones = jnp.ones(amo["npsr"])

    # --- the realisation, exactly as run_cell built it ---
    import stagec_anchor_a2 as A2
    from cw_helpers import compute_mode_spacing
    ent = Vf["ent_psrs"]
    L0 = np.array([pe.pdist[0] for pe in ent])
    dL = np.array([min(compute_mode_spacing(s[GP.I_COSGT], s[GP.I_GWPHI], s[GP.I_FGW],
                                            Vf["disco_psrs"][i].pos) for s in slots)
                   for i in range(len(ent))])
    Pd = dict(npsr=len(ent), L0=L0,
              priors=A2.load_real_priors([pe.name for pe in ent]), ent_psrs=ent)
    EV = C.build_EV(Pd, dL)
    noise_seed = GL.SEED_NOISE_BASE + 1000 * a.sky
    L_true, _ = IG.draw_true_distances_tier(Pd, dL, EV, seed=noise_seed + 10_000_000,
                                            tier=GL.TIER)
    theta_true = np.asarray(amo["theta_truth"], float).copy()
    theta_true[np.arange(nd)] = L_true
    clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
    theta_truth0 = np.asarray(amo["theta_truth"], float)
    data_zero = amo["inject_delay"](jnp.asarray(theta_truth0), ones)
    sp = C.B1Split(amo, theta_truth0, data_zero, np.ones(amo["npsr"]))
    ndraw = C.NoiseDrawer(sp, amo, lgwb_path=GL.LGWB_BANKS[a.t], strict=True)
    noise = ndraw.draw(noise_seed, components=("white", "rn"))
    data_full = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                      for c, n in zip(clean, noise))
    print(f"[T2] realisation built (noise_seed {noise_seed})", flush=True)

    fed_idx_all = np.where(fed_mask)[0]
    if len(fed_idx_all) > a.max_fed:
        print(f"[T2] SCOPE: {int(fed_mask.sum())} fed slots; scoring the first "
              f"{a.max_fed} in slot order (brightest-first). The rest are NOT scored "
              f"and the counts below say so.", flush=True)
    fed_idx = fed_idx_all[:a.max_fed]
    carried = np.where(~fed_mask)[0]
    states = {"anchored": theta_true, "wandered": theta_wander}

    # keep only host arrays; nothing jax-side from the full venue survives this point
    clean_np = tuple(np.asarray(c) for c in clean)
    data_np = tuple(np.asarray(d) for d in data_full)
    theta_keys_full = list(amo["theta_keys"])
    theta_truth0 = np.asarray(theta_truth0, float)
    npsr_full = int(amo["npsr"])

    # --- FREE THE FULL VENUE ---------------------------------------------------------
    # THE CPU-LANE MAP-COUNT HAZARD, hit twice for real. Job 12834458: three live
    # `build_b1_amortised` venues at ncw=287 x 116 pulsars exhausted vm.max_map_count
    # (65530 on these nodes) and the third build died in XLA-CPU with "INTERNAL: Failed
    # to materialize symbols". Job 12834779: freeing the full venue was not enough,
    # because COMPILING its logL on top of a live NoiseDrawer (a 227 MB banked L_gwb +
    # the refactored GWB block) already exhausted the maps -- LLVM
    # `allocateMappedMemory failed: Cannot allocate memory`. Each JIT section is an mmap,
    # so compiles cost maps, not just RAM.
    #
    # THE DESIGN IS THEREFORE STRICTLY ONE VENUE LIVE AT A TIME, and the full venue's
    # logL is NEVER compiled -- it is only ever asked for `inject_delay` (needed for the
    # clean signal and the gate reference) and the noise draw. The consequence, declared:
    # the `dlnL_full` reference column is NOT produced. It was a nice-to-have; the
    # pre-registered test is on the HALVES and is unaffected.
    import gc
    del Vf, amo, sp, ndraw, clean, data_full, EV
    gc.collect()
    jax.clear_caches()
    print("[T2] full venue released WITHOUT compiling its logL "
          "(CPU map-count hazard: peak 1 live venue)", flush=True)

    scale = TE.phi_scale({"ncw": 1})
    rows = []
    gate = {}
    norefit = {}          # (state, slot, half) -> dlnL
    fits = {}             # (state, slot, half) -> (src_hat, widths)
    heldout = {}          # (state, slot, fit_half) -> dlnL on the OTHER half

    def open_half(nm, parity):
        t = time.time()
        V = build_venue(a.t, pop_s, parity=parity)
        print(f"[T2] {nm} venue {time.time()-t:.0f}s  "
              f"TOAs {sum(int(m.sum()) for m in V['masks'])}", flush=True)
        Vh = V["amo"]
        if list(Vh["theta_keys"]) != theta_keys_full:
            raise RuntimeError(f"{nm}: theta_keys differ from the full venue.")
        if not np.allclose(np.asarray(Vh["theta_truth"]), theta_truth0, rtol=0, atol=0):
            raise RuntimeError(f"{nm}: theta_truth differs from the full venue.")
        ch = Vh["inject_delay"](jnp.asarray(theta_true), jnp.ones(Vh["npsr"]))
        d = max(float(np.max(np.abs(np.asarray(ch[i]) - clean_np[i][V["masks"][i]])))
                for i in range(Vh["npsr"]))
        gate[nm] = d
        print(f"[T2] GATE G-T2a {nm}: max |half - full[mask]| = {d:.3e}"
              f"  {'PASS' if d == 0.0 else 'FAIL'}", flush=True)
        dat = tuple(jnp.asarray(data_np[i][V["masks"][i]]) for i in range(npsr_full))
        return V, Vh, dat, jnp.ones(Vh["npsr"])

    def close_half():
        """Callers MUST drop their own references first (`V = Vh = dat = one_h = None`)
        -- `del` inside this function would only unbind the local name and the venue
        would stay alive, which is precisely the failure this whole restructure exists
        to avoid."""
        gc.collect()
        jax.clear_caches()

    def do_fit(Vh, dat, one_h, th, k):
        th_eval = GL.theta_with_absent(th, nd, np.setdiff1d(carried, [k]))

        def marg_fn(srcs):
            ths = np.tile(th_eval, (len(srcs), 1))
            ths[:, nd:] = srcs
            return np.asarray(Vh["logL_batch_theta"](jnp.asarray(ths), dat, one_h))

        sh, w, _ = GL.mstep_quadratic(marg_fn, th_eval[nd:].copy(), [int(k)], scale)
        return sh, np.asarray(w, float).ravel()

    # ---- PASS 1: EVEN. no-refit column + the even fits. --------------------------
    print("\n[T2] PASS 1 -- even half")
    V, Vh, dat, one_h = open_half("even", 0)
    for sname, th in states.items():
        for k in fed_idx:
            norefit[(sname, int(k), "even")] = dlnl_member(Vh, jnp, th, nd, carried, k,
                                                           dat, one_h)
            fits[(sname, int(k), "even")] = do_fit(Vh, dat, one_h, th, k)
    print(f"[T2] even: {len(fed_idx)*len(states)} no-refit scores + fits done",
          flush=True)
    V = Vh = dat = one_h = None
    close_half()

    # ---- PASS 2: ODD. no-refit column, the odd fits, and the even->odd held-out. --
    print("\n[T2] PASS 2 -- odd half (held-out for the even fits)")
    V, Vh, dat, one_h = open_half("odd", 1)
    for sname, th in states.items():
        for k in fed_idx:
            norefit[(sname, int(k), "odd")] = dlnl_member(Vh, jnp, th, nd, carried, k,
                                                          dat, one_h)
            fits[(sname, int(k), "odd")] = do_fit(Vh, dat, one_h, th, k)
            th_fit = th.copy(); th_fit[nd:] = fits[(sname, int(k), "even")][0]
            heldout[(sname, int(k), "even")] = dlnl_member(Vh, jnp, th_fit, nd, carried,
                                                           k, dat, one_h)
    V = Vh = dat = one_h = None
    close_half()

    # ---- PASS 3: EVEN again, for the reverse direction (odd fit -> even test). ----
    # A fourth build, ~10 min, so that the verdict does not rest on one direction.
    print("\n[T2] PASS 3 -- even half rebuilt (held-out for the odd fits)")
    try:
        V, Vh, dat, one_h = open_half("even2", 0)
        for sname, th in states.items():
            for k in fed_idx:
                th_fit = th.copy(); th_fit[nd:] = fits[(sname, int(k), "odd")][0]
                heldout[(sname, int(k), "odd")] = dlnl_member(Vh, jnp, th_fit, nd,
                                                              carried, k, dat, one_h)
        V = Vh = dat = one_h = None
        close_half()
    except Exception as ex:
        print(f"[T2] PASS 3 FAILED ({type(ex).__name__}: {ex}). The reverse direction "
              f"is NOT reported; passes 1-2 stand and the counts below say so.",
              flush=True)

    gate_max = max(gate.values()) if gate else np.nan
    if not np.isfinite(gate_max) or gate_max > 0.0:
        print(f"   GATE FAIL: {gate_max}. The halves are not the same venue; "
              f"no downstream number is quotable.")
        return 1
    print(f"\n[T2] GATE G-T2a PASS: bit-exact on every half built "
          f"({', '.join(gate)})\n")

    # ---- readout ----
    print("[T2] (B) NO-REFIT: the state scored directly on each half")
    print("   state      slot    dlnL_even     dlnL_odd   both>0")
    for sname in states:
        for k in fed_idx:
            de = norefit[(sname, int(k), "even")]
            do = norefit[(sname, int(k), "odd")]
            ok = de > 0 and do > 0
            print(f"   {sname:9s} {int(k):5d} {de:12.3f} {do:12.3f}   "
                  f"{'YES' if ok else 'no'}")
            rows.append(["norefit", sname, int(k), np.nan, de, do,
                         np.nan, np.nan, np.nan, np.nan, bool(ok)])

    print("\n[T2] (A) HELD-OUT WITH REFIT: fit one half, score the other")
    print("   state      slot  fit->test   dlnL_heldout   d(fgw)/w   d(mc)/w   verdict")
    for sname in states:
        for k in fed_idx:
            se, we = fits[(sname, int(k), "even")]
            so, wo = fits[(sname, int(k), "odd")]
            b = GP.NP_SRC * int(k)
            dfgw = float(se[b + GP.I_FGW] - so[b + GP.I_FGW])
            dmc = float(se[b + GP.I_MC] - so[b + GP.I_MC])
            wf = float(we[0]) if len(we) > 0 and np.isfinite(we[0]) else np.nan
            wm = float(we[1]) if len(we) > 1 and np.isfinite(we[1]) else np.nan
            for fitnm, testnm in (("even", "odd"), ("odd", "even")):
                if (sname, int(k), fitnm) not in heldout:
                    continue
                dh = heldout[(sname, int(k), fitnm)]
                rf = dfgw / wf if wf and np.isfinite(wf) else np.nan
                rm = dmc / wm if wm and np.isfinite(wm) else np.nan
                print(f"   {sname:9s} {int(k):5d} {fitnm}->{testnm:5s} {dh:13.3f} "
                      f"{rf:10.2f} {rm:9.2f}   {'PASS' if dh > 0 else 'FAIL'}")
                rows.append(["refit", sname, int(k), np.nan,
                             dh if testnm == "even" else np.nan,
                             dh if testnm == "odd" else np.nan,
                             dfgw, dmc, wf, wm, bool(dh > 0)])

    # --- verdict ---
    print("\n[T2] VERDICT TABLE")
    for mode in ("norefit", "refit"):
        for sname in states:
            sel = [r for r in rows if r[0] == mode and r[1] == sname]
            npass = sum(1 for r in sel if r[10])
            print(f"   {mode:8s} {sname:9s}: {npass}/{len(sel)} PASS "
                  f"(held-out dlnL > 0)")

    path = GP.bank_npz(
        f"sieve_t2_xval_{a.rung}_{a.arm}_s{a.sky}_i{a.iter}",
        note=("SIEVE T2: even/odd TOA cross-validation of the fed template on a banked "
              "GLACIER cell. Halves are stock builds with the discovery pulsars subset "
              "to one TOA parity after the T-extension; the realisation is built ONCE on "
              "the full venue and restricted. GATE G-T2a (half inject_delay == full "
              "restricted, bit-exact) passed. REPORT-ONLY."),
        gate_gt2a_max=gate_max, gate_pass=True,
        gate_per_half=np.array([[k, v] for k, v in gate.items()],
                               dtype=object),
        dlnL_full_available=False,
        cell_stem=stem, iter=a.iter, rung=a.rung, arm=a.arm, sky=a.sky, t_label=a.t,
        n_fed_total=int(fed_mask.sum()), n_fed_scored=len(fed_idx),
        fed_scored=np.array(fed_idx),
        noise_seed=noise_seed,
        rows=np.array(rows, dtype=object),
        rows_keys=np.array(["mode", "state", "slot", "dlnL_full", "dlnL_even",
                            "dlnL_odd", "d_fgw", "d_mc", "w_fgw", "w_mc", "pass"]))
    print(f"\nbanked -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
