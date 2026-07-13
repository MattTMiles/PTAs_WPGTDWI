"""SAMPLER — THE SUB-FRINGE OFFSET TEST, re-spec'd (§9.4): the u-grid and its three reductions.

WHY (SAMPLER_dev.md §9)
-----------------------
`A2.eval_grid` returns FRINGE CENTRES (one point per fringe) whenever the +-3 sigma window holds
more than DENSE_FRINGE_MAX = 64 fringes -- 111/116 pulsars; four more get centres + `L0` PADDING
DUPLICATES; one gets a dense grid. So the sub-fringe offset `u` is **PINNED AT u = 0**, not
profiled and not marginalised, in every banked number. And because the fringe spacing is DEFINED by
dPhi_p(dL) = 2*pi, every fringe centre carries the SAME pulsar-term phase (the phase at L0) --
so in ARM B, where truth sits at L0 + (n_true + u)*dL with u ~ U(-1/2,1/2), the likelihood is
evaluated with a pulsar-term phase wrong by 2*pi*u AT EVERY FRINGE, IDENTICALLY.

THE TEST
--------
Evaluate lnL on a U-GRID -- L = L0 + (n + u_m) * dL for the pulsar's retained fringes n and a set of
offsets u_m -- and score THREE reductions off the SAME evaluations:

    CENTRE    m_n = lnL(n, u=0)                      + lnprior(centre_n)     [the banked object]
    PROFILE   m_n = max_u lnL(n, u)                  + lnprior(centre_n)     [what snap_to_peak did]
    MARGINAL  m_n = LSE_u[ lnL(n,u) + lnprior(L_nu) ] + log(dL * du)         [the honest object]

Implementation note that matters: the u-axis is swept by REPEATED CALLS to the banked split E-step
(one call per offset, each on the pulsar's UNIQUE fringe centres shifted by u*dL). This (a) reuses
the gated `B1Split.estep` verbatim, (b) is memory-safe (never a K x M stack), and (c) **de-duplicates
the `L0` padding by construction** -- the grid is built from `np.unique` fringe indices, never from
the padded EV columns. A naive logsumexp over the banked EV would have counted the padding up to
374x on four pulsars (§9.3).

GATES (pre-registered)
----------------------
  G-REPRO  CENTRE reproduces a BANKED FORGE Arm-B flat (`reports/flat_b10_g0_B_n9000xx.npz`):
           q_max, map_k bit-comparable, and the banked `dlnL_det` definition is IDENTIFIED by
           reproduction (FORGE's scorer is not in this repo, so the definition is disambiguated
           against the bank rather than assumed).
  G-ARM-A  In ARM A (u_true = 0 by construction) CENTRE == PROFILE: the within-fringe likelihood
           peak must SIT at u = 0. If it does not, the phase argument of §9.2 is wrong.
  G-ARM-B  In ARM B the within-fringe argmax must sit at u ~ u_true, THE SAME u FOR EVERY FRINGE.
           That is the direct, falsifiable form of the "common phase error" claim.

Env: jug-gpu, x64. Reads banked npz; writes only to SAMPLER_results/. Matt commits; never here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time, argparse, glob
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import sampler_core as S
import trackB_b1_core as C
import trackB_estimator as TE

RES = S.RES
REPO = "/home/mattm/projects/HSYMT"


def _lse(x, axis=-1):
    m = np.max(x, axis=axis, keepdims=True)
    return np.squeeze(m + np.log(np.sum(np.exp(x - m), axis=axis, keepdims=True)), axis=axis)


class UGrid:
    """Per-pulsar unique fringe set + offset grid; sweeps u by repeated split-E-step calls."""

    def __init__(self, P, prior_key="lit", M=16):
        self.P = P
        self.npsr = P["npsr"]
        self.L0 = np.asarray(P["L0"], float)
        self.dL = np.asarray(P["dL_truth"], float)
        self.sig = np.asarray(P["priors"][prior_key], float)
        EV = np.asarray(P["EV_truth"], float)
        # UNIQUE fringe indices actually sampled by the banked grid (this DROPS the L0 padding)
        self.uk = []
        for a in range(self.npsr):
            k = np.round((EV[a] - self.L0[a]) / self.dL[a]).astype(int)
            self.uk.append(np.unique(k))
        self.K = np.array([len(u) for u in self.uk])
        self.Kmax = int(self.K.max())
        # padded fringe-index table; padded slots carry lnprior = -inf and never contribute
        self.kk = np.zeros((self.npsr, self.Kmax), int)
        self.valid = np.zeros((self.npsr, self.Kmax), bool)
        for a in range(self.npsr):
            self.kk[a, :self.K[a]] = self.uk[a]
            self.valid[a, :self.K[a]] = True
        # offsets: u = 0 (the banked CENTRE) + M midpoints spanning (-1/2, +1/2), no endpoint dup
        self.u_mid = (np.arange(M) + 0.5) / M - 0.5
        self.u_all = np.concatenate([[0.0], self.u_mid])
        self.M = M
        self.du = 1.0 / M
        # prior at fringe centres (the banked convention) and pointwise on the u-grid
        self.lnprior_c = np.where(self.valid,
                                  -(self.kk * self.dL[:, None]) ** 2 / (2 * self.sig[:, None] ** 2),
                                  -np.inf)

    def lnprior_at(self, u):
        off = (self.kk + u) * self.dL[:, None]                       # L - L0
        return np.where(self.valid, -off ** 2 / (2 * self.sig[:, None] ** 2), -np.inf)

    def L_at(self, u):
        return self.L0[:, None] + (self.kk + u) * self.dL[:, None]   # (npsr, Kmax)

    def sweep(self, sp, theta_base, data, mask, AI, verbose=True):
        """LNL[m, a, k] over the offset grid. One split-E-step call per offset."""
        t0 = time.time()
        out = np.empty((len(self.u_all), self.npsr, self.Kmax))
        for m, u in enumerate(self.u_all):
            EV_u = np.maximum(self.L_at(u), 1e-6)                    # keep distances physical
            L = sp.estep(theta_base, EV_u, AI, data, mask)
            C._finite(f"ugrid estep u={u:+.3f}", L)
            out[m] = L
            if verbose and (m == 0 or m == len(self.u_all) - 1):
                print(f"     u = {u:+.4f}: E-step done ({time.time()-t0:.0f}s cum)", flush=True)
        return out

    # ---------------- the three reductions ----------------
    def reduce_all(self, LNL):
        """LNL (M+1, npsr, Kmax) -> per-fringe log-weights m_n under each reduction."""
        cen = LNL[0] + self.lnprior_c                                       # (npsr, Kmax)
        prof = LNL[1:].max(axis=0) + self.lnprior_c
        lp = np.stack([self.lnprior_at(u) for u in self.u_mid])             # (M, npsr, Kmax)
        marg = _lse(LNL[1:] + lp, axis=0) + np.log(self.dL * self.du)[:, None]
        marg = np.where(self.valid, marg, -np.inf)
        return dict(CENTRE=cen, PROFILE=prof, MARGINAL=marg)

    def argmax_u(self, LNL):
        """within-fringe argmax offset per (pulsar, fringe) -- the phase test."""
        j = np.argmax(LNL[1:], axis=0)
        return self.u_mid[j]                                                # (npsr, Kmax)

    def qpost(self, m_n):
        """q_p(n), MAP fringe, and the fringe-gap statistic, from per-fringe log-weights."""
        npsr = self.npsr
        q = np.zeros((npsr, self.Kmax))
        qmax = np.zeros(npsr); mapk = np.zeros(npsr, int)
        for a in range(npsr):
            w = m_n[a, :self.K[a]]
            w = w - w.max()
            e = np.exp(w); e /= e.sum()
            q[a, :self.K[a]] = e
            j = int(np.argmax(e))
            qmax[a] = e[j]; mapk[a] = self.uk[a][j]
        return q, qmax, mapk

    def ptrue(self, q, n_true):
        out = np.zeros(self.npsr)
        for a in range(self.npsr):
            hit = np.where(self.uk[a] == n_true[a])[0]
            out[a] = q[a, hit[0]] if len(hit) else 0.0
        return out

    def gap(self, LNL_like):
        """dlnL = best - runner-up over fringes, on a LIKELIHOOD-only per-fringe column."""
        g = np.zeros(self.npsr)
        for a in range(self.npsr):
            v = np.sort(LNL_like[a, :self.K[a]])[::-1]
            g[a] = v[0] - v[1] if len(v) > 1 else np.inf
        return g


# ============================================================
# realisation construction (the FORGE Arm-B T=15 recipe, deterministic from banked seeds)
# ============================================================
def make_realisation(P, arm, dist_seed, noise_seed):
    tt = P["theta_truth"].copy()
    npsr = P["npsr"]
    if arm == "B":
        L_true, n_true, u_true = C.draw_true_distances(P, seed=dist_seed)
        tt[:P["n_dist"]] = L_true
    else:
        n_true = np.zeros(npsr, int); u_true = np.zeros(npsr)
    mask1 = P["mask_one"]
    clean = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(mask1))
    if noise_seed is None:
        data = clean
    else:
        nz = P["nd"].draw(seed=noise_seed)
        data = tuple(jnp.asarray(np.asarray(d) + np.asarray(n)) for d, n in zip(clean, nz))
    C._finite("realisation data", *[np.asarray(d) for d in data])
    # the CONDITIONAL pipeline scores at the TRUE source with BASE (prior-mean) distances
    theta_base = np.concatenate([P["L0"], tt[P["n_dist"]:]])
    return dict(theta_true=tt, theta_base=theta_base, data=data, mask=mask1,
                n_true=n_true, u_true=u_true)


def mode_gate(a):
    print("=== U-GRID GATES ===", flush=True)
    P = C.build_b1_problem()
    sp = P["sp"]; AI = P["AI"]
    G = UGrid(P, "lit", M=a.M)
    print(f"  offsets: u = 0 plus {a.M} midpoints; fringes/pulsar K in "
          f"[{G.K.min()}, {G.K.max()}]; L0-padding DROPPED by construction", flush=True)

    # ---------- G-REPRO: CENTRE vs a banked FORGE Arm-B flat ----------
    f = f"{REPO}/reports/flat_b10_g0_B_n{a.flat}.npz"
    Z = np.load(f, allow_pickle=True)
    print(f"\n-- G-REPRO: CENTRE vs {os.path.basename(f)} "
          f"(arm {Z['arm']}, dist_seed {int(Z['dist_seed'])}, noise_seed {int(Z['noise_seed'])}) --",
          flush=True)
    R = make_realisation(P, str(Z["arm"]), int(Z["dist_seed"]), int(Z["noise_seed"]))
    nt_bank = np.asarray(Z["n_true_grid"])
    d_nt = int(np.sum(R["n_true"] != nt_bank))
    print(f"   n_true reproduced: {P['npsr']-d_nt}/{P['npsr']} pulsars match the bank", flush=True)

    LNL = G.sweep(sp, R["theta_base"], R["data"], R["mask"], AI)
    red = G.reduce_all(LNL)
    q_c, qmax_c, mapk_c = G.qpost(red["CENTRE"])
    # The <=64-fringe pulsars get a DENSE uniform grid from eval_grid, so B1Marg's segment_max
    # there really IS a within-fringe profile -- my CENTRE (u = 0 exactly) is NOT expected to
    # reproduce the bank on those. Classify rather than conflate.
    nf = np.array([int((6.0 * G.sig[p]) / G.dL[p]) + 1 for p in range(P["npsr"])])
    dense = nf <= 64
    dq_all = np.abs(qmax_c - np.asarray(Z["qmax"]))
    mm = mapk_c != np.asarray(Z["mapk"])
    print(f"   grid regime: {int((~dense).sum())} centre-grid pulsars (CENTRE must reproduce), "
          f"{int(dense.sum())} dense-grid (bank profiles there -- disagreement EXPECTED)",
          flush=True)
    print(f"   q_max : max|diff| = {float(dq_all[~dense].max()):.3e} (centre-grid) | "
          f"{float(dq_all[dense].max()) if dense.any() else 0.0:.3e} (dense-grid)", flush=True)
    print(f"   map_k : {int(mm[~dense].sum())}/{int((~dense).sum())} disagree (centre-grid) | "
          f"{int(mm[dense].sum())}/{int(dense.sum())} (dense-grid)", flush=True)
    dq = float(dq_all[~dense].max())
    dm = int(mm[~dense].sum())
    # identify FORGE's dlnL_det by reproduction: gap on the prior-weighted vs likelihood-only column
    cand = dict(gap_post=G.gap(red["CENTRE"]), gap_lik=G.gap(LNL[0]))
    bank_d = np.asarray(Z["dlnL_det"])
    for nm, v in cand.items():
        print(f"   dlnL candidate '{nm}': max|diff| vs banked dlnL_det = "
              f"{float(np.max(np.abs(v - bank_d))):.3e}", flush=True)
    ok_repro = (dq < 1e-6) and (dm == 0) and (d_nt == 0)
    print(f"   [G-REPRO] -> {'PASS' if ok_repro else 'FAIL'}", flush=True)

    # ---------- G-ARM-A / G-ARM-B: where does the within-fringe peak sit? ----------
    for arm, dseed in (("A", None), ("B", int(Z["dist_seed"]))):
        print(f"\n-- G-ARM-{arm}: within-fringe argmax offset u* (the PHASE test) --", flush=True)
        RR = make_realisation(P, arm, dseed if dseed else 0, None)   # ZERO-NOISE, isolate the phase
        L2 = G.sweep(sp, RR["theta_base"], RR["data"], RR["mask"], AI, verbose=False)
        us = G.argmax_u(L2)
        ut = RR["u_true"]
        # THE SHARP TEST: u* AT THE TRUE FRINGE. On WRONG fringes the pulsar-term model is
        # decorrelated from the data, the likelihood is nearly flat in u, and its argmax is
        # meaningless noise (measured: std_n u* ~ 0.42, i.e. ~uniform) -- averaging that in
        # swamps the one fringe that carries the offset. Scored separately, and the
        # true-fringe peak is the one the CENTRE convention is supposed to be sitting on.
        ustar_true = np.full(P["npsr"], np.nan)
        for p in range(P["npsr"]):
            hit = np.where(G.uk[p] == RR["n_true"][p])[0]
            if len(hit):
                ustar_true[p] = us[p, hit[0]]
        good = np.isfinite(ustar_true)
        err_t = np.abs(((ustar_true[good] - ut[good] + 0.5) % 1.0) - 0.5)
        print(f"   [TRUE-FRINGE] |u*(n_true) - u_true| : median {np.median(err_t):.4f}, "
              f"90th pct {np.percentile(err_t, 90):.4f}  (grid resolution {G.du/2:.4f}; "
              f"{int(good.sum())}/{P['npsr']} pulsars have n_true on the grid)", flush=True)
        print(f"   [TRUE-FRINGE] the CENTRE convention evaluates at u = 0, i.e. off the true "
              f"peak by |u_true| = {np.median(np.abs(ut)):.4f} (median)", flush=True)
        np.savez(f"{RES}/ugrid_phase_true_arm{arm}.npz", ustar_true=ustar_true,
                 u_true=ut, err=err_t)
        # the OLD (mis-specified) all-fringe statistic, kept for the record
        med = np.array([np.median(us[p, :G.K[p]]) for p in range(P["npsr"])])
        spread = np.array([np.std(us[p, :G.K[p]]) for p in range(P["npsr"])])
        print(f"   median over pulsars of (median_n u*)        = {np.median(med):+.4f}", flush=True)
        print(f"   median over pulsars of (std_n u*)           = {np.median(spread):.4f}  "
              f"(0 => the SAME offset at every fringe)", flush=True)
        if arm == "B":
            err = np.abs(((med - ut + 0.5) % 1.0) - 0.5)
            print(f"   |median_n u*  -  u_true| : median {np.median(err):.4f}, "
                  f"90th pct {np.percentile(err, 90):.4f}   (0 => peak sits AT the true offset)",
                  flush=True)
            print(f"   -> the phase claim of §9.2 is {'CONFIRMED' if np.median(err) < 0.1 else 'REFUTED'}",
                  flush=True)
        else:
            print(f"   |median_n u*| = {np.median(np.abs(med)):.4f}  "
                  f"(0 => CENTRE is exact in arm A, as claimed)", flush=True)
        np.savez(f"{RES}/ugrid_phase_arm{arm}.npz", u_star=us, u_true=RR["u_true"],
                 med=med, spread=spread, K=G.K)
    return 0


def mode_impact(a):
    """Impact rows: CENTRE -> MARGINAL shifts on >=10 banked Arm-B T=15 realisations."""
    print(f"=== U-GRID IMPACT: CENTRE vs PROFILE vs MARGINAL, {a.n} Arm-B T=15 realisations ===",
          flush=True)
    P = C.build_b1_problem()
    sp = P["sp"]; AI = P["AI"]
    G = UGrid(P, "lit", M=a.M)
    flats = sorted(glob.glob(f"{REPO}/reports/flat_b10_g*_B_n*.npz"))[:a.n]
    print(f"  banked flats: {len(flats)}", flush=True)
    rows = []
    for i, f in enumerate(flats):
        Z = np.load(f, allow_pickle=True)
        R = make_realisation(P, "B", int(Z["dist_seed"]), int(Z["noise_seed"]))
        t0 = time.time()
        LNL = G.sweep(sp, R["theta_base"], R["data"], R["mask"], AI, verbose=False)
        red = G.reduce_all(LNL)
        rec = dict(file=os.path.basename(f))
        for nm in ("CENTRE", "PROFILE", "MARGINAL"):
            q, qmax, mapk = G.qpost(red[nm])
            pt = G.ptrue(q, R["n_true"])
            on_true = (mapk == R["n_true"])
            gap = G.gap(red[nm])
            rec[nm] = dict(qmax=qmax, mapk=mapk, ptrue=pt, on_true=on_true, gap=gap)
        rows.append((rec, R["n_true"], R["u_true"]))
        c, m = rec["CENTRE"], rec["MARGINAL"]
        print(f"  [{i}] {os.path.basename(f)} ({time.time()-t0:.0f}s): "
              f"median dq_max {np.median(m['qmax']-c['qmax']):+.4f}; "
              f"median dP_true {np.median(m['ptrue']-c['ptrue']):+.4f}; "
              f"on_true {int(c['on_true'].sum())} -> {int(m['on_true'].sum())}; "
              f"median dgap {np.median(m['gap']-c['gap']):+.3f} nat", flush=True)
        np.savez(f"{RES}/ugrid_impact.npz",
                 **{f"{nm}_{k}": np.array([r[0][nm][k] for r in rows])
                    for nm in ("CENTRE", "PROFILE", "MARGINAL")
                    for k in ("qmax", "mapk", "ptrue", "on_true", "gap")},
                 n_true=np.array([r[1] for r in rows]),
                 u_true=np.array([r[2] for r in rows]),
                 files=np.array([r[0]["file"] for r in rows]))
    print(f"\n  saved {RES}/ugrid_impact.npz", flush=True)
    return 0


def gumbel_floor(offenders, alpha=0.05):
    """D2's estimator: Gumbel-MLE (1-alpha) quantile. floor = mu + beta*z; err = 2.80*beta/sqrt(N)."""
    x = np.asarray([o for o in offenders if np.isfinite(o)], float)
    if len(x) < 3:
        return np.nan, np.nan, np.nan, np.nan
    # method-of-moments seed then a few Newton steps on the Gumbel MLE
    beta = np.std(x) * np.sqrt(6) / np.pi
    for _ in range(200):
        e = np.exp(-x / max(beta, 1e-12))
        beta_new = np.mean(x) - np.sum(x * e) / max(np.sum(e), 1e-300)
        if not np.isfinite(beta_new) or beta_new <= 0:
            break
        if abs(beta_new - beta) < 1e-10:
            beta = beta_new
            break
        beta = beta_new
    mu = -beta * np.log(np.mean(np.exp(-x / max(beta, 1e-12))))
    z = 2.9702                                    # (1 - 0.05) Gumbel quantile factor
    floor = mu + beta * z
    err = 2.80 * beta / np.sqrt(len(x))
    return float(floor), float(err), float(mu), float(beta)


def mode_floors(a):
    """(ii)+(iii): refit the floor from NULLS under each reduction, then re-cut the counts."""
    print(f"=== U-GRID FLOORS + criterion-v2 CUT: {a.nulls} nulls, {a.n} signal reals ===",
          flush=True)
    P = C.build_b1_problem()
    sp = P["sp"]; AI = P["AI"]
    G = UGrid(P, "lit", M=a.M)
    lnK = np.log(np.maximum(np.array(
        [2 * int(np.floor(3.0 * G.sig[p] / G.dL[p] + 1e-9)) for p in range(P["npsr"])]), 2.0))
    RED = ("CENTRE", "PROFILE", "MARGINAL")

    def score(theta_base, data, mask):
        LNL = G.sweep(sp, theta_base, data, mask, AI, verbose=False)
        red = G.reduce_all(LNL)
        out = {}
        for nm in RED:
            q, qmax, mapk = G.qpost(red[nm])
            out[nm] = dict(qmax=qmax, mapk=mapk, gap=G.gap(red[nm]))
        return out

    # ---------- NULLS: pure noise, NO CW in the data (the fN family) ----------
    off = {nm: [] for nm in RED}
    for i in range(a.nulls):
        seed = 950000 + i
        nz = P["nd"].draw(seed=seed)
        data = tuple(jnp.asarray(np.asarray(n)) for n in nz)          # NO CW
        theta_base = np.concatenate([P["L0"], P["theta_truth"][P["n_dist"]:]])
        s = score(theta_base, data, P["mask_one"])
        for nm in RED:
            m = (s[nm]["gap"] > lnK) & (s[nm]["qmax"] > 0.9)          # layers 1 AND 3
            off[nm].append(float(s[nm]["gap"][m].max()) if m.any() else 0.0)
        if i % 5 == 0:
            print(f"  null {i}: offenders " +
                  ", ".join(f"{nm} {off[nm][-1]:.2f}" for nm in RED), flush=True)
    floors = {}
    print(f"\n-- FLOORS (Gumbel-MLE, alpha = 0.05, N = {a.nulls}) --", flush=True)
    for nm in RED:
        f, e, mu, b = gumbel_floor(off[nm])
        floors[nm] = (f, e)
        nz_frac = float(np.mean(np.asarray(off[nm]) == 0.0))
        print(f"   {nm:>9s}: floor = {f:7.2f} +- {e:5.2f} nat  (mu {mu:6.2f}, beta {b:5.2f}; "
              f"zero-offender fraction {nz_frac:.2f})", flush=True)

    # ---------- SIGNAL: criterion-v2 cut under each reduction's OWN floor ----------
    flats = sorted(glob.glob(f"{REPO}/reports/flat_b10_g*_B_n*.npz"))[:a.n]
    tot = {nm: dict(det=0, cert=0, corr=0, wrong=0) for nm in RED}
    for f in flats:
        Z = np.load(f, allow_pickle=True)
        R = make_realisation(P, "B", int(Z["dist_seed"]), int(Z["noise_seed"]))
        s = score(R["theta_base"], R["data"], R["mask"])
        for nm in RED:
            det = s[nm]["gap"] > np.maximum(lnK, floors[nm][0])
            cert = det & (s[nm]["qmax"] > 0.9)
            on_t = s[nm]["mapk"] == R["n_true"]
            tot[nm]["det"] += int(det.sum()); tot[nm]["cert"] += int(cert.sum())
            tot[nm]["corr"] += int((cert & on_t).sum()); tot[nm]["wrong"] += int((cert & ~on_t).sum())
    print(f"\n-- criterion-v2 CUT ({a.n} Arm-B realisations, each reduction under ITS OWN floor) --",
          flush=True)
    print(f"   {'reduction':>9s} {'floor':>8s} {'det':>5s} {'cert':>5s} {'correct':>8s} "
          f"{'wrong':>6s} {'certs/real':>11s}", flush=True)
    for nm in RED:
        t = tot[nm]
        print(f"   {nm:>9s} {floors[nm][0]:8.2f} {t['det']:5d} {t['cert']:5d} {t['corr']:8d} "
              f"{t['wrong']:6d} {t['cert']/a.n:11.2f}", flush=True)
    np.savez(f"{RES}/ugrid_floors.npz",
             **{f"off_{nm}": np.array(off[nm]) for nm in RED},
             **{f"floor_{nm}": np.array(floors[nm]) for nm in RED},
             **{f"tot_{nm}": np.array([tot[nm][k] for k in ("det", "cert", "corr", "wrong")])
                for nm in RED},
             nulls=a.nulls, n=a.n)
    print(f"\n  saved {RES}/ugrid_floors.npz", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "impact", "floors"])
    ap.add_argument("--nulls", type=int, default=30)
    ap.add_argument("--M", type=int, default=16)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--flat", default="900001")
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    return {"gate": mode_gate, "impact": mode_impact, "floors": mode_floors}[a.mode](a)


if __name__ == "__main__":
    sys.exit(main())
