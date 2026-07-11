"""Track B — B1 STEP 1D: the SOFT-CASCADE PROBE. The last door.

WHY THIS EXISTS. STEP 1 measured the registration radius R_a = the box-shrink factor at which
pulsar a's pulsar-term phase wander falls to 1 rad, and found ZERO pulsars with R >= 1 at every
conditioning tier (A 0.0214, B 0.0220, C 0.271). That is the HARD criterion: it asks when a
fringe can be LOCKED. It does not ask what a SOFT, posterior-weighted fringe mixture leaks.
At R ~ 0.27 a soft fix spreads weight over ~1/R ~ 4 fringes rather than locking one, and a
4-fringe mixture still carries chirp-mass information. Whether that leak compounds is a
measurable question, and it is the only door left before the targeted-scenario STOP is written.

THE LOOP (pre-registered, Matt 2026-07-09; bounded, no tuning):
  Tier-C conditioning, Asimov. Start from a float draw inside Tier C's box. Then iterate <= 5x:
    1. SOFT E-step over ALL 116 pulsars -- q_p(n) from the B_CERT=512 fringe posterior.
       NO hard fixes anywhere, ever. Every pulsar stays a mixture.
    2. Mc-marginal update FROM the posterior-weighted fringe mixture. The mixture-marginalised
       objective IS `lnL_marg` (trackB_b1.B1Marg): its curvature over the 6 (f, mc) dims, with
       the tier's own Gaussian prior added, gives sigma(log10_mc) per loud source.
    3. Damped Newton step on lnL_marg, clipped to the tier box. Re-E-step.

REPORTED PER ITERATION (all three mandated):
  * sigma(log10_mc) per loud source (dex).
  * S = sum_p [ max_n q_p(n) ] - S_prior, the aggregate SOFT-REGISTRATION signal, where S_prior
    is the same sum under the fringe PRIOR alone (no likelihood). S > 0 = the data are
    concentrating fringe posteriors beyond what the EM distance priors already do.
  * The WRONG-FIX COLUMN ANALOGUE: W = sum_p [1 - q_p(n_true)] = total posterior mass on FALSE
    fringes, and dW_grew = the number of pulsars whose false-fringe mass GREW since the last
    iteration. A leak that works must shrink sigma_mc WITHOUT growing W: a loop that sharpens
    by concentrating on wrong fringes is the GPS wrong-fix failure, softened but not cured.

PRE-REGISTERED VERDICT (no mid-run softening):
  PASS: cumulative sigma_mc shrinkage >= 1.3x by iter 5 AND monotone -> the loop leaks.
        Extrapolate iterations-to-9.5x (the Tier-C deficit factor) and report the price.
  FAIL: shrinkage < 1.1x OR non-monotone -> the door is checked and closed.
  (1.1x-1.3x monotone: reported as INCONCLUSIVE-WEAK, no PASS claim.)

Truth use: n_true = 0 (Asimov; the injection puts every true distance at its EM prior mean) is
used ONLY to SCORE W, never to steer the loop. The start point and every update are computed
from the data + tier conditioning alone.

Env: jug-gpu, autotune off, x64. Matt commits; never commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_b1_core as C
import trackB_b1 as B1
import trackB_b1_referendum as REF
import trackB_b1_ladder as LAD
import trackB_estimator as TE

CWT = "/home/mattm/projects/HSYMT/CW_transition"
I_MC, I_FGW = 3, 4
DEFICIT_TARGET = 9.5          # Tier-C break-even deficit (STEP 2): sigma_mc must improve this much


def prior_max_sum(FT):
    """S_prior = sum_p max_n q_p(n) under the fringe PRIOR alone (likelihood flat)."""
    tot = 0.0
    for a in range(FT.npsr):
        lp = FT.lnprior[a]
        w = np.exp(lp - lp.max()); w /= w.sum()
        tot += float(w.max())
    return tot


def q_diagnostics(FT, LNL, n_true):
    """max_n q_p, and false-fringe mass 1 - q_p(n_true), per pulsar."""
    qmax = np.zeros(FT.npsr); qfalse = np.zeros(FT.npsr)
    for a in range(FT.npsr):
        pk = np.maximum.reduceat(LNL[a][FT.order[a]], FT.redidx[a])
        logw = (pk - pk.max()) + FT.lnprior[a]
        w = np.exp(logw - logw.max()); w /= w.sum()
        qmax[a] = w.max()
        hit = np.where(FT.uk[a] == n_true[a])[0]
        qfalse[a] = 1.0 - (w[hit[0]] if len(hit) else 0.0)
    return qmax, qfalse


def sigma_profile(G, x, i, sig_prior_i, ngrid=33, target=0.5):
    """Conditional 1-sigma half-width of the LOG-POSTERIOR along active dim i, by profiling.

    NOT a curvature. The Hessian of lnL_marg at the 3e-6 razor scale measures the micro-structure
    of the shell, not the marginal width over the tier prior's scale (sigma_mc ~ 0.073 scaled), and
    it is not positive-definite (the non-positive eigendirection), so inv(Hpos + Pi) returns
    NEGATIVE variances -> sigma = 0 -> nan ratios. Profiling is well-defined regardless: scan the
    log-posterior along e_i, find where it falls `target` nat (0.5 = 1 sigma for a Gaussian), and
    interpolate the crossing in log |s|. Returns the mean of the +/- half-widths, in SCALED units.
    """
    g0 = float(G(x[None])[0])

    def logpost(S):
        X = np.tile(x, (len(S), 1)); X[:, i] = x[i] + S
        gl = G(X)
        pri = -0.5 * ((x[i] + S) ** 2 - x[i] ** 2) / sig_prior_i ** 2
        return gl + pri - g0

    hw = []
    for sign in (-1.0, +1.0):
        S = sign * np.geomspace(1e-7 * sig_prior_i, 3.0 * sig_prior_i, ngrid)
        d = logpost(S)
        below = np.where(d <= -target)[0]
        if len(below) == 0:
            hw.append(3.0 * sig_prior_i)                    # never falls: prior-limited
            continue
        j = below[0]
        if j == 0:
            hw.append(abs(S[0]))
        else:
            x0, x1 = np.log(abs(S[j - 1])), np.log(abs(S[j]))
            y0, y1 = d[j - 1], d[j]
            hw.append(float(np.exp(x0 + ((-target) - y0) * (x1 - x0) / (y1 - y0))))
    return float(np.mean(hw))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters", type=int, default=5)
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--damp", type=float, default=0.5)
    a = ap.parse_args()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    print("=== B1 STEP 1D: SOFT-CASCADE PROBE (Tier-C conditioning, Asimov) ===", flush=True)

    P = C.build_b1_problem()
    tt = P["theta_truth"]; nd = P["n_dist"]
    mask1 = P["mask_one"]
    data = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(mask1))
    G = REF.TargetedMarg(P, data, mask1)
    G.E.gate_fast_full(G.src_from_x(np.zeros((1, G.D))))

    # ---- Tier-C conditioning: box + Gaussian prior precision on the 6 active dims ----
    sig_f, sig_mc, prov = LAD.tier_sigmas(P, "C")
    print(f"  Tier-C float: {prov}", flush=True)
    sig_tier = np.empty(G.D)
    for k in range(C.N_LOUD):
        sig_tier[REF.active_index(k, I_FGW)] = sig_f[k]
        sig_tier[REF.active_index(k, I_MC)] = sig_mc[k]
    Pi = np.diag(1.0 / sig_tier ** 2)
    boxhalf = 3.0 * sig_tier

    # ---- start point: ONE draw from the tier float (truth-blind; the loop never sees truth) ----
    rng = np.random.default_rng(a.seed)
    x = rng.normal(size=G.D) * sig_tier
    x = np.clip(x, -boxhalf, boxhalf)
    print(f"  start x (scaled, from Tier-C float draw seed {a.seed}): "
          f"{np.array2string(x, precision=4)}", flush=True)

    FT = G.E.FT
    n_true = np.zeros(P["npsr"], int)          # Asimov: truth distance == EM prior mean
    S_prior = prior_max_sum(FT)
    print(f"  S_prior (sum_p max_n q_p under the fringe prior alone) = {S_prior:.3f}", flush=True)

    scale = TE.phi_scale(P)
    hstep = np.full(G.D, 3e-6)                 # for the Newton step only (razor-scale curvature)

    hist = []
    W_prev = None
    sig_mc_hist = []
    for it in range(a.iters + 1):
        t0 = time.time()
        # --- 1. soft E-step over ALL pulsars (no hard fixes) ---
        src = G.src_from_x(x[None])[0]
        LNL, lnref = G.E.estep_batch(G.E.theta_of(src[None]))
        C._finite("softcascade LNL", LNL[0])
        qmax, qfalse = q_diagnostics(FT, LNL[0], n_true)
        S = float(qmax.sum()) - S_prior
        W = float(qfalse.sum())
        grew = -1 if W_prev is None else int(np.sum(qfalse > W_prev + 1e-12))
        # --- 2. Mc-marginal update from the mixture-marginalised objective ---
        g0 = float(G(x[None])[0])
        # sigma(log10_mc) = PROFILE half-width of the log-posterior (lnL_marg + tier prior) along
        # each mc dim. See sigma_profile(): the curvature route is ill-defined here (Hpos is not PD).
        smc = np.array([sigma_profile(G, x, REF.active_index(k, I_MC),
                                      sig_tier[REF.active_index(k, I_MC)])
                        for k in range(C.N_LOUD)]) * scale[I_MC]
        Hpos = REF._fd_hessian(lambda X: G(x[None, :] + np.atleast_2d(X)), g0, hstep, G.D)
        sig_mc_hist.append(smc)
        hist.append(dict(it=it, x=x.copy(), smc=smc, S=S, W=W, grew=grew, g0=g0,
                         qmax_med=float(np.median(qmax))))
        print(f"\n  iter {it}: lnL_marg = {g0:.3f} ({time.time()-t0:.0f}s)", flush=True)
        print(f"    sigma(log10_mc) dex = {np.round(smc,5)}   (median q_max {np.median(qmax):.4f})",
              flush=True)
        print(f"    S = sum_p max_n q_p - S_prior = {S:+.4f}  (soft-registration signal)", flush=True)
        print(f"    W = sum_p [1 - q_p(n_true)] = {W:.4f}  (false-fringe mass); "
              f"pulsars whose W GREW = {grew}", flush=True)
        W_prev = qfalse.copy()
        if it == a.iters:
            break
        # --- 3. damped Newton on the log-POSTERIOR (lnL_marg + tier prior), clipped to the box ---
        A_ = Hpos + Pi
        pts = []
        for i in range(G.D):
            e = np.zeros(G.D); e[i] = hstep[i]
            pts.append(x + e); pts.append(x - e)
        gv = G(np.array(pts))                       # batched: 12 points, one padded call
        grad = np.array([(gv[2 * i] - gv[2 * i + 1]) / (2 * hstep[i]) for i in range(G.D)])
        grad_post = grad - Pi @ x                   # prior is N(0, sig_tier^2) about the float centre
        try:
            step = np.linalg.solve(A_, grad_post)
        except np.linalg.LinAlgError:
            step = grad_post / np.maximum(np.diag(A_), 1e-30)
        x = np.clip(x + a.damp * step, -boxhalf, boxhalf)

    # ---- verdict ----
    smc0 = np.array(sig_mc_hist[0]); smcN = np.array(sig_mc_hist[-1])
    shrink = float(np.median(smc0 / smcN))
    per_iter = np.array([np.median(sig_mc_hist[i] / sig_mc_hist[i + 1])
                         for i in range(len(sig_mc_hist) - 1)])
    monotone = bool(np.all(per_iter >= 1.0 - 1e-9))
    Wser = np.array([h["W"] for h in hist])
    print(f"\n=== SOFT-CASCADE VERDICT ===", flush=True)
    print(f"  sigma_mc median: {np.median(smc0):.5f} -> {np.median(smcN):.5f} dex; "
          f"cumulative shrink {shrink:.3f}x over {a.iters} iters", flush=True)
    print(f"  per-iter shrink factors: {np.round(per_iter,4)}  monotone={monotone}", flush=True)
    print(f"  false-fringe mass W: {np.round(Wser,3)}  (grew at any iter: "
          f"{bool(np.any(np.diff(Wser) > 0))})", flush=True)
    W_grew = bool(np.any(np.diff(Wser) > 0))
    if shrink >= 1.3 and monotone and not W_grew:
        r = float(np.median(per_iter))
        n_need = np.log(DEFICIT_TARGET) / np.log(max(r, 1.0 + 1e-9))
        verdict = (f"PASS -- the loop LEAKS. Extrapolated iterations to close the {DEFICIT_TARGET}x "
                   f"Tier-C deficit at the measured per-iter factor {r:.4f}: N = {n_need:.1f}. "
                   f"DOOR OPEN, price = {n_need:.1f} iterations of soft re-registration.")
    elif shrink < 1.1 or not monotone or W_grew:
        why = []
        if shrink < 1.1: why.append(f"shrinkage {shrink:.3f}x < 1.1x")
        if not monotone: why.append("non-monotone")
        if W_grew: why.append("false-fringe mass W GREW (soft wrong-fix)")
        verdict = ("FAIL -- " + "; ".join(why) + ". The soft mixture does not compound. "
                   "THE DOOR IS CHECKED AND CLOSED.")
    else:
        verdict = ("INCONCLUSIVE-WEAK (1.1x <= shrink < 1.3x, monotone). No PASS claim; "
                   "reported as a marginal leak that does not reach the pre-registered bar.")
    print(f"  {verdict}", flush=True)
    np.savez(f"{CWT}/b1_softcascade.npz",
             sig_mc=np.array(sig_mc_hist), S=np.array([h["S"] for h in hist]), W=Wser,
             grew=np.array([h["grew"] for h in hist]), x=np.array([h["x"] for h in hist]),
             g0=np.array([h["g0"] for h in hist]), shrink=shrink, per_iter=per_iter,
             monotone=monotone, verdict=verdict, seed=a.seed, S_prior=S_prior)
    print(f"  saved b1_softcascade.npz", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
