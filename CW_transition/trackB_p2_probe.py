"""Track B — P2 probe: is the flat 0.05-2 deg GAP a HARD-profiling artifact?

P1 found the profiled joint lnL (per-pulsar MAP-fringe snap) is a smooth funnel > ~2 deg,
a near-flat noisy plateau ~2 deg -> ~0.05 deg, and DEAD FLAT + a 0.006 deg delta cusp
below ~0.05 deg. If a cold-start search must cross that flat gap it stalls -> the needle
is unreachable. HYPOTHESIS (Matt): the plateau is an artifact of the HARD max_n reduction;
the SOFT marginal surface (logsumexp_n, the validated A2 soft E-step) has gradient across
the gap.

Pre-registered diagnostic (decision rule set BEFORE running):
  * SOFT gradient present + coherent across 0.05-2 deg where HARD is flat -> P2 = SOFT-
    surface descent (logsumexp objective), no RTK build needed for search.
  * No soft gradient either -> the gap is PHYSICAL -> FAIL verdict or RTK build.

Three measurements along a log-spaced line from +/-2 deg down to the cusp, both sky axes
(everything else at truth = the P1 conditioning; stated in the writeup):
  1. HARD_surf = sum_p max_n [lnL_p(n) - offs_n^2/2sig^2]      (per-pulsar MAP + EM prior)
     SOFT_surf = sum_p logsumexp_n [same]                       (per-pulsar marginal)
  2. QTRUE     = sum_p q_p(k=0)  (softmax weight on the TRUE fringe; the hidden gradient)
  3. SLICE CHECK: repeat the gwphi line with loud0 log10_fgw polished to its conditional
     best SOFT per trial (5-pt freq line search) -> the sky-only slice isn't hiding joint
     sky-freq gradient.

Run (jug-gpu, background): python trackB_p2_probe.py
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
from trackB_p1_map import batched_scan, _local_sky_metric, ang_deg

CWT = "/home/mattm/projects/HSYMT/CW_transition"


def _lse(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def surfaces(P, LNL, lit):
    """HARD (sum_p max_n), SOFT (sum_p logsumexp_n), QTRUE (sum_p q_p(k=0)) from one
    trial's per-pulsar fringe scan LNL (npsr,B). Uses the A2 EM prior (Gaussian tail) and
    the per-fringe peak reduction, exactly as classify_full / fringe_posterior."""
    post = TE.fringe_posterior(P, LNL, None, lit, 1.0)
    hard = soft = qtrue = 0.0
    for a in range(P["npsr"]):
        uk, offs_u, lnL_u, w = post["qlist"][a]
        sig = lit[a]
        lp = lnL_u - offs_u ** 2 / (2.0 * sig ** 2)     # log unnormalised posterior per fringe
        hard += float(np.max(lp))
        soft += float(_lse(lp))
        qtrue += float(w[uk == 0].sum()) if (uk == 0).any() else 0.0
    return hard, soft, qtrue


# log-spaced signed offsets (deg) spanning the whole gap: 2 deg -> cusp core -> 0
OFF = np.array([2.0, 1.0, 0.5, 0.25, 0.1, 0.05, 0.025, 0.012, 0.006, 0.003, 0.0015, 0.0007])
LINE_OFFSETS = np.concatenate([-OFF[::-1], [0.0], OFF])       # 25 points


def make_line(P, axis):
    """Trials along one sky axis at LINE_OFFSETS (deg) about truth; everything else truth."""
    tt = P["theta_truth"]; ndist = P["n_dist"]
    cg0, gp0 = tt[ndist + 0], tt[ndist + 1]
    m_cg, m_gp = _local_sky_metric(cg0)
    trials = []
    for off in LINE_OFFSETS:
        th = tt.copy()
        if axis == 0:
            th[ndist + 0] = np.clip(cg0 + np.radians(off) / m_cg, -1, 1)
        else:
            th[ndist + 1] = gp0 + np.radians(off) / m_gp
        trials.append(th)
    return trials


def run_line(P, axis, lit, tag):
    trials = make_line(P, axis)
    LNL, t = batched_scan(P, trials)
    H = np.zeros(len(trials)); S = np.zeros(len(trials)); Q = np.zeros(len(trials))
    for i in range(len(trials)):
        H[i], S[i], Q[i] = surfaces(P, LNL[i], lit)
    axname = "cos_gwtheta" if axis == 0 else "gwphi"
    print(f"\n=== {tag} ({axname} line, {t:.0f}s) ===", flush=True)
    print(f"  {'offset_deg':>11s} {'HARD-max':>10s} {'SOFT-max':>10s} {'QTRUE':>8s}", flush=True)
    H0, S0 = H.max(), S.max()
    for i, off in enumerate(LINE_OFFSETS):
        print(f"  {off:11.4f} {H[i]-H0:10.2f} {S[i]-S0:10.2f} {Q[i]:8.2f}", flush=True)
    np.savez(f"{CWT}/trackB_p2probe_{tag}.npz", offsets_deg=LINE_OFFSETS, axis=axis,
             hard=H, soft=S, qtrue=Q)
    return dict(off=LINE_OFFSETS, H=H, S=S, Q=Q)


def run_freq_polished_line(P, axis, lit, tag, nf=5, half_bins=1.0):
    """Slice check: same line but loud0 log10_fgw polished to conditional-best SOFT per
    trial (nf-pt freq line search). Confirms the sky-only slice hides no sky-freq gradient."""
    trials = make_line(P, axis)
    tt = P["theta_truth"]; ndist = P["n_dist"]
    f0 = tt[ndist + 4]
    T = 6.992e8; df = 1.0 / (3.0 * T); dlf = half_bins * df / (10 ** f0 * np.log(10))
    lfs = f0 + np.linspace(-dlf, dlf, nf)
    Sp = np.zeros(len(trials)); Qp = np.zeros(len(trials)); bestf = np.zeros(len(trials))
    t0 = time.time()
    for i, th in enumerate(trials):
        variants = []
        for lf in lfs:
            v = th.copy(); v[ndist + 4] = lf; variants.append(v)
        LNL, _ = batched_scan(P, variants)
        best = (-np.inf, f0, None)
        for k in range(nf):
            _, s, q = surfaces(P, LNL[k], lit)
            if s > best[0]:
                best = (s, lfs[k], q)
        Sp[i] = best[0]; bestf[i] = best[1]; Qp[i] = best[2]
    print(f"\n=== {tag} (freq-polished {('cos_gwtheta' if axis==0 else 'gwphi')} line, "
          f"{time.time()-t0:.0f}s) ===", flush=True)
    print(f"  {'offset_deg':>11s} {'SOFT-max':>10s} {'QTRUE':>8s} {'df_best(bins)':>13s}", flush=True)
    S0 = Sp.max()
    for i, off in enumerate(LINE_OFFSETS):
        dfb = (bestf[i] - f0) * (10 ** f0 * np.log(10)) / df
        print(f"  {off:11.4f} {Sp[i]-S0:10.2f} {Qp[i]:8.2f} {dfb:13.3f}", flush=True)
    np.savez(f"{CWT}/trackB_p2probe_{tag}.npz", offsets_deg=LINE_OFFSETS, axis=axis,
             soft=Sp, qtrue=Qp, best_dfbins=(bestf - f0) * (10 ** f0 * np.log(10)) / df)
    return dict(off=LINE_OFFSETS, S=Sp, Q=Qp)


def verdict(lines):
    """Pre-registered decision: is there SOFT (or QTRUE) gradient across the 0.05-2 deg gap
    where HARD is flat? Measure slope of each on the positive side in the gap window."""
    print("\n=== PRE-REGISTERED DIAGNOSTIC (gap window 0.05-2 deg) ===", flush=True)
    off = LINE_OFFSETS
    gapmask = (off >= 0.05) & (off <= 2.0)          # positive-side gap points
    for tag, d in lines.items():
        if "H" in d:
            dH = d["H"][gapmask].max() - d["H"][gapmask].min()
        else:
            dH = np.nan
        dS = d["S"][gapmask].max() - d["S"][gapmask].min()
        dQ = d["Q"][gapmask].max() - d["Q"][gapmask].min()
        # monotone inward? compare gap-edge (2 deg) to gap-inner (0.05 deg)
        i2 = np.argmin(np.abs(off - 2.0)); i05 = np.argmin(np.abs(off - 0.05))
        soft_rise = d["S"][i05] - d["S"][i2]; q_rise = d["Q"][i05] - d["Q"][i2]
        print(f"  {tag:20s}: HARD span {dH:7.2f} | SOFT span {dS:7.2f} (inward rise "
              f"{soft_rise:+7.2f}) | QTRUE span {dQ:6.2f} (inward rise {q_rise:+6.2f})", flush=True)
    print("\n  RULE: SOFT/QTRUE inward-rise >> HARD span AND monotone -> GAP IS ARTIFACT "
          "(P2 = soft descent). Both flat -> GAP PHYSICAL (FAIL verdict / RTK).", flush=True)


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop")
    lit = P["priors"]["lit"]
    print(f"  build_problem {time.time()-t0:.0f}s", flush=True)
    lines = {}
    lines["cos_line"] = run_line(P, 0, lit, "cos_line")
    lines["gwphi_line"] = run_line(P, 1, lit, "gwphi_line")
    lines["gwphi_freqpolished"] = run_freq_polished_line(P, 1, lit, "gwphi_freqpolished")
    verdict(lines)
    print(f"\n  total {time.time()-t0:.0f}s", flush=True)
