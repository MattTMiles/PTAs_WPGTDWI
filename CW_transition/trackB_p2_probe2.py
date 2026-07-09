"""Track B — P2 PROBE 2 (decisive): the ANCHOR-FRACTION / loop-gain sweep.

PROBE 1 showed the flat P1 plateau is a GAUGE CONSPIRACY (simultaneous re-registration
of all 116 combs), not absent physics: with the array fully anchored (others at truth) the
source-sky gradient is steep and QTRUE spans the gap smoothly. The decisive question: how
MANY anchors break the conspiracy enough to leave a FOLLOWABLE gradient on the surface the
real algorithm sees?

Objective = the ESTIMATOR's surface, self-consistent E-step (NOT others-at-truth):
  * hold k pulsars at TRUTH-anchored fringes (k=0 fringe, distance = L0);
  * the remaining 116-k FREE pulsars snap to their MAP fringe with all others (anchors at
    truth + free at their current MAP) held fixed -- iterated `PASSES` times to self-
    consistency;
  * HARD_k(theta) = fisher_logL(theta, anchors=truth, free=MAP)   [P1's all-snapped surface
    at k=0, which is FLAT];
  * QTRUE_k(theta) = sum over FREE pulsars of q_p(k=0)  (posterior weight on the true fringe).

Sweep k in {0,1,3,6,12,24} along the radial gwphi line (2 deg -> cusp). Anchor ordering =
by true-fringe registration STRENGTH at truth (the pulsars that carry the needle). Plus an
explicit k=3 = CENSUS SET (J0711/J1713/J1909) variant -- the real pipeline's anchors are
prior/data-concentrated, not arbitrary.

Deliverable: the minimum anchor set that restores a followable gradient across the gap ->
"loop closes at k=..". Pre-registered rule:
  * loop closes at k <= 3 census-class -> P2 = soft-EM alternation (build + gate).
  * needs k >> 3 / high-q anchors priors can't supply -> gap physical -> stall-confirm FAIL
    + RTK/LAMBDA integer solve.

Run (jug-gpu, background): python trackB_p2_probe2.py
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
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
from trackB_p1_map import batched_scan, _local_sky_metric

CWT = "/home/mattm/projects/HSYMT/CW_transition"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]
PASSES = 3                          # self-consistent E-step passes
OFFSETS = np.array([0.0, 0.006, 0.012, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0])   # deg, +gwphi
K_LIST = [0, 1, 3, 6, 12, 24]


def _lse(x):
    m = np.max(x); return m + np.log(np.sum(np.exp(x - m)))


def estep_anchored(P, theta_src, anchor_mask, lit):
    """Self-consistent anchored E-step. anchor_mask (npsr,) bool: those held at truth fringe.
    Returns HARD (joint lnL at anchors=truth, free=MAP), QTRUE_free (sum_free q_p(k=0)),
    and soft_free (sum_free logsumexp over fringes -- the marginal surface)."""
    tt = P["theta_truth"]; ndist = P["n_dist"]; npsr = P["npsr"]
    fl = P["amo_full"]["fisher_logL"]; data = P["data_obs"]
    dist = tt[:ndist].copy()                      # start all at truth
    free = np.where(~anchor_mask)[0]
    post = None
    for _ in range(PASSES):
        th = theta_src.copy(); th[:ndist] = dist
        LNL, _ = batched_scan(P, [th])            # (1,npsr,B): each pulsar swept, others at dist
        post = TE.fringe_posterior(P, LNL[0], None, lit, 1.0)
        newd = dist.copy()
        newd[free] = post["map_evalL"][free]      # free -> MAP
        newd[anchor_mask] = tt[:ndist][anchor_mask]   # anchors -> truth
        if np.max(np.abs(newd - dist)) < 1e-9:
            dist = newd; break
        dist = newd
    th = theta_src.copy(); th[:ndist] = dist
    hard = float(fl(jnp.asarray(th), data))
    # QTRUE + soft over FREE pulsars only, from the last posterior
    qtrue = 0.0; soft = 0.0
    for p in free:
        uk, offs_u, lnL_u, w = post["qlist"][p]
        sig = lit[p]
        lp = lnL_u - offs_u ** 2 / (2.0 * sig ** 2)
        soft += float(_lse(lp))
        qtrue += float(w[uk == 0].sum()) if (uk == 0).any() else 0.0
    return hard, qtrue, soft


def anchor_order_by_registration(P, lit):
    """Order pulsars by true-fringe posterior weight q_p(k=0) at TRUTH (strongest first)."""
    LNL, _ = batched_scan(P, [P["theta_truth"]])
    post = TE.fringe_posterior(P, LNL[0], None, lit, 1.0)
    q0 = np.array([float(post["qlist"][p][3][post["qlist"][p][0] == 0].sum())
                   if (post["qlist"][p][0] == 0).any() else 0.0 for p in range(P["npsr"])])
    return np.argsort(-q0), q0


def gwphi_theta(P, off_deg):
    tt = P["theta_truth"].copy(); ndist = P["n_dist"]
    cg0 = tt[ndist + 0]; _, m_gp = _local_sky_metric(cg0)
    tt[ndist + 1] = tt[ndist + 1] + np.radians(off_deg) / m_gp
    return tt


def run_sweep(P, lit, order, anchor_sets):
    """anchor_sets: dict label -> boolean mask (npsr,). Sweep offsets, report HARD & QTRUE."""
    npsr = P["npsr"]
    res = {}
    for label, mask in anchor_sets.items():
        H = np.zeros(len(OFFSETS)); Q = np.zeros(len(OFFSETS)); S = np.zeros(len(OFFSETS))
        t0 = time.time()
        for i, off in enumerate(OFFSETS):
            H[i], Q[i], S[i] = estep_anchored(P, gwphi_theta(P, off), mask, lit)
        H = H - H[0]                               # relative to truth (off=0)
        res[label] = dict(H=H, Q=Q, S=S - S[0])
        # gap-gradient diagnostic: HARD drop 0.05->2 deg, QTRUE rise 2->0.05 deg
        i05 = np.argmin(np.abs(OFFSETS - 0.05)); i2 = np.argmin(np.abs(OFFSETS - 2.0))
        gap_drop = H[i05] - H[i2]                  # >0 & monotone -> gradient across gap
        q_rise = Q[i05] - Q[i2]
        mono = all(H[j] >= H[j + 1] - 1e-6 for j in range(len(OFFSETS) - 1))   # monotone outward
        k = int(mask.sum())
        print(f"\n  [{label}] k={k}: HARD gap-drop(0.05->2deg)={gap_drop:.1f} nat, "
              f"QTRUE rise(2->0.05deg)={q_rise:.2f}, monotone={mono}, {time.time()-t0:.0f}s", flush=True)
        print(f"    {'off_deg':>8s} {'HARD-truth':>12s} {'QTRUE_free':>11s}", flush=True)
        for i, off in enumerate(OFFSETS):
            print(f"    {off:8.3f} {H[i]:12.1f} {Q[i]:11.2f}", flush=True)
    return res


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop")
    lit = P["priors"]["lit"]
    npsr = P["npsr"]
    print(f"  build_problem {time.time()-t0:.0f}s", flush=True)

    order, q0 = anchor_order_by_registration(P, lit)
    census_idx = np.array([P["names"].index(n) for n in CENSUS])
    print(f"  top-8 registrars (idx:q0): "
          + " ".join(f"{int(i)}:{q0[i]:.2f}" for i in order[:8]), flush=True)
    print(f"  census idx {census_idx.tolist()} q0 "
          + " ".join(f"{q0[i]:.2f}" for i in census_idx), flush=True)

    anchor_sets = {}
    for k in K_LIST:
        m = np.zeros(npsr, bool)
        if k > 0:
            m[order[:k]] = True
        anchor_sets[f"reg_top{k}"] = m
    # explicit census-3 variant
    mc = np.zeros(npsr, bool); mc[census_idx] = True
    anchor_sets["census3"] = mc

    res = run_sweep(P, lit, order, anchor_sets)

    # loop-closure summary
    print("\n=== LOOP-GAIN SUMMARY (does the anchored surface give a followable gap gradient?) ===", flush=True)
    for label, d in res.items():
        i05 = np.argmin(np.abs(OFFSETS - 0.05)); i2 = np.argmin(np.abs(OFFSETS - 2.0))
        print(f"  {label:12s}: HARD gap-drop {d['H'][i05]-d['H'][i2]:9.1f} nat | "
              f"QTRUE {d['Q'][i2]:.2f}(2deg)->{d['Q'][i05]:.2f}(0.05deg)", flush=True)
    np.savez(f"{CWT}/trackB_p2probe2_anchorsweep.npz",
             offsets_deg=OFFSETS, q0=q0, order=order, census_idx=census_idx,
             **{f"{lab}_H": d["H"] for lab, d in res.items()},
             **{f"{lab}_Q": d["Q"] for lab, d in res.items()})
    print(f"  saved trackB_p2probe2_anchorsweep.npz ; total {time.time()-t0:.0f}s", flush=True)
