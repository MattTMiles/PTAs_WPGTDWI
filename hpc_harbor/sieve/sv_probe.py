#!/usr/bin/env python3
"""SIEVE-A cost probe: what does a GLACIER venue build cost on the CPU lane, and can
an EVEN/ODD half-venue be built from the same pulsars?  No banks, no verdicts.

Times, in order:
  1. load_pulsars(116) + the T=30 extension
  2. build_glacier_problem at the cascade cell's census (ncw = 256 + harmonic slots)
  3. one logL_batch_theta evaluation (2 rows)
  4. the same build with the pulsars SUBSET to even TOA epochs (the T2 half-venue)

The half-venue is the object T2 needs: an exact held-out likelihood is the full data
vector restricted to the held-out TOA indices, evaluated under a venue built on those
same indices.  This probe only checks that the subset survives the build path and what
it costs; it evaluates no cross-validation statistic.
"""
import os, sys, time
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
          "hpc_harbor/atlas", "hpc_harbor/spark", "hpc_harbor/glacier"):
    sys.path.insert(0, f"{HSYMT}/{p}")

RUNG, ARM, SKY, T = "r13p25", "e07", 0, 30


def main():
    t0 = time.time()
    import glacier_loop as GL
    import glacier_pop as GP
    print(f"[probe] imports {time.time()-t0:.1f}s", flush=True)

    t = time.time()
    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()
    print(f"[probe] _stack {time.time()-t:.1f}s", flush=True)

    t = time.time()
    pop, cond = GP.draw_population_conditional(SKY, RUNG, n_src=256,
                                               band_log10f=GL.BAND_CAMPAIGN)
    slots, n_harm, active, chan, n_clip = GL.embed_igniter(
        pop, GL.E_IGNITER[ARM], GL.VENUE_SPAN_S[T])
    print(f"[probe] census+embed {time.time()-t:.1f}s  n_slot={len(slots)} "
          f"(+{n_harm} harmonic)", flush=True)

    pop_s = dict(pop); pop_s["src"] = slots; pop_s["n_src"] = len(slots)
    t = time.time()
    G = GP.build_glacier_problem(T, pop_s, verbose=True)
    t_build = time.time() - t
    print(f"[probe] build_glacier_problem FULL {t_build:.1f}s", flush=True)

    amo = G["amo"]; nd = amo["n_dist"]
    ones = jnp.ones(amo["npsr"])
    ntoa = [len(p.toas) for p in G["disco_psrs"]]
    print(f"[probe] npsr={amo['npsr']} n_dist={nd} n_theta={amo['n_theta']} "
          f"TOAs total={sum(ntoa)} min/med/max={min(ntoa)}/{int(np.median(ntoa))}/"
          f"{max(ntoa)}", flush=True)

    theta = np.asarray(amo["theta_truth"], float)
    zero = tuple(jnp.zeros(n) for n in ntoa)
    t = time.time()
    ll = np.asarray(amo["logL_batch_theta"](jnp.asarray(np.stack([theta, theta])),
                                            zero, ones))
    print(f"[probe] logL_batch_theta(2) first call (incl. compile) "
          f"{time.time()-t:.1f}s  -> {ll}", flush=True)
    t = time.time()
    _ = np.asarray(amo["logL_batch_theta"](jnp.asarray(np.stack([theta, theta])),
                                           zero, ones))
    print(f"[probe] logL_batch_theta(2) warm {time.time()-t:.2f}s", flush=True)

    # ---- the T2 half-venue: same pulsars, EVEN TOA epochs only ----
    print("[probe] --- half-venue (even TOA epochs) ---", flush=True)
    t = time.time()
    try:
        import discovery as ds
        from cw_helpers import load_pulsars
        ent2, dis2 = load_pulsars(116)
        for pe in dis2:
            IG.extend_pulsar(pe, float(T - 15))
        masks = []
        for pe in dis2:
            n = len(pe.toas)
            m = np.zeros(n, bool); m[::2] = True
            masks.append(m)
        print(f"[probe]   reload+extend {time.time()-t:.1f}s; would keep "
              f"{sum(m.sum() for m in masks)} of {sum(len(p.toas) for p in dis2)} TOAs",
              flush=True)
        # what attributes carry per-TOA length? report, do not mutate blindly
        pe = dis2[0]
        n = len(pe.toas)
        per_toa = sorted(k for k, v in vars(pe).items()
                         if isinstance(v, np.ndarray) and v.shape and v.shape[0] == n)
        print(f"[probe]   per-TOA attributes on disco pulsar: {per_toa}", flush=True)
        pen = ent2[0]
        ne = len(pen.toas)
        per_toa_e = sorted(k for k, v in vars(pen).items()
                           if isinstance(v, np.ndarray) and v.shape
                           and v.shape[0] == ne)
        print(f"[probe]   per-TOA attributes on enterprise pulsar "
              f"({type(pen).__name__}): {per_toa_e}", flush=True)
    except Exception as ex:
        print(f"[probe]   half-venue recon FAILED: {type(ex).__name__}: {ex}",
              flush=True)

    print(f"[probe] TOTAL {time.time()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
