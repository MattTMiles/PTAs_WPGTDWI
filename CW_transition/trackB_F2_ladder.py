"""F2 — the lever-arm ladder spectrum. Per (pulsar, loud-source) registration tolerance in
SCALED source units. A source-sky error d(theta_s) shifts the pulsar-term phase
Phi_p = 2*pi*f_gw*L_p*(1-cos mu_ps)/c by dPhi ~ (dPhi/dtheta_s) d(theta_s); the fringe
re-registers when |dPhi| ~ 1 rad. So the per-baseline sky tolerance (rad) ~ 1 / |grad_theta Phi_p|,
and in SCALED units we divide by the sky param scale. Frequency analogue: dPhi/df = 2*pi*L_p*
(1-cos mu)/c ... plus the Earth-term T*... ~ 2*pi*(L_p(1-cos mu)/c) for the pterm; tol_f(scaled)
= 1/(2*pi*f*L_p(1-cos mu)/c ... ) mapped through the log10_fgw scale.

QUESTION: does the sorted tolerance spectrum span CONTINUOUSLY from ~0.05 scaled (de-risk float
floor) down to ~1e-4 (L2c pull-in), with no gap wider than ~one rung? That is the wide-lane
cascade feasibility criterion. Pure geometry -- no likelihood, no GPU.
Run: python trackB_F2_ladder.py
"""
import os
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
import sys
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
CWT = "/home/mattm/projects/HSYMT/CW_transition"
C_KPC_S = 1.0293e11   # speed of light in kpc/s  (c = 9.716e-12 kpc/s^-1 -> 1/c below); use SI-ish


def main():
    print("=== F2 lever-arm ladder (per pulsar x loud source registration tolerance) ===", flush=True)
    P = TE.build_problem("pop", build_earth=False, verbose=False)
    ndist = P["n_dist"]; names = P["names"]; npsr = P["npsr"]
    src_t = np.asarray(P["theta_truth"][ndist:]); scale = TE.phi_scale(P)
    n_loud = TE.CONFIGS["pop"]["population"][0]
    L0 = np.asarray(P["L0"])                       # kpc
    disco = P["disco_psrs"]
    # pulsar unit vectors
    ppos = np.array([disco[a].pos for a in range(npsr)])   # (npsr,3)
    import jax.numpy as jnp

    def phi_of(theta_s_loud, a):
        """pulsar-term phase Phi_p(source) = 2*pi f L_p (1 - cos mu)/c ; return as fn for autograd."""
        cosgt, gwphi, lfgw = theta_s_loud[0], theta_s_loud[1], theta_s_loud[4]
        st = jnp.sqrt(jnp.clip(1 - cosgt**2, 0, 1))
        khat = jnp.array([st*jnp.cos(gwphi), st*jnp.sin(gwphi), cosgt])   # source sky unit (GW prop ~ -khat)
        cosmu = jnp.dot(khat, ppos[a])
        f = 10.0**lfgw
        c_kpc_per_s = 9.71561e-12                                        # kpc / s  (c)
        L_s = L0[a] / c_kpc_per_s                                        # light-travel time (s)
        return 2*np.pi * f * L_s * (1 - cosmu)                           # radians

    # gradient of Phi wrt the 8 loud params (only 0,1,4 matter); tolerance = 1/|d Phi/d p| in that param,
    # then scaled. Take the sky tolerance = 1/|grad over (cosgt,gwphi)| combined, and freq separately.
    grad_phi = jax.grad(lambda ts, a: phi_of(ts, a), argnums=0)
    rows = []
    for k in range(n_loud):
        ts = jnp.asarray(src_t[8*k:8*k+8])
        for a in range(npsr):
            g = np.asarray(grad_phi(ts, a))                             # d Phi / d param (rad per raw unit)
            # sky: combine cosgt,gwphi gradient magnitude in SCALED units (raw*scale per scaled unit)
            gs = np.sqrt((g[0]*scale[0])**2 + (g[1]*scale[1])**2)       # rad per scaled-sky
            gf = abs(g[4]*scale[4])                                     # rad per scaled-freq
            tol_sky = 1.0/gs if gs > 0 else np.inf                      # scaled units for 1 rad shift
            tol_frq = 1.0/gf if gf > 0 else np.inf
            rows.append((k, a, tol_sky, tol_frq))
    rows = np.array([(r[0], r[1], r[2], r[3]) for r in rows])
    tol_sky = rows[:, 2]; tol_frq = rows[:, 3]
    finite = np.isfinite(tol_sky)
    ts_sorted = np.sort(tol_sky[finite])[::-1]                          # loosest (largest tol) first
    print(f"  {int(finite.sum())}/{len(rows)} finite sky tolerances (scaled units for 1-rad reg shift)")
    print(f"  sky tol spectrum (scaled): max(loosest)={ts_sorted[0]:.3e} median={np.median(ts_sorted):.3e} "
          f"min(tightest)={ts_sorted[-1]:.3e}")
    for thr, lab in [(5e-2,'>=0.05 (float floor)'), (1e-2,'>=1e-2'), (1e-3,'>=1e-3'), (1e-4,'>=1e-4 (pull-in)')]:
        print(f"    # (psr,src) pairs with sky tol {lab}: {int(np.sum(tol_sky>=thr))}")
    # gap analysis on the sorted spectrum between 0.05 and 1e-4
    band = ts_sorted[(ts_sorted <= 5e-2) & (ts_sorted >= 1e-4)]
    if len(band) > 1:
        lr = np.log10(band)
        gaps = -np.diff(lr)                                            # dex gaps descending
        print(f"  in-band [1e-4, 0.05]: {len(band)} rungs; max log10 gap = {gaps.max():.3f} dex "
              f"at ~{10**lr[gaps.argmax()]:.2e} scaled")
        spans = band.min() <= 1.5e-4 and band.max() >= 4e-2 and gaps.max() < 0.7
        print(f"  LADDER SPANS 0.05 -> 1e-4 with no >0.7dex gap: {spans}")
    else:
        print("  in-band rungs < 2 -> ladder does NOT span")
    np.savez(f"{CWT}/trackB_F2_ladder.npz", rows=rows, tol_sky=tol_sky, tol_frq=tol_frq, names=np.array(names))
    print("  saved trackB_F2_ladder.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
