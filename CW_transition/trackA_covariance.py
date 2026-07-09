"""TRACK A — anisotropic distance information via the cross-pulsar GWB covariance.

Referendum on Farr's isotropic zero-information theorem for pulsar distances: in the
confusion regime the background is a finite-population shot-noise sky, anisotropic at low
multipoles. Distances enter the CROSS-pulsar covariance through pulsar-term phase factors
that cancel under the isotropic sky integral but may survive for anisotropic multipoles.

Frequency-domain cross-pulsar covariance INCLUDING pulsar terms with finite distances:
  C_ab(f) = Integral_sky P(f,nhat) [Fp_a Fp_b + Fx_a Fx_b]
                                   (1 - e^{-2 pi i f tau_a})(1 - e^{+2 pi i f tau_b})
  tau_p = (L_p/c)(1 - cos mu_p(nhat)) ;  Fp,Fx from prong2's antenna() (incl. 1/(1-cos mu)).
The polarisation sum Fp_a Fp_b + Fx_a Fx_b is invariant under the common sky-frame rotation
psi, so psi=0 is exact.

F0 = model + isotropic-null GATES:
 (a) isotropic null: P=monopole -> dC_ab/dL_a -> 0 for a!=b (short-wavelength f*tau>>1);
     residual converges down with sky resolution and -> 0 as f*tau grows.
 (b) auto-term: dC_aa/dL_a survives isotropy (nonzero); the comparison baseline.
 (c) quadrature convergence: double resolution, well-resolved integrals stable to 1e-8.

Run (jug-gpu):  python trackA_covariance.py f0
"""
import os
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
import sys, argparse
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
from prong2_transition import antenna, KPC_OVER_C, build_array


# ---------------- sky quadrature ----------------
def sky_grid(n_mu, n_phi):
    """Gauss-Legendre in x=cos(theta) x uniform phi. Returns theta,phi,W (flattened),
    with sum(W) = 4 pi (solid angle)."""
    x, wx = np.polynomial.legendre.leggauss(n_mu)      # nodes in [-1,1], weights sum to 2
    theta = np.arccos(x)                                # (n_mu,)
    phi = (np.arange(n_phi) + 0.5) * (2 * np.pi / n_phi)  # midpoint uniform
    TH, PH = np.meshgrid(theta, phi, indexing="ij")     # (n_mu,n_phi)
    W = np.outer(wx, np.full(n_phi, 2 * np.pi / n_phi))  # (n_mu,n_phi), sum = 2*2pi = 4pi
    return (jnp.asarray(TH.ravel()), jnp.asarray(PH.ravel()), jnp.asarray(W.ravel()))


def pulsar_response(theta, phi, p_hats):
    """Fp,Fx,cos_mu for every (sky point, pulsar). psi=0 (pol-sum invariant).
    Returns Fp,Fx,cos_mu each (n_psr, n_pix)."""
    def one(ph):
        return jax.vmap(lambda t, p: antenna(t, p, 0.0, ph))(theta, phi)  # (n_pix,) each
    Fp, Fx, cmu = jax.vmap(one)(p_hats)                 # (n_psr, n_pix)
    return Fp, Fx, cmu


def covariance(theta, phi, W, p_hats, L_kpc, f, P=None):
    """C_ab(f) (n_psr,n_psr complex) and dC/dL_a for each a (n_psr,n_psr,n_psr) where the
    last axis is d/dL_a (only a=row or a=col nonzero). P: (n_pix,) sky power, default monopole=1.
    Normalised so the monopole sky-average is 1: divide by 4 pi."""
    Fp, Fx, cmu = pulsar_response(theta, phi, p_hats)   # (n_psr,n_pix)
    if P is None:
        P = jnp.ones_like(theta)
    tau = L_kpc[:, None] * KPC_OVER_C * (1.0 - cmu)     # (n_psr,n_pix) seconds
    psi = 2.0 * jnp.pi * f * tau                        # phase (n_psr,n_pix)
    g = 1.0 - jnp.exp(-1j * psi)                        # (n_psr,n_pix)
    dpsi_dL = 2.0 * jnp.pi * f * KPC_OVER_C * (1.0 - cmu)   # d psi_a / d L_a  (n_psr,n_pix)
    dg = 1j * dpsi_dL * jnp.exp(-1j * psi)             # dg_a/dL_a (n_psr,n_pix)
    R = Fp[:, None, :] * Fp[None, :, :] + Fx[:, None, :] * Fx[None, :, :]   # pol sum (n,n,n_pix)
    WP = (W * P) / (4.0 * jnp.pi)                       # (n_pix,)
    C = jnp.sum(WP[None, None, :] * R * g[:, None, :] * jnp.conj(g)[None, :, :], axis=-1)
    # dC_ab/dL_a : replace g_a by dg_a
    dC_da = jnp.sum(WP[None, None, :] * R * dg[:, None, :] * jnp.conj(g)[None, :, :], axis=-1)
    return C, dC_da


# ---------------- F0 gates ----------------
def f0(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--npsr", type=int, default=5)
    ap.add_argument("--L", type=float, default=1.0)             # kpc, all pulsars
    a = ap.parse_args(argv)

    p_hats = build_array(a.npsr, seed=777)
    L = jnp.full(a.npsr, a.L)

    print(f"jax {jax.__version__}, {jax.devices()}; npsr={a.npsr}, L={a.L} kpc\n")

    # f parametrised by ftau0 = f*L*(kpc/c) = phase at (1-cos mu)=1 (max phase = 2*ftau0)
    def f_of(ftau0):
        return ftau0 / (a.L * KPC_OVER_C)

    print("=== GATE (a) isotropic null + (b) auto-term, vs sky resolution ===")
    print("    (P=monopole; report |dC_ab/dL_a| a!=b  vs  |dC_aa/dL_a|)\n")
    for ftau0 in (10.0, 30.0, 100.0):
        f = f_of(ftau0)
        print(f"  ftau0 = {ftau0:6.1f}  (f = {f:.3e} Hz)")
        prev_auto = None
        for n_mu in (200, 400, 800, 1600):
            n_phi = 2 * n_mu
            th, ph, W = sky_grid(n_mu, n_phi)
            C, dC = covariance(th, ph, W, p_hats, L, f)
            adC = np.abs(np.asarray(dC))                            # (n,n): dC_ab/dL_a
            auto = np.diag(adC)                                      # dC_aa/dL_a
            off = adC.copy()
            np.fill_diagonal(off, 0.0)
            cross_max = off.max()                                    # max_{a!=b} |dC_ab/dL_a|
            auto_med = np.median(auto)
            ratio = cross_max / auto_med
            conv = "" if prev_auto is None else f" auto d(2x)={abs(auto_med-prev_auto)/auto_med:.2e}"
            print(f"     n_mu={n_mu:4d} n_phi={n_phi:4d} pix={n_mu*n_phi:>8d} | "
                  f"cross_max={cross_max:.3e}  auto_med={auto_med:.3e}  cross/auto={ratio:.2e}{conv}")
            prev_auto = auto_med
        print()

    print("=== GATE (c) quadrature convergence of C_ab (auto + a well-separated cross) ===")
    ftau0 = 30.0; f = f_of(ftau0)
    Cs = {}
    for n_mu in (400, 800, 1600):
        th, ph, W = sky_grid(n_mu, 2 * n_mu)
        C, _ = covariance(th, ph, W, p_hats, L, f)
        Cs[n_mu] = np.asarray(C)
    for n_mu in (800, 1600):
        prev = Cs[n_mu // 2]; cur = Cs[n_mu]
        d_auto = np.abs(cur[0, 0] - prev[0, 0]) / np.abs(cur[0, 0])
        # a well-separated pair (smallest |cos_mu diff|); just use (0,1)
        d_cross = np.abs(cur[0, 1] - prev[0, 1]) / (np.abs(cur[0, 1]) + 1e-30)
        print(f"  n_mu {n_mu//2}->{n_mu}: rel d|C_00|={d_auto:.2e}  |C_00|={np.abs(cur[0,0]):.4e} ; "
              f"rel d|C_01|={d_cross:.2e}  |C_01|={np.abs(cur[0,1]):.4e}")

    # save a small record
    np.savez("/home/mattm/projects/HSYMT/CW_transition/trackA_f0.npz",
             p_hats=np.asarray(p_hats), L=np.asarray(L))
    print("\n[F0] saved trackA_f0.npz")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("mode", choices=["f0"], nargs="?", default="f0")
    a, rest = ap.parse_known_args()
    sys.exit(f0(rest))
