"""Prong-2 close-out P2-B (v2, corrected): coherence-axis physical grounding.

Coherence Y-axis unit: SNR*sigma_phi (Stage A.2: distance info halves at
SNR*sigma_phi = 1). Pulsar-term lag tau_p = (L/c)(1-cos mu) ~ 1-10 kyr for L~kpc.

CORRECTIONS vs the first pass:
 (1) Pulsar-side phase uncertainty is the NOISE-TO-SIGNAL ratio, NOT a time-base shift.
     A red/DM timing fluctuation sigma_res(f) in the CW band perturbs the pulsar-term
     PHASE by sigma_phi = sigma_res / A_CW, where A_CW = h/(2 pi f) is the CW pulsar-term
     timing amplitude. First pass wrote sigma_phi = 2 pi f sigma_res (treating the noise
     as a time shift), which is a factor 1/A_CW*(1/(2 pi f))^-1 ... = 1/h ~ 5.6e13 (~13.7
     orders) too SMALL. So sigma_phi = sigma_res * 2 pi f / h = 1/SNR_pterm ~ O(0.1-few).
 (2) The deterministic GW chirp accumulates Delta_phi = pi * fdot * tau_p^2 relative to a
     monochromatic extrapolation. It is MODELLED (evolve=True) -> recoverable info, OFF
     the coherence axis; shown as a table across (Mc,f), NOT as a band.
 (3) Environmental drift: report the REQUIRED df/f to cross the knee as a family of
     SNR*sigma_phi(df/f) curves for SNR=5,20,100; astrophysical magnitude left OPEN.

Run in jug-gpu:  python stagec_p2b.py
"""

import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import load_pulsars
from discovery.deterministic import make_phase_connected_binary

YR = 365.25 * 86400.0
TSUN = 4.925490947e-6   # G Msun / c^3 [s]
F0 = 10e-9              # 10 nHz
H_REF = -13.75          # equal-strain reference (D4)
TAU_P_REF = 3e3 * YR    # 3 kyr reference lag

_FULL = jax.vmap(make_phase_connected_binary(pulsarterm=True),
                 in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0, None))
_EARTH = jax.vmap(make_phase_connected_binary(pulsarterm=False),
                  in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0))
_KEYS = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
         "log10_fgw", "log10_h", "phase0", "psi")


def chirp_fdot(f, mc_msun):
    m = mc_msun * TSUN
    return (96.0 / 5.0) * np.pi ** (8.0 / 3.0) * m ** (5.0 / 3.0) * f ** (11.0 / 3.0)


def pterm_amplitude(disco_psrs, src):
    """Per-pulsar pulsar-term timing RMS A_CW [s] for one reference source."""
    st = {k: jnp.array([src[k]]) for k in _KEYS}
    A = []
    for p in disco_psrs:
        ba = (p.toas, p.pos, st["cos_gwtheta"], st["gwphi"], st["cos_inc"],
              st["log10_mc"], st["log10_fgw"], st["log10_h"], st["phase0"], st["psi"])
        r_pterm = np.array(jnp.sum(_FULL(*ba, p.pdist[0]), axis=0)
                           - jnp.sum(_EARTH(*ba), axis=0))
        A.append(float(np.sqrt(np.mean(r_pterm ** 2))))
    return np.array(A)


def white_sigma_res_bin(disco_psrs):
    """Per-pulsar white-noise timing fluctuation in the 1/T bin at 10 nHz:
    sigma_res = sigma_TOA * sqrt(2/N_TOA)  (real toaerrs from the feathers)."""
    return np.array([float(np.sqrt(np.mean(np.asarray(p.toaerrs) ** 2))
                           * np.sqrt(2.0 / len(p.toas))) for p in disco_psrs])


def main():
    print("=== P2-B (v2): coherence-axis physical grounding ===\n")
    ent, disco = load_pulsars(116)
    names = [p.name for p in disco]

    # reference source (equal strain h=-13.75, 10 nHz, mid geometry)
    src = {"cos_gwtheta": 0.3, "gwphi": 2.0, "cos_inc": 0.4, "log10_mc": 9.0,
           "log10_fgw": np.log10(F0), "log10_h": H_REF, "phase0": 1.0, "psi": 0.7}
    A_CW = pterm_amplitude(disco, src)          # s
    sig_res = white_sigma_res_bin(disco)        # s (real white)
    sigma_phi = sig_res / A_CW                   # rad  (= 1/SNR_pterm)
    snr_pterm = A_CW / sig_res

    h = 10.0 ** H_REF
    # ---------- Item 2: full arithmetic chain, 3 named pulsars ----------
    print("--- Item 2 CORRECTED: pulsar-side red/DM -> pulsar-term phase uncertainty ---")
    print("chain: sigma_res(bin) = sigma_TOA*sqrt(2/N) [s] ;  A_CW = RMS pulsar-term [s]")
    print(f"       (A_CW ~ h/(2pi f) * antenna, h=1e{H_REF}, f=10nHz -> h/(2pi f)="
          f"{h/(2*np.pi*F0)*1e9:.1f} ns) ;  sigma_phi = sigma_res/A_CW = 1/SNR_pterm\n")
    named = ["J1713+0747", "J0437-4715", "J1909-3744"]
    print(f"{'pulsar':13s} {'N_TOA':>6s} {'sigTOA_ns':>10s} {'sigres_ns':>10s} "
          f"{'A_CW_ns':>9s} {'sigma_phi':>10s} {'SNR_pterm':>10s}")
    for nm in named:
        a = names.index(nm)
        print(f"{nm:13s} {len(disco[a].toas):>6d} "
              f"{np.sqrt(np.mean(np.asarray(disco[a].toaerrs)**2))*1e9:>10.1f} "
              f"{sig_res[a]*1e9:>10.1f} {A_CW[a]*1e9:>9.1f} "
              f"{sigma_phi[a]:>10.3f} {snr_pterm[a]:>10.3f}")
    print(f"\n  FIRST-PASS ERROR: wrote sigma_phi = 2*pi*f*sigma_res (time-shift) = "
          f"{2*np.pi*F0*np.median(sig_res):.2e} rad;")
    print(f"  correct sigma_phi = sigma_res/A_CW = (2*pi*f/h)*sigma_res = "
          f"{np.median(sigma_phi):.3f} rad -> factor 1/h = {1/h:.2e} (~{np.log10(1/h):.1f} orders) larger.\n")
    print(f"  ARRAY (116): sigma_phi median {np.median(sigma_phi):.2f} rad "
          f"(16/84 {np.percentile(sigma_phi,16):.2f}/{np.percentile(sigma_phi,84):.2f}); "
          f"SNR_pterm median {np.median(snr_pterm):.2f}.")
    print(f"  -> sigma_phi ~ O(1) rad = 1/SNR_pterm: this IS the SNR/measurement floor,\n"
          "     already inside the marginalised-noise Fisher (D2-D5). Red/DM noise is an\n"
          "     SNR-axis effect that sits AT the coherence knee BY CONSTRUCTION "
          "(SNR_pterm*sigma_phi=1);\n     it is not an INDEPENDENT extra decoherence term.\n")

    # ---------- Item 1b: GW chirp accumulated phase (Mc, f) plane ----------
    print("--- Item 1b CORRECTED: deterministic GW chirp Delta_phi = pi*fdot*tau_p^2 (rad), "
          "tau_p=3 kyr ---")
    fs = [3e-9, 10e-9, 30e-9]
    print(f"{'Mc(Msun)':>10s} " + "".join(f"f={f*1e9:.0f}nHz{'':>4s}" for f in fs))
    for logmc in (8.0, 8.5, 9.0, 9.5):
        mc = 10 ** logmc
        row = [np.pi * chirp_fdot(f, mc) * TAU_P_REF ** 2 for f in fs]
        print(f"  1e{logmc:<7.1f} " + "".join(f"{v:>11.2e}" for v in row))
    print("  -> Delta_phi ~ 1e-3 .. 1e4 rad (large), but DETERMINISTIC and MODELLED by\n"
          "     evolve=True -> recoverable info, OFF the coherence axis. Not plotted as a band.\n")

    # ---------- Item 1a: required df/f family ----------
    print("--- Item 1a CORRECTED: REQUIRED stochastic df/f to cross the knee "
          "(SNR*sigma_phi=1) ---")
    print("  sigma_phi_env = 2*pi*(df/f)*f*tau_p ;  knee at df/f = 1/(2*pi*f*tau_p*SNR)")
    print(f"  (f=10nHz, tau_p=3 kyr -> 2*pi*f*tau_p = {2*np.pi*F0*TAU_P_REF:.0f})")
    for snr in (5, 20, 100):
        dff = 1.0 / (2 * np.pi * F0 * TAU_P_REF * snr)
        print(f"    SNR={snr:<4d} -> required df/f = {dff:.2e}")
    print("  astrophysical df/f magnitude over kyr lags: OPEN (handoff).\n")

    np.savez("/home/mattm/projects/HSYMT/CW_transition/stagec_p2b_results.npz",
             names=np.array(names), A_CW=A_CW, sig_res=sig_res,
             sigma_phi=sigma_phi, snr_pterm=snr_pterm)

    # ---------- figure: 2 panels, clipped ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # Panel 1: environmental threshold family SNR*sigma_phi vs df/f
    dff = np.logspace(-7, -2, 200)
    C = 2 * np.pi * F0 * TAU_P_REF
    for snr, col in [(5, "#2c7fb8"), (20, "#238b45"), (100, "#d95f0e")]:
        ax1.plot(dff, snr * C * dff, color=col, label=f"SNR={snr}")
        cross = 1.0 / (snr * C)
        ax1.plot([cross], [1.0], "o", color=col, ms=7)
    ax1.axhline(1.0, color="k", ls="--", lw=1, label="knee (info halves)")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("required stochastic df/f over tau_p=3 kyr")
    ax1.set_ylabel("SNR * sigma_phi")
    ax1.set_ylim(1e-3, 1e2); ax1.set_xlim(1e-7, 1e-2)
    ax1.set_title("Item 1a: required df/f to decohere (astro magnitude OPEN)")
    ax1.legend(fontsize=8)

    # Panel 2: 116 per-pulsar pulsar-side sigma_phi = 1/SNR_pterm
    ax2.hist(np.clip(sigma_phi, 0, 3), bins=30, color="#2c7fb8", alpha=0.8)
    ax2.axvline(1.0, color="k", ls="--", lw=1)
    ax2.text(1.02, ax2.get_ylim()[1]*0.9, "sigma_phi = 1 (knee)", fontsize=8)
    ax2.text(0.5, ax2.get_ylim()[1]*0.6,
             f"median {np.median(sigma_phi):.2f} rad\n= 1/SNR_pterm\n(measurement floor,\nin D2-D5 already)",
             fontsize=8)
    ax2.set_xlabel("pulsar-side sigma_phi = sigma_res/A_CW  (rad)")
    ax2.set_ylabel("# pulsars (of 116)")
    ax2.set_title(f"Item 2: per-pulsar pulsar-term phase uncertainty (h={H_REF})")
    fig.tight_layout()
    fig.savefig("/home/mattm/projects/HSYMT/CW_transition/p2b_coherence_grounding.png", dpi=130)
    print("saved p2b_coherence_grounding.png, stagec_p2b_results.npz")

    print("=== VERDICT (corrected) ===")
    print("Pulsar-side red/DM noise gives sigma_phi = sigma_res/A_CW = 1/SNR_pterm ~ O(1) rad")
    print("-- it IS the measurement/SNR floor, already inside the marginalised-noise Fisher")
    print("(D2-D5) and sitting at the knee by construction; NOT an independent coherence term.")
    print("So the coherence axis is not a SEPARATE limit beyond SNR on the pulsar side.")
    print("The genuinely independent coherence question is SOURCE-side stochastic df/f: it")
    print("crosses the knee at df/f ~ 1e-6..1e-5 over kyr lags (SNR 100..5). Whether real SMBHB")
    print("environments reach that is an astrophysics question (gas/stellar coupling over kyr)")
    print("this project cannot answer -> HANDOFF (Taylor/Farr). The deterministic chirp")
    print("(Delta_phi up to ~1e4 rad) is modelled and off-axis.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
