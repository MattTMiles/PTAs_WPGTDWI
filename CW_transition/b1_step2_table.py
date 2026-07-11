"""B1 STEP 2 -- assemble the three-tier referendum table from the saved npz files.

Recomputes, from banked quantities only (no GPU):
  f_tier +- sigma_f          (delta method through the logistic; sigma from the SMC s.e.m.)
  the pre-registered +-2sigma PASS / FAIL / STRADDLE verdict
  break-even mc box at f = 0.95 and the DEFICIT FACTOR vs the tier's own mc box

The deficit factor is the number the paper quotes: how much external chirp-mass information --
eccentricity's harmonic lever (kappa), a better host D_L, or a louder source's sigma(log10_h) --
must be supplied ON TOP of that tier's conditioning for the targeted posterior to concentrate.

Usage:  python b1_step2_table.py
"""
import sys
import numpy as np

CWT = "/home/mattm/projects/HSYMT/CW_transition"
N_LOUD = 3
SCALE_MC = 0.5                      # phi_scale for log10_mc
LOGIT95 = np.log(0.95 / 0.05)


def verdict_of(f, sig_f, gate_ok):
    if not gate_ok or not np.isfinite(sig_f):
        return "NO VERDICT (SMC gate failed)"
    if f - 2 * sig_f >= 0.95:
        return "PASS"
    if f + 2 * sig_f < 0.95:
        return "FAIL"
    return "STRADDLE"


def main():
    rows = []
    for tier in ("A", "B", "C"):
        try:
            z = np.load(f"{CWT}/b1_referendum_tier{tier}.npz", allow_pickle=True)
        except FileNotFoundError:
            print(f"  tier {tier}: not run yet")
            continue
        lnZn = float(z["lnZn_quad"]); lnZb = float(z["lnZbox"])
        lnZb_e = float(z["lnZbox_err"])
        d = lnZn - lnZb
        f = 1.0 / (1.0 + np.exp(-d))
        sig_f = f * (1 - f) * lnZb_e
        gate_ok = bool(z["gate_ok"]) if "gate_ok" in z.files else False
        half_mc = np.asarray(z["half_mc"]) if "half_mc" in z.files else np.full(N_LOUD, np.nan)
        half_f = np.asarray(z["half_f"]) if "half_f" in z.files else np.full(N_LOUD, np.nan)
        lam_mc = float(np.exp((d - LOGIT95) / N_LOUD))
        be_mc_dex = half_mc * lam_mc * SCALE_MC
        rows.append(dict(tier=tier, lnZn=lnZn, lnZb=lnZb, lnZb_e=lnZb_e, d=d, f=f, sig_f=sig_f,
                         gate_ok=gate_ok, verdict=verdict_of(f, sig_f, gate_ok),
                         half_f=half_f, half_mc=half_mc, lam_mc=lam_mc, be_mc_dex=be_mc_dex,
                         deficit=1.0 / lam_mc,
                         doubling=float(z["doubling"]) if "doubling" in z.files else np.nan,
                         npos=int(z["npos"]) if "npos" in z.files else -1,
                         nseed=int(z["nseed"]) if "nseed" in z.files else -1,
                         spread=float(z["seed_spread"]) if "seed_spread" in z.files else np.nan,
                         acc_hi=float(z["acc_hi"]) if "acc_hi" in z.files else np.nan))

    if not rows:
        print("no tiers available"); return 1

    print("\n=== B1 STEP 2: the targeted (f, mc) referendum ===\n")
    print(f"{'tier':>5s} {'lnZ_needle':>12s} {'lnZ_box':>12s} {'d (nat)':>9s} "
          f"{'f':>9s} {'+-2sig':>10s} {'verdict':>10s}")
    for r in rows:
        band = f"+-{2*r['sig_f']:.3g}"
        print(f"{r['tier']:>5s} {r['lnZn']:12.3f} {r['lnZb']:12.3f} {r['d']:+9.3f} "
              f"{r['f']:9.4g} {band:>10s} {r['verdict']:>10s}")

    print(f"\n{'tier':>5s} {'mc box (dex)':>14s} {'break-even mc (dex)':>21s} {'DEFICIT':>10s}")
    for r in rows:
        mc_dex = np.round(r['half_mc'] * SCALE_MC, 5)
        be = np.round(r['be_mc_dex'], 6)
        print(f"{r['tier']:>5s} {str(mc_dex):>14s} {str(be):>21s} {r['deficit']:9.3g}x")

    print(f"\n{'tier':>5s} {'seeds':>6s} {'spread(nat)':>12s} {'acc_hi':>8s} "
          f"{'doubling':>9s} {'pos-curv':>9s} {'gate':>6s}")
    for r in rows:
        print(f"{r['tier']:>5s} {r['nseed']:6d} {r['spread']:12.3f} {r['acc_hi']:8.3f} "
              f"{r['doubling']:9.3f} {r['npos']:6d}/6 {str(r['gate_ok']):>6s}")

    print("\nGATES: SMC seed spread <= 0.3 nat AND high-beta acceptance >= 0.25;")
    print("       quadrature doubling |dlnJ| < 0.2. A tier with either failing gets NO VERDICT.")
    print("PRE-REGISTERED: PASS iff f - 2sig >= 0.95; FAIL iff f + 2sig < 0.95; else STRADDLE")
    print("       (a STRADDLE names no branch -- it is reported as a measured near-miss with the")
    print("        deficit factor as its price tag).")
    np.savez(f"{CWT}/b1_step2_table.npz", rows=np.array([rows], dtype=object))
    print(f"\nsaved b1_step2_table.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
