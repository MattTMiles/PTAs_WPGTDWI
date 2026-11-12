#!/usr/bin/env python3
"""GLACIER-LITE harvest -- pure numpy/filesystem, no GPU, no jax. Reads the gl_lite_*
banks and emits the report's tables as VALUES (the brief: "values, not just plots").

  * floor re-cut summary (per lap, per cell: floor, estimator, zero-fraction)
  * the per-lap table for every realisation
  * certifications-vs-lap and residual-vs-lap curves
  * the anchor cell vs the two published papers
  * the pre-logged prediction quoted next to outcomes
"""
import os, sys, glob
import numpy as np

OUT = os.environ.get("GLACIER_LITE_OUT",
                     "/data/taylor_group/matt_miles/PTAs_WPGTDWI/GLACIER_LITE_results")
KPC_PC = 1000.0

# Published anchors (fetched 2026-07-27; setups recorded so the comparison is honest)
LIT = {
    "2503.23017": dict(
        title="Yu & Pan, 'Sub-parsec precision measurement of pulsar distances with "
              "nanohertz gravitational waves'",
        npsr=20, tspan_yr=30, sigma_n_ns=20, n_src="multiple resolved SMBHBs",
        priors="SKA-era timing-parallax",
        claims=[("D ~ 1.0 kpc", 0.4, "better than 0.4 pc"),
                ("D ~ 2.2 kpc", 1.0, "better than 1 pc")]),
    "2512.10729": dict(
        title="Xiao et al., 'Two-Dimensional Pulsar Distance Inference from Nanohertz "
              "Gravitational Waves'",
        npsr=20, tspan_yr=30, sigma_n_ns=20, n_src="N = 2,3,4,5 (+7,10)",
        priors="SKA-era timing-parallax (Lee et al. relation)",
        # per-pulsar table, 4 sources at 40 ns, 30 yr, isotropic
        table=[("J0030+0451", 0.323, 0.04), ("J0613-0200", 0.99, 1.02),
               ("J0751+1807", 1.17, 0.10), ("J1012+5307", 1.07, 0.09),
               ("J1022+1001", 0.85, 0.06), ("J1024-0719", 0.98, 0.09),
               ("J1455-3330", 0.76, 0.05), ("J1600-3053", 1.39, 0.13),
               ("J1640+2224", 1.08, 0.09), ("J1713+0747", 1.136, 0.09),
               ("J1730-2304", 0.48, 0.03), ("J1744-1134", 0.388, 0.03),
               ("J1751-2857", 0.79, 0.07), ("J1801-1417", 1.0, 0.09),
               ("J1804-2717", 0.8, 0.07), ("J1857+0943", 1.11, 0.24),
               ("J1909-3744", 1.06, 0.23), ("J1911+1347", 2.2, 4.32),
               ("J1918-0642", 1.3, 0.41), ("J2124-3358", 0.47, 0.05)],
        claims=[("most D <~ 1.4 kpc, ~15 yr", 1.0, "below 1 pc")]),
}


def _load(pat):
    return sorted(glob.glob(os.path.join(OUT, pat)))


def laps_for(wkey, real, renorm=False):
    r = "_ng15" if renorm else ""
    rows = [np.load(f, allow_pickle=True)
            for f in _load(f"gl_lite_loop_{wkey}_r{real}_T15_lit{r}_lap*__*.npz")]
    return sorted(rows, key=lambda z: int(z["lap"]))


def laps_complete(pw_tag, real):
    rows = [np.load(f, allow_pickle=True)
            for f in _load(f"gl_lite_cmpl_w20_pw{pw_tag}_r{real}_T15_lit_lap*__*.npz")]
    return sorted(rows, key=lambda z: int(z["lap"]))


def completeness_block():
    """PART 3 -- C vs lap, saturation, K distribution, and the never-certified autopsy."""
    print("\n" + "=" * 100)
    print("PART 3 -- THE COMPLETENESS AXIS (20 ns; 3 loudest at T2; prior widths scaled)")
    print("=" * 100)
    any_ = False
    print(f"\n{'pw':>6} {'real':>4} {'satLap':>6} {'C':>7} {'cumTrue':>7} "
          f"{'K_min':>7} {'K_med':>7} {'K_max':>8} {'floored':>7} "
          f"{'nAudib':>6} {'nAmbig':>6}")
    for tag in ("1", "3", "10", "30"):
        for real in (0, 1):
            zs = laps_complete(tag, real)
            if not zs:
                continue
            any_ = True
            z = zs[-1]
            print(f"{'1/'+tag:>6} {real:>4} {int(z['lap']):>6} "
                  f"{float(z['completeness']):>7.4f} {int(z['cum_true_cert']):>7} "
                  f"{float(z['K_min']):>7.0f} {float(z['K_med']):>7.0f} "
                  f"{float(z['K_max']):>8.0f} {int(z['pw_n_floored']):>7} "
                  f"{int(z['n_audibility_limited']):>6} "
                  f"{int(z['n_ambiguity_limited']):>6}")
    if not any_:
        print("   [no completeness banks present]")
        return
    print("\nC vs lap (cumulative TRUE certifications / 116, DISTINCT pulsars):")
    for tag in ("1", "3", "10", "30"):
        for real in (0, 1):
            zs = laps_complete(tag, real)
            if not zs:
                continue
            print(f"   pw 1/{tag:<3} r{real}: " +
                  " ".join(f"{float(z['completeness']):.4f}" for z in zs))
    print("\nAUTOPSY of the never-certified set (final lap of each cell):")
    for tag in ("1", "3", "10", "30"):
        for real in (0, 1):
            zs = laps_complete(tag, real)
            if not zs:
                continue
            z = zs[-1]
            ev = np.asarray(z["ever_true_mask"])
            aud = np.asarray(z["autopsy_audibility_limited"])
            amb = np.asarray(z["autopsy_ambiguity_limited"])
            snr = np.asarray(z["autopsy_snr_p"])
            K = np.asarray(z["K_counted"])
            n_un = int((~ev).sum())
            print(f"   pw 1/{tag:<3} r{real}: uncertified {n_un}/116 -> "
                  f"audibility-limited {int(aud.sum())}, "
                  f"ambiguity-limited {int(amb.sum())}")
            if aud.any():
                print(f"        audibility set: SNR_p median {np.median(snr[aud]):.2f}, "
                      f"K median {np.median(K[aud]):.0f}")
            if amb.any():
                print(f"        ambiguity  set: SNR_p median {np.median(snr[amb]):.2f}, "
                      f"K median {np.median(K[amb]):.0f}")
    print(f"\n   autopsy rule (pre-stated): {str(zs[-1]['autopsy_rule'])}")


def secondary_block():
    print("\n" + "=" * 100)
    print("SECONDARY ARM -- NG15-renormalised contrast (A_equivalent -14.60), 20 ns")
    print("=" * 100)
    any_ = False
    for real in (0, 1):
        zs = laps_for("w20", real, renorm=True)
        if not zs:
            continue
        any_ = True
        print(f"\n--- w20 r{real} [renorm; A_eq {float(zs[0]['a_equivalent']):.3f}; "
              f"shift {float(zs[0]['renorm_dex']):+.3f} dex] ---")
        print(f"{'lap':>3} {'newTrue':>7} {'cumTrue':>7} {'nCert':>5} {'wrong':>5} "
              f"{'nRes':>4} {'A_bg':>8} {'floor':>7} {'zf':>5} {'wall_s':>7}")
        for z in zs:
            print(f"{int(z['lap']):>3} {int(z['n_true_cert']):>7} "
                  f"{int(z['cum_true_cert']):>7} {int(z['n_cert']):>5} "
                  f"{int(z['wrong_cert']):>5} {int(z['n_resolved']):>4} "
                  f"{float(z['a_bg']):>8.3f} {float(z['floor']):>7.3f} "
                  f"{float(z['zero_fraction']):>5.2f} {float(z['wall_s']):>7.0f}")
    if not any_:
        print("   [no secondary-arm banks present]")


def per_lap_table():
    print("=" * 100)
    print("PART 1 -- THE PER-LAP TABLE (one block per realisation)")
    print("=" * 100)
    curves = {}
    for wkey in ("w50", "w20"):
        for real in range(4):
            zs = laps_for(wkey, real)
            if not zs:
                continue
            ns = float(zs[0]["white_ns"])
            print(f"\n--- {wkey} ({ns:.0f} ns white), realisation r{real} "
                  f"[{len(zs)} laps banked] ---")
            print(f"{'lap':>3} {'newTrue':>7} {'cumTrue':>7} {'nCert':>5} {'wrong':>5} "
                  f"{'newRes':>6} {'nRes':>4} {'A_bg':>8} {'+-':>6} {'floor':>7} "
                  f"{'zf':>5} {'estimator':>20} {'kill':>10} {'wall_s':>7}")
            c_true, c_abg, c_lap = [], [], []
            for z in zs:
                kill = ("n/a" if not bool(z["kill_scored"])
                        else ("would-KILL" if bool(z["kill_would_kill"]) else "RETAIN"))
                print(f"{int(z['lap']):>3} {int(z['n_true_cert']):>7} "
                      f"{int(z['cum_true_cert']):>7} {int(z['n_cert']):>5} "
                      f"{int(z['wrong_cert']):>5} {int(z['n_new_resolved']):>6} "
                      f"{int(z['n_resolved']):>4} {float(z['a_bg']):>8.3f} "
                      f"{float(z['a_bg_sig']):>6.3f} {float(z['floor']):>7.3f} "
                      f"{float(z['zero_fraction']):>5.2f} {str(z['floor_est']):>20} "
                      f"{kill:>10} {float(z['wall_s']):>7.0f}")
                c_lap.append(int(z["lap"])); c_true.append(int(z["cum_true_cert"]))
                c_abg.append(float(z["a_bg"]))
            curves[(wkey, real)] = (c_lap, c_true, c_abg)
    return curves


def curves_block(curves):
    print("\n" + "=" * 100)
    print("THE CORE PRODUCT -- certifications vs lap, and residual power vs lap (VALUES)")
    print("=" * 100)
    for key in ("w50", "w20"):
        got = {r: v for (w, r), v in curves.items() if w == key}
        if not got:
            continue
        print(f"\n[{key}] cumulative TRUE certifications by lap")
        for r, (lp, ct, ab) in sorted(got.items()):
            print(f"   r{r}: " + " ".join(f"{v:>3d}" for v in ct))
        print(f"[{key}] fitted background amplitude log10 A_bg by lap (the drain)")
        for r, (lp, ct, ab) in sorted(got.items()):
            print(f"   r{r}: " + " ".join(f"{v:>8.3f}" for v in ab))
            d = np.diff(ab)
            mono = "monotone-decreasing" if np.all(d < 0) else (
                "monotone-increasing" if np.all(d > 0) else "non-monotone")
            print(f"        deltas: " + " ".join(f"{v:+8.4f}" for v in d) + f"   [{mono}]")


def anchor_block():
    fs = _load("gl_lite_anchor_*__*.npz")
    print("\n" + "=" * 100)
    print("PART 2 -- THE ANCHOR CELL vs THE PUBLISHED SINGLE-SHOT CLAIMS")
    print("=" * 100)
    for k, v in LIT.items():
        print(f"\narXiv:{k} -- {v['title']}")
        print(f"   setup: {v['npsr']} pulsars, {v['tspan_yr']} yr, "
              f"sigma_n = {v['sigma_n_ns']} ns, sources: {v['n_src']}")
        print(f"   distance priors: {v['priors']}")
        for c in v.get("claims", []):
            print(f"   claim: {c[0]:>28}  ->  {c[2]}")
        if "table" in v:
            t = np.array([(d, s) for _, d, s in v["table"]])
            w = t[:, 0] < 1.4
            print(f"   per-pulsar table (4 src, 40 ns, 30 yr): median sigma_L "
                  f"{np.median(t[:,1]):.3f} pc; D<1.4 kpc median "
                  f"{np.median(t[w,1]):.3f} pc; sub-pc fraction "
                  f"{float((t[:,1]<1.0).mean()):.2f}")
    if not fs:
        print("\n   [anchor bank not present]")
        return
    for f in fs:
        z = np.load(f, allow_pickle=True)
        sc = np.asarray(z["sigma_L_cond_kpc"]) * KPC_PC
        sm = np.asarray(z["sigma_L_marg_kpc"]) * KPC_PC
        L0 = np.asarray(z["L0_kpc"])
        w1 = np.asarray(z["within_1kpc"])
        w14 = L0 < 1.4
        on_true = np.asarray(z["on_true"])
        q = np.asarray(z["q_max"])
        print(f"\nTHIS CAMPAIGN -- {os.path.basename(f)}")
        print(f"   setup: 116 pulsars, 22.15 yr, sigma_n = {float(z['white_ns']):.0f} ns, "
              f"{int(z['n_src_fed'])} sources at FULL truth (T3)")
        print(f"   distance priors: REAL published (tier 'lit') -- NOT idealized")
        print(f"   sigma_L_cond (within-fringe, CW known): median {np.median(sc):.4f} pc; "
              f"D<1 kpc median {np.median(sc[w1]):.4f} pc; "
              f"D<1.4 kpc median {np.median(sc[w14]):.4f} pc")
        print(f"        sub-pc fraction: all {float((sc<1).mean()):.3f}, "
              f"D<1 kpc {float((sc[w1]<1).mean()):.3f}, "
              f"D<1.4 kpc {float((sc[w14]<1).mean()):.3f}")
        print(f"   sigma_L_marg (fringe-comb, ambiguity kept): median {np.median(sm):.4f} pc; "
              f"D<1 kpc median {np.median(sm[w1]):.4f} pc")
        print(f"        sub-pc fraction: all {float((sm<1).mean()):.3f}, "
              f"D<1 kpc {float((sm[w1]<1).mean()):.3f}")
        print(f"   fringe: MAP == true {int(on_true.sum())}/{len(on_true)}; "
              f"q_max median {np.median(q):.3f}; K_counted median "
              f"{np.median(np.asarray(z['K_counted'])):.0f}; "
              f"dL median {np.median(np.asarray(z['dL_kpc']))*KPC_PC:.2f} pc")
        print(f"   step-size stability: median {np.median(np.asarray(z['step_agree'])):.2e}")


def floor_block():
    print("\n" + "=" * 100)
    print("FLOOR RE-CUT SUMMARY (mandatory; banked floors do not apply at these levels)")
    print("=" * 100)
    print(f"{'cell':>22} {'lap':>3} {'floor':>8} {'err':>7} {'zf':>5} {'estimator':>22} "
          f"{'n_null':>6} {'off_max':>8}")
    for wkey in ("w50", "w20"):
        for real in range(4):
            for z in laps_for(wkey, real):
                off = np.asarray(z["null_offenders"])
                print(f"{wkey+'_r'+str(real):>22} {int(z['lap']):>3} "
                      f"{float(z['floor']):>8.3f} {float(z['floor_err']):>7.3f} "
                      f"{float(z['zero_fraction']):>5.2f} {str(z['floor_est']):>22} "
                      f"{int(z['n_null']):>6} {off.max():>8.3f}")


def main():
    fs = _load("gl_lite_*__*.npz")
    print(f"GLACIER-LITE harvest: {len(fs)} banks in {OUT}\n")
    if fs:
        z = np.load(fs[0], allow_pickle=True)
        if "prediction" in z:
            print("PRE-LOGGED PREDICTION (verbatim, banked before any scoring):")
            print("  " + str(z["prediction"]))
            print(f"  cksum_prediction = {str(z.get('cksum_prediction',''))}  "
                  f"cksum_driver = {str(z.get('cksum_driver',''))}\n")
    zz = [np.load(f, allow_pickle=True) for f in _load("gl_lite_loop_*_lap0__*.npz")]
    if zz:
        z = zz[0]
        print("LOUDNESS LABEL (the phrase 'NG15-consistent' is struck from this campaign):")
        print(f"  primary arm A_equivalent = {float(z['a_equivalent']):.3f}  "
              f"(band-power ratio {float(z['band_power_ratio']):.1f}x NG15)\n")
    floor_block()
    curves = per_lap_table()
    curves_block(curves)
    secondary_block()
    completeness_block()
    anchor_block()
    # GPU-hr ledger
    tot = 0.0
    for f in _load("gl_lite_loop_*__*.npz") + _load("gl_lite_anchor_*__*.npz"):
        z = np.load(f, allow_pickle=True)
        if "wall_s" in z:
            tot += float(z["wall_s"])
    print("\n" + "=" * 100)
    print(f"LEDGER: banked in-cell wall = {tot/3600:.2f} GPU-hr "
          f"(excludes venue build + JIT; see the job logs for the billed total)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
