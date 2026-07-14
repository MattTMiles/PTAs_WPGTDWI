"""KINDLE pilot analysis — turn the k_pilot_*.npz banks into the §4 table.
The load-bearing question: is the loop inert ONLY because it starts at truth?
Reads KINDLE_results/. CPU-only."""
import glob, os
import numpy as np

KRES = "/data/taylor_group/matt_miles/PTAs_WPGTDWI/KINDLE_results"


def main():
    files = sorted(glob.glob(f"{KRES}/k_pilot*.npz"))
    if not files:
        print("no pilot banks yet")
        return 1
    rows = []
    for f in files:
        d = np.load(f, allow_pickle=True)
        stepn = np.asarray(d["traj_step_norm"], float)
        moved = bool(np.any(stepn > 0))
        smo = np.asarray(d["traj_src_mc_off"], float)
        rows.append(dict(
            sm=str(d["start_mode"]), sc=float(d["start_scale"]), scr=bool(d["scrambled"]),
            gid=int(d["geo_id"]), ns=int(d["noise_seed"]),
            off_mc=float(d["off_mc_sig"]), C0s=int(d["C0_start"]), C0t=int(d["C0_truth"]),
            Cf=int(d["C_fix"]), dN=int(d["dN_final"]), moved=moved,
            nstep=int(np.sum(stepn > 0)), niter=int(d["n_iter"]),
            mc_off0=float(smo[0]), mc_offf=float(smo[-1]),
            wrong=int(np.asarray(d["traj_wrong"], float)[-1]),
            in_basin=bool(d["start_in_basin"]),
            mflow=float(np.asarray(d["traj_margin_flow"], float).sum())))

    def grp(scr, sm):
        return [r for r in rows if r["scr"] == scr and r["sm"] == sm]

    print("=" * 96)
    print("SIGNAL — arm (a) truth vs arm (b-FIX) 1sig(mc); does a known displacement get "
          "repaired?")
    print("=" * 96)
    for sm, lab in [("truth", "arm(a) TRUTH"), ("fix_mc", "arm(b-FIX) 1sig(mc)")]:
        g = grp(False, sm)
        if not g:
            continue
        moved = sum(r["moved"] for r in g)
        dN = np.array([r["dN"] for r in g])
        print(f"\n{lab}: {len(g)} reals")
        print(f"  M-step accepted >=1 step : {moved}/{len(g)}   "
              f"(g1: truth-start loops are 25/30 FROZEN)")
        print(f"  additive dN              : {dN.tolist()}   "
              f"(ignitions {int((dN>0).sum())}, losses {int((dN<0).sum())})")
        print(f"  |C_0(start)|             : {[r['C0s'] for r in g]}")
        if sm == "fix_mc":
            print(f"  mc offset (sigma units)  : {[round(r['off_mc'], 2) for r in g]}")
            trail = [f"{r['mc_off0']:.4f}->{r['mc_offf']:.4f}" for r in g]
            print(f"  src_mc_off  start->final : {trail}")
    print("\n" + "=" * 96)
    print("SCRAMBLED — must stay SILENT (any detection = STOP + anatomy)")
    print("=" * 96)
    for sm in ("truth", "fix_mc"):
        g = grp(True, sm)
        if not g:
            continue
        cf = [r["Cf"] for r in g]
        wr = [r["wrong"] for r in g]
        print(f"  {sm:8s}: {len(g)} reals  final |C_fix| {cf}  wrong {wr}  "
              f"-> {'SILENT' if sum(cf)==0 else 'DETECTION — STOP'}")

    # basin sanity
    print("\n" + "-" * 96)
    sig = [r for r in rows if not r["scr"]]
    print(f"start-in-basin (|C_0(start)|>=1): "
          f"{sum(r['in_basin'] for r in sig)}/{len(sig)} signal reals")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
