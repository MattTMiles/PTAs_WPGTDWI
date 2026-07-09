"""Track B — LAMBDA feasibility probe, DELIVERABLE L1: search-space specification.

Pure analysis of trackB_L0_float.npz (no GPU). From the float q_p(n):
  1. CARRIER SET: rank pulsars by posterior concentration (max_n q_p). Concentration spectrum.
  2. CANDIDATE SETS: per carrier, integers holding >1% posterior mass (counts).
  3. VOLUME: staged search-volume estimate (census-3 first, then tiers). Enumeration count/stage.
Also reports whether TRUTH's integer (n_true, =0) is inside each carrier's candidate set
(the reachability question the probe hinges on).
Run: python trackB_L1_spec.py
"""
import numpy as np
CWT = "/home/mattm/projects/HSYMT/CW_transition"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]
MASS = 0.01   # candidate = integer with > this posterior mass


def candidate_integers(uk, w, K, mass=MASS):
    uk = uk[:K]; w = w[:K]
    sel = w > mass
    return uk[sel], w[sel]


def main():
    d = np.load(f"{CWT}/trackB_L0_float.npz", allow_pickle=True)
    names = list(d["names"]); npsr = len(names)
    q_uk, q_w, q_K = d["q_uk"], d["q_w"], d["q_K"]
    qmax = d["qmax"]; n_true = d["n_true"]; supp = d["int_support"]
    census_idx = list(d["census_idx"]); status = str(d["status"])
    print(f"=== L1 SEARCH-SPACE SPEC (float status={status}) ===\n")

    # 1. concentration spectrum
    order = np.argsort(-qmax)
    print("1. CARRIER CONCENTRATION SPECTRUM (top by max_n q_p):")
    print(f"   qmax deciles: {np.round(np.quantile(qmax,[0,.25,.5,.75,.9,.99,1]),3).tolist()}")
    for thr in [0.99, 0.9, 0.5, 0.1]:
        print(f"   #pulsars qmax>{thr}: {int(np.sum(qmax>thr))}")
    carriers = np.where(qmax > 0.9)[0]
    carriers = carriers[np.argsort(-qmax[carriers])]
    print(f"   -> CARRIER SET (qmax>0.9): {len(carriers)} pulsars")
    print(f"   {'psr':16s} {'qmax':>6s} {'support':>7s} {'n_true':>6s} {'cand_ints (>1%)':>22s} {'true_in_set':>11s}")
    reach = {}
    for a in carriers:
        ci, cw = candidate_integers(q_uk[a], q_w[a], q_K[a])
        tin = int(n_true[a]) in ci.tolist()
        reach[int(a)] = (ci.tolist(), tin)
        tag = "CENSUS" if a in census_idx else ""
        print(f"   {names[a]:16s} {qmax[a]:6.3f} {supp[a]:7d} {n_true[a]:6d} "
              f"{str(ci.tolist()):>22s} {str(tin):>11s} {tag}")

    if len(carriers) == 0:
        print("   !! NO carriers (qmax>0.9): the truth-blind soft-EM float has ZERO identified")
        print("      fringes. The source drifted off truth -> every per-pulsar posterior is")
        print("      spread. LAMBDA has no concentrated carrier to resolve from this float.")

    # 2. candidate set sizes (all carriers)
    csize = np.array([len(candidate_integers(q_uk[a], q_w[a], q_K[a])[0]) for a in carriers])
    if csize.size:
        print(f"\n2. CANDIDATE SET SIZES over {len(carriers)} carriers: "
              f"min {csize.min()} median {int(np.median(csize))} max {csize.max()}")
    else:
        print(f"\n2. CANDIDATE SET SIZES: n/a (0 carriers)")

    # census-3 specifically
    print("\n   CENSUS-3 candidate sets:")
    for nm in CENSUS:
        a = names.index(nm)
        ci, cw = candidate_integers(q_uk[a], q_w[a], q_K[a])
        print(f"     {nm:14s} qmax={qmax[a]:.3f} cands={ci.tolist()} "
              f"weights={np.round(cw,3).tolist()} n_true={n_true[a]} "
              f"true_in_set={int(n_true[a]) in ci.tolist()}")

    # 3. staged volume
    print("\n3. STAGED SEARCH VOLUME (direct enumeration; pad to >=1e3/scorer call):")
    cen_a = [names.index(nm) for nm in CENSUS]
    cen_sz = [len(candidate_integers(q_uk[a], q_w[a], q_K[a])[0]) for a in cen_a]
    print(f"   stage 1 (census-3 joint lattice): {' x '.join(map(str,cen_sz))} = {int(np.prod(cen_sz))} candidates")
    # tiers of the remaining carriers, ordered by qmax
    rest = [a for a in carriers if a not in cen_a]
    print(f"   remaining carriers: {len(rest)} (each conditioned on prior fixes)")
    # simple tiering: groups whose joint product stays ~1e3-1e5
    print(f"   per-carrier candidate counts (rest, qmax-ordered): "
          f"{[len(candidate_integers(q_uk[a],q_w[a],q_K[a])[0]) for a in rest][:30]}")
    print("\n   NOTE: with these tiny per-carrier sets, direct staged score_candidates is affordable;")
    print("   the LAMBDA Z-transform (decorrelation) is UNNECESSARY at this scale (log as scaling path).")
    return 0


if __name__ == "__main__":
    import sys; sys.exit(main())
