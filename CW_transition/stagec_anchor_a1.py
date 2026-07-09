"""ANCHOR CENSUS — Deliverable A1: dual-prior corrected K table.

K_a = 6 * sigma_EM / dL_a   (candidate fringes in the +/-3 sigma_EM window;
                             dL_a = fringe spacing = c/[f_gw (1-cos mu)]).
dL from stagec_d3_results.npz (10 single-CW geometry draws = the D4 N_CW=1 ensemble).

Two prior columns per pulsar:
  K_canonical : sigma_EM exactly as best_distances.txt delivers (script priority).
  K_optimal   : sigma_EM = MIN sigma among geometric-quality methods:
                {all Cornell techniques, legacy VLBI, Timing_PX, secure PBDOT_Shk}.
                PBDOT_Shk included ONLY if Pbdot SNR>10 AND proper motion measured
                AND masses known (m1,m2 -> GR subtracted).  Flagged when it wins.
Plus an inverse-variance combined sigma (informational 3rd column) where >=2
independent geometric measurements exist; K uses MIN sigma (conservative).

Run:  jug-gpu python stagec_anchor_a1.py
"""
import os, glob, importlib.util
import numpy as np

REPO = "/home/mattm/projects/HSYMT"
CWT = REPO + "/CW_transition"
GET_DIST = "/home/mattm/soft/scripts/MISC/get_distance.py"
PAR_DIRS = [
    "/home/mattm/soft/JUG/data/pulsars/MPTA_data/DR3_partim",
    "/home/mattm/projects/NG/chimecheck/partim",
    "/home/mattm/projects/MPTA/github/mpta-6yr/data/sixth_pass",
]
KPC_TO_PC = 1000.0
B2J_ALIAS = {"B1855+09": "J1857+0943", "B1937+21": "J1939+2134", "B1953+29": "J1955+2908"}

spec = importlib.util.spec_from_file_location("gd", GET_DIST)
gd = importlib.util.module_from_spec(spec); spec.loader.exec_module(gd)

# par index
def scan_psrj(path):
    j = None
    try:
        with open(path) as f:
            for line in f:
                p = line.split()
                if p and p[0] == "PSRJ":
                    j = p[1]; break
    except OSError:
        pass
    return j
jname2par = {}
for d in PAR_DIRS:
    for path in sorted(glob.glob(d + "/*.par")):
        if path.endswith(".test.par"):
            continue
        j = scan_psrj(path)
        if j and j not in jname2par:
            jname2par[j] = path

# A0 priors (canonical) + dL from D3
a0 = np.load(CWT + "/anchor_a0_priors.npz", allow_pickle=True)
fname = a0["fname"]; jname = a0["jname"]
sig_canon_kpc = a0["sig_real_kpc"]; sig_feath_kpc = a0["sig_feather_kpc"]
dist_canon = a0["dist_real_kpc"]; method_canon = a0["method"]
d3 = np.load(CWT + "/stagec_d3_results.npz", allow_pickle=True)
d3names = list(d3["psr_names"]); dL_draws = d3["dL"] * KPC_TO_PC   # (10,116) kpc->pc
# map dL to A0 ordering
idx = [d3names.index(n) for n in fname]
dL = dL_draws[:, idx]                                        # (10,116) pc, A0 order
dL_med = np.median(dL, axis=0); dL_min = dL.min(0); dL_max = dL.max(0)

def cornell_min_sigma(jn):
    """Min distance-sigma (pc) over ALL Cornell techniques (SNR>=2), + (d,ref,tech)."""
    best = None
    for px, err, tech, ref in gd.PARALLAX_CATALOG.get(jn, []):
        if px > 0 and err > 0 and px / err >= 2.0:
            d = 1.0 / px; sd = d * err / px * KPC_TO_PC
            if best is None or sd < best[0]:
                best = (sd, d, tech, ref)
    return best

def secure_pbdot(psr):
    """PBDOT_Shk distance (kpc, sigma_pc) only if secure: Pbdot SNR>10, PM present,
    masses known (m1&m2 -> GR subtracted). Returns (sigma_pc,d_kpc,note) or None."""
    if psr is None or psr.pbdot is None or psr.pbdot_err is None or psr.pb is None:
        return None
    if psr.pmra is None or psr.pmdec is None:
        return None
    if abs(psr.pbdot / psr.pbdot_err) < 10.0:
        return None
    if psr.m1 is None or psr.m2 is None or psr.ecc is None:
        return None          # cannot subtract GR with known masses -> not secure
    est = gd.get_pbdot_distance(psr)
    if est is None:
        return None
    return (est.error_kpc * KPC_TO_PC, est.value_kpc,
            f"Pbdot SNR={abs(psr.pbdot/psr.pbdot_err):.0f}; {est.notes}")

rows = []
for i, fn in enumerate(fname):
    jn = str(jname[i])
    par = jname2par.get(jn)
    psr = gd.parse_par_file(par) if par else None
    if psr is not None:
        psr.name = jn
    cands = []   # (sigma_pc, dist_kpc, method, ref)
    cm = cornell_min_sigma(jn)
    if cm:
        cands.append((cm[0], cm[1], f"Catalog_{cm[2]}", cm[3]))
    if jn in gd.PUBLISHED_VLBI_DISTANCES:
        dd, ee, src = gd.PUBLISHED_VLBI_DISTANCES[jn]
        cands.append((ee * KPC_TO_PC, dd, "VLBI", src))
    if psr is not None:
        tp = gd.get_timing_parallax_distance(psr)
        if tp:
            cands.append((tp.error_kpc * KPC_TO_PC, tp.value_kpc, "Timing_PX", tp.notes))
    pb = secure_pbdot(psr)
    if pb:
        cands.append((pb[0], pb[1], "PBDOT_Shk", pb[2]))

    if cands:
        cands.sort(key=lambda c: c[0])
        sig_opt_pc, dist_opt, meth_opt, ref_opt = cands[0]
        pbdot_wins = meth_opt == "PBDOT_Shk"
        # inverse-variance combo across DISTINCT method-classes (informational)
        seen, combo_terms = set(), []
        for s, dd, m, r in cands:
            cls = m.split("_")[0] if m.startswith(("Catalog", "Timing", "PBDOT")) else m
            if cls not in seen:
                seen.add(cls); combo_terms.append((dd, s))
        if len(combo_terms) >= 2:
            w = np.array([1.0 / s**2 for _, s in combo_terms])
            sig_combo = float(1.0 / np.sqrt(w.sum()))
        else:
            sig_combo = np.nan
    else:
        # no geometric method -> keep feather (parallax-class absent)
        sig_opt_pc = sig_feath_kpc[i] * KPC_TO_PC
        dist_opt = dist_canon[i]; meth_opt = "feather"; ref_opt = "no geometric method"
        pbdot_wins = False; sig_combo = np.nan

    sig_canon_pc = sig_canon_kpc[i] * KPC_TO_PC
    sig_feath_pc = sig_feath_kpc[i] * KPC_TO_PC
    rows.append(dict(
        fn=fn, jn=jn, dL_med=dL_med[i], dL_min=dL_min[i], dL_max=dL_max[i],
        sig_canon=sig_canon_pc, sig_opt=sig_opt_pc, sig_feath=sig_feath_pc,
        sig_combo=sig_combo, meth_canon=str(method_canon[i]), meth_opt=meth_opt,
        ref_opt=ref_opt, pbdot_wins=pbdot_wins,
        K_canon=6 * sig_canon_pc / dL_med[i], K_opt=6 * sig_opt_pc / dL_med[i],
        K_feath=6 * sig_feath_pc / dL_med[i],
        Kopt_lo=6 * sig_opt_pc / dL_max[i], Kopt_hi=6 * sig_opt_pc / dL_min[i],
        Kcan_lo=6 * sig_canon_pc / dL_max[i], Kcan_hi=6 * sig_canon_pc / dL_min[i],
    ))

def counts(key):
    return {t: sum(1 for r in rows if r[key] <= t) for t in (3, 10, 30)}
print("[A1] fiducial dL (median over 10 draws) used for headline K.")
print("[A1] COUNTS (of 116):")
print(f"    K_canonical <=3/10/30 : {counts('K_canon')}")
print(f"    K_optimal   <=3/10/30 : {counts('K_opt')}")
print(f"    K_feather   <=3/10/30 : {counts('K_feath')}")

for key, lab in [("K_opt", "K_OPTIMAL"), ("K_canon", "K_CANONICAL")]:
    sub = sorted([r for r in rows if r[key] <= 10], key=lambda r: r[key])
    print(f"\n[A1] pulsars with {lab} <= 10  (n={len(sub)}):")
    print(f"    {'pulsar':<12s} {'K':>6s} {'sigEM_pc':>9s} {'dL_pc':>7s} {'Krange':>13s}  method / ref")
    for r in sub:
        s = r["sig_opt"] if key == "K_opt" else r["sig_canon"]
        m = (r["meth_opt"] + " " + r["ref_opt"]) if key == "K_opt" else r["meth_canon"]
        lo, hi = (r["Kopt_lo"], r["Kopt_hi"]) if key == "K_opt" else (r["Kcan_lo"], r["Kcan_hi"])
        star = " *PBDOT" if (key == "K_opt" and r["pbdot_wins"]) else ""
        print(f"    {r['fn']:<12s} {r[key]:6.2f} {s:9.3f} {r['dL_med']:7.3f} "
              f"{lo:5.1f}-{hi:5.1f}  {m[:52]}{star}")

# hand checks
for nm in ["J0437-4715", "J1713+0747"]:
    r = next(x for x in rows if x["fn"] == nm)
    print(f"\n[A1] HAND-CHECK {nm}:")
    print(f"    dL_med={r['dL_med']:.4f} pc (range {r['dL_min']:.4f}-{r['dL_max']:.4f})")
    print(f"    sig_canon={r['sig_canon']:.4f} pc [{r['meth_canon']}] -> "
          f"K_canon = 6*{r['sig_canon']:.4f}/{r['dL_med']:.4f} = {r['K_canon']:.3f}")
    print(f"    sig_opt  ={r['sig_opt']:.4f} pc [{r['meth_opt']}: {r['ref_opt'][:50]}] -> "
          f"K_opt   = 6*{r['sig_opt']:.4f}/{r['dL_med']:.4f} = {r['K_opt']:.3f}")
    print(f"    sig_feather={r['sig_feath']:.4f} pc -> K_feather={r['K_feath']:.3f}")
    if not np.isnan(r["sig_combo"]):
        print(f"    inv-var combined sigma={r['sig_combo']:.4f} pc (informational)")

np.savez(CWT + "/anchor_a1_Ktable.npz",
         fn=np.array([r["fn"] for r in rows]), jn=np.array([r["jn"] for r in rows]),
         sig_canon_pc=np.array([r["sig_canon"] for r in rows]),
         sig_opt_pc=np.array([r["sig_opt"] for r in rows]),
         sig_feath_pc=np.array([r["sig_feath"] for r in rows]),
         sig_combo_pc=np.array([r["sig_combo"] for r in rows]),
         dL_med=np.array([r["dL_med"] for r in rows]),
         dL_min=np.array([r["dL_min"] for r in rows]),
         dL_max=np.array([r["dL_max"] for r in rows]),
         K_canon=np.array([r["K_canon"] for r in rows]),
         K_opt=np.array([r["K_opt"] for r in rows]),
         K_feath=np.array([r["K_feath"] for r in rows]),
         meth_opt=np.array([r["meth_opt"] for r in rows]),
         ref_opt=np.array([r["ref_opt"] for r in rows]),
         pbdot_wins=np.array([r["pbdot_wins"] for r in rows]))
print(f"\n[A1] wrote {CWT}/anchor_a1_Ktable.npz")
