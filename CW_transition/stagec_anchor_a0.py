"""ANCHOR CENSUS — Deliverable A0: real EM priors from the canonical get_distance.py.

Live-fetches the Cornell parallax catalogue (via importing get_distance.py, whose
module top-level runs load_parallax_catalog(), cached ~/.cache), matches to the 116
MPTA+NG array pulsars (the data_products/*.feather set), and classifies each pulsar's
best REAL distance by method.

Parallax-class methods (used as a real sigma_EM): Catalog_PX, VLBI, Timing_PX.
Everything else (DM_YMW16/NE2001, PBDOT_Shk, Kinematic, LK_Parallax) => KEEP the
feather prior and flag (do NOT substitute a DM distance as a Gaussian sigma_EM).

Outputs (written into CW_transition/):
  best_distances.txt        canonical get_distance.py output for the 116 array pulsars
  anchor_a0_priors.npz      per-pulsar arrays for A1/A2 (name, dist, sig_real, sig_feather,
                            method, is_parallax_class)
Run:  jug-gpu python stagec_anchor_a0.py
"""
import os, sys, glob, re, importlib.util
import numpy as np

REPO = "/home/mattm/projects/HSYMT"
FEATHER_DIR = REPO + "/data_products/"
OUT_DIR = REPO + "/CW_transition"
GET_DIST = "/home/mattm/soft/scripts/MISC/get_distance.py"
PAR_DIRS = [
    "/home/mattm/soft/JUG/data/pulsars/MPTA_data/DR3_partim",   # MPTA DR3 (J-named)
    "/home/mattm/projects/NG/chimecheck/partim",                # NANOGrav CHIME pars
    "/home/mattm/projects/MPTA/github/mpta-6yr/data/sixth_pass",  # MPTA fallback
]
KPC_TO_PC = 1000.0
PARALLAX_CLASS = {"Catalog_PX", "VLBI", "Timing_PX"}

# ---- import the canonical get_distance.py as a module (triggers live Cornell fetch) ----
spec = importlib.util.spec_from_file_location("get_distance", GET_DIST)
gd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gd)          # <-- fetches + caches the parallax catalogue
print(f"[A0] Cornell parallax catalogue entries: {len(gd.PARALLAX_CATALOG)}")

# ---- build a par-file index: PSRJ -> par path, PSRB -> PSRJ ----
def scan_par_names(path):
    psrj = psrb = None
    try:
        with open(path) as f:
            for line in f:
                p = line.split()
                if not p:
                    continue
                if p[0] == "PSRJ":
                    psrj = p[1]
                elif p[0] in ("PSRB", "PSR") and psrj is None:
                    # some pars only have 'PSR'; treat J* as J, B* as B
                    if p[1].startswith("B"):
                        psrb = p[1]
                    elif p[1].startswith("J"):
                        psrj = p[1]
                elif p[0] == "PSRB":
                    psrb = p[1]
    except OSError:
        pass
    return psrj, psrb

# Standard B->J aliases for the 3 B-name feathers (PINT pars carry only PSRJ, so
# the par scan cannot recover these; Cornell/VLBI are keyed by J-name).
B2J_ALIAS = {
    "B1855+09": "J1857+0943",
    "B1937+21": "J1939+2134",
    "B1953+29": "J1955+2908",
}

jname2par, bname2jname = {}, dict(B2J_ALIAS)
for d in PAR_DIRS:
    for path in sorted(glob.glob(d + "/*.par")):
        if path.endswith(".test.par"):
            continue
        j, b = scan_par_names(path)
        if j and j not in jname2par:      # first dir wins (DR3 preferred)
            jname2par[j] = path
        if b and j:
            bname2jname.setdefault(b, j)
print(f"[A0] indexed {len(jname2par)} PSRJ par files, {len(bname2jname)} B->J maps")

# ---- load the 116 array pulsars + feather priors ----
from enterprise_extensions import load_feathers
feath = load_feathers.load_feathers_from_folder(FEATHER_DIR)
feath = sorted(feath, key=lambda p: p.name)
print(f"[A0] loaded {len(feath)} array pulsars")

rows = []            # per array pulsar
best_for_writer = []  # (PulsarData, DistanceEstimate) for the canonical writer
unmatched_par = []
no_catalog = []

for p in feath:
    fname = p.name
    dist_f, sig_f = float(p.pdist[0]), float(p.pdist[1])   # feather prior (kpc)
    # resolve J-name (for Cornell/VLBI lookup, keyed by PSRJ)
    if fname.startswith("B"):
        jname = bname2jname.get(fname)
        if jname is None:
            # fall back: no B->J map found
            jname = fname
    else:
        jname = fname
    par = jname2par.get(jname)
    if par is not None:
        psr = gd.parse_par_file(par)
        if psr is None:
            psr = gd.PulsarData(name=jname, raj="00:00:00", decj="00:00:00")
        # force the catalogue key to the resolved J-name
        psr.name = jname
    else:
        unmatched_par.append((fname, jname))
        psr = gd.PulsarData(name=jname, raj="00:00:00", decj="00:00:00")

    ests = gd.get_all_distances(psr)
    best = gd.select_best_distance(ests) if ests else None
    if best is not None:
        best_for_writer.append((psr, best))

    has_cat = any(e.method == "Catalog_PX" for e in ests)
    if not has_cat and jname not in gd.PARALLAX_CATALOG:
        no_catalog.append(fname)

    if best is not None and best.method in PARALLAX_CLASS:
        sig_real = best.error_kpc          # kpc, real parallax-class
        dist_real = best.value_kpc
        method = best.method
        note = best.notes
        is_px = True
    else:
        sig_real = sig_f                   # KEEP feather (flagged)
        dist_real = dist_f
        method = "feather(" + (best.method if best else "none") + ")"
        note = "kept feather prior; no parallax-class measurement"
        is_px = False

    rows.append(dict(fname=fname, jname=jname, dist_real=dist_real, sig_real=sig_real,
                     dist_f=dist_f, sig_f=sig_f, method=method, is_px=is_px,
                     note=note, has_par=par is not None))

# ---- write canonical best_distances.txt (frozen copy pinned in CW_transition/) ----
gd.write_best_distances(best_for_writer, os.path.join(OUT_DIR, "best_distances.txt"))

# ---- coverage by method ----
from collections import Counter
cov = Counter()
for r in rows:
    m = r["method"].split("(")[0]
    cov[m] += 1
print("\n[A0] COVERAGE by method (of 116):")
for m in ["Catalog_PX", "VLBI", "Timing_PX", "feather"]:
    print(f"    {m:<12s}: {cov.get(m,0):3d}")
print(f"    (parallax-class total: {sum(1 for r in rows if r['is_px'])})")

print(f"\n[A0] par-file unmatched (used name-only lookup): {len(unmatched_par)}")
for fn, jn in unmatched_par:
    print(f"      {fn:<12s} (J={jn})")
print(f"[A0] no Cornell catalogue entry at all: {len(no_catalog)}")
print("      " + " ".join(no_catalog))

# ---- 15 tightest REAL priors (parallax-class only) ----
px_rows = [r for r in rows if r["is_px"]]
px_rows.sort(key=lambda r: r["sig_real"])
print("\n[A0] 15 tightest-prior pulsars (parallax-class):")
print(f"    {'pulsar':<12s} {'d_kpc':>8s} {'sigR_pc':>9s} {'sigF_pc':>9s} {'R=sigR/sigF':>11s}  method")
for r in px_rows[:15]:
    sr, sf = r["sig_real"]*KPC_TO_PC, r["sig_f"]*KPC_TO_PC
    print(f"    {r['fname']:<12s} {r['dist_real']:8.4f} {sr:9.3f} {sf:9.3f} {sr/sf:11.3f}  {r['method']}")

# ---- GATE: J0437-4715 ----
g = next((r for r in rows if r["fname"] == "J0437-4715"), None)
print("\n[A0] GATE J0437-4715:")
if g:
    print(f"    real distance = {g['dist_real']:.4f} kpc  (expect ~0.157)")
    print(f"    real sigma_EM = {g['sig_real']*KPC_TO_PC:.3f} pc   method={g['method']}")
    print(f"    reference/notes: {g['note']}")

# ---- save arrays for A1/A2 ----
np.savez(os.path.join(OUT_DIR, "anchor_a0_priors.npz"),
         fname=np.array([r["fname"] for r in rows]),
         jname=np.array([r["jname"] for r in rows]),
         dist_real_kpc=np.array([r["dist_real"] for r in rows]),
         sig_real_kpc=np.array([r["sig_real"] for r in rows]),
         dist_feather_kpc=np.array([r["dist_f"] for r in rows]),
         sig_feather_kpc=np.array([r["sig_f"] for r in rows]),
         method=np.array([r["method"] for r in rows]),
         is_parallax_class=np.array([r["is_px"] for r in rows]),
         note=np.array([r["note"] for r in rows]))
print(f"\n[A0] wrote {OUT_DIR}/best_distances.txt and anchor_a0_priors.npz")
