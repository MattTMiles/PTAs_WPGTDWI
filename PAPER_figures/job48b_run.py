"""JOB 48b -- Job 48 recomputed with fisher_nested_fast (closed-form Dirichlet Fisher,
nested draws, one Cholesky per sky-arm), 200 skies, plotted in Figure 1's style.
Estimator: booked prong2 'marg' (source phase/amplitude marginalised, other pulsars'
distances held fixed) throughout."""
import sys, time, json, csv
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import jax, jax.numpy as jnp
import prong2_transition as P
import fisher_nested_fast as FN

OUT = "/home/mattm/projects/HSYMT/PAPER_figures/"
T = 15*P.YR; SIG = 1e-7; NMAX = 1000; NSKY = int(sys.argv[1]) if len(sys.argv) > 1 else 200
FB = 1.0/T
FIG1 = [1, 2, 3, 4, 6, 10, 16, 25, 40, 63, 100, 158, 251, 398, 631, 1000]
GRID = np.unique(np.concatenate([np.arange(1, 31), np.round(np.geomspace(30, 1000, 92)).astype(int), FIG1]))
ARMS = [("B8", "8 frequency bins (as Figure 1)", ("log", 3e-9, 2e-8), T*(2e-8-3e-9)),
        ("B4", "4 frequency bins", ("lin", 2*FB, 6*FB), 4.0),
        ("B2", "2 frequency bins", ("lin", 3*FB, 5*FB), 2.0),
        ("B1", "1 frequency bin", ("lin", 3*FB, 4*FB), 1.0),
        ("SF", "all at one frequency", ("single", 3.5*FB), 1.0)]

def nested_draw(s):
    """one set of NMAX sources per sky; every arm uses the same underlying variates."""
    rng = np.random.default_rng([4802, s])
    th = np.arccos(rng.uniform(-1, 1, NMAX)); ph = rng.uniform(0, 2*np.pi, NMAX)
    ps = rng.uniform(0, np.pi, NMAX); ci = rng.uniform(-1, 1, NMAX)
    u = rng.uniform(0, 1, NMAX); ph0 = rng.uniform(0, 2*np.pi, NMAX)
    return th, ph, ps, ci, u, ph0

def arm_srcs(base, fspec):
    th, ph, ps, ci, u, ph0 = base
    if fspec[0] == "log":
        l10f = np.log10(fspec[1]) + (np.log10(fspec[2]) - np.log10(fspec[1]))*u
    elif fspec[0] == "lin":
        l10f = np.log10(fspec[1] + u*(fspec[2] - fspec[1]))
    else:
        l10f = np.full(NMAX, np.log10(fspec[1]))
    return np.stack([th, ph, ps, ci, l10f, np.full(NMAX, -14.3), ph0], 1)

t_all = time.time()
nf = FN.NestedFisher(P.build_array(116), T, 7.0, SIG, NMAX, GRID)
print(f"grid: {len(nf.grid)} N values (every integer 1-30, ~60/decade to 1000, incl. Figure 1 points)")
# warm-up / compile (not counted in the per-arm timing)
t0 = time.time(); nf.run(arm_srcs(nested_draw(0), ARMS[0][2])); t_compile = time.time() - t0
nf.t_asm = nf.t_fac = 0.0
print(f"compile + first call: {t_compile:.1f} s")
RES = {}
for key, lab, fspec, B in ARMS:
    CM = np.zeros((NSKY, NMAX)); MM = np.full((NSKY, len(nf.grid)), np.nan); nok = []
    t0 = time.time()
    for s in range(NSKY):
        r = nf.run(arm_srcs(nested_draw(s), fspec))
        CM[s] = np.median(r["cond"], 0)
        MM[s] = np.median(r["marg"], 1)         # NaN rows stay NaN if Cholesky broke before them
        nok.append(r["n_chol_ok"])
    nok = np.array(nok)
    RES[key] = dict(CM=CM, MM=MM, nok=nok, B=B, lab=lab)
    print(f"  arm {key}: {NSKY} skies in {time.time()-t0:.1f} s; Cholesky factors all 1000 sources in"
          f" {np.sum(nok == NMAX)}/{NSKY} skies; min/median failure N = {nok.min()}/{int(np.median(nok))}", flush=True)
t_compute = time.time() - t_all
print(f"COMPUTE wall {t_compute:.1f} s = compile {t_compile:.1f} + assembly {nf.t_asm:.1f}"
      f" + factorisation/solve {nf.t_fac:.1f} + python/draw/transfer {t_compute - t_compile - nf.t_asm - nf.t_fac:.1f}")
np.savez(OUT + "job48_data/job48b_results.npz", grid=nf.grid,
         **{f"{k}_{f}": v[f] for k, v in RES.items() for f in ("CM", "MM", "nok")})
json.dump(dict(t_compute=t_compute, t_compile=t_compile, t_asm=nf.t_asm, t_fac=nf.t_fac),
          open(OUT + "job48_data/job48b_timing.json", "w"))
