"""JOB 48 -- the turnover per frequency bin: confine the sources to fewer bins.

Machinery: prong2_transition (Prong 2, Stage A.1) VERBATIM -- analytic per-pulsar
gradients, joint-array Fisher, eigh pseudo-inverse at rcond 1e-10. The only new
ingredient is WHERE IN FREQUENCY the sources are placed. Assembly is chunked over
pulsars (F_ss is a sum over pulsars, so the chunk sum is the same matrix) so that
N up to 4200 fits on the 24 GB card; equivalence to the booked fisher_blocks is
asserted below before any new arm runs.

TWO ESTIMATORS, both reported (the job's S0 and the booked code disagree):
  booked : I_i = [F_LL - F_Ls F_ss^+ F_sL]_ii  -- sources marginalised, OTHER
           pulsars' distances HELD FIXED. This is what Figure 1 plots and what
           reproduces 0.99/0.87/0.53/0.08.
  full   : I_i = 1/[S^-1]_ii with S the matrix above -- every other parameter
           marginalised, including the other 115 distances (the job's S0 wording).

PRE-REGISTERED (written before any Job 48 arm ran):
  booked estimator -> job's H1: ratio ~ 1 - n/N_psr, knee n ~ N_psr/2 = 58 for
      every B; same-frequency: linear to 0 at N = N_psr.
  full estimator   -> H2: each pulsar spends one data number on its own L, spread
      over its B bins, so spare data/bin = 2N_psr - N_psr/B and
      ratio ~ 1 - 2n/(N_psr(2 - 1/B)), knee n ~ N_psr(2-1/B)/4 = 29/44/51/54 at
      B = 1/2/4/8.05; same-frequency: linear to 0 at N = N_psr/2.
  C.3 (6 params/source: phase, log h, theta, phi, log f, cos i):
      booked -> knee n ~ N_psr/6; full -> ratio ~1e-5 already at N = 1 (sky
      position enters the pulsar-term phase with lever arm 2 pi f L/c ~ 4800 rad;
      a precision problem, not a count; amplitude-invariant).
"""
import os, sys, time, json
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import prong2_transition as P

OUT = "/home/mattm/projects/HSYMT/PAPER_figures/job48_data"
T = 15 * P.YR; CAD = 7.0; SIG = 1e-7; L10H = -14.3; HTOT = 5e-14
RC = (1e-10, 1e-13)
FB = 1.0 / T
toas = P.build_toas(T, CAD)
t_wall = time.time()

def hr(msg=""):
    print(f"\n{'=' * 78}\n{msg}\n{'=' * 78}", flush=True)

# ----------------------------------------------------------------------------
# 0. CONFIG + REPRODUCE THE BOOKED RATIOS (nothing changes before this passes)
# ----------------------------------------------------------------------------
hr("JOB 48 -- STEP 0: booked config and reproduction")
d = np.load("/home/mattm/projects/HSYMT/CW_transition/prong2_results.npz", allow_pickle=True)
print(f"  N_psr={int(d['ps_n_psr'])}  T={float(d['ps_T'])/P.YR:.1f} yr  cadence={float(d['ps_cadence_days'])} d"
      f"  sigma={float(d['ps_sigma'])*1e9:.0f} ns  band={d['ps_f_band']*1e9} nHz  skies={int(d['ps_n_seeds'])}"
      f"  rcond={float(d['prod_rcond']):.0e}  log10 h={L10H}")
print(f"  n_toa per pulsar = {toas.shape[0]}   1/T = {FB*1e9:.4f} nHz   T*Df(3-20 nHz) = {T*(2e-8-3e-9):.3f} bins")
print(f"  booked knee = {float(d['knee_N']):.1f}  (0.45 N_psr T*Df -> n = {float(d['knee_N'])/(T*(2e-8-3e-9)):.1f} per bin)")
print(f"  Fisher params: 116 distances + (phase0, log10 h) per source; sky, f, psi, cos i fixed")
print(f"  pulsars: build_array(116, seed=777) isotropic synthetic, all L = 1 kpc")
# ----------------------------------------------------------------------------
# machinery: chunked assembly (same sums as P.fisher_blocks) + both estimators
# ----------------------------------------------------------------------------
@jax.jit
def _chunk(toas_, p_c, L_c, srcs, mask, inv_s2):
    g, H = jax.vmap(lambda ph, L: P.analytic_psr_grads(toas_, ph, L, srcs, mask))(p_c, L_c)
    return (jnp.einsum("pt,pt->p", g, g) * inv_s2,
            jnp.einsum("pt,pta->pa", g, H) * inv_s2,
            jnp.einsum("pta,ptb->ab", H, H) * inv_s2)

def chunk_size(n_psr, tier, budget=60000):
    ds = [k for k in range(1, n_psr + 1) if n_psr % k == 0 and k * tier <= budget]
    return max(ds) if ds else 1

def blocks(p_hats, srcs_pad, mask):
    n_psr = p_hats.shape[0]; cs = chunk_size(n_psr, srcs_pad.shape[0])
    L = jnp.ones(n_psr); c_, f_, F = [], [], None
    for a in range(0, n_psr, cs):
        c, f, s = _chunk(toas, p_hats[a:a+cs], L[a:a+cs], srcs_pad, mask, 1.0 / SIG**2)
        c_.append(c); f_.append(f); F = s if F is None else F + s
    return jnp.concatenate(c_), jnp.concatenate(f_), F

_eigh = jax.jit(jnp.linalg.eigh)

def estimators(cond, FLs, Fss):
    """booked marg at each rcond, and the full marginal (all other params incl.
    other distances marginalised) at rcond[0]."""
    w, V = _eigh(Fss)
    U = V.T @ FLs.T
    wmax = jnp.max(w); scale = jnp.where(wmax > 0, wmax, 1.0)
    cond_np = np.asarray(cond); out = {}
    for rc in RC:
        keep = w > rc * scale
        inv = jnp.where(keep, 1.0 / jnp.where(keep, w, 1.0), 0.0)
        Q = np.asarray(U.T @ (inv[:, None] * U))
        S = np.diag(cond_np) - Q
        out[rc] = np.clip(np.diag(S), 0.0, None)
        if rc == RC[0]:
            S = 0.5 * (S + S.T)
            lam, W = np.linalg.eigh(S)
            kp = lam > 1e-10 * max(lam.max(), 1e-300)
            v = (W[:, kp] ** 2 / lam[kp]).sum(1)
            nullw = (W[:, ~kp] ** 2).sum(1)
            full = np.where((nullw > 1e-8) | (v <= 0), 0.0, 1.0 / np.where(v > 0, v, 1.0))
    return cond_np, out[RC[0]], out[RC[1]], full

def draw(n, rng, fspec, mode):
    """Same draw ORDER as P.draw_sources (theta, phi, psi, cos i, frequency, phase0),
    so every arm consumes the random stream identically: arms share geometry,
    phases and within-band quantile u; only the frequency map differs."""
    theta = np.arccos(rng.uniform(-1, 1, n)); phi = rng.uniform(0, 2*np.pi, n)
    psi = rng.uniform(0, np.pi, n); cinc = rng.uniform(-1, 1, n)
    if fspec[0] == "log":
        l10f = rng.uniform(np.log10(fspec[1]), np.log10(fspec[2]), n)
    else:
        u = rng.uniform(0, 1, n)
        f = fspec[1] + u * (fspec[2] - fspec[1]) if fspec[0] == "lin" else np.full(n, fspec[1])
        l10f = np.log10(f)
    phase0 = rng.uniform(0, 2*np.pi, n)
    l10h = np.full(n, L10H) if mode == "fixed_persource" else np.full(n, np.log10(HTOT/np.sqrt(n)))
    return np.stack([theta, phi, psi, cinc, l10f, l10h, phase0], axis=1)

TIERS = [8, 16, 32, 64, 128, 256, 512, 768, 1024, 1536, 2048]
def tier(n):
    return next(t for t in TIERS if t >= n)

# ---- equivalence gates ----------------------------------------------------------
hr("STEP 0b: equivalence of the chunked assembly and the draw to the booked code")
rng_a = np.random.default_rng(10 + 1000*3 + 300); rng_b = np.random.default_rng(10 + 1000*3 + 300)
sa = P.draw_sources(300, rng_a, (3e-9, 2e-8), L10H); sb = draw(300, rng_b, ("log", 3e-9, 2e-8), "fixed_persource")
print(f"  draw('log') vs P.draw_sources, booked seed: max |diff| = {np.max(np.abs(sa-sb)):.1e}")
assert np.array_equal(sa, sb)
ph116 = P.build_array(116)
pad, mask = P.pad_sources(sa, 512)
c0, f0, F0 = P.fisher_blocks(toas, ph116, jnp.ones(116), jnp.array(pad), jnp.array(mask), 116, SIG)
m0 = np.asarray(P.marg_from_blocks(c0, f0, F0, 1e-10))
c1, mb, _, _ = estimators(*blocks(ph116, jnp.array(pad), jnp.array(mask)))
e1 = np.max(np.abs(c1 - np.asarray(c0)) / np.asarray(c0)); e2 = np.max(np.abs(mb - m0) / np.maximum(m0, 1e-300))
print(f"  chunked vs booked fisher_blocks+marg_from_blocks (N=300): cond rel {e1:.1e}, marg rel {e2:.1e}")
assert e1 < 1e-9 and e2 < 1e-6

hr("STEP 0c: reproduce the booked ratios (booked seeds, 30 skies, tiered path)")
ok = True; bk = list(d["ps_n"])
for n in [1, 100, 398, 1000]:
    C, M = [], []
    for sd in range(30):
        rng = np.random.default_rng(10 + 1000 * sd + n)
        srcs = draw(n, rng, ("log", 3e-9, 2e-8), "fixed_persource")
        pad, mask = P.pad_sources(srcs, tier(n))
        c, m10, _, _ = estimators(*blocks(ph116, jnp.array(pad), jnp.array(mask)))
        C.append(np.median(c)); M.append(np.median(m10))
    new = np.median(M) / np.median(C); j = bk.index(n); old = d["ps_marg"][j] / d["ps_cond"][j]
    ok &= abs(new - old) < 1e-4
    print(f"  N={n:5d}  ratio reproduced {new:.4f}  booked {old:.4f}  |diff| {abs(new-old):.1e}", flush=True)
print(f"  reproduction {'PASS' if ok else 'FAIL -- STOP'}")
if not ok:
    sys.exit(3)

# ----------------------------------------------------------------------------
# arms
# ----------------------------------------------------------------------------
def ngrid(k_job):
    """log grid 1 -> 20x the job's predicted knee (5/decade), plus 10/decade
    within [k/6, 3k] so both pre-registered knees (k and ~k/2) sit in the dense
    part. Capped at 2048 for speed (stated in the report)."""
    top = min(20 * k_job, 2048)
    base = np.geomspace(1, top, int(np.ceil(5 * np.log10(top))) + 1)
    lo, hi = max(1.0, k_job / 6), min(top, 3 * k_job)
    dense = np.geomspace(lo, hi, int(np.ceil(8 * np.log10(hi / lo))) + 1)
    return np.unique(np.round(np.concatenate([base, dense])).astype(int))

def nsky(N):
    return 20

def run_arm(name, n_psr, fspec, mode, k_job, B_label, fc, grid=None):
    fp = os.path.join(OUT, f"{name}.npz")
    if os.path.exists(fp):
        print(f"  [{name}] checkpoint exists -- skipped"); return
    ph = P.build_array(n_psr); Ng = ngrid(k_job) if grid is None else np.array(grid)
    hr(f"ARM {name}: N_psr={n_psr}  B={B_label}  f={fspec}  mode={mode}\n  N grid ({len(Ng)} pts): {Ng.tolist()}")
    rec = dict(N=[], s=[], cond=[], mb10=[], mb13=[], mf=[])
    t_arm = time.time()
    for N in Ng:
        t0 = time.time()
        for s in range(nsky(N)):
            rng = np.random.default_rng([48, s, int(N)])
            srcs = draw(int(N), rng, fspec, mode)
            pad, mask = P.pad_sources(srcs, tier(N))
            c, m10, m13, mf = estimators(*blocks(ph, jnp.array(pad), jnp.array(mask)))
            rec["N"].append(N); rec["s"].append(s); rec["cond"].append(c)
            rec["mb10"].append(m10); rec["mb13"].append(m13); rec["mf"].append(mf)
        k = np.array(rec["N"]) == N
        C = np.array(rec["cond"])[k]; M = np.array(rec["mb10"])[k]; MF = np.array(rec["mf"])[k]
        cm = np.median(C, 1); r_b = np.median(M, 1) / cm; r_f = np.median(MF, 1) / cm
        rb_ = np.median(np.median(M, 1)) / np.median(cm)
        print(f"  N={N:5d} n={N/B_label:8.2f} skies={k.sum():2d}  booked {rb_:.4f} "
              f"[{np.percentile(r_b,16):.3f},{np.percentile(r_b,84):.3f}]  full {np.median(r_f):.4f}"
              f"  ({time.time()-t0:5.1f} s)", flush=True)
    np.savez(fp, **{k_: np.array(v) for k_, v in rec.items()}, n_psr=n_psr, B=B_label, fc=fc,
             fspec=json.dumps([str(x) for x in fspec]), mode=mode, k_job=k_job,
             rcond=np.array(RC), runtime=time.time() - t_arm)
    print(f"  [{name}] saved ({(time.time()-t_arm)/60:.1f} min; wall {(time.time()-t_wall)/60:.1f} min)", flush=True)

TB8 = T * (2e-8 - 3e-9)
ARMS = [
    ("A_B1",    116, ("lin", 3*FB, 4*FB),  "fixed_persource",  52, 1.0,  3.5*FB),
    ("A_SF",    116, ("single", 3.5*FB),   "fixed_persource",  52, 1.0,  3.5*FB),
    ("A_B2",    116, ("lin", 3*FB, 5*FB),  "fixed_persource", 104, 2.0,  4.0*FB),
    ("A_B4",    116, ("lin", 2*FB, 6*FB),  "fixed_persource", 208, 4.0,  4.0*FB),
    ("A_B8",    116, ("log", 3e-9, 2e-8),  "fixed_persource", 423, TB8,  np.sqrt(3e-9*2e-8)),
    ("A_B1hi",  116, ("lin", 7*FB, 8*FB),  "fixed_persource",  52, 1.0,  7.5*FB),
    ("B_58",     58, ("lin", 3*FB, 4*FB),  "fixed_persource",  26, 1.0,  3.5*FB),
    ("B_232",   232, ("lin", 3*FB, 4*FB),  "fixed_persource", 104, 1.0,  3.5*FB),
    ("C_tot",   116, ("lin", 3*FB, 4*FB),  "fixed_total",      52, 1.0,  3.5*FB),
]
hr("BANDS (bin edges k/T, 1/T = %.4f nHz)" % (FB*1e9))
for a in ARMS:
    fs = a[2]
    rng_txt = (f"{fs[1]*1e9:.3f}-{fs[2]*1e9:.3f} nHz" if fs[0] != "single" else f"all at {fs[1]*1e9:.3f} nHz")
    print(f"  {a[0]:7s} N_psr={a[1]:3d} B={a[5]:5.2f} {fs[0]:6s} {rng_txt}  mode={a[3]}")
for a in ARMS:
    run_arm(*a)

# ----------------------------------------------------------------------------
# C.3 optional: 6 free parameters per source (autodiff per source; block-sparse)
# ----------------------------------------------------------------------------
FREE6 = jnp.array([6, 5, 0, 1, 4, 3])     # phase0, log10 h, theta, phi, log10 f, cos i
def _src_jac(toas_, p_hat, L, src):
    f = lambda L_, s_: P.single_source_residual(toas_, p_hat, L_, s_)
    return jax.jacfwd(f, argnums=0)(L, src), jax.jacfwd(f, argnums=1)(L, src)

def make_blocks_free(free):
    @jax.jit
    def _b(toas_, p_c, L_c, srcs, mask, inv_s2):
        def one(ph, L):
            jL, js = jax.vmap(lambda s: _src_jac(toas_, ph, L, s))(srcs)
            g = jnp.sum(mask[:, None] * jL, 0)
            H = (mask[:, None, None] * js[:, :, free]).transpose(1, 0, 2).reshape(toas_.shape[0], -1)
            return g, H
        g, H = jax.vmap(one)(p_c, L_c)
        return (jnp.einsum("pt,pt->p", g, g) * inv_s2,
                jnp.einsum("pt,pta->pa", g, H) * inv_s2,
                jnp.einsum("pta,ptb->ab", H, H) * inv_s2)
    return _b

def blocks_free(bfun, p_hats, srcs_pad, mask, cs=29):
    n_psr = p_hats.shape[0]; L = jnp.ones(n_psr); c_, f_, F = [], [], None
    for a in range(0, n_psr, cs):
        c, f, s = bfun(toas, p_hats[a:a+cs], L[a:a+cs], srcs_pad, mask, 1.0/SIG**2)
        c_.append(c); f_.append(f); F = s if F is None else F + s
    return jnp.concatenate(c_), jnp.concatenate(f_), F

fp = os.path.join(OUT, "C3_6p.npz")
if not os.path.exists(fp):
    hr("C.3 (optional): N_psr=58, B=1, 6 free params per source")
    ph58 = P.build_array(58)
    b2 = make_blocks_free(jnp.array([6, 5]))
    rng = np.random.default_rng([48, 0, 20]); srcs = draw(20, rng, ("lin", 3*FB, 4*FB), "fixed_persource")
    pad, mask = P.pad_sources(srcs, 32)
    ca, _, _, _ = estimators(*blocks(ph58, jnp.array(pad), jnp.array(mask)))
    _, ma, _, _ = estimators(*blocks(ph58, jnp.array(pad), jnp.array(mask)))
    cb, mbb, _, _ = estimators(*blocks_free(b2, ph58, jnp.array(pad), jnp.array(mask)))
    print(f"  autodiff 2-param route vs analytic (N=20): cond rel {np.max(np.abs(cb-ca)/ca):.1e}, "
          f"marg rel {np.max(np.abs(mbb-ma)/np.maximum(ma,1e-300)):.1e}")
    b6 = make_blocks_free(FREE6)
    Ng = np.array([1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 17, 20, 25, 30, 40, 60, 100, 200])
    rec = dict(N=[], s=[], cond=[], mb10=[], mb13=[], mf=[])
    t_arm = time.time()
    for N in Ng:
        t0 = time.time()
        for s in range(20):
            rng = np.random.default_rng([48, s, int(N)])
            srcs = draw(int(N), rng, ("lin", 3*FB, 4*FB), "fixed_persource")
            pad, mask = P.pad_sources(srcs, tier(N))
            c, m10, m13, mf = estimators(*blocks_free(b6, ph58, jnp.array(pad), jnp.array(mask)))
            rec["N"].append(N); rec["s"].append(s); rec["cond"].append(c)
            rec["mb10"].append(m10); rec["mb13"].append(m13); rec["mf"].append(mf)
        k = np.array(rec["N"]) == N
        C = np.array(rec["cond"])[k]; cm = np.median(C, 1)
        rb_ = np.median(np.median(np.array(rec["mb10"])[k], 1)) / np.median(cm)
        rb13 = np.median(np.median(np.array(rec["mb13"])[k], 1)) / np.median(cm)
        rf = np.median(np.median(np.array(rec["mf"])[k], 1) / cm)
        print(f"  N={N:4d}  booked {rb_:.4f} (rcond 1e-13: {rb13:.4f})  full {rf:.2e}  ({time.time()-t0:.1f} s)", flush=True)
    np.savez(fp, **{k_: np.array(v) for k_, v in rec.items()}, n_psr=58, B=1.0, fc=3.5*FB,
             fspec=json.dumps(["lin", str(3*FB), str(4*FB)]), mode="fixed_persource_6param",
             k_job=58/6, rcond=np.array(RC), runtime=time.time()-t_arm)

print(f"\n  TOTAL COMPUTE RUNTIME: {(time.time()-t_wall)/60:.1f} min", flush=True)
