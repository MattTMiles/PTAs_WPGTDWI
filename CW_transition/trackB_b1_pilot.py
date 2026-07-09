"""Track B — B1 PILOT: measure the objects that size (and possibly re-scope) B1.

Motivation. B1's brief specifies the targeted pipeline as "sky FIXED at truth per loud
source, FREQUENCY free (1-D scan to needle resolution), extrinsics fit". That presumes the
only registration axis left once the sky is given is log10_fgw. But the model's pulsar term
is CHIRPED: the phase-connected pulsar-term phase carries `pi * fdot * tau_p^2` with
tau_p = L_p (1 - cos mu_p)/c ~ kyr, and fdot ∝ Mc^{5/3}. F2's lever-arm ladder
(`trackB_F2_ladder.py:47`) used the NON-chirped phase `2 pi f L (1-cos mu)/c`, whose
gradient wrt log10_mc is IDENTICALLY ZERO -- so F2 could not have seen an Mc lane even if
one exists. This pilot tests, on the real likelihood, whether log10_mc is a THIRD
registration axis of tolerance comparable to sky/freq, and whether the Earth term (which is
all that constrains it before the pulsar terms register) can supply it.

Deliverables of the pilot (all zero-noise Asimov, pop draw, truth-placed = LABELLED):
  M1  chirp-phase census: pi fdot tau_p^2 per (pulsar, loud source), and d(phase)/dlog10_mc.
  M2  DECISIVE, likelihood-side: E-step census P(true fringe) vs a perturbation of ONE
      registration axis at a time (sky-only / freq-only / mc-only). Gives each axis its own
      certification tolerance, directly comparable to B0.2's ~1e-4 scaled.
  M3  Earth-term-only Fisher: sigma(log10_fgw) and sigma(log10_mc) from the Earth term alone
      (the information a targeted pipeline actually has before registration). Ratio to M2's
      tolerance = the per-axis wall, the targeted-scenario analogue of F2's 27x.
  M4  fringe-shift structure: at a small mc offset, does each pulsar's MAP fringe move
      proportional to tau_p^2 (registration) rather than randomly (noise)?
  M5  cost model: wall time of one B1Split E-step; sizes every B1 scan.

Env: jug-gpu, autotune off, x64. Matt commits; never commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_b1_core as C
import trackB_estimator as TE

CWT = "/home/mattm/projects/HSYMT/CW_transition"
C_KPC_PER_S = 9.71561e-12          # speed of light, kpc/s  (matches trackB_F2_ladder)
MSUN_S = 4.925490947e-6            # G M_sun / c^3, seconds

# theta layout inside a source block
I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)


def fdot_of(log10_mc, log10_fgw):
    """Quadrupole chirp rate, Hz/s.  fdot = (96/5) pi^{8/3} (G Mc/c^3)^{5/3} f^{11/3}.
    Diagnostic/reporting only; every TOLERANCE below is autodiffed off the EXACT phase."""
    mc = 10.0 ** log10_mc * MSUN_S
    f = 10.0 ** log10_fgw
    return (96.0 / 5.0) * np.pi ** (8.0 / 3.0) * mc ** (5.0 / 3.0) * f ** (11.0 / 3.0)


def tau_of(P, k, a):
    """Pulsar-term light-travel lag tau_p = L_p (1 - cos mu_pa)/c, seconds."""
    src = P["theta_truth"][P["n_dist"]:]
    cg = src[8 * k + I_COSGT]; gp = src[8 * k + I_GWPHI]
    st = np.sqrt(max(1 - cg ** 2, 0.0))
    khat = np.array([st * np.cos(gp), st * np.sin(gp), cg])
    cosmu = float(np.dot(khat, P["disco_psrs"][a].pos))
    return (P["L0"][a] / C_KPC_PER_S) * (1.0 - cosmu)


# ---- EXACT pulsar-term phase, mirroring discovery.deterministic.make_phase_connected_binary ----
from discovery import const as _dsc


def _khat(cg, gp):
    st = jnp.sqrt(jnp.clip(1 - cg ** 2, 0, 1))
    return jnp.array([st * jnp.cos(gp), st * jnp.sin(gp), cg])


def dphi_exact(params4, pos, L_kpc, t_rel):
    """GW phase of the PULSAR term MINUS that of the EARTH term, radians, at time t_rel.

    Exactly the phase discovery uses: phase(t) = phase0_orb + (1/(32 mc^{5/3}))
    (w0^{-5/3} - omega(t)^{-5/3}), omega(t) = w0 (1 - (256/5) mc^{5/3} w0^{8/3} t)^{-3/8},
    evaluated at t = t_rel (Earth) and t = t_rel - (kpc/c) L (1-cos mu) (pulsar).
    GW phase = 2 * orbital phase. phase0 cancels in the difference.
    """
    cg, gp, lmc, lf = params4
    mc = (10.0 ** lmc) * _dsc.Tsun
    w0 = jnp.pi * (10.0 ** lf)
    cosmu = jnp.dot(_khat(cg, gp), pos)
    tp = t_rel - (_dsc.kpc / _dsc.c) * L_kpc * (1.0 - cosmu)

    def orb_phase(t):
        term = 1.0 - (256.0 / 5.0) * mc ** (5.0 / 3.0) * w0 ** (8.0 / 3.0) * t
        omega = w0 * jnp.power(term, -3.0 / 8.0)
        return (1.0 / (32.0 * mc ** (5.0 / 3.0))) * (w0 ** (-5.0 / 3.0) - omega ** (-5.0 / 3.0))

    return 2.0 * (orb_phase(tp) - orb_phase(t_rel))


_grad_dphi = jax.jit(jax.grad(dphi_exact, argnums=0))


def dphi_linear(params4, pos, L_kpc):
    """F2's non-chirped object: Phi = 2 pi f L (1-cos mu)/c. No log10_mc dependence."""
    cg, gp, _, lf = params4
    cosmu = float(np.dot(np.asarray(_khat(cg, gp)), pos))
    return 2 * np.pi * (10.0 ** lf) * (_dsc.kpc / _dsc.c) * L_kpc * (1.0 - cosmu)


# ============================================================
def m1_chirp_census(P, t_rel):
    """Per (pulsar, loud) 1-rad registration tolerance in EVERY source param that shifts the
    pulsar-term phase: sky (2), log10_mc, log10_fgw. Autodiff of the EXACT phase."""
    print("\n=== M1: exact pulsar-term phase + per-axis registration ladder ===", flush=True)
    src = P["theta_truth"][P["n_dist"]:]
    scale = TE.phi_scale(P)
    npsr = P["npsr"]
    rows = []
    for k in range(C.N_LOUD):
        lmc = src[8 * k + I_MC]; lf = src[8 * k + I_FGW]
        p4 = jnp.array([src[8 * k + I_COSGT], src[8 * k + I_GWPHI], lmc, lf])
        fd = fdot_of(lmc, lf)
        print(f"  loud{k}: log10_mc={lmc:.4f} log10_fgw={lf:.4f} (f={10**lf*1e9:.2f} nHz) "
              f"fdot={fd:.3e} Hz/s  log10_h={src[8*k+I_H]:.3f}", flush=True)
        for a in range(npsr):
            pos = jnp.asarray(P["disco_psrs"][a].pos)
            L = float(P["L0"][a])
            dp = float(dphi_exact(p4, pos, L, t_rel))              # exact pterm-earth GW phase
            g = np.asarray(_grad_dphi(p4, pos, L, t_rel))          # d(dphi)/d(cg,gp,lmc,lf)
            lin = dphi_linear(np.asarray(p4), np.asarray(pos), L)  # F2's non-chirped phase
            gs = np.hypot(g[0] * scale[I_COSGT], g[1] * scale[I_GWPHI])   # rad per scaled sky
            g_mc = abs(g[2] * scale[I_MC])                                # rad per scaled mc
            g_f = abs(g[3] * scale[I_FGW])                                # rad per scaled freq
            rows.append((k, a, tau_of(P, k, a), dp, lin, gs, g_mc, g_f))
    R = np.array(rows)
    tau, dphi, lin = R[:, 2], R[:, 3], R[:, 4]
    with np.errstate(divide="ignore"):
        tol_sky = 1.0 / R[:, 5]; tol_mc = 1.0 / R[:, 6]; tol_f = 1.0 / R[:, 7]
    C._finite("M1 tolerances", tol_sky[np.isfinite(tol_sky)])
    print(f"\n  tau_p (kyr): median {np.median(tau)/3.156e10:.2f}  min {tau.min()/3.156e10:.3f}  "
          f"max {tau.max()/3.156e10:.1f}", flush=True)
    print(f"  |exact pterm-earth GW phase| (rad): median {np.median(np.abs(dphi)):.3e}  "
          f"max {np.abs(dphi).max():.3e}", flush=True)
    print(f"  |non-chirped 2 pi f tau|     (rad): median {np.median(lin):.3e}  max {lin.max():.3e}",
          flush=True)
    print(f"  chirp correction |exact|/|linear|: median {np.median(np.abs(dphi)/lin):.3f}", flush=True)
    print(f"\n  1-rad tolerance ladders, SCALED units (loosest / median / tightest):", flush=True)
    for lab, t in [("sky", tol_sky), ("log10_mc", tol_mc), ("log10_fgw", tol_f)]:
        tt_ = t[np.isfinite(t)]
        print(f"    {lab:10s} {tt_.max():.3e} / {np.median(tt_):.3e} / {tt_.min():.3e}   "
              f"(# >= 1e-2: {int(np.sum(tt_>=1e-2))}, # >= 1e-3: {int(np.sum(tt_>=1e-3))})", flush=True)
    tol_f_nochirp = 1.0 / (np.log(10.0) * lin * scale[I_FGW])
    print(f"    {'freq(F2)':10s} {tol_f_nochirp.max():.3e} / {np.median(tol_f_nochirp):.3e} / "
          f"{tol_f_nochirp.min():.3e}   <- F2's chirp-free object "
          f"(reported 2.52e-2 / 1.04e-4)", flush=True)
    print(f"  => including the chirp TIGHTENS the freq ladder by a median factor "
          f"{np.median(tol_f_nochirp)/np.median(tol_f[np.isfinite(tol_f)]):.2f}x", flush=True)
    np.savez(f"{CWT}/b1_pilot_m1.npz", rows=R, tol_sky=tol_sky, tol_mc=tol_mc, tol_f=tol_f,
             tol_f_nochirp=tol_f_nochirp, names=np.array(P["names"]), t_rel=t_rel)
    return dict(tol_sky=tol_sky, tol_mc=tol_mc, tol_f=tol_f)


# ============================================================
def _estep_ceilings(P, sp, theta, data, FT, n_true):
    LNL = sp.estep(theta, P["EV_truth"], P["AI"], data, jnp.asarray(P["mask_one"]))
    C._finite("estep LNL", LNL)
    post = FT.posterior(LNL)
    pt, rpost, rlik = FT.p_of_fringe(n_true)
    return post, pt, LNL


def m2_axis_tolerances(P, sp, FT):
    print("\n=== M2: per-AXIS certification tolerance (E-step only, zero noise, truth-placed) ===",
          flush=True)
    tt = P["theta_truth"]; nd = P["n_dist"]
    data0 = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(P["mask_one"]))
    scale = TE.phi_scale(P)
    n_true = np.zeros(P["npsr"], int)
    rng = np.random.default_rng(1234)
    # fixed unit perturbation directions (same across deltas so the scan is a clean ray)
    dirs = {}
    for lab, idxs in [("sky", [I_COSGT, I_GWPHI]), ("freq", [I_FGW]), ("mc", [I_MC]),
                      ("extrinsic", [I_COSINC, I_H, I_PH0, I_PSI])]:
        v = np.zeros(8 * C.N_LOUD)
        for k in range(C.N_LOUD):
            for j in idxs:
                v[8 * k + j] = rng.normal()
        dirs[lab] = v

    print(f"  {'axis':10s} {'delta(scaled)':>13s} {'J0711':>7s} {'J1713':>7s} {'J1909':>7s} "
          f"{'J0437':>7s} {'ncertT':>7s} {'ncertF':>7s}", flush=True)
    out = {}
    deltas = [0.0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
    for lab, v in dirs.items():
        rec = []
        for d in deltas:
            th = tt.copy()
            th[nd:nd + 8 * C.N_LOUD] += v * scale[:8 * C.N_LOUD] * d
            post, pt, _ = _estep_ceilings(P, sp, th, data0, FT, n_true)
            cert = post["q_max"] > 0.9
            nT = int(np.sum(cert & (post["map_k"] == 0)))
            nF = int(np.sum(cert & (post["map_k"] != 0)))
            vals = [pt[P["names"].index(nm)] for nm in C.CENSUS] + [pt[P["anchor_idx"]]]
            rec.append((d, *vals, nT, nF))
            print(f"  {lab:10s} {d:13.1e} {vals[0]:7.3f} {vals[1]:7.3f} {vals[2]:7.3f} "
                  f"{vals[3]:7.3f} {nT:7d} {nF:7d}", flush=True)
            if lab == "extrinsic" and d >= 1e-3 and vals[1] > 0.9:
                pass
        out[lab] = np.array(rec)
    np.savez(f"{CWT}/b1_pilot_m2.npz", **{k: v for k, v in out.items()}, deltas=np.array(deltas))
    return out


# ============================================================
def m3_earth_fisher(P):
    print("\n=== M3: Earth-term-only Fisher on the loud sources (what a targeted pipeline has) ===",
          flush=True)
    amo = P["amo"]; tt = P["theta_truth"]; nd = P["n_dist"]
    zero = jnp.asarray(P["mask_zero"])
    # Asimov data for the EARTH-ONLY model -> residual(truth)=0 -> H = -Fisher exactly
    data_e = amo["inject_delay"](jnp.asarray(tt), zero)
    ll = float(amo["logL"](jnp.asarray(tt), data_e, zero))
    print(f"  Earth-only Asimov lnL(truth) = {ll:.4f}", flush=True)

    nloud = 8 * C.N_LOUD
    sel = np.arange(nd, nd + nloud)

    def f_of_x(x):
        th = jnp.asarray(tt).at[jnp.asarray(sel)].set(x)
        return amo["logL"](th, data_e, zero)

    # HVP assembly (D1's production path): jax.hessian's forward-over-reverse transient blows
    # up the compile on the 116-pulsar masked graph (observed: 22 GiB RSS, no progress).
    g_of_x = jax.grad(f_of_x)

    def _hvp(x, v):
        return jax.jvp(g_of_x, (x,), (v,))[1]

    hvp = jax.jit(jax.vmap(_hvp, in_axes=(None, 0)))
    x0 = jnp.asarray(tt[sel])
    t0 = time.time()
    cols = []
    CH = 6
    for c0 in range(0, nloud, CH):
        V = jnp.asarray(np.eye(nloud)[c0:c0 + CH])
        cols.append(np.asarray(hvp(x0, V)))
        print(f"    HVP cols {c0}-{c0+V.shape[0]-1} ({time.time()-t0:.0f}s)", flush=True)
    H = np.concatenate(cols, axis=0)
    C._finite("earth Hessian", H)
    # symmetry residual is the honest accuracy check on the assembled Hessian
    print(f"  Hessian symmetry residual max|H-H^T|/max|H| = "
          f"{np.max(np.abs(H-H.T))/np.max(np.abs(H)):.2e}", flush=True)
    F = -0.5 * (H + H.T)
    scale = TE.phi_scale(P)[:nloud]
    # scale to the standard scaled units so numbers are directly comparable to M1/M2
    Fs = F * scale[:, None] * scale[None, :]
    w = np.linalg.eigvalsh(Fs)
    print(f"  scaled Earth Fisher eigenvalues: min {w.min():.3e} max {w.max():.3e} "
          f"(cond {w.max()/max(w.min(),1e-30):.2e})", flush=True)

    # Generative-prior precision, in the SAME scaled units. NEVER pinv here: an axis the
    # Earth term does not constrain has a ZERO Fisher eigenvalue, and pinv would silently
    # report sigma = 0 for it (the exact failure mode for log10_mc). The honest quantity is
    # the POSTERIOR sigma under the model's own uniform priors, whose variance is
    # (range/scale)^2 / 12.
    lo, hi = TE.phi_bounds(P)
    var_prior = ((hi[:nloud] - lo[:nloud]) / scale) ** 2 / 12.0
    Pi = np.diag(1.0 / var_prior)
    Cov = np.linalg.inv(Fs + Pi)
    sig = np.sqrt(np.clip(np.diag(Cov), 0, None))            # posterior sigma, SCALED
    sig_prior = np.sqrt(var_prior)                            # what the prior alone gives
    sig_cond = 1.0 / np.sqrt(np.clip(np.diag(Fs), 1e-300, None))   # conditional (others fixed)
    lab8 = ["cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h", "phase0", "psi"]
    print(f"\n  1-sigma from the EARTH TERM + generative prior (scaled units; 24 loud params free).")
    print(f"  'prior' = prior alone; a posterior sigma == prior sigma means the Earth term "
          f"carries NO information on that axis.", flush=True)
    print(f"  {'param':13s} {'prior':>10s} " +
          " ".join(f"{'loud'+str(k)+'(post)':>14s}" for k in range(C.N_LOUD)), flush=True)
    for j, nm in enumerate(lab8):
        print(f"  {nm:13s} {sig_prior[j]:10.3e} " +
              " ".join(f"{sig[8*k+j]:14.3e}" for k in range(C.N_LOUD)), flush=True)
    print(f"\n  information gain over the prior (sigma_prior / sigma_post):", flush=True)
    for j, nm in enumerate(lab8):
        print(f"  {nm:13s} " + " ".join(f"{sig_prior[j]/sig[8*k+j]:14.2f}" for k in range(C.N_LOUD)),
              flush=True)
    print(f"\n  in PHYSICAL units (posterior):")
    for k in range(C.N_LOUD):
        print(f"    loud{k}: sigma(log10_fgw) = {sig[8*k+I_FGW]*scale[I_FGW]:.3e} dex ; "
              f"sigma(log10_mc) = {sig[8*k+I_MC]*scale[I_MC]:.3e} dex", flush=True)
    np.savez(f"{CWT}/b1_pilot_m3.npz", H=H, Fs=Fs, sigma_scaled=sig, sigma_prior=sig_prior,
             sigma_cond=sig_cond, scale=scale, var_prior=var_prior)
    return sig


# ============================================================
def m4_fringe_shift_structure(P, sp, FT, tau_of_psr, t_rel, delta_mc=3e-4):
    print(f"\n=== M4: fringe-shift structure under an mc-only offset (delta={delta_mc:.0e} scaled) ===",
          flush=True)
    tt = P["theta_truth"]; nd = P["n_dist"]
    scale = TE.phi_scale(P)
    data0 = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(P["mask_one"]))
    n_true = np.zeros(P["npsr"], int)
    th = tt.copy()
    for k in range(C.N_LOUD):
        th[nd + 8 * k + I_MC] += scale[I_MC] * delta_mc
    post, pt, _ = _estep_ceilings(P, sp, th, data0, FT, n_true)
    k_map = post["map_k"]
    # PREDICTED fringe shift: the mc offset moves pulsar p's pterm phase by d(dphi)/dlog10_mc
    # x delta_mc x scale; the fringe that restores it is n = -dphi/(2 pi). Take the source
    # with the tightest dL (the one that sets the fringe grid), via the largest |dphi|.
    src = tt[nd:]
    pred = np.zeros(P["npsr"])
    for a in range(P["npsr"]):
        pos = jnp.asarray(P["disco_psrs"][a].pos); L = float(P["L0"][a])
        best = 0.0
        for k in range(C.N_LOUD):
            p4 = jnp.array([src[8*k+I_COSGT], src[8*k+I_GWPHI], src[8*k+I_MC], src[8*k+I_FGW]])
            g = float(np.asarray(_grad_dphi(p4, pos, L, t_rel))[2])
            dphi = g * scale[I_MC] * delta_mc
            if abs(dphi) > abs(best):
                best = dphi
        pred[a] = -best / (2 * np.pi)
    good = np.abs(k_map) > 0
    if good.sum() >= 5:
        r = float(np.corrcoef(pred[good], k_map[good])[0, 1])
        sl = float(np.polyfit(pred[good], k_map[good], 1)[0])
    else:
        r, sl = np.nan, np.nan
    print(f"  pulsars with MAP fringe != true: {int(np.sum(k_map!=0))}/{P['npsr']}", flush=True)
    print(f"  corr(predicted chirp-driven shift, observed MAP fringe) = {r:.3f}; slope {sl:.3f} "
          f"(expect ~1 if the mc error re-registers via the chirp)", flush=True)
    print(f"  {'pulsar':12s} {'tau(kyr)':>9s} {'pred n':>9s} {'MAP n':>7s} {'q_max':>7s} {'P(true)':>8s}",
          flush=True)
    for nm in C.CENSUS + [C.ANCHOR, "J2222-0137"]:
        a = P["names"].index(nm)
        print(f"  {nm:12s} {max(tau_of_psr[k][a] for k in range(C.N_LOUD))/3.156e10:9.2f} "
              f"{pred[a]:9.2f} {k_map[a]:7d} {post['q_max'][a]:7.3f} {pt[a]:8.4f}", flush=True)
    np.savez(f"{CWT}/b1_pilot_m4.npz", pred=pred, k_map=k_map, q_max=post["q_max"], p_true=pt,
             corr=r, slope=sl, delta_mc=delta_mc)
    return r


# ============================================================
def m5_cost(P, sp):
    print("\n=== M5: cost model ===", flush=True)
    tt = P["theta_truth"]
    data0 = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(P["mask_one"]))
    one = jnp.asarray(P["mask_one"])
    _ = sp.estep(tt, P["EV_truth"], P["AI"], data0, one)          # warm
    t0 = time.time(); _ = sp.estep(tt, P["EV_truth"], P["AI"], data0, one); t_es = time.time() - t0
    t0 = time.time()
    for _ in range(5):
        _ = float(P["amo"]["logL"](jnp.asarray(tt), data0, one))
    t_ln = (time.time() - t0) / 5
    ntoa = [len(p.toas) for p in P["disco_psrs"]]
    print(f"  n_toa: total {sum(ntoa)}, median {int(np.median(ntoa))}, max {max(ntoa)}", flush=True)
    print(f"  full-array logL       : {t_ln*1e3:.2f} ms", flush=True)
    print(f"  B1Split E-step (116x{C.B_CERT}) : {t_es:.3f} s "
          f"({t_es/(P['npsr']*C.B_CERT)*1e6:.2f} us per pulsar-fringe)", flush=True)
    print(f"  -> a 1-D lnL_marg scan of N points costs N x {t_es:.2f} s", flush=True)
    return t_es


def main():
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    P = C.build_b1_problem()
    sp = P["sp"]
    FT = C.FringeTables(P, P["EV_truth"], P["dL_truth"])
    tref = 86400.0 * 51544.5
    t_rel = float(np.mean(np.concatenate([np.asarray(p.toas) for p in P["disco_psrs"]])) - tref)
    print(f"  reference epoch t_rel = {t_rel:.4e} s from J2000 (array mean TOA)", flush=True)
    t_es = m5_cost(P, sp)
    m1 = m1_chirp_census(P, t_rel)
    sig = m3_earth_fisher(P)
    m2 = m2_axis_tolerances(P, sp, FT)
    tau_of_psr = [np.array([tau_of(P, k, a) for a in range(P["npsr"])]) for k in range(C.N_LOUD)]
    m4_fringe_shift_structure(P, sp, FT, tau_of_psr, t_rel)

    # ---- the targeted-scenario wall, per axis ----
    print("\n=== PILOT SUMMARY: the targeted-scenario wall, per axis ===", flush=True)
    scale = TE.phi_scale(P)
    print(f"  (all in SCALED units; wall = Earth-term sigma / certification tolerance)", flush=True)
    for j, nm, arr in [(I_FGW, "log10_fgw", m2["freq"]), (I_MC, "log10_mc", m2["mc"])]:
        # certification tolerance = largest delta keeping >=2 census pulsars above P(true)>0.9
        keep = arr[(arr[:, 1] > 0.9).astype(int) + (arr[:, 2] > 0.9).astype(int)
                   + (arr[:, 3] > 0.9).astype(int) >= 2]
        tol = keep[:, 0].max() if len(keep) else 0.0
        se = float(np.median([sig[8 * k + j] for k in range(C.N_LOUD)]))
        print(f"  {nm:11s} cert tol {tol:.1e} scaled ; Earth sigma {se:.3e} scaled ; "
              f"WALL {se/max(tol,1e-12):.3g}x", flush=True)
    print(f"\n  1-D scan feasibility: a needle-resolution scan of ONE axis over +-4 Earth-sigma "
          f"costs ~{8*np.median([sig[8*k+I_FGW] for k in range(C.N_LOUD)])/1e-5:.0f} x {t_es:.2f} s "
          f"per source (at 1e-5-scaled steps).", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
