"""JOB 54 tests: residual vs prong2, closed-form likelihood vs brute force, evolution size."""
import numpy as np, jax, jax.numpy as jnp, json, os
import job54_loc_core as C
p2 = C.p2
arr = C.build_all(); src = C.draw_sources()
toas = jnp.asarray(arr["toas"]); out = {}

def test_residual_vs_prong2(n=50):
    rng = np.random.default_rng(1); worst = 0.0
    for k in range(n):
        th = np.arccos(rng.uniform(-1, 1)); ph = rng.uniform(0, 2*np.pi); psi = rng.uniform(0, np.pi)
        ci = rng.uniform(-1, 1); f = 10**rng.uniform(np.log10(3e-9), np.log10(2e-8)); h = 10**rng.uniform(-15, -13)
        p0 = rng.uniform(0, 2*np.pi); a = rng.integers(0, C.N_PSR); L = rng.uniform(100, 5000)
        ref = p2.single_source_residual(toas, jnp.asarray(arr["p_hats"][a]), L/1000.0,
                                        jnp.array([th, ph, psi, ci, np.log10(f), np.log10(h), p0]))
        mine = C.residual(toas, jnp.asarray(arr["p_hats"][a]), L, th, ph, psi, ci, f, h, p0, 1e-6)
        worst = max(worst, float(jnp.max(jnp.abs(ref-mine))/jnp.max(jnp.abs(ref))))
    return worst

def test_closed_form(f=1e-8, snr=10.0, isrc=3, npts=40):
    s = C.src_at(src, isrc); mc = 1e9
    h = snr*C.snr_unit_h(arr, s, f, mc)
    sky = C.make_sky(arr, s, f, h, mc)
    rng = np.random.default_rng(2)
    u = jnp.asarray(np.concatenate([[0.0], rng.normal(0, 0.05, npts)])); v = jnp.asarray(np.concatenate([[0.0], rng.normal(0, 0.05, npts)]))
    coh, inc, a, alpha, sphi1, base = sky["parts"](u, v)
    nv = np.asarray(sky["uv2n"](u, v))
    th = np.arccos(nv[:, 2]); ph = np.arctan2(nv[:, 1], nv[:, 0])
    Pj = jnp.asarray(arr["p_hats"]); Lj = jnp.asarray(arr["L_pc"])
    d = C._res_arr(toas, Pj, Lj, s["theta"], s["phi"], s["psi"], s["cinc"], f, h, s["phase0"], mc)
    e_coh = 0; e_inc = 0; e_inc_approx = 0; e_part = 0
    phis = np.linspace(0, 2*np.pi, 256, endpoint=False)
    for k in range(npts+1):
        t = C._res_arr(toas, Pj, Lj, float(th[k]), float(ph[k]), s["psi"], s["cinc"], f, h, s["phase0"], mc)
        bf = -0.5*np.asarray(jnp.sum((d-t)**2, 1))/C.SIGMA**2
        e_coh = max(e_coh, np.max(np.abs(bf-np.asarray(coh[k]))))
        if k % 8 == 0:
            # exact phase marginal by quadrature: shift pulsar-term phase by varying L by dphi/sphi1
            sp = np.asarray(sphi1[k])
            acc = []
            for p in phis:
                tt = C._res_arr(toas, Pj, Lj + p/jnp.asarray(sp), float(th[k]), float(ph[k]), s["psi"], s["cinc"], f, h, s["phase0"], mc)
                acc.append(-0.5*np.asarray(jnp.sum((d-tt)**2, 1))/C.SIGMA**2)
            acc = np.array(acc); mx = acc.max(0)
            ex = mx + np.log(np.mean(np.exp(acc-mx), 0))
            good = sp > 0.5      # need a finite distance shift for a 2 pi phase turn
            e_inc = max(e_inc, np.max(np.abs(ex-np.asarray(inc[k]))[good]))
            # partial: Gaussian in L with sigma_L = 0.3 pc, numeric
            sl = 0.3; xs = np.linspace(-6, 6, 241); w = np.exp(-0.5*xs**2); w /= w.sum(); acc = []
            for xx in xs:
                tt = C._res_arr(toas, Pj, Lj + xx*sl, float(th[k]), float(ph[k]), s["psi"], s["cinc"], f, h, s["phase0"], mc)
                acc.append(-0.5*np.asarray(jnp.sum((d-tt)**2, 1))/C.SIGMA**2)
            acc = np.array(acc); mx = acc.max(0)
            ex = mx + np.log(np.sum(w[:, None]*np.exp(acc-mx), 0))
            pl = np.asarray(C.partial_lnl(base[k], a[k], alpha[k], sl*sphi1[k]))
            e_part = max(e_part, np.max(np.abs(ex-pl)))
    return dict(max_abs_err_coh=float(e_coh), max_abs_err_inc_vs_exact_phase_marginal=float(e_inc),
                max_abs_err_partial_0p3pc=float(e_part), coh_at_truth_max=float(np.max(np.abs(np.asarray(coh[0])))),
                max_a=float(jnp.max(a)))

def evolution_size():
    r = {}
    for f in (1e-8, 2e-8):
        tau = 1000.0*C.PC_OVER_C
        fr, dPhi = C.evol(f, 1e9*C.GMSUN_C3, tau)
        fdot = (96/5)*np.pi**(8/3)*(1e9*C.GMSUN_C3)**(5/3)*f**(11/3)
        r[f"{f*1e9:.0f}nHz"] = dict(fdot_Hz_per_s=float(fdot), fp_over_f_at_tau_1kpc=float(fr),
            df_nHz=float(f*(1-fr)*1e9), df_times_T_cycles=float(f*(1-fr)*C.T_SPAN),
            extra_phase_rad=float(dPhi-2*np.pi*f*tau), earth_drift_phase_over_span_rad=float(np.pi*fdot*C.T_SPAN**2))
    return r

out["residual_vs_prong2_max_rel"] = test_residual_vs_prong2()
print(out, flush=True)
out["evolution"] = evolution_size(); print(out["evolution"], flush=True)
for f, snr in ((1e-8, 10.0), (2e-8, 10.0), (2e-8, 20.0)):
    out[f"closed_form_f{f*1e9:.0f}_snr{snr:.0f}"] = test_closed_form(f, snr); print(out, flush=True)
json.dump(out, open(os.path.join(C.DATA, "tests.json"), "w"), indent=1)
