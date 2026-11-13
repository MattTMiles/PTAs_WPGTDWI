# ================= JOB 54 Part 1A-1C: one-pulsar examples for the talk ============================
# natural sky, census priors, reference configuration (retardation on, drawn masses, timing x1, 45 yr).
# Pulsar 1 (rule): first natural sky (100); among pulsars with 100 <= K <= 300 (finest tooth of the N = 12
#   known set) that the N = 12 joint fit gets wrong, the one with the median margin.
# Pulsar 2 (rule): first natural sky in order that has a pulsar with prior width 1-3 pc which the N = 12 fit gets
#   right; among those, the median margin.
import os, sys, json, pickle
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job54_common as H, job54_run as R
os.makedirs("talk_data", exist_ok=True)
P53 = lambda j: pickle.load(open(f"job54_data/parts/mainr_{j:03d}.pkl", "rb"))[0]      # corrected run (random pulsar-term phases)
def median_pick(c):          # median margin; lower middle when the count is even
    c = sorted(c, key=lambda r: r["m|ref|12|t1"]); return c[(len(c)-1)//2]
r1 = [r for r in P53(100) if 100 <= r["Kfin|12"] <= 300 and r["w|ref|12|t1"]]
p1 = median_pick(r1); sel = [(100, p1["psr"], f"sky 100: {len(r1)} pulsars with 100 <= K <= 300 that the N = 12 fit gets wrong; median margin")]
for j in range(100, 300):
    r2 = [r for r in P53(j) if 1.0 <= r["pri"] <= 3.0 and not r["w|ref|12|t1"]]
    if r2: p2 = median_pick(r2); sel.append((j, p2["psr"], f"sky {j} (first sky with one): {len(r2)} pulsars with prior width 1-3 pc that the N = 12 fit gets right; median margin")); break
print("selected:", sel, flush=True)
C.J.load(verbose=False); CC = C.NS["C_PC_S"]; f8 = lambda x: [float(v) for v in np.asarray(x).ravel()]
def cx(z): return [[float(np.real(v)), float(np.imag(v))] for v in np.asarray(z).ravel()]


class Export:
    def __init__(self, rank, targets): self.rank = rank; self.targets = targets; self.out = {}
    def sky(self, loc): return {}
    def psr(self, loc):
        i = loc["i"]; j = loc["j"]
        if i not in self.targets: return {}
        Lt = loc["Lt"]; sd = float(loc["pri"][i]); B = H.build(loc, self.rank[j], {"on", "off"}); ST = B["ST"]; D = ST["on"]; order = B["order"]
        B["phi"] = np.random.default_rng([54, 7, loc["seed0"], j, i]).uniform(0, 2*np.pi, B["nC"]); cPc = D["cP"]*np.exp(1j*B["phi"])[:, None]      # as Joint54(randph=True)
        Gd = H.grid(D["kap"], sd); G = H.phasors(B, Gd, Lt); Lg = Gd["Lg"]; M = Gd["M"]; x = Lg-Lt; pr = -Lg**2/(2*sd**2)
        kq = G["on"]["kq"]; eP = G["on"]["eP"]; nb = B["nb"]
        def joint(n, d, Pun):
            S = np.empty(M)
            for c0 in range(0, M, 4096):
                E = np.exp(1j*kq[:n, None]*Lg[None, c0:c0+4096]); T = cPc[:n].T @ E
                S[c0:c0+4096] = ((2.0*np.real(np.conj(T)*d[:, None])-np.abs(T)**2)/Pun[:, None]).sum(0)
            return S
        def single(jj, r, Pun):       # one source's comb: 2 Re(w e^{-i kappa L}) - const, the other known pulsar terms removed from the data
            w = np.sum(r*np.conj(cPc[jj])/Pun); return 2.0*np.real(w*np.exp(-1j*kq[jj]*Lg))-np.sum(D["cP"][jj]**2/Pun)
        cur = {}
        for n in (3, 12):
            d, Pun, PEu, PPu = H.setup(B, G, "ref", n, "t1"); dc = G["on"]["cumeP"][n]
            Sr = joint(n, d, Pun)+pr; Sc = joint(n, dc, Pun)+pr
            Sf = H.stat_fft(B, Gd, G, "ref", n, d, Pun)+pr; assert np.max(np.abs(Sf-Sr)) < 1e-6*max(1.0, np.max(np.abs(Sr)))
            cur[f"joint_N{n}_real"] = Sr; cur[f"joint_N{n}_clean"] = Sc
            if n == 12:
                for jj in range(3):
                    oth = eP[:12].sum(0)-eP[jj]
                    cur[f"single_src{jj+1}_real"] = single(jj, d-oth, Pun); cur[f"single_src{jj+1}_clean"] = single(jj, dc-oth, Pun)
                d12, Pun12, PEu12, PPu12 = d, Pun, PEu, PPu
        cur["prior"] = pr
        win = np.abs(x) <= 15.0; nblk = int(np.ceil(M/4000)); mb = (M//nblk)*nblk if nblk > 1 else M
        red = lambda y: y[:mb].reshape(-1, nblk).max(1) if nblk > 1 else y
        xr = x[:mb].reshape(-1, nblk).mean(1) if nblk > 1 else x
        cp2 = D["cP"][:12]**2; own = cp2.sum(1); w_ = cp2/own[:, None]; bun_on = own/(own+w_ @ Pun12)
        _, Pof, _, _ = H.setup(B, G, "off", 12, "t1"); co2 = ST["off"]["cP"][:12]**2; oo = co2.sum(1); bun_off = oo/(oo+(co2/oo[:, None]) @ Pof)
        mu = loc["mu"][order][:12]
        srcs = [dict(rank=k+1, f_earth_nHz=float(B["fE"][k]*1e9), f_pulsar_nHz=float(D["fP"][k]*1e9), angle_deg=float(np.degrees(mu[k])),
                     tooth_spacing_pc=float(2*np.pi/kq[k]), beta_unmod=float(bun_on[k]), beta_unmod_retardation_off=float(bun_off[k]),
                     power_earth=float(loc["Pc"][order][k]*B["sc"]), power_pulsar=float(loc["Pc"][order][k]*D["pf"][k]*B["sc"])) for k in range(12)]
        ex = dict(sky=int(j), pulsar_index=int(i), true_distance_pc=float(loc["Lpc"]), prior_centre_pc=float(loc["Lpc"]-Lt), prior_width_pc=sd,
                  truth_offset_from_prior_centre_pc=float(Lt), grid_step_pc=float(Gd["dL"]), n_grid=int(M), block_size_whole_range=int(nblk),
                  K_finest_tooth_N12=float(6*sd*kq[:12].max()/(2*np.pi)), K_coarsest_tooth_N12=float(6*sd*kq[:12].min()/(2*np.pi)),
                  winner_offset_pc={k: float(x[int(np.argmax(v))]) for k, v in cur.items() if k.startswith("joint")},
                  sources=srcs, window=dict(offset_pc=f8(x[win]), **{k: f8(v[win]) for k, v in cur.items()}),
                  whole_range=dict(offset_pc=f8(xr), **{k: f8(red(v)) for k, v in cur.items()}),
                  units="offset_pc = trial distance minus true distance [pc]; statistics in natural-log likelihood units (joint_* include the prior; single_* do not)")
        self.out[(j, i)] = dict(combs=ex)
        # ---------------- 1B: one bin, one number (loudest known source) -------------------------------
        bs = int(np.argmax(D["cP"][0]**2)); ref = eP[0, bs]; fb = B["fb"]; Ts = B["Ts"]; non = B["non"]; f2 = loc["f2"]
        items = []
        for k in range(B["nC"]):
            kn = k < 12
            items.append((B["eE"][k, bs], "earth", kn, k));
            if k != 0: items.append((eP[k, bs], "pulsar", kn, k))
        sEc = np.sinc((f2[non]-fb[bs])*Ts); sPc = np.sinc((D["fPn"]-fb[bs])*Ts)
        zE = B["aE"]*np.exp(1j*B["thE"])*sEc; zP = D["aP"]*np.exp(1j*B["thP"])*sPc
        allz = np.concatenate([np.array([t[0] for t in items]), zE, zP])
        tags = [(t[1], t[2]) for t in items]+[("earth", False)]*len(zE)+[("pulsar", False)]*len(zP)
        o = np.argsort(-np.abs(allz)); top = o[:30]; rest = allz[o[30:]].sum()
        raw = ref+allz.sum()+B["d_sea"][bs]+B["d_tim1"][bs]
        fit = d12[bs]; knE = B["cumE"][12][bs]
        resid = abs(raw-knE-fit)/abs(fit)
        self.out[(j, i)]["phasor"] = dict(sky=int(j), pulsar_index=int(i), bin_index=int(bs+1), bin_frequency_nHz=float(fb[bs]*1e9), source_rank=1,
            own_pulsar_term=[1.0, 0.0], own_power=float(abs(ref)**2),
            largest_30=[dict(re=float((allz[k]/ref).real), im=float((allz[k]/ref).imag), term=tags[k][0], known=bool(tags[k][1]),
                             subtracted_in_fit=bool(tags[k][1] and tags[k][0] == "earth")) for k in top],
            remainder_of_discrete_sources=cx(rest/ref)[0], n_in_remainder=int(len(allz)-30), sea_draw=cx(B["d_sea"][bs]/ref)[0], timing_noise_draw=cx(B["d_tim1"][bs]/ref)[0],
            bin_datum_all_terms=cx(raw/ref)[0], known_earth_terms_sum=cx(knE/ref)[0], bin_datum_used_by_fit=cx(fit/ref)[0],
            residual_relative=float(resid), beta_unmod_this_source=float(bun_on[0]),
            units="complex numbers as [real, imaginary], in units of the source's own pulsar-term phasor in this bin; bin_datum_all_terms = own + all listed + remainder + sea + timing; the fit uses that minus known_earth_terms_sum")
        assert resid < 1e-9
        # ---------------- 1C: spectrum ------------------------------------------------------------------
        kP = D["cumPP"][12]; kE = B["cumPE"][12]
        o3 = np.argsort(-loc["P0"])[:3000]
        self.out[(j, i)]["spectrum"] = dict(sky=int(j), pulsar_index=int(i),
            bins=dict(frequency_nHz=f8(fb*1e9), known_pulsar_terms=f8(kP), known_earth_terms=f8(kE), unknown_earth_terms=f8(PEu12), unknown_pulsar_terms=f8(PPu12),
                      sea_earth_plus_pulsar=f8(B["seab"]), timing_noise=f8(B["Nb"]), total_gw=f8(kP+kE+PEu12+PPu12+B["seab"])),
            known_sources=srcs,
            loudest_3000=dict(f_earth_nHz=f8(np.round(f2[o3]*1e9, 5)), f_pulsar_nHz=f8(np.round(loc["fP_all"][o3]*1e9, 5)),
                              power_earth=f8(loc["P0"][o3]*B["sc"]), power_pulsar=f8(loc["P0"][o3]*loc["pf_all"][o3]*B["sc"])),
            units="power in the machine's units (strain^2 x frequency x span, per 45-yr bin; same units for every entry); frequencies in nHz; bins are k/T, T = 45 yr")
        return {}


res = {}
for (j, i, why) in sel:
    src = pickle.load(open(f"job53_data/parts/stat_nat_{j:03d}.pkl", "rb"))[1]
    hk = Export({d["sky"]: (d["cand"], d["arr_opt"]) for d in src}, {i})
    C.collect50([j], hook=hk, **R.RUNS["main"][0]); res[(j, i)] = (hk.out[(j, i)], why)
(o1, w1), (o2, w2) = res[(sel[0][0], sel[0][1])], res[(sel[1][0], sel[1][1])]
meta = dict(made_by="job54_export.py (uncommitted, 2026-11-12); fit job54_common.build/setup; sky type: natural", selection_rule_1=w1, selection_rule_2=w2)
json.dump(dict(meta=meta, pulsar_wrong=o1["combs"], pulsar_right=o2["combs"]), open("talk_data/combs_example.json", "w"))
json.dump(dict(meta=meta, **o1["phasor"]), open("talk_data/phasor_example.json", "w"))
json.dump(dict(meta=meta, **o1["spectrum"]), open("talk_data/spectrum_example.json", "w"))
for nm, o in (("wrong", o1), ("right", o2)):
    c = o["combs"]; print(f"pulsar_{nm}: sky {c['sky']} pulsar {c['pulsar_index']}: prior width {c['prior_width_pc']:.2f} pc, K finest {c['K_finest_tooth_N12']:.0f}, grid {c['n_grid']} points step {c['grid_step_pc']:.4f} pc, winners (offset from truth, pc): {c['winner_offset_pc']}")
p = o1["phasor"]; print(f"phasor: bin {p['bin_index']} at {p['bin_frequency_nHz']:.2f} nHz; residual {p['residual_relative']:.1e}; |datum used by fit| {np.hypot(*p['bin_datum_used_by_fit']):.2f} x own; beta_unmod {p['beta_unmod_this_source']:.3f}")
print("sizes:", {f: os.path.getsize(f"talk_data/{f}") for f in ("combs_example.json", "phasor_example.json", "spectrum_example.json")})
