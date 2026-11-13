# Loads every file in talk_data/ and prints its shape and one sanity number each.   python talk_data/check.py
import os, csv, json
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); P = lambda f: os.path.join(D, f)
def rd(f):
    r = list(csv.DictReader(open(P(f)))); print(f"{f:<26s} {len(r):6d} rows x {len(r[0]):2d} columns", end="  | "); return r
c = json.load(open(P("combs_example.json")))
for k in ("pulsar_wrong", "pulsar_right"):
    e = c[k]; w = e["whole_range"]; x = np.array(w["offset_pc"]); s = np.array(w["joint_N12_real"]); sc = np.array(w["joint_N12_clean"])
    print(f"combs_example.json[{k}]: {len(x)} points over the prior, {len(e['window']['offset_pc'])} in the window; N = 12 winner at {x[s.argmax()]:+.2f} pc from truth (stored {e['winner_offset_pc']['joint_N12_real']:+.2f}); clean winner at {x[sc.argmax()]:+.3f} pc")
p = json.load(open(P("phasor_example.json"))); z = lambda v: complex(*v)
tot = 1+sum(complex(t["re"], t["im"]) for t in p["largest_30"])+z(p["remainder_of_discrete_sources"])+z(p["sea_draw"])+z(p["timing_noise_draw"])
print(f"phasor_example.json: 30 terms + remainder + sea + timing; |sum - bin datum| = {abs(tot-z(p['bin_datum_all_terms'])):.1e}; |sum - known Earth - fit datum| = {abs(tot-z(p['known_earth_terms_sum'])-z(p['bin_datum_used_by_fit'])):.1e}")
s = json.load(open(P("spectrum_example.json"))); b = s["bins"]
parts = sum(np.array(b[k]) for k in ("known_pulsar_terms", "known_earth_terms", "unknown_earth_terms", "unknown_pulsar_terms", "sea_earth_plus_pulsar"))
print(f"spectrum_example.json: {len(b['frequency_nHz'])} bins, {len(s['known_sources'])} known sources, {len(s['loudest_3000']['f_earth_nHz'])} loudest; max |parts - total|/total = {np.max(np.abs(parts-np.array(b['total_gw']))/np.array(b['total_gw'])):.1e}")
r = rd("curve.csv"); v = [x for x in r if x["arm"] == "reference" and x["timing"] == "x1" and x["N_known"] == "every_candidate"][0]; print(f"reference, every candidate: {v['frac_all']} (Job 53, before the phase correction: 0.260); N = 12: {[x for x in r if x['arm'] == 'reference' and x['timing'] == 'x1' and x['N_known'] == '12'][0]['frac_all']} (Job 53: 0.046)")
r = rd("prior_slice.csv"); print("N = 12 reference at 1 pc:", [x for x in r if float(x["prior_width_pc"]) == 1.0][0]["reference_N12"])
for n in (2, 4, 8, 12):
    r = rd(f"heatmap_Ns{n}.csv"); g = np.array([float(x["fraction_right"]) for x in r]); print(f"min {g.min():.3f} max {g.max():.3f}")
r = rd("heatmap_contours.csv"); print("levels:", sorted({x["level_fraction_right"] for x in r}))
r = rd("jitter_points.csv"); d = max(abs(float(x["cos_dphi_measured"])-float(x["cos_dphi_exact_rice"]))/max(float(x["cos_dphi_se"]), 1e-9) for x in r); print(f"largest |measured - exact <cos>| in standard errors: {d:.1f}")
r = rd("beta_hist.csv"); print("pairs:", sum(int(x["n_pairs_retardation_on"]) for x in r), sum(int(x["n_pairs_retardation_off"]) for x in r))
r = rd("ladder.csv"); v = [x for x in r if x["N_known"] == "12" and x["order"] == "sea_first" and x["skies"] == "all 200"]; print("N = 12 ladder:", " -> ".join(x["frac_all"] for x in v))
r = rd("answers_example.csv"); print(f"right: {sum(int(x['right_within_0.25pc']) for x in r)} of {len(r)}; marked {sum(int(x['marked_for_animation_every_4th_row']) for x in r)}")
r = rd("fisher_curves.csv"); print("modes:", sorted({x["mode"] for x in r}))
for f in ("sensitivity_vs_nph.csv", "sources_strain.csv", "localisation.csv"):
    if os.path.exists(P(f)): r = rd(f); print("columns:", list(r[0].keys())[:6])
if os.path.exists(P("skymaps.npz")):
    z = np.load(P("skymaps.npz"), allow_pickle=True); print("skymaps.npz:", {k: z[k].shape for k in z.files})
tot = sum(os.path.getsize(P(f)) for f in os.listdir(D)); print(f"total size of talk_data/: {tot/1e6:.2f} MB")
