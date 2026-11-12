# ================= JOB 45 -- T1 (alignment lottery) + T2.4 (calendar revision) ==
# Runs ONLY after job45_preflight.py has PASSED (anchor reproduced). Loads the
# machine namespace snapshot (falls back to a full cell rebuild), then:
#
# T1 PRE-REGISTRATIONS (before running -- standing rule):
#   1. Golden pairs (K_pair = 6 sigma_d/DeltaL(mu) < 2, DeltaL = (c/f_E)/(1-cos mu))
#      concentrate where sigma_d is small: expected order class A ~ few per sky,
#      class B ~ 0.1-0.3, class C ~ 0.01 (a-priori estimate: per detected source
#      P(golden) ~ lambda/(3 sigma_d)/2 under isotropy).
#   2. Golden pairs: f_P/f_E -> 1 (tau -> 0 at alignment) and the target's own-
#      element share s ABOVE the detected-population median.
#   3. Rank variance across skies correlates positively with P(>=1 golden pair);
#      the alignment-attributable fraction is measured by an arm with the
#      geometry draw FROZEN across skies (loudness lottery only) -- no
#      pre-registered magnitude, the mechanism was never identified (booked).
#   4. The FINAL corner's certifications contain a NONZERO count of golden-pair
#      certifications, and the corner value is UNCHANGED (nothing here touches
#      the gates). NOTE, stated up front: the banked machine draws cos mu in
#      (-0.9, 0.9), so near-alignment (cos mu > 0.9) is structurally excluded
#      there -- machine golden pairs can only arise from large-lambda combs on
#      tight priors.
# T2.4 PRE-REGISTRATION: with the parallax channel at Yu & Pan's fortnightly
#   normalisation (px -> px * sqrt(100/(26.09*30)) = 0.3573 * px, i.e. 2.80x
#   tighter than the machine's banked px), T(50%) moves EARLIER for near
#   pulsars, unchanged for distant ones; the crossover distance is reported.
#   WHICH-TERM CHECK: px enters the PRIOR WIDTH pri = min(px, pb) only -- it
#   moves K (the ln K term / tooth count), never sigma_phi and never N_eff.
import os, sys, json, time, re, csv, traceback
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
SCRATCH = "/tmp/claude-1005/-home-mattm-soft/986e71f8-eff1-4a6d-b07f-676b52e4ab76/scratchpad"
os.chdir(BASE)
t_wall = time.time()

# ---- machine namespace: snapshot, else rebuild ----------------------------------
NS = {"__name__": "__main__"}
snap = os.path.join(SCRATCH, "job45_banked_state.pkl")
loaded = False
if os.path.exists(snap):
    try:
        import dill
        with open(snap, "rb") as f:
            NS.update(dill.load(f))
        _sd, _f2, _h2, _P0, _sea = NS["_pop18"](0, 0.35)
        assert len(_f2) > 100
        loaded = True
        print(f"  machine namespace: snapshot loaded ({len(NS)} objects), _pop18 OK")
    except Exception as e:
        print(f"  snapshot load failed ({e!r}) -- rebuilding from cells")
        NS = {"__name__": "__main__"}
if not loaded:
    nb = json.load(open("active_wheel_demo_v4.ipynb"))
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] != "code" or i in (12, 43):
            continue
        try:
            exec(compile("".join(c["source"]), f"<v4-cell-{i}>", "exec"), NS)
        except Exception as e:
            print(f"  [cell {i}] FAIL {e!r}")
    print("  machine namespace rebuilt from cells")

for k in ("_pop18", "rn_draw", "depro_sigma", "alias_score24", "TOOTH_W",
          "C_PC_S", "YR_S", "T_OBS", "bin_edges", "bin_cent", "draw_population"):
    assert k in NS, f"machine name missing: {k}"
_pop18 = NS["_pop18"]; depro_sigma = NS["depro_sigma"]
alias_score24 = NS["alias_score24"]; rn_draw = NS["rn_draw"]
C_PC_S = NS["C_PC_S"]; YR_S = NS["YR_S"]; TOOTH_W = NS["TOOTH_W"]
draw_population = NS["draw_population"]

# job43c fragments, verbatim by line span (inv_rate; book; collect36m defined below
# via span 710-756 only -- the STEPS block is re-run here for T1.4)
lines43 = open("job43c_rebooking.py").read().splitlines()
def run_span(lo, hi):
    exec(compile("\n".join(lines43[lo - 1:hi]), f"<43c:{lo}-{hi}>", "exec"), NS)
run_span(485, 509)      # T45 + inv_rate
run_span(661, 689)      # book
run_span(712, 756)      # collect36m def only

# ---- census ---------------------------------------------------------------------
rows = []
for ln in open("/home/mattm/projects/HSYMT/CW_transition/best_distances.txt"):
    if ln.startswith("#") or not ln.strip():
        continue
    p = ln.split()
    if len(p) < 5 or not p[0].startswith("J"):
        continue
    try:
        d = float(p[1]); e = float(p[2])
    except ValueError:
        continue
    rows.append((p[0], d, e, float(p[3].rstrip("%"))/100.0, p[4]))

def jpos(nm):
    m = re.match(r"J(\d{2})(\d{2})([+-])(\d{2})(\d{2})", nm)
    if not m:
        return None
    ra = (int(m.group(1))+int(m.group(2))/60.0)*15.0*np.pi/180
    dec = (int(m.group(4))+int(m.group(5))/60.0)*(1 if m.group(3) == "+" else -1)*np.pi/180
    return np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])

def cls_abs(sd_pc, rel):
    return "A" if sd_pc < 5.0 else ("B" if rel <= 0.10 else "C")

cen = []
for nm, d, e, rel, meth in rows:
    v = jpos(nm)
    if v is None:
        continue
    cen.append(dict(name=nm, d=d, sd=max(e*1000.0, 0.05), rel=rel, meth=meth,
                    cls=cls_abs(e*1000.0, rel), vec=v))
print(f"  census: {len(cen)} pulsars, classes "
      f"{[sum(1 for x in cen if x['cls'] == c) for c in 'ABC']} (booked [2,16,68])")

# ================= T1.1 + T1.2 -- golden-pair counting ===========================
# DETECTION THRESHOLD (stated): the machine's own gate -- Earth-term array SNR
# arr = rho_E * sqrt(0.2*120) >= 3 with rho_E = h/(2 pi f) * sqrt(n_toa/2)/100ns
# at T = 45 yr, 14 d cadence, AND 2e-9 <= f_E <= 3e-8 (the comb band).
print("\n" + "=" * 74)
print("T1 -- THE ALIGNMENT LOTTERY")
print("=" * 74)
NDRAW = 250
T45y = 45.0
nt = T45y*365.25/14.0
TSUN = 4.9255e-6

per = {x["name"]: dict(n=[], nc=[]) for x in cen}
ndet_all = []; ncomb_all = []
gold_srcs = []          # (sky, source idx, pulsar, cosmu, D_kpc) for T1.2
t0 = time.time()
for j in range(NDRAW):
    sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
    rp = (h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100e-9)
    arr = rp*np.sqrt(0.2*120)
    det = (arr >= 3.0) & (f2 >= 2e-9) & (f2 <= 3e-8)
    di = np.where(det)[0]
    combgrade = rp[di] >= 3.0        # the machine's per-pulsar comb-candidate bar
    ndet_all.append(len(di)); ncomb_all.append(int(combgrade.sum()))
    rg = np.random.default_rng(4500 + j)
    cth = rg.uniform(-1, 1, len(di)); phi = rg.uniform(0, 2*np.pi, len(di))
    sth = np.sqrt(1-cth**2)
    S = np.stack([sth*np.cos(phi), sth*np.sin(phi), cth], 1)
    lampc = C_PC_S/f2[di]                       # lambda = c/f_E in pc
    for x in cen:
        cmu = S @ x["vec"]
        dL = lampc/np.maximum(1-cmu, 1e-12)     # DeltaL = lambda/(1-cos mu)
        K = 6*x["sd"]/dL
        g = K < 2.0
        per[x["name"]]["n"].append(int(g.sum()))
        per[x["name"]]["nc"].append(int((g & combgrade).sum()))
        if g.any() and x["cls"] in "AB":
            for a in np.where(g)[0]:
                gold_srcs.append((j, int(di[a]), x["name"], float(cmu[a]),
                                  float(x["d"])))
print(f"  golden-pair sweep: {NDRAW} skies, {time.time()-t0:.0f} s")
print(f"  DETECTION THRESHOLD (stated): array Earth-term SNR arr = rho_E*sqrt(0.2*120)"
      f" >= 3,\n  2-30 nHz, T = 45 yr / 100 ns / 14 d. N_detected median"
      f" {np.median(ndet_all):.0f} [10/90 {np.quantile(ndet_all,.1):.0f},"
      f" {np.quantile(ndet_all,.9):.0f}];\n  of these, comb-candidate grade"
      f" (per-pulsar rho_E >= 3): median {np.median(ncomb_all):.0f}")

tab = []
for x in cen:
    n = np.array(per[x["name"]]["n"]); ncg = np.array(per[x["name"]]["nc"])
    tab.append((x["name"], x["cls"], x["sd"], float(n.mean()),
                float((n >= 1).mean()), float(ncg.mean()), float((ncg >= 1).mean())))
tab.sort(key=lambda r: -r[3])
print(f"\n  T1.1 -- expected golden pairs per sky, per pulsar (top 15 of {len(tab)}):")
print(f"  {'pulsar':>13} {'cls':>3} {'sigma_d pc':>10} {'E[#gold]':>9} {'P(>=1)':>7}"
      f" {'E[gold,rho>=3]':>14} {'P(>=1)':>7}")
for r in tab[:15]:
    print(f"  {r[0]:>13} {r[1]:>3} {r[2]:10.2f} {r[3]:9.2f} {r[4]:7.3f}"
          f" {r[5]:14.3f} {r[6]:7.3f}")
print(f"\n  by class (mean over pulsars):")
for c in "ABC":
    v = [r for r in tab if r[1] == c]
    print(f"    class {c} (n={len(v)}): E[#golden] {np.mean([r[3] for r in v]):.3f}"
          f" [range {min(r[3] for r in v):.3f}-{max(r[3] for r in v):.3f}],"
          f" P(>=1) {np.mean([r[4] for r in v]):.3f};"
          f" comb-grade E {np.mean([r[5] for r in v]):.3f},"
          f" P(>=1) {np.mean([r[6] for r in v]):.3f}")
with open("job45_t1_golden_pairs.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["pulsar", "class", "sigma_d_pc", "E_n_golden_per_sky",
                 "P_ge1_golden", "E_n_golden_combgrade", "P_ge1_combgrade",
                 "n_skies", "median_N_detected", "median_N_combgrade"])
    for r in tab:
        cw.writerow([r[0], r[1], f"{r[2]:.2f}", f"{r[3]:.4f}", f"{r[4]:.4f}",
                     f"{r[5]:.4f}", f"{r[6]:.4f}", NDRAW,
                     int(np.median(ndet_all)), int(np.median(ncomb_all))])
print("  wrote job45_t1_golden_pairs.csv")

# ---- T1.2 pile-up-escape check --------------------------------------------------
# f_P/f_E via the toy's retarded map on the population's own chirp masses
# (draw_population(seed) re-drawn; MCx is the source-frame chirp mass -- the
# redshift correction is omitted, which UNDERSTATES the displacement, i.e. is
# conservative for the escape claim). s = the source's share of its own
# element's total power (discrete window sum + sea), machine bookkeeping.
print("\n  T1.2 -- pile-up escape for golden pairs:")
bin_edges = NS["bin_edges"]; bin_cent = NS["bin_cent"]; T_OBS = NS["T_OBS"]
popcache = {}
def pop_full(j):
    if j not in popcache:
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        F, H, MC, DX, Sf, _ = draw_population(seed=sd)
        if len(F) != len(f2) or not np.allclose(F, f2):
            popcache[j] = None                 # SEED_DATA entry: Mc not recoverable
        else:
            popcache[j] = (f2, h2, P0, sea, MC)
    return popcache[j]

def share_of(f2, P0, sea, k, T=45.0):
    T_s = T*YR_S
    o = np.argsort(f2); fs = f2[o]; Ps = P0[o]; half = 100.0/T_s
    lo = np.searchsorted(fs, f2[k]-half); hi = np.searchsorted(fs, f2[k]+half)
    w = np.sinc((fs[lo:hi]-f2[k])*T_s)**2
    c = (w*Ps[lo:hi]*(T_s/T_OBS)).sum()
    b = int(np.clip(np.digitize(f2[k], bin_edges)-1, 0, len(bin_cent)-1))
    tot = c + sea[b]*(T_OBS/T_s)
    Pt = P0[k]*(T_s/T_OBS)
    return Pt/max(tot, 1e-300)

fr, sh, shmed = [], [], []
nskip = 0
for (j, k, nm, cmu, Dkpc) in gold_srcs[:4000]:
    pc = pop_full(j)
    if pc is None:
        nskip += 1
        continue
    f2, h2, P0, sea, MC = pc
    tau = Dkpc*1000.0*(1-cmu)/C_PC_S
    tc = (5/256)*(TSUN*MC[k])**(-5/3)*(np.pi*f2[k])**(-8/3)
    fr.append(float((1+tau/max(tc, 1e-9))**(-3/8)))
    sh.append(share_of(f2, P0, sea, k))
if fr:
    # population-median share among DETECTED sources, same bookkeeping (20 skies)
    for j in range(0, 20):
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        rp = (h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100e-9)
        det = (rp*np.sqrt(0.2*120) >= 3.0) & (f2 >= 2e-9) & (f2 <= 3e-8)
        for k in np.where(det)[0]:
            shmed.append(share_of(f2, P0, sea, k))
    print(f"    golden pairs sampled {len(fr)} (skipped {nskip} on SEED_DATA skies)")
    print(f"    f_P/f_E: median {np.median(fr):.4f} [10/90 {np.quantile(fr,.1):.4f},"
          f" {np.quantile(fr,.9):.4f}]  (pre-reg: -> 1)")
    print(f"    own-element share s: golden median {np.median(sh):.3f} vs detected-"
          f"population median {np.median(shmed):.3f}  (pre-reg: golden larger)")
    print(f"    sigma_phi = sqrt((1-s)/s): golden {np.sqrt((1-np.median(sh))/max(np.median(sh),1e-9)):.2f}"
          f" vs population {np.sqrt((1-np.median(shmed))/max(np.median(shmed),1e-9)):.2f} rad")

# ================= T1.3 -- attribution: alignment vs loudness ====================
print("\n  T1.3 -- attribution of the sky-sensitivity (rank variance):")
NDRAW3 = 200
inventory20 = NS["inventory20"]
def score_matrix(freeze_geo):
    """Per (pulsar, sky) geometry-gate score, job40-T6 machinery. freeze_geo=True
    pins the source-direction draw to sky 0's seed -- the alignment lottery is
    frozen; the population (loudness) lottery still runs."""
    M = np.full((len(cen), NDRAW3), np.nan)
    for j in range(NDRAW3):
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        NS["_CAP32"] = 1300
        ks, sig, fk, rp, nats = inventory20(f2, h2, P0, P0, sea, 45.0, 100.0, 14.0, 120)
        if len(ks) < 2:
            continue
        rng = np.random.default_rng(7700 + (0 if freeze_geo else j))
        sth = np.arccos(rng.uniform(-1, 1, len(ks)))
        sph = rng.uniform(0, 2*np.pi, len(ks))
        S = np.stack([np.sin(sth)*np.cos(sph), np.sin(sth)*np.sin(sph),
                      np.cos(sth)], 1)
        sgv0 = np.array([depro_sigma(np.sqrt(max(sig[a]**2-1/max(rp[a], 1e-9)**2, 0.0)),
                                     rp[a], fk[a], 0, "loo") for a in range(len(ks))])
        for xi, x in enumerate(cen):
            cmu = S @ x["vec"]; okc = cmu < 0.9
            if okc.sum() < 2:
                continue
            mu = np.arccos(cmu[okc])
            sgv = sgv0[okc]*np.sqrt((1/3)/np.maximum((1+np.cos(mu))**2/4, 1e-3))
            keep = sgv < np.pi
            if keep.sum() < 2:
                continue
            sgv = sgv[keep]
            kap = 2*np.pi*fk[okc][keep]*(1-np.cos(mu[keep]))/C_PC_S
            sigL = 1.0/np.sqrt(max(float(np.sum(kap**2/sgv**2)), 1e-30))
            K = 6*x["sd"]/max(2*np.pi/max(kap.min(), 1e-30), 1e-9)
            M[xi, j] = (1.0 if sigL < TOOTH_W else 0.0)/max(K, 1.0)
    return M

t0 = time.time()
Mfull = score_matrix(False)
Mfroz = score_matrix(True)
print(f"    two score arms done ({time.time()-t0:.0f} s)")

def rank_var(M):
    R = np.full_like(M, np.nan)
    for j in range(M.shape[1]):
        col = M[:, j]
        ok = np.isfinite(col)
        if ok.sum() < 10:
            continue
        r = np.empty(ok.sum())
        r[np.argsort(-col[ok])] = np.arange(ok.sum())
        R[ok, j] = r
    return np.nanvar(R, axis=1), R

vfull, Rfull = rank_var(Mfull)
vfroz, _ = rank_var(Mfroz)
pge1 = np.array([dict((r[0], r[4]) for r in tab)[x["name"]] for x in cen])
ok = np.isfinite(vfull) & np.isfinite(vfroz)
from scipy.stats import spearmanr
rho_g, p_g = spearmanr(vfull[ok], pge1[ok])
frac = 1.0 - np.nansum(vfroz[ok])/max(np.nansum(vfull[ok]), 1e-12)
medr = np.nanmedian(Rfull, axis=1)
mid = ok & (medr > np.nanquantile(medr, 0.25)) & (medr < np.nanquantile(medr, 0.75))
frac_mid = 1.0 - np.nansum(vfroz[mid])/max(np.nansum(vfull[mid]), 1e-12)
print(f"    Spearman(rank variance, P(>=1 golden)) = {rho_g:+.3f} (p = {p_g:.2e})")
print(f"    variance fraction attributable to the geometry draw (1 - frozen/full):")
print(f"      all pulsars   : {frac:.3f}")
print(f"      middle half   : {frac_mid:.3f}   <- the booked '100% sky-sensitive' tier")
print(f"    residual (loudness-lottery) fraction: all {1-frac:.3f}, middle {1-frac_mid:.3f}")

# ================= T1.4 -- the corner does not move ==============================
print("\n  T1.4 -- golden-pair certifications inside the FINAL corner:")
collect36m = NS["collect36m"]; book = NS["book"]
rec = collect36m(fuzzy=True, full_dmu=True, marg=True)
r = book(rec)
print(f"    corner re-run: {r['y']:.2f} [{r['yb'][0]:.2f},{r['yb'][1]:.2f}] @"
      f" {r['pur']:.4f} (gate {r['gate']:.1f}%)  -- must equal the pre-flight row")
sk = np.array([q["sky"] for q in rec]); wr = np.array([q["wrong"] for q in rec], bool)
al = np.array([alias_score24(q) for q in rec])
cal = sk < 100; tst = sk >= 100
qs = np.linspace(0.5, 60, 120); kker = np.ones(7)/7.0
th = [np.quantile(al[cal], q/100) for q in qs]
pu = []
for t in th:
    s = cal & (al <= t); pu.append(1-wr[s].mean() if s.sum() >= 10 else np.nan)
pu = np.array(pu)
ps = np.convolve(np.nan_to_num(pu, nan=np.nanmean(pu[np.isfinite(pu)])), kker, mode="same")
thr = np.nan
for a in range(len(ps)-1):
    if ps[a] >= 0.90 >= ps[a+1]:
        fct = (ps[a]-0.90)/max(ps[a]-ps[a+1], 1e-9)
        thr = th[a]+fct*(th[a+1]-th[a]); break
sel = tst & (al <= thr)
Kmin = np.array([6*q["pri"]*q["kap"].min()/(2*np.pi) for q in rec])
gold = Kmin < 2.0
print(f"    certified (test half): {int(sel.sum())}; of these GOLDEN (K_pair<2):"
      f" {int((sel & gold).sum())} ({(sel & gold).sum()/max(sel.sum(),1):.1%})")
print(f"    golden certified correct: {int((sel & gold & ~wr).sum())}"
      f" of {int((sel & gold).sum())}")
print(f"    ungated golden rate in the record set: {gold.mean():.3f}")
print(f"    corner value unchanged: nothing above touches the gates -- the same")
print(f"    record set and threshold reproduce the pre-flight row printed above.")

# ================= T2.4 -- calendar revision =====================================
print("\n" + "=" * 74)
print("T2.4 -- TIER CALENDAR, PARALLAX AT YU&PAN FORTNIGHTLY NORMALISATION")
print("=" * 74)
print("  CAVEAT: the revised px is the WHITE-NOISE-ONLY optimistic bound (see T2).")
print("  WHICH-TERM: px enters pri = min(px, pb) -> the K / prior-window term only.")
FORT = float(np.sqrt(100.0/(365.25/14.0*30.0)))    # = 0.3573: machine px -> Yu&Pan
print(f"  px multiplier (fortnightly N_obs vs banked normalisation): {FORT:.4f}")

def tier_run45(sig_ns, T, nsky=18, px_fac=1.0):
    NS["_CAP32"] = 1300
    tot = right = 0
    inv20 = NS["inventory20"]
    for j in range(nsky):
        sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
        rng = np.random.default_rng(9300+sd+977*j); Np = 120
        beta_ = np.arcsin(rng.uniform(-1, 1, Np)); Dk = rng.uniform(0.5, 3.0, Np)
        pb = 3.0*(T/30.0)**(-2.5)*(sig_ns/100.0)
        px = 2.34/np.cos(beta_)**2*Dk**2*10.0*np.sqrt(30.0/T)*(sig_ns/100.0)*px_fac
        pri = np.maximum(np.minimum(px, pb), 0.05)
        Arn, _g = rn_draw(Np); rncl = Arn == 0.0
        ks, sig, fk, rp, nats = inv20(f2, h2, P0, P0, sea, T, sig_ns, 14.0, Np)
        if len(ks) < 2:
            continue
        for i in np.where(rncl)[0]:
            sg = np.array([depro_sigma(np.sqrt(max(sig[a]**2-1/max(rp[a], 1e-9)**2, 0.0)),
                           rp[a], fk[a], i, "loo") for a in range(len(ks))])
            okk = sg < np.pi
            if okk.sum() < 2:
                continue
            sgu = sg[okk]; mu = np.arccos(rng.uniform(-0.9, 0.9, len(sgu)))
            sgu = sgu*np.sqrt((1/3)/np.maximum((1+np.cos(mu))**2/4, 1e-3))
            ok2 = sgu < np.pi
            if ok2.sum() < 2:
                continue
            sgu = sgu[ok2]; kap = 2*np.pi*fk[okk][ok2]*(1-np.cos(mu[ok2]))/C_PC_S
            Lt = float(np.clip(rng.normal(0, pri[i]), -2.9*pri[i], 2.9*pri[i]))
            Lg = np.linspace(-3*pri[i], 3*pri[i], 401)
            eps = kap*Lt+rng.normal(0, sgu)
            perr = (kap[:, None]*Lg[None, :]-eps[:, None]+np.pi) % (2*np.pi)-np.pi
            Ssur = -(perr**2/(2*sgu[:, None]**2)).sum(0)-Lg**2/(2*pri[i]**2)
            if 1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)), 1e-30)) >= TOOTH_W:
                continue
            tot += 1
            right += int(abs(float(Lg[int(np.argmax(Ssur))])-Lt) <= TOOTH_W)
    return right/max(tot, 1), tot

Ts = [12, 15, 18, 21, 25, 30, 35, 45]
print(f"\n  {'tier':>8} {'arm':>9} " + " ".join(f"{'T='+str(t):>7}" for t in Ts)
      + f" {'T(50%)':>8}")
T50 = {}
for sn in (100.0, 30.0, 10.0):
    for lab, fac in (("banked", 1.0), ("fortnight", FORT)):
        ps_ = []
        for T in Ts:
            p, n = tier_run45(sn, T)
            ps_.append(p)
        t50 = np.nan
        for a in range(len(ps_)-1):
            if ps_[a] <= 0.5 <= ps_[a+1]:
                t50 = Ts[a]+(0.5-ps_[a])/max(ps_[a+1]-ps_[a], 1e-9)*(Ts[a+1]-Ts[a])
                break
        if np.isnan(t50) and ps_[0] > 0.5:
            t50 = float(Ts[0])
        T50[(sn, lab)] = t50
        print(f"  {int(sn):6d}ns {lab:>9} "
              + " ".join(f"{p:7.3f}" for p in ps_) + f" {t50:8.1f}")
print(f"\n  T(50%) old -> new (anchor class, parallax channel at fortnightly"
      f" normalisation):")
for sn in (100.0, 30.0, 10.0):
    print(f"    {int(sn):4d} ns: {T50[(sn,'banked')]:.1f} -> {T50[(sn,'fortnight')]:.1f} yr"
          f"  [booked {dict([(100.0,30.1),(30.0,23.2),(10.0,17.8)])[sn]}]")
# model-side crossover: px_new < pb at the model's Dk in [0.5, 3] kpc?
for T in (20.0, 30.0):
    D2 = 3.0*(T/30.0)**(-2.5)/100.0/(2.34*np.sqrt(30.0/T)/10.0*FORT/1.0)
    print(f"    model crossover (cos beta = 1) at T = {T:.0f} yr:"
          f" D < {np.sqrt(max(D2, 0))*1000:.0f} pc -> parallax wins")
with open("job45_t2_tier_revision.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["tier_ns", "arm", "T50_yr", "caveat"])
    for (sn, lab), v in T50.items():
        cw.writerow([int(sn), lab, f"{v:.2f}", "white-noise-only optimistic px"])
print("  wrote job45_t2_tier_revision.csv")
print(f"\n  TOTAL RUNTIME: {(time.time()-t_wall)/60:.1f} min (CPU)")
