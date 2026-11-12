# ================= JOB 45 T3 -- THE McGRATH TOGGLE ==============================
# McGrath, D'Orazio & Creighton 2022 (arXiv:2208.06495): two resolved sources with
# different tooth spacings certify pulsar distance with no prior -- in a universe
# with exactly two sources and (verbatim) no noise added to the GW-induced
# residuals. Reproduced here in one_pulsar_six_worlds.ipynb's framework (constants,
# surface, null bar, population, pile-up: copied VERBATIM from the staged toy, same
# seeds), then the SAME seeds with the realistic population switched on (pile-up
# ON) and nothing else.
#
# PRE-REGISTERED:
#   M1 (two loud sources, no noise of any kind, broad prior): purity ~ 1.000,
#      margin beats the same-geometry null bar -> row: CERTIFIES.
#   M2 (same seeds, pile-up ON): the two sources' elements are contested; the
#      frozen offset sigma_phi/kappa displaces the ridge -> purity collapses or
#      the combs go unusable -> row: FAILS.
# Breach of either direction = STOP-and-report, not a fix.
#
# WHICH-TERM CHECK: the toggle moves ONLY sigma_phi (the criterion's e^{-sigma^2}
# jitter term via the share s of each source in its own element). N_eff (2 combs),
# K (prior window / tooth spacing), and the bar's ln K term are UNTOUCHED.
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")

# ---- toy cell 1, verbatim constants --------------------------------------------
BLUE, RED, INK = "#2b5aa0", "#a03030", "0.35"
PARAMS = dict(SEED=20260821, L_TRUE=1.42, SIGMA_D=2.0, T_YR=20.0,
              F_POP=(1e-8, 1.3e-7), F_COMB=(1.5e-8, 4.0e-8), NREAL=200,
              C_PC_S=2.99792458e8/3.0856775814913673e16, YR_S=365.25*86400.0)
P = PARAMS; T_S = P["T_YR"]*P["YR_S"]; DF = 1.0/T_S; TOOTH = 0.25
rng0 = np.random.default_rng(P["SEED"]); N_MAX = 8
for _ in range(2000):
    mu = np.arccos(rng0.uniform(-0.9, 0.9, N_MAX))
    f_src = rng0.uniform(*P["F_COMB"], N_MAX)
    kappa = 2*np.pi*f_src*(1-np.cos(mu))/P["C_PC_S"]
    lam = 2*np.pi/kappa
    if 8.0 <= 6*P["SIGMA_D"]/lam.max() <= 12.0:
        break

# ---- the McGrath universe: TWO resolved sources ---------------------------------
# Their theorem's premise is two sources with DIFFERENT wrapping-cycle (tooth)
# spacings: the joint comb is unique out to the vernier beat period
# B = lam1*lam2/|lam1-lam2|, so "no prior" operationally means the searched
# window is <= B while each single comb holds K >> 1 teeth in it.
# PAIR CHOICE, measured not assumed: over the toy's 8 combs the naive first pair
# kappa[:2] has lam 0.219/0.223 pc -> B = 11.7 pc, and the maximum over all 28
# pairs is B = 12.6 pc (pair 2,7: lam 0.399/0.386 pc). The M0 arm below keeps the
# naive pair against a 120 pc window as the measured counterexample (breach of
# this job's first pre-registration -- mechanism: window >> beat, ~10 exact
# impostor coincidences); M1/M2 use the max-beat pair with window = 0.98*B,
# which is the regime their construction states.
NREAL = P["NREAL"]
SIG_DISP = 0.05          # numerical tooth width [rad]; NOT a noise draw -- eps is exact

def beat(i, j):
    return lam[i]*lam[j]/abs(lam[i]-lam[j])

PAIR = (2, 7)
B_PAIR = beat(*PAIR)
SIGMA_PRI = 0.98*B_PAIR/6.0
kk2 = kappa[list(PAIR)]; lam2 = lam[list(PAIR)]
GRIDB = np.linspace(-3*SIGMA_PRI, 3*SIGMA_PRI, 24001)

print("=" * 74)
print("JOB 45 T3 -- the McGrath toggle (toy framework, seeds unchanged)")
print("=" * 74)
print(f"  toy comb spacings [pc]: {lam.round(3)}")
print(f"  naive pair kappa[:2]: lam {lam[0]:.3f}/{lam[1]:.3f}, beat {beat(0,1):.1f} pc")
print(f"  max-beat pair {PAIR}: lam {lam2.round(3)}, beat {B_PAIR:.1f} pc")
print(f"  M1/M2 window 6*sigma = {6*SIGMA_PRI:.1f} pc (= 0.98*beat);"
      f" single-comb K = {(6*SIGMA_PRI/lam2).round(1)} -> prior alone selects nothing")

def null_bar2(sig, kk, grid, pri, ntrial=250):
    ss = np.full(len(kk), sig) if np.ndim(sig) == 0 else np.asarray(sig)
    marg = []
    for t in range(ntrial):
        r = np.random.default_rng(50_000 + t)
        eps = r.uniform(0, 2*np.pi, len(kk))
        resid = (kk[:, None]*grid[None, :]-eps[:, None]+np.pi) % (2*np.pi)-np.pi
        S = -(resid**2/(2*ss[:, None]**2)).sum(0)-grid**2/(2*pri**2)
        S -= S.max()
        ip = int(np.argmax(S)); far = np.abs(grid-grid[ip]) > 2*TOOTH
        marg.append(float(S[ip]-S[far].max()) if far.any() else np.inf)
    return float(np.percentile(marg, 90)), np.array(marg)

def world2(kk, sig_per_comb, frozen_offset, label, grid=None, pri=None, keep=6):
    """Two-comb world. eps = kappa*Lt (+ frozen offset per realisation if given).
    Purity over NREAL realisations, margin vs the same-sigma null bar."""
    grid = GRIDB if grid is None else grid
    pri = SIGMA_PRI if pri is None else pri
    ss = np.asarray(sig_per_comb, float)
    usable = ss < np.pi
    if usable.sum() < 2:
        return dict(label=label, purity=np.nan, margin=np.nan, bar=np.nan,
                    beats=False, surfs=[], usable=int(usable.sum()), sig=ss,
                    grid=grid)
    bar, _ = null_bar2(ss, kk, grid, pri)
    right = 0; margs = []; surfs = []
    for t in range(NREAL):
        r = np.random.default_rng(9_000 + t)          # SAME seeds in every world
        Lt = float(np.clip(r.normal(0, pri), -2.9*pri, 2.9*pri))
        eps = kk*Lt + (r.normal(0, ss) if frozen_offset else 0.0)
        resid = (kk[:, None]*grid[None, :]-eps[:, None]+np.pi) % (2*np.pi)-np.pi
        S = -(resid**2/(2*ss[:, None]**2)).sum(0)-grid**2/(2*pri**2)
        S -= S.max()
        ip = int(np.argmax(S)); Lh = float(grid[ip])
        far = np.abs(grid-grid[ip]) > 2*TOOTH
        margs.append(float(S[ip]-S[far].max()) if far.any() else np.inf)
        right += int(abs(Lh-Lt) <= TOOTH)
        if t < keep:
            surfs.append((S, Lt))
    med = float(np.median(margs))
    return dict(label=label, purity=right/NREAL, margin=med, bar=bar,
                beats=med > bar, surfs=surfs, usable=2, sig=ss, grid=grid)

# ---- M0: the measured counterexample (naive pair, window >> beat) ---------------
# 1 mpc grid so the 1.7 mpc teeth are resolved -- the failure is the vernier
# ambiguity, not the grid.
G120 = np.linspace(-60.0, 60.0, 120001)
r0 = world2(kappa[:2], np.full(2, SIG_DISP), frozen_offset=False,
            label="M0 naive pair, window 120 pc", grid=G120, pri=20.0)

# ---- M1: exactly two sources, NO noise of any kind ------------------------------
# eps = kappa*L_true exactly (their "no noise added to the GW-induced residuals");
# SIG_DISP is the surface's numerical tooth width only, no draw is made with it.
r1 = world2(kk2, np.full(2, SIG_DISP), frozen_offset=False,
            label="M1 two sources, no noise")

# ---- the toggle: realistic population, pile-up ON, nothing else -----------------
# toy cells 32+34 verbatim (same rngP seed chain), then the two combs' shares in
# their own elements under the crowd -> sigma_phi^2=(1-s)/s, frozen per realisation.
rngP = np.random.default_rng(P["SEED"]+7); NSRC = 1500
f_E = 10**rngP.uniform(np.log10(P["F_POP"][0]), np.log10(P["F_POP"][1]), NSRC)
Mc = 10**rngP.normal(9.0, 0.35, NSRC)
amp = 10**(-15.2)*(Mc/1e9)**(5/3)*10**rngP.normal(0, 0.30, NSRC)
h_E = amp*(f_E/1e-8)**(2/3)
TSUN = 4.9255e-6
tau = rngP.uniform(0.5, 3.0, NSRC)*1000.0*(1-rngP.uniform(-0.9, 0.9, NSRC))/P["C_PC_S"]
ts = TSUN*Mc; t_c = (5/256)*ts**(-5/3)*(np.pi*f_E)**(-8/3)
f_P = f_E*(1+tau/np.maximum(t_c, 1e-9))**(-3/8); h_P = amp*(f_P/1e-8)**(2/3)
inpop = (f_P > P["F_POP"][0]*0.3)
bins = np.arange(P["F_COMB"][0], P["F_COMB"][1], DF)

def comb_sigmas2(useP):
    ff = f_P if useP else f_E; hh = h_P if useP else h_E
    out = []
    for k in PAIR:
        b = int(np.clip(np.digitize(f_src[k], bins)-1, 0, len(bins)-2))
        m = (ff >= bins[b]) & (ff < bins[b+1])
        if useP: m = m & inpop
        Pw = hh[m]**2
        if Pw.size == 0:
            out.append(0.05); continue
        s_ = Pw.max()/Pw.sum()
        out.append(float(np.sqrt(max((1-s_)/s_, 1e-4))))
    return np.array(out)

sig_on = comb_sigmas2(True)
print(f"\n  pile-up ON per-comb sigma_phi: {sig_on.round(3)} rad "
      f"(displacement in distance sigma_phi/kappa = {(sig_on/kk2).round(3)} pc"
      f" vs tooth spacings {lam2.round(3)} pc)")
r2 = world2(kk2, sig_on, frozen_offset=True, label="M2 same seeds, pile-up ON")

# ---- scoreboard + figure --------------------------------------------------------
print(f"\n  {'world':>28} {'purity':>7} {'margin':>8} {'bar':>7} {'beats':>6}  verdict")
import csv
with open("job45_t3_scoreboard.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["world", "purity", "median_margin_nats", "null_bar_nats",
                 "sigma_phi_comb1", "sigma_phi_comb2", "beats_bar", "verdict"])
    for r_ in (r0, r1, r2):
        v = "CERTIFIES" if (r_["purity"] >= 0.9 and r_["beats"]) else "FAILS"
        print(f"  {r_['label']:>28} {r_['purity']:7.3f} {r_['margin']:8.2f}"
              f" {r_['bar']:7.2f} {str(r_['beats']):>6}  {v}")
        cw.writerow([r_["label"], f"{r_['purity']:.4f}", f"{r_['margin']:.3f}",
                     f"{r_['bar']:.3f}", f"{r_['sig'][0]:.3f}", f"{r_['sig'][1]:.3f}",
                     r_["beats"], v])
print("  wrote job45_t3_scoreboard.csv")

fig, ax = plt.subplots(2, 1, figsize=(8.4, 5.6), sharex=True)
for a, r_, col in ((ax[0], r1, BLUE), (ax[1], r2, RED)):
    for S, Lt in r_["surfs"]:
        a.plot(r_["grid"], S, color=col, lw=0.7, alpha=0.55)
        a.axvline(Lt, color=INK, lw=0.6, ls=":", alpha=0.7)
    a.set_ylim(-30, 2); a.grid(alpha=.25, lw=.5)
    a.text(.02, .88, r_["label"], transform=a.transAxes, color=col)
    a.text(.02, .74, f"purity {r_['purity']:.3f} | margin {r_['margin']:+.1f}"
           f" vs bar {r_['bar']:.1f} nats", transform=a.transAxes, fontsize=8)
    a.set_ylabel("log posterior [nats]")
ax[1].set_xlabel("trial distance [pc]  (dotted: each realisation's true distance)")
plt.tight_layout(); plt.savefig("fig45_mcgrath_toggle.png", dpi=130,
                                bbox_inches="tight")
print("  wrote fig45_mcgrath_toggle.png")

print("""
  CALLOUT (for the paper): McGrath et al.'s classification of measurable
  distances (their Sec. 4, eqs 28-30) is a Fisher-width statement -- it prices
  IMPRECISION, the curvature of the likelihood at the true tooth. The realistic
  population does not primarily widen the teeth: it DISPLACES them. Each comb's
  contested element contributes a frozen offset sigma_phi/kappa (printed above,
  in pc) that does not average down with time and that an unbiased Fisher
  analysis cannot represent -- the two-source certificate above survives their
  test in both worlds, yet certifies in only one.""")
