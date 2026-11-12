"""JOB 48 -- analysis of job48_compute.py outputs: knees, fits, tail, CSV, figures.
Primary estimator = booked (Figure 1 quantity): ratio of medians, i.e.
median_skies(median_pulsars I_marg) / median_skies(median_pulsars I_cond).
Band = 16/84 % over skies of the per-sky ratio. 'full' = every other parameter
(including the other distances) marginalised, same summary."""
import os, json, csv, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/home/mattm/projects/HSYMT/PAPER_figures/job48_data"
OUTF = "/home/mattm/projects/HSYMT/PAPER_figures"
rng = np.random.default_rng(4848)

def load(name):
    fp = os.path.join(D, f"{name}.npz")
    if not os.path.exists(fp):
        return None
    z = np.load(fp, allow_pickle=True)
    a = {k: z[k] for k in z.files}
    a["name"] = name
    return a

def summarise(a, key="mb10"):
    Ns = np.unique(a["N"]); out = dict(N=Ns, r=[], lo=[], hi=[], marg=[], mlo=[], mhi=[],
                                        cond=[], nsk=[], per=[])
    for N in Ns:
        k = a["N"] == N
        cm = np.median(a["cond"][k], 1); mm = np.median(a[key][k], 1)
        rs = mm / cm
        out["r"].append(np.median(mm) / np.median(cm))
        out["lo"].append(np.percentile(rs, 16)); out["hi"].append(np.percentile(rs, 84))
        out["marg"].append(np.median(mm)); out["mlo"].append(np.percentile(mm, 16)); out["mhi"].append(np.percentile(mm, 84))
        out["cond"].append(np.median(cm)); out["nsk"].append(k.sum()); out["per"].append((mm, cm))
    for k_ in ("r", "lo", "hi", "marg", "mlo", "mhi", "cond", "nsk"):
        out[k_] = np.array(out[k_])
    return out

def cross(x, y, level=0.5):
    """first downward crossing of `level`, interpolated in log x"""
    for i in range(len(x) - 1):
        if y[i] >= level > y[i + 1]:
            t = (y[i] - level) / max(y[i] - y[i + 1], 1e-12)
            return float(np.exp(np.log(x[i]) + t * (np.log(x[i + 1]) - np.log(x[i]))))
    return np.nan

def knee(S, nboot=400):
    k = cross(S["N"], S["r"])
    klo, khi = cross(S["N"], S["lo"]), cross(S["N"], S["hi"])      # sky-scatter interval
    bs = []
    for _ in range(nboot):
        r = []
        for mm, cm in S["per"]:
            ii = rng.integers(0, len(mm), len(mm))
            r.append(np.median(mm[ii]) / np.median(cm[ii]))
        bs.append(cross(S["N"], np.array(r)))
    bs = np.array(bs)
    return k, klo, khi, bs

ARMS = ["A_B1", "A_SF", "A_B2", "A_B4", "A_B8", "A_B1hi", "B_58", "B_232", "C_tot"]
A = {n: load(n) for n in ARMS}
A = {n: a for n, a in A.items() if a is not None}
S = {n: summarise(a) for n, a in A.items()}
SF = {n: summarise(a, "mf") for n, a in A.items()}

print("=" * 78)
print("JOB 48 -- RESULTS")
print("=" * 78)
print("  arms present:", list(A))

# ---- knee table ------------------------------------------------------------
print("\n  KNEE TABLE (0.5 crossing, interpolated in log N)")
print("  interval = N where the 16th / 84th percentile sky curves cross 0.5 (sky scatter);")
print("  +/- = bootstrap 16-84 % of the median-curve knee (uncertainty of the median)")
KT = {}
for est, SS, pred_lab in (("booked", S, "H1: 0.50 N_psr"), ("full", SF, "H2: N_psr(2-1/B)/4")):
    print(f"\n  estimator: {est}   (pre-registered {pred_lab})")
    print(f"  {'arm':8s} {'N_psr':>5s} {'B':>5s} {'N_knee':>8s} {'sky 16-84':>17s} {'boot 16-84':>17s}"
          f" {'n_knee':>7s} {'n/N_psr':>8s} {'pred n':>7s}")
    for nm, s in SS.items():
        a = A[nm]; B = float(a["B"]); npsr = int(a["n_psr"])
        k, klo, khi, bs = knee(s)
        pred = 0.5 * npsr if est == "booked" else npsr * (2 - 1 / B) / 4
        b16, b84 = np.nanpercentile(bs, [16, 84]) if np.isfinite(bs).any() else (np.nan, np.nan)
        KT[(est, nm)] = dict(k=k, klo=klo, khi=khi, bs=bs, B=B, npsr=npsr)
        print(f"  {nm:8s} {npsr:5d} {B:5.2f} {k:8.1f} [{klo:7.1f},{khi:7.1f}] [{b16:7.1f},{b84:7.1f}]"
              f" {k/B:7.1f} {k/B/npsr:8.3f} {pred:7.1f}")

# ---- fit N_knee = a B^beta ------------------------------------------------------
print("\n  FIT N_knee = a * B^beta over B = 1, 2, 4, 8.05 (N_psr = 116); errors: bootstrap")
FIT = {}
for est in ("booked", "full"):
    arms = [n for n in ("A_B1", "A_B2", "A_B4", "A_B8") if (est, n) in KT]
    if len(arms) < 3:
        continue
    Bv = np.array([KT[(est, n)]["B"] for n in arms]); kv = np.array([KT[(est, n)]["k"] for n in arms])
    p = np.polyfit(np.log(Bv), np.log(kv), 1)
    boots = []
    for i in range(len(KT[(est, arms[0])]["bs"])):
        kb = np.array([KT[(est, n)]["bs"][i] for n in arms])
        if np.all(np.isfinite(kb)):
            boots.append(np.polyfit(np.log(Bv), np.log(kb), 1))
    boots = np.array(boots)
    be16, be84 = np.percentile(boots[:, 0], [16, 84]); a16, a84 = np.percentile(np.exp(boots[:, 1]), [16, 84])
    FIT[est] = (np.exp(p[1]), p[0])
    print(f"    {est:6s}: a = {np.exp(p[1]):6.1f} [{a16:.1f}, {a84:.1f}]   beta = {p[0]:.3f} [{be16:.3f}, {be84:.3f}]"
          f"   (pred beta 1, a 52-58)")

# ---- same frequency, zero point -----------------------------------------------
for est, SS in (("booked", S), ("full", SF)):
    if "A_SF" in SS:
        s = SS["A_SF"]; z = s["N"][np.argmax(s["r"] < 1e-3)] if (s["r"] < 1e-3).any() else np.nan
        print(f"\n  same-frequency ({est}): ratio first < 1e-3 at N = {z}"
              f"  (pre-reg: {'116' if est == 'booked' else '58'})")
        print("    N    : " + " ".join(f"{int(n):6d}" for n in s["N"]))
        print("    ratio: " + " ".join(f"{v:6.3f}" for v in s["r"]))

# ---- conditional is the same for all B --------------------------------------------
if all(n in S for n in ("A_B1", "A_B2", "A_B4", "A_B8")):
    common = set(S["A_B1"]["N"])
    for n in ("A_B2", "A_B4", "A_B8"):
        common &= set(S[n]["N"])
    common = sorted(common)
    spread = []
    for N in common:
        cs = [S[n]["cond"][list(S[n]["N"]).index(N)] for n in ("A_B1", "A_B2", "A_B4", "A_B8")]
        spread.append(np.ptp(cs) / np.mean(cs))
    print(f"\n  conditional info across B at common N {common}: max relative spread {max(spread):.3f}")

# ---- shape: deviation from the linear law, tail ------------------------------------
print("\n  SHAPE on the collapsed axis x = n/N_psr")
TAIL = {}
for est, SS in (("booked", S), ("full", SF)):
    print(f"\n  estimator: {est}")
    for nm in ("A_B1", "A_B2", "A_B4", "A_B8", "A_SF", "B_58", "B_232"):
        if nm not in SS:
            continue
        s = SS[nm]; a = A[nm]; B = float(a["B"]); npsr = int(a["n_psr"])
        n = s["N"] / B; x = n / npsr
        lin = 1 - x if est == "booked" else 1 - 2 * n / (npsr * (2 - 1 / B))
        kk = KT[(est, nm)]["k"] / B
        below = n < kk
        dev = (s["r"] - lin)[below]
        md = dev[np.argmax(np.abs(dev))] if below.any() else np.nan
        tail = (n > kk) & (s["r"] > 1e-4) & (s["r"] < 0.5)
        g = np.nan
        if tail.sum() >= 3:
            g = -np.polyfit(np.log(n[tail] / kk), np.log(s["r"][tail]), 1)[0]
        TAIL[(est, nm)] = g
        ls = np.diff(np.log(np.maximum(s["r"], 1e-12))) / np.diff(np.log(n))
        lsl = [f"{v:.1f}" for v, m in zip(ls, (n[1:] > kk) & (s["r"][1:] > 1e-4)) if m]
        print(f"    {nm:7s} largest deviation below knee (measured - linear law) = {md:+.3f};"
              f" tail gamma = {g:.2f} ({tail.sum()} pts); local slopes past knee: {', '.join(lsl[:8])}")
    for nm in ("A_B1", "A_B2", "A_B4", "A_B8"):
        pass
    pool_x, pool_y = [], []
    for nm in ("A_B1", "A_B2", "A_B4", "A_B8"):
        if nm in SS:
            s = SS[nm]; B = float(A[nm]["B"]); kk = KT[(est, nm)]["k"] / B; n = s["N"] / B
            m = (n > kk) & (s["r"] > 1e-4) & (s["r"] < 0.5)
            pool_x += list(np.log(n[m] / kk)); pool_y += list(np.log(s["r"][m]))
    if len(pool_x) >= 3:
        p = np.polyfit(pool_x, pool_y, 1); res = np.array(pool_y) - np.polyval(p, pool_x)
        print(f"    POOLED (B = 1,2,4,8.05): gamma = {-p[0]:.2f}, rms residual in ln(ratio) = {res.std():.2f}")
        TAIL[(est, "pooled")] = -p[0]

# ---- A.2 frequency dependence -----------------------------------------------------
if "A_B1" in S and "A_B1hi" in S:
    for est in ("booked", "full"):
        k1, k2 = KT[(est, "A_B1")]["k"], KT[(est, "A_B1hi")]["k"]
        b1, b2 = KT[(est, "A_B1")]["bs"], KT[(est, "A_B1hi")]["bs"]
        dd = b2 - b1
        print(f"\n  A.2 ({est}): B=1 knee at 7.39 nHz = {k1:.1f}, at 15.84 nHz = {k2:.1f};"
              f" difference {k2-k1:+.1f} (bootstrap 16-84 [{np.nanpercentile(dd,16):+.1f},{np.nanpercentile(dd,84):+.1f}])")

# ---- C.2 normalisation --------------------------------------------------------------
if "C_tot" in S and "A_B1" in S:
    c = S["C_tot"]; b = S["A_B1"]
    diffs = [abs(c["r"][i] - b["r"][list(b["N"]).index(N)]) for i, N in enumerate(c["N"]) if N in b["N"]]
    print(f"\n  C.2: fixed total power vs fixed per-source strain at B=1, {len(diffs)} common N:"
          f" max |ratio difference| = {max(diffs):.1e}")

# ---- C.3 ---------------------------------------------------------------------------
c3 = load("C3_6p")
if c3 is not None:
    s6 = summarise(c3); s6f = summarise(c3, "mf"); s6b = summarise(c3, "mb13")
    k6 = cross(s6["N"], s6["r"])
    print(f"\n  C.3 (N_psr=58, B=1, 6 params/source): booked knee N = {k6:.2f}"
          f" (n/N_psr = {k6/58:.3f}; pred N_psr/6 = {58/6:.2f})")
    print("    N        : " + " ".join(f"{int(n):8d}" for n in s6["N"]))
    print("    booked   : " + " ".join(f"{v:8.4f}" for v in s6["r"]))
    print("    booked13 : " + " ".join(f"{v:8.4f}" for v in s6b["r"]))
    print("    full     : " + " ".join(f"{v:8.1e}" for v in s6f["r"]))

# ---- CSV --------------------------------------------------------------------------
with open(os.path.join(OUTF, "job48_turnover_per_bin.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["arm", "B", "centre_frequency_nHz", "N_psr", "N", "normalisation", "n_real",
                "median_ratio", "p16", "p84", "full_marg_median_ratio", "full_p16", "full_p84",
                "I_cond_median_kpc-2", "I_marg_median_kpc-2"])
    for nm, s in S.items():
        a = A[nm]; sf = SF[nm]
        norm = "fixed_total_power" if "total" in str(a["mode"]) else "fixed_per_source_strain"
        fcs = "same-frequency " + f"{float(a['fc'])*1e9:.3f}" if nm == "A_SF" else f"{float(a['fc'])*1e9:.3f}"
        for i, N in enumerate(s["N"]):
            w.writerow([nm, f"{float(a['B']):.3f}", fcs, int(a["n_psr"]), int(N), norm, int(s["nsk"][i]),
                        f"{s['r'][i]:.5f}", f"{s['lo'][i]:.5f}", f"{s['hi'][i]:.5f}",
                        f"{sf['r'][i]:.5f}", f"{sf['lo'][i]:.5f}", f"{sf['hi'][i]:.5f}",
                        f"{s['cond'][i]:.5e}", f"{s['marg'][i]:.5e}"])
    if c3 is not None:
        for i, N in enumerate(s6["N"]):
            w.writerow(["C3_6param", "1.000", f"{float(c3['fc'])*1e9:.3f}", 58, int(N), "fixed_per_source_strain_6param",
                        int(s6["nsk"][i]), f"{s6['r'][i]:.5f}", f"{s6['lo'][i]:.5f}", f"{s6['hi'][i]:.5f}",
                        f"{s6f['r'][i]:.3e}", f"{s6f['lo'][i]:.3e}", f"{s6f['hi'][i]:.3e}",
                        f"{s6['cond'][i]:.5e}", f"{s6['marg'][i]:.5e}"])
print("\n  wrote job48_turnover_per_bin.csv")

# ---- figures ---------------------------------------------------------------------
cmap = plt.get_cmap("Reds")
BARMS = [("A_B8", 0.35), ("A_B4", 0.55), ("A_B2", 0.75), ("A_B1", 0.95)]
def blab(nm):
    B = float(A[nm]["B"])
    return f"sources spread over B = {round(B, 2):g} frequency bin{'s' if B != 1 else ''}"

def figA(SS, fn, ylab_extra=""):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ref = SS.get("A_B1") or next(iter(SS.values()))
    ax.plot(ref["N"], ref["cond"], color="#2b5aa0", lw=1.6, label="conditional (other parameters fixed; same for all B)")
    for nm, c in BARMS:
        if nm in SS:
            s = SS[nm]; m = s["marg"] > 0
            ax.plot(s["N"][m], s["marg"][m], color=cmap(c), lw=1.6, label=blab(nm))
            ax.fill_between(s["N"], np.maximum(s["mlo"], 1e-30), s["mhi"], color=cmap(c), alpha=0.15, lw=0)
    if "A_SF" in SS:
        s = SS["A_SF"]; m = s["marg"] > 0
        ax.plot(s["N"][m], s["marg"][m], color="k", lw=1.4, label="all sources at one frequency")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("number of CW sources on the sky, N (array fixed at N$_{psr}$ = 116)")
    ax.set_ylabel("distance information per pulsar, I$_{LL}$ = 1/$\\sigma_L^2$ [kpc$^{-2}$]" + ylab_extra)
    lo = min(np.min(s_["marg"][s_["marg"] > 0]) for s_ in SS.values() if (s_["marg"] > 0).any())
    ax.set_ylim(max(lo * 0.5, ref["cond"][0] * 1e-3), ref["cond"].max() * 2)
    ax.legend(fontsize=7.5, loc="lower left")
    plt.tight_layout(); plt.savefig(fn, dpi=180); plt.close()

def figB(SS, fn, pred="H1"):
    two = all(n in SS for n in ("B_58", "A_B1", "B_232"))
    fig, axs = plt.subplots(1, 2 if two else 1, figsize=(12 if two else 6.5, 4.8), squeeze=False)
    ax = axs[0, 0]
    for nm, c in BARMS:
        if nm in SS:
            s = SS[nm]; n = s["N"] / float(A[nm]["B"])
            ax.plot(n, s["r"], color=cmap(c), lw=1.6, label=blab(nm))
            ax.fill_between(n, s["lo"], s["hi"], color=cmap(c), alpha=0.2, lw=0)
    if "A_SF" in SS:
        s = SS["A_SF"]; ax.plot(s["N"], s["r"], color="k", lw=1.3, label="all sources at one frequency")
        ax.fill_between(s["N"], s["lo"], s["hi"], color="k", alpha=0.2, lw=0)
    xx = np.geomspace(0.5, 3000, 400)
    ax.plot(xx, np.clip(1 - xx / 116, 0, None), color="0.5", lw=0.9, ls="--", label="1 − n/N$_{psr}$")
    ax.axhline(0.5, color="0.3", lw=0.8)
    ax.set_xscale("log"); ax.set_ylim(0, 1.05); ax.set_xlim(0.8, 2500)
    ax.set_xlabel("sources per frequency bin, n = N/B")
    ax.set_ylabel("fraction of distance information that survives, I$_{marg}$/I$_{cond}$")
    ax.legend(fontsize=7.5, loc="lower left")
    if two:
        ax = axs[0, 1]
        for nm, col, lab in (("B_58", "#9ecae1", "N$_{psr}$ = 58"), ("A_B1", "#4292c6", "N$_{psr}$ = 116"),
                             ("B_232", "#08306b", "N$_{psr}$ = 232")):
            s = SS[nm]; x = s["N"] / float(A[nm]["n_psr"])
            ax.plot(x, s["r"], color=col, lw=1.6, label=lab + ", B = 1")
            ax.fill_between(x, s["lo"], s["hi"], color=col, alpha=0.2, lw=0)
        xx = np.geomspace(0.005, 30, 400)
        ax.plot(xx, np.clip(1 - xx, 0, None), color="0.5", lw=0.9, ls="--", label="1 − n/N$_{psr}$")
        ax.axhline(0.5, color="0.3", lw=0.8)
        ax.set_xscale("log"); ax.set_ylim(0, 1.05); ax.set_xlim(0.008, 20)
        ax.set_xlabel("sources per frequency bin per pulsar, n/N$_{psr}$")
        ax.set_ylabel("fraction of distance information that survives, I$_{marg}$/I$_{cond}$")
        ax.legend(fontsize=7.5, loc="lower left")
    plt.tight_layout(); plt.savefig(fn, dpi=180); plt.close()

figA(S, os.path.join(OUTF, "job48_figA_turnover_moves.png"))
figB(S, os.path.join(OUTF, "job48_figB_collapse.png"))
figA(SF, os.path.join(OUTF, "job48_figA_turnover_moves_FULLMARG.png"))
figB(SF, os.path.join(OUTF, "job48_figB_collapse_FULLMARG.png"))
print("  wrote job48_figA_turnover_moves.png, job48_figB_collapse.png (+ _FULLMARG supplementary versions)")
