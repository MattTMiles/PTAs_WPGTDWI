"""SPARK readback: re-derive every reported number FROM THE BANK, never from the log.
Pure numpy/scipy -- no jax, no GPU (the point of lean-npz: every floor is re-cuttable)."""
import numpy as np
from scipy.stats import gumbel_r

OUT = "/home/mattm/projects/HSYMT/SPARK_results"
ALPHA, ZF_GATE, C_SD, Z_ALPHA, BOOT, BOOT_SEED = 0.05, 0.20, 2.80, 2.9701952521018403, 4000, 20260716

def adopt(x):
    x = np.asarray(x, float); n = len(x); zf = float(np.mean(x == 0.0))
    if n < 2 or np.ptp(x) <= 0:
        gu = gu_sd = np.nan
    else:
        mu, beta = gumbel_r.fit(x); gu = mu + beta*Z_ALPHA; gu_sd = C_SD*beta/np.sqrt(n)
    emp = float(np.quantile(x, 1-ALPHA))
    rng = np.random.default_rng(BOOT_SEED)
    esd = float(np.std([np.quantile(rng.choice(x, n, replace=True), 1-ALPHA) for _ in range(BOOT)]))
    if zf <= ZF_GATE: return dict(floor=float(gu), err=float(gu_sd), est="gumbel", zf=zf, n=n)
    return dict(floor=emp, err=esd, est="emp_q95", zf=zf, n=n)

g0 = np.load(f"{OUT}/spark_g0.npz", allow_pickle=True)
z  = np.load(f"{OUT}/spark_s2c.npz", allow_pickle=True)

print("=== S1 g0 gate chain (banked) ===")
for k in ["g0a_maxdev", "g0b_maxdev", "g0b_free_maxdev", "g0c_ncap", "passed"]:
    print(f"  {k:20s} {g0[k]}")

states = ["s0", "sC_m1e07", "sC_m1e05", "sC_m2e03", "sMAX"]
print("\n=== S2c ledger, RE-DERIVED from raw banked columns ===")
print(f"{'state':12s} {'Ncert':>6s} {'2F med':>8s} {'2F max':>8s} {'floor':>8s} {'err':>6s} "
      f"{'est':>8s} {'zf':>6s} {'bar':>7s} {'clears':>7s} {'gain':>8s}")
res = {}
s0v = z["twoF_s0"]
for s in states:
    v = z[f"twoF_{s}"]; pm = z[f"pmask_{s}"]
    a = adopt(np.maximum(z[f"null_twoF_{s}"], 0.0))
    dev = abs(a["floor"] - float(z[f"floor_{s}_floor"]))
    bar = a["floor"] + a["err"]; cl = int((v > bar).sum())
    res[s] = dict(cl=cl, med=float(np.median(v)), mx=float(v.max()), bar=bar,
                  fl=a["floor"], dev=dev, n=int((pm > 0).sum()))
    print(f"{s:12s} {int((pm>0).sum()):6d} {np.median(v):8.4f} {v.max():8.4f} {a['floor']:8.4f} "
          f"{a['err']:6.4f} {a['est']:>8s} {a['zf']:6.3f} {bar:7.4f} {cl:7d} "
          f"{float(np.median(v-s0v)):+8.4f}   [floor re-cut |dev|={dev:.2e} "
          f"{'OK' if dev < 1e-9 else '*** MISMATCH'}]")

print("\n=== RECRUITMENT (the cascade readout) ===")
b = res["s0"]["cl"]
for s in states:
    print(f"  {s:12s} clears {res[s]['cl']:2d}/13   vs s0: {res[s]['cl']-b:+d}")

print("\n=== certified sets, re-derived from CHORUS raw ===")
for g in ["m1e07", "m1e05", "m2e03"]:
    n, nt, nw = z[f"ncert_{g}"], z[f"ncert_true_{g}"], z[f"ncert_wrong_{g}"]
    cb = float(z[f"corr_banked_{g}"])
    print(f"  {g}: N_cert mean {n.mean():.3f} (min {n.min()}, max {n.max()}, used={res['sC_'+g]['n']}) "
          f"| correct {nt.mean():.3f} vs banked corr {cb:.3f} | wrong {nw.mean():.3f} "
          f"| wrong-rate {nw.sum()/max(n.sum(),1):.4f}")

print("\n=== COST: lnK growth (exact, CPU) ===")
ns, ks, lm = z["lnK_n_src"], z["lnK_K_sum"], z["lnK_med"]
for i in [0, 2, 3, 15]:
    print(f"  n_src={ns[i]:2d}  K_sum={ks[i]:10.0f}  lnK_med={lm[i]:.4f}")
print(f"  lnK cost 3->16: +{float(z['lnK_cost_3_to_16']):.4f} nat | K_sum ratio {ks[15]/ks[2]:.3f}x")
print(f"\nVERDICT (banked): {z['verdict']}")
