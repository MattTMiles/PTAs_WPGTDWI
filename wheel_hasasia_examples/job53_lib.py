# JOB 53 analysis helpers
import pickle, numpy as np
TW = 0.25
def ld(n): return pickle.load(open(f"job53_data/run_{n}.pkl", "rb"))
def popid(rec, sky_type): return np.array([(r["sky"] % 50) if sky_type == "constructed" else r["sky"] for r in rec])
def frac(ok, m, pop, nboot=2000, seed=1):
    """fraction of ok among m, 95% interval from resampling populations"""
    ok = np.asarray(ok, float); m = np.asarray(m, bool)
    if m.sum() == 0: return (np.nan, np.nan, np.nan, 0)
    ids = np.unique(pop); num = np.array([(ok*m)[pop == p].sum() for p in ids]); den = np.array([m[pop == p].sum() for p in ids])
    rg = np.random.default_rng(seed); pick = rg.integers(0, len(ids), (nboot, len(ids)))
    bs = num[pick].sum(1)/np.maximum(den[pick].sum(1), 1)
    return (float(num.sum()/den.sum()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), int(m.sum()))
def classes(rec):
    K = np.array([r["Kc"] for r in rec]); return {"all": np.ones(len(rec), bool), "K60": K >= 60, "K600": K >= 600}
def table_row(rec, key, sky_type, ci=True):
    ok = ~np.array([r[key] for r in rec]); pop = popid(rec, sky_type); out = []
    for cn, m in classes(rec).items():
        f = frac(ok, m, pop); out.append(f"{f[0]:.3f} [{f[1]:.3f}, {f[2]:.3f}]" if ci else f"{f[0]:.3f}")
    return out, float(ok.sum()/len(np.unique([r["sky"] for r in rec])))
def prior_only(rec, sky_type):
    ok = np.abs(np.array([r["Ltrue"] for r in rec])) <= TW; pop = popid(rec, sky_type)
    return [f"{frac(ok, m, pop)[0]:.3f}" for m in classes(rec).values()]
def fmt3(x): return " / ".join(x)
