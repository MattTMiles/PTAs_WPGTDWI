"""LOTTERY tier-1 (CPU) -- P(switch-on) over candidate eccentricity distributions.

The mixed-population census is CHORUS's: 3 LOUD members + 13 faint (3+13), all at
census loudness. Each of the 3 loud members draws an eccentricity from the candidate
population distribution; the 13 faint members are circular. Per draw we score TWO
switch-on statistics and compare them:

  * CHANNEL-BUDGET (the mechanism, S7.6.4a-c / MAGPIE J1):
        n_active = 13 (faint, 1 channel each) + sum_{loud i} c(e_i)
        switch-on  iff  n_active >= 30
    c(e) = # active harmonics of one member = #{ n in 1..32 : g(n,e) >= 1e-3 max_n g }
    computed from the REAL ATLAS harmonic algebra (atlas_harmonics.g_n), which reproduces
    the ch_README census mapping (e=0.3->8, 0.5->17, 0.7->32, 0.8->32; circular->1) and
    every banked n_active anchor (m0=16, m1e03=23, m2e03=30, m3e03=37, m1e05=32, m1e07=47)
    EXACTLY. N_HARM=32 supplies the e=0.8 magnitude truncation intrinsically (c saturates
    at 32 by e~0.65).

  * THRESHOLD FORM (RECUT_floors.md 2.1, the externally-quoted rule):
        switch-on  iff  (any loud member e_i >= 0.5)  OR  (>= 2 loud members with e_i >= 0.3)

Where the two DISAGREE on a draw is where the proxy (channel budget) and the quoted rule
part company -- the band tier-2 spends GPU on.

Distributions (parameterised so Taylor/Farr posteriors plug straight in):
  (a) circ   : GW-circularised, delta at e=0
  (b) lnN    : log-normal in e, mode e_peak in {0.1,0.3,0.5}, log-space width sigma in {0.1,0.2}
               (X~LogNormal(mu,sigma), mu = ln(e_peak) + sigma^2 so the MODE is e_peak;
                mu,sigma are the pluggable shape params) clipped to [0, EMAX]
  (c) unif   : uniform on [0, 0.9]
  (d) mix    : two-population -- each loud member is eccentric (e=e_char) w.p. f_ecc, else
               circular; sweep f_ecc x e_char  (the observer-facing panel)

CPU only, no GPU, no discovery import. Pure numpy/scipy via atlas_harmonics.
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "atlas"))
import atlas_harmonics as H

# ---- fixed conventions (verbatim from chorus.py) ----------------------------
N_HARM      = 32          # ATLAS production stack depth (validity exact through e=0.7;
                          # e=0.8 magnitude-truncated -- this cap IS the truncation)
G_ACTIVE_REL = 1e-3       # (C3) a harmonic is ACTIVE iff g(n,e) >= 1e-3 * max_n g(n,e)
N_FAINT     = 13          # 3+13 census: 13 faint members, circular, 1 channel each
N_LOUD      = 3
N_BASE      = N_FAINT     # faint contribution to n_active (13 * 1 channel)
CHAN_THRESH = 30          # channel-budget switch: n_active >= 30
E_SINGLE    = 0.5         # threshold-form: any single member >= 0.5
E_PAIR      = 0.3         # threshold-form: any PAIR (>=2 members) >= 0.3
EMAX        = 0.9         # physical cap on eccentricity we consider
SEED        = 20260717
N_DRAW_NAMED = 200_000    # draws per named distribution (a/b/c)  (>> 1e4)
N_DRAW_MIX   = 50_000     # draws per (f_ecc, e_char) mix cell

_NS = np.arange(1, N_HARM + 1)


# ---- c(e): active-channel count for ONE member ------------------------------
def _c_scalar(e):
    if e <= 0.0:
        return 1                                   # circular: only n=2 active
    g = H.g_n(_NS, float(e))
    return int(np.sum(g >= G_ACTIVE_REL * g.max()))

# grid lookup (c(e) is an integer step function; 0.0005 spacing resolves every jump)
_EGRID = np.arange(0.0, EMAX + 1e-9, 0.0005)
_CGRID = np.array([_c_scalar(e) for e in _EGRID], dtype=np.int32)

def c_of_e(e):
    """Vectorised active-channel count per member. e any shape; clipped to [0,EMAX]."""
    e = np.clip(np.asarray(e, float), 0.0, EMAX)
    idx = np.rint(e / 0.0005).astype(int)
    return _CGRID[idx]


# ---- the two switch-on statistics on a (Ndraw x 3) eccentricity matrix ------
def score(e_loud):
    """e_loud: (N,3) eccentricities of the 3 loud members. Returns dict of per-draw
    booleans + aggregate P for both statistics."""
    e_loud = np.asarray(e_loud, float)
    # channel budget
    n_active = N_BASE + c_of_e(e_loud).sum(axis=1)
    chan_on = n_active >= CHAN_THRESH
    # threshold form
    single_on = (e_loud >= E_SINGLE).any(axis=1)
    pair_on   = (e_loud >= E_PAIR).sum(axis=1) >= 2
    thr_on    = single_on | pair_on
    return dict(n_active=n_active, chan_on=chan_on, thr_on=thr_on,
                single_on=single_on, pair_on=pair_on)


def _binom_err(k, n):
    p = k / n
    return np.sqrt(p * (1.0 - p) / n)

def summarise(sc):
    n = len(sc["chan_on"]); ch = int(sc["chan_on"].sum()); th = int(sc["thr_on"].sum())
    disagree = sc["chan_on"] != sc["thr_on"]
    ch_not_th = int((sc["chan_on"] & ~sc["thr_on"]).sum())
    th_not_ch = int((~sc["chan_on"] & sc["thr_on"]).sum())
    return dict(
        n=n,
        P_chan=ch / n, P_chan_err=_binom_err(ch, n),
        P_thr=th / n,  P_thr_err=_binom_err(th, n),
        P_disagree=int(disagree.sum()) / n,
        ch_not_th=ch_not_th / n, th_not_ch=th_not_ch / n,
        n_active_mean=float(sc["n_active"].mean()),
    )


# ---- distribution samplers (return (N,3) loud eccentricities) ---------------
def draw_circ(rng, n):
    return np.zeros((n, N_LOUD))

def draw_lognormal(rng, n, e_peak, sigma):
    """MODE = e_peak, log-space width sigma. mu = ln(e_peak) + sigma^2."""
    mu = np.log(e_peak) + sigma**2
    e = rng.lognormal(mean=mu, sigma=sigma, size=(n, N_LOUD))
    return np.clip(e, 0.0, EMAX)

def draw_uniform(rng, n, lo=0.0, hi=0.9):
    return rng.uniform(lo, hi, size=(n, N_LOUD))

def draw_mix(rng, n, f_ecc, e_char):
    """each loud member eccentric (e=e_char) w.p. f_ecc, else circular."""
    is_ecc = rng.random((n, N_LOUD)) < f_ecc
    return np.where(is_ecc, e_char, 0.0)


# ============================================================================
def main():
    rng = np.random.default_rng(SEED)

    print("=" * 78)
    print("LOTTERY tier-1  --  P(switch-on) over candidate eccentricity distributions")
    print("=" * 78)
    print(f"census 3+13 | channel switch n_active>={CHAN_THRESH} | "
          f"threshold single>={E_SINGLE} or pair>={E_PAIR}")
    print(f"c(e) sanity: e=0->{_c_scalar(0.0)} 0.3->{_c_scalar(0.3)} 0.5->{_c_scalar(0.5)} "
          f"0.7->{_c_scalar(0.7)} 0.8->{_c_scalar(0.8)}  (want 1/8/17/32/32)")
    print(f"n_active sanity m1e03={N_BASE+_c_scalar(0.3)+2} (23) "
          f"m2e03={N_BASE+2*_c_scalar(0.3)+1} (30) m3e03={N_BASE+3*_c_scalar(0.3)} (37) "
          f"m1e05={N_BASE+_c_scalar(0.5)+2} (32) m1e07={N_BASE+_c_scalar(0.7)+2} (47)")
    print()

    named = []   # list of dict rows

    # (a) circularised
    sc = score(draw_circ(rng, N_DRAW_NAMED)); s = summarise(sc)
    named.append(dict(name="circ (delta e=0)", family="a", **s)); del sc

    # (b) log-normal families
    for e_peak in (0.1, 0.3, 0.5):
        for sigma in (0.1, 0.2):
            sc = score(draw_lognormal(rng, N_DRAW_NAMED, e_peak, sigma)); s = summarise(sc)
            named.append(dict(name=f"lnN peak={e_peak} w={sigma}", family="b",
                              e_peak=e_peak, sigma=sigma, **s)); del sc

    # (c) uniform
    sc = score(draw_uniform(rng, N_DRAW_NAMED)); s = summarise(sc)
    named.append(dict(name="unif [0,0.9]", family="c", **s)); del sc

    # ---- print named table -------------------------------------------------
    print("NAMED DISTRIBUTIONS  (P +- binomial sd, N=%d)" % N_DRAW_NAMED)
    print(f"{'distribution':<22} {'P_chan':>14} {'P_thr':>14} {'|dis|':>7} "
          f"{'chNOTth':>8} {'thNOTch':>8} {'<nact>':>7}")
    for r in named:
        print(f"{r['name']:<22} {r['P_chan']:.4f}+-{r['P_chan_err']:.4f} "
              f"{r['P_thr']:.4f}+-{r['P_thr_err']:.4f} {r['P_disagree']:7.4f} "
              f"{r['ch_not_th']:8.4f} {r['th_not_ch']:8.4f} {r['n_active_mean']:7.2f}")
    print()

    # ---- (d) mix grid ------------------------------------------------------
    f_grid = np.round(np.arange(0.0, 1.0001, 0.05), 3)     # 21
    e_grid = np.round(np.arange(0.1, 0.9001, 0.05), 3)     # 17
    P_chan = np.zeros((len(e_grid), len(f_grid)))
    P_thr  = np.zeros_like(P_chan)
    P_dis  = np.zeros_like(P_chan)
    E_chan = np.zeros_like(P_chan)
    E_thr  = np.zeros_like(P_chan)
    for i, ec in enumerate(e_grid):
        for j, fe in enumerate(f_grid):
            sc = score(draw_mix(rng, N_DRAW_MIX, fe, ec)); s = summarise(sc)
            P_chan[i, j] = s["P_chan"]; P_thr[i, j] = s["P_thr"]
            P_dis[i, j]  = s["P_chan"] - s["P_thr"]
            E_chan[i, j] = s["P_chan_err"]; E_thr[i, j] = s["P_thr_err"]
            del sc

    print("MIX GRID computed: f_ecc x e_char = %dx%d cells, N=%d each"
          % (len(f_grid), len(e_grid), N_DRAW_MIX))
    print(f"  max |P_chan - P_thr| in mix = {np.abs(P_dis).max():.3f}")

    # ---- save bank ---------------------------------------------------------
    out = os.path.join(HERE, "..", "..", "reports", "lottery_tier1.npz")
    out = os.path.abspath(out)
    np.savez(out,
             # named
             named_json=json.dumps(named),
             # channel mapping provenance
             egrid_c=_EGRID, cgrid_c=_CGRID,
             # mix surface
             f_grid=f_grid, e_grid=e_grid,
             P_chan=P_chan, P_thr=P_thr, P_dis=P_dis, E_chan=E_chan, E_thr=E_thr,
             # conventions
             conventions=json.dumps(dict(
                 N_HARM=N_HARM, G_ACTIVE_REL=G_ACTIVE_REL, N_FAINT=N_FAINT, N_LOUD=N_LOUD,
                 CHAN_THRESH=CHAN_THRESH, E_SINGLE=E_SINGLE, E_PAIR=E_PAIR, EMAX=EMAX,
                 SEED=SEED, N_DRAW_NAMED=N_DRAW_NAMED, N_DRAW_MIX=N_DRAW_MIX)))
    print(f"\nbanked -> {out}")

    return named, dict(f_grid=f_grid, e_grid=e_grid, P_chan=P_chan, P_thr=P_thr,
                       P_dis=P_dis, E_chan=E_chan, E_thr=E_thr)


if __name__ == "__main__":
    main()
