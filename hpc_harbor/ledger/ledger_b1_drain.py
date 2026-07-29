#!/usr/bin/env python3
"""LEDGER B1 -- THE DRAIN, RE-REFERENCED AGAINST THE SCRAMBLED ARM (CPU, banked npz only).

THE DEFECT UNDER AUDIT. Every drain number the campaign has quoted is a BITE

    bite(cell, i) = a_bg(cell, i) - a_eff_drawn(cell)                      [dex]

i.e. the fitted background amplitude minus the DEFINITIONAL projection of the drawn
census onto the 13/3 reference. Its implicit null is ZERO: "if the loop resolved
nothing, the fitted background would sit at the drawn population's own projection".
That null is wrong in two independent ways, and both push the same direction:

  (1) FEED-FREE OFFSET. a_bg is a FITTED amplitude of a band-limited GP against a
      finite noisy realisation; a_eff_drawn is an exact projection. They differ at
      iteration 0 BEFORE anything is fed, by an amount that is a property of the
      venue and the realisation, not of the loop. Quoting bite against zero charges
      that offset to the loop.
  (2) SELECTION. The headline numbers are the LARGEST bites across a fan (S4.15's
      "-0.22 dex" is e07_s7, the best of 16; S4.17's table takes a mean over 4 skies
      per rung but the campaign's spoken headline is the max). The reference
      distribution for a max over N cells is the null's max over N cells, not the
      null's mean.

THE CONTROL THAT ALREADY EXISTS. The scrambled arm (`gl1_null{r}_*`, `gl2_*_null{r}_*`)
is the counterpart null: SAME census seed, SAME noise seed at real=0, SAME venue, the
igniter's sky scrambled in the RECOVERY TEMPLATE only (data at truth). It runs the FULL
loop, and -- since S4.15.1's amended null semantics -- it FEEDS real census members too.
So it produces its own feed-drain. That is the missing reference distribution.

WHAT THIS SCRIPT COMPUTES (pre-registered here, before any number is read):

  A. PAIRED EXCESS (the exact contrast). real=0 nulls carry noise_seed
     9_500_000 + 1000*sky, which is EXACTLY the signal cell's noise seed at that sky.
     Same census, same noise, same venue; only the template's igniter sky differs.
         excess_paired(sky, i) = bite_signal(sky, i) - bite_null0(sky, i)
     Band: the two profiles' own curvature widths added in quadrature.

  B. UNPAIRED EXCESS (the distributional contrast). Over the null arm's full
     (sky in {0,1}) x (real 0..4) set at the matched venue:
         excess(sky, i) = bite_signal(sky, i) - mean_null(i)
     Band: sd_null(i) (the reference SPREAD -- the question is whether a signal-arm
     bite is beyond what the null feed produces, not whether the null mean is known).

  C. MAX-OVER-FAN, the winner's-curse statistic proper. max over signal cells of
     bite(i0) vs the null arm's own max over an equal-sized draw. Reported with the
     null max's mean and sd over all C(n_null, n_sig) subsets when small, else all
     null cells as the reference set with the count declared.

  D. FEED-CONDITIONED. bite split by n_fed, because the arms do not feed the same
     number of members: an excess at unmatched n_fed is a feed-count difference, not
     a drain difference. Declared as a CONFOUND wherever the matched cell count is 0.

NOTHING HERE RETRO-CHANGES A BANKED VERDICT. The banked a_bg / a_eff columns are read,
never rewritten. Output bank: reports/ledger_b1_drain.npz. Report-only.

Usage:  python3 hpc_harbor/ledger/ledger_b1_drain.py [--repo <path>]
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import numpy as np

STEM_RE = re.compile(
    r"^(?P<pre>gl1|gl2)_"
    r"(?:(?P<rung>r13p\d+)(?:_w(?P<w>[0-9p]+))?_)?"
    r"(?P<kind>cell|null\d+)_"
    r"(?P<arm>e07|none)_s(?P<sky>\d+)_T(?P<T>\d+)_(?P<tier>\w+?)_i(?P<iter>\d+)$"
)


def scan(results_dir):
    """Read every per-iteration cell/null bank into flat records. Raw columns only."""
    rows = []
    for fn in sorted(os.listdir(results_dir)):
        if not fn.endswith(".npz") or "__" not in fn:
            continue
        stem = fn.split("__")[0]
        m = STEM_RE.match(stem)
        if m is None:
            continue
        try:
            z = np.load(os.path.join(results_dir, fn), allow_pickle=True)
        except Exception as e:                                  # noqa: BLE001
            print(f"  !! unreadable {fn}: {e}", file=sys.stderr)
            continue
        if "a_bg" not in z.files or "a_eff_drawn" not in z.files:
            continue
        g = m.groupdict()
        kind = g["kind"]
        rows.append(dict(
            stem=stem, stage=g["pre"], rung=g["rung"] or "stage1",
            wscale=float((g["w"] or "1").replace("p", ".")),
            arm=g["arm"], sky=int(g["sky"]), T=int(g["T"]), it=int(g["iter"]),
            scrambled=(kind != "cell"),
            real=(0 if kind == "cell" else int(kind[4:])),
            a_bg=float(z["a_bg"]), a_bg_sig=float(z["a_bg_sig"]),
            a_eff=float(z["a_eff_drawn"]),
            n_fed=int(np.asarray(z["fed_mask"]).sum()),
            n_res=int(z["n_resolved"]), n_cert=int(z["n_cert"]),
        ))
    return rows


def venue(r):
    return (r["stage"], r["rung"], r["wscale"], r["T"], r["arm"])


def fmt(x, n=3):
    return "  nan" if not np.isfinite(x) else f"{x:+.{n}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/taylor_group/matt_miles/PTAs_WPGTDWI")
    ap.add_argument("--iter", type=int, default=0,
                    help="the FIRST-BITE iteration the campaign quotes (default 0)")
    a = ap.parse_args()
    rdir = os.path.join(a.repo, "GLACIER_results")
    rows = scan(rdir)
    print(f"[LEDGER-B1] {len(rows)} per-iteration banks read from {rdir}\n")

    for r in rows:
        r["bite"] = r["a_bg"] - r["a_eff"]

    # ---------------- feed-free offset: the i0 bite AT ZERO FED ----------------
    # (1) above, measured. Any cell whose iteration-0 feed is empty is a direct read
    # of "what does bite say when the loop has done nothing".
    z0 = [r["bite"] for r in rows if r["it"] == 0 and r["n_fed"] == 0]
    print("-- (1) FEED-FREE OFFSET: bite at iteration 0 with ZERO fed members --")
    if z0:
        print(f"   n={len(z0)}  mean {np.mean(z0):+.3f}  sd {np.std(z0, ddof=1):.3f}  "
              f"min {np.min(z0):+.3f}  max {np.max(z0):+.3f} dex")
        print("   -> this is the zero-point the campaign's bite is quoted AGAINST ZERO.")
    else:
        print("   no zero-fed iteration-0 cell in the bank")
    print()

    # ---------------- A. paired excess ----------------
    key = lambda r: (venue(r), r["sky"], r["it"])                      # noqa: E731
    sig = {key(r): r for r in rows if not r["scrambled"]}
    nul0 = {key(r): r for r in rows if r["scrambled"] and r["real"] == 0}
    paired = []
    for k, s in sig.items():
        n = nul0.get(k)
        if n is None:
            continue
        band = float(np.hypot(s["a_bg_sig"], n["a_bg_sig"]))
        paired.append(dict(venue=k[0], sky=k[1], it=k[2],
                           bite_sig=s["bite"], bite_null=n["bite"],
                           excess=s["bite"] - n["bite"], band=band,
                           nfed_sig=s["n_fed"], nfed_null=n["n_fed"]))
    print(f"-- (A) PAIRED EXCESS (identical census AND noise seed; template scramble only) --")
    print(f"   {len(paired)} exactly-paired (cell, null real=0) iterations\n")
    hdr = ("   venue                              sky it  bite_sig bite_null  excess"
           "   band  nfed s/n")
    print(hdr)
    pa = [p for p in paired if p["it"] == a.iter]
    for p in sorted(pa, key=lambda d: (str(d["venue"]), d["sky"])):
        v = "/".join(str(x) for x in p["venue"])
        print(f"   {v:34s} {p['sky']:3d} {p['it']:2d}  {fmt(p['bite_sig'])}  "
              f"{fmt(p['bite_null'])}  {fmt(p['excess'])} {p['band']:6.3f}  "
              f"{p['nfed_sig']:3d}/{p['nfed_null']:d}")
    if pa:
        ex = np.array([p["excess"] for p in pa])
        bd = np.array([p["band"] for p in pa])
        print(f"\n   iter {a.iter} paired excess: mean {ex.mean():+.3f} dex, "
              f"sd {ex.std(ddof=1):.3f}, n={len(ex)}; "
              f"|excess| > 2*band in {int(np.sum(np.abs(ex) > 2*bd))}/{len(ex)} cells")
    print()

    # ---------------- B. unpaired excess vs the null distribution ----------------
    print("-- (B) UNPAIRED EXCESS vs the matched-venue scrambled-arm feed-drain --")
    nulls = defaultdict(list)
    for r in rows:
        if r["scrambled"]:
            nulls[(venue(r), r["it"])].append(r)
    unp = []
    print("   venue                              it  n_null  null_mean  null_sd   "
          "sig_mean   excess    z")
    for (v, it), ns in sorted(nulls.items(), key=lambda kv: (str(kv[0][0]), kv[0][1])):
        if it != a.iter:
            continue
        ss = [r for r in rows if venue(r) == v and not r["scrambled"] and r["it"] == it]
        if not ss:
            continue
        nb = np.array([r["bite"] for r in ns])
        sb = np.array([r["bite"] for r in ss])
        sd = float(np.std(nb, ddof=1)) if len(nb) > 1 else np.nan
        exc = float(sb.mean() - nb.mean())
        z = exc / sd if np.isfinite(sd) and sd > 0 else np.nan
        vs = "/".join(str(x) for x in v)
        print(f"   {vs:34s} {it:2d}  {len(nb):5d}  {fmt(nb.mean())}  {sd:7.3f}  "
              f"{fmt(sb.mean())}  {fmt(exc)}  {fmt(z, 2)}")
        unp.append(dict(venue=v, it=it, n_null=len(nb), null_mean=float(nb.mean()),
                        null_sd=sd, sig_mean=float(sb.mean()), n_sig=len(sb),
                        excess=exc, z=z))
    print()

    # ---------------- C. the winner's-curse statistic ----------------
    print("-- (C) MAX-OVER-FAN (the selected headline) vs the null arm's own max --")
    print("   The campaign quotes the LARGEST first bite over a fan. The reference for")
    print("   a max over n cells is the null's max over n cells, not the null's mean.\n")
    print("   venue                              n_sig  max_sig  n_null  max_null  "
          "E[max_null|n_sig]  sd   excess")
    curse = []
    rng = np.random.default_rng(20260729)
    for (v, it), ns in sorted(nulls.items(), key=lambda kv: (str(kv[0][0]), kv[0][1])):
        if it != a.iter:
            continue
        ss = [r for r in rows if venue(r) == v and not r["scrambled"] and r["it"] == it]
        if len(ss) < 2 or len(ns) < 2:
            continue
        sb = np.array([r["bite"] for r in ss])
        nb = np.array([r["bite"] for r in ns])
        # the drain is a FALL, so "largest bite" is the MOST NEGATIVE
        n_sig = len(sb)
        draws = np.array([nb[rng.choice(len(nb), size=min(n_sig, len(nb)),
                                        replace=False)].min()
                          for _ in range(20000)])
        vs = "/".join(str(x) for x in v)
        print(f"   {vs:34s} {n_sig:5d}  {fmt(sb.min())}  {len(nb):5d}  {fmt(nb.min())}"
              f"      {fmt(draws.mean())}     {draws.std(ddof=1):.3f} "
              f"{fmt(sb.min() - draws.mean())}")
        curse.append(dict(venue=v, it=it, n_sig=n_sig, max_sig=float(sb.min()),
                          n_null=len(nb), max_null=float(nb.min()),
                          exp_max_null=float(draws.mean()),
                          sd_max_null=float(draws.std(ddof=1)),
                          excess=float(sb.min() - draws.mean())))
    print()

    # ---------------- D. feed-conditioned ----------------
    print("-- (D) FEED-CONDITIONED bite (the arms do not feed the same count) --")
    print("   arm-side   n_fed   n_cells   mean bite    sd")
    fc = []
    for scr in (False, True):
        by = defaultdict(list)
        for r in rows:
            if r["scrambled"] == scr and r["it"] == a.iter:
                by[r["n_fed"]].append(r["bite"])
        for nf in sorted(by):
            b = np.array(by[nf])
            sd = float(np.std(b, ddof=1)) if len(b) > 1 else np.nan
            print(f"   {'NULL   ' if scr else 'SIGNAL ':10s} {nf:5d}   {len(b):6d}   "
                  f"{fmt(b.mean())}     {sd:.3f}")
            fc.append(dict(scrambled=scr, n_fed=nf, n=len(b),
                           mean=float(b.mean()), sd=sd))
    print()

    # ---------------- E. per-venue zero-fed baseline + the identity of the maxima --------
    print("-- (E) PER-VENUE ZERO-FED BASELINE (the venue's own zero-point for bite) --")
    print("   venue                              n0  zero_fed_bite  sd   | max-bite cell(s)")
    base = {}
    for v in sorted({venue(r) for r in rows}, key=str):
        b0 = [r["bite"] for r in rows if venue(r) == v and r["it"] == a.iter
              and r["n_fed"] == 0]
        cand = [r for r in rows if venue(r) == v and r["it"] == a.iter]
        if not cand:
            continue
        mx = min(cand, key=lambda r: r["bite"])
        vs = "/".join(str(x) for x in v)
        sd = float(np.std(b0, ddof=1)) if len(b0) > 1 else np.nan
        mu = float(np.mean(b0)) if b0 else np.nan
        base[v] = (mu, sd, len(b0))
        tag = f"{'null' + str(mx['real']) if mx['scrambled'] else 'cell'} s{mx['sky']}"
        print(f"   {vs:34s} {len(b0):3d}  {fmt(mu)}       {sd:6.3f} | "
              f"{fmt(mx['bite'])} {tag} (nfed {mx['n_fed']})")
    print()

    out = os.path.join(a.repo, "reports", "ledger_b1_drain.npz")
    np.savez(
        out,
        note=("LEDGER B1: drain re-referenced against the scrambled (counterpart-null) "
              "arm's own feed-drain. REPORT-ONLY: no banked verdict is changed. bite = "
              "a_bg - a_eff_drawn (dex); excess = bite_signal - bite_null. Paired rows "
              "share census AND noise seed (null real=0)."),
        iter_quoted=a.iter,
        zero_fed_i0=np.array(z0),
        paired=np.array([(str(p["venue"]), p["sky"], p["it"], p["bite_sig"],
                          p["bite_null"], p["excess"], p["band"], p["nfed_sig"],
                          p["nfed_null"]) for p in paired], dtype=object),
        unpaired=np.array([tuple(str(x) if isinstance(x, tuple) else x
                                 for x in d.values()) for d in unp], dtype=object),
        unpaired_keys=np.array(list(unp[0].keys()) if unp else []),
        curse=np.array([tuple(str(x) if isinstance(x, tuple) else x
                              for x in d.values()) for d in curse], dtype=object),
        curse_keys=np.array(list(curse[0].keys()) if curse else []),
        feed_cond=np.array([tuple(d.values()) for d in fc], dtype=object),
        feed_cond_keys=np.array(list(fc[0].keys()) if fc else []),
        all_rows_keys=np.array(list(rows[0].keys()) if rows else []),
        all_rows=np.array([tuple(r.values()) for r in rows], dtype=object),
    )
    print(f"[LEDGER-B1] banked -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
