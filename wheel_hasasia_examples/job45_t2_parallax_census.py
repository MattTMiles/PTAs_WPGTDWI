# ================= JOB 45 T2 -- TIMING PARALLAX AS THE THIRD K-CHANNEL ==========
# Yu & Pan 2025 timing-parallax bound, evaluated over the real census
# (best_distances.txt) at the calendar's tiers. EVERY table below carries the
# caveat line -- print it, never drop it:
CAVEAT = ("  CAVEAT: sigma_D here is Yu & Pan's WHITE-NOISE-ONLY timing-parallax "
          "bound\n  (sigma_D = 2.34/cos^2(beta) * (N_obs/100)^(-1/2) * (D/kpc)^2 "
          "* (sigma_n/10ns) pc,\n  fortnightly cadence) -- an OPTIMISTIC floor: "
          "red noise, DM noise and covariances\n  with position/PM/binary terms "
          "are all priced at zero.")
# Channels: timing parallax (all pulsars, NEW); Pbdot via the booked machine law
# sigma_d(Pbdot) = 3.0*(T/30)^(-2.5)*(sigma_n/100) pc (binaries; here only the
# five PBDOT_Shk rows have a measured Pbdot -- binary flags for the rest await
# the psrqpy join, booked limitation). VLBI = external event, excluded.
# Class rule (Job 40 amendment, adopted): A = sigma_d < 5 pc; B = rel <= 10 %
# (not A); C = remainder. Corner entry (booked Job 35 T1b / Job 40): sub-pc,
# sigma_d < 1 pc.
import numpy as np, re, csv

EPS = np.deg2rad(23.4392911)
CAD_PER_YR = 365.25/14.0

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

def jpos_radec(nm):
    m = re.match(r"J(\d{2})(\d{2})([+-])(\d{2})(\d{2})", nm)
    if not m:
        return None
    ra = (int(m.group(1))+int(m.group(2))/60.0)*15.0*np.pi/180
    dec = (int(m.group(4))+int(m.group(5))/60.0)*(1 if m.group(3) == "+" else -1)*np.pi/180
    return ra, dec

def ecl_beta(ra, dec):
    sb = np.sin(dec)*np.cos(EPS) - np.cos(dec)*np.sin(EPS)*np.sin(ra)
    return np.arcsin(np.clip(sb, -1, 1))

def cls_abs(sd_pc, rel):
    if sd_pc < 5.0:
        return "A"
    if rel <= 0.10:
        return "B"
    return "C"

cen = []
for nm, d, e, rel, meth in rows:
    rd = jpos_radec(nm)
    if rd is None:
        continue
    beta = ecl_beta(*rd)
    cen.append(dict(name=nm, d=d, sd=max(e*1000.0, 0.05), rel=rel, meth=meth,
                    cls=cls_abs(e*1000.0, rel), beta=beta))
from collections import Counter
cc = Counter(x["cls"] for x in cen)
print("=" * 74)
print("JOB 45 T2 -- timing parallax over the real census")
print("=" * 74)
print(f"  parsed {len(cen)} pulsars; classes (adopted absolute rule)"
      f" A {cc['A']} / B {cc['B']} / C {cc['C']}  [booked 2/16/68]")

def sig_px(D_kpc, beta, T_yr, sig_ns):
    N = CAD_PER_YR*T_yr
    return 2.34/np.cos(beta)**2 * (N/100.0)**(-0.5) * D_kpc**2 * (sig_ns/10.0)

def sig_pb(T_yr, sig_ns):
    return 3.0*(T_yr/30.0)**(-2.5)*(sig_ns/100.0)

TIERS = (100.0, 30.0, 10.0)
TGRID = np.arange(5.0, 40.0+1e-9, 0.5)

def entry_T(fn, thresh):
    """Earliest span [yr] (<=40) at which fn(T) < thresh; NaN if never."""
    v = np.array([fn(T) for T in TGRID])
    ok = np.where(v < thresh)[0]
    return float(TGRID[ok[0]]) if len(ok) else np.nan

# ---- (2) per-pulsar class-entry dates under PARALLAX ALONE ----------------------
print("\n  (2) CLASS-ENTRY DATES, PARALLAX ALONE (span in yr at which the pulsar")
print("      would enter class A [sigma_D < 5 pc] / sub-pc corner entry"
      " [sigma_D < 1 pc]; '--' = beyond 40 yr)")
print(CAVEAT)
tab2 = []
for x in cen:
    r = dict(name=x["name"], cls=x["cls"], d=x["d"], sd_now=x["sd"],
             beta=np.degrees(x["beta"]))
    for sn in TIERS:
        fn = lambda T, sn=sn: sig_px(x["d"], x["beta"], T, sn)
        r[f"A_{int(sn)}"] = entry_T(fn, 5.0)
        r[f"pc_{int(sn)}"] = entry_T(fn, 1.0)
        r[f"sig40_{int(sn)}"] = fn(40.0)
    tab2.append(r)

def fmt(t):
    return f"{t:5.1f}" if np.isfinite(t) else "   --"

print(f"\n  {'pulsar':>13} {'cls':>3} {'D kpc':>6} {'|b_ecl|':>7} "
      f"{'A@100':>6} {'A@30':>6} {'A@10':>6} {'pc@100':>7} {'pc@30':>6} {'pc@10':>6}")
for r in sorted(tab2, key=lambda r: r["d"]):
    print(f"  {r['name']:>13} {r['cls']:>3} {r['d']:6.2f} {abs(r['beta']):7.1f} "
          f"{fmt(r['A_100']):>6} {fmt(r['A_30']):>6} {fmt(r['A_10']):>6} "
          f"{fmt(r['pc_100']):>7} {fmt(r['pc_30']):>6} {fmt(r['pc_10']):>6}")

for tag, key in (("class-A entry (sigma_D<5pc)", "A"),
                 ("corner entry (sigma_D<1pc)", "pc")):
    print(f"\n  census summary -- {tag}, parallax alone:")
    print(CAVEAT)
    print(f"    {'tier':>7} {'n<=15yr':>8} {'n<=25yr':>8} {'n<=40yr':>8} {'never':>7}")
    for sn in TIERS:
        v = np.array([r[f"{key}_{int(sn)}"] for r in tab2])
        print(f"    {int(sn):5d}ns {np.sum(v <= 15):8d} {np.sum(v <= 25):8d}"
              f" {np.sum(np.isfinite(v)):8d} {np.sum(~np.isfinite(v)):7d}")

# ---- (3) combined calendar: earliest span at which ANY channel delivers ---------
# corner membership's booked requirement: sigma_d < 1 pc (sub-pc rows). Pbdot
# channel only where a Pbdot is measured (the 5 PBDOT_Shk rows).
print("\n  (3) COMBINED CALENDAR -- earliest span where ANY channel reaches the")
print("      sub-pc corner requirement (parallax: all; Pbdot: PBDOT_Shk rows only;")
print("      VLBI external, excluded). Per-pulsar, tiers 100/30/10 ns.")
print(CAVEAT)
tab3 = []
for x in cen:
    haspb = x["meth"] == "PBDOT_Shk"
    r = dict(name=x["name"], cls=x["cls"], d=x["d"], haspb=haspb)
    for sn in TIERS:
        f_px = lambda T, sn=sn: sig_px(x["d"], x["beta"], T, sn)
        f_pb = lambda T, sn=sn: sig_pb(T, sn)
        f_min = (lambda T, a=f_px, b=f_pb: min(a(T), b(T))) if haspb else f_px
        r[f"T_{int(sn)}"] = entry_T(f_min, 1.0)
        r[f"chan_{int(sn)}"] = ("pb" if haspb and np.isfinite(entry_T(f_pb, 1.0))
                                and (not np.isfinite(entry_T(f_px, 1.0))
                                     or entry_T(f_pb, 1.0) <= entry_T(f_px, 1.0))
                                else "px")
    tab3.append(r)
print(f"\n  {'pulsar':>13} {'cls':>3} {'D kpc':>6} {'Pbdot?':>6} "
      f"{'T@100':>6} {'T@30':>6} {'T@10':>6} {'chan@10':>8}")
for r in sorted(tab3, key=lambda r: (r["T_10"] if np.isfinite(r["T_10"]) else 99)):
    if np.isfinite(r["T_10"]) or np.isfinite(r["T_100"]):
        print(f"  {r['name']:>13} {r['cls']:>3} {r['d']:6.2f} {str(r['haspb']):>6} "
              f"{fmt(r['T_100']):>6} {fmt(r['T_30']):>6} {fmt(r['T_10']):>6}"
              f" {r['chan_10']:>8}")
nfin = {int(sn): int(np.sum([np.isfinite(r[f"T_{int(sn)}"]) for r in tab3]))
        for sn in TIERS}
print(f"\n  pulsars reaching sub-pc by 40 yr: {nfin[100]} @100ns,"
      f" {nfin[30]} @30ns, {nfin[10]} @10ns  (of {len(tab3)})")

# ---- crossover distance: parallax vs Pbdot --------------------------------------
# sig_px < sig_pb  <=>  D^2 < 3.0*(T/30)^(-2.5)*(sigma_n/100) /
#                              [2.34/cos^2 beta*(N/100)^(-1/2)*(sigma_n/10)]
# The tier cancels EXCEPT the 10x unit mismatch (pb is quoted per 100 ns, px per
# 10 ns), leaving D_cross INDEPENDENT of sigma_n.
print("\n  (4) CROSSOVER: parallax overtakes Pbdot below")
for T in (20.0, 30.0, 40.0):
    N = CAD_PER_YR*T
    for cb in (1.0, 0.79):     # cos^2 beta = 1 (ecliptic pole-ward best) and census median
        D2 = 3.0*(T/30.0)**(-2.5)/100.0 / (2.34/cb*(N/100.0)**(-0.5)/10.0)
        print(f"    T={T:4.0f} yr, cos^2beta={cb:4.2f}: D_cross ="
              f" {np.sqrt(max(D2,0))*1000:6.0f} pc")
cb_med = float(np.median([np.cos(x['beta'])**2 for x in cen]))
print(f"    census median cos^2(beta) = {cb_med:.3f}")

with open("job45_t2_parallax_census.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["pulsar", "class", "D_kpc", "sigma_d_now_pc", "beta_ecl_deg",
                 "entryA_100ns_yr", "entryA_30ns_yr", "entryA_10ns_yr",
                 "subpc_100ns_yr", "subpc_30ns_yr", "subpc_10ns_yr",
                 "combined_subpc_10ns_yr", "caveat"])
    t3d = {r["name"]: r for r in tab3}
    for r in tab2:
        cw.writerow([r["name"], r["cls"], f"{r['d']:.3f}", f"{r['sd_now']:.2f}",
                     f"{r['beta']:.1f}",
                     *(f"{r[k]:.1f}" if np.isfinite(r[k]) else ""
                       for k in ("A_100", "A_30", "A_10", "pc_100", "pc_30", "pc_10")),
                     (f"{t3d[r['name']]['T_10']:.1f}"
                      if np.isfinite(t3d[r['name']]['T_10']) else ""),
                     "white-noise-only optimistic bound"])
print("  wrote job45_t2_parallax_census.csv")
