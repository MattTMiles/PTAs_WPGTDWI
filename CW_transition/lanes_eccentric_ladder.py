"""LANES -- eccentric harmonic-comb fringe/lane ladder (E-track scoping).

TOY / analytic. CPU only (numpy/scipy). No likelihood machinery, no GPU.
Scopes whether an ECCENTRIC binary's Peters-Mathews harmonic comb populates the
F2 "lane gap": the circular Track-B float ceiling ~0.05 scaled vs the loosest
fringe-registration rung 1.85e-3 scaled (a 27x span with ZERO rungs between).

CONVENTIONS / APPROXIMATIONS (state every one):
 A1. Fringe registration tolerance in the loud-source sky param ("scaled" units)
     scales as 1/f_gw (F2: tol = 1/|dPhi/dtheta|, Phi = 2*pi f L(1-cos mu)/c, so
     dPhi/dtheta ~ f). Harmonic n has f_n = n*f_orb -> tol_n = tol_ref * (f_ref/f_n).
 A2. FAIR-DETECTABILITY ANCHOR: the eccentric source is tuned so its POWER-dominant
     harmonic n_peak sits at the same f_gw as the circular F2 reference. Then
     tol at n_peak = 1.85e-3 (matches F2 loosest rung), and rung_n = 1.85e-3*n_peak/n.
     The fundamental (n=1) is the widest lane = 1.85e-3 * n_peak.
 A3. Peters-Mathews (1963) power per harmonic g(n,e) via Bessel J_n(ne). Power
     fraction p_n = g(n,e)/F(e), F(e)=(1+73/24 e^2+37/96 e^4)/(1-e^2)^3.5.
 A4. RESIDUAL SNR weighting: a PTA measures timing residuals r ~ h/(2 pi f).
     Strain h_n ~ sqrt(g)/n (since GW luminosity ~ f^2 h^2 = L_circ*g). So residual
     amplitude R_n ~ h_n/f_n ~ sqrt(g_n)/n^2, and (white-noise) SNR_n^2 ~ g_n/n^4.
     This is the SNR the fringe-fixing actually earns; it REWEIGHTS toward low n.
     [Noise color (red noise / GWB rising to low f) would push usable band even
      lower -- so A4 is an OPTIMISTIC-for-high-n bound; flagged.]
 A5. Fringe suppression (certification budget) at rung n ~ 0.5 * SNR_n^2 * (1-M),
     M = match between adjacent fringes. Toy: adjacent-fringe match at a harmonic is
     ~0 for a resolved rung (orthogonal comb teeth), so (1-M)~1; we report the raw
     0.5*SNR_n^2 per-lane budget and its distribution across usable rungs.
 A6. Evolution envelope: closed-form Peters (1964) e(t),a(t) inspiral. Over the
     pulsar-term lag tau_p = L(1-cos mu)/c ~ kyr the (f_orb,e) differ Earth<->pulsar;
     match decay of the FUNDAMENTAL over one of its own cycles estimates when the
     fundamental lane's K_eff -> 1 (adjacent fundamental fringes de-cohere).

Outputs: LANES_eccentric_ladder.md (written by hand from these numbers),
 lanes_eccentric_ladder.npz, lanes_ladder.png, lanes_kcontour.png.
"""
import numpy as np
from scipy.special import jv

# ---------------------------------------------------------------- F2 anchors
TOL_LOOSEST = 1.85e-3     # scaled: F2 loosest registration rung (at f_gw, circular)
FLOAT_CEIL  = 0.05        # scaled: blind-search float ceiling (de-risk floor)
GAP_SPAN    = FLOAT_CEIL / TOL_LOOSEST   # ~27x
THRESH_FRAC = 0.05        # a harmonic is "usable" if its weight >= 5% of the peak weight

# ---------------------------------------------------------------- Peters-Mathews
def g_n(n, e):
    """Peters & Mathews 1963 Eq(20): power in n-th harmonic / (circular power).
    sum_n g_n = F(e). n integer array, e scalar."""
    n = np.asarray(n, dtype=float)
    ne = n * e
    Jm2 = jv(n - 2, ne); Jm1 = jv(n - 1, ne); J0 = jv(n, ne)
    Jp1 = jv(n + 1, ne); Jp2 = jv(n + 2, ne)
    t1 = Jm2 - 2*e*Jm1 + (2.0/n)*J0 + 2*e*Jp1 - Jp2
    t2 = Jm2 - 2*J0 + Jp2
    return (n**4 / 32.0) * (t1**2 + (1 - e**2)*t2**2 + (4.0/(3*n**2))*J0**2)

def F_e(e):
    return (1 + 73/24*e**2 + 37/96*e**4) / (1 - e**2)**3.5

def harmonic_spectrum(e, nmax=400):
    n = np.arange(1, nmax + 1)
    g = g_n(n, e)
    g = np.clip(g, 0, None)
    p = g / F_e(e)                       # power fraction (A3); sum ~ 1
    snr2 = g / n**4.0                    # unnormalised residual SNR^2 per harmonic (A4)
    w = snr2 / snr2.sum()               # residual SNR^2 fraction
    return n, g, p, w

# ---------------------------------------------------------------- lane ladder
def lane_ladder(e, weight="power"):
    n, g, p, w = harmonic_spectrum(e)
    wt = p if weight == "power" else w
    n_peak = n[np.argmax(wt)]                        # power(or snr)-dominant harmonic
    usable = wt >= THRESH_FRAC * wt.max()            # >=5% of peak (A2/threshold)
    nu = n[usable]
    n_min, n_max = nu.min(), nu.max()
    # rung tolerance in scaled units, anchored so n_peak -> 1.85e-3 (A2)
    rung = TOL_LOOSEST * (n_peak / nu.astype(float))
    widest = rung.max()                              # at n_min
    tightest = rung.min()                            # at n_max
    reaches_ceiling = widest >= FLOAT_CEIL
    # gap analysis: consecutive usable rungs, largest ratio jump in n (in log)
    if len(nu) > 1:
        rungs_sorted = np.sort(rung)[::-1]
        max_gap_dex = np.max(-np.diff(np.log10(rungs_sorted)))
    else:
        max_gap_dex = np.inf
    return dict(e=e, weight=weight, n_peak=int(n_peak), n_min=int(n_min), n_max=int(n_max),
                span=n_max / n_min, usable_n=nu, rung=rung, widest=widest, tightest=tightest,
                reaches_ceiling=bool(reaches_ceiling), max_gap_dex=float(max_gap_dex),
                fund_usable=bool(usable[0]), p=p, w=w, n=n)

# ---------------------------------------------------------------- min-e bridge scan
def min_e_bridge(weight="power", egrid=None):
    if egrid is None:
        egrid = np.linspace(0.01, 0.97, 400)
    reached = []
    for e in egrid:
        L = lane_ladder(e, weight)
        reached.append(L["widest"] >= FLOAT_CEIL)
    reached = np.array(reached)
    if reached.any():
        return float(egrid[np.argmax(reached)]), egrid, reached
    return None, egrid, reached

# ---------------------------------------------------------------- SNR tax / optimal e
def cert_budget(e, SNR_tot=15.0, weight="snr"):
    """Per-lane certification budget 0.5*SNR_n^2*(1-M) ~ 0.5*SNR_tot^2*w_n (A5, (1-M)~1).
    Report: number of usable lanes, min per-lane budget (the weakest rung on the
    vernier -- the bottleneck), and whether every usable rung clears a threshold."""
    L = lane_ladder(e, weight)
    nu = L["usable_n"]
    n, g, p, w = harmonic_spectrum(e)
    w_usable = w[nu - 1]
    budget = 0.5 * SNR_tot**2 * w_usable             # nat per lane
    return dict(e=e, n_lanes=len(nu), min_budget=float(budget.min()),
                max_budget=float(budget.max()), budget=budget, usable_n=nu,
                widest=L["widest"], reaches=L["reaches_ceiling"], span=L["span"])

# ---------------------------------------------------------------- Peters evolution
# closed-form: de/dt, da/dt (Peters 1964). Work in geometric-ish scaled units; we only
# need the FRACTIONAL change of f_orb over the pulsar lag, so absolute G,c drop out via Mc.
G = 6.674e-11; c = 2.998e8; Msun = 1.989e30; yr = 3.156e7; pc = 3.086e16

def peters_rhs(a, e, Mc_kg, m1, m2):
    # a in m; standard Peters 1964 Eqs (5.6),(5.7)
    beta = 64.0/5.0 * G**3 * m1*m2*(m1+m2) / c**5
    dadt = -beta/(a**3*(1-e**2)**3.5) * (1 + 73/24*e**2 + 37/96*e**4)
    dedt = -19.0/12.0 * beta/(a**4*(1-e**2)**2.5) * e * (1 + 121/304*e**2)
    return dadt, dedt

def evolve_back(f_orb0, e0, Mc_solar, tau_s, q=1.0, nsteps=4000):
    """Evolve the binary BACKWARD by tau_s (pulsar term lags Earth term). q=m2/m1.
    Return (f_orb, e) at the pulsar epoch. Mc = chirp mass; split into m1,m2 via q."""
    m1 = Mc_solar*Msun * (1+q)**0.2 / q**0.6      # from Mc=(m1 m2)^.6/(m1+m2)^.2
    m2 = q*m1
    Mtot = m1+m2
    # a from Kepler: (2 pi f_orb)^2 = G Mtot / a^3
    a = (G*Mtot/(2*np.pi*f_orb0)**2)**(1.0/3.0)
    e = e0
    dt = -tau_s/nsteps                              # backward
    for _ in range(nsteps):
        dadt, dedt = peters_rhs(a, e, None, m1, m2)
        a = a + dadt*dt
        e = np.clip(e + dedt*dt, 0.0, 0.999)
        if a <= 0: break
    f_orb = np.sqrt(G*Mtot/a**3)/(2*np.pi)
    return f_orb, e

def fundamental_match_decay(Mc_solar, f_orb0, e0, tau_s):
    """Phase drift of the FUNDAMENTAL over ONE fundamental cycle between Earth and
    pulsar epochs. K_eff (# resolvable fringes) -> 1 when the accumulated fundamental
    phase difference over the observation exceeds ~pi (adjacent fringes de-cohere).
    TOY match: M = |cos(dphi/2)|-like; report dphi_cycle = 2*pi*(f_p - f_e)/f_e
    (fractional freq shift x 2pi = phase slip per cycle)."""
    f_p, e_p = evolve_back(f_orb0, e0, Mc_solar, tau_s)
    dfrac = (f_p - f_orb0) / f_orb0                 # fractional f_orb shift over the lag
    de = e_p - e0
    dphi_cycle = 2*np.pi*abs(dfrac)                 # rad slip accumulated per fundamental cycle
    return dict(f_p=f_p, e_p=e_p, dfrac=float(dfrac), de=float(de),
                dphi_cycle=float(dphi_cycle))

# ================================================================ MAIN
if __name__ == "__main__":
    ECC = [0.1, 0.3, 0.5, 0.7, 0.9]
    print(f"=== LANES eccentric lane ladder ===  GAP span = {GAP_SPAN:.1f}x "
          f"({TOL_LOOSEST:.2e} -> {FLOAT_CEIL:.2e} scaled)\n")

    print("--- Part 2: LANE SPECTRUM (power-weighted usable band, A2 anchor) ---")
    print(f"{'e':>4} {'n_peak':>6} {'n_min':>5} {'n_max':>5} {'span':>6} {'fund?':>5} "
          f"{'widest':>9} {'tight':>9} {'reach.05?':>9} {'maxgap_dex':>10}")
    rows_pw = {}
    for e in ECC:
        L = lane_ladder(e, "power"); rows_pw[e] = L
        print(f"{e:>4} {L['n_peak']:>6} {L['n_min']:>5} {L['n_max']:>5} {L['span']:>6.1f} "
              f"{str(L['fund_usable']):>5} {L['widest']:>9.2e} {L['tightest']:>9.2e} "
              f"{str(L['reaches_ceiling']):>9} {L['max_gap_dex']:>10.3f}")

    print("\n--- Part 2b: same, RESIDUAL-SNR-weighted usable band (A4) ---")
    print(f"{'e':>4} {'n_peak':>6} {'n_min':>5} {'n_max':>5} {'span':>6} {'fund?':>5} "
          f"{'widest':>9} {'reach.05?':>9}")
    rows_sw = {}
    for e in ECC:
        L = lane_ladder(e, "snr"); rows_sw[e] = L
        print(f"{e:>4} {L['n_peak']:>6} {L['n_min']:>5} {L['n_max']:>5} {L['span']:>6.1f} "
              f"{str(L['fund_usable']):>5} {L['widest']:>9.2e} {str(L['reaches_ceiling']):>9}")

    me_pw, eg, rp = min_e_bridge("power")
    me_sw, _, rs = min_e_bridge("snr")
    print(f"\n  MIN e to bridge (widest usable lane >= 0.05):")
    print(f"    power-weighted band : e_min = {me_pw}")
    print(f"    snr-weighted   band : e_min = {me_sw}")

    print("\n--- Part 3: SNR TAX (per-lane cert budget 0.5*SNR_n^2, SNR_tot=15) ---")
    print(f"{'e':>4} {'n_lanes':>7} {'min_budget':>11} {'max_budget':>11} {'widest':>9} {'reach?':>6}")
    budg = {}
    for e in np.round(np.arange(0.1, 0.96, 0.05), 2):
        cb = cert_budget(float(e), SNR_tot=15.0, weight="snr"); budg[float(e)] = cb
        print(f"{e:>4} {cb['n_lanes']:>7} {cb['min_budget']:>11.3f} {cb['max_budget']:>11.3f} "
              f"{cb['widest']:>9.2e} {str(cb['reaches']):>6}")

    print("\n--- Part 4: EVOLUTION ENVELOPE (pulsar-lag fundamental de-coherence) ---")
    # grid over Mc, f_orb; tau_p ~ kyr
    Mc_grid = [1e8, 1e9, 5e9]          # solar masses (SMBHB chirp mass)
    forb_grid = [1e-9, 3e-9, 1e-8]     # Hz (orbital fundamental; f_gw ~ n_peak*forb)
    tau_grid = [1e3*yr, 3e3*yr, 1e4*yr]
    print(f"{'Mc(Msun)':>9} {'f_orb(Hz)':>10} {'tau(kyr)':>9} {'e0':>4} "
          f"{'dfrac_f':>10} {'de':>9} {'dphi/cyc':>9} {'K_eff~1?':>8}")
    evo = {}
    for Mc in Mc_grid:
        for fo in forb_grid:
            for tau in tau_grid:
                for e0 in [0.3, 0.7]:
                    d = fundamental_match_decay(Mc, fo, e0, tau)
                    kc = d['dphi_cycle'] > np.pi     # fundamental lanes de-cohere
                    key = (Mc, fo, tau, e0)
                    evo[str(key)] = d
                    print(f"{Mc:>9.0e} {fo:>10.0e} {tau/yr/1e3:>9.1f} {e0:>4} "
                          f"{d['dfrac']:>10.2e} {d['de']:>9.2e} {d['dphi_cycle']:>9.2e} "
                          f"{str(kc):>8}")

    # ---------------- save
    np.savez("lanes_eccentric_ladder.npz",
             ECC=np.array(ECC), TOL_LOOSEST=TOL_LOOSEST, FLOAT_CEIL=FLOAT_CEIL,
             GAP_SPAN=GAP_SPAN, THRESH_FRAC=THRESH_FRAC,
             min_e_power=(me_pw if me_pw else np.nan),
             min_e_snr=(me_sw if me_sw else np.nan),
             egrid=eg, reached_power=rp, reached_snr=rs,
             **{f"pw_{e}_n": rows_pw[e]["n"] for e in ECC},
             **{f"pw_{e}_p": rows_pw[e]["p"] for e in ECC},
             **{f"pw_{e}_w": rows_pw[e]["w"] for e in ECC},
             **{f"pw_{e}_usable": rows_pw[e]["usable_n"] for e in ECC},
             **{f"pw_{e}_rung": rows_pw[e]["rung"] for e in ECC})
    print("\nsaved lanes_eccentric_ladder.npz")
