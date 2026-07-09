"""PENCIL — crossing-time of the registration wall (pure arithmetic, no sim/likelihood).

Q: the F2 lever-arm ladder's LOOSEST registration baseline tolerates 1.85e-3 scaled
sky error; blind Earth-term source localisation floors at 0.05 scaled (6.4 deg) on the
15-yr dataset. The wall closes from ABOVE as data accumulate (sigma_sky shrinks with SNR).
When does the float cross onto the first rung?

CPU only. No GPU, no discovery calls. All arithmetic below.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- ANCHOR (measured)
T0      = 15.0        # yr, current dataset span
SIG0    = 0.05        # scaled, blind F-stat sky float floor at T0 (== 6.4 deg)
SIG0_DEG= 6.4         # deg, same anchor (for readout only)
SNR0    = 10.7        # per-source optimal matched-filter SNR at T0, h=-13.75
LOGH0   = -13.75      # log10 strain of the anchor source
DEG_PER_SCALED = SIG0_DEG / SIG0   # 128 deg per scaled unit (sky param prior width)

# ---------------------------------------------------------------- RUNGS (targets)
RUNG_LOOSE = 1.85e-3  # F2 loosest registration baseline (top of the ladder)
RUNG_L2C   = 1e-4     # L2c conditional-re-solve pull-in radius (bottom, banked)

# ---------------------------------------------------------------- ASSUMPTIONS (named)
# A1  Persistent monochromatic CW, always-on. Coherent matched filter over the full span.
# A2  White-noise-dominated band, FIXED cadence, NO new pulsars, NO noise re-weighting.
#     => per-sinusoid  SNR^2 = h^2 T / (2 S_n)  is LINEAR in T   ->  SNR ~ T^(1/2).
# A3  Blind sky localisation precision is SNR-limited (geometric array response):
#     sigma_sky ~ 1/SNR x (fixed baseline geometry).  Geometry frozen under A2
#     (same pulsars, same sky) -> only the SNR factor moves with T.
#          sigma_sky(T) = SIG0 * (T/T0)^(-1/2).
# A4  Frequency resolution improves SEPARATELY as sigma_f ~ 1/(SNR*T) ~ T^(-3/2)
#     (faster than sky). IRRELEVANT to the wall: F2 showed SKY BINDS -- a 0.05 scaled
#     sky error alone wraps the loosest pulsar-term phase by 0.05/1.85e-3 ~ 27 rad
#     regardless of frequency precision. So the binding axis scales as T^(-1/2).

def sigma_sky(T, p=0.5):
    """Blind sky precision (scaled) at span T. p=0.5 baseline (SNR^2~T); p=1.5 stress test."""
    return SIG0 * (T / T0) ** (-p)

def T_cross(sigma_target, p=0.5):
    """Span T (yr) at which sigma_sky(T) == sigma_target.  T = T0 * (SIG0/target)^(1/p)."""
    return T0 * (SIG0 / sigma_target) ** (1.0 / p)

# =============================================================== (2) crossing times
print("="*70)
print("STEP 2 — crossing times, baseline scaling SNR^2 ~ T  (sigma_sky ~ T^-1/2)")
print("="*70)
for name, rung in [("loosest rung 1.85e-3", RUNG_LOOSE), ("L2c pull-in 1e-4", RUNG_L2C)]:
    ratio = SIG0 / rung
    T = T_cross(rung, p=0.5)
    print(f"\n {name}:")
    print(f"   need sigma_sky shrink factor  SIG0/target = {SIG0:.3g}/{rung:.3g} = {ratio:.4g}")
    print(f"   SNR^2 ~ T  =>  T/T0 = (SIG0/target)^2 = {ratio**2:.6g}")
    print(f"   T = {T0} * {ratio**2:.6g}  =  {T:,.0f} yr   (= {T/1e3:.4g} kyr)")

# =============================================================== (3) vs strain
print("\n" + "="*70)
print("STEP 3 — vs STRAIN instead of T (fix T=T0=15yr)")
print("="*70)
# sigma_sky ~ 1/SNR ; SNR ~ h (linear in strain at fixed noise & span).
# Reach a rung by raising SNR (hence h) instead of waiting.
for name, rung in [("loosest rung 1.85e-3", RUNG_LOOSE), ("L2c pull-in 1e-4", RUNG_L2C)]:
    fac  = SIG0 / rung                 # required SNR boost factor
    snr  = SNR0 * fac                  # threshold per-source SNR
    dlogh= np.log10(fac)               # h ~ SNR linearly -> log10_h shifts by log10(fac)
    logh = LOGH0 + dlogh
    print(f"\n {name}:")
    print(f"   SNR threshold = SNR0 * (SIG0/target) = {SNR0} * {fac:.4g} = {snr:,.1f}")
    print(f"   h ~ SNR (fixed T) -> log10_h = {LOGH0} + log10({fac:.4g}) = {LOGH0}+{dlogh:.3f} = {logh:.3f}")
    print(f"   => h = 10^{logh:.3f} = {10**logh:.3g}   (anchor h = 10^{LOGH0} = {10**LOGH0:.3g})")

# =============================================================== (4) sensitivity SNR^2~T^3
print("\n" + "="*70)
print("STEP 4 — sensitivity: SNR^2 ~ T^3  (sigma_sky ~ T^-3/2)")
print("="*70)
print(" Justification test: for a SINUSOID in WHITE noise, SNR^2 = h^2 T/(2 S_n) is")
print(" strictly LINEAR in T (A2). T^3 requires per-sample sensitivity ALSO improving ~T,")
print(" i.e. red-noise / spin-down-limited spectra where the growing 1/T frequency bin")
print(" suppresses in-band noise. Not our regime (Stage C freezes noise at truth, white-")
print(" dominated 3-20 nHz band). => T^3 is an OPTIMISTIC bound, not the baseline. Reported")
print(" only to bracket; even it does not rescue the wall.")
for name, rung in [("loosest rung 1.85e-3", RUNG_LOOSE), ("L2c pull-in 1e-4", RUNG_L2C)]:
    T = T_cross(rung, p=1.5)
    print(f"   {name}:  T = T0*(SIG0/target)^(2/3) = {T0}*{(SIG0/rung)**(2/3):.4g} = {T:,.0f} yr")

# =============================================================== (5) figure
T_grid = np.logspace(np.log10(10), np.log10(2e4), 400)      # 10 yr .. 20 kyr
fig, ax = plt.subplots(figsize=(8.4, 5.6))
ax.loglog(T_grid, sigma_sky(T_grid, 0.5), color="#1f77b4", lw=2.2,
          label=r"$\sigma_{\rm sky}(T)=0.05\,(T/15)^{-1/2}$  (SNR$^2\!\sim T$, baseline)")
ax.loglog(T_grid, sigma_sky(T_grid, 1.5), color="#1f77b4", lw=1.4, ls=":",
          label=r"$\sigma_{\rm sky}(T)=0.05\,(T/15)^{-3/2}$  (SNR$^2\!\sim T^3$, optimistic)")

# rungs as horizontal lines
ax.axhline(RUNG_LOOSE, color="#d62728", lw=1.6, ls="--",
           label=f"loosest F2 rung  1.85e-3  (T={T_cross(RUNG_LOOSE):,.0f} yr)")
ax.axhline(RUNG_L2C,   color="#9467bd", lw=1.6, ls="--",
           label=f"L2c pull-in  1e-4  (T={T_cross(RUNG_L2C):,.0f} yr)")

# anchor + float-floor band
ax.plot(T0, SIG0, "o", color="black", ms=8, zorder=5)
ax.annotate(f"anchor: 0.05 scaled (6.4$^\\circ$)\n@ T=15 yr, SNR=10.7",
            (T0, SIG0), textcoords="offset points", xytext=(10, 12), fontsize=9)
ax.axhspan(SIG0, 1.0, color="grey", alpha=0.10)
ax.text(11, 0.16, "achievable blind-float region\n(sky-localisation floor)",
        fontsize=8, color="grey")

# crossing markers on baseline curve
for rung, c in [(RUNG_LOOSE, "#d62728"), (RUNG_L2C, "#9467bd")]:
    Tc = T_cross(rung, 0.5)
    ax.plot(Tc, rung, "v", color=c, ms=9, zorder=6)

ax.set_xlabel("dataset span  T  [yr]")
ax.set_ylabel(r"blind sky precision  $\sigma_{\rm sky}$  [scaled]")
ax.set_title("Registration wall: when does the Earth-term float reach the first rung?")
ax.grid(True, which="both", alpha=0.25)
ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
ax.set_xlim(10, 2e4); ax.set_ylim(5e-5, 1.0)
fig.tight_layout()
fig.savefig("CW_transition/PENCIL_t_crossing.png", dpi=140)
print("\n[fig] CW_transition/PENCIL_t_crossing.png")

# =============================================================== save npz
np.savez("CW_transition/PENCIL_t_crossing.npz",
    T0=T0, sigma0_scaled=SIG0, sigma0_deg=SIG0_DEG, snr0=SNR0, logh0=LOGH0,
    deg_per_scaled=DEG_PER_SCALED,
    rung_loose=RUNG_LOOSE, rung_l2c=RUNG_L2C,
    # step2 baseline SNR^2~T
    Tcross_loose_p05=T_cross(RUNG_LOOSE, 0.5), Tcross_l2c_p05=T_cross(RUNG_L2C, 0.5),
    # step4 stress SNR^2~T^3
    Tcross_loose_p15=T_cross(RUNG_LOOSE, 1.5), Tcross_l2c_p15=T_cross(RUNG_L2C, 1.5),
    # step3 strain thresholds
    snr_thresh_loose=SNR0*(SIG0/RUNG_LOOSE), snr_thresh_l2c=SNR0*(SIG0/RUNG_L2C),
    logh_thresh_loose=LOGH0+np.log10(SIG0/RUNG_LOOSE),
    logh_thresh_l2c=LOGH0+np.log10(SIG0/RUNG_L2C),
    assumptions="A1 persistent monochromatic CW; A2 white-noise-dom, fixed cadence, no "
        "new pulsars => SNR^2~T; A3 sigma_sky~1/SNR x fixed geometry; A4 sky binds "
        "(freq ~T^-3/2 irrelevant); p=0.5 baseline, p=1.5 optimistic red/cadence bound.")
print("[npz] CW_transition/PENCIL_t_crossing.npz")
