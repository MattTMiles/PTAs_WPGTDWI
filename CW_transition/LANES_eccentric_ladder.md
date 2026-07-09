# LANES — Eccentric harmonic-comb lane ladder (E-track scoping)

**Agent:** LANES · **Date:** 2026-07-07 · CPU-only toy (numpy/scipy). No likelihood
machinery, no GPU. This SCOPES the full E-track run; it does not replace it.

**Question.** Track B's NO-GO rests on a lane gap: the loosest fringe-registration
rung on the array is **1.85e-3 scaled** (F2) vs the blind-search float ceiling
**~0.05 scaled** — a **27× span with zero rungs between**. E-track hypothesis: an
ECCENTRIC binary emits a phase-locked Peters–Mathews harmonic comb at n·f_orb; the
registration tolerance scales as 1/f, so the **orbital fundamental (n=1)** gives a lane
n_peak× wider than the highest harmonic — a *vernier* that might populate the gap.

## BRIDGE VERDICT — **NO-GO** (does not span the 27× gap; no optimal-e window)

- **Minimum e to reach the float ceiling: none.** No eccentricity brings the widest
  *usable* lane to 0.05 scaled, in either weighting.
- Widest reachable lane peaks at **~8.6e-3 scaled at e=0.9** (power-anchor) — still
  **~5.8× short** of 0.05. Physically (residual-SNR anchor) the widest lane is
  **~1.85e-3–3.7e-3** — essentially the F2 rung itself, ~13–27× short.
- **No optimal-e window** bridges. e≈0.5 optimally maximises lane *count* (5 usable
  rungs) but every one lands in the *already-populated fine band* (≤1.85e-3), not the
  empty coarse gap.
- **What eccentricity buys:** finer rungs (precision, the L2c pull-in end), **not**
  wider rungs (capture, the float-gap end). It densifies the ladder Track B already
  has; it does not extend it upward.

The negative is a **scissors** — two independent mechanisms, either alone fatal:

1. **Power-anchor geometry.** To sit the power-dominant harmonic at 27× (n_peak≈27)
   needs e≳0.85. But at high e the fundamental n=1 falls **below 5% of the peak**, so
   it is not usable; the >5%-of-peak band spans only ~5–14× in n and is centred high.
   Widest usable lane = 1.85e-3·(n_peak/n_min) tops out at ~8.6e-3.
2. **Residual-SNR reality.** A PTA detects via *timing residuals* r ~ h/(2πf), so
   per-harmonic SNR² ∝ g(n,e)/n⁴. The SNR-dominant harmonic is therefore **always
   n=1–2**, regardless of e. The detection band *is* the fundamental, so F2's 1.85e-3
   anchor already sits at the fundamental — there is no wider lane above it to climb to.

Both point the same way: the lever the E-track hoped for (detect high, register low)
**does not exist**, because the same 1/n⁴ that would make the fundamental loud also
makes it the detection band.

## Literature — the angle is NOVEL / UNCLAIMED (but the toy says it doesn't bridge)

- **Taylor, Huerta, Gair & McWilliams 2016** (arXiv:1505.06208) — the eccentric
  resolvable-source paper — **explicitly excludes the pulsar term** and defers it to
  future work; distances are "searched or marginalized," likelihood "highly sensitive"
  to distance. Harmonics used for **detection SNR / PE only**, never distance/fringe.
- **Eccentric CW searches** (Taylor+2016 SNR arXiv:1504.00928; Huerta; Sesana): comb
  used for sensitivity; pulsar distance marginalized. No comb-as-vernier. UNCLAIMED.
- **EOB paper arXiv:2511.19611** (Manzini & Babak, 2025-11): maps the "rich harmonic
  structure," treats the pulsar term, but states it cannot determine the pulsar phase
  "given poor knowledge of pulsar distances." Does **not** link comb→distance. UNCLAIMED.
- **CPTA / others:** no eccentric-comb distance claim. Closest is Xiao+2025
  (arXiv:2512.10729) — L_p from the pulsar term of **multiple** CGW sources (different
  mechanism). Zhu+2016 (arXiv:1606.04539) is the canonical circular result: pulsar term
  refines L_p only if prior σ_L < λ_GW/(1−cosθ) (~0.66 pc J0437 @10nHz; still ~3× short).
- **Verdict:** the eccentric harmonic-comb *vernier for PTA pulsar-distance
  certification* is **novel/unclaimed** — the comb (detection) and pulsar-term distance
  (circular) exist as separate mature threads; no one joins them. The idea is worth a
  footnote of priority, but this toy shows it does **not** close the F2 gap.

## Part 2 — Lane spectrum (Peters–Mathews g(n,e), 5%-of-peak usable band)

Anchor A2: power-dominant harmonic n_peak placed at f_gw = F2's 1.85e-3 rung; rung
tolerance_n = 1.85e-3·(n_peak/n). "widest" = loosest usable rung (at n_min).

| e | n_peak | usable n | span | fund. usable? | widest (scaled) | reaches 0.05? |
|----|----|----|----|----|----|----|
| 0.1 | 2 | 2–3 | 1.5× | no | 1.85e-3 | no |
| 0.3 | 3 | 2–6 | 3.0× | no | 2.78e-3 | no |
| 0.5 | 4 | 2–11 | 5.5× | no | 3.70e-3 | no |
| 0.7 | 10 | 3–27 | 9.0× | no | 6.17e-3 | no |
| 0.9 | 51 | 11–153 | 13.9× | no | 8.58e-3 | no |

Residual-SNR-weighted band (the physical one): n_peak = 1–2 for all e; fundamental
usable for e≥0.3; **widest lane 1.85e-3–3.7e-3** — never reaches 0.05.

## Part 3 — SNR tax (per-lane cert budget 0.5·SNR_n²·(1−M), SNR_tot=15, (1−M)~1)

| e | usable lanes | min per-lane budget (nat) | widest (scaled) | reaches? |
|----|----|----|----|----|
| 0.1 | 1 | 107 | 1.85e-3 | no |
| 0.3 | 3 | 15.0 | 3.70e-3 | no |
| **0.5** | **5** | **2.9** | 1.85e-3 | no |
| 0.65 | 5 | 4.1 | 1.85e-3 | no |
| 0.8 | 1 | 98.5 | 1.85e-3 | no |
| 0.9 | 2 | 8.0 | 1.85e-3 | no |

SNR is **not** the binding constraint: even at the 5-lane optimum (e≈0.5) the weakest
rung keeps ~3 nat — comfortably fixable. The optimum e for lane *count* is ~0.5–0.65,
but those lanes tile the fine band (≤1.85e-3), so the healthy budget buys precision the
ladder already has, not reach. Above e≈0.7 power concentrates in a narrow high-n band →
lanes collapse to 1–2. **There is an SNR-optimal e (~0.5); there is no bridge-optimal e.**

## Part 4 — Evolution envelope (pulsar-lag fundamental de-coherence)

Closed-form Peters (1964) e(t),a(t) over the pulsar-term lag τ_p; K_eff ≈ π/Δφ_cycle
(resolvable fundamental fringes before Earth↔pulsar de-coherence). Over τ_p = 3 kyr the
fundamental lane collapses (**K_eff → 1**, Δφ_cycle > π) **only in the high-Mc /
high-f_orb / high-e corner** (Mc ≳ 1e9 M⊙, f_orb ~ 1e-8 Hz, e ≳ 0.7). For the bulk of
(Mc ≤ 1e8, or f_orb ≤ 3e-9) the fundamental stays coherent (K_eff ≫ 1). So evolution
does **not** independently rescue *or* kill the idea across most of parameter space: the
fundamental lane exists (good for a vernier in principle) but, per Parts 2–3, is not
wide enough to bridge; evolution only removes it in the loud/fast/massive corner.

## Figures

- `lanes_ladder.png` — lane rungs vs e over the F2 gap (shaded 1.85e-3→0.05). Even
  e=0.9 tops out ~8.6e-3; physical (★) lanes hug 1.85e-3. The gap stays empty.
- `lanes_kcontour.png` — K_eff(fundamental) over (e, f_orb) for Mc∈{1e8,1e9,5e9};
  K_eff=1 (red) and K_eff=3 (white) contours; collapse only in the top-right corner.

## Approximations (all stated; this scopes, not settles)

- **A1** tol ∝ 1/f (F2 lever-arm geometry, verified against trackB_F2_ladder.py).
- **A2** fair-detectability anchor: power-dominant harmonic at f_gw = F2 rung. *This
  anchor is the crux*; the residual-SNR framing (A4) is shown alongside and agrees.
- **A3** Peters–Mathews 1963 g(n,e) via Bessel J_n(ne); F(e) enhancement.
- **A4** residual SNR² ∝ g(n,e)/n⁴ (white noise). **Optimistic for high n** — real red
  noise / GWB rising to low f pushes the usable band *lower*, strengthening the NO-GO.
- **A5** cert budget 0.5·SNR²·(1−M), (1−M)~1 for resolved teeth (toy).
- **A6** closed-form Peters inspiral for the evolution channel; K_eff ≈ π/Δφ_cycle.

**Escape hatches for the full E-track run (flagged, not pursued here):** (i) an
anomalously aligned pulsar (1−cosμ→0) with an unusually loose fundamental — array-
specific, check the real ladder tail; (ii) a strain-limited (non-residual) detection
regime — not the PTA case. Absent these, eccentricity aids **precision, not capture**;
it does not reopen Track B.

*Files:* `lanes_eccentric_ladder.py`, `lanes_figures.py`,
`lanes_eccentric_ladder.npz`, `lanes_ladder.png`, `lanes_kcontour.png`,
`logs/lanes_run.log`. (Agent LANES)
