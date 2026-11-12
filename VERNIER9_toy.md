# VERNIER-9 — raise the ceiling

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~6 min.
Clean room: imports only `vernier_toy` / `vernier6b` / `vernier7` / `vernier8`.

Code `hpc_harbor/vernier/vernier9.py` · data `reports/vernier9_wheel.npz`.
Staged, not committed.

---

## Lead — the ceiling does rise, the wall does move, and it buys 1–3 ranks

The key algebra, which VERNIER-8's confirmed budget makes exact:

> r = (f_P/f_E)^(−1/3) and f_P/f_E = (1 + τ/t_c)^(−3/8), so **r = (1 + τ/t_c)^(1/8)**,
> and τ = L(1−cos μ)/c scales with **pulsar distance**.

So r is a per-pulsar quantity. Giving the array a realistic 0.3–5 kpc log-uniform spread:

| f_E | t_c (yr) | med τ/t_c | r median | r min–max | **array ceiling** | vs single-distance | **final rank** |
|---|---|---|---|---|---|---|---|
| 10 nHz | 6484 | 0.49 | 1.051 | 1.00–1.23 | **0.1660** | −0.027 | — |
| 15 | 2199 | 1.43 | 1.117 | 1.00–1.39 | 0.1830 | −0.010 | **9** |
| **20** | 1021 | 3.24 | 1.198 | 1.01–1.52 | **0.1996** | +0.007 | **10** |
| 25 | 563 | 5.82 | 1.271 | 1.01–1.63 | 0.2157 | +0.023 | **10** |
| 32 | 292 | 11.1 | 1.366 | 1.03–1.77 | 0.2364 | +0.044 | **10** |
| 50 | 89 | 36.5 | 1.573 | 1.07–2.06 | **0.2779** | +0.085 | **12** |

*(VERNIER-8 single-distance reference: r = 1.196, ceiling 0.1929 dex, final rank 9.
Step to rank 10 = 0.1980 dex.)*

**The ceiling crosses the rank-10 step at ~20 nHz.** Below that the distance spread makes
things *worse* than the single-distance figure (0.166 at 10 nHz vs 0.1929), because most of a
log-uniform array sits at small τ/t_c where r → 1. Above it, the ceiling climbs steeply and
the wheel reaches rank 10 at 20–32 nHz and rank 12 at 50 nHz.

**The formula is validated by the simulation, not just asserted**: the wheel's measured max
gain tracks the analytic array ceiling to <2% at every frequency (0.1986 vs 0.1999 at 20 nHz;
0.2165 vs 0.2156 at 25; 0.2846 vs 0.2810 at 50).

**So VERNIER-8's 0.005 dex wall is real but not robust — it is an artefact of evaluating r at
one distance and one frequency. It moves. It just does not move far: +1 to +3 ranks.**

## The array-design question: the answer is no

**The most distant quartile contributes 28–32% of the weighted gain, against its 25%
population share.** That is a ~1.2× enrichment, not a concentration.

| f_E | 10 | 20 | 32 | 50 nHz |
|---|---|---|---|---|
| far-quartile share of Σ r²ρ_E² | 27.9% | 30.4% | 31.4% | 31.9% |

**This is a negative worth stating on its own: r ∝ (1+τ/t_c)^(1/8) is a one-eighth power, so a
15× spread in τ produces only a 15^(1/8) = 1.4× spread in r.** Distant pulsars are better, but
the exponent is too weak for them to dominate. **Preferentially timing distant pulsars is not
an effective lever** — the measured r range across the whole array is only 1.00–2.06.

## Physical consistency — the achievable range stops around 25–32 nHz

Large evolution requires being near merger, and the required proximity is severe:

| f_E | t_c | sources per bin vs 10 nHz | chirp over 30 yr |
|---|---|---|---|
| 15 nHz | 2199 yr | **2.9× fewer** | 0.5% |
| 20 | 1021 yr | 6.3× fewer | 1.1% |
| 25 | 563 yr | 11.5× fewer | 2.1% |
| 32 | 292 yr | 22× fewer | 4.2% |
| **50** | **89 yr** | **73× fewer** | **16.7%** |

Sources per frequency bin scale as the residence time t_c ∝ f^(−8/3) Mc^(−5/3), so the
frequencies where the ceiling is highest are exactly where the sources are rarest.

**Flagging plainly: the 50 nHz row requires catching a 10^9.5 M☉ binary 89 years from merger,
chirping 17% during the observation. That is not an achievable operating point** — such a
source is both 73× rarer and actively leaving the band. **The achievable range is ~20–32 nHz**,
where t_c is 300–1000 yr, the chirp is 1–4%, and the ceiling is 0.20–0.24 dex. That buys
**rank 10 instead of rank 9.**

The favourable correlation is worth noting: h ∝ Mc^(5/3) and t_c ∝ Mc^(−5/3), so the *loudest*
sources are also the *most evolved* — the large-r regime and the detectable regime point the
same way. But both are the rare corner.

## Null check — done in w_true, per VERNIER-8

VERNIER-8 found the scrambled null itself becomes confident at high SNR, so K_eff is not a
safety metric. Comparing in w_true instead, at ρ_P = 6 (near-perfect coherence, R = 1.000):

| f_E | 15 | 20 | 25 | 32 | 50 nHz |
|---|---|---|---|---|---|
| w_true signal | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| w_true null | 0.0011 | 0.0017 | 0.0021 | 0.0017 | 0.0012 |

**Separation is ~10³× at every frequency.** The weights are on the true fringe and the null is
not; K_eff would not have shown this.

*(Bug found and fixed: the source-level scramble offset could push a nearby pulsar's trial
distance negative, giving `(1+τ/t_c)^(5/8)` a negative base and a nan null at 32/50 nHz. Trial
distances are now floored at 50 pc.)*

---

## What this does not settle

1. **The gain is still bounded by ~0.28 dex even at 50 nHz**, and the census step grows
   without limit with depth (0.20 at rank 10, 0.34 at rank 16). The wheel reaches a *deeper*
   wall, not no wall. Ranks 9 → 10 in the achievable band.
2. **The ceiling is computed at Mc = 10^9.5 throughout.** At the band-median Mc ≈ 10^9 the
   t_c values are 46× longer, τ/t_c 46× smaller, and r collapses toward 1 — the whole effect
   is a heavy-chirp-mass phenomenon.
3. **The distance distribution is log-uniform 0.3–5 kpc**, not the actual IPTA spread. The
   1/8 power makes the answer insensitive to this, which is itself the reason the
   array-design lever is weak.
4. **Still no GWB, sources pre-resolved, white noise, known sky positions.** Same optimistic
   frame as VERNIER-6 through 8.
5. **The DOF effect from VERNIER-8 (0.4–0.5 dex) remains unsimulated and remains larger than
   everything measured here.** Raising the ceiling by 0.085 dex is a second-order fix next to
   it.

## One-line verdict

**The ceiling is a per-pulsar quantity and it does rise — from 0.166 dex at 10 nHz to 0.278
at 50 nHz, crossing VERNIER-8's rank-10 wall at ~20 nHz and validated against the simulation
to better than 2% — but r = (1+τ/t_c)^(1/8) is a one-eighth power, so a 15× distance spread
buys only 1.4× in r, the distant quartile carries 32% rather than a dominant share, and the
frequencies where the ceiling is highest need sources 89–292 years from merger that are
22–73× rarer; in the achievable 20–32 nHz band the wheel goes from rank 9 to rank 10.**
