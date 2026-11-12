# VERNIER-2 — walk the wheel to failure

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~12 min compute.
Clean room preserved: imports only `vernier_toy` (my own VERNIER-1 implementation,
anchored to 0.42 % on the banked comb). Nothing from the campaign.

Code `hpc_harbor/vernier/vernier2.py` · data `reports/vernier2.npz` ·
plots `reports/VERNIER2_{fcert,dials,gearing,topk}.png`. Staged, not committed.

---

## Lead

**The brightness gap binds first, and it binds alone.** It is the only dial that stops
the wheel *turning* anywhere inside its realistic range: 5 turns at 0 dex, 4 at 0.15,
2 at 0.3, **stalled at 0.5 dex**. Every other dial either never stalls it (N_psr, prior)
or stalls only at an absurd endpoint (1 comparable source).

**The wheel engages for Δ ≤ 0.5 dex on top-3 spread, and it is cheap.** Rejection cost
is **4.0 seeds per accept at Δ = 0.15 dex** and 24.6 at Δ = 0.05. The natural top-3
spread of a 256-member, σ_logh = 0.6 census is **0.258 dex** (16/84 = 0.111/0.493), so
top-3 conditioning is buying a factor of a few in rarity, **not orders of magnitude**.

**And the wheel does buy sources — about 3 of them.** That is the biggest change from
VERNIER-1 and it is not because the gearing got bigger. Gearing is ≤ 0.19 dex. The yield
comes from **how many sources sit within that 0.19 dex of threshold**, which is a property
of the luminosity function, not of the wheel.

> **The wheel's output is set by the density of sources near threshold, not by the size of
> the gearing.** Under top-3 conditioning ~3 sources land in the window and the wheel buys
> them. On A-SKY as-built, nothing lies within 0.19 dex of threshold and the wheel buys
> **zero**, despite certifying 111 of 116 pulsars.

### Dial ranking, most binding first

| rank | dial | stalls *turning* at | stalls *certification* at | note |
|---|---|---|---|---|
| **1** | **brightness gap** | **0.5 dex** (5→4→2→1 turns) | halves it by 0.15 dex (100→44) | the only real binder |
| 2 | distance prior | never (5 turns at every tier) | **1 % of L** (2/116), 0 at 10 % | kills the *output*, not the turning |
| 3 | noise | never (<2 not reached; GWB→2) | GWB: 98→**8** certified | GWB costs 3 turns and 92 % of certifications |
| 4 | n comparable | 1 source | — | endpoint only, not informative |
| 5 | N_psr | never (50/116/250 all 5 turns) | scales ∝ N (48/96/140) | **not binding over 50–250** |

The prior and noise dials do not stop the wheel turning because at gap 0 every source is
individually detectable — **the turning is supplied by the sky, and the dials only remove
the wheel's product.** For a conditions-of-validity box those are two different failures
and both belong in it.

---

## Gearing: measured, and it saturates

**Verdict: it saturates. It does not compound.** Confirmed, and now with the mechanism.

Gearing is not a constant. It is

> **gain(f_cert) = log₁₀ √(1 + r² f_cert)**,  r = ρ_P/ρ_E, f_cert = fraction of the array's
> pulsar-term **power** certified.

| f_cert | equiv. pulsars of 116 | gain, r = 1.02 | gain, r = 1.20 |
|---|---|---|---|
| 0.0086 | **1** | 0.0019 | **0.0027** |
| 0.05 | 5.8 | 0.0109 | 0.0150 |
| 0.25 | 29 | 0.0499 | 0.0664 |
| 0.50 | 58 | 0.0904 | 0.1172 |
| 0.70 | 81 | 0.1181 | **0.1507** |
| 1.00 | 116 | 0.1540 | **0.1929** |

**Addendum item 1 confirmed.** Certifying one pulsar of 116 buys **0.0027 dex**, not 0.15.
The 0.15 dex of VERNIER-1 was the f_cert = 1 corner misread as a per-turn gain.

**f_cert required to close a gap:** 0.05 dex → 0.18; **0.15 dex → 0.70** (r = 1.20) or 0.96
(r = 1.02); **0.30 dex and beyond → impossible at any f_cert.** The ceiling is 0.19 dex, so
one natural census rung (0.15–0.2 dex) is *just* reachable and two rungs never are.

**Measured per-turn, against the honest baseline** (8 sources at 0.15 dex, 1 pc, 116 psr):

| turn | certified | f_cert | gain vs **discard** | gain vs **marginalise** |
|---|---|---|---|---|
| 1 | 0 | 0.000 | 0.0889 | **0.0000** |
| 2 | 24 | 0.238 | 0.0707 | 0.0417 |
| 3 | 26 | 0.243 | 0.0601 | 0.0425 |
| 4 | 39 | 0.343 | 0.0752 | 0.0547 |
| 5 | 53 | 0.446 | 0.0907 | 0.0805 |

**Addendum item 3 matters a lot.** At turn 1, with *zero* pulsars certified, the discard
baseline already shows 0.0889 dex of "gain" — that is pure artefact of throwing the
uncertified pulsar terms away. Marginalising them instead (ρ²_P,marg = max(0, ρ_P² − 2 ln K_a))
gives the correct 0.0000. **The honest marginal gain is roughly half the discard figure**
and it simply tracks f_cert; there is no compounding term anywhere.

Why it cannot compound: certification unlocks *one* reservoir — pulsar-term coherence.
Once a pulsar's term is in the coherent sum, certifying it again buys nothing, and there is
no second reservoir. **This is a mechanism for LEDGER-C3's banked "reach does not
compound", not merely a restatement of it.**

**Addendum item 2 confirmed and my VERNIER-1 exponent corrected.** PTA sensitivity lives in
*residual* space: residual ∝ h/f and h ∝ f^{2/3}, so residual ∝ **f^(−1/3)**. I had used the
strain exponent f^(+2/3). At f_p/f0 = 0.584 the pulsar term is therefore **1.196× the Earth
term, not 0.70×**, and the full-array ceiling is **0.1929 dex, not 0.1505**. Measured over a
realistic band/Mc population r = 1.016 (16/84 = 1.002/1.057), ceiling 0.1533 dex
[0.1511, 0.1579] — barely above √2 because most draws barely evolve.

### Target selection (addendum item 4) — the μ effect is real but does not bind at 1 pc

Power per fringe is analytic: w ∝ sin²μ and K_a ∝ (1−cos μ), so **w/K_a ∝ (1 + cos μ)**,
monotonically falling. Best targets are low-μ. But:

| K_a cap | sources needed | f_cert | pulsars | gain |
|---|---|---|---|---|
| < 1 | prior alone | **0.032** | 13/116 | 0.0097 |
| < 12 | 3 | 1.000 | 116/116 | 0.1929 |
| < 250 | 5–6 | 1.000 | 116/116 | 0.1929 |

**I have to report this against my own expectation.** At σ = 1 pc the K_a cap is *not*
binding: dL spans 0.3–15 pc so K_a = 6σ/dL < 12 for essentially the whole array, and f_cert
→ 1 for every cap ≥ 12. The μ-selection tension I predicted — "certifiable pulsars are the
low-weight ones" — **only bites in the sub-fringe regime**, where K_a < 1 admits just 3.2 %
of pulsar-term power (13 pulsars, 0.0097 dex). And K_a < cap is **necessary, not
sufficient**: certification also needs the joint comb intersection to collapse. The wheel's
**measured** f_cert of 0.24–0.45 is the real reservoir; the K_a proxy is optimistic.

---

## What changed once the pulsar term was allowed to evolve

**It eases the ambiguity, but by a quarter of a source, not by 1.71×.**

| tier | K_trials (K=1) static → evolving | modes at K=1 | K₉₀ (fractional) | shift |
|---|---|---|---|---|
| vlbi (1 pc) | 4.5 → 3.9 | 5 → 3 | 2.59 → **2.32** | **+0.27 sources easier** |
| 20 % of L | 1368 → 1189 | 1367 → 1040 | 5.56 → **5.05** | **+0.51 sources easier** |

**Correction to the brief's premise.** The 1.71× coarser comb holds at f_p/f0 = 0.584, which
is the *extreme corner* (20 nHz, Mc = 10⁹·⁵). Over a realistic population (3–20 nHz,
log10 Mc 8.5–9.5) the comb is coarser by only **1.15×** median, because most draws have
τ/t_c ≪ 1 and barely evolve. So the easing is real but sub-integer in K — VERNIER-1's
K₉₀ = 3 at 1 pc stands, now at 2.32 fractionally. The clearest single effect is at
K = 4 with 20 % priors: P(unique) **0.15 → 0.28**, nearly doubled.

### Is the separate spectral line usable? Partly — as a prior, not as an identifier

Because f_p ≠ f0 the pulsar term is a distinct line, resolvable when |f0 − f_p| > 1/T_obs.
That happens in **6 of 20** band cells — only the high-f, high-Mc corner. Where resolvable,
the line *position* is a **monotonic, wrap-free** function of L, so unlike the comb it has no
fringe ambiguity. Precision σ_L = σ_f/|df_p/dL| with df_p/dL = (3/8) f_p(1−cos μ)/(c(t_c+τ)):

- σ_L spans **10.1 – 6.8×10⁴ pc**, median 102.7 pc; best cell (20 nHz, Mc = 10⁹·⁵) **10.1 pc**
- in comb units that is **9.7 – 2×10⁴ fringes**, median 78

**So it is a coarse unambiguous prior, not a fringe identifier** — but a *data-derived* one,
which is exactly what the programme lacks: it could replace a 235 pc feather prior with a
~10–100 pc one at no EM cost.

**The catch, and it is severe.** f_p depends on t_c ∝ Mc^(−5/3), and the banked chirp-mass
posterior is prior-dominated (gain 1.00–1.73×, σ = 0.50–0.87 dex). A 0.5 dex Mc error is
0.83 dex in t_c, which swamps the L dependence — **the handle is Mc–L degenerate for a
single pulsar.** It breaks only with the array: one shared Mc against N_psr distinct τ_a,
i.e. N_psr line positions for N_psr + 1 unknowns. **Worth a dedicated look; this is the one
genuinely new distance handle the toy turned up.**

---

## The two skies

**(a) A-SKY as built** — 5 loud at log10 h = −12.75 + 11 faint at −14.25, 30 realisations:

| turn | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| median certified | 0 | 44 | 53 | 81 | 100 | 111 |
| median gain (dex) | 0.106 | 0.122 | 0.124 | 0.134 | 0.143 | 0.163 |

Turns 5/5 in every realisation. Peak certified **111/116**. **Geared turns: 0, max 0.**

This is the sharpest test and it splits cleanly: **the vernier works — 111 of 116 pulsars
certified — and the wheel still buys nothing.** All 5 loud members clear threshold on their
own, so the sky supplies every turn; the 1.5 dex cliff to the faint tier is 8× the 0.19 dex
ceiling. A-SKY demonstrates certification at scale and gearing yield of exactly zero.

**(b) Top-3 conditioning at fixed band power** — 256 members, σ_logh = 0.6 dex,
400 000 skies drawn:

| Δ dex | P(accept) | **seeds/accept** | turns | geared turns |
|---|---|---|---|---|
| 0.05 | 0.0407 | **24.6** | 16 | 3.5 |
| 0.10 | 0.1345 | **7.4** | 16 | 3.5 |
| 0.15 | 0.2502 | **4.0** | 13 | 3.0 |
| 0.20 | 0.3698 | **2.7** | 14 | 3.0 |
| 0.30 | 0.5829 | 1.7 | 13 | 3.0 |
| 0.50 | 0.8463 | 1.2 | 7 | 2.0 |

**The wheel engages at every Δ ≤ 0.5, buying 2–3.5 geared turns, and the rarity is mild.**
Natural top-3 spread is 0.258 dex median, so demanding 0.15 dex costs 4 seeds and demanding
0.05 dex costs 25. **Plug straight into LOTTERY-2: top-3 conditioning is a ~4–25× calendar
penalty, not a 10³× one.** Contrast the 0.77–1.17 dex that rank-0-only conditioning
produces — conditioning on rank 0 actively *pushes the top three apart*, which is precisely
the wrong target. **Change the conditioning variable.**

---

## What this cannot settle

1. **The gearing yield depends on a luminosity function I made up.** "~3 geared turns" is
   the count of sources within 0.19 dex of threshold in my ladder. That number is a
   *census* measurement, and it is now the single most load-bearing input. GLACIER/CHORUS
   must supply dN/dlog h near threshold.
2. **The marginalisation penalty is a caricature.** ρ²_P,marg = max(0, ρ_P² − 2 ln K_a) is a
   trials-factor argument, not a computed marginal likelihood — and the banked R referendum
   showed real fringe-marginalised posteriors *smear* in ways this does not capture.
3. **Noise is a two-parameter caricature** (r_red, r_gwb set so each is 3× white at 10 nHz),
   with no cross-pulsar GWB correlation and no timing-model marginalisation. The GWB dial's
   teeth (98→8 certified) are qualitatively right and quantitatively soft.
4. **Detection is a threshold, not a search.** No trials factor over sky/frequency for
   *detection*, so "geared turn" is optimistic.
5. **The f_p distance handle is unvalidated.** σ_L came from a Fisher-style σ_f ≈ 1/(2πTρ),
   not from a likelihood scan, and the Mc–L degeneracy was argued rather than fitted.
6. **Still Asimov.** No noise realisations, so nothing here has a false-alarm rate and no
   count is a certification in the programme's sense. Flicker remains invisible.

## One-line verdict

Gearing is small (≤ 0.19 dex), saturating, and now mechanistically explained; the brightness
gap is the only dial that binds; and the wheel's yield is set by how many sources sit within
0.19 dex of threshold — so **stop trying to make the gearing bigger and start conditioning
the sky on its top three, which costs only ~4 seeds at 0.15 dex.**
