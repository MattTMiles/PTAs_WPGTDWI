# EXPLAINER — six figures, for a scientist outside the field

**Agent:** EASEL-2 · ACCRE · tag `criterion-v2.1+` (`git rev-parse HEAD` → `6bec3d6`) ·
**Date:** 2026-07-13 · **READ-ONLY.** No likelihood evaluation, no new realisations. Every
number below is read out of a banked `.npz`, or re-cut from banked raw statistics using the
campaign's own definitions. Script: `EXPLAINER_results/explainer_figures.py`.

## The one-paragraph story

A pulsar timing array watches ~116 millisecond pulsars and looks for a single supermassive
black-hole binary. The signal reaches us twice: once now (the "Earth term"), and once as it
passed each pulsar, thousands of years ago (the "pulsar term"). The old copy is enormously
valuable — it gives a kiloyear-long lever arm on the binary's chirp mass, which is what
converts the source into a *standard siren* that measures cosmological distance. To use it you
must know how far away each pulsar is to within a fraction of a gravitational wavelength
(~0.1 pc). Nobody does. What follows is: the problem (F1), why the obvious fix fails (F2), how
we learned to tell a real detection from a fluctuation (F3), the one thing that turns the
measurement on (F4), the map of when it works (F5), and what it is worth (F6).

**Everything here is simulated** — a mock 116-pulsar array with injected sources. No real
pulsar data is analysed anywhere in this repository.

---

### F1 — THE PROBLEM · `F1_the_problem.png`

> The gravitational wave arriving from the pulsar's direction is periodic, so the data fit the
> pulsar equally well at *any* distance that shifts the wave by a whole number of cycles: the
> likelihood is a comb of hundreds of near-identical peaks (a). Astronomers already measure
> pulsar distances by other means, and that prior knowledge (red band, ±1 parsec) narrows the
> choice — in the easy case shown in (b) it nearly isolates a single tooth. But in the real
> array the comb is far finer than the prior is narrow: the typical pulsar has ~6,000 candidate
> teeth inside its own distance error bar (c), and only the nearby anchor J0437−4715, with just
> 6, is close to unambiguous. *"Which tooth?" is the whole problem.*
>
> Source: banked distance scan `lnLs_GWAmp_phase_connected/runD_3CW_test/J0437-4715_cw_p_dist/`
> (5,000-point likelihood scan, 0.1–5.0 kpc); distance priors `CW_transition/anchor_a0_priors.npz`;
> candidate-slot counts `K_counted` from `reports/geo_dlnl_bank.npz` (median over 40 sky draws).

### F2 — THE WALL · `F2_the_wall.png`

> The obvious fix is to find the source first, then use its position to work out which tooth
> each pulsar sits on. This figure shows why that fails. Each blue bar counts pulsar–source
> pairs by how precisely you would need to already know the source's sky position for the tooth
> identification to survive: all 348 pairs need better than 1.9×10⁻³ (green), while a blind
> search over the same data delivers only 0.05–0.21 (red). **Between the two lies a gap 27× to
> 112× wide containing zero pairs** — there is no partially-good pulsar you could solve first and
> bootstrap from, which is why the cold-start version of this measurement was ruled out.
>
> Source: `CW_transition/trackB_F2_ladder.npz` (348 pulsar–source registration tolerances, pure
> geometry). Blind-search bracket 0.05 → 0.21 and the pull-in basin < 10⁻⁴ are the banked Track-B
> close-out numbers.

### F3 — THE HONEST CRITERION · `F3_the_honest_criterion.png`

> Certifying a pulsar means claiming its evidence is too strong to be a fluke — which requires
> knowing how strong a fluke can get. Grey: 100 datasets containing **no source at all**, scored
> identically; noise alone routinely manufactures 15–20 nats of apparent evidence, because we
> take the best of ~116 pulsars × thousands of teeth. The earlier bar (red dashed) was set by the
> loudest of just *ten* noise trials, landed at 13.5, and yielded 1.03 "certifications" per
> dataset — it looked like a detection. Fitting the bar properly, so that pure noise clears it
> only 5% of the time (green, 19.5 ± 1.4 from 100 trials), the same signal data yields **0.37 —
> and the detection is retracted.** *Confidence never lied; the comparison did.*
>
> Source: re-cut from `SURFACE_results/sf_{nullN,sig}_h1300_T30_lit_k3_*.npz` (100 noise-only +
> 30 signal realisations) at cell h = 10⁻¹³, T = 30 yr. The refit reproduces SURFACE's banked
> floor (19.46 ± 1.40) and its retraction count (0.37) exactly.

### F4 — THE SWITCH · `F4_the_switch.png`

> The census population is three loud sources on circular orbits, and it does not work: 0.70
> certified pulsars per dataset, below the bar of 1 (a, leftmost). Now make **one** of those three
> orbits oval — same loudness, same array, same observing time. A slightly oval orbit (e = 0.3)
> crosses the bar; a strongly oval one (e = 0.7) delivers 7.83, eleven times the circular count.
> An eccentric orbit emits on many harmonics instead of one, and that harmonic comb is a clock
> the analysis can lock onto. Panel (b) is the honesty check: the noise bar *rises* at the same
> time (4.3 → 8.5 nats), so the gain is a real gain and not a lowered standard.
>
> Source: `CHORUS_results/ch_analysis.npz` (counts) and `ch_floors.npz` (floors), 26 mixture cells
> × 30 realisations, floors refit per mixture. Caveat: the extra certifications are attributable
> to the eccentric member's *own* template — the clock is not shared with its circular neighbours.

### F5 — THE MAP · `F5_the_map.png`

> This is the capability map: how loud the sources are (left→right) against how long you watch
> (bottom→top), with the number of certified pulsars in each cell. The red line is the on/off
> boundary, drawn only where the count survives its own floor's calibration error. On the
> population we think we have (a), the measurement switches on only in the loud, long-baseline
> corner — and note it is *not* monotone in time: 40 years beats 50, because the longer you watch
> the higher the noise bar climbs too. Promote just two more of the sixteen sources from faint to
> loud (b) and the boundary sweeps to the faint edge of the map. **How many loud sources exist
> matters more than how loud the loudest one is** — and nobody has measured that number.
>
> Source: `SURFACE_results/surface_analysis.npz`, 108 cells × 30 signal realisations, each with
> its own 200-null floor. Panels show the "today's distance knowledge" (lit) tier.

### F6 — THE PAYOFF · `F6_the_payoff.png`

> This is what certification is *for*. Each line is one black-hole pair; the vertical axis is how
> well you could measure its distance. With zero certified pulsars the distance is essentially
> unmeasured (332% — the Earth term alone cannot separate a heavy far source from a light near
> one). Each certified pulsar term adds a kiloyear lever arm on the chirp mass and breaks that
> degeneracy: for the heavier pairs, three to five certified pulsars land inside the 10–30% band
> that dark-siren cosmology already treats as useful. The red star is F4's punchline — **an
> eccentric source reaches the same 12% from its own harmonic comb, with zero certified pulsars.**
>
> Source: `SIREN_results/siren_summary.npz` (arm B: seed distances known to 0.1 pc; 189 exact
> Fisher matrices on zero-noise data, so these are best-case bounds) and the self-clocking point
> from `reports/atlas_m2m4_summary.npz` (f_orb = 10⁻⁸ Hz, M_c = 10⁹·⁵ M_☉, e ≳ 0.58).

---

## What travels with these figures

1. **Simulated throughout.** Mock array, injected sources, stated observing-extension convention.
2. **F6 is conditional.** It prices *what certification buys*, not whether you can achieve it —
   and F2 is the standing negative on achieving it cold.
3. **Certified ≠ correct.** Every count in F4 and F5 carries a small rate of *wrong* certifications
   (a pulsar locked onto the wrong tooth), which no layer of the criterion currently removes.
4. **The two figures that point forward are F4 and F5**, and they point the same way: certification
   is a property of the *population*, not of any single source.
