# EXPLAINER — six figures, for a scientist outside the field

**Agent:** EASEL-2 · ACCRE · tag `criterion-v2.1+` (`git rev-parse HEAD` → `6bec3d6`) ·
**Date:** 2026-07-13 · **READ-ONLY.** No likelihood evaluation, no new realisations. Every
number below is read out of a banked `.npz`, or re-cut from banked raw statistics using the
campaign's own definitions. Script: `EXPLAINER_results/explainer_figures.py`.

> ## ⚠ FOUR FIGURES ARE STALE — THE FLOOR FIX (RECUT, 2026-07-13) MOVED THEIR NUMBERS.
>
> **`F3_the_honest_criterion.png`, `F4_the_switch.png`, `F5_the_map.png`, `F9_the_growth_path.png`
> must be REGENERATED from the `_recut` banks before they are shown to anyone.** The captions below
> have been corrected in place and now disagree with the PNGs they describe. **The other ten
> figures are unaffected.**
>
> **The one that must not be shown as it stands is F4/F9.** They teach *"a slightly oval orbit
> (e = 0.3) crosses the bar"* — **that is refuted.** The single-member switch-on is **e = 0.5**.
>
> **And F3 carries an irony worth stating:** the figure that teaches *"fit the bar properly"*
> illustrates it with the one cell whose properly-fitted Gumbel is **invalid** (its null is 27 %
> silent). Its lesson survives — the detection is still retracted — but its bar and its count both
> moved, and *the figure argues for an estimator the programme has since had to bound.*
>
> Regenerate with `reports/explainer_figures.py` against `recut_surface.npz` / `recut_chorus.npz`.
> [FLAGGED FOR MATT — figure regeneration is compute + design, not a doc edit, and is not done here.]

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
> **The comb in (b) is a banked legacy single-source scan and its teeth are coarser than the
> census's** — it is shown because it is the only banked likelihood-vs-distance curve in the repo;
> panel (c) is the real distribution and it is the one that carries the argument.
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
> only 5% of the time (**16.6 ± 1.6** from 100 trials), the same signal data yields **0.60 — and
> the detection is retracted.** *Confidence never lied; the comparison did.*
>
> **⚠ FIGURE STALE — and instructively so.** The PNG draws the green bar at **19.5 ± 1.4** (a
> Gumbel fit) and the count at **0.37**. Under the floor convention adopted after this figure was
> made, **that Gumbel is INVALID at this very cell**: 27 % of its noise trials produce no candidate
> at all, so the fit is describing a point mass at zero rather than a tail. The honest bar is the
> **empirical 95th percentile, 16.60 ± 1.60**, and the count is **0.60**. **The retraction — the
> entire point of the figure — stands, and stands more firmly** (0.60 is still far below 1). But
> the figure teaches the reader to trust an estimator that the programme has since had to bound,
> and it must be redrawn on the empirical bar.
>
> Source: re-cut from `SURFACE_results/sf_{nullN,sig}_h1300_T30_lit_k3_*.npz` (100 noise-only +
> 30 signal realisations) at cell h = 10⁻¹³, T = 30 yr; corrected numbers from
> `reports/recut_surface.npz` (zero-fraction 0.27 → empirical q95 adopted).

### F4 — THE SWITCH · `F4_the_switch.png`

> The census population is three loud sources on circular orbits, and it does not work: **0.37**
> certified pulsars per dataset, below the bar of 1 (a, leftmost). Now make **one** of those three
> orbits oval — same loudness, same array, same observing time. **A slightly oval orbit (e = 0.3)
> does NOT cross the bar** (0.70 — it fails). **A moderately oval one (e = 0.5) does: 3.13.** A
> strongly oval one (e = 0.7) delivers **5.43, roughly fifteen times the circular count.**
> An eccentric orbit emits on many harmonics instead of one, and that harmonic comb is a clock
> the analysis can lock onto. Panel (b) is the honesty check: the noise bar *rises* at the same
> time (**7.0 → 11.7 nats**), so the gain is a real gain and not a lowered standard.
>
> **⚠ FIGURE STALE — AND ITS HEADLINE IS REFUTED.** The PNG shows the e = 0.3 bin clearing the bar
> at 1.57. **Under the corrected floors it reads 0.70 and fails.** *"One slightly-oval orbit is
> enough"* is **not true**; **"one moderately-oval orbit (e ≳ 0.5) is enough, OR two mildly-oval
> ones (e ≳ 0.3)"** is. **Do not present this figure until it is redrawn.** The honesty check in
> panel (b) is the reason: the noise bar rose *more* than the figure knew — the e = 0.3 cell's
> bar was understated by 53 %, and it was the bar, not the signal, that was wrong.
>
> Source: `recut_chorus.npz` (counts + adopted floors), 26 mixture cells × 30 realisations, floors
> refit per mixture **and re-cut against the empirical q95** (all 26 cells fail the Gumbel validity
> gate). Caveat: the extra certifications are attributable
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
> **⚠ FIGURE STALE (mildly).** The re-cut moves **4 of the 108 cells** across the boundary — two
> onsets died at the faint edge and two were born — so **the red on/off contour is drawn in the
> wrong place at h = 10⁻¹³·²⁵**, though the total (59 cells on) is unchanged. Both headline
> readings survive exactly: 40 yr still beats 50 in the loud columns (12 of 12), and 30 yr is still
> optimal nowhere (0 of 36). **The lesson holds; the contour must be redrawn.**
>
> Source: `reports/surface_analysis_recut.npz`, 108 cells × 30 signal realisations, each with
> its own floor (Gumbel where the null's zero-fraction ≤ 20 %, empirical 95th percentile with a
> bootstrap error above it). Panels show the "today's distance knowledge" (lit) tier.

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

# The machine half — F7 to F10

F1–F6 say what the problem is and what solving it is worth. F7–F10 say how the machine works:
why a certified pulsar has any value (F7), whether certifications breed more certifications (F8),
how many you actually need (F9), and what the certified set is worth *afterwards* (F10).
Script: `EXPLAINER_results/explainer_figures2.py`.

### F7 — THE TWO CLOCKS · `F7_the_two_clocks.png`

> The wave that reaches Earth today passed the pulsar centuries ago — for J1909−3744, 692 years
> ago — so a pulsar is a second, delayed sample of the *same* waveform, and the size of that delay
> is set by how far away the pulsar is. **The delay is the distance.** That is why certification is
> worth anything: one sample of a wave cannot separate a heavy distant binary from a light nearby
> one, but two samples a kiloyear apart measure how fast the binary is spiralling in, and therefore
> its mass. Panel (c) is the payoff: the uncertainty on the binary's mass falls from 0.865 dex
> (the prior — i.e. unmeasured) to 0.326 with one certified pulsar and 0.025 with five, a 35×
> sharpening.
>
> Source: `SIREN_results/siren_geom.npz` (the 18 real lags) and `siren_summary.npz` (arm B,
> f = 10⁻⁸ Hz, M_c = 10⁹ M_☉). Zero-noise Fisher bound — a best case.

### F8 — THE LOOP · `F8_the_loop.png`

> If each certified pulsar sharpens your model of the source, the sharper source should certify the
> *next* pulsar — a cycle that ought to snowball. Panel (a) shows the cycle with its measured
> numbers: the first seed does sharpen the binary's mass (0.87 → 0.33 dex), but the noise bar the
> next pulsar must clear does not move (38.9 nats, fixed), so the only place a gain can come from is
> pulsars newly clearing a *fixed* bar. Panel (b) is every one of the 70 banked loop runs: **64 hold
> flat, 2 gain a real pulsar, and 4 spontaneously shed a false alarm** — the loop is safe,
> self-cleaning, and consolidating, with a measured gain of exactly 1.000 and never above it.
> **It does not yet snowball at today's sensitivity.**
>
> Source: all 70 banked soft-loop runs (`IGNITE2_results/ig_sloop*`, `CHORUS_results/ch_sloop*`),
> 50 clean + 20 sky-scrambled, per-iteration `traj_n_cert_true` / `traj_wrong` / `traj_gain`.

### F9 — THE GROWTH PATH · `F9_the_growth_path.png`

> The whole array is 116 pulsars, but **the science does not need the array — it needs the first
> handful.** Three certified pulsars make the source a useful standard siren *and* pin it to a
> single galaxy; five saturates (more buy nothing). Today's all-circular population delivers
> **0.37** per dataset — below even one. Making a single loud orbit **slightly** oval (e = 0.3)
> **is not enough — it reaches 0.70 and still fails.** Making it **moderately** oval (e = 0.5)
> reaches **3.13 — the three-pulsar siren threshold, exactly**; making it strongly oval (e = 0.7)
> reaches **5.43, at saturation.** **One eccentric source still supplies the entire handful the
> science needs — but it must be moderately eccentric, not slightly**, and it does so at unchanged
> loudness.
>
> **⚠ FIGURE STALE — its e = 0.3 rung is refuted** (see F4). The growth path is real; **its first
> rung is one step further up the eccentricity axis than the figure draws it.** The corrected
> reading is arguably the *cleaner* story: e = 0.5 lands the count at **3.13**, which is precisely
> the three-pulsar threshold the figure is built around — *the switch-on and the science threshold
> coincide.*
>
> Source: `recut_chorus.npz`, CHORUS mixture cells (census loudness, 30 yr, lit tier); SURFACE's best cell over the
> 108-cell box (7.93); GEO's 40 sky draws (18 pulsars certify in ≥1 sky — the same 18 as SIREN's
> geometry union); SURFACE's `percert` (70 distinct pulsars certify *somewhere* in the box, so 46
> never certify at any modelled loudness).

### F10 — THE UPGRADED ARRAY · `F10_the_upgraded_array.png`

> Certification pays twice. Panel (a): the certified set does not just measure the source's
> distance, it *localises* it — the error box shrinks from ~1.6 square degrees (with none certified,
> "somewhere over there") to ~5×10⁻⁴ square degrees with five, a ~3,000× shrink that crosses the
> threshold where the box contains a single candidate host galaxy. Panel (b) is the second payment:
> **once a pulsar's distance is pinned to a fraction of a wavelength, it stays pinned.** It has
> become a permanent phase-connected baseline — available to the *next* source, and to mapping the
> background's anisotropy. The array phases itself up, one clock at a time.
>
> Source: `SIREN_results/siren_summary.npz`, arm B sky areas across the 9 (frequency × mass) cells.
> Zero-noise Fisher bound. Panel (b) is a schematic, not a measurement.

---

# The rest of the walls, and where this lands — F2.5, F2.7, F11–F14

F1–F6 pose the problem; F7–F10 explain the machine. These six close the argument: why you
cannot buy your way past the wall with telescopes (F2.5), what the one physical escape route
actually is (F2.7), where in the universe this puts the first measurable binary (F11), why you
do not get to pick your pulsars (F12), what happens if you ignore all of it (F13), and a check
that the whole capacity calculation agrees with theory (F14).
Script: `EXPLAINER_results/explainer_figures3.py`.

### F2.5 — THE SECOND WALL · `F2p5_the_second_wall.png`

> The obvious escape from F2's wall is to point a telescope at the source and let astronomy
> supply what the array cannot. It does not work. Conditioning on the source's sky position gets
> the posterior mass on the truth to 0.085; adding the orbital period takes it to 0.048; adding a
> 1%-accurate host-galaxy distance takes it to 0.043 — **against the 0.95 you would need, and the
> trend is flat and very slightly downward.** Each of those measurements is a serious observational
> campaign, and together they buy nothing, because the quantity that actually sets the pulsar
> clock's tick rate is the binary's **mass** — and no telescope supplies it. Panel (b) prices what
> would: the mass would have to be known **more than 20× better than today**, and even a 20×
> improvement only reaches 0.67, so this is a bound, not a target.
>
> **PROVENANCE (important):** these are the canonical values (STORY.md S4.2.8 / S4.2.10). Their npz
> — `b1_step2_table.npz`, `b1_breakeven_curve.npz` — are **not in this checkout**; they were produced
> on cronus and never pushed. The two referendum npz that *are* here (`b1_referendum_tier{A,C}.npz`)
> are the superseded 2-seed runs (f = 0.769 / 0.096) and are **not** plotted. Stated on the figure.

### F2.7 — THE SELF-CLOCKING CORNER · `F2p7_self_clocking_corner.png`

> If no telescope can supply the mass, the binary has to supply it itself — and an eccentric one
> can. A circular orbit emits a single tone that carries no timing information; an oval one emits a
> comb of harmonics whose spacing *is* its orbital clock. Panel (a) is the price of admission: the
> minimum eccentricity at which the binary measures its own mass 20× better, across mass and orbital
> frequency — roughly e ≳ 0.5–0.65, easiest for fast orbits. Panel (b) is the subtlety, and it is
> not in the textbook formula: **the same orbit self-clocks wildly differently depending on where its
> harmonics land.** At 10⁻⁹ Hz they fall below the band the array can hear and the clock is silent;
> half a decade up at 10⁻⁸·⁵ Hz they land in-band and it rings, beating the textbook prediction.
>
> Source: `reports/atlas_m2m4_summary.npz` (banked κ ≥ 20 contour) and `ATLAS_results/atlas_kappa_forb{0,1}.npz`.
> **Carried, not reconciled:** STORY flags a DISPUTED inconsistency — ATLAS's markdown labels a
> *different* column (0.58–0.70) as e*(>20×) while the npz holds the κ ≥ 20 values (0.50–0.65). I plot
> the npz and say so.

### F11 — WHERE IN THE UNIVERSE · `F11_where_in_the_universe.png`

> Putting the pieces together: how far away can the first measurable binary actually be? On the
> population we think we have, a 10⁹·⁵ M☉ binary at the fastest orbit must sit within **~9 Mpc** —
> *inside* the Virgo cluster, which is a small volume to bet a discovery on. Add two more loud
> sources to the population and the same binary works out to **~53 Mpc**, roughly a 180× larger
> volume. The red stars mark ATLAS's self-clocking corners with the eccentricity each demands.
> **Marked PROVISIONAL:** the onset floors these curves rest on are pending the RECUT re-fit, and
> the faint end of that surface is exactly where the floors are known to be permissive.
>
> Source: SURFACE's own amplitude relation, evaluated at its banked onsets (h* = 10⁻¹²·⁵ census,
> 10⁻¹³·²⁵ at 5-loud); reproduces SURFACE §10's 9.4 / 53.1 Mpc exactly.

### F12 — THE GEOMETRY LOTTERY · `F12_the_geometry_lottery.png`

> You do not get to choose which pulsars pay off — the source does, by where it happens to sit.
> Across 40 random source skies, **only 18 of the 116 pulsars ever certify at all**, and not the
> same ones twice: J1909−3744 pays off in 95% of skies, J0437−4715 in 80%, and the tail is a
> scatter of one-off winners. The other 98 never certify in any sky. Panel (b) shows the mechanism —
> a cliff. Once the source sits on the *far* side of the sky from a pulsar, that pulsar's chance of
> certifying collapses from ~30% to a few percent. Even the best pulsar in the array fails in 2 of
> 40 skies, and in both, the source had landed on the wrong side of it.
>
> Source: `reports/geo_summary.npz` (40 isotropic sky draws, census population, zero noise).

### F13 — THE WARNING · `F13_the_warning.png`

> This is what happens if you ignore all of the above and analyse a loud source with today's pulsar
> distances. The error does not shrink — the answer sits a stubborn 4–6° from the truth however loud
> the source gets (a). But the error box you *claim* shrinks by 160× (b). The two together are the
> failure: the fraction of the time the truth is actually inside your 90% box collapses **90% → 50%
> → 0%** as the source gets louder (c). **Wrong distances do not blur the answer, they move it — and
> a louder source makes it worse, because a louder source makes you more confident about a location
> that is wrong.** With parsec-class distances (green) it stays honest at every loudness. That is the
> case for a precision distance campaign stated as *prevention*, not enhancement.
>
> Source: `RING_results/npz/cell_fgw-8_A_snr{5,10,20}.npz`, scenario A, 10 realisations per cell;
> 'feather' = today's real distance priors, 'tier3' = near-exact distances.

### F14 — THE THEORY CHECK · `F14_the_theory_check.png`

> A sanity check that the whole capacity picture is right. How many sources a pulsar array can tell
> apart before they blur together was measured here at **0.399 · N^1.03** — linear in the number of
> pulsars, with no saturation out to 200. Boyle & Pen derived the theoretical bound in 2012 as
> **2N/7**. Ours is **2N/5**, and the ratio is 7/5 = 1.400, matched to under 1%. The two are not in
> conflict: Boyle & Pen count 7 unknowns per source, and our measurement *hands* the analysis the
> source's sky position, removing 2 of them. **The gap is exactly the price of not knowing where the
> sources are — worth 40% of the array's source-separating capacity.**
>
> Source: `CW_transition/stagec_p2c_results.npz`; the 2/5-vs-2/7 identification is STORY S2.1.3.

---

## What travels with these figures

1. **Simulated throughout.** Mock array, injected sources, stated observing-extension convention.
2. **F6 is conditional.** It prices *what certification buys*, not whether you can achieve it —
   and F2 is the standing negative on achieving it cold.
3. **Certified ≠ correct.** Every count in F4 and F5 carries a small rate of *wrong* certifications
   (a pulsar locked onto the wrong tooth), which no layer of the criterion currently removes.
4. **The two figures that point forward are F4 and F5**, and they point the same way: certification
   is a property of the *population*, not of any single source.
5. **F7, F9 and F10 are Fisher bounds on zero-noise data** — best cases. A real noisy analysis can
   only be wider.
6. **F9's ladder uses banked counts only.** The brief's "readable sub-array ≈ 30" is not a number
   any bank contains; the two quantities that do exist are **18** (pulsars certifying in ≥1 of GEO's
   40 sky draws — identical to SIREN's geometry-permitted union) and **70** (pulsars certifying
   somewhere in SURFACE's 108-cell box). Both are plotted; neither is 30.
7. **F2.5 is plotted from STORY, not from npz** — its two source banks are not in this checkout
   (see its caption). This is the only figure in the set not built from a banked array.
8. **F11 is PROVISIONAL** and watermarked as such: its onset floors are pending the RECUT re-fit.
9. **F2.7 carries a DISPUTED flag from STORY** — ATLAS's npz and markdown disagree on which
   eccentricity threshold the M4 payoff numbers are quoted at. Not reconciled here.
10. **F5's faintest column is estimator-limited.** Its floor is fitted where 60% of noise-only
   realisations produce no offender at all, which ANCHOR measured as biased permissive (its
   convention: the Gumbel floor is valid only below a ~20% zero-fraction). The count there is an
   upper bound. Flagging this *on* F5 is offered and not yet applied.
