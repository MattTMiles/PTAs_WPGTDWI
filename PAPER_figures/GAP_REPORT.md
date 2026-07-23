# GAP_REPORT — claims in STORY / the brief that could NOT be sourced to a bank on this machine

**Agent EASEL-5 · ACCRE · 2026-07-23.** One page. Everything below was verified by repo-wide
search or direct npz readback before being declared missing; nothing was improvised around a
gap — each gap is either displayed on its figure or the panel is reserved as PENDING.

## 1. P2(c) — the tier referendum and break-even curve npz are cronus-only
`b1_step2_table.npz` and `b1_breakeven_curve.npz` do not exist anywhere in this checkout
(repo-wide `find` empty; only the generator .py files exist). The on-disk
`b1_referendum_tier{A,C}.npz` are the SUPERSEDED 2-seed runs (f = 0.769/0.096) and are not
plotted. P2(c) therefore plots the canonical STORY values (S4.2.8 frozen 4-seed primary per
adjudicated D-2: 0.0847/0.0481/0.0323, with the 5-seed 0.0431 as a trail marker; S4.2.10
curve 0.032/0.107/0.289/0.673), with the provenance stated on the figure. **Pre-publication
flag (STORY D-2): re-run the frozen protocol at a stated seed count and adopt one number.**

## 2. P13 — the honest-posterior figure is unbuildable from banks anywhere in the repo
SAMPLER_dev.md §3: the posterior gates (g1/g2/g3/SBC) "did not run to completion on cronus"
(~7 GPU-hr/chain; host-RAM OOM). No NUTS chain, fringe-marginalised distance posterior, or SBC
bank exists on this machine (searched: `*chain*`, `*post*`, `*sbc*` npz — none). Panel (c)'s
"σ_MLE vs σ_sampler, 4.5–17×" does not exist in the SAMPLER record at all; the repo's only
"4–17×" is a different quantity (ATLAS sky-area shrink per SNR doubling, STORY:1399). P13 is
delivered as a reserved layout with per-panel missing-artifact statements. This is not
cronus-vs-ACCRE — the artifacts were never produced.

## 3. P15 — partially banked; bank contradicts report prose
- Banked coverage realisations: **N = 3** (scenario B: c = 0.053, 0.753; scenario C:
  c = 0.943), not the brief's "N = 5". The SAMPLER_dev.md g4 prose table quotes C = "0.903,
  0.813" — neither value is in `s4_scenC.npz` (which holds one realisation, 0.943). The npz is
  plotted; the discrepancy is stated here and N on the figure.
- The saddle spectrum ("3/8 negative eigenvalues", eig range) was read off `sampler_s4b.log`,
  which is absent (CW_transition/logs/ has no sampler logs) and was never banked
  (`sampler_s4.py` saves only `cond`). Only cond = 3.198e8 (C) and 4.9e7/1.3e9 (B) are banked
  and drawn; the spectrum is carried as provenance text on the panel.

## 4. P3(c) — ANCHOR's realised per-rung noise rms is report-prose, not npz
The 1.993 → 2.634 μs (×1.32) ladder rms values exist only in the tracked report table
(`ANCHOR_realdata_null.md` §2); no rms key exists in `reports/anchor_ladder.npz` or
`ANCHOR_results/anchor_floors.npz`. The floor shifts (the panel's data) are fully banked
(`d_emp`, `mw_p`); the rms row is annotated with an asterisk on the figure.

## 5. P14(b) — the MK9 corner is a notebook-only artifact
The MK9 delta-on-truth demonstration lives in
`CW_node_sampling/data_likelihood_sandbox_MK9_executed.ipynb` (parameters quoted via
SAMPLER_dev.md §8.2); there is no npz bank. It is cited as an annotation on a fully-banked
panel (archaeology_snap.npz), not plotted as data. The VLBI/T30 variant of the q-vs-purity
ladder (q 0.31→1.00, purity 0.78→0.37) is text-only in SAMPLER_dev.md; the banked
`archaeology_snap.npz` is the lit/T30 arm and is what P14(a) draws.

## 6. P11(d) — GLACIER's drain curve does not exist yet (by design)
The Stage-1/2b fan is queued (two jobs observed PENDING on dependency; nothing touched). Only
the zero-noise Stage-0 precursor is banked (`reports/glacier_g2_population.npz`, whose
`drain_Abg` is a zero placeholder). The half-panel is reserved and watermarked PENDING, per the
brief.

## 7. Small brief-vs-bank discrepancies (bank wins, per ARTIFACT READBACK)
- **P9(c) ladder**: brief says "7.83"; the corrected bank's best cell is **7.93**
  (`surface_analysis_recut.npz` max corr), matching STORY S6.1.3. 7.93 is drawn.
- **P7(c) motion specificity**: EMBER report rounds to 0.72; the banked contingency table
  (`cont_motion` = [5,27,0,80]) gives 80/107 = **0.75**. Derived-from-bank values are drawn.
- **P12 sampler-coverage cell**: brief says "TESTED at N=5"; the bank holds N = 3 (see §3).
  The cell reads "N=3 banked".
- **P1**: the only banked likelihood-vs-distance curve is a legacy single-source scan with
  coarser teeth than the census (stated on the panel; panel (c) carries the claim).
- **P2(a) blind-float band** (0.05–0.21) and pull-in (10⁻⁴): banked Track-B close-out numbers
  carried from the spec (trackB_L0/L2b/L2c npz exist but the band values are quoted as the
  spec's close-out constants, as EASEL-4 did).

## 8. Not gaps, but worth restating
- The 116-psr feather "residuals" are the injected CW+CURN mock (ANCHOR §0); nothing plotted
  here treats them as data.
- All counts/floors are from the `_recut` banks (criterion-v2.2); no superseded number
  (9.01-nat floor, single-member e = 0.3 switch, retired ratio-gain statistic, IGNITE h*) is
  drawn anywhere in the set. The retired objects appear only as explicitly-labelled retired
  bars (P3a old bar) or refuted points (P9c's 0.70).
