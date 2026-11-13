# JOB 54 status (refreshed 2026-11-12 18:45)

HEAD 74cdf93 (MM_playground). Nothing committed; Jobs 52-54 files staged. ALL PARTS DONE.

- Part 1: talk_data/ ready (7.0 MB; README.md; check.py passes) -> JOB54_DATA.md
  ladder N=12: 0.999 -> 0.744 (sea) -> 0.071 (other binaries) -> 0.041 (retardation) -> 0.042
  (timing) -> 0.023 (25 yr) -> 0.015 (15 yr); all 200 skies, none skipped
- DEFECT FOUND AND CORRECTED: unknown candidates' pulsar terms had no random phase
  (job53_common.py:119). Census results drop slightly (N=12 0.046 -> 0.041; all 0.260 -> 0.238);
  Job 53 prior-width slice was too optimistic (N=12 at 1 pc 0.72 -> 0.50; 3 pc 0.39 -> 0.29);
  certification at 1 pc: 3.9 per sky (was 45.7); at 3 pc: none. All talk files use corrected runs.
- Part 3: per-bin model with real shares = corrected skies within 0.007 at every prior width (W7 holds)
- Part 4: amplitude-only arm explains 82-103% of the retardation loss, frequency-only 48-55%;
  pulsar terms 47% of the competition; destination bins not louder (0.97) (W8: 2 of 4)
- Side finding: candidate cut (single-pulsar SNR >= 3) hides most background-detectable sources
  (12.6 / 30.5 per sky, not 4.6 / 10.1); N = 7, 12 results unchanged with them included
- Part 2: JOB54_LOCALISATION.md, Figures 1-3. Anchors pass (1.0e-5 and 330 deg^2). At SNR 10,
  3 phased pulsars shrink the 90% area x1.1, 10 by x1.5; collapse only at 40-80 phased.
  Fisher is wrong by 1e2-1e6 for few phased pulsars. Detection reach: all phased x1.43 (A), x1.39 (B)
  W1 pass, W2 fail, W3 half, W4 pass, W5 pass, W6 holds
- GPU: XLA_PYTHON_CLIENT_PREALLOCATE=false now set in job54_loc_core.py
