# JOB 52 status (refreshed 2026-11-12 12:12)

HEAD 74cdf93 (MM_playground, in sync with origin). wheel_hasasia_examples/ IS committed.
Job 52 files uncommitted, staged (run_*.pkl and parts/ not staged).

- Parts 1+2: DONE -> JOB52_MAIN.md, job52_fig1_fraction_vs_N.png
  census N=7..12 known: 4.2-5.3% of pulsars right (prior 0.7%); flat with timing
  Arm A (unknown Earth terms removed) doubles everything; B, C change K>=60 by <= 0.02
  Arm D every candidate, census: 0.540 all / 0.472 K>=60 / 0.290 K>=600 (x1 timing)
  Arm D falls at zero timing noise: bins above 33 nHz treated as independent (diag_bins.log)
  P1 half wrong, P7 wrong, P2-P6 hold
- Parts 3-6: DONE -> JOB52_REST.md, job52_fig2_phase_error.png, job52_fig3_information.png
  P8 mixed, P9 holds, P10 holds (x2.7), P11 wrong (C - ln K not universal)
- Part 7: STOPPED. Machine draws pulsar-source angles independently (job50_common.py:181);
  no real sky matches, so Fisher code cannot share skies without major work.
- Part 8: DONE. Full rerun 2179 s (cap lifted by Matt). All 103 arrays match June npz,
  worst curve 5.8e-7 (one knee-diagnostic ratio), main curves <= 1.7e-10. June file untouched.

Cost: ranking 7 min, production 54 min (22 workers), Part 4/5 runs, analysis 2 min, Fisher 30 + 36 min.
All work done. Nothing committed.
