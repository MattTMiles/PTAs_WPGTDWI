# JOB 50 status (refreshed 2026-11-11, after MAIN)

DONE
- JOB50_A.md (+ addendum A6: booked 601-pt grid too coarse for wide priors).
- JOB50_MAIN.md: Parts B-E. tau=0 reproduces p2_mf exactly; 232/44 reproduced (notebook masses).
  12 arms x harness + census priors; sensitivity x1/x3/x10; re-book + amplitude table.
  Q1,Q2,Q4,Q5,Q6,Q7 hold; Q3 holds without retardation (0.54), fails with it (0.24).
- Headline with retardation: no 0.90 catalogue (top 12); 17.5 @ 0.881 (every candidate).

PART F
- Part F done (JOB50_EXTRA.md): F1 figure, F2 sqrt2 correction. F3 not run.

AT RISK
- p2r uses Job 41 (redshifted) masses; notebook 232/44 uses source-frame masses (gives 44 vs 20).
  Retardation-on numbers not yet rerun with notebook masses.
- "Every candidate" = incoherent stack with independent noise per candidate, not a joint fit.
- First production attempt (40 min) lost: no incremental saves. Fixed (per-chunk saves).
- run_*.pkl and parts/ not staged (regenerate with job50_run.py, ~37 min on 22 cores).
