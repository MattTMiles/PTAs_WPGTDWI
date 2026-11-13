# JOB 53 status (refreshed 2026-11-12)

HEAD 74cdf93 (MM_playground). Nothing committed; Job 52 + 53 files staged.

- Part 0: DONE -> JOB53_SKIES.md (job53_test.py, job53_skies.py, job53_diag_mc.py)
  natural switch in job50_common.collect50(sky=..., mass=...); Job 52 reproduces 79/79
  natural N=7..12: 2.7-4.2% right (constructed 4.2-5.3%); every candidate 25.8% (27.1%)
  new defect: chirp-mass lookup wrong for 5% of sources + all 3 raised sources
- Parts 1-3: DONE -> JOB53_MAIN.md (job53_main.py, job53_heatmap.py; Figures 1-4)
  natural, 200 populations: N=7 0.029, N=12 0.046, every candidate 0.260; flat with timing
  retardation costs x1.6-1.7, all from the known source's own pulsar term (V3 holds)
  heatmap N_s=4: 90% contour at 2.6 pc (beta 0.8), 0.34 pc (0.5); Job 49 panel reproduced exactly
  prior-width slice N=12: 0.96 / 0.72 / 0.39 / 0.17 at 0.3 / 1 / 3 / 10 pc
  V1 half wrong, V5 2 of 3, V7 half wrong; V2, V3, V4, V6 hold
- Part 4 span: DONE -> JOB53_REST.md (job53_rest.py, job53_diag_span.py)
  N=12 all pulsars: 0.022 (15 yr) / 0.032 (25 yr) / 0.046 (45 yr); V8 holds (x2.1)
  collector silently skips skies with < 2 usable booked sources: 43 of 200 at 15 yr, 3 at 25 yr
- Part 5 six worlds: DONE -> JOB53_REST.md (job53_sixworlds.py)
  notebook stores 1.000/0.965/0.515/0.000 at K 8.2, not the quoted 0.985/0.540 at K 10.1
  joint fit on complex data: 1.000/0.99/0.79/0.00; final zero = one source + prior, not a statistic

All three reports on disk and staged. Nothing committed. Run pkl files and parts/ not staged.
Cost: Part 0 runs 11 min; production 88 min (20 workers); span 13 min; heatmap 8 min; analysis < 5 min.
