# JOB 51 status (refreshed 2026-11-11 19:55)

ALL PARTS DONE. Nothing committed; reports, code, logs, analysis.json staged (run_*.pkl, parts/ not staged).

- Part 0 -> JOB51_PROVENANCE.md: only version on this machine; wheel_hasasia_examples/ never committed.
  Missing: fig_phase_jitter_from_s.png, any trace of Job 44, a Job 45 report.
- Part 1 -> JOB51_MAIN.md: observed f, redshifted-mass strain, notebook FP_ALL uses source-frame mass
  (inconsistent). Both conventions run. arr = timing noise only.
  Headline every candidate: 17.5 @ 0.881 (redshifted) vs 36.9 @ 0.890 (notebook masses).
- Part 2: joint oracle fit K>=60 0.238 (x1), 0.270 (x10); Job 50 stack 0.236. R1 WRONG.
  Same stack on shared bins: 0.020. Limit = GW power outside candidate list (56% of band).
- Part 3: arr_bg booked-12 median 3.9 at any timing (need 6931 for 1 rad). 7 of 1161 candidates/sky
  have arr_bg >= 3. Full-template decode = prior at every timing. R2, R3 HOLD.
- Part 4: 0.5 sources/sky with >=4 pulsars at share>0.5; 3.3 at >0.2; census K<1: none.

Not run: Parts 2-4 with source-frame masses; census joint fit on 20 skies only.
Files: job51_common.py (hooks), job51_run.py, job51_analysis.py, job51_test.py, job50_frames.py;
       job51_logs/analysis.log, job51_data/analysis.json; job50_logs/analysis_src.log.
