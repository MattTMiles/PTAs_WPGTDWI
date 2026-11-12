# JOB 49 status (refreshed 2026-11-11, end of session)

DONE
- Part A/B/C -> JOB49_AUDIT.md. H1 HOLDS (Gaussian draw, var (1-s)/s; job43c:749, job47c:124).
  Ledger finding: s counts each other source once; machine's own f_E split implies 2 terms ->
  incumbent small-angle variance is already correct; error is distribution SHAPE only.
- All booked numbers reproduce exactly (36.27@0.9009; 42.35/40.99/31.57; no crossing 8.3/10.6;
  47c grid; 47d grid with CSV-rounded seeds).
- Part D/E/F -> JOB49_PHYSICAL.md. P1 holds, P2 partly fails, P3 fails for D0 arm (holds for
  one-term arm only), P4 HOLDS (K>=60 <=13%, K>=600 <=1.6%), P5 holds.
- Part F: pile-up 232/44 unchanged.

NOT DONE
- Part G (real-likelihood cross-check): not started; needs GPU campaign + sqrt2 fix in
  CW_transition/stagec_p2b.py:86 first.

AT RISK / FLAGGED
- book() gate edge artifact (zero-padded convolve): booked 58.6% gates are artifacts; edge-safe
  headline crossing gives 39.69 @ 0.876 out-of-sample (< 0.90).
- Run order deviated: Part D launched before JOB49_AUDIT.md was written (stated in the audit).
- rec_*.pkl (120 MB) not staged; regenerate with job49_D2D3.py (~4 min).
