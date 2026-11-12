# JOB 50 — Part F

**Five-line summary.**
1. The share-versus-frequency figure is redrawn with both terms of every source and with true shares. In 15–40 nHz the 20 largest pulsar-term shares have a median of 0.07 (old figure: 0.19 with pulsar terms only).
2. Only 6 pulsar terms in that band hold more than 10% of their bin; about 100 hold more than 1%.
3. The real-likelihood √2 is corrected in a copy: median pulsar-term SNR 0.47 → 0.67, small-angle phase error 2.13 → 1.51 rad.
4. With exact statistics the median pulsar keeps a coherence of 0.40, against 0.10 from the booked formula. 47 of 116 pulsars have SNR above 1 (booked: 30).
5. F3 (same skies through both codes) was not run.

## F1 — share versus frequency, both terms

New file: `PAPER_figures/job50_share_vs_frequency_two_term.png`. The old figure is untouched. Frozen sky (draw_population seed 3), pulsar 41 (L = 2105 pc). Each point is one pulsar term at its retarded frequency with its retarded amplitude. The vertical axis is its true share: its power over all gravitational-wave power in its 1/T bin for this pulsar (Earth and pulsar terms of every source, plus twice the sea). No readability line.

The two panels differ only in the chirp-mass convention (see JOB50_MAIN.md, B2).

| 15–40 nHz | notebook masses | Job 41 (redshifted) masses |
|---|---|---|
| pulsar terms in band | 4153 | 3801 |
| median true share | 0.00032 | 0.00027 |
| median of the 20 largest | 0.070 | 0.052 |
| terms with share > 0.1 | 6 | 6 |
| terms with share > 0.01 | 104 | 67 |

For comparison, the old figure's top-20 median was 0.188 (pulsar terms only in the denominator).

## F2 — real-likelihood √2

`CW_transition/stagec_p2b.py:57` takes A_CW as the RMS of the pulsar-term residual, and `:86` divides the per-quadrature amplitude error by it. The corrected lines are in `job50_partF.py`:

```python
amp       = np.sqrt(2.0) * A_CW      # sinusoid amplitude, not RMS
rho       = amp / sig_res            # standard pulsar-term SNR
sigma_phi = 1.0 / rho                # small-angle phase error
```

They are applied to the saved per-pulsar arrays in `stagec_p2b_results.npz` (A_CW and sig_res are stored, so the correction is exact and nothing booked is rerun or overwritten). Output: `job50_data/F2_stagec_p2b_corrected.npz`.

| Array median (116 pulsars) | Booked | Corrected |
|---|---|---|
| pulsar-term SNR | 0.470 | 0.665 |
| small-angle phase error | 2.128 rad | 1.505 rad |
| coherence ⟨cos δφ⟩ | 0.104 (e^{−σ²/2}) | **0.395** (exact) |
| pulsars with SNR > 1 | 30 | 47 |

## F3

Not run. It needs a GPU campaign through `CW_transition` and, to be meaningful now, a decision on the template term and the mass convention first.
