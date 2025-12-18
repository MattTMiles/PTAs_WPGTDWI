# Q1 scaffolding: MAP-based ring PTA experiments

Files:
- `sim.py`: load pulsars from `data_products`, build ring geometries, assemble enterprise/discovery likelihood pieces, and simulate residuals.
- `optimize.py`: small JAX MAP helpers (momentum ascent + optional LBFGS refinement).
- `run_q1_map.py`: minimal end-to-end example for one scenario/tier.

Quick start:
```bash
cd MAIN_PROJECT_QUESTIONS/Q1
python run_q1_map.py --scenario B --tier-scale 1.0 --tier-name real --opt-mode svi
```
Flags:
- `--scenario` A (WN), B (WN+RN), C (WN+RN+GWB)
- `--tier-scale` multiply pulsar distance σ (e.g., 3.0 loosens, 0.6 tightens, 0.0 ~ exact)
- `--tier-name` label for bookkeeping
- `--n-pulsars` size of ring (default 30)
- `--opt-mode` choose `svi` (NumPyro AutoDelta, recommended) or `map` (simple gradient ascent)
- `--steps`, `--lr` tweak the `map` optimiser (ignored for `svi`)
- `--seed` to reproduce CW/noise draws

How it works:
1. Draw CW(s) and noise params, set ring geometry (radius expanded if sources are separated), and force TOA uncertainties to 1e-6 before simulation to avoid outliers.
2. Build enterprise PTA, simulate a dataset, and mirror it into a discovery likelihood (pulsar-term phase-connected by default).
3. Run JAX MAP ascent over all free parameters; prints a small summary.

Next steps (intended):
- Hook up SNR targeting (scale `log10_h` to hit desired SNRs).
- Add NumPyro AutoDelta SVI option (see `/home/mattm/soft/pulsar-map-noise-estimates`).
- Add sweeps over SNR and distance tiers with result logging/plotting for detection/localisation metrics.
