# NG15 free-spectrum KDE posteriors

Source: Zenodo record 21844115, "KDE Representations of the Gravitational Wave
Background Free Spectra Present in the NANOGrav 15-Year Dataset" (ceffyl data
product; la-forge core `model_crn_30f_ptswapfix`; FFTKDE, epanechnikov kernel).
Downloaded 2026-08-13. Underlying dataset: NANOGrav 15 yr (10.5281/zenodo.7967584).

- `hd/`     : 30-frequency HD-correlated free spectrum  (primary for V4u)
- `cp/`     : 30-frequency CURN (auto-only) free spectrum (robustness)
- `hdonly/` : HD process from the multi-process fs{hd+mp+dp+cp} run (robustness)

Each dir: freqs.npy (30 bins, f_i = i/16.03 yr), log10rhogrid.npy (10000 pts,
[-9,-4]), density.npy (1,30,10000) per-bin KDE marginals. NOTE: marginals only —
no inter-bin covariance is published in this product.
