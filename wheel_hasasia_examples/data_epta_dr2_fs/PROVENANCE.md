# EPTA DR2new free-spectrum posterior samples (compact)

Source: Zenodo record 8091568, "The second data release from the European Pulsar
Timing Array III. Search for gravitational wave signals" — chains.zip →
chains/DR2new/{hd_fs,crn_fs}.tar.gz (PTMCMC chain_1.txt + pars.txt). Downloaded
2026-08-13. DR2new span 10.33 yr; 20 log10_rho bins at f_i = i/T.

These .npz files hold ONLY the gw log10_rho columns (float32, 25% burn-in
removed) — the full joint posterior over bins, so inter-bin covariance IS
available here (unlike the NG15 KDE product, which publishes marginals only).
