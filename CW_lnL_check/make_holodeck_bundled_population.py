"""Build notebook-ready CW population arrays from bundled Holodeck sample data.

This does not modify the active Python environment.  It only reads the
Holodeck package data file and writes a compact NPZ consumed by
lnL_distance_scans/01c_singleCW_evolution_heatmap.ipynb.

The bundled file has axes (mtot, mrat, redz, freq_bin) but does not store the
frequency-bin centers.  We therefore use a conservative PTA-band log grid from
1e-9 to 1e-7 Hz and store that assumption in the output metadata.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_HOLODECK_DATA = (
    HERE
    / ".venv-holodeck"
    / "lib"
    / "python3.11"
    / "site-packages"
    / "holodeck"
    / "data"
    / "mbhb-pops-continuous_casey-clyde_mingarelli_2021-02-17.npz"
)
DEFAULT_OUTPUT = HERE / "holodeck_bundled_population_sample.npz"

C_LIGHT = 299792458.0
G_NEWTON = 6.67430e-11
MSUN_KG = 1.98847e30
MSUN_GEOM = G_NEWTON * MSUN_KG / C_LIGHT**3
MPC_M = 3.0856775814913673e22
YR = 365.25 * 24 * 3600


def luminosity_distance_seconds(z: np.ndarray) -> np.ndarray:
    """Flat LCDM luminosity distance, returned as light-seconds."""
    z = np.maximum(z, 1.0e-5)
    h0 = 70.0
    omega_m = 0.3
    omega_l = 0.7
    z_grid = np.linspace(0.0, max(3.0, float(np.nanmax(z)) + 0.01), 10000)
    ez = np.sqrt(omega_m * (1.0 + z_grid) ** 3 + omega_l)
    integ = 1.0 / ez
    dz = np.diff(z_grid)
    dc = np.zeros_like(z_grid)
    dc[1:] = (C_LIGHT / 1000.0 / h0) * np.cumsum(
        0.5 * (integ[1:] + integ[:-1]) * dz
    )
    dl_mpc = (1.0 + z) * np.interp(z, z_grid, dc)
    return dl_mpc * MPC_M / C_LIGHT


def make_sample(
    holodeck_data: Path = DEFAULT_HOLODECK_DATA,
    output: Path = DEFAULT_OUTPUT,
    n_sample: int = 250_000,
    tspan_years: float = 16.0,
    seed: int = 12345,
) -> Path:
    data = np.load(holodeck_data)
    weights = np.asarray(data["pops"], dtype=float)
    mtot_grid = np.asarray(data["mtot"], dtype=float)
    mrat_grid = np.asarray(data["mrat"], dtype=float)
    redz_grid = np.asarray(data["redz"], dtype=float)

    if weights.ndim != 4:
        raise ValueError(f"Expected 4D pops grid, got shape {weights.shape}")
    if weights.shape[:3] != (mtot_grid.size, mrat_grid.size, redz_grid.size):
        raise ValueError(
            "Expected pops axes (mtot, mrat, redz, freq); got "
            f"{weights.shape[:3]} versus "
            f"{(mtot_grid.size, mrat_grid.size, redz_grid.size)}"
        )

    nf = weights.shape[-1]
    fgw_grid = np.geomspace(1e-9, 1e-7, nf)
    prob = weights.ravel()
    prob = prob / np.sum(prob)

    rng = np.random.default_rng(seed)
    flat_idx = rng.choice(prob.size, size=n_sample, replace=True, p=prob)
    im, iq, iz, iff = np.unravel_index(flat_idx, weights.shape)

    mtot = mtot_grid[im]
    q = mrat_grid[iq]
    z = redz_grid[iz]
    fgw = fgw_grid[iff]

    mc_src = mtot * q ** (3.0 / 5.0) / (1.0 + q) ** (6.0 / 5.0)
    mc_z = mc_src * (1.0 + z)
    mc_geom_z = MSUN_GEOM * mc_z
    dl_sec = luminosity_distance_seconds(z)
    tspan = tspan_years * YR

    fdot = (96.0 / 5.0) * np.pi ** (8.0 / 3.0) * mc_geom_z ** (5.0 / 3.0) * fgw ** (11.0 / 3.0)
    ncorr = fgw / (np.abs(fdot) * tspan)
    h0 = 4.0 * mc_geom_z ** (5.0 / 3.0) * (np.pi * fgw) ** (2.0 / 3.0) / dl_sec

    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output,
        z=z,
        mtot=mtot,
        q=q,
        mc_src=mc_src,
        mc_z=mc_z,
        fgw=fgw,
        h0=h0,
        ncorr=ncorr,
        source_file=str(holodeck_data),
        model_name="holodeck bundled Mingarelli-style continuous population",
        frequency_assumption="last axis mapped to log-spaced fgw centers from 1e-9 to 1e-7 Hz",
        tspan_years=tspan_years,
    )
    return output


if __name__ == "__main__":
    path = make_sample()
    print(path)
