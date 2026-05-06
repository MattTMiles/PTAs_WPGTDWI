#!/usr/bin/env python3
"""
Pulsar Distance Estimation Script
==================================
Computes distance estimates using multiple methods and outputs:
1. best_distances.txt - Most precise distance estimate per pulsar
2. distance_comparison.txt - All distance estimates for comparison

Methods implemented:
- Cornell parallax catalog (live fetch: hosting.astro.cornell.edu/research/parallax/)
- VLBI parallax (legacy hand-curated list)
- Timing parallax (from PX in .par files)
- Lutz-Kelker corrected parallax
- DM-based distance (YMW16 model)
- DM-based distance (NE2001 model)
- PBDOT-derived distance (from Shklovskii effect inversion)
- Kinematic distance (from proper motion + Galactic rotation)
"""

import numpy as np
import glob
import os
import re
import json
import time
import urllib.request
from html.parser import HTMLParser
from astropy.coordinates import SkyCoord, Galactocentric
from astropy import units as u
from astropy import constants as const
from dataclasses import dataclass
from typing import Optional, Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# Cornell Parallax Catalog — fetched at runtime, cached locally
# https://hosting.astro.cornell.edu/research/parallax/
# Format: {jname: [[px_mas, px_err_mas, technique, reference], ...]}
# =============================================================================

_CATALOG_URL = 'https://hosting.astro.cornell.edu/research/parallax/'
_CATALOG_CACHE = os.path.expanduser('~/.cache/pulsar_parallax_catalog.json')
_CATALOG_CACHE_AGE_DAYS = 30


class _CellParser(HTMLParser):
    """Extract text from table cells, converting <br> to newlines."""

    def __init__(self):
        super().__init__()
        self.rows = []
        self._in_table = self._in_tr = self._in_cell = False
        self._cur_row: list = []
        self._cur_cell: list = []

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self._in_table = True
        elif tag == 'tr' and self._in_table:
            self._in_tr = True
            self._cur_row = []
        elif tag in ('td', 'th') and self._in_tr:
            self._in_cell = True
            self._cur_cell = []
        elif tag == 'br' and self._in_cell:
            self._cur_cell.append('\n')

    def handle_endtag(self, tag):
        if tag == 'table':
            self._in_table = False
        elif tag == 'tr' and self._in_tr:
            self._in_tr = False
            if self._cur_row:
                self.rows.append(self._cur_row)
        elif tag in ('td', 'th') and self._in_cell:
            self._in_cell = False
            self._cur_row.append(''.join(self._cur_cell).strip())

    def handle_data(self, data):
        if self._in_cell:
            self._cur_cell.append(data)


def _parse_px_cell(text):
    """
    Parse a parallax cell into list of (px_mas, err_mas).
    Handles symmetric (3.02 ± 0.07), asymmetric (0.93+0.08 / -0.07),
    and the case where the lower bound wraps to a new line after <br>.
    """
    # Merge continuation lines: a bare negative number on its own line
    # is the lower error bound of the previous asymmetric measurement.
    lines = text.replace('−', '-').replace('–', '-').split('\n')
    merged = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if re.match(r'^-\d+\.?\d*([eE][+-]?\d+)?$', line) and merged:
            merged[-1] += line  # "0.93+0.08" + "-0.07" → "0.93+0.08-0.07"
        else:
            merged.append(line)

    results = []
    for line in merged:
        line = line.replace('±', '+-').strip()
        # Symmetric: val +- err
        m = re.match(r'^([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\+-\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)$', line)
        if m:
            results.append((float(m.group(1)), float(m.group(2))))
            continue
        # Asymmetric: val+upper-lower
        m = re.match(r'^([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\+(\d+\.?\d*(?:[eE][+-]?\d+)?)-(\d+\.?\d*(?:[eE][+-]?\d+)?)$', line)
        if m:
            px = float(m.group(1))
            err = (float(m.group(2)) + float(m.group(3))) / 2.0
            results.append((px, err))
            continue
        # Plain value
        m = re.match(r'^([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)$', line)
        if m:
            results.append((float(m.group(1)), 0.0))
    return results


def _norm_technique(s):
    u = s.upper()
    if 'VLBI' in u:
        return 'VLBI'
    if any(k in u for k in ('GAIA', 'OPTICAL', 'HST', 'CHANDRA', 'XMM')):
        return 'Optical'
    return 'Timing'


def _parse_catalog_html(html):
    parser = _CellParser()
    parser.feed(html)

    catalog = {}
    header_found = False

    for row in parser.rows:
        if len(row) < 5:
            continue
        if not header_found:
            if any('parallax' in c.lower() or 'jname' in c.lower() for c in row):
                header_found = True
            continue

        jname = row[0].replace('–', '-').replace('—', '-').strip()
        if not jname:
            continue

        px_values = _parse_px_cell(row[2])
        refs  = [r.strip() for r in row[3].split('\n') if r.strip()]
        techs = [t.strip() for t in row[4].split('\n') if t.strip()]

        entries = []
        for i, (px, err) in enumerate(px_values):
            ref  = refs[i]  if i < len(refs)  else (refs[-1]  if refs  else 'Unknown')
            tech = _norm_technique(techs[i] if i < len(techs) else (techs[-1] if techs else ''))
            entries.append([px, err, tech, ref])

        if entries:
            catalog.setdefault(jname, []).extend(entries)

    return catalog


def load_parallax_catalog(force_refresh=False):
    """
    Load Cornell pulsar parallax catalog. Fetches from the web on first call
    or when the local cache is older than 30 days. Falls back to a stale cache
    if the network is unavailable.
    """
    cache_fresh = False
    try:
        age = time.time() - os.path.getmtime(_CATALOG_CACHE)
        cache_fresh = age < _CATALOG_CACHE_AGE_DAYS * 86400
    except FileNotFoundError:
        pass

    if not force_refresh and cache_fresh:
        try:
            with open(_CATALOG_CACHE) as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError):
            pass

    try:
        print(f"Fetching parallax catalog from {_CATALOG_URL} ...")
        with urllib.request.urlopen(_CATALOG_URL, timeout=15) as resp:
            html = resp.read().decode('utf-8', errors='replace')
        catalog = _parse_catalog_html(html)
        if len(catalog) < 50:
            raise ValueError(f"Only {len(catalog)} pulsars parsed — page structure may have changed")
        os.makedirs(os.path.dirname(_CATALOG_CACHE), exist_ok=True)
        with open(_CATALOG_CACHE, 'w') as f:
            json.dump(catalog, f)
        print(f"Cached {len(catalog)} pulsars to {_CATALOG_CACHE}")
        return catalog
    except Exception as e:
        print(f"Warning: could not fetch catalog ({e}), trying cache ...")
        try:
            with open(_CATALOG_CACHE) as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            print("Warning: no parallax catalog available. Catalog_PX method disabled.")
            return {}


PARALLAX_CATALOG = load_parallax_catalog()



# =============================================================================
# Published VLBI Distances — legacy hand-curated list, superseded by catalog.
# Units: kpc (distances, not parallaxes). Used directly by get_vlbi_distance().
# =============================================================================
PUBLISHED_VLBI_DISTANCES = {
    # name: (dist_kpc, err_kpc, source)
    'J0030+0451': (0.325, 0.009, 'NICER'),      # NICER geometric
    'J0437-4715': (0.157, 0.001, 'VLBI'),       # VLBI parallax
    'J1909-3744': (1.14, 0.03, 'VLBI'),         # VLBI, Reardon et al. 2016
    'J0613-0200': (0.78, 0.06, 'VLBI'),         # Deller et al. 2019
    'J1012+5307': (0.83, 0.06, 'VLBI'),         # Deller et al. 2019
    'J1024-0719': (1.08, 0.08, 'VLBI'),         # Deller et al. 2019
    'J1600-3053': (3.0, 0.5, 'VLBI'),           # Deller et al. 2019
    'J1713+0747': (1.18, 0.04, 'VLBI'),         # Chatterjee et al. 2009
    'J1738+0333': (1.47, 0.10, 'VLBI'),         # Deller et al. 2019
    'J2145-0750': (0.50, 0.03, 'VLBI'),         # Deller et al. 2019
}


@dataclass
class PulsarData:
    """Container for pulsar parameters from .par file"""
    name: str
    raj: str
    decj: str
    ra_deg: float = 0.0
    dec_deg: float = 0.0
    gl: float = 0.0
    gb: float = 0.0
    pb: Optional[float] = None          # Orbital period (days)
    pbdot: Optional[float] = None       # Observed Pbdot
    pbdot_err: Optional[float] = None
    px: Optional[float] = None          # Parallax (mas)
    px_err: Optional[float] = None
    dm: Optional[float] = None          # Dispersion measure
    pmra: Optional[float] = None        # Proper motion RA (mas/yr)
    pmdec: Optional[float] = None       # Proper motion Dec (mas/yr)
    m1: Optional[float] = None          # Pulsar mass (Msun)
    m2: Optional[float] = None          # Companion mass (Msun)
    mtot: Optional[float] = None        # Total mass (Msun)
    ecc: Optional[float] = None         # Eccentricity
    used_xpbdot: bool = False           # Whether XPBDOT was used


@dataclass
class DistanceEstimate:
    """Container for a single distance estimate"""
    value_kpc: float
    error_kpc: float
    method: str
    notes: str = ""

    @property
    def relative_error(self) -> float:
        if self.value_kpc > 0:
            return self.error_kpc / self.value_kpc
        return float('inf')


def parse_par_file(filename: str) -> Optional[PulsarData]:
    """Parse a TEMPO2 .par file to extract pulsar parameters."""
    data = {}
    pbdot_data = {}
    xpbdot_data = {}

    with open(filename, 'r') as f:
        for line in f:
            parts = line.split()
            if not parts or parts[0].startswith('#'):
                continue

            param = parts[0]

            if param == 'PSRJ':
                data['name'] = parts[1]
            elif param == 'RAJ':
                data['raj'] = parts[1]
            elif param == 'DECJ':
                data['decj'] = parts[1]
            elif param == 'PB':
                data['pb'] = float(parts[1])
            elif param == 'ECC':
                data['ecc'] = float(parts[1])
            elif param == 'M2':
                data['m2'] = float(parts[1])
            elif param == 'MTOT':
                data['mtot'] = float(parts[1])
            elif param == 'PX':
                data['px'] = float(parts[1])
                if len(parts) > 3 and parts[2] == '1':
                    data['px_err'] = float(parts[3])
            elif param == 'DM':
                data['dm'] = float(parts[1])
            elif param == 'PMRA':
                data['pmra'] = float(parts[1])
            elif param == 'PMDEC':
                data['pmdec'] = float(parts[1])
            elif param == 'PBDOT':
                pbdot_data['value'] = float(parts[1])
                if len(parts) > 3 and parts[2] == '1':
                    pbdot_data['err'] = float(parts[3])
            elif param == 'XPBDOT':
                xpbdot_data['value'] = float(parts[1])
                if len(parts) > 3 and parts[2] == '1':
                    xpbdot_data['err'] = float(parts[3])

    # Prefer XPBDOT over PBDOT
    if xpbdot_data.get('value') is not None:
        data['pbdot'] = xpbdot_data['value']
        data['pbdot_err'] = xpbdot_data.get('err')
        data['used_xpbdot'] = True
    elif pbdot_data.get('value') is not None:
        data['pbdot'] = pbdot_data['value']
        data['pbdot_err'] = pbdot_data.get('err')
        data['used_xpbdot'] = False

    # Estimate m1 if we have m2 and mtot
    if 'm2' in data and 'mtot' in data:
        data['m1'] = data['mtot'] - data['m2']

    # Check minimum required data
    if 'name' not in data or 'raj' not in data or 'decj' not in data:
        return None

    # Create PulsarData object
    psr = PulsarData(
        name=data['name'],
        raj=data['raj'],
        decj=data['decj'],
        pb=data.get('pb'),
        pbdot=data.get('pbdot'),
        pbdot_err=data.get('pbdot_err'),
        px=data.get('px'),
        px_err=data.get('px_err'),
        dm=data.get('dm'),
        pmra=data.get('pmra'),
        pmdec=data.get('pmdec'),
        m1=data.get('m1'),
        m2=data.get('m2'),
        mtot=data.get('mtot'),
        ecc=data.get('ecc'),
        used_xpbdot=data.get('used_xpbdot', False)
    )

    # Compute coordinates
    try:
        coord = SkyCoord(psr.raj, psr.decj, unit=(u.hourangle, u.deg))
        psr.ra_deg = coord.ra.deg
        psr.dec_deg = coord.dec.deg
        psr.gl = coord.galactic.l.deg
        psr.gb = coord.galactic.b.deg
    except:
        pass

    return psr


# =============================================================================
# Distance Estimation Methods
# =============================================================================

def get_catalog_parallax_distance(psr: PulsarData) -> Optional[DistanceEstimate]:
    """
    Distance from Cornell parallax catalog (March 2025).
    Selects highest-SNR measurement, preferring VLBI > Optical > Timing.
    Requires px > 0 and SNR >= 2.
    """
    entries = PARALLAX_CATALOG.get(psr.name)
    if not entries:
        return None

    tech_priority = {'VLBI': 0, 'Optical': 1, 'Timing': 2}

    valid = []
    for px, err, tech, ref in entries:
        if px > 0 and err > 0 and px / err >= 2.0:
            valid.append((px, err, tech, ref, px / err))

    if not valid:
        return None

    valid.sort(key=lambda x: (tech_priority.get(x[2], 3), -x[4]))
    px, err, tech, ref, _ = valid[0]

    dist_kpc = 1.0 / px          # d[kpc] = 1 / π[mas]  (exact by definition of parsec)
    err_kpc = dist_kpc * (err / px)

    return DistanceEstimate(dist_kpc, err_kpc, 'Catalog_PX',
                           f'{ref} ({tech}), PX={px:.3f}±{err:.3f} mas')


def get_vlbi_distance(psr: PulsarData) -> Optional[DistanceEstimate]:
    """Check for published VLBI distance (legacy list)."""
    if psr.name in PUBLISHED_VLBI_DISTANCES:
        dist, err, source = PUBLISHED_VLBI_DISTANCES[psr.name]
        return DistanceEstimate(dist, err, 'VLBI', f'Published {source}')
    return None


def get_timing_parallax_distance(psr: PulsarData) -> Optional[DistanceEstimate]:
    """Distance from timing parallax (PX in .par file)."""
    if psr.px is None or psr.px <= 0:
        return None

    # Require positive parallax with reasonable S/N
    px_err = psr.px_err if psr.px_err else psr.px * 0.5
    if psr.px < 0.1 or px_err / psr.px > 1.0:
        return None

    dist_kpc = 1.0 / psr.px  # mas to kpc
    err_kpc = dist_kpc * (px_err / psr.px)

    return DistanceEstimate(dist_kpc, err_kpc, 'Timing_PX',
                           f'PX={psr.px:.3f}±{px_err:.3f} mas')


def get_lutz_kelker_distance(psr: PulsarData) -> Optional[DistanceEstimate]:
    """
    Lutz-Kelker corrected parallax distance.

    Applies statistical correction for parallax bias when
    selecting objects based on measured parallax.

    Reference: Lutz & Kelker (1973), PASP 85, 573
    """
    if psr.px is None or psr.px <= 0 or psr.px_err is None:
        return None

    # Require reasonable S/N for LK correction to be meaningful
    snr = psr.px / psr.px_err
    if snr < 2.0:
        return None

    # Fractional parallax error
    f = psr.px_err / psr.px

    # LK correction factor (approximate for uniform space density)
    # For small f: correction ~ 1 + 1.5*f^2
    # More accurate: use the full LK bias formula
    if f < 0.175:
        # Low error regime - apply standard LK correction
        lk_factor = 1.0 + 1.5 * f**2 + 3.75 * f**4
    else:
        # Higher error - use modified correction to avoid overcorrection
        lk_factor = 1.0 + 1.5 * f**2

    # Corrected parallax
    px_corrected = psr.px / lk_factor

    # Distance and error
    dist_kpc = 1.0 / px_corrected
    # Error increases due to correction uncertainty
    err_kpc = dist_kpc * f * np.sqrt(1 + 2 * f**2)

    return DistanceEstimate(dist_kpc, err_kpc, 'LK_Parallax',
                           f'LK correction factor={lk_factor:.3f}')


def get_dm_distance_ymw16(psr: PulsarData) -> Optional[DistanceEstimate]:
    """
    Distance from DM using simplified YMW16 electron density model.

    Reference: Yao, Manchester & Wang (2017)
    """
    if psr.dm is None or psr.dm <= 0:
        return None

    gl, gb = psr.gl, psr.gb

    # Simplified YMW16-like model
    # Mean electron density depends on Galactic latitude
    abs_b = abs(gb)

    if abs_b < 5:
        # Near Galactic plane - higher density, more uncertain
        ne_mean = 0.03  # cm^-3
        ne_scatter = 0.015
    elif abs_b < 20:
        # Intermediate latitude
        ne_mean = 0.025
        ne_scatter = 0.01
    else:
        # High latitude - lower density
        ne_mean = 0.02
        ne_scatter = 0.008

    # Correct for scale height
    scale_height = 1.0  # kpc
    z_factor = 1.0 / np.cos(np.radians(gb)) if abs_b < 85 else 10.0

    # Distance estimate
    dist_kpc = psr.dm / (ne_mean * 1000 * z_factor)  # DM in pc/cm^3

    # Uncertainty (dominated by ne scatter, typically 30-50%)
    rel_err = ne_scatter / ne_mean + 0.2  # Add 20% systematic
    err_kpc = dist_kpc * rel_err

    # Sanity bounds
    dist_kpc = max(0.1, min(dist_kpc, 30.0))

    return DistanceEstimate(dist_kpc, err_kpc, 'DM_YMW16',
                           f'DM={psr.dm:.2f}, ne={ne_mean:.3f}')


def get_dm_distance_ne2001(psr: PulsarData) -> Optional[DistanceEstimate]:
    """
    Distance from DM using simplified NE2001 electron density model.

    Reference: Cordes & Lazio (2002)
    """
    if psr.dm is None or psr.dm <= 0:
        return None

    gl, gb = psr.gl, psr.gb

    # NE2001 uses different scale heights and spiral arm contributions
    abs_b = abs(gb)

    # Thick disk component
    n1 = 0.033  # cm^-3 midplane density
    h1 = 0.97   # kpc scale height

    # Thin disk component
    n2 = 0.090
    h2 = 0.14

    # Effective mean density along line of sight
    sin_b = np.sin(np.radians(abs_b)) if abs_b > 0 else 0.01

    # Approximate integral through disk
    if abs_b < 5:
        ne_eff = n1 + n2 * 0.5  # Both components contribute
    elif abs_b < 30:
        ne_eff = n1 * np.exp(-abs_b/30) + n2 * 0.1
    else:
        ne_eff = n1 * np.exp(-1) * 0.5  # Far from plane

    # Distance
    path_factor = 1.0 / max(sin_b, 0.1) if abs_b > 5 else 1.0
    dist_kpc = psr.dm / (ne_eff * 1000) / min(path_factor, 3.0)

    # NE2001 typically has ~20-40% uncertainty
    err_kpc = dist_kpc * 0.35

    # Sanity bounds
    dist_kpc = max(0.1, min(dist_kpc, 30.0))

    return DistanceEstimate(dist_kpc, err_kpc, 'DM_NE2001',
                           f'DM={psr.dm:.2f}, ne_eff={ne_eff:.3f}')


def compute_pbdot_gr(pb: float, ecc: float, m1: float, m2: float) -> float:
    """Compute GR orbital decay contribution to Pbdot."""
    pb_sec = pb * 86400.0
    m1_kg = m1 * const.M_sun.value
    m2_kg = m2 * const.M_sun.value

    # Chirp mass
    m_chirp = (m1_kg * m2_kg)**(3/5) / (m1_kg + m2_kg)**(1/5)

    # Eccentricity function
    f_e = (1 + 73/24 * ecc**2 + 37/96 * ecc**4) / (1 - ecc**2)**(7/2)

    # GR formula (Peters & Mathews 1963)
    pbdot_gr = -(192 * np.pi / 5) * (pb_sec / (2*np.pi))**(-5/3) * \
               (const.G.value * m_chirp / const.c.value**3)**(5/3) * f_e

    return pbdot_gr


def get_pbdot_distance(psr: PulsarData) -> Optional[DistanceEstimate]:
    """
    Distance derived from PBDOT via Shklovskii effect inversion.

    Shklovskii effect: Pbdot_shk = Pb * d * μ² / c

    Rearranging: d = Pbdot_shk * c / (Pb * μ²)

    Where Pbdot_shk = Pbdot_obs - Pbdot_GR (for binaries with known masses)
    """
    # Need PBDOT, proper motion, and orbital period
    if psr.pbdot is None or psr.pb is None:
        return None
    if psr.pmra is None or psr.pmdec is None:
        return None

    # Total proper motion
    mu_total = np.sqrt(psr.pmra**2 + psr.pmdec**2)  # mas/yr
    if mu_total < 1.0:  # Need significant proper motion
        return None

    # Convert proper motion to rad/s
    mu_rad_s = mu_total * (1e-3 / 3600) * (np.pi / 180) / (365.25 * 86400)

    # Get Pbdot_shk by subtracting GR contribution if masses known
    pbdot_shk = psr.pbdot

    if psr.m1 is not None and psr.m2 is not None and psr.ecc is not None:
        pbdot_gr = compute_pbdot_gr(psr.pb, psr.ecc, psr.m1, psr.m2)
        pbdot_shk = psr.pbdot - pbdot_gr
        notes = f'Pbdot_obs={psr.pbdot:.2e}, Pbdot_GR={pbdot_gr:.2e}'
    else:
        # Assume GR contribution is small for wide binaries
        notes = f'Pbdot_obs={psr.pbdot:.2e} (GR not subtracted)'

    # Shklovskii must be positive (apparent spindown)
    if pbdot_shk <= 0:
        return None

    # Distance from inverted Shklovskii formula
    pb_sec = psr.pb * 86400.0
    d_m = pbdot_shk * const.c.value / (pb_sec * mu_rad_s**2)
    dist_kpc = d_m / (1000 * const.pc.value)

    # Error estimation (propagate uncertainties)
    # Dominated by proper motion and Pbdot uncertainties
    rel_err_mu = 0.05  # Assume 5% PM uncertainty if not given
    rel_err_pbdot = abs(psr.pbdot_err / psr.pbdot) if psr.pbdot_err else 0.1

    rel_err = np.sqrt(rel_err_pbdot**2 + 4 * rel_err_mu**2)  # μ² so 2x
    err_kpc = dist_kpc * rel_err

    # Sanity check
    if dist_kpc < 0.01 or dist_kpc > 50:
        return None

    return DistanceEstimate(dist_kpc, err_kpc, 'PBDOT_Shk', notes)


def get_kinematic_distance(psr: PulsarData) -> Optional[DistanceEstimate]:
    """
    Kinematic distance from proper motion and Galactic rotation model.

    Uses the fact that objects at different distances have different
    expected proper motions due to differential Galactic rotation.
    """
    if psr.pmra is None or psr.pmdec is None:
        return None

    # Need Galactic coordinates
    gl, gb = psr.gl, psr.gb

    # Only reliable near Galactic plane
    if abs(gb) > 15:
        return None

    # Total proper motion
    mu_total = np.sqrt(psr.pmra**2 + psr.pmdec**2)  # mas/yr

    if mu_total < 2.0:  # Need measurable proper motion
        return None

    # Galactic rotation parameters
    R0 = 8.122  # kpc, distance to Galactic center
    V0 = 220.0  # km/s, circular velocity at Sun

    # Convert to Galactic proper motion components (approximate)
    coord = SkyCoord(ra=psr.ra_deg*u.deg, dec=psr.dec_deg*u.deg)

    # For objects in the disk, proper motion is dominated by:
    # 1. Solar motion (reflex of Sun's motion)
    # 2. Galactic rotation (differential rotation)
    # 3. Peculiar velocity of the pulsar

    # Solar motion contribution (dominant for nearby objects)
    # V_sun ~ 220 km/s toward l=90°
    # At distance d, this gives μ ~ V_sun / (4.74 * d)

    # Expected proper motion from solar motion alone
    l_rad = np.radians(gl)

    # Approximate: μ_l ~ V0 * sin(l) / (4.74 * d)
    # Solve for d: d ~ V0 * |sin(l)| / (4.74 * μ)

    sin_l = np.sin(l_rad)
    if abs(sin_l) < 0.3:  # Near l=0 or l=180, method unreliable
        return None

    # Convert proper motion to km/s/kpc
    # μ (mas/yr) = v (km/s) / (4.74047 * d (kpc))

    # Estimate distance assuming motion dominated by Galactic rotation
    # This is a rough approximation
    v_expected = V0 * abs(sin_l) * 0.5  # Factor 0.5 for average effect

    dist_kpc = v_expected / (4.74047 * mu_total)

    # Large uncertainty due to assumptions
    err_kpc = dist_kpc * 0.5  # 50% uncertainty typical

    # Sanity bounds
    if dist_kpc < 0.1 or dist_kpc > 20:
        return None

    return DistanceEstimate(dist_kpc, err_kpc, 'Kinematic',
                           f'μ={mu_total:.1f} mas/yr, l={gl:.1f}°')


# =============================================================================
# Main Processing
# =============================================================================

def get_all_distances(psr: PulsarData) -> List[DistanceEstimate]:
    """Get all available distance estimates for a pulsar."""
    estimates = []

    # Try each method
    methods = [
        get_catalog_parallax_distance,
        get_vlbi_distance,
        get_timing_parallax_distance,
        get_lutz_kelker_distance,
        get_dm_distance_ymw16,
        get_dm_distance_ne2001,
        get_pbdot_distance,
        get_kinematic_distance,
    ]

    for method in methods:
        try:
            est = method(psr)
            if est is not None:
                estimates.append(est)
        except Exception as e:
            pass

    return estimates


def select_best_distance(estimates: List[DistanceEstimate]) -> Optional[DistanceEstimate]:
    """Select the best distance estimate based on method reliability and precision."""
    if not estimates:
        return None

    # Priority ranking (lower = better)
    method_priority = {
        'Catalog_PX': 1,
        'VLBI': 2,
        'Timing_PX': 3,
        'LK_Parallax': 4,
        'PBDOT_Shk': 5,
        'DM_YMW16': 6,
        'DM_NE2001': 7,
        'Kinematic': 8,
    }

    # Sort by: (1) method priority, (2) relative error
    def sort_key(est):
        priority = method_priority.get(est.method, 10)
        return (priority, est.relative_error)

    estimates_sorted = sorted(estimates, key=sort_key)
    return estimates_sorted[0]


def write_best_distances(pulsars: List[Tuple[PulsarData, DistanceEstimate]],
                         filename: str):
    """Write the best distance estimates to a file."""
    with open(filename, 'w') as f:
        f.write("# Best Distance Estimates for Pulsars\n")
        f.write("# Generated by get_distance.py\n")
        f.write("#\n")
        f.write("# Methods (in priority order):\n")
        f.write("#   Catalog_PX - Cornell parallax catalog (VLBI/Optical/Timing, March 2025)\n")
        f.write("#   VLBI       - Legacy published VLBI parallax (geometric)\n")
        f.write("#   Timing_PX  - Timing parallax from .par file\n")
        f.write("#   LK_Parallax- Lutz-Kelker corrected parallax\n")
        f.write("#   PBDOT_Shk  - Distance from Shklovskii effect inversion\n")
        f.write("#   DM_YMW16   - DM-based distance (YMW16 model)\n")
        f.write("#   DM_NE2001  - DM-based distance (NE2001 model)\n")
        f.write("#   Kinematic  - Kinematic distance from proper motion\n")
        f.write("#\n")
        f.write(f"# {'Pulsar':<14} {'Dist(kpc)':>10} {'Err(kpc)':>10} {'Rel.Err':>8} {'Method':<12} Notes\n")
        f.write("#" + "-"*90 + "\n")

        for psr, est in sorted(pulsars, key=lambda x: x[0].name):
            f.write(f"  {psr.name:<14} {est.value_kpc:>10.4f} {est.error_kpc:>10.4f} "
                   f"{est.relative_error:>8.1%} {est.method:<12} {est.notes}\n")

    print(f"Wrote best distances to: {filename}")


def write_distance_comparison(pulsars: List[Tuple[PulsarData, List[DistanceEstimate]]],
                              filename: str):
    """Write comparison of all distance estimates to a file."""
    # Get all unique methods
    all_methods = ['Catalog_PX', 'VLBI', 'Timing_PX', 'LK_Parallax',
                   'PBDOT_Shk', 'DM_YMW16', 'DM_NE2001', 'Kinematic']

    with open(filename, 'w') as f:
        f.write("# Distance Estimate Comparison for Pulsars\n")
        f.write("# Generated by get_distance.py\n")
        f.write("#\n")
        f.write("# All values in kpc. Format: distance ± error\n")
        f.write("# '-' indicates method not applicable or failed\n")
        f.write("#\n")

        # Header
        header = f"{'Pulsar':<14}"
        for method in all_methods:
            header += f" {method:>18}"
        f.write(header + "\n")
        f.write("#" + "-"*145 + "\n")

        for psr, estimates in sorted(pulsars, key=lambda x: x[0].name):
            # Build dict for easy lookup
            est_dict = {e.method: e for e in estimates}

            line = f"{psr.name:<14}"
            for method in all_methods:
                if method in est_dict:
                    e = est_dict[method]
                    line += f" {e.value_kpc:>7.3f}±{e.error_kpc:<7.3f}"
                else:
                    line += f" {'-':>18}"
            f.write(line + "\n")

        # Summary statistics
        f.write("\n\n# Summary Statistics\n")
        f.write("#" + "="*60 + "\n")

        method_counts = {m: 0 for m in all_methods}
        for psr, estimates in pulsars:
            for e in estimates:
                if e.method in method_counts:
                    method_counts[e.method] += 1

        f.write(f"# Total pulsars: {len(pulsars)}\n")
        f.write("# Methods available per pulsar:\n")
        for method in all_methods:
            f.write(f"#   {method:<14}: {method_counts[method]:3d} pulsars "
                   f"({100*method_counts[method]/len(pulsars):5.1f}%)\n")

    print(f"Wrote distance comparison to: {filename}")


def main():
    """Main function to process all pulsars and generate output files."""
    # Find all .par files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    par_pattern = os.path.join(script_dir, 'partim', '*_tdb.par')
    par_files = sorted(glob.glob(par_pattern))

    if not par_files:
        # Try without _tdb suffix
        par_pattern = os.path.join(script_dir, 'partim', '*.par')
        par_files = sorted(glob.glob(par_pattern))

    print(f"Found {len(par_files)} .par files")

    # Parse all pulsars
    pulsars_data = []
    for pf in par_files:
        psr = parse_par_file(pf)
        if psr:
            pulsars_data.append(psr)

    print(f"Successfully parsed {len(pulsars_data)} pulsars")

    # Get all distance estimates
    pulsars_with_best = []
    pulsars_with_all = []

    for psr in pulsars_data:
        estimates = get_all_distances(psr)
        if estimates:
            best = select_best_distance(estimates)
            if best:
                pulsars_with_best.append((psr, best))
            pulsars_with_all.append((psr, estimates))

    print(f"Distance estimates available for {len(pulsars_with_best)} pulsars")

    # Write output files
    output_dir = script_dir
    write_best_distances(pulsars_with_best,
                        os.path.join(output_dir, 'best_distances.txt'))
    write_distance_comparison(pulsars_with_all,
                             os.path.join(output_dir, 'distance_comparison.txt'))

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    method_counts = {}
    for psr, est in pulsars_with_best:
        method_counts[est.method] = method_counts.get(est.method, 0) + 1

    print("\nBest method distribution:")
    for method, count in sorted(method_counts.items(), key=lambda x: -x[1]):
        print(f"  {method:<14}: {count:3d} pulsars ({100*count/len(pulsars_with_best):5.1f}%)")

    # Distance statistics
    distances = [est.value_kpc for _, est in pulsars_with_best]
    print(f"\nDistance range: {min(distances):.2f} - {max(distances):.2f} kpc")
    print(f"Median distance: {np.median(distances):.2f} kpc")
    print(f"Mean distance: {np.mean(distances):.2f} ± {np.std(distances):.2f} kpc")


if __name__ == '__main__':
    main()
