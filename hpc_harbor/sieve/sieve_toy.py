"""SIEVE-C toy spine -- the deterministic half of the SIEVE battery (agent SIEVE-C, cronus).

HARD SCOPE (declared, travels with every number produced here): this box is
RECOMPUTED-UNSAFE for noisy realisations.  EVERYTHING in this module is
ZERO-NOISE / ASIMOV.  No `np.random` draw ever enters a datum: the data is the
noiseless injected model and the likelihood is the exact quadratic mismatch

    lnL(theta) = -0.5 * sum_p || r_p(theta) - r_p(truth) ||^2 / sigma_p^2

so lnL(truth) = 0 by construction and every dlnL quoted is a pure
model-mismatch cost, not a noise realisation.  Any test that needs a fresh
noise draw DEFERS to SIEVE-A on ACCRE.

THE TOY SPINE (brief): 12-20 pulsars from the top of the census by noise rank,
REAL sky positions and REAL distances; 2-4 sources; T = 15-20 yr; zero noise.

  positions : `pos_t` (row 0) out of the tracked feathers in data_products/ --
              the enterprise unit vector, so the sky geometry is the array's own.
  distances : CW_transition/best_distances.txt (the canonical prior file).
  noise rank: median TOA error out of reports/anchor_data_forensics.npz
              (ascending = best-timed first).  The measured per-pulsar SPANS are
              heterogeneous (6.5-20 yr); the toy imposes a UNIFORM T and cadence
              on every pulsar, which is a TOY CHOICE, declared, not a claim
              about the real array's sampling.

CW MODEL: circular, non-evolving (fdot = 0 unless `fdot` is passed), Earth term
+ pulsar term, standard PTA antenna pattern.  The pulsar-term retardation is

    dt_p = (L_p / c) * (1 - shat . phat)                      [PLANE WAVE]

and the curved-front (T1) replacement is the EXACT source-to-pulsar path

    dt_p = (1/c) * ( |D shat - L_p phat| + L_p - D )          [SPHERICAL]

whose leading term is the brief's Fresnel phase pi L_perp^2 / (lambda D).
"""

import os
import numpy as np

C_LIGHT = 299792458.0
PC = 3.0856775814913673e16          # m
KPC = 1e3 * PC
MPC = 1e6 * PC
YR = 365.25 * 86400.0

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FEATHER_DIR = os.path.join(ROOT, "data_products")
DIST_FILE = os.path.join(ROOT, "CW_transition", "best_distances.txt")
FORENSIC = os.path.join(ROOT, "reports", "anchor_data_forensics.npz")


# ------------------------------------------------------------------
# 1. the spine
# ------------------------------------------------------------------
def _read_distances():
    """canonical prior file -> {name: (dist_kpc, err_kpc, method)}"""
    out = {}
    with open(DIST_FILE) as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            tok = line.split()
            if len(tok) < 5:
                continue
            out[tok[0]] = (float(tok[1]), float(tok[2]), tok[4])
    return out


def _psr_pos(name):
    """real sky unit vector from the tracked feather (`pos_t`, row 0)."""
    import pyarrow.feather as pf
    t = pf.read_table(os.path.join(FEATHER_DIR, name + ".feather"),
                      columns=["pos_t_0", "pos_t_1", "pos_t_2"])
    v = np.array([t["pos_t_0"][0].as_py(),
                  t["pos_t_1"][0].as_py(),
                  t["pos_t_2"][0].as_py()], float)
    return v / np.linalg.norm(v)


def build_spine(npsr=16, T_yr=18.0, cadence_per_yr=20.0, verbose=True):
    """Top-`npsr` census pulsars by noise rank, with REAL positions/distances.

    Returns a dict; every array is float64 and deterministic (no RNG anywhere).
    """
    fo = np.load(FORENSIC, allow_pickle=True)
    names_all = [str(x) for x in fo["names"]]
    toaerr_all = np.asarray(fo["median_toaerr_s"], float)
    dists = _read_distances()

    # noise rank = median TOA error ascending, restricted to pulsars that have a
    # canonical distance prior AND a tracked feather.
    cand = [(toaerr_all[i], names_all[i]) for i in range(len(names_all))
            if names_all[i] in dists
            and os.path.exists(os.path.join(FEATHER_DIR, names_all[i] + ".feather"))]
    cand.sort()
    sel = cand[:npsr]

    names = [n for _, n in sel]
    sigma = np.array([e for e, _ in sel], float)              # s
    phat = np.array([_psr_pos(n) for n in names], float)      # (npsr, 3)
    L_kpc = np.array([dists[n][0] for n in names], float)
    Lerr_kpc = np.array([dists[n][1] for n in names], float)
    meth = [dists[n][2] for n in names]

    ntoa = int(round(T_yr * cadence_per_yr))
    toas = np.linspace(0.0, T_yr * YR, ntoa)                  # uniform toy cadence

    spine = dict(names=names, phat=phat, sigma=sigma, L_m=L_kpc * KPC,
                 Lerr_m=Lerr_kpc * KPC, L_kpc=L_kpc, Lerr_kpc=Lerr_kpc,
                 method=meth, toas=toas, T_yr=T_yr, npsr=len(names))
    if verbose:
        print(f"[spine] npsr={len(names)} T={T_yr}yr ntoa/psr={ntoa} "
              f"sigma {sigma.min():.2e}-{sigma.max():.2e}s "
              f"L {L_kpc.min():.3f}-{L_kpc.max():.3f}kpc")
    return spine


def sky_unit(costheta, phi):
    """source unit vector shat (Earth -> source) from (cos theta, phi)."""
    st = np.sqrt(max(0.0, 1.0 - costheta ** 2))
    return np.array([st * np.cos(phi), st * np.sin(phi), costheta], float)


# ------------------------------------------------------------------
# 2. antenna pattern + waveform
# ------------------------------------------------------------------
def antenna(shat, phat):
    """(F_plus, F_cross) for source direction shat, pulsar directions phat (n,3).

    Standard PTA pattern with the GW propagating along Omega = -shat.
    The polarisation basis (m, n) is built from shat; a rotation by psi is
    applied later in the waveform (F rotates as cos2psi/sin2psi).
    """
    z = np.array([0.0, 0.0, 1.0])
    m = np.cross(z, shat)
    nm = np.linalg.norm(m)
    if nm < 1e-12:                       # source at a pole -- pick any basis
        m = np.array([1.0, 0.0, 0.0])
        nm = 1.0
    m = m / nm
    n = np.cross(shat, m)
    Om = -shat
    denom = 1.0 + phat @ Om              # (npsr,)
    mp = phat @ m
    np_ = phat @ n
    Fp = 0.5 * (mp ** 2 - np_ ** 2) / denom
    Fx = (mp * np_) / denom
    return Fp, Fx


def _retard_plane(shat, phat, L_m):
    """plane-wave pulsar-term retardation dt_p (s)."""
    return (L_m / C_LIGHT) * (1.0 - phat @ shat)


def _retard_sphere(shat, phat, L_m, D_m):
    """EXACT spherical-front retardation: (|D shat - L phat| + L - D)/c."""
    cosmu = phat @ shat
    R = np.sqrt(D_m ** 2 + L_m ** 2 - 2.0 * D_m * L_m * cosmu)
    return (R + L_m - D_m) / C_LIGHT


def _retard_fresnel(shat, phat, L_m, D_m):
    """Fresnel (leading-curvature) retardation -- the brief's pi L_perp^2/(lam D)
    written as a time delay: plane + L^2 sin^2 mu / (2 c D)."""
    cosmu = phat @ shat
    return (L_m / C_LIGHT) * (1.0 - cosmu) + \
           (L_m ** 2 * (1.0 - cosmu ** 2)) / (2.0 * C_LIGHT * D_m)


RETARD = {"plane": None, "sphere": _retard_sphere, "fresnel": _retard_fresnel}


def residuals(spine, src, front="plane", earth_only=False, pulsar_only=False):
    """Zero-noise CW residuals, (npsr, ntoa).

    src: dict with costheta, phi, log10_h, log10_f, cosi, psi, phi0
         (+ D_m if front != 'plane'; + fdot optional; + L_m to override the
         spine's pulsar distances -- that is the fringe/nuisance axis).
    Sums over a LIST of source dicts if given a list.
    """
    if isinstance(src, (list, tuple)):
        out = 0.0
        for s in src:
            out = out + residuals(spine, s, front=front, earth_only=earth_only,
                                  pulsar_only=pulsar_only)
        return out

    t = spine["toas"]
    phat = spine["phat"]
    L_m = np.asarray(src.get("L_m", spine["L_m"]), float)
    shat = sky_unit(src["costheta"], src["phi"])
    Fp, Fx = antenna(shat, phat)

    h0 = 10.0 ** src["log10_h"]
    f = 10.0 ** src["log10_f"]
    cosi = src["cosi"]
    psi = src["psi"]
    phi0 = src["phi0"]
    fdot = src.get("fdot", 0.0)

    Ap = h0 * (1.0 + cosi ** 2) / 2.0
    Ax = h0 * cosi
    c2, s2 = np.cos(2.0 * psi), np.sin(2.0 * psi)
    Fp_r = Fp * c2 + Fx * s2
    Fx_r = -Fp * s2 + Fx * c2

    if front == "plane":
        dt = _retard_plane(shat, phat, L_m)
    else:
        dt = RETARD[front](shat, phat, L_m, src["D_m"])

    def _wave(tt):
        # residual (not strain): the 1/(2 pi f) integral of the sinusoid.
        ph = 2.0 * np.pi * (f * tt + 0.5 * fdot * tt ** 2) + phi0
        fe = f + fdot * tt
        sp = (Ap / (2.0 * np.pi * fe)) * np.sin(ph)
        sx = -(Ax / (2.0 * np.pi * fe)) * np.cos(ph)
        return sp, sx

    tE = t[None, :]
    tP = t[None, :] - dt[:, None]
    spE, sxE = _wave(tE)
    spP, sxP = _wave(tP)
    if earth_only:
        spP = sxP = 0.0
    if pulsar_only:
        spE = sxE = 0.0
    return Fp_r[:, None] * (spE - spP) + Fx_r[:, None] * (sxE - sxP)


def lnL(spine, data, src, front="plane"):
    """Exact zero-noise quadratic log-likelihood; lnL = 0 at a perfect template."""
    r = residuals(spine, src, front=front)
    d = (data - r) / spine["sigma"][:, None]
    return -0.5 * float(np.sum(d * d))


def snr2(spine, data):
    """optimal (zero-noise) SNR^2 of a data vector = 2 * |lnL(no signal)|."""
    d = data / spine["sigma"][:, None]
    return float(np.sum(d * d))


# ------------------------------------------------------------------
# 3. the fringe geometry (used by T1/T5)
# ------------------------------------------------------------------
def fringe_spacing(spine, src):
    """dL_p = c / (f (1 - cos mu_p))  -- one pulsar-term phase cycle, in metres."""
    shat = sky_unit(src["costheta"], src["phi"])
    f = 10.0 ** src["log10_f"]
    return C_LIGHT / (f * (1.0 - spine["phat"] @ shat))


def curvature_dL_eff(spine, src, D_m):
    """The T1 identity: to Fresnel order the curved front is EXACTLY a per-pulsar
    effective distance offset

        dL_eff,p = L_p^2 (1 + cos mu_p) / (2 D)

    because both curvature and a distance error act on the observables ONLY
    through the retardation dt_p.  Returned in metres.
    """
    shat = sky_unit(src["costheta"], src["phi"])
    cosmu = spine["phat"] @ shat
    return spine["L_m"] ** 2 * (1.0 + cosmu) / (2.0 * D_m)
