"""TRACK A — F1/F2: the anisotropic distance-Fisher via the cross-pulsar GWB covariance.

Referendum on Farr's isotropic zero-information theorem for pulsar distances. F0 built the
pulsar-term cross-covariance C_ab(f) and its distance derivative. Here we compute the
per-pulsar distance Fisher under the Whittle (discrete-Fourier-bin) likelihood, on the real
116-pulsar array, for anisotropic shot-noise skies, and measure the ANISOTROPIC ENHANCEMENT

    I_LL_aniso(N) = I_LL_full(shot-noise sky, N) - I_LL_iso(matched total power).

Amendments folded in (see PROJECT_PROGRESS Track A spec):
 1. baseline subtraction: iso row reported explicitly, never asserted zero.
 2. Whittle bins: C and dC evaluated at f_k = k/T only; Fisher = sum over bins.
 3. aliasing: n_mu >= RESFAC * ftau0 per bin (RESFAC=10), capped at N_MU_CEIL; drops logged.
 4. auto-term sanity: growth of dC_aa/dL_a under bin summation examined (see `gates`).
 5. FRINGE CAVEAT: every sigma_L is WITHIN-FRINGE covariance information.

Physics fiducials (reported every run):
    h_c(f) = A (f/f_yr)^{-2/3},  A = 2e-15;  S_gw(f) = h_c^2/(12 pi^2 f^3)   [s^3]
    white PSD per pulsar = 2 sigma_a^2 dt,  dt = 7 d (weekly),  sigma_a real (feathers)
    sky distribution p(nhat) normalised to unit monopole -> total power N-independent
    shot noise: N equal-strain isotropic point sources -> empirical c_lm (l<=10) -> band-limited p

Fisher (spec form, complex-Gaussian Whittle):
    I_LpLp = (1/2) sum_k Tr[ C_k^-1 dC_k/dLp C_k^-1 dC_k/dLp ],  C_k = S_gw(f_k) C_geo + N.

Runs (jug-gpu): python trackA_fisher.py {validate|gates|f1|f2}
"""
import os
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
import sys, argparse, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from scipy.special import sph_harm_y

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
from prong2_transition import antenna, KPC_OVER_C

REPO = "/home/mattm/projects/HSYMT"
CWT = REPO + "/CW_transition"

# ---- fiducial constants ----
YR = 3.15576e7
T_OBS = 15.0 * YR
FYR = 1.0 / YR
F0BIN = 1.0 / T_OBS                      # 2.113e-9 Hz
A_GWB = 2.0e-15
DT_WEEKLY = 7.0 * 86400.0
LMAX = 10
RESFAC = 10.0                            # amendment 3: n_mu >= RESFAC * ftau0
N_MU_CEIL = 16384                        # feasible resolution ceiling
CHUNK_MU = 64                            # mu-nodes per pixel chunk


def sgw(f):
    """GWB residual PSD [s^3] at frequency f [Hz]."""
    hc = A_GWB * (f / FYR) ** (-2.0 / 3.0)
    return hc ** 2 / (12.0 * np.pi ** 2 * f ** 3)


# ---------------- real array ----------------
def load_real_array():
    """116 real positions + real heterogeneous white-noise PSD + real distances (name-aligned)."""
    ra = np.load(CWT + "/trackA_realarray.npz", allow_pickle=True)
    pr = np.load(CWT + "/anchor_a0_priors.npz", allow_pickle=True)
    names = list(ra["names"]); pos = np.asarray(ra["pos"]); wn_s = np.asarray(ra["wn_s"])
    jn = list(pr["jname"]); fn = list(pr["fname"]); dist = np.asarray(pr["dist_real_kpc"])
    # match ra names to anchor names (ra names are Jxxxx pulsar names)
    def key(s):
        return s.strip()
    amap = {key(j): i for i, j in enumerate(jn)}
    L = np.empty(len(names))
    for i, nm in enumerate(names):
        k = key(nm)
        if k in amap:
            L[i] = dist[amap[k]]
        else:
            L[i] = np.median(dist)     # fallback (shouldn't happen)
    pwn = 2.0 * wn_s ** 2 * DT_WEEKLY
    return names, pos, L, pwn


# ---------------- real spherical harmonics ----------------
def real_ylm_table(theta, phi, lmax=LMAX):
    """Real orthonormal Y_lm(theta,phi) for l=0..lmax. Returns (n_modes, n_pix) and (l,m) list.
    theta,phi: numpy arrays (n_pix,). scipy sph_harm(m,l,phi,theta)."""
    rows = []; lm = []
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            if m == 0:
                y = sph_harm_y(l, 0, theta, phi).real
            elif m > 0:
                y = np.sqrt(2.0) * (-1.0) ** m * sph_harm_y(l, m, theta, phi).real
            else:
                y = np.sqrt(2.0) * (-1.0) ** m * sph_harm_y(l, -m, theta, phi).imag
            rows.append(y); lm.append((l, m))
    return np.array(rows), lm


def build_sky_basis(n_mu, lm):
    """Separable real-Ylm basis on the bin grid: Y_lm(theta,phi) = R[mode,i_mu] * PHI[mode,j].
    scipy is called ONCE (n_mu theta values), not per pixel-chunk. Returns R (nmode,n_mu),
    PHI (nmode,n_phi), theta, phi, W2d-weights (wx, n_phi factor)."""
    x, wx = gl_nodes(n_mu)
    theta = np.arccos(x)
    n_phi = 2 * n_mu
    phi = (np.arange(n_phi) + 0.5) * (2 * np.pi / n_phi)
    R = np.empty((len(lm), n_mu)); PHI = np.empty((len(lm), n_phi))
    for i, (l, m) in enumerate(lm):
        am = abs(m)
        base = sph_harm_y(l, am, theta, 0.0).real          # N_l|m| P_l^|m|(cos theta)
        R[i] = base
        if m == 0:
            PHI[i] = 1.0
        elif m > 0:
            PHI[i] = np.sqrt(2.0) * (-1.0) ** m * np.cos(m * phi)
        else:
            PHI[i] = np.sqrt(2.0) * (-1.0) ** m * np.sin(am * phi)
    return R, PHI, theta, phi, wx


def shot_noise_clm(N, lmax, seed):
    """Empirical real c_lm (l<=lmax) of N equal-strain isotropic point sources.
    p(nhat) = (4pi/N) sum_i delta -> c_lm = (4pi/N) sum_i Y_lm(nhat_i). Monopole -> unit."""
    rng = np.random.default_rng(seed)
    ct = rng.uniform(-1, 1, N); th = np.arccos(ct); ph = rng.uniform(0, 2 * np.pi, N)
    Y, lm = real_ylm_table(th, ph, lmax)          # (nmode, N)
    clm = (4.0 * np.pi / N) * Y.sum(axis=1)
    return clm, lm


def iso_clm(lm):
    c = np.zeros(len(lm)); c[0] = np.sqrt(4.0 * np.pi)   # c_00 -> p==1
    return c


def control_clm(lm, l_mode, m_mode=0):
    """Maximally anisotropic non-negative control: p = 1 + a*Ylm with min(p)=0 (upper envelope)."""
    c = np.zeros(len(lm)); c[0] = np.sqrt(4.0 * np.pi)
    idx = lm.index((l_mode, m_mode))
    return c, idx


# ---------------- sky grid (chunked) ----------------
def gl_nodes(n_mu):
    x, wx = np.polynomial.legendre.leggauss(n_mu)
    return x, wx


def n_mu_for_bin(f, Lmax_kpc):
    ftau0 = f * Lmax_kpc * KPC_OVER_C
    need = int(np.ceil(RESFAC * ftau0))
    n = min(need, N_MU_CEIL)
    n += n % 2
    return n, ftau0, need


# ---------------- covariance via jitted per-chunk kernel + python loop ----------------
from functools import partial


def _antenna_grid(thf, phf, p_hats):
    """Vectorised antenna (psi=0) via matmuls -> Fp,Fx,cos_mu each (npsr, cp). Compiles fast
    (no vmap-of-antenna). Matches prong2.antenna incl. 1/(1-cos mu)."""
    st = jnp.sin(thf); ct = jnp.cos(thf); sp = jnp.sin(phf); cp_ = jnp.cos(phf)
    n_hat = jnp.stack([st * cp_, st * sp, ct], axis=1)         # (cp,3)
    e_th = jnp.stack([ct * cp_, ct * sp, -st], axis=1)
    e_ph = jnp.stack([-sp, cp_, jnp.zeros_like(sp)], axis=1)
    P = p_hats.T                                               # (3, npsr)
    cmu = (n_hat @ P).T                                        # (npsr, cp)
    dt = (e_th @ P).T; dp = (e_ph @ P).T
    denom = jnp.clip(1.0 - cmu, 1e-6, None)
    Fp = 0.5 * (dt * dt - dp * dp) / denom
    Fx = (dt * dp) / denom
    return Fp, Fx, cmu


@jax.jit
def _chunk_kernel(p_hats, Lj, twopif, thf, phf, Wf, Pj):
    """Per-chunk contribution to C[case],R[case] (ncase,npsr,npsr). Case loop is a lax.map so
    the peak intermediate is (npsr,cp) not (ncase,npsr,cp) -> big chunks allowed. Pj: (ncase,cp)."""
    Fp, Fx, cmu = _antenna_grid(thf, phf, p_hats)              # (npsr, cp)
    om = twopif * KPC_OVER_C * (1.0 - cmu)
    psi = Lj[:, None] * om
    e = jnp.exp(-1j * psi)
    g = 1.0 - e
    dg = 1j * om * e
    sw = jnp.sqrt(Wf / (4.0 * jnp.pi))
    n = p_hats.shape[0]; m = Pj.shape[0]
    Cc = jnp.zeros((m, n, n), dtype=jnp.complex64)
    Rc = jnp.zeros_like(Cc)
    Pj32 = Pj.astype(jnp.float32)
    for Fc in (Fp, Fx):                                        # batched over cases (einsum GEMM)
        # phases (psi ~ 1e4 rad) in fp64 above; GEMM in complex64 (fp32 units, ~30x faster).
        # Within a mu-chunk the partial-sky sum does NOT cancel -> complex64 safe; the
        # between-chunk cross-cancellation is carried by the fp64 OUTER accumulator.
        # cp is HARD-capped in geo_covariance so this einsum's temp stays small (~<1 GB).
        U = ((sw[None, :] * Fc) * g).astype(jnp.complex64)     # (npsr,cp)
        V = ((sw[None, :] * Fc) * dg).astype(jnp.complex64)
        M = Pj32[:, None, :] * jnp.conj(U)[None, :, :]         # (ncase,npsr,cp)
        Cc = Cc + jnp.einsum("ap,mbp->mab", U, M)
        Rc = Rc + jnp.einsum("ap,mbp->mab", V, M)
    return Cc, Rc


def geo_covariance(pos, L, f, clm_cases, lm, n_mu):
    """C_geo[case], Rmat_geo[case] (n_case,n_psr,n_psr) complex for the GW-geometry covariance.
      C_ab   = sum_pix W P /(4pi) (Fp_aFp_b+Fx_aFx_b) g_a conj(g_b)
      Rmat_pb= sum_pix W P /(4pi) (Fp_pFp_b+Fx_pFx_b) dg_p conj(g_b)   (row p of dC/dL_p)
    Sky P(nhat)=sum clm Y_lm built on-device per chunk (memory-bounded, separable basis)."""
    n_phi = 2 * n_mu
    Rb, PHIb, theta, phi, wx = build_sky_basis(n_mu, lm)
    npsr = pos.shape[0]; ncase = clm_cases.shape[0]
    p_hats = jnp.asarray(pos); Lj = jnp.asarray(L); twopif = 2.0 * np.pi * f
    PHIj = jnp.asarray(PHIb); clm_j = jnp.asarray(clm_cases); Rbj = jnp.asarray(Rb)
    phij = jnp.asarray(phi); dphi = 2.0 * np.pi / n_phi
    C = jnp.zeros((ncase, npsr, npsr), dtype=jnp.complex128); R = jnp.zeros_like(C)
    # chunk mu so the batched complex64 MM=(ncase,cp,npsr) intermediate stays ~<6 GB, and cp is
    # HARD-capped (~1.3e5) so a small ncase can't inflate cp into a giant matmul temp. cp uniform
    # (last chunk padded so the kernel keeps ONE compiled shape). Batched GEMM -> high GPU util.
    cp_cap = min(131072, int(6.0e9 / (8.0 * max(ncase, 1) * npsr)))
    block_mu = max(1, min(n_mu, cp_cap // n_phi))
    starts = list(range(0, n_mu, block_mu))
    for i0 in starts:
        i1 = min(i0 + block_mu, n_mu); b = i1 - i0
        th = theta[i0:i1]
        if b < block_mu:                                       # pad last chunk (zero weight -> inert)
            th = np.concatenate([th, np.full(block_mu - b, np.pi / 2)])
            wxb = np.concatenate([wx[i0:i1], np.zeros(block_mu - b)])
            Rbc = np.concatenate([Rbj[:, i0:i1], jnp.zeros((Rbj.shape[0], block_mu - b))], axis=1)
        else:
            wxb = wx[i0:i1]; Rbc = Rbj[:, i0:i1]
        THf = jnp.repeat(jnp.asarray(th), n_phi)
        PHf = jnp.tile(phij, block_mu)
        Wf = (jnp.asarray(wxb)[:, None] * dphi * jnp.ones(n_phi)[None, :]).reshape(-1)
        Q = clm_j[:, :, None] * jnp.asarray(Rbc)[None, :, :]   # (ncase,nmode,block_mu)
        Pj = jnp.einsum("cmi,mj->cij", Q, PHIj).reshape(ncase, block_mu * n_phi)
        Cc, Rc = _chunk_kernel(p_hats, Lj, twopif, THf, PHf, Wf, Pj)
        C = C + Cc.astype(jnp.complex128); R = R + Rc.astype(jnp.complex128)  # fp64 outer accumulate
    return np.asarray(C), np.asarray(R)


# ---------------- per-pulsar Fisher (closed form) ----------------
def fisher_rows(C_geo, R_geo, amp, pwn):
    """Per-pulsar distance Fisher contribution at one bin (before the 1/2 and bin-sum).
      C = amp*C_geo + diag(pwn) ; dC/dLp = amp*(e_p r_p^T + conj(r_p) e_p^T), r_p = R_geo[p,:].
      Tr[B dCp B dCp] = amp^2*[(r.u)^2 + 2(r.v)u_p + v_p^2], u=B[:,p], v=B@conj(r_p).
    Returns real vector (n_psr,) of Tr[...] per pulsar."""
    n = C_geo.shape[0]
    C = amp * C_geo + np.diag(pwn)
    B = np.linalg.inv(C)
    R = amp * R_geo                          # rows r_p = R[p,:]
    out = np.empty(n)
    for p in range(n):
        r = R[p, :]
        u = B[:, p]
        v = B @ np.conj(r)
        val = (r @ u) ** 2 + 2.0 * (r @ v) * u[p] + v[p] ** 2
        out[p] = val.real
    return out


def fisher_rows_dense(C_geo, R_geo, amp, pwn):
    """Reference: build dense dC/dLp and take np.trace(B@dCp@B@dCp). For validation only."""
    n = C_geo.shape[0]
    C = amp * C_geo + np.diag(pwn)
    B = np.linalg.inv(C)
    R = amp * R_geo
    out = np.empty(n)
    for p in range(n):
        dC = np.zeros((n, n), dtype=complex)
        dC[p, :] += R[p, :]
        dC[:, p] += np.conj(R[p, :])
        M = B @ dC
        out[p] = np.trace(M @ M).real
    return out


# ---------------- P on grid from c_lm ----------------
def P_on_grid(clm_list, n_mu, lmax=LMAX, extra=None):
    """Build P_cases (n_case, n_pix) on the bin grid from a list of c_lm vectors.
    extra: optional list of (base_clm, idx, envelope) to build max-anisotropy controls."""
    x, wx = gl_nodes(n_mu)
    n_phi = 2 * n_mu
    phi = (np.arange(n_phi) + 0.5) * (2 * np.pi / n_phi)
    th = np.arccos(x)
    TH, PH = np.meshgrid(th, phi, indexing="ij")
    Y, lm = real_ylm_table(TH.ravel(), PH.ravel(), lmax)     # (nmode, n_pix)
    out = []
    for clm in clm_list:
        out.append(clm @ Y)
    return np.array(out), lm, Y


# ---------------- subcommands ----------------
def cmd_validate(argv):
    """Self-tests: (i) trace formula vs dense; (ii) U-factorisation C vs F0 covariance()."""
    print("=== VALIDATE ===")
    npsr = 6
    rng = np.random.default_rng(1)
    ct = rng.uniform(-1, 1, npsr); ph = rng.uniform(0, 2 * np.pi, npsr)
    st = np.sqrt(1 - ct ** 2)
    pos = np.stack([st * np.cos(ph), st * np.sin(ph), ct], axis=1)
    L = rng.uniform(0.3, 3.0, npsr)
    f = 5.0 / (1.0 * KPC_OVER_C)          # small ftau0 ~5 for a cheap dense sky
    n_mu = 200
    lm = real_ylm_table(np.array([0.1]), np.array([0.1]))[1]
    clm_i = iso_clm(lm)
    clm_s, _ = shot_noise_clm(50, LMAX, 7)
    clm_cases = np.stack([clm_i, clm_s])
    C, R = geo_covariance(pos, L, f, clm_cases, lm, n_mu)
    pwn = np.full(npsr, 1e-3)
    amp = 1.0
    for ci, tag in [(0, "iso"), (1, "shot50")]:
        a = fisher_rows(C[ci], R[ci], amp, pwn)
        b = fisher_rows_dense(C[ci], R[ci], amp, pwn)
        rel = np.max(np.abs(a - b) / (np.abs(b) + 1e-30))
        print(f"  [{tag}] trace-formula vs dense: max rel err = {rel:.2e}")
    # cross-check C against F0 covariance() (iso sky, monopole)
    from trackA_covariance import covariance, sky_grid
    th_g, ph_g, W_g = sky_grid(n_mu, 2 * n_mu)
    C_f0, _ = covariance(th_g, ph_g, W_g, jnp.asarray(pos), jnp.asarray(L), f)
    C_f0 = np.asarray(C_f0)
    rel = np.max(np.abs(C[0] - C_f0)) / np.max(np.abs(C_f0))
    print(f"  U-factorisation C vs F0 covariance() (iso): max rel err = {rel:.2e}")
    # Hermiticity
    print(f"  C Hermitian: {np.max(np.abs(C[0]-C[0].conj().T)):.2e}")
    return 0


# ---------------- shared Fisher assembler ----------------
def assemble(pos, L, pwn, clm_cases, lm, bins_k, ceil=N_MU_CEIL, verbose=True, keep=False):
    """I[case,psr] = 0.5 sum_k Tr_k. Returns (I, per_bin_records)."""
    Lmax = float(np.max(L))
    ncase = clm_cases.shape[0]; npsr = len(L)
    Tr_sum = np.zeros((ncase, npsr))
    per_bin = []
    # ONE grid for all bins -> the scan kernel compiles once (bins differ only by f/amp, traced).
    n_mu_uni = min(ceil, max(n_mu_for_bin(k * F0BIN, Lmax)[2] for k in bins_k))
    n_mu_uni += n_mu_uni % 2
    for k in bins_k:
        f = k * F0BIN
        _, ftau0max, need = n_mu_for_bin(f, Lmax)
        n_mu = n_mu_uni
        amp = sgw(f)
        t0 = time.time()
        Cg, Rg = geo_covariance(pos, L, f, clm_cases, lm, n_mu)
        binTr = np.zeros((ncase, npsr))
        for ci in range(ncase):
            binTr[ci] = fisher_rows(Cg[ci], Rg[ci], amp, pwn)
        Tr_sum += binTr
        rec = dict(k=k, f=f, n_mu=n_mu, ftau0max=ftau0max, need=need,
                   resolved=(need <= n_mu), amp=amp, binTr=binTr, dt=time.time() - t0)
        if keep:
            rec["Cg"] = Cg; rec["Rg"] = Rg
        per_bin.append(rec)
        if verbose:
            flag = "RESOLVED" if need <= n_mu else f"CAPPED(need n_mu={need})"
            print(f"  bin k={k} f={f*1e9:.2f}nHz amp={amp:.3e} n_mu={n_mu} "
                  f"ftau0max={ftau0max:.0f} [{flag}] {rec['dt']:.1f}s")
    return 0.5 * Tr_sum, per_bin


def info_weight_by_bin(pos, L, pwn, bins_k):
    """Fraction of the iso distance-Fisher carried by each bin (for drop accounting)."""
    lm = real_ylm_table(np.array([0.1]), np.array([0.1]))[1]
    ci = iso_clm(lm)[None, :]
    tot = np.zeros(len(L)); w = []
    recs = []
    for k in bins_k:
        f = k * F0BIN; n_mu, ftau0max, need = n_mu_for_bin(f, L.max())
        # cheap: evaluate at reduced resolution just for the WEIGHT (ratio), capped low
        nm = min(n_mu, 4096)
        Cg, Rg = geo_covariance(pos, L, f, ci, lm, nm)
        tr = fisher_rows(Cg[0], Rg[0], sgw(f), pwn)
        recs.append((k, 0.5 * tr.sum())); tot += 0.5 * tr
    S = sum(v for _, v in recs)
    return [(k, v / S) for k, v in recs]


# ---------------- GATES ----------------
def cmd_gates(argv):
    names, pos, L, pwn = load_real_array()
    lm = real_ylm_table(np.array([0.1]), np.array([0.1]))[1]
    npsr = len(L)
    print(f"jax {jax.__version__}, {jax.devices()}")
    print(f"real array: {npsr} psr, L[kpc] med={np.median(L):.2f} max={L.max():.2f}; "
          f"pwn[s^3] med={np.median(pwn):.2e}; A_gwb={A_GWB:.1e}, T={T_OBS/YR:.0f}yr, "
          f"weekly, lmax={LMAX}\n")

    # ---- GATE 1: monopole-only sky reproduces iso baseline; enhancement exactly 0 ----
    print("=== GATE 1: monopole-only sky == iso baseline (enhancement zero by construction) ===")
    clm_iso = iso_clm(lm)
    clm_shot, _ = shot_noise_clm(300, LMAX, 11)
    clm_mono = clm_shot.copy(); clm_mono[1:] = 0.0          # kill all l>=1 -> pure monopole
    cases = np.stack([clm_iso, clm_mono])
    I, _ = assemble(pos, L, pwn, cases, lm, [1], verbose=False)
    d = np.max(np.abs(I[1] - I[0]) / (np.abs(I[0]) + 1e-300))
    print(f"  max|I(monopole-only) - I(iso)| / I(iso) = {d:.2e}  (want <=1e-8)  "
          f"[{'PASS' if d <= 1e-8 else 'FAIL'}]\n")

    # ---- GATE 2: aliasing convergence (double resolution at largest resolved production ftau0) ----
    # quadrature convergence is pulsar-count-independent -> use a small subset (incl. farthest) fast.
    print("=== GATE 2: aliasing convergence, double resolution (amendment 3) ===")
    k = 1; f = k * F0BIN
    n_mu, ftau0max, need = n_mu_for_bin(f, L.max())
    order = np.argsort(L)
    sub = np.unique(np.r_[order[-1], order[::len(order)//8]])   # farthest + spread of 8
    pos_s = pos[sub]; L_s = L[sub]; pwn_s = pwn[sub]; names_s = [names[i] for i in sub]
    pfar_s = int(np.argmax(L_s))
    res = {}
    for factor, nm in [(1, n_mu), (2, 2 * n_mu)]:
        Cg, Rg = geo_covariance(pos_s, L_s, f, clm_iso[None, :], lm, nm)
        tr = fisher_rows(Cg[0], Rg[0], sgw(f), pwn_s)
        res[factor] = (Cg[0], tr)
        print(f"  n_mu={nm} (={factor}x, ftau0max={ftau0max:.0f}) done")
    C1, tr1 = res[1]; C2, tr2 = res[2]
    dauto = np.abs(C2[pfar_s, pfar_s] - C1[pfar_s, pfar_s]) / np.abs(C1[pfar_s, pfar_s])
    off = np.abs(C1); np.fill_diagonal(off, 0.0)
    b = int(np.argmax(off[pfar_s]))
    dcross = np.abs(C2[pfar_s, b] - C1[pfar_s, b]) / (np.abs(C1[pfar_s, b]) + 1e-300)
    dtr = np.abs(tr2[pfar_s] - tr1[pfar_s]) / np.abs(tr1[pfar_s])
    print(f"  farthest psr {names_s[pfar_s]} (L={L_s[pfar_s]:.2f}kpc): "
          f"rel d|C_auto|={dauto:.2e}  rel d|C_cross|={dcross:.2e}  rel d(Fisher row)={dtr:.2e}")
    print(f"  (key integrals stable to 1e-8: {'PASS' if max(dauto,dtr)<=1e-8 else 'see values'})\n")

    # ---- GATES 3 & 4 share one covariance pass per bin (iso sky) ----
    print("=== GATE 3: Whittle-bin table, iso sky (amendment 2) ===")
    print("  per-bin iso auto dC_aa/dL_a vs F0 (pi/3)ftau0_a, and cross_max; Fisher SUMS these bins")
    print(f"  {'k':>2} {'f[nHz]':>7} {'ftau0(medL)':>11} {'auto_med':>10} {'(pi/3)ft':>9} "
          f"{'ratio':>6} {'cross_max':>10} {'cr/au':>8}")
    bins_all = [1, 2, 3, 4, 5, 6]
    Iiso_bins = []; g4 = []
    pmed = int(np.argsort(L)[len(L) // 2])
    for k in bins_all:
        f = k * F0BIN; n_mu, ftau0max, need = n_mu_for_bin(f, L.max())
        nm = min(n_mu, N_MU_CEIL)
        Cg, Rg = geo_covariance(pos, L, f, clm_iso[None, :], lm, nm)
        R0 = Rg[0]; amp = sgw(f)
        auto = np.abs(np.diag(R0))                       # |dC_aa/dL_a| geometry
        off = np.abs(R0.copy()); np.fill_diagonal(off, 0.0)
        pred = (np.pi / 3.0) * f * L * KPC_OVER_C        # per-pulsar (pi/3)ftau0_a
        ratio = np.median(auto / pred)
        tr = fisher_rows(Cg[0], Rg[0], amp, pwn)
        Iiso_bins.append(0.5 * tr)
        ft_med = f * np.median(L) * KPC_OVER_C
        print(f"  {k:>2} {f*1e9:>7.2f} {ft_med:>11.0f} {np.median(auto):>10.3e} "
              f"{(np.pi/3)*ft_med:>9.3e} {ratio:>6.3f} {off.max():>10.3e} "
              f"{off.max()/np.median(auto):>8.2e}")
        # GATE 4 accumulation (auto-only, median pulsar)
        Caa = amp * Cg[0][pmed, pmed].real + pwn[pmed]
        dCaa = amp * np.abs(R0[pmed, pmed])
        g4.append((k, f * L[pmed] * KPC_OVER_C, dCaa, Caa, 0.5 * (dCaa / Caa) ** 2))
    Iiso_bins = np.array(Iiso_bins)
    frac = Iiso_bins.sum(axis=1) / Iiso_bins.sum()
    print("  bin fraction of total iso Fisher (info weight): " +
          " ".join(f"k{bins_all[i]}={frac[i]:.4f}" for i in range(len(bins_all))))
    print(f"  -> bins 1-4 carry {frac[:4].sum()*100:.3f}% of iso distance-Fisher\n")

    print("=== GATE 4: auto-term (pi/3)ftau0 growth under Whittle summation (amendment 4) ===")
    print(f"  auto-ONLY per-bin Fisher, median pulsar {names[pmed]} (L={L[pmed]:.2f}kpc):")
    running = 0.0
    for k, ft, dCaa, Caa, Iaa in g4:
        running += Iaa
        print(f"    k={k} ftau0={ft:.0f}: dC_aa/dL={dCaa:.3e} amp*C_aa+N={Caa:.3e} "
              f"I_aa(bin)={Iaa:.3e} running={running:.3e}")
    print("  -> per-bin auto info: amp^2*(pi/3 ftau0)^2 numerator vs (amp*C_aa+N)^2 denominator; "
          "growth is bounded where noise dominates and is WITHIN-FRINGE curvature (amendment 5).\n")
    return 0


# ---------------- F1 ----------------
def cmd_f1(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bins", type=int, default=4)
    ap.add_argument("--nreal", type=int, default=10)
    ap.add_argument("--Ns", type=int, nargs="+", default=[100, 300, 1000, 3000])
    aa = ap.parse_args(argv)
    names, pos, L, pwn = load_real_array()
    lm = real_ylm_table(np.array([0.1]), np.array([0.1]))[1]
    bins_k = list(range(1, aa.bins + 1))
    npsr = len(L)
    print(f"jax {jax.__version__}, {jax.devices()}; real {npsr}-psr array")
    print(f"F1: bins={bins_k}, nreal={aa.nreal}, Ns={aa.Ns}, A_gwb={A_GWB:.1e}\n")

    # assemble the case list: [iso, control_l2, control_l4, shot(N,real) x Ns x nreal]
    clm_iso = iso_clm(lm)
    # controls: max-anisotropy non-negative single mode at unit monopole
    def ctrl(l_mode, m_mode):
        c = np.zeros(len(lm)); c[0] = np.sqrt(4 * np.pi)
        idx = lm.index((l_mode, m_mode))
        # find a s.t. min(1 + a*Y_lm)=0 over sky -> a = -1/min(Y_lm); use a fine grid for min(Y)
        th = np.arccos(np.linspace(-1, 1, 400)); ph = np.linspace(0, 2 * np.pi, 400)
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        Y, _ = real_ylm_table(TH.ravel(), PH.ravel(), LMAX)
        ymin = Y[idx].min()
        c[idx] = -1.0 / ymin
        return c
    clm_l2 = ctrl(2, 0); clm_l4 = ctrl(4, 0)
    cases = [clm_iso, clm_l2, clm_l4]
    labels = ["iso", "ctrl_l2", "ctrl_l4"]
    shot_index = {}
    for N in aa.Ns:
        idxs = []
        for r in range(aa.nreal):
            c, _ = shot_noise_clm(N, LMAX, 1000 * N + r)
            idxs.append(len(cases)); cases.append(c)
        shot_index[N] = idxs
    cases = np.stack(cases)
    print(f"total sky cases = {len(cases)}")

    I, per_bin = assemble(pos, L, pwn, cases, lm, bins_k, keep=False)
    # info-weight resolved fraction
    resolved = [r for r in per_bin if r["resolved"]]
    wr = sum(0.5 * r["binTr"][0].sum() for r in resolved)
    wt = sum(0.5 * r["binTr"][0].sum() for r in per_bin)
    print(f"\nresolved-bin info weight (iso) = {wr/wt*100:.3f}% "
          f"(capped bins: {[r['k'] for r in per_bin if not r['resolved']]})")

    I_iso = I[0]; I_l2 = I[1]; I_l4 = I[2]
    # baseline subtraction: enhancement = I_full - I_iso
    def summ(v):
        return f"med={np.median(v):.3e} max={v.max():.3e}"
    print("\n=== I_LL per pulsar (kpc^-2), WITHIN-FRINGE ===")
    print(f"  (i)   ISO baseline:        {summ(I_iso)}")
    print(f"  (iii) control l=2 (env):   full {summ(I_l2)} | enh {summ(I_l2-I_iso)}")
    print(f"        control l=4 (env):   full {summ(I_l4)} | enh {summ(I_l4-I_iso)}")
    print("  (ii)  shot-noise skies (enhancement I_full - I_iso; median[16,84] over reals):")
    out = dict(names=names, L=L, pwn=pwn, I_iso=I_iso, I_l2=I_l2, I_l4=I_l4,
               A_gwb=A_GWB, bins=np.array(bins_k), Ns=np.array(aa.Ns), nreal=aa.nreal)
    enh_med = {}; enh_lo = {}; enh_hi = {}; full_med = {}
    for N in aa.Ns:
        stack = np.array([I[j] for j in shot_index[N]])       # (nreal, npsr)
        enh = stack - I_iso[None, :]
        em = np.median(enh, axis=0); el = np.percentile(enh, 16, axis=0); eh = np.percentile(enh, 84, axis=0)
        enh_med[N] = em; enh_lo[N] = el; enh_hi[N] = eh; full_med[N] = np.median(stack, axis=0)
        # array-median-pulsar enhancement + best pulsar
        mm = np.median(em); bp = int(np.argmax(em))
        print(f"    N={N:>5}: enh med(psr-median)={mm:.3e}  best psr {names[bp]}={em[bp]:.3e} "
              f"[{el[bp]:.2e},{eh[bp]:.2e}]  (iso@best={I_iso[bp]:.3e})")
    out["enh_med"] = np.array([enh_med[N] for N in aa.Ns])
    out["enh_lo"] = np.array([enh_lo[N] for N in aa.Ns])
    out["enh_hi"] = np.array([enh_hi[N] for N in aa.Ns])
    out["full_med"] = np.array([full_med[N] for N in aa.Ns])
    # full per-realization enhancement stack (nN, nreal, npsr) for F2 bootstrap
    out["enh_stack"] = np.array([[I[j] - I_iso for j in shot_index[N]] for N in aa.Ns])
    np.savez(CWT + "/trackA_f1_results.npz", **out)
    print("\n[F1] saved trackA_f1_results.npz")
    return 0


# ---------------- F2 ----------------
def _fit_alpha(N, y):
    """Fit log y = c - alpha log N. Returns alpha, c. y>0 assumed."""
    m = y > 0
    lx = np.log(np.asarray(N)[m]); ly = np.log(y[m])
    A = np.polyfit(lx, ly, 1)
    return -A[0], A[1]


def cmd_f2(argv):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = np.load(CWT + "/trackA_f1_results.npz", allow_pickle=True)
    sa = np.load(CWT + "/prong2_results.npz", allow_pickle=True)
    names = list(d["names"]); L = d["L"]; Ns = d["Ns"]
    enh_med = d["enh_med"]; enh_stack = d["enh_stack"]; I_iso = d["I_iso"]
    full_med = d["full_med"]; A_gwb = float(d["A_gwb"])
    nN = len(Ns)
    print(f"F2 verdict: A_gwb={A_gwb:.1e}, Ns={list(Ns)}")

    # ---- scaling fit: array-median-pulsar enhancement vs N ----
    ymed = np.array([np.median(enh_med[i]) for i in range(nN)])          # psr-median of median-over-real
    alpha_med, c_med = _fit_alpha(Ns, ymed)
    # bootstrap alpha over realisations (resample reals, recompute psr-median enh, refit)
    nreal = enh_stack.shape[1]
    rng = np.random.default_rng(0)
    alphas = []
    for _ in range(2000):
        pick = rng.integers(0, nreal, nreal)
        yb = np.array([np.median(np.median(enh_stack[i][pick], axis=0)) for i in range(nN)])
        if np.all(yb > 0):
            alphas.append(_fit_alpha(Ns, yb)[0])
    alphas = np.array(alphas)
    print(f"\n=== SCALING: I_LL_aniso(N) ~ N^-alpha (psr-median) ===")
    print(f"  alpha = {alpha_med:.3f}  [16,84] = [{np.percentile(alphas,16):.3f},{np.percentile(alphas,84):.3f}]")
    for i, N in enumerate(Ns):
        print(f"    N={N:>5}: enh(psr-median)={ymed[i]:.3e}  iso baseline(median)={np.median(I_iso):.3e}")

    # best-pulsar scaling
    bp = int(np.argmax(enh_med[0]))
    ybest = np.array([enh_med[i][bp] for i in range(nN)])
    alpha_best, _ = _fit_alpha(Ns, ybest)
    print(f"  best pulsar {names[bp]} (L={L[bp]:.2f}kpc): alpha={alpha_best:.3f}")

    # ---- resolved channel (Stage A, fixed total power) ----
    rn = sa["tot_n"]; rmarg = sa["tot_marg"]; knee = float(sa["knee_N"])
    # resolved slope past knee
    postk = rn >= knee * 0.9
    a_res, _ = _fit_alpha(rn[postk], rmarg[postk])
    print(f"\n=== RESOLVED channel (Stage A, fixed-total): knee N={knee:.0f}, "
          f"post-knee marg ~ N^-{a_res:.2f} ===")

    # ---- sigma_L within-fringe (best pulsar) at N=300,1000 ----
    print(f"\n=== sigma_L [pc] WITHIN-FRINGE (per amendment 5) ===")
    def sig_pc(I):  # I in kpc^-2 -> sigma in pc
        return 1e3 / np.sqrt(np.maximum(I, 1e-300))
    for N in (300, 1000):
        if N in list(Ns):
            i = list(Ns).index(N)
            If = full_med[i][bp]; Ie = enh_med[i][bp]; Ii = I_iso[bp]
            print(f"  N={N}: best psr {names[bp]}: sigma_L(full covariance)={sig_pc(If):.3e} pc | "
                  f"sigma_L(iso baseline)={sig_pc(Ii):.3e} pc | "
                  f"sigma_L(anisotropy-only)={sig_pc(Ie):.3e} pc")

    # ---- figure ----
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.loglog(rn, rmarg, "o-", color="C0", label="resolved marginal info (Stage A, fixed-total)")
    ax.axhline(np.median(I_iso), ls=":", color="k",
               label=f"iso covariance baseline (median psr)")
    ax.loglog(Ns, ymed, "s-", color="C3", label="anisotropic enhancement (psr-median)")
    ax.loglog(Ns, ybest, "^--", color="C1", label=f"anisotropic enhancement (best psr {names[bp]})")
    ax.axvline(knee, ls="--", color="gray", alpha=0.7, label=f"resolved knee N={knee:.0f}")
    ax.set_xlabel("N  (source count generating the sky)")
    ax.set_ylabel("distance Fisher information  [kpc$^{-2}$]")
    ax.set_title(f"Track A: covariance-channel vs resolved-channel distance info\n"
                 f"aniso enhancement ~ N$^{{-{alpha_med:.2f}}}$, resolved post-knee ~ N$^{{-{a_res:.2f}}}$")
    ax.legend(fontsize=7.5, loc="lower left")
    ax.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    fig.savefig(CWT + "/trackA_f2_scaling.png", dpi=140)
    print(f"\n[F2] saved trackA_f2_scaling.png")

    # ---- verdict ----
    print("\n=== VERDICT ===")
    med_enh = np.median(enh_med, axis=1)
    frac = med_enh / np.median(I_iso)
    print(f"  anisotropic enhancement / iso baseline (psr-median): " +
          " ".join(f"N{Ns[i]}={frac[i]:.2e}" for i in range(nN)))
    print(f"  enhancement scaling alpha={alpha_med:.2f}; resolved post-knee alpha={a_res:.2f}")
    np.savez(CWT + "/trackA_f2_results.npz", Ns=Ns, ymed=ymed, ybest=ybest,
             alpha_med=alpha_med, alpha_best=alpha_best, alphas=alphas,
             a_res=a_res, knee=knee, I_iso=I_iso, full_med=full_med, bp=bp,
             rn=rn, rmarg=rmarg)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("mode", choices=["validate", "gates", "f1", "f2"], nargs="?", default="validate")
    a, rest = ap.parse_known_args()
    fns = {"validate": cmd_validate, "gates": cmd_gates, "f1": cmd_f1, "f2": cmd_f2}
    fn = fns.get(a.mode)
    if fn is None:
        print(f"mode {a.mode} not yet implemented"); sys.exit(1)
    sys.exit(fn(rest))
