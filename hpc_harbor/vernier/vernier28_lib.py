#!/usr/bin/env python3
"""
VERNIER-28 library -- the three recovery ladders, on ONE universe.

EXPLORATORY, NOT pre-registered.

Clean room: imports vernier27_lib / vernier25 / vernier24 / vernier16 /
vernier6 / vernier_toy and edits none of them. Everything new is here.

WHAT THIS FILE ADDS
  R1  BASELINE   -- V27's chain (= V25's chain + V27's four corrections),
                    re-implemented VERBATIM here for one reason only: the
                    control needs to report the GRID INDEX of every component
                    so that all three ladders can be scored by one accounting
                    function. `validate_R1()` asserts bit-equality with
                    vernier27_lib.chain on the same universe. The control is
                    not improved, only instrumented.
  R2  JOINT-BIN  -- all sources sharing a bin fitted SIMULTANEOUSLY. Positions
                    are found by greedy CLEAN and then refined by space-
                    alternating conditional maximisation (SAGE): each source's
                    sky position is re-solved against the residual of all the
                    OTHERS, with their amplitudes re-fitted, until the set stops
                    moving. M is scanned 1..M_MAX and chosen by a null-
                    CALIBRATED admission floor -- the recut floor of the same
                    statistic's INCREMENT, measured on pure-noise realisations
                    of the same bin -- so "the M the data supports" is a
                    measured threshold, not an information criterion asserted
                    by me, and it has the same shape as R3's.
  R3  + PULSAR   -- R2 plus the pulsar-term block as a SECOND observable, with
                    the pulsar distances MARGINALISED over their 10 pc priors.

THE PULSAR BLOCK, AND WHY IT IS A SEPARATE OBSERVABLE
  VERNIER-24's `dvec` already treats a source's pulsar term as a distinct data
  block from its Earth term (its column 6 -- the source's own frequency -- is
  what makes it independent information). This file makes that block explicit
  in the toy's per-bin data:

      y_a = sum_k E_ak R_AMP env_ak exp(i phi_ak(L_a))  +  n_a ,   n_a ~ CN(0,1)

  E_ak is exactly the Earth-term contribution of source k to pulsar a that
  `inject` already puts in the Earth block, R_AMP = 1.196 is VERNIER-8's
  pulsar/Earth residual ratio, phi is the EXACT retarded circular-chirp phase
  (vernier_toy._retarded, the same function certify() uses through
  vernier6.phase_env) and env is its sinc envelope. Every source in the bin
  contributes -- all 300 -- so the pulsar block carries the same confusion the
  Earth block does.

  R1 and R2 never look at y. That is the standing VERNIER-16/22 Earth-only
  limitation and it keeps R1 bit-identical to V25/V27. R3 uses both blocks and
  marginalises L_a over the SKA-era 10 pc parallax prior:

      ln Z_a = ln \\int dL pi(L) exp( 2 Re(y_a^* yhat_a(L)) - |yhat_a(L)|^2 )

  The M sources in a bin share ONE unknown L_a per pulsar, so for M >= 2 the
  block is over-determined and informative: this is VERNIER-23 Part 1's degree-
  of-freedom surplus, spent on RECOVERY instead of on certification. Whether it
  pays is the measurement.

  R3's score is T(m) = F_joint(m) + sum_a ln Z_a(m), both in the same
  2 x log-likelihood units, and it is null-calibrated exactly as R2's is.

TWO DELIBERATE DEPARTURES FROM VERNIER-27, BOTH CONSERVATIVE
  1. SUB_SPAN = 0. V27 gives packed seeds sub-bin frequency offsets. The toy's
     F-statistic searches at the bin centre only, so a sub-bin offset is a
     frequency the search cannot measure -- and R3's comb model is exquisitely
     sensitive to it (2 pi df tau ~ 35 rad at df = 0.35/T). Handing R3 the true
     sub-bin frequency would be handing it information the toy's search cannot
     have. Here every source in a bin sits at the bin centre, so the model
     frequency is exactly the searched frequency for every ladder. The cost is
     that packed sources share a comb-spacing family, which makes the vernier's
     commensurability lever WEAKER than V27's construction, not stronger.
  2. ONE CHIRP MASS. V25 draws log10 Mc uniformly in 8.5-9.5 per source. The
     retarded pulsar-term phase carries a term (3 pi / 8) f tau^2 / tc worth
     ~90 rad, and tc ~ Mc^(-5/3) varies by 46x across that prior -- so a model
     that used each source's TRUE Mc would be using an oracle no fit has, while
     one that used a common Mc would break on real sources and not on phantoms,
     which is the same oracle wearing a hat. Here every source in the sky and in
     every model carries log10 Mc = 9.0, so the comb model is exact for
     genuine and phantom components alike and no ladder is told anything. The
     cost is that all combs in a bin share one chirp evolution, which makes them
     MORE alike and the vernier's commensurability lever WEAKER -- conservative
     for R3, which is the ladder under test.
  3. The genuine pool is the AVAILABLE set. Conversion efficiency is
     recovered/available, so "available" and "recovered" have to name the same
     set: the pool a component may be matched against is exactly the set of
     sources above the detection floor in the unsubtracted lambda_eff budget
     (V25.k_available's criterion). The chance-match floor of that pool is
     reported per cell.
"""
import numpy as np

import vernier27_lib as L27
import vernier25 as V25
import vernier24 as W
import vernier16 as S
import vernier6 as X
import vernier_toy as V

YR, DF, T_OBS, C_PC_S = W.YR, W.DF, W.T_OBS, V.C_PC_S
KAPPA, R_AMP = W.KAPPA, W.R_AMP
NGRID, TOL_RAD = W.NGRID, W.TOL_RAD
PRIOR_PC = W.PRIOR_PC
N_BRIGHT = V25.N_BRIGHT
BAND = V25.BAND
RHO_SEED = V25.RHO_SEED
N_NULL_C = V25.N_NULL_C
NPTS_CERT = V25.NPTS

# ---------------------------------------------------------------- V28 controls
M_MAX      = 12        # brief: M scanned over 1..min(N_bin, 12)
N_SWEEP    = 2         # SAGE sweeps per M in the Earth-only joint fit (R2)
N_SWEEP_P  = 1         # combined-objective sweeps on top of R2's solution (R3)
N_CAND     = 32        # spread candidate positions offered to R3 per source
CAND_SEP   = 0.15      # rad; candidates are separated by at least this (= TOL_RAD)
GUARD_SEP  = 0.02      # rad; two fitted components may not collapse closer
NPTS_L     = 512       # EVEN, sub-fringe rule (V24): distance grid for ln Z
NJ_NULL    = 12        # pure-noise fits calibrating R2's model-selection penalty
NJ_NULL_P  = 6         # pure-noise band fits calibrating R3's admission floor
MC_FIX     = 9.0       # log10 Mc, IDENTICAL for every source (see below)
COMB_CAP   = 16        # combs carried into a certification (V25 used 8, V27 32)
MAX_TURN   = 36
RIDGE      = 1e-10


# ============================================================ sky construction
def seed_plan(K, band):
    """
    V27's seed packing with SUB_SPAN = 0: K seeds spread ceil/floor(K/nb) per
    bin, every one of them AT ITS BIN CENTRE (departure 1 of the docstring).
    """
    nb = len(band)
    base, extra = divmod(int(K), nb)
    return [tuple([float(f)] * (base + (1 if b < extra else 0)))
            for b, f in enumerate(band)]


def available(rng, arr, bins, top=40):
    """
    The AVAILABLE set, per bin: source indices above the detection floor in the
    lambda_eff-consistent confusion budget with nothing subtracted. V25's
    k_available criterion, but returning WHICH sources rather than how many.
    """
    out = []
    for B in bins:
        fl = V25.floor_for(rng, arr, B)
        thr = np.sqrt(max(fl["floor"] + fl["err"], 0.0))
        idx = [i for i in range(min(B["M"], top))
               if np.sqrt(np.sum(W.rho_budget(B, i, set()) ** 2)) >= thr]
        out.append(set(idx))
    return out


def pterm_block(rng, arr, B, null=False):
    """
    The pulsar-term observable of one bin (see the file docstring).

    Built from EVERY source in the bin, at the exact retarded phase of the true
    pulsar distances, with the sinc envelope of the retarded-frequency offset.
    """
    n = B["Tn"].shape[0]
    y = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
    if null:
        return y
    rho_arr = np.sqrt((B["rho0"] ** 2).sum(axis=0))
    E = B["Tn"] * rho_arr[None, :]                       # (n_psr, M) Earth block
    Lst = arr["L"]
    for k in range(B["M"]):
        omc = B["omc"][:, k]
        fk = float(B["fk"][k]) if "fk" in B else float(B["f"])
        tc = float(V.t_merge(fk, float(B["mc"][k])))
        ph, env = _abs_phase(Lst, omc, fk, tc)
        y = y + R_AMP * E[:, k] * env * np.exp(1j * ph)
    return y


def _abs_phase(L_pc, omc, f, tc):
    """
    ABSOLUTE retarded pulsar-term phase and its band envelope, for distance
    L_pc (scalar or array broadcast against omc) and geometry omc.

    Uses vernier_toy._retarded -- the same exact circular-chirp propagator that
    vernier6.phase_env (and hence certify) uses. The envelope is the sinc of the
    offset between the retarded carrier and the source's own frequency, which is
    the frequency the bin is analysed at.
    """
    tau = V.tau_of(L_pc, omc)
    ph, fo = V._retarded(tau, f, tc, V.PTERM_STATIC)
    return ph, np.sinc((fo - f) * T_OBS)


def universe(seed, arr, gamma, K_seed, noise_seed, band=BAND, null=False,
             mc_fix=MC_FIX):
    """
    ONE realisation, built in exactly vernier27_lib.chain's RNG order so that
    R1 run through L27.chain(seed, ..., noise_seed=noise_seed) sees bit-identical
    bins, floors and Earth data. The pulsar block is drawn from a SEPARATE
    stream so that adding it perturbs nothing R1 or R2 sees.
    """
    sky_rng = np.random.default_rng(seed)
    plan = seed_plan(K_seed, band)
    bins = [L27.build_bin(sky_rng, arr, f, gamma, seed_fs=plan[b])
            for b, f in enumerate(band)]
    if mc_fix is not None:
        for B in bins:
            B["mc"][:] = mc_fix
    # WARM the floor cache off a throwaway stream first. floor_for() consumes rng
    # only on a cache MISS, so an un-warmed cache would put the ladders' noise
    # streams at different offsets and break the RNG matching between them.
    wrng = np.random.default_rng(noise_seed + 700000)
    for B in bins:
        V25.floor_for(wrng, arr, B)
    rng = np.random.default_rng(noise_seed)
    floors = [V25.floor_for(rng, arr, B) for B in bins]
    data = [W.inject(rng, B, null=null) for B in bins]
    yrng = np.random.default_rng(noise_seed + 500000)
    ydata = [pterm_block(yrng, arr, B, null=null) for B in bins]
    arng = np.random.default_rng(noise_seed + 900000)
    avail = available(arng, arr, bins)
    return dict(bins=bins, floors=floors, data=data, ydata=ydata, avail=avail,
                K_av=int(sum(len(a) for a in avail)),
                PG=arr["P"] @ arr["G"].T, seed=seed, noise_seed=noise_seed,
                gamma=gamma, K_seed=K_seed, null=null)


def null_twin(U, arr):
    """
    The PAIRED pure-noise null: the same sky, the same array, the same floors,
    the same noise draws -- with the signal removed from BOTH blocks. Because
    inject() and pterm_block() draw their noise before adding any source, the
    null twin's data is bit-for-bit the noise part of the signal arm's, so a
    phantom in the null arm and a phantom in the signal arm are the same noise.
    """
    rng = np.random.default_rng(U["noise_seed"])
    for B in U["bins"]:
        V25.floor_for(rng, arr, B)
    data = [W.inject(rng, B, null=True) for B in U["bins"]]
    yrng = np.random.default_rng(U["noise_seed"] + 500000)
    ydata = [pterm_block(yrng, arr, B, null=True) for B in U["bins"]]
    Un = dict(U)
    Un["data"], Un["ydata"], Un["null"] = data, ydata, True
    return Un


# ====================================================== joint Earth-block fit
def fstat(r, B):
    """
    VERNIER-16's F-statistic over the whole sky grid, with the three template
    Gram sums cached on the bin.

    They depend only on the whitening, not on the data, and recomputing them on
    every call was ~60 % of the joint fit's cost (the M-scan makes ~250 calls
    per bin and the null calibration makes ~4000). The arithmetic is V16's,
    unchanged and in the same order; the control ladder still calls V16's own
    function so that its bit-equality check means what it says.
    """
    m = B.get("_fs")
    if m is None:
        Fp, Fc = B["Fpw"], B["Fcw"]
        Muu = np.sum(Fp * Fp, axis=0)
        Muv = np.sum(Fp * Fc, axis=0)
        Mvv = np.sum(Fc * Fc, axis=0)
        det = Muu * Mvv - Muv ** 2
        det = np.where(np.abs(det) < 1e-12, 1e-12, det)
        # contiguous REAL transposes: `Fpw.T @ r` with complex r promotes the
        # whole (6000, 116) template block to complex on every call, which was
        # 17 ms a call and ~90 % of the joint fit. Two real mat-vecs are 1 ms.
        m = (Muu, Muv, Mvv, det, np.ascontiguousarray(Fp.T),
             np.ascontiguousarray(Fc.T))
        B["_fs"] = m
    Muu, Muv, Mvv, det, FpT, FcT = m
    rr, ri = r.real, r.imag
    bu = FpT @ rr + 1j * (FpT @ ri)
    bv = FcT @ rr + 1j * (FcT @ ri)
    return (Mvv * np.abs(bu) ** 2 - 2 * Muv * np.real(np.conj(bu) * bv)
            + Muu * np.abs(bv) ** 2) / det


def _design(B, pos):
    """Real design matrix of the Earth block for grid positions `pos`."""
    return np.concatenate([np.stack([B["Fpw"][:, j], B["Fcw"][:, j]], 1)
                           for j in pos], axis=1)


def joint_solve(d, B, pos):
    """
    Simultaneous least-squares fit of all sources at `pos`.

    Returns (coef, F_joint, residual). F_joint = b^H M^-1 b is the maximised
    2 x log-likelihood gain, the same quantity and the same units as V16's
    F-statistic at a single position.
    """
    if not pos:
        return np.zeros(0, complex), 0.0, d
    A = _design(B, pos)
    Mm = A.T @ A
    b = A.T @ d
    try:
        coef = np.linalg.solve(Mm + RIDGE * np.eye(Mm.shape[0]) * np.trace(Mm), b)
    except np.linalg.LinAlgError:
        return np.zeros(2 * len(pos), complex), 0.0, d
    return coef, float(np.real(np.vdot(b, coef))), d - A @ coef


def _blocked(G, pos, i, sep=GUARD_SEP):
    """Grid points too close to any CURRENT position other than source i."""
    if len(pos) < 2:
        return None
    others = [pos[q] for q in range(len(pos)) if q != i]
    ca = G @ G[others].T
    return (ca > np.cos(sep)).any(axis=1)


def sage_refine(d, B, pos, n_sweep=N_SWEEP):
    """
    Space-alternating conditional maximisation: re-solve every source's sky
    position against the residual of all the others (whose amplitudes are
    re-fitted each time). This is what makes the fit JOINT rather than greedy --
    a component that greedy dropped between two sources can migrate onto one of
    them once the other is being modelled.
    """
    pos = list(pos)
    G = B["G"]
    for _ in range(n_sweep):
        moved = False
        for i in range(len(pos)):
            rest = [pos[q] for q in range(len(pos)) if q != i]
            _, _, r = joint_solve(d, B, rest)
            Fv = fstat(r, B)
            bl = _blocked(G, pos, i)
            if bl is not None:
                Fv = np.where(bl, -np.inf, Fv)
            j = int(np.argmax(Fv))
            if j != pos[i]:
                pos[i] = j
                moved = True
        if not moved:
            break
    return pos


def mscan_earth(d, B, m_max=M_MAX, n_sweep=N_SWEEP):
    """
    Nested M = 1..m_max joint fits. Each M starts from the previous solution
    plus the peak of its residual, then refines the WHOLE set.

    Returns (list of position lists, array of F_joint).
    """
    pos, POS, FJ = [], [], []
    for m in range(m_max):
        _, _, r = joint_solve(d, B, pos)
        Fv = fstat(r, B)
        if pos:
            ca = B["G"] @ B["G"][pos].T
            Fv = np.where((ca > np.cos(GUARD_SEP)).any(axis=1), -np.inf, Fv)
        pos = pos + [int(np.argmax(Fv))]
        pos = sage_refine(d, B, pos, n_sweep)
        _, F, _ = joint_solve(d, B, pos)
        POS.append(list(pos))
        FJ.append(F)
    return POS, np.array(FJ)


def recut_curve(out, rng):
    """
    Column-wise VERNIER-24 recut floor (Gumbel where the null zero-fraction is
    low, empirical q95 with a bootstrap error above the gate) of an
    (n_null, m_max) table of null statistics. The same estimator the detection
    and certification floors use everywhere else in the series, so the model-
    selection penalty is set by the same rule as every other threshold here.
    """
    return np.array([sum(W.recut_floor(out[:, m], rng).get(k, 0.0)
                         for k in ("floor", "err"))
                     for m in range(out.shape[1])])


def null_curve_earth(rng, B, m_max=M_MAX, n_null=NJ_NULL, n_sweep=N_SWEEP):
    """
    The model-selection penalty, MEASURED: the recut floor of the INCREMENT
    F_joint(m) - F_joint(m-1) over pure-noise realisations of this bin. A
    component is admitted iff its own increment clears its own null floor, which
    is a calibrated 5 %-per-component test rather than an information criterion
    asserted by me -- and is the same shape of rule as R3's Delta T floor.
    """
    out = np.zeros((n_null, m_max))
    n = B["Tn"].shape[0]
    for q in range(n_null):
        dn = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
        _, FJ = mscan_earth(dn, B, m_max, n_sweep)
        out[q] = np.diff(np.concatenate([[0.0], FJ]))
    return recut_curve(out, rng)


def select_M(FJ, bar):
    """Largest m whose joint statistic clears its own null floor (0 if none)."""
    ok = np.where(np.asarray(FJ) > np.asarray(bar))[0]
    return int(ok[-1] + 1) if ok.size else 0


# ================================================== pulsar block, marginalised
def _lgrid(Lstar):
    """EVEN-count distance grid centred on the prior mean -- V24's sub-fringe rule."""
    return Lstar + np.linspace(-3 * PRIOR_PC, 3 * PRIOR_PC, NPTS_L)


class BandCtx:
    """
    The pulsar block of the WHOLE BAND, with the pulsar distances marginalised
    ONCE over their 10 pc priors.

    This is the point of R3 and it is why the marginalisation cannot be done
    bin by bin: pulsar a has ONE distance, and every source in every bin sees
    the same L_a. Ten bins x M sources therefore impose 10M phase constraints
    on one unknown per pulsar. Marginalising L_a separately in each bin throws
    that lock away and leaves each bin's comb phase almost free -- which is a
    model that can eat noise instead of a vernier that can resolve blends. The
    per-bin version was tried first and behaved exactly like that (a single
    source's ln Z was as large 1.5 rad from the truth as on it).

        ll_a(L)  = sum_b [ 2 Re(y_ab^* yhat_ab(L)) - |yhat_ab(L)|^2 ]
        ln Z     = sum_a ln \\int dL pi(L) exp(ll_a(L))

    `ll` is kept per bin so that changing one component touches one bin's term
    and the band total is an array add.
    """

    def __init__(self, arr, U, cache=128):
        self.arr, self.U = arr, U
        self.bins = U["bins"]
        self.nb = len(self.bins)
        n = arr["n_psr"]
        self.n = n
        Lst = arr["L"]
        self.Lg = np.stack([_lgrid(Lst[a]) for a in range(n)])
        self.npts = self.Lg.shape[1]
        self.w = np.exp(-0.5 * ((self.Lg - Lst[:, None]) / PRIOR_PC) ** 2)
        self.lognorm = np.log(np.trapezoid(self.w, self.Lg, axis=1))
        self.tc = [float(V.t_merge(float(B["f"]), MC_FIX)) for B in self.bins]
        self.ll = [np.zeros((n, self.npts)) for _ in range(self.nb)]
        self.acc = np.zeros((n, self.npts))
        self._c = {}
        self._cap = cache

    # ---- comb kernels (data independent, so they survive across realisations)
    def kernel(self, b, j):
        k = self._c.get((b, j))
        if k is None:
            omc = np.clip(1.0 - self.U["PG"][:, j], 5e-2, None)
            ph, env = _abs_phase(self.Lg, omc[:, None],
                                 float(self.bins[b]["f"]), self.tc[b])
            k = env * np.exp(1j * ph)
            if len(self._c) >= self._cap:
                self._c.pop(next(iter(self._c)))
            self._c[(b, j)] = k
        return k

    def yhat(self, b, pos, coef):
        B = self.bins[b]
        out = np.zeros((self.n, self.npts), complex)
        for m, j in enumerate(pos):
            E = B["Fpw"][:, j] * coef[2 * m] + B["Fcw"][:, j] * coef[2 * m + 1]
            out += R_AMP * E[:, None] * self.kernel(b, j)
        return out

    def ll_of(self, b, y, yh):
        return 2.0 * np.real(np.conj(y)[:, None] * yh) - np.abs(yh) ** 2

    def _marg(self, acc):
        mx = acc.max(axis=1)
        z = np.trapezoid(self.w * np.exp(acc - mx[:, None]), self.Lg, axis=1)
        return float(np.sum(mx + np.log(np.maximum(z, 1e-300)) - self.lognorm))

    # ---- state
    def reset(self):
        for b in range(self.nb):
            self.ll[b][:] = 0.0
        self.acc[:] = 0.0

    def set_bin(self, b, y, pos, coef):
        new = (self.ll_of(b, y, self.yhat(b, pos, coef)) if pos
               else np.zeros((self.n, self.npts)))
        self.acc += new - self.ll[b]
        self.ll[b] = new

    def lnZ(self):
        return self._marg(self.acc)

    def lnZ_if(self, b, ll_b):
        """ln Z if bin b's term were `ll_b` instead of its current one."""
        return self._marg(self.acc + ll_b - self.ll[b])


def spread_candidates(Fv, G, n=N_CAND, sep=CAND_SEP):
    """
    n candidate positions separated by at least `sep`: the F peak, then the best
    point outside sep of everything already taken, and so on.

    The top-n grid points of an F surface are all one blob; a blended component
    has to be offered positions FAR from where greedy put it, which is what this
    returns.
    """
    order = np.argsort(Fv)[::-1]
    take, keep = [], np.ones(G.shape[0], bool)
    cs = np.cos(sep)
    for j in order:
        if not keep[j]:
            continue
        take.append(int(j))
        if len(take) >= n:
            break
        keep &= (G @ G[j]) <= cs
    return take


# ------------------------------------------------------- band-level R3 engine
class R3State:
    """
    The band-wide model under R3: positions per bin, their joint Earth fits, and
    the shared marginalised pulsar block.

        T = sum_b F_joint,b  +  ln Z(all components)

    Every operation below (move / add / drop) is scored by the change in T.
    """

    def __init__(self, U, ctx):
        self.U, self.ctx = U, ctx
        self.pos = {b: [] for b in range(ctx.nb)}
        self.F = np.zeros(ctx.nb)
        ctx.reset()

    def set(self, b, pos):
        d, y, B = self.U["data"][b], self.U["ydata"][b], self.U["bins"][b]
        coef, F, _ = joint_solve(d, B, pos)
        self.pos[b] = list(pos)
        self.F[b] = F
        self.ctx.set_bin(b, y, pos, coef)

    def T(self):
        return float(self.F.sum()) + self.ctx.lnZ()

    def trial(self, b, pos):
        """T if bin b held `pos`, without committing to it."""
        d, y, B = self.U["data"][b], self.U["ydata"][b], self.U["bins"][b]
        coef, F, _ = joint_solve(d, B, pos)
        ll = (self.ctx.ll_of(b, y, self.ctx.yhat(b, pos, coef)) if pos
              else np.zeros((self.ctx.n, self.ctx.npts)))
        return float(self.F.sum() - self.F[b] + F) + self.ctx.lnZ_if(b, ll)

    def cand_scan(self, b, rest, cands):
        """
        Best of `cands` as the one component to put in bin b beside `rest`.

        The others' comb contribution is built ONCE and the candidate only adds
        its own term -- the naive version rebuilt every component's kernel for
        every candidate and cost 6x more. The candidate's amplitude is its
        conditional least-squares fit against the others' residual (the same
        SAGE convention R2's refinement uses); the committing `set` call then
        does the exact joint solve.
        """
        d, y, B = self.U["data"][b], self.U["ydata"][b], self.U["bins"][b]
        crest, Frest, r = joint_solve(d, B, rest)
        yh_rest = (self.ctx.yhat(b, rest, crest) if rest
                   else np.zeros((self.ctx.n, self.ctx.npts), complex))
        Fbase = float(self.F.sum() - self.F[b] + Frest)
        best, bj = -np.inf, None
        for j in cands:
            a = np.stack([B["Fpw"][:, j], B["Fcw"][:, j]], 1)
            Mm = a.T @ a
            bj_ = a.T @ r
            cj = np.linalg.solve(Mm + RIDGE * np.eye(2) * np.trace(Mm), bj_)
            F = Fbase + float(np.real(np.vdot(bj_, cj)))
            E = B["Fpw"][:, j] * cj[0] + B["Fcw"][:, j] * cj[1]
            yh = yh_rest + R_AMP * E[:, None] * self.ctx.kernel(b, j)
            T = F + self.ctx.lnZ_if(b, self.ctx.ll_of(b, y, yh))
            if T > best:
                best, bj = T, int(j)
        return bj, best, r

    # ---- the three moves
    def refine(self, n_sweep=N_SWEEP_P, n_cand=N_CAND):
        """Move each component to the candidate position maximising T."""
        for _ in range(n_sweep):
            moved = False
            for b in range(self.ctx.nb):
                for i in range(len(self.pos[b])):
                    B = self.U["bins"][b]
                    rest = [self.pos[b][q] for q in range(len(self.pos[b]))
                            if q != i]
                    _, _, r = joint_solve(self.U["data"][b], B, rest)
                    Fv = fstat(r, B)
                    bl = _blocked(B["G"], self.pos[b], i)
                    if bl is not None:
                        Fv = np.where(bl, -np.inf, Fv)
                    cands = spread_candidates(Fv, B["G"], n_cand)
                    cur = self.pos[b][i]
                    if cur not in cands:
                        cands.append(cur)
                    bj, _, _ = self.cand_scan(b, rest, cands)
                    if bj is not None and bj != cur:
                        self.set(b, rest + [bj])
                        moved = True
            if not moved:
                break

    def prune(self, thr):
        """
        Drop any component whose contribution to T is below the null floor for
        a model of that size. `thr` is the per-step floor curve.
        """
        dropped = 0
        changed = True
        while changed:
            changed = False
            for b in range(self.ctx.nb):
                for i in range(len(self.pos[b])):
                    rest = [self.pos[b][q] for q in range(len(self.pos[b]))
                            if q != i]
                    bar = thr[min(len(self.pos[b]) - 1, len(thr) - 1)]
                    if self.T() - self.trial(b, rest) < bar:
                        self.set(b, rest)
                        dropped += 1
                        changed = True
                        break
                if changed:
                    break
        return dropped

    def grow(self, thr, m_max=M_MAX, n_cand=N_CAND):
        """
        Offer each bin one more component at a time; keep it only if its own
        Delta T clears the null floor. This is where the pulsar block can
        support a component the Earth block alone could not -- or fail to.
        """
        added = 0
        for b in range(self.ctx.nb):
            while len(self.pos[b]) < m_max:
                B = self.U["bins"][b]
                _, _, r = joint_solve(self.U["data"][b], B, self.pos[b])
                Fv = fstat(r, B)
                if self.pos[b]:
                    ca = B["G"] @ B["G"][self.pos[b]].T
                    Fv = np.where((ca > np.cos(GUARD_SEP)).any(axis=1), -np.inf, Fv)
                cands = spread_candidates(Fv, B["G"], max(n_cand // 4, 4))
                base = self.T()
                bj, best, _ = self.cand_scan(b, self.pos[b], cands)
                bar = thr[min(len(self.pos[b]), len(thr) - 1)]
                if bj is None or best - base < bar:
                    break
                self.set(b, self.pos[b] + [bj])
                added += 1
        return added


def null_dT(rng, U, arr, ctx, n_null=NJ_NULL_P, m_max=M_MAX, n_cand=N_CAND):
    """
    The R3 admission threshold, MEASURED: the distribution of Delta T for a
    component fitted to PURE NOISE in both blocks, over the same greedy adds the
    real fit makes. The recut floor of that distribution is what a component has
    to clear to be admitted -- so R3's model selection is calibrated by the same
    estimator as every other threshold in the series, and the over-fitting bias
    the marginalisation itself carries is inside the calibration rather than
    outside it.

    PER STEP, and ROUND-ROBIN, and why both matter. Delta T for one more
    component depends on (a) how many components the bin already holds -- by
    m = 12 on pure noise the bin is saturated and every further add is worth
    almost nothing -- and (b) how well L is already pinned, and L is shared
    across the WHOLE BAND.

    Getting either wrong sets a floor the real fit walks straight through, and
    both earlier versions of this function did: pooling every step into one
    scalar floor let the many tiny high-m increments drag it down, and the
    pure-noise arm then claimed ~7 components a run against R1/R2's ~0.3. So
    the floor is a CURVE in m, exactly like R2's, and the adds go round the
    band -- every bin gets its m-th component before any bin gets its
    (m+1)-th -- so the band is populated as the threshold is measured.
    """
    vals = np.zeros((n_null * ctx.nb, m_max))
    Un = dict(U)
    n = arr["n_psr"]
    for q in range(n_null):
        Un["data"] = [(rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
                      for _ in U["bins"]]
        Un["ydata"] = [(rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
                       for _ in U["bins"]]
        st = R3State(Un, ctx)
        for b in range(ctx.nb):
            st.set(b, [])
        for m in range(m_max):
            for b in range(ctx.nb):
                B = Un["bins"][b]
                _, _, r = joint_solve(Un["data"][b], B, st.pos[b])
                Fv = fstat(r, B)
                if st.pos[b]:
                    ca = B["G"] @ B["G"][st.pos[b]].T
                    Fv = np.where((ca > np.cos(GUARD_SEP)).any(axis=1),
                                  -np.inf, Fv)
                cands = spread_candidates(Fv, B["G"], max(n_cand // 4, 4))
                base = st.T()
                bj, best, _ = st.cand_scan(b, st.pos[b], cands)
                vals[q * ctx.nb + b, m] = best - base
                st.set(b, st.pos[b] + [bj])
    ctx.reset()
    return recut_curve(vals, rng)


# =============================================================== bin drivers
def fit_R2(d, B, dbar, m_max=M_MAX, n_sweep=N_SWEEP):
    """
    One bin under ladder R2.

    Stepwise admission on the joint statistic: component m is admitted iff its
    own increment F_joint(m) - F_joint(m-1) clears that increment's null floor.
    Matched, term for term, with R3's Delta T rule, so "does the pulsar block
    change M_hat" is a comparison between two thresholds of the same shape.
    """
    POS, FJ = mscan_earth(d, B, m_max, n_sweep)
    dF = np.diff(np.concatenate([[0.0], FJ]))
    m = 0
    while m < m_max and dF[m] > dbar[m]:
        m += 1
    return (POS[m - 1] if m else []), m, FJ, POS


# ================================================================= accounting
def account(B, pos, avail, claimed=None):
    """
    VERNIER-18 genuine / split / blend accounting, one rule for all three
    ladders and no credit for fitting one phantom to two truths.

    Every fitted component is matched to its nearest AVAILABLE true source. It
    is GENUINE if that source is inside TOL_RAD and has not already been claimed
    by another component of this fit; SPLIT if it is inside TOL_RAD but the
    source is already claimed; FAINT if it lands on a true source that was not
    available (a real object, but not one of the ones we are converting);
    PHANTOM if it lands on nothing at all.

    Returns (list of (kind, true_index_or_-1), claimed_set).
    """
    claimed = set() if claimed is None else set(claimed)
    av = sorted(avail)
    out = []
    for j in pos:
        kind, which = "phantom", -1
        if av:
            ca = B["Gt"][av] @ B["G"][j]
            t = int(np.argmax(ca))
            if np.arccos(np.clip(ca[t], -1, 1)) <= TOL_RAD:
                i = av[t]
                if i in claimed:
                    kind = "split"
                else:
                    kind, which = "genuine", i
                    claimed.add(i)
        if kind == "phantom":
            ca = B["Gt"] @ B["G"][j]
            t = int(np.argmax(ca))
            if np.arccos(np.clip(ca[t], -1, 1)) <= TOL_RAD:
                kind = "faint"
        out.append((kind, which))
    return out, claimed


def chance_floor(n_pool):
    """Probability a random sky point lands within TOL_RAD of one of n_pool sources."""
    p1 = (1.0 - np.cos(TOL_RAD)) / 2.0
    return 1.0 - (1.0 - p1) ** max(int(n_pool), 0)


# ============================================================ certification
def certify_pass(rng, arr, U, rec, comb_cap=COMB_CAP, soft=True,
                 npts=NPTS_CERT, n_null=N_NULL_C):
    """
    V24's certification, unmodified, run over a ladder's recovered set.

    `rec` is a list of (bin index, true source index) in the order the ladder
    would have banked them (loudest first). Each addition is certified against
    the accumulated comb, and the SOURCE-LEVEL SCRAMBLED NULL is run at every
    step -- the brief's "scrambled null on every joint fit".

    The confusion budget of a certified source excludes the OTHER recovered
    sources in its own bin, because a joint fit models them simultaneously. For
    R1 that set is the ones already subtracted, which is V25/V27's rule; for
    R2/R3 it is all of the bin's recoveries. This is a real part of what joint
    fitting buys and it is stated rather than hidden.
    """
    bins = U["bins"]
    per_bin = {}
    for (b, k) in rec:
        per_bin.setdefault(b, []).append(k)
    det, hist = [], []
    prev_cert = np.zeros(arr["n_psr"], bool)
    seen = {}
    for (b, k) in rec:
        B = bins[b]
        sub = set(per_bin[b]) - {k}
        det.append((B, k, sub))
        others = sorted(det[:-1], key=lambda t: -float(t[0]["R2"][t[1]]))
        use = others[:comb_cap - 1] + [det[-1]]
        C = W.certify(rng, arr, use, soft=soft, npts=npts, n_null=n_null)
        Cn = W.certify(rng, arr, use, scramble=True, soft=soft, npts=npts,
                       n_null=n_null)
        new_cert = int(np.sum(C["cert"] & ~prev_cert))
        prev_cert = prev_cert | C["cert"]
        seen[b] = seen.get(b, 0) + 1
        hist.append(dict(bin=b, which=k, n_combs=len(det),
                         n_cert=C["n_cert"], n_true=C["n_true"],
                         n_wrong=C["n_cert"] - C["n_true"],
                         n_cert_null=Cn["n_cert"], new_cert=new_cert,
                         lnK=C["lnK"], S_max=C["S_max"], S_med=C["S_med"],
                         rho_arr=C["rho_arr"], rho_P_p90=C["rho_P_p90"],
                         f_soft=C["f_soft"], floor_inc=C["floor_inc"]))
    return hist


def rec_order(U, per_bin):
    """
    Bank order for the certification pass: recovered sources loudest first,
    which is the order any real pipeline would carry them in.
    `per_bin` maps bin index -> list of true source indices.
    """
    rec = [(b, k) for b, ks in per_bin.items() for k in ks]
    return sorted(rec, key=lambda t: -float(U["bins"][t[0]]["R2"][t[1]]))


# ==================================================================== ladders
def ladder_R1(U, arr, soft=True, comb_cap=COMB_CAP, max_turn=MAX_TURN,
              npts=NPTS_CERT, n_null=N_NULL_C):
    """
    THE CONTROL. vernier27_lib.chain, copied verbatim, with one addition: the
    grid index of every component is recorded so that account() can score all
    three ladders by one rule. Nothing about the search, the floors, the
    subtraction, the certification or the deepening logic is changed --
    `validate_R1` checks that against the imported original.
    """
    bins, floors, data = U["bins"], U["floors"], U["data"]
    rng = np.random.default_rng(U["noise_seed"])
    # the imported chain draws floors and data from this stream before the loop;
    # replay those draws so the search's own randomness matches turn for turn.
    for B in bins:
        V25.floor_for(rng, arr, B)
    for B in bins:
        W.inject(rng, B, null=U["null"])
    found = [[] for _ in bins]
    sub = [set() for _ in bins]
    det, hist = [], []
    gen = spl = blend = 0
    prevF = [None] * len(bins)
    F0 = [None] * len(bins)
    prev_cert = np.zeros(arr["n_psr"], bool)
    for t in range(max_turn):
        best, Fall = None, [None] * len(bins)
        for b, B in enumerate(bins):
            d = data[b]
            if found[b]:
                A = _design(B, found[b])
                Mm = A.T @ A
                try:
                    coef = np.linalg.solve(
                        Mm + RIDGE * np.eye(Mm.shape[0]) * np.trace(Mm), A.T @ d)
                    r = d - A @ coef
                except np.linalg.LinAlgError:
                    continue
            else:
                r = d
            Fv = S.fstat(r, B["Fpw"], B["Fcw"])
            Fall[b] = Fv
            j = int(np.argmax(Fv))
            exc = float(Fv[j]) - (floors[b]["floor"] + floors[b]["err"])
            if exc > 0 and (best is None or exc > best[0]):
                best = (exc, b, j, float(Fv[j]))
        if best is None:
            break
        if F0[0] is None:
            F0 = [None if x is None else x.copy() for x in Fall]
        exc, b, j, Fmax = best
        B = bins[b]
        bar_b = floors[b]["floor"] + floors[b]["err"]
        newly = bool(prevF[b] is not None and prevF[b][j] <= bar_b)
        newly0 = bool(F0[b] is not None and F0[b][j] <= bar_b)
        which, kind = -1, "blend"
        ca = B["Gt"] @ B["G"][j]
        i = int(np.argmax(ca))
        n_true_pool = B["n_seed"] + W.TOP_TRUE
        if np.arccos(np.clip(ca[i], -1, 1)) <= TOL_RAD and i < n_true_pool:
            if i in sub[b]:
                kind, spl = "split", spl + 1
            else:
                kind, which, gen = "genuine", i, gen + 1
        else:
            blend += 1
        found[b].append(j)
        prevF = Fall
        if which < 0:
            hist.append(dict(turn=t + 1, bin=b, grid=j, kind=kind, which=-1,
                             newly0=newly0, deepen=False, new_cert=0,
                             n_cert=int(prev_cert.sum()), n_true=0, n_wrong=0,
                             n_cert_null=0, lnK=np.nan, S_max=0.0))
            continue
        det.append((B, which, set(sub[b])))
        others = sorted(det[:-1], key=lambda t_: -float(t_[0]["R2"][t_[1]]))
        use = others[:comb_cap - 1] + [det[-1]]
        C = W.certify(rng, arr, use, soft=soft, npts=npts, n_null=n_null)
        new_cert = int(np.sum(C["cert"] & ~prev_cert))
        if C["n_cert"] > 0 or newly:
            Cn = W.certify(rng, arr, use, scramble=True, soft=soft, npts=npts,
                           n_null=n_null)
        else:
            Cn = dict(n_cert=0)
        sub[b].add(which)
        prev_cert = prev_cert | C["cert"]
        deepen = bool((newly or new_cert > 0) and Cn["n_cert"] == 0)
        hist.append(dict(turn=t + 1, bin=b, grid=j, kind=kind, which=which,
                         newly0=newly0, deepen=deepen, new_cert=new_cert,
                         n_cert=C["n_cert"], n_true=C["n_true"],
                         n_wrong=C["n_cert"] - C["n_true"],
                         n_cert_null=Cn["n_cert"], lnK=C["lnK"],
                         S_max=C["S_max"]))
    dp = [h["deepen"] for h in hist]
    per_bin_pos = {}
    for h in hist:
        per_bin_pos.setdefault(h["bin"], []).append(h["grid"])
    return dict(hist=hist, pos=per_bin_pos, turns=len(hist),
                native_genuine=gen, native_split=spl, native_blend=blend,
                ignite=any(dp[i] and dp[i + 1] for i in range(len(dp) - 1)),
                n_cert_tot=sum(h["new_cert"] for h in hist),
                n_wrong_tot=sum(h["n_wrong"] for h in hist),
                n_null_tot=sum(h["n_cert_null"] for h in hist))


def ladder_R2(U, arr, bars, m_max=M_MAX, n_sweep=N_SWEEP):
    """Joint-bin ladder: every bin fitted simultaneously, M chosen per bin."""
    pos, Mhat, POSN = {}, {}, {}
    for b, B in enumerate(U["bins"]):
        p, m, FJ, POS = fit_R2(U["data"][b], B, bars[b], m_max, n_sweep)
        if p:
            pos[b] = p
        Mhat[b] = m
        POSN[b] = POS
    return dict(pos=pos, Mhat=Mhat, POS=POSN)


def ladder_R3(U, arr, ctx, R2, thrT, m_max=M_MAX, n_cand=N_CAND):
    """
    Joint-bin ladder WITH the marginalised pulsar block, band-wide.

    Starts from R2's selected model -- the two ladders therefore differ only in
    what the pulsar block does to it -- then, scoring every move by the change
    in T = sum_b F_joint,b + ln Z:

      REFINE  every component may move to any of N_CAND well-separated candidate
              positions. This is the blend-resolution step: a component that
              greedy dropped between two sources can be pulled onto one of them
              by the comb, which only locks at a true geometry.
      PRUNE   any component not earning the null-calibrated admission floor is
              dropped. This is the anti-phantom step.
      GROW    each bin is offered further components, each admitted only on the
              same floor. This is where the surplus either buys a detection the
              Earth block could not make, or does not.
    """
    st = R3State(U, ctx)
    for b in range(ctx.nb):
        st.set(b, R2["pos"].get(b, []))
    st.refine(N_SWEEP_P, n_cand)
    n_drop = st.prune(thrT)
    n_add = st.grow(thrT, m_max, n_cand)
    if n_add:
        st.refine(1, n_cand)
    pos = {b: p for b, p in st.pos.items() if p}
    return dict(pos=pos, Mhat={b: len(p) for b, p in st.pos.items()},
                n_drop=n_drop, n_add=n_add, T=st.T())


def score(U, pos):
    """
    Uniform V18 accounting of a ladder's fitted positions, bin by bin.

    Returns counts and the per-bin genuine true-source indices, which are what
    the certification pass and the conversion efficiency are both built on.
    """
    n = dict(genuine=0, split=0, blend=0, faint=0, phantom=0, claims=0)
    per_bin, perbin_kinds = {}, {}
    for b, ps in pos.items():
        kinds, _ = account(U["bins"][b], ps, U["avail"][b])
        perbin_kinds[b] = kinds
        for kind, which in kinds:
            n["claims"] += 1
            n[kind if kind != "phantom" else "phantom"] += 1
            if kind == "genuine":
                per_bin.setdefault(b, []).append(which)
    n["blend"] = n["faint"] + n["phantom"]      # V18 "blend" = anything not a
    return n, per_bin, perbin_kinds             # clean match to an available source


def end_state(hist, comb_cap=COMB_CAP):
    """End-state fringe contrast, its trials bar and the accumulated comb count."""
    if not hist:
        return dict(S_end=0.0, bar=np.nan, combs=0.0, cross=False,
                    cert=0, true=0, wrong=0, null=0)
    S = max(h["S_max"] for h in hist)
    bar = float(np.median([h["lnK"] for h in hist]))
    return dict(S_end=float(S), bar=bar, combs=float(min(len(hist), comb_cap)),
                cross=bool(S >= bar),
                cert=int(sum(h["new_cert"] for h in hist)),
                true=int(sum(h["n_true"] for h in hist)),
                wrong=int(sum(h["n_wrong"] for h in hist)),
                null=int(sum(h["n_cert_null"] for h in hist)))
