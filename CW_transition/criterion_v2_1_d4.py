"""criterion-v2.1 / D4 — the REALISATION-LEVEL DISCORDANCE GATE.

Machinery: the D3 co-registration statistic, re-implemented on cronus from banked arrays
(`reports/ig2_levers.npz` + the per-realisation banks) and GATED bit-comparably against
IGNITE-2's banked per-candidate values (`reports/ig2_purity.npz`: sig_R_all / sig_R_det,
809 candidates; null_R_all, 280 null detections).  See `d4_gate()`.

D3 vetoed PER PULSAR and failed its pre-registration on (b): above onset the leave-one-out
reference u_R is itself poisoned, so a TRUE cert fails concordance with its own reference.
D3's one perfect number was (c): 42/42 wrong-counterpart DETECTIONS rejected.  D4 promotes
that to the REALISATION level -- flag the whole realisation as WRONG-COUNTERPART and veto
every certification in it -- and asks whether an aggregate exists that separates the two
populations where the per-pulsar form could not.

Definitions (all in the source parameters theta = (log10 f, log10 mc), levers from SIREN):

    dk_a    = mapk_a - n_true_grid_a                    fringe error of pulsar a
    u_a     = 2*pi*dk_a * g_a / |g_a|^2                 min-norm source shift a demands
    s2_a    = (2*pi)^2 * [ 1/6 + (1-qmax_a)*(sig_EM_a/dL_a)^2 ]     phase variance
    I_a     = g_a g_a^T / s2_a                          (rank-1) information
    R_a     = (u_a-u_R)^T (Sigma_a + Sigma_R)^-1 (u_a-u_R)          leave-one-out, chi2(dof)

Realisation-level aggregates over the DETECTED set (D):

    max_R / min_R / mean_R / frac_R   of R_all over a in D
    S_ref(R)  = J^T I^-1 J   with  I = sum_b I_b,  J = sum_b g_b (2*pi*dk_b)/s2_b
                the significance of the reference solution's DISPLACEMENT from the assumed
                counterpart; algebraically chi2(u=0) - min_u chi2(u), i.e. the GLRT for
                "the array's consensus source is not the one we targeted".   ~ chi2(2)

CAVEAT INHERITED FROM D3 (stated, not resolved): dk uses the true fringe index
n_true_grid.  An implementable form replaces k* by the EM-prior fringe and inherits an
extra per-pulsar variance (sig_EM/dL)^2 that s2_a budgets only through its (1-qmax)
mixture term.  Every D4 number is therefore an UPPER BOUND on the implementable form's
power, exactly as every D3 number was.
"""
import numpy as np

CRIT2 = 9.210340371976182   # chi2 crit, dof 2, p = 0.01  (inherited from D3, NOT retuned)
CRIT1 = 6.634896601021217   # chi2 crit, dof 1, p = 0.01
NPSR = 116


class Levers:
    """Analytic lag-levers + fringe spacings, from IGNITE-2's banked `ig2_levers.npz`."""

    def __init__(self, path='reports/ig2_levers.npz'):
        L = np.load(path, allow_pickle=True)
        self.names = L['names']
        self.skies = L['skies']
        self.T_grid = L['T_grid']
        self.g_lev = L['g_lev']          # (T, sky, 116, 2)
        self.dL = L['dL']                # (sky, 116)
        self.sig_lit = L['sig_lit']
        self.sig_vlbi = L['sig_vlbi']

    def get(self, T, geo, tier):
        ti = int(np.where(self.T_grid == int(T))[0][0])
        si = int(np.where(self.skies == int(geo))[0][0])
        g = self.g_lev[ti, si]                                  # (116,2)
        sg = self.sig_lit if tier == 'lit' else self.sig_vlbi
        r2 = (sg / self.dL[si]) ** 2                            # EM prior width in fringes^2
        return g, r2


def phase_var(qmax, r2):
    """s2_a: within-fringe U(+-1/2) + MAP-snap variance + two-point fringe-posterior mixture."""
    return (2 * np.pi) ** 2 * (1.0 / 6.0 + (1.0 - qmax) * r2)


def _pinv_sym(M, rcond=1e-10):
    w, V = np.linalg.eigh(M)
    keep = w > rcond * max(w.max(), 1e-300)
    winv = np.zeros_like(w)
    winv[keep] = 1.0 / w[keep]
    return (V * winv) @ V.T, int(keep.sum())


def _basis(M, rcond=1e-10):
    """Eigenvectors spanning the range of symmetric PSD M, and the rank."""
    w, V = np.linalg.eigh(M)
    keep = w > rcond * max(w.max(), 1e-300)
    return V[:, keep], int(keep.sum())


def _quad_in(S, d, Vk):
    """d^T S^-1 d evaluated INSIDE the reference's information subspace Vk.

    D3's rank handling: dof = rank of the REFERENCE set's information (rank-1 arises when
    the reference is a single pulsar -- the R_det control).  A rank-2 reference reduces to
    the plain 2-d form, which is why R_all is unaffected.
    """
    if Vk.shape[1] == 0:
        return np.nan, 0
    Sp = Vk.T @ S @ Vk
    dp = Vk.T @ d
    w, V = np.linalg.eigh(Sp)
    keep = w > 1e-10 * max(w.max(), 1e-300)
    if not keep.any():
        return np.nan, 0
    q = V[:, keep].T @ dp
    return float(np.sum(q ** 2 / w[keep])), int(keep.sum())


def R_leave_one_out(dk, qmax, g, s2, ref_mask, cand_idx):
    """R_a for each candidate a, against the reference set ref_mask \\ {a}."""
    gg = np.sum(g * g, axis=1)
    ginv2 = 1.0 / gg
    u = (2 * np.pi) * dk[:, None] * g * ginv2[:, None]
    I_b = np.einsum('bi,bj->bij', g, g) / s2[:, None, None]
    J_b = g * ((2 * np.pi) * dk / s2)[:, None]
    R = np.full(NPSR, np.nan)
    dof = np.zeros(NPSR, int)
    nref = np.zeros(NPSR, int)
    for a in cand_idx:
        m = ref_mask.copy()
        m[a] = False
        n = int(m.sum())
        nref[a] = n
        if n == 0:
            continue
        IR = I_b[m].sum(0)
        SigR, rk = _pinv_sym(IR)
        if rk == 0:
            continue
        uR = SigR @ J_b[m].sum(0)
        Sa = np.outer(g[a], g[a]) * s2[a] * ginv2[a] ** 2
        Vk, _ = _basis(IR)
        R[a], dof[a] = _quad_in(Sa + SigR, u[a] - uR, Vk)
    return R, dof, nref


def S_ref(dk, g, s2, mask):
    """Displacement significance of the reference solution built from `mask`.

    S = J^T I^-1 J ~ chi2(rank).  Equals chi2(u=0) - min_u chi2(u): the GLRT statistic for
    'the pulsars in `mask` co-register at a source OTHER than the assumed counterpart'.
    Also returns the reduced residual of that best common solution (the co-registration
    goodness-of-fit: small => the array agrees on SOME source; large => no source at all
    explains the fringe field).
    """
    n = int(mask.sum())
    if n == 0:
        return np.nan, 0, np.nan
    gm, s2m, dkm = g[mask], s2[mask], dk[mask]
    I = np.einsum('bi,bj->bij', gm, gm) / s2m[:, None, None]
    I = I.sum(0)
    J = (gm * ((2 * np.pi) * dkm / s2m)[:, None]).sum(0)
    Iinv, rk = _pinv_sym(I)
    if rk == 0:
        return np.nan, 0, np.nan
    S = float(J @ Iinv @ J)
    chi2_0 = float(np.sum(((2 * np.pi) * dkm) ** 2 / s2m))
    dofr = max(n - rk, 1)
    return S, rk, (chi2_0 - S) / dofr


def realisation_stats(lev, dk, qmax, dlnL, lnK, tier, T, geo, floor):
    """Every realisation-level aggregate D4 considers, from one realisation's banked arrays."""
    g, r2 = lev.get(T, geo, tier)
    s2 = phase_var(np.asarray(qmax, float), r2)
    dk = np.asarray(dk, int)
    det = np.asarray(dlnL) > np.maximum(np.asarray(lnK), floor)
    cert = det & (np.asarray(qmax) > 0.9)
    d_idx = np.where(det)[0]
    R_all, dof_all, _ = R_leave_one_out(dk, qmax, g, s2, np.ones(NPSR, bool), d_idx)
    Rv = R_all[d_idx]
    Rv = Rv[np.isfinite(Rv)]
    S_a, _, gof_a = S_ref(dk, g, s2, np.ones(NPSR, bool))
    S_d, _, gof_d = S_ref(dk, g, s2, det)
    out = dict(
        n_det=int(det.sum()), n_cert=int(cert.sum()),
        n_cert_wrong=int((cert & (dk != 0)).sum()), ndk0=int((dk == 0).sum()),
        max_R=float(np.max(Rv)) if len(Rv) else np.nan,
        min_R=float(np.min(Rv)) if len(Rv) else np.nan,
        mean_R=float(np.mean(Rv)) if len(Rv) else np.nan,
        frac_R=float(np.mean(Rv > CRIT2)) if len(Rv) else np.nan,
        S_ref_all=S_a, gof_all=gof_a, S_ref_det=S_d, gof_det=gof_d,
        det=det, cert=cert, dk=dk, R_all=R_all,
    )
    return out


# ---------------------------------------------------------------- value gate
def d4_gate(verbose=True):
    """Reproduce IGNITE-2's banked per-candidate R_all / R_det. Must pass before any scoring."""
    lev = Levers()
    P = np.load('reports/ig2_purity.npz', allow_pickle=True)
    B = np.load('reports/ignite_bank.npz', allow_pickle=True)
    ok = True
    for tag, cellk, realk, ak, Rallk, Rdetk in [
        ('sig', 'sig_cell', 'sig_real', 'sig_a', 'sig_R_all', 'sig_R_det'),
        ('null', 'null_cell', 'null_real', 'null_a', 'null_R_all', 'null_R_det')]:
        rows = np.arange(len(P[realk]))
        e_all, e_det, dkbad = [], [], 0
        floors = {}
        for ci in range(24):
            h, T, _ = P['fN_cells'][ci]
            key = "(%s, %d, '%s')" % (repr(round(float(h), 2)), int(T), str(P['fN_tier'][ci]))
            floors[key] = float(P['fN'][ci])
        for r in np.unique(P[realk]):
            i = int(r)
            sel = rows[P[realk] == r]
            key = str(P[cellk][sel[0]])
            tier = str(B['tier'][i]); T = int(B['T_label'][i]); geo = int(B['geo_id'][i])
            g, r2 = lev.get(T, geo, tier)
            qm = B['qmax'][i]
            s2 = phase_var(qm, r2)
            dk = B['mapk'][i] - B['n_true_grid'][i]
            det = B['dlnL_det'][i] > np.maximum(B['lnK'][i], floors[key])
            cands = P[ak][sel]
            Ra, _, _ = R_leave_one_out(dk, qm, g, s2, np.ones(NPSR, bool), cands)
            Rd, _, _ = R_leave_one_out(dk, qm, g, s2, det, cands)
            for j, a in zip(sel, cands):
                dkbad += int(dk[a] != P[tag + '_dk'][j])
                for got, want, acc in [(Ra[a], P[Rallk][j], e_all), (Rd[a], P[Rdetk][j], e_det)]:
                    if np.isfinite(got) and np.isfinite(want):
                        acc.append(abs(np.log10(max(got, 1e-300)) - np.log10(max(want, 1e-300))))
                    elif np.isnan(got) != np.isnan(want):
                        acc.append(np.inf)
        m_all = max(e_all) if e_all else 0.0
        m_det = max(e_det) if e_det else 0.0
        good = (dkbad == 0) and m_all < 1e-8 and m_det < 1e-8
        ok &= good
        if verbose:
            print('  D4-GATE %-5s n=%4d  dk mismatch=%d  max|dlog10 R_all|=%.2e  max|dlog10 R_det|=%.2e  %s'
                  % (tag, len(P[realk]), dkbad, m_all, m_det, 'PASS' if good else 'FAIL'))
    return bool(ok)


# ================================================================= THE D4 PRE-REGISTRATION
PREREG = """
D4 -- THE REALISATION-LEVEL DISCORDANCE GATE.  PRE-REGISTERED, BINDING, FIXED BEFORE SCORING.

THE STATISTIC (chosen on the banked distributions, before any condition was scored):

    S_det = J^T I^-1 J   over the DETECTED set D,   I = sum_{b in D} g_b g_b^T / s2_b,
                                                    J = sum_{b in D} g_b (2 pi dk_b) / s2_b

  = chi2(u=0) - min_u chi2(u) on D: the GLRT for "the pulsars the data actually registered
  co-register at a source OTHER than the assumed counterpart".  FLAG => veto every
  certification in the realisation.  Undefined (no detections) => NO FLAG, and the
  realisation leaves the denominator (a gate with nothing to veto is vacuous, not passing).

WHY THIS AGGREGATE, AND WHY NOT THE OTHERS (justified on the banked distributions):
  - max/mean/min/frac of the per-candidate R_all: at the (-12.75,30,lit) onset cell TRUE-signal
    realisations have median min_R = 1.4e4 against a chi2(2) bar of 9.21 -- even the MOST
    concordant detection of a true realisation is discordant.  These aggregates inherit D3's
    poisoned leave-one-out reference in full and cannot separate.  REJECTED at design time.
  - S_ref over ALL 116 pulsars: TRUE median 2.2e4 vs wrong-counterpart 4e5 -- but PURE NOISE
    sits at 4.4e5, ABOVE the wrong-counterpart population.  It measures "this fringe field is
    poisoned", which is true of every realisation above onset (IGNITE-2 sec 1.2).  REJECTED.
  - S_det: TRUE-signal median 0.0 at BOTH onset cells -- the detected pulsars are the ones the
    data registered, so their consensus sits AT the assumed counterpart -- against 3e4 (lit) /
    2e2 (vlbi) for wrong counterparts.  The ONLY aggregate whose true-signal distribution
    concentrates at the null value.  SELECTED.  It is also the cheapest: one 2x2 solve per
    realisation, no leave-one-out loop.
  This INVERTS D3's variant ranking (R_all carried the per-pulsar test, R_det was the control)
  and IGNITE-2 sec 1.4 says why: "the R_det control kills 6/6 -- but for the degenerate reason
  that every detection is discordant with every other under a scrambled source, which is
  exactly the realisation-level (not per-pulsar) signal".  D4 is that sentence, made a gate.

TWO dk CONVENTIONS, BOTH SCORED (the grid is indexed -256..255 about the EM-PRIOR MEAN fringe,
measured: per-pulsar mean(n_true_grid) ~ 0, std ~ 150):
  - D4-ORACLE  dk = mapk - n_true_grid,  s2 = (2pi)^2[1/6 + (1-qmax)(sig_EM/dL)^2]
      D3's exact convention.  Referenced to the TRUE fringe, which a real analysis does not
      know.  Its numbers are an UPPER BOUND on any implementable form's power.
  - D4-OBS     dk = mapk,                s2 = (2pi)^2[1/6 + (sig_EM/dL)^2]
      Referenced to the EM-PRIOR MEAN fringe (k=0), which a real analysis DOES know.  The
      prior's own distance error is then present in every dk regardless of how confident the
      fringe posterior is, so the (1-qmax) factor is dropped -- forced, not tuned.

TWO THRESHOLDS, BOTH SCORED:
  - T-CHI2:  the analytic bar chi2(2, p=0.01) = 9.2103, inherited from D3, no new dial.
  - T-CAL:   the (1-alpha) quantile of S_det on a DISJOINT true-signal CALIBRATION set
             (IGNITE Stage-1, 50 reals/cell), alpha = 0.10, per this project's own convention
             ("a calibrated threshold must state its false-alarm rate and its sampling
             scatter").  The chi2(2) null is PROVABLY false above onset -- that is what killed
             D3 -- so an empirical bar at a stated FPR is the honest form.  Condition (ii) is
             then scored OUT OF SAMPLE on IGNITE-2's sloop realisations (disjoint seeds).

CONDITIONS (adopt only if (i) AND (ii) pass; (iii) reported either way):
  (i)   flag >= 95% of wrong-counterpart/scrambled realisations that DETECT
  (ii)  flag <= 10% of true-signal realisations that DETECT, at BOTH onset cells
  (iii) the small-|dk| keepers (J0437 dk=-4 class): does the realisation-level form catch what
        the per-pulsar form missed?  The D1 hole's closure test.  Reported honestly either way.
  (iv)  zero flags on pure-noise realisations with NO detections (vacuous, by construction --
        reported as such, not as evidence); the substantive number is the flag rate on
        pure-noise realisations that DO detect, which is reported beside it.
NO PARTIAL ADOPTION.  NO POST-HOC THRESHOLD TUNING.
"""


def S_det_of(lev, dk_true, mapk, qmax, dlnL, lnK, tier, T, geo, floor, form):
    """S_det under one of the two pre-registered dk conventions."""
    g, r2 = lev.get(T, geo, tier)
    if form == 'oracle':
        dk = np.asarray(dk_true, int)
        s2 = (2 * np.pi) ** 2 * (1.0 / 6.0 + (1.0 - np.asarray(qmax, float)) * r2)
    elif form == 'obs':
        dk = np.asarray(mapk, int)
        s2 = (2 * np.pi) ** 2 * (1.0 / 6.0 + r2)
    else:
        raise ValueError(form)
    det = np.asarray(dlnL) > np.maximum(np.asarray(lnK), floor)
    S, rk, gof = S_ref(dk, g, s2, det)
    cert = det & (np.asarray(qmax) > 0.9)
    return dict(S=S, rank=rk, gof=gof, n_det=int(det.sum()), n_cert=int(cert.sum()),
                n_cert_wrong=int((cert & (np.asarray(dk_true) != 0)).sum()), det=det, cert=cert)


if __name__ == '__main__':
    print('D4 value gate (reproduce IGNITE-2 banked co-registration statistic):')
    print('  ALL PASS' if d4_gate() else '  FAILED')
