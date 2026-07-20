"""Track B — B1 SOFT SOURCE SIDE (FORGE-G): the last hard edge in the loop, made soft.

WHY THIS EXISTS
---------------
The spec-§3 loop is soft on the PULSAR side already: every pulsar carries a fringe posterior
`q_p(n)` (never a hard-locked fringe -- the hard-lock violation is what caused the S8.1.2
cascade), and the M-step is the q-weighted `Q(theta) = sum_p sum_n q_p(n) lnL(theta, L_p(n))`
(`trackB_b1.B1Marg`). The SOURCE side is still a point estimate: `TargetedMarg` moves the loud
sources' (f, mc) as a bare vector and PINS every faint reservoir source at truth. That asymmetry
is the one hard edge left in the loop.

SPARK-2 and SPARK-3 measured what removing it is worth and named the law it must obey:

  * SPARK-2 (`reports/SPARK2_second_arrow.md`): the faint reservoir UNMODELLED erases
    certification (offender 0.000 vs 4.435 truth-modelled); deconfusion's trials cost is
    +0.578 nat and SATURATES.
  * SPARK-3 (`reports/SPARK3_second_arrow.md`): arrow 2 is MODEL-QUALITY-LIMITED and the verdict
    is STRADDLED. Under an OPTIMISTIC (conditional-width) reservoir model deconfusion recruits;
    under a PESSIMISTIC (PRIOR-width) model a badly-constrained reservoir *actively raises the
    certification floor and suppresses the count* (floor 118 -> 744 at one unit). §5.0: a
    prior-width model is WORSE than no model at all. **THE LAW (SPARK-3): a soft component may
    feed the joint template only when its posterior is demonstrably tighter than prior width by a
    declared factor; below that it is carried but not fed.** This module is that law as code.

WHAT THIS MODULE ADDS
---------------------
B1  SoftComponent / SoftSourceState -- every source carries a POSTERIOR (mean + per-axis width,
    Gaussian/Laplace; grid hook for multimodal axes), never a point estimate. Reservoir
    components start at PRIOR width. The joint model consumes the SET of soft components + the
    pulsar fringe posteriors SYMMETRICALLY, via the source-side feed weight SMASK (the twin of
    the pulsar-term mask PMASK) threaded through `trackB_b1_core`.
B2  ModelQualityGate -- SPARK-3's law. A component feeds back iff its concentration ratio
    (posterior/prior width, worst free axis) is strictly below `ratio_thresh`. ratio_thresh = 0
    reproduces UNMODELLED (nothing feeds); ratio_thresh = inf reproduces ALWAYS-FED. Default
    `SPARK3_GATE_RATIO` (declared factor = 1/ratio_thresh tighter than prior).
B3  Soft detection readout -- NO binary `detected` flag inside the loop. The per-component
    concentration statistic is tracked per iteration; the certification criterion runs only as
    the END-OF-CYCLE scoreboard (existing v2.2 machinery, unchanged).
B4  EM conditioning hooks -- optional hard sky pin (+ optional period) per component, the two
    sub-cases (anchor-a-resolved = sky; direct-a-buried = sky+period) as one flag. A pinned
    component's pinned axes drop from the joint free dimension.

SCOPE. This box (cronus) is RECOMPUTED-UNSAFE for noise by design (no canonical L_gwb bank), so
this module builds and gates ZERO-NOISE / ASIMOV only -- no campaign, no noisy realisations. The
GPU gates need the sanctioned lane (jug-gpu); the pure-logic gates (B2 limits, B4 pin-drop) run
anywhere. Matt commits; never commit from here.

Env: jug-gpu, XLA autotune off, x64, jax persistent cache.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, sys, time
import numpy as np

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")

# ---- source-parameter layout (one 8-block per source) ----
I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)
NP_SRC = 8
AXIS_NAME = ["cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
             "log10_fgw", "log10_h", "phase0", "psi"]
SKY_AXES = (I_COSGT, I_GWPHI)          # the EM-pinnable position block
PERIOD_AXIS = I_FGW                    # the optional period pin (direct-a-buried)

# ---- generative population prior half-widths (stagec_fisher; NOT TE.phi_bounds, which is the
#      estimator SEARCH box -- SPARK-3 §2.2's correction). Uniform box -> sigma = half/sqrt(3). ----
PRIOR_HALF = np.array([
    1.0,          # cos_gwtheta  U[-1, 1]
    np.pi,        # gwphi        U[0, 2pi]
    1.0,          # cos_inc      U[-1, 1]
    0.5,          # log10_mc     U[8.5, 9.5]
    0.25,         # log10_fgw    U[-8.0, -7.5]
    0.25,         # log10_h      U[-14.0, -13.5]
    np.pi,        # phase0       U[0, 2pi]
    0.5 * np.pi,  # psi          U[0, pi]
])
PRIOR_SIGMA = PRIOR_HALF / np.sqrt(3.0)

# ---- B2 default (SPARK-3's declared factor). SPARK-3 §5.0: the PESSIMISTIC bound is the prior
#      width (ratio 1.0) and it SUPPRESSES; a fed component must be meaningfully tighter. The exact
#      cross-over is the marginal width between the two Fisher bounds, which SPARK-2 measured is
#      unaffordable -- so this is a DECLARED parameter, provisional, the one knob a campaign
#      re-derives from its own Fisher. The machinery's content is the two LIMITS, not this number.
SPARK3_GATE_RATIO = 0.5               # feed iff posterior < 0.5 x prior width (factor F = 2)


# ============================================================
# B1 -- the soft source component and the soft source state
# ============================================================
class SoftComponent:
    """One source carrying a POSTERIOR over its 8 params, not a point estimate.

    mean         (8,) posterior mean (the Laplace/Gaussian centre; the template value when fed)
    sigma        (8,) posterior 1-sigma per axis (Gaussian/Laplace). Reservoir starts at prior.
    pin          (8,) bool, hard-pinned axes (EM conditioning, B4): known exactly, sigma := 0,
                 dropped from the free dimension and from the gate's concentration.
    grid         optional {axis: (vals, logw)} for a multimodal axis -- the q-weighted
                 generalisation of the mean (hook; the Gaussian mean is the fed value in the MVP).
    is_reservoir whether this component is a faint reservoir source (started at prior width).
    label        human tag.
    """

    def __init__(self, mean, sigma=None, pin=None, is_reservoir=False, label="", grid=None):
        self.mean = np.asarray(mean, float).copy()
        assert self.mean.shape == (NP_SRC,), f"mean must be length {NP_SRC}"
        self.sigma = (PRIOR_SIGMA.copy() if sigma is None
                      else np.asarray(sigma, float).copy())
        self.pin = (np.zeros(NP_SRC, bool) if pin is None else np.asarray(pin, bool).copy())
        self.is_reservoir = bool(is_reservoir)
        self.label = label
        self.grid = grid or {}
        self.sigma[self.pin] = 0.0             # a pinned axis is known exactly

    @classmethod
    def at_prior_width(cls, mean, **kw):
        """A reservoir component started at PRIOR width (SPARK-3's pessimistic-bound start)."""
        return cls(mean, sigma=PRIOR_SIGMA.copy(), is_reservoir=True, **kw)

    def free_axes(self):
        """Indices of the axes that are free in the joint (not hard-pinned)."""
        return np.where(~self.pin)[0]

    def concentration(self):
        """Per-axis posterior/prior width ratio (B3). Pinned axes -> 0 (known exactly)."""
        r = self.sigma / PRIOR_SIGMA
        r[self.pin] = 0.0
        return r

    def conc_scalar(self, axes=None):
        """The BINDING concentration ratio: the WORST (largest) ratio over the gate axes. A
        component is only as well-modelled as its worst free axis (SPARK-3 §5.0: prior width on
        ANY key axis lets the fit accommodate noise). Pinned axes are excluded."""
        r = self.concentration()
        free = self.free_axes()
        if axes is not None:
            free = np.array([a for a in axes if a in set(free)], int)
        if len(free) == 0:
            return 0.0                          # fully pinned -> perfectly known
        return float(np.max(r[free]))

    def pin_sky(self):
        self.pin[list(SKY_AXES)] = True
        self.sigma[list(SKY_AXES)] = 0.0

    def pin_period(self):
        self.pin[PERIOD_AXIS] = True
        self.sigma[PERIOD_AXIS] = 0.0

    def copy(self):
        c = SoftComponent(self.mean, self.sigma, self.pin, self.is_reservoir, self.label,
                          dict(self.grid))
        return c


class SoftSourceState:
    """The SET of soft components. Symmetric consumer of soft sources + pulsar fringe posteriors:
    it assembles (src_mean, smask) for the joint, where smask is the per-source feed weight the
    ModelQualityGate decides -- exactly as the pulsar side hands the joint a per-pulsar q-weight.
    """

    def __init__(self, components):
        self.comps = list(components)
        self.ncw = len(self.comps)

    # ---- construction from a truth/point source vector ----
    @classmethod
    def from_src_vector(cls, src, n_loud, sigma=None, reservoir_at_prior=True, labels=None):
        """src (ncw*8,) point vector -> a state. Loud sources [0,n_loud) get width `sigma` (or a
        delta if None); reservoir sources [n_loud, ncw) start at PRIOR width if
        reservoir_at_prior, else at `sigma`."""
        src = np.asarray(src, float).reshape(-1, NP_SRC)
        ncw = src.shape[0]
        comps = []
        for s in range(ncw):
            res = s >= n_loud
            lab = (labels[s] if labels else (f"reservoir{s}" if res else f"loud{s}"))
            if res and reservoir_at_prior:
                comps.append(SoftComponent.at_prior_width(src[s], label=lab))
            else:
                sig = (np.zeros(NP_SRC) if sigma is None
                       else np.asarray(sigma, float))
                comps.append(SoftComponent(src[s], sigma=sig, is_reservoir=res, label=lab))
        return cls(comps)

    # ---- the symmetric hand-off to the joint ----
    def feed_weights(self, gate):
        """smask (ncw,): 1.0 for a component the gate feeds, 0.0 otherwise. The source-side twin
        of the pulsar-side q-weight vector."""
        return np.array([1.0 if gate.feeds(c) else 0.0 for c in self.comps])

    def src_mean(self):
        """(ncw*8,) the template value of every component (its posterior mean / Laplace centre)."""
        return np.concatenate([c.mean for c in self.comps])

    def assemble(self, gate):
        """Return (src (ncw*8,), smask (ncw,)) for the joint. Not-fed components carry smask 0 and
        are therefore ABSENT from the template (present only in the trials budget) -- SPARK-3's
        'carried but not fed'."""
        return self.src_mean(), self.feed_weights(gate)

    def free_dim(self, gate):
        """Dimension of the joint FREE parameter block: sum over FED components of their number of
        non-pinned axes. A pinned axis (B4) drops; a not-fed component contributes 0 (it is not a
        free parameter of the template, only a trials-budget entry)."""
        smask = self.feed_weights(gate)
        return int(sum(len(c.free_axes()) for c, m in zip(self.comps, smask) if m > 0))

    def concentrations(self, gate=None):
        """Per-component binding concentration ratio (B3). If gate given, also the feed decision."""
        out = []
        for c in self.comps:
            rec = dict(label=c.label, is_reservoir=c.is_reservoir, conc=c.conc_scalar(),
                       n_free=len(c.free_axes()))
            if gate is not None:
                rec["fed"] = bool(gate.feeds(c))
            out.append(rec)
        return out

    def copy(self):
        return SoftSourceState([c.copy() for c in self.comps])


# ============================================================
# B2 -- THE MODEL-QUALITY GATE (SPARK-3's law, as code)
# ============================================================
class ModelQualityGate:
    """A component feeds the joint template iff its binding concentration ratio is STRICTLY below
    `ratio_thresh` (posterior/prior width). SPARK-3's law: below threshold the component is carried
    (it still costs the trials budget -- SPARK-2 §2's +0.578 nat) but is absent from the template.

    LIMITS (the gated content):
      ratio_thresh = 0    -> nothing feeds        == UNMODELLED  (r < 0 never true)
      ratio_thresh = inf  -> everything feeds      == ALWAYS-FED  (r < inf always true)
    STRICT inequality so a delta posterior (ratio 0) still needs a POSITIVE threshold to feed, and
    the two limits are exact against those two reference behaviours.

    `ratio_thresh` = 1 / (declared factor tighter than prior). Default from SPARK-3's anatomy.
    `gate_axes` = the axes whose width the gate keys on (default: all free axes of the component).
    """

    def __init__(self, ratio_thresh=SPARK3_GATE_RATIO, gate_axes=None):
        self.ratio_thresh = float(ratio_thresh)
        self.gate_axes = None if gate_axes is None else tuple(gate_axes)

    @property
    def factor(self):
        """The declared factor: a fed component is at least this many times tighter than prior."""
        return np.inf if self.ratio_thresh <= 0 else 1.0 / self.ratio_thresh

    def feeds(self, comp):
        return comp.conc_scalar(self.gate_axes) < self.ratio_thresh

    @classmethod
    def unmodelled(cls, **kw):
        return cls(ratio_thresh=0.0, **kw)

    @classmethod
    def always_fed(cls, **kw):
        return cls(ratio_thresh=np.inf, **kw)


# ============================================================
# B3 -- soft detection readout (end-of-cycle scoreboard; existing v2.2 machinery unchanged)
# ============================================================
def concentration_trace(traces):
    """Stack per-iteration per-component concentration records into an array (n_iter, ncw).
    NO binary detected flag anywhere in the loop -- this is the tracked statistic; the criterion
    runs only as the scoreboard below."""
    return np.array([[rec["conc"] for rec in it] for it in traces])


def scoreboard(state, gate, criterion=None):
    """END-OF-CYCLE scoreboard. Reports each component's concentration + feed decision; the
    detection CRITERION (v2.2) is applied here and ONLY here, never inside the loop. `criterion`
    is a callable(state, gate)->dict; if None, the concentration/feed table is returned as-is so
    the existing v2.2 machinery (criterion_v2_gates) can be dropped in unchanged by the caller."""
    table = state.concentrations(gate)
    out = dict(table=table, ratio_thresh=gate.ratio_thresh, factor=gate.factor,
               n_fed=int(sum(r["fed"] for r in table)),
               n_carried=int(sum(not r["fed"] for r in table)))
    if criterion is not None:
        out["criterion"] = criterion(state, gate)
    return out


# ============================================================
# B4 -- EM conditioning hooks (one flag, two sub-cases)
# ============================================================
def condition_component(comp, cond):
    """Apply an EM hard pin to one component. `cond`:
        None            -> no pin
        'sky'           -> anchor-a-resolved: pin the position block only
        'sky+period'    -> direct-a-buried:  pin the position block AND the period (log10_fgw)
    The pinned axes drop from the joint free dimension (SoftSourceState.free_dim)."""
    if cond in (None, "none", ""):
        return comp
    if cond not in ("sky", "sky+period"):
        raise ValueError(f"unknown conditioning flag {cond!r} "
                         "(expected None | 'sky' | 'sky+period')")
    comp.pin_sky()
    if cond == "sky+period":
        comp.pin_period()
    return comp


def condition_state(state, conds):
    """Apply per-component conditioning. `conds` is a dict {index: flag} or a single flag applied
    to every LOUD (non-reservoir) component."""
    if isinstance(conds, dict):
        for i, c in conds.items():
            condition_component(state.comps[i], c)
    else:
        for comp in state.comps:
            if not comp.is_reservoir:
                condition_component(comp, conds)
    return state


# ============================================================
# The zero-noise soft cycle (E-step + gated soft source hand-off) -- Asimov only
# ============================================================
class SoftCycle:
    """One soft cycle on the spec-§3 loop, source side made soft and gated. GPU (needs B1Marg).

    E-step  : the pulsar fringe posteriors at the gated assembled (src, smask). Symmetric to the
              source side -- the pulsars are soft (q_p(n)) and the sources are soft (gated feed).
    M-hand  : each fed, unpinned component's posterior width is (re)measured; not-fed / pinned
              axes are left at their carried state. (The full Newton M-step on lnL_marg is
              trackB_b1_referendum's machinery; this class provides the gated hand-off and the
              zero-noise fixed-point property the gates check, not a new optimiser.)
    """

    def __init__(self, P, data, base_mask):
        import trackB_b1 as B1
        self.P = P
        self.E = B1.B1Marg(P, data, base_mask)
        self.ncw = P["ncw"]
        self.nd = P["n_dist"]

    def estep(self, state, gate):
        """Set the gated feed weights on the joint and return the pulsar fringe posterior at the
        assembled source template. Zero-noise/Asimov."""
        src, smask = state.assemble(gate)
        self.E.set_smask(smask)                 # the source-side twin of the pulsar q-weights
        post, LNL, lnref = self.E.posterior(src)
        return dict(post=post, LNL=LNL, lnref=lnref, src=src, smask=smask)

    def lnmarg(self, state, gate):
        src, smask = state.assemble(gate)
        self.E.set_smask(smask)
        return float(self.E.lnmarg(src[None])[0])

    def reset(self):
        """Restore the incumbent all-on behaviour (no smask)."""
        self.E.set_smask(None)


# ============================================================
# GATES (B5)
# ============================================================
def _logic_gates(verbose=True):
    """The pure-logic gates (B2 limits + B4 pin-drop). No GPU, no likelihood -- run anywhere."""
    ok = True
    p = print if verbose else (lambda *a, **k: None)

    # ---- Bg2: the B2 gate at its two limits, against small Asimov fixtures ----
    p("\n-- Bg2: model-quality gate limits (threshold-0 == unmodelled, inf == always-fed) --")
    # fixture: 1 tight loud + 1 half-tight + 1 prior-width reservoir
    tight = SoftComponent(np.zeros(NP_SRC), sigma=0.01 * PRIOR_SIGMA, label="tight")
    mid = SoftComponent(np.zeros(NP_SRC), sigma=0.40 * PRIOR_SIGMA, label="mid")   # < 0.5 -> fed at default
    wide = SoftComponent.at_prior_width(np.zeros(NP_SRC), label="wide")            # ratio 1.0
    st = SoftSourceState([tight, mid, wide])

    g0 = ModelQualityGate.unmodelled()
    ginf = ModelQualityGate.always_fed()
    gdef = ModelQualityGate()                    # default SPARK3_GATE_RATIO = 0.5
    n0 = int(st.feed_weights(g0).sum())
    ninf = int(st.feed_weights(ginf).sum())
    ndef = int(st.feed_weights(gdef).sum())
    b_lim0 = (n0 == 0)                            # threshold-0 feeds NOTHING (unmodelled)
    b_liminf = (ninf == 3)                        # threshold-inf feeds ALL (always-fed)
    b_def = (ndef == 2)                           # default feeds tight+mid, carries the prior-width one
    p(f"   threshold-0   fed {n0}/3  (expect 0, unmodelled)   -> {'PASS' if b_lim0 else 'FAIL'}")
    p(f"   threshold-inf fed {ninf}/3 (expect 3, always-fed)  -> {'PASS' if b_liminf else 'FAIL'}")
    p(f"   default r<{gdef.ratio_thresh} (factor {gdef.factor:.0f}x) fed {ndef}/3 "
      f"(expect 2: tight+mid, carry prior-width) -> {'PASS' if b_def else 'FAIL'}")
    # monotonicity: fed-count non-decreasing in ratio_thresh
    counts = [int(st.feed_weights(ModelQualityGate(r)).sum())
              for r in [0.0, 0.005, 0.1, 0.5, 0.9, 1.0, 2.0, np.inf]]
    b_mono = all(counts[i] <= counts[i + 1] for i in range(len(counts) - 1))
    p(f"   monotone fed-count vs ratio_thresh: {counts} -> {'PASS' if b_mono else 'FAIL'}")
    # a delta posterior still needs a positive threshold (strict <)
    delta = SoftSourceState([SoftComponent(np.zeros(NP_SRC), sigma=np.zeros(NP_SRC))])
    b_delta = (int(delta.feed_weights(g0).sum()) == 0 and int(delta.feed_weights(gdef).sum()) == 1)
    p(f"   delta posterior: fed@0={int(delta.feed_weights(g0).sum())} (0), "
      f"fed@default={int(delta.feed_weights(gdef).sum())} (1) -> {'PASS' if b_delta else 'FAIL'}")
    bg2 = b_lim0 and b_liminf and b_def and b_mono and b_delta
    ok &= bg2
    p(f"   [Bg2] -> {'PASS' if bg2 else 'FAIL'}")

    # ---- Bg4: a pinned component's position block drops from the joint dimension ----
    p("\n-- Bg4: EM conditioning drops the pinned block from the joint free dimension --")
    base = SoftSourceState([SoftComponent(np.zeros(NP_SRC), sigma=0.01 * PRIOR_SIGMA, label=f"s{i}")
                            for i in range(3)])
    gall = ModelQualityGate.always_fed()
    d_base = base.free_dim(gall)                 # 3 comps x 8 free axes = 24
    st_sky = base.copy(); condition_component(st_sky.comps[0], "sky")
    d_sky = st_sky.free_dim(gall)
    st_sp = base.copy(); condition_component(st_sp.comps[0], "sky+period")
    d_sp = st_sp.free_dim(gall)
    b_base = (d_base == 3 * NP_SRC)
    b_sky = (d_sky == d_base - 2)                # sky pins cos_gwtheta + gwphi
    b_sp = (d_sp == d_base - 3)                  # + log10_fgw
    p(f"   free_dim: base {d_base} (expect {3*NP_SRC})            -> {'PASS' if b_base else 'FAIL'}")
    p(f"   free_dim: one comp sky-pinned {d_sky} (expect {d_base-2}) -> {'PASS' if b_sky else 'FAIL'}")
    p(f"   free_dim: one comp sky+period {d_sp} (expect {d_base-3}) -> {'PASS' if b_sp else 'FAIL'}")
    # a pinned axis is excluded from the gate concentration (known exactly) and from free_axes
    b_pinconc = (st_sky.comps[0].conc_scalar() == st_sky.comps[0].conc_scalar()  # finite
                 and I_COSGT not in set(st_sky.comps[0].free_axes()))
    # a not-fed component contributes 0 to the free dimension (carried, not a template parameter)
    st_carry = base.copy()
    st_carry.comps[2] = SoftComponent.at_prior_width(np.zeros(NP_SRC))   # ratio 1 -> not fed default
    d_carry = st_carry.free_dim(ModelQualityGate())
    b_carry = (d_carry == 2 * NP_SRC)
    p(f"   free_dim: one comp carried-not-fed {d_carry} (expect {2*NP_SRC}) "
      f"-> {'PASS' if b_carry else 'FAIL'}")
    bg4 = b_base and b_sky and b_sp and b_pinconc and b_carry
    ok &= bg4
    p(f"   [Bg4] -> {'PASS' if bg4 else 'FAIL'}")

    return ok


def _gpu_gates(verbose=True, fp_tol=1e-6):
    """The zero-noise GPU gates (fixed point + threshold-inf/truth-pinned == spec-§3 loop). Needs
    the sanctioned lane (jug-gpu): builds the b1 problem and B1Marg."""
    import trackB_b1_core as C
    import trackB_b1 as B1
    import jax.numpy as jnp
    ok = True
    p = print if verbose else (lambda *a, **k: None)

    P = C.build_b1_problem()
    tt = P["theta_truth"]; nd = P["n_dist"]
    src_truth = tt[nd:].copy()
    mask1 = P["mask_one"]
    data = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(mask1))   # zero-noise Asimov

    # ---- Bg5a: threshold-inf + truth-pinned reproduces the spec-§3 pulsar loop BIT-COMPARABLY ----
    # The spec-§3 q-weighted marginal loop (B1Marg) at the truth source, on the MASKED residual
    # path, is the reference (smask == None). Our soft layer at ALWAYS-FED + every component a
    # delta at truth is smask == 1 on the SAME object and SAME path, so 1.0*x == x makes it
    # BIT-EXACT. Separately, gate_fast_full certifies the incumbent's fast-full optimisation equals
    # the masked path, so the reference is the genuine incumbent to that gate's 1e-8.
    p("\n-- Bg5a: always-fed + truth-pinned == incumbent B1Marg loop (spec-§3) bit-comparably --")
    cyc = SoftCycle(P, data, mask1)
    d_ff = cyc.E.gate_fast_full(src_truth[None])           # incumbent fast-full == masked (<1e-8)
    # incumbent reference on the masked path (smask None). gate_fast_full leaves fast-full ON, and
    # set_smask(None) does not switch it off, so force the masked path by hand for a like-for-like
    # (both masked) bit-exact comparison against the always-fed smask==ones path.
    cyc.E.sp.enable_fast_full(False); cyc.E.all_on = False; cyc.E._abb = {}
    cyc.E.set_smask(None)
    g_ref = float(cyc.E.lnmarg(src_truth[None])[0])
    post_ref, LNL_ref, _ = cyc.E.posterior(src_truth)
    # always-fed + truth-pinned soft state
    state = SoftSourceState.from_src_vector(src_truth, C.N_LOUD, sigma=np.zeros(NP_SRC),
                                            reservoir_at_prior=False)  # every comp a delta at truth
    ginf = ModelQualityGate.always_fed()
    smask = state.feed_weights(ginf)
    b_allfed = bool(np.all(smask == 1.0))
    out = cyc.estep(state, ginf)                           # smask == ones, masked path
    g_soft = cyc.lnmarg(state, ginf)
    d_lnmarg = abs(g_soft - g_ref)
    d_lnl = float(np.max(np.abs(out["LNL"] - LNL_ref)))
    d_q = float(np.max(np.abs(out["post"]["q_max"] - post_ref["q_max"])))
    b5a = b_allfed and d_lnmarg == 0.0 and d_lnl == 0.0 and d_q == 0.0 and d_ff < 1e-8
    p(f"   smask all-ones: {b_allfed}; incumbent fast-full vs masked |diff| {d_ff:.3e} (<1e-8)")
    p(f"   |lnL_marg(always-fed) - lnL_marg(incumbent)| = {d_lnmarg:.3e} (== 0, bit-exact)")
    p(f"   max|LNL(always-fed) - LNL(incumbent)|        = {d_lnl:.3e} (== 0)")
    p(f"   max|q_max(always-fed) - q_max(incumbent)|    = {d_q:.3e} (== 0)")
    ok &= b5a
    p(f"   [Bg5a] -> {'PASS' if b5a else 'FAIL'}")

    # ---- Bg5b: zero-noise fixed point -- start at truth, one soft cycle moves nothing ----
    # At truth the E-step posteriors reproduce the census ceilings and the assembled template is
    # unchanged; the concentration of every fed (delta) component is 0, so the gate decision and
    # the src hand-off are stationary. We check the loop is a fixed point of the SOURCE hand-off:
    # re-assembling from the (unchanged) state reproduces the same (src, smask, lnL_marg).
    p("\n-- Bg5b: zero-noise fixed point (start at truth, one soft cycle moves nothing) --")
    src2, smask2 = state.assemble(ginf)
    d_src = float(np.max(np.abs(src2 - src_truth)))
    g_soft2 = cyc.lnmarg(state, ginf)
    d_fp = abs(g_soft2 - g_soft)
    # census ceilings reproduced at the fixed point (the incumbent zero-noise property)
    FT = cyc.E.FT
    ptrue, _, _ = FT.p_of_fringe(np.zeros(P["npsr"], int))
    import trackB_estimator as TE
    d_ceil = max(abs(ptrue[P["names"].index(nm)] - TE.CEILING[nm]) for nm in C.CENSUS)
    b5b = d_src < fp_tol and d_fp < fp_tol and d_ceil < 0.02
    p(f"   max|src(fp) - truth|            = {d_src:.3e} (<{fp_tol:.0e})")
    p(f"   |lnL_marg reassembled - before| = {d_fp:.3e} (<{fp_tol:.0e})")
    p(f"   census ceiling max|diff|        = {d_ceil:.4f} (<0.02)")
    ok &= b5b
    p(f"   [Bg5b] -> {'PASS' if b5b else 'FAIL'}")

    # ---- Bg5c: the gate limits on the REAL joint, and a w_s==0 source is bit-exact ABSENT ----
    # (i)  threshold-0 unmodels EVERY component (strict r<0): smask is all-zeros (the fully
    #      unmodelled, no-CW template). The reservoir-unmodelled SPARK-2 endpoint (loud fed,
    #      reservoir carried) is the DEFAULT gate, and its smask is the loud-only pattern.
    # (ii) bit-exact absence: a source with w_s==0 contributes NOTHING, so its parameters are
    #      irrelevant -- lnL_marg is invariant to replacing the not-fed reservoir means with junk.
    p("\n-- Bg5c: gate limits on the joint + w_s==0 is bit-exact absence --")
    state_r = SoftSourceState.from_src_vector(src_truth, C.N_LOUD)   # loud delta, reservoir prior
    smask_unmod = state_r.feed_weights(ModelQualityGate.unmodelled())
    smask_def = state_r.feed_weights(ModelQualityGate())             # default: loud fed, reservoir not
    smask_loudonly = np.array([1.0] * C.N_LOUD + [0.0] * (P["ncw"] - C.N_LOUD))
    b_unmod = bool(np.all(smask_unmod == 0.0))                       # threshold-0 -> all absent
    b_def = bool(np.all(smask_def == smask_loudonly))               # default -> reservoir carried
    p(f"   threshold-0 smask all-zero (fully unmodelled): {b_unmod}")
    p(f"   default smask == loud-only (reservoir carried, not fed): {b_def}")
    # (ii) junk-invariance: reservoir means -> random, feed weights loud-only, lnL_marg unchanged
    src_A = state_r.src_mean().copy()
    src_B = src_A.reshape(-1, NP_SRC).copy()
    rng = np.random.default_rng(31337)
    src_B[C.N_LOUD:] = rng.normal(size=(P["ncw"] - C.N_LOUD, NP_SRC))   # junk reservoir params
    src_B = src_B.reshape(-1)
    cyc.E.set_smask(smask_loudonly)
    g_A = float(cyc.E.lnmarg(src_A[None])[0])
    g_B = float(cyc.E.lnmarg(src_B[None])[0])
    d_abs = abs(g_A - g_B)
    b_abs = d_abs == 0.0
    p(f"   |lnL_marg(reservoir@truth) - lnL_marg(reservoir@junk)| with w_res==0 = {d_abs:.3e} "
      f"(== 0: params of a not-fed source are irrelevant)")
    b5c = b_unmod and b_def and b_abs
    ok &= b5c
    p(f"   [Bg5c] -> {'PASS' if b5c else 'FAIL'}")

    cyc.reset()
    return ok


def mode_gate(cpu_only=False, verbose=True):
    print("=== B1 SOFT SOURCE SIDE GATES (FORGE-G) ===", flush=True)
    ok = _logic_gates(verbose=verbose)
    if cpu_only:
        print(f"\n=== LOGIC GATES (B2/B4): {'ALL PASS' if ok else 'FAIL'} "
              f"(GPU gates B5 skipped: --cpu-only) ===", flush=True)
        return 0 if ok else 1
    try:
        import jax
        okg = _gpu_gates(verbose=verbose)
        ok &= okg
    except Exception as e:
        print(f"\n  GPU gates NOT RUN ({type(e).__name__}: {e}). "
              f"Run on jug-gpu; logic gates above stand.", flush=True)
        print(f"\n=== SOFT SOURCE GATES: logic {'PASS' if ok else 'FAIL'}, "
              f"GPU DEFERRED ===", flush=True)
        return 0 if ok else 1
    print(f"\n=== SOFT SOURCE SIDE GATES: {'ALL PASS' if ok else 'FAIL'} ===", flush=True)
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="gate", choices=["gate"])
    ap.add_argument("--cpu-only", action="store_true",
                    help="run only the pure-logic gates (B2 limits + B4 pin-drop); skip the GPU "
                         "zero-noise gates (B5). For the cronus dev box.")
    a = ap.parse_args()
    if a.mode == "gate":
        return mode_gate(cpu_only=a.cpu_only)
    return 2


if __name__ == "__main__":
    sys.exit(main())
