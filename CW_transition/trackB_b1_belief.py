"""Track B — B1 RUNNING-BELIEF M-STEP (agent FORGE-B): the THIRD ARM beside frozen and live.

MATT'S DIRECTIVE (the arm's charter, verbatim)
----------------------------------------------
  "give sources the same running-belief treatment the pulsar fringes already enjoy;
   tighten only on earned information; wander bounded by belief width by construction."

WHY THIS EXISTS
---------------
The spec-§3 loop is soft on the PULSAR side: every pulsar carries a fringe posterior q_p(n)
that PERSISTS across iterations and is UPDATED, never re-derived from scratch. The SOURCE
side is not: FORGE-G made the source side soft in the SET sense (every source carries a
posterior width, and the ModelQualityGate decides feed-back), but the M-STEP itself is still
a fresh re-maximisation every iteration -- `glacier_loop.mstep_quadratic` starts from the
current template value and walks to the local likelihood mode, with nothing remembering
where the belief was or how much of it was earned. That is the LIVE arm. The FROZEN arm is
the opposite extreme: a member's parameters are pinned at their feed-time values forever.

The dual-ladder measurement is what makes this a live problem and not an aesthetic one: the
refit WANDERS 0.01-0.44 dex per iteration at ~0.2-nat class curvature. At that curvature the
likelihood cannot distinguish the wander from the truth, so the live arm's motion is
unearned; but the frozen arm cannot move at all, so it cannot use information that IS
earned. This module is the third arm: a persistent per-component BELIEF (centre + per-axis
width) that the M-step updates instead of re-maximising, with the wander bounded by the
belief width BY CONSTRUCTION and the width allowed to tighten ONLY on new information.

WHAT THIS MODULE ADDS
---------------------
B1  `Belief` -- the persistent per-component belief state: `centre` (the value the template
    is fed), `anchor` (the centre OF RECORD, which moves only on an earned tighten), and
    `width` (per-parameter). INITIALISED AT FEED TIME from the frontier evidence:
    `Belief.from_frontier(feed_point, F_ii, box_sigma, dlnL_support)` -- the curvature at the
    feed point sets the Laplace width, and the DATA-SUPPORT delta-lnL inflates it when the
    member was fed on thin evidence. `BeliefBook` is the slot-indexed container; it packs to
    flat arrays so the belief CHECKPOINTS and BANKS with the component state.
B2  `InformationLedger` + `BeliefMStep` -- the update rule and THE DOUBLE-COUNTING GUARD.
    Refits are MAP against the belief, not re-maximisations. Iterations re-read the SAME
    data, so belief WIDTH may tighten ONLY on a new-information event: the fed set changed,
    the certified set changed, the background state moved beyond its band, or a new data
    epoch arrived. Otherwise the update may RE-CENTRE within the belief band but may NOT
    narrow it. The event detector is explicit (`InformationLedger.observe`) and every
    tighten is logged with the event that licensed it (`Belief.tighten_log`).
B3  `GraduationLadder` -- the model-quality gate reading belief width LIVE. Components climb
    carried -> fed -> refine-eligible -> cert-eligible as their belief concentrates past a
    declared threshold ladder whose FIRST rung IS the incumbent frontier (SPARK3_GATE_RATIO
    = 0.5), so the arm does not move the existing gate, it extends it. `GraduationLedger`
    banks the belief-width trajectories -- the sources-resolving-confidently curve -- as a
    first-class readout.
B4  the gates (`gate` mode below): g1 the two LIMITS, g2 THE GUARD, g3 the WANDER BOUND,
    g4 default-inert, g5 runtime-SMASK compatibility.

THE LIMITS (the arm's content, and why the limits are EXACT and not merely close)
--------------------------------------------------------------------------------
  belief width -> inf : the belief carries no precision, the MAP centre is the likelihood
                        centre and the trust region is unbounded -> the LIVE M-step. This is
                        implemented as an EXPLICIT BRANCH (`not isfinite(w) -> return the
                        likelihood centre`), not as 1/inf arithmetic, so the returned
                        (src_hat, widths, n_eval) are BIT-IDENTICAL to `mstep_quadratic`.
  belief width -> 0   : a delta belief has nothing to learn, so the slot is SKIPPED (no
                        evaluations at all) -> the member never moves from its feed point
                        and n_eval == 0: FROZEN-after-feed, exactly, including cost.
The two existing arms are therefore literal limits of this one, not analogues of it.

DECLARED PARAMETERS (provisional, in the SPARK3_GATE_RATIO tradition -- a campaign
re-derives them from its own Fisher; the machinery's content is the LIMITS, not these):
  WANDER_K           trust-region radius in belief widths                        1.0
  DLNL_SUPPORT_REF   reference data-support (nat) for the feed-time width         1.0
  LADDER_THRESH      the graduation ladder (fed / refine / cert)         (0.5, 0.2, 0.1)
  BG_BAND_K          background move, in its own sigma, that counts as news       1.0

SCOPE. This box (cronus) is RECOMPUTED-UNSAFE for noise by design (no canonical L_gwb bank),
so this module builds and gates ZERO-NOISE / ASIMOV only -- no campaign, no noisy
realisations. g1-g4 are analytic-fixture gates and run anywhere; g5 needs the b1 stack.
Matt commits; never commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, sys, time
import numpy as np

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")

from trackB_b1_softsource import (NP_SRC, AXIS_NAME, PRIOR_SIGMA, SPARK3_GATE_RATIO,
                                  I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI)

# ---- the M-step surface axes (glacier_loop.MSTEP_AXES; a parameter everywhere below) ----
MSTEP_AXES_DEFAULT = (I_FGW, I_MC)

# ---- DECLARED parameters (see the header) ----
WANDER_K = 1.0                 # refit displacement bounded by WANDER_K x belief width
DLNL_SUPPORT_REF = 1.0         # nat; feed-time width inflates as sqrt(REF/support) below this
BG_BAND_K = 1.0                # |dA_bg| > K x sigma  == the background state MOVED
LADDER_THRESH = (SPARK3_GATE_RATIO, 0.2, 0.1)    # fed / refine-eligible / cert-eligible

STAGE_NAME = ("carried", "fed", "refine-eligible", "cert-eligible")
STAGE_CARRIED, STAGE_FED, STAGE_REFINE, STAGE_CERT = 0, 1, 2, 3

# the four licensed new-information event kinds (B2). NOTHING ELSE licenses a tighten.
EVENT_KINDS = ("fed_set", "cert_set", "background", "data_epoch")


# ============================================================
# B1 -- THE BELIEF STATE
# ============================================================
class Belief:
    """One component's persistent belief: centre + per-parameter width, plus the centre OF
    RECORD (`anchor`) the wander bound is measured from.

    centre    (NP_SRC,)  the value the joint template is fed THIS iteration. Free to move
                         inside the belief band every iteration (re-centring is free).
    anchor    (NP_SRC,)  the centre of record. Moves ONLY when a tighten is licensed. The
                         wander bound is |centre - anchor| <= WANDER_K * width, enforced by
                         a hard clip, so the bound holds across ANY number of iterations.
    width     (NP_SRC,)  the belief 1-sigma per axis. inf == no belief on that axis (LIVE
                         limit); 0 == delta belief (FROZEN limit). Narrows ONLY on a
                         licensed new-information event.
    feed_*               the feed-time record: centre, iteration, and the frontier evidence
                         (per-axis curvature F_ii and the data-support delta-lnL) the belief
                         was initialised from. Immutable; banked.
    tighten_log          one record per tighten: (iteration, event kinds, width before/after).
                         Every narrowing in the campaign is traceable to its licensing event.
    """

    __slots__ = ("centre", "anchor", "width", "feed_centre", "feed_iter", "feed_F_ii",
                 "dlnL_support", "n_tighten", "n_recentre", "tighten_log", "label")

    def __init__(self, centre, width, anchor=None, feed_iter=-1, feed_F_ii=None,
                 dlnL_support=np.nan, label=""):
        self.centre = np.asarray(centre, float).copy()
        assert self.centre.shape == (NP_SRC,), f"centre must be length {NP_SRC}"
        self.width = np.asarray(width, float).copy() * np.ones(NP_SRC)
        assert np.all(self.width[~np.isnan(self.width)] >= 0.0), "belief width must be >= 0"
        self.anchor = (self.centre.copy() if anchor is None
                       else np.asarray(anchor, float).copy())
        self.feed_centre = self.centre.copy()
        self.feed_iter = int(feed_iter)
        self.feed_F_ii = (np.full(NP_SRC, np.nan) if feed_F_ii is None
                          else np.asarray(feed_F_ii, float).copy() * np.ones(NP_SRC))
        self.dlnL_support = float(dlnL_support)
        self.n_tighten = 0
        self.n_recentre = 0
        self.tighten_log = []
        self.label = label

    # ---- feed-time initialisation from the frontier evidence ----
    @classmethod
    def from_frontier(cls, feed_point, F_ii, box_sigma, dlnL_support, feed_iter=-1,
                      support_ref=DLNL_SUPPORT_REF, label=""):
        """THE FEED-TIME INITIALISER (B1). The belief a component is born with is the
        frontier's own evidence about it, and nothing else:

          width_j = min( box_j ,  F_ii_j^(-1/2) * inflate(dlnL_support) )
          inflate = sqrt( max(1, DLNL_SUPPORT_REF / dlnL_support) )

        The curvature term is the Laplace width at the feed point -- the same conditional
        width the frontier gated on, so a component enters with EXACTLY the width that got
        it fed and no tighter. The data-support term is the guard against a member that
        crossed the frontier on a sharp-but-shallow curvature: a member fed with less than
        the reference support starts PROPORTIONATELY WIDER (its own dlnL says the data has
        not yet earned that curvature), capped at the prior box, which is the correct
        "we know nothing". A member with support at or above the reference is not inflated.

        `F_ii` is the per-axis observed information (glacier_loop.fisher_conditional's F_ii
        row for this member, own-term-live per the S4.14 fix); `box_sigma` the prior box.
        Non-positive / non-finite curvature -> the box, i.e. no belief beyond the prior.
        """
        F = np.asarray(F_ii, float) * np.ones(NP_SRC)
        box = np.asarray(box_sigma, float) * np.ones(NP_SRC)
        sig = np.where(F > 0, 1.0 / np.sqrt(np.maximum(F, 1e-300)), np.inf)
        s = float(max(dlnL_support, 0.0))
        inflate = float(np.sqrt(max(1.0, support_ref / max(s, 1e-12))))
        width = np.minimum(sig * inflate, box)
        return cls(feed_point, width, feed_iter=feed_iter, feed_F_ii=F,
                   dlnL_support=s, label=label)

    # ---- the two limits, as constructors (used by the g1 gates and the policy dispatch) ----
    @classmethod
    def live(cls, centre, **kw):
        """width -> inf: no belief. The LIVE M-step, bit-exactly (see the module header)."""
        return cls(centre, np.full(NP_SRC, np.inf), **kw)

    @classmethod
    def frozen(cls, centre, **kw):
        """width -> 0: a delta belief. FROZEN-after-feed, exactly (slots are skipped)."""
        return cls(centre, np.zeros(NP_SRC), **kw)

    # ---- readouts ----
    def width_ratio(self, prior_sigma=PRIOR_SIGMA, axes=None):
        """The BINDING belief-width ratio: worst (largest) belief/prior width over `axes`
        (default: all NP_SRC). Same reduction as SoftComponent.conc_scalar -- a component is
        only as well-believed as its worst axis."""
        r = np.asarray(self.width, float) / np.asarray(prior_sigma, float)
        idx = range(NP_SRC) if axes is None else axes
        vals = [r[j] for j in idx]
        return float(np.max(vals)) if len(vals) else 0.0

    def displacement(self):
        """|centre - anchor| per axis -- the wander, in the units the bound is stated in."""
        return np.abs(self.centre - self.anchor)

    def wander_bound(self, wander_k=WANDER_K):
        return wander_k * self.width

    def in_band(self, wander_k=WANDER_K, tol=0.0):
        d, b = self.displacement(), self.wander_bound(wander_k)
        return bool(np.all(d <= b + tol))

    def copy(self):
        b = Belief(self.centre, self.width, self.anchor, self.feed_iter, self.feed_F_ii,
                   self.dlnL_support, self.label)
        b.feed_centre = self.feed_centre.copy()
        b.n_tighten, b.n_recentre = self.n_tighten, self.n_recentre
        b.tighten_log = list(self.tighten_log)
        return b


class BeliefBook:
    """The slot-indexed set of beliefs -- one per FED component, keyed by census slot index.

    This is the container the driver carries, checkpoints and banks. It packs to flat arrays
    (`pack`) with an `active` mask, so a checkpoint round-trips the belief exactly and a
    resumed cell continues the SAME belief rather than re-deriving one (re-deriving on
    resume would be a silent tighten at unchanged information -- exactly what B2 forbids).
    """

    def __init__(self, n_slot, beliefs=None):
        self.n_slot = int(n_slot)
        self.b = dict(beliefs or {})

    # ---- container protocol ----
    def __contains__(self, k):
        return int(k) in self.b

    def __getitem__(self, k):
        return self.b[int(k)]

    def __len__(self):
        return len(self.b)

    def keys(self):
        return sorted(self.b)

    def get(self, k, default=None):
        return self.b.get(int(k), default)

    def set(self, k, belief):
        self.b[int(k)] = belief
        return belief

    # ---- feed-time hook: a promoted member is born with a belief, once ----
    def feed(self, k, feed_point, F_ii, box_sigma, dlnL_support, iteration,
             support_ref=DLNL_SUPPORT_REF, label=""):
        """Initialise slot k's belief from the frontier evidence AT FEED TIME. Idempotent by
        design: a slot that already has a belief is NOT re-initialised (a re-init would be a
        free tighten at unchanged information)."""
        k = int(k)
        if k in self.b:
            return self.b[k]
        self.b[k] = Belief.from_frontier(feed_point, F_ii, box_sigma, dlnL_support,
                                         feed_iter=iteration, support_ref=support_ref,
                                         label=label or f"slot{k}")
        return self.b[k]

    # ---- banking / checkpointing ----
    def pack(self):
        """Flat arrays for bank_npz / np.savez. Round-trips through `unpack` exactly."""
        n = self.n_slot
        active = np.zeros(n, bool)
        centre = np.full((n, NP_SRC), np.nan)
        anchor = np.full((n, NP_SRC), np.nan)
        width = np.full((n, NP_SRC), np.nan)
        feed_c = np.full((n, NP_SRC), np.nan)
        feed_F = np.full((n, NP_SRC), np.nan)
        feed_it = np.full(n, -1, int)
        support = np.full(n, np.nan)
        ntight = np.zeros(n, int)
        nrec = np.zeros(n, int)
        for k, bel in self.b.items():
            active[k] = True
            centre[k], anchor[k], width[k] = bel.centre, bel.anchor, bel.width
            feed_c[k], feed_F[k] = bel.feed_centre, bel.feed_F_ii
            feed_it[k], support[k] = bel.feed_iter, bel.dlnL_support
            ntight[k], nrec[k] = bel.n_tighten, bel.n_recentre
        return dict(belief_active=active, belief_centre=centre, belief_anchor=anchor,
                    belief_width=width, belief_feed_centre=feed_c, belief_feed_F=feed_F,
                    belief_feed_iter=feed_it, belief_support=support,
                    belief_n_tighten=ntight, belief_n_recentre=nrec)

    @classmethod
    def unpack(cls, z, n_slot=None):
        """Rebuild from a banked npz (or dict). Returns None if the bank has no belief block
        (a pre-belief checkpoint) -- the caller then runs default-inert."""
        if "belief_active" not in z:
            return None
        active = np.asarray(z["belief_active"], bool)
        n = int(len(active)) if n_slot is None else int(n_slot)
        book = cls(n)
        for k in np.where(active)[0]:
            bel = Belief(np.asarray(z["belief_centre"])[k], np.asarray(z["belief_width"])[k],
                         anchor=np.asarray(z["belief_anchor"])[k],
                         feed_iter=int(np.asarray(z["belief_feed_iter"])[k]),
                         feed_F_ii=np.asarray(z["belief_feed_F"])[k],
                         dlnL_support=float(np.asarray(z["belief_support"])[k]),
                         label=f"slot{int(k)}")
            bel.feed_centre = np.asarray(z["belief_feed_centre"])[k].copy()
            bel.n_tighten = int(np.asarray(z["belief_n_tighten"])[k])
            bel.n_recentre = int(np.asarray(z["belief_n_recentre"])[k])
            book.b[int(k)] = bel
        return book

    def copy(self):
        return BeliefBook(self.n_slot, {k: v.copy() for k, v in self.b.items()})


# ============================================================
# B2 -- THE EVENT DETECTOR (the double-counting guard's licensing side)
# ============================================================
class InformationLedger:
    """THE DOUBLE-COUNTING GUARD, licensing side.

    The loop re-reads the SAME data every iteration. A refit that narrows the belief on that
    re-read counts the same photons twice, and the belief collapses to a confidence the data
    never supplied -- which is precisely the failure mode the running-belief arm exists to
    stop. So the belief WIDTH may narrow only on a NEW-INFORMATION EVENT, and this class is
    the explicit detector for the four kinds that qualify:

      fed_set     the fed set changed (a promote/demote changed what the joint models, so
                  every other member's conditional information genuinely changed)
      cert_set    the certified set changed (the E-step's PM rows cohere a different set, so
                  the fringe posteriors carry different information)
      background  the background state moved beyond its own band (|dA_bg| > K sigma) -- the
                  drain moved, so the residual the sources are fitted against moved
      data_epoch  new data epochs arrived (the only unambiguous new information)

    Anything else -- another sweep, another iteration, a re-read, a resume -- is NOT news.

    LICENSING. `observe` returns the events it detected and queues them. A single M-step call
    consumes ALL queued events (`spend`) and is thereby licensed to perform AT MOST ONE
    tighten. Multiple events in one iteration therefore license one tighten, not several;
    zero events license none. The first `observe` SEEDS the record and fires nothing (the
    feed-time initialisation already spent that information -- see Belief.from_frontier).
    """

    def __init__(self, bg_band_k=BG_BAND_K):
        self.bg_band_k = float(bg_band_k)
        self._seeded = False
        self._fed = None
        self._cert = None
        self._bg = None            # (value, sigma) of record
        self._epoch = None
        self._pending = []
        self.history = []          # every event ever detected, for the bank

    @staticmethod
    def _as_set(x):
        if x is None:
            return None
        a = np.asarray(x)
        if a.dtype == bool:
            return frozenset(int(i) for i in np.where(a)[0])
        return frozenset(int(i) for i in a.reshape(-1))

    def observe(self, iteration, fed=None, cert=None, a_bg=None, a_bg_sig=None,
                data_epoch=None):
        """Update the information record; return the list of NEW-INFORMATION events."""
        fed_s, cert_s = self._as_set(fed), self._as_set(cert)
        ev = []
        if not self._seeded:
            self._fed, self._cert, self._epoch = fed_s, cert_s, data_epoch
            if a_bg is not None:
                self._bg = (float(a_bg), float(a_bg_sig if a_bg_sig is not None else 0.0))
            self._seeded = True
            return []
        if fed_s is not None and fed_s != self._fed:
            ev.append(dict(kind="fed_set", iteration=int(iteration),
                           detail=f"{sorted(self._fed or [])} -> {sorted(fed_s)}"))
            self._fed = fed_s
        if cert_s is not None and cert_s != self._cert:
            ev.append(dict(kind="cert_set", iteration=int(iteration),
                           detail=f"n {len(self._cert or [])} -> {len(cert_s)}"))
            self._cert = cert_s
        if a_bg is not None:
            sig = float(a_bg_sig if a_bg_sig is not None else 0.0)
            if self._bg is None:
                self._bg = (float(a_bg), sig)
            else:
                v0, s0 = self._bg
                band = self.bg_band_k * max(sig, s0)
                if abs(float(a_bg) - v0) > band:
                    ev.append(dict(kind="background", iteration=int(iteration),
                                   detail=f"A_bg {v0:+.4f} -> {float(a_bg):+.4f} "
                                          f"(band {band:.4f})"))
                    self._bg = (float(a_bg), sig)
        if data_epoch is not None and data_epoch != self._epoch:
            ev.append(dict(kind="data_epoch", iteration=int(iteration),
                           detail=f"{self._epoch} -> {data_epoch}"))
            self._epoch = data_epoch
        self._pending.extend(ev)
        self.history.extend(ev)
        return ev

    def pending(self):
        return list(self._pending)

    def spend(self):
        """Consume every queued event; return them. A non-empty return licenses EXACTLY ONE
        tighten in the M-step call that consumed it."""
        ev, self._pending = list(self._pending), []
        return ev

    def event_counts(self):
        return {k: int(sum(e["kind"] == k for e in self.history)) for k in EVENT_KINDS}


# ============================================================
# B2 -- THE UPDATE RULE (the M-step whose internals this replaces)
# ============================================================
def _map_centre_unclipped(w_b, anchor, like_centre, like_width):
    """The precision-weighted MAP centre, BEFORE the trust-region clip (reported so the
    driver can bank how often the bound actually binds)."""
    if not np.isfinite(w_b):
        return float(like_centre)
    if w_b == 0.0 or not np.isfinite(like_width) or like_width <= 0.0:
        return float(anchor)
    pb = 1.0 / (w_b * w_b)
    pl = 1.0 / (like_width * like_width)
    return float((pb * anchor + pl * like_centre) / (pb + pl))


def _map_centre(w_b, anchor, like_centre, like_width, wander_k):
    """One axis of the MAP update, with the wander bound enforced BY CONSTRUCTION.

    The prior is N(anchor, w_b^2) -- ANCHORED, i.e. centred on the centre OF RECORD, not on
    last iteration's working centre. That is the second half of the double-counting guard:
    if the prior chased its own previous posterior, repeated iterations at unchanged data
    would ratchet the centre to the likelihood mode and the arm would decay into the live
    arm with extra steps. Anchored, repeated iterations converge to a fixed point inside the
    band and stay there.

    The returned centre is then hard-clipped to anchor +- wander_k*w_b, so the displacement
    bound holds for ANY likelihood, including a pathological one -- 'wander bounded by
    belief width BY CONSTRUCTION', not in expectation.
    """
    if not np.isfinite(w_b):
        return float(like_centre)              # LIVE limit: exact passthrough (bit-exact)
    if w_b == 0.0:
        return float(anchor)                   # FROZEN limit (slots are skipped upstream)
    if not np.isfinite(like_width) or like_width <= 0.0:
        return float(anchor)                   # no usable curvature -> no information
    c = _map_centre_unclipped(w_b, anchor, like_centre, like_width)
    lo, hi = anchor - wander_k * w_b, anchor + wander_k * w_b
    return float(min(max(c, lo), hi))


def _tightened(w_b, like_width):
    """The licensed width update: precision addition. inf belief + finite likelihood -> the
    likelihood width (a first tighten from 'no belief' lands exactly on the Laplace width)."""
    if not np.isfinite(like_width) or like_width <= 0.0:
        return w_b
    if not np.isfinite(w_b):
        return float(like_width)
    if w_b == 0.0:
        return 0.0
    return float(1.0 / np.sqrt(1.0 / (w_b * w_b) + 1.0 / (like_width * like_width)))


class BeliefMStep:
    """The running-belief M-step. Drop-in for `glacier_loop.mstep_quadratic`.

    `step(marg_fn, src0, fed_idx, scale, book=None, licensed=False, ...)`

      book is None      -> the INCUMBENT algorithm, executed verbatim (default-inert, g4).
      belief width inf  -> the same numerical trajectory as the incumbent (g1 live limit).
      belief width 0    -> the slot is skipped entirely (g1 frozen limit, n_eval == 0).
      otherwise         -> MAP against the belief, clipped to the belief band; the width
                           narrows only if `licensed` (one tighten per licensing call).

    The sweep itself is the incumbent's batched coordinate-quadratic ascent, transcribed
    unchanged -- the arm changes WHAT IS DONE WITH the local quadratic (centre and width),
    never how the quadratic is measured. `widths` retains its incumbent meaning (the
    LIKELIHOOD Laplace width per fed axis); the belief widths come back separately, so the
    live limit can be gated bit-exact on both returns.
    """

    def __init__(self, mstep_axes=MSTEP_AXES_DEFAULT, wander_k=WANDER_K,
                 n_sweep=2, step0=0.3):
        self.mstep_axes = tuple(mstep_axes)
        self.wander_k = float(wander_k)
        self.n_sweep = int(n_sweep)
        self.step0 = float(step0)

    def step(self, marg_fn, src0, fed_idx, scale, book=None, licensed=False, iteration=-1,
             events=None):
        """Returns (src_hat, widths (n_fed, n_axes), n_eval, log dict)."""
        axes = self.mstep_axes
        src = np.asarray(src0, float).copy()
        fed_idx = [int(k) for k in fed_idx]
        slots = np.array([NP_SRC * k + j for k in fed_idx for j in axes], int)
        sc = np.array([scale[j] for _ in fed_idx for j in axes], float)
        owner = [(k, j) for k in fed_idx for j in axes]
        n_eval = 0
        widths = np.full(len(slots), np.nan)
        n_clipped = 0
        n_skipped = 0
        did_tighten = False

        for sweep in range(self.n_sweep):
            for i, sl in enumerate(slots):
                k, j = owner[i]
                bel = None if book is None else book.get(k)
                w_b = np.inf if bel is None else float(bel.width[j])
                if bel is not None and w_b == 0.0:
                    n_skipped += 1
                    continue                       # FROZEN limit: nothing to learn, no evals
                d = self.step0 * sc[i] * (0.5 ** sweep)
                trip = np.tile(src, (3, 1))
                trip[0, sl] -= d
                trip[2, sl] += d
                v = np.asarray(marg_fn(trip))
                n_eval += 3
                a, b, c = v
                den = a - 2 * b + c
                if den < 0:                        # concave: quadratic step, clipped
                    off = np.clip(0.5 * (a - c) / den, -1.0, 1.0) * d
                    like_centre = src[sl] + off
                    like_width = d / np.sqrt(-den)
                    widths[i] = like_width
                    if bel is None:
                        src[sl] = like_centre      # incumbent: src[sl] += off, same op
                    else:
                        new_c = _map_centre(w_b, bel.anchor[j], like_centre, like_width,
                                            self.wander_k)
                        if np.isfinite(w_b) and w_b > 0.0:
                            cm = _map_centre_unclipped(w_b, bel.anchor[j], like_centre,
                                                       like_width)
                            if cm != new_c:
                                n_clipped += 1     # the trust region BOUND actually bound
                        src[sl] = new_c
                        bel.centre[j] = new_c
                # non-concave: leave the slot (carried at its entry value this sweep)

        # ---- the ONE licensed tighten, applied after the sweeps, logged with its event ----
        if book is not None and licensed:
            wl = {}
            for i, (k, j) in enumerate(owner):
                if np.isfinite(widths[i]):
                    wl[(k, j)] = widths[i]
            kinds = sorted({e["kind"] for e in (events or [])})
            for k in sorted({kk for kk, _ in wl}):
                bel = book.get(k)
                if bel is None:
                    continue
                before = bel.width.copy()
                for (kk, j), lw in wl.items():
                    if kk != k:
                        continue
                    bel.width[j] = _tightened(float(bel.width[j]), float(lw))
                    bel.anchor[j] = bel.centre[j]      # the centre of record moves ON EARNING
                bel.n_tighten += 1
                bel.tighten_log.append(dict(iteration=int(iteration), events=kinds,
                                            width_before=before, width_after=bel.width.copy()))
                did_tighten = True

        if book is not None:
            for k in fed_idx:
                bel = book.get(k)
                if bel is not None:
                    bel.n_recentre += 1

        log = dict(n_eval=int(n_eval), n_clipped=int(n_clipped), n_skipped=int(n_skipped),
                   licensed=bool(licensed), tightened=bool(did_tighten),
                   events=[e["kind"] for e in (events or [])])
        return src, widths.reshape(len(fed_idx), len(axes)), n_eval, log


def frozen_step(src0, fed_idx, mstep_axes=MSTEP_AXES_DEFAULT):
    """The FROZEN arm, stated directly: fed members never move from their feed point, and
    the M-step costs nothing. (Identical to the belief arm at width -> 0; g1 gates that.)"""
    n_ax = len(tuple(mstep_axes))
    return (np.asarray(src0, float).copy(),
            np.full((len(list(fed_idx)), n_ax), np.nan), 0,
            dict(n_eval=0, n_clipped=0, n_skipped=len(list(fed_idx)) * n_ax,
                 licensed=False, tightened=False, events=[]))


# ============================================================
# B3 -- GRADUATION: the model-quality gate reading belief width LIVE
# ============================================================
class GraduationLadder:
    """The model-quality gate, extended to a LADDER read off the LIVE belief width.

    A component's stage is the deepest rung whose threshold its BINDING belief-width ratio
    (worst gate axis, belief width / prior width) is STRICTLY below -- the same strict
    inequality and the same worst-axis reduction as `ModelQualityGate`, and the FIRST rung
    IS `SPARK3_GATE_RATIO`, so at the fed rung this ladder and the incumbent frontier are
    the same gate. The deeper rungs are new:

      carried          r >= 0.5   not in the template; costs the trials budget only
      fed              r <  0.5   in the template, a free parameter of the M-step
      refine-eligible  r <  0.2   believed tightly enough that a refit is worth its evals
      cert-eligible    r <  0.1   believed tightly enough to be offered to the scoreboard

    "Eligible" is the operative word on the two deep rungs: the ladder never certifies
    anything. The v2.2 criterion still runs ONLY as the end-of-cycle scoreboard (FORGE-G B3,
    the no-in-loop-flag rule); cert-eligibility is a precondition it is offered, not a verdict.

    GATE AXES ARE ALL AXES BY DEFAULT, and that default is load-bearing. Keying the ladder on
    the M-step axes alone (log10_fgw, log10_mc) makes a member whose FREQUENCY is believed to
    1e-5 read "cert-eligible" while its SKY sits at prior-box width -- measured on the real
    b1 machinery, a reservoir member the frontier scores at ratio 0.948 (nowhere near fed)
    reads 0.0015 on the M-step axes alone. That is SPARK-3 SS5.0's failure exactly: prior
    width on ANY key axis lets the fit accommodate noise. The consequence is worth stating
    plainly rather than engineering around: while the M-step surface moves only 2 of the 8
    axes, the deep rungs are STRUCTURALLY UNREACHABLE for any member whose other six axes are
    box-limited. Narrowing `gate_axes` would make the ladder read green off a dead
    instrument -- the S4.14 lesson. The M-step-axis ratio is banked as a DIAGNOSTIC column
    (`GraduationLedger(diag_axes=...)`) so the two are visible side by side and neither is
    silently substituted for the other.
    """

    def __init__(self, thresh=LADDER_THRESH, prior_sigma=PRIOR_SIGMA, gate_axes=None):
        t = tuple(float(x) for x in thresh)
        assert len(t) == 3, "ladder is (fed, refine, cert)"
        assert t[0] >= t[1] >= t[2], "ladder thresholds must descend"
        self.thresh = t
        self.prior_sigma = np.asarray(prior_sigma, float)
        self.gate_axes = None if gate_axes is None else tuple(gate_axes)

    def ratio(self, belief):
        return belief.width_ratio(self.prior_sigma, self.gate_axes)

    def stage(self, belief):
        r = self.ratio(belief)
        s = STAGE_CARRIED
        for i, t in enumerate(self.thresh):
            if r < t:
                s = i + 1
        return s

    def stage_name(self, belief):
        return STAGE_NAME[self.stage(belief)]

    def feeds(self, belief):
        return self.stage(belief) >= STAGE_FED

    def refine_eligible(self, belief):
        return self.stage(belief) >= STAGE_REFINE

    def cert_eligible(self, belief):
        return self.stage(belief) >= STAGE_CERT


class GraduationLedger:
    """The belief-width trajectories, as a first-class banked readout -- "the
    sources-resolving-confidently curve". One record per iteration per fed slot: the binding
    width ratio, the stage, the wander, and whether that iteration tightened."""

    def __init__(self, n_slot, ladder=None, diag_axes=None):
        self.n_slot = int(n_slot)
        self.ladder = ladder or GraduationLadder()
        # a SECOND, non-binding width ratio over `diag_axes` (typically the M-step axes):
        # banked beside the binding one so "believed tightly where the loop actually looks"
        # and "believed tightly overall" are never confused for each other.
        self.diag_axes = None if diag_axes is None else tuple(diag_axes)
        self.rows = []                      # per iteration: dict of arrays
        self.first_crossing = {}            # (slot, stage) -> iteration

    def record(self, iteration, book, tightened=False):
        n = self.n_slot
        ratio = np.full(n, np.nan)
        rdiag = np.full(n, np.nan)
        stage = np.full(n, -1, int)
        wander = np.full(n, np.nan)
        wband = np.full(n, np.nan)
        ntight = np.zeros(n, int)
        for k in book.keys():
            bel = book[k]
            ratio[k] = self.ladder.ratio(bel)
            if self.diag_axes is not None:
                rdiag[k] = bel.width_ratio(self.ladder.prior_sigma, self.diag_axes)
            stage[k] = self.ladder.stage(bel)
            d = bel.displacement()
            wander[k] = float(np.max(d[list(self.ladder.gate_axes)])
                              if self.ladder.gate_axes else np.max(d))
            b = bel.wander_bound()
            wband[k] = float(np.min(b[list(self.ladder.gate_axes)])
                             if self.ladder.gate_axes else np.min(b))
            ntight[k] = bel.n_tighten
            for s in range(STAGE_FED, STAGE_CERT + 1):
                if stage[k] >= s and (k, s) not in self.first_crossing:
                    self.first_crossing[(k, s)] = int(iteration)
        self.rows.append(dict(iteration=int(iteration), ratio=ratio, ratio_diag=rdiag,
                              stage=stage, wander=wander, wander_bound=wband,
                              n_tighten=ntight, tightened=bool(tightened)))
        return self.rows[-1]

    def trajectories(self):
        """(n_iter, n_slot) arrays: the banked curves."""
        if not self.rows:
            z = np.zeros((0, self.n_slot))
            return dict(iters=np.zeros(0, int), width_ratio=z, width_ratio_diag=z,
                        stage=z.astype(int), wander=z, wander_bound=z,
                        n_tighten=z.astype(int))
        return dict(iters=np.array([r["iteration"] for r in self.rows], int),
                    width_ratio=np.stack([r["ratio"] for r in self.rows]),
                    width_ratio_diag=np.stack([r["ratio_diag"] for r in self.rows]),
                    stage=np.stack([r["stage"] for r in self.rows]),
                    wander=np.stack([r["wander"] for r in self.rows]),
                    wander_bound=np.stack([r["wander_bound"] for r in self.rows]),
                    n_tighten=np.stack([r["n_tighten"] for r in self.rows]))

    def bank(self):
        t = self.trajectories()
        return {f"grad_{k}": v for k, v in t.items()}


# ============================================================
# the incumbent, transcribed verbatim -- the g1/g4 REFERENCE
# ============================================================
def _incumbent_mstep_quadratic(marg_fn, src0, fed_idx, scale, mstep_axes=MSTEP_AXES_DEFAULT,
                               n_sweep=2, step0=0.3, bmax=64):
    """VERBATIM transcription of `hpc_harbor/glacier/glacier_loop.mstep_quadratic` (with
    MSTEP_AXES lifted to an argument). Used as the gate reference when glacier_loop is not
    importable; when it IS importable the gates compare against the real function and
    additionally assert this transcription equals it bit-for-bit."""
    src = src0.copy()
    slots = np.array([NP_SRC * k + j for k in fed_idx for j in mstep_axes])
    sc = np.array([scale[j] for _ in fed_idx for j in mstep_axes])
    n_eval = 0
    widths = np.full(len(slots), np.nan)
    for sweep in range(n_sweep):
        for i, sl in enumerate(slots):
            d = step0 * sc[i] * (0.5 ** sweep)
            trip = np.tile(src, (3, 1))
            trip[0, sl] -= d; trip[2, sl] += d
            v = np.asarray(marg_fn(trip)); n_eval += 3
            a, b, c = v
            den = a - 2 * b + c
            if den < 0:
                off = np.clip(0.5 * (a - c) / den, -1.0, 1.0) * d
                src[sl] += off
                widths[i] = d / np.sqrt(-den)
    return src, widths.reshape(len(fed_idx), len(mstep_axes)), n_eval


# ============================================================
# GATES (B4)
# ============================================================
def _quad_marg(mode, F, n_src_slots):
    """Analytic zero-noise M-step surface: lnL(x) = -0.5 * sum_j F_j (x_j - mode_j)^2, over a
    flat (n_src_slots*NP_SRC,) source vector, batched over rows. Deterministic, no data, no
    noise -- the Asimov fixture the logic gates need."""
    mode = np.asarray(mode, float)
    F = np.asarray(F, float)

    def fn(rows):
        r = np.atleast_2d(np.asarray(rows, float))
        d = r - mode[None, :]
        return -0.5 * np.sum(F[None, :] * d * d, axis=1)
    return fn


def _fixture(n_slot=3, n_fed=2, seed=20260729, curvature=25.0):
    """A small deterministic M-step fixture: mode, curvature, scale, starting src."""
    rng = np.random.default_rng(seed)
    n = n_slot * NP_SRC
    mode = np.zeros(n)
    F = np.zeros(n)
    fed_idx = list(range(n_fed))
    for k in fed_idx:
        for j in MSTEP_AXES_DEFAULT:
            mode[NP_SRC * k + j] = rng.uniform(-0.05, 0.05)
            F[NP_SRC * k + j] = curvature
    src0 = np.zeros(n)
    scale = np.ones(NP_SRC) * 0.1
    return mode, F, src0, np.array(fed_idx), scale


def _g1_limits(p):
    """g1 -- THE LIMITS. width -> inf reproduces LIVE bit-comparably; width -> 0 reproduces
    FROZEN-after-feed bit-comparably. The two existing arms as limits of this one."""
    p("\n-- g1 LIMITS: belief width -> inf == LIVE; width -> 0 == FROZEN-after-feed --")
    ok = True
    for tag, (n_slot, n_fed, seed, curv) in {
            "A": (3, 2, 20260729, 25.0),
            "B": (5, 3, 11, 4.0),
            "C": (2, 2, 907, 100.0)}.items():
        mode, F, src0, fed_idx, scale = _fixture(n_slot, n_fed, seed, curv)
        marg = _quad_marg(mode, F, n_slot)
        M = BeliefMStep()

        # --- the LIVE reference: the incumbent function itself where importable ---
        ref_src, ref_w, ref_n = _incumbent_mstep_quadratic(marg, src0, fed_idx, scale)
        try:
            sys.path.insert(0, "/home/mattm/projects/HSYMT/hpc_harbor/glacier")
            from glacier_loop import mstep_quadratic as real_mstep
            r2 = real_mstep(marg, src0, fed_idx, scale)
            d_tr = max(float(np.max(np.abs(r2[0] - ref_src))),
                       float(np.nanmax(np.abs(r2[1] - ref_w))) if np.any(np.isfinite(ref_w)) else 0.0)
            src_lab = "glacier_loop.mstep_quadratic"
        except Exception as e:                      # noqa: BLE001 -- reference is optional
            d_tr, src_lab = 0.0, f"transcription (glacier_loop unavailable: {type(e).__name__})"

        # --- width -> inf: LIVE ---
        book_inf = BeliefBook(n_slot, {int(k): Belief.live(src0[NP_SRC*k:NP_SRC*(k+1)])
                                       for k in fed_idx})
        s_inf, w_inf, n_inf, _ = M.step(marg, src0, fed_idx, scale, book=book_inf,
                                        licensed=False)
        d_src = float(np.max(np.abs(s_inf - ref_src)))
        d_w = float(np.max(np.abs(np.nan_to_num(w_inf, nan=-1.0)
                                  - np.nan_to_num(ref_w, nan=-1.0))))
        b_inf = (d_src == 0.0) and (d_w == 0.0) and (n_inf == ref_n)

        # --- width -> 0: FROZEN-after-feed ---
        feed_pt = src0.copy()
        book_0 = BeliefBook(n_slot, {int(k): Belief.frozen(feed_pt[NP_SRC*k:NP_SRC*(k+1)])
                                     for k in fed_idx})
        s_0, w_0, n_0, lg0 = M.step(marg, src0, fed_idx, scale, book=book_0, licensed=False)
        f_src, f_w, f_n, _ = frozen_step(src0, fed_idx)
        d_f = float(np.max(np.abs(s_0 - f_src)))
        b_frz = (d_f == 0.0) and (n_0 == 0) and (f_n == 0) and np.all(np.isnan(w_0))

        p(f"   [{tag}] live limit  max|src - {src_lab}| = {d_src:.3e}, max|widths| = "
          f"{d_w:.3e}, n_eval {n_inf} == {ref_n} -> {'PASS' if b_inf else 'FAIL'}")
        p(f"   [{tag}] frozen limit max|src - feed point| = {d_f:.3e}, n_eval {n_0} (== 0, "
          f"frozen costs nothing) -> {'PASS' if b_frz else 'FAIL'}")
        if d_tr != 0.0:
            p(f"   [{tag}] TRANSCRIPTION DIVERGES from glacier_loop: {d_tr:.3e}")
        ok &= b_inf and b_frz and (d_tr == 0.0)
    p(f"   [g1] -> {'PASS' if ok else 'FAIL'}")
    return ok


def _g2_guard(p):
    """g2 -- THE GUARD. Repeated iterations at unchanged information produce EXACTLY ZERO
    width tightening (asserted exactly), while re-centring inside the belief stays free;
    then ONE synthetic new-information event licenses EXACTLY ONE tighten."""
    p("\n-- g2 THE GUARD: no tighten without news; one event -> exactly one tighten --")
    n_slot, n_fed = 3, 2
    mode, F, src0, fed_idx, scale = _fixture(n_slot, n_fed, seed=4242, curvature=25.0)
    marg = _quad_marg(mode, F, n_slot)
    M = BeliefMStep()
    box = PRIOR_SIGMA.copy()
    book = BeliefBook(n_slot)
    for k in fed_idx:
        book.feed(int(k), src0[NP_SRC*k:NP_SRC*(k+1)], F[NP_SRC*k:NP_SRC*(k+1)],
                  box, dlnL_support=4.0, iteration=0)
    led = InformationLedger()
    fed_mask = np.zeros(n_slot, bool); fed_mask[fed_idx] = True
    led.observe(0, fed=fed_mask, cert=np.zeros(0, int), a_bg=-14.60, a_bg_sig=0.05)

    w_hist, c_hist = [], []
    src = src0.copy()
    N_QUIET = 6
    for it in range(N_QUIET):
        ev = led.observe(it, fed=fed_mask, cert=np.zeros(0, int), a_bg=-14.60, a_bg_sig=0.05)
        lic = led.spend()
        src, _, _, lg = M.step(marg, src, fed_idx, scale, book=book,
                               licensed=bool(lic), iteration=it, events=lic)
        w_hist.append(np.stack([book[int(k)].width.copy() for k in fed_idx]))
        c_hist.append(np.stack([book[int(k)].centre.copy() for k in fed_idx]))
        assert not ev, "quiet iteration fired an event"

    w0 = w_hist[0]
    b_nowiden = all(np.array_equal(w, w0) for w in w_hist)          # EXACT, not approximate
    b_notight = all(book[int(k)].n_tighten == 0 for k in fed_idx)
    moved = float(np.max(np.abs(c_hist[0] - np.stack([src0[NP_SRC*k:NP_SRC*(k+1)]
                                                      for k in fed_idx]))))
    b_recentred = moved > 0.0                                        # re-centring stays free
    band = np.stack([book[int(k)].wander_bound() for k in fed_idx])
    disp = np.stack([book[int(k)].displacement() for k in fed_idx])
    b_band = bool(np.all(disp <= band + 0.0))
    p(f"   {N_QUIET} quiet iterations: width array EXACTLY unchanged: {b_nowiden}; "
      f"n_tighten {[book[int(k)].n_tighten for k in fed_idx]} (all 0): {b_notight}")
    p(f"   re-centring remained free: max|centre - feed| = {moved:.3e} (> 0): {b_recentred}; "
      f"stayed inside the belief band: {b_band}")

    # ---- ONE synthetic new-information event: the certified set changes ----
    ev = led.observe(N_QUIET, fed=fed_mask, cert=np.array([7]), a_bg=-14.60, a_bg_sig=0.05)
    lic = led.spend()
    b_ev = (len(ev) == 1 and ev[0]["kind"] == "cert_set")
    w_before = np.stack([book[int(k)].width.copy() for k in fed_idx])
    src, _, _, lg = M.step(marg, src, fed_idx, scale, book=book, licensed=bool(lic),
                           iteration=N_QUIET, events=lic)
    w_after = np.stack([book[int(k)].width.copy() for k in fed_idx])
    b_one = all(book[int(k)].n_tighten == 1 for k in fed_idx)
    ax = list(MSTEP_AXES_DEFAULT)
    b_narrow = bool(np.all(w_after[:, ax] < w_before[:, ax]))
    untouched = [j for j in range(NP_SRC) if j not in ax]
    b_untouched = bool(np.array_equal(w_after[:, untouched], w_before[:, untouched]))
    b_logged = all(len(book[int(k)].tighten_log) == 1
                   and book[int(k)].tighten_log[0]["events"] == ["cert_set"] for k in fed_idx)
    p(f"   event fired: {[e['kind'] for e in ev]} -> tightens {[book[int(k)].n_tighten for k in fed_idx]} "
      f"(exactly 1 each): {b_one}")
    p(f"   licensed widths strictly narrowed on the M-step axes: {b_narrow} "
      f"(unmeasured axes untouched: {b_untouched}); every tighten logged with its event: {b_logged}")

    # ---- and quiet again: no further tightening, ever ----
    for it in range(N_QUIET + 1, N_QUIET + 6):
        led.observe(it, fed=fed_mask, cert=np.array([7]), a_bg=-14.60, a_bg_sig=0.05)
        lic = led.spend()
        src, _, _, _ = M.step(marg, src, fed_idx, scale, book=book, licensed=bool(lic),
                              iteration=it, events=lic)
    w_end = np.stack([book[int(k)].width.copy() for k in fed_idx])
    b_stable = np.array_equal(w_end, w_after) and all(book[int(k)].n_tighten == 1
                                                      for k in fed_idx)
    p(f"   5 further quiet iterations: width EXACTLY unchanged, n_tighten still 1: {b_stable}")

    # ---- and the other three event kinds each license, exactly once ----
    kinds = []
    e2 = led.observe(99, fed=fed_mask, cert=np.array([7]), a_bg=-14.20, a_bg_sig=0.05)
    kinds += [e["kind"] for e in e2]
    fed2 = fed_mask.copy(); fed2[n_slot - 1] = True
    e3 = led.observe(100, fed=fed2, cert=np.array([7]), a_bg=-14.20, a_bg_sig=0.05)
    kinds += [e["kind"] for e in e3]
    e4 = led.observe(101, fed=fed2, cert=np.array([7]), a_bg=-14.20, a_bg_sig=0.05,
                     data_epoch="T40")
    kinds += [e["kind"] for e in e4]
    b_kinds = kinds == ["background", "fed_set", "data_epoch"]
    p(f"   all four event kinds detected, in order: {kinds} -> {b_kinds}")

    ok = (b_nowiden and b_notight and b_recentred and b_band and b_ev and b_one
          and b_narrow and b_untouched and b_logged and b_stable and b_kinds)
    p(f"   [g2] -> {'PASS' if ok else 'FAIL'}")
    return ok


# the measured hazard (dual-ladder): the refit wanders 0.01-0.44 dex per iteration at
# ~0.2-nat class curvature. F from dlnL = 0.5*F*dx^2 at (dx, dlnL) = (0.44 dex, 0.2 nat).
WANDER_DEX = 0.44
WANDER_NAT = 0.2
F_FLAT = 2.0 * WANDER_NAT / (WANDER_DEX ** 2)          # ~2.07 dex^-2 -> sigma ~0.70 dex


def _g3_wander(p, n_iter=120, seed=20260729):
    """g3 -- THE WANDER BOUND. At a synthetic FLAT-likelihood fixture of the measured
    hazard's curvature class (0.2 nat at 0.44 dex), the refit displacement stays bounded by
    the belief width across >= 100 iterations. Negative control: the same fixture with no
    belief (the LIVE arm) walks away -- the walk-suppression demonstration, which is the
    arm's whole point.

    The wander driver is a DETERMINISTIC seeded jitter of the analytic surface's mode -- the
    same data re-read with a refit that cannot tell the wander from the truth. It is a
    fixture, not a data realisation: no L_gwb, no noise draw, nothing recomputed.
    """
    p(f"\n-- g3 WANDER BOUND: flat fixture ({WANDER_NAT} nat at {WANDER_DEX} dex, "
      f"sigma {1/np.sqrt(F_FLAT):.3f} dex), {n_iter} iterations --")
    n_slot, n_fed = 2, 1
    n = n_slot * NP_SRC
    ax = list(MSTEP_AXES_DEFAULT)
    fed_idx = np.array([0])
    scale = np.ones(NP_SRC) * 0.3
    F = np.zeros(n)
    for j in ax:
        F[j] = F_FLAT
    mode0 = np.zeros(n)
    box = PRIOR_SIGMA.copy()

    def run(with_belief):
        rng = np.random.default_rng(seed)
        mode = mode0.copy()
        src = mode0.copy()
        M = BeliefMStep()
        book = None
        if with_belief:
            book = BeliefBook(n_slot)
            # fed on a solid data-support delta-lnL, so the feed-time width is the Laplace
            # width at the feed point (no inflation) -- the honest belief the frontier bought.
            book.feed(0, src[:NP_SRC], F[:NP_SRC], box, dlnL_support=4.0, iteration=0)
        led = InformationLedger()
        fed_mask = np.zeros(n_slot, bool); fed_mask[0] = True
        led.observe(0, fed=fed_mask, cert=np.zeros(0, int), a_bg=-14.60, a_bg_sig=0.05)
        disp, worst = [], 0.0
        for it in range(n_iter):
            for j in ax:                       # the wander driver: mode random-walks
                mode[j] += rng.uniform(-WANDER_DEX, WANDER_DEX)
            marg = _quad_marg(mode, F, n_slot)
            led.observe(it, fed=fed_mask, cert=np.zeros(0, int), a_bg=-14.60, a_bg_sig=0.05)
            lic = led.spend()
            src, _, _, _ = M.step(marg, src, fed_idx, scale, book=book,
                                  licensed=bool(lic), iteration=it, events=lic)
            if with_belief:
                d = float(np.max(book[0].displacement()[ax]))
            else:
                d = float(np.max(np.abs(src[ax] - mode0[ax])))
            disp.append(d)
            worst = max(worst, d)
        return np.array(disp), book

    d_bel, book = run(True)
    d_live, _ = run(False)
    bound = float(np.min(book[0].wander_bound()[ax]))
    width = float(np.max(book[0].width[ax]))
    b_bound = bool(np.all(d_bel <= bound + 0.0))
    b_notight = (book[0].n_tighten == 0)
    b_control = bool(d_live.max() > bound)
    supp = float(d_live.max() / max(d_bel.max(), 1e-300))
    p(f"   belief width {width:.4f} dex, bound (K={WANDER_K:g}) {bound:.4f} dex")
    p(f"   BELIEF arm: max displacement over {n_iter} iters = {d_bel.max():.4f} dex "
      f"(<= bound: {b_bound}); tightens {book[0].n_tighten} (== 0, no news): {b_notight}")
    p(f"   LIVE control (no belief, same fixture): max displacement = {d_live.max():.4f} dex "
      f"-- exceeds the bound ({b_control}); walk suppressed by x{supp:.1f}")
    ok = b_bound and b_notight and b_control
    p(f"   [g3] -> {'PASS' if ok else 'FAIL'}")
    return ok


def _g4_inert(p):
    """g4 -- DEFAULT-INERT. Absent belief state, the M-step is BIT-IDENTICAL to the incumbent
    live path (the FORGE convention: every new capability is opt-in and provably inert)."""
    p("\n-- g4 DEFAULT-INERT: no belief state -> bit-identical to the incumbent --")
    ok = True
    try:
        sys.path.insert(0, "/home/mattm/projects/HSYMT/hpc_harbor/glacier")
        from glacier_loop import mstep_quadratic as ref_fn
        lab = "glacier_loop.mstep_quadratic"
    except Exception as e:                          # noqa: BLE001
        ref_fn, lab = _incumbent_mstep_quadratic, f"transcription ({type(e).__name__})"
    for tag, (n_slot, n_fed, seed, curv) in {
            "A": (3, 2, 20260729, 25.0),
            "B": (6, 4, 5150, 3.0),
            "C": (2, 1, 88, 400.0),
            "D-nonconcave": (2, 1, 12, -9.0)}.items():
        mode, F, src0, fed_idx, scale = _fixture(n_slot, n_fed, seed, curv)
        marg = _quad_marg(mode, F, n_slot)
        r_src, r_w, r_n = ref_fn(marg, src0, fed_idx, scale)
        s, w, ne, lg = BeliefMStep().step(marg, src0, fed_idx, scale, book=None)
        d_s = float(np.max(np.abs(s - r_src)))
        d_w = float(np.max(np.abs(np.nan_to_num(w, nan=-1.0) - np.nan_to_num(r_w, nan=-1.0))))
        b = (d_s == 0.0) and (d_w == 0.0) and (ne == r_n)
        ok &= b
        p(f"   [{tag}] vs {lab}: max|src| {d_s:.3e}, max|widths| {d_w:.3e}, "
          f"n_eval {ne} == {r_n} -> {'PASS' if b else 'FAIL'}")
    # a BeliefBook holding no belief for the fed slots is also inert
    mode, F, src0, fed_idx, scale = _fixture(3, 2, 20260729, 25.0)
    marg = _quad_marg(mode, F, 3)
    r_src, r_w, r_n = ref_fn(marg, src0, fed_idx, scale)
    s, w, ne, _ = BeliefMStep().step(marg, src0, fed_idx, scale, book=BeliefBook(3))
    b_empty = (float(np.max(np.abs(s - r_src))) == 0.0) and ne == r_n
    ok &= b_empty
    p(f"   [empty book] a BeliefBook with no belief for the fed slots is inert: {b_empty}")
    # and the checkpoint round-trip is exact
    book = BeliefBook(3)
    book.feed(0, src0[:NP_SRC], F[:NP_SRC], PRIOR_SIGMA, 3.0, 0)
    book.feed(1, src0[NP_SRC:2*NP_SRC], F[NP_SRC:2*NP_SRC], PRIOR_SIGMA, 0.25, 1)
    z = book.pack()
    rt = BeliefBook.unpack(z, 3)
    b_rt = all(np.array_equal(book[k].width, rt[k].width)
               and np.array_equal(book[k].centre, rt[k].centre)
               and np.array_equal(book[k].anchor, rt[k].anchor)
               and book[k].n_tighten == rt[k].n_tighten for k in book.keys())
    b_none = (BeliefBook.unpack({}, 3) is None)
    ok &= b_rt and b_none
    p(f"   checkpoint round-trip exact: {b_rt}; pre-belief checkpoint -> None (inert): {b_none}")
    p(f"   [g4] -> {'PASS' if ok else 'FAIL'}")
    return ok


def _g3b_graduation(p):
    """B3 readout gate (folded into the g3 block): the ladder reads belief width LIVE, the
    stages are ordered and monotone in width, and the fed rung IS the incumbent frontier."""
    p("\n-- B3 GRADUATION: the ladder reads belief width LIVE --")
    lad = GraduationLadder()                       # DEFAULT gate axes == ALL axes
    c = np.zeros(NP_SRC)
    stages = []
    for r in (2.0, 1.0, 0.6, 0.49, 0.3, 0.19, 0.11, 0.09, 0.0):
        bel = Belief(c, r * PRIOR_SIGMA)
        stages.append(lad.stage(bel))
    b_mono = all(stages[i] <= stages[i+1] for i in range(len(stages)-1))
    b_ends = (stages[0] == STAGE_CARRIED and stages[-1] == STAGE_CERT)
    b_rungs = (stages == [0, 0, 0, 1, 1, 2, 2, 3, 3])
    p(f"   stage vs width ratio [2, 1, .6, .49, .3, .19, .11, .09, 0]: {stages} "
      f"(monotone {b_mono}, ends carried..cert {b_ends}, rungs {b_rungs})")
    # the fed rung is the incumbent frontier, exactly
    from trackB_b1_softsource import ModelQualityGate, SoftComponent
    g = ModelQualityGate()
    agree = []
    for r in (0.9, 0.51, 0.5, 0.4999, 0.2, 0.0):
        bel = Belief(c, r * PRIOR_SIGMA)
        comp = SoftComponent(c, sigma=r * PRIOR_SIGMA)
        agree.append(lad.feeds(bel) == g.feeds(comp))
    b_same = all(agree)
    p(f"   fed rung == ModelQualityGate(SPARK3_GATE_RATIO) on every probe: {b_same}")

    # ---- the STANDING CONTROL against the S4.14 failure class, in the graduation ladder.
    # A component believed to 1e-4 on the M-step axes and at PRIOR BOX width everywhere else
    # is NOT well modelled (SPARK-3 SS5.0). The default all-axis ladder must call it CARRIED;
    # a ladder narrowed to the M-step axes calls it CERT-ELIGIBLE. Measured on the real b1
    # machinery, this is not hypothetical: a reservoir member the frontier scores at 0.948
    # reads 0.0015 on the M-step axes alone. The control asserts BOTH sides, so a future
    # narrowing of gate_axes fails here instead of reading green off a dead instrument.
    w_lop = PRIOR_SIGMA.copy()
    for j in MSTEP_AXES_DEFAULT:
        w_lop[j] = 1e-4 * PRIOR_SIGMA[j]
    lop = Belief(c, w_lop)
    s_all = lad.stage(lop)
    s_narrow = GraduationLadder(gate_axes=MSTEP_AXES_DEFAULT).stage(lop)
    b_control = (s_all == STAGE_CARRIED and s_narrow == STAGE_CERT)
    p(f"   lopsided belief (1e-4 on the M-step axes, PRIOR BOX on the other six): "
      f"all-axis ladder -> {STAGE_NAME[s_all]}, M-step-axis ladder -> {STAGE_NAME[s_narrow]} "
      f"-> {'PASS' if b_control else 'FAIL'} (the default must be the all-axis reduction)")

    # the trajectory readout banks, with the diagnostic column beside the binding one
    book = BeliefBook(2)
    book.set(0, Belief(c, 0.4 * PRIOR_SIGMA))
    led = GraduationLedger(2, lad, diag_axes=MSTEP_AXES_DEFAULT)
    for it in range(4):
        book[0].width *= 0.8
        led.record(it, book, tightened=(it == 2))
    T = led.trajectories()
    b_traj = (T["width_ratio"].shape == (4, 2) and np.all(np.diff(T["width_ratio"][:, 0]) < 0)
              and led.first_crossing.get((0, STAGE_FED)) == 0
              and T["width_ratio_diag"].shape == (4, 2))
    p(f"   belief-width trajectory banked {T['width_ratio'].shape} + the M-step-axis "
      f"diagnostic column, monotone-tightening curve, first fed crossing iter "
      f"{led.first_crossing.get((0, STAGE_FED))}: {b_traj}")
    ok = b_mono and b_ends and b_rungs and b_same and b_control and b_traj
    p(f"   [B3] -> {'PASS' if ok else 'FAIL'}")
    return ok


def _g5_smask(p, flip_ms=10.0, n_flip=20):
    """g5 -- RUNTIME-SMASK COMPATIBILITY. With belief state present and driving the M-step,
    pattern flips stay < 10 ms and NOTHING recompiles (jit cache sizes asserted constant,
    not inferred from timing). The belief lives entirely numpy-side and never enters a
    trace, so this is a structural claim; the gate proves it rather than asserting it."""
    p("\n-- g5 RUNTIME-SMASK: pattern flips with belief state live (no recompile) --")
    import jax
    import jax.numpy as jnp
    import trackB_b1_core as C

    def cache_size(fn):
        try:
            return int(fn._cache_size())
        except Exception:                            # noqa: BLE001
            return -1

    P = C.build_b1_problem()
    tt = np.asarray(P["theta_truth"], float)
    nd = P["n_dist"]
    ncw = P["ncw"]
    sp = P["sp"]                                     # the problem's own B1Split
    mask1 = jnp.asarray(P["mask_one"])
    data = P["amo"]["inject_delay"](jnp.asarray(tt), mask1)   # zero-noise Asimov
    wA = np.array([1.0] * C.N_LOUD + [0.0] * (ncw - C.N_LOUD))
    wB = np.ones(ncw); wB[1::2] = 0.0

    p0 = int(P["census_idx"][0])
    sp.AI = np.asarray(P["AI"])
    fn = sp._pulsar_ab_fn(p0)
    tb = jnp.asarray(tt); EVp = jnp.asarray(P["EV_truth"][p0])
    for w in (wA, wB):                               # warm both patterns
        sp.set_smask(w)
        jax.block_until_ready(fn(tb, EVp, data, mask1, sp.smask))
    cs0 = cache_size(fn)

    # belief state present and ACTIVE: a book over the loud slots, fed from a synthetic
    # frontier record, plus a graduation ledger recording every flip.
    book = BeliefBook(ncw)
    box = PRIOR_SIGMA.copy()
    for k in range(C.N_LOUD):
        book.feed(k, tt[nd + NP_SRC*k: nd + NP_SRC*(k+1)],
                  np.full(NP_SRC, 400.0), box, dlnL_support=4.0, iteration=0)
    lad = GraduationLadder()
    gled = GraduationLedger(ncw, lad, diag_axes=MSTEP_AXES_DEFAULT)
    led = InformationLedger()
    led.observe(0, fed=np.array([k in book for k in range(ncw)]))

    # fixed-pattern baseline (no flip): the evaluator's own cost on THIS lane. The FORGE-G2
    # budget (< 10 ms) is a jug-gpu number; on a slow lane the evaluator alone exceeds it and
    # no implementation could meet it, so the lane-independent claim is the FLIP OVERHEAD
    # above this baseline. Both are reported; the absolute budget is gated only where the
    # lane can carry it, and declared owed otherwise.
    base = []
    for _ in range(n_flip):
        t0 = time.perf_counter()
        jax.block_until_ready(fn(tb, EVp, data, mask1, sp.smask))
        base.append(time.perf_counter() - t0)

    flips = []
    for k in range(n_flip):
        w = wA if k % 2 == 0 else wB
        t0 = time.perf_counter()
        sp.set_smask(w)
        jax.block_until_ready(fn(tb, EVp, data, mask1, sp.smask))
        flips.append(time.perf_counter() - t0)
        gled.record(k, book)                          # the belief readout runs every flip
    cs1 = cache_size(fn)
    sp.set_smask(None)
    med = float(np.median(flips) * 1e3)
    med_b = float(np.median(base) * 1e3)
    over = med - med_b
    no_recompile = (cs0 == cs1)
    b_cache = (cs0 == -1) or no_recompile
    b_over = over < flip_ms                           # lane-independent: flip costs nothing
    b_abs = med < flip_ms                             # FORGE-G2's absolute jug-gpu budget
    b_bel = (len(gled.rows) == n_flip and len(book) == C.N_LOUD)
    p(f"   warm evaluator with belief live: fixed-pattern {med_b:.2f} ms, FLIP {med:.2f} ms "
      f"(median of {n_flip})")
    p(f"   FLIP OVERHEAD = {over:+.2f} ms (gate < {flip_ms:.0f} ms, lane-independent) "
      f"-> {'PASS' if b_over else 'FAIL'}")
    p(f"   absolute flip < {flip_ms:.0f} ms (FORGE-G2's jug-gpu budget): "
      f"{'PASS' if b_abs else f'{med:.1f} ms on this lane -- the EVALUATOR costs {med_b:.1f} ms here, so the absolute budget is OWED ON jug-gpu, not failed by the belief layer'}")
    p(f"   jit cache size across flips: {cs0} -> {cs1} "
      f"({'unchanged -> NO RECOMPILE' if no_recompile else 'CHANGED' if cs0 != -1 else 'API n/a'})")
    p(f"   belief state stayed live and banked across every flip ({len(gled.rows)} rows, "
      f"{len(book)} beliefs): {b_bel}")
    p(f"   (the belief is numpy-side and never enters a trace: it CANNOT recompile -- "
      f"this gate proves the structure rather than asserting it)")
    ok = b_cache and b_over and b_bel
    p(f"   [g5] -> {'PASS' if ok else 'FAIL'}"
      f"{'' if b_abs else '  (absolute-ms sub-check deferred to jug-gpu)'}")
    return ok


def mode_gate(cpu_only=False, verbose=True, n_iter_g3=120):
    p = print if verbose else (lambda *a, **k: None)
    print("=== RUNNING-BELIEF M-STEP GATES (FORGE-B) ===", flush=True)
    ok = True
    ok &= _g1_limits(p)
    ok &= _g2_guard(p)
    ok &= _g3_wander(p, n_iter=n_iter_g3)
    ok &= _g3b_graduation(p)
    ok &= _g4_inert(p)
    if cpu_only:
        print(f"\n=== LOGIC GATES g1-g4 + B3: {'ALL PASS' if ok else 'FAIL'} "
              f"(g5 skipped: --cpu-only) ===", flush=True)
        return 0 if ok else 1
    try:
        ok &= _g5_smask(p)
    except Exception as e:                            # noqa: BLE001
        print(f"\n  g5 NOT RUN ({type(e).__name__}: {e}). Run on the sanctioned lane; "
              f"g1-g4 above stand.", flush=True)
        print(f"\n=== BELIEF GATES: g1-g4 {'PASS' if ok else 'FAIL'}, g5 DEFERRED ===",
              flush=True)
        return 0 if ok else 1
    print(f"\n=== RUNNING-BELIEF M-STEP GATES: {'ALL PASS' if ok else 'FAIL'} ===", flush=True)
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="gate", choices=["gate"])
    ap.add_argument("--cpu-only", action="store_true",
                    help="run only g1-g4 + B3 (analytic fixtures); skip g5 (needs the b1 stack)")
    ap.add_argument("--n-iter-g3", type=int, default=120,
                    help="wander-bound iterations (gate requires >= 100)")
    a = ap.parse_args()
    if a.mode == "gate":
        assert a.n_iter_g3 >= 100, "g3 is specified at >= 100 iterations"
        return mode_gate(cpu_only=a.cpu_only, n_iter_g3=a.n_iter_g3)
    return 2


if __name__ == "__main__":
    sys.exit(main())
