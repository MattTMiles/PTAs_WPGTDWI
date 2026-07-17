# SPARK-2 — the second arrow (certification-side feedback): **VERDICT NOT REACHED**

**Status: INCOMPLETE AND STOPPED. Neither EDGE-POSITIVE nor EDGE-ZERO is claimed.** The
(realisation × rung) grid is **2 of 12 units** banked; 3 units died on a numerical pathology
diagnosed but not yet fixed, and the rest were cancelled rather than left running past the
point where I could analyse them. **Do not read the indicative numbers in §3 as the ledger.**

What SPARK-2 *did* deliver is three measured, banked, reusable results — two of which
change how arrow 2 must be built — plus the trials-factor honesty item the addendum asked
for. They are in §1–§2 and they stand on their own.

**Tree.** `MM_playground` @ `d87db93`. Arrow 1 = `reports/SPARK_launch_criterion.md`
(CASCADE-ALIVE). Cell: h = −13.25, T = 30, lit, fiducial POP sky. L_gwb fingerprint
**`9fd547b39b02c705`** (the banked CHORUS T=30 value) on every job; `cpus-per-task=8`.
**MOCK spine** (AXIS 1440 MHz, §10.15(a)) — the residuals *are* the injected CW+CURN.

**SCOPE OF INFERENCE.** §1 and §2 are statements about the machinery and about exact,
CPU-checkable counting facts. §3 is one realisation at two rungs and is **not** a measurement
of arrow 2.

---

## 1. THE STRUCTURAL FINDING — arrow 2's "rigidity" question cannot be asked the way arrow 1's was

The addendum asked whether SPARK's **selectivity channel** (the null floor falling because the
coherent template is rigid — arrow 1's second and larger mechanism) also operates on the
certification side. **It cannot be tested by masking the certification E-step, and this is now
measured rather than argued** (`SPARK2_results/spark2_mask.npz`, mode `mask`):

| pmask rung | LIVE fringe rows | offender | median qmax |
|---|---|---|---|
| **0** | **0 / 116** | 0.000 | 0.0020 |
| **8** | **8 / 116** | 4.447 | 0.0020 |
| **116** (pmask = ONE) | **116 / 116** | 4.435 | 0.3506 |

`B1Split.estep(theta_base, EV, AI, data, pmask)` sweeps **every** pulsar's distance in **one
call under a global pmask**. A pulsar at `m_p = 0` has its *own* pulsar term switched off, so
its distance is inert, its fringe row is **flat**, and it cannot certify at all. `n_live` tracks
the rung **exactly** — 0, 8, 116. So an E-step at the rung's global pmask kills the uncertified
rows outright, and the offender statistic (a **max over pulsars**) then falls for a **counting**
reason, not template rigidity. **At rung 0 the certification statistic does not exist.**

> **The certified-coherent E-step arrow 2 needs — target pulsar's own term ON, other certified
> pulsars at q, uncertified OFF — is one call PER TARGET, not one call.** That is a ~116× cost
> and a real build, not a parameter. It is the single largest thing standing between this
> addendum and its verdict.

Note also `offender(rung 8) = 4.447 ≈ offender(rung 116) = 4.435`: the offender is a max, and
here it is set by one pulsar that is already in the rung-8 set. **A max-over-pulsars statistic is
nearly blind to the rung** — a second reason the naive test is uninformative.

(b) therefore keeps **pmask = ONE**, the convention every banked IGNITE/CHORUS/SURFACE number
uses (`ignite.py:349`).

---

## 2. THE TRIALS-FACTOR HONESTY ITEM — the cost column, de-oracled

From SPARK §5(a): arrow 1's cost column was oracle-flattered. Measured exactly on CPU
(`spark2_trials.npz`, mode `trials`). `dL` = **min** over the **modelled** source list of the
mode spacing (`forge_b1.apply_geometry:91-94`), so modelling the reservoir shrinks `dL` → more
fringes in the prior window → `K_counted` and the `lnK` bar **grow**. This is what a
*deconfusing* loop pays for the privilege:

| n_src modelled | K_sum | lnK_med | Δ vs 3 loud |
|---|---|---|---|
| 1 | 348 330 | 7.164 | −1.019 |
| 3 (the loud) | 581 900 | 8.183 | — |
| 4 (+1 recruited) | 622 702 | 8.291 | **+0.108** |
| 5 | 681 814 | 8.341 | +0.158 |
| 16 (full census) | 993 132 | 8.761 | **+0.578** |

**The cost of modelling the whole reservoir is +0.578 nat on the median pulsar's trials bar, and
it SATURATES** — the curve is flat from n_src ≈ 13 onward (8.761 at 13, 14, 15 and 16). The
first recruit costs +0.108 nat; the thirteenth costs ~0. **Deconfusion's trials cost is
sublinear and bounded**, which is favourable to arrow 2 and is now a measured number rather
than an oracle anchor.

---

## 3. WHAT THE TWO SURVIVING UNITS INDICATE — **not a verdict**

`s2b_r0_k0` and `s2b_r0_k8`: **one** realisation, **two** rungs, against a pre-registered
3 × 4. Quoted only to say where the effect is going and to justify §4's fix list.

| unit | rung | coherent | med σ_cond(faint) | offender: unmodelled | modelled(soft) | floor_un (zf) | floor_mo (zf) | N_cert un | N_cert mo |
|---|---|---|---|---|---|---|---|---|---|
| r0_k0 | 0 | 0 | 0.1708 | 0.000 | 0.000 | 0.000 (0.98) | 6.681 (0.66) | 0 | 0 |
| r0_k8 | 8 | 8 | **0.1293** | 0.000 | 0.000 | 0.000 (0.98) | **5.426** (0.78) | 0 | 0 |

Three things are visible, none of them a verdict:

**(a) tightening is real but small.** Median conditional width **0.1708 → 0.1293 (1.32×)** for
8 coherent pulsars. Consistent with arrow 1's finding that the per-pulsar lever is roughly
linear and unspectacular at reachable N_cert.

**(b) the deconfusion ladder is the right object, and it is steep.** With the reservoir
**unmodelled**, the signal offender is **0.000** and *nothing certifies* — while the pilot's
**truth**-modelled control on the same machinery gave **4.435**. The 13 faint sources are not a
perturbation on certification; unmodelled, they erase it. That is a large deconfusion effect and
it is arrow 2's premise looking healthy. **But** a soft model at rung-0 *or* rung-8 widths also
gives **0.000** — i.e. at reachable rungs the reservoir cannot yet be modelled well enough to
recover any of it. The ladder to measure is `unmodelled (0) → soft@rung k (?) → truth (4.435)`.

**(c) the floor does fall with the rung** — 6.681 → 5.426 in the modelled state — but `N_cert`
stays **0** because the signal offender is 0. On this evidence arrow 2 would read EDGE-ZERO, **and
one realisation at two rungs is not evidence.** The zero-fractions (0.66–0.98) are far above
RECUT's `ZF_GATE = 0.20`, so these floors are emp-q95 + bootstrap, correctly — and a 98%
zero-fraction floor of 0.000 is a floor in name only.

---

## 4. WHAT BROKE, AND THE FIX LIST

**(a) The joint/marginal Fisher is not affordable — measured twice, not assumed.** (a) was
specified as *joint* constraints. The marginal width needs the full 78×78 Hessian of lnL over
the faint block. **Two builds failed to return on an A100-80GB inside a 1 h walltime:** a
monolithic `jax.hessian` (job 12583525, >27 min on the Hessian, cancelled) and a chunked
batched-JVP-of-grad Hessian at CH=8 (job 12583899, >17 min, cancelled). The likelihood carries a
2668-dim correlated-GP solve per evaluation; 78 JVPs of its gradient is not a few-GPU-hr object.
**SPARK-2 fell back to the CONDITIONAL width** (1-D curvature scan per axis, one batched `logL`
call per chunk, seconds) **and labels it as such on every number.** The bias direction is stated
and is safe: `σ_cond ≤ σ_marg` always, so the conditional width **overstates** how well the
reservoir is known and therefore **favours** arrow 2 — an EDGE-ZERO measured this way is a
fortiori EDGE-ZERO with true joint widths. **An EDGE-POSITIVE measured this way would not be
safe** and would have to be re-cut against the joint width before being quoted.

**(b) The soft draw is numerically pathological — partially fixed, not closed.** The MODELLED
state draws the faint sources from `N(truth, σ(rung))`. Unclipped, this put `cos_inc` outside
[−1, 1] (its width at loose rungs is comparable to its whole box), the waveform went NaN, and the
E-step returned **NON-FINITE on all 59 392 entries — all 12 units died** (12584400–11). Clipping
the draw to `TE.phi_bounds` fixed rungs 0 and 8 but **rungs 2 and 5 still die all-NaN**
(12584896/97, 12584901). This is IGNITE-2/KINDLE's known double-perturbation pathology
(`kindle-gain-diagnosis`: *"scrambled + fix_mc is a double perturbation … degenerate covariance →
all-non-finite E-step"*) arriving by a new route: a wide draw plus the fringe sweep. Remaining
suspect: `log10_mc` high in its 0.87-dex box × the faint sources' own `fgw` → coalescence inside
the observation. **EMBER's `mc_scan` measured "non-evaluable 0% across the entire mc box" for the
LOUD sources — that result does not transfer to the faint, and SPARK-2 is the counter-example.**
**Required fix: reject-and-redraw on non-finite, with the rejection count BANKED** — the fraction
of the faint posterior that is non-evaluable is itself a measurement about how soft "soft" can be,
not a nuisance to suppress.

**(c) Per-job JAX cache isolation costs a cold build every time.** `HARBOR_JAXCACHE` keyed on
`SLURM_JOB_NAME` gives each job an empty cache: builds ran 302 s (warm) to 606 s (cold), and the
first E-step compile is 315–500 s. Correct (it defuses EMBER §7(a)'s race) but expensive at 12
one-shot units. **Fix: one task per REALISATION looping the 4 rungs in-process** — 3 builds
instead of 12, which is what the addendum's "many units per process" asked for and which I did
not do.

---

## 5. TO FINISH ARROW 2 (the successor's list, in order)

1. **Reject-and-redraw the soft faint draw**, bank the non-evaluable fraction per rung (§4b).
2. **One task per realisation, 4 rungs in-process**, skip-on-exist per unit (§4c).
3. **Add the `truth`-modelled control to every unit** — the ladder `unmodelled → soft@rung →
   truth` is the measurement, and §3(b) shows the endpoints are 0.000 and 4.435. Without the
   truth rung the ledger has no scale.
4. **Then, and only if the ledger reads EDGE-POSITIVE**, pay for the per-target E-step (§1, ~116×)
   to test the rigidity channel honestly, and re-cut (a) against true joint widths (§4a).
5. `N_cert = 0` in **both** states at both rungs means the current ladder may be **entirely below
   the certification bar at this cell** — check a louder cell or the truth control first, before
   spending on a grid that can only return zeros.

**Arrow 1's verdict is unaffected.** CASCADE-ALIVE stands, and with arrow 2 unmeasured the
EM-mediated **EPOCH** route (SPARK §6) remains the only pre-registered campaign. **GLACIER — the
soft-joint fixed-list loop — is NOT pre-registered here**, because its licence was an
EDGE-POSITIVE verdict and no verdict was reached.

---

## 6. FILES

```
hpc_harbor/spark/spark2.py      modes: pilot | mask | unit | ledger | trials
hpc_harbor/spark/sp2_one.sbatch parameterised runner (SP2_ARGS)
hpc_harbor/spark/rb2.py         readback of the surviving units + mask + trials
SPARK2_results/spark2_mask.npz     §1, the global-pmask confound (COMPLETE)
SPARK2_results/spark2_trials.npz   §2, the trials cost column (COMPLETE)
SPARK2_results/s2b_r0_k0.npz       §3, indicative only
SPARK2_results/s2b_r0_k8.npz       §3, indicative only
```

`SPARK2_results/` is gitignored (`*_results/`, the repo convention — no peer result dir is
tracked); the banks live on ACCRE. Lean-npz throughout: raw per-pulsar `dlnL`/`lnK`/`qmax` and
raw null offender vectors are banked per unit, never verdicts.

## STOP

**Arrow 2 is NOT measured. No verdict — not EDGE-POSITIVE, not EDGE-ZERO.** §1 (the E-step's
global pmask makes the rigidity question a per-target, ~116× build — and at rung 0 the
certification statistic does not exist) and §2 (deconfusion's trials cost is +0.578 nat and
saturates) are complete and stand. §3 is one realisation and is not the ledger. §4 lists two
defects I diagnosed but did not close and one I introduced. Awaiting Matt's decision before any
further spend.
