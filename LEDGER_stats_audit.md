# LEDGER_stats_audit.md — agent LEDGER, ACCRE, 2026-07-29

Statistical-treatment audit of the loop's inference layer. **Report-only: nothing here
arms a protocol step, moves a banked verdict, or enters a closure claim.** Every new rule
is **pre-registered forward**; the two banked re-cuts (B1, B2) are re-*readings* of raw
banked columns and are labelled as such wherever a campaign number moves.

Charter posture: operational items proceeded; the one budget/lane-class event (B3's GPU
leg) is **POSTED AND PARKED**, not queued behind another session. `LEDGER` stages;
**Matt commits**.

---

## STEP 0 — SYNC, AND WHY `git pull --ff-only` DID NOT RUN

```
branch            glacier_lite     (no upstream configured)
origin/master     10 behind / 49 ahead of HEAD
```

`git pull --ff-only` **cannot execute on this branch**: `glacier_lite` has no upstream, and
against `origin/master` the histories have diverged 49/10, so no fast-forward exists.
`git fetch origin` ran (clean). **Nothing was merged, rebased, reset, or forced** — that is
charter #4 territory and it is Matt's. The 10 unmerged upstream commits are notebook /
README / sampler work (`c0c8294 cleanup`, `5e0de21 CW likelihood scan notebooks`,
`0d3ac6d …sampler to find pulsar distances`, …) and touch nothing this audit reads.

**Concurrency note.** Two other sessions are live and were read before acting:
PHOENIX (frozen-vs-live M-step arm, pinned dgx03 A100-80) and BASELINE (field baseline,
pinned dgx01 A100-40). A third, **SIEVE-A**, has since written
`hpc_harbor/ledger/sieve_t7_eprocess.py` into the directory this audit created and cites
"LEDGER B2 / criterion-v2.3" — i.e. it is already consuming §B2 below. That file is
**not** in LEDGER's add-list.

---

## THE VERDICTS, IN ONE TABLE

| item | verdict | one line |
|---|---|---|
| **A1** E-step source-conditioning | **REAL** | the fringe posterior is evaluated at ONE source point; the loop computes a belief and discards it (`glacier_loop.py:417,496` → never reaches `:673`) |
| **A2** noise-leak / RN refit | **REAL (frozen)** — and the data shares the model's RN exactly | RN hyperparameters land in `frozen` (`trackB_b1_core.py:256`) and can never move; `NoiseDrawer` draws from the *same* `Phi_rn_diag` |
| **B1** drain winner's curse | **PARTIAL, and two-sided** | the e07 arm's bites are ~entirely a feed-free fit offset; the *none* arm's bites are **understated**, not overstated |
| **B2** iteration multiplicity | **REAL** | 6 scored looks, priced nowhere; the campaign's own true certs are non-consecutive flickers |
| **B3** feed-order dependence | **NOT where it looks** | feed *decisions* are order-invariant; the M-step's coordinate ascent is not. Venue leg **PARKED** on lane; CPU leg done and it refutes the obvious hypothesis |
| **C1** belief multimodality | spec delivered (§C1) | posted here, not in the capstone — that file is frozen pending Matt's GLACIER-LITE review |
| **C3** q-tail payoff | **REAL** | SIREN's 6–12 % is the width of the all-right component of a mixture the criterion itself implies |

---

## A1 — E-STEP SOURCE-CONDITIONING: **REAL**

### The finding, at the lines

The E-step is evaluated at a single source point and nowhere else:

```
glacier_loop.py:668   theta_base = theta_with_absent(theta_rec, nd, carried)
glacier_loop.py:673   LNL = estep_per_target(self.sp, theta_base, EV, AI, data, PM, jnp)
spark3.py:269         tb = jnp.asarray(theta_base)      <- one point
spark3.py:273-289     the only loop is over TARGET PULSARS
trackB_b1_core.py:700 FringeTables.posterior(LNL)       <- q_max, map_k off that one surface
```

The loop **does** carry a source belief, and throws it away at exactly this step:

| where the belief is made | what happens to it |
|---|---|
| `glacier_loop.py:417` `fisher_conditional` → `sig_opt` | used only for the frontier ratio (`:420`), banked (`:555`) |
| `glacier_loop.py:496` `mstep_quadratic` → `widths` | **not even banked** |
| `forge_b1.py:311-316` `Cov = inv(Fs + Pi)` → `sig_f, sig_mc` | `forge_b1.py:331` re-E-steps at `theta_src`, the point |

**Consequence.** `q_max` is `P(fringe | this exact source)` and is charged against the bar as
if it were `P(fringe | data)`. This is the same error class the programme has already
retired once: IGNITE-2's "hard lock" was over-confidence in a **pinned fringe**; this is
over-confidence in a **pinned source**.

### How much it matters — measured, not asserted

`ledger_a1_sigma_estep.py gate`, leg G-L3b, on a fringe toy whose peak position moves with
the source (which is the actual coupling — pulsar-term phase depends on `(fgw, mc)`):

| peak motion per 1σ of belief | q_max point → marginal |
|---|---|
| 0.00 fringe | 0.9094 → 0.9094 |
| 0.10 | 0.9094 → 0.8944 |
| 0.25 | 0.9094 → 0.8113 |
| 0.50 | 0.9094 → **0.5443** |
| 1.00 | 0.9094 → 0.3224 |

**At half a fringe of belief-induced peak motion a q = 0.91 certification becomes q = 0.54 —
it crosses the QBAR = 0.9 bar from the wrong side.** How much motion the real venue has is
the venue's own number and is not claimed here.

### The upgrade (pre-registered FORWARD, for the belief arm's fan only)

Replace the point evaluation with a quadrature of the **likelihood** over the belief:

```
LNL_marg[p,b] = ln ∫ dθ_s N(θ_s; θ̂, Σ) exp(LNL(θ_s)[p,b])  ≈  logsumexp_i(ln w_i + LNL(χ_i)[p,b])
```

Note the order: average in **likelihood** space, then take the log. Averaging
log-likelihoods gives the geometric mean and cannot widen a posterior at all.

Rule: scaled unscented transform, `w0 = 1/3` (Julier–Uhlmann's Gaussian value) so **every
weight is strictly positive** — a likelihood quadrature requires it, since a negative weight
can drive the integrand negative and the log undefined. `m = 2·n_fed` over `MSTEP_AXES`,
which is exactly FORGE-B's `(sig_f, sig_mc)` block.

### Gates — **ALL PASS** (`hpc_harbor/ledger/ledger_a1_sigma_estep.py gate`)

| gate | result |
|---|---|
| **G-L1 belief-width → 0 reproduces the incumbent** | `max\|Δ\| = 0.000e+00` **exactly** |
| G-L2 degenerate rule (`w0 = 1`) at full belief width | `0.000e+00` exactly |
| G-L2b weights strictly positive, sum to 1 | min w 0.0833, sum−1 = 0 |
| G-L2c point-set moments match (θ̂, Σ) | 4.4e-16 / 1.4e-17 |
| G-L3 vs the analytic Gaussian marginal | captures 99.5 % of the shift |

**One honesty note on G-L1.** The bit-exactness is *structural*, not numerical: `sigma_points`
collapses to the single incumbent point when `Σ = 0`, the same way `embed_igniter` returns
the census bit-exactly at `e = 0` (CHORUS C1). Computed *through* the logsumexp with
unit-sum weights the same quantity lands **1 ulp (2.2e-16)** off — which is why the collapse
is written explicitly rather than "bit-comparable" being a tolerance in disguise.

### The cost, and the cap this forces

`n_eval = 4·n_fed + 1` per E-step: **5× at n_fed 1, 21× at 5, 81× at 20.**
The r13p25-class cells are **not affordable**. The fan therefore caps the sigma set at the
top `n_belief` fed members by (belief width / prior box width) — the members the loop is
least sure of, i.e. the ones that can actually move a fringe — with the rest held at point.
`n_belief` is banked per cell and **`n_belief = 0` IS the incumbent.**

### Wiring for the belief arm (**NOT applied — see the caveat**)

```python
# in the belief arm's driver, replacing CertScoreboard.columns' single estep call:
from ledger_a1_sigma_estep import estep_sigma, rank_belief_members
sel   = rank_belief_members(mstep_widths, box_sigma[list(MSTEP_AXES)], fed_idx, n_belief)
axes  = np.array([nd + NP_SRC*k + j for k in sel for j in MSTEP_AXES])
sigs  = np.array([w for k, w in zip(sel, mstep_widths[sel]) for _ in MSTEP_AXES]).ravel()
LNL, n_eval = estep_sigma(lambda th: SP3.estep_per_target(sp, th, EV, AI, data, PM, jnp),
                          theta_base, axes, sigs)
```

**Caveat, stated plainly: I did not apply this, and I did not apply the A2 monitor wiring
either.** `glacier_loop.py` is the HELD driver and **PHOENIX's frozen arm (`frzgate`,
job 12833771) is running against it right now**; the frozen arm's whole value is
bit-comparability with the banked live arm, and editing the file underneath a live fan
risks splitting that. The modules and gates are delivered; the two wiring diffs are Matt's
to apply between fans.

---

## A2 — NOISE-LEAK: **REAL (frozen)**, and the idealisation is stronger than "frozen"

### The finding, at the lines

```
glacier_pop.py:530-535 (docstring)  "RN frozen at the seed-0 reference draw … fixed noise
                                     covariance across skies"
trackB_b1_core.py:252   theta_keys = dist_keys + cw_param_keys(num_cw)
trackB_b1_core.py:256   frozen = {k: temp_dict[k] for k in truth_full if k not in theta_keys}
                        -> every *_rednoise_log10_A / *_rednoise_gamma lands in `frozen`
trackB_b1_core.py:357-8 Pinv, ldP = vsm.P_var.make_inv()(self.frozen);  Phi_rn_diag = 1/Pinv
```

RN is **not in θ**, so no stage — Fisher, frontier, drain, E-step, M-step, scoreboard — can
move it, and the per-pulsar Cholesky built from it is cut once per cell.

**And it is stronger than frozen.** `NoiseDrawer` draws the data's RN from the *same*
object (`trackB_b1_core.py:590,609`: `sig_rn = sqrt(sp.Phi_rn_diag)`). Model RN and data RN
share hyperparameters **exactly**. Named plainly: **the campaign never pays for not knowing
the noise.**

**Already-handled, and credited:** ANCHOR §4 measured that deliberate RN mis-specification
does **not** move the certification floor (the tail is template-dominated). So this is *not*
a floor defect, and A2 is not a retraction of anything. The exposure it leaves is different:

> Feeding a source removes power from the residual; the M-step then wanders the fed member
> off truth by **0.01–0.44 dex** (S4.15, measured). That mismatched power has to go
> somewhere. `BackgroundFit` absorbs the common-mode part into the fitted GWB amplitude —
> but the **per-pulsar** part has nowhere to go except the residual, where it inflates the
> effective noise the fringe likelihood is scored against. Frozen RN cannot show that as a
> parameter moving. It shows up only as posterior RN power drifting from its prior.

### The monitor (added; additive, and it cannot change a banked verdict)

Exact Gaussian conditional posterior on the per-pulsar RN GP coefficients — **every object
it needs is already held by `B1Split`**, so it costs one 1-D solve + one matvec + one
triangular solve per pulsar per iteration. **No new factorisation, no extra likelihood call,
and no RNG draw** (so it cannot perturb a realisation and a resumed checkpoint is unaffected).

```
S_p = (1/n_k) Σ_k ( E[c_pk]² + Cov[c_p]_kk ) / Φ_pk        the whitened RN occupancy
Var(S_p) = (1/n_k²) Σ_k ( 2 Cov_kk² + 4 E[c_pk]² Cov_kk ) / Φ_pk²      (analytic; no MC)
```

**The band and flag, pre-registered here, forward:** reference `S_p(0)` in the same cell
(self-referenced — cross-venue RN occupancy is not comparable); band
`sqrt(Var(0)+Var(i))`; **FLAG** at `|S_p(i) − S_p(0)| > 3·band` per pulsar; **CELL FLAG** if
any pulsar flags **or** the array-median drift exceeds 3× its own band — the common-mode leg,
which is precisely what a per-pulsar 3σ test misses by construction.

**Report-only at introduction.** It does not STOP a cell and does not enter any
certification. Promotion to a STOP is a pre-registered bar change = charter #1 = Matt's.

### Gates — **ALL PASS** (`ledger_a2_rn_monitor.py gate`)

| gate | result |
|---|---|
| G-R1 zero-residual occupancy finite, below prior (shrinkage) | S ∈ [0.0210, 0.0401] |
| G-R2 single-pulsar RN-band leak flags **that pulsar only** | +25.9σ, flagged `[2]` |
| G-R3 unchanged state raises **no** flag (false-positive leg) | 0 flagged, cell flag False |
| G-R4 common-mode leak caught by the **median** leg | cell flag True |

### Venue replay — **SUBMITTED** (CPU lane, job `12834356`, `ldg_a2_replay.sbatch`)

Rather than edit the live driver, the monitor is run **post hoc** on banked cells:
`theta_rec` and `fed_mask` are banked per iteration (`glacier_loop.py:566,557`), so the whole
history can be re-scored after the fact. Deliberately CPU: `rn_occupancy` never calls
`estep_per_target`, so the CPU-lane `vm.max_map_count` hazard does not apply, and requesting
no GRES means it **cannot contend with PHOENIX or BASELINE**. Leg 0 is the synthetic gate
(hard fail); leg 2 is Stage-1 `none/s1` at T=30/n=256 and is **attempted, not promised** — a
CPU build at that venue is unmeasured in this programme.

### RESULT — the monitor RAN on a production fan cell, and it is QUIET

Job `12834356` (`cn1461`, `lane_tag = cn1461_noGPU`). Leg 0 passed on-lane (all four gates
green **inside the container**, not just on the login node). **Leg 2 completed** — the
T=30/n=256 venue build, which this programme had never measured on CPU, works — and produced
the full 6-iteration monitor for `gl1_cell_none_s1_T30_lit`:

| iter | median S_p | max S_p | median drift vs i0 | flagged | CELL FLAG |
|---|---|---|---|---|---|
| 0 | 1.02465 | 1.336 | — (reference) | 0/116 | False |
| 1 | 1.02144 | 1.655 | +0.00059 | 0/116 | False |
| 2 | 1.02569 | 1.569 | +0.00206 | 0/116 | False |
| 3 | 1.02843 | 1.342 | +0.00268 | 0/116 | False |
| 4 | 1.02860 | 1.555 | +0.00155 | 0/116 | False |
| 5 | 1.03218 | 1.429 | +0.00110 | 0/116 | False |

**Two readings.**

1. **Median `S_p` ≈ 1.02–1.03 — occupancy sits at the prior.** That is A2's finding shown
   directly rather than argued from source: the data's RN and the model's RN come from the
   same `Phi_rn_diag`, so the whitened occupancy is 1 by construction. It also makes the
   G-R1 synthetic result (`S < 1`, pure shrinkage at zero residual) read correctly — that
   was a *no-RN-in-the-data* case, and real data lifts it to 1.
2. **No leak in this cell.** 0/116 flagged and CELL FLAG False at every iteration; drift is
   +0.0006 to +0.0027 and **non-monotonic** (peaks at i3, falls at i4–i5), which is the
   signature of noise rather than accumulation. This is the G-R3 false-positive leg passing
   on real data — the monitor is quiet when nothing has happened, which is what makes a
   future flag worth reading.

**Scope, and it matters: this is a scope-limited null.** `none/s1` feeds **2** members and
carries a −0.176 dex bite. The cells where the leak hypothesis has teeth are the r13p25/e07
cascade cells — **n_fed 20–23 with M-step wander up to 0.44 dex** — and those are exactly the
cells this leg has not touched. The honest statement is *"the frozen-RN idealisation is not
leaking measurably where the loop barely moves"*, not *"there is no leak"*. Running the same
replay on one cascade cell is the obvious next step and is cheap (same CPU lane, one build).

**One defect found and fixed in my own code, recorded rather than quietly patched:** the
first submission computed all six iterations correctly and then died on the bank write —
`glacier_pop.bank_npz`'s own first parameter is named `stem`, so passing the cell's stem as
`stem=` collided (`TypeError: bank_npz() got multiple values for argument 'stem'`). Renamed
to `cell_stem=`; re-gated (PASS) and resubmitted as job **`12843673`**. The numbers in the
table above are from the completed first run's log and will be reproduced bit-for-bit by the
resubmission — the fix touches only the write, not one line of the arithmetic.

### A2-SPEC — refit-in-loop (**spec only, not implemented**, per the brief)

Refitting per-pulsar RN inside the loop is a design change, and the reason it is not a free
improvement is worth recording: **it would move every banked floor.** Every `recut_surface`
floor in the programme was cut against a fixed noise covariance; a fitted RN makes the null
offender distribution a function of the template state, which is exactly the circularity the
floor exists to break. If it is ever done, the pre-registration must be: (i) refit RN in the
**null draws too**, or the floor is cut against a different model than the data; (ii) profile,
not sample, the per-pulsar amplitude via the same band-selective diagonal-rescale identity
`BackgroundFit` already uses for the GWB (`glacier_pop.py:393-430`) — one extra
`ncomp × ncomp` Cholesky per pulsar per grid point, affordable; (iii) a fresh floor re-cut
for every cell, with the old floors retired rather than reused. **Cost is dominated by (iii),
not by the fit.**

---

## B1 — THE DRAIN, RE-REFERENCED AGAINST THE SCRAMBLED ARM

`hpc_harbor/ledger/ledger_b1_drain.py` → `reports/ledger_b1_drain.npz`. 706 per-iteration
banks, pure numpy, banked columns read and never rewritten.

### What the campaign has been quoting, and against what

`bite(cell,i) = a_bg(i) − a_eff_drawn`. Its implicit null is **zero**. That null is wrong two
ways, and both push the same direction: (1) `a_bg` is a *fitted* amplitude of a band-limited
GP against a finite noisy realisation while `a_eff_drawn` is an *exact* projection — they
differ **before anything is fed**; (2) the spoken headline is a **max over a fan**, whose
reference is the null's max over an equal-sized fan, not zero.

### The control that already existed

The scrambled arm shares the census **and, at `real = 0`, the noise seed exactly**
(`9_500_000 + 1000·sky`). Only the recovery template's igniter sky differs. Since S4.15.1's
amended null semantics it feeds real members too — so it produces its **own feed-drain**.
That is the missing reference distribution, and it was already banked.

### (1) The feed-free offset is large

**Bite at iteration 0 with ZERO fed members, n = 51: mean −0.389, sd 0.374 dex.**
Per venue:

| venue | zero-fed baseline (n) | max bite in fan | at n_fed |
|---|---|---|---|
| gl1 stage-1 / e07 | **−0.024** ± 0.015 (3) | −0.222 (s7) | 4 |
| gl1 stage-1 / none | **+0.017** ± 0.095 (4) | −0.176 (s1) | 2 |
| gl2 r13p9 / e07 / T30 | **−0.501** ± 0.099 (4) | −0.586 (s3) | **0** |
| gl2 r13p9 / e07 / w0.25 / T40 | **−0.704** ± 0.071 (10) | −0.788 (s1) | **0** |
| gl2 r13p25 / none | **+0.193** ± 0.032 (5) | −0.768 (s0) | 1 |
| gl2 r13p9 / none / T30 | **+0.235** (1) | −0.573 (s2) | 1 |

### (2) Paired excess — the exact contrast (same census, same noise, template scramble only)

| venue | sky | bite_sig | bite_null | **excess** | band | n_fed s/n |
|---|---|---|---|---|---|---|
| gl1 stage-1 / e07 | 0 | +0.061 | +0.061 | **+0.000** | 0.033 | 1/1 |
| gl1 stage-1 / e07 | 1 | −0.217 | −0.217 | **+0.000** | 0.056 | 1/1 |
| gl1 stage-1 / none | 1 | −0.176 | −0.035 | **−0.141** | 0.048 | 2/1 |
| gl2 r13p25 / e07 | 0 | −0.656 | −0.553 | −0.103 | 0.031 | 20/18 |
| gl2 r13p25 / none | 0 | −0.768 | +0.162 | **−0.929** | 0.029 | **1/0** |
| gl2 r13p25 / none | 1 | −0.635 | +0.194 | **−0.829** | 0.032 | **1/0** |
| gl2 r13p9 / e07 / w0.25 / T40 | 0,1 | −0.711/−0.788 | identical | **+0.000** | — | 0/0 |
| gl2 r13p9 / none / w0.25 / T40 | 0,1 | −0.707/−0.877 | −0.231/−0.342 | **−0.48/−0.54** | 0.06 | 1/0 |

7/12 exceed 2× their own band.

### (3) The winner's-curse statistic proper

| venue | n_sig | max_sig | E[max_null \| n_sig] | sd | **excess** |
|---|---|---|---|---|---|
| gl1 stage-1 / e07 | 8 | −0.222 | **−0.220** | 0.005 | **−0.002** |
| gl1 stage-1 / none | 8 | −0.176 | −0.047 | 0.009 | −0.129 |
| gl2 r13p25 / e07 | 4 | −0.677 | −0.553 | 0.000 | −0.124 |
| gl2 r13p25 / none | 4 | −0.768 | +0.165 | 0.008 | **−0.933** |
| gl2 r13p9 / e07 / w0.25 / T40 | 4 | −0.788 | −0.780 | 0.019 | −0.008 |

### THE FINDING — the curse is real, and it runs in **opposite directions** in the two arms

**e07 arm.** Its biggest "bites" occur in cells that fed **nothing** (`r13p9`: max bite
−0.586 to −0.788 at **n_fed = 0**), and its Stage-1 headline is reproduced by the null to
0.002 dex. S4.17 already footnotes the mechanism — `a_eff_drawn` is computed *before* the
comb redistributes the bright member's power across 31 harmonics, partly out of band — and
that footnote is **correct and sufficient**; what this re-reference adds is the *size*
(−0.50 to −0.70 dex) and the fact that it swamps everything the loop does in that arm.
**S4.15's headline "e07_s7 … a real −0.22 dex first bite" does not survive**: the null arm's
max over an equal fan is −0.220 ± 0.005.

**none arm.** The opposite. Its zero-fed baseline is **positive** (+0.017 to +0.235: the
in-band GP legitimately absorbs some white+RN), so quoting against zero **understates** the
drain. Re-referenced, feeding a *single* member at r13p25 drains **−0.93 ± 0.03 dex** — the
cleanest measurement of "the array eats the background" the campaign has, and it is bigger
than what was reported, not smaller.

### Numbers that move beyond their stated errors — proposed §S4.17 amendment (Matt's call)

S4.17 (iii) reads *"The drain's first bite at bright rungs is LARGE and real (−0.30..−0.67 dex
truth-anchored)"*. Re-referenced to each venue's own zero-fed baseline that range becomes
**−0.81 to −0.96 dex for the *none* arm**, and **is not a drain measurement at all for the
e07 arm**. Proposed replacement sentence, for Matt to accept or reject:

> The drain's first bite at bright rungs, measured as **excess over the counterpart-null
> arm's own feed-drain**, is **−0.83 to −0.93 dex (band 0.03)** in the circular arm at
> r13p25, from a **single** fed member. In the eccentric arm the quantity is not
> interpretable as a drain: `a_eff_drawn` is cut pre-comb, and the arm's zero-fed offset
> (−0.50 to −0.70 dex) exceeds any loop effect in it.

**This retro-changes nothing on its own.** The banked `a_bg` / `a_eff_drawn` columns are
untouched; the amendment is a change of *reference*, and it needs Matt's sign-off precisely
because it is a re-reading of a spoken headline.

### Confound, declared

The arms do not feed the same count (feed-conditioned table, §D of the script output): at
`n_fed = 0` signal cells mean −0.489 vs null −0.234. Wherever the matched-`n_fed` cell count
is 0 the excess mixes a feed-count difference with a drain difference, and the paired table
above marks the `n_fed s/n` for exactly that reason.

---

## B2 — ITERATION MULTIPLICITY: THE PERSISTENCE RULE (criterion **v2.3**)

`hpc_harbor/ledger/ledger_b2_persistence.py` → `reports/ledger_b2_persistence.npz`.

### The defect

The scoreboard runs at the end of **every** iteration (`glacier_loop.py` step (f)) and a cell
counts as certified if **any** scored iteration certifies. With `n_iter = 6` that is six
looks at a moving template. The trials term the criterion charges (`lnK` + `TRIALS_NAT`
0.578) prices the **fringe and census** multiplicity and prices **no part of the iteration
axis**.

### The companion term — PRE-REGISTERED FORWARD, 2026-07-29

```
criterion v2.3  =  criterion v2.2  AND  PERSISTENCE
PERSISTENCE(p,i) := cert_v22(p,i) AND cert_v22(p,i−1)
```

Declared before scoring: `i = 0` can never certify (no predecessor) — deliberate, since i0 is
the truth-anchored feed, the most optimistic state the loop occupies. v2.3 ⊂ v2.2 always, so
it can only **remove** certifications and can never create a wrong one. It is a
**conjunction on the existing bar** — no floor is re-cut, no estimator switched.

**Why not an additive `ln(n_iter) = 1.79` nat charge?** (a) the six looks are not independent
— consecutive iterations share the data and most of the template, so it over-charges by an
unknown factor; (b) an additive bar change would retro-alter the meaning of every banked
floor, which is forbidden.

### GATE — **PASS**

Re-derived v2.2 mask reproduces the banked `n_cert` column on **706/706** banks exactly.
Nothing downstream would have been quotable otherwise.

### The re-cut

| | events | wrong | on_true |
|---|---|---|---|
| **v2.2** (as banked) | 18 | 16 | 2 |
| **v2.3** (persistence) | **2** | **2** | **0** |
| killed | 16 | 14 | **2** |

Flicker anatomy (iterations certified → consecutive-run lengths):

```
r13p25/e07/s0 psr62   [0,1,4,5]  runs [2,2]   SURVIVES   on_true=False
r13p25/e07/s0 psr36   [1,4]      runs [1,1]   dies
r13p5 /e07/s0 psr62   [0,2]      runs [1,1]   dies
r13p5 /none/s3 psr8   [0,4]      runs [1,1]   dies       on_true=TRUE
… 8 further single-iteration flickers, all wrong, all die
```

### THE FINDING — **the pre-registered expectation is REFUTED by the record**

The brief's expectation was *"nothing true dies; flickers die"*. Measured: persistence kills
**14 of 16 wrong certs — and both true ones.** The campaign's only true certifications
(§S4.24 Finding 1, `r13p5/none/s3` psr 8 at **i0 and i4**) are themselves flickers. The
result is robust to the rule's shape: no 2-consecutive rule and no 2-of-3-window rule
contains `[0,4]`.

And the two survivors are **wrong** — `psr62` in the r13p25 cascade persists with runs
`[2,2]`. So:

> **Persistence is not a manufacturing discriminant.** It is a *stability* test, and the
> manufactured cascade certs are stable while the true ones are not. The thing that
> separates true from manufactured in this campaign remains the **D2 rigidity gate**
> (20/20 kills, R2 the killer every time) — v2.3 does not substitute for it and must not be
> presented as if it did.

**Honest reading, and it cuts against the campaign:** the true certs sit at dlnL 2.35 / 1.60
against a 1.27 bar. A marginal cert that does not survive one iteration of the joint template
moving is telling us something about the **loop's stability**, not about the pulsar. That is
the most valuable output of B2 and it is a *negative* one.

**Status: v2.3 is pre-registered FORWARD only.** It does **not** retro-apply. §S4.24's
verdicts stand as banked under v2.2 with the §S4.24.1/.2 conditions unchanged. Adopting v2.3
for future cells is a pre-registered bar change = charter #1 = **Matt's**.

---

## B3 — FEED-ORDER DEPENDENCE

### Finding 0: it is not where the brief points

Within an iteration the **feed decisions are order-invariant**, and that is worth stating
because it is the first place one looks: `carried` is fixed *before* the candidate loop
(`glacier_loop.py:416`), so `feed_dlnl[k]` is a function of `k` alone (`:434-438`), and
`fed_idx` is re-sorted afterwards (`:465`). **Permuting the feed order cannot change who is
fed.**

The order dependence is real one stage later, in `mstep_quadratic`
(`glacier_loop.py:275,279-291`): **in-place sequential coordinate ascent**, on the
programme's most correlated pair `(log10_fgw, log10_mc)`, ordered by *census rank* — an
arbitrary labelling convention. Two amplifiers: `n_sweep = 2` (truncated, so it is not at a
fixed point at all), and the non-concave branch (`:291`) **silently leaves a slot at its entry
value**, so order changes *which coordinates move at all*.

### Finding 1: **B3 as specified cannot be run at an unquarantined rung**

Measured over the whole bank (n_fed per iteration, every stage-2 venue):

```
r13p9  (all wscale, T30/T40, e07)   n_fed = 0 everywhere
r13p9  (all wscale, T30/T40, none)  n_fed = 0-1
r13p5  / none                       n_fed = 1-2
r13p25 / none                       n_fed = 1-3
r13p5  / e07                        n_fed = 3-19     <- quarantined arm
r13p25 / e07                        n_fed = 3-23     <- quarantined arm
```

**A one-member feed has no feed order to permute.** Every rung outside the eccentric cascade
feeds 0–1 members. So the probe necessarily runs in the e07 arm; `r13p5/e07` skies {0,1,2,3}
(n_fed 16/10/3/3) is chosen over r13p25 because **skies 2 and 3 there carry no certification
claim at all**. Cells stay QUARANTINED from closure claims; this is a diagnostic of the
**M-step operator**, not a claim about those cells.

### Finding 2: the CPU leg **refutes the obvious hypothesis**

`ledger_b3_feedorder.py surrogate` — the operator on a surface whose answer is known
(quadratic, curvature calibrated to S4.15's ~0.2 nat, `SCALE_MC = 0.5`, `step0 = 0.3`,
`n_sweep = 2` as wired), 16 cells of (n_fed ∈ {2,3,5,10}) × (ρ ∈ {0, 0.5, 0.9, 0.99}):

**Max 2-sweep spread over feed orders: 0.0022 dex. 0/16 of the grid exceeds the 0.02 dex
materiality bar — even at ρ = 0.99, n_fed = 10.**

So truncated coordinate ascent is **not**, by itself, order-sensitive at this curvature. This
sharpens B3 into a pre-registration rather than a fishing trip:

> If the venue leg finds order dependence above 0.02 dex, **it cannot be blamed on the sweep
> count.** It would have to come from the two non-quadratic features the operator actually
> has on the comb surface: (i) the non-concave branch silently skipping slots, and (ii)
> genuine multimodality of the harmonic-comb surface — which is exactly the object §C1 below
> asks FORGE-B to detect. **B3 and C1 would then be one finding seen from two sides.**
> A NULL venue result is equally informative: the M-step's 0.01–0.44 dex wander would be a
> property of the *surface*, not the ordering, and no canonicalisation would have prevented it.

**Pre-registered bar (stated before the venue run):** feed order is IMMATERIAL if the
max-over-permutations **converged** spread is below **0.02 dex** (the drain's own 1σ) on every
fed axis in every sky. Above that, the M-step ordering must be made canonical before any
per-member drain number is quoted.

### GPU leg: **POSTED AND PARKED** (charter — a stopped arm never idles another lane)

`hpc_harbor/ledger/ldg_b3.sbatch` is complete and submit-ready; **it was not submitted.**
`sbatch --test-only` returns **`QOSGrpGRES`** on *both* live entitlements:

| lane | state |
|---|---|
| `dgx_iacc` A100-80 (dgx03) | 8/8 allocated; PHOENIX's pinned frozen arm |
| `dgx_iacc` A100-40 (dgx01) | BASELINE's `%4` share |
| `taylor_group_account_batch_gpu` H200 | behind the `nodeupgrade` drain to 2026-11-04 |

There is no GPU headroom for a third session without displacing a claim that is not
LEDGER's to take. Cost when it can run: **< 1 GPU-hr** for all four skies (replay from banked
checkpoints — one venue build per sky plus ~4·n_fed·3 likelihood evals per order). LEDGER's
claim, for the moment it can be honoured, is in the sbatch header: `%2`, `ldg*` prefix,
`ledger_b3_*` banks only, yields immediately, and **host-agnostic by design** — B3 is a
within-run permutation contrast, so unlike the three-way belief arm it does **not** need the
dgx03 pin.

---

## C1 — COMB MULTIMODALITY: SPEC FOR FORGE-B *(for FORGE-B / Matt)*

*Posted here rather than in `GLACIER_capstone.md`: that file is frozen pending Matt's review
of the unstaged GLACIER-LITE addendum (PHOENIX order #2). Fold in with the rest when the
review lands. **No cronus-owned file was edited** — `CW_transition/` is untouched by this
audit, per the harbor hard rule.*

### The problem

FORGE-B's belief is a **Laplace/Fisher Gaussian**: `Cov = inv(Fs + Pi)`, `sig_f`, `sig_mc`
(`forge_b1.py:311-316`). A Gaussian width is only meaningful if the surface has one mode. On
the eccentric comb it very plausibly does not: the WEAVE-3.3 tie makes 32 harmonics of one
source, and `tie_log10_fgw_n(f_orb, n)` places templates at `n·f_orb`. A recovery template
sitting one harmonic off is a **distinct near-degenerate mode**, not a tail of the central
one. The campaign has already seen the consequence in every other layer: S4.17's cascade
("~20 fed comb-harmonic slots, M-step wander among them, coherent joint template
manufacturing maximally-confident fringes at WRONG distances") and D1's verdict that "the
wander built every false fringe". **What has never been checked is whether the belief
FORGE-B hands downstream is unimodal at all.** If it is not, `sig_mc` is not a width — it is
the curvature of one mode among several, and A1's sigma-point average, fed that number,
would smear over the *wrong* set.

### The check (cheap; runs inside the existing source-fit step)

At each source fit, after the damped Newton step, evaluate `lnL_marg` on a 1-D comb probe in
`log10_fgw` at the harmonic offsets the tie itself defines:

```
f_probe(δ) = log10( 10**log10_fgw_hat * (n+δ)/n )    for  δ ∈ {±1, ±2}, n the fitted harmonic
```
plus the same in `log10_mc` along the degeneracy direction `d(log10 mc)/d(log10 f) = −2/5·…`
(the amplitude relation's own slope, so the probe follows the ridge instead of crossing it).

**Report three numbers per fit:**
1. `Δ_multi = max_δ≠0 lnL(δ) − lnL(0)` — how close the best rival mode is, in nat;
2. `w_rival = softmax weight` of the best rival under the fitted prior;
3. `n_modes` = count of probe points with `Δ > −1` nat.

**PRE-REGISTERED THRESHOLD (state it before running):** the belief is **UNIMODAL-CERTIFIED**
if `Δ_multi < −3` nat (rival ≥ 3 nat down) at every probe point. Otherwise it is
**MULTIMODAL** and the Gaussian is refused.

### The fallback, when the check fails

Two options, in cost order:

**(F1) GRID belief.** Replace the Gaussian on `(log10_fgw, log10_mc)` with a normalised
weight over the probe grid, widened to `±3` harmonics and `±3σ_mc`. Downstream this is a
drop-in for A1's sigma-point rule: `estep_sigma` takes `(points, weights)` and does not care
that the weights came from a grid rather than a UT — **the same `logsumexp` reduction, the
same width→0 collapse, the same G-L1 identity.** Cost: `n_grid` E-step evaluations instead of
`2m+1`, so cap `n_grid` the same way `n_belief` is capped.

**(F2) MIXTURE belief.** Fit a Laplace Gaussian at each surviving mode and carry
`Σ_j π_j N(θ_j, Σ_j)`. Sigma points are then the union of each component's, with weights
`π_j · w_i`. Strictly better than F1 where modes are well-separated and sharp; strictly worse
where they merge. Recommended only once F1 has *measured* how many modes there actually are.

**What must NOT happen:** carrying a Gaussian width from a multimodal surface into the
E-step and calling the resulting `q_max` a marginal confidence. That would combine A1's
defect with a wrong belief and produce a number more over-confident than the incumbent —
i.e. the upgrade would make things worse. **A1's sigma-point E-step should therefore be
gated on C1's unimodality check passing, or on F1/F2 being in place.** That gating is the
one hard dependency this audit creates between two items.

---

## C2 — SBC: REPORT-ONLY RESTATEMENT OF SCOPE

C2 remains the **queued production item**; nothing was run for it here. Its scope, restated
now that A1 exists, is strictly larger than when it was queued:

Simulation-based calibration was queued to check the **fringe posterior's** coverage —
"does `q_max` mean what it says". A1 shows `q_max` is a *conditional* quantity, so SBC on the
incumbent pipeline would calibrate `P(fringe | θ_s = θ̂)`, which is not the quantity the
criterion uses it as. **SBC must therefore now cover the belief widths too:** draw
`θ_s ~ prior`, draw the belief, run the E-step (point *and* sigma-point), and check rank
uniformity of the true fringe under **both**. Two ranks, not one. The point-rule ranks are
predicted (by A1's G-L3b mechanism) to be **over-concentrated** — the true fringe landing in
the extreme ranks more often than uniform — with the effect growing in belief width; the
sigma-point ranks are predicted uniform if the belief is unimodal. **If the sigma-point
ranks are *also* non-uniform, that is C1's multimodality showing up in the calibration, and
the two items diagnose each other.** Cost scales as C1's/A1's cap `n_belief`, so C2 should be
queued *behind* A1's cap being fixed, not before.

---

## C3 — Q-TAIL PAYOFF PROPAGATION

`hpc_harbor/ledger/ledger_c3_qtail.py` → `reports/ledger_c3_qtail.npz`,
`reports/LEDGER_C3_qtail.png`.

### The paragraph, for the payoff section

> **SIREN §4's `σ(D_L)/D_L` is the width of the all-right component of a mixture, quoted as
> if it were the width of the posterior.** The criterion that supplies the seeds admits them
> at `q_max > 0.9` (`glacier_loop.py:104`), so each "certified" seed carries up to 10 %
> posterior mass on a *wrong* fringe, and the payoff table propagates none of it. Convolving
> that mass in — as the Gaussian mixture the criterion itself implies, with
> `k ~ Binom(N, 1−q)` seeds wrong and two declared bounds on what a wrong seed does
> (**benign**: it is merely uninformative, so the analysis is an `(N−k)`-seed analysis off the
> banked ladder; **adversarial**: it is confidently wrong, displacing `log10 Mc` by the
> one-fringe reabsorption `2π·σ_mc^(1)` and thence `D_L` through the amplitude relation) —
> gives, at the headline cell `(log10 f, log10 Mc) = (−8.0, +9.5)`, `N_seed = 3`, comparing
> **like for like** on the 90 % half-width: **16.1 % → 23.6 % (benign) or 56.1 %
> (adversarial) at q = 0.90**, `16.6 %`/`17.6 %` at `q = 0.99`. The probability that the
> quoted core width is the right answer at all is `q^N` = **0.729** at `q = 0.90, N = 3` and
> **0.590** at `N = 5`. **The dominant effect is not the widening of any component — it is
> that the mixture has a heavy component at all**: one wrong seed out of three drops the
> analysis onto the two-seed rung, which at this cell is 2.4× wider. Two structural
> consequences follow. First, **more seeds is not monotonically better**: at `q = 0.90`,
> `N = 5` has `P(all right) = 0.59` against `N = 3`'s `0.73`, so the fifth seed buys core
> width while spending assignment reliability. Second, **cells that saturate early are
> robust and cells that need every seed are fragile** — at `(−7.7, +9.5)`, which SIREN §4.1
> notes saturates at `N_seed = 1`, the q-tail costs 1.03× at `q = 0.90`, while at
> `(−8.0, +9.5)`, which needs all three, it costs 1.47×. **The honest form of the payoff
> claim is therefore conditional on `q`, not just on `N_seed`**: the 6–12 % class survives
> only at `q ≳ 0.99`, and at the criterion's own `q > 0.9` admission bar the same
> configuration delivers **17–24 %** on the benign bound — still inside the dark-siren-useful
> band, but no longer the headline number.

### Table (90 % half-width in `ln D_L`; ×N = inflation over the core's own 90 % half-width)

| cell (log10 f, log10 Mc) | N | q | P(all right) | core hw90 | benign | adversarial |
|---|---|---|---|---|---|---|
| (−8.0, +9.5) | 3 | 0.90 | 0.729 | 0.1607 | 0.2364 (**1.47×**) | 0.5607 (3.49×) |
| (−8.0, +9.5) | 3 | 0.95 | 0.857 | 0.1607 | 0.1941 (1.21×) | 0.4613 (2.87×) |
| (−8.0, +9.5) | 3 | 0.99 | 0.970 | 0.1607 | 0.1666 (1.04×) | 0.1755 (1.09×) |
| (−8.0, +9.5) | 5 | 0.90 | **0.590** | 0.1024 | 0.1176 (1.15×) | 0.4036 (3.94×) |
| (−7.7, +9.0) | 3 | 0.90 | 0.729 | 0.2029 | 0.2729 (1.34×) | 0.5567 (2.74×) |
| (−7.7, +9.5) | 3 | 0.90 | 0.729 | 0.1672 | 0.1728 (**1.03×**) | 0.1964 (1.17×) |

### Declared limitations (so the number is not over-read)

`q` is treated as **independent per seed** — correlated fringe errors (a biased template
mis-registering several pulsars the same way, which is precisely the S4.17/D1 wander
mechanism) would be **worse**, not better. Dropping a wrong seed is modelled as dropping the
*last* of the nested freq-ranked ladder, which understates the loss when the wrong seed is
the informative one. `q_max` is a criterion-side confidence, not calibrated frequentist
coverage — and per ANCHOR's D2 finding the Gumbel floor is biased **permissive** at high
zero-fraction, so the true `(1−q)` is plausibly **larger** than nominal. **Every one of these
biases runs the same way: the honest number is worse than the one above.** Asimov/Fisher
cores throughout, inherited from SIREN.

**Figure:** `reports/LEDGER_C3_qtail.png` — left, the implied mixture at `N = 3` for
`q ∈ {0.90, 0.95, 0.99}` against SIREN's core (log density, so the heavy component is
visible); right, 90 % half-width vs `q` for both bounds and both `N`, with the three `q`
values marked. Single linear axis per panel, CVD-safe fixed hues, `N` also carried by line
style.

---

## HOW THE ITEMS BIND TO EACH OTHER

Three of these are not independent, and shipping them in the wrong order would be worse than
shipping none:

1. **C1 gates A1.** A sigma-point average fed a Gaussian width from a multimodal surface is
   *more* over-confident than the incumbent point rule. A1 must not fan before C1's
   unimodality check passes or F1/F2 is in place.
2. **B3 tests C1.** If B3's venue leg finds material order dependence, the CPU leg has already
   excluded the sweep count, leaving comb multimodality as the live explanation.
3. **B2 does not replace the D2 rigidity gate**, and the re-cut says so with data: the certs
   that survive persistence are the manufactured ones.

---

## ADD-LIST (staged, **not committed** — charter #4)

New, all LEDGER's own:

```
LEDGER_stats_audit.md                              44 K   this report
hpc_harbor/ledger/ledger_a1_sigma_estep.py         17 K   A1 machinery + CPU gates (PASS)
hpc_harbor/ledger/ledger_a2_rn_monitor.py          17 K   A2 monitor + gates (PASS) + replay
hpc_harbor/ledger/ledger_b1_drain.py               15 K   B1 re-reference (CPU, banked)
hpc_harbor/ledger/ledger_b2_persistence.py         10 K   B2 v2.3 re-cut (CPU, banked)
hpc_harbor/ledger/ledger_b3_feedorder.py           19 K   B3 replay driver + CPU surrogate
hpc_harbor/ledger/ledger_c3_qtail.py               14 K   C3 mixture calc + figure
hpc_harbor/ledger/ldg_b3.sbatch                   2.7 K   B3 GPU leg — PARKED, submit-ready
hpc_harbor/ledger/ldg_a2_replay.sbatch            2.4 K   A2 venue replay, CPU lane (submitted)
reports/ledger_b1_drain.npz                        92 K   B1 bank
reports/ledger_b2_persistence.npz                 6.6 K   B2 bank
reports/ledger_c3_qtail.npz                       5.1 K   C3 bank
reports/LEDGER_C3_qtail.png                       165 K   C3 figure
```

**Total: ~410 K.** No bank over 100 K; `*_results/` is gitignored, so the three `reports/`
npz files and the png are the only data artifacts and they are deliberately small and lean
(raw columns, not verdicts — per the lean-npz discipline).

**NOT LEDGER's, flagged so authorship is not misattributed at commit time:**
`hpc_harbor/ledger/sieve_t7_eprocess.py` (18 K) — **agent SIEVE-A's**, written into this
directory concurrently. LEDGER never staged it; **SIEVE-A's own session has since staged it**
(it shows as `A` in `git status` alongside LEDGER's files, index touched 11:34). It is in the
tree legitimately and should ride along — but it is **their** work, not part of this audit,
and the commit message should say so.

**Suggested `git add`:**

```bash
git add LEDGER_stats_audit.md \
        hpc_harbor/ledger/ledger_a1_sigma_estep.py \
        hpc_harbor/ledger/ledger_a2_rn_monitor.py \
        hpc_harbor/ledger/ledger_b1_drain.py \
        hpc_harbor/ledger/ledger_b2_persistence.py \
        hpc_harbor/ledger/ledger_b3_feedorder.py \
        hpc_harbor/ledger/ledger_c3_qtail.py \
        hpc_harbor/ledger/ldg_b3.sbatch \
        hpc_harbor/ledger/ldg_a2_replay.sbatch \
        reports/ledger_b1_drain.npz \
        reports/ledger_b2_persistence.npz \
        reports/ledger_c3_qtail.npz \
        reports/LEDGER_C3_qtail.png
```

---

## STOP — FOR MATT

Everything above is delivered and gated. Five things are **not** mine to take:

1. **§S4.17's drain sentence.** B1's re-reference moves the *none*-arm bites beyond their
   stated errors (−0.30..−0.67 → **−0.81..−0.96 dex**) and empties the e07-arm ones. The
   proposed replacement sentence is in §B1; **adopting it is a change to a spoken campaign
   headline** and needs your sign-off. Nothing is changed until you say so.
2. **Adopting criterion v2.3.** Pre-registered forward here, re-cut here, **not adopted** —
   a bar change is charter #1. Note before deciding: it kills 14/16 wrong certs **and both
   true ones**, and the two survivors are wrong. That is not the trade the brief expected.
3. **Applying the A1 and A2 wiring diffs.** Both are written and gated; neither is applied,
   because `glacier_loop.py` is the HELD driver and **PHOENIX's frozen arm is live against
   it right now** (job 12833771). These want to go in **between** fans, not under one.
4. **B3's GPU leg** — `ldg_b3.sbatch`, POSTED AND PARKED on `QOSGrpGRES`. Submit unchanged
   whenever either entitlement frees; < 1 GPU-hr. Not queued behind PHOENIX or BASELINE,
   per the charter.
5. **Git.** `glacier_lite` has no upstream and has diverged 49/10 from `origin/master`.
   Nothing was merged or forced. Staging only; **you commit.**

**In flight:** CPU job `12843673` (`ldga2`) — the A2 monitor replay, resubmitted after a
bank-write fix. The science already landed on the first run (§A2: monitor quiet, 0/116
flagged at every iteration, occupancy at prior); the resubmission only writes the bank.

**One thing I'd put on your list that wasn't in the brief:** the A2 replay proves a
**T=30/n=256 venue build runs on the free `batch` CPU partition**. That is a lane the
programme has been treating as GPU-only, and several read-only post-hoc analyses
(monitors, replays, anything that does not call `estep_per_target`) could run there without
touching either GPU entitlement.
