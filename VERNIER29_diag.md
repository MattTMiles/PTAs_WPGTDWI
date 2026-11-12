# VERNIER-29 DIAGNOSTIC — the K = 16 certification null leak

**DIAGNOSIS ONLY. No new science claims, no ranking changes.** Every cell is RNG-matched
to the banked VERNIER-29 deliverable-6 cells; the diagnostic kernel is verified to
reproduce `vernier29.cert_rows` bit-for-bit before any number is read off it.

Code `hpc_harbor/vernier/vernier29_diag.py`, `hpc_harbor/vernier/vernier29_diag2.py` ·
logs `hpc_harbor/vernier/logs/vernier29_diag{,2}.log` ·
data `reports/vernier29_diag{,2}.npz`. **Staged, not committed.**

---

## Lead — all four hypotheses are dead, and so is the premise

**H1, H2, H3 and H4 are all falsified.** So is the framing that produced them: the leak
is not a V29 configuration effect, the real priors are not the cause, and V25 was never
clean at K = 16 in the first place.

**The surviving mechanism is a shared defect in `vernier24.certify`, inherited unchanged
by V25, V27, V28 and V29: the certification bar is `max(ln K_i, floor_i + err_i)`, and
`ln K_i` exceeds the null-calibrated floor for 95–100 % of pulsars, in every arm, at
every comb count, by 5–8 nat. The `recut_floor` α = 0.05 calibration is dead code. The
operative gate is a bare Bonferroni-style trials term whose realised false-alarm rate is
calibrated to nothing** — empirically 0.02–0.46 % per pulsar-certification, varying with
array geometry and *non-monotonically* with comb count.

**The V29 ranking is untouched.** It is built on Δf_marg, the soft channel, which uses no
threshold, no floor and no `ln K`. Products (c) localisation and (d) chirp mass are
likewise threshold-free. Only the (a) certification column touches the gate, and that
column needs the correction in §7.

---

## 0 — RNG identity

The diagnostic kernel reproduces `vernier29.cert_rows` exactly (`stat`, `bar`, `Kfr`,
`Rcoh`, `cert`, `cert_true`, `pinned` all bit-identical) at K = 8 and 16, scrambled and
not. The banked D6 row reproduces:

| K | f_cert | f_cert_true | f_cert SCRAM | n_cert | n_SCRAM |
|---|---|---|---|---|---|
| 8 | 0.02101 | 0.02101 | 0.00000 | 0.75 | 0.00 |
| 16 | 0.03574 | 0.03318 | 0.01574 | 1.12 | 0.50 |

V28's C20-R1 cell also replays exactly: 6 scrambled certifications over 12 realisations
against the banked 0.50 per run.

---

## 1 — H3 is dead (the brief's metric), and so is my own variant of it

**Metric A, the brief's:** the scramble displaces each comb by U(−500, 500) pc. Measured
across the array that is a **median 299 fringes**, range **196–445**. **Zero of 116
pulsars** are displaced by fewer than 3 fringes. The control is valid by the brief's own
criterion, everywhere, and H3 as written cannot fire.

**Metric B, which the brief did not name and which I expected to be the answer:** the
scramble is ±500 pc but the searched window is ±min(3σ_a, 0.9 L_a), and for **76 of 116**
real-prior pulsars that window is *wider* than 500 pc — so 83 % of displaced combs land
back inside the grid. That looks decisive and it is wrong. Two independent falsifications:

1. **The peak is never on a displaced comb.** For every leaker, the winning grid point
   sits more than half a fringe from all K displaced comb centres (`peak@Ls` = 0.000).
   The null is not certifying a real comb it merely moved.
2. **The arm comparison runs the wrong way.** C1 (real array, *uniform* 10 pc priors,
   in-window fraction 0.063) leaks **more** than C0 (real array, real priors, in-window
   fraction 0.827). In-window probability anti-correlates with the leak.

C20 seals it: its in-window probability is 0.304 and its FPR is 0.043 %, the *lowest* of
any high-comb cell measured.

---

## 2 — H1 is dead

At K = 16 exactly **three** pulsars fire in the scrambled null across 8 realisations:

| pulsar | category | null fire | signal fire | K_i | σ_d (pc) | window (pc) | disp (fringes) |
|---|---|---|---|---|---|---|---|
| J0437−4715 | VLBI | 0.250 | 0.875 | **2.5** | 0.11 | 0.33 | 378 |
| J1630+3734 | DM default | 0.125 | 0.000 | **531** | 24.00 | 72.0 | 263 |
| J2124−3358 | timing PX | 0.125 | 0.125 | **1225** | 57.23 | 171.7 | 334 |

Two of the three are high-K_i, not sub-fringe. Only **one** pulsar in the whole array is
at K_i ≤ 3. Median K_i of leakers 531 vs 5434 for non-leakers; Spearman(fire, K_i) =
**−0.138**. H1 predicted leakers *exclusively* at low K_i; they are not.

**C20 kills H1 outright and without a run.** C20's priors are uniform 10 pc for every
pulsar and its array is V24 synthetic — no J0437, no anchor, no heterogeneous prior. Its
K_i range is **146–192** (pure geometry, since the prior is constant) and **zero** pulsars
sit at K_i ≤ 3. H1 has nothing to bite on, and C20 leaks anyway.

**Consequence for the brief's item 6:** the sub-fringe rule (exclude or floor K_i ≥ 3) is
**not the repair**. In the V29 cell it would remove one pulsar; in C20 it would remove
none, and C20 leaks.

---

## 3 — H4 is dead

By construction: `cert_rows` builds `n_inc` inside the per-pulsar loop from the comps in
scope and calls `recut_floor` in the same iteration. There is no floor cache and no
cross-K state — `vernier25.floor_for`'s `_FLOOR_CACHE` holds the *detection* F-floor and
is never touched by the certification path. Empirically, in both codebases:

| n_null | V29 K=16 f_SCRAM | V29 median floor | C20 null certs |
|---|---|---|---|
| 32 | 0.02359 | 0.431 | 1 |
| 128 | 0.01498 | 0.421 | 0 |
| 384 | 0.02359 | 0.402 | — |

Flat. (The bounce between 0.015 and 0.024 on the *same* cell is itself a measurement of
how noise-dominated the banked number was.)

---

## 4 — H2's mechanism is dead; the leak is K-dependent but NOT monotonic

`ln K_i` is computed from each pulsar's own window and its own `dL_min`; it is not global
and not a median. And the bar *rises* with K (7.51 → 8.56 from K = 2 to 16, because
`dL_min` is a minimum over more combs), which is the opposite of a bar failing to keep up.

The high-statistics sweep (32 realisations per K, 3712 pulsar-certifications per cell):

| K | signal ev | null ev | null FPR % | 95 % CI % | null/signal | med lnK |
|---|---|---|---|---|---|---|
| 2 | 23 | **17** | **0.458** | [0.267, 0.733] | 0.739 | 7.51 |
| 4 | 14 | 3 | 0.081 | [0.017, 0.236] | 0.214 | 7.97 |
| 8 | 26 | 3 | 0.081 | [0.017, 0.236] | 0.115 | 8.37 |
| 12 | 40 | 4 | 0.108 | [0.029, 0.276] | 0.100 | 8.51 |
| 16 | 42 | 9 | 0.242 | [0.111, 0.460] | 0.214 | 8.56 |

χ² for a constant rate = **20.11** on 4 dof, **p = 0.000** — the rate is genuinely not
constant. But **Spearman(null events, K) = −0.100**: it is not monotonic, and **K = 2 is
the worst cell in the sweep**, at nearly twice the K = 16 rate. Signal events, by
contrast, rise cleanly with K (Spearman **+0.900**). So "the leak grows with K" is false,
and a K-dependent defect in the ordinary sense is not what this is.

---

## 5 — What survives: the gate is not null-calibrated at all

| arm | K | bar == lnK | med lnK | med floor | gap (nat) |
|---|---|---|---|---|---|
| C0 real array + real priors | 2 | 99.7 % | 7.61 | 0.00 | +7.61 |
| C0 | 8 | 98.9 % | 8.35 | 0.23 | +8.11 |
| C0 | 16 | 96.0 % | 8.57 | 0.47 | +8.10 |
| C1 real array + uniform 10 pc | 2 | 100.0 % | 4.57 | 0.00 | +4.57 |
| C1 | 8 | 98.9 % | 5.21 | 0.22 | +4.99 |
| C1 | 16 | 94.7 % | 5.42 | 0.46 | +4.96 |
| C2 synthetic + uniform 10 pc | 2 | 100.0 % | 4.55 | 0.00 | +4.55 |
| C2 | 8 | 99.1 % | 5.25 | 0.08 | +5.18 |
| C2 | 16 | 97.3 % | 5.44 | 0.21 | +5.22 |

**The incoherent-comb floor never binds.** It sits at 0.0–0.5 nat while `ln K_i` sits at
4.5–8.6. So `recut_floor`'s α = 0.05 — the whole CRITERION-v2.2 apparatus, the Gumbel
fit, the zero-fraction gate, the bootstrap error — is bypassed in 95–100 % of pulsars.
Certification is decided by `ln K_i` alone.

That is a Bonferroni-style trials correction, and it is exactly right *only if* the
contrast statistic's tail is exponential with unit scale. Nothing guarantees that. The
realised false-alarm rate is therefore whatever the tail gives, and it varies with the
things that shape the tail:

### The decisive arm test — the leak tracks the ARRAY, not the priors

40 realisations per arm, 4640 pulsar-certifications per cell:

| arm | configuration | K | sig ev | null ev | null FPR % | 95 % CI % |
|---|---|---|---|---|---|---|
| C0 | real array + REAL priors | 8 | 34 | 3 | 0.065 | [0.013, 0.189] |
| C0 | real array + REAL priors | 16 | 51 | 9 | 0.194 | [0.089, 0.368] |
| C1 | real array + uniform 10 pc | 8 | 4 | 2 | 0.043 | [0.005, 0.156] |
| C1 | real array + uniform 10 pc | 16 | 184 | **14** | **0.302** | [0.165, 0.506] |
| C2 | synthetic array + uniform 10 pc | 8 | 3 | 1 | 0.022 | [0.001, 0.120] |
| C2 | synthetic array + uniform 10 pc | 16 | 110 | 3 | 0.065 | [0.013, 0.189] |

Exact conditional binomial, K = 16: **C2 vs C1 p = 0.013 (significant)**; C1 vs C0
p = 0.405; C2 vs C0 p = 0.146.

**C1 has uniform 10 pc priors and leaks at the highest rate of any arm.** Swapping the
real priors *out* does not close the leak — it widens it. What separates the leaking arms
from the clean one is the **array**: real sky positions and real rms weights versus V24's
isotropic synthetic array. The brief's premise — that V29's real heterogeneous priors and
the J0437 anchor opened this — is falsified by its own control.

---

## 6 — The rider: C20, and were V25's nulls clean?

### C20 (V28, γ = 0.87 + 10 seeds, K_available ≈ 19, ladder R1), replayed on V28's seeds

12 realisations, **121 certification steps, 5.9 combs per step** (COMB_CAP = 16 never
binds — the ladder only banks ~10), exposure 14036 pulsar-certifications:

- signal certifications **36** (all 36 true-fringe), scrambled-null certifications **6**
- null FPR **0.043 %**, 95 % CI [0.016, 0.093] %
- null/signal **0.167**

Steps 1–3 on C20: leakers sit at K_i ≈ 173 against an array range of 146–192, i.e.
indistinguishable from the bulk; **zero** pulsars at K_i ≤ 3; scramble displacement
**714 fringes** (range 607–798), so the control is valid by metric A; and the floor recut
at 4× nulls does not close it. **H1 and H3 cannot explain C20, and neither can H2 or H4.**
C20's rate is consistent with arm C2 — same synthetic array, same uniform priors — which
is exactly what the surviving mechanism predicts.

**Label correction worth carrying:** C20 is filed under "K_av = 18", but the quantity that
enters a certification is the *comb count*, and that is **5.9**, not 18. K_available and
combs-in-certification are different numbers and the cell labels conflate them.

### Were V25's K = 16 nulls clean, or starved? — **Both, and it never ran a 16-comb null**

V25 Arm K at K = 16, as banked: CERT 0.1/run, true 0.1, wrong 0.0, null 0.0, over 10
realisations.

- **Starved on signal**: 0.1 per run is **one certification in the entire cell**. V25's
  "every one of them on the correct fringe" is a statement about a single event.
- **The null arm's exposure is not what it looks like.** It ran on every certifying turn,
  giving ~9164 pulsar-certifications — but V25 is a *chain*: turn *t* certifies against
  *t* accumulated combs, capped at **8**. The exposure is spread over comb counts 1–8 and
  **none of it is at 16 combs. V25 never ran a 16-comb null.**
- **Is V25's zero even in tension with this run?** Against the matched arm C2:
  expected 2.0 events at C2's 8-comb rate → **P(observe 0) = 0.139, consistent**;
  expected 5.9 at C2's 16-comb rate → P = 0.003.

**So "V25 was clean at K = 16 and V29 opened it" does not survive.** V25 has no 16-comb
null to be clean, and its low-comb zero agrees with the same machinery measured here.

---

## 7 — What stands, what is corrected, what is exposed

**The V29 ranking is NOT touched.** Δf_marg (the soft channel, VERNIER-23 Part 3) applies
no threshold, consults no floor and never evaluates `ln K`. Products (c) sky-area and (d)
chirp mass are Fisher quantities with no gate. The ranked table, the driver decomposition,
the diminishing-returns curves and the J0437 accounting all stand as published.

**One VERNIER-29 statement is corrected.** V29 D6 reported the certification channel as
"clean at K = 8, scrambled null = 0.000". At 8 realisations that was 0 events. At **40**
realisations C0 K = 8 gives **3 null events, FPR 0.065 %**, i.e. ~0.075 false
certifications per realisation against 0.85 real ones — **≈ 9 % contamination, not zero**.
The K = 8 (a) column should be quoted with that contamination, not as clean. The direction
of V29's conclusion — use (a) at K = 8, not at K = 16 — is unchanged; its absolute
cleanliness claim was an artefact of 8 realisations.

**The quarantine list, itemised.**

| item | exposed? | why |
|---|---|---|
| V29 ranking (soft channel) | **No** | threshold-free |
| V29 (a) column, K = 8 | **Partly** | ~9 % contamination, was reported as 0 % |
| V29 (a) column, K = 16 | **Yes** | already withdrawn in V29 D6; FPR 0.194 % |
| **combs\* = 7.8** | **Yes, most of all** | an *onset* threshold |
| **K\*₅₀ ≈ 14–16** | **Yes** | an *onset* threshold |
| V27 K ≥ 16 columns | **Measured, not predicted** | V27 banks its own null column |

**Onset thresholds are the most exposed quantities in the series, and for a structural
reason: contamination is worst where the signal rate is lowest, which is precisely at the
onset.** At C2-like rates (0.022–0.065 %) the machinery produces 0.026–0.075 false
certifications per step, against onset-level signal rates of ~0.1 per step. **combs\* and
K\*₅₀ should therefore be read as *lower* bounds on the true threshold** — chance
certifications make the channel look like it switches on earlier than it does.

**V27's K ≥ 16 columns do not need a quarantine, because V27 already measured this.**

| cell | CER/step | TRUE | WRONG | NULL | NULL/CER | **WRONG/CER** |
|---|---|---|---|---|---|---|
| P12 | 0.075 | 0.081 | 0.006 | 0.000 | 0.000 | 0.083 |
| P16 | 0.181 | 0.113 | 0.087 | 0.006 | 0.034 | **0.483** |
| P20 | 0.181 | 0.181 | 0.019 | 0.000 | 0.000 | 0.103 |
| P24 | 0.094 | 0.106 | 0.019 | 0.013 | 0.133 | 0.200 |
| P28 | 0.106 | 0.131 | 0.000 | 0.000 | 0.000 | 0.000 |
| P32 | 0.075 | 0.094 | 0.013 | 0.000 | 0.000 | 0.167 |

Predicted contamination under the surviving mechanism is 15–100 %; V27's own NULL column
sits at or below the low end (0–13 %). The action is to **quote CER with NULL beside it**,
not to withdraw the columns.

**And the number that dwarfs all of this is already banked:** at P16, **WRONG/CER = 0.48**
— roughly half of V27's certifications in that cell are on the *wrong fringe*. That is an
order of magnitude larger than any scrambled-null leak, it is independent of this
diagnosis, and it is the thing that should actually govern how the certification counts
are read.

---

## 8 — The repair, and what it is not

**The sub-fringe rule (brief item 6) is not the repair** — §2. Neither is the
window-aware scramble I built as a candidate control fix: it moves K = 8 from 0.000 to
0.00862 and K = 16 from 0.01574 to 0.00395, i.e. it moves the two cells in *opposite*
directions, which is what noise does.

**The repair is to make the gate actually calibrated.** The bar should be a floor recut at
the *trials-corrected* level — calibrate `recut_floor` at α/K_i rather than taking
`max(ln K_i, floor at α)`. That replaces an uncalibrated Bonferroni term with a null
distribution measured on the same statistic, the same window and the same comb set, which
is what the incoherent-comb null was built for and is currently prevented from doing.
Under that change `combs*` and `K*₅₀` must be re-derived, and both are expected to move
*up*.

This diagnostic does not implement it — that is a machinery change with consequences for
five instalments' banked columns, and it is a separate decision.

---

## 9 — Distance-file coupling (scope check only)

Of the three V29 leakers, **one** has its K_i move by more than 2× between the tracked
feather `pdist` and `get_distance.py`:

| pulsar | σ_d used | σ_d get_distance | K_i used | K_i get_distance | ratio |
|---|---|---|---|---|---|
| J0437−4715 | 0.11 | 1.30 | 2.5 | 22.6 | **8.90** |
| J1630+3734 | 24.00 | (absent from file) | 531 | — | — |
| J2124−3358 | 57.23 | 57.20 | 1225 | 785 | 1.56 |

**Coupled for J0437 only**, and only because the brief's 0.11 pc and the repo's Deller-2008
1.32 pc differ — the same discrepancy VERNIER-29 D0 already flagged. It moves J0437 across
the K_i ≤ 3 boundary, so any rule keyed on that boundary inherits the dispute. It does not
affect the verdict, because the verdict does not rest on K_i.

---

## One-line verdict

**All four hypotheses are dead — H3 because the scramble displaces by a median 299 fringes
and never fewer than 196, so the control is valid by the brief's own metric everywhere;
H1 because two of the three V29 leakers are at K_i = 531 and 1225 rather than sub-fringe,
and because C20 leaks with uniform 10 pc priors, zero pulsars at K_i ≤ 3 and a K_i range
of 146–192 that is pure geometry; H2 because `ln K_i` is per-pulsar and the bar *rises*
with comb count while the leak does not (χ² p = 0.000 for non-constancy but Spearman
−0.100 in K, with K = 2 the worst cell in the sweep at 0.458 %); and H4 because the floor
is recut inside every call and the rate is flat from 32 to 384 nulls in both codebases —
and the premise dies with them, because the decisive 40-realisation arm test shows the
leak tracks the ARRAY and not the priors (C1, real array with *uniform* 10 pc priors,
leaks hardest at 0.302 % against synthetic-array C2's 0.065 %, p = 0.013, while real
priors add nothing significant at p = 0.405), and because V25 never ran a 16-comb null at
all (COMB_CAP = 8, one signal certification in the whole cell) so its "clean" is a
low-comb bound that this run's matched arm reproduces at P(observe 0) = 0.139; what
survives is a shared defect in `vernier24.certify` inherited unchanged by V25/V27/V28/V29
— the bar is `max(ln K_i, floor + err)` and `ln K_i` wins for 95–100 % of pulsars by 5–8
nat, so CRITERION-v2.2's α = 0.05 null calibration is dead code and certification is
gated by a bare trials term whose realised false-alarm rate is calibrated to nothing and
runs 0.02–0.46 % depending on array geometry and comb count; the V29 RANKING IS UNTOUCHED
because it is built on the threshold-free soft channel, but V29's own "K = 8 is clean"
must be corrected to ~9 % contamination (3 null events in 4640 at 40 realisations, not 0
in 928 at 8), and the quantities genuinely exposed are the ONSET thresholds — combs\* =
7.8 and K\*₅₀ ≈ 14–16 — which are exposed hardest precisely because contamination bites
worst where the signal rate is lowest, and which should now be read as LOWER bounds;
V27's K ≥ 16 columns need their own banked NULL column quoted beside CER rather than
quarantine, and the already-banked WRONG/CER = 0.48 at P16 is an order of magnitude
larger a problem than anything in this diagnosis.**
