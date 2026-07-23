# LOTTERY re-cut — the CHORUS m1 e=0.3 (lit) grade, quoted as a BAND

**Agent LOTTERY, ACCRE, 2026-07-20.** General-queue H200 (`batch_gpu`, `taylor_group_acc` /
`taylor_group_account_batch_gpu`). Reserved `dgx_iacc` share and every SPARK-3 job untouched —
the dgx numbers below are **read from the existing bank**, never re-run. Staged, not committed.

> **Budget (measured).** 3 basis states × 2 cells × (100 nullN + 20 sig) = **720 realisations,
> 3 190 GPU-s ≈ 0.89 GPU-hr** on general-queue H200 (3 concurrent jobs, 3 of 8 grantable h200).
> Campaign running total (tier 1 + tier 2 + this re-cut) ≈ **1.65 GPU-hr** against the ~60 GPU-hr
> cap — **~2.8 %**. No trim needed.
>
> Scripts: `hpc_harbor/lottery/lottery_recut_band.py` (+ `lot_band.sbatch`,
> `lottery_recut_band_figure.py`). Bank: `reports/lottery_recut_band.npz`. Figure:
> `LOTTERY_recut_band.png`. Provenance sidecars: `LOTTERY_results/basis_B{2,3,4}.json`.

---

## 0. THE DECISION, AND WHAT THE ENSEMBLE ACTUALLY IS

LOTTERY tier-2 §2 found that `m1e03_lit` — RECUT's headline *"a single e = 0.3 member does NOT
switch the count on"*, called *"the single most expensive consequence of the floor defect in the
whole repo"* — **flipped `below` → `CONFIRMED`** when its floor was refit from a fresh 100-null
draw on a different host. This re-cut replaces that point grade with a **band over an ensemble of
GWB-basis states**, under a **pre-registered two-sided reading**.

**What was reachable, stated plainly.** The brief asked for ≥3 hosts if the general queue allows,
else host + thread-count variations. **Three hosts are not reachable**, for two independent
reasons, both verified rather than assumed:

- the QOS `taylor_group_account_batch_gpu` grants **only** `gres/gpu:nvidia_h200=8` and
  `nvidia_gh200=2` (DenyOnLimit) — and the *only* h200 node on the cluster is `hgx03`;
- the two `gracehopper` (gh200) nodes are **aarch64**, while the harbor container's `jug-gpu` env
  is **x86-64** (`python3.12` ELF x86-64; `jaxlib` + `nvidia-*-cu12` wheels). It cannot run there.

So the **fallback branch was taken**: one host varied by BLAS thread count, plus the banked dgx
cut supplying the one genuine second host. **Five basis states:**

| basis | host | BLAS thr | `L_gwb` fingerprint | source |
|---|---|---|---|---|
| **B0** | dgx | 8 | *(banked)* | CHORUS `recut_chorus.npz` — **read only** |
| **B1** | hgx03 | 8 | `f92c9e36b460d6f5` | LOTTERY tier-2 |
| **B2** | hgx03 | 4 | `18f7afd5c9607a44` | this re-cut |
| **B3** | hgx03 | 16 | `92c61cc47e801d78` | this re-cut |
| **B4** | hgx03 | 2 | `3d93e3637b5d5300` | this re-cut |

**The four fingerprints are all distinct** — the thread count really did rotate the
eigenbasis, so B1–B4 are four genuinely different bases and not four repeats. Seeds are
**identical** across states (the work plan is imported verbatim from `lottery_tier2._plan()`), so
the basis is the only free variable.

`m2e03_lit` is carried throughout as a **far-from-bar contrast**. Without it, a wide band on the
target would be unreadable — it could mean "this cell is fragile" *or* "the estimator is noisy
everywhere". The contrast separates those.

### 0.1 The pre-registered rule (fixed before the ensemble was scored)

The standing RECUT convention is **one-sided** — `CONFIRMED` iff count(floor+err) > 1, `MARGINAL`
iff count(floor) > 1, else `below` — so it has no way to say *"the refutation is inside the
error"*, and a cell whose floor moves by ~1 σ flips its headline. The band rule reads the floor
error in **both** directions:

> **CONFIRMED** iff `count(floor + err) > 1` — survives the pessimistic floor
> **below** iff `count(floor − err) ≤ 1` — fails even the optimistic floor
> **MARGINAL** otherwise — **the bar lies inside 1 σ; it neither confirms nor refutes**
> **Ensemble grade** = the unanimous label if all states agree, else **MARGINAL (band)**, quoting
> [min, max] of floor and count.

**The external quote is fixed in advance and does not depend on the outcome:** it stays
**e = 0.5 single / e = 0.3 pair** (conservative, band-safe) whichever side this cell's central
value lands.

---

## 1. THE RESULT — the grade is a band, and it straddles

Figure `LOTTERY_recut_band.png`.

| cell | basis | host/thr | floor ± bs | cnt(fl−e) | count | cnt(fl+e) | **BAND** |
|---|---|---|---|---|---|---|---|
| **m1 e=0.30** | B0 | dgx / 8 | 11.30 ± 1.02 | 0.97~ | 0.70 | 0.43 | **below** |
| | B1 | hgx03 / 8 | 9.91 ± 0.86 | 1.45 | 1.20 | 1.05 | **CONFIRMED** |
| | B2 | hgx03 / 4 | 10.21 ± 0.64 | 1.20 | 1.20 | 1.05 | **CONFIRMED** |
| | B3 | hgx03 / 16 | 10.02 ± 0.67 | 1.25 | 1.20 | 1.10 | **CONFIRMED** |
| | B4 | hgx03 / 2 | 9.93 ± 0.84 | 1.45 | 1.20 | 1.00 | **MARGINAL** |
| **m2 e=0.30** | B0 | dgx / 8 | 10.05 ± 1.41 | 3.47~ | 2.77 | 2.07 | CONFIRMED |
| | B1 | hgx03 / 8 | 9.93 ± 0.49 | 3.35 | 2.80 | 2.40 | CONFIRMED |
| | B2 | hgx03 / 4 | 10.27 ± 0.56 | 2.85 | 2.45 | 2.35 | CONFIRMED |
| | B3 | hgx03 / 16 | 9.89 ± 1.11 | 3.65 | 2.75 | 2.15 | CONFIRMED |
| | B4 | hgx03 / 2 | 10.25 ± 0.36 | 2.75 | 2.60 | 2.40 | CONFIRMED |

*(~ = optimistic-side count linearised from the banked `corr`/`corr_lo` pair; B0 is the only row
where the raw signal columns were not re-scored here.)*

| cell | states | floor band | spread / mean err | count band | **ENSEMBLE GRADE** | pointwise labels |
|---|---|---|---|---|---|---|
| **m1 e=0.30 lit** | 5 | **9.91 – 11.30** | **1.74** | **0.70 – 1.20** | **MARGINAL (band)** | CONFIRMED / MARGINAL / below ← **unstable** |
| m2 e=0.30 lit | 5 | 9.89 – 10.27 | 0.48 | 2.45 – 2.80 | **CONFIRMED** | CONFIRMED (all) |

> ### **THE VERDICT — `m1 e = 0.3, lit` is MARGINAL (band): floor 9.91–11.30, count 0.70–1.20.**
> **It neither confirms nor refutes the single-member e = 0.3 switch-on.** Its across-basis floor
> spread is **1.74× its own mean bootstrap error** — the systematic is *larger than the quoted
> error*. The contrast cell `m2 e = 0.3` sits at **0.48×** and is **CONFIRMED in all five states**:
> the estimator is *not* noisy everywhere. **This cell specifically is on the knife-edge.**

---

## 2. THE SHARPEST FINDING — the flip is NOT a two-sided-error effect, it is a HOST effect

Two decompositions, both of which change what the convention has to require.

### 2.1 At fixed basis, the band rule does not rescue the cell

Applying the two-sided rule to the **banked dgx columns alone** leaves `m1e03_lit`
**`below`** — its optimistic-side count is **0.97**, just under the bar. Repo-wide, at fixed
basis (banked columns, no GPU):

| bank | cells | grades that move | CONFIRMED touched |
|---|---|---|---|
| CHORUS `recut_chorus.npz` | 26 | **0** | 0 |
| SURFACE `recut_surface.npz` | 108 | **5** (`below` → `MARGINAL`) | **0** |

So the reversal is **not** produced by widening the error bar. **It is produced by changing the
basis.** The across-basis spread is a systematic that the bootstrap error *does not model at all*.

*(SURFACE already banks `corr_hi` — the true count at floor−err — so the rule costs nothing to
compute there. The 5 movers are all near-bar `below` cells at h = −13.25/−13.00.)*

### 2.2 Thread-count rotation is benign; the HOST is where the systematic lives

| cell | within-host (B1–B4: 4 distinct bases) | cross-host (hgx03 − dgx) | vs mean bs err |
|---|---|---|---|
| **m1 e=0.30** | floor 10.02, **sd 0.136**, range 0.30; counts **1.20 / 1.20 / 1.20 / 1.20** | **Δfloor −1.29 nat (−11.4 %)**, Δcount +0.50 | within **0.18σ** · cross **1.72σ** |
| m2 e=0.30 | floor 10.08, sd 0.203, range 0.38; counts 2.80 / 2.45 / 2.75 / 2.60 | Δfloor **+0.03 nat (+0.3 %)**, Δcount −0.12 | within 0.32σ · cross **0.05σ** |

Two things worth reading twice:

1. **Across four genuinely different eigenbases on one host, the certified count for `m1e03` is
   IDENTICAL — 1.20, four times.** The floor moves by 0.136 nat (0.18 σ). **A pure basis rotation
   is benign.** This is the first direct measurement of RECUT §0.1's rotation-invariance claim at
   **T = 30**; it was previously established only at T = 15.
2. **The host step is 1.72 σ for `m1e03` and 0.05 σ for `m2e03` — same host pair, same code,
   same seeds.** Whatever the host does, it is **not a common-mode bar offset**: it moved one cell
   by −11.4 % and left the other at +0.3 %.

**The honest caveat, stated because it bounds the claim:** there is only **one** dgx estimate per
cell (n = 1, its own error 1.02 / 1.41). The `m1e03` gap is **1.27 σ of that single estimate's own
error**. So *"a cell-dependent host systematic"* and *"one unlucky dgx null draw"* are **not
separable with the data in hand** — distinguishing them needs replicate dgx null blocks, which I
have not run and will not (reserved share). **Either way the operational conclusion is identical:
the point grade is not reproducible, and the band is the honest quote.**

---

## 3. CONVERGENCE WITH SPARK-3 ARM (b) — same hazard class, found twice this week, and one correction

`reports/SPARK3_second_arrow.md` §4.4 measured the **same dgx-vs-hgx03 floor offset** on a
different problem, independently, within the same week:

| | SPARK-3 §4.4 (venue-A floor) | LOTTERY (this re-cut) |
|---|---|---|
| dgx / cronus lane | 122.461 ± 4.895 | m1e03 **11.30 ± 1.02** · m2e03 **10.05 ± 1.41** |
| hgx03 lane | 138.66 (3 blocks: 133.8 / 141.2 / 140.9) | m1e03 **10.02** (4 bases) · m2e03 **10.08** |
| offset | **+16.19 = +13 % = 3.0 σ** | m1e03 **−11.4 % = 1.7 σ** · m2e03 **+0.3 % = 0.05 σ** |
| replicates | 3 independent null blocks on hgx03 | 4 independent bases on hgx03 |

**Both campaigns found the same thing: the null floor carries a host systematic of order 10 %,
comparable to or larger than its own quoted error, and it is a systematic rather than a draw.**
Neither was looking for it.

**LOTTERY independently corroborates SPARK-3's mechanism.** SPARK-3 §4.4 argued the culprit is
`load_or_build_L_gwb`'s `w = np.clip(w, 0.0, None)` on a **near-degenerate** spectrum: a different
host's `eigh` returns slightly different tiny/negative eigenvalues, and the clip turns that into a
slightly different **effective covariance** — so a host change is **not a pure rotation**. This
re-cut supplies the missing control: **rotating the basis while holding the host fixed (four
thread counts, four fingerprints) moves the floor by 0.18 σ and the count not at all.** The
rotation part is demonstrably benign, so the residual host effect must come from the
non-rotation part — exactly the clip, exactly as SPARK-3 argued.

> **ONE CORRECTION THAT TRAVELS BACK TO SPARK-3.** SPARK-3 treats its offset as a **common-mode
> bar offset** ("+13 % bar → a fortiori for an EDGE-POSITIVE"). That reasoning is sound *for
> SPARK-3's own problem*, where the offset was measured directly with three replicate blocks —
> **it is not safe to generalise.** Here the *same host pair* produced **−11.4 % on one cell and
> +0.3 % on its neighbour**, and with the **opposite sign** to SPARK-3's. **The offset is
> problem- and cell-dependent; it must be measured per problem, never inherited.** SPARK-3's
> a-fortiori verdict stands on its own measurement; the *number* should not be reused elsewhere.

### 3.1 The mechanism has an un-closed hole — and it is exactly the one this cell fell into

RECUT §6 records: *"`b1_L_gwb` is canonical and gated — `--cpus-per-task` is no longer an
undeclared input."* **That is true at T = 15 only.** `reports/b1_L_gwb_manifest.npz` reads:

```
shape       [3248 3248]        <- T = 15 baseline
cpus        8
gate        ANCHOR g1: 80 banked ig_nullN_* T=15 realisations replayed THROUGH this bank
g1_passed   True
```

At **T = 30** the extended GWB is **5336²**. There is **no banked square root for it** — the bank
shape-mismatches and `NoiseDrawer` raises, so every T = 30 realisation goes down the
**recompute** path, which prints its own provenance as `RECOMPUTED-UNSAFE`. **All 26 CHORUS cells
are T = 30.** So for CHORUS — including the headline cell — **`--cpus-per-task` and the host are
still undeclared inputs.** That is not a new defect; it is the *known* defect with an uncovered
range, and it is the direct cause of the flip reported here. **Banking `L_gwb` at T = 30 and
gating it the way ANCHOR g1 gated T = 15 would close it** — the same action SPARK-3 §4.4 called
for ("this is a further argument for banking `L_gwb`"), now with a specific shape and a specific
set of cells attached.

---

## 4. THE STANDING CONVENTION PROPOSAL — for the next doc pass

Offered as a proposal, not an adopted change; the grades below are Matt's to cut.

1. **Every floor carries host/basis provenance in its bank metadata** — `host`, affinity CPU
   count, the printed `L_gwb` fingerprint, and `T`. *Demonstrated here*: `lottery_recut_band.py`
   writes `LOTTERY_results/basis_B*.json` per state, and the driver **refuses to run** if the
   observed CPU affinity disagrees with the declared basis label, so a mislabelled basis cannot
   enter a bank. (It caught a real error: `harbor_py` runs apptainer with `--cleanenv`, so
   `SLURM_CPUS_PER_TASK` is *not* visible inside the container — the first submission refused
   rather than banking three states all silently labelled "8 threads".)
2. **Near-bar grades are quoted MARGINAL with the band, never flipped by a re-cut on a different
   host.** "Near-bar" is the pre-registered §0.1 test: `count(floor−err) > 1 ≥ count(floor+err)`.
3. **A grade may only be *changed* by a re-cut that either reports ≥ 2 basis states, or runs on
   the same `L_gwb` fingerprint as the original.** A single-host re-cut may *report* a
   discrepancy; it may not *settle* one. (This re-cut is bound by its own rule: it reports a
   band, and does not overturn RECUT.)
4. **Bank `L_gwb` at T = 30 and gate it as ANCHOR g1 gated T = 15** (§3.1). Until then, T = 30
   floors should be quoted with basis provenance attached.

**The cost of adopting rules 1–3 is bounded and measured, not speculative:**

- **No CONFIRMED grade anywhere can move.** `CONFIRMED` is defined *identically* in both rules,
  and `count(floor−err) ≥ count(floor)`, so `below`(band) ⊆ `below`(pointwise). **The band rule
  can only ever convert `below` → `MARGINAL`.**
- Measured: **0 of 26** CHORUS grades and **5 of 108** SURFACE grades move. **SURFACE's headline
  — onset exists, 59/108 cells — is preserved exactly**; all 59 CONFIRMED cells are untouched.
- The only *headline* affected repo-wide is this one cell.

---

## 5. WHAT THIS CHANGES IN THE DOCS — and what it deliberately does not

**Changes.** RECUT §2's wording *"the e = 0.3 switch-on **does not survive** … **REFUTED**"* and
§2.2's *"a published CHORUS sentence that is now FALSE"* both rest on `m1e03_lit = 0.70, below`.
Under the band that cell is **MARGINAL — undecided**. The right wording is **"not confirmed"**,
not "refuted": **CHORUS's original claim is not restored, and RECUT's correction is not
overturned — the cell simply does not settle either way at the precision available.**

**Does NOT change — the externally quotable number is untouched, as pre-registered:**

> ### **With ONE eccentric member: switch-on at e = 0.5. With TWO OR MORE: at e = 0.3.**

This survives *because* of how the band lands, on both sides:
- `m1 e=0.3 lit` → **MARGINAL**, so a single e = 0.3 member is **still not confirmed** → the
  single-member threshold stays at **e = 0.5** ✓
- `m2 e=0.3 lit` → **CONFIRMED in all five basis states** → the pair threshold at **e = 0.3**
  stands on a cell that is demonstrably basis-robust ✓

The companion `m1 e=0.3 vlbi` (banked 1.03 [0.60], already **MARGINAL** pointwise) **needs no GPU**:
by the monotonicity in §4, a pointwise-MARGINAL cell is MARGINAL under the band rule too. Its
label is unchanged.

RECUT's structural results are untouched: the floor-defect diagnosis, the 26-cell re-cut, the
Gumbel-vs-zero-fraction gate, the clock-not-shared result, and the two-regime channel statement
all stand — none is a near-bar grade.

---

## 6. WHAT IS HANDED BACK

1. **The re-cut grade: `m1 e = 0.3 lit` = MARGINAL (band), floor 9.91–11.30, count 0.70–1.20.**
   Neither confirms nor refutes. External quote unchanged (§5).
2. **The convention proposal (§4)** — provenance in bank metadata; near-bar grades quoted as
   bands; re-cuts require ≥ 2 basis states to settle a grade. Cost measured: 0/26 CHORUS,
   5/108 SURFACE, **0 CONFIRMED grades ever**.
3. **Bank `L_gwb` at T = 30** (§3.1) — the concrete hole behind this flip, and behind SPARK-3's.
4. **Whether the host offset is cell-dependent or a single unlucky dgx draw is NOT resolved**
   (§2.2, n = 1 on dgx). Resolving it needs replicate dgx null blocks — **reserved share, not
   mine to run.** Flagged, not attempted.
5. **A number that should stop travelling:** SPARK-3's "+13 % bar" is problem-specific (§3).

```
add-list (stage, do not commit):
  LOTTERY_recut_band.md
  LOTTERY_recut_band.png
  reports/lottery_recut_band.npz
  hpc_harbor/lottery/lottery_recut_band.py
  hpc_harbor/lottery/lottery_recut_band_figure.py
  hpc_harbor/lottery/lot_band.sbatch
  (raw realisation banks: LOTTERY_results/lotb{02,04,16}_*.npz — 720 lean npz — plus the
   basis_B*.json provenance sidecars, on ACCRE disk, GITIGNORED by the repo's `*_results/`
   rule exactly as CHORUS_results; the committed bank is reports/lottery_recut_band.npz)
```

**STOP.** The decision is executed: multi-basis ensemble run, grade quoted as a band under the
pre-registered reading, convention proposal folded in, SPARK-3 convergence noted (with one
correction travelling back). The tier-1/tier-2 staging stands as handed off. **No further submits
without a fresh brief.**
