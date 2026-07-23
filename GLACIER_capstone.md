# GLACIER — the capstone: does the array EAT THE BACKGROUND?

*Agent GLACIER, ACCRE. Dev/doc authority. Code STAGED (Matt commits). No git commit/push/tag
from here.*

**STATUS: STAGE 0 — GATES. `g0` CLOSED (H200, all pass). `g1` CLOSED (independently re-derived).
`g2` SCIENCE gates PASS; fit gate owed on a smoke-N mismatch (§g2.2). Blocking the wide launch:
(1) a CONCURRENT second GLACIER session (§g2.4), (2) the g2 fit-gate remedy, (3) BUDGET 1.8–3.3×
the STOP (§B). NO CAMPAIGN GPU SPENT; nothing committed/moved/deleted.**

---

## HEADER — the campaign's declared conditions

| item | value | authority |
|---|---|---|
| **Floor host, PINNED** | `hgx03` / NVIDIA **H200**, `batch_gpu`, `cpus-per-task=8` | g1 below; SPARK-3 §4.5 cut its verdict on this lane |
| **Lane** | general queue (`batch_gpu`, `taylor_group_acc`) | brief. The reserved `dgx_iacc` share is SPARK-3C's and is **NOT touched** |
| **Licence** | SPARK-3 §5.3, **EDGE-POSITIVE** | verified in source, not inherited — see §L |
| **GWB basis** | **`RECOMPUTED-UNSAFE`, host-pinned** — there is no banked `L_gwb` at T≠15 | §g0.2 |
| **Budget STOP** | ~150 GPU-hr | brief. Current estimate of the literal Stage 1+2 is **~350 GPU-hr** — see §B |

---

## L. THE LICENCE — verified, because two tracked files say the opposite

Two tracked documents state GLACIER is **UNLICENSED**:

- `reports/FORGE_G_build.md:35-36` — *"GLACIER stays UNLICENSED (SPARK-3 §5.3: STRADDLED is not
  EDGE-POSITIVE)"*
- `reports/SPARK2_second_arrow.md:166-168` — *"GLACIER … is NOT pre-registered here, because its
  licence was an EDGE-POSITIVE verdict and no verdict was reached."*

**Both predate the split.** `reports/SPARK3_second_arrow.md:1001` is the operative text:

> ### 5.3 GLACIER — **LICENSED** by the EDGE-POSITIVE split, WITH the model-quality law as a design gate

and §5.1a records the split itself — *"4 of 5 units survive at the true marginal width. ≥2
required. Scrambled-clean. → EDGE-POSITIVE."* The licence condition is met. **Checked in the
source rather than taken from the brief**, because a campaign that spends 150 GPU-hr on an
inherited adjective is the failure mode this programme keeps writing reports about.

### L.1 The pre-registration GLACIER is held to, and where this brief DEVIATES

SPARK-3 §5.3 pre-registers GLACIER with specific conditions. The brief differs in four places.
Declared here rather than discovered later:

| §5.3 requires | this brief | disposition |
|---|---|---|
| **fixed list** — *"the source list **never grows** — no recruitment, no self-found step"* | a drawn population with a moving **frontier** and `N_resolved(iter)` rising | **COMPATIBLE, and the reason is load-bearing — see L.2** |
| ≥8 realisations **at both venues** | one venue (`lit`) | **DEVIATION.** Brief overrides; scoped in every number |
| **both Fisher bounds** | not specified | **RESTORED** — the model-quality gate needs the marginal width anyway, so both bounds are free on the scoreboard |
| scrambled arm **at every rung** | scrambled at 2 skies, ≥5 reals | **DEVIATION.** Brief overrides; scoped |
| per-iteration floors, never inherited; trials cost restored (+0.578 nat, saturating) | same | **HELD** |

### L.2 The frontier is not a search — and this is the campaign's whole safety case

§5.3's safety argument is that *"the list is frozen and every source enters at its true
position"*, which closes EMBER's manufacturing regime (motion under a wrong counterpart,
sensitivity 1.00, p = 0.002, S8.5.3) **by construction**, so B0.2's search gap (S8.5.5 — self-found
source error ~0.5 scaled against a ~1e-4 certification tolerance, 3–4 orders, confident-wrong)
never enters. That gap is what **stopped AVALANCHE pre-flight at zero GPU**.

A campaign whose `N_resolved` grows could re-open it. It does not, **provided one construction
rule holds**, which is therefore pre-registered as a hard rule and gated, not assumed:

> **THE FROZEN-CENSUS RULE.** The census is drawn once per sky and never grows. A frontier update
> changes only a source's **feed weight** `w_s` (FORGE-G's SMASK) — `carried → fed`. Every
> component enters at its **true drawn position**. There is no search step, no self-found source,
> and no counterpart to get wrong. `N_resolved` rises because the *fed set* grows inside a frozen
> census, never because the census does.

Under this rule "sources → distances → more sources → more distances" is a re-partition of known
sources between the discrete template and the collective GP — exactly SMASK's `fed / carried /
pinned` trichotomy (FORGE_G §2.2) — and §5.3's safety case survives intact. **If any frontier
update ever promotes a source at a fitted rather than a drawn position, B0.2 applies and the
campaign STOPs.** Gate `g2c` below measures this rather than trusting it.

---

## g1 — THE C2 DISCREPANCY: **CLOSED**. It is a provenance bug, not a stale read and not (only) a host effect.

**The inherited gate.** Five rungs read `UNREACHABLE, N_cert < rung` on H200 (Jul 17) but computed
`N_cert = 10` on dgx (Jul 18): `s3r_A_g3_r5_k8`, `s3r_B_g3_r3_k5`, `s3r_B_g3_r3_k8`,
`s3r_B_g3_r4_k8`, `s3r_B_g3_r5_k8`.

**THE MECHANISM — one missing flag, and it is neither hypothesis the brief offered.**

`certified_of` (`hpc_harbor/spark/spark3.py:838`) branches on `self_lane` and returns a different
floor, which sets the certified set through `cert = (dlnL > max(lnK, floor)) & (qmax > QBAR)`:

| launcher | flag | lane | floor A | floor B |
|---|---|---|---|---|
| `sp3_h200fan.sbatch` (Jul 17) | `--self` | **SELF** (`spark3_venueself_*`) | **133.844** | **49.564** |
| `sp3_fan.sbatch` (Jul 18) | *(absent)* | **SURFACE** (`spark3_venue.npz`) | **122.461** | **44.397** |

Verified three independent ways, none of them the report's narrative:
1. **Source** — the branch at `spark3.py:845`; `sp3_h200fan.sbatch:61` ends its `spark3.py unit`
   line with `--self`, and `sp3_fan.sbatch` contains no `--self` anywhere.
2. **Logs** — the two runs print their own floors, at *identical seeds*:
   `Jul-17 … noise 20560400 dist 30560400  N_cert(own) 9  floor 133.844` versus
   `Jul-18 … noise 20560400 dist 30560400  N_cert(own) 10  floor 122.461`.
3. **Re-derivation from raw columns** — `N_cert` recomputed from `dlnL_det`/`lnK`/`qmax` under
   both floors reproduces the published 2×6 table **in all 12 cells, zero discrepancies**, and
   accounts for **5 of 5** named files (each unreachable under SELF, reachable under SURFACE).

**Disposition per the brief's fork — BOTH branches apply.**

- **Driver bug → named.** It is a *wrong-bank-path* bug, delivered by naming: `_rung_path`
  (`spark3.py:1070`) emits `s3r_{key}_g{geo}_r{real}_k{k}.npz` with **no lane tag**, while
  `_lane_tag()` (`spark3.py:535`) — *"An artifact whose name does not distinguish the conditions
  it was produced under is not an artifact, it is a race"* — already guards the `venueself` and
  `replay` artifacts one file away. Two provenances share one directory.
  **GLACIER's driver does not share this path**: a separate `GLACIER_results/`, and `_lane_tag()`
  applied to *every* checkpoint and floor bank from iteration 0. Gated by `g2d`.
- **Host-dependent floors → host-pinning ADOPTED anyway**, because the underlying offset is real:
  three same-host null blocks land 133.844 / 141.231 / 140.891 (mean 138.66, scatter 4.17 against
  a per-estimate error ~4.79), **3.0σ from SURFACE's 122.461**. GLACIER therefore declares its
  floor host in this header, carries host/GPU/`fp=`/`cpus` metadata on every floor bank, and
  quotes near-bar grades **MARGINAL-with-band**.

### g1.1 TWO REFINEMENTS TO §4.5, found in the re-derivation — one is a live defect

**(a) The gap is floor *and* basis, not floor alone.** §4.5 presents this as a floor story. The
SELF bank carries its **own** `dlnL`/`lnK`/`qmax`, re-derived on the H200's own GWB basis. The
floor-swap-only control — SURFACE columns re-thresholded at the SELF floor — gives
`A 8,11,14,9,10,7 / B 10,8,11,5,7,7`, which differs from the true SELF column at **3 of 12 cells**
(A r0, B r2, B r3). The floor is the dominant term; it is not the only one. Both published columns
are individually exact against their own lane's raw data.

**(b) LIVE DEFECT — the venue-A SELF floor is selected by filename sort order.** `certified_of`
takes `sorted(glob(spark3_venueself_{key}_g3_*.npz))[0]`. **Three** co-resident venue-A banks
exist, and they are the same three blocks §4.5 quotes as evidence the offset is systematic:

| bank | floor | N_cert A r0..r5 |
|---|---|---|
| `…_hgx03_NVIDIAH200.npz` *(selected — `.` sorts before `_`)* | 133.844 | 9, 10, 14, 9, 10, 7 |
| `…_hgx03_NVIDIAH200_b1.npz` | 141.231 | 8, 9, 13, 9, 9, 5 |
| `…_hgx03_NVIDIAH200_b2.npz` | 140.891 | 8, 9, 13, 9, 9, 5 |

§4.5 treats the three as a *scatter measurement*; it does not note that the code **silently picks
one of them by lexical accident**, and that the choice moves `N_cert` and hence rung reachability.
This is the same provenance-hazard family a fourth time. **Reported, not fixed** (HARD RULE) —
and it is a second, independent reason GLACIER writes its own banks with `_lane_tag()`.

**Consequence for the five files: they are STILL UNQUARANTINED.** All 48 `s3r_*` are on disk under
their original names; §4.5's recommended `mv … .SURFACElane.npz` was explicitly *"not done here"*.
Anyone running `mode_ledger` (which globs `s3r_*`, `spark3.py:1710`) still folds two lanes and
still flips STRADDLED → EDGE-POSITIVE on 8 cross-lane crossings. **GLACIER has not touched them**
(they sit under SPARK-3C's pre-registered STOP and are not GLACIER's to move); flagged for Matt.

---

## g0 — FORGE-G's MACHINERY ON THIS CLUSTER'S GPUs: **CLOSED, ALL PASS**

FORGE_G_build.md §3.3 ran the full B5 suite on cronus CPU (all `0.000e+00`) and auto-deferred the
GPU path (that box's GPU DNN library fails to init). The brief makes the GPU re-run a **blocking**
gate. Ran on the pinned floor host (`hgx03`, H200, `batch_gpu`, `cpus=8`), job `12635733`, wall
1267 s:

```
[Bg2] model-quality gate limits  -> PASS   (fed 0/3 @thr0, 3/3 @inf, 2/3 @default; monotone; delta)
[Bg4] EM conditioning free_dim   -> PASS   (24 / 22 / 21 / 16)
[Bg5a] always-fed+truth-pinned == incumbent B1Marg loop:
        |lnL_marg(always-fed) - incumbent| = 0.000e+00   (bit-exact)
        max|LNL - incumbent| = 0.000e+00 ; max|q_max - incumbent| = 0.000e+00
        fast-full vs masked = 0.000e+00 (<1e-8)                 -> PASS
[Bg5b] zero-noise fixed point: max|src(fp)-truth| 0.000e+00 ; census ceiling 0.0004 (<0.02) -> PASS
[Bg5c] gate limits + w_s==0 bit-exact absence: |lnL(res@truth)-lnL(res@junk)| = 0.000e+00 -> PASS
=== SOFT SOURCE SIDE GATES: ALL PASS ===
```

**Every load-bearing difference is `0.000e+00` on the H200, reproducing cronus.** As pre-declared
in the sbatch (a scope line, not a hedge): these are within-machine *contrasts*, so they are
bit-exact by construction of the gate; the absolute lnL values are not compared cross-host (CPU/GPU
reduce in different orders — arithmetic, not a defect, SPARK3 §4.3). All gates are zero-noise, so
none touches the `L_gwb` basis hazard — which is why they are lane-portable and the campaign's
noisy floors are not. **SMASK is proven live on this cluster's GPUs.**

---

## g2 — POPULATION BUILD: science gates PASS; fit gate owed on a smoke-N mismatch; a CONCURRENCY hazard surfaced

### g2.0 TWO IMPLEMENTATIONS, TWO READINGS — and a correction to my own first framing

There are **two g2 modules on disk**, from two GLACIER work-sessions, and they encode the two
readings of "the background." This section states the relationship plainly, and corrects a framing
I initially wrote and then found too strong.

- **`glacier_pop.py` (the brief's reading, the primary build)** — draws N~200–500 binaries,
  **NORMALISES** the drawn set so its band-integrated characteristic-strain power equals the NG15
  target (`Σ h_i² = A² ∫ (f/f_yr)^{-4/3} df/f`, `log10 A=−14.6`), injects the whole population as
  **deterministic delays with NO stochastic GWB in the data** ("the background IS the unresolved
  sum"), and **FITS the `gwb` GP amplitude** by profiling `lnL(log10_A)` (re-factoring `Pinv(A)`
  per grid point). `A_background(iter)` is the drain. It carries a `PromoteLedger` that enforces
  the frozen-census rule (§L.2) *in code* — a promote at a fitted position raises `CampaignStop` —
  and a `bank_npz` that lane-tags every artifact and refuses cross-lane collisions (the g1 lesson,
  encoded). **This is a sound, self-consistent realisation of the brief's actual question**, and it
  is the primary g2.
- **`glacier_population.py` (a §5.3-reading cross-check I wrote)** — the same draw + a
  reading-independent power-bookkeeping (`A_background := sqrt(Σ carried P)`), gated on CPU. Its
  value is the *incumbent-census* view (§g2.3) and the Stage-0 figure; it is subordinate to
  `glacier_pop.py`, not a competing primary.

**The correction.** I first framed the brief's reading as a "premise conflict" that severed
comparability with the fixed-`A=−14.6` floors and should be deferred to §5.3. **That was too
strong.** `glacier_pop.py`'s construction — inject the sources as the data, put *no* stochastic GWB
in it, and fit the GP amplitude to soak up the unresolved sum — is internally self-consistent and
is gated on its **own** terms (recover the drawn `A_eff`), in its own `GLACIER_results/` with its
own lane tags. It does not claim comparability with the banked floors; it is a different,
well-posed experiment, and it is what the brief asks for. The genuine open item is the **budget**
(§B), not the generative model.

### g2.1 WHAT RAN, AND THE VERDICTS (job `12645842`, hgx03 H200; and `12645876`)

`glacier_pop.py`'s gates ran on the H200. The **science gates PASS**, at the production draw:

```
g2a-i  DRAW NORMALISATION (N=256): Σ h_i²/target − 1 = 2.2e-16 (exact by construction);
       |log10 A_eff − A_target| = 0.0048 (<0.01)                                    PASS
g2b    CONSERVATION across all 257 frontier positions: max rel dev 3.3e-16 (<1e-12)  PASS
g2c    FROZEN-CENSUS LEDGER: promote@truth bit-exact; promote@fitted → CampaignStop  PASS
g2d    PROVENANCE: metadata in bank, readback identical, cross-lane overwrite refused PASS
```

**g2a-i is the brief's iteration-0 gate, met at production N:** the drawn population reproduces the
NG15-class summed power to 0.005 dex, and total (resolved + residual) is conserved across the whole
frontier by construction. **g2c enforces §L.2's safety case in code** — the frozen-census rule is
not a promise, it is a `CampaignStop`.

### g2.2 WHAT DID NOT CLOSE — the GPU FIT GATE (g2a-ii), and why: a smoke-rung N mismatch

The sbatch runs a **T=15, N=32 smoke** before the T=30, N=256 fit. At **N=32**, g2a-i **FAILS**:
`|log10 A_eff − A_target| = 0.064 > 0.01`. This is **not a defect** — it is single-draw
discreteness: 32 sources cannot tile a smooth γ=13/3 powerlaw over 24 bins to 0.01 dex, whereas 256
can (0.005). The `0.01` normalisation gate is **N-agnostic and unachievable at the smoke's N**, so
the smoke fails and **blocks the GPU fit gate from ever running**. The load-bearing g2a-ii —
*does the fitted background amplitude recover the summed power on the real 116×N_gp joint?* —
therefore has **no verdict yet**. The one-line remedy (raise the smoke N, or scale the g2a-i
tolerance as ~1/√N — discreteness is a √N effect) is **named, not applied**: see §g2.4.

### g2.3 THE §5.3-READING VIEW, and the number it delivers — the incumbent reservoir is 4%

Under the fixed-list §5.3 reading (`glacier_population.py`, CPU, ALL PASS), on the **incumbent
3+13 census** the 3 loud members carry **95.85%** of the residual power; the faint reservoir is
**4.15%** (= **20.4%** in amplitude, `A_background/A_total` after the loud tier resolves —
`GLACIER_g2_drain.png`). **This is the crux of why the reading matters:** the §5.3 frozen census has
little to drain (a 4% reservoir — the arithmetic face of SPARK-3's kill-shot, *"+1 at the
optimistic bound → practically irrelevant"*), whereas the brief's reading normalises N~256 sources
to the **full** NG15 background, so the drainable object is the entire unresolved sum. **The two
readings differ not in method but in how much background there is to eat** — and Stage 2b (into the
knee) is where even the §5.3 reservoir grows. This is a real Stage-0 result under either reading.

### g2.4 THE CONCURRENCY HAZARD — a SECOND GLACIER session is acting on this repo, and it must be reconciled

While closing g2 I found **jobs submitted in `milesmt` that this session did not launch**
(`glg2` 12645842 at 16:33, 12645876 at 16:37) and **`GLACIER_results/` cleared and rewritten**
(dir mtime 16:37) outside my actions — i.e. a **second GLACIER work-session is running
concurrently** on the same tree. It authored `glacier_pop.py` + `gl_g2.sbatch` and is resubmitting
the fit gate. **Per this programme's own core lesson (g1: "two provenances must not share one
artifact"), two agents editing one tree and one queue is that hazard at the session level.** I have
therefore **stopped acting on the shared queue and shared files**: I cancelled my one duplicate
submit (`12645878`), and I did **not** edit `glacier_pop.py`/`gl_g2.sbatch` (in active use by the
other session), delete any bank, or resubmit. **Reconciling the two sessions is Matt's to do**
(§DECISION) — including which of `glacier_pop.py` / `glacier_population.py` is kept. I flag it rather
than race it.

---

## B. BUDGET NOTE — from SPARK-3's measured per-unit cost, before the wide launch

Measured SPARK-3 H200 rung cost (`sp3hf_*.out`, `exit=0 wall=`): mean **1532 s = 0.43 GPU-hr** per
`(venue, real, rung)` at `n-null=100`. A GLACIER iteration is one such certify-and-refloor pass
plus a soft M-step over the census (~0.12 GPU-hr) → **~0.55 GPU-hr/iteration**.

| arm | cells × iter | GPU-hr |
|---|---|---|
| S1 primary (8 sky × 4 igniter) | 32 × 6 | 105 |
| S1 scrambled null (2 sky × 5 real) | 10 × 6 | 33 |
| S2a time-ladder (2 sky × 3 T × 4 igniter, T=40 ×1.5) | 24 × 5 | 98 |
| S2b knee (2 sky × 4 igniter, denser ×1.4) | 8 × 6 | 37 |
| **LITERAL TOTAL (mean 6 iter)** | | **272** |
| **… to the 12-iteration cap (loops run long)** | | **~490** |

**The literal campaign is 1.8×–3.3× the brief's own ~150 GPU-hr STOP.** Per the brief
(*"STOP at ~150 GPU-hr for re-scope"*), the wide launch does not proceed as written.

A re-scope holding the **full primary arm** (the verdict-bearing 8 skies × 4 igniters) and the
scrambled null, trimming the iteration cap to 5/4 and the time-ladder to 2 igniters, lands at
**148 GPU-hr** — but it defers S2b, the knee, which §g2.2 just argued is the arm most likely to
make the drain non-trivial. **The re-scope is therefore itself a science decision, not a mechanical
trim**, and is surfaced rather than chosen.

---

## THE CAMPAIGN AS IT WOULD RUN (pre-registered; behind the DECISION)

Under the §5.3 reading, on the built + gated machinery, the loop is: per `(sky, igniter)`, iterate
`[E] certify → [D] certified-coherent re-score of the frozen reservoir → [M] soft-joint tighten
(MargJax, SMASK per the model-quality gate) → frontier update → repeat` to fixed point or the
iteration cap. Per `(sky, igniter, iteration)` BANK — with the frozen-census rule (§L.2) enforced
and floors **refit per iteration-state, never inherited, host-pinned (§g1), zero-fraction column
required**: `N_cert`, `N_resolved`, `A_background` (THE DRAIN CURVE), per-component concentration,
per-source localisation (host-unique crossings flagged EPOCH, report-only), channel budget,
**wrong-certs (must stay 0)**, trials cost restored (+0.578 nat, saturating — SPARK-2 §2). The
**scrambled-counterpart** arm runs the full loop; any loop-grown wrong cert or spurious resolution
→ STOP + anatomy.

**Pre-registered verdicts (per arm):** COMPOUNDING (per-iteration gains grow/hold ≥3 iter; report
what saturates) / CONVERGENT (gains shrink to a fixed point; report its height above ignition and
the shortfall to the next recruit) / INERT (nothing beyond iteration 0; report nearest-miss
margins). Plus the **DRAIN verdict**: does `A_background` measurably fall as sources resolve,
quoted with its error. **Kill-shot (SPARK-3 §5.3, adopted):** a fixed point of `+1` at the
optimistic bound = arithmetically alive, practically irrelevant — that sentence closes GLACIER
rather than a campaign.

---

## DECISION — surfaced for Matt (I am not proceeding past Stage 0 autonomously)

Three items block the wide launch; none is mine to resolve by fiat:

1. **CONCURRENCY (do first).** A second GLACIER session is acting on this tree and queue (§g2.4).
   Two agents on one repo is the g1 hazard at the session level. **Reconcile the sessions** —
   decide which is live, and which of `glacier_pop.py` (the brief-reading primary) /
   `glacier_population.py` (my §5.3 cross-check) is kept — before any further submits.
2. **Close the g2 fit gate.** The science gates PASS (§g2.1); only the GPU fit gate (g2a-ii) is
   owed, blocked by an **N=32 smoke rung** that cannot meet the N-agnostic 0.01 normalisation gate
   (§g2.2). One-line remedy: raise the smoke N (or scale the g2a-i tolerance ~1/√N). ~1 GPU-hr,
   zero-noise. **Not applied here** — it edits a file the other session is actively running.
3. **Budget.** The literal Stage 1+2 is **272 GPU-hr (up to ~490 at the iteration cap)** — 1.8–3.3×
   the brief's own ~150 STOP (§B). Release the full campaign, *or* the **148-GPU-hr re-scope**
   (full primary arm + null; ladder trimmed; **knee deferred** — which §g2.3 argues is the arm
   most likely to make the drain non-trivial, so the trim is a science call), *or* a middle scope.

**What is DONE and needs nothing further:** g0 (H200, all pass), g1 (closed + a live defect
reported), g2 **science** gates (draw-normalisation at production N, conservation, frozen-census,
provenance — all PASS), the budget analysis, and the §5.3-reading Stage-0 drain figure.
**Owed before the wide loop:** the concurrency reconciliation, the g2 fit-gate verdict, and the
budget release.

---

## STOP

Stage 0 stands on its own terms: **g0 and g1 — the two blocking GPU gates — are CLOSED with file
evidence; g2's SCIENCE gates PASS** (draw sums to the NG15 band power to 0.005 dex; conservation
exact across the frontier; the frozen-census safety rule enforced as a `CampaignStop`). In the
course of it I established (a) a **live provenance defect** in `certified_of`'s SELF-lane bank
selection (§g1.1b, reported not fixed — HARD RULE); (b) that the brief's drainable-background g2 is
**built and sound** (`glacier_pop.py`, prior session) — I **retract** my first "premise conflict"
framing (§g2.0); (c) the g2 **fit gate is owed** on a smoke-N gate mismatch (§g2.2); (d) a
**second GLACIER session is running concurrently** on this tree and queue (§g2.4); and (e) the
literal campaign is **272 GPU-hr, 1.8–3.3× the STOP** (§B).

**No campaign GPU spent** (g0 and the g2 gates are zero-noise Asimov, not campaign realisations).
**No noisy realisation drawn. Nothing committed. I moved/deleted nothing, and — on finding a
concurrent session — I stopped acting on the shared queue and shared files rather than race it.**
The wide launch waits on Matt's three decisions above — as the convention requires (you stage, Matt
commits), and as this programme's pre-flight discipline requires when the gates, the budget, and now
a second session all say *stop and reconcile*.

---

## ADD-LIST (staged; Matt commits)

```
GLACIER_capstone.md                        new   this report (Stage 0: gates + decisions)
GLACIER_g2_drain.png                       new   Stage-0 §5.3-reading: census power hierarchy + drain
hpc_harbor/glacier/gl_g0.sbatch            new   g0 — FORGE-G B5 suite on H200 (ALL PASS)
hpc_harbor/glacier/glacier_population.py   new   g2 §5.3-reading cross-check (CPU): drain bookkeeping
hpc_harbor/glacier/glacier_g2_figure.py    new   the Stage-0 drain figure generator
reports/glacier_g2_population.npz          new   lean bank: incumbent census, per-source power, drain
```

**Authored by the concurrent session, NOT me — do not double-stage; reconcile first (§g2.4):**
`hpc_harbor/glacier/glacier_pop.py`, `hpc_harbor/glacier/gl_g2.sbatch`, and any `GLACIER_results/*`
(gitignored). `glacier_pop.py` is the brief-reading **primary** g2 and is the better base for the
campaign; my `glacier_population.py` is a subordinate §5.3 cross-check kept for its figure.

Nothing in `SPARK3_results/` was touched (the five cross-lane rungs remain under SPARK-3C's STOP;
§g1.1 flags them and the `certified_of` defect for Matt). No tracked code was edited.

## STILL OWED / named (not dropped)

- The **three decisions** above — concurrency, the g2 fit-gate remedy, budget — before any wide GPU.
- The two headline figures owed at campaign end: the DRAIN curve `A_background(iter)`, and from
  Stage 2b `N_resolved(iter)` + `A_background(iter)` **against the knee line**. `GLACIER_g2_drain.png`
  is their zero-noise precursor.
- **`L_gwb` at T=20/30/40 is host-pinned `RECOMPUTED-UNSAFE`** (fp `8548f148b50a5b44` at T=40, per
  SPARK-3 §5.4 / EMBER) — no banked basis off the incumbent T=15. Any time-ladder (S2a) is portable
  only to hosts reproducing that fingerprint; declared so it is not rediscovered mid-run.
- The `certified_of` SELF-lane filename-sort defect (§g1.1b) and the five unquarantined rungs
  (§g1.1) — both SPARK-3C's code, reported for Matt.

---
---

# SESSION 2 ADDENDUM — the brief-carrying resume (written after the sections above went quiescent)

*Agent GLACIER session 2 — the session the other author detected in §g2.4. Everything above this
line is the concurrent session's text and is left intact; this addendum reconciles rather than
edits it. I am the author of `glacier_pop.py`, `gl_g2.sbatch`, `glacier_loop.py`, `gl_fan.sbatch`,
and jobs `12645816/42/76/12645909`.*

## R. RECONCILIATION — two of the three DECISION items are already answered BY THE BRIEF

Session 2 was launched by Matt with a resume brief written **after** the Stage-0 readback, and the
brief is the paper trail for DECISION items 2 and 3:

1. **CONCURRENCY — the one live item, Matt's to close.** The brief opens *"You are agent GLACIER
   (session 2) on ACCRE, resuming the capstone campaign after a session loss."* Session 1 was
   believed lost mid-Stage-0; it turned out alive and flushed its report + artifacts at 16:30–16:42
   during this session's g2 runs. Both sessions independently detected the collision and both
   applied the same discipline (it stood down from the shared queue; I paused capstone writes until
   quiescence). **Recommend: confirm session 2 as the live session; retire session 1.** Its work
   product is preserved above and its `glacier_population.py` / drain figure stay staged as the
   §5.3 cross-check.
2. **THE g2 READING AND THE FIT GATE — resolved by the brief's own g2 clause.** The brief
   re-specifies g2 verbatim as *"a GP whose amplitude is a FITTED parameter tied to the unresolved
   sum (the drainable background)"* with the iteration-0 reproduction gate — i.e. Matt confirmed
   the fitted-GP reading after the premise-conflict readback (the §g2.0 conflict was raised, and
   §STOP above already retracts it). The smoke-N mismatch is **fixed the principled way** (the
   projection-vs-target check is a population-level statement: gated at n ≥ 200, report-only
   below; the fit is always cut against `A_eff`, well-defined at any N) and the fit gate is
   **running as job `12645909`** (T=30, ncw=256, this lane). Verdict recorded in §R.2 below when
   the artifact lands.
3. **BUDGET — resolved by Option 1.** The 272-GPU-hr table above prices a **4-igniter** Stage 1.
   The brief records Matt's decision as **TWO igniters** (e=0.7 single + all-circular NONE, the
   escalation to e=0.5 / e=0.3×2 behind a first-sky readback authorisation), the held floor line
   as **full floors on the verdict-carrying subset only** (both scrambled cells + 2 calibration
   skies per igniter; others refit at {0,3,6} behind a drift gate), and the TURBINE adoptions
   (2-way packing, shared warm cache, runtime-SMASK). Under that scope: 36 cells (16 signal + 20
   null) × 6 iters, full-floor iteration ~0.3–0.45 GPU-hr (TURBINE rung class at T=30), carried
   ~0.05–0.1, ÷1.33 packing → **~53 GPU-hr Stage 1, ~75–85 with extension riders + Stage 2** —
   the brief's own number, reproduced independently. **The 150 STOP stands with ~2× headroom.**

## R.1 STAGE-1 BUILD STATE (held; nothing submitted beyond gates)

`glacier_loop.py` + `gl_fan.sbatch` (from TURBINE's T3 template): the soft-joint iteration as
pre-registered — Fisher both bounds (columns, never merged) → ModelQualityGate (ratio 0.5,
declared) → PromoteLedger **at drawn truth** (fitted-position promote = `CampaignStop` = B0.2) →
BackgroundFit drain refit at an **explicit SMASK always** (an absent key means all-fed and silently
zeroes the residual — caught before any GPU ran) → E-step → M-step (dynamic fed slots, numpy-side:
no recompile on frontier change) → end-of-cycle v2.2 scoreboard (+0.578 nat trials restored) →
lane-tagged bank + checkpoint resume. Igniter arms paired on census + noise seeds (CHORUS C6);
e=0.7 comb per CHORUS C1–C3; **scope line: the igniter is the only eccentric waveform in data or
template — the drawn population e (LOTTERY realistic point, parameterised) is banked, report-only.**
Scrambled-null STOPs are code: null-arm promote (spurious resolution) and null-arm wrong cert.
2-way packing keeps the 8-cpu convention by `taskset`-splitting a 16-cpu allocation — each drawing
process sees affinity exactly 8, and the module **refuses to run** otherwise (LOTTERY-recut's
affinity refusal, adopted).

**THE HOLDS, all enforced in code (`cell`/`null` modes refuse to run):**
(i) FORGE-G2 runtime-SMASK pull + re-gate on this lane (expect all Bg5 PASS + bit-exact flip 0.0 +
warm flip <10 ms); (ii) **driver gate G-d2** — the per-target E-step scoreboard wiring
(`_cert_columns` is a **declared stub**; a fan launched today would bank empty scoreboards, and the
refusal exists so that cannot happen silently); (iii) the resume drill.

## R.2 FINDING FOR THE FIRST-SKY READBACK — the realistic point is FAINT (report-only)

At the NG15-summed normalisation the drawn brightest member is **`log10_h ≈ −14.73`** — ~30× in
amplitude (~10³ in power) below the `−13.25` census loudness of every prior onset result, which
STORY S7.3.3 itself puts ~0.35 dex *above* the NG15 sky-averaged limit; S7.3.4's population clock
(N̄ ≲ 0.01 detectable) says the same from the literature. **The iteration-0 frontier may therefore
resolve zero members and Stage 1 may legitimately read INERT on `N_cert`** — a pre-registered
verdict class, with the DRAIN CURVE and nearest-miss margins as the live measurement. This is the
same physics as the concurrent session's 4%-reservoir finding under the incumbent census: both say
the drain, not the count, is where the realistic-point information lives — and both point at
**Stage 2b (into the knee)** as the arm that can make it non-trivial. A population re-inflated to
the old census loudness would measure a louder universe than NG15 permits; not done.

*(g2a-ii verdict from job `12645909` is appended below when the bank lands.)*

## R.3 ADD-LIST ADDITIONS (session 2; Matt commits)

```
hpc_harbor/glacier/glacier_pop.py     new   g2 PRIMARY (brief reading): population synthesis,
                                            fitted-GP BackgroundFit, PromoteLedger, gates g2a-d
hpc_harbor/glacier/gl_g2.sbatch       new   the g2 gate job (smoke T=15/n=32 + main T=30/n=256)
hpc_harbor/glacier/glacier_loop.py    new   Stage-1 soft-joint loop driver (HELD: G-d2 stub +
                                            runtime-SMASK re-gate + resume drill before any fan)
hpc_harbor/glacier/gl_fan.sbatch      new   Stage-1 fan (2-way taskset packing, 8-cpu convention,
                                            cache/quota prologue) -- DO NOT SUBMIT until holds clear
GLACIER_capstone.md                   mod   this addendum
```

---
---

# SESSION 3 ADDENDUM — the resume that closed g2a-ii's anatomy (2026-07-22)

*Agent GLACIER session 3, resuming per the brief ("inventory before building, from disk and
squeue"). Everything above is prior sessions' text, left intact. HEAD at resume: `2b7ce2f`
(pull was already-up-to-date). Queue at resume: EMPTY — no prior session's jobs survived.*

## I. INVENTORY (the brief's four questions, answered from disk)

| item | state found |
|---|---|
| (a) Bg5a / g0 | **CLOSED, artifact readable** — `glg0_12635733.out` re-read end-to-end: `[Bg5a/b/c] PASS`, `=== SOFT SOURCE SIDE GATES: ALL PASS ===`. No re-run needed. |
| (b) g2 | Build exists (`glacier_pop.py`, primary). Science gates PASS (§g2.1). **The fit gate g2a-ii: FAIL, twice** — jobs `12645909` AND the post-addendum re-run `12646058` (band-matched) both read `Ahat = −15.6 ± inf [EDGE HIT]` vs `A_eff = −14.54`. Session 2's §R.2 slot ("appended when the bank lands") was never filled: **the bank landed FAIL at 17:00 Jul 20** (`g2gate_fit_T15_s9000000_n32__hgx03_NVIDIAH200.npz`, verdict column FAIL) and the T=30/n=256 main fit never ran (fail-fast). |
| (c) Stage-1 remnants | `glacier_loop.py` + `gl_fan.sbatch` (HELD, refusal-as-code confirmed in `main()`); no checkpoints, no `gl1_*` banks, nothing queued. |
| (d) FORGE-G2 runtime-SMASK | **NOT LANDED.** Read in code, per the brief: `trackB_b1_core.py:372-390` — `set_smask` still closes smask over `_params` inside the jitted evaluators and **invalidates every compiled per-pulsar evaluator on change** (the docstring says so in its own words). Last commit touching the core is `f4e2a85`. The g0 HARD HOLD therefore STANDS: build everything else, wait for Matt's "pushed". |

## II. g2a-ii ANATOMY — CLOSED: it is neither a plumbing bug nor a tolerance mis-set. THE VENUE IS DEAF.

The banked FAIL profile is **pure Occam**: `lnl(−15.6) − lnl(−14.6) = +3.8` nats, monotone
across the whole 2-dex grid, no visible quadratic gain anywhere. Two hypotheses produce that
curve — H-BUG (the `__smask` key dropped between `BackgroundFit._data_terms` and `MaskedDelay`;
template full; zero-noise residual identically zero) and H-FAINT (residual real but below the
in-band noise). They demand opposite remedies, so a forensic measured which
(`glacier_g2_forensic.py`, job `12706051`, hgx03/H200, zero-noise gate-class):

**F1/F2 — H-BUG REFUTED at every layer.** `max|template(smask=0)| = 0.000e+00` exactly;
`‖FtNmy‖ = 4.045e+07` with the population under `smask=0` versus `0.0` under all-fed/absent-key
(the `1.0*x == x` identity bit-exact at `max|diff| = 0`). The dict flows unfiltered through
`kernelterms → yprod → MaskedDelay`; the machinery is SOUND and the population (median residual
31.2 ns, max 58.3 ns) reaches the fit.

**F3/F4 — H-FAINT CONFIRMED, and it is worse than "the population is faint": the ARRAY is deaf
in the generative band.** Decomposed profile (quad and Occam separately), floor→A0:

| venue | population quad rise | control quad rise | Occam rise | ahat (pop) | ahat (ctrl) |
|---|---|---|---|---|---|
| T=15, n=32 | +0.56 | +0.39 | +8.32 | −17.6 [EDGE] | −17.6 [EDGE] |
| T=30, n=256 | +2.34 | +1.38 | +17.10 | −17.6 [EDGE] | −17.6 [EDGE] |

The **control** is the load-bearing row: a synthetic GP realisation drawn AT the NG15 reference
amplitude from the exact banded prior (correlated, in-band, seed 4242) — i.e. a TRUE background
of exactly the amplitude the gate asks the fit to recover — **also pegs at the extended grid's
edge at both venues**. The venue's honest amplitude precision is unbounded (`sig_ctrl = inf`).
No tolerance, no grid, no estimator fixes that: the 116-psr mock at T=15–30 has essentially zero
Fisher information on a 13/3 GP amplitude confined to the generative box **[10, 32] nHz** — PTA
GWB sensitivity lives at ~2–9 nHz, and the box (inherited from the CW estimator population
convention, `stagec_fisher`) sits entirely above it. An NG15-total population confined there is
invisible AS A BACKGROUND by construction. Banks: `g2forensic_T15_n32__hgx03_NVIDIAH200.npz`,
`g2forensic_T30_n256__hgx03_NVIDIAH200.npz` (grids, both decompositions, both controls).

**Consequence, stated plainly: the pre-registered g2a-ii gate CANNOT PASS at any campaign venue
under the incumbent band, and the DRAIN CURVE `A_background(iter)` as pre-registered would be a
measurement with infinite error bars at every iteration.** This is a Stage-0 blocking finding of
the same family as session 2's §R.2 (the realistic point is faint on the COUNT observable) — the
two together say: at the NG15-normalised realistic point in the incumbent box, there is nothing
measurable to eat, on either axis.

## III. THE REMEDY FORK (Matt's call; pre-analysed, one arm pre-measured)

- **(a) EXTEND THE GENERATIVE BAND down to where the array listens** (flo → −8.7, ~2 nHz; the
  estimator's own search box already reaches −8.7, so no estimator premise moves). The
  unresolved sum then has low-frequency power, the iteration-0 fit recovers it, g2a-ii becomes
  cuttable, and the drain question becomes well-posed. The band is now a DECLARED PARAMETER of
  `draw_population` (banked per draw) rather than an inheritance. Pre-measured at the main-gate
  venue by job `12706074` (readback in §V below). Drain arithmetic from the banked n=256 draw:
  top-1/3/10 members carry 19.7% / 44.8% / 77.6% of population power → A_bg drops of
  0.048 / 0.129 / 0.324 dex if resolved — so a fit σ of ~0.03–0.05 dex makes a top-3 resolution
  a ≥3σ drain detection. SCOPE: population realism of the mock band is a modelling choice,
  declared, not a claim about the real sky.
- **(b) Re-observable the drain as band-limited excess POWER** (quad-based, not 13/3-amplitude):
  REJECTED BY THE SAME NUMBERS — F3/F4 show the total in-band excess is itself sub-unity in
  SNR²; there is no information to re-parameterise.
- **(c) Run Stage 1 with the drain as LEDGER bookkeeping only** (resolved+carried conservation,
  which is exact by g2b): honest, but it answers a weaker question than the capstone's headline
  ("does the array MEASURABLY eat the background") and reduces the DRAIN verdict to arithmetic
  that was true at iteration 0 by construction.

**Recommendation: (a)**, with g2a-ii re-cut at the extended band (same 0.15-dex tolerance IF
§V's measured σ supports it; else tolerance = max(0.15, 3σ_ctrl), declared).

## IV. G-d2 LANDED — the per-target E-step scoreboard is wired (driver hold (ii) cleared, pending its gate)

`_cert_columns` (the declared stub that refused the fan) is replaced by **`CertScoreboard`**
(`glacier_loop.py`): spark3's AUDITED pieces — `estep_per_target` (SPARK-2 §1's structural fix:
target's own term live under its OWN mask; no 0/116 dead rows at rung 0), `rung_masks`
(previous iteration's certified set coherent at their q; uncertified DECOHERED — IGNITE-2's
hard-lock retirement), `score_from_LNL` (ignite.py:349-358 banked certification semantics) —
constructed on GLACIER's own venue (own census fringe spacing `dL` = min mode spacing over ALL
census+harmonic slots; `build_EV`/`FringeTables` at B_CERT=512; tier prior column; Arm-B
distance truth drawn per cell at `dist_seed = noise_seed + 10_000_000`). Carried sources are
stated on this path via H_ABSENT theta slots (gated == SMASK by G-d3), so the scoreboard never
recompiles on a frontier change even on the baked core.

**Three defects found and fixed in the held driver while wiring** (none had run on GPU):
1. `NoiseDrawer(amo)` — wrong constructor arity (it takes `(sp, amo)`); would have
   TypeError'd on first noisy cell. Now built once per cell from the scoreboard's `sp`, and
   `_null_offenders` REUSES it (a per-refit rebuild would re-factor the GWB block per floor).
2. **T=30 L_gwb path**: `run_cell` never took the `force_recompute_lgwb` branch — at T=30 there
   is no canonical `b1_L_gwb` bank (shape 3248 vs 5336) and the loader RAISES. Now forced, with
   the RECOMPUTED-UNSAFE/host-pinned declaration in the header standing.
3. **M-step double-count**: `marg_fn` tiled `theta_rec` directly, leaving CARRIED sources in the
   M-step template while they are simultaneously in the fitted GP — the soft-joint state's one
   forbidden configuration. Now evaluated at `theta_with_absent(...)` (recomputed AFTER the
   frontier promotes — the stage-(a) carried set is stale by then), with only fed-slot values
   written back.
4. (Bonus, session-2 latent) distance truth was never drawn — data were injected at
   `theta_truth`'s default distances, making `on_true` vacuous. Now drawn per cell (Arm-B tier
   convention).

Driver gates re-run on H200 as job `12706069` (mode `gate`: G-d1 plumbing, G-d2a rung-mask
structure, G-d2b rung-0 liveness with a fed member + estep wall time, G-d2c coherent path,
G-d3 identity, drain-profile finiteness — the profile EDGE at the T=15 gate venue is now
REPORT-ONLY, per §II it is physics, not plumbing). Verdict in §V.

**The RESUME DRILL (hold iii) is BUILT, not run**: `gl_drill.sbatch` — interrupted leg (2
iters), resumed leg (3 iters, i0/i1 must resume cold), uninterrupted control leg in its own
subdirectory, then `drillcmp`: the resumed leg's final bank must equal the control's
COLUMN-EXACT (`theta_rec`, `fed_mask`, `a_bg`, `n_cert`, raw columns, floor, zf). DO NOT SUBMIT
until FORGE-G2 is pulled — the drill's bit-identity baseline must be cut on the final tree.

## V. READBACKS LANDED AFTER THE SECTIONS ABOVE WERE WRITTEN

*(filled in this session, same day, before staging)*

### V.1 REMEDY-A MEASURED: at the band the array hears, the fit gate is not just passable — it is PRECISE (job `12706074`, FIT-OK)

Same forensic, same lane, T=30/n=256, generative + fit band extended to `log10 f in
[-8.7, -7.5]` (~2–32 nHz; the estimator's own search box already reaches −8.7, so no estimator
premise moves — only the population's band, now a declared parameter):

| quantity | incumbent band | REMEDY A |
|---|---|---|
| population residual RMS (median) | 31 ns | **570 ns** |
| quad rise floor→A0 | +2.3 | **+1538** |
| Occam rise floor→A0 | +17.1 | +387 |
| `Ahat` (population) | −17.6 ± inf [EDGE] | **−14.629 ± 0.021** |
| `A_eff` (drawn) | −14.595 | −14.583 |
| |Ahat − A_eff| | 3.0 dex, unbounded | **0.045 dex** |
| honest gate precision (`sig_ctrl`) | inf | **~0.03 dex** |

The population recovers its own drawn effective amplitude to **0.045 dex (2.2σ of the profile's
own width)** — the pre-stated 0.15-dex tolerance is comfortably supported, with a peak inside the
grid and finite curvature. One scope line: the TRUE-A0 control lands at −14.81 ± 0.029, ~0.21 dex
from its injected amplitude — that is single-realisation scatter of ONE steep-spectrum GP draw
(the information sits in the few loudest modes; the quoted ± is the likelihood width, not the
draw scatter). The GATE object is the deterministic population Asimov, which has no such scatter.

**The drain becomes a real measurement under remedy A.** With σ(log10 A_bg) ≈ 0.02–0.03 dex per
refit and the drawn hierarchy carrying 19.7% / 44.8% / 77.6% of population power in its top
1 / 3 / 10 members: a top-1 resolution moves A_bg by 0.048 dex (~2σ), top-3 by 0.129 dex (~6σ),
top-10 by 0.324 dex (~16σ). The capstone's headline curve is measurable exactly when the loop
does anything at all.

Figure: `GLACIER_g2aii_forensic.png` (three panels: the coincident monotone profiles at both
incumbent venues; the T=30 quad/Occam decomposition; the remedy-A recovery). Lean bank:
`reports/glacier_g2aii_forensic.npz` (all three venues' profiles, decompositions, controls,
prefixed `t15_/t30_/bandA_`).

### V.2 What remedy A changes in the pre-registration (proposed, NOT enacted — Matt authorises)

1. `draw_population(band_log10f=(-8.7, -7.5))` for every campaign draw (the parameter is
   threaded and banked; the incumbent default is untouched, so every existing bank re-reads
   under its own convention). `BackgroundFit(band_log10f=...)` matched to the generative band.
2. g2a-ii re-gated at the extended band, T=30/n=256, tolerance 0.15 dex as pre-stated (measured
   headroom: 0.045 observed, σ≈0.02). The smoke rung stays plumbing-only (edge REPORT-ONLY at
   T=15/n=32 — §II shows the smoke venue is information-poor even band-extended... to be
   measured by the re-gate's own smoke print, not assumed).
3. The two banked FAIL fit-gate banks and all three forensic banks stay on disk untouched — the
   FAIL is evidence, not an accident (stems are band-tagged so the re-gate cannot overwrite them).
4. Population realism scope line (travels with every campaign number): the band is a modelling
   choice of the mock universe, declared; the e-mix, hierarchy width, and NG15 normalisation are
   unchanged (LOTTERY realistic point).

### V.3 DRIVER-GATE VERDICTS (G-d1..G-d3) — landed, with two findings of their own

Three runs of `gl_dgate.sbatch` (all hgx03 H200, batch_gpu, cpus=8 — the pinned floor lane):

**Run 1 — job `12706069`: NaN DATA, crashed.** First diagnosis was the missing merger guard
(embed_igniter's docstring promised the CHORUS `tie_member` cap and did not implement it; the
WEAVE tie shortens t_coal by 2^{8/3}·F(e) ≈ 173× at e=0.7, so an uncapped comb CAN inject
coalesced kernels — SPARK-3 §2's `jnp.power(negative, -3/8)` NaN pathology, in the data itself).
The guard was implemented per `chorus.tie_member` verbatim (TERM_FLOOR=0.2, chirp-mass cap,
`n_clip` banked per cell). **Run 2 then printed `n_clip 0` and `min t_coal/span = 11.109` for
this same census — no member was anywhere near the cap, so the guard did not fix the NaN; it
was never the cause.** The NaN is attributed to the EMBER-class shared-JAX-cache race: 12706069
ran concurrent with the forensic/band-A jobs on the same `HARBOR_JAXCACHE` dir. The guard STAYS
(the arithmetic requires it — an e=0.7 igniter at a tighter draw WILL cross it), but the honest
record is: guard = required by arithmetic, NaN = cache race. Consequence for the fan: per-job
cache dirs or serialized first-compile, already the standing EMBER mitigation.

**Run 2 — job `12706109` (still on H_ABSENT=-18 code): 6 of 7 gates PASS, G-d3 FAIL — and the
FAIL is a measurement.** Verdicts:

| gate | check | result |
|---|---|---|
| G-d1a | embed: 32→63 slots, chan 63, n_clip 0 | PASS |
| G-d1b | merger guard: min t_coal/span 11.109 > 0.8 | PASS |
| G-d1c | e=0 structural collapse (census unchanged) | PASS |
| drain | profile finite; ahat −15.600±inf, edge (REPORT-ONLY at this venue, per §II) | PASS |
| G-d2a | rung-mask structure: own-term live, certified at q, rest decohered | PASS |
| G-d2b | rung-0 columns, brightest fed: **live rows 100.0%** (need 100), dlnL finite, q∈[0,1]; **E-step wall 300.0 s** | PASS |
| G-d2c | rung-2 coherent path finite; max\|LNL shift\| 1.335e−04 (report-only, faint venue) | PASS |
| G-d3 | carried-absent ≡ smask-zero: \|d p0\| 9.127e−03, max\|d FtNmy\| 1.376e+04 vs 1e−6 | **FAIL** |

The G-d2b number is the one SPARK-2's global-pmask confound killed: 100% live rows at rung 0
with the per-target E-step, at a 300 s wall for the full column set — the scoreboard is wired
and affordable (the fan's per-iteration E-step cost is this, not the ~116× per-target rebuild).

**The G-d3 FAIL is the H_ABSENT=-18 remnant, measured.** Zeroing a carried source via
theta (h→1e−18) leaves a 3.5e−4-relative template remnant against true smask-zero absence —
invisible in single-source arithmetic, but G-d3 compares whole-basis projections (`FtNmy` norms
~4e7, so 1.4e4 absolute = exactly that relative remnant). This is a real (tiny) bias the
carried-via-theta design would have injected into every M-step. Fix: `H_ABSENT_GL = −30.0`
(strain 1e−30 is fp-exact absence at float64; verified the evaluator returns bit-zero template).
The -18 convention is spark3's own H_ABSENT — fine there (scores, ratios), not fine here
(basis-level identity gate). One-line diff, compile-checked.

**Run 3 — job `12707708` (clean run on the −30 code): === DRIVER GATES: PASS ===** (hgx03
H200, wall 474 s). G-d3 now: **\|d p0\| = 0.000e+00, max\|d FtNmy\| = 1.397e−08** (vs 1e−6) —
carried-via-theta at −30 IS smask-zero to fp. All other gates re-PASS (G-d2b 100% live rows,
E-step wall 163.4 s this run; G-d2c shift 7.9e−05). The banked gate record for the campaign is
run 3; runs 1–2 are kept as the forensic trail. The job's own final print stands as the state
of the holds: `PASS (fan HELD on FORGE-G2 runtime-SMASK + resume drill)`.

### V.4 DECISION — what Matt authorises before anything wide runs (nothing below fires by default)

The holds as they stand tonight, in dependency order:

1. **REMEDY A (the band): authorise or refuse.** §III/§V.1 is the case: the pre-registered
   g2a-ii cannot pass at ANY venue in the incumbent 10–32 nHz box (even a TRUE-A0 control pegs);
   extended to flo=−8.7 the same gate object lands 0.045 dex from truth with σ≈0.02, and the
   drain becomes a 2σ/6σ/16σ measurement at top-1/3/10. Mechanics are ready and band-tagged:
   `harbor_py hpc_harbor/glacier/glacier_pop.py fitgate --t 30 --n 256 --flo -8.7` re-gates
   g2a-ii under the pre-stated 0.15 dex tolerance. Without remedy A the campaign's headline
   observable has infinite error bars — refusing it means re-scoping to the ledger-only drain
   (§III arm c, a weaker question), not running the fan as pre-registered.
2. **FORGE-G2 runtime-SMASK: push.** Still NOT landed (`trackB_b1_core.py:372-390` bakes smask
   via `_params`; re-checked in code this session). The Stage-1 driver refuses `cell`/`null`
   modes in code until this lands — the M-step refit path pays the 17–32 GPU-hr recompile tax
   without it. Waiting on Matt's "pushed."
3. **After (2): the resume drill** (`gl_drill.sbatch`, built, DO NOT SUBMIT before (2)) — resume
   leg vs straight leg, column-exact `drillcmp`. Then the driver-gate job ONCE more on the
   post-push tree (same sbatch; the SMASK landing touches the evaluators G-d3 gates).
4. **After (1)+(2)+(3): the fan** — 2 igniters × 8 skies at T=30, 6 iterations, per §R budget
   (~75–85 GPU-hr, 150 STOP stands). Scrambled-null cells and full floors on the verdict subset
   as pre-registered. Nothing about the fan's shape changed this session except: band per
   remedy A, `n_clip` banked per cell, distance truth drawn per cell (Arm B), M-step carried
   handling via `theta_with_absent` at −30.

What does NOT need authorisation (already done, evidence-class): the forensic + remedy-A study
banks, the driver-gate record, the figure, this report.

### V.5 ADD-LIST (session 3; staged with `git add`, Matt commits — nothing committed by me)

New:
- `hpc_harbor/glacier/glacier_g2_forensic.py` — the F1–F4 forensic (H-BUG/H-FAINT/FIT-OK)
- `hpc_harbor/glacier/gl_forensic.sbatch`, `gl_frn_banda.sbatch`, `gl_dgate.sbatch` — the three
  gate/study jobs (H200 lane, cpus=8)
- `hpc_harbor/glacier/gl_drill.sbatch` — resume drill (HELD: post-FORGE-G2 only)
- `hpc_harbor/glacier/glacier_g2ii_forensic_figure.py` + `GLACIER_g2aii_forensic.png` — the
  3-panel deaf-venue figure
- `reports/glacier_g2aii_forensic.npz` — lean bank (t15_/t30_/bandA_ profiles, decompositions,
  controls)

Modified:
- `hpc_harbor/glacier/glacier_pop.py` — `band_log10f` threaded end-to-end (drawn, banked,
  band-tagged stems), `--flo/--fhi/--grid-lo/--grid-hi`, deaf-venue finding recorded in the
  gate docstring
- `hpc_harbor/glacier/glacier_loop.py` — CertScoreboard (G-d2 landed), merger guard, NoiseDrawer
  arity, force_recompute_lgwb at T=30, M-step `carried_now`/`theta_with_absent`, distance-truth
  draw, `H_ABSENT_GL=−30`, drill modes, gate mode reworked
- `GLACIER_capstone.md` — this addendum

Held on disk, NOT staged (working state, not evidence): `GLACIER_results/*` banks (the two FAIL
fit-gate banks, three forensic banks, band-A bank — per repo convention results dirs stay
untracked; the LEAN bank above is the tracked evidence).

### V.6 STOP

Stopping here, per the brief. Not run, by hold: the fan (`cell`/`null` — refused in code),
`gl_drill.sbatch` (post-FORGE-G2 only), the remedy-A re-gate (behind authorisation (1)).
Budget spent this session: ~1.9 GPU-hr (forensic 2-leg ~0.6, band-A study ~0.5, three
driver-gate runs ~0.8) against the 150 STOP — cumulative campaign spend remains under 5 GPU-hr.
The next session fires on Matt's readback of §V.4; items (1) and (2) are independent and can
land in either order.
