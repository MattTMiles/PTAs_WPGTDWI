# SUMMIT — the closure campaign: where does the self-calibrating loop CLOSE on a realistic sky?

*Agent SUMMIT, ACCRE. You stage, Matt commits; no git commit/push/tag from here. Operating
under the STANDING AUTONOMY CHARTER (2026-07-25): operational decisions proceed and are
posted here; STOPs are posted-and-parked, never idled on.*

**STATUS (refreshed by agent PHOENIX, 2026-07-29 — recovery session; the SUMMIT session
died 2026-07-25 with Phase-1 fans in flight):**

**Phase 0 harvest COMPLETE and RE-CUT (96 cells; the array ladder has fully drained and
both null arms have landed). Phase 1 COMPLETE on both dials — D2 (priors) written at
§1.3; D3 (N_psr) ran to completion 2026-07-26 and is harvested here at §1.4. Mission
unchanged: BOUNDARY-MAPPING centered on the feasible rung (r13p9, none, T=30, w=1).
THE §0.2a SIGHTING CALL IS NOW ANSWERED — the same-rung null baseline landed (§0.2c);
the control kills one of the two strong candidates and leaves the other standing.
Phase 2 cell list PROPOSED at §2.0, partly blocked by a new STOP.**

**STOPs on the board: #S1 (wrong certs at super-NG15 rungs — ruled by Matt at S4.18,
response ran, 20/20 killed). #S2 (Stage-1 control manufactured at i2 — ruled at S4.20,
frontier-v2 the response). >>> #S3 NEW AND UNREAD (§1.5): the D3 +100 control
manufactured at ITERATION 0, ratio 0.182 against a 0.5 bar, under the v1 frontier.
Posted and parked; the +100 rung's feed claim is quarantined pending Matt's ruling.**

---

## §0 — PHASE-0 HARVEST READBACK (the Stage-2a plane, from the GLACIER ladder banks)

### §0.0 What was on disk when SUMMIT arrived (2026-07-25 ~12:00)

| arm | job | state |
|---|---|---|
| Sky ladder signal (3 rungs × 2 arms × 4 skies × 6 iters) | 12763562 (A100) + 12763563 (H200) | **ALL 24 CELLS BANKED** |
| Sky-ladder scrambled nulls (12 cells @ r13p25) | 12763563 + rerun 12766445 | 5 banked (the no-promote survivors of the pre-amendment STOP); **7 rerunning** under §S4.15.1 semantics |
| Array ladder (5 (w,T) blocks × 2 arms × 4 skies + 12 nulls) | 12763564 (A100/dgx03) | **(w=1,T=40) e07 block banked (4 cells) + none s0/s1 at i0**; 46 tasks pending behind the dgx reservation |
| Stage-1 rerun signal | 12763565 | complete (context for §S4.15; not re-harvested here) |
| Stage-1 null rerun (9 tasks) | 12766444 | running (hgx03) |

No prior harvest existed; SUMMIT performed it (CPU only, zero GPU):
`hpc_harbor/summit/summit_harvest.py` → bank `reports/summit_s0_harvest.npz` (35 cells).
HOLD is scored exactly as pre-registered: ≥1 fed member and ≥2 **consecutive** refit
intervals in which every commonly-fed member moves ≤ 0.05 dex in both (log10_fgw,
log10_mc). Bite/regress = `a_bg − a_eff_drawn` at i0 / final iter.

### §0.1 THE PLANE — three contours over (sky loudness × capability), as banked so far

Loudness axis: r13p9 = NG15-feasible ceiling (tension 0σ), r13p5 (+7.9σ), r13p25 (+12.5σ).
Capability axis so far: (w=1,T=30) fully sampled; (w=1,T=40) e07-only; w=1/2, 1/4 pending.

**FIRST-BITE contour (N_resolved > 0, any iteration):**
- ON at r13p25 and r13p5 in **8/8 cells each** (both arms, every sky), feeding 1–23 slots.
- At the feasible rung r13p9 it is ON **only in the circular arm**: none 3/4 skies
  (s0,s1,s2; s3 no), **e07 0/4 skies at T=30 AND 0/4 at T=40**.
- **NEW DIAL READING — eccentricity KILLS first-bite at the feasible loudness.** The e0.7
  comb spreads the conditioned brightest member over 32 harmonic channels; min
  concentration ratio across all 287 slots is 0.981 (nothing near the 0.5 feed bar) vs
  0.266 (slot 0 crosses) in the paired circular cell. Structure assistance at cert onset
  (SURFACE) inverts to structure DILUTION at the feed frontier when the total is
  NG15-bounded. Igniter choice at the realistic point should be circular.

**M-STEP TRUST (HOLD) contour:**
- **The FEASIBLE rung HOLDS: r13p9/none/T30 skies 0 and 1.** Sky 0 is the strong form —
  the fed member converges to an exact M-step **fixed point**: Δ(fgw,mc) = (0.0000,
  0.0000) dex over i3→i5 (series +0.023/−0.044 → −0.008/−0.022 → −0.004/0.000 → 0/0 →
  0/0). Not wander-below-bar; a fixed point.
- Also HOLD: r13p5 e07 s3, r13p5 none s0/s1, r13p25 none s3.
- NOT holding: every bright-rung cell with multi-slot feeds (wander 0.10–0.23 dex/interval
  — the §S4.15 ACTIVE-REGRESSIVE wander, reproduced on the ladder).
- **The M-step trust boundary is INVERTED in loudness: the loop is trustworthy at the
  faint edge (1 fed member, sharp fixed point) and regressive at the bright rungs (many
  fed members, mutual confusion).** The "rung where the wander stops" (Matt's decision-2
  readout) is the BOTTOM of the ladder, not the top.

**CERTIFICATION contour:** N_cert > 0 only at r13p25 (s0: 3, s1: 1) and r13p5 (s0: 2,
none-s3: 1). **All but one of these cells carry wrong certs (STOP #S1 below); the single
clean-cert cell (r13p5/none/s3, N_cert=1 at i0 and i4, wrong=0) has a flickering cert at
a degenerate floor (zf 0.97–1.00, floor 0.0) — not attributable.** At the feasible rung
certification is ZERO everywhere, as §S4.14.1 predicted (conditioned brightest lands at
~2.3 nHz, fringe-starved). The cert contour does not reach the NG15-consistent region on
the (w=1) column.

**DRAIN:** first bites are large and real at every rung where feed fires (i0 bite −0.36
to −0.77 dex vs a_eff). End-state: bright rungs regress (+0.10..+0.48 above baseline, the
ACTIVE-REGRESSIVE signature), but the r13p9 HOLD cells regress once at i1 and then park
~0.03–0.30 dex BELOW or at baseline (s0: −14.91 bite → parks −14.51; s2 retains −0.31 at
i5). At the feasible rung the loop keeps part of its bite.

### §0.2 CANDIDATE CLOSURE EVENTS — posted for Matt's read, NOT declared (charter §5)

The sighting bar: iteration N+1 feeds/certifies a member iteration N could not,
attributable to iteration N's improvement, scrambled-clean, safety intact. The banks
contain four **clean** (wrong_cert = 0) late feeds of census members (promote_events is
(iteration, slot); slots < n_src are census members, ≥ n_src are igniter harmonics):

| cell | event | safety | scrambled arm |
|---|---|---|---|
| **r13p25/none/s0** | census slot 1 fed at **i1** (i0 could not) | wrong=0 | rung has nulls (7 rerunning) |
| **r13p25/none/s1** | census slot 8 fed at **i3** | wrong=0 | same |
| r13p5/e07/s2 | census slot 3 fed at i4 | wrong=0 | **no nulls at r13p5** |
| r13p5/none/s1 | census slot 3 fed at i4 | wrong=0 | **no nulls at r13p5** |

(Two further late census feeds exist in wrong-cert cells — r13p25/e07/s0, r13p5/e07/s0 —
and are quarantined with #S1.) Attribution analysis (did i0's drain/fed term tighten the
late member's conditional widths?) and the control-arm late-feed baseline (do the
scrambled cells late-feed at the same rate?) are OWED before any of these can face the
bar — the null reruns supply the baseline; Phase-2 adds scrambled cells at r13p5. These
are CLOSES-VIA-DRAIN-shaped; no CLOSES-VIA-CERT candidate exists anywhere yet.

### §0.3 >>> STOP #S1 (SAFETY-CLASS, posted-and-parked): WRONG CERTS AT THE SUPER-NG15 RUNGS

At r13p25/e07/s0 the scoreboard certifies 1–3 pulsars per iteration (11 iteration-events
total), **q_max = 1.0, dlnL 17–28 nat over a 15.3-nat floor, on_true = FALSE** — evidence
not anchored at the true parameters. Also r13p25/e07/s1 (2), r13p5/e07/s0 (3). All in the
**eccentric arm at tension ≥ 7.9σ**; zero wrong certs in any circular cell, any r13p9
cell, or any banked null. This is the IGNITE loop-wrong-cert-cascade class arriving in
GLACIER's loop at high loudness. **Parked:** the three cells' cert columns are quarantined
from every claim in this campaign (their feed/drain columns stand, labelled); nothing
about Phase 1/2 depends on them. **For Matt:** whether this becomes a campaign finding
("the loop's certification layer breaks first in the eccentric arm above +8σ" — arguably
the M-step wander feeding template-mismatched power to the scoreboard) or a machinery
audit is your call; the per-pulsar anatomy is one script away.

Instrument notes riding with the harvest (report-only): (i) the drain fit PEGS at its
grid ceiling (a_bg = −13.6, σ=inf) at r13p25/e07/s0 from i1 — the +0.25-dex grid cap
above a_eff is too tight one rung above NG15; Phase-2 cells at bright rungs should widen
`a_bg_grid`. (ii) zf = 1.00 → floor 0.0 (emp-q95-degenerate) in most faint-state cells —
correct per ANCHOR S4 but it makes the cert bar at the feasible rung entirely the lnK
term; stated so no one reads N_cert=0 there as a floor artefact. (iii) *(added on
reconciliation, from GLACIER S4.17)* every e07 `bite_i0`/`regress_end` column in the
harvest carries the a_eff SCOPE ARTIFACT — `a_eff_drawn` is computed before the comb
embedding redistributes the bright member's power across 31 harmonics (partly out of
band), so e07 drain offsets are NOT comparable to none-arm ones as printed.

### §0.3a RECONCILIATION — #S1 is GLACIER S4.17's STOP, and Matt has ALREADY RULED (S4.18)

A concurrent GLACIER session (same charter, §S4.16) independently found and posted the
identical cascade STOP (S4.17) hours before SUMMIT's harvest, and **Matt's S4.18
response is on the record**: the interpretation is BANKED as the capstone's safety
headline ("the loop's first certifications are false, they are structure-made, and the
criterion layer that catches everything else does not catch this composition" — the
comb adverse at BOTH ends: thins information below the feed bar when faint,
manufactures when bright). Diagnostics D1–D3 are authorized to GLACIER (≤8 GPU-hr; D1
mechanism cut queued as `12766749`). **Standing consequence for SUMMIT: a
rigidity-gate v2.3-candidate is pre-registered for Phase 2/3** (threshold defined from
SPARK's banked selectivity statistics before scoring, never tuned on the cascade cells,
never retro-applied) — any certification claim in a Phase-2 composed cell with
eccentric structure WAITS on that gate landing green; circular-arm cells are unaffected
(the cascade is e07-only, none-arm clean at every rung). #S1 therefore converts from
"parked awaiting ruling" to "ruled; response running in the GLACIER lane."

CONCURRENCY LEDGER (the g2.4 discipline, kept explicit): GLACIER owns
`glacier_loop.py` (its 12:01 staged edit = the S4.15.1 amended null semantics + G-d4c
rework — verified signature-compatible with SUMMIT's wrappers; signal-cell path
untouched), the cascade diagnostics, the array ladder, and the null reruns. SUMMIT
owns `hpc_harbor/summit/*`, `SUMMIT_results/`, and this report. Shared banks are
read-only to SUMMIT. Neither edits the other's files; both queues share fair-share
only. The array-ladder driver-version split (tasks 0-5 pre-12:01, 6+ post) is BENIGN:
the diff touches only scrambled-arm STOP logic (all array-ladder nulls are pending →
uniform semantics) and `mode gate`.

### §0.4 DECISION RULE — fired as pre-registered

**HOLD exists on the ladder (strongest at the feasible rung itself) → SUMMIT's mission is
BOUNDARY-MAPPING, Phase 2 centered on (r13p9, none, T=30, w=1).** Operational reading
(charter: rung selection under a stated rule): the D2/D3 single-dial ladders are pulled
INTO the boundary-mapping mission at the holding rung — they are the two unmeasured axes
through the center cell, and Phase-2 composition requires their boundaries. D1/D4
boundaries INHERIT from the array ladder as it drains (nothing re-run).

### §0.5 THE DSA ANSWER + control baseline — OWED, watcher armed

The DSA question (does rms/4 bring the loudest NG15-consistent sky across first-bite at
T=30–40?) needs the w0p25 blocks of 12763564 — not yet run (46 tasks pending behind the
dgx reservation; ETA unknowable from queue state, watcher posts on drain). Note the w=1
harvest ALREADY gives a partial answer: first-bite is ON at r13p9 in the circular arm at
TODAY'S rms — the DSA rungs will measure margin, not existence. The control-arm
first-bite baseline (12766444/45) completes the §0.2 table when it lands.

### §0.6 Ledger

SUMMIT GPU spend: **0.0 GPU-hr** (harvest is CPU). Ladder spend is GLACIER's, excluded
per the brief. Hard STOP: 150 GPU-hr cumulative for SUMMIT's own submits.

### §0.2a ATTRIBUTION ANATOMY (CPU, from the banks) — two strong candidates, two weak

Concentration-ratio trajectory of each late-fed slot (feed bar 0.5):

| cell / slot | ratio trajectory to feed | reading |
|---|---|---|
| **r13p25/none/s1, slot 8** | 1.000 → 0.622 → 0.563 → **0.480** (fed i3) | **monotone three-iteration descent — improvement-attributable in shape**; the frontier tightens every refit until the member crosses |
| **r13p25/none/s0, slot 1** | 0.626 → **0.460** (fed i1) | a real 0.17 one-step improvement after i0's feed+drain — attributable |
| r13p5/e07/s2, slot 3 | 1.0, 1.0, 0.571, 1.0, **0.496** (fed i4) | flicker across the bar — threshold noise class |
| r13p5/none/s1, slot 3 | 0.538, 0.546, 0.524, 0.538, **0.491** (fed i4) | four-iteration hover at the bar — the bar sits inside the jitter |

Both STRONG candidates live at r13p25 in the circular arm — the one rung that HAS a
scrambled arm. The completed baseline test is: do the rerunning null cells (12766445)
late-feed their real members at the same rate/shape? A matching rate kills the
attribution; a null arm that only ever i0-feeds (or never feeds) leaves the two
candidates standing. Post-feed both cells REGRESS in drain (a_bg −14.6 → −13.7/−13.8,
the bright-rung M-step signature), so these are candidate CLOSES-VIA-DRAIN events whose
drain benefit is subsequently eaten — if they survive the null baseline, the loop
CLOSES at +12.5σ but does not COMPOUND there. **Held for Matt with the sighting call.**

**First control baseline landed (Stage-1 null rerun, 12 cells banked, amended
semantics):** control cells feed real members at i0 in 8/12 (the live frontier working
in both arms, as S4.15.1 predicted) and LATE-feed in exactly one (gl1_null0/e07/s0,
slot 7 at i5) — and that one is a **single-step discontinuous crossing** (ratio 1.000
flat for five iterations, then 0.475 in one step: frontier jitter, no build-up). The
control arm so far contains NO analog of the monotone multi-refit descent that marks
the strong candidate (1.000 → 0.622 → 0.563 → 0.480). Shape taxonomy as measured:
single-step crossings and bar-hover occur in both arms (jitter class); gradual descents
have appeared only in the signal arm. The same-rung (r13p25) null baseline from
12766445 remains the decisive comparison. Scope: Stage-1 nulls sit at the realistic
normalisation, one loudness class below the candidates' rung.

### §0.3b >>> STOP #S2 (SAFETY-CLASS, posted-and-parked): FIRST TRUE MANUFACTURING EVENT — the control promoted the SCRAMBLED member

2026-07-25 14:09, Stage-1 null rerun task `12766444_15`: cell `gl1_null0_none_s1`
raised the **surviving** (post-amendment) manufacturing STOP — **NULL-ARM
SCRAMBLED-MEMBER PROMOTE at iter 2**: the control arm fed slot 0 *from its wrong-sky
template*. This is manufacturing proper under §S4.15.1's narrowed semantics — the STOP
class that was deliberately kept. Trail: ratio[0] = 1.000 (i0), 1.000 (i1), crossed
the 0.5 bar at i2; meanwhile the real member (slot 1) was fed at i0 and the cell's
drain was RISING (−14.63 → −14.50) — the crossing rode the fed-member wander + drain
refit, the same (structure × wander) composition as the cascade, here in the
manufacturing direction. The task's packed partner cell completed cleanly through i5.

**Consequences for the closure evidence (why this matters to SUMMIT):** the frontier
has now been *demonstrated* to cross the feed bar on a structurally wrong template at
realistic SNR. Single-step bar crossings are therefore NOT evidence of real structure
on their own — the null just produced one onto a wrong-sky template. Of the §0.2
taxonomy, only the **monotone multi-refit descent** class remains as candidate closure
evidence, and the same-rung scrambled baseline (12766445, queued) is now load-bearing
for the sighting call. Nothing new of SUMMIT's is parked (the sighting call already
was); the bar for it just measurably rose.

**Instrument defect, report-only (GLACIER's file, not fixed per house rule + the
no-edit constraint):** the STOP message promises anatomy "banked at …_i2__…npz" but
**no i2 bank exists** — `CampaignStop` fires at promote time (step b), the bank writes
at iteration end, so the manufacturing STOP loses its own anatomy. The i2 state is
recoverable deterministically from the i1 checkpoint (checkpoints are continuations,
per the drill). Flagged for GLACIER/Matt.

### §0.2c THE SIGHTING CALL — the same-rung null baseline LANDED; one candidate dies, one stands

*(PHOENIX, 2026-07-29. Instrument: `hpc_harbor/summit/smt_sighting_baseline.py`, CPU,
zero GPU. §0.2a named the decisive test; `12766445` has since banked the r13p25 none-arm
control 6/6 across null0/null1/null2 × s0/s1, so the test is now runnable and was run.)*

The test as SUMMIT pre-registered it, verbatim: *"do the rerunning null cells (12766445)
late-feed their real members at the same rate/shape? A matching rate kills the
attribution; a null arm that only ever i0-feeds (or never feeds) leaves the two
candidates standing."*

**RESULT — the control never late-feeds at all.** Across the 6 same-rung control cells:
**0 late census feeds**, of any shape. Across an independently re-scored Stage-1 control
(20 cells, at the realistic normalisation): **2 late census feeds, both SINGLE-step**.
Combined control exposure: **26 cells, 2 late census feeds, ZERO of the GRADUAL class.**

| arm | cell | slot | fed@ | class | trajectory to the 0.5 bar |
|---|---|---|---|---|---|
| SIGNAL | r13p25/none/s1 | 8 | i3 | **GRADUAL** | 1.000 → 0.622 → 0.563 → **0.480** |
| SIGNAL | r13p25/none/s0 | 1 | i1 | SINGLE | 0.626 → **0.460** |
| CONTROL | — (6 cells) | — | — | — | *no late feeds* |

**THE CALL, in two parts:**

1. **Candidate 2 (r13p25/none/s0, slot 1) is DEMOTED — the control already killed this
   shape.** §0.2a read its one-step 0.17 improvement as "attributable". STOP #S2 then
   demonstrated a single-step bar crossing onto a *wrong-sky* template, and the Stage-1
   control has now produced two more single-step crossings. A SINGLE-step crossing is a
   measured control behaviour. It is not evidence. This candidate comes off the board.
2. **Candidate 1 (r13p25/none/s1, slot 8) STANDS.** The monotone three-refit descent has
   no analog anywhere in 26 control cells. The control is not merely quiet — in Stage-1
   it feeds at i0 in 16/20 cells, so it has the *exposure* to produce late feeds and
   produces only the jitter class.

**Honest scope, stated because it bounds the claim:**
- **n = 1.** GRADUAL occurs exactly once in the entire campaign. A class defined by one
  member is a hypothesis, not a rate.
- **Same-rung control exposure is weak on its own** — only 1 of 6 r13p25 control cells
  i0-feeds at all, so its zero is partly an exposure effect. The load-bearing leg is the
  Stage-1 control (16/20 cells feeding, still 0 GRADUAL), which is one loudness class
  below the candidate.
- **Cross-lane:** r13p25 signal cells are dgx03/A100, the null reruns hgx03/H200.
  `conc_ratio` cross-lane wobble is 1.66e-02 (SG-D2a) — 14× under the feed-bar margin and
  1–2 orders under the descents classified, so the SHAPE taxonomy is readable cross-lane.
  A *rate* comparison at the bar would not be, and none is made.
- The event is at **+12.5σ**, three rungs above anything NG15-consistent, and it is a
  CLOSES-VIA-DRAIN candidate whose drain benefit is subsequently eaten (§0.2a).

**>>> DECLARING A SIGHTING IS CHARTER STOP #5. NOT DECLARED. This is the baseline
result Matt asked to read on return; the call is his.**

### §0.2d Correction to §0.2a's control count

§0.2a reported the Stage-1 control late-feeding "in exactly one" cell. The independent
re-score finds **two**: `gl1_null0/e07/s0` slot 7 at i5 (1.000×5 → 0.475) and
`gl1_null2/none/s0` slot 2 at i1 (0.504 → 0.460). Both SINGLE. The conclusion §0.2a drew
is unchanged and slightly strengthened — the control's late-feed repertoire is jitter-only
across a larger sample than was quoted.

---

## §1 — PHASE 1 (as boundary-mapping instruments): the D2/D3 dial builds

Both dials run at the boundary-mapping center (r13p9, circular arm, T=30, w=1) — 
the D2/D3 axes through the HOLD cell. Builds are wrappers over the STOCK driver
(`hpc_harbor/summit/summit_loop.py`, `summit_d3.py`): `glacier_loop.py` is NOT edited —
the 46 pending array-ladder tasks import it at start time, and editing it mid-array is
the g1 hazard. Runtime patches only (chorus-C7/COMPASS pattern), banks in
`SUMMIT_results/` (never GLACIER's dir), evidence gates COPIED with a copy-record line.

**D2 (priors).** Tiers: lit (inherited from the banked ladder) / vlbi = min(lit, 1 pc)
on the union-18 (ignite.py IG4 verbatim; truth draw follows the tier) / subpc =
min(lit, 0.1 pc) same rule (DECLARED, the extension rung). Gates: SG-D2b/c (CPU: tier
columns off-union-identical + monotone; tiered truth-draw stream alignment) + SG-D2a
(GPU: wrapper inertness — the stock lit path through the wrapper, i0, column-compared
to the banked ladder cell). Gate job `12766491`. Fan (11 cells: 2 new tiers × 4 skies +
3 scrambled nulls at subpc, the manufacture-riskiest tier): 
**header estimate ~5.5 GPU-hr (11 × ~0.5 h, the measured feasible-rung H200 cell wall);
wallclock ~1.7 h at 4-wide, ~6 h at 1-wide backfill.** The gate job's `gatea` leg warms
the exact campaign compile shapes, so the fan shares one cache race-free.

**D3 (N_psr).** Rungs: 116 (inherited) / +30 / +100. Additions EXTEND the real array
(real 116 untouched): clones of the top-N pulsars by the COMPASS quality rank at
seeded golden-angle uniform positions (SEED_D3 = 9300), THE BUNDLE INTACT — TOAs,
design matrix, EM distance, K-table prior row, frozen reference RN all the donor's;
only name+position differ (COMPASS rank-map convention for additions; its S1 shuffle
prices the pairing choice at 0.14 sky-σ — cited, not re-measured). Name-keyed lookups
strip the clone suffix. Per-rung array-keyed L_gwb banks
(`smt_L_gwb_T30_ext{N}_s9300.npz`, cpus=8, machine-local — the COMPASS venue-bank
hazard: shape alone does not key an array). Per-rung OUT subdirs = the array key on
disk. Gates: SG-D3b TWIN (a +1 clone at its donor's own position must reproduce the
donor's scoreboard row — every name-keyed inheritance verified in one shot), SG-D3c
position audit, SG-D3e bank assert. Prep+gate array `12766511` (bank cuts + batteries).
Fan (11 cells: 2 rungs × 4 skies + 3 nulls at +100): **header estimate ~8 GPU-hr
(cost scales ~npsr: ×1.26 at +30, ×1.86 at +100 on the 0.5-h base); wallclock ~2.5 h
at 4-wide AFTER a serialized first-compile pass per rung (tasks 0 and 4 run first —
the EMBER cache-race mitigation), ~9 h at 1-wide.**

Fans fire on gate green (charter: phase transitions are operational). Both fans +
both gate jobs ≈ **15 GPU-hr against the 150 STOP; cumulative SUMMIT spend after
Phase 1 ≈ 15.5 GPU-hr.**

### §1.1 D2 GATES: ALL PASS (with one instrument-calibration amendment, trail below) — FAN LAUNCHED

Job `12766491` (hgx03 H200): **SG-D2b PASS** (tier columns off-union-identical,
monotone, caps respected; the vlbi dial moves 16 of 18 union pulsars — 2 already sit
at/below 1 pc under lit). **SG-D2c PASS** (off-union truth draws bit-identical across
tiers; union |n_true| tightens 14.2 → 2.4 → 0.2 lit→vlbi→subpc — the tier dial works
end-to-end through the truth draw). **SG-D2a PASS after an amendment to MY OWN compare
bar** (not a pre-registered campaign bar — a gate authored and mis-calibrated by SUMMIT
this session): the stock-lit wrapper reproduction matched every deterministic column
EXACTLY cross-lane (fed_mask ==, n_res ==, |Δa_bg| = 0.0, |Δfloor| = 0.0, zf ==) but
showed max|Δconc_ratio| = 1.66e-02 against a 1e-6 tolerance. conc_ratio is
Fisher-derived — second differences at ~0.2-nat curvature — where cross-arch
(A100-banked vs H200-gate) reduction noise amplifies to ~1e-2; 1e-6 was never
achievable for that column cross-lane. Amended bar: **decision invariance** (zero
slots cross the 0.5 feed bar between banks; measured margin 0.233 vs wobble 0.0166,
a 14× separation) + a 0.1 sanity cap; all non-Fisher columns stay exact. Same family
as SPARK-3 §4.3's "arithmetic, not a defect" scope line. **D2 FAN LAUNCHED:
`12768821` --array=0-10%4** (8 signal + 3 subpc nulls; warm cache from the gatea leg).

D3 iteration trail (instrument-class, ~15 min GPU total): (1) missing `CW_lnL_check`
sys.path → 3-s import failure, fixed; (2) twin gate scored the all-carried (empty)
template — non-finite by construction, never run in production — fixed by feeding
slot 0 first (G-d2b convention); (3) **a real venue property found by the gate: an
EXACT positional twin makes the HD-ORF GWB block singular** (duplicate sky position →
Phi_gwb rank-deficient → every E-step row non-finite; CPU-container forensic confirmed
theta and data fully finite, isolating the GP block). Fix: the twin sits 0.1° off the
donor (inheritance breaks are O(1–100) nat, 0.1° of geometry is O(1e-3) — tolerances
recalibrated to <0.5 nat), plus a production guard refusing any extension whose min
pairwise separation < 0.05° (chance near-coincidence would hit the same singularity).
En route both **array-keyed L_gwb banks cut clean**: ext+30 `fp=5c8750ef9de488d9`
(6716², recon 3.3e-26), ext+100 `fp=585df47ec7d550ca` (9936², recon 4.3e-26), both
strict-loading; SG-D3b structural half PASSED twice (clone L0/prior == donor exact).
Gate battery: `12768840_0/1`.

### §1.2 D3 GATES: ALL PASS ×2 — FAN LAUNCHED

`12768840_0/1` (hgx03 H200, 295–299 s each): **SG-D3b twin** — L0/prior inheritance
EXACT, |Δ dlnL| = 5.9e-04 nat and |Δ q| = 5.0e-03 at the 0.1° twin offset against an
O(1–100)-nat break scale, lnK bit-exact; **SG-D3c** — +30 additions NN sep 32.7°/35.2°
min/med (dipole 0.005), +100 additions 17.8°/19.4° (dipole 0.000); **SG-D3e** — both
array-keyed banks strict-load at their extended shapes. **D3 FAN LAUNCHED:** warm pass
`12769018` (tasks 0, 4 — one first-compile per rung, the EMBER cache-race mitigation),
main pass `12769019` (afterany, %4). 8 signal (2 rungs × 4 skies) + 3 nulls at +100.

### §1.3 D2 COMPLETE — THE PRIORS DIAL MOVES TRUST, NOT REACH, AND SATURATES AT 1 pc

All 11 cells banked (fan `12768821`, zero failures; 8 signal + 3 subpc nulls, all
clean: zero feeds, zero certs, zero wrong certs in the scrambled arm; no late feeds
anywhere in D2). The dial readout at the feasible rung, circular arm:

| axis | lit (banked ladder) | vlbi-1pc | subpc-0.1pc |
|---|---|---|---|
| FEED (skies fed / 4) | 3 (s0,s1,s2) | 3 (same skies) | 3 (same skies) |
| HOLD (of fed) | 2/3 | **3/3** | **3/3** (s0 = exact fixed point from i0, wander 0.000 every interval) |
| drain retained at i5 (s1, s2) | −0.19, −0.31 | −0.21, −0.31 | −0.23, −0.25 |
| CERT | 0 | 0 | 0 |

> **D2 VERDICT (single-dial hold boundary):** the distance-prior dial does NOT move
> the feed boundary (the same 3/4 skies feed at every tier — feed is source-geometry-
> gated) and does NOT un-starve certification (0 at every tier — cert onset is
> f-position, which priors cannot touch). What it moves is **M-step trust**: HOLD goes
> 2/3 → 3/3 fed skies from lit → vlbi, with wander collapsing toward exact fixed
> points, and the loop retains −0.2..−0.3 dex of drain at i5 where lit regressed. **The
> gain SATURATES at 1 pc: subpc-0.1pc is statistically indistinguishable from vlbi-1pc
> on every column.** The D2 hold boundary therefore sits between lit and vlbi-1pc, and
> the useful engineering statement is: *one parsec of VLBI-class distance information
> on the union-18 converts the feasible rung's HOLD from a 2-of-3-sky property to a
> universal one — and more precision buys nothing further at this rung.*

### §1.4 D3 COMPLETE — THE ARRAY-SIZE DIAL MOVES REACH, NOT TRUST, AND NOT CERTIFICATION

*(PHOENIX, 2026-07-29. The D3 fan `12769215` (warm) + `12769216` (main) ran to completion
2026-07-25/26; the session died before harvesting it. Harvester
`hpc_harbor/summit/smt_d3_harvest.py`, bank `reports/smt_d3_ladder.npz`. CPU, zero GPU.
11 of 12 fan tasks banked 6/6; task `_10` STOPPED — see §1.5.)*

| axis | 116 (inherited) | 146 (+30) | 216 (+100) |
|---|---|---|---|
| FEED (skies fed / 4) | 3 (s0,s1,s2) | 3 (same skies) | **4 — sky 3 crosses** |
| HOLD (of 4) | 2 | **3** | 2 |
| CERT | 0 | 0 | 0 |
| wrong certs | 0 | 0 | 0 |
| mean bite at i0 | −0.302 | −0.303 | **−0.492** |
| mean regress at i5 | −0.056 | −0.046 | −0.047 |
| mean max wander | 0.081 | 0.111 | 0.117 |

> **D3 VERDICT (single-dial reach boundary):** the array-size dial **moves the FEED
> boundary** — the one sky that never feeds at 116 or 146 (sky 3) crosses first-bite at
> +100, and the mean first bite deepens 0.30 → 0.49 dex. It does **not** un-starve
> certification (0 at every rung, as with priors — cert onset is f-position, and neither
> knowledge nor array size touches it). Its effect on **M-step trust is non-monotone and
> small**: HOLD goes 2 → 3 → 2 while mean wander *rises* 0.081 → 0.117, i.e. more
> pulsars buy reach and pay a little trust back. **The useful engineering statement:
> +100 pulsars converts the feasible rung's feed frontier from a 3-of-4-sky property to
> a universal one, buys ~0.19 dex of extra first bite, and buys nothing at the
> certification layer.**

**THE CONVERGENT READING (this is the campaign-level finding, and it is new):** *sky 3 —
the one sky that will not feed at baseline — is turned ON by **two independent levers**:
lower white noise (rms/2 and rms/4 at T=30, array ladder, §0.1) and more pulsars (+100,
here). Both levers move FEED. **Neither lever, alone or anywhere on the measured plane,
moves CERT off zero at the feasible rung.** Reach is purchasable with sensitivity or
aperture; certification at NG15-consistent loudness is not — it is structure/venue-gated,
exactly as §S4.14.1 predicted and as the Stage-2a plane now shows in one panel.*

**>>> DRAIN RE-REFERENCE CAVEAT (LEDGER B1, 2026-07-29 — concurrent audit, read after this
table was cut).** LEDGER measures the feed-free offset at i0 (n_fed = 0, n = 51) at
**−0.389 ± 0.374 dex**, so `bite = a_bg − a_eff` referenced to zero is the wrong
reference, *and the error is two-sided*: for the **none arm the zero-fed baseline is
POSITIVE (+0.02..+0.24)**, so none-arm bites are **UNDERSTATED**. Consequences for this
table, stated precisely:

- The **bite column above is understated**, not overstated. LEDGER's proposed §S4.17
  amendment moves the none-arm range −0.30..−0.67 → −0.81..−0.96 (awaiting Matt).
- The **rung-to-rung comparison survives** — all three rungs are the same arm at the same
  loudness, so they share one reference and the *difference* (116 → +100 deepening by
  0.19 dex) is unaffected by a common offset. The absolute numbers are not quotable until
  the amendment lands.
- **It explains sky 3.** Sky 3's `bite_i0 = +0.235` (116) and `+0.197` (146) sit exactly
  in LEDGER's measured positive zero-fed baseline band (+0.02..+0.24) — because sky 3
  **feeds nothing** at those rungs. That is not a regression; it is the no-feed baseline.
  When sky 3 does feed (+100 here, or rms/2–rms/4 on the array ladder) it drops to
  −0.61 / −0.74. The convergent reading above is strengthened, not weakened.
- **The e07 arm's drain columns are the ones actually in danger**, not these: LEDGER finds
  the eccentric arm's biggest "bites" occur in cells that fed NOTHING (r13p9 e07 max −0.79
  at n_fed = 0) — and every r13p9 e07 cell on the Stage-2a plane is FB = 0 with bites
  −0.36..−0.97. Those are winner's-curse artifacts and must not be read as drain. **No
  e07 drain number is quoted in §1.4, §0.1, or on the Stage-2a figure** — the figure
  carries feed/hold/certify/manufacture only, no drain layer, so it is unaffected.

**SCOPE LINE, load-bearing:** every D3 bank was cut on **2026-07-25/26 under the v1
frontier** (frontier-v2 was wired at commit `ed78a1c`, 2026-07-26T23:34, ~4 h after the
last D3 task finished). The whole of §1.4 therefore carries the v1 scope line of
S4.19/S4.20, and the +100 rung additionally carries §1.5.

### §1.5 >>> STOP #S3 (SAFETY-CLASS, charter #3), POSTED AND PARKED: THE D3 +100 CONTROL MANUFACTURED AT ITERATION ZERO

Task `12769216_10`, cell `gl2_r13p9_null0_none_s1` at the **+100 rung**, banked
`SUMMIT_results/d3_ext100/..._STOPANAT_i0__hgx03_NVIDIAH200.npz`:

```
stop_kind  = scrambled_promote
iter       = 0
promoted_k = 0
conc_ratio[0] = 0.1822        (feed bar 0.5; every other slot 1.000)
```

Driver message: *"NULL-ARM SCRAMBLED-MEMBER PROMOTE at iter 0: the control resolved the
sky-scrambled igniter — manufacturing."*

**Why this is worse than STOP #S2, and why it is nevertheless not a new criterion
failure:**

- **Worse:** #S2 crossed the bar at **i2**, from 1.000 → 1.000 → cross, and SUMMIT read it
  as riding the fed-member wander plus the drain refit. **This one is at i0** — the
  honest starting state, before a single M-step has run. There is no wander for it to
  ride. And it is not marginal: **0.182 against a 0.5 bar**, a factor of 2.7 inside.
  At 216 pulsars the v1 frontier mistook a *sky-scrambled* igniter for a resolvable
  member on first look, with a large margin.
- **Not a new failure:** the run predates frontier-v2 by ~4 hours. This is a **v1-era
  event of exactly the class Matt already ruled on** at S4.20, whose pre-registered
  response (v2's data-support term) was validated at S4.20.1 against — notably — the
  *same scramble family* (`null0`, sky 1) at Stage-1, where it refused the wrong-sky copy
  at −6.5 nat with 0 false refusals.
- **Suggestive, not established:** the same (null0, s1) scramble manufactured at Stage-1
  (116 psr, i2) and again here (216 psr, i0) — earlier and deeper as the array grows.
  If that is real it says array growth *widens* the v1 frontier's blind spot rather than
  narrowing it. One pair of events is not a trend; stated as a hypothesis only.

**PARKED, and what is quarantined:** the **+100 rung's FEED claim in §1.4** (the sky-3
crossing and the −0.49 bite) is quarantined from every closure claim and from the
Stage-2a figure until ruled. The 116 and 146 rungs are untouched (their controls are
clean). Nothing in §0, §1.3 or the plane depends on +100.

**THE UNPARKER, costed and NOT run (charter #3 — Matt's call):** re-score this banked i0
state under the frozen frontier-v2 criterion and ask whether v2's data-support term
refuses the promote. Direct scoring from the existing bank, no replay, **~0.3 GPU-hr**.
Two outcomes, both decisive: *v2 refuses* → #S3 is a closed v1-era artifact, §1.4's +100
row unparks under the v2 scope line, and the D3 ladder needs **no** re-running; *v2
promotes* → v2 has a hole that S4.20.1's validation missed, which is an immediate
escalation and would put ~8 GPU-hr of D3 re-runs on the table. I have not fired it: it is
a diagnostic attached to a safety-class STOP, and the charter parks those for your read.
One word and it goes.

### §1.6 FROZEN-vs-LIVE M-STEP — COMPLETE. FREEZING DOES NOT RESCUE THE DRAIN, AND THE REASON IS THE *FIRST* FIT

*(PHOENIX, 2026-07-29. Discharges GLACIER **S4.15.1 item 2**, whose gate condition — "IF the
ladder shows a rung where frozen templates would have held the drain" — was met once the
ladder readback completed. Wrapper `hpc_harbor/frozen/frozen_mstep.py` (stock driver,
runtime-patched, per-slot freeze-after-first-fit); gate **SG-F1 PASS bit-exact**
(job 12834518); fan `12834892` 4/4 COMPLETED, 24/24 banked; readback
`hpc_harbor/frozen/frozen_compare.py`, bank `reports/frozen_vs_live.npz`. ~2.6 GPU-hr.
Live arm INHERITED, never re-run. Same lane (dgx03 A100-80GB) as the live banks.)*

| sky | give-back LIVE | give-back FROZEN | RESCUE | wander LIVE | wander FROZEN |
|---|---|---|---|---|---|
| 0 | +0.400 | +0.371 | +0.029 | 0.044 | 0.000 |
| 1 | +0.320 | +0.273 | +0.047 | 0.061 | 0.000 |
| 2 | +0.266 | +0.343 | **−0.077** | 0.139 | 0.000 |
| 3 | 0.000 (never feeds) | 0.000 | 0.000 | — | — |

**mean give-back LIVE +0.247, FROZEN +0.247, mean RESCUE −0.000 dex.**

> **VERDICT: FREEZING DOES NOT RESCUE THE DRAIN.** |rescue| is 15% of the give-back, both
> signs, signed mean **−0.000 dex**. Freezing removes template wander completely
> (0.044–0.139 → 0.000) and recovers **none** of the drain.

**WHY — and this is the useful part.** The give-back is **complete at the first post-feed
refit (i0 → i1)**: 92.8% / 85.3% / 128.8% of the total in skies 0/1/2. And at i0 and i1 the
two arms are **identical to +0.000 dex by construction** — freeze-after-*first*-fit means
nothing differs until i2. So the M-step's *later* wander cannot be the cause of
ACTIVE-REGRESSIVE; it is second-order jitter that freezing cleans up for free.

Reading the give-back against the driver's step order (`glacier_loop.run_cell`:
(b) promote → (c) drain refit → (d)(e) E-step + M-step) pins the mechanism:

> **The i0 bite is measured against a TRUTH-ANCHORED template** — the promote writes the
> drawn truth (`led.promote(...)  # AT DRAWN TRUTH`) and the drain is refit at (c) *before*
> the M-step runs at (e). **The first M-step fit then moves that template off truth to the
> data's peak, and the very next drain refit prices the move.** That single step is the
> whole regression. Freezing after the first fit cannot recover it, because **the first fit
> is what costs.**

**Consequences, stated carefully:**

- **The M-step is EXONERATED as a wandering process** — S4.15's ACTIVE-REGRESSIVE verdict
  should not be read as "the loop degrades itself by drifting." It takes its loss once, at
  the first fit, and thereafter (live) jitters ±0.03–0.08 dex with no trend.
- **The i0 bite is optimistic by construction** — it is a truth-anchored number that no
  self-calibrating loop can hold, because the loop does not know the truth. This is an
  independent route to the same caution LEDGER's B1 reaches from the reference side, and
  the two should be reconciled before any drain number enters the paper.
- **HOLD IS DEGENERATE IN THE FROZEN ARM — do not read `holdF`.** A frozen template has
  wander 0.000 by construction, so HOLD is trivially TRUE for every fed frozen cell (this
  is why sky 2 shows holdL=0 → holdF=1; it is arithmetic, not physics). The M-step-trust
  contour is meaningful in the **live arm alone**, and the Stage-2a figure uses live only.
- **Safety clean:** frozen arm 0 wrong certs, 0 certs; live 0 certs. No STOP fired.
- **Internal control passed:** sky 3 feeds in neither arm and its paired difference is
  **exactly 0.000 at every iteration** (a_bg flat at −14.397) — the wrapper is provably
  inert where it has nothing to freeze, independently of SG-F1.
- **Immune to LEDGER B1:** this is a paired contrast (same census seed, noise seed, venue,
  driver; only the write-back differs), so the feed-free reference offset cancels exactly in
  the difference. The verdict uses only paired differences.

**>>> THIRD ARM (running-belief M-step / FORGE-B): its pre-registered condition is NOT MET.**
Matt's rider fans the belief arm "if frozen rescues the drain, and FORGE-B has landed."
Frozen does not rescue the drain, so **PHOENIX does not fan it.** *But the mechanism above
argues the belief arm is now MORE interesting, not less:* the loss is entirely in the
**first fit's bias**, which is exactly what a belief-propagating M-step
(`Cov = inv(Fs+Pi)`, LEDGER A1's sigma-point route) would change — whereas freezing, which
only suppresses later motion, provably cannot. A one-line re-authorisation would point the
third arm at the first fit rather than at the wander. **Matt's call; nothing fired.**

---

## §2 — PHASE 2: COMPOSED CELLS (proposed, per charter — fan launch follows the posting)

**The composition question the single dials have set up.** Three dials have now been
measured through the boundary-mapping centre, and they separate cleanly:

| dial | moves FEED? | moves HOLD (trust)? | moves CERT? |
|---|---|---|---|
| D2 priors (lit → vlbi-1pc) | no | **yes** (2/3 → 3/3, saturates at 1 pc) | no |
| D3 array size (116 → +100) | **yes** (3/4 → 4/4) | non-monotone, ~none | no |
| Array ladder rms (→ rms/4) | **yes** (3/4 → 4/4 at T=30) | no clear move | no |

**Every dial alone leaves CERT at exactly zero at the feasible rung.** Composition is
therefore the only untested route to certification at NG15-consistent loudness, and that
— not more single dials — is what §2 is for.

**PROPOSED CELL LIST** (r13p9, circular arm, T=30, the boundary centre; 4 skies each):

| # | cell | dials composed | asks | cost |
|---|---|---|---|---|
| **C3** | `vlbi × w0p5` | trust × sensitivity | does 1 pc of distance info hold the drain while rms/2 widens the feed set? | 4 × ~0.5 h ≈ **2 GPU-hr** |
| **C1** | `w0p5 × +100` | the two feed-movers | do two reach levers compose, or do they redundantly turn on the same sky 3? | 4 × ~0.93 h ≈ **3.7 GPU-hr** |
| **C2** | `w0p25 × +100 × vlbi` | all three | the maximum-capability corner — the only cell on the map where CERT could plausibly leave zero | 4 × ~0.93 h ≈ **3.7 GPU-hr** |
| **N** | 3 scrambled nulls at C2 | — | manufacture control at the riskiest composed cell | 3 × ~0.93 h ≈ **2.8 GPU-hr** |

**Total ≈ 12.2 GPU-hr.** SUMMIT cumulative after Phase 2 ≈ **27 GPU-hr against the 150
STOP** — no cap issue, so this is charter-operational and does not need a ruling.

**LAUNCH DEPENDENCY, stated plainly:** **C1, C2 and N all use the +100 rung, which §1.5
has just parked.** They do not fire until #S3 is ruled. **C3 uses only clean dials
(vlbi ✓ §1.3, w0p5 ✓ §0.1, both with clean controls) and is launchable now** — it is the
cell I will fire first when the frozen-M-step arm clears the dgx03 lane, unless you say
otherwise.

**Drain deadline:** every Phase-2 fan must complete by **2026-08-08** —
`AugustDowntime2026` takes the whole cluster 08-10 → 08-13. At ~3 h wallclock per block
on a 1-free-GPU dgx03 there is ample margin, but the list must not be started after 08-06.

---

*(§3 closure test follows as it cuts.)*
