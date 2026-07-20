# SCRIBE — coherence audit of the canonical record

**Agent:** SCRIBE (cronus, read-only). **Date:** 2026-07-15. **Deliverable:** defect list, not fixes.
**Scope read in full:** `STORY.md` (2117 ln), `project_progress.md` (3199 ln),
`CW_transition/trackB_estimator_spec.md` (1962 ln), and all 16 `reports/*.md` (via fan-out extract).
**Step 0:** `git pull --ff-only` → 10 files, +1479 (KINDLE landed). Clean tree at HEAD `d87db93`.

---

## THE ONE-LINE VERDICT

The record is **near-paper-clean on the two master docs** — `STORY.md` and `project_progress.md §10`
agree number-for-number, carry disciplined supersession trails, and their `[MEASURED]` tags survive
spot-checking against the banked reports. **The incoherence is almost entirely a DOC-ORDERING artifact:
two writes happened after the master docs froze and neither was folded in** — (1) the floor-fix doc
session updated STORY + project_progress to `criterion-v2.2 / N_onset=59` but left the **spec a full
version behind**, and (2) the **KINDLE** campaign committed *after* the freeze and closed items the
masters still show open. Fix those two seams and STORY is draftable.

**The timeline that explains everything (git):**
| commit | time (07-13) | what it wrote |
|---|---|---|
| `d587b29` HPC sync | 13:11 | **spec** frozen here — criterion-v2.1 + PROVISIONAL zero-fraction convention |
| `db2075a` floor adopted, onset recut to 59 | 15:34 | **STORY + project_progress** → criterion-v2.2, N_onset=59, e-switch refuted |
| `d87db93` KINDLE | 19:29 | **KINDLE report** → gain statistic retired, D-7 closed — masters NOT re-touched |

The spec was never brought forward from 13:11; the masters were never brought forward from 15:34.

---

## RANKED DEFECT LIST

### D1 — THE SPEC IS A FULL CRITERION VERSION BEHIND (blocking) — cross-doc
`trackB_estimator_spec.md` predates the `criterion-v2.2` adoption and directly contradicts both masters
on the two headline numbers a paper would lead with.

- **D1a — the onset count. HARD CONTRADICTION.**
  `trackB_estimator_spec.md:1897`: *"The corrected onset count is **bounded 57 ≤ N ≤ 67**. … **Do not
  quote 59.**"* (and the whole zero-fraction section is tagged *"consequences PROVISIONAL"*, `:1878`,
  `:1830-1832`). Against `STORY.md:1126` (S6.3.2) *"**N_onset = 59** of 108 cells"* and
  `project_progress.md:3054` (§10.16c) *"SURFACE: N_onset = 59."* The spec actively forbids the number
  the masters feature. Source `reports/RECUT_floors.md:41` confirms the masters (*"N_onset = 59. The
  provisional bound (57 ≤ N ≤ 67) is resolved"*); `reports/FLOOR_FIX_provisional.md:77` is the SUPERSEDED
  origin of the spec's language and is itself banner-marked superseded (`:3`).

- **D1b — no `criterion-v2.2` exists in the spec.** Spec's top canonical criterion is `criterion-v2.1`
  (`trackB_estimator_spec.md:1666`); the floor-validity gate is *"adopted; consequences PROVISIONAL"*
  (`:1878`). STORY treats v2.2 as binding and names the spec as its authority: `STORY.md:20-21`
  (*"full text in `…trackB_estimator_spec.md`"*) then `STORY.md:28` cites *"criterion-v2.2, adopted from
  ANCHOR §4"*. `project_progress.md:3017,3140` declare v2.2 binding. **A reader following STORY's own
  citation lands on a doc that holds the convention provisionally and says "do not quote 59."**
  (KINDLE independently flags the tag is absent from git: `reports/KINDLE_gain_contour.md` §0.)

- **D1c — the e-switch STATUS contradicts.** `trackB_estimator_spec.md:1899-1900`: *"The e = 0.3
  switch-on is **INDETERMINATE** … **Not refuted; not established.**"* Against `STORY.md:1534` (S7.6.4)
  and `project_progress.md:2877` (§10.14b), both of which say the single-member e = 0.3 switch is
  **REFUTED** (collapses to 0.70, below bar). The *bottom-line* number agrees (single member → e = 0.5)
  but the record disagrees on whether it is settled or open. Source: `reports/RECUT_floors.md:122-134`
  (REFUTED) vs `reports/FLOOR_FIX_provisional.md:143-167` (INDETERMINATE, superseded).

- **D1d — a completed campaign is still "queued" in the spec.**
  `trackB_estimator_spec.md:1794` *"CHORUS — PRE-REGISTRATION … queued, no compute yet."* CHORUS has
  fully delivered (`STORY.md:1493` marks the tag struck; `project_progress.md:2848` §10.14 is its full
  write-up; `reports/CHORUS_mixed_pop.md` exists). The spec's record says an executed campaign has not
  run.

- **D1e — the spec contradicts ITSELF on 59.** `trackB_estimator_spec.md:1872` (*"SURFACE then
  measured the box open (**59 onset cells**)"*) vs `:1897` (*"**Do not quote 59**"*) — the SCOPE-OF-
  INFERENCE convention paragraph quotes the number the provisional paragraph 25 lines later forbids.

> **D1 is one fix: fold criterion-v2.2 / RECUT into the spec (retire the PROVISIONAL section, promote
> N_onset=59, mark CHORUS delivered).** Until then, every STORY sentence that cites the spec for a
> convention inherits a provisional or contradictory source.

### D2 — KINDLE CLOSED OPEN ITEMS THE MASTERS STILL SHOW OPEN (moderate) — cross-doc / tag audit
KINDLE (`reports/KINDLE_gain_contour.md`, committed `d87db93`, *after* the masters froze) resolves two
items the masters carry as OPEN. The record is internally inconsistent *right now*.

- **D2a — D-7 (`fALL` re-cut) is CLOSED, not open.** `reports/KINDLE_gain_contour.md` §2 (*"D-7 …
  CLOSED (RECOVERED) … 21 → 21, zero verdict changes, gates 0.000e+00"*). Against `STORY.md:2117`
  (Appendix B, D-7 listed OPEN, *"stands on the estimator that criterion-v2.2 bounds"*) and
  `project_progress.md:3185` (§10.17 open-item 2) and `STORY.md:2104` (ACCRE-CPU pending row). The
  masters' *"probably is not a re-cut"* is now falsified — it WAS re-cut.

- **D2b — the loop-gain "1.000" statistic is RETIRED.** `reports/KINDLE_gain_contour.md` §1.4
  (*"gain = 1.000 is a measured property of a **degenerate** statistic … Ratio gain RETIRED"*). The
  masters don't record the retirement. **Mitigating:** STORY does not actually quote loop-gain 1.000 as
  a live property (the only `gain` hits — `STORY.md:1603` *"next-cycle gain 0.00"* and SIREN's
  *"gain 1.000–1.051"* at `STORY.md:586` — are different statistics), so this is an *unrecorded
  retirement*, not a live stale claim. Low urgency, but the F8 "measured gain 1.000" figure box KINDLE
  flags for retirement is still notionally live.

- **D2c — KINDLE does NOT close D-6, and the masters correctly keep it open.** Verified
  (`reports/KINDLE_gain_contour.md` §1.4, L126-128: the onset-bar `>`/`≥` convention *"stands as an open
  convention question"*). So `STORY.md:2116` (D-6) and the KINDLE(i) pending note (`STORY.md:2096`) are
  **correct** — no defect. Flagged here only to pre-empt a false conflation of the two "1.000"s.

> **Design note, not a fault:** STORY's single-writer rule (`STORY.md:3`) means fenced-agent reports
> (KINDLE) are *expected* to lead the masters until a cronus doc session folds them in. KINDLE stages
> exactly those edits (its §6 add-list). D2 is the un-folded seam, not a contradiction anyone introduced
> carelessly.

### D3 — STORY'S STANDING-CONVENTIONS LIST OMITS FOUR BINDING CONVENTIONS (minor) — convention drift
`STORY.md:22-39` enumerates the conventions that *"govern every number below"* and says their full text
is in the spec. The spec adopted **four more binding conventions on 2026-07-13**
(`trackB_estimator_spec.md:1828-1876, 1940`) that STORY's list does not carry:
1. **ARTIFACT READBACK** (`:1834`) — *"a number is not a result until it is read back off the artifact."*
2. **CODE-BEHAVIOUR CLAIMS REQUIRE A UNIT TEST** (`:1855`).
3. **SCOPE-OF-INFERENCE LINE ON EVERY VERDICT** (`:1866`).
4. **AGENTS PREPARE COMMITS; MATT EXECUTES** (`:1940`).
These are *practiced* in the recent reports (SURFACE/CHORUS/ANCHOR all carry scope/caveat sections; the
anchor-transpose catch at `reports/RECUT_floors.md:232` is artifact-readback in action) but never
promoted into STORY's governing list. Convention practiced, under-declared in the master.

### D4 — STORY'S TWO-CALIBRATED-CELLS TABLE SHOWS A SUPERSEDED FLOOR INLINE (cosmetic) — STORY internal
`STORY.md:1094` (S6.2.1 table) lists the vlbi cell floor as **7.59 ± 0.48 nat**; the very next
paragraph S6.2.2 (`STORY.md:1109-1111`) restates it to the empirical q95 **7.06 ± 0.40** (Gumbel invalid
at zero-fraction 0.45). Trailed at paragraph level, so not a stale claim — but the table itself still
displays the retired Gumbel number, and a table is what gets lifted into a paper. Source
`reports/IGNITE2_softloop.md` §2 (7.59/±0.48) + spec restatement `trackB_estimator_spec.md:1902-1903`.

### D5 — project_progress §10.0 CARRIES criterion-v1's "PERFECTLY PURIFIES" WITHOUT AN INLINE POINTER (cosmetic)
`project_progress.md:2046-2047` (*"The gate does not merely thin the count, it perfectly purifies what is
left"*) is the exact claim `STORY.md:848` (S5.4.1) calls a **census-loudness artifact, dead above onset**.
§10.0 sits under a section-level SUPERSEDED-by-v2 banner (`project_progress.md:2002-2008`) scoping it to
census loudness, so it is trailed — but the specific purity sentence has no inline pointer to §10.8.3
where it dies. Low risk because §10.0 is explicitly a "fixed cell of the map."

---

## TAG AUDIT (STORY.md `[MEASURED]` / `[PENDING]`)

**`[MEASURED]` spot-checks — all PASS** (20+ checked, incl. every tag touched since criterion-v2.1):
- S3.2.1 `4.50 ± 1.48` (range 1–9), strict `1.57 ± 0.98` ✔ `GEO §2` (report L84-85).
- S3.3.2 per-pulsar freqs J1909 38/40, J0437 32/40, J1713 27/40, J1744 20/40, J0711 16/40 ✔ `GEO §3` (L110-127).
- S5.4.5 D4 scorecard: (a) 23/23, (b) 3%/67% (2/77, 39/58), (d) J1909 Δk=−4 R_all 4.65 vs 9.21 ✔ `IGNITE2 §1` / `D4`.
- S5.4.5 D4 (i) fails 8/8; best 90.3% vs 95%; 97.5%-setting flags 44% ✔ `D4 §3`.
- S6.1.2 floors 38.86±1.47 (μ19.82 β6.41) / 7.59±0.48 (μ1.35 β2.10); counts 1.54→0.92, 1.16→0.54; wrong 23/50→14/50, 5/50→3/50 ✔ `IGNITE2 §2`.
- S6.3.2 census onset h*=−12.50, 1.13 [1.10,1.23], floor 106.04±4.62, zero-frac 0.00 ✔ `SURFACE §0/§3/§7`.
- S6.4.3 ANCHOR ladder: 18/40 inflated, median Δq95 −0.18, 1/40 sig, rms 1.993→2.634 ×1.32 ✔ `ANCHOR §3`.
- S6.5.1 defect: 0.845 fitted vs 2.395 empirical, 24σ/12σ, zero-fracs 57/80/93% ✔ `ANCHOR §4`.
- S7.6.2 CHORUS 14.8× lit / 12.4× vlbi (was 11.2×/14.2×) ✔ `CHORUS §1`.
- S9.1.1 SIREN σ(D_L)/D_L 6–12% vs 332%; σ(log10 Mc) 0.866→7e-4–0.03 dex ✔ `SIREN §4.2`.
- S4.1.10 R: f=6.9e-7, ln Z_needle 405617.64/405617.84, ln Z_plateau 405631.83±0.053 ✔ `spec DELIVERABLE R` (L650-667).
- **ATLAS D-1 dispute is REAL and correctly carried:** direct npz read of `reports/atlas_m2m4_summary.npz`
  confirms the stored min-e column holds **0.501/0.526/0.516-class** values while the markdown labels
  them **0.59/0.58/0.66** — exactly the split `STORY.md:1832-1835` (S9.4) and `:2111` (D-1) flag. No defect;
  DISPUTED handled correctly. (Note: the ATLAS *report* itself does not record the dispute — only STORY/
  project_progress do, which is the right place for a cross-source dispute.)

No `[MEASURED]` tag was found citing a number its source does not contain.

**`[PENDING]` tags — quietly closed?**
- **KINDLE** (`STORY.md:2096`) — PARTIALLY closed post-freeze: its D-7 sub-item is done (D2a); its
  gain-diagnosis stage is done (D2b); its *ladder/contour* deliverable is genuinely still pending. Tag
  should split, not simply flip.
- **D-7** (`STORY.md:2104,2117`) — CLOSED by KINDLE (D2a). Should flip.
- ROBUST, SAMPLER, PIPELINE, QUILL, EOB, IMAGING, REAL-ARRAY — all verified still genuinely open (no
  committed report addresses them).
- D-1, D-2, D-6 — verified still genuinely open/unresolved.

---

## CROSS-DOC CONSISTENCY (the triples)

| quantity | STORY | project_progress | spec | verdict |
|---|---|---|---|---|
| **onset count** | 59 (S6.3.2, `:1126`) | 59 (§10.16, `:3054`) | **"do not quote 59"; 57≤N≤67** (`:1897`) | **spec CONFLICTS (D1a)** |
| **criterion label** | v2.2 binding (`:28`) | v2.2 binding (`:3140`) | **v2.1 + PROVISIONAL** (`:1666,1878`) | **spec CONFLICTS (D1b)** |
| **e=0.3 single switch** | REFUTED (`:1534`) | REFUTED (`:2877`) | **INDETERMINATE** (`:1899`) | **spec CONFLICTS (D1c); number agrees** |
| **CHORUS status** | delivered (`:1493`) | delivered (§10.14) | **queued, no compute** (`:1794`) | **spec CONFLICTS (D1d)** |
| **D-7** | OPEN (`:2117`) | OPEN (`:3185`) | n/a | **CLOSED by KINDLE (D2a)** |
| **h*=−12.50 census onset** | S6.3.2 | §10.16c | n/a (predates) | STORY↔pp AGREE |
| **switch thresholds (1/e=0.5, 2+/e=0.3)** | S7.6.4 table | §10.14b table | e=0.5 defensible only | STORY↔pp AGREE; spec partial |
| **loop verdicts (soft=ref impl, hard retired)** | S8.2.4 | §10.10c | criterion-v2.1 (`:1666`) | ALL AGREE |
| **D1 wrong-counterpart hole OPEN** | S5.6.5 | §10.12(3) | `:1756` | ALL AGREE |

Where STORY and project_progress are compared alone, **no triple disagrees**. Every disagreement above
involves the spec (D1) or a post-freeze report (D2).

---

## WRITING-READINESS VERDICT (one page)

**Can the paper be drafted from `STORY.md` alone? — YES for structure and nearly all numbers, with a
short must-fix list first.** STORY is an unusually disciplined document: numbered claims, every number
tagged to a banked source, supersession trails kept rather than deleted, and — verified here — its tags
hold against the reports. `project_progress.md §10` is a faithful long-form mirror; the two never
disagree. A drafter working from STORY will not be misled on any *measured* result.

**What must be fixed first (ranked):**
1. **Fold criterion-v2.2 / RECUT into the spec (D1).** This is the only blocking item. STORY cites the
   spec as its conventions authority, and the spec currently (a) holds those conventions provisionally,
   (b) instructs "do not quote 59," (c) calls the e-switch indeterminate, (d) lists CHORUS as un-run.
   Anyone who checks a STORY convention against its cited source hits a contradiction. One editing pass.
2. **Fold KINDLE into the masters (D2).** Flip D-7 to CLOSED (`STORY.md:2117`, `:2104`;
   `project_progress.md:3185`); record the loop-gain retirement; split the KINDLE tag into
   done-vs-pending. Small, mechanical, but the record is self-inconsistent until done.
3. **Resolve or explicitly defer the three live disputes** a paper cannot straddle: **D-1** (ATLAS M4
   e-column, npz 0.50-class vs markdown 0.59/0.58/0.66 — a real number a siren figure would quote),
   **D-2** (Tier-C f, 0.0323 vs 0.0431 — verdict-invariant, so a footnote suffices), **D-6** (onset
   bar `>` vs `≥`, which moves N_onset 59↔61 — must be *declared*, per KINDLE(i)).

**Cosmetic (does not block a draft):**
- STORY S6.2.1 table shows the retired 7.59 floor inline though S6.2.2 restates it to 7.06 (D4) — swap
  in the table if that table is lifted.
- STORY's standing-conventions header omits the four 07-13 binding conventions (D3) — add them for
  completeness of the "governing conventions" list.
- project_progress §10.0's "perfectly purifies" sentence wants an inline pointer to S5.4.1/§10.8.3 (D5).

**Bottom line:** the *analysis* is paper-ready; the *record* has two un-folded seams (spec, KINDLE) and
three open disputes. None touches a measured verdict. Fix D1 and D2 — an afternoon of doc surgery, no
compute — and STORY alone is a sound drafting base.

---

## APPENDIX — HAZARD-CLASS SWEEP (mission item 1)

Each named hazard, checked for stale quotation *without* a supersession trail:
- **e=0.3 single-member switch** — present only with REFUTED trail (`STORY.md:1534`, S7.6.4). Spec is
  stale-but-labelled-provisional (D1c). No untrailed instance.
- **4.5±1.5 census count** — quoted only as GEO's Bayesian count with the full trail table
  (`STORY.md:335,345`). ✔ clean.
- **0.275/draw** — every instance carries RETIRED/category-error trail (`STORY.md:347,358,2113`). ✔.
- **N_onset ≠ 59** — the only non-59 value ("57≤N≤67") lives in the spec (D1a) and the superseded
  FLOOR_FIX report. Masters quote 59 only. ✔ (defect is spec-side).
- **"gain = 1.000" as a loop property** — NOT present in STORY as a live claim (D2b); retirement
  unrecorded but no stale live quotation.
- **IGNITE h* (−12.75/−13.25)** — always carries [SUPERSEDED] (`STORY.md:1024,1983`, S6.1.1/S11.1.2). ✔.
- **max-of-N floors** — discussed only as retired (`STORY.md:807-813`, S5.3.4). ✔.
- **"no onset exists"** — every instance struck/superseded (`STORY.md:1044,1057,1072`). ✔.
- **IGNITE h* values / pre-v2.2 floors quoted as current** — none found except the D4 table cosmetic.

*Untrailed stale claims in the two master docs: NONE found.* The hazard is confined to the un-folded
spec (D1) and the post-freeze KINDLE items (D2).
