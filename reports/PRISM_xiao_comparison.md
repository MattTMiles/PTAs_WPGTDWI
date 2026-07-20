# PRISM — the positioning analysis: Xiao et al. (2D pulsar-distance inference) against our record

**Agent:** PRISM (cronus, CPU-only, read-and-map). **Date:** 2026-07-16.
**Status:** untracked deliverable. No canonical doc edited. Nothing staged.

---

## P0 — SOURCES, AND WHAT I ACTUALLY READ

| paper | id | status |
|---|---|---|
| **Xiao, Song, Shao, Wang, Zhang & Zhang**, "Two-Dimensional Pulsar Distance Inference from Nanohertz Gravitational Waves" | [arXiv:2512.10729**v3**](https://arxiv.org/abs/2512.10729v3) (11 Dec 2025 v1; **rev. 9 Jul 2026 v3**) | **v3 read in full** — main text + all 8 Supplements, from raw LaTeXML HTML. **v2 also read**; delta reported below. |
| **Yu & Pan**, "Sub-parsec precision measurement of pulsar distances with nanohertz gravitational waves" | [arXiv:2503.23017](https://arxiv.org/abs/2503.23017) | **v1 read in full** (only version; 33 eqs, Tables 1–4). |
| **Wen, Chen, Zhao, Ding & Zhu**, "From Detection to Host Galaxy Identification: Precision CGW Localization with a Few Anchor Pulsars" | [arXiv:2603.28897](https://arxiv.org/abs/2603.28897) (30 Mar 2026; rev. 28 Apr 2026) · also *Chin. Phys. Lett.* **43** 061102 | **secondary scan** — see Appendix B. **Fidelity caveat inside.** |

**The brief said v2; I read v3 and use v3 as canonical.** Two v2→v3 changes are load-bearing for positioning
(§P0.2). The brief's original reading of the method — *the 2D object is the pulsar-pair distance–distance
joint posterior* — is **CORRECT**, and was confirmed independently from the raw HTML before the amendment
arrived.

### P0.1 Our record, as read

STORY.md (whole, 2135 lines); criterion-v2.2 (**authoritative text = `project_progress.md` §10.16**, *not*
the spec — see the caveat below); `SCRIBE_coherence_audit.md`.

> **CAVEAT ON THE BRIEF'S CITATION.** The brief named "the spec's CRITERION v2.2". **The spec does not hold
> v2.2 canonically.** `CW_transition/trackB_estimator_spec.md:1666` tops out at **criterion-v2.1**; v2.2
> appears only as a forward-reference note (`:1671–1676`) pointing *back* at the masters, and the spec holds
> the floor-validity gate as *"adopted; consequences PROVISIONAL"* (`:1878`) while saying **"do not quote
> 59"** (`:1897`) — the number both masters feature. This is **SCRIBE D1b/D1e**, already logged and unfixed.
> I therefore took criterion-v2.2 from `project_progress.md` §10.16 + STORY's standing conventions. *A reader
> following STORY's own citation lands on a doc that contradicts it.*

### P0.2 The v2→v3 delta — two items that matter for positioning

1. **🚩 The headline was silently downgraded; the abstract was not.** v2: *"four GW sources yield sub-parsec
   precision in about **90 %** of the realizations."* v3: *"**~72 %** for N=4 and ~87 % for N=5."* The
   σ_n = 100 ns case went *"roughly half"* → **47 %**. **The abstract's "only a few CGW sources" was not
   revised.** Any spec built against v2's N=4 number is stale by 18 points.
2. **🚩 v3 ADDED a supplement that undercuts its own core mechanism — and it converges on our frontier
   statement.** Supplement III is entirely new. Its Monte Carlo finding, verbatim:
   > *"For every pulsar pair, ~85 % of the realizations satisfy |S_eff − L_q/L_p| < 1.0 … implying that
   > **degeneracy bands from different sources are generally nearly parallel**."*
   > *"In summary, **breaking the degeneracy requires sources in the chirp-mass-dominated regime**
   > (R_sys ≲ 0.25)."*

   The entire method rests on *mismatched* band slopes. Their own new analysis says the slopes are nearly
   parallel in ~85 % of realizations, and that degeneracy-breaking works **only in the chirp-mass-dominated
   corner**. **This is our S11.2.1 frontier statement, in their coordinates, derived independently, under
   referee pressure.** Not reflected in their abstract or conclusions.

Also new in v3: the optimality claim was **retracted** (*"is statistically optimal"* → *"retains all
inter-parameter correlations without information loss from marginalization"*); Supp VI (pairing sensitivity,
factor ~12); Supp VII (prior-domination cliff beyond ~1–1.5 kpc); Supp VIII (localization payoff — which
candidly admits the fiducial 20-yr setup **saturates**, *"making a quantitative localization comparison
between the 2D combination and 1D combination methods uninformative"*, forcing a weaker 10-yr/100-ns config).
**No null test, false-alarm rate, trials factor, or accuracy/coverage check was added in v3.** The referees
pushed on geometry, pairing, reach and payoff — **not on noise calibration.** That gap is identical in v2 and v3.

### P0.3 What their method actually is

**Two-stage, not one likelihood.** *Stage 1:* MCMC (Eryn) over 8 source parameters
(θ, φ, d_L, ι, M, ψ, Φ₀, f₀) **+ Φ_p sampled free in (0, 2π]** — a deliberate reparameterisation to dodge
the oscillatory likelihood (*"we sample them as free parameters to avoid the highly oscillatory likelihood
structure that makes direct sampling in L_p inefficient"*). **Sky (θ, φ) is FIXED AT TRUTH.** *Stage 2:*
post-hoc analytic remap Φ_p → L_p summing over fringe index k (Eq. S2), then the pair joint, Eq. (3):

> PDF_pq(L_p, L_q) ∝ π_pq(L_p, L_q) · Π_{n=1..N} P_{Φ_pq,n}(L_p, L_q)

marginalised over the partner (Eq. 4). **The reported ΔL_p is the 68 % credible-interval half-width of the
best partner, selected by minimising over all 84 candidate partners.** The likelihood is never written down;
priors are inherited by reference to Kato & Takahashi (2026).

---

## P1 — THE TERM-BY-TERM MAP

### P1.1 The structural claim, stated once

Their pair-2D object is a **rank-2 slice of the joint fringe posterior**. Our machinery holds the two
neighbouring ranks and **not** that one:

| rank | object | ours | status |
|---|---|---|---|
| **1** | per-pulsar fringe posterior `q_p(n)` at fixed source θ | **E-step** (star topology, spec §1) | **BUILT, banked.** Exact *conditionally* on θ; **drops cross-pulsar fringe correlation under source uncertainty** |
| **2** | **pair joint `PDF_pq(L_p, L_q)`** | **— none —** | **THEIRS. Built and run.** |
| **N_psr** | full joint (source × integer × distance) | **SAMPLER** | **[PENDING]** — specced (S10.1), never built |

**So the honest form of "ours contains theirs" is: ours contains theirs IN SPECIFICATION, NOT IN DELIVERY.**
SAMPLER's joint marginalised to any pair (p, q) *is* their PDF_pq (up to their two-stage remap
approximation). But SAMPLER is a PENDING tag and their pair-2D is a run result. **We cannot claim containment
by a campaign we have not executed.**

### P1.2 THE VERDICT TABLE

| # | their component | verdict | our correspondent, and the exact relation |
|---|---|---|---|
| 1 | **Signal model / pulsar-term delay.** Eq. (1) Δs = s(t) − s(t_p), t_p = t − L_p(1+Ω̂·p̂) | **SAME** | S1.1.1. Identical object. Their *"the CGW source generates a series of distance-degenerate modes, producing an intrinsically multimodal likelihood"* ≡ our fringe periodicity, dL = λ_gw/(1−cos μ). |
| 2 | **Multi-source degeneracy-breaking.** *"mismatched bands are suppressed, and the posteriors remain mutually consistent only near the true pulsar distances"* | **SAME (mechanism)** | S1.1.1: *"only the true distance phases all sources simultaneously."* **The mechanism agreement is real and it strengthens both records.** |
| 3 | **Chirp mass is the clock.** v3 Supp III: degeneracy-breaking *"requires sources in the chirp-mass-dominated regime (R_sys ≲ 0.25)"*; in the frequency-dominated limit S^(f₀) ≈ L_q/L_p, *"nearly identical slopes"* | **SAME — and CONVERGENT, independently** | **S4.2.16 / S11.2.1:** *"The pulsar term is a kyr-baseline TIMESTAMP of the source's phase. It cannot be read without the CLOCK RATE. The clock rate is fdot, i.e. Mc."* **They reached our frontier statement from the opposite direction and did not notice.** Their R_sys ≲ 0.25 is our clock requirement in band-slope coordinates. **This is the single strongest citation in the positioning.** |
| 4 | **Reparameterise Φ_p free, remap post-hoc** (Eq. S2; Yu & Pan Eq. 26 is the origin) | **SAME idea; THEIRS-ONLY in one property** | Our E-step evaluates lnL **at fringe centres**. Theirs samples Φ_p ∈ (0,2π] free, so the within-fringe offset u is **inferred, not pinned**. → **THEIRS-ONLY #2**, below. |
| 5 | **The pair-2D joint** PDF_pq(L_p, L_q), Eqs. (3)–(4) | **THEIRS-ONLY** | See P1.1. We hold rank-1 and rank-N; **the rank-2 rung is theirs**. Their Supp I states its value: *"the multimodal structure is blurred in the marginalized posterior. In contrast, multidimensional posteriors retain more of the structure from the k-indexed solutions."* → **P5.** |
| 6 | **Multi-source count** ("~4 signals suffice": N=4 → 72 %, N=5 → 87 %) | **SUBSET** | Our **N_CW axis** (S2.4.1) + **channel budget** (MAGPIE J1 ⚠ *readback-pending*). We measured the *threshold and its mechanism*; they have *one instance*. **But the bars are not the same object** — see P2(b). **Channel arithmetic:** a circular source contributes **one** harmonic channel, so their N=5 ≈ **5 channels** against our population threshold of **≳30**. *Their bar is a 1-pc CI half-width; ours is a null-calibrated certification. The 6× gap is a statement about the bars, not (yet) about the physics.* |
| 7 | **Priors.** Gaussian on L_p, σ_p ~ O(10 pc) at ~1 kpc, from **forecast SKA timing parallax** (Lee et al. 2011) | **SUBSET** | Our **tier system** (lit / vlbi, canonical `best_distances.txt`). Their O(10 pc) sits **between our tiers**. We measured what a tier **BUYS** (ΣK 88 454 → 470, S7.2.1) **and what it COSTS** (+2.9 ± 1.0 nat of null floor at h = −13.25, S7.2.8). They have one prior and **no cost accounting, because they have no floor.** **Shared lineage:** both records inherit the parallax-prior de-aliasing idea from **Lee et al. (2011)** — Yu & Pan say so explicitly. |
| 8 | **Multimodality "recovery"** (*"explicitly recovers the multimodal distance information lost in single-pulsar analyses"*) | **SUBSET** | Our E-step's `q_p(n)` **is** the explicit fringe posterior. Ours additionally carries `q_max` as a **certification quantity with a measured reliability curve** (claimed 0.51→0.96 vs realised true-fringe fraction) and **BH-FDR@0.05 → realised true fraction 1.000** (S10.1.3). Theirs is a width; **a width is not a probability of being on the right fringe.** |
| 9 | **Cross-pulsar structure used CONSTRUCTIVELY** (the pair product sharpens L_p) | **THEIRS-ONLY in form — and see the warning** | Ours used cross-pulsar consistency as a **VETO** (D3/D4) and **both were REJECTED** (S5.4). The constructive form **does not evade the failure — it inverts it.** See P1.3. |
| 10 | **Partner selection:** *"We adopt the one with the smallest half-width of the 68 % credible interval as the final constraint"* — **min over 84 partners** | **THEIRS-ONLY — and it is a DEFECT by our own convention** | Our binding convention: *"A calibrated threshold states its false-alarm rate α and its sampling scatter. **An order statistic is not a threshold**."* Their ΔL_p is a **min-of-84 order statistic applied to the reported precision**, uncorrected. **Their own Supp VI measures the spread it selects over: 0.4 pc vs 4.8 pc, factor ~12.** This is **D2.2's lesson (S5.3.4) in the opposite tail**: *"max-of-N has no fixed false-alarm rate; its stringency was an accident of how many nulls happened to be banked."* Theirs is min-of-84 on a 12×-spread quantity, and it is the headline. |
| 11 | **Analytic band-slope diagnostic** S^(f₀), S^(M), R_i, R_sys, S_eff (Eqs. S3–S5) | **THEIRS-ONLY** | We have **no analytic geometric predictor** of fringe-breaking. Our **GEO §4** selection function is *empirical* — stratified ρ(1−cos μ, dlnL) = −0.251 over 4640 (pulsar, draw) cells. **Their R_sys is the closed-form version of what GEO measured by brute force.** → adopt candidate, P4. |
| 12 | **Delivered scale** — 85 pulsars, 1000 realizations/config, full-array count (median 47/85 at N=5, 50 ns) | **THEIRS-ONLY** | **We have never run our stack on a real or SKA-class array.** REAL-ARRAY is [PENDING]; every prior in our criterion is keyed to the 116-pulsar **mock**, whose "residuals" column is *itself an injected CW* (S6.4.2, ANCHOR §0). Their forecast is a forecast — **but it is delivered, and ours is not.** |
| 13 | **Iteration / feedback** | **SUBSET — and it is our strongest generous citation** | v3 *adds*: *"more precise pulsar-distance measurements improve the detectability of additional CGW sources, which in turn enables further distance refinement. **We leave for future work** a more realistic analysis of such cyclic refinement."* **That future work is our S8, and it is done.** We can tell them which implementation cascades (hard lock: 3 → 116, **356 of 359 wrong**, S8.1.2) and which is safe (soft, spec-§3: **no cascade in 40/40**, S8.2.1). |
| 14 | **The two-stage tie-break** (Φ_p sampled free in stage 1, re-tied to L_p in stage 2) | **NEITHER — an unquantified approximation in both papers** | L_p enters the pulsar-term **frequency** as well as the phase (their Eq. 1), so stage 1 samples L_p *and* Φ_p and **breaks the very tie that makes the measurement**; stage 2 re-imposes it post-hoc. Yu & Pan state the over-parameterisation explicitly (*"There are actually only 8+N_psr independent parameters"*) and neither paper bounds the error of the two-stage remap against the correctly-tied joint. **This is precisely what SAMPLER (i) — "no conditioning on truth anywhere" — would do correctly, and precisely what our E-step's star topology makes exact only *given* θ.** |

### P1.3 THE AMENDMENT'S SHARPEST QUESTION: does the constructive form evade the poisoned reference?

**Answer: NO. It does not evade it — it converts a detectable veto failure into a silent bias, and their
selection rule then actively rewards the poison.**

The chain, each link sourced:

1. **D3 failed because the reference is poisoned** (S5.4.4): above onset the leave-one-out reference `u_R` is
   dragged by confident wrong fringes everywhere, so a TRUE cert fails concordance with its own poisoned
   reference. **D4 failed one level up** (S5.4.7): *"D3 failed because the REFERENCE was poisoned; D4 fails
   because the POPULATION IT MUST PROTECT is itself poisoned. Same disease, one level up."*
2. **A veto that cannot distinguish at least FIRES.** Their product does not fire — it **multiplies**. If
   partner q carries a confident wrong-fringe posterior, Eq. (3) places mass at the wrong (L_p, L_q). The
   poison propagates **multiplicatively and silently**. Constructive use is not immune; it is **unprotected**.
3. **And their selection rule is biased *toward* poisoned partners.** They select the partner **minimising the
   68 % half-width**. A confidently noise-locked wrong-fringe partner produces a **narrow** PDF_pq — hence a
   narrow marginal — hence **it wins the selection**. Our S8.1.2 named this object exactly: ***"Tight local
   width + wrong global registration = confident nonsense."*** Their rule selects for the first half of that
   sentence and has no instrument that can see the second.
4. **They could not have caught it**: no null, no wrong-fringe rate, no bias/coverage statistic anywhere
   (P2a). **Their ΔL_p is a width, so a fabricated narrow peak scores as a SUCCESS.**

**Scope line (binding, per our own convention).** This is a *structural* argument from our D3/D4 anatomy and
their stated selection rule. **It is not a measurement of their pipeline.** Our poisoning was measured
above *our* onset, at census loudness, on a 116-pulsar mock with drawn RN+GWB noise. Whether a 50-ns
white-noise SKA array at their loudness populates confident wrong fringes at a rate that matters **is
unmeasured, by us and by them** — and it is exactly what P3 would settle.

---

## P2 — THE CRITICAL DELTAS

*Every item verified against their body text, not their abstract. Absence claims state what was searched.*

### (a) NULL CALIBRATION — **ABSENT IN BOTH. This is the largest single delta.**

**Searched** (Xiao v3, full main text + Supp I–VIII + all figure captions): `noise-only`, `null
test/hypothesis/distribution`, `false alarm`, `background`, `p-value`, `scrambl`, `shuffl`, `sky scramble`,
`phase shift`, `trials factor`, `look-elsewhere`, `Bayes factor`, `odds ratio`, `model selection`, `detection
statistic`, `signal-absent`, `zero-signal`, `coverage`, `P-P plot`, `bias`, `unbias`. **Every `background`
hit is "stochastic gravitational-wave background" or a reference title. Zero otherwise.** Yu & Pan: same
search, same result.

**Their closest passage** — a *prior-only* baseline, which is not a null:
> *"With parallax priors alone, only 4 (1) pulsars meet this criterion at σ_n = 50 (100) ns."*

**And the nearest thing to a sanity limit** (Supp VII):
> *"Beyond ~1–1.5 kpc, the posterior uncertainty approaches the prior limit, indicating that the CGW data no
> longer provide meaningful constraints."*

That is the method correctly returning the prior at large distance. **It is a limit, not a calibration:** no
FPR, no null distribution, no statement of how often a confidently-narrow-but-wrong peak occurs. **No accuracy
statistic of any kind is reported in either paper — no bias, no coverage, no P-P plot, no
fraction-recovering-correct-k. The results sections report widths only.**

**Our two-layer lesson, which is exactly the thing they have no instrument for** (S5.2.1–S5.2.2, FORGE §3):
> **`nullN` — pure noise, NO CW in the data — certified 0.8 pulsars/realisation at the Bayesian bar, and 2 of
> them at P > 0.99.** And the scrambled null *fired harder*: **2.8/realisation with all 16 sources
> scrambled.** *Scrambling all 16 does not reduce the count* — the floor is **intrinsic to the Bayesian
> criterion**, not a coherence leak.

**The one-sentence delta:** *a confident multimodal posterior can be manufactured by noise, and neither paper
has run the experiment that would show it.*

**Fairness, stated because it matters.** **Xiao DRAWS noise** ("1000 realizations", σ_n = 50/100 ns) — a real
advance over its predecessor. **Yu & Pan's likelihood is zero-noise Asimov**: their Eq. (24) sets the data
vector to the **noiseless injected signal** μ(θ₀, t_j); no noise realization is drawn anywhere. **Their entire
Table 4 headline never saw a noise draw.** Our own zero-noise ceiling was **1.350 ± 0.82/draw** (GEO, flat
gate); under real drawn noise with truth off the prior mean, honestly gated, **Arm B detects 0.000**
(S3.2.4). ***That transition — zero-noise to noise — is where our count died. Yu & Pan's headline sits on the
near side of it.***

### (b) DETECTION vs CERTIFICATION, AND TRIALS — **ABSENT**

Their detection is a **strain threshold only**: *"a source is considered detectable if its strain exceeds the
PTA sensitivity curve … for our initial array configuration."* **No detection statistic is ever computed.**
No trials factor, no Bayes factor, no model-selection bar (searched as above). **Their bar is "68 % CI
half-width < 1 pc."**

**Our criterion (v2.2, `project_progress.md` §10.16 / STORY S5.1.1):**
```
DETECTION      dlnL_a > max( ln K_counted,a , floor(h, T, tol) )
CERTIFICATION  q_max,a > 0.9   (strict: > 0.99)   applied ONLY within detections
PURITY         NONE. Both candidate layers tested; BOTH REJECTED.
```
plus the floor-validity gate (Gumbel valid only where the nullN zero-fraction ≲ 20 %; above it, empirical
(1−α) quantile with a **bootstrap** error; **the zero-fraction is a REQUIRED column**).

**The trials number they never apply, estimated in their own units.** Their showcase pulsar J0613−0200 sits
at ~0.99 kpc with σ_p ~ O(10 pc); at their band (3–15 nHz) the fringe spacing dL = λ_gw/(1−cos μ) is ~1 pc
for (1−cos μ) ~ 1. So **K ≈ 6σ_p/dL ≈ 60 fringes in the prior window — our ln K_counted ≈ 4.1 nat per
pulsar, uncorrected.** *Their prior is tens of fringes wide; that is why they need the multi-source
intersection; and there is no bar anywhere in the paper against which the resulting confidence is scored.*
**Add the min-over-84 partner selection (P1.2 #10) and the look-elsewhere is applied twice, in the same
direction, on the headline.**

**Our convention, which is the whole of this delta** (S5.2.4):
> ***"Confidence without a detection statistic is prior-pinning in disguise. Every confidence bar must sit
> downstream of a detection statistic that CAN RETURN ZERO. A criterion that cannot fire on a null is not a
> criterion."***

**A width cannot return zero.** Their criterion cannot fire on a null, because widths always exist.

### (c) THE WALLS — **ASSUMED, and they say so cleanly. But we RAN their assumption.**

They are explicit and, to their credit, pre-emptive:
> *"we assume that the sky locations of the CGW sources used for pulsar-distance inference are known and
> therefore treat their polar and azimuthal angles as fixed parameters. This does not rely on precise PTA
> localization as an input to our method. Such sky positions can be provided either by targeted CGW searches
> toward SMBHB candidates suggested by electromagnetic observations … or by all-sky CGW searches followed by
> host-galaxy filtering and ranking within the PTA credible sky region, which can in favorable cases reduce
> the candidates to a single host."*

**We agree the sky is free — we measured it.** S4.1.11: break-even sky prior **θ\* = 0.188°/source**, against
the F-stat blind floor ~6° — a **32× gap** — and *"an EM host supplies exactly this."* S7.3.4: **EM astrometry
beats the 22-arcsec tolerance by ≥10⁴×.** ***"The counterpart doesn't cost you; the loudness does."***

**But their assumption is licensed only inside our targeted arm — and we ran that arm.** S4.2.8, the targeted
referendum, sky **exact**, frozen protocol:

| tier | conditioning | **f** | gate |
|---|---|---|---|
| A | sky only | 0.0847 ± 0.131 | FAILED |
| B | + EM period | 0.0481 ± 0.0227 | FAILED |
| C | + host D_L | **0.0323 ± 0.0134** *(frozen 4-seed, PRIMARY per D-2)* | FAILED by **16.1×** |

> ***~97 % of the posterior's evidence sits on the wrong-fringe plateau EVEN WITH THE SKY EXACT. The tier
> gradient is FLAT and mildly INVERTED — f falls as conditioning tightens. Counterpart information does not
> help.***

**And there is a circularity in their v3 worth naming.** Their Supp VIII sells improved **sky localization**
as the payoff of certified distances — while the method **requires** sky positions as fixed external input.
They pre-empt it (*"This does not rely on precise PTA localization as an input"*), leaning on EM hosts. Fair.
But Supp VIII then admits the fiducial config **saturates**: *"already far exceeds 9 across nearly all
realizations, saturating at the diffraction limit for both methods and making a quantitative localization
comparison between the 2D combination and 1D combination methods uninformative."* **They had to weaken to
10 yr / 100 ns to show any 2D-vs-1D difference in the payoff at all.**

**Scope line.** Our f = 0.0323 is at **census loudness, 116-psr mock, N_CW = 1 star topology, circular**.
Theirs is 85 psr, 50 ns, N = 2–5. **Our number is not a refutation of theirs — it is a different cell.** The
*structure* — sky-exact does not rescue fringe identification — is ours, is measured, and is what P3 exists
to test in their cell.

### (d) WRONG COUNTERPART / WRONG ASSOCIATION — **ABSENT**

Searched Xiao v3 (main + Supp I–VIII + captions): `wrong`, `mis-identif`, `misassoc`, `incorrect`,
`spurious`, `offset`, `accura`. **`spurious` appears only as the thing the method claims to suppress, never
as a measured rate.** Their host-galaxy conditioning carries the hedge *"in favorable cases"* — **no rate is
attached**.

**Ours, measured:** the **D1 hole is OPEN and structurally un-closable by co-registration** (S5.4.7,
S5.6.4–5). IGNITE's scrambled-source loop **detects in 2 of 5** realisations under the `fN` floor (dlnL up to
~15 nat vs floor 5.46); IGNITE-2's soft loop: **6 of 10** scrambled realisations certify at some iteration,
**3 keep one to the fixed point**, and **one clears even the blind-robust `fALL` floor by 0.15 nat** — *even
the blind-robust column leaks this event.* The implementable D4 form catches only **25–52 %** (S5.5.2),
because *"σ_EM/dL is O(150–800) fringes in the lit tier: a wrong counterpart does not look displaced relative
to a prior that was never going to localise the fringe anyway."*

**And the price of the fork is on our record, permanently** (S5.6.3): under `fALL` — the
wrong-counterpart-robust null family — **there is NO onset anywhere in the modelled grid.** *That is the
honest cost of the targeted framing, and Xiao et al. carry the same exposure with no column for it.*

### (e) DYNAMICS — **ONE-SHOT. And v3 names our loop as future work.**

Verbatim, new in v3:
> *"In principle, more precise pulsar-distance measurements improve the detectability of additional CGW
> sources, which in turn enables further distance refinement. **We leave for future work** a more realistic
> analysis of such cyclic refinement within a full multiple-source search."*

**That is our S8, and it is measured — including the safety boundary they would need.**

| | hard lock (retired) | soft loop (spec §3, reference impl.) |
|---|---|---|
| cascade | **3 → 116 in three iterations; 356 of 359 certs WRONG** | **none in 40/40 trajectories** |
| src_mc_off | up to **1.6 dex in a single step** | **< 1e-4 dex** at truth |
| wrong-cert count | grows | **never grows** |
| scrambled source | 2/5 detect, all wrong | **self-cleans 3 of 6** |

**The mechanism** (S8.1.3): the hard M-step **pins each certified pulsar at its MAP fringe centre** — up to
half a fringe off — biasing the (f, mc) gradient; one damped Newton step moves the source materially; the
re-E-step re-registers; the poison compounds. ***A GPS wrong-fix failure at loop level.*** The spec's own
soft, fringe-posterior-weighted M-step **predicted exactly this failure, and the first implementation did not
follow it.**

**This is the most generous and most useful thing we can tell them:** their proposed future work has a known
failure mode, a known safe implementation, and a known ceiling — ***"The loop works given seeds. The modelled
box supplies none"*** (S8.3.2), and ***"THE LOOP IS NOT THE PROBLEM. THE CRITERION IS"*** (S8.2.5).

### (f) REGIME HONESTY — **SKA-era, white-noise-only, given-sky, given-detection**

| axis | theirs (Xiao v3) | ours | delta |
|---|---|---|---|
| noise | **50 or 100 ns white Gaussian**. *"Our analysis is performed in a simplified PTA setting that includes a CGW signal and white noise."* No red noise, no SGWB, no overlapping CGWs, no DM, no timing-model marginalisation | drawn heterogeneous white + per-pulsar RN GP + **HD-correlated GWB GP** | **S10.2.1: *"Noise is NOT a perturbation"* — median per-pulsar CW rms 458 ns vs drawn noise 2005 ns (white 1414, RN 805).** They omit the two components carrying ~1000 ns. |
| array | 85 psr, 20 yr, 2-wk cadence | 116-psr mock, T = 15 native | comparable in size; **theirs is SKA-era in precision** |
| parallax priors | σ_p ~ O(10 pc), **forecast SKA** (15 μas VLBI vs current 45; 10–70 μas timing vs current 100) | lit / vlbi tiers from real catalogue | **theirs is a forecast capability, ours a measured one** |
| sky | **fixed at truth** | free (blind) or supplied (targeted arm) | P2(c) |
| detection | **assumed** (strain > sensitivity curve) | `dlnL > max(ln K, floor)` | P2(b) |
| **source count** | **N = 2–5 simultaneously resolved** | — | **see below** |

**THE REGIME NUMBER THAT MATTERS MOST, AND NEITHER PAPER PRICES IT.** Their headline needs **N = 4–5
simultaneously resolved, sky-localized CGW sources**. **Zero are currently detected.** Our **SCOUT** population
clock (S7.3.4): the mean number of individually detectable SMBHBs is
> **N̄ ≲ 0.01 (current), ≲ 0.1 (any SMBH distribution consistent with the SGWB), and O(0.1–1) in the SKA era.**

***Their headline input is 5–50× above the best population estimate for the very era they assume.*** Their
catalog is **monochromatic — every one of the 100 sources at Mc = 5×10⁹ M_⊙, d_L = 1 Gpc** — and **no
per-source SNR is reported anywhere in the paper** (SNR appears only in Supp VIII, as an *assumed* SNR = 10
plugged into an analytic scaling). Add **SCOUT's structural scissors** (S7.3.3): our −13.75 operating point
already sits **~0.35 dex ABOVE the NG15 sky-averaged CW upper limit** — *"any source genuinely at the
operating point, in band, would already be an individual-source DETECTION."*

**Given-seeds vs achievability.** Their paper is structurally **SIREN** (S9.3.5): a given-resource forecast.
Our own version carries the caveat verbatim (S9.3.4): *"These are Cramér–Rao bounds on zero-noise Asimov data
with the fringe integers GIVEN … **Every σ SIREN quotes is a LOWER BOUND**"* — and the matching honesty
(S9.3.5): ***"SIREN's N_seed = 3–5 columns currently price a resource the noisy pipeline CANNOT DELIVER AT
ALL. The payoff is real; the road to it is not through cold-start certification."*** **Xiao et al. price
N = 4–5 resolved sources with no equivalent line.**

---

## P3 — THE CROSS-CHECK CELL SPEC *(for a future ACCRE task — DO NOT RUN)*

> ⚠ **READBACK-PENDING.** This spec is written against SPARK's staged results as **relayed in-session**
> (`reports/SPARK_launch_criterion.md`, ACCRE, not yet pushed — absent from cronus after `git fetch`).
> **Our binding ARTIFACT READBACK convention (adopted 2026-07-13) states: *"a number relayed in-session
> before its npz exists is a guess with a decimal point."*** Every SPARK-derived number below is tagged ⚠ and
> **must be read off its artifact before this spec is executed or quoted.**

### P3.1 The cell

**Their headline scenario, read from their prose** (there are **no tables in the paper** — 3 main + 6
supplemental figures, all numbers in prose; the brief's "read from their tables" cannot be satisfied):

| quantity | value | our nearest banked analogue |
|---|---|---|
| pulsars | **85**, positions from current PTA datasets | 116-psr mock |
| timespan | **20 yr**, 2-wk cadence | T = 20 in the SURFACE box |
| noise | **50 ns white only** | drawn white+RN+GWB, ~2005 ns rms |
| sources | **N = 5**, from a 100-SMBHB catalog | N_CW axis |
| chirp mass | **fixed 5×10⁹ M_⊙** (all sources) | census structure 3+13 |
| frequency | **3–15 nHz** | fiducial band |
| d_L | **1 Gpc** | — |
| distance prior | Gaussian σ_p ~ O(10 pc) | **between lit and vlbi** |
| **headline** | **~87 % of realizations reach sub-pc on J0613−0200; median 47/85 pulsars sub-0.5 pc** | — |

### P3.2 What to compute

Re-score **their claim, not their pipeline**. The deliverable is: *what verdict would criterion-v2.2 return
for a precision claim of this shape?*

1. **Build the cell**: 85-psr subset, T = 20, σ_n = 50 ns **white-only** (deliberately theirs, not ours — the
   point is to score their regime), N_CW = 5 circular at Mc = 5×10⁹, d_L = 1 Gpc, f ∈ [3,15] nHz, sky **fixed
   at truth** (their assumption), Gaussian σ_p = 10 pc priors (a **new tier — call it `ska10`**, between lit
   and vlbi).
2. **Floors** (the whole point): fit `floor(h, T, tol)` at this cell from **N ≥ 100 nullN** (S5.3.5 sizing:
   scale-free sd < 10 % of the floor; **N ≥ 150 at onset cells**). **BANK THE ZERO-FRACTION** — it is a
   required column, not a caveat (criterion-v2.2). **If zero-fraction > 20 %, the Gumbel is void: quote the
   empirical q95 with a bootstrap error.**
3. **Both null families.** `fN` operative; **`fALL` travels beside it** (S5.6.1) — non-negotiable, and it is
   where their wrong-association exposure (P2d) would surface.
4. **Score**: correct certs/realisation at floor and at floor+error; **wrong-cert count beside it** — *an
   above-onset count without its purity number is meaningless* (S5.4.1).
5. **Their statistic, in our units**: reproduce their ΔL_p (68 % half-width, **min over 84 partners**) and
   report it **alongside** `q_max` and `dlnL`. **The comparison of interest is not the width — it is the
   fraction of narrow widths sitting on the WRONG fringe.**
6. **The look-elsewhere term they omit**: report ΔL_p under (i) their min-over-84 rule and (ii) a
   **pre-registered fixed partner**. The difference **is** the uncorrected trials factor on their headline.
7. ⚠ **SPARK's certified-coherent detection statistic** (MaskedDelay wired into an ncw=1 detector,
   bit-identical gate) — *if* it generalises to ncw=5, it is the natural instrument for step 5.
   **READBACK REQUIRED before this step is specced further.**

### P3.3 Pre-registered expectations — **BOTH WAYS, before compute**

**Expectation A — their claim survives our bar (the honest case for it):**
- Their σ_p = 10 pc at dL ~ 1 pc gives **K ~ 60 → ln K ≈ 4.1 nat**: a *low* trials bar.
- 50 ns white-only across 85 pulsars is **~40× quieter than our drawn noise**; SNR_pterm would be far above
  our array's median 0.47 (S2.3.1), so the readable sub-array is much larger than our 30/116.
- **T = 20 yr with 5 sources sits in a regime our box never sampled**: our onset is baseline-driven
  (`T^{5/2}` beats the `h^{1.66}` floor race) and our N_CW axis tops out at 16 with the *census* loudness
  structure. **SURFACE's structure lever is the relevant precedent: 3+13 → 5+11 raised the count up to 6.1×
  and moved the frontier ≥ 0.75 dex** (S7.6.2). Five loud sources is *closer to 5+11 than to 3+13.*
- **If the floor at 50 ns white-only is small, their claim could clear our bar.** *This is a live possibility
  and the spec must not be written to foreclose it.*

**Expectation B — their claim does not survive (the case against):**
- **The floor is loudness-relative and template-driven**: `floor_fN ∝ h^1.66`, and the mechanism *"runs on
  data containing NO CW at all"* — the E-step's matched-filter cross term is linear in the **model**
  amplitude (S5.3.1). At Mc = 5×10⁹ / d_L = 1 Gpc the hypothesis is **loud**, so the floor is **high**, and
  their quiet noise does not lower it. **ANCHOR reproduced the exponent from scratch at 1.67 at a baseline
  where it was never fitted** (S6.4.4) — *"the tail is template-dominated; the body is noise-dominated."*
  ***Their 40× noise advantage buys them less than it looks like, and this is the crux.***
- **Louder is non-monotone**: S5.3.2 — *"a 10× louder source LOWERS the honest count."*
- **Their 5 circular sources ≈ 5 channels against our ≳30-channel population threshold** (⚠ MAGPIE J1,
  readback-pending).
- **Their own Supp III**: bands nearly parallel in ~85 % of realizations; degeneracy-breaking needs
  R_sys ≲ 0.25. **Their monochromatic Mc catalog cannot populate the chirp-mass-dominated regime their own
  mechanism requires** — every source has the *same* Mc, so the S^(M) lever has no diversity to exploit.
  ***This is the internal tension most likely to decide the cell, and it is theirs, not ours.***

**Pre-registered STOP:** if the cell's nullN zero-fraction exceeds 20 %, **the Gumbel is void** — quote the
bootstrap empirical quantile and say so. *That defect cost us CHORUS's loudest headline; it is not to be
rediscovered.*

**Cost note (binding, and it decides where this runs).** **`b1_L_gwb` is ACCRE-only and gitignored**
(§10.16e). ***A machine without it — cronus — CANNOT reproduce any banked noisy number*** and will hit
`RECOMPUTED-UNSAFE` and draw a rotated GWB realisation. **This cell is ACCRE, not cronus.** *A cronus session
that sees that warning and proceeds anyway is fabricating a number.*

**Honest scope of what P3 can deliver.** It scores **a claim of their shape in a cell of their construction
under our criterion.** It is **not** a re-run of their pipeline, does not use their code, and cannot
adjudicate their two-stage remap (P1.2 #14). **It answers exactly one question: does a 1-pc CI half-width,
obtained this way, survive a bar that a null cannot clear?**

---

## P4 — THE POSITIONING PARAGRAPH *(drafted for the paper)*

> ⚠ **The final sentence is READBACK-PENDING** (SPARK, relayed not banked). Draft below marks it.

**Draft (4 sentences + the adoption sentence):**

> Pulsar-distance inference from the phases of multiple resolved continuous-wave sources has developed
> rapidly and independently: Yu & Pan (2025) proposed combining phase information across individually
> resolved SMBHBs and, building on the parallax-prior de-aliasing of Lee et al. (2011), demonstrated
> sub-parsec forecasts for a 20-pulsar SKA-era array; Xiao et al. (2026) extended this to the two-dimensional
> pulsar-pair distance–distance joint posterior, showing that the degeneracy bands of different sources
> intersect to suppress spurious modes, and that the multimodal structure a marginalised one-dimensional
> posterior blurs is retained in the higher-dimensional object. **We agree on the mechanism, and the
> agreement is worth stating plainly: their band-slope analysis concludes that breaking the distance
> degeneracy requires sources in the chirp-mass-dominated regime, which is — in different coordinates — the
> conclusion this work reaches from the pulsar-term side, namely that the kyr-baseline pulsar term is a
> timestamp that cannot be read without the clock rate, and that the clock rate is the chirp mass.** Our
> contribution is not a competing estimator but the surrounding structure within which any such estimator,
> including theirs, must operate: the feasibility walls and the onset surface that say *where* in
> (loudness, baseline, prior, population) a distance claim can be made at all; a criterion calibrated against
> noise-only and scrambled data, which is necessary because a confident multimodal distance posterior can be
> manufactured by noise alone — our Bayesian bar certified 0.8 pulsars per realisation on data containing no
> continuous wave at all, and 2.8 per realisation under a scrambled source; the dynamics and safety boundary
> of the cyclic refinement that Xiao et al. leave to future work, where we find that the natural hard-locked
> implementation cascades to 356 wrong certifications in 359 while the fringe-marginalised implementation is
> stable in 40 of 40 trajectories; the transition-region limits that bound where distance information
> survives at all; and an honest accounting of sampling, in which every quoted precision carries the
> population it generalises to and the population it does not.
>
> ⚠ **[READBACK-PENDING — SPARK]** What we would adopt from them is concrete: the pair-level
> distance–distance joint is a rank-2 object sitting between our per-pulsar fringe posterior and the full
> joint we specify but have not built, and their closed-form band-slope diagnostic (R_sys) is the analytic
> form of a selection function we measured only empirically; **conversely, the mechanism our record adds and
> theirs cannot express is selectivity — that the rigid two-component Earth-plus-pulsar template is
> unforgeable by noise, so that coherence lowers the null floor rather than merely raising the signal.**

**Notes for Matt on the draft.**
- **It cites generously and the mechanism agreement leads.** Their Supp III convergence is the strongest card
  in the whole analysis and it is *theirs*, which is precisely why it is credible.
- **It does not attack their numbers.** Every differentiator is a *thing we have and they lack*, not a claim
  their result is wrong. **We have not measured their pipeline and must not imply we have.**
- **The two sharpest available attacks are deliberately NOT in the paragraph** — the min-over-84
  look-elsewhere on their headline (P1.2 #10) and the poisoned-partner selection (P1.3). **They are correct
  and they are ours to make, but they belong in a referee report or a discussion section, not in a
  positioning paragraph that opens by citing them generously.** *Recommend holding both.*
- **The "SUBSET" framing is softened deliberately.** P1.1's honest finding is that our containment is
  **specification, not delivery** (SAMPLER is PENDING). *"Not a competing estimator but the surrounding
  structure"* is true and defensible; *"ours contains theirs"* would not survive a referee who asks to see
  the joint sampler.
- **The final sentence must not ship until SPARK's npz is read back.** As written it makes a mechanism claim
  ("coherence lowers the null floor", ⚠ relayed as 19.5 → 7.5 nat) that is exactly the class our ARTIFACT
  READBACK convention was adopted to stop.

---

## P5 — IS THE PAIR-2D WORTH ADOPTING AS A CERTIFICATION AID?

### P5.1 The verdict: **PLAUSIBLY YES, for a reason that is NOT theirs — and the cheap test is specced below.**

**Not as inference.** Adopting their *constructive* use would import the poisoned-partner failure (P1.3) into
a pipeline whose entire value is that it does not have one. **Rejected.**

**As a diagnostic — the pair-level `q` statistic the amendment asks for.** Our E-step's `q_p(n)` is
**conditionally exact given θ** (star topology, spec §1) and **drops cross-pulsar fringe correlation under
source uncertainty**. Their PDF_pq is the object that *carries* the pair term. **The difference between them
is not noise — it is exactly the quantity our rank-1 approximation discards.** So define

```
I_pq  =  KL[ q_pq(n_p, n_q)  ‖  q_p(n_p) ⊗ q_q(n_q) ]
```

— the fringe-integer mutual information between pulsars p and q induced by marginalising the source. **`I_pq`
is a SENSITIVITY statistic, not a CONSISTENCY statistic**, and that distinction is the whole argument:

- **D3/D4 asked "do p and q AGREE?"** — a consistency question, which **needs a reference**, and the reference
  is poisoned above onset (S5.4.4). **D4 failed one level up because the population itself is poisoned**
  (S5.4.7), so *"NO co-registration statistic can close the D1 hole in general."*
- **`I_pq` asks "does p's fringe answer DEPEND on q's, through θ?"** — it **has no reference**. It compares two
  objects computed from the *same* realisation, one of which (the product) is our own approximation. **The
  poisoned-reference mechanism does not obviously apply.** *That is a genuinely new mechanism, not a re-run of
  a rejected one — which is the only thing that justifies spending compute after two pre-registered
  rejections.*

**Why it could certify.** A per-pulsar `q_max` is confident **because it conditioned on a θ it does not
know**. If `I_pq ≈ 0` across p's partners, that conditioning is harmless and the E-step certification is
sound. If `I_pq` is large, `q_max` is **confidence borrowed from a fixed θ** — which is the *exact* failure
S8.1.2 named: ***"Tight local width + wrong global registration = confident nonsense."***

**The honest worry, stated before the test.** A true pair at a poorly-determined θ **also** shows large
`I_pq` → false flags → **D4's condition (ii)** (≤10 % false-flag) is where this most likely dies. **Our
prediction is FAIL on that condition** — recorded here, before compute, so it cannot be claimed afterwards.

### P5.2 The cheapest test on banked E-step outputs

**Inputs — all already banked, no new realisations:**
- The **lit onset cell (−12.75, 30, lit)**: 50 realisations, **14/50 carry a wrong certification** under the
  fresh floors (S5.4.9) — a **labelled** wrong/right population, which is what makes this cheap.
- The **vlbi cell (−13.25, 30, vlbi)**: 3/50 wrong — the low-impurity contrast.
- **SURFACE's reserved Pair B** — (−12.75, 40, vlbi, 5+11) at 4.07/real and (−13.00, 40, vlbi, 5+11) at 3.57,
  **wrong-cert 0.07–0.13**: *the first genuinely above-onset, low-impurity seed set the programme has ever
  had, and both untouched by the floor fix* (STORY App. A, KINDLE). **The natural high-signal arm.**

**Method (CPU, hours not days):**
1. Read banked `q_p(n)` and the E-step columns. **Do not re-run the E-step** — that would need `b1_L_gwb`
   and is therefore **ACCRE-only** (§10.16e). *If the banked q-columns suffice, this runs on cronus; if not,
   it is an ACCRE task. Determine this first — it decides the whole cost.*
2. Propagate source uncertainty by **linear response / Laplace** about the banked source, using the B1 pilot's
   HVP-assembled Earth-only Fisher **inverted WITH the generative priors, never pinv** (S4.2.3 — *pinv reports
   σ = 0 for an unconstrained axis*). **This is the approximation that makes it cheap**, and it is the first
   thing a referee will attack: it assumes the pair correlation is captured to first order in δθ. **State it.**
3. Form `q_pq` to first order, compute `I_pq` for each certified p over a **pre-registered** partner set
   (e.g. the 5 largest-lever partners — *fixed before scoring, never selected*; **their min-over-84 rule is
   the error we are diagnosing, not one to repeat**).
4. **Score BOTH the ORACLE and the IMPLEMENTABLE forms** (S5.5.1, binding). *Every D3 number was oracle-form;
   the implementable form is 2–4× weaker. This caveat travels backward onto anything new.*

**Pre-registration (fixed now, inherited from D4 so it cannot be tuned):**
- **PASS** iff `I_pq` catches **≥95 %** of wrong certs at a **≤10 %** false-flag rate on true-signal
  realisations, **in the implementable form**, at **both** cells. **No partial adoption. No threshold
  tuning.** *(D3 and D4 were both rejected by their own pre-registration; this inherits the bar.)*
- **PREDICTED: FAIL on the false-flag condition**, for the S5.4.6 reason — at the lit onset cell **36 % of
  true-signal realisations already have an impure detected set**, and *"the gate faithfully measures the
  cell's own impurity and cannot beat it."* **If `I_pq` false-flags at ~36 % it has rediscovered D4, and the
  answer is that the pair-2D adds nothing a veto can use.**
- **The result that would make it worth it:** `I_pq` false-flags **well below** the cell's own impurity —
  which would show it is *not* measuring impurity but source-conditioning sensitivity, i.e. **a genuinely
  different axis**, and the first statistic in the programme that is.

**⚠ Artefact caveat that must be checked FIRST, or the test is void.** Our banked E-step evaluates lnL **at
fringe centres**, pinning the within-fringe offset u = 0 — which **fabricates fringe evidence on pure
noise** (null fires ~30 % vs 3 %; floor 3.59 vs 0.27 nat). The floor refit **absorbs** it, so banked counts
survive, but the floor/FPR *meaning* is off. **`I_pq` computed on u=0-pinned posteriors inherits the
artefact.** *Their Φ_p-free reparameterisation does not have it (P1.2 #4) — which is the irony of this whole
section: **the thing most worth adopting from them may be the sampling reparameterisation, not the pair
object.*** **Decide this before spending anything on `I_pq`.**

---

## Appendix A — FINDINGS AGAINST OUR OWN RECORD

**A1. S3.1.1 needs a frequency scope line, by our own convention.** STORY S3.1.1 states: *"The
prior-alone-anchor assumption of arXiv:2603.28897 **FAILS** on the real array."* **That is true as measured
and under-scoped as written.** Our census `K = 6σ_d/dL` was evaluated at our fiducial `dL ≈ 0.36–0.49 pc`.
**Wen et al. work at f_GW = 4 nHz, where λ_GW ≈ 2.4 pc** — fringes ~5× wider, so K is ~5× smaller at fixed
σ_d. Their D_err = 1 pc anchors sit near **K ~ 2.5**, not K ≤ 1 — and their threshold (`D_err < λ_GW`) omits
both our 6× factor and the (1−cos μ) geometric factor. **Separately, their 1 pc is an assumed forecast**
(*"Forecasting the exact astrometric capabilities of future facilities … is beyond the scope of this
paper"*), while our census measures **current** priors. **The two claims are not in contradiction; they are
about different frequencies and different epochs.** Our own **SCOPE-OF-INFERENCE convention** (*"every verdict
states the population it generalises to and the population it does not"*) applies to us here.
**RECOMMEND** (not applied — PRISM edits no canonical doc): scope S3.1.1 to *"at our fiducial dL, under
current priors."*

**A2. The two named sources in the brief do not exist on cronus.** `reports/SPARK_launch_criterion.md` and
`MAGPIE_audit.md` are **absent after `git fetch --all`** — no file, no commit, no branch, zero mentions in any
tracked file. Matt confirms both are **staged on ACCRE, unpushed**. Their content was **relayed in-session**
and is used above **under the ⚠ readback tag**, per the binding convention: *"a number relayed in-session
before its npz exists is a guess with a decimal point."* **`ch_README`, MAGPIE J1's cited authority for the
per-member channel counts (e=0.3→8, e=0.5→17, e=0.7→32), is likewise not in the repo.**

**A3. The brief's "read their tables" is unsatisfiable.** **Xiao et al. contains no data tables** — the 10
`<table>` tags in its HTML are LaTeXML equation wrappers. It is a Letter: 3 main + 6 supplemental figures,
**all scenario numbers in prose**. (Yu & Pan *does* have Tables 1–4; those are read and reproduced in the
source material.)

---

## Appendix B — arXiv 2603.28897 (Wen et al.) vs RING and our anchor/census work

**Fidelity caveat, stated first.** The Wen et al. body text was reachable only through a summarising fetch.
Load-bearing quotes were cross-checked across three fetches, **but character-level verbatim fidelity is NOT
certified** for body quotes (the abstract is). Two symptoms: the λ criterion returned in two phrasings, and
the sky-area range differs between abstract (**~0.1–9.2 deg²**) and Discussion (**~0.1–7.6 deg²**) —
**unresolved**; possibly a v1/v2 revision. The paper uses **named, unnumbered sections**. **Verify against the
PDF before quoting onward.**

**What they do.** 25-pulsar array, **20 yr**, white+red noise from real datasets, **SNR = 20** primary; a
minimal subset of **3 anchors {J0030+0451, J0437−4715, J2222−0137}** (or 6, adding J0636−3044, J1744−1134,
J1909−3744) at **D_err = 1 pc**, with the **other 22 at 20 % fractional error**. Threshold: **D_err < λ_GW ≈
2.4 pc at 4 nHz**. Claim: a few anchors *"is sufficient to phase-lock the array and drastically shrink the
sky-localization error"* — **factor ~30** in favourable directions. **Distances are PRIORS, not inferred**
(`D_p = D_p,0 + D_err · D_prior`). **They cite neither Yu & Pan nor Xiao et al.**

**RING's bias result lands on them squarely, and this is the paragraph:**

> **RING (S7.2.4) measured that bad distance priors BIAS the sky — they do not merely broaden it.** Every
> non-exact prior drives the sky MAP **3–6° off truth, INDEPENDENT of SNR**, while the 90 % area shrinks
> 4–17× per SNR doubling, so **coverage DEGRADES as the signal gets louder**: `inside90` = **0.90 → 0.50 →
> 0.00** at SNR 5/10/20. ***This was proven, not inferred***: a zero-noise control gives bias **2.73–5.28°**
> for every imperfect tier and **exactly 0.0000°** for the exact tier, in all four configurations. The
> mechanism was isolated: **bias ∝ un-modelled pulsar-term power, and nothing else** — it collapses to 0.033°
> only when κ̄ reaches **0.290**.
>
> **Wen et al. sit below that collapse point by construction.** With 22 of 25 pulsars at 20 % fractional
> error (~200 pc at 1 kpc ≫ their own 2.4 pc threshold), those 22 contribute κ ≈ 0, so κ̄ ≈ **3/25 = 0.12**
> (or **6/25 = 0.24** for the six-anchor set) — **both below 0.290.** RING therefore predicts a **residual,
> SNR-independent sky bias** in exactly their configuration. Their reported 90 % areas (~0.1–9.2 deg², i.e.
> radii ~0.2–1.7°) are **smaller than the bias RING measures**, which is RING's undercoverage signature.
>
> **And their own paper contains the symptom, explained away as noise:** *"We also note a slight deviation
> between the inferred Right Ascension (RA) and the injected value in Figure 1. This marginal offset arises
> from random noise realization, which shifts the maximum likelihood peak slightly away from the true injected
> parameters."* ***RING's zero-noise control exists precisely to catch that misdiagnosis: at zero noise, with
> imperfect distances, the offset is 2.73–5.28° and there is no noise to blame it on.***
>
> **Worse, theirs is RING's named worst regime.** S7.2.5: *"The worst regime is PARTIAL distance
> knowledge"* — tier1 → tier2_k3 at SNR 10 shrinks area90 **6.5×** while the zero-noise bias barely moves
> (2.73° → 3.07°), so `inside90` falls 0.60 → 0.30: **the credible region contracts around an offset that
> does not.** *A few anchors among many bad distances is the definition of partial knowledge, and their
> headline **factor ~30 shrink** is what that contraction looks like from the inside.*

**Caveats that bound the above — it is a PREDICTION, not a refutation:**
- **RING ran fgw = −8 (10 nHz); Wen et al. are at 4 nHz (−8.4).** **RING's own stop-point S-3 says its harness
  is internally inconsistent in the timing-model prior (factor ≈ 19) — harmless at −8 (KS p = 0.73) but
  broken at −9 (KS p = 0.0000): *"treat every fgw ≲ −9 NOISY result from this harness as uncalibrated."*
  4 nHz is between. **RING is not certified at their frequency.**
- At 4 nHz the coherence threshold relaxes (**S7.2.3**: σ_d < 3.02 pc at −8; σ_d < 30 pc at −9), so **their
  1-pc anchors genuinely do cross it** and the ladder there is **graded, not binary**. *Their anchor choice is
  sound; it is the other 22 pulsars that RING indicts.*
- **RING's prior mean is ALWAYS exactly the true distance** (S7.2.7) — a "bad prior" there is **wide**, never
  **mis-centred**. **Wen et al.'s `D_p = D_p,0 + D_err·D_prior` has the same property.** So both share the
  optimism, and *"a mis-centred prior can only be worse."*
- RING used 30 ring pulsars and feather σ_d, not their array.

**Their gaps, against our record:**
- **Bias analysis: effectively ABSENT.** They vary distance *precision* (0.2–10.0 pc) but **never test a
  wrong or offset anchor**. Their prior is centred at truth by construction. **Our RING S-5 names the
  ready-to-run fix we also have not run:** the canonical `best_distances.txt` means differ from the feather
  means by **1.40σ on J0437 — a 0.55 rad pulsar-term phase error** — *"a ready-to-run mis-centred-prior arm
  and the single highest-value follow-up RING points to."* **Neither they nor we have run it. Say so.**
- **Null / FAR / trials: ABSENT.** Parameter recovery on injected signals only; detection assumed. Same delta
  as P2(a)/(b).
- **The anchor overlap is a genuine convergence worth citing:** their anchor set includes **J0437−4715**, and
  our Anchor Census independently finds **J0437 is the array's smallest-K pulsar** (K_lit = 3.07, the sole
  K ≤ 3, via the Reardon+2016 composite) — and **GEO found it certifies in 32/40 sky redraws**, more often
  than any census name, *omitted from our census triple by bad luck alone* (S3.3.3). **Two independent
  routes to the same pulsar.**

**Also surfaced, not scanned** — flagging only: **arXiv:2606.28721**, *"VLBI-Enabled Precision Localization of
CGW Sources with PTAs in the SKA Era"* (Jun 2026), same problem space, plausibly the same group. **Directly
adjacent to RING's VLBI-tier result (S7.2) and to the +2.9 ± 1.0 nat VLBI floor cost (S7.2.8). Recommend a
scan before the paper's related-work is frozen.**

---

## SUMMARY — the three things to take away

1. **The mechanism agreement is real, and their v3 strengthens us more than it threatens us.** Supp III —
   added under referee pressure, absent from their abstract — concludes that degeneracy-breaking **requires
   the chirp-mass-dominated regime**. *That is our frontier statement, reached independently from the
   opposite side.* **Cite it generously; it is the best card in the deck precisely because it is theirs.**
2. **The THEIRS-ONLY column is non-empty and matters.** The **pair-2D rank-2 object** (we hold rank-1 and an
   unbuilt rank-N); the **Φ_p-free reparameterisation** (structurally immune to the sub-fringe-offset artefact
   our banked E-step carries); the **analytic R_sys band-slope diagnostic** (the closed form of what GEO
   measured empirically); and **delivered SKA-scale results** (we have never touched a real TOA). **Our
   "containment" is specification, not delivery — SAMPLER is PENDING, and a referee will check.**
3. **The differentiator is the surrounding structure, and every piece of it is measured.** They have no null,
   no trials factor, no detection statistic that can return zero, no wrong-association rate, no
   accuracy/coverage statistic of any kind, a **min-over-84 order statistic on their headline precision**
   (12× pairing spread, their own Supp VI), and a source count (N = 4–5 resolved) that is **5–50× above
   SCOUT's SKA-era population clock**. *We have all of those, and the record of what they cost us.*
