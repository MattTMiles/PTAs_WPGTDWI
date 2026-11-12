# VERNIER-6 — everything on, maximum optimism

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~8 min.
Clean room: imports only `vernier_toy` (my own anchored code).
**Primary metric is CONTRAST. Every cell is reported against a scrambled null.**

Code `hpc_harbor/vernier/vernier6.py` · data `reports/vernier6.npz` ·
plot `reports/VERNIER6_contrast.png`. Staged, not committed. **No realism added, as
instructed.**

---

## Lead — the loop closes

| | SIGNAL | SCRAMBLED NULL |
|---|---|---|
| contrast (median) | **136.2 nats** | −6.0 |
| contrast 16/84 | 122.9 / 152.9 | 95th: **+0.37** |
| q_top | **1.0000** | 0.447 |
| bar = max(ln K, 2) | 3.45 (K = 31.4) | — |
| **P(clear bar)** | **1.000** | **0.025** |
| **P(top = true fringe)** | **1.000** | — |

**Excess over null: 142.3 nats. Null false-alarm rate at the bar: 0.025.**

### Settings that closed it — recorded verbatim

```
n_psr = 50            T = 30 yr             sigma_toa = 20 ns      GWB = False
f_gw_top = 25 nHz     log10_mc = 9.5        L = 1 kpc              prior = 1 pc
rho_P = 5.0 per pulsar (LOUD, deliberately)  n_src = 8
ecc = 0.7             n_harm = 30           sigma_f_rel = 1e-4
FLOOR on the bar = 2.0 nats
```

The reframe is what did it. Under the old binary metric this configuration reports
"P(unique) < 1" and looks like a partial failure; under contrast it is 136 nats clear of the
best false fringe with the true fringe on top in 100% of trials.

---

## Ablation — ranked by contrast lost, matched trials

All arms use **the same seeds**, so they see the same sky and the same noise.

| variant | contrast | loss | pass bar | q_top | null 95th |
|---|---|---|---|---|---|
| **all on** | **143.4** | — | 1.00 | 1.000 | 0.83 |
| T5 seeding OFF (1 source) | 0.04 | **−143.4** | 0.19 | 0.643 | 0.00 |
| T2 joint/pairwise OFF | 38.9 | **−104.5** | 1.00 | 1.000 | 0.39 |
| T3 eccentric clock OFF | 99.2 | **−44.3** | 1.00 | 1.000 | **+30.4** |

| rank | technique | contribution |
|---|---|---|
| **1** | **T5 source multiplicity** | **−143.4 nats** — nothing works with one source |
| **2** | **T2 joint/pairwise** | **−104.5 nats** (at σ_f/f = 1e-4; see sensitivity) |
| **3** | **T3 eccentric clock** | **−44.3 nats**, *and* a 5× worse null — see below |
| **4** | **T1 posterior product** | **0.0** in power, **5× in safety** — see below |
| **5** | **T4 weight not gate** | **0.0 here**; matters only where q ∈ [0.5, 0.9] |

**T1 and T4 are estimators, not signal generators**, so they are scored on the *same*
trials rather than as contrast losses.

**T1 — posterior product vs survivor counting, identical data:**

| | P(success) | P(success) on the NULL |
|---|---|---|
| posterior product | 1.000 | **0.040** |
| survivor counting (VERNIER-1…5) | 1.000 | **0.200** |

At this optimistic point both estimators succeed on every trial, so T1 buys **zero power** —
which does *not* support the conjecture that a weak estimator explains VERNIER-5's
reproduce-failure. What T1 buys is **safety**: survivor counting false-alarms on 20% of
scrambled nulls, the product on 4%. **A 5× lower false-alarm rate for free.**

**T2 sensitivity** — its contribution is entirely set by the shared-parameter error:

| σ_f/f | 3e-5 | **1e-4** | 3e-4 | 1e-3 |
|---|---|---|---|---|
| contrast T2 on | 150.4 | **146.7** | 121.8 | 15.0 |
| contrast T2 off | 134.9 | **49.2** | 0.004 | 0.000 |
| **gain from T2** | 15.5 | **97.5** | 121.8 | 15.0 |

I chose 1e-4 as the fiducial so T2-on is nearly unpenalised and T2-off is degraded but not
annihilated. **At 1e-3 the T2-off envelope is e^−131 — the ablation would be a degenerate
limit, not a measurement.** T2's headline number is therefore a *choice*, and I have shown
the whole curve rather than one point.

**T4** — at this working point q = 1.000 for essentially every pulsar, so gate and weight
agree exactly (both 0.1929 dex). **T4 contributes nothing here by construction.** It can only
matter where q sits at 0.5–0.9, which is where a realistic array will sit and this one does
not. Deliberately optimistic settings hide T4's value.

---

## Scaling

**Contrast grows LINEARLY in N — log-likelihoods add, exactly as the addendum predicted.**

| N sources | 1 | 2 | 3 | 4 | 6 | 8 | 12 |
|---|---|---|---|---|---|---|---|
| contrast | 0.07 | 16.9 | 37.6 | 54.6 | 97.8 | 142.6 | **228.9** |
| P(clear bar) | 0.26 | 0.92 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

Linear fit **21.0 nats per source** (R² by inspection excellent); the log fit is visibly
worse. This is the single most important structural result: **under contrast the techniques
compound, whereas under the old binary metric each looked useless alone.**

**The bar grows only logarithmically with the prior — D2's saturation, mechanised.**

| prior | 0.1 pc | 0.3 | 1.0 | 3.0 | 10 | 100 |
|---|---|---|---|---|---|---|
| K fringes | 3.9 | 10.0 | 43.4 | 120 | 391 | **3188** |
| **bar = ln K** | 2.00 | 2.31 | 3.77 | 4.79 | 5.97 | **8.07** |
| contrast | 172.2 | 162.7 | 135.2 | 127.4 | 126.1 | 106.2 |

**Confirmed.** A **1000× wider prior costs 6.1 nats of bar** and about 66 nats of contrast
(the latter because a wider grid offers more chances for a high false peak — so contrast is
weakly falling, not flat as I stated in the code comment). Both effects are logarithmic and
both are dwarfed by the 21 nats/source that multiplicity supplies.

> **Tightening the distance prior buys you a logarithm. Adding a source buys you a
> constant, per source. That is why D2 saturated at 1 pc, and it says the prior was never
> the lever.**

## T3 — the eccentric clock removes regions, and that is a safety result

| | excluded volume | contrast | **null FPR** | null 95th |
|---|---|---|---|---|
| circular (n=2 only) | 0.994 | 97.2 | **0.150** | **+15.0** |
| eccentric e=0.7, n ≤ 30 | 0.999 | 137.2 | **0.033** | **+0.32** |

Contrast gain **+40.0 nats**; extra prior volume excluded only **+0.5 percentage points**.

**The answer to your question is "removes regions", but the evidence for it is in the null,
not in the excluded fraction.** Without eccentricity the scrambled null clears the bar **15%**
of the time with a 95th percentile of **+15 nats**; with harmonics to n = 30 it clears it
**3%** with a 95th of **+0.3 nats**. The harmonics' sinc envelopes narrow as 1/n and veto
stretches of the comb that noise would otherwise be free to populate.

**This is the ORACLE failure mode (q = 0.986 on a wrong fringe) being suppressed by
construction rather than by a threshold.** VERNIER-1 tested this handle for *circular*
sources (n = 2) and found it met in 0/60 band cells; at n = 30 the frequency shift is 15×
larger, which is why it turns on here.

---

## Two bugs found, both of which would have inverted conclusions

1. **Ablation arms were not RNG-matched.** The first run reported "T1 off" as *better* than
   "all on" — impossible. Contrast has a broad distribution (16/84 spanning 123–153), so
   independent draws per arm are dominated by trial scatter. Fixed with shared seeds.
2. **The scramble displaced each *component* independently.** A source's 30 harmonics share
   one true distance; scrambling them separately destroyed an eccentric source far more
   thoroughly than a circular one, making the T3 null unfair and producing a spurious
   **null95 = +60 nats**. Fixed to scramble at the *source* level. This is exactly the class
   of error the null exists to catch, and it was in the null itself.

## What this does not establish

1. **These settings are deliberately optimistic and several are unphysical in combination**:
   no GWB, sources pre-resolved with known sky positions, ρ_P = 5 *per pulsar* (VERNIER-4's
   realistic census gives ~0.1), and e = 0.7 with 30 usable harmonics. **This is a
   demonstration that the mechanism exists, not evidence that it operates.**
2. **T2's 104 nats is a choice of σ_f/f**, spanning 15→122 nats across a factor 30. The real
   number needs the actual source-parameter posterior.
3. **T1 and T4 both contributed zero power at this working point.** Their value is claimed
   for the realistic regime, which this run does not touch — so the ablation ranks them
   lowest *here* while their whole purpose is to matter *elsewhere*. That is a limitation of
   ablating at an optimistic point, not evidence they are useless.
4. **The null is a source-level scramble only.** It destroys cross-source consistency but
   preserves each source's comb; other failure modes (mis-specified Mc, wrong sky position)
   are not tested.
5. **Still no GWB, no timing-model marginalisation, no realistic pulsar noise, sources not
   found but given.** All four are the next dials, and per instruction none were touched.

## One-line verdict

**The loop closes — 136 nats of contrast against a null 95th of +0.4, true fringe on top in
100% of trials — and the ranking for the walk back is: source multiplicity (21 nats each,
compounding linearly) ≫ joint source-parameter treatment (choice-dependent, 15–122 nats) >
eccentric harmonics (44 nats *and* a 5× better null) > posterior product (no power, 5× better
null) > weight-not-gate (nothing here, everything in the regime this run avoided).**
