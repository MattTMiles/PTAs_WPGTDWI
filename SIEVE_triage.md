# SIEVE_triage — the deterministic half (agent SIEVE-C, cronus)

**Session 2026-07-29. Venue: cronus, CPU only. `git pull --ff-only` clean (up-to-date at
`5461142`). Nothing committed, nothing pushed, no tracked file edited — new files staged
only; Matt commits.**

**HARD CONSTRAINT HONOURED.** cronus is RECOMPUTED-UNSAFE for noisy realisations.
**Every number in this report is zero-noise / Asimov.** No `rng` call touches a datum in
any of the three new compute modules; the two RNG calls that exist draw *nuisance
parameters* (T4's algebra gate) and *scrambled sky directions* (T4's manufacture ensemble),
neither of which is a noise realisation. The likelihood is the exact quadratic mismatch
`lnL = -0.5 Σ_p ||r_p(θ) − r_p(truth)||² / σ_p²`, so `lnL(truth) ≡ 0` and every ΔlnL quoted
is pure model-mismatch cost. The one test that genuinely needs fresh noise —
calibration-against-truth for T5 — is **deferred to SIEVE-A** and is not attempted here.

**BUDGET LAW — no cap approached.** T1 cap 2 GPU-hr, T5 cap 1 GPU-hr, T4 "CPU + light
4090": **0.00 GPU-hr spent on all four tests.** All three tests reduced to float64 numpy on
16×360 arrays (~18 min wall total) once the right structure was exploited — the pulsar-term
residual of pulsar *p* depends on `L_p` alone, so the distance block is **exactly separable**
and no 23-D optimiser is ever needed. This is reported as a *method* finding: the toy-scale
SIEVE questions are not GPU questions.

**THE TOY SPINE** (`hpc_harbor/sieve/sieve_toy.py`). 16 pulsars, top of the census by noise
rank (median TOA error ascending, restricted to pulsars holding both a canonical distance
prior and a tracked feather); **real** sky positions from `pos_t` in the tracked
`data_products/*.feather`; **real** distances and prior widths from
`CW_transition/best_distances.txt` (the canonical prior file). T = 18 yr, 20 epochs/yr.
*Toy choice, declared:* the measured per-pulsar spans are heterogeneous (6.5–20 yr); the toy
imposes a **uniform** T and cadence on every pulsar. That is a statement about the toy, not
about the real array's sampling.

Spine: J1909−3744, J1713+0747, J2234+0944, J0437−4715, J2241−5236, J1216−6410, J0125−2327,
J0613−0200, J1640+2224, J1614−2230, J2017+0603, J1811−2405, J1843−1113, J1853+1303,
J2010−1323, J0030+0451 (σ 7.4e−8 – 2.5e−7 s; L 0.156 – 2.499 kpc).

**A POINTER THAT DOES NOT RESOLVE (reported, not worked around).** The brief cites
`GLACIER_capstone.md §S4.24`. **No such section exists** — `GLACIER_capstone.md` (1154 lines)
runs to §S4.14 and the string `S4.2[0-9]` does not occur anywhere in the tracked `.md` set.
The "wander widths" the brief asks T5 to use are **`STORY.md §S4.2.4`** (+ §S4.2.3 for the
per-axis σ), which I have used and named explicitly. Likewise no file named as "the criterion
spec" exists; I took the criterion convention from `reports/GENERALISE.md`'s header
(criterion-v2.2: Gumbel valid only at zero-fraction ≤ 20 %, else empirical q95 ± bootstrap;
floors refit per cell). **None of T1/T4/T5/T6a is a floor-vs-count test, so criterion-v2.2 is
context here, not machinery.**

---

## T1 — WAVEFRONT CURVATURE

> ### VERDICT: **NEGLIGIBLE-EVERYWHERE** for the pipeline as it actually stands. The BEDROCK curvature rung is **KILLED**. Curvature is not a small effect that happens to be small — it is *exactly degenerate*, to Fresnel order, with the per-pulsar distance nuisance the programme already carries. It re-enters only if pulsar distances are ever pinned tighter than **7×10⁻⁴ of their canonical prior width**, which no pulsar in the census is.

**The mechanism, and it is an identity, not a coincidence.** Curvature and a pulsar-distance
error act on the observables through the *same single scalar* — the pulsar-term retardation
`dt_p`. Expanding the exact source-to-pulsar path
`dt_p = (|D ŝ − L_p p̂| + L_p − D)/c` to leading order in `L/D` gives the brief's Fresnel
phase `π L_⊥²/(λD)`, and that phase is reproduced **identically** by shifting the pulsar
distance by

```
    dL_eff,p  =  L_p² (1 + cos μ_p) / (2 D)
```

Because *any* observable that depends on the curvature only through `dt_p` cannot separate
it from a distance error, the degeneracy is exact — and it survives chirp, because a
frequency-evolving source still enters the pulsar term only via the retarded epoch.

**G0, the identity gate (measured).** Injecting the exact spherical front and scoring it
against a flat template with distances shifted by `dL_eff` costs, at SNR ≈ 97:

| D | SNR | max \|dL_eff\|/fringe | ΔlnL(sphere vs dL_eff-shifted flat) | ΔlnL(sphere vs Fresnel) |
|---|---|---|---|---|
| 16 Mpc | 96.7 | 0.0846 | **1.35e−07** | 1.35e−07 |
| 50 Mpc | 98.0 | 0.0271 | **1.42e−09** | 1.42e−09 |
| 200 Mpc | 98.4 | 0.0068 | **1.37e−11** | 1.37e−11 |

The two right-hand columns agree to every digit printed: the 1e−7 nat residual **is** the
Fresnel truncation error (the next order in `L/D`), not a failure of the identity. Against
`SNR²/2 ≈ 4.7e3` nat of signal, the identity holds to 3×10⁻¹¹ in relative terms.

**(a) + (b) — recovering curved data with FLAT templates.** ΔlnL cost of flatness and the
recovered source bias, in three distance regimes:

| D | regime | ΔlnL cost | Δsky | Δlog₁₀h | Δlog₁₀f | max\|ΔL\| |
|---|---|---|---|---|---|---|
| 16 Mpc | pinned | **1.9206** | 0.0010° | +6.1e−3 | +7.5e−6 | — |
| 16 Mpc | prior | **9.50e−3** | 0.0009° | +2.4e−3 | +8.6e−6 | 1.90 pc |
| 16 Mpc | free | **1.25e−3** | 0.0009° | +8.1e−4 | +8.6e−6 | 15.2 pc |
| 50 Mpc | pinned | 0.2020 | 0.0003° | +2.3e−3 | +2.4e−6 | — |
| 50 Mpc | prior | 7.65e−4 | 0.0003° | +6.7e−4 | +2.8e−6 | 0.013 pc |
| 50 Mpc | free | 8.18e−5 | 0.0003° | +1.9e−4 | +2.8e−6 | 15.2 pc |
| 200 Mpc | pinned | 0.01274 | 0.0001° | +6.0e−4 | +6.1e−7 | — |
| 200 Mpc | prior | 4.45e−5 | 0.0001° | +1.6e−4 | +6.9e−7 | 0.003 pc |
| 200 Mpc | free | 4.59e−6 | 0.0001° | +4.4e−5 | +6.9e−7 | 15.2 pc |

The pinned column scales as `D⁻²` exactly (1.9206/0.2020 = 9.51 vs (50/16)² = 9.77;
0.2020/0.01274 = 15.9 vs (200/50)² = 16.0), which is the Fresnel scaling and confirms the
fit is converged rather than optimiser-limited. **Even with distances pinned at truth,
curvature buys > 1 nat only inside D ≲ 22 Mpc at SNR ≈ 97** (cost = 1.92 × (16 Mpc/D)²,
∝ SNR²). Letting distances move under their *canonical* priors suppresses the 16 Mpc cost
**202×** to 9.5e−3 nat; unconstrained, **1540×** to 1.2e−3 nat. **Sky bias never exceeds
0.0010°** in any cell — three orders below the 0.188°/src break-even sky prior that the
R-posterior referendum identified, so curvature does not move that trade at all.

**Why it can never bite, per pulsar (16 Mpc, the worst case):**

| pulsar | L (kpc) | fringe (pc) | σ_L (pc) | dL_eff (pc) | cycles | **dL_eff/σ_L** | σ_L/fringe |
|---|---|---|---|---|---|---|---|
| J1909−3744 | 1.136 | 0.497 | 25.8 | 0.0018 | 0.0037 | 7.1e−05 | 51.9 |
| J1713+0747 | 1.053 | 0.554 | 60.9 | 0.0086 | 0.0154 | 1.4e−04 | 109.9 |
| J2234+0944 | 1.921 | 0.632 | 309.0 | 0.0535 | 0.0846 | 1.7e−04 | 488.6 |
| **J0437−4715** | 0.156 | 1.266 | **1.30** | 0.0009 | 0.0007 | **7.2e−04** | **1.03** |
| J0125−2327 | 1.348 | 0.880 | 91.2 | 0.0509 | 0.0579 | 5.6e−04 | 103.6 |
| J0030+0451 | 0.329 | 0.843 | 5.40 | 0.0029 | 0.0034 | 5.3e−04 | 6.41 |

**max dL_eff/σ_L over the whole spine = 7.2×10⁻⁴.** The curvature offset is at worst 0.07 %
of the tightest canonical distance prior in the census, and **every** prior spans ≥ 1 fringe
(min σ_L/fringe = 1.03, J0437 — the array's usual suspect, again the only pulsar anywhere
near a bar). Curvature is buried inside a nuisance the pipeline already marginalises.

**(c) — the GW-parallax check at 16 Mpc, curvature IN the template.** Profiling lnL over D:

| regime | peak | lnL(16 Mpc) − lnL(10⁵ Mpc) |
|---|---|---|
| pinned | **16 Mpc** | **+1.917 nat** |
| prior | 10⁵ Mpc | −1.05e−3 nat |
| free | 10³ Mpc | −5.62e−4 nat |

**GW parallax is real but unreachable.** With distances pinned at truth the profile peaks at
the true 16 Mpc and D is measured at 1.92 nat — the curvature phase *is* physically
recoverable, and this is the one place the effect exists. The moment distances are given
their canonical priors (let alone freedom) the profile is flat to ~10⁻³ nat and the truth is
not even favoured. The reason is structural: D enters through **one** parameter but its
signature is a per-pulsar pattern `∝ L_p²(1+cos μ_p)` living in a **16-dimensional** distance
nuisance space whose priors are 10³–10⁶× wider than the pattern. One parameter cannot be
separated from sixteen that already swamp it.

**KILL.** No promotion. The BEDROCK curvature rung should not be built.
*Scope lines that travel:* single non-evolving circular source, f = 10 nHz (λ = 0.972 pc),
SNR ≈ 97, 16-pulsar spine, local fits (no fringe hopping in the flat-recovery arm). Costs
scale as SNR² and as D⁻², so a source 10× louder at 16 Mpc reaches ~192 nat *pinned* — but
the pinned column is not the pipeline's regime, and the prior/free columns scale identically,
leaving the 202×/1540× suppression intact.

---

## T4 — CLOSURE OBSERVABLES

> ### VERDICT: **CLOSURE KILLS THE MANUFACTURE CLASS — but it is not the only thing that does.** The C3 double difference cancels the per-source phase and the per-pulsar achromatic clock exactly (to machine precision), and it drops the wrong-sky manufacture rate from **2/16 to 0/16**. However an explicitly modelled per-pulsar phase nuisance achieves the identical 0/16 at a statistically indistinguishable cost. Closure's value is **robustness and model-freedom, not new information** — and it has a real hole: it does **not** cancel a per-pulsar *chromatic* clock.

**Toy:** 6 pulsars (the spine's top 6) × 3 sources at f = {10.0, 15.1, 21.9} nHz, real
positions/distances, zero noise.

**PART 1 — the algebra, verified numerically.** Observable
`φ_pk = g_pk(sky, L_p, f_k) + Φ_k + c_p [+ 2π f_k τ_p]`, where `Φ_k` is the per-source phase
offset, `c_p` a per-pulsar **achromatic** phase offset, and `τ_p` a per-pulsar **time**
offset (a genuine clock error → chromatic phase `2π f_k τ_p`). Max deviation of each closure
class from its nuisance-free value:

| nuisance injected | C1 source-pair | C2 pulsar-pair | **C3 closure** | C4 chromatic |
|---|---|---|---|---|
| Φ_k only | 4.222e+00 | 3.64e−12 | **3.64e−12** | 3.66e−04 |
| c_p only | 3.64e−12 | 4.116e+00 | **7.28e−12** | 1.397e+08 |
| Φ_k + c_p | 4.222e+00 | 4.116e+00 | **3.64e−12** | 1.397e+08 |
| τ_p only | 3.117e−01 | 1.314e+00 | **4.048e−01** | 4.88e−04 |
| Φ_k + c_p + τ_p | 4.534e+00 | 4.204e+00 | **4.048e−01** | 1.397e+08 |

The geometric phases here are O(10⁴) rad and the C4 column is in units of seconds (O(10⁸)),
so **every "zero" above is machine precision in relative terms** (3.6e−12/1e4 ≈ 4e−16;
3.7e−4/1.4e8 ≈ 3e−12). Reading the kernels off the table:

- **C1** (source-pair, fixed p) kills `c_p`, keeps `Φ_k`.
- **C2** (pulsar-pair, fixed k) kills `Φ_k`, keeps `c_p`.
- **C3** (double difference) kills **both** — as advertised.
- **C3 does NOT kill a chromatic clock**: `τ_p` leaves a 0.405 rad residual. This is a
  *finding*, not a footnote — "per-pulsar clock-like term" is two different nuisances and
  closure only kills one of them. A real timing-clock error is chromatic.
- **C4** (the same double difference formed on `φ/f` rather than `φ`) kills `τ_p` **and**
  `Φ_k`, but **not** `c_p`. The two closure flavours have *complementary*, not nested,
  kernels; neither alone covers both clock classes.

**PART 2 — the money test (the STOP-#2 wrong-sky manufacture mode).** Data at truth; the
recovery template's source-0 sky scrambled isotropically (16 deterministic draws). Because
the observables are phases, everything runs in the **wrapped** (mod 2π) domain with pulsar
distances free over their **full canonical ±3σ prior**, which spans **8.8 to 2087 fringes**
per pulsar — i.e. **fringe hopping is allowed**, as it is in the real problem. This matters:
a first pass restricted to ±1 fringe rejected every scramble trivially and measured nothing.
Three statistics, identical freedom:

- `raw_blind` — wrapped residual with only `Φ_k` profiled out (**the pipeline's actual
  freedom — it carries no per-pulsar phase term**)
- `raw_nuis` — both `Φ_k` and `c_p` profiled out (the fair comparison)
- `closure` — C3, immune to both by construction

| arm | statistic | truth | scram min | scram median | **manufacture (scram beats truth)** |
|---|---|---|---|---|---|
| clean (no clock) | raw_blind | 0.0000 | 1.9474 | 18.288 | 0/16 |
| clean | raw_nuis | 0.0000 | 0.0366 | 0.1920 | 0/16 |
| clean | closure | 0.0000 | 0.0331 | 1.2458 | 0/16 |
| **clock injected** | **raw_blind** | **6.0990** | **4.9558** | 14.938 | **2/16** |
| clock injected | raw_nuis | 1.50e−13 | 0.0366 | 0.1920 | **0/16** |
| clock injected | **closure** | **0.0000** | 0.0331 | 1.2458 | **0/16** |

**The manufacture mode fires exactly when a per-pulsar systematic is present and unmodelled,
and only then.** In the clean arm all three statistics put truth on top and no scrambled sky
comes near — closure buys nothing. Inject an unmodelled achromatic per-pulsar clock and the
c_p-blind statistic is displaced from truth (χ² = 6.099 instead of 0) by enough that **2 of
16 scrambled skies fit the data better than the true sky** — a 12.5 % wrong-sky manufacture
rate, manufactured *not* by the sky freedom but by a systematic the model cannot represent.
Closure returns truth to exactly 0.0 and rejects all 16.

**But `raw_nuis` does the same job (0/16, truth 1.5e−13).** And the two agree numerically
where it counts: identical scram minima to four digits (0.0366 vs 0.0331 are different
statistics on different normalisations; the *ordering and the rejection* match). That is not
a coincidence — the closure is the *nuisance-profiled* raw statistic. Both remove 8 of the
18 phase dof and leave 4 after the 6 distances, so they are information-equivalent by
construction. **Closure is therefore a reparametrisation, not a new lever.**

**PROMOTE — narrowly, as an audit instrument, not as a detection statistic.**
→ *Destination:* a **C3+C4 closure-consistency diagnostic** run beside the existing
certification layer (never in place of it), whose job is to answer "is there an unmodelled
per-pulsar phase systematic?" — the one question that provably manufactures wrong-sky
certifications. Its value over `raw_nuis` is that it needs no nuisance model and costs 10
numbers instead of a 9-parameter profile; its cost is that C3 alone is blind to chromatic
clocks, so **both** flavours must be carried.
*Scope lines:* 6×3 zero-noise, phase-domain, one source's sky scrambled (not all three);
n = 16 scrambles, so 2/16 carries a ±1-event counting error and is a *demonstration that the
mode exists*, not a rate. **The rate itself, and its behaviour under noise, is SIEVE-A's.**

---

## T5 — SIGMA-POINT E-STEP, GATE HALF

> ### VERDICT: **GATE PASSES.** The sigma-point E-step reproduces the point E-step bit-comparably as the belief width → 0 (max TV = **1.11e−16** at exactly zero width, **7.92e−11** at 1e−8, converging as width² once the perturbation is sub-radian). At realistic-class source uncertainty the fringe posterior shifts by TV median **0.069–0.138** and the MAP fringe flips for **15 of 16** pulsars — but the underlying posterior is **near-uniform** (p(true fringe) 0.081 vs 0.059 uniform), so the flips are cheap and the *decision-relevant* shift is small. **PROMOTE the machinery, KILL the expectation that it rescues anything.**

**Why the toy had to chirp.** §S4.2.4 says 95–100 % of the pulsar-term phase wander is
`log10_mc`; a non-evolving toy has no `mc` at all. T5 therefore runs a GR-driven circular
binary with the exact closed form `f(t) = [f₀^(−8/3) − k t]^(−3/8)`, evaluated at the Earth
time and at the retarded epoch. Retardations are **171–16 010 yr**, so the pulsar term is a
genuinely earlier, lower-frequency epoch of the same binary. SNR = 128.5, 16 pulsars.

**Cross-check against §S4.2.4 (independent reproduction).** Per-dex gradients are
`|∂φ/∂log₁₀mc| = 0.94–6209 rad` and `|∂φ/∂log₁₀f| = 778–58 920 rad` — i.e. **per dex the
wander is dominated by f, not mc**. The banked "95–100 % is mc" is a *Fisher-weighted*
statement (gradient × that axis's posterior σ), and once weighted by σ(log₁₀Mc) = 0.501 dex
against the Earth-term-pinned σ(log₁₀f) = 0.144/389 dex, the mc share of the wander variance
is **72.8 – 100.0 %**, bracketing the banked 95–100 %. The low end (72.8 %) is a
different-pulsar-subset / fixed-Mc = 10⁹ M⊙ effect, not a contradiction. **The banked
statement is confirmed and its normalisation is now explicit: it is per-σ, not per-dex.**

**Sub-fringe artefact guard — the banked hazard is large here.** Scoring each fringe cell at
its **profiled** sub-fringe offset (ladder anchored on the *prior mean*, never on the
injected truth) versus the centre-pinned grid changes the fringe posterior by
**TV = 0.144 – 0.941, median 0.718.** That is a near-total rewrite of the posterior, and it
is a direct, independent confirmation of the banked sub-fringe-offset artefact at toy scale.
**Every T5 number here uses the guarded version.**

**GATE — width → 0:**

| belief-width scale | max TV vs point E-step | max \|Δw\| |
|---|---|---|
| 1e−0 | 1.667e−01 | 1.666e−01 |
| 1e−2 | 8.074e−02 | 4.039e−02 |
| 1e−4 | 3.784e−02 | 1.923e−02 |
| 1e−6 | 7.925e−07 | 3.982e−07 |
| 1e−8 | 7.925e−11 | 4.054e−11 |
| **0** | **1.110e−16** | **5.551e−17** |

**PASS.** At exactly zero width the sigma points collapse and the two E-steps agree to
machine epsilon. The 1e−6 → 1e−8 leg falls as width² — correct unscented behaviour. The
non-quadratic middle (1e−2, 1e−4) is *physical, not a gate failure*: at scale 1e−2 the
effective σ(log₁₀mc) = 5e−3 dex still drives up to 31 rad ≈ 5 fringes of pulsar-term wander,
so the perturbation is not small until scale ≲ 1e−6. Quoting this leg as a failure would be
a misreading; the gate is the width → 0 limit, and it passes.

**READOUT — the shift at §S4.2.3/§S4.2.4 widths:**

| σ(log₁₀Mc) | TV median (min–max) | MAP-fringe flips | mean p(true fringe) point → σ-point | wander at this width |
|---|---|---|---|---|
| 0.501 dex | 0.0689 (0.0021–0.1667) | 15/16 | 0.0814 → 0.0751 | 0.47–3111 rad = 0.075–495 fringes |
| 0.864 dex | 0.1352 (0.0086–0.2122) | 15/16 | 0.0814 → 0.0742 | 0.81–5360 rad = 0.13–854 fringes |
| 0.866 dex | 0.1376 (0.0066–0.2123) | 15/16 | 0.0814 → 0.0734 | 0.81–5380 rad = 0.13–856 fringes |

**Read this carefully, because the two headline numbers point opposite ways.** The MAP
fringe flips for 15 of 16 pulsars at every realistic width — which sounds decisive — but the
guarded fringe posterior is **already almost flat**: p(true fringe) = 0.0814 against a
uniform 1/17 = 0.0588. The argmax of a near-uniform distribution is not a meaningful
quantity, and the flips cost almost nothing. The honest metric is TV (0.069–0.138) and the
change in p(true fringe) (−8 % relative at the widest width), both **modest**. The
mechanism: at σ(log₁₀Mc) = 0.501 dex the pulsar-term wander is **up to 495 fringes**, so
averaging over the belief does not blur the posterior *between* neighbouring fringes — it
scrambles the ladder wholesale, and a wholesale scramble of an already-flat ladder returns
something flat. **Sigma-point averaging is the correct E-step, and it is correct for the
reason that makes it nearly pointless here: the fringe posterior carries so little
information at realistic source uncertainty that there is very little for the averaging to
change.** This is fully consistent with the banked D3 result that a single CW resolves
fringes for zero pulsars.

**PROMOTE (machinery), KILL (the hope).**
→ *Destination:* land the 5-point unscented E-step in the M-step/E-step driver as the
**default** — it is gate-clean, costs 5× one E-step, and removes an approximation that is
indefensible at σ(log₁₀Mc) ≈ 0.5–0.87 dex. Do **not** budget for it to move certified
counts; the readout says it will not.
**HANDOFF TO SIEVE-A (explicit):** calibration-against-truth — does the sigma-point E-step
produce *better-calibrated* fringe posteriors than the point E-step over noisy realisations?
— is **not attempted here** and cannot be, on banked data or otherwise, without noise draws
this box cannot make safely. SIEVE-A owns it. What SIEVE-A should inherit: the guarded
fringe ladder (prior-mean-anchored, sub-fringe-profiled), the 5-point unscented weights, and
the finding that TV is small *because the posterior is flat* — so the calibration test must
be powered against a near-uniform null, not against a peaked one.

---

## T6a — CRN AUDIT (code read, zero compute)

> ### VERDICT: **THE ARMS SHARE SEEDS. Both contrasts are fully PAIRED (common random numbers) by construction — census seed, noise seed, distance-truth seed, and scramble seed are all arm-independent.** This is deliberate and documented, not accidental. The paired-vs-independent variance demo is therefore a measurement of *how much* the pairing buys, not a check of *whether* it exists.

The design is stated in the module docstring itself —
`hpc_harbor/glacier/glacier_loop.py:10-12`: *"TWO IGNITERS (Matt's budget decision, Option
1), **PAIRED on census + noise seeds** (the CHORUS C6 discipline — the igniter contrast is
exact)"* — and the seed convention at `glacier_loop.py:69-70`: *"census 9_0xx_xxx (9_000_000
+ sky); noise 9_5xx_xxx (9_500_000 + 1000*sky + real); dist = noise + 10_000_000; scrambles
9_8xx_xxx."* **No term in any of those formulas is a function of the arm.**

**Line cites — every seed, and what it does not depend on:**

| seed | line | expression | arm-dependent? |
|---|---|---|---|
| census (Stage-1) | `glacier_loop.py:394` | `draw_population(SEED_POP_BASE + sky, ...)` | **no** |
| noise | `glacier_loop.py:421` | `noise_seed = SEED_NOISE_BASE + 1000*sky + real` | **no** |
| distance truth | `glacier_loop.py:422` | `dist_seed = noise_seed + 10_000_000` | **no** |
| scramble | `glacier_loop.py:443` | `default_rng(SEED_SCRAM_BASE + 1000*sky + real)` | **no** |
| bases | `glacier_loop.py:123-124` | `SEED_NOISE_BASE = 9_500_000`, `SEED_SCRAM_BASE = 9_800_000` | — |

**e07-vs-none.** The arm enters at `glacier_loop.py:404`,
`embed_igniter(pop, E_IGNITER[arm], VENUE_SPAN_S[T])` — it changes the *source slot list*
(e07 appends 31 harmonic slots; `none` is the same census all-circular,
`glacier_loop.py:120`), never a seed. The noise vector is therefore **bit-identical across
arms**, and this is provable rather than merely plausible: `NoiseDrawer.draw`
(`CW_transition/trackB_b1_core.py:592-613`) consumes its RNG stream as `c_gwb` (skipped
entirely when `"gwb"` is not in `components` — and GLACIER passes
`components=("white","rn")`, `glacier_loop.py:437`), then per-pulsar white, then per-pulsar
red. **The stream depends only on the pulsar set and TOA layout, never on the source list.**
Same pulsars + same T + same seed ⇒ same draw, whatever the arm did to the census.

**frozen-vs-live.** `mstep_policy` is a *venue axis of the output stem only* —
`glacier_loop.py:387-388` appends `_m<policy>` to the stem, with the explicit comment that
`"live"` appends nothing so *"incumbent stems"* stay bit-for-bit. The belief machinery is
constructed at `glacier_loop.py:466-471`, i.e. **after** the noise draw at line 437, so it
cannot perturb the RNG stream even incidentally. The seeds above are untouched by it.

> **CONCURRENCY, recorded (the §g2.4 hazard, live again).** `hpc_harbor/glacier/glacier_loop.py`
> **changed while this audit was running** — a concurrent session is staging FORGE-B belief
> work on this tree (`+177/−11` on `glacier_loop.py`, plus `CW_transition/trackB_b1_belief.py`,
> `reports/FORGE_B_belief.md`, `CW_transition/logs/forge_b_belief_gates.log`, and an `M` on
> `CW_transition/trackB_b1_softsource.py`, all of which appeared **in the index** and are
> **not mine**). I re-verified every T6a cite against the file as it now stands: lines
> 120/123/124/387/394/404/421/422/437/443 are **unchanged and still say what is quoted above**;
> the belief-book cite moved 456 → 466 and is corrected here. **The CRN verdict is unaffected** —
> no seed expression was touched. I did not stage, revert, or edit any of the concurrent
> session's files, and I did not race it.

**Two facts that travel beside the pairing, both from the same read:**

1. **Signal and scrambled-null cells at the same `(sky, real)` share the noise draw too**, since
   `noise_seed` does not depend on `scrambled` — and the scramble is applied to the *recovery
   template only*, never the data (`glacier_loop.py:440`, comment; `442-445`). The null is
   therefore *the same data with a scrambled recovery sky*, which is the strongest possible
   form of the control. Worth stating explicitly because it means the null and signal arms
   are paired as well.
2. **The pairing is exact within a `(T, wscale)` cell and NOT across T.** `white_scale`
   preserves the stream by construction (`trackB_b1_core.py:593-598`: the same normals are
   drawn at any scale, and `wscale = 1.0` is bit-identical). But across different `T` the
   same seed produces a *different-length* draw and hence a different realisation — so any
   T = 30 vs T = 40 comparison is **unpaired**, unlike every arm comparison. If SIEVE-A's
   variance demo touches the time ladder, that leg needs different bookkeeping.

**PROMOTE (as a measurement, already de-risked).**
→ *Destination:* SIEVE-A's paired-vs-independent variance demo can proceed directly on the
banked GLACIER arm pairs with **no plumbing change** — the CRN is already there and exact.
The demo's deliverable is the variance-reduction factor on ΔN and on the drain, per cell;
its scope must exclude the T axis per finding (2).

---

## SPEND

**0.00 GPU-hr** against caps of 2 (T1) + 1 (T5) + "light 4090" (T4). ~18 min of single-core
CPU wall across all three compute modules, plus one code read. No cap approached; no test
was trimmed, deferred, or reduced in scope for budget. The only deferral in this report is
the T5 noisy-calibration handoff, which is a *venue-safety* deferral (cronus is
RECOMPUTED-UNSAFE), not a budget one.

## ADD-LIST (staged with `git add`; **Matt commits — I have committed and pushed nothing**)

```
SIEVE_triage.md                              new   this report (T1/T4/T5/T6a)
hpc_harbor/sieve/sieve_toy.py                new   the zero-noise toy spine + CW model
hpc_harbor/sieve/sieve_t1_curvature.py       new   T1 wavefront curvature
hpc_harbor/sieve/sieve_t4_closure.py         new   T4 closure observables + money test
hpc_harbor/sieve/sieve_t5_sigmapoint.py      new   T5 sigma-point E-step, gate half
reports/sieve_t1_curvature.npz               new   T1 bank (gate, rows, D-profiles)
reports/sieve_t4_closure.npz                 new   T4 bank (algebra table, money ensemble)
reports/sieve_t5_sigmapoint.npz              new   T5 bank (gate, readout, wander shares)
```

No tracked file was edited. `hpc_harbor/glacier/glacier_loop.py` was **read only** for T6a.

**DO NOT DOUBLE-STAGE — the index is shared.** `git add` was run on the eight paths above and
**only** those. The index also currently holds, authored by a **concurrent session, not me**:
`hpc_harbor/glacier/glacier_loop.py` (M, +177/−11), `CW_transition/trackB_b1_belief.py`,
`CW_transition/trackB_b1_softsource.py` (M), `CW_transition/logs/forge_b_belief_gates.log`,
`reports/FORGE_B_belief.md`, `LOTTERY2.md`, `hpc_harbor/lottery/lottery2_realistic_pop.py`,
`reports/lottery2_realistic_pop.npz`, `reports/OBSERVERS_TABLE.md` (M). **Reconcile those with
their author before committing**; the remaining dirty working tree (CW_transition/*.npz,
SAMPLER_dev.md, notebooks, logs) is likewise not mine and not staged by me.

## STILL OWED / named (not dropped)

- **The brief's `§S4.24` pointer does not resolve** (§ above). I proceeded on `STORY.md
  §S4.2.4 + §S4.2.3` and named it everywhere. If that substitution is wrong, T5's readout
  widths change and the readout must be re-cut; the gate is unaffected.
- **T4's 2/16 is a demonstration, not a rate.** The wrong-sky manufacture rate under an
  unmodelled per-pulsar clock, and its noise behaviour, is SIEVE-A's.
- **T4's chromatic-clock hole (C3 leaves 0.405 rad on `τ_p`)** is unclosed by C3 alone. If
  the closure diagnostic is built, C4 must be built with it.
- **T5 calibration-against-truth** — SIEVE-A, on banked data. Power it against a
  near-uniform null.
- **T1's `free`-regime `max|ΔL| = 15.2 pc`** is a genuine walk along the exactly-degenerate
  direction inside the ±1-fringe local bracket, not a convergence failure — but it is
  bracket-limited, so the *unconstrained* suppression factor (1540×) is a lower bound on the
  true one.

---

# CLOSEOUT DELIVERABLES

## C1 — T1 paper text → *the draft's systematics discussion*

Drop-in paragraph. Numbers are the toy's (16 pulsars, 18 yr, SNR ≈ 97, `f = 10⁻⁸ Hz`);
re-cut against the real array before submission — the *scalings* and the identity are
array-independent, the absolute ΔlnL is not.

> **Wavefront curvature.** Continuous-wave searches model the source as a plane wave, but a
> source at finite distance `D` emits a spherical front, and the pulsar term samples it a
> transverse baseline `L_⊥ = L_p sin μ_p` off the Earth–source axis. Expanding the exact path
> `(|D ŝ − L_p p̂| + L_p − D)/c` to first order in `L_p/D` shows that the correction enters
> the likelihood *only* through the pulsar-term retardation `Δt_p`, in exactly the place the
> unknown pulsar distance already enters. Curvature is therefore not an independent
> systematic but an **exact reparametrisation of the per-pulsar distance nuisance**, with
> effective offset
>
> `dL_eff,p = L_p² (1 + cos μ_p) / (2 D)`,
>
> equivalently the Fresnel phase `π L_⊥²/(λ D)`. We verify the identity numerically: injecting
> a spherical front and recovering with a plane-wave template whose distances are shifted by
> `dL_eff,p` leaves a residual of 1.3 × 10⁻⁷ nat at `D` = 16 Mpc, falling as `D⁻²` to
> 1.4 × 10⁻¹¹ nat at 200 Mpc — identical, to all printed digits, to the truncation error of
> the first-order expansion itself, i.e. the shifted plane-wave template is the spherical
> template to the order we work at. Two consequences follow. First, the effect is
> **correctable at zero cost**: any pipeline that already fits pulsar distances absorbs it
> automatically. With distances pinned at their catalogue values the flatness penalty is
> 1.92 nat at 16 Mpc, but simply opening the distances under their standard Gaussian priors
> suppresses it by 202× to 9.5 × 10⁻³ nat, and the induced sky bias is ≤ 0.001° in every
> configuration tested. Second, the effect is **negligible for all but the nearest plausible
> source**: `dL_eff,p ≤ 0.053 pc` at 16 Mpc, which is 7 × 10⁻⁴ of the smallest distance
> uncertainty in the array and 0.08 of a fringe, and the uncorrected penalty falls below
> 1 nat beyond `D` ≈ 22 Mpc at this SNR (and further out still at lower SNR, since the penalty
> scales as `h² D⁻²`). The same identity forecloses the optimistic reading: because curvature
> acts *only* along the already-free distance direction, it carries no independent distance
> information. Profiling the likelihood over `D` with curvature in the template recovers the
> true 16 Mpc at +1.92 nat when distances are pinned, but that peak collapses to ≲ 3 × 10⁻³
> nat — flat, with the truth not even preferred — the moment the distances are allowed their
> priors. **There is no GW parallax to be had from a single continuous wave at these
> baselines.** We therefore adopt plane-wave templates throughout and quote the curvature
> systematic as bounded by the table below.

**One-line correction table** (16-pulsar toy, `f = 10⁻⁸ Hz`, SNR ≈ 97; `Δ` = flat-template
penalty, `resid` = what survives after applying `dL_eff,p`):

| `D` (Mpc) | max `dL_eff` (pc) | median | max `dL_eff/σ_L` | max cycles | ΔlnL pinned | ΔlnL prior-free | resid after correction |
|---|---|---|---|---|---|---|---|
| 16 | 0.0535 | 0.0077 | 7.2e−4 | 0.085 | 1.92 | 9.50e−3 | 1.35e−7 |
| 50 | 0.0171 | 0.0025 | 2.3e−4 | 0.027 | 2.02e−1 | 7.65e−4 | 1.42e−9 |
| 200 | 0.0043 | 0.0006 | 5.8e−5 | 0.0068 | 1.27e−2 | 4.45e−5 | 1.37e−11 |

Every column scales as `D⁻¹` (offsets) or `D⁻²` (ΔlnL) to 2–3 %, which is the convergence
check as well as the physics. Source: `reports/sieve_t1_curvature.npz`,
`hpc_harbor/sieve/sieve_t1_curvature.py`.

## C2 — T5 erratum note → **flagged for the ACCRE fold; capstone NOT edited here**

> **Erratum / clarification, §S4.2.4 — "95–100 % of the pulsar-term phase wander is
> `log10_mc`".** The statement is correct but **normalisation-dependent, and the normalisation
> is not stated in the text.** It is a *per-σ* (Fisher-weighted) decomposition: each axis's
> phase gradient multiplied by that axis's posterior width. Under that weighting the SIEVE-C
> toy independently reproduces an `mc` share of the wander variance of **72.8 – 100.0 %**,
> bracketing the banked 95–100 % (the low tail is a different-pulsar-subset and fixed
> `Mc = 10⁹ M⊙` effect, not a disagreement). Under the *per-dex* normalisation the ranking
> **inverts**: bare gradients are `|∂φ/∂log₁₀mc| = 0.94 – 6.21e3` rad against
> `|∂φ/∂log₁₀f| = 778 – 5.89e4` rad, giving an `mc` share of only **1.5e−6 – 1.1e−2**, i.e.
> per dex the wander is overwhelmingly **`f`-dominated**. Both statements are true of the same
> system; they differ only in whether the axes are compared per dex or per posterior sigma.
> The `mc` dominance in §S4.2.4 is bought entirely by `σ(log₁₀Mc) ≈ 0.5 dex` against an
> Earth-term-pinned `σ(log₁₀f) = 0.144/389 dex ≈ 3.7e−4 dex`. **Recommended fix:** add "per
> posterior σ" to the §S4.2.4 sentence, and note that the conclusion does not survive a
> reader who assumes per-dex. Anything downstream that inherits the "`mc` dominates" claim
> while using an *unpinned* or wider `f` posterior must re-derive the share — the crossover
> is at `σ(log₁₀f)/σ(log₁₀Mc) ≈ 10⁻²`, not far below the pinned value.
>
> Source: `reports/sieve_t5_sigmapoint.npz` keys `mc_share` (per-σ) and `mc_share_perdex`,
> `hpc_harbor/sieve/sieve_t5_sigmapoint.py`. Zero-noise toy; the *shares* are geometric and
> noise-free, so this does not defer to SIEVE-A.

## C3 — STAGE HYGIENE: files authored by SIEVE-C **only**

Per-author commits — the eight paths below are mine and complete. Everything else in the
index or the working tree belongs to FORGE-B or another concurrent session; reconcile with
their add-list separately.

```
SIEVE_triage.md
hpc_harbor/sieve/sieve_toy.py
hpc_harbor/sieve/sieve_t1_curvature.py
hpc_harbor/sieve/sieve_t4_closure.py
hpc_harbor/sieve/sieve_t5_sigmapoint.py
reports/sieve_t1_curvature.npz
reports/sieve_t4_closure.npz
reports/sieve_t5_sigmapoint.npz
```

I edited no tracked file, committed nothing, pushed nothing, and staged nothing outside this
list.
