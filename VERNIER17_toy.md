# VERNIER-17 — turn the foreground dial

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~8 min.
Clean room: imports only `vernier6` / `vernier7` / `vernier8` / `vernier9`.

Data `reports/vernier17.npz`, `reports/vernier17_wheel.npz`. Staged, not committed.

---

## Lead — our threshold is 0.70, exactly the field's, and most of the census sits below it

**Where certification stops working: λ_fore = 0.70 by the contrast criterion, 0.80 by the
stricter per-pulsar certification bar.** The field's stated 0.7 for reliable recovery is
reproduced *exactly* — not looser, not stricter.

| band | λ_fore | ρ_eff | contrast | bar | **P(clear bar)** | null 95th | w_true sig/null |
|---|---|---|---|---|---|---|---|
| 2.00–2.85 nHz (156/bin) | 0.263 | 0.593 | **−1.30** | 3.69 | **0.00** | −0.01 | 0.125 / 0.068 |
| field: detection | 0.500 | 0.981 | 0.05 | 3.69 | **0.06** | 0.07 | 0.533 / 0.065 |
| **census median** | **0.580** | 1.144 | 1.29 | 3.69 | **0.24** | 0.12 | 0.833 / 0.064 |
| 2.85–3.71 (46/bin) | 0.579 | 1.142 | 1.27 | 3.69 | **0.22** | 0.12 | 0.830 / 0.064 |
| 4.56–5.41 (11/bin) | 0.636 | 1.278 | 2.67 | 3.69 | **0.38** | 0.17 | 0.958 / 0.064 |
| **field: reliable** | **0.700** | 1.461 | 5.07 | 3.69 | **0.68** | 0.24 | 0.997 / 0.063 |
| 3.71–4.56 (22/bin) | 0.722 | 1.534 | 6.19 | 3.69 | **0.76** | 0.27 | 0.999 / 0.062 |
| 5.41–6.27 (5/bin) | 0.811 | 1.914 | 12.99 | 3.69 | **0.94** | 0.47 | 1.000 / 0.057 |
| above 8 nHz | 0.999 | 4.939 | 137.5 | 3.69 | 1.00 | 1.15 | 1.000 / 0.008 |
| VERNIER-6 (λ=1) | 1.000 | 5.000 | 141.2 | 3.69 | 1.00 | 1.18 | 1.000 / 0.007 |

**Of PHOENIX's measured bands, only 5.41–6.27 nHz (λ = 0.811) works cleanly. 3.71–4.56
(0.722) is marginal. Everything at or below the census median 0.58 fails.**

## The ceiling, which is the mechanism

λ_fore = h²_CW/(h²_CW + h²_conf) means h²_conf = h²_CW(1−λ)/λ, so the residual foreground's
power scales *with the source's own power*. Treating it as an added noise component:

> ρ_eff² = ρ₀²/(1 + ρ₀²(1−λ)/λ), which for a loud source saturates at
> **ρ_eff → √(λ/(1−λ))**

**A hard ceiling independent of source brightness.** Brightness cannot buy past it, because
the confusion scales with the population that produced the source.

| λ | 0.263 | 0.50 | 0.58 | 0.70 | 0.722 | 0.811 |
|---|---|---|---|---|---|---|
| ceiling on ρ_P | 0.60 | 1.00 | 1.17 | 1.53 | 1.61 | **2.07** |

**Certification needs ρ_P ≥ 2, hence λ/(1−λ) ≥ 4, hence λ ≥ 0.80.** Only the 5.41–6.27 nHz
band clears that. The contrast criterion — clearing ln K rather than certifying every pulsar
— is more forgiving and turns on at 0.70.

**Two thresholds, and they are not the same test:** contrast > ln K at **λ ≥ 0.70** (the
vernier resolves the fringe), per-pulsar ρ_P ≥ 2 at **λ ≥ 0.80** (every pulsar's distance
certifies). The field's 0.7 corresponds to the first.

## The wheel

| band | λ | ρ_eff | **final rank** | w_true sig | w_true null |
|---|---|---|---|---|---|
| 2.00–2.85 nHz | 0.263 | 0.59 | **5** | 0.0080 | 0.0064 |
| census median | 0.580 | 1.14 | **5** | 0.0113 | 0.0063 |
| 4.56–5.41 | 0.636 | 1.28 | 5 | 0.0121 | 0.0063 |
| field 0.70 | 0.700 | 1.46 | 6 | 0.0146 | 0.0056 |
| 3.71–4.56 | 0.722 | 1.53 | 6 | 0.0224 | 0.0051 |
| 5.41–6.27 | 0.811 | 1.91 | **7** | 0.0206 | 0.0026 |
| above 8 nHz | 0.999 | 4.94 | **10** | 1.0000 | 0.0000 |

**Against the standing 12–13: the wheel drops to rank 5–7 across the real census, and only
recovers rank 10 in the λ → 1 bins above 8 nHz.** Note the w_true column: at λ ≤ 0.64 the
signal sits at 0.008–0.012 against a null of 0.006 — barely above chance, so those
"certifications" are not certifications at all. Only λ = 0.999 gives clean separation
(1.0000 vs 0.0000).

**The standing rank 12–13 was obtained at λ_fore = 1, which the banked cells also assume by
construction.** On PHOENIX's measured foreground it does not survive.

---

## The Gaussian approximation, and what it costs

**I approximated the residual foreground as an added Gaussian noise component.** It is not —
it is other CW sources: coherent, sky-localised, and in principle separable. That
approximation is **conservative in one direction and optimistic in another**, and I cannot
sign the net:

- **Conservative:** a real analysis can fit or subtract the brightest interlopers, recovering
  some of the lost SNR. VERNIER-16 measured how well that works with a greedy search — 70–80%
  spurious — so the recoverable fraction is small in practice, but a better search would do
  better.
- **Optimistic:** Gaussian noise averages down across pulsars; a coherent foreground does
  not, and can bias parameters rather than merely inflating their errors. VERNIER-4 hit the
  same fork (uncorrelated vs correlated confusion) and found the two branches differed by
  "N = 10⁴ pulsars" versus "never".

**A proper treatment injects the actual interloper population and searches jointly.** That is
VERNIER-16's territory and it failed its safety gate, so this instalment could not go there.

## Other things this does not settle

1. **ρ₀ = 5 per pulsar is VERNIER-6's optimistic working point.** The ceiling √(λ/(1−λ)) is
   independent of ρ₀, so the *threshold* λ ≥ 0.70/0.80 is robust; the contrast values are not.
2. **The bar is ln K with K from a 1 pc prior.** A wider prior raises the bar logarithmically
   (VERNIER-6 measured 6.1 nats per 1000×), which would push the λ threshold up slightly.
3. **PHOENIX's λ values are taken as given.** I have not re-derived them from the GLACIER
   census.
4. **Single frequency (20 nHz) for the wheel**, with the λ values applied as a uniform
   degradation rather than band-by-band.

## One-line verdict

**Turning the foreground dial to PHOENIX's measured values collapses the wheel from rank
12–13 to rank 5–7, because λ_fore imposes a brightness-independent ceiling ρ_P ≤ √(λ/(1−λ))
— so certification needs λ ≥ 0.80 and the contrast criterion needs λ ≥ 0.70, the latter
reproducing the field's stated 0.7 exactly rather than beating it; of the measured census only
the 5.41–6.27 nHz band (λ = 0.811) clears both, the median 0.58 clears neither, and the
standing rank 12–13 turns out to have been a λ = 1 result that the banked cells share by
construction.**
