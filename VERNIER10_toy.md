# VERNIER-10 — the targeted search

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~5 min.
Clean room: imports only `vernier_toy` / `vernier7` / `vernier8` / `vernier9`.

Code + data `reports/vernier10.npz`. Staged, not committed.

---

## Lead — refuted, and the reason is instructive

**D over A is 0.278 dex, not 0.49. The wheel reaches rank 13, not rank 21.**

| config | N_trials | ρ_thresh | signal | **net vs A** | **rank reached** |
|---|---|---|---|---|---|
| **A** blind sky, Earth terms only | 5.57e5 | 4.425 | 1.000 | 0.0000 | 4 |
| **B** blind sky, + pulsar terms | 7.10e13 | 5.210 | 1.584 | **0.1288** | 7 |
| **C** sky known, + pulsar terms | 1.35e3 | 4.031 | 1.584 | **0.2402** | 11 |
| **D** sky **and** frequency known | 30 | 3.691 | 1.584 | **0.2784** | **13** |

*(f_E = 20 nHz, L_max = 5 kpc, ⟨r²⟩_w = 1.508 from VERNIER-9, α = 0.01, k = 4,
N_f = 45 independent frequency bins over 2–50 nHz at T = 30 yr, N_Mc = 30.)*

Rank 21 needs **0.4272 dex**. Rank 20 needs 0.4121. **D falls 0.15 dex short. The wheel does
not sustain.**

**Why the estimate was high: the trials penalty enters logarithmically, so it was never worth
0.3 dex.** Decompose D − A:

- **signal** √(1+r²) = 1.584 → **+0.200 dex** (this is the VERNIER-8/9 ceiling, unchanged)
- **trials** 4.425 → 3.691 → **+0.079 dex**

**Going from 5.6×10⁵ trials to 30 — a factor 1.9×10⁴ — buys only 0.079 dex**, because
ρ²_thresh = k + √(2k)·Φ⁻¹(1−α/N) and the quantile grows like √(2 ln N). The whole
sky-knowledge benefit (B → C) is **0.111 dex**; making frequency free too (C → D) adds
**0.038 dex**.

**So fixing the sky removes a penalty that was only ~0.11 dex to begin with.** The dominant
term remains the capped pulsar-term power gain, exactly as VERNIER-8 and 9 found.

**Wheel confirmation** (per-pulsar distances, ρ_P = 6, null every turn): ranks 10 / 11 / 12
for A / C / D against the analytic 4 / 11 / 13. Signal w_true = 1.0000 vs null 0.0000 in
every configuration — separation is complete, and K_eff would not have shown it (VERNIER-8's
finding). The wheel's "A" row still uses pulsar terms internally, so its rank 10 is an upper
bound on a true Earth-only baseline; the analytic ladder is the cleaner statement.

## The trials factor, derived rather than assumed

Pulsar-term phase φ_a = 2π f_P τ_a with τ_a = L_a(1−cos μ_a)/c. Moving the source by δθ gives
δ(cos μ_a) = sin μ_a δθ, so

> **δφ_a = 2π f_P L_a sin μ_a δθ / c**, and the template decorrelates at
> **δθ = c / (2π f_P L_max sin μ)** — set by the *most distant* pulsar, since the product
> over pulsars decorrelates as soon as any one term does.

| f_P | L_max | δθ | **N_sky** |
|---|---|---|---|
| 10 nHz | 1 kpc | 31.9″ | 5.3e8 |
| 20 nHz | 1 kpc | 15.9″ | 2.1e9 |
| **20 nHz** | **5 kpc** | **3.19″** | **5.26e10** |
| 32 nHz | 5 kpc | 1.99″ | 1.35e11 |

**Earth-term only**: the Earth term carries no τ, so sky enters only through the antenna
pattern F₊(θ,φ), which varies on radian scales — N_sky ~ O(10²–10³) (413 for a 100 deg²
localisation).

> **N_sky(pulsar-term)/N_sky(Earth-only) = 1.3×10⁸ at 20 nHz — and my VERNIER-8 working
> figure of 1e8 sky pixels was low by ~500× in absolute terms.** Correcting it changes the
> answer by only ~0.02 dex, which is itself the point: the logarithm absorbs it.

## How well must the sky be known? Arcseconds — and that is fine

| sky prior | N_cells | ρ_thresh | net vs A |
|---|---|---|---|
| 1° | 1.27e6 | 4.816 | 0.163 |
| 10′ | 3.54e4 | 4.655 | 0.178 |
| 1′ | 354 | 4.417 | 0.201 |
| 10″ | 9.8 | 4.196 | 0.223 |
| **3.2″** | **1.0** | **4.031** | **0.240** |
| 0.1″ | 1.0 | 4.031 | 0.240 |

**The threshold stops improving at ~3″ — one decorrelation cell — so the answer to your
question is arcseconds, not arcminutes.** An arcminute prior recovers 0.201 of the 0.240 dex;
you need ~10″ for 93% of it and ~3″ to saturate.

**But that is not a practical obstacle.** An EM-identified AGN has an optical or radio
position good to ≲0.1″ — far inside the 3″ requirement. **So for an LSST-identified SMBHB
candidate the sky penalty genuinely vanishes and configuration C/D is free.** The precision
demanded is arcsecond-level; the precision available from EM identification is two orders
better.

---

## What this does not settle

1. **N_Mc = 30 and N_f = 45 are stated inputs**, not derived. Both enter logarithmically, so
   the answer is insensitive: changing either by 10× moves D−A by <0.02 dex.
2. **The chirp mass may not be free even for an EM candidate.** An orbital period fixes f, but
   Mc still has to be searched unless the EM data constrains the masses — I have kept N_Mc = 30
   in configuration D for that reason. Setting N_Mc = 1 would add only ~0.03 dex.
3. **⟨r²⟩ = 1.508 is the 20 nHz array-weighted value from VERNIER-9.** At 32 nHz it is 1.97,
   giving signal 1.723 and D−A ≈ 0.315 dex → rank 14. Still far short of rank 21, and it
   requires the near-merger sources VERNIER-9 flagged as 22× rarer.
4. **The 2N-nuisance state from VERNIER-8 is not in this ladder** — all four configurations
   assume distances are already known (k = 4). The 0.4–0.5 dex DOF effect is a *different*
   comparison (phases free vs certified) and remains the largest unsimulated term in the
   series.
5. **Still no GWB, white noise only, sources pre-resolved.** Same optimistic frame as
   VERNIER-6 through 9.

## One-line verdict

**Refuted: the targeted search buys 0.278 dex over the blind Earth-only baseline, reaching
rank 13 rather than the rank 21 that 0.49 dex would have given — because the trials factor
enters as √(2 ln N) and so a 1.9×10⁴ reduction in trials is worth only 0.079 dex, leaving the
capped 0.200 dex pulsar-term power gain as the dominant and still-binding term; the sky must
be known to ~3″ for the penalty to vanish, which any EM identification supplies with two
orders to spare, so the benefit is real and free but small.**
