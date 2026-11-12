# JOB 50 — Part A: what the machine contains

**Five-line summary.**
1. **The booked headline (36.27 at 0.9009) does not contain the pile-up. Neither does the heatmap.** Every source, target and interferer, sits at its Earth-term frequency, once. Forcing τ = 0 changes nothing in either path except through one frozen lookup table (the Earth fraction f_E), which moves the median share from 0.440 to 0.454.
2. **The pile-up is real but lives in other code:** the single-sky notebook audit (V4o), the share-versus-frequency figure and the Job 41 script. Those are where "232 against 44", "28 against 42 usable", "2.08 against 1.37 rad" and "18.7 against 1.66" come from. Conflict 1 is a mix-up between paths, not a contradiction.
3. **Every earlier "weak sources carry nothing" test used the Gaussian phase model or its σ < π cut.** None used more than 12 sources in a decode. The N-scan used at most 8, on one pulsar. Conflict 2 is untested, not settled.
4. **The template term has no derivation.** It first appears in Job 30 labelled "illustrative scaling, not calibrated". It uses the prior width where the true distance belongs: median 0.15 rad as coded, 142 rad with the true distance. Only 1.1% of decode pairs would have a template good to 1 rad.
5. Job 49's "0.303" and the 0.261 implied by 0.414 are both right; they differ because 0.414 includes the template term and 0.303 does not.

**Things in the prompt that the code shows to be wrong or unfindable**
- "D12" (per-pulsar pulsar-term SNR ≈ 0.2) is not in the repo under that name or wording. Not assessed.
- The retardation toggle "28 against 42 usable, 2.08 against 1.37 rad" is not a machine run. It is the share-versus-frequency figure (`PAPER_figures/job46_partB_figure.py:15-21`, one sky, one pulsar). I reproduce its 2.08 / 1.36 below.
- "3.4 usable per sky" was measured under the old top-40 cap and is a per-sky count of decoder combs, not per pulsar.

Evidence: `job50_partA.py`, `job50_logs/partA.log`, `job50_data/partA.json`. Notebook lines are cited as `v4 cell N:line`.

---

## A1 — Retardation, path by path

`f2` is the source's Earth-frame gravitational-wave frequency from the population draw: `SEED50[j][1]`, built at v4 cell 14:1072-1085 from `draw_population(seed=_sd)`, passed through unchanged by `_pop18` (v4 cell 18:12-17).

| Path | Target's pulsar term placed at | Each other source's pulsar term placed at | Retarded amplitude | τ from this pulsar's L, μ | Earth and pulsar terms separate | Booked numbers from this path |
|---|---|---|---|---|---|---|
| **(i-a)** notebook `sigma_phi` (v4 cell 3:262-276) | f_E (its 1/T bin) | not placed; other power is the bin total at f_E | no | no | no, one term per source | V2–V4 wheel, the V3 exhibit |
| **(i-b)** notebook V4o `sigma_audited` (v4 cell 14:1598-1617) | **retarded** `fpk = FP_ALL[i, k]` | **retarded** `w_p = np.sinc((FP_ALL[i] - fpk)*T_OBS)**2` | **no** (`P_A = h_src**2*f_src*T_OBS`, Earth power) | **yes** (`_dt_all = L_AUD[:, None]*(1 - np.cos(_mu_all))/C_PC_S`) | **yes** (`w_e` at `f_src`, `w_p` at `FP_ALL`) | V4o audit, the 232/44 table, V4o wheel verdicts; one sky |
| **(ii)** `inventory20` / `inv_rate` (v4 cell 22:24-28; job43c:498-505) | f_E (`f2[k]`) | f_E, once (`w=np.sinc((fs[lo:hi]-f2[k])*T)**2`) | no | no (no L or μ argument) | no | every SEED50 result |
| **(iii)** Job 43C headline (`collect36m`, job43c:712-756) | f_E (`fk` from `inv_rate`) | as (ii) | no | κ uses μ; the share does not | no; f_E table only (see test) | **36.27 at 0.9009; 42.35 / 40.99 / 31.57** |
| **(iv)** Job 47c heatmap (job47c:108-128) | f_E | none: σ_φ set by hand (`YS = np.linspace(0.05, 3.0, 20)`) | no | no | no | **47c / 47d grids** |
| **(v)** `job46_partB_figure.py:11-16` | **retarded** (`sON = shares(fP, hP)`) | **retarded**, pulsar terms only | **yes** (`hP = hE*(fP/fE)**(2/3)`) | yes, pulsar 41 | pulsar terms only; Earth terms left out of the denominator | pile-up figure; 28 / 42 usable; 2.08 / 1.37 rad |
| **(vi)** Job 41 (`job41_displaced_excess.py:209-214`) | n/a (population spectrum) | **retarded** (`fP=f2[None,:]*(1.0+tau/...)**(-3.0/8.0)`) | **yes** (`hP=amp[None,:]*np.maximum(fP,1e-30)**(2.0/3.0)`) | yes, 60 random pulsars | pulsar terms only, weighted by (1 − f_E) | 18.7 against 1.66 |
| **(vii)** Job 49 `p1` / `p2` (job49_common.py `contributors`) | f_E | f_E (two phasors each in `p2`, same frequency) | no | no | `p2`: yes in number, no in frequency | 36.56 at 0.9083 |

**Numerical test: force τ = 0** (set `FP_ALL` to the Earth frequencies). 10 headline skies; shares are 10th / 50th / 90th percentile.

| Path | Unit | Booked | τ = 0, nothing re-derived | τ = 0, f_E table re-derived |
|---|---|---|---|---|
| (ii) inventory share | per (sky, source), n = 106 | 0.134 / 0.405 / 0.752 | **identical** | **identical** |
| (iii) headline decode share | per (pulsar, source) pair, n = 6359 | 0.172 / 0.440 / 0.799 | **identical, bit for bit** | 0.180 / 0.454 / 0.804 (n = 6371) |
| (vii) `p2` true share | same pairs | 0.096 / 0.327 / 0.688 | **identical** | 0.103 / 0.342 / 0.689 |
| (iv) heatmap comb sets | per sky | — | **identical** | identical |
| (i-b) V4o audited share | 60 loudest sources × 3 pulsars, one sky | 0.007 / 0.039 / 0.232 | 0.008 / 0.055 / 0.619 | — |
| (v) figure, 15–40 nHz | per pulsar term, one sky, pulsar 41 | median of top-20 shares 0.188 (σ = 2.08 rad) | 0.350 (σ = 1.36 rad) | — |
| (vi) Job 41 ratio | 10 skies | 21.6 | 2.19 | — |

The one place the headline touches retardation: `depro_sigma` reads a seven-point table f_E(f), computed once from the presentation sky with retarded frequencies (v4 cell 17:96-104). With τ = 0 the table becomes 0.500 everywhere (booked: 0.495 / 0.477 / 0.455 at 3 / 10 / 25 nHz). Ungated purity moves 0.694 → 0.700 (incumbent) and 0.725 → 0.719 (`p2_mf`): no pile-up effect.

Note for (i-b): my test treats every source as undetected, because the detected-set bookkeeping needs the hasasia curve. The retardation logic is the notebook's own.

**Verdict.**
- Booked headline: **no pile-up.**
- Heatmap: **no pile-up.**
- Share-versus-frequency figure: **contains the pile-up**, pulsar terms only, retarded amplitudes, one pulsar of one sky.

## A2 — The Job 49 D4 counts

There are 120 pulsars per sky; 79 are red-noise-clean (the same 79 every sky, `rn_draw` seed 202); 69.6 per sky yield a decode record.

| Number | Unit | What it is |
|---|---|---|
| 1160.7 | sources per sky | pass ρ ≥ 3, array SNR ≥ 3, 2–30 nHz, cap 1300 |
| 16.6 | sources per sky | of those, σ² = res/Pt + 1/ρ² < π² (before any per-pulsar step) |
| 10.9 | sources per sky | the 16.6 cut to the best 12 |
| 1102.7 | sources per sky | pass the exact-statistics cut, two-term ledger |
| 9.5 | sources per (pulsar) decode record | what one pulsar actually uses |

How 16.6 and 9.5 fit: 16.6 → at most 12 per sky (mean 10.9) → each pulsar then drops sources whose template term pushes σ past π → 9.5 on average (median 12; 53% of records use all 12). All pulsars in a sky share the same list of at most 12 sources.

## A3 — The 0.303

Medians over the 132,514 decode pairs of the 200 headline skies:

| Definition | Median | s/(2 − s) |
|---|---|---|
| machine share 1/(1 + σ²), σ² = confusion + 1/ρ² + template term | 0.414 | 0.261 |
| 1/(1 + confusion + 1/ρ²), no template term | 0.465 | **0.303** |
| Job 49 "true two-term share" 1/(1 + 2·confusion + 2/ρ²) | **0.303** | — |

Job 49's 0.303 leaves the template term out of the share, because a template error is not power in the bin. The 0.414 includes it. The template term is 10% of σ² on average.

## A4 — What the earlier tests did

**Job 27 complex matched filter (+3.7%).** The data were a complex sum: `z = A*_np.exp(1j*(kap*Lt)) + (rg.normal(0,1/_np.sqrt(2),len(A)) + 1j*rg.normal(...))` with `A = 1.0/np.maximum(sgu,1e-6)` (job27_dm_oneshot_mf.py:329-330). Both decoders ran on that same datum. Depends on the Gaussian phase model: **no**. It does use A = 1/σ, which is the one-term ledger (Job 49's `p1`), and the same top-12 source list.

**Job 40 T5 N-scan, rebuilt (job40_final_batch.py:381-408).** Sources were added to the **decoder** only; the sky was fixed. One hand-picked pulsar, at most **8** sources (`NMAX=min(8,len(ks))`), 40 trials. Data: `eps=kap_all[:N]*p["Lt"]+r.normal(0,sgv)`. Weight: 1/σ² in the wrapped quadratic. Depends on the Gaussian model: **yes** (that line). The "peak at N ≈ 5" is 0.16 nats against −0.20 at N = 7, 8 (`fig40_nscan.csv`), on 40 trials of one pulsar.

**Census law, 3.4 usable per sky at ×1, ×3, ×10 (job29_joint_bin_gate.py:291-301).** Only the **timing noise** was scaled: `sig_ns = 100.0/boost`. Strains were untouched. "Usable" is the number of sources per sky returned by `inventory20`: top 40 by per-pulsar SNR, then σ < π, then the best `nmax`. Depends on the Gaussian model: **yes**, through `good = sig < np.pi` (v4 cell 22:29). The count cannot move with timing noise by construction: the top-40 cap picks the same sources, and their σ is set by res/Pt, which has no timing noise in it.

**D12.** Not found (see top).

## A5 — The template term

**Derivation:** none exists. First use: job30_middle_tier.py:293-296,
`dmu=1.0/_np.maximum(arrk,1e-6); dfq=1.0/(T*_np.maximum(arrk,1e-6))` … `sg=_np.sqrt(sg**2+(dkap*pri[i])**2)`.
The Job 30 report calls it "ILLUSTRATIVE scaling, not calibrated" (ACTIVE_WHEEL_REPORT.md:3141). Job 43C describes it as "adding template phase error δκ·L" (report :6253), but the code multiplies by `pri[i]`, the prior width.

**Does it assume a calibrated offset?** Implicitly, yes. The decoder's distance grid is an offset from the prior centre (`Lg=np.linspace(-3*pri[i],3*pri[i],601)`), and the template is κ × offset. That is only valid if the phase κ × (prior-centre distance) is known for each (source, pulsar) pair. No code estimates it.

**In the `p2` arms:** the datum is rotated by η ~ N(0, δκ·σ_d) per pair; the matched filter down-weights by e^{−(δκσ_d)²/2}; the incumbent decoder carries it inside σ.

**Sizes over the 132,514 decode pairs** (10th / 50th / 90th percentile, radians):

| Quantity | 10% | 50% | 90% |
|---|---|---|---|
| δκ·σ_d (as coded) | 0.006 | 0.147 | 0.687 |
| δκ·L (true distance, 0.5–3 kpc) | 5.3 | **142** | 868 |
| … frequency part | 1.1 | 12.0 | 57 |
| … sky-position part | 4.9 | 140 | 867 |

Fraction of pairs with δκ·L below 1 rad: **1.1%**; below π: 6.2%. The sky-position part dominates. The median source in the decode set has array SNR 37.

---

## How this changes the reading of Parts B–E

- Part B is needed: `p2_mf` has no retardation, so `p2r_mf` must be built.
- The "full" template arm in Part C will, on these numbers, scramble 99% of pairs. Q5 is close to decided before the run; the run will confirm or refute.
- The honest open question is Part C's oracle and machine-template arms with every candidate: nothing in the record tests them.

---

## Addendum A6 — the booked decode grid (found while building Part B; changes how Job 49 D5 must be read)

Every booked decoder searches distance on `Lg=np.linspace(-3*pri[i],3*pri[i],601)` (job43c:748) and counts a decode correct if `abs(Lh-Lt)<=TOOTH_W` with TOOTH_W = 0.25 pc. The grid step is 6σ_d/600. For pulsars with K ≥ 60 the median step is **1.47 pc**; for K ≥ 600 it is **7.5 pc**. A 0.25 pc window cannot be hit reliably on that grid, and the finest teeth (about 0.2–0.3 pc) are not resolved at all.

Test (100 held-out skies, same data, only the grid changes; fine grid step ≤ 0.05 pc):

| Arm | K ≥ 60, booked grid | K ≥ 60, fine grid | K ≥ 600, booked | K ≥ 600, fine |
|---|---|---|---|---|
| top 12, oracle templates, no retardation | 0.073 | 0.252 | 0.013 | 0.133 |
| every candidate, oracle templates, no retardation | 0.154 | 0.543 | 0.031 | 0.389 |
| `p2_mf` (top 12, machine template) | 0.013 | 0.017 | 0.000 | 0.000 |

**Reading:** Job 49's "K ≥ 60: 2–3%, K ≥ 600: 0 of 79" is correct for the machine as booked, but with good templates most of that failure is the grid, not physics. With the machine's template term the grid makes no difference, because the template term has already removed the information. Details in JOB50_MAIN.md.
