# JOB 49 — Audit of the phase-jitter model (Parts A, B, C)

**Plain-language summary.**
1. **Verdict on H1: H1 holds as stated.** On every booked path, the machine generates each source's phase error by drawing a Gaussian with variance σ² = res/Pt (+ 1/ρ² + fuzzy) = (1 − s)/s. It never forms a complex sum. (Part A: job43c_rebooking.py:749, job47c_heatmap.py:124. Part C: matches a wrapped normal with variance (1 − β)/β to within 2 s.e. at all 7 β, and misses the Rice law by 100–180 s.e.)
2. **Each booked number reproduces exactly** with untouched code and the same seeds (36.27 @ 0.9009; 42.35 / 40.99 / 31.57; no crossing at 8.3 or 10.6; both the 47c and 47d grids). Part B reproduces both §1 tables exactly.
3. **H1 does not by itself make the machine pessimistic.** The machine's share s counts each other source **once**. Its own Earth/pulsar split (notebook cell 17, f_E ≈ 0.485) says each other source puts **two** terms of equal power in the bin. With both terms, the true confusion is 2·res, and the correct small-angle phase variance is 2·res/(2Pt) = res/Pt. That is what the machine draws. Two factors of 2 cancel. What remains is the *shape*: Gaussian versus Rice.
4. **Two procedural defects in booked code** surfaced while reproducing, unrelated to jitter:
   - `book()`'s gate smoothing zero-pads the end of its 0.5–60% search range. Train purity on the headline is still 0.919 at 60%, so the booked crossing at 58.59% is a padding artifact.
   - The real-likelihood σ_φ (stagec_p2b.py:86) is √2 too large, because it divides by an RMS rather than an amplitude.
5. **Not found:** the "exact-Bayes 0.252 / measured 0.247" numbers appear nowhere in the repo. The 0.255 toy value is the H1 model itself, so any Gaussian-draw path at σ = 1.81, N = 8, K = 11 gives about 0.25.

**Process note (deviation, stated plainly).** The brief says not to start Part D until this file is on disk. I launched the Part D computations before writing this file, once Parts A–C were complete. Every statement in Parts A–C rests on code reads and on the Part B/C runs. The ledger observation (item 3) came from the A11 factor sweep *before* any Part D number existed. The only addition after Part D is a cross-reference to the gate-artifact diagnostic in JOB49_PHYSICAL.md. Wall time: Part A ≈ 45 min reading; Part B 32 s; Part C 0.5 s.

Repo: `PTAs_WPGTDWI` (`/home/mattm/projects/HSYMT`), branch `MM_playground`. Machine files are in `wheel_hasasia_examples/`. Notebook cells are cited as `v4 cell N:line` (code cell N of `active_wheel_demo_v4.ipynb`, md5 0dca167a…, line within that cell).

---

## Part A — the code

### A-table

| Q | Answer | Evidence (file:line — quoted) | Conf. |
|---|---|---|---|
| **A1** | **(a) drawn from a Gaussian** with σ = the per-comb σ_φ, in every booked path. No path forms a complex sum of simulated contributions. | headline `job43c_rebooking.py:749` — `eps=kap*Lt+rng.normal(0,sgu)`; T3 amplitude rows `job43c:822` (identical line in `_fixed`); heatmap `job47c_heatmap.py:124` — `eps = kap*Lt + rg.normal(0, sigphi, len(kap))`; 47d `job47d_contour.py:78` (same); notebook: v4 cell 22:98 `else rng.normal(0, sgu))`, cell 27:50 `eps = kap*Ltrue + rng.normal(0, sgu)`, cell 11:20 `return kap, rho, rng.normal(0, sigphi)` | high |
| **A2** | s = Pt/(Pt + res). **Numerator:** the source's full power P0 = h²fT_OBS at its Earth frequency f2, scaled by T/T_OBS. **Denominator:** includes the source. Other sources are each counted **once**, at f2, sinc²-weighted over ±100/T (no separate Earth and pulsar terms; no retardation). The faint sea is included (scaled by T_OBS/T). **Measurement noise is not** included; it is a separate 1/ρ² term. Power is **total** (|z|², both quadratures). The one-/two-sided convention cancels because only ratios are used. | `job43c:500-505` — `c=w*Pb[lo:hi]*(T/T_OBS)` … `res=(c.sum()-P0[k]*(T/T_OBS)+sea[b]*(T_OBS/T))*rate` … `Pt=P0[k]*(T/T_OBS)` … `sg[a]=np.sqrt(max(res,0)/max(Pt,1e-300)+1/max(rp[k],1e-9)**2)`; share as used by 47b/c `job47c:83` — `s_el = 1.0/(1.0 + np.maximum(sg**2 - 1.0/np.maximum(rp, 1e-9)**2, 0.0))`; P definition v4 cell 3:233 `SRC_P = h_src**2*f_src*T_OBS`; original σ_φ v4 cell 3:275-276 | high |
| **A3** | σ² = res/Pt + 1/ρ² (inventory). It then goes through **depro**: (res/Pt)(1 − f_E E_K)(1 + K/(N−1)) + 1/ρ². Then **fuzzy**: += (δκ σ_d)². σ < π is applied three times (inventory, after depro, after fuzzy). Inventory ranks by nats = ln(π/σ) and keeps the top 12. The decoder weight is 1/σ². The gate score uses the same σ. The criterion's e^{−σ²} appears in the V18 margin. | `job43c:506` `good=sg<np.pi`; `:508` `nats=np.log(np.pi/sg); ordr=np.argsort(nats)[::-1][:nmax]`; `:730-732` `sg=np.array([depro_sigma(...)])` / `ok=sg<np.pi`; v4 cell 17:122-124 `res = sig_conf**2*(1 - fE*_EKUSE)` … `return np.sqrt(res*_LOONOISE + 1.0/max(rho, 1e-9)**2)`; `job43c:743-744` `sgu=np.sqrt(sgu**2+(dkap*pri[i])**2)` / `o2=sgu<np.pi`; v4 cell 15:10 `_PI16 = float(np.pi)  # admit combs with POSITIVE information only`; v4 cell 42:47 `marg = float(np.sum(np.exp(-sgu**2))/2) - lnKi` | high |
| **A4** | **Per source:** −wrap(κL − ε)²/(2σ²). **Stacked:** sum over sources, minus L²/(2σ_d²), argmax on a 601-point grid over ±3σ_d; correct if \|L̂ − L\| ≤ 0.25 pc. **Phase only**, weights 1/σ², no amplitude. Certification (σ_L = (Σκ²/σ²)^{−1/2} < 0.25 pc) and the gate score (alias_score24, evaluated about the true L) are **data-independent**: they depend on κ, σ, σ_d only. | `job43c:750-754` `r_=(kap[:,None]*Lg[None,:]-eps[:,None]+np.pi)%(2*np.pi)-np.pi` / `S=-(r_**2/(2*sgu[:,None]**2)).sum(0)-Lg**2/(2*pri[i]**2)` / `if 1.0/np.sqrt(max(float(np.sum(kap**2/sgu**2)),1e-30))>=TOOTH_W: continue`; v4 cell 27:67-78 `alias_score24` (`S += -rr**2/(2*sg**2)` with `rr` measured from `r["Ltrue"]`) | high |
| **A5** | ρ_k = (h/2πf)·√(n_t/2)/σ_TOA is the standard ρ² = (s\|s) for a sinusoid in white noise. In the booked paths ρ enters only as 1/ρ² in σ². The decoder −Δ²/(2σ²) ≈ −(1/σ²)(1 − cos Δ), so in the noise-only limit it is −ρ²(1 − cos): **standard**. The −(ρ²/2)(1 − cos) form lives in **v4 cell 11:27**, the V3 teaching exhibit. There ρ is a "DEMO link" (4·snr·U(0.3, 1)), not (s\|s). If ρ were (s\|s), that form would be half the standard curvature. It is on no booked path. I did not find the ρ²/2 form in any note text. | `job43c:491` `rp=(h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100.0*1e-9)`; v4 cell 11:27 `return -np.sum((rho[:n_cw, None]**2/2) *` … `(1 - np.cos(kap[:n_cw, None]*dL[None, :] - eps[:n_cw, None])), axis=0)`; cell 11 `rho = np.clip(4.0 * snr_src * rng.uniform(0.3, 1.0, K_MAX), 0, 12)  # comb SNR; 4.0 = DEMO link` | high |
| **A6** | K = 6σ_d/ΔL. There are **two conventions**. ΔL = 2π/κ_min (**coarsest** tooth) in job43B/36 records and in the 47c heatmap. ΔL = λ_min (**finest** tooth, larger K) in every notebook wheel (cells 17–42). The headline booking does not use K; only the 43B law-fit diagnostic does. | `job43c:551` `KK.append(6*pri[i]/max(2*np.pi/max(kapu.min(),1e-30),1e-9))`; `job47c:121` `dl_coarse.append(2*np.pi/kap.min())`; v4 cell 22:95 `lnKi = float(np.log(max(6*pri[i]/max(lam.min(),1e-9), 1.001)))` (`lam = 2*np.pi/kap`) | high |
| **A7** | δκ² = (2π(1−cos μ)/c)²δf² + (2πf sin μ/c)²δμ², with δμ ~ Rayleigh(mean 1/arr) and δf = 1/(T·arr). It adds **(δκ·σ_d)²** to σ², so it scales with the **prior width σ_d**, not with L. The datum receives it as part of the same single Gaussian draw. | `job43c:738-743` `md=1.0/np.maximum(arru,1e-6)` / `dmu=rng.rayleigh(md*np.sqrt(2/np.pi)) if full_dmu else md` / `dfq=1.0/(T*YR_S*np.maximum(arru,1e-6))` / `dkap=np.sqrt(...)` / `sgu=np.sqrt(sgu**2+(dkap*pri[i])**2)` | high |
| **A8** | **No timing-model projection exists** in the headline path (absent, not switched off). The only projection is the Earth-term deprojection (`depro_sigma`, LOO, K = 15). Low-frequency combs therefore carry no quadratic-spindown absorption. | searched v4 cells and job43c/40/47c for `timing model`, `tm_proj`, `spindown`, `F0F1`, `design matrix`: no hits on the booked path (the only "design matrix" is v4 cell 14:3067, an Earth-term fit) | medium (absence is proven by search, not by a positive line) |
| **A9** | The heatmap sets σ_φ directly (y-grid `YS = linspace(0.05, 3.0, 20)`). s = 1/(1 + σ_φ²) is **only an axis label**: the same σ_φ is used for every one of the top-4 combs, and data are `eps = κL + N(0, σ_φ)`. The comb frequencies come from the machine's `inv_rate` (top 4 by nats). The per-comb shares `s_el` computed at :83 are used only for the strain-invariance check, never for data. | `job47c:10` `#   s <-> sigma_phi: s = 1/(1+sigma_phi^2), applied to every comb`; `:132` `YS = np.linspace(0.05, 3.0, 20)`; `:124-126`; `:231` `sax.set_yticks(np.sqrt((1-sticks)/sticks))` | high |
| **A10** | `CW_transition/stagec_p2b.py:86-87`: `sigma_phi = sig_res / A_CW  # rad (= 1/SNR_pterm)` and `snr_pterm = A_CW / sig_res`. Here σ_res = σ_TOA·√(2/N), the per-quadrature error on a Fourier amplitude (:62-66). A_CW is the **RMS** of the pulsar-term residual, i.e. a/√2 (:57 `np.sqrt(np.mean(r_pterm ** 2))`). So SNR_pterm = (a/√2)/σ_res = ρ/√2, and σ_φ = √2/ρ: **√2 larger than the correct 1/ρ**. The ratio is an amplitude ratio against **one quadrature**, with an RMS in the numerator. Corrected array medians: σ_φ ≈ 2.13/√2 = 1.51 rad and ρ ≈ 0.66. Both are still small-angle formulas applied at ρ < 1, where the phase is far from Gaussian. | quoted lines above; `GALLERY_steve_update.md:56` | high |
| **A11** | See the sweep below. **Three factors flagged as unjustifiable**: (i) the power ratio used as a phase variance (missing ½); (ii) the one-term-per-source ledger (missing ×2); (iii) the sea's T_OBS/T scaling. (i) and (ii) cancel at small angles. **Two flagged elsewhere:** the V3 demo's ρ²/2 and the real-likelihood √2. | below | high for (i)/(ii)/real-likelihood √2; medium for (iii) |

### A1 — every code path, and which booked number used which

| Path | Generator | Booked numbers it produced |
|---|---|---|
| `job43c_rebooking.py:712-756` `collect36m` (T2) | Gaussian draw (:749) | **36.27 @ 0.9009** (step 3, marg) and steps 1–2 |
| `job43c:790-829` `_fixed` (T3) | Gaussian draw (:822) | **42.35 / 40.99 / 31.57**; no crossing at 8.3 / 10.6 |
| `job43c:618-660` `collect36f` (T1 toggle) | Gaussian (:652) | 28.49 / 34.91 toggle rows |
| `job47c_heatmap.py:108-128` `purity_cell` | Gaussian (:124) | **the 47c grid** and its caption numbers |
| `job47d_contour.py:64-82` | Gaussian (:78) | **the 47d refined grid and contour** |
| v4 cells 15, 17–25, 27, 28, 30, 32 (:47, :169), 36 | Gaussian `rng.normal(0, sg…)` | historical wheel / ladder / gate results (V5–V17) |
| v4 cell 11 (V3 exhibit) | Gaussian, ρ²/2 cosine likelihood | teaching figures only |
| **v4 cell 32:84-95** (Job 25 T2 `cmp_channels`), **cell 32:223-230**, **cell 42:52-56** (V18), `job43c:63-69` (`soft_sky`, Job 26) | **complex datum** `z = A e^{iθ} + CN(0,1)`, `A = 1/σ` | Job 25 T2 "complex MF 0.620 vs phase-only 0.598"; Job 26 soft wheel; V18 b(margin). **Not** used by the 36.27 headline, the T3 rows or the heatmap. |

The complex-datum path builds A²/P = 1/σ². It is the physical datum *at the machine's own share* (variant P1 in Part D), with the confusion modelled as a single Gaussian. It is the one place the machine has ever been "physical", and the booked headline does not go through it.

### A2/A11 — the ledger point (central finding)

The machine contradicts itself on how much confusion power a pulsar sees:

- `inventory20` / `inv_rate` sum **each other source once**, at its Earth-term frequency: `job43c:500-503`; v4 cell 22:24-26 `w = np.sinc((f2 - f2[k])*T)**2` / `c = w*Pbkg*(T/T_OBS)` / `res = c.sum() - c[k] + sea[b]*(T_OBS/T)`.
- `depro_sigma` treats that same res as **Earth + pulsar** and removes an Earth fraction f_E·E_K from it (v4 cell 17:122). f_E was measured in v4 cell 17 Task 2 from a ledger in which **every source contributes an Earth term at f_src and a pulsar term at FP_ALL, each of power SRC_P**. The saved output says so: *"f[nHz] Earth pulsar f_E … 2.50 6.93e-28 7.76e-28 0.472 … low-band (2-10 nHz) mean f_E = 0.485"*. The Earth and pulsar powers per element are about equal.
- For one pulsar, the faint sea likewise contributes both terms. This is the standard self-term: Γ(0) = 1 versus HD's ½ at zero separation.

So the per-pulsar confusion consistent with the machine's own physics is about **2·res**, and the true share is β_true = 1/(1 + 2·res/Pt) = s/(2 − s), not s. The physically correct small-angle phase variance is then 2·res·D/(2Pt) = res·D/Pt, **exactly what the machine draws**. The brief's correction (variance (1 − β)/(2β)) applied to the machine's s would halve the variance a second time. On the headline decode combs the median machine share is s = 0.414 and the median true two-term share is 0.303 (`job49_data/analysis.json` → D4).

**Consequence.** H1 is true, but it is a *shape* error, not a *scale* error, at the machine's own β. The Gaussian draw has the right small-angle variance and the wrong distribution: no Rice floor for weak sources, light tails for moderate ones. Part D therefore runs two physical ledgers. **P1** keeps one phasor per other source and isolates the generator change at the machine's own share. **P2** uses two phasors per other source, as D0 specifies, and is the physically consistent arm.

### A11 — factor-of-two (and related) sweep, source strain → decode statistic

| # | Factor | Where | Why it is there | Verdict |
|---|---|---|---|---|
| 1 | P = h²fT_OBS | v4 cell 3:233 | per-1/T-bin power in h_c² units; used only in ratios | OK (cancels) |
| 2 | ×T/T_OBS on source and discrete terms | job43c:500, :504 | coherent power ∝ T | OK |
| 3 | ×T_OBS/T on the sea | job43c:503 | In these units a dense discrete sum is T-independent (Σ sinc² ≈ n/T cancels the ×T), and a smooth background must be too. The sea should therefore carry ×1, not ×T_OBS/T. At T = 45 yr it is low by 1/1.5 relative to the discrete confusion. In sky 0 the sea is 4–30% of res per comb (`contributors()` smoke test), so res is low by ≈ 1–10%. | **FLAG: cannot justify** (not changed; it does not touch the jitter model) |
| 4 | sinc² weight | job43c:500 | power response of a length-T DFT to an off-bin sinusoid | OK |
| 5 | σ² = res/Pt | job43c:505 | power ratio used as phase variance. Correct small-angle variance is (noise power)/(2·source power): only the perpendicular quadrature rotates the phase. | **FLAG: missing ½** |
| 6 | each other source counted once | job43c:500-503; v4 cell 22:24-26 | the machine's own f_E model (cell 17 Task 2) has two equal-power terms per source; the per-pulsar self-term of a background is 2× the cross term | **FLAG: missing ×2. Cancels #5 at small angle.** |
| 7 | ρ = a√(n_t/2)/σ, a = h/(2πf) | job43c:491 | ρ² = N a²/(2σ²) is the standard (s\|s); phase variance 1/ρ² is correct | OK. (Antenna pattern omitted in collect36m; 43B's `collect_banked` had √((1/3)/R). An approximation, not a factor of 2.) |
| 8 | √(0.2·N_p) array SNR | job43c:492 | effective-fraction convention | convention-tier, not audited further |
| 9 | (1 − f_E E_K)(1 + K/(N−1)) | v4 cell 17:122-124 | Earth deprojection with LOO noise, derived in Job 17 | OK only if res is the E+P total (see #6) |
| 10 | Rayleigh scale md·√(2/π) | job43c:739 | gives mean δμ = md exactly | OK. δf = 1/(T·arr) is about 1.8× the CRB √3/(πTρ): mildly conservative. |
| 11 | κ = 2πf(1 − cos μ)/c | job43c:736 | one 2π (verified in Job 42A differential pack) | OK |
| 12 | −wrap(Δ)²/(2σ²) | job43c:751 | single-image wrapped-normal log density, matched to the generator; → −ρ²(1 − cos) | OK (single-image approximation degrades for σ ≳ 1.5) |
| 13 | −(ρ²/2)(1 − cos) | v4 cell 11:27 | V3 demo; ρ is a demo link | FLAG (exhibit only) |
| 14 | σ_φ = σ_res/A_CW, A_CW = RMS | CW_transition/stagec_p2b.py:57, :86 | RMS = a/√2 in the denominator | **FLAG: σ_φ is √2 too large** |
| 15 | K = 6σ_d/ΔL | job43c:551, cell 22:95 | ±3σ prior window | OK; two ΔL conventions (A6) |
| 16 | criterion N_eff e^{−σ²}/2 | cell 42:47 | e^{−σ²} = (Gaussian coherence)²; the /2 is absorbed by the fitted c(p) | law-tier, not derivable from code |

### Other things found while reproducing (not jitter, but they touch booked numbers)

- **G1. The booked gate is an edge artifact.** `book()` (`job43c:666-676`) smooths train purity with `np.convolve(..., np.ones(7)/7, mode="same")` over q ∈ [0.5, 60]%. That call zero-pads, so the last three smoothed points fall to 6/7, 5/7 and 4/7 of their true values. On the headline, raw train purity at q = 58 / 60% is **0.924 / 0.919**. The smoothed curve "crosses" 0.90 only because of the padding, at **58.59%**. With edge-safe smoothing and a range to 100%, the real train crossing is at **66.1%**, giving 39.69 [37.52, 41.94] @ **0.876** [0.863, 0.889] out-of-sample. So the true 0.90 train crossing does **not** generalise to 0.90, and the padding accidentally made the booked gate more conservative. The same artifact sets the 58.8% gates at A = 3.7 and 5.0 (purity 0.997 / 0.971: those rows are **not** at a 0.90 crossing). The A = 6.4 row (gate 50.3%) is a genuine crossing. (`job49_data/analysis.json` → `diag_gate`.)
- **G2. 47d reproduction depends on CSV rounding.** `job47d_contour.py:86-87` reads the grid back from `job47c_grid.csv` (4 decimals), and the seed formula `int(1e6*sigd)` then sees the rounded values. Exact XS gives different realisations: refined cells differ by rms 0.016, about binomial noise. With rounded values the 47d grid reproduces exactly.
- **G3. V4o (0e) empirical coherence check** (v4 cell 14:1575-1700, saved output): this is the one earlier time-domain test of the jitter formula. It used 200 realisations at 6 frequencies for one pulsar. The noise floor of |⟨e^{iφ}⟩| at n = 200 is about 0.06. The empirical R exceeded the formula's R at 3.0, 16.2 and 20.7 nHz (0.38 vs 0.09; 0.09 vs 0.00; 0.12 vs 0.00). Its verdict line ("Empirical R ≤ audited R ≤ old R") is supported in only 3 of 6 rows. It was too weak to detect H1.

---

## Part B — reference physics (scripts run unchanged)

`phase_check.py` (scripts and outputs in `job49_logs/partB/`): every row matches §1a to the printed precision. For example β = 0.5: rms 0.871, ⟨cos⟩ 0.710 (simulation) and 0.710 (Rice). β = 0.1: 1.481 and 0.287. Uniform rms is 1.8138.

`decode_toy.py`: every §1b row reproduces exactly (same seeds). (N, K, β) = (8, 11, 0.234): 0.687 / 0.627 / **0.255**; (4, 6, 0.2): 0.557 / 0.525 / 0.260; (4, 60, 0.2): 0.152 / 0.143 / 0.036; (4, 600, 0.2): 0.025 / 0.022 / 0.003; (4, 60, 0.5): 0.486 / 0.450 / 0.270; (4, 60, 0.8): 0.920 / 0.908 / 0.803.

The script's extra column, a Gaussian draw at σ/√2 (the correct small-angle variance at the same β), is the comparison that matters for the ledger point. At β = 0.2: 0.448 (N4, K11), 0.473 (K6), 0.101 (K60), 0.015 (K600). The physical MF gets 0.486 / 0.557 / 0.152 / 0.025. At β ≥ 0.5 the σ/√2 Gaussian slightly *beats* the physical data (e.g. K60: 0.544 vs 0.486). At equal small-angle variance, physical data win only at low β (the Rice floor) and lose a little at moderate β (heavier tails).

**Part B: PASS.**

## Part C — the machine's own generator, in isolation

`job49_partC.py` execs the **booked source lines read from the files**: share→σ is `job43c:505`, generator `job43c:749`; heatmap generator `job47c:124` with its s↔σ map. Measurement noise is off (ρ = 10¹²) and κ = 0, L = 0, so ε is the phase error. 2×10⁵ draws per β. Output in `job49_partC_table.csv`.

| β | machine rms (wrapped / raw) | machine ⟨cos⟩ | 47c path ⟨cos⟩ | Rice rms, ⟨cos⟩ | WN var (1−β)/β: rms, ⟨cos⟩ | WN var (1−β)/(2β): rms, ⟨cos⟩ |
|---|---|---|---|---|---|---|
| 0.99 | 0.100 / 0.100 | 0.995 | 0.995 | 0.071, 0.997 | 0.101, 0.995 | 0.071, 0.997 |
| 0.90 | 0.333 / 0.333 | 0.946 | 0.946 | 0.244, 0.971 | 0.333, 0.946 | 0.236, 0.973 |
| 0.70 | 0.656 / 0.656 | 0.807 | 0.807 | 0.552, 0.868 | 0.655, 0.807 | 0.463, 0.898 |
| 0.50 | 0.996 / 0.999 | 0.607 | 0.607 | 0.871, 0.710 | 0.997, 0.607 | 0.707, 0.779 |
| 0.30 | 1.430 / 1.523 | 0.313 | 0.311 | 1.166, 0.524 | 1.433, 0.311 | 1.073, 0.558 |
| 0.20 | 1.654 / 1.995 | 0.138 | 0.137 | 1.314, 0.417 | 1.658, 0.135 | 1.355, 0.368 |
| 0.10 | 1.802 / 2.998 | 0.010 | 0.014 | 1.480, 0.287 | 1.802, 0.011 | 1.694, 0.105 |

(Machine ⟨cos⟩ s.e. ≤ 0.002.) The z-scores of machine ⟨cos⟩ against each reference:

- **WN var (1−β)/β: within ±2 s.e. at every β.**
- Rice: −104 to −180 s.e.
- WN var (1−β)/(2β): −60 to −173 s.e.

Histograms at β = 0.5 and 0.2 are in `job49_partC_hist.png`. The machine's histogram lies on the red wrapped normal of variance (1−β)/β. The Rice curve is taller and narrower at β = 0.5, and at β = 0.2 it has a pronounced peak over a floor that the machine lacks.

**In situ:** in the full headline run (132,514 comb draws through `collect36m` itself), (ε − κL)/σ has mean −0.004, s.d. 1.003, excess kurtosis −0.006, and KS p = 0.31 against N(0, 1).

**Verdict line: the machine matches the wrapped normal with variance (1 − β)/β, not the Rice law and not a physical complex sum. H1 holds.**

### The 0.247 / 0.252 / 0.255 coincidence

I could not find "0.247" or "0.252" with N = 8, K = 11 or σ = 1.81 anywhere in the repo (`.md`, `.py`, `.ipynb` sources and outputs searched). The 0.255 in §1b is, by construction, the H1 generator plus its matched wrapped-normal decoder. Any machine path that draws ε ~ N(0, 1.81²) and decodes with the matching likelihood is that same calculation and must give about 0.25. That includes an "exact-Bayes" computation that uses the wrapped normal as the likelihood. The coincidence is therefore what H1 predicts. It is not independent support for the machine.

## What I could not find

- The 0.247 / 0.252 numbers (above).
- Any note text containing −(ρ_k²/2)(1 − cos(…)); the form exists only in v4 cell 11:27.
- A derivation in code of the "/2" in N_eff·e^{−σ²}/2 (it is a fitted law).
- Any timing-model projection (A8).
- Census prior widths wired into the SEED50 harness. `job45_t2_parallax_census.csv` exists, but `collect36m` draws synthetic px/pb priors.
