# THE OBSERVER'S TABLE — what to measure, in what order, and why

*Assembled by HERALD (cronus, 2026-07-24; K-vs-frequency columns backfilled from the
GENERALISE commit a8d3d97 and the sky-marginalisation slot filled from the Arm A-SKY commit
943dcf3, same day). Audience: the distance-measurement community —
PSRpi/MSPSRpi VLBI programmes, Pbdot/orbital-parallax timing groups, Verbiest-lineage
distance-catalogue readers. Every number below is an artifact readback from a banked file in
this repository — with one flagged exception, the COMPASS isotropic-control attribution in
framing item 3, whose bank is staged on ACCRE and not yet pushed here; per-column provenance
is in the footer. No number is new analysis.*

---

## The half page for an observer

If you can measure one pulsar's distance to a parsec, measure the pulsars in this table, in
this order. Here is why the order is what it is, and why it is not the order you would guess.

A continuous-wave (CW) source in the PTA band makes each pulsar's likelihood **periodic in its
distance**, with fringe spacing `dL = λ_gw/(1 − cos μ)` ≈ 0.35–0.7 pc for the best-coupled
pulsars at the top of the band. A distance prior is useful to this problem in proportion to
**K = 6σ_d/dL** — roughly, how many fringes your error bar spans. K ≤ 1 means your catalogue
value alone picks the fringe (a true "anchor"); K ≤ 3 is the census's "free anchor" bar;
K ≤ 10 means the array's data can plausibly finish the job. The fringe spacing scales as
1/f_gw, so K falls linearly as the source frequency drops — and that scaling is now banked,
exactly, per pulsar (GENERALISE Arm B1, gate g2: formula reproduces the banked census columns
to 0.0).

**The headline for an observing programme is one sentence: at the fiducial band no pulsar in
the 116-pulsar array is a free anchor — and at ~4.6 nHz exactly one (J0437−4715, K = 2.97
under its optimal prior) nearly is; and the field's loudest current candidate (the
PKS 2131-class periods) lives at that band.** That is the single best motivation for a
pc-class distance campaign, and every number in it is banked. The only other pulsar within
reach of the bar at 4.6 nHz is J2222−0137 (K_canon 13.56 → 3.39). Under canonical priors
alone, "anchors do not exist" (STORY S3.1) survives the frequency axis everywhere in
{fid, f/2, f/4}.

Three consequences an observing programme should internalise:

1. **The priority order is NOT timing precision.** Certification is statistically decoupled
   from RMS (ρ = +0.036, p = 0.7): on a real array the **loudest-coupled** pulsars, not the
   best-timed ones, carry the first certifications (STORY S3.3.9 keeper). The order below is
   the measured per-pulsar certification frequency across a 40-sky ensemble — how often that
   pulsar's fringe evidence certifies when the source sky is treated as the random variable it
   is. The top of the list (J1909−3744, 38/40 skies) is robust; the tail is a lottery. The
   sky-ensemble frequency is the right statistic here: GENERALISE Arm C has since shown that
   single-sky certification counts are not population numbers (zero of six CONFIRMED
   switch-table cells survived sky marginalisation), which the 40-sky column is immune to by
   construction.

2. **A parsec-class campaign is double-edged, and both edges are banked.** Tightening σ_d
   collapses K — the VLBI tier converts the union-18's trials mass ΣK = 88 454 → 470, which is
   detections at fixed evidence (S7.2.1), and it extends the certification onset to fainter
   sources (the re-cut floors gained the `vlbi, 2+14, T=50, h=−13.25` onset cell:
   count 0.90 → 1.17, below → ONSET; RECUT §1.1 — single-sky-set number). **But the same
   tightening lowers the trials bar that pure noise must clear, and the null floor rises to
   meet it: the banked price of the VLBI tier is +2.9 ± 1.0 nat of absolute floor at
   h = −13.25** (ANCHOR §7; STORY S7.2.8). The frequency axis runs the same race at population
   level: halving f widens every fringe and the null floor rises 3.7×; quartering it, 16× —
   **wider fringes buy anchors and pay for them in floor** (GENERALISE Arm B2; single-member
   counts, quotable as sky-median with band only). Better distances buy detections on the
   signal side and pay a floor tax on the null side. Budget both.

3. **Where the switch to eccentric-source certification sits — the sky-marginalised
   statement (Arm A-SKY, banked):** *"In a 5+11 population at h ≥ −12.75, a single e = 0.3
   member switches the certified count on in ≥ 6 of 8 independent skies (sky median 1.77,
   mean 3.18 ± 3.02). Every single-member 3+13 switch claim — including the old e = 0.5 edge
   at h = −13.0 — is a sky-lottery statement (2–4 of 8 skies), quotable only as a
   sky-ensemble median with its band, never as a per-sky guarantee."* (GENERALISE §4b.) The
   eccentricity and structure levers survive as levers; what sky marginalisation removes is
   every threshold claim that rested on a single-member 3+13 population. Loudness buys no
   floor calm: the per-sky floor swing stays ×2.5–4.3 at every h (the `floor ∝ h^1.66` race,
   confirmed along the sky axis).

   **What the band is made of — the isotropic control (COMPASS):** *"COMPASS (isotropic
   control array): the surviving cell's sky band is INTRINSIC SOURCE GEOMETRY (scatter ratio
   iso/real = 0.89, persists; strengthens to 8/8 skies, mean 4.82 under isotropy); the
   sparse-population/single-member lottery is ARRAY ANISOTROPY (ratio 0.23, collapses)."*
   Read for the table: the ±3.02 band on the surviving 5+11 / e = 0.3 cell is **not** an
   artifact of where this array's pulsars happen to sit — it survives an isotropic array
   almost undiminished, so it is a property of the source population and must be carried
   whatever the array. The 3+13 single-member lottery is the opposite: it is **this array's**
   anisotropy, and it does not survive isotropy. Neither result changes the quoting
   convention — sky-median with band, never a per-sky guarantee — but it changes what an
   observer can hope to buy: adding pulsars can retire the lottery, and cannot retire the
   band. [ACCRE: `reports/COMPASS.md`, pending push — re-verify citation at Matt's next push.]

**One flag before the table. J0437−4715 has a dual role.** It is the array's best anchor
candidate — smallest K at every frequency (11.88 → 5.94 → **2.97** under the optimal prior at
fid → f/2 → 4.6 nHz; K = 3.07 at fiducial under the Reardon+2016 0.25-pc joint-timing
composite), the only pulsar whose fringe information survives collapse of the joint sky (89%
true-fringe MAP on the failure plateau, S3.1.7), and the first name above onset (20/50 correct
certifications, 0 wrong — IGNITE §4.4, single-sky-set). **It
is also the dominant false-alarm channel:** its trials bar is the array minimum (ln K = 1.39),
so noise the model does not know about clears it first — J0437 posts the largest offender
shift in the realistic-noise null ladder (+31 under R3, ANCHOR §6) — and the VLBI floor tax
and the low-frequency floor explosion above are this same mechanism at tier and population
level. Robustness to source error and vulnerability to noise are the same property viewed from
two sides (S5.2.4). A campaign that tightens J0437 further must tighten its null calibration
with it.

---

## The table

Sorted by certification frequency across the 40-sky GEO ensemble (strict P > 0.99 frequency in
parentheses). Union-18 = every pulsar certifying in ≥ 1 of 40 skies; the ≥ 10%-of-skies set is
rows 1–10. K columns are `K = 6σ_d/dL` with per-pulsar median fiducial fringe spacing `dL`
from the frozen census bank (fiducial = D3/D4 ensemble draw, log₁₀f_gw ~ U(−8.0, −7.5), median
≈ 18 nHz; f/4 ≈ 4.6 nHz); "(opt)" uses the census's min-σ optimal prior where it differs from
canonical. σ_d(current) is the canonical catalogue entry (`best_distances.txt`); for the three
feather-kept rows the census K uses the feather prior, noted per row.

| # | pulsar | current σ_d (method) | K @ fid. f (opt) | K @ f/2 (opt) | K @ ~4.6 nHz (opt) | cert. freq (strict) | onset tier, first cert. | pc-class campaign buys |
|---|---|---|---|---|---|---|---|---|
| 1 | J1909−3744 | 25.8 pc (timing px, Verbiest+2009) | 365 (249) | 183 (125) | 91 (62) | **38/40** (27/40) | first name at census loudness (IGNITE §4.4) | K → 14.2; the workhorse gets a real prior |
| 2 | J0437−4715 | 1.3 pc (VLBI, Deller+2008); 0.97 pc timing (Reardon+2024); 0.25 pc composite (Reardon+2016) | 16.2 (11.9; **3.07** composite) | 8.1 (5.94) | 4.1 (**2.97**) | 32/40 (14/40) | first above onset; 20/50 correct, 0 wrong (IGNITE §4.4) | K → 12.3; already at the pc-class ceiling — **dual role, see flag** |
| 3 | J1713+0747 | 60.9 pc (VLBI, Chatterjee+2009) | 816 (536) | 408 (268) | 204 (134) | 27/40 (12/40) | T = 30 or h ≥ −12.75 (IGNITE §4.4) | K → 13.4 (61× collapse) |
| 4 | J1744−1134 | 17.4 pc (timing px, Verbiest+2009) | 289 (62) | 144 (31) | 72 (15.5) | 20/40 (3/40) | — (not banked) | K → 16.6 |
| 5 | J0711−6830 | 507 pc (DM YMW16; census: feather 21.2 pc — no geometric method) | 250 | 125 | 62 | 16/40 (3/40) | T = 30 or h ≥ −12.75 (IGNITE §4.4) | K → 11.8; **first geometric distance ever** for a census-triple name |
| 6 | J1045−4509 | 219 pc (timing px) | 2474 | 1237 | 618 | 12/40 (1/40) | T = 30 or h ≥ −12.75 (IGNITE §4.4) | K → 11.3 (219× collapse) |
| 7 | J0030+0451 | 5.4 pc (VLBI, Ding+2023) | 67 | 33 | 16.7 | 8/40 (1/40) | — (not banked) | K → 12.4 |
| 8 | J1603−7202 | 307 pc (timing px) | 4440 | 2220 | 1110 | 4/40 (1/40) | — (not banked) | K → 14.4 |
| 9 | J1640+2224 | 113 pc (VLBI, Ding+2023) | 1231 | 615 | 308 | 4/40 (0/40) | — (not banked) | K → 10.9 |
| 10 | J1824−2452A | — (no entry in canonical file; census: feather 1074 pc — no geometric method) | 18307 | 9153 | 4577 | 4/40 (0/40) | — (not banked) | K → 17.0; largest single collapse in the array (1074×) |
| 11 | J1125−6014 | 178 pc (Pbdot–Shklovskii; census: feather 197.6 pc) | 2396 | 1198 | 599 | 3/40 (0/40) | — (not banked) | K → 12.1 |
| 12 | J2241−5236 | 30.4 pc (timing px) | 364 | 182 | 91 | 3/40 (0/40) | — (not banked) | K → 12.0 |
| 13 | B1937+21 | 245 pc (VLBI, Ding+2023) | 3115 | 1557 | 779 | 2/40 (0/40) | — (not banked) | K → 12.7 |
| 14 | J0613−0200 | 136 pc (timing px, Hotan+2006) | 1189 (509) | 595 (255) | 297 (127) | 2/40 (0/40) | — (not banked) | K → 8.7 |
| 15 | J2317+1439 | 185 pc (timing px) | 2317 | 1158 | 579 | 2/40 (1/40) | — (not banked) | K → 12.5 |
| 16 | J1012+5307 | 30.8 pc (VLBI, Ding+2023) | 119 | 60 | 30 | 1/40 (0/40) | — (not banked) | K → 3.9 (widest fringe in the union, dL = 1.55 pc) |
| 17 | J1545−4550 | 432 pc (timing px) | 6271 | 3135 | 1568 | 1/40 (0/40) | T = 30 or h ≥ −12.75 (IGNITE §4.4) | K → 14.5 |
| 18 | J2222−0137 | 0.97 pc (VLBI, Ding+2024) | **13.6** | 6.8 | **3.39** | 1/40 (0/40) | — (not banked) | already pc-class: K → 13.9 at fid. f — its payoff is at low f, where it is the only other pulsar within reach of the K ≤ 3 bar |

**Reading notes.**
- *K @ f/2 and K @ ~4.6 nHz* are exact banked re-cuts, not simulations: at fixed geometry
  `dL ∝ 1/f_gw`, so `K(f) = K_fid·(f/f_fid)` pulsar-by-pulsar (GENERALISE Arm B1; gate g2
  proved the formula reproduces the banked fiducial table to 0.0). Population summary at the
  census bars — K_canon ≤3/≤10/≤30: fiducial 0/0/2, f/2 0/2/2, f/4 0/2/4 (K_opt at f/4:
  **1**/2/5, the one being J0437).
- *"K → x" in the last column* is `6 × (1 pc)/dL_med` — the K this pulsar would have at the
  fiducial frequency with a 1-pc distance, computed from the banked per-pulsar `dL_med` with
  the census's own banked formula. It is a derived target, not a bank readback, and it makes
  the structural point: **at the fiducial frequency, pc-class distances land at K ≈ 9–17 —
  detection-grade trials reduction, not anchors.** The anchor regime opens at low frequency:
  a 1-pc pulsar at 4.6 nHz sits near K ≈ 3 (J0437's banked 2.97 is exactly this point reached
  with today's best prior).
- *The frequency race has a null side*, and it is banked: at f/2 the certification floor rises
  3.7× (count 1.57, CONFIRMED); at f/4 it rises 16× and the count dies (0.27, below) —
  GENERALISE Arm B2. These are single-member (m1e05, 3+13) counts: per the A-SKY convention,
  sky-median with band, not a per-sky guarantee. Wider fringes buy anchors and pay for them
  in floor. A low-frequency anchor campaign is a campaign *against this floor*, which is
  precisely what a K ≈ 3 prior is for: it removes the trials mass the floor feeds on.
- *Onset-tier caveat:* the per-pulsar arrival ordering is IGNITE §4.4 prose (banked); the
  onset-cell **verdicts** were subsequently re-cut (SURFACE/RECUT) — the confirmed onset cells
  under corrected floors are (−12.50, 30, lit) and the RECUT frontier set (7/18 columns). All
  per-cell onset counts quoted anywhere in this document are single-sky(-set) numbers, and the
  A-SKY precedent (boundary-cell counts swing across the >1 bar with sky) applies: quote them
  as sky-median-with-band-class statements, not per-sky guarantees. Per-pulsar onset tiers for
  the unnamed rows were never banked (the GENERALISE bank carries K columns, not per-pulsar
  onset tiers).

---

## LaTeX block (AASTeX, makecell convention, no p{} columns)

```latex
% Requires \usepackage{makecell} — AASTeX deluxetable*, all-l columns per repo convention.
\begin{deluxetable*}{llccccccc}
\tablecaption{The observer's table: union-18 CW-certifying pulsars, sorted by certification
frequency across the 40-sky geometry ensemble. $K \equiv 6\sigma_d/d_L$ counts prior-spanned
fringes at the per-pulsar median fiducial fringe spacing $d_L$ (fiducial
$\log_{10}f_{\rm gw}\sim U(-8.0,-7.5)$, median $\approx18$\,nHz; $f/4\approx4.6$\,nHz);
$K\leq3$ is the census free-anchor bar, $K\leq10$ the data-completable regime. Parenthesised
values use the census min-$\sigma$ optimal prior where it differs from canonical.
\label{tab:observers}}
\tablehead{
  \colhead{\#} &
  \colhead{Pulsar} &
  \colhead{\makecell{$\sigma_d$ (method)}} &
  \colhead{\makecell{$K$ \\ fid.~$f$}} &
  \colhead{\makecell{$K$ \\ $f/2$}} &
  \colhead{\makecell{$K$ \\ $4.6$\,nHz}} &
  \colhead{\makecell{cert.\ freq.\ \\ (strict)}} &
  \colhead{\makecell{onset tier, \\ first cert.}} &
  \colhead{\makecell{1-pc campaign: \\ $K\rightarrow$}}
}
\startdata
 1 & J1909$-$3744 & \makecell[l]{25.8\,pc (timing px)} & 365 (249) & 183 (125) & 91 (62) & $38/40$ ($27/40$) & \makecell[l]{first, census loudness} & 14.2 \\
 2 & J0437$-$4715\tablenotemark{a} & \makecell[l]{1.3\,pc (VLBI); \\ 0.25\,pc (composite)} & 16.2 (11.9; 3.07) & 8.1 (5.94) & 4.1 (\textbf{2.97}) & $32/40$ ($14/40$) & \makecell[l]{first above onset; \\ 20/50 correct, 0 wrong\tablenotemark{b}} & 12.3 \\
 3 & J1713$+$0747 & \makecell[l]{60.9\,pc (VLBI)} & 816 (536) & 408 (268) & 204 (134) & $27/40$ ($12/40$) & \makecell[l]{$T=30$ or $h\geq-12.75$} & 13.4 \\
 4 & J1744$-$1134 & \makecell[l]{17.4\,pc (timing px)} & 289 (62) & 144 (31) & 72 (15.5) & $20/40$ ($3/40$) & \nodata & 16.6 \\
 5 & J0711$-$6830 & \makecell[l]{507\,pc (DM); \\ feather 21.2\,pc} & 250 & 125 & 62 & $16/40$ ($3/40$) & \makecell[l]{$T=30$ or $h\geq-12.75$} & 11.8 \\
 6 & J1045$-$4509 & \makecell[l]{219\,pc (timing px)} & 2474 & 1237 & 618 & $12/40$ ($1/40$) & \makecell[l]{$T=30$ or $h\geq-12.75$} & 11.3 \\
 7 & J0030$+$0451 & \makecell[l]{5.4\,pc (VLBI)} & 67 & 33 & 16.7 & $8/40$ ($1/40$) & \nodata & 12.4 \\
 8 & J1603$-$7202 & \makecell[l]{307\,pc (timing px)} & 4440 & 2220 & 1110 & $4/40$ ($1/40$) & \nodata & 14.4 \\
 9 & J1640$+$2224 & \makecell[l]{113\,pc (VLBI)} & 1231 & 615 & 308 & $4/40$ ($0/40$) & \nodata & 10.9 \\
10 & J1824$-$2452A & \makecell[l]{--- (feather 1074\,pc)} & 18307 & 9153 & 4577 & $4/40$ ($0/40$) & \nodata & 17.0 \\
11 & J1125$-$6014 & \makecell[l]{178\,pc ($\dot{P}_b$--Shk); \\ feather 197.6\,pc} & 2396 & 1198 & 599 & $3/40$ ($0/40$) & \nodata & 12.1 \\
12 & J2241$-$5236 & \makecell[l]{30.4\,pc (timing px)} & 364 & 182 & 91 & $3/40$ ($0/40$) & \nodata & 12.0 \\
13 & B1937$+$21 & \makecell[l]{245\,pc (VLBI)} & 3115 & 1557 & 779 & $2/40$ ($0/40$) & \nodata & 12.7 \\
14 & J0613$-$0200 & \makecell[l]{136\,pc (timing px)} & 1189 (509) & 595 (255) & 297 (127) & $2/40$ ($0/40$) & \nodata & 8.7 \\
15 & J2317$+$1439 & \makecell[l]{185\,pc (timing px)} & 2317 & 1158 & 579 & $2/40$ ($1/40$) & \nodata & 12.5 \\
16 & J1012$+$5307 & \makecell[l]{30.8\,pc (VLBI)} & 119 & 60 & 30 & $1/40$ ($0/40$) & \nodata & 3.9 \\
17 & J1545$-$4550 & \makecell[l]{432\,pc (timing px)} & 6271 & 3135 & 1568 & $1/40$ ($0/40$) & \makecell[l]{$T=30$ or $h\geq-12.75$} & 14.5 \\
18 & J2222$-$0137 & \makecell[l]{0.97\,pc (VLBI)} & 13.6 & 6.8 & \textbf{3.39} & $1/40$ ($0/40$) & \nodata & 13.9 \\
\enddata
\tablenotetext{a}{Dual role: best anchor candidate at every frequency (array-minimum $K$;
$K_{\rm opt}=2.97$ at 4.6\,nHz — the array's only free-anchor crossing, and the band where
the field's loudest current candidate lives) and dominant false-alarm channel (array-minimum
trials bar $\ln K = 1.39$; largest offender shift, $+31$, under realistic noise; the VLBI
tier's $+2.9\pm1.0$\,nat floor tax is this mechanism at tier level).}
\tablenotetext{b}{Per-cell onset counts are single-sky(-set) statements, quotable as
sky-ensemble medians with their band, never as per-sky guarantees (GENERALISE Arm~A-SKY
convention). The 40-sky certification-frequency column is unaffected (it is the
sky-marginalised statistic).}
\tablecomments{Sorted by certification frequency (40-sky ensemble, Bayesian gate; strict
$P>0.99$ in parentheses). $K(f) = K_{\rm fid}\,(f/f_{\rm fid})$ exactly at fixed geometry
(banked re-cut, gate-verified). ``$K\rightarrow$'' is the fiducial-$f$ $K$ under a 1\,pc
distance: detection-grade trials reduction, not anchors --- the $K\leq3$ anchor regime opens
only at low frequency. A pc-class campaign is double-edged: it extends onset (the re-cut
floors gained the vlbi $2{+}14$, $T{=}50$, $h{=}{-}13.25$ cell) and pays $+2.9\pm1.0$\,nat of
null floor at $h=-13.25$; on the frequency axis the null floor rises $3.7\times$ at $f/2$ and
$16\times$ at $f/4$ --- wider fringes buy anchors and pay for them in floor.}
\end{deluxetable*}
```

---

## Per-column provenance (one line each)

- **pulsar / row set** — union-18 (certifies in ≥1 of 40 skies) read from `reports/geo_summary.npz::names,freq`; the ≥10%-of-skies set is rows 1–10 of the same array. B1937+21 = J1939+2134 in the K-table banks.
- **current σ_d + method** — canonical catalogue `CW_transition/best_distances.txt` (frozen, 116 rows = 69 parallax-class + 47 feather-kept; method strings verbatim); feather values for the 3 no-geometric-method rows from `CW_transition/anchor_a1_Ktable.npz::sig_canon_pc`.
- **K @ fiducial f (opt)** — `reports/gen_b1_kcensus.npz::K_canon_f1 (K_opt_f1)` ≡ `anchor_a1_Ktable.npz::K_canon (K_opt)` (gate g2: max|diff| = 0.0), at per-pulsar `dL_med` (median over 10 D3 geometry draws; fiducial log₁₀f_gw ~ U(−8.0, −7.5), median ≈ 1.8e−8 Hz); J0437 composite K = 3.07 from STORY S3.1.3 / `stagec_anchor_a2.py::LIT_INJECT`.
- **K @ f/2, K @ ~4.6 nHz (opt)** — `reports/gen_b1_kcensus.npz::K_canon_f2/K_opt_f2, K_canon_f4/K_opt_f4` (GENERALISE Arm B1, exact banked re-cut, `K(f) = K_fid/div`, fdivs = [1,2,4]; f/4 ≈ 4.6 nHz); population bars and the J0437/J2222 crossings from `reports/GENERALISE.md` §1.
- **certification frequency (strict)** — `reports/geo_summary.npz::freq (freq_strict)`, 40 isotropic sky redraws, zero-noise Asimov, literature priors (GEO §3 / STORY S3.3.2). **This is the sky-marginalised statistic and is unaffected by GENERALISE Arm C's finding that single-sky certification counts are not population numbers — the 40-sky ensemble is exactly the marginalisation Arm C demands.**
- **onset tier, first cert. (named rows)** — IGNITE §4.4 per-pulsar arrival prose (`reports/IGNITE_onset_map.md`); onset-cell verdicts as re-cut by `reports/RECUT_floors.md` §1.1 and SURFACE (confirmed cell (−12.50, 30, lit)). All per-cell counts are single-sky(-set), sky-median-with-band class per the A-SKY convention. Per-pulsar onset tiers for unnamed rows: never banked.
- **pc-class campaign column** — derived: 6 × (1 pc)/`dL_med` (banked formula K = 6σ_d/dL, STORY S3.1.1; banked dL); double-edge citations: onset extension RECUT §1.1 (vlbi 2+14 T=50: 0.90 → 1.17, below → ONSET; single-sky-set), floor tax +2.9 ± 1.0 nat at h = −13.25 from ANCHOR §7 (`reports/ANCHOR_realdata_null.md`) / STORY S7.2.8; frequency-axis floor race (×3.7 at f/2, ×16 at f/4; single-member, sky-median-with-band class) from `reports/GENERALISE.md` §3 (Arm B2).
- **J0437 dual-role flag** — anchor side: `gen_b1_kcensus.npz` (K_opt 11.88 → 5.94 → 2.97), STORY S3.1.3 (K_lit = 3.07), S3.1.7 (89% true-fringe MAP on the R plateau), IGNITE §4.4 (20/50 correct, 0 wrong; single-sky-set); false-alarm side: ANCHOR §6 (+31 offenders under R3, table row J0437), ln K = 1.39 array minimum (quoted in ANCHOR), convention STORY S5.2.4.
- **sky-marginalised switch statement (framing item 3)** — GENERALISE §4b (Arm A-SKY, commit 943dcf3; bank `reports/gen_armAS_sky.npz`, figure `GENERALISE_AS_sky.png`): sole surviving boundary cell (5+11, e = 0.3, h = −12.75) at 6/8 skies, sky median 1.77, mean 3.18 ± 3.02; all single-member 3+13 threshold claims sky-lottery (2–4/8); floor sky-swing ×2.5–4.3 at every h. Single-member counts throughout this document carry the A-SKY convention: sky-median with band, never per-sky guarantee.
- **band-attribution sentence (framing item 3, isotropic control)** — COMPASS, isotropic control array: surviving-cell sky band = intrinsic source geometry (scatter ratio iso/real = 0.89, persists; strengthens to 8/8 skies, mean 4.82 under isotropy); sparse-population/single-member lottery = array anisotropy (ratio 0.23, collapses). **[ACCRE: `reports/COMPASS.md`, pending push — re-verify citation at Matt's next push.]** This is the one entry in this footer not yet readable from a bank in this repository; every other number above is a local artifact readback.
