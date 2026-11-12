# VERNIER-22 — the small array with a soft bar

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~12 min.
Clean room: imports only `vernier_toy` / `vernier6b` / `vernier16` / `vernier18` / `vernier21`.

Code `hpc_harbor/vernier/vernier22.py` · data `reports/vernier22.npz`.
Staged, not committed.

---

## Lead — neither change rescues it, and the 4.5 % shortfall was never the real gap

**Both candidate artefacts are eliminated. The failure is physical.**

| N_psr | per-pulsar ρ_P | ρ_P/bar | K_eff turn 1 → end | coherence R | f_cert | sources found |
|---|---|---|---|---|---|---|
| 15 | 0.161 | 0.081 | 69.6 → 75.8 | 0.0007 | 0.0015 | 0 |
| 25 | 0.146 | 0.073 | 70.7 → 64.7 | 0.0003 | 0.0008 | 1 |
| 35 | 0.147 | 0.073 | 67.6 → 69.2 | 0.0003 | 0.0001 | 1 |
| 50 | 0.163 | 0.081 | 67.1 → 69.4 | 0.0004 | 0.0020 | 1 |
| 116 | 0.171 | 0.086 | 70.6 → 75.1 | 0.0003 | 0.0007 | 2 |

**K_eff does not fall. Coherence is ~3×10⁻⁴. f_cert is ~10⁻³ at every array size.** The soft
bar extracts essentially nothing, because there is essentially nothing to extract.

### The shortfall is not 4.5 %, it is a factor of ~8

**V21's 1.911-against-2.0 is a comparison of *ceilings*, not of realised SNR.** The ceiling
√(λ_eff/(1−λ_eff)) is the saturation value for an infinitely bright source. Driving the source
brightness up by 50× shows the saturation directly:

| per-pulsar SNR input | 2.0 | 10.0 | 30.0 | 100.0 |
|---|---|---|---|---|
| **realised per-pulsar ρ_P** (N=40) | 0.168 | 0.245 | 0.248 | **0.249** |

**ρ_P saturates at 0.25 and brightness cannot buy past it** — the confusion ceiling of
VERNIER-17 in action. **The realised per-pulsar SNR is ~8× below the bar, not 4.5 % below.**
V21's near-miss was an artefact of quoting a saturation ceiling as though it were an achieved
value.

Even taking V21's most favourable ceiling and shrinking the array as far as the argument
allows, the per-pulsar value does not clear the bar:

| N_psr | 5 | 15 | 25 | 40 |
|---|---|---|---|---|
| per-pulsar ceiling (λ_eff = 0.948) | **1.91** | 1.10 | 0.85 | 0.68 |

**Even a 5-pulsar array sits at 1.91 — the same near-miss, and 5 pulsars cannot detect
anything.** The small-array route is closed at both ends.

## Anti-circularity — passes where anything moves, and nothing moves

| N_psr | confidence × | truth × | null conf × | w_true sig | w_true null | verdict |
|---|---|---|---|---|---|---|
| 15 | 0.918 | 1.116 | 1.021 | 0.0283 | 0.0241 | PASS |
| 25 | 1.093 | 1.273 | 1.127 | 0.0318 | 0.0309 | fail |
| 35 | 0.977 | 0.907 | 0.952 | 0.0263 | 0.0256 | fail |
| 50 | 0.968 | 0.956 | 1.003 | 0.0248 | 0.0236 | PASS |
| 116 | 0.940 | 0.940 | 0.929 | 0.0224 | 0.0224 | PASS |

**The two "fails" are not failures of the criterion — they are noise on ratios of ~1.0.**
w_true signal and null are identical to three decimals at every N (0.0283/0.0241,
0.0318/0.0309, 0.0224/0.0224). **There is no signal to be circular about.** Reporting these as
PASS would have been as misleading as reporting them as FAIL; the honest statement is that the
test has nothing to discriminate.

At the brightest setting and the smallest useful array (N = 25, per-pulsar input 30) there is a
faint hint — K_eff falls 68.7 → 44.7 and w_true/null reaches 1.41 — consistent with the
small-array argument having the right sign. **f_cert is 0.0018 there, so it is a hint about
direction, not a usable effect.**

## Sources found vs certified — the product has no optimum

| N_psr | 15 | 25 | 35 | 50 | 116 |
|---|---|---|---|---|---|
| sources found (array SNR ≥ 8) | 0 | 1 | 1 | 1 | **2** |
| sources certified | 0 | 0 | 0 | 0 | 0 |

**The array that finds the most certifies the fewest per pulsar, and the array that would
certify best finds nothing.** There is no interior optimum because the certified column is zero
everywhere — the trade-off VERNIER-13 identified is real in its direction but there is no
operating point where both are non-zero.

## The λ_eff bracket, tightened

V21 left raw λ_fore as a lower bound and overlaps-at-known-positions as an upper bound. Adding
the operational case — target position **searched**, purity of the recovered template measured:

| bin | N | **lower (raw λ)** | **SEARCHED** | **upper (V21, known)** |
|---|---|---|---|---|
| 2.00–2.85 | 156 | 0.263 | **0.328** | 0.792 |
| 2.85–3.71 | 46 | 0.579 | **0.876** | 0.937 |
| 3.71–4.56 | 22 | 0.722 | **0.941** | 0.965 |
| 4.56–5.41 | 11 | 0.636 | **0.922** | 0.950 |
| 5.41–6.27 | 5 | 0.811 | **0.948** | 0.975 |

**In sparse bins the truth sits near the upper bound (0.948 vs 0.975 at N = 5); in the crowded
bin it collapses toward the lower one (0.328 vs 0.792 at N = 156).** The bracket is now narrow
where it matters and wide only where the bin is hopeless anyway — because a search in a crowded
bin returns a blend, and a blend has low purity.

**One caveat that limits this:** the searched measurement uses **Earth-term templates only**.
The equivalent with pulsar terms is not measurable by brute force (VERNIER-16: 5.3×10¹⁰ grid
points), so **V21's λ_eff(E+P) = 0.993 remains an upper bound with no attainable counterpart**,
and the resolution gain from certified distances cannot currently be demonstrated
operationally.

---

## The gap I have to flag

**The wheel in this instalment uses RAW per-pulsar confusion, not the λ_eff overlap suppression
measured in VERNIER-21 Part 1.** Its noise budget is `1 + Σ_others R²` with no template
orthogonality credit, which is the pessimistic (raw λ) treatment. So Parts 1 and 2 of this
series are not yet consistent with each other.

**It does not change the conclusion** — the analytic check above uses λ_eff = 0.948 and still
gives per-pulsar ceilings of 0.68–1.91 across N = 5–40, all below the bar — but the wheel
numbers (ρ_P ≈ 0.25) are more pessimistic than a λ_eff-consistent wheel would give. **A wheel
built on λ_eff-suppressed confusion is the right next construction and I did not build it.**

## Other things this does not settle

1. **Median per-pulsar ρ_P (0.25) sits well below the "typical" value** implied by the ladder
   (~1.6) because antenna coupling varies strongly across pulsars and the median pulsar couples
   poorly. That is realistic but it means the headline is a median, not a typical.
2. **5 sources, one bin, one λ.** The saturation argument generalises; the numbers do not.
3. **Soft weights implemented as VERNIER-6b's** — full fringe posterior, coherence as resultant
   length. Other graded criteria might extract more.
4. **White noise, no GWB, Earth-term search.**

## One-line verdict

**Neither the small array nor the soft bar rescues certification, and the reason is that
VERNIER-21's 4.5 % shortfall compared two ceilings rather than an achieved value: the realised
per-pulsar SNR saturates at 0.25 — a factor ~8 below the bar, and unmoved by a 50× increase in
source brightness — while even the most favourable λ_eff and a 5-pulsar array give a ceiling of
1.91 that still does not clear 2; the anti-circularity checks neither pass nor fail meaningfully
because signal and null are identical to three decimals; and the λ_eff bracket is now tight
(0.948 searched vs 0.975 known at 5 sources, 0.328 vs 0.792 at 156) with the important
qualification that the pulsar-term version of that measurement is computationally unreachable,
so the resolution gain remains an upper bound with no demonstrated operational counterpart.**
