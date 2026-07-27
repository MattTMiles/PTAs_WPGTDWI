# ACCRE EXPORT — cronus, 2026-07-26

Pure export: no artifact was modified, regenerated, or approximated. All files in this
directory are byte-copies (`cp -p`) of the cronus originals listed below; `SHA256SUMS`
covers every copy and was verified with `sha256sum -c` after copying. Transfer path:
this repo is the established cronus→ACCRE channel (commit → push to origin → pull on
ACCRE; `hpc_harbor/` is the in-repo staging convention).

---

## Item 1 (REVISED) — SAMPLER g3 driver / configs: **DOES NOT EXIST ON CRONUS**

**Finding: no g3 driver, config, chain output, log, or checkpoint exists on this
machine.** The ACCRE working hypothesis is confirmed: the g3 driver was never written.

Evidence, exhaustively:

- `CW_transition/sampler_run.py` (working tree, clean vs git): argparse accepts
  **`mode choices=["g12"]` only** (line 265). The docstring's g3/sbc mentions are
  aspirational. The g3 *machinery* (`fringe_posterior_at_draws` line 171,
  `cert_from_qbar` line 195) exists as library functions, called only by
  `sampler_sbc.py` — there is no g3 mode or g3 entry point anywhere.
- Git history of `sampler_run.py`: **one commit** (`27dbcab`), byte-identical argparse
  (`choices=["g12"]`). No newer version exists here or in any branch.
- `SAMPLER_dev.md` §3 gate ladder records g3 as **"NOT RUN — machinery built"**; §2
  records the measured reason (1.09 s/gradient ⇒ even single-realisation gates belong
  on ACCRE).
- Searched: `CW_transition/`, `CW_node_sampling/`, `SAMPLER_results/`,
  `CW_transition/logs/`, whole project tree. Patterns: `*g3*` filenames, `g3` in all
  `sampler_*.py`, `g3` in all `CW_transition/logs/*.log`, `chain_*.npz` outputs
  (none banked). Hits were all name collisions, excluded per spec:
  - `reports/ig_*_g3_*.npz`, `reports/flat_b10_g3_A_*.npz` — IGNITE null banks
    (IGNITE's own generation tag, not SAMPLER g3);
  - CHORUS's internal gate named g3 (`reports/chorus_ch_README.md`);
  - `eig3` lines in `b1_ref_tierA*.log` (eigenvalue index).
- Per instruction, the g12-only `sampler_run.py` is **not** exported as a substitute.

## Item 1b (NEW) — IGNITE GEO union-18 membership: **FOUND** → `geo_summary.npz`

- Original: `reports/geo_summary.npz` — 145184 bytes, mtime 2026-07-12 13:55:07 −0500
- sha256 `d020911ba753f4d04c2fc256ae711f36b19a13cc7da4a61b20f6e00b7792dcf0`
- This is the file `IGNITE_onset_map.md` §2 names as the build-time verification source
  for the VLBI tier (σ_d → min(lit, 1 pc) on the GEO union-18) at the onset cells,
  including h = −13.25, T = 30 yr.
- Membership rule verified on this copy: `names[freq > 0]` yields **exactly 18
  pulsars**, matching `GEO_geometry_ensemble.md` §9.2 name-by-name and count-by-count
  (J1909-3744 38/40, J0437-4715 32/40, J1713+0747 27/40, J1744-1134 20/40,
  J0711-6830 16/40, J1045-4509 12/40, J0030+0451 8/40, J1603-7202 4/40,
  J1640+2224 4/40, J1824-2452A 4/40, J1125-6014 3/40, J2241-5236 3/40, B1937+21 2/40,
  J0613-0200 2/40, J2317+1439 2/40, J1012+5307 1/40, J1545-4550 1/40, J2222-0137 1/40).

## Item 2 — `sampler_s4b.log`: **FOUND, claim survives verbatim**

- Original: `CW_transition/logs/sampler_s4b.log` — 1008146 bytes,
  mtime 2026-07-14 10:15:51 −0500, git-untracked (exists ONLY via this export)
- sha256 `e133bdf739402afcb3a5615bf19626b821fb5e504d4f103ec70bcbfa00b0b002`
- The claim's source line is **7206**:
  `Hessian at init: 193s; eig(-H) in [-1.710e+09, 8.916e+08]; cond = 3.20e+08; neg eigs 3/8`
  — scenario C realisation [C0], seed 770001; `sampler_s4.py` initialises at truth
  (`x_tr = T.x_truth()`), so "at init" = "at truth". Companion realisations in the same
  log: C1 4/8, C2 4/8, C3 2/8, C4 2/8, B0 3/8.
- **Caveat the project lead needs:** the log records only the eigenvalue RANGE, the
  condition number, and the negative COUNT. The **full 8-eigenvalue spectrum is not
  banked anywhere** — not in this log, not in `SAMPLER_results/s4_scenB.npz` /
  `s4_scenC.npz` (which carry only `cond`). "3/8 negative" is recoverable; the
  spectrum itself is not.

## Item 3 — STORY frozen 4-seed protocol npz (D-2): **PARTIAL — standalone file no longer exists**

The adjudicated-primary Tier-C number (STORY S4.2.8, D-2: f = 0.0323 ± 0.0134) was
banked in `CW_transition/b1_referendum_tierC.npz` at the 4-seed freeze on 2026-07-09,
then **overwritten in place at 19:30:01 the same evening** when seed s4 completed and
the 5-seed aggregate was written. The 4-seed aggregate file is gone:

- On-disk `b1_referendum_tierC.npz`: `nseed = 5`, `f = 0.0431` (the 5-seed
  auto-ingest = the TRAIL value, not the D-2 primary).
- Git-committed version (`634aab8`): an **older 2-seed schema** (`N=192, seeds=2,
  f = 0.0963`) — also not the frozen 4-seed bank.

**What survives (all exported here):**

1. `b1_referendum_tierC.npz` (5-seed, sha256 `4d4c6d7d…0478`, 10109 B,
   mtime 2026-07-09 19:30:01): its `lnZbox_all[:4]` are the four frozen-protocol
   seeds. Verified reconstruction on this copy:
   mean = 405633.0345 (std 0.8586, s.e.m. 0.4293), d = lnZn_quad − mean = −3.4008 nat,
   f = e^d/(1+e^d) = **0.03227**, σ_f = f(1−f)·s.e.m. = **0.01341** — reproducing the
   adjudicated 0.0323 ± 0.0134 exactly (and matching project_progress.md 2026-07-09:
   lnZbox 405633.035, std 0.859, s.e.m. 0.429, d −3.4008).
2. The four per-seed SMC checkpoints of the frozen run, untouched since the run
   (mtimes 17:08–18:53 on 2026-07-09): `b1_ref_tC_s0_ckpt.npz` (`d5a3ed6d…9478`),
   `s1` (`07a9685e…48a3`), `s2` (`f88e2194…04e4`), `s3` (`dae5257b…63a2`). s2/s3 are
   git-untracked; s0/s1 differ from the `634aab8` blobs (those belong to the older
   2-seed run, not to any clobber — cronus mtimes predate the ACCRE D-2 keying fix
   and are consistent with the frozen run itself).
3. `b1_ref_tierC.log` (`e9226786…b18c`, 1016227 B, mtime 2026-07-09 19:30:01) — the
   frozen run's log.
4. `b1_step2_table.npz` (`af065512…a868`, 1613 B, mtime 2026-07-10 03:26:09) — the
   5-seed auto-ingest table (tier C `nseed=5, f=0.0431`), exported as the TRAIL
   record, git-untracked.

**No 4-seed aggregate was regenerated** — the reconstruction above is arithmetic on
the surviving bank, documented for verification only. STORY's own pre-publication flag
(re-run the frozen protocol at a stated seed count) stands.

## Item 4 — `data_likelihood_sandbox_MK9_executed.ipynb`: **FOUND, intact**

- Original: `CW_node_sampling/data_likelihood_sandbox_MK9_executed.ipynb` —
  1374428 bytes, mtime 2026-05-04 15:07:02 −0500, git-tracked and clean
  (working tree == committed copy).
- sha256 `770c893e35f233ded638bd283349debfc04eb002194b7f99f4358f749ff2b6b8`

---

## SHA256 (originals; copies verified identical)

```
e133bdf739402afcb3a5615bf19626b821fb5e504d4f103ec70bcbfa00b0b002  CW_transition/logs/sampler_s4b.log
4d4c6d7d7a0db55bf4517ba18adecc2069fac753429b57d2176b1efba16d0478  CW_transition/b1_referendum_tierC.npz
d5a3ed6db3d5b5f73963784d4dd9510bf9dbd648026e311fc34de7a5a0b29478  CW_transition/b1_ref_tC_s0_ckpt.npz
07a9685eef72bb401a6044d0a0bb88364addf84eed61ff121268756383b048a3  CW_transition/b1_ref_tC_s1_ckpt.npz
f88e2194afc2504a4803616181f96f488b7b3d467f4d7122b41fa35626d704e4  CW_transition/b1_ref_tC_s2_ckpt.npz
dae5257bea32b3fbb8d264ad6a4033302873d5584c0d298d6dd725b1184563a2  CW_transition/b1_ref_tC_s3_ckpt.npz
af06551256395533580c3a322e4b588ff0c4ee21be4572405ca03bd3d86aa868  CW_transition/b1_step2_table.npz
e9226786b1089be727ce22c0d37babb99ef3861189c6a4a0d9aa8942b4b7b18c  CW_transition/logs/b1_ref_tierC.log
770c893e35f233ded638bd283349debfc04eb002194b7f99f4358f749ff2b6b8  CW_node_sampling/data_likelihood_sandbox_MK9_executed.ipynb
d020911ba753f4d04c2fc256ae711f36b19a13cc7da4a61b20f6e00b7792dcf0  reports/geo_summary.npz
```
