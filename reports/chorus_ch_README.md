# CHORUS launch sequence (agent notes, untracked)

Tag `criterion-v2.1` @ 6bec3d6. All jobs via `ch_job.sbatch` / `ch_gate.sbatch`,
lane `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1`.

Order (STOP at any gate failure):
1. `sbatch ch_gate.sbatch`                      # g1, g2, g3 (g3a CPU part already PASS on login)
2. `sbatch --job-name=ch_warm --time=05:00:00 ch_job.sbatch warm 1`   # 4 graph shapes, sizes arrays
3. RESUME DRILL: submit nulls array, `scancel` one task mid-run after >=1 npz banked,
   resubmit same command, verify "already banked: N; to run: M" line. Production license.
4. `sbatch --job-name=ch_nulls --array=0-59%12 --time=<sized> ch_job.sbatch nulls 60`
5. `harbor_py chorus.py floors` (login node, CPU) -> ch_floors.npz (g3 convention: N + fit
   error; m0e00_vlbi reuses the IGNITE-2 banked 150/270-null cell, provenance recorded)
6. `sbatch --job-name=ch_stage1 --array=0-31%12 --time=<sized> ch_job.sbatch stage1 32`
7. `sbatch --job-name=ch_pairs --time=<sized> ch_job.sbatch pairs 1`
8. `sbatch --job-name=ch_sloop --array=0-3 --time=<sized> ch_job.sbatch sloop 4`
9. `harbor_py chorus_analysis.py` (login, CPU) -> ch_analysis.npz + printed tables
10. report CHORUS_mixed_pop.md

Budget ledger (fill from warm):
- per-realisation wall by shape (16/47/78/109 slots): __ / __ / __ / __ s
- nulls: 25 cells x 160 = 4000; stage1: 26 x 30 = 780; pairs: 40; sloops: 30
- disk: lean npz ~15-40 KB x ~4900 ~ 200 MB (fine; group headroom 3.1 TiB)

e-assignments (C8, pre-registered, drawn once):
- eU n_ecc=1: [0.5801]
- eU n_ecc=2: [0.5631, 0.6753]
- eU n_ecc=3: [0.3522, 0.2957, 0.6792]

Active-harmonic counts at census parameters (mc=9, fgw=-7.9): e=0.3 -> 8, e=0.5 -> 17,
e=0.7 -> 32, e=0.8 -> 32 (magnitude-truncated per ATLAS validity note). n_clip = 0.

## REPRODUCIBILITY HAZARD FOUND BY THE g1 GATE (report to cronus)
Banked noisy realisations reproduce bit-identically ONLY at --cpus-per-task=8.
NoiseDrawer's GWB square root (np.linalg.eigh of the 5336x5336 near-degenerate Phi_gwb, CPU
LAPACK) rotates its eigenvector basis with the BLAS thread count; L_gwb @ z at the same seed
is then a DIFFERENT equal-distribution GWB realisation (dlnL shifts O(0.1-1) nat everywhere,
certified/detected sets unchanged in the case measured). Verified: 8 cpus -> bit-identical to
banked (jobs 12507594/12507755, 3 cards, 2 nodes); 16 cpus -> reproducibly shifted (jobs
12507521/12507900, full anatomy in their logs). ALL chorus jobs pinned to 8 cpus.
Durable fix (cronus's call): deterministic single-threaded factorisation in NoiseDrawer, or
bank L_gwb.
