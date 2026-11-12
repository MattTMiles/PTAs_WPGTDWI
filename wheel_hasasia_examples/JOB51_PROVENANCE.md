# JOB 51 Part 0 — Provenance

**Verdict: this directory cannot be shown to be the latest version by version control, but nothing newer exists anywhere I can read.**
`wheel_hasasia_examples/` has never been committed on any branch, so git cannot vouch for it. Against that: a search of every readable disk location found exactly one copy of the machine notebook and of every key script, all in this repository, and the three saved notebook values I checked reproduce today.

## Five-line summary

1. **No git history for the machine.** All 383 tracked files in `wheel_hasasia_examples/` are staged-but-never-committed. HEAD is dated 2026-07-26; every job since (26 to 50) lives only in the working tree and index of one machine.
2. **No other copy on this machine.** One `active_wheel_demo_v4.ipynb`, one of each key script. The only duplicate is an older PDF export in `~/soft/`. No v5, "final" or "paper" notebook exists.
3. **No branch or stash holds a different version.** No branch has any commit touching `wheel_hasasia_examples/`. The remote had 4 commits not yet fetched (July 29–30); none touch the three directories.
4. **Notebook outputs agree with the code as it runs now** (4 of 4 values checked, and three notebook figures regenerate byte-identical).
5. **Gaps Matt should know about:** `fig_phase_jitter_from_s.png` does not exist on disk; no trace of Job 44 anywhere; Job 45 has scripts and tables but no report; two figures depend on data that is gone or older than its script; other machines (ACCRE) and other users' home directories could not be checked.

## Things in the prompt that the disk contradicts

- **`inventory20` and `SEED50` are not files.** `inventory20` is a function in notebook cell 22; `SEED50` is a variable in notebook cell 14. No file or directory with either name exists (the only `seed50` hits are unrelated SCULPT outputs under `~/soft/SCULPT/`).
- **`job43c_rebooking.py` is not a standalone script.** Its first line is `# ================= JOB 26 -- THE SOFT WHEEL`; it is a 45 kB accumulation of Jobs 26–43C meant to run inside the notebook's namespace. Jobs 47c/47d/49/50 execute selected line spans of it (see 0.5).
- **`fig_phase_jitter_from_s.png`: not found** anywhere on this machine, and no script, notebook or report mentions that filename.
- **Corner versus amplitude (Job 43C T3): no plot exists.** `job43c_rebooking.py` saves only `fig1_comb.png` and `fig5_gate.png`.

Everything below was produced by read-only shell commands on 2026-11-11 (`git`, `find`, `sha256sum`, `stat`), plus one `git fetch origin` (updates remote-tracking refs only; working tree untouched).

---

## 0.1 The repository

| Item | Value |
|---|---|
| Top level | `/home/mattm/projects/HSYMT` |
| Remote | `origin` = github.com/MattTMiles/PTAs_WPGTDWI |
| Branch | `MM_playground` |
| HEAD | `5461142`, 2026-07-26 23:37 −0500, "EXPORT: ACCRE artifact package 2026-07-26 (hpc_harbor/export_20260726)" |
| Versus remote | local is 1 commit ahead and 4 behind `origin/MM_playground` (`86b3302`, 2026-07-30) |
| Changed tracked/staged paths | 559 (546 of them staged) |
| Untracked paths | 104 |
| Stashes | 2, both "WIP on master: dfbd6ec", dated 2026-04-14 and 2026-04-17; each touches only `data_likelihood_sandbox_sampling.ipynb` |
| Worktrees | 1 (the repository itself) |

The 4 remote-only commits (ed78a1c, c234d62, 657db13, 86b3302; 2026-07-26 to 07-30; "GLACIER…", "update", "update", "toy problem updates") touch none of the three directories (`git log MM_playground..origin/MM_playground -- <dirs>` is empty; `git diff --stat` between the two over the three directories is empty). They were pushed from somewhere other than this working tree.

## 0.2 History of the three directories

| Directory | Commits touching it (HEAD branch) | Last commit | Files on disk / in index | Most recent uncommitted change |
|---|---|---|---|---|
| `wheel_hasasia_examples/` | **0** | none, ever | 858 / 383 (all 383 are new-file staged) | 2026-11-11 (Job 50 files); last pre-Job-49 change 2026-09-25 (`one_pulsar_six_worlds.ipynb`) |
| `PAPER_figures/` | 1 | `863d5cd` 2026-07-23 "pub plots potentially" | 88 / 59 (11 staged new, 13 untracked) | 2026-11-11 16:48 (`job50_share_vs_frequency_two_term.png`); before that 2026-09-25 12:35 (Job 48b) |
| `CW_transition/` | many (15 listed below) | `a8d3d97` 2026-07-24 "GENERALISE check" | 291 / 217 (13 modified, 7 staged new, 48 untracked) | 2026-09-25 12:10 (`fisher_nested_fast.py`, **untracked**) |

Last 15 commits touching `CW_transition/`:

```
a8d3d97 2026-07-24 GENERALISE check
a2fd4df 2026-07-22 FORGE-G2: runtime-SMASK — smask as jitted argument
f4e2a85 2026-07-20 Forge run, HPC sync to accre
3a2ed4b 2026-07-16 update for HPC sync
423ca12 2026-07-13 reassess floor
d587b29 2026-07-13 HPC sync
27dbcab 2026-07-13 machinery: MargJax tracked, BLAS thread hazard closed, floor fix provisional
6bec3d6 2026-07-12 criterion-v2.1: soft loop adopted, onset retracted, both purity layers rejected
fb277db 2026-07-12 criterion-v2: null family, floor-as-function, purity layer pre-registered
db6ff05 2026-07-12 criterion-v1: three-layer certification criterion, fitted floor, campaign consolidation
6243ae7 2026-07-11 accre branch push
634aab8 2026-07-09 backlog of investigationns
06826cd 2026-06-25 prong2: Stage A.3 robustness
c23a42e 2026-06-25 prong2: N x t_c map (Stage A.3)
49a5ab0 2026-06-25 prong2: probe Farr's coherence transition (Stage A.2)
```

Consequence: for `wheel_hasasia_examples/` there is no earlier version to compare with and no way to recover a file if it is overwritten. The index holds one snapshot (whatever was last `git add`ed).

## 0.3 Other branches

Date of the last commit touching each directory, per branch (after `git fetch`):

| Branch | wheel_hasasia_examples | PAPER_figures | CW_transition |
|---|---|---|---|
| `MM_playground` (local) | none | 2026-07-23 `863d5cd` | 2026-07-24 `a8d3d97` |
| `origin/MM_playground` | none | 2026-07-23 `863d5cd` | 2026-07-24 `a8d3d97` |
| `master`, `origin/master` | none | none | none |
| `add-sampling-notebook` (local and remote) | none | none | none |
| `stash@{0}`, `stash@{1}` | none | none | none |

No branch holds a version of any key file that is newer than, or different from, the working tree. For the key files that are in git at all (`stagec_p2b.py`, `best_distances.txt`: last commit `634aab8`, 2026-07-09), the working tree equals HEAD.

## 0.4 Other copies on this machine

Searched by name: `/home` (only `/home/mattm` is readable; 12 other users' home directories returned "Permission denied"), `/opt`, `/srv`, `/mnt`, `/media`, `/tmp`, `/var/tmp`. One disk (`/`, 1.9 TB); no separate scratch or project filesystem is mounted. Checked for `.ipynb_checkpoints` copies and VS Code local-history entries of the notebook and key scripts: none.

| File | Path (under `/home/mattm/projects/HSYMT/`) | Size | Modified | sha256 (first 16) | Other copies |
|---|---|---|---|---|---|
| active_wheel_demo_v4.ipynb | wheel_hasasia_examples/ | 3,313,397 | 2026-08-17 02:18 | d2f2ee6376d9f7d9 | none |
| active_wheel_demo_v4.pdf | wheel_hasasia_examples/ | 3,203,812 | 2026-08-18 20:04 | f0cc133bc69b0527 | **`/home/mattm/soft/active_wheel_demo_v4.pdf`**, 1,897,354 bytes, 2026-08-11 01:03, d7907086ecf9be6b — differs, **older** |
| ACTIVE_WHEEL_REPORT.md | wheel_hasasia_examples/ | 338,328 | 2026-08-21 01:54 | 862cda2c5cc7c0ab | none |
| job43c_rebooking.py | wheel_hasasia_examples/ | 45,170 | 2026-08-21 00:56 | 597c23b010a3252c | none |
| job47c_heatmap.py | wheel_hasasia_examples/ | 10,738 | 2026-09-16 22:51 | e622c21034f1ea44 | none (untracked) |
| job47d_contour.py | wheel_hasasia_examples/ | 8,315 | 2026-09-16 23:05 | 52a7e0541bb4f119 | none (untracked) |
| job46_partB_figure.py | PAPER_figures/ | 6,871 | 2026-09-16 01:02 | b258bfd6ec2e7a8b | none |
| job41_displaced_excess.py | wheel_hasasia_examples/ | 13,417 | 2026-08-19 00:19 | 390c46830df1a3ec | none |
| job29_joint_bin_gate.py | wheel_hasasia_examples/ | 17,864 | 2026-08-17 15:40 | 9c1c5aafcc188ae2 | none |
| job30_middle_tier.py | wheel_hasasia_examples/ | 22,512 | 2026-08-17 15:54 | f3db0b610d16fa8b | none |
| job40_final_batch.py | wheel_hasasia_examples/ | 21,606 | 2026-08-19 00:02 | 28af404f1001fef8 | none |
| fisher_nested_fast.py | CW_transition/ | 8,750 | 2026-09-25 12:10 | a78560070a11ed95 | none (**untracked**, never committed) |
| stagec_p2b.py | CW_transition/ | 10,013 | 2026-07-02 21:26 | f6f3a98fecf67797 | none; equals HEAD |
| job45_t2_parallax_census.csv | wheel_hasasia_examples/ | 7,194 | 2026-09-03 21:02 | 465408329ea54965 | none (untracked) |
| best_distances.txt | CW_transition/ | 8,652 | 2026-07-03 00:15 | c55d5e08e4f06040 | none; equals HEAD |
| `*inventory20*`, `*SEED50*` | — | — | — | — | no such file (see top) |

For every key file that is staged, the working-tree copy equals the staged copy (git status shows `A `, not `AM`). So the files Jobs 49 and 50 read are the only versions that exist here.

**Not checkable from here:** the ACCRE cluster (HEAD's subject line and several commits mention HPC syncs; the 4 unfetched remote commits came from another clone), any laptop, Google Drive or Overleaf copy, and other users' home directories.

## 0.5 Notebook versions

Notebooks named `*wheel_demo*` anywhere on the machine (all in `wheel_hasasia_examples/`):

| File | Modified | Size | Code cells | Code characters |
|---|---|---|---|---|
| wheel_demo_v2.ipynb | 2026-08-06 23:21 | 2.4 MB | 10 | 32,744 |
| wheel_demo_v3.ipynb | 2026-08-10 16:09 | 3.1 MB | 10 | 53,806 |
| wheel_demo_v4.ipynb | 2026-08-10 19:15 | 3.6 MB | 12 | 77,735 |
| **active_wheel_demo_v4.ipynb** | 2026-08-17 02:18 | 3.3 MB | 32 | 543,136 |

- v4 is the highest. Nothing is named v5, final or paper. There are **two notebooks called v4**: `wheel_demo_v4.ipynb` (small, a week older) and `active_wheel_demo_v4.ipynb` (seven times more code). Only the `active_` one contains `inventory20`, `SEED50`, `FP_ALL`; it is the one every job script reads.
- Other notebooks in the directory (`one_pulsar_six_worlds.ipynb` 2026-09-25, `step1…step6_*.ipynb` 2026-09-03) are teaching notebooks, not the machine.
- The notebook's last save (08-17) is **earlier** than `job43c_rebooking.py` (08-21) and the report (08-21). `collect36m`, `inv_rate` and `book` are **not in the notebook at all**; they exist only in `job43c_rebooking.py`. The headline machine is therefore notebook cells plus that script, not the notebook alone.

**How the job scripts get notebook code: by executing cell source, not importing.** `job49_common.py:49-81` (same pattern in `job47c_heatmap.py:29-62`, `job47d_contour.py:26-48`):

```
nb = json.load(open("active_wheel_demo_v4.ipynb"))
...
        exec(compile(src, f"<v4-cell-{i}{':'+tag if tag else ''}>", "exec"), NS)
...
    for lo, hi in ((485, 509), (661, 689), (712, 756)):
        exec(compile("\n".join(lines43[lo-1:hi]), f"<43c:{lo}-{hi}>", "exec"), NS)
```

So the scripts run whatever the cell source is today; the saved outputs play no part. The line spans of `job43c_rebooking.py` are hard-coded: if that file is ever edited above line 756, the loaders silently pick up the wrong lines.

**Do saved outputs agree with what the cells produce now?** Four checks, all pass:

| Saved in notebook | Produced now | Where |
|---|---|---|
| Cell 3 output: "18,727 discrete sources above the split" | `N_DISC` = 18,727 | `job50_frames.py` (uncommitted, 2026-11-11), log `job50_logs/frames.log` line (3) |
| Cell 3 output: "CW upper-limit rejections during the draw: 4" | `N_REJ` = 4 | `job49_common.load()` run 2026-11-11 (uncommitted) |
| Cell 14 output: "50- 100 nHz: Earth 232 pulsar 44" | 232 / 44 | `job50_data/B2_pileup.json` (Job 50 B2, uncommitted, 2026-11-11) |
| `fig_v1_population.png`, `fig_v2_curves.png`, `fig_seed_robustness.png` (staged copies) | regenerated by today's loader run at 18:39–18:40; git still shows them unchanged from the staged copies, i.e. byte-identical | `job49_common.load()` |

Side effect to know about: every `load()` rewrites those three figure files in place. Today the bytes are identical, so nothing was lost, but the modification times now read 2026-11-11.

## 0.6 Figures

| Figure | File | Modified | Script | Script modified | Similar names elsewhere | Flag |
|---|---|---|---|---|---|---|
| Fisher, fixed strain | PAPER_figures/job46_info_fixed_per_source_strain.png | 09-16 01:08 | PAPER_figures/job46_partA_figure.py | 09-16 01:08 | none | input `CW_transition/prong2_results.npz` is dated 06-25, its generator `prong2_transition.py` 08-06: **data older than its script** |
| Fisher, fixed total power | PAPER_figures/job46_info_fixed_total_power.png | 09-16 01:08 | same | 09-16 01:08 | none | same flag |
| (also) surviving fraction | PAPER_figures/job46_info_surviving_fraction.png | 09-16 01:08 | same | 09-16 01:08 | none | same flag |
| Job 48b confined bins | PAPER_figures/job48b_turnover_per_bin_fig1style.png, job48b_collapse_Beff.png | 09-25 12:35 | PAPER_figures/job48b_plot.py (data: job48b_run.py 12:10 → job48_data/job48b_results.npz 12:12) | 09-25 12:35 | Job 48 (not b) versions: job48_figA_turnover_moves[_FULLMARG].png, job48_figB_collapse[_FULLMARG].png, 09-25 11:28, from job48_analyse.py | consistent; four older "48" figures sit beside the two "48b" ones |
| fig_phase_jitter_from_s.png | **not found** | — | **none found** | — | none | **missing figure, no script** |
| Job 47c heatmap | wheel_hasasia_examples/fig47c_criterion_heatmap.png | 09-16 23:00 | job47c_heatmap.py | 09-16 22:51 | fig47b_criterion_heatmap.png (22:25, job47b_heatmap.py 22:20) | consistent |
| Job 47d contour | wheel_hasasia_examples/fig47d_criterion_heatmap.png | 09-16 23:06 | job47d_contour.py | 09-16 23:05 | as above | consistent. Job 49 found 47d reproduces only with seeds taken from the rounded CSV grid |
| Share vs frequency (Job 46 B) | PAPER_figures/job46_pileup_histogram.png, job46_pileup_histogram_arcs15.png | 09-16 01:02 | PAPER_figures/job46_partB_figure.py | 09-16 01:02 | none | **its input `partB_sky.npz` no longer exists**: `job46_partB_extract.py:37` wrote it to a session scratch directory under `/tmp` that has been deleted. The figure cannot be remade without rerunning the extract |
| Share vs frequency, two-term | PAPER_figures/job50_share_vs_frequency_two_term.png | 11-11 16:48 | wheel_hasasia_examples/job50_partF.py | 11-11 16:47 | none | consistent |
| Corner vs amplitude (43C T3) | **no plot exists** | — | — | — | — | numbers only, in the report |
| job49_D1_heatmaps.png | wheel_hasasia_examples/ | 11-11 13:03 | job49_D1.py | 11-11 13:01 | none | consistent |
| job49_partC_hist.png | wheel_hasasia_examples/ | 11-11 12:52 | job49_partC.py | 11-11 12:52 | none | consistent |
| job50_fraction_vs_K.png | wheel_hasasia_examples/ | 11-11 16:43 | job50_analysis.py | **11-11 18:41** | none | **figure older than script.** The 18:41 edit only added an output-name suffix for the second mass convention; the plotting code is unchanged |

All dates 2026. No figure has a same-named copy in another directory.

## 0.7 Missing reports

- **Job 44 (clean-world phase-up): nothing on disk.** No file named `job44*` or `JOB44*`; the strings "JOB 44", "Job 44", "job44" appear in no `.md`, `.py`, `.txt` or `.log` file in the repository; `ACTIVE_WHEEL_REPORT.md` goes from Job 43C (line 6234) to Job 42B sections and ends at line 6559 with no Job 44 heading.
- **Job 45: scripts and tables, no report.** In `wheel_hasasia_examples/`, all dated 2026-09-03, all untracked: `job45_preflight.py` (20:52), `job45_t3_mcgrath.py` (21:01), `job45_t2_parallax_census.py` and `.csv` (21:02), `job45_tasks.py` (21:06), `job45_t3_scoreboard.csv` (21:01), `job45_t1_golden_pairs.csv` (22:17), `job45_t2_tier_revision.csv` (22:19), `fig45_mcgrath_toggle.png` (21:01). No `JOB45*.md`, no Job 45 section in the report, no log file. Job 45 also saved its machine state to a scratch directory that no longer exists.
- Jobs 46, 47 and 48 likewise have no section in `ACTIVE_WHEEL_REPORT.md` (last modified 2026-08-21) and no standalone report file; only scripts, tables and figures.

## 0.8 What Matt should look at by eye

1. **Commit `wheel_hasasia_examples/`.** Three months of work, including the booked headline code, exists as one uncommitted copy on one disk that is 92% full.
2. **Is there a copy on ACCRE or another machine?** I cannot see it. The remote received 4 commits from elsewhere after this clone's last fetch.
3. **`fig_phase_jitter_from_s.png`**: where was it made? It is not here.
4. **Job 44**: was it run in a session whose output was never saved here?
5. **`wheel_demo_v4.ipynb` versus `active_wheel_demo_v4.ipynb`**: confirm the small one is a dead ancestor.
6. **`/home/mattm/soft/active_wheel_demo_v4.pdf`** (08-11) is an older export than the one in the repository (08-18); neither was exported from the final notebook-plus-`job43c` state.
7. **Job 46 Part B figure** cannot be regenerated as is (input file deleted with its scratch directory).
8. **Job 46 Part A figures** read `prong2_results.npz` from June; its generating script changed in August.

## Decision on "stop there"

The brief says to stop if this directory is not the latest version. Nothing found says it is not: there is no newer copy, branch, stash or notebook. I am therefore continuing with Parts 1–4, on the understanding that "latest" here means "the only version on this machine", not "verified against version control".
