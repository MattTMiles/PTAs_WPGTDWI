"""LOTTERY tier-2 -- GPU validation of the channel-budget switch-on proxy.

Reuses CHORUS's exact machinery (chorus.build_chorus_problem / run_realisation) and the RECUT
floor+score+grade convention (CW_transition/recut_chorus.py, verbatim) on NEW low-eccentricity
mixture cells that tier-1 flagged in the DISAGREEMENT BAND (channel budget ON while the quoted
threshold rule is OFF) and at the channel KNEE.

The proxy under test:  n_active >= 30  <=>  the certified count clears the >1 bar (CONFIRMED).
Each cell has a DETERMINISTIC n_active (fixed n_ecc, fixed e), so the analytic channel verdict is
ON/OFF, not a probability; the GPU measures the grade. A distribution point's analytic
P(switch-on) is the n_ecc-weighted sum of its constituent cells' ON indicators, so validating the
cells validates the whole P(switch-on) surface.

VALIDATION CELLS (tier=lit, T=30 census; lit is where the CHORUS headline lives):
  NEW (GPU here):
    m3 e=0.25  n_active=34  channel ON   threshold OFF   <- decisive disagreement
    m3 e=0.20  n_active=31  channel ON   threshold OFF   <- decisive (barely over knee)
    m2 e=0.28  n_active=30  channel ON   threshold OFF   <- decisive at the knee
    m2 e=0.25  n_active=28  channel OFF  threshold OFF   <- knee, below
    m1 e=0.25  n_active=22  channel OFF  threshold OFF   <- concordant-low + Point-A constituent
    m1 e=0.45  n_active=29  channel OFF  threshold OFF   <- single-member knee, below
  CARRIED from CHORUS banks (no GPU; recut_chorus.npz):
    m3 e=0.70 (=m3e07)  n_active=109  ON  ON   CONFIRMED  <- concordant HIGH anchor
    m2 e=0.30 (=m2e03)  n_active=30   ON  ON   CONFIRMED  <- concordance x-check of the refit path
    m1 e=0.30 (=m1e03)  n_active=23   OFF OFF  below      <- concordant-low anchor

Distribution point A (the observer panel's steep ridge): mix f_ecc, e_char=0.25. A draw is
channel-ON iff n_ecc=3 (n_active=16+6 n_ecc), so P_chan = f_ecc^3, threshold always OFF. Its
constituent cells m1/m2/m3 e=0.25 (+ banked m0) let us read measured-P off the cell grades.

MODES:
  run       : GPU. 100 nullN + 20 sig per NEW cell, one node/host, cpus-per-task=8 (basis
              consistency). Skip-on-exist => checkpoint per draw, free resume.
  analyze   : CPU. refit floors (RECUT adopt), score, grade, calibration curve + validated table.

Usage (inside harbor container):
  harbor_py hpc_harbor/lottery/lottery_tier2.py run       # on a GPU node
  harbor_py hpc_harbor/lottery/lottery_tier2.py analyze    # login/CPU
"""
import os, sys, glob, time, socket, json
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
sys.path.insert(0, f"{HSYMT}/hpc_harbor/chorus")
sys.path.insert(0, f"{HSYMT}/hpc_harbor/atlas")
sys.path.insert(0, f"{HSYMT}/CW_transition")

import chorus as CH
import atlas_harmonics as H

LOT_OUT = f"{HSYMT}/LOTTERY_results"
REPORTS = f"{HSYMT}/reports"
PREFIX = "lot"

# ---- channel map c(e) (identical to tier 1, from the real harmonic algebra) --
_NS = np.arange(1, CH.N_HARM + 1)
def c_of_e(e):
    if e <= 0.0:
        return 1
    g = H.g_n(_NS, float(e))
    return int(np.sum(g >= CH.G_ACTIVE_REL * g.max()))
def n_active_of(n_ecc, e):
    return 13 + n_ecc * c_of_e(e) + (3 - n_ecc) * 1

# ---- the NEW validation cells: (n_ecc, e, edist_label) ----------------------
#   edist_label feeds CH.E_FIXED so CH.e_values(n_ecc,label) returns full(n_ecc, e).
NEW_CELLS = [
    (3, 0.25, "L025"),
    (3, 0.20, "L020"),
    (2, 0.28, "L028"),
    (2, 0.25, "L025"),
    (1, 0.25, "L025"),
    (1, 0.45, "L045"),
]
E_LABELS = {"L020": 0.20, "L025": 0.25, "L028": 0.28, "L045": 0.45}
# BASIS CONTROL: CHORUS's OWN cells re-run on THIS host/basis. If their refit floor+count match
# the banked (dgx) values in recut_chorus.npz, the GWB basis is robust across hosts and every
# cross-basis comparison below is licensed; if they diverge ~nat, the near-bar grade is
# host-sensitive and that systematic -- not the proxy -- dominates. ("e03" is CHORUS's own label.)
CONTROL_CELLS = [(1, 0.30, "e03"), (2, 0.30, "e03")]
ALL_CELLS = NEW_CELLS + CONTROL_CELLS
TIER = "lit"

N_NULLN_LOT = 100         # match CHORUS N_NULLN for floor comparability
N_SIG_SKY   = 4           # 4 reps x 5 skies = 20 sig/cell (>= 20 required)
CI_BASE     = 900         # distinct cell index base so seeds never collide with CHORUS
CHAN_THRESH = 30          # channel-budget switch: n_active >= 30

# ---- RECUT floor/score/grade convention (recut_chorus.py, verbatim) ---------
from scipy.stats import gumbel_r
ALPHA = 0.05
Z_ALPHA = -np.log(-np.log(1.0 - ALPHA))
C_SD = 2.80
QBAR = 0.9
ZF_GATE = 0.20
BOOT = 4000
SEED = 20260717

def _offender(dlnL, lnK, qmax):
    dlnL = np.nan_to_num(np.asarray(dlnL, float), posinf=1e30)
    m = (dlnL > np.asarray(lnK)) & (np.asarray(qmax) > QBAR)
    return float(dlnL[m].max()) if m.any() else 0.0

def gumbel_floor(x):
    x = np.asarray(x, float); mu, beta = gumbel_r.fit(x)
    return float(mu + beta * Z_ALPHA), float(C_SD * beta / np.sqrt(len(x)))

def emp_quantile(off):
    return float(np.quantile(np.asarray(off, float), 1.0 - ALPHA))

def boot_sd(off):
    rng = np.random.default_rng(SEED); off = np.asarray(off, float); n = len(off)
    return float(np.std([np.quantile(rng.choice(off, n, replace=True), 1.0 - ALPHA)
                         for _ in range(BOOT)]))

def adopt(off, zf):
    gu, gsd = gumbel_floor(off)
    if zf <= ZF_GATE:
        return float(gu), float(gsd), "gumbel"
    return emp_quantile(off), boot_sd(off), "emp_q95"


# ============================================================================
# GPU run mode
# ============================================================================
def _plan():
    """Deterministic work list: nullN (no_cw) + sig per NEW cell."""
    plan = []
    for idx, (k, e, lab) in enumerate(ALL_CELLS):
        ci = CI_BASE + idx
        # nulls (nullN only -- the RECUT floor uses oN alone)
        for rep in range(N_NULLN_LOT):
            gid = CH.SKIES[rep % len(CH.SKIES)]
            ns = 9_000_000 + 10_000 * ci + rep
            plan.append(dict(kind="fnullN", n_ecc=k, edist=lab, tier=TIER, geo_id=gid,
                             noise_seed=ns, dist_seed=ns + 10_000_000, no_cw=True))
        # signal
        for si, gid in enumerate(CH.SKIES):
            for rep in range(N_SIG_SKY):
                ns = 9_500_000 + 10_000 * ci + 100 * si + rep
                plan.append(dict(kind="sig", n_ecc=k, edist=lab, tier=TIER, geo_id=gid,
                                 noise_seed=ns, dist_seed=ns + 10_000_000))
    return plan

def _path_of(e):
    return f"{LOT_OUT}/{PREFIX}_{e['kind']}_{CH.cell_tag(e['n_ecc'], e['edist'], e['tier'])}_g{e['geo_id']}_n{e['noise_seed']}.npz"

def mode_run():
    # register the lottery eccentricity labels + redirect output into LOTTERY_results
    CH.E_FIXED.update(E_LABELS)
    CH.OUT = LOT_OUT
    os.makedirs(LOT_OUT, exist_ok=True)
    host = socket.gethostname()
    print(f"[LOTTERY tier-2 run] host={host} cpus={os.environ.get('SLURM_CPUS_PER_TASK','?')} "
          f"start={time.strftime('%FT%T')}", flush=True)
    jax, jnp, C, B1, TE, F, IG = CH._import_stack()
    print(f"jax {jax.__version__} {jax.devices()}", flush=True)

    # --- GWB basis convention: match CHORUS's T=30 realisations EXACTLY ---------
    # The frozen b1_L_gwb bank is the T=15 baseline (3248^2). At T=30 the extended GWB is
    # 5336^2, so the bank shape-mismatches and NoiseDrawer would RAISE. CHORUS's own T=30
    # nulls/signals were drawn on 2026-07-13 BEFORE that bank existed, via the recompute
    # fallback at cpus-per-task=8 (OpenBLAS 8 threads). We reproduce that convention: force
    # the recompute path so every LOTTERY cell shares ONE host's 8-thread GWB basis (nulls and
    # signals together). The certified COUNT is robust to the basis rotation (RECUT §0.1: the
    # certified/detected sets are unchanged under the thread-count rotation) — which is exactly
    # why counts from this host are comparable to CHORUS's banked counts on the calibration.
    _orig_lgwb = C.load_or_build_L_gwb
    def _force_recompute(Phi_gwb, path=None, strict=False):
        return _orig_lgwb(Phi_gwb, path="/nonexistent/lottery_forces_recompute.npz", strict=False)
    C.load_or_build_L_gwb = _force_recompute
    print("[LOTTERY] NoiseDrawer forced to recompute L_gwb (T=30 extended; matches CHORUS "
          "T=30 convention, 8-thread basis)", flush=True)

    plan = _plan()
    todo = [e for e in plan if not os.path.exists(_path_of(e))]
    print(f"plan: {len(plan)} entries; already banked: {len(plan)-len(todo)}; to run: {len(todo)}",
          flush=True)
    if not todo:
        print("NOTHING TO DO (all banked)", flush=True)
        return 0
    skies = {g: (None if g < 0 else F.load_geo_skies([g])[0]) for g in set(CH.SKIES)}
    by_shape = {}
    for e in todo:
        by_shape.setdefault(e["n_ecc"], []).append(e)
    for n_ecc in sorted(by_shape):
        P = CH.build_chorus_problem(n_ecc)
        print(f"  [NoiseDrawer] {P['nd'].lgwb_prov}", flush=True)   # basis provenance -> log
        for e in by_shape[n_ecc]:
            t0 = time.time()
            tag, msg, _ = CH.run_realisation(
                P, C, F, jnp, e["kind"], e["n_ecc"], e["edist"], e["tier"],
                e["geo_id"], skies[e["geo_id"]], e["noise_seed"], e["dist_seed"],
                no_cw=e.get("no_cw", False), prefix=PREFIX)
            print(f"  {tag:40s} {msg}  ({time.time()-t0:.1f}s)", flush=True)
        del P
        import gc; gc.collect()
    print(f"RUN COMPLETE {time.strftime('%FT%T')}", flush=True)
    return 0


# ============================================================================
# CPU analyze mode
# ============================================================================
def _load_sig(tag):
    fs = sorted(glob.glob(f"{LOT_OUT}/{PREFIX}_sig_{tag}_g*_n*.npz"))
    D, K, Q, M, N = [], [], [], [], []
    for f in fs:
        d = np.load(f, allow_pickle=True)
        D.append(np.nan_to_num(np.asarray(d["dlnL_det"]), posinf=1e30))
        K.append(np.asarray(d["lnK"])); Q.append(np.asarray(d["qmax"]))
        M.append(np.asarray(d["mapk"])); N.append(np.asarray(d["n_true_grid"]))
    return dict(dlnL=np.array(D), lnK=np.array(K), qmax=np.array(Q),
                mapk=np.array(M), ntg=np.array(N), n=len(fs))

def _load_nulls(tag):
    fs = sorted(glob.glob(f"{LOT_OUT}/{PREFIX}_fnullN_{tag}_g*_n*.npz"))
    return np.array([_offender(np.load(f)["dlnL_det"], np.load(f)["lnK"],
                               np.load(f)["qmax"]) for f in fs])

def _score(S, floor):
    """recut_chorus.score, verbatim. Returns (corr_per_real, wrong_per_real, per-draw corr counts)."""
    cert = (S["dlnL"] > np.maximum(S["lnK"], floor)) & (S["qmax"] > QBAR)
    corr = cert & (S["mapk"] == S["ntg"])
    wrong = cert & (S["mapk"] != S["ntg"])
    per_draw = corr.sum(axis=1)                      # certified-correct count per realisation
    return float(corr.sum()) / S["n"], float(wrong.sum()) / S["n"], per_draw

def _load_carried():
    """Carry banked CHORUS cells from recut_chorus.npz (same convention, self-consistent basis)."""
    z = np.load(f"{REPORTS}/recut_chorus.npz", allow_pickle=True)
    cols = [str(c) for c in z["cols"]]; tab = z["table"]; tags = [str(t) for t in z["tags"]]
    grade = [str(g) for g in z["grade"]]
    want = {"m3e07_lit": (3, 0.70), "m2e03_lit": (2, 0.30), "m1e03_lit": (1, 0.30)}
    out = []
    for tg, (k, e) in want.items():
        i = tags.index(tg)
        row = {c: float(tab[i][j]) for j, c in enumerate(cols)}
        out.append(dict(tag=tg, n_ecc=k, e=e, n_active=n_active_of(k, e),
                        floor=row["floor"], err=row["err"], corr=row["corr"],
                        corr_lo=row["corr_lo"], n=int(row["n"]), zf=row["zf"],
                        grade=grade[i], carried=True))
    return out

def mode_analyze():
    print("=" * 100)
    print("LOTTERY tier-2 analyze -- refit floors (RECUT adopt), score, grade, calibrate")
    print("=" * 100)
    rows = []
    for (k, e, lab) in ALL_CELLS:
        tag = CH.cell_tag(k, lab, TIER)
        oN = _load_nulls(tag); S = _load_sig(tag)
        if len(oN) == 0 or S["n"] == 0:
            print(f"  {tag}: INCOMPLETE (nulls={len(oN)}, sig={S['n']}) -- skipped", flush=True)
            continue
        zf = float(np.mean(oN == 0.0))
        floor, err, est = adopt(oN, zf)
        cc, ww, per_draw = _score(S, floor)
        clo, _, _ = _score(S, floor + err)
        grade = "CONFIRMED" if clo > 1.0 else ("MARGINAL" if cc > 1.0 else "below")
        nact = n_active_of(k, e)
        chan_pred = "ON" if nact >= CHAN_THRESH else "OFF"
        thr_pred = "ON" if (e >= 0.5 or (k >= 2 and e >= 0.3)) else "OFF"
        # per-draw switch-on rate (measured P that a single realisation clears >1)
        p_meas = float(np.mean(per_draw > 1.0))
        p_err = float(np.sqrt(p_meas * (1 - p_meas) / S["n"]))
        is_control = (k, e, lab) in CONTROL_CELLS
        rows.append(dict(tag=(tag + " (hgx)" if is_control else tag), n_ecc=k, e=e,
                         n_active=nact, zf=zf,
                         floor=floor, err=err, est=est, corr=cc, corr_lo=clo, wrong=ww,
                         n=S["n"], grade=grade, chan_pred=chan_pred, thr_pred=thr_pred,
                         p_meas=p_meas, p_err=p_err, n_nulls=len(oN), carried=False,
                         is_control=is_control, base_tag=tag))
        print(f"  {tag:12s} nact={nact:3d} zf={zf:.2f} floor={floor:6.2f}±{err:.2f} [{est}] "
              f"count={cc:.2f}[{clo:.2f}] grade={grade:9s} chan_pred={chan_pred} thr_pred={thr_pred}",
              flush=True)

    carried = _load_carried()
    for r in carried:
        r["chan_pred"] = "ON" if r["n_active"] >= CHAN_THRESH else "OFF"
        r["thr_pred"] = "ON" if (r["e"] >= 0.5 or (r["n_ecc"] >= 2 and r["e"] >= 0.3)) else "OFF"

    # ---- BASIS CONTROL: my-host m1e03/m2e03 vs CHORUS banked (dgx) ----------
    print("\n" + "=" * 100)
    print("BASIS CONTROL -- CHORUS's own cells re-run on THIS host vs banked (dgx). If floor+count")
    print("match, cross-host counts are comparable; if they diverge, the near-bar grade is host-sensitive.")
    carried_by_tag = {r["tag"]: r for r in carried}
    basis = []
    for r in rows:
        if not r.get("is_control"):
            continue
        b = carried_by_tag.get(r["base_tag"])
        if not b:
            continue
        dfloor = r["floor"] - b["floor"]; dcount = r["corr"] - b["corr"]
        basis.append(dict(tag=r["base_tag"], mine_floor=r["floor"], bank_floor=b["floor"],
                          dfloor=dfloor, mine_count=r["corr"], bank_count=b["corr"], dcount=dcount,
                          mine_grade=r["grade"], bank_grade=b["grade"]))
        print(f"  {r['base_tag']:12s}  floor hgx={r['floor']:6.2f} vs dgx={b['floor']:6.2f} "
              f"(Δ={dfloor:+.2f})  count hgx={r['corr']:.2f} vs dgx={b['corr']:.2f} "
              f"(Δ={dcount:+.2f})  grade hgx={r['grade']} vs dgx={b['grade']}")
    maxdf = max((abs(x["dfloor"]) for x in basis), default=0.0)
    flips = [x["tag"] for x in basis if x["mine_grade"] != x["bank_grade"]]
    print(f"  => max |Δfloor| across hosts = {maxdf:.2f} nat; grade flips: {flips or 'none'}")
    print(f"  => VERDICT: {'basis-ROBUST (cross-host counts comparable)' if maxdf < 1.0 and not flips else 'basis-SENSITIVE (host systematic dominates near-bar grades)'}")

    all_rows = rows + carried
    # save bank
    os.makedirs(REPORTS, exist_ok=True)
    np.savez(f"{REPORTS}/lottery_tier2.npz",
             rows_json=json.dumps(all_rows),
             basis_json=json.dumps(basis),
             conventions=json.dumps(dict(N_NULLN=N_NULLN_LOT, N_SIG_SKY=N_SIG_SKY,
                                         ALPHA=ALPHA, ZF_GATE=ZF_GATE, QBAR=QBAR, TIER=TIER,
                                         chan_thresh=CHAN_THRESH)))
    print(f"\nbanked -> {REPORTS}/lottery_tier2.npz")

    # ---- validated table -------------------------------------------------
    print("\n" + "=" * 100)
    print(f"{'cell':12s} {'n_act':>5} {'chan':>4} {'thr':>4} | {'floor':>7} {'±':>5} "
          f"{'count':>6} {'@+bs':>6} {'grade':>10} | {'agrees channel?':>16}")
    for r in sorted(all_rows, key=lambda x: x["n_active"]):
        graded_on = r["grade"] in ("CONFIRMED",)
        agree = "YES" if (graded_on == (r["chan_pred"] == "ON")) else "*** NO ***"
        tail = " (carried)" if r.get("carried") else ""
        print(f"{r['tag']:12s} {r['n_active']:5d} {r['chan_pred']:>4} {r['thr_pred']:>4} | "
              f"{r['floor']:7.2f} {r['err']:5.2f} {r['corr']:6.2f} {r['corr_lo']:6.2f} "
              f"{r['grade']:>10} | {agree:>16}{tail}")
    return 0


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    sys.exit({"run": mode_run, "analyze": mode_analyze}[mode]())
