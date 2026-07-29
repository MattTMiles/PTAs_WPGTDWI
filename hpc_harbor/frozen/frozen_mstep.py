#!/usr/bin/env python
"""PHOENIX -- FROZEN-vs-LIVE M-step at the feasible rung.

Discharges the pre-registered follow-up authorized at GLACIER capstone S4.15.1 item 2:
  "PRE-REGISTERED FOLLOW-UP (authorized now, gated behind the ladder readback):
   frozen-vs-live M-step comparison at ONE rung (cheapest: refit-once-then-freeze at
   the feasible sky rung) IF the ladder shows a rung where frozen templates would have
   held the drain."

GATE CONDITION MET (2026-07-29, from the now-complete ladder readback S4.24 + SUMMIT
S0.1): at r13p9/none/T30 the HOLD cells take a real regression at i1 and then park
0.03-0.30 dex below or at baseline. There IS a regression for freezing to act on, and
it is at the rung the boundary-mapping mission is centered on. Condition satisfied.

WRAPPER, NOT A FORK (the chorus-C7 / COMPASS / SUMMIT-D2 pattern). glacier_loop.py is
imported STOCK and runtime-patched. Nothing in hpc_harbor/glacier/ is edited.

THE DIAL (declared): PER-SLOT FREEZE-AFTER-FIRST-FIT.
  A fed slot is fit by the stock `mstep_quadratic` on the iteration it is FIRST fed;
  on every later iteration its M-step values are HELD at that first fit. Everything
  else stays live -- the frontier still promotes, the drain still refits, floors still
  refit on policy, the scoreboard still scores. Only the fed members' (fgw, mc) stop
  moving.
  At the feasible rung the circular arm feeds exactly one member at i0, so this
  reduces EXACTLY to the pre-registered "refit-once-then-freeze". The per-slot form is
  the strictly-more-correct generalisation for any cell that late-feeds, and it is
  what is coded so a late feed cannot silently produce a half-frozen template.

LIVE ARM: INHERITED, NEVER RE-RUN.
  GLACIER_results/gl2_r13p9_cell_none_s{0..3}_T30_lit_i{0..5}__dgx03_NVIDIAA100-SXM4-80GB.npz
  24/24 banked, complete. The frozen arm runs on the SAME LANE (dgx03, A100-80GB,
  cpus-per-task=8). This is not a convenience: the host systematic measured in the
  noise-draw hazard is 1.72 sigma and is NOT common-mode, while the wander this
  experiment measures is 0.01-0.44 dex. A cross-lane frozen-vs-live comparison would
  not be readable. If dgx03 is unavailable the run WAITS; it does not migrate.

BANKS: FROZEN_results/ (own directory -- two provenances never share one, the g1
hazard). Stems are the STOCK stems; the directory IS the arm key, exactly as SUMMIT's
D3 rungs use d3_ext30/ d3_ext100/ with stock stems.

FRESH-DIR DISCIPLINE: this driver REFUSES to start a cell whose banks already exist.
It does not skip-on-exist. That is the SPARK-3C C3 lesson recorded verbatim in
SPARK3_second_arrow.md S5.4 -- "the Jul-18 array was NOT that re-cut; skip-on-exist
reduced it to a five-rung gap-fill on the wrong floor lane". A frozen cell resumed
from a live checkpoint would silently re-fit at the resume iteration and produce a
half-frozen template with no marker. Refuse, don't resume.

MODES:
  gate  -- SG-F1 wrapper inertness. Freeze DISABLED, i0 only, column-compared to the
           banked live cell on the same lane. Bar = SG-D2a's amended bar (deterministic
           columns exact; Fisher-derived columns decision-invariant + 0.1 sanity cap).
  cell  -- the frozen cell (--sky 0..3, circular arm, r13p9, T=30, w=1).
"""
import os, sys, glob, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "glacier"))
sys.path.insert(0, os.path.join(HERE, "..", "ignite"))
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")

import numpy as np

FROZEN_OUT_REQUIRED = "FROZEN_results"
LIVE_DIR = "/home/mattm/projects/HSYMT/GLACIER_results"
LIVE_LANE = "dgx03_NVIDIAA100-SXM4-80GB"


def _check_out():
    out = os.environ.get("GLACIER_OUT", "")
    if FROZEN_OUT_REQUIRED not in out:
        print(f"REFUSED: GLACIER_OUT={out!r} does not point at {FROZEN_OUT_REQUIRED}. "
              "The frozen arm never banks into GLACIER_results (two provenances, one "
              "dir = the g1 hazard). Set APPTAINERENV_GLACIER_OUT.", flush=True)
        sys.exit(3)
    os.makedirs(out, exist_ok=True)
    return out


def _check_lane():
    """The frozen arm is only readable against the live arm on the live arm's lane."""
    import glacier_pop as GPM
    lane = GPM.lane_tag() if hasattr(GPM, "lane_tag") else None
    if lane is None:
        import glacier_loop as GL
        lane = GL.lane_tag()
    if lane != LIVE_LANE:
        print(f"REFUSED: lane={lane!r} but the banked live arm is on {LIVE_LANE!r}. "
              "The host systematic (1.72 sigma, noise-draw hazard) is not common-mode "
              "and the wander under measurement is 0.01-0.44 dex; a cross-lane "
              "frozen-vs-live comparison is not readable. Run on dgx03 A100-80GB.",
              flush=True)
        sys.exit(3)
    return lane


def _copy_evidence(out):
    """Evidence gate COPIED with a copy-record line (the SUMMIT convention). The holds
    are the same green record on the same tree -- copied, never re-declared."""
    rec = []
    for src in [f"{LIVE_DIR}/HOLDS_CLEARED"] + sorted(glob.glob(f"{LIVE_DIR}/gl2_ladder_gates__*.npz")):
        dst = os.path.join(out, os.path.basename(src))
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.copy2(src, dst)
            rec.append(os.path.basename(src))
    if rec:
        with open(os.path.join(out, "EVIDENCE_COPY_RECORD"), "a") as fh:
            fh.write(f"copied from {LIVE_DIR}: {', '.join(rec)}\n")
    return rec


# ============================================================
# THE PATCH
# ============================================================
_FROZEN = {"on": False, "slots": {}, "log": []}


def reset_freeze(on=True):
    _FROZEN["on"] = bool(on)
    _FROZEN["slots"] = {}
    _FROZEN["log"] = []


def _patch_mstep():
    """Wrap glacier_loop.mstep_quadratic with per-slot freeze-after-first-fit."""
    import glacier_loop as GL
    if getattr(GL.mstep_quadratic, "_phoenix_frozen", False):
        return
    orig = GL.mstep_quadratic
    NP = GL.NP_SRC
    NAX = len(GL.MSTEP_AXES)

    def mstep_frozen(marg_fn, src0, fed_idx, scale, **kw):
        if not _FROZEN["on"]:
            return orig(marg_fn, src0, fed_idx, scale, **kw)
        fed_idx = list(fed_idx)
        new = [k for k in fed_idx if k not in _FROZEN["slots"]]
        n_eval = 0
        widths = np.full((len(fed_idx), NAX), np.nan)
        if new:
            # fit ONLY the slots being fed for the first time -- stock code, stock bar
            src_hat, w_new, n_eval = orig(marg_fn, src0, new, scale, **kw)
            for j, k in enumerate(new):
                _FROZEN["slots"][k] = src_hat[NP * k: NP * (k + 1)].copy()
                _FROZEN["log"].append((int(k), "fit"))
                widths[fed_idx.index(k)] = w_new[j]
        for k in fed_idx:
            if k not in new:
                _FROZEN["log"].append((int(k), "held"))
        out = src0.copy()
        for k in fed_idx:
            out[NP * k: NP * (k + 1)] = _FROZEN["slots"][k]
        return out, widths, n_eval

    mstep_frozen._phoenix_frozen = True
    GL.mstep_quadratic = mstep_frozen
    print("[FROZEN] mstep_quadratic patched: per-slot freeze-after-first-fit", flush=True)


# ============================================================
# MODES
# ============================================================
def mode_cell(sky, arm="none", rung="r13p9", t=30, wscale=1.0, n_iter=None):
    out = _check_out()
    import glacier_loop as GL
    import glacier_pop as GPM
    _check_lane()
    _copy_evidence(out)

    marker = f"{GPM.OUT}/HOLDS_CLEARED"
    gates2 = glob.glob(f"{GPM.OUT}/gl2_ladder_gates__*.npz")
    if not (os.path.exists(marker) and gates2):
        print(f"REFUSED: need {marker} AND a gl2_ladder_gates bank.", flush=True)
        return 3

    stem = f"gl2_{rung}_cell_{arm}_s{sky}_T{t}_{GL.TIER}"
    pre = glob.glob(f"{GPM.OUT}/{stem}_i*__*.npz")
    if pre:
        print(f"REFUSED (fresh-dir discipline): {len(pre)} bank(s) for {stem} already "
              f"exist in {GPM.OUT}. The frozen arm never resumes from an existing "
              "checkpoint -- a resumed cell re-fits at the resume iteration and yields "
              "a half-frozen template with no marker (SPARK-3C S5.4, skip-on-exist). "
              "Move or delete them and relaunch.", flush=True)
        return 3

    _patch_mstep()
    reset_freeze(on=True)
    n_iter = GL.N_ITER if n_iter is None else int(n_iter)
    print(f"[FROZEN] cell {stem} arm=FROZEN n_iter={n_iter} lane={GL.lane_tag()}", flush=True)
    rc = GL.run_cell(arm, sky, rung=rung, t_label=t, wscale=wscale, n_iter=n_iter)
    fits = [k for k, w in _FROZEN["log"] if w == "fit"]
    held = [k for k, w in _FROZEN["log"] if w == "held"]
    print(f"[FROZEN] freeze ledger: {len(fits)} first-fits (slots {sorted(set(fits))}), "
          f"{len(held)} held slot-iterations", flush=True)
    return rc


def mode_gate(sky=0, arm="none", rung="r13p9", t=30):
    """SG-F1 -- wrapper inertness with the freeze DISABLED. One iteration, columns
    compared to the banked live cell on the same lane."""
    out = _check_out()
    import glacier_loop as GL
    import glacier_pop as GPM
    lane = _check_lane()
    _copy_evidence(out)

    sub = os.path.join(out, "gate_f1")
    os.makedirs(sub, exist_ok=True)
    GPM.OUT = sub
    GL.OUT = sub
    for f in glob.glob(f"{sub}/*.npz"):
        os.remove(f)                      # the gate is always a fresh cut
    for src in [f"{LIVE_DIR}/HOLDS_CLEARED"] + sorted(glob.glob(f"{LIVE_DIR}/gl2_ladder_gates__*.npz")):
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(sub, os.path.basename(src)))

    _patch_mstep()
    reset_freeze(on=False)                # INERT: stock path through the wrapper
    print(f"[FROZEN-GATE] SG-F1 inertness, freeze OFF, i0 only, lane={lane}", flush=True)
    GL.run_cell(arm, sky, rung=rung, t_label=t, n_iter=1)

    stem = f"gl2_{rung}_cell_{arm}_s{sky}_T{t}_{GL.TIER}"
    got = f"{sub}/{stem}_i0__{lane}.npz"
    ref = f"{LIVE_DIR}/{stem}_i0__{lane}.npz"
    if not (os.path.exists(got) and os.path.exists(ref)):
        print(f"SG-F1 FAIL: missing bank\n  got={got} {os.path.exists(got)}"
              f"\n  ref={ref} {os.path.exists(ref)}", flush=True)
        return 1

    A, B = np.load(got, allow_pickle=True), np.load(ref, allow_pickle=True)
    EXACT = ["fed_mask", "n_res", "a_bg", "floor", "zf", "n_cert", "cert_idx",
             "wrong_cert", "promote_events", "on_true"]
    FISHER = ["conc_ratio", "widths", "F_ii", "loc_area_deg2"]
    ok, lines = True, []
    keys = [k for k in A.files if k in B.files]
    for k in keys:
        a, b = A[k], B[k]
        if a.dtype == object or b.dtype == object:
            same = str(a) == str(b)
            lines.append(f"  {k:20s} object  {'==' if same else 'DIFF'}")
            if not same and k in EXACT:
                ok = False
            continue
        try:
            a = np.asarray(a, float); b = np.asarray(b, float)
        except (TypeError, ValueError):
            continue
        if a.shape != b.shape:
            lines.append(f"  {k:20s} SHAPE {a.shape} vs {b.shape}")
            ok = False
            continue
        d = np.abs(a - b)
        d = np.nanmax(d) if d.size and np.isfinite(d).any() else 0.0
        tag = "EXACT" if k in EXACT else ("FISHER" if k in FISHER else "other")
        lines.append(f"  {k:20s} {tag:6s} max|d| = {d:.6e}")
        if k in EXACT and d != 0.0:
            ok = False
        if k in FISHER and d > 0.1:
            ok = False
    print("SG-F1 column compare (frozen wrapper, freeze OFF) vs banked live cell:", flush=True)
    print("\n".join(lines), flush=True)

    # decision invariance on the feed bar -- the SG-D2a amended bar, verbatim in form
    if "conc_ratio" in keys:
        ca = np.asarray(A["conc_ratio"], float)
        cb = np.asarray(B["conc_ratio"], float)
        m = np.isfinite(ca) & np.isfinite(cb)
        crossed = int(np.sum((ca[m] < 0.5) != (cb[m] < 0.5)))
        marg = float(np.nanmin(np.abs(cb[m] - 0.5))) if m.any() else float("nan")
        wob = float(np.nanmax(np.abs(ca[m] - cb[m]))) if m.any() else 0.0
        print(f"  decision invariance: slots crossing the 0.5 feed bar = {crossed} "
              f"(margin {marg:.4f} vs wobble {wob:.4f})", flush=True)
        if crossed:
            ok = False
    print(f"SG-F1 {'PASS' if ok else 'FAIL'}", flush=True)
    return 0 if ok else 1


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "cell"])
    ap.add_argument("--sky", type=int, default=0)
    ap.add_argument("--arm", default="none")
    ap.add_argument("--rung", default="r13p9")
    ap.add_argument("--t", type=int, default=30)
    ap.add_argument("--wscale", type=float, default=1.0)
    ap.add_argument("--iters", type=int, default=None)
    a = ap.parse_args()
    if a.mode == "gate":
        return mode_gate(sky=a.sky, arm=a.arm, rung=a.rung, t=a.t)
    return mode_cell(sky=a.sky, arm=a.arm, rung=a.rung, t=a.t,
                     wscale=a.wscale, n_iter=a.iters)


if __name__ == "__main__":
    sys.exit(main() or 0)
