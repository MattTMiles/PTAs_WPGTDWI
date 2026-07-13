"""Emit SCHEMA-COMPATIBLE corrected twins of the four campaign banks.

WHY TWINS AND NOT OVERWRITES. `surface_floors.npz`, `surface_analysis.npz`, `ch_floors.npz` and
`ch_analysis.npz` are the artifacts the PUBLISHED reports (SURFACE_onset.md, CHORUS_mixed_pop.md,
FLOOR_FIX_provisional.md) read back against. Clobbering them would destroy the pre-fix record
that those reports cite, and the artifact-readback convention exists precisely to stop that. So
the corrected banks ship ALONGSIDE, under `*_recut.npz`, with the SAME SCHEMA as their originals
— so `surface_figures.py` / `chorus_figures.py` and any downstream reader can consume them by
changing one path, not their code.

Every corrected file carries the zero-fraction as a REQUIRED column (the convention), plus the
estimator actually used per cell, plus the pre-fix columns for side-by-side comparison.

Source of truth: reports/recut_surface.npz, reports/recut_chorus.npz (both readback-gated).

Run: python CW_transition/stage_recut_banks.py       (CPU, seconds)
"""
import os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda f: os.path.join(ROOT, "reports", f)


def main():
    S = np.load(R("recut_surface.npz"), allow_pickle=True)
    C = np.load(R("recut_chorus.npz"), allow_pickle=True)
    sc = [str(x) for x in S["cols"]]
    st = S["table"]
    sg = lambda k: st[:, sc.index(k)]
    cc = [str(x) for x in C["cols"]]
    ct = C["table"]
    cg = lambda k: ct[:, cc.index(k)]

    # ---------------- surface_floors_recut.npz (mirrors surface_floors.npz) ----------------
    O = np.load(R("surface_floors.npz"), allow_pickle=True)
    ocols = [str(x) for x in O["tab_cols"]]
    otab = O["table"]
    og = lambda k: otab[:, ocols.index(k)]
    # cells() order is identical in both (verified by the recut gates), so rows align 1:1
    assert np.allclose(og("h"), sg("h")) and np.allclose(og("T"), sg("T")), "row order drift"

    tab_cols = ["h", "T", "k", "fN", "fN_sd", "fN_emp", "fN_emp_sd", "fN_zerofrac",
                "floor_adopted", "floor_adopted_sd", "fALL", "fALL_sd"]
    rows = np.column_stack([
        sg("h"), sg("T"), sg("k"),
        sg("gumbel"), sg("gumbel_sd"),          # the Gumbel, kept for comparison
        sg("emp"), sg("emp_sd"),                # the empirical q95 + BOOTSTRAP sd (new)
        sg("zf"),                               # REQUIRED column
        sg("floor"), sg("err"),                 # what the convention actually adopts
        og("fALL"), og("fALL_sd"),
    ])
    np.savez(R("surface_floors_recut.npz"),
             alpha=float(S["alpha"]), zf_gate=float(S["zf_gate"]),
             boot=int(S["boot"]), seed=int(S["seed"]),
             tab_cols=np.array(tab_cols), table=rows,
             tiers=S["tiers"], struct=S["struct"], estimator=S["estimator"],
             cols=np.array(tab_cols),
             note=("ADOPTED floor per ANCHOR §4: floor_adopted = fN (Gumbel) where "
                   "fN_zerofrac <= 0.20, else fN_emp (empirical q95) with a BOOTSTRAP sd. "
                   "Schema mirrors surface_floors.npz. Originals untouched."))
    print(f"  surface_floors_recut.npz    {rows.shape[0]} cells")

    # ---------------- surface_analysis_recut.npz (mirrors surface_analysis.npz) ------------
    A = np.load(R("surface_analysis.npz"), allow_pickle=True)
    acols = [str(x) for x in A["cols"]]
    atab = A["table"]
    ag = lambda k: atab[:, acols.index(k)]

    keys = ["h", "T", "k", "n", "fN", "fN_sd", "fN_zerofrac", "fN_emp", "fN_emp_sd",
            "floor_adopted", "floor_adopted_sd", "corr", "corr_lo", "corr_hi", "wrong",
            "corr_se", "corr_prefix", "corr_lo_prefix"]
    tab = np.column_stack([
        sg("h"), sg("T"), sg("k"), sg("n"),
        sg("gumbel"), sg("gumbel_sd"), sg("zf"), sg("emp"), sg("emp_sd"),
        sg("floor"), sg("err"),
        sg("corr"), sg("corr_lo"), sg("corr_hi"), sg("wrong"), sg("corr_se"),
        sg("old_corr"), sg("old_corr_lo"),
    ])
    np.savez(R("surface_analysis_recut.npz"),
             cols=np.array(keys), table=tab,
             tiers=S["tiers"], struct=S["struct"],
             verdicts=S["verdict"], verdicts_prefix=S["old_verdict"],
             estimator=S["estimator"], touched=S["touched"],
             n_onset=int(S["n_onset"]), n_onset_prefix=int(S["n_onset_banked"]),
             n_marginal=int((S["verdict"] == "MARGINAL").sum()),
             n_lost=int(S["n_lost"]), n_gained=int(S["n_gained"]),
             hstar_edge_cols=int(S["hstar_edge_cols_new"]),
             hstar_edge_cols_prefix=int(S["hstar_edge_cols_old"]),
             geo_ceiling=float(A["geo_ceiling"]),
             gateA_maxdev=float(S["gateA_maxdev"]), gateB_maxdev=float(S["gateB_maxdev"]),
             note=("counts RE-CUT from the per-realisation sf_sig banks against the ADOPTED "
                   "floor. ONSET iff corr at (floor_adopted + floor_adopted_sd) > 1. "
                   "*_prefix columns are the pre-fix banked values. Schema mirrors "
                   "surface_analysis.npz."))
    print(f"  surface_analysis_recut.npz  {tab.shape[0]} cells, N_onset {int(S['n_onset'])}")

    # ---------------- ch_floors_recut.npz (mirrors ch_floors.npz) --------------------------
    F = np.load(R("ch_floors.npz"), allow_pickle=True)
    cells = [str(x) for x in C["cells"]]
    out = {}
    for i in range(len(cells)):
        out[f"offN_{i}"] = F[f"offN_{i}"]                  # the raw offenders, carried through
        out[f"fN_{i}"] = cg("gumbel")[i]
        out[f"fN_sd_{i}"] = cg("gumbel_sd")[i]
        out[f"fN_emp_{i}"] = cg("emp")[i]
        out[f"fN_emp_sd_{i}"] = cg("emp_sd")[i]            # BOOTSTRAP sd (new)
        out[f"floor_adopted_{i}"] = cg("floor")[i]
        out[f"floor_adopted_sd_{i}"] = cg("err")[i]
        out[f"zero_frac_{i}"] = cg("zf")[i]                # REQUIRED column
        out[f"fN_n_{i}"] = int(F[f"fN_n_{i}"])
        out[f"fALL_{i}"] = float(F[f"fALL_{i}"])
        out[f"fALL_sd_{i}"] = float(F[f"fALL_sd_{i}"])
        out[f"prov_{i}"] = str(F[f"prov_{i}"])
    np.savez(R("ch_floors_recut.npz"),
             alpha=float(C["alpha"]), zf_gate=float(C["zf_gate"]),
             boot=int(C["boot"]), seed=int(C["seed"]),
             cells=C["cells"], estimator=C["estimator"],
             n_invalid=int(C["n_invalid"]), n_floor_up=int(C["n_floor_up"]),
             note=("ALL 26 cells have zero_frac > 0.20, so floor_adopted = fN_emp (empirical "
                   "q95) with a BOOTSTRAP sd for every cell. fN/fN_sd (the Gumbel) are kept "
                   "for comparison ONLY and are INVALID as floors here. Schema mirrors "
                   "ch_floors.npz."),
             **out)
    print(f"  ch_floors_recut.npz         {len(cells)} cells "
          f"({int(C['n_invalid'])} Gumbel-invalid)")

    # ---------------- ch_analysis_recut.npz (mirrors ch_analysis.npz surface_* block) -------
    AN = np.load(R("ch_analysis.npz"), allow_pickle=True)
    carry = {k: AN[k] for k in AN.files if not k.startswith("surface_")}
    np.savez(R("ch_analysis_recut.npz"),
             surface_tags=C["tags"],
             surface_corr=cg("corr"),
             surface_corr_lo=cg("corr_lo"),
             surface_corr_hi=cg("corr_hi"),
             surface_wrong=cg("wrong"),
             surface_corr_prefix=cg("old_corr"),
             surface_wrong_prefix=cg("old_wrong"),
             surface_fN=cg("gumbel"), surface_fN_sd=cg("gumbel_sd"),
             surface_fN_emp=cg("emp"), surface_fN_emp_sd=cg("emp_sd"),
             surface_floor_adopted=cg("floor"), surface_floor_adopted_sd=cg("err"),
             surface_zero_frac=cg("zf"),
             surface_estimator=C["estimator"], surface_grade=C["grade"],
             surface_nreal=cg("n"),
             gateA_maxdev=float(C["gateA_maxdev"]), gateB_maxdev=float(C["gateB_maxdev"]),
             note=("counts RE-CUT from the per-realisation ch_sig banks against the ADOPTED "
                   "floor. grade CONFIRMED iff count at (floor + bootstrap err) > 1; MARGINAL "
                   "iff it only clears the >1 bar at the floor. *_prefix = pre-fix banked. "
                   "The clock-sharing / pair / loop blocks are carried through UNCHANGED from "
                   "ch_analysis.npz — they do not depend on the floor."),
             **carry)
    print(f"  ch_analysis_recut.npz       {len(cells)} cells "
          f"(clock-sharing blocks carried through unchanged)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
