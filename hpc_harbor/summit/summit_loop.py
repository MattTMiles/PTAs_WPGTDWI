#!/usr/bin/env python
"""SUMMIT D2 — the prior-tier single-dial ladder at the feasible rung (r13p9).

WRAPPER, NOT A FORK: the GLACIER driver (glacier_loop.py) is imported STOCK and
runtime-patched (the chorus-C7 / COMPASS pattern). glacier_loop.py itself is NOT
edited -- the pending array-ladder tasks (12763564) import it at start time, and
editing it mid-array is the g1 provenance hazard.

THE DIAL (declared):
  lit    -- the banked ladder cells (INHERITED, never re-run here)
  vlbi   -- sigma_d -> min(lit, 0.001 kpc = 1 pc) on the UNION18   [ignite.py IG4,
            verbatim: improve, never degrade; truth draw follows the tier]
  subpc  -- sigma_d -> min(lit, 0.0001 kpc = 0.1 pc) on the UNION18 (same rule, 10x)

BANKS: $GLACIER_OUT must point at SUMMIT_results (set via APPTAINERENV_GLACIER_OUT --
plain exports do not survive harbor_py --cleanenv, the 12740157 lesson). Two
provenances never share one directory. Stems inherit glacier_loop's cell2 form and
carry the tier: gl2_r13p9_cell_none_s{sky}_T30_{tier}.
Evidence gate: HOLDS_CLEARED + gl2_ladder_gates are COPIED into SUMMIT_results with a
copy-record line -- the holds are the same green record on the same tree.

Modes:
  gate  -- SG-D2a/b/c (b/c CPU; a = one-iteration lit reproduction vs the banked cell)
  cell2 / null2 -- the D2 ladder cell (tier from --tier), signal / scrambled arm
"""
import os, sys, glob

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "glacier"))
sys.path.insert(0, os.path.join(HERE, "..", "ignite"))
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")

SUMMIT_OUT_REQUIRED = "SUMMIT_results"


def _check_out():
    out = os.environ.get("GLACIER_OUT", "")
    if SUMMIT_OUT_REQUIRED not in out:
        print(f"REFUSED: GLACIER_OUT={out!r} does not point at {SUMMIT_OUT_REQUIRED}. "
              "SUMMIT banks never land in GLACIER_results (two provenances, one dir = "
              "the g1 hazard). Set APPTAINERENV_GLACIER_OUT.", flush=True)
        sys.exit(3)
    os.makedirs(out, exist_ok=True)
    return out


VLBI_SIG_KPC = 0.001     # 1 pc      (ignite.py convention, verbatim)
SUBPC_SIG_KPC = 0.0001   # 0.1 pc    (same rule, the D2 extension rung -- DECLARED)


def _patch_priors():
    """Wrap A2.load_real_priors to append the vlbi/subpc columns (union-18 min rule)."""
    import stagec_anchor_a2 as A2
    from ignite import UNION18
    if getattr(A2.load_real_priors, "_summit_patched", False):
        return
    orig = A2.load_real_priors

    def load_real_priors_summit(names):
        pr = orig(names)
        for key, sig_cap in (("vlbi", VLBI_SIG_KPC), ("subpc", SUBPC_SIG_KPC)):
            col = pr["lit"].copy()
            for n in UNION18:
                if n in names:
                    i = list(names).index(n)
                    col[i] = min(col[i], sig_cap)
            pr[key] = col
        return pr

    load_real_priors_summit._summit_patched = True
    A2.load_real_priors = load_real_priors_summit


def _gates(GL):
    """SG-D2b/c (CPU): prior-column structure + tier truth-draw stream alignment."""
    import numpy as np
    import stagec_anchor_a2 as A2
    from ignite import UNION18
    # names via the pop loader is heavy; use the K-table bank's name column (the same
    # source load_real_priors itself keys on)
    a1 = np.load(A2.CWT + "/anchor_a1_Ktable.npz", allow_pickle=True)
    names = [str(n) for n in a1["fn"]]
    pr = A2.load_real_priors(names)
    assert set(("vlbi", "subpc")) <= set(pr), "patch did not land"
    lit, vl, sp = pr["lit"], pr["vlbi"], pr["subpc"]
    u = np.array([n in UNION18 for n in names])
    ok_b = (np.all(vl[~u] == lit[~u]) and np.all(sp[~u] == lit[~u])
            and np.all(vl[u] <= lit[u]) and np.all(sp[u] <= vl[u])
            and np.all(vl[u] <= VLBI_SIG_KPC + 1e-15)
            and np.all(sp[u] <= SUBPC_SIG_KPC + 1e-15))
    print(f"[SG-D2b] tier columns: off-union identical to lit; union monotone "
          f"subpc<=vlbi<=lit; caps respected -> {'PASS' if ok_b else 'FAIL'}")
    n_lit_below = int((lit[u] <= VLBI_SIG_KPC).sum())
    print(f"         (union-18 rows already at/below 1 pc under lit: {n_lit_below} "
          f"-- the vlbi dial moves {int(u.sum()) - n_lit_below} pulsars)")
    # SG-D2c: the tiered truth draw consumes an identical RNG stream per pulsar
    # (one choice + one uniform), so off-union pulsars must draw IDENTICAL n_true.
    import ignite as IG
    rng_probe = np.random.default_rng(0)
    npsr = len(names)
    L0 = 1.0 + rng_probe.uniform(0, 2, npsr)
    dL = 10 ** rng_probe.uniform(-4, -2, npsr)
    EV = np.array([L0[a] + dL[a] * np.arange(-40, 41) for a in range(npsr)])
    P = dict(npsr=npsr, L0=L0, priors=pr)
    out = {}
    for tier in ("lit", "vlbi", "subpc"):
        out[tier] = IG.draw_true_distances_tier(P, dL, EV, seed=12345, tier=tier)[1]
    same_off = (np.all(out["lit"][~u] == out["vlbi"][~u])
                and np.all(out["lit"][~u] == out["subpc"][~u]))
    tighter = (np.abs(out["subpc"][u]).mean() <= np.abs(out["vlbi"][u]).mean() + 1e-9
               <= np.abs(out["lit"][u]).mean() + 2.0)
    print(f"[SG-D2c] tier truth draw: off-union n_true identical across tiers "
          f"({'PASS' if same_off else 'FAIL'}); union |n_true| means "
          f"lit/vlbi/subpc = {np.abs(out['lit'][u]).mean():.2f}/"
          f"{np.abs(out['vlbi'][u]).mean():.2f}/{np.abs(out['subpc'][u]).mean():.2f} "
          f"({'PASS' if tighter else 'FAIL'})")
    if not (ok_b and same_off and tighter):
        sys.exit(2)
    print("=== SUMMIT D2 CPU GATES: PASS ===")


def _gatea_compare():
    """CPU: compare the d2gate lit i0 bank against the banked GLACIER ladder cell.
    Cross-lane tolerant (draws are host-free through the banked L_gwb, S4.8 item 2;
    lnL reductions differ at ~1e-12 cross-arch)."""
    import numpy as np
    out = os.environ["GLACIER_OUT"]
    a = glob.glob(f"{out}/d2gate/gl2_r13p9_cell_none_s0_T30_lit_i0__*.npz")
    b = glob.glob("/home/mattm/projects/HSYMT/GLACIER_results/"
                  "gl2_r13p9_cell_none_s0_T30_lit_i0__*.npz")
    if not (a and b):
        print(f"[SG-D2a] MISSING bank: gate={a} ref={b}"); sys.exit(2)
    za, zb = np.load(a[0], allow_pickle=True), np.load(b[0], allow_pickle=True)
    same_fed = np.array_equal(za["fed_mask"], zb["fed_mask"])
    same_nres = int(za["n_resolved"]) == int(zb["n_resolved"])
    d_abg = abs(float(za["a_bg"]) - float(zb["a_bg"]))
    ra, rb = za["conc_ratio"], zb["conc_ratio"]
    d_conc = float(np.nanmax(np.abs(ra - rb)))
    # conc_ratio is Fisher-derived (second differences at ~0.2-nat curvature): cross-arch
    # reduction noise amplifies to ~1e-2 there. The gate requirement is DECISION
    # INVARIANCE (no slot's ratio crosses the 0.5 feed bar between banks) + a 0.1 sanity
    # cap; measured cross-lane wobble 0.0166 vs a 0.233 smallest margin (trail: SUMMIT
    # report SS1). All non-Fisher columns stay exact.
    crossings = int(((ra < 0.5) != (rb < 0.5)).sum())
    d_floor = abs(float(za["floor"]) - float(zb["floor"]))
    same_zf = float(za["zero_fraction"]) == float(zb["zero_fraction"])
    ok = same_fed and same_nres and d_abg < 1e-3 and crossings == 0 and \
        d_conc < 0.1 and d_floor < 1e-6 and same_zf
    print(f"[SG-D2a] wrapper inertness at tier=lit vs banked "
          f"({str(zb['_lane'])} -> {str(za['_lane'])}): fed_mask "
          f"{'==' if same_fed else 'DIFF'}, n_res {'==' if same_nres else 'DIFF'}, "
          f"|d a_bg|={d_abg:.2e}, max|d conc|={d_conc:.2e} (bar-crossings {crossings}), "
          f"|d floor|={d_floor:.2e}, zf {'==' if same_zf else 'DIFF'} -> "
          f"{'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 2)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "gatea", "gateacmp", "cell2", "null2"])
    ap.add_argument("--tier", choices=["lit", "vlbi", "subpc"], default="vlbi")
    ap.add_argument("--rung", default="r13p9")
    ap.add_argument("--arm", default="none")
    ap.add_argument("--sky", type=int, default=0)
    ap.add_argument("--real", type=int, default=0)
    ap.add_argument("--iters", type=int, default=6)
    # PHASE-2 COMPOSITION (PHOENIX 2026-07-29): forward the ARRAY-side venue axes so a
    # composed cell (prior tier x array capability) can be cut. Defaults 1.0 / 30
    # reproduce every banked D2 stem bit-for-bit -- this is an argument plumb, not a
    # behaviour change. Both dials are separately validated: the tier dial by SG-D2a/b/c
    # (S1.1) and the wscale/T dial by the banked array ladder (12763564).
    ap.add_argument("--wscale", type=float, default=1.0, choices=[1.0, 0.5, 0.25])
    ap.add_argument("--t", type=int, default=30, choices=[30, 40])
    a = ap.parse_args()

    out = _check_out()
    _patch_priors()

    if a.mode == "gate":
        return _gates(None)
    if a.mode == "gateacmp":
        return _gatea_compare()
    if a.mode == "gatea":
        # SG-D2a GPU leg: the STOCK lit path through the wrapper, first iteration
        # only, own subdir (never mistakable for a campaign cell). No marker needed
        # -- this IS a gate (glacier_loop mode_gate precedent).
        import glacier_loop as GL
        GL._redirect_out("d2gate")
        print(f"[SUMMIT-D2] gatea: lit s0 i0 -> {out}/d2gate", flush=True)
        return GL.run_cell("none", 0, n_iter=1, rung="r13p9")

    # evidence gate, honored not bypassed: the green record travels with the banks
    marker = f"{out}/HOLDS_CLEARED"
    gates2 = glob.glob(f"{out}/gl2_ladder_gates__*.npz")
    if not (os.path.exists(marker) and gates2):
        print(f"REFUSED: {marker} + gl2_ladder_gates bank required in SUMMIT_results "
              "(copy from GLACIER_results with a copy-record line; same tree, same "
              "green jobs).", flush=True)
        return 3
    with open(marker) as fh:
        print(f"[SUMMIT-D2] holds cleared: {fh.read().strip()}", flush=True)

    import glacier_loop as GL
    GL.TIER = a.tier                       # call-time global: scoreboard prior column,
    #                                        tier truth draw, and the bank stem
    print(f"[SUMMIT-D2] tier={a.tier} rung={a.rung} arm={a.arm} sky={a.sky} "
          f"wscale={a.wscale:g} T={a.t} OUT={out}", flush=True)
    return GL.run_cell(a.arm, a.sky, scrambled=(a.mode == "null2"), real=a.real,
                       n_iter=a.iters, rung=a.rung, t_label=a.t, wscale=a.wscale)


if __name__ == "__main__":
    sys.exit(main() or 0)
