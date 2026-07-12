"""Track B certification criterion (criterion-v1) -- fit and apply.

THE CRITERION (three-layer, pre-registered 2026-07-12):

    DETECTION     dlnL_a > max(ln K_counted,a , DLNL_FLOOR)
    CERTIFICATION q_max,a > 0.9   (strict: > 0.99)   *within detections*

Layer 1 (ln K) is the per-pulsar trials factor: the fringe-breaking evidence must
beat the number of candidate fringes. Layer 2 (DLNL_FLOOR) is an ABSOLUTE floor,
fitted here against the banked null realisations; it exists because a pulsar with a
tight EM prior has so low a trials bar (J0437-4715: ln K = 1.39 nat, the array
minimum) that a pure-noise likelihood fluctuation clears it. Layer 3 is the
Bayesian bar, applied only to what already passed detection.

Provenance: FORGE_b1_loop.md sec.9 ran the two-layer gate (layers 1+3) and found a
residual null floor of 0.2-0.4 certs/realisation, all from the smallest-K pulsars;
it flagged the absolute floor as the fix and explicitly did NOT impose it ("it would
be a new criterion, not a re-score"). This script imposes it, fitted to the nulls.

Banks (reports/, from PORTER-2; no new realisations are ever drawn here):
  flat_*.npz (89)      per-realisation scores: b10 (60 signal), null/nullL/nullA/nullN
                       (27 null), wrongpos (2, tol_scale=5)
  geo_dlnl_bank.npz    GEO 40 zero-noise sky draws x 116 pulsars
  b1_dlnl_bank_{A,B}.npz  FORGE B1.0, 30 noisy realisations x 116 pulsars per arm

BINDING INVARIANTS (from the bank-builder; asserted below, do not weaken):
  1. Certification is defined on q_max, NOT on P_true. cert90 == (qmax > 0.9)
     bit-exact. The cells where qmax > 0.9 but P_true <= 0.9 are exactly the
     wrong-certifications (2 in Arm A, 8 in Arm B). Within detections qmax == P_true.
  2. GEO carries NO wrong-cert field by construction -- its Bayesian criterion is
     defined on the true fringe (P_true), so "wrong-cert" is not a notion that
     exists there. Do not synthesise one. GEO's canonical distance column is `lit`;
     the feather/script raw dlnL columns are carried in the bank but are not primary.

Usage:  python trackB_criterion.py [--reports DIR]
"""

import argparse
import glob
import os

import numpy as np

# The fitted floor. See fit_floor(): smallest floor zeroing ALL null certifications,
# = the largest dlnL among null cells that pass layers 1+3 (9.0094 nat, a pure-noise
# nullN J1713+0747 fluctuation), rounded UP so the zeroing property is preserved.
DLNL_FLOOR = 9.01

NULL_KINDS = ("null", "nullL", "nullA", "nullN")


def load_flats(reports):
    out = []
    for f in sorted(glob.glob(os.path.join(reports, "flat_*.npz"))):
        d = np.load(f, allow_pickle=True)
        out.append(
            dict(
                path=f,
                kind=str(d["kind"]),
                arm=str(d["arm"]),
                tol=float(d["tol_scale"]),
                dlnL=d["dlnL_det"],
                lnK=d["lnK"],
                qmax=d["qmax"],
                on_true=d["on_true"],
                names=d["names"],
            )
        )
    return out


def detect(dlnL, lnK, floor):
    """Layer 1 + layer 2. Strict `>`: a floor set AT an observed dlnL kills that cell."""
    return dlnL > np.maximum(lnK, floor)


def score(dlnL, lnK, qmax, floor, bar=0.9):
    det = detect(dlnL, lnK, floor)
    return det, det & (qmax > bar)


def check_invariants(reports):
    """Assert the two binding invariants. Raises AssertionError if a bank drifts."""
    for arm in "AB":
        b = np.load(os.path.join(reports, f"b1_dlnl_bank_{arm}.npz"), allow_pickle=True)
        assert np.array_equal(b["cert90"], b["qmax"] > 0.9), f"arm {arm}: cert90 != (qmax>0.9)"
        assert np.array_equal(b["cert99"], b["qmax"] > 0.99), f"arm {arm}: cert99 != (qmax>0.99)"
        # P_true does NOT reproduce the certifications -- that gap IS the wrong-certs.
        n_wrong = int(b["wrong90"].sum())
        assert not np.array_equal(b["cert90"], b["P_true"] > 0.9), f"arm {arm}: ptrue reproduced cert90"
        assert np.allclose(b["qmax"][b["detect"]], b["P_true"][b["detect"]]), "qmax != ptrue within detections"
        expect = {"A": 2, "B": 8}[arm]
        assert n_wrong == expect, f"arm {arm}: wrong90 {n_wrong} != {expect}"
    g = np.load(os.path.join(reports, "geo_dlnl_bank.npz"), allow_pickle=True)
    assert "wrong90" not in g.files, "GEO gained a wrong-cert field -- criterion assumption broken"
    assert str(g["variant"]) == "lit", "GEO canonical column is not `lit`"
    return True


def fit_floor(flats):
    """Smallest absolute floor that zeroes ALL null certifications.

    A null cell certifies under layers 1+3 iff dlnL > lnK and qmax > 0.9. Because
    detection uses a strict `>`, setting floor = max(dlnL over those cells) kills
    every one of them, and nothing smaller does. Returns (floor_min, offending cells).
    """
    nulls = [r for r in flats if r["kind"] in NULL_KINDS]
    off = []
    for r in nulls:
        m = (r["dlnL"] > r["lnK"]) & (r["qmax"] > 0.9)
        for i in np.where(m)[0]:
            off.append(
                dict(
                    kind=r["kind"],
                    name=str(r["names"][i]),
                    dlnL=float(r["dlnL"][i]),
                    lnK=float(r["lnK"][i]),
                    qmax=float(r["qmax"][i]),
                    on_true=bool(r["on_true"][i]),
                )
            )
    off.sort(key=lambda o: -o["dlnL"])
    return (max(o["dlnL"] for o in off) if off else 0.0), off, nulls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reports", default=os.path.join(os.path.dirname(__file__), "..", "reports"))
    a = ap.parse_args()
    R = os.path.abspath(a.reports)

    check_invariants(R)
    print("binding invariants: PASS (cert on q_max; GEO has no wrong-cert field)\n")

    flats = load_flats(R)
    floor_min, off, nulls = fit_floor(flats)

    print(f"=== K1 FLOOR FIT -- {len(nulls)} null realisations, {len(off)} null certs at floor=0 ===")
    print(f"{'kind':7s} {'pulsar':12s} {'dlnL':>7s} {'lnK':>6s} {'qmax':>6s}")
    for o in off:
        print(f"{o['kind']:7s} {o['name']:12s} {o['dlnL']:7.3f} {o['lnK']:6.3f} {o['qmax']:6.3f}")
    print(f"\nfloor_min = {floor_min:.4f} nat   ->  adopted DLNL_FLOOR = {DLNL_FLOOR} nat")
    assert DLNL_FLOOR >= floor_min, "adopted floor does not zero the nulls"

    g = np.load(os.path.join(R, "geo_dlnl_bank.npz"), allow_pickle=True)
    BA = np.load(os.path.join(R, "b1_dlnl_bank_A.npz"), allow_pickle=True)
    BB = np.load(os.path.join(R, "b1_dlnl_bank_B.npz"), allow_pickle=True)

    print(f"\n=== K2 FINAL TABLE (floor = {DLNL_FLOOR} nat) ===")
    print(f"{'population':24s} {'N':>3s} {'Bayes':>6s} {'2-layer':>7s} {'FINAL':>6s} {'strict':>6s} {'wrong':>5s}")

    # GEO: certification bar is on P_true (its Bayesian criterion is defined on the
    # true fringe). No wrong-cert notion exists here -- invariant 2.
    _, c90 = score(g["dlnL"], g["lnK"], g["P_true"], DLNL_FLOOR)
    _, c99 = score(g["dlnL"], g["lnK"], g["P_true"], DLNL_FLOOR, bar=0.99)
    _, two = score(g["dlnL"], g["lnK"], g["P_true"], 0.0)
    print(f"{'GEO zero-noise / draw':24s} {40:3d} {g['bayes'].sum()/40:6.2f} {two.sum()/40:7.2f} "
          f"{c90.sum()/40:6.3f} {c99.sum()/40:6.3f} {'n/a':>5s}")

    for tag, b in (("FORGE B1.0 Arm A", BA), ("FORGE B1.0 Arm B", BB)):
        _, c90 = score(b["dlnL"], b["lnK"], b["qmax"], DLNL_FLOOR)
        _, c99 = score(b["dlnL"], b["lnK"], b["qmax"], DLNL_FLOOR, bar=0.99)
        _, two = score(b["dlnL"], b["lnK"], b["qmax"], 0.0)
        wrong = int((c90 & ~b["on_true"]).sum())
        print(f"{tag+' / real':24s} {30:3d} {b['cert90'].sum()/30:6.2f} {two.sum()/30:7.2f} "
              f"{c90.sum()/30:6.3f} {c99.sum()/30:6.3f} {wrong:5d}")

    print()
    total = 0
    for kind in NULL_KINDS:
        rs = [r for r in flats if r["kind"] == kind]
        bay = sum(int((r["qmax"] > 0.9).sum()) for r in rs)
        two = sum(int(((r["dlnL"] > r["lnK"]) & (r["qmax"] > 0.9)).sum()) for r in rs)
        fin = sum(int(score(r["dlnL"], r["lnK"], r["qmax"], DLNL_FLOOR)[1].sum()) for r in rs)
        total += fin
        print(f"{'null: '+kind:24s} {len(rs):3d} {bay/len(rs):6.2f} {two/len(rs):7.2f} "
              f"{fin/len(rs):6.3f} {0:6.3f} {0:5d}")
    print(f"{'ALL NULLS':24s} {sum(1 for r in flats if r['kind'] in NULL_KINDS):3d} "
          f"{'':6s} {'':7s} {total:6d}   <== must be 0")
    assert total == 0, "floor does not zero the nulls"


if __name__ == "__main__":
    main()
