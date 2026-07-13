"""Re-bank anchor_ladder.npz's `offenders` in METADATA-ROW order.

WHAT FLOOR_FIX_provisional §5 SAID, AND WHAT IS ACTUALLY TRUE
------------------------------------------------------------
§5 reported that `offenders` is stored (rung-major, cell-minor) while every sibling column
is (cell-major, rung-minor), and concluded:

    "The array is the TRANSPOSE of its own metadata, and NOTHING IN THE FILE SAYS SO."
    "FIX ON ACCRE: re-bank `offenders` in metadata-row order, OR ADD AN EXPLICIT INDEX KEY."

The first half is right; the second half is not. **The file already carries an explicit index
key.** `anchor_ladder.npz` banks `offender_index`, a 48-element string column whose entry j
is the label of `offenders[j]`:

    offender_index[0] = 'R0|-13.25|lit'      offenders[0] belongs to that (rung, h, tier)
    offender_index[1] = 'R0|-13.25|vlbi'     ... i.e. rung-major, cell-minor. Declared.

Re-cutting by that key reproduces the banked `emp_q95` AND `zero_frac` to exactly 0.000e+00
(verified below), and the permutation §5 reverse-engineered, perm[j] = (j%6)*8 + (j//6), is
*identically* the mapping `offender_index` already encodes. So:

  * ANCHOR's published §3 ladder table is CORRECT.                       (§5 got this right)
  * The bank was SELF-DESCRIBING, not undeclared.                        (§5 got this wrong)
  * The trap is real but narrower: a re-cut that ignores `offender_index` and assumes
    offenders[j] pairs with metadata row j marries the wrong cell's label, silently, and is
    off by up to 79.9 nat.

THE FIX APPLIED HERE: fix the ARRAY, not the metadata.
`offenders` is permuted into metadata-row order, so offenders[j] now pairs with rung[j],
h[j], tier[j], emp_q95[j] — the naive re-cut becomes the CORRECT one. `offender_index` is
rewritten in the same order (so it still labels each row, and now agrees with its siblings
rather than cross-cutting them), and `offender_index_orientation` states the convention in
words. Nothing else in the file is touched: every published column keeps its banked value.

Why fix the array rather than relabel the metadata: the metadata order (cell-major) is the
order of the PUBLISHED §3 table and of every downstream consumer. Permuting the one column
that disagrees is the change that leaves every existing reader correct.

Run: python CW_transition/rebank_anchor_ladder.py       (CPU, seconds)
"""
import os, sys, shutil
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORTS = os.path.join(ROOT, "reports")
SRC = os.path.join(REPORTS, "anchor_ladder.npz")
ALPHA = 0.05
NR, NC = 6, 8              # 6 rungs (R0,R1,R2,R3,R2o,R3o) x 8 cells (4 h x 2 tier)


def label(rung, h, tier):
    return f"{rung}|{h}|{tier}"


def key(rung, h, tier):
    return (str(rung), round(float(h), 4), str(tier))


def main():
    z = np.load(SRC, allow_pickle=True)
    offs = z["offenders"]
    oi = [str(x) for x in z["offender_index"]]
    rung = [str(x) for x in z["rung"]]
    h = z["h"]
    tier = [str(x) for x in z["tier"]]
    empb, zfb = z["emp_q95"], z["zero_frac"]
    n = len(rung)

    print("=" * 98)
    print("RE-BANK anchor_ladder.npz — `offenders` into metadata-row order")
    print("=" * 98)
    print(f"  {n} (rung, cell) rows, {offs.shape[1]} offenders each")

    # ---- 1. the declared index is authoritative: prove it ----
    pos = {}
    for i, lab in enumerate(oi):
        r, hv, t = lab.split("|")
        pos[key(r, hv, t)] = i
    src = np.array([pos[key(rung[j], h[j], tier[j])] for j in range(n)])

    rec = np.array([np.quantile(offs[src[j]], 1 - ALPHA) for j in range(n)])
    zrec = np.array([float((offs[src[j]] == 0).mean()) for j in range(n)])
    naive = np.array([np.quantile(offs[j], 1 - ALPHA) for j in range(n)])
    PERM = np.array([(j % NR) * NC + (j // NR) for j in range(n)])

    print("\n  A. IS THE DECLARED INDEX AUTHORITATIVE?")
    print(f"     re-cut by offender_index   : max|dev| emp_q95 = {np.abs(rec-empb).max():.3e}   "
          f"zero_frac = {np.abs(zrec-zfb).max():.3e}")
    print(f"     naive offenders[j] <-> [j] : max|dev| emp_q95 = {np.abs(naive-empb).max():8.3f}"
          f"   <- the 79.9-nat trap")
    print(f"     FLOOR_FIX §5's PERM equals the index mapping? "
          f"{bool(np.array_equal(PERM, src))}")
    if np.abs(rec - empb).max() != 0.0 or np.abs(zrec - zfb).max() != 0.0:
        print("     *** STOP: the declared index does NOT reproduce the banked ladder.")
        return 1
    print("     => YES. The bank was self-describing. §5's 'nothing in the file says so' is "
          "incorrect;\n        the trap is only sprung by a re-cut that IGNORES offender_index.")

    # ---- 2. permute the array into metadata-row order ----
    # The index strings are carried over verbatim (just reordered), so the labels keep the
    # exact spelling the bank already used -- no reformatting of h can creep in.
    new_offs = offs[src]
    new_index = np.array([oi[src[j]] for j in range(n)])

    # ---- 3. GATE: the naive re-cut must now be the correct one ----
    rec2 = np.array([np.quantile(new_offs[j], 1 - ALPHA) for j in range(n)])
    zrec2 = np.array([float((new_offs[j] == 0).mean()) for j in range(n)])
    d_q, d_z = np.abs(rec2 - empb).max(), np.abs(zrec2 - zfb).max()
    print("\n  B. AFTER THE PERMUTATION — the NAIVE re-cut on the re-banked file:")
    print(f"     max|dev| vs banked emp_q95   = {d_q:.3e}")
    print(f"     max|dev| vs banked zero_frac = {d_z:.3e}")
    if d_q != 0.0 or d_z != 0.0:
        print("     *** STOP: permuted array does not reproduce the banked ladder.")
        return 1
    print("     => the trap is closed: offenders[j] now pairs with metadata row j, exactly.")

    # the reordered index must now label its own row
    consistent = all(new_index[j].split("|")[0] == rung[j]
                     and new_index[j].split("|")[2] == tier[j]
                     and abs(float(new_index[j].split("|")[1]) - float(h[j])) < 1e-9
                     for j in range(n))
    print(f"     offender_index[j] now agrees with (rung[j], h[j], tier[j]): {consistent}")
    if not consistent:
        print("     *** STOP: re-ordered index does not match its sibling metadata.")
        return 1

    # ---- 4. verify the PUBLISHED §3 table reproduces from the re-banked file ----
    print("\n  C. ANCHOR's PUBLISHED §3 TABLE, reproduced from the RE-BANKED offenders")
    print(f"     (R0 control rows; empirical q95 is the estimator ANCHOR issues the verdict on)")
    print(f"     {'rung':>5} {'cell':>16} | {'zero%':>6} | {'emp q95 (re-cut)':>17} "
          f"{'banked':>9} {'dev':>9}")
    shown = 0
    for j in range(n):
        if rung[j] != "R0":
            continue
        print(f"     {rung[j]:>5} {f'{h[j]:.2f} {tier[j]}':>16} | {100*zrec2[j]:6.1f} | "
              f"{rec2[j]:17.3f} {empb[j]:9.3f} {abs(rec2[j]-empb[j]):9.1e}")
        shown += 1
    print(f"     all {n} rows (not just the {shown} shown) reproduce to "
          f"{max(d_q, d_z):.1e}")

    # ---- 5. write ----
    bak = SRC + ".preRECUT.bak"
    if not os.path.exists(bak):
        shutil.copy2(SRC, bak)
        print(f"\n  original preserved -> {os.path.basename(bak)}")

    payload = {k: z[k] for k in z.files}
    payload["offenders"] = new_offs
    payload["offender_index"] = new_index
    payload["offender_index_orientation"] = (
        "row j of `offenders` belongs to metadata row j: (rung[j], h[j], tier[j]), and "
        "offender_index[j] is its label. Cell-major, rung-minor -- the SAME order as every "
        "other column in this file and as the published §3 table. RE-BANKED by "
        "CW_transition/rebank_anchor_ladder.py; the original stored `offenders` rung-major "
        "(declared then only via offender_index), which made a naive offenders[j]<->row-j "
        "re-cut wrong by up to 79.9 nat.")
    payload["rebanked_by"] = "RECUT 2026-07-13 (CW_transition/rebank_anchor_ladder.py)"
    np.savez(SRC, **payload)
    print(f"  RE-BANKED -> reports/anchor_ladder.npz")
    print("  (only `offenders` was permuted and `offender_index` reordered to match; every "
          "published\n   column keeps its banked value, and the §3 table above proves it.)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
