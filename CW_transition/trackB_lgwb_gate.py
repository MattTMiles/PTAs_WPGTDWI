"""GATE — the GWB square root is thread-count invariant (the BLAS hazard, CHORUS §0.1).

WHAT IS BEING GATED, and why it needed a gate at all
----------------------------------------------------
`NoiseDrawer` used to build the GWB noise square root by `np.linalg.eigh(Phi_gwb)`. Phi_gwb
is HD-correlated and NEAR-DEGENERATE, so LAPACK's eigenvector basis inside a degenerate
subspace depends on the BLAS THREAD COUNT. Same seed + different `--cpus-per-task` = a
different (rotated, equal-in-distribution) GWB realisation. CHORUS's g1 gate caught this as
an O(0.1-1) nat dlnL shift on essentially every pulsar. Until the fix, cpus-per-task was an
undeclared input to every banked noisy number in the repo.

The fix banks `L_gwb` to disk, so no BLAS call sits between the seed and the draw.

THE GATE (run: python trackB_lgwb_gate.py)
------------------------------------------
G1  HAZARD IS REAL      — the eigh path returns a DIFFERENT L_gwb at different thread counts
                          on a near-degenerate Phi. (If this ever stops failing, the hazard
                          is gone and the bank is belt-and-braces. It is not gone today.)
G2  BANKED PATH IS SAFE — with a banked L_gwb, one `NoiseDrawer.draw(seed)` is BIT-IDENTICAL
                          across BLAS thread counts {1, 2, 8, 16, 24}. This is the property
                          every banked realisation needs and never had.
G3  CORRUPTION IS CAUGHT— a tampered bank raises rather than silently drawing wrong noise.

This is a UNIT test: it drives the real `load_or_build_L_gwb` / `NoiseDrawer.draw` code paths
against a small synthetic near-degenerate Phi and a stubbed Split, so it runs on CPU in
seconds and needs no GPU, no 116-pulsar build, and no banked npz.

SCOPE OF INFERENCE. This gate proves the DRAW IS DETERMINISTIC given a bank. It does NOT
prove the bank is the right one. The canonical `b1_L_gwb.npz` must be generated on ACCRE at
the banked convention `--cpus-per-task=8` and validated against ANCHOR's g1 replay (80
banked `ig_nullN_*` realisations, bit-identical) before any banked result is re-derived
through it. cronus (24 cores) cannot produce the canonical basis.
"""
import os, sys, json, subprocess, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
THREAD_COUNTS = [1, 2, 8, 16, 24]
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")

# ---------------------------------------------------------------- synthetic problem
NPSR, NGP = 24, 14          # 336x336 — same STRUCTURE as the real 5336x5336, small enough to unit-test


def make_phi(npsr=NPSR, ngp=NGP):
    """A near-degenerate HD-correlated GWB prior covariance, in the real one's image."""
    rng = np.random.default_rng(20260713)
    pos = rng.normal(size=(npsr, 3))
    pos /= np.linalg.norm(pos, axis=1, keepdims=True)
    c = np.clip(pos @ pos.T, -1.0, 1.0)
    x = 0.5 * (1.0 - c)
    with np.errstate(divide="ignore", invalid="ignore"):        # Hellings-Downs
        hd = 1.5 * x * np.log(x) - 0.25 * x + 0.5
    hd[np.isnan(hd)] = 1.0
    np.fill_diagonal(hd, 1.0)
    hd = 0.5 * (hd + hd.T)
    # a power-law GP spectrum, repeated across pulsars -> the DEGENERACY that bites
    f = np.arange(1, ngp + 1, dtype=float)
    rho = f ** (-13.0 / 3.0)
    Phi = np.kron(hd, np.diag(rho))
    return 0.5 * (Phi + Phi.T)


class StubSplit:
    """The only surface NoiseDrawer touches."""
    def __init__(self, npsr=NPSR, ngp=NGP, ntoa=7, nrn=5):
        rng = np.random.default_rng(11)
        self.npsr, self.ngp_gwb = npsr, ngp
        self.N_diag = [np.full(ntoa, 1e-13) for _ in range(npsr)]
        self.Phi_rn_diag = np.full((npsr, nrn), 1e-15)
        self.Fs_rn = [rng.normal(size=(ntoa, nrn)) for _ in range(npsr)]
        self.Ts_gwb = [rng.normal(size=(ntoa, ngp)) for _ in range(npsr)]


def child():
    """Run inside a subprocess at a pinned thread count; print a JSON fingerprint."""
    sys.path.insert(0, HERE)
    import trackB_b1_core as C
    Phi = make_phi()
    bank = os.environ.get("GATE_BANK", "")
    sp = StubSplit()
    amo = {"internals": {"Pinv_gwb": np.linalg.inv(Phi)}}
    if bank:
        nd = C.NoiseDrawer(sp, amo, lgwb_path=bank)
    else:                      # force the legacy eigh path: point at a file that cannot exist
        nd = C.NoiseDrawer(sp, amo, lgwb_path=os.path.join(tempfile.gettempdir(), "__no_such_lgwb__.npz"))
    r = np.concatenate([np.asarray(x) for x in nd.draw(4242)])
    print("@@" + json.dumps({
        "L_fp": C.lgwb_fingerprint(nd.L_gwb),
        "draw_fp": C.lgwb_fingerprint(r),
        "prov": nd.lgwb_prov,
    }))


def run_at(nthreads, bank=""):
    env = dict(os.environ)
    for v in THREAD_VARS:
        env[v] = str(nthreads)
    env["GATE_BANK"] = bank
    env["JAX_PLATFORMS"] = "cpu"
    env["GATE_CHILD"] = "1"
    out = subprocess.run([sys.executable, os.path.abspath(__file__)],
                         capture_output=True, text=True, env=env)
    line = [l for l in out.stdout.splitlines() if l.startswith("@@")]
    if not line:
        raise RuntimeError(f"child failed at {nthreads} threads:\n{out.stdout}\n{out.stderr}")
    return json.loads(line[-1][2:])


def main():
    print(__doc__.split("THE GATE")[0].strip()[:0] or "", end="")
    print("=" * 78)
    print("GATE — L_gwb thread-count invariance   (BLAS hazard, CHORUS 2026-07-13 §0.1)")
    print("=" * 78)
    Phi = make_phi()
    ev = np.linalg.eigvalsh(Phi)
    print(f"synthetic Phi: {Phi.shape[0]}x{Phi.shape[0]}, cond = {ev.max()/max(ev.min(),1e-300):.2e}, "
          f"min eig = {ev.min():.2e}   (near-degenerate, as the real one is)")
    print()

    # ---- G1: the hazard is real (legacy eigh path) ----
    print("G1  HAZARD — legacy eigh path, no bank:")
    legacy = {}
    for n in THREAD_COUNTS:
        r = run_at(n, bank="")
        legacy[n] = r
        print(f"      threads={n:<3d}  L_gwb fp = {r['L_fp']}   draw fp = {r['draw_fp']}")
    uniq_L = {r["L_fp"] for r in legacy.values()}
    uniq_d = {r["draw_fp"] for r in legacy.values()}
    g1 = len(uniq_L) > 1 or len(uniq_d) > 1
    print(f"      -> {len(uniq_L)} distinct L_gwb, {len(uniq_d)} distinct draws across thread counts")
    print(f"      -> G1 {'CONFIRMED: the draw MOVES with the thread count' if g1 else 'not reproduced on this BLAS today'}")
    print()

    # ---- bank it ----
    bank = os.path.join(tempfile.gettempdir(), "gate_L_gwb.npz")
    w, V = np.linalg.eigh(Phi)
    L = V * np.sqrt(np.clip(w, 0.0, None))
    sys.path.insert(0, HERE)
    import trackB_b1_core as C
    np.savez(bank, L_gwb=L, fingerprint=C.lgwb_fingerprint(L), cpus=os.cpu_count())
    print(f"    banked a test L_gwb -> {bank}  fp={C.lgwb_fingerprint(L)}")
    print()

    # ---- G2: the banked path is invariant ----
    print("G2  FIX — banked L_gwb, BLAS out of the draw path:")
    banked = {}
    for n in THREAD_COUNTS:
        r = run_at(n, bank=bank)
        banked[n] = r
        print(f"      threads={n:<3d}  L_gwb fp = {r['L_fp']}   draw fp = {r['draw_fp']}")
    uL = {r["L_fp"] for r in banked.values()}
    uD = {r["draw_fp"] for r in banked.values()}
    g2 = len(uL) == 1 and len(uD) == 1
    print(f"      -> {len(uL)} distinct L_gwb, {len(uD)} distinct draws")
    print(f"      -> G2 {'PASS — one draw, BIT-IDENTICAL at 1/2/8/16/24 threads' if g2 else 'FAIL'}")
    print()

    # ---- G3: corruption is caught ----
    print("G3  CORRUPTION:")
    bad = os.path.join(tempfile.gettempdir(), "gate_L_gwb_bad.npz")
    Lb = L.copy(); Lb[0, 0] += 1e-9
    np.savez(bad, L_gwb=Lb, fingerprint=C.lgwb_fingerprint(L), cpus=8)   # hash of the ORIGINAL
    try:
        C.load_or_build_L_gwb(Phi, bad)
        g3 = False; print("      -> G3 FAIL: tampered bank loaded silently")
    except RuntimeError as e:
        g3 = True; print(f"      -> G3 PASS: {str(e)[:70]}...")
    print()

    ok = g2 and g3
    print("=" * 78)
    print(f"G1 hazard-confirmed={g1}   G2 invariance={'PASS' if g2 else 'FAIL'}   "
          f"G3 corruption={'PASS' if g3 else 'FAIL'}")
    print("VERDICT:", "PASS — the draw no longer depends on the thread count." if ok else "FAIL")
    print("=" * 78)
    print("\nSCOPE: proves the draw is deterministic GIVEN a bank. Does NOT prove the bank is\n"
          "the one the repo's banked realisations were drawn under. Canonical b1_L_gwb.npz must\n"
          "be generated on ACCRE at --cpus-per-task=8 and validated against ANCHOR's g1 replay\n"
          "(80 banked ig_nullN_* realisations, bit-identical) before re-deriving banked results.")
    return 0 if ok else 1


if __name__ == "__main__":
    if os.environ.get("GATE_CHILD"):
        child()
    else:
        sys.exit(main())
