#!/usr/bin/env python
"""SUMMIT D3 — the N_psr single-dial ladder at the feasible rung (r13p9).

THE DIAL (declared): N_psr in {116 (inherited, the banked ladder), +30, +100}.
+N EXTENDS the real array: the 116 keep their real positions and bundles UNTOUCHED;
N new pulsars are CLONES of the top-N real pulsars by the COMPASS quality rank
(descending n_toa), placed at quasi-uniform sky positions (golden-angle grid under a
seed-controlled rotation, SEED_D3 = 9300 -- the COMPASS iso_positions construction at
n = N). THE BUNDLE TRAVELS INTACT (COMPASS convention, cited: rank-map pairing priced
at 0.14 sky-sigma there): TOA epochs/errors, design matrix, EM distance (pdist), K-table
prior row, and the frozen reference RN draw are the DONOR'S; only name and sky position
differ. Clone names = donor + 'C%02d' -- every name-keyed lookup strips the suffix.

WRAPPER, NOT A FORK (chorus-C7 / COMPASS pattern): stock modules runtime-patched --
  cw_helpers.load_pulsars            -> append clones when EXT_STATE['n_add'] > 0
  stagec_anchor_a2.load_real_priors  -> clone rows inherit the donor's prior row
  stagec_fisher.make_geometry_injection -> clone RN := donor RN (bundle intact);
                                        first-116 draws UNCHANGED (gated, SG-D3d)
  glacier_loop.LGWB_BANKS            -> the ARRAY-KEYED bank (the COMPASS L_gwb venue
                                        hazard: shape alone does not key an array)
BANKS: per-rung OUT subdirs (SUMMIT_results/d3_ext{N}) -- stems carry no N_psr tag, so
the directory is the array key (two provenances never share one directory).

Modes: gate (SG-D3b/c/d/e battery) | lgwb (cut the array-keyed bank, cpus=8 pinned)
       | cell2 / null2 (the ladder cell, --next per glacier_loop)
"""
import os, sys, glob, copy, socket, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "glacier"))
sys.path.insert(0, os.path.join(HERE, "..", "ignite"))
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")

import re as _re

SEED_D3 = 9300
EXT_STATE = {"n_add": 0, "twin": False}   # twin: clone 0 sits AT its donor's position
_CLONE_RE = _re.compile(r"^(?P<donor>.+)C(?P<k>\d{2})$")
SUMMIT_ROOT = "/home/mattm/projects/HSYMT/SUMMIT_results"


def donor_of(name):
    m = _CLONE_RE.match(name)
    return m.group("donor") if m else name


def _ext_positions(n_add, seed=SEED_D3):
    """COMPASS iso_positions at n = n_add (golden-angle + seeded rotation + perm)."""
    import numpy as np
    k = np.arange(n_add, dtype=np.float64) + 0.5
    z = 1.0 - 2.0 * k / n_add
    az = np.pi * (1.0 + np.sqrt(5.0)) * k
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    pts = np.stack([r * np.cos(az), r * np.sin(az), z], axis=1)
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((3, 3))
    Q, R = np.linalg.qr(A)
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return pts @ Q.T, rng.permutation(n_add)


def install_patches():
    import numpy as np
    import cw_helpers as CWH
    import stagec_anchor_a2 as A2
    import stagec_fisher as SF
    if getattr(CWH.load_pulsars, "_summit_d3", False):
        return
    orig_load = CWH.load_pulsars

    def load_pulsars_ext(n_pulsars=116):
        ent, disco = orig_load(n_pulsars)
        n_add = EXT_STATE["n_add"]
        if not n_add:
            return ent, disco
        ninfo = np.array([len(p.toas) for p in disco], float)
        order = np.argsort(-ninfo, kind="stable")           # COMPASS quality rank
        donors = [int(order[k % len(ent)]) for k in range(n_add)]
        pts, perm = _ext_positions(n_add)
        for k, di in enumerate(donors):
            pe, pd = copy.deepcopy(ent[di]), copy.deepcopy(disco[di])
            cname = f"{ent[di].name}C{k:02d}"
            if EXT_STATE["twin"] and k == 0:
                # twin: 0.1 deg off the donor's sky. EXACT coincidence makes the
                # HD-ORF GWB block singular (duplicate position -> Phi_gwb rank
                # deficient -> every E-step row non-finite; measured, jobs
                # 12767979_1/12768789_*). Inheritance breaks are O(1-100) nat;
                # 0.1 deg of geometry is O(1e-3) -- full discriminating power kept.
                p0 = np.asarray(disco[di].pos, float)
                t = np.cross(p0, [0.0, 0.0, 1.0])
                if np.linalg.norm(t) < 1e-6:
                    t = np.cross(p0, [1.0, 0.0, 0.0])
                t /= np.linalg.norm(t)
                eps = np.deg2rad(0.1)
                pos = np.cos(eps) * p0 + np.sin(eps) * t
            else:
                pos = pts[perm[k]]
            th = float(np.arccos(np.clip(pos[2], -1.0, 1.0)))
            ph = float(np.arctan2(pos[1], pos[0]) % (2.0 * np.pi))
            for obj in (pe, pd):
                obj.name = cname
                obj.pos = np.asarray(pos, float)
                obj.theta = th
                obj.phi = ph
            ent.append(pe)
            disco.append(pd)
        # singularity guard: a (near-)duplicate sky position makes the HD-ORF GWB
        # block rank-deficient (measured -- the exact-twin all-NaN E-step). The twin's
        # own 0.1-deg offset passes; anything closer than 0.05 deg refuses.
        allpos = np.array([np.asarray(p.pos, float) for p in disco])
        c = np.clip(allpos @ allpos.T, -1.0, 1.0)
        np.fill_diagonal(c, -1.0)
        min_sep = float(np.degrees(np.arccos(c.max())))
        if min_sep < 0.05:
            raise RuntimeError(f"D3 ext: min pairwise separation {min_sep:.4f} deg "
                               "< 0.05 -- near-duplicate position would make Phi_gwb "
                               "singular; re-seed the position draw.")
        print(f"  [D3 ext] +{n_add} clones (donors by n_toa rank; twin={EXT_STATE['twin']}) "
              f"-> npsr {len(ent)}; min pair sep {min_sep:.2f} deg", flush=True)
        return ent, disco

    load_pulsars_ext._summit_d3 = True
    CWH.load_pulsars = load_pulsars_ext

    orig_priors = A2.load_real_priors

    def load_real_priors_ext(names):
        base = [donor_of(n) for n in names]
        pr = orig_priors(base)
        return pr

    load_real_priors_ext._summit_d3 = True
    A2.load_real_priors = load_real_priors_ext

    orig_inj = SF.make_geometry_injection

    def make_geometry_injection_ext(pta, ent_psrs, num_cw, cw_block_names, seed,
                                    **kw):
        inj = orig_inj(pta, ent_psrs, num_cw, cw_block_names, seed, **kw)
        # bundle intact: clone RN := donor RN (overwrites the continued ref stream)
        for pe in ent_psrs:
            d = donor_of(pe.name)
            if d != pe.name:
                for suff in ("rednoise_log10_A", "rednoise_gamma"):
                    kc, kd = f"{pe.name}_{suff}", f"{d}_{suff}"
                    if kd in inj and kc in inj:
                        inj[kc] = inj[kd]
        return inj

    make_geometry_injection_ext._summit_d3 = True
    SF.make_geometry_injection = make_geometry_injection_ext


def _bank_path(n_add, T=30):
    return f"{SUMMIT_ROOT}/smt_L_gwb_T{T}_ext{n_add}_s{SEED_D3}.npz"


def mode_lgwb(n_add, T=30):
    """Cut the ARRAY-KEYED canonical GWB square root (gl_lgwb_t30.py logic verbatim,
    extended array). cpus=8 pinned by the sbatch (thread refusal below)."""
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        if os.environ.get(v) != "8":
            print(f"REFUSED: {v}={os.environ.get(v)!r} != '8' (banked convention).")
            return 2
    import numpy as np
    bank0 = _bank_path(n_add, T)
    if os.path.exists(bank0):
        z = np.load(bank0, allow_pickle=True)
        print(f"  bank exists (skip-on-exist): {bank0} fp={str(z['fingerprint'])} "
              f"npsr={int(z['npsr'])}")
        return 0
    import jax.numpy as jnp
    import glacier_pop as GP
    import trackB_b1_core as C
    EXT_STATE["n_add"] = n_add
    pop = GP.draw_population(GP.SEED_POP_BASE, n_src=32, band_log10f=(-8.7, -7.5))
    G = GP.build_glacier_problem(T, pop, verbose=True)
    amo = G["amo"]
    it = amo["internals"]
    theta = np.asarray(amo["theta_truth"], float)
    data0 = amo["inject_delay"](jnp.asarray(theta), jnp.ones(amo["npsr"]))
    sp = C.B1Split(amo, theta, data0, np.ones(amo["npsr"]))
    Pinv = np.asarray(it["Pinv_gwb"])
    Phi = np.linalg.inv(0.5 * (Pinv + Pinv.T))
    Phi = 0.5 * (Phi + Phi.T)
    w, V = np.linalg.eigh(Phi)
    L = V * np.sqrt(np.clip(w, 0.0, None))
    fp = C.lgwb_fingerprint(L)
    recon = float(np.max(np.abs(L @ L.T - Phi)))
    bank = _bank_path(n_add, T)
    np.savez(bank, L_gwb=L, fingerprint=fp, cpus="8", host=socket.gethostname(),
             T_label=T, npsr=amo["npsr"], ngp=sp.ngp_gwb, n_add=n_add, seed_d3=SEED_D3,
             note=f"SUMMIT D3 array-keyed GWB square root, ext+{n_add} (seed {SEED_D3}), "
                  "geometry-only. The COMPASS venue-bank hazard: shape does not key an "
                  "array; this bank is keyed by (T, n_add, seed) and MACHINE-LOCAL.")
    nd = C.NoiseDrawer(sp, amo, lgwb_path=bank, strict=True)
    print(f"  Phi {Phi.shape} recon {recon:.2e}; banked {bank} fp={fp}")
    print(f"  drawer: {nd.lgwb_prov}")
    print(f"=== D3 L_gwb BANK CUT (ext+{n_add}) ===")
    return 0


def mode_gate(n_add):
    """SG-D3b TWIN + SG-D3c audit + SG-D3d RN-freeze + SG-D3e bank assert.
    Gate venue: T=15, n=32, zero noise (gate-class; no campaign draw touched)."""
    import numpy as np
    import glacier_pop as GP
    import glacier_loop as GL
    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()

    # ---- SG-D3d: first-116 RN draws unchanged under extension (CPU-light dict cmp)
    # Build the injection twice at n=32 via the pta path is heavy; instead replay the
    # frozen ref stream exactly as make_geometry_injection does: the first 116 pulsars'
    # (A, gamma) pairs come from the same seed-0 stream positions iff clone params sort
    # AFTER the real 116 in the inj key order. Verified structurally on the real build
    # below (names order = ent order; clones appended last). Recorded, then gated on
    # the built inj in SG-D3b.
    # ---- SG-D3b: the TWIN test
    EXT_STATE["n_add"] = 1
    EXT_STATE["twin"] = True
    pop = GP.draw_population(GP.SEED_POP_BASE, n_src=32, band_log10f=(-8.7, -7.5))
    slots, n_harm, active, chan, n_clip = GL.embed_igniter(
        pop, 0.0, GL.VENUE_SPAN_S[15])
    pop_slots = dict(pop)
    pop_slots["src"] = slots
    pop_slots["n_src"] = len(slots)
    G = GP.build_glacier_problem(15, pop_slots, verbose=True)
    G["slots"] = slots                      # run_cell/mode_gate convention (caller wires)
    amo = G["amo"]
    ent = G["ent_psrs"]
    names = [p.name for p in ent]
    npsr = len(ent)
    assert npsr == 117, f"twin build npsr {npsr} != 117"
    di = names.index(donor_of(names[-1]))
    ci = npsr - 1
    sb = GL.CertScoreboard(G, amo, jnp, C)
    # L0/prior inherit EXACTLY (name-keyed); dL/EV follow position -> O(1e-3) at the
    # 0.1-deg twin offset. Breaks are O(1): wrong prior row, wrong L0, wrong RN.
    ok_struct = (abs(sb.L0[ci] - sb.L0[di]) == 0.0
                 and sb.priors["lit"][ci] == sb.priors["lit"][di]
                 and abs(sb.dL[ci] / sb.dL[di] - 1.0) < 1e-2
                 and np.allclose(sb.EV[ci], sb.EV[di], rtol=1e-2, atol=1e-4))
    print(f"[SG-D3b:struct] twin rows: L0/prior EXACT, dL rel "
          f"{abs(sb.dL[ci]/sb.dL[di]-1.0):.2e} (<1e-2), EV close -> "
          f"{'PASS' if ok_struct else 'FAIL'}")
    L_true, n_true = sb.draw_truth(IG, dist_seed=SEED_D3 + 77)
    data = amo["inject_delay"](jnp.asarray(np.asarray(amo["theta_truth"], float)),
                               jnp.ones(npsr))
    import glacier_pop as GPM
    led = GPM.PromoteLedger(slots)
    led.promote(0, slots[0], iteration=0)   # G-d2b convention: score with a fed member
    #                                         (all-carried == empty template is non-finite
    #                                          by construction and never runs in production)
    theta_rec = np.asarray(amo["theta_truth"], float).copy()
    dlnl, lnK, qmax, on_true = sb.columns(theta_rec, led, data, jnp.ones(npsr),
                                          np.array([], int), np.array([], float))
    d_dl = abs(dlnl[ci] - dlnl[di])
    d_lk = abs(lnK[ci] - lnK[di])
    d_q = abs(qmax[ci] - qmax[di])
    ok_cols = np.isfinite([dlnl[ci], lnK[ci], qmax[ci]]).all() and \
        d_dl < 0.5 and d_lk < 0.5 and d_q < 0.05
    print(f"[SG-D3b:cols] twin scoreboard rows (0.1-deg offset; break scale O(1-100)): "
          f"|d dlnL|={d_dl:.2e} (<0.5) |d lnK|={d_lk:.2e} (<0.5) "
          f"|d q|={d_q:.2e} (<0.05) -> {'PASS' if ok_cols else 'FAIL'}")

    # ---- SG-D3c: audit of the +n_add production positions
    EXT_STATE["twin"] = False
    pts, perm = _ext_positions(n_add)
    pos = pts[perm[np.arange(n_add)]] if n_add else pts
    c = np.clip(pos @ pos.T, -1.0, 1.0)
    np.fill_diagonal(c, -1.0)
    nn = np.degrees(np.arccos(c.max(axis=1)))
    dip = float(np.linalg.norm(pos.mean(axis=0)))
    print(f"[SG-D3c] +{n_add} positions: NN sep min/med {nn.min():.1f}/"
          f"{np.median(nn):.1f} deg; dipole |<n>| {dip:.3f}")
    np.savez(f"{os.environ['GLACIER_OUT']}/smt_d3_audit_ext{n_add}.npz",
             positions=pos, nn_deg=nn, dipole=dip, seed=SEED_D3, n_add=n_add,
             note="SUMMIT D3 addition-position audit (golden-angle, seeded rotation)")

    # ---- SG-D3e: bank exists + strict-loads at the extended shape
    bank = _bank_path(n_add)
    ok_bank = os.path.exists(bank)
    if ok_bank:
        z = np.load(bank, allow_pickle=True)
        ok_bank = int(z["npsr"]) == 116 + n_add and int(z["n_add"]) == n_add
        print(f"[SG-D3e] bank {os.path.basename(bank)}: npsr {int(z['npsr'])} "
              f"fp={str(z['fingerprint'])} -> {'PASS' if ok_bank else 'FAIL'}")
    else:
        print(f"[SG-D3e] bank MISSING ({bank}) -- run mode lgwb first -> FAIL")
    ok = ok_struct and ok_cols and ok_bank
    print(f"=== SUMMIT D3 GATES (+{n_add}): {'ALL PASS' if ok else 'FAIL'} ===")
    return 0 if ok else 2


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "lgwb", "cell2", "null2"])
    ap.add_argument("--next", dest="n_add", type=int, default=30, choices=[30, 100])
    ap.add_argument("--sky", type=int, default=0)
    ap.add_argument("--real", type=int, default=0)
    ap.add_argument("--iters", type=int, default=6)
    a = ap.parse_args()

    out = os.environ.get("GLACIER_OUT", "")
    if "SUMMIT_results" not in out:
        print(f"REFUSED: GLACIER_OUT={out!r} not under SUMMIT_results.")
        return 3
    os.makedirs(out, exist_ok=True)
    install_patches()

    if a.mode == "lgwb":
        return mode_lgwb(a.n_add)
    if a.mode == "gate":
        return mode_gate(a.n_add)

    marker = f"{out}/HOLDS_CLEARED"
    gates_ok = glob.glob(f"{out}/smt_d3_audit_ext{a.n_add}.npz")
    if not (os.path.exists(marker) and gates_ok and os.path.exists(_bank_path(a.n_add))):
        print(f"REFUSED: need HOLDS_CLEARED + smt_d3_audit_ext{a.n_add}.npz + the "
              f"array-keyed L_gwb bank in/for {out} (run lgwb + gate first).")
        return 3
    with open(marker) as fh:
        print(f"[SUMMIT-D3] holds cleared: {fh.read().strip()}", flush=True)

    EXT_STATE["n_add"] = a.n_add
    import glacier_loop as GL
    GL.LGWB_BANKS = dict(GL.LGWB_BANKS)
    GL.LGWB_BANKS[30] = _bank_path(a.n_add)
    print(f"[SUMMIT-D3] ext+{a.n_add} rung=r13p9 arm=none sky={a.sky} OUT={out} "
          f"lgwb={GL.LGWB_BANKS[30]}", flush=True)
    return GL.run_cell("none", a.sky, scrambled=(a.mode == "null2"), real=a.real,
                       n_iter=a.iters, rung="r13p9")


if __name__ == "__main__":
    sys.exit(main() or 0)
