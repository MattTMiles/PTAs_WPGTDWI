"""COMPASS -- THE ISOTROPIC CONTROL ARRAY: separating array-geometry effects from
source-sky and physics effects, for the first time.

PURE EXECUTION (agent COMPASS, ACCRE, general H200 lane only -- the reserved dgx share,
the GLACIER pools and their spillover are untouched). No tracked-file edits, no commits
(stage-only). Results -> COMPASS_results/ (lean keyed npz, raw statistics ALWAYS banked).
Report reports/COMPASS.md.

THE FAKE ARRAY (the BUILD, gated before any spend)
---------------------------------------------------
116 pulsars at quasi-uniform sky positions (Fibonacci-sphere golden-angle grid, a
seed-controlled random rotation + a seed-controlled random enumeration of the grid
points), carrying the REAL array's noise properties by a RANK-PRESERVING map: the k-th
best real pulsar's noise bundle -> the k-th enumerated position. The BUNDLE travels
INTACT -- TOA epochs/errors, timing-model design matrix, the frozen per-pulsar red-noise
assignment, the EM distance (L0 = pdist) and the lit/vlbi prior widths all stay attached
to their pulsar; the sky position (pos/theta/phi, on BOTH the enterprise and discovery
objects) is the ONLY thing changed. Consequences, all automatic and all position-
consistent because the patch runs BEFORE any model build: CW antenna patterns +
pulsar-term geometry (injection AND recovery), the HD-ORF GWB correlation (Phi_gwb),
and every fringe spacing dL_a = c/[f_gw(1 - cos mu)].

  - QUALITY RANK: descending TOA count. The harness flattens the discovery-side TOA
    errors to 1 us (cw_helpers.load_pulsars), so the per-pulsar white-noise information
    that actually enters the likelihood is proportional to n_toa; the real (feather)
    median TOA error is banked beside it for the record.
  - "Distances drawn from the real distance distribution" holds EXACTLY: the empirical
    distance/prior census of the iso array IS the real array's (bundles intact).
  - SCOPE (stated here and in the report): the rank-map assignment is ONE choice of
    noise-geometry pairing; arm S1 re-runs the survivor cell under a second,
    shuffled assignment (ISO_SEED_SHUF) to price that choice. The timing-model design
    matrix (incl. its astrometry columns, fit at the ORIGINAL position) travels with
    the bundle -- the marginalisation structure is treated as a noise property, not a
    sky property. Second-order by construction; stated, not hidden.

THE L_gwb VENUE-BANK HAZARD, CLOSED BEFORE IT FIRES: the GENERALISE banks
gen_L_gwb_n{5336,6960}.npz are keyed by SHAPE ONLY, and the iso array's Phi_gwb has the
SAME shape with DIFFERENT (position-dependent) values -- the stock loader would silently
draw GWB noise with the REAL array's correlations. The compass loader routes by iso
state: iso OFF -> the GENERALISE venue bank (bit-compatible with the banked campaign;
required by gate gA); iso ON -> COMPASS_results/comp_L_gwb_n{n}_iso{seed}.npz, keyed by
the assignment seed (arm S1's shuffled pairing changes Phi_gwb and gets its OWN bank).
Every realisation banks lgwb_prov; the run loop ASSERTS the loaded bank matches the
unit's iso state before any realisation is written.

MACHINERY REUSED, NOT REIMPLEMENTED: hpc_harbor/generalise/generalise.py's
run_gen_realisation / ledger_unit / floor_v22 / _offender / _counts_at / grade_of are
used VERBATIM (module-level path+plan hooks repointed at COMPASS_results) -- the scoring
path is the one GENERALISE gates g1/g3 already licensed. Criterion-v2.2 floors
throughout (Gumbel valid only at zero-fraction <= 20%, else empirical q95 +/- bootstrap;
zero-fraction a REQUIRED column; floors REFIT per unit -- per sky in arm I1/S1, never
carried). Grades use the LOTTERY two-sided band rule. Host/thread provenance banked per
npz; floors + signal of a unit ride the same task/process (8-BLAS-thread convention,
2-way GPU packing).

ARMS
----
  gA (gate, iso OFF): the compass path re-runs 5 nulls + 5 signal realisations of the
     banked A-SKY survivor unit (AS_e03_h1275_k5_s4) with the banked seeds. BINDING:
     certified + detected sets equal at the unit's banked adopted floor and total certs
     within +/-2 (the g3 count convention); bit-identity REPORTED beside it (same lane
     -> expected, but cross-host class if the node differs; the g3 amended near-bar
     benign class applies then). Proves the wrapper adds nothing.
  g0 (gate, iso ON, THE MISSION GATE): at zero signal the iso array's floor machinery
     runs one cell x N=100 nulls (the survivor cell at GEO sky 4 -- these ARE arm I1's
     unit nulls, reused by skip-on-exist). BINDING: zero-fraction <= 20% (the estimator
     branch of 30/32 real A-SKY units at these loudnesses) and a finite, sane Gumbel
     MLE. REPORTED: the floor value against the real cell's per-sky floor band
     (59.9-259.8 nat) -- a level SHIFT is a measurement, not a machinery failure.
     Plus the array audit: nearest-neighbour separation stats, dipole moment, noise
     census identity (sorted n_toa / L0 / prior multisets EQUAL to the real array's).
  I1 -- THE ATTRIBUTION EXPERIMENT (the headline): the 4 A-SKY boundary cells
     (incl. the survivor) x 8 skies (GEO draws 4-11, sky-paired with the real
     ensemble; placements via the SAME as_placement map) x (100 nulls + 15 signal),
     per-sky v2.2 floors. PRE-REGISTERED READING, per cell: the sky-to-sky scatter
     ratio rho = sd_sky(iso count)/sd_sky(real count), real side from the banked
     gen_ledger_AS_* files. rho < 0.5 -> COLLAPSES ("a uniform array de-lotteries the
     switch": the lottery is array anisotropy); rho >= 0.75 -> PERSISTS ("no array
     design escapes it": intrinsic source-geometry variance); between -> MARGINAL
     (band). Primary = raw ratio (identical protocol both sides; the 15-real sampling
     noise is common-mode); a sampling-noise-decomposed ratio is REPORTED beside it,
     never substituted. The survivor's 6/8-sky verdict is re-graded on iso in Arm-C
     vocabulary, and the per-cell floor sky-swing (max/min) rides beside the counts.
  I2 -- THE CONSTANTS RE-MEASURED (cheap, piggybacked; pooled CHORUS protocol --
     5 skies x 6 noise + 100 pooled nulls -- so the census cells are read against the
     banked RECUT table):
     (a) channel budget: m1/m2/m3 e=0.3 lit 3+13 at census h=-13.25 (the MAGPIE-J1
         equal-kappa contrast, kappa fixed at 2.65, channels 23->30->37: does the
         ON flip between n_ecc=1 and 2 survive isotropy? real: 0.70/2.77/2.50) and
         m2/m3 e=0.3 at the survivor operating point (h=-12.75, 5+11; m1 = the I1
         survivor cell itself, 23 channels).
     (b) union-ceiling analogue (CPU, from banks): union certified set + per-pulsar
         certification frequency over the I1 ensemble, iso vs real (SAME protocol both
         sides -- NOT comparable to GEO's zero-noise union-18 or GALLERY's 30-readable
         SNR_pterm census, which are different machineries; scope line in the report).
     (c) one floor-law column (CPU, from I1 ledgers): mean per-sky floor vs h over
         {(-13.0, e05, 3+13), (-12.75, e03, 3+13), (-12.5, e03, 3+13)} -- the SAME
         e-mixture column as the real A-SKY floor readout -- does 1.66 survive?
  S1 -- SHUFFLED-ASSIGNMENT ROBUSTNESS (budget-gated): the survivor cell x 8 skies
     under ISO_SEED_SHUF (same grid, different pulsar->position pairing, own L_gwb
     bank). Reading: |count/floor shift| vs the I1 survivor row, in units of its sky
     scatter -- prices the rank-map choice.

SEEDS (disjoint from every banked range -- IGNITE <= 4.xM, IGNITE-2 5/6.xM, CHORUS
7.0-8.5M, FORGE/ATLAS 9.xM, SURFACE 20/40M (+10M -> 30/50M), GENERALISE 61-68.xM
(+10M -> 71-78.xM), GLACIER 500M/900M):
  I1 signal 81.xM, nulls 82.xM; I2 signal 83.xM, nulls 84.xM; S1 signal 85.xM,
  nulls 86.xM; dist = noise + 10_000_000 (house convention) -> dists in 91-96.xM.
ISO_SEED_MAIN = 8100 (array build), ISO_SEED_SHUF = 8101 (S1). Not noise seeds; they
key the position assignment and the L_gwb banks only.

BUDGET: header estimate printed by the gate job from measured warm walls, BEFORE any
wide submit. Mission class ~A-SKY x2 (est 3-6 GPU-hr); HARD STOP at 20 GPU-hr.
Trim ladder if the header breaches: drop S1, then the I2 operating-point pair, then
I2 entirely (I1 + g0 are the mission and are never trimmed below 8 skies).
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time, argparse, glob, socket
import numpy as np

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
sys.path.insert(0, "/home/mattm/projects/HSYMT/hpc_harbor/forge")
sys.path.insert(0, "/home/mattm/projects/HSYMT/hpc_harbor/ignite")
sys.path.insert(0, "/home/mattm/projects/HSYMT/hpc_harbor/atlas")
sys.path.insert(0, "/home/mattm/projects/HSYMT/hpc_harbor/chorus")
sys.path.insert(0, "/home/mattm/projects/HSYMT/hpc_harbor/generalise")

HSYMT = "/home/mattm/projects/HSYMT"
OUT = f"{HSYMT}/COMPASS_results"
GEN_OUT = f"{HSYMT}/GENERALISE_results"

import chorus as CH                      # jax-free at import
import generalise as GEN                 # the gated scoring machinery, reused verbatim

I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)

NPSR = 116
ISO_SEED_MAIN = 8100
ISO_SEED_SHUF = 8101
ISO_STATE = {"seed": None}               # None = real array (gate gA); else iso seed

H_CENSUS = GEN.H_CENSUS                  # -13.25
SKIES_AB = GEN.SKIES_AB                  # [-1, 0, 1, 2, 3]
C_SKIES = GEN.C_SKIES                    # GEO draws 4-11 (sky-paired with Arm C / A-SKY)
N_NULLN = GEN.N_NULLN                    # 100
N_SIG_SKY = GEN.N_SIG_C                  # 15 signal / sky (I1, S1)
N_NOISE_AB = GEN.N_NOISE_AB              # 6 (I2 pooled protocol)
D_OFF = GEN.D_OFF                        # +10M dist offset

# the A-SKY boundary cells, SAME order as GEN.AS_CELLS -> ci-paired with the real bank
I1_CELLS = list(GEN.AS_CELLS)            # [(-12.75,.3,3), (-12.75,.3,5), (-12.5,.3,3), (-13.0,.5,3)]
SURV = (-12.75, 0.3, 5)                  # the A-SKY survivor (6/8 skies on the real array)
SURV_CI = I1_CELLS.index(SURV)

# I2 pooled cells: (n_ecc, h, k_loud); e=0.3, lit, T=30 throughout
I2_CELLS = [(1, H_CENSUS, 3), (2, H_CENSUS, 3), (3, H_CENSUS, 3),
            (2, -12.75, 5), (3, -12.75, 5)]

# seeds
I1_SIG, I1_NUL = 81_000_000, 82_000_000
I2_SIG, I2_NUL = 83_000_000, 84_000_000
S1_SIG, S1_NUL = 85_000_000, 86_000_000

# real-side comparator: the banked A-SKY ledgers (per-sky floors + counts)
REAL_AS_LEDGER = GEN_OUT + "/gen_ledger_{tag}.npz"      # tag arrives WITH its AS_ prefix

STOP_GPU_HR = 20.0


def hkey(h):
    return GEN.hkey(h)


# ============================================================
# THE ISOTROPIC ARRAY (build + audit)
# ============================================================
def fibonacci_sphere(n=NPSR):
    """Golden-angle spiral: quasi-uniform unit vectors, index 0 at the +z cap."""
    k = np.arange(n, dtype=np.float64) + 0.5
    z = 1.0 - 2.0 * k / n
    az = np.pi * (1.0 + np.sqrt(5.0)) * k
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    return np.stack([r * np.cos(az), r * np.sin(az), z], axis=1)


def iso_positions(seed):
    """The grid + the seed-controlled assignment order. The rotation kills any
    alignment between the spiral axis and equatorial coordinates; the enumeration
    permutation makes 'the k-th position' a random point of the uniform grid (a
    rank->grid-index map would pile the best pulsars into one polar cap -- exactly
    the anisotropy this array exists to remove)."""
    pts = fibonacci_sphere()
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((3, 3))
    Q, R = np.linalg.qr(A)
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    pts = pts @ Q.T
    perm = rng.permutation(len(pts))
    return pts, perm


def quality_rank(disco_psrs):
    """Descending TOA count (see module docstring). Stable tie-break by load order."""
    ninfo = np.array([len(p.toas) for p in disco_psrs], float)
    order = np.argsort(-ninfo, kind="stable")
    rank = np.empty(len(ninfo), int)
    rank[order] = np.arange(len(ninfo))
    return rank, ninfo


def isotropize(ent_psrs, disco_psrs, seed, verbose=True):
    """Reposition the loaded pulsars onto the iso grid. Bundle-preserving: ONLY
    pos/theta/phi change, on BOTH object families. Returns the audit record."""
    pts, perm = iso_positions(seed)
    rank, ninfo = quality_rank(disco_psrs)
    old_pos = np.array([np.asarray(p.pos, float) for p in disco_psrs])
    new_pos = np.empty_like(old_pos)
    for i, (pe, pd) in enumerate(zip(ent_psrs, disco_psrs)):
        pos = pts[perm[rank[i]]]
        th = float(np.arccos(np.clip(pos[2], -1.0, 1.0)))
        ph = float(np.arctan2(pos[1], pos[0]) % (2.0 * np.pi))
        for obj in (pe, pd):
            obj.pos = np.asarray(pos, float)
            obj.theta = th
            obj.phi = ph
        new_pos[i] = pos
    med_err = np.array([float(np.median(np.asarray(p.toaerrs, float)))
                        for p in ent_psrs])
    rec = dict(seed=seed, positions=new_pos, old_positions=old_pos, perm=perm,
               rank=rank, n_toa=ninfo.astype(int), med_toaerr_real=med_err,
               names=np.array([p.name for p in ent_psrs]))
    if verbose:
        print(f"  [ISO seed={seed}] repositioned {len(ent_psrs)} pulsars "
              f"(rank metric: n_toa desc; best={rec['names'][np.argmax(ninfo)]} "
              f"n_toa={int(ninfo.max())})", flush=True)
    return rec


def nn_sep_deg(pos):
    """Nearest-neighbour angular separations (deg) -- the uniformity audit."""
    c = np.clip(pos @ pos.T, -1.0, 1.0)
    np.fill_diagonal(c, -1.0)
    return np.degrees(np.arccos(c.max(axis=1)))


# ============================================================
# runtime patches (chorus-C7 pattern: this process only, no tracked file touched)
# ============================================================
_PATCHED = {"done": False}
_ISO_REC = {}                            # seed -> audit record of the LAST build


def _install_patches(C):
    """(1) cw_helpers.load_pulsars -> isotropize when ISO_STATE['seed'] is set.
    (2) trackB_b1_core.load_or_build_L_gwb -> route by iso state (docstring)."""
    if _PATCHED["done"]:
        return
    import cw_helpers as CWH
    orig_load = CWH.load_pulsars

    def _comp_load(n_pulsars=10):
        ent, disco = orig_load(n_pulsars)
        if ISO_STATE["seed"] is not None:
            _ISO_REC[ISO_STATE["seed"]] = isotropize(ent, disco, ISO_STATE["seed"])
        return ent, disco

    CWH.load_pulsars = _comp_load

    orig_lgwb = C.load_or_build_L_gwb

    def _comp_lgwb(Phi_gwb, path=None, strict=False):
        n = int(Phi_gwb.shape[0])
        if ISO_STATE["seed"] is None:
            bank = f"{GEN_OUT}/gen_L_gwb_n{n}.npz"
            if not os.path.exists(bank):
                raise RuntimeError(f"real-array venue bank missing: {bank} -- "
                                   "gate gA needs the GENERALISE bank, never a recompute")
        else:
            bank = f"{OUT}/comp_L_gwb_n{n}_iso{ISO_STATE['seed']}.npz"
            if not os.path.exists(bank):
                nth = (os.environ.get("OPENBLAS_NUM_THREADS")
                       or os.environ.get("OMP_NUM_THREADS")
                       or str(len(os.sched_getaffinity(0))))
                t0 = time.time()
                w, V = np.linalg.eigh(Phi_gwb)
                L = V * np.sqrt(np.clip(w, 0.0, None))
                fp = C.lgwb_fingerprint(L)
                tmp = f"{bank}.{os.getpid()}.tmp.npz"
                np.savez(tmp, L_gwb=L, fingerprint=fp, cpus=nth,
                         host=socket.gethostname(), iso_seed=ISO_STATE["seed"])
                os.replace(tmp, bank)
                print(f"[COMP lgwb] banked {bank} fp={fp} threads={nth} "
                      f"host={socket.gethostname()} ({time.time()-t0:.0f}s)", flush=True)
        return orig_lgwb(Phi_gwb, path=bank, strict=strict)

    C.load_or_build_L_gwb = _comp_lgwb
    C._gen_lgwb_patched = True           # GEN._import_stack must NOT re-patch on top
    _PATCHED["done"] = True


def _import_stack():
    jax, jnp, C, B1, TE, F, IG = CH._import_stack()
    _install_patches(C)
    return jax, jnp, C, B1, TE, F, IG


def build_problem(n_ecc, T_label, iso_seed, verbose=True):
    """CH.build_chorus_problem under the compass patches, with the lgwb-prov guard."""
    _import_stack()
    ISO_STATE["seed"] = iso_seed
    P = CH.build_chorus_problem(n_ecc, T_label=T_label, verbose=verbose)
    prov = str(P["nd"].lgwb_prov)
    want = "comp_L_gwb" if iso_seed is not None else "gen_L_gwb"
    assert want in prov, (f"L_gwb routing violated: iso_seed={iso_seed} but "
                          f"lgwb_prov={prov} -- STOP (the venue-bank hazard)")
    if iso_seed is not None:
        P["iso_seed"] = iso_seed
        P["iso_rec"] = _ISO_REC[iso_seed]
    else:
        P["iso_seed"] = -1
    return P


# ============================================================
# units + plans (GEN's granule convention: a unit = one floor-refit granule)
# ============================================================
def i1_units():
    out = []
    for ci, (h, e, k) in enumerate(I1_CELLS):
        ed = f"e{int(round(e * 10)):02d}"
        for si, sky in enumerate(C_SKIES):
            out.append(dict(arm="I1", e=e, edist=ed, h=h, T=30, k=k, tier="lit",
                            fdiv=1, n_ecc=1, sky=sky, ci=ci, si=si,
                            placement=GEN.as_placement(k, sky),
                            iso=ISO_SEED_MAIN,
                            tag=f"I1_{ed}_h{hkey(h)}_k{k}_s{sky}"))
    return out


def i2_units():
    out = []
    for ci, (n_ecc, h, k) in enumerate(I2_CELLS):
        out.append(dict(arm="I2", e=0.3, edist="e03", h=h, T=30, k=k, tier="lit",
                        fdiv=1, n_ecc=n_ecc, ci=ci,
                        placement=tuple(range(n_ecc)),     # the banked CHORUS convention
                        iso=ISO_SEED_MAIN,
                        tag=f"I2_m{n_ecc}e03_h{hkey(h)}_k{k}"))
    return out


def s1_units():
    h, e, k = SURV
    ed = f"e{int(round(e * 10)):02d}"
    out = []
    for si, sky in enumerate(C_SKIES):
        out.append(dict(arm="S1", e=e, edist=ed, h=h, T=30, k=k, tier="lit",
                        fdiv=1, n_ecc=1, sky=sky, ci=0, si=si,
                        placement=GEN.as_placement(k, sky),
                        iso=ISO_SEED_SHUF,
                        tag=f"S1_{ed}_h{hkey(h)}_k{k}_s{sky}_shuf"))
    return out


def units_for(arm):
    return {"I1": i1_units(), "I2": i2_units(), "S1": s1_units()}[arm]


def comp_unit_entries(arm, ui, u):
    """Nulls FIRST (the floors are the verdicts), then signal. Idempotent."""
    ent = []
    if arm == "I2":                       # pooled CHORUS protocol (5 skies x 6 + 100)
        for rep in range(N_NULLN):
            ns = I2_NUL + 10_000 * ui + rep
            ent.append(dict(kind="fnullN", geo_id=SKIES_AB[rep % len(SKIES_AB)],
                            noise_seed=ns, dist_seed=ns + D_OFF, no_cw=True))
        for si, gid in enumerate(SKIES_AB):
            for rep in range(N_NOISE_AB):
                ns = I2_SIG + 10_000 * ui + 100 * si + rep
                ent.append(dict(kind="sig", geo_id=gid, noise_seed=ns,
                                dist_seed=ns + D_OFF, no_cw=False))
    else:                                 # I1 / S1: per-sky units
        nul0, sig0 = (I1_NUL, I1_SIG) if arm == "I1" else (S1_NUL, S1_SIG)
        for rep in range(N_NULLN):
            ns = nul0 + 10_000 * ui + rep
            ent.append(dict(kind="fnullN", geo_id=u["sky"], noise_seed=ns,
                            dist_seed=ns + D_OFF, no_cw=True))
        for rep in range(N_SIG_SKY):
            ns = sig0 + 10_000 * ui + rep
            ent.append(dict(kind="sig", geo_id=u["sky"], noise_seed=ns,
                            dist_seed=ns + D_OFF, no_cw=False))
    return ent


def comp_real_path(u, e):
    return f"{OUT}/comp_{e['kind']}_{u['tag']}_g{e['geo_id']}_n{e['noise_seed']}.npz"


def comp_ledger_path(u):
    return f"{OUT}/comp_ledger_{u['tag']}.npz"


def _repoint_gen():
    """Repoint GEN's module-level hooks at COMPASS_results + compass plans. Everything
    downstream (run_gen_realisation, ledger_unit, floor_v22, ...) is then GEN verbatim."""
    GEN.OUT = OUT
    GEN.real_path = comp_real_path
    GEN.ledger_path = comp_ledger_path
    GEN.unit_entries = comp_unit_entries
    GEN.units_for = units_for


_repoint_gen()


# ============================================================
# run mode (the fan)
# ============================================================
def mode_run(args):
    jax, jnp, C, B1, TE, F, IG = _import_stack()
    print(f"jax {jax.__version__}, {jax.devices()} host={socket.gethostname()} "
          f"threads={os.environ.get('OPENBLAS_NUM_THREADS', 'auto')}", flush=True)
    os.makedirs(OUT, exist_ok=True)
    units = units_for(args.arm)
    blocks = np.array_split(np.arange(len(units)), args.nshards)
    mine = [int(i) for i in blocks[args.shard]]
    n_todo = 0
    for ui in mine:
        u = units[ui]
        ent = comp_unit_entries(args.arm, ui, u)
        n_todo += sum(0 if os.path.exists(comp_real_path(u, e)) else 1 for e in ent)
    n_tot = sum(len(comp_unit_entries(args.arm, ui, units[ui])) for ui in mine)
    print(f"arm {args.arm}: {len(units)} units, shard {args.shard}/{args.nshards} -> "
          f"units {mine}; already banked: {n_tot - n_todo}; to run: {n_todo}", flush=True)
    if not n_todo:
        for ui in mine:
            GEN.ledger_unit(units[ui], ui)
        print("NOTHING TO DO (all banked)", flush=True)
        return 0
    by_shape = {}
    for ui in mine:
        u = units[ui]
        by_shape.setdefault((u["n_ecc"], u["T"], u["iso"]), []).append(ui)
    gids = sorted(set(e["geo_id"] for ui in mine
                      for e in comp_unit_entries(args.arm, ui, units[ui])))
    skies = {g: (None if g < 0 else F.load_geo_skies([g])[0]) for g in gids}
    for (n_ecc, T, iso), uis in sorted(by_shape.items()):
        P = build_problem(n_ecc, T, iso)
        print(f"  [shape n_ecc={n_ecc} T={T} iso={iso}] lgwb: {P['nd'].lgwb_prov}",
              flush=True)
        for ui in uis:
            u = units[ui]
            t0 = time.time(); nran = 0
            for e in comp_unit_entries(args.arm, ui, u):
                te = time.time()
                tag, msg, _ = GEN.run_gen_realisation(P, C, jnp, u, e, skies[e["geo_id"]])
                if msg != "skip":
                    nran += 1
                    if nran <= 3 or nran % 25 == 0:
                        print(f"    [{nran}] {tag:52s} {msg} ({time.time()-te:.1f}s)",
                              flush=True)
            print(f"  unit {u['tag']}: {nran} run ({time.time()-t0:.0f}s)", flush=True)
            GEN.ledger_unit(u, ui)
        del P
        import gc
        gc.collect()
    print("SHARD COMPLETE", flush=True)
    return 0


def mode_ledgers(args):
    os.makedirs(OUT, exist_ok=True)
    arms = [args.arm] if args.arm else ["I1", "I2", "S1"]
    bad = 0
    for arm in arms:
        print(f"=== ledgers: arm {arm} ===", flush=True)
        for ui, u in enumerate(units_for(arm)):
            if GEN.ledger_unit(u, ui) is None:
                bad += 1
    return 1 if bad else 0


def mode_status(args):
    for arm in ["I1", "I2", "S1"]:
        units = units_for(arm)
        done_r, tot_r, done_u = 0, 0, 0
        for ui, u in enumerate(units):
            ent = comp_unit_entries(arm, ui, u)
            have = sum(1 for e in ent if os.path.exists(comp_real_path(u, e)))
            done_r += have; tot_r += len(ent)
            if os.path.exists(comp_ledger_path(u)):
                done_u += 1
        print(f"arm {arm}: {done_r}/{tot_r} realisations, {done_u}/{len(units)} ledgers",
              flush=True)
    return 0


# ============================================================
# GATES (one job, STOP on failure, before any wide submit)
# ============================================================
def gate_gA(C, F, jnp):
    """iso OFF: the compass path re-runs banked A-SKY survivor realisations."""
    print("\n-- gA: compass path (iso OFF) vs banked AS survivor unit "
          "(AS_e03_h1275_k5_s4) --", flush=True)
    h, e, k = SURV
    sky = 4
    as_ui = SURV_CI * len(C_SKIES) + 0          # GEN's AS unit index for (ci=1, s4)
    u = dict(arm="GA", e=e, edist="e03", h=h, T=30, k=k, tier="lit", fdiv=1,
             n_ecc=1, sky=sky, placement=GEN.as_placement(k, sky), iso=None,
             tag="gA_e03_h1275_k5_s4")
    led_p = REAL_AS_LEDGER.format(tag=f"AS_e03_h{hkey(h)}_k{k}_s{sky}")
    if not os.path.exists(led_p):
        print(f"   FAIL: banked ledger missing: {led_p}", flush=True)
        return False, np.nan
    led = np.load(led_p, allow_pickle=True)
    fl = float(led["floor_adopted"])
    banked = []
    for kind, base, nrep in (("fnullN", GEN.AS_NUL, 5), ("sig", GEN.AS_SIG, 5)):
        for rep in range(nrep):
            ns = base + 10_000 * as_ui + rep
            p = f"{GEN_OUT}/gen_{kind}_AS_e03_h{hkey(h)}_k{k}_s{sky}_g{sky}_n{ns}.npz"
            if not os.path.exists(p):
                print(f"   FAIL: banked realisation missing: {p}", flush=True)
                return False, np.nan
            banked.append((kind, ns, np.load(p, allow_pickle=True)))
    sky_src = F.load_geo_skies([sky])[0]
    ISO_STATE["seed"] = None
    P = build_problem(1, 30, None)
    times = []
    strict_bad, tot_b, tot_n, sets_ok = 0, 0, 0, True
    for kind, ns, b in banked:
        ent = dict(kind="ga", geo_id=sky, noise_seed=ns,
                   dist_seed=ns + D_OFF, no_cw=(kind == "fnullN"))
        t0 = time.time()
        tag, msg, res = GEN.run_gen_realisation(P, C, jnp, u, ent, sky_src, keep=False)
        times.append(time.time() - t0)
        dl_b = np.nan_to_num(np.asarray(b["dlnL_det"], float), posinf=1e30)
        dl_n = np.nan_to_num(np.asarray(res["dlnL"], float), posinf=1e30)
        d = float(np.max(np.abs(dl_n - dl_b)))
        if d > 1e-9:
            strict_bad += 1
        lnKb = np.asarray(b["lnK"], float)
        cb = (dl_b > np.maximum(lnKb, fl)) & (np.asarray(b["qmax"]) > GEN.QBAR)
        cn = (dl_n > np.maximum(res["lnK"], fl)) & (res["qmax"] > GEN.QBAR)
        db_ = dl_b > np.maximum(lnKb, fl); dn_ = dl_n > np.maximum(res["lnK"], fl)
        tot_b += int(cb.sum()); tot_n += int(cn.sum())
        same = bool(np.array_equal(cb, cn) and np.array_equal(db_, dn_))
        if not same:
            bad = np.where((cb != cn) | (db_ != dn_))[0]
            benign = True
            for a in bad:
                bar = max(lnKb[a], fl)
                near = abs(dl_b[a] - bar) < 2.0
                deep = (dl_b[a] < bar - 10.0) and (dl_n[a] < bar - 10.0)
                benign &= (near or deep)
                print(f"     flip psr[{a}]: dlnL {dl_b[a]:.3f} -> {dl_n[a]:.3f} "
                      f"(bar {bar:.2f}) -> {'NEAR-BAR benign' if (near or deep) else 'NOT benign'}",
                      flush=True)
            sets_ok &= benign
        print(f"   {kind} n{ns}: max|d dlnL|={d:.3e} cert {int(cb.sum())}->{int(cn.sum())}"
              f" sets_equal={same} ({times[-1]:.1f}s)", flush=True)
    count_ok = abs(tot_n - tot_b) <= 2
    prov = str(P["nd"].lgwb_prov)
    print(f"   lgwb prov (must be the GENERALISE venue bank): {prov}", flush=True)
    print(f"   strict bit-identity {10 - strict_bad}/10 (REPORTED, not required); "
          f"total certs banked {tot_b} vs new {tot_n} (|d|<=2: {count_ok})", flush=True)
    ok = sets_ok and count_ok and ("gen_L_gwb" in prov)
    print(f"   gA: {'PASS' if ok else 'FAIL'}", flush=True)
    del P
    return ok, float(np.median(times[2:]))


def gate_g0(C, F, jnp):
    """iso ON: the mission gate -- one cell x N=100 nulls, floor machinery in class."""
    print("\n-- g0: iso array (seed 8100), survivor cell, sky 4, N=100 nulls --",
          flush=True)
    units = i1_units()
    ui = SURV_CI * len(C_SKIES) + 0
    u = units[ui]
    assert u["sky"] == 4 and u["k"] == 5, u
    P = build_problem(1, 30, ISO_SEED_MAIN)
    rec = P["iso_rec"]
    # ---- array audit ----
    nn_new = nn_sep_deg(rec["positions"]); nn_old = nn_sep_deg(rec["old_positions"])
    dip_new = float(np.linalg.norm(rec["positions"].mean(axis=0)))
    dip_old = float(np.linalg.norm(rec["old_positions"].mean(axis=0)))
    print(f"   NN sep (deg) iso: min/med/max = {nn_new.min():.2f}/"
          f"{np.median(nn_new):.2f}/{nn_new.max():.2f} | real: {nn_old.min():.2f}/"
          f"{np.median(nn_old):.2f}/{nn_old.max():.2f}", flush=True)
    print(f"   dipole |<pos>|: iso {dip_new:.4f} vs real {dip_old:.4f}", flush=True)
    audit_ok = (np.median(nn_new) > np.median(nn_old)) and (dip_new < dip_old) \
        and (nn_new.min() > 2.0)
    np.savez(f"{OUT}/comp_isoarray_{ISO_SEED_MAIN}.npz", **rec,
             nn_sep_deg=nn_new, nn_sep_deg_real=nn_old,
             dipole=dip_new, dipole_real=dip_old,
             rank_metric=np.array("n_toa descending (disco toaerrs flattened to 1e-6; "
                                  "info ~ n_toa)"))
    # ---- the 100 nulls (these ARE arm I1's unit nulls; skip-on-exist) ----
    sky_src = F.load_geo_skies([4])[0]
    times = []
    for e in comp_unit_entries("I1", ui, u):
        if e["kind"] != "fnullN":
            continue
        t0 = time.time()
        GEN.run_gen_realisation(P, C, jnp, u, e, sky_src)
        times.append(time.time() - t0)
    off = []
    for e in comp_unit_entries("I1", ui, u):
        if e["kind"] != "fnullN":
            continue
        d = np.load(comp_real_path(u, e), allow_pickle=True)
        off.append(GEN._offender(np.asarray(d["dlnL_det"], float), d["lnK"], d["qmax"]))
    fl = GEN.floor_v22(np.array(off))
    led_real = [np.load(REAL_AS_LEDGER.format(tag=f"AS_e03_h1275_k5_s{s}"),
                        allow_pickle=True) for s in C_SKIES]
    rf = np.array([float(l["floor_adopted"]) for l in led_real])
    rzf = np.array([float(l["floor_zero_frac"]) for l in led_real])
    print(f"   iso floor: {fl['adopted']:.2f} +/- {fl['adopted_sd']:.2f} "
          f"[{fl['estimator']}, zf {fl['zero_frac']:.2f}]", flush=True)
    print(f"   real survivor per-sky floors: {rf.min():.1f}-{rf.max():.1f} "
          f"(zf {rzf.min():.2f}-{rzf.max():.2f}) -- level shift REPORTED, not gated",
          flush=True)
    ok = (fl["zero_frac"] <= GEN.ZF_GATE) and np.isfinite(fl["gumbel"]) \
        and np.isfinite(fl["gumbel_sd"]) and audit_ok
    print(f"   g0: {'PASS' if ok else 'FAIL'} (binding: zf <= 0.20, finite Gumbel, "
          f"uniformity audit)", flush=True)
    wall = float(np.median(times[3:])) if len(times) > 4 else np.nan
    del P
    return ok, wall


def mode_gate(args):
    os.makedirs(OUT, exist_ok=True)
    jax, jnp, C, B1, TE, F, IG = _import_stack()
    print(f"jax {jax.__version__}, {jax.devices()} host={socket.gethostname()}",
          flush=True)
    okA, t_real_off = gate_gA(C, F, jnp)
    ok0, t_real_iso = gate_g0(C, F, jnp)
    ok = okA and ok0
    # ---- header estimate (BEFORE any wide submit; measured warm walls) ----
    tr = t_real_iso if np.isfinite(t_real_iso) else 0.55
    # GEN banked warm walls for the heavier shapes (gen_warm.npz), scaled by this
    # lane's measured/banked ratio at (1, 30)
    scale = tr / 0.51
    t2, t3 = 0.66 * scale, 0.83 * scale
    build = 300.0
    est_I1 = (32 * 115 * tr + 16 * build) / 3600.0
    est_I2 = ((130 * tr) + (2 * 130 * t2) + (2 * 130 * t3) + 3 * build * 2) / 3600.0
    est_S1 = (8 * 115 * tr + 2 * (build + 300)) / 3600.0
    tot = est_I1 + est_I2 + est_S1
    print(f"\n=== HEADER ESTIMATE (GPU-hr of compute, warm wall {tr:.2f}s/real): "
          f"I1 {est_I1:.1f} + I2 {est_I2:.1f} + S1 {est_S1:.1f} = {tot:.1f} "
          f"vs STOP {STOP_GPU_HR} ===", flush=True)
    if tot > STOP_GPU_HR:
        print("OVER THE STOP BAR -- apply the trim ladder (S1, then I2 op-point pair, "
              "then I2); do NOT submit the full fan", flush=True)
    np.savez(f"{OUT}/comp_gates.npz", gA=okA, g0=ok0, t_real_iso=t_real_iso,
             t_real_off=t_real_off, est_I1=est_I1, est_I2=est_I2, est_S1=est_S1,
             est_total=tot, host=socket.gethostname())
    print(f"\n=== COMPASS gates: {'ALL PASS' if ok else 'FAIL -- STOP'} ===", flush=True)
    return 0 if ok else 1


# ============================================================
# analysis (CPU, from banks; iso side COMPASS_results, real side GENERALISE_results)
# ============================================================
def _unit_rows(arm):
    rows = []
    for ui, u in enumerate(units_for(arm)):
        p = comp_ledger_path(u)
        if not os.path.exists(p):
            print(f"  MISSING ledger: {u['tag']}", flush=True)
            continue
        z = np.load(p, allow_pickle=True)
        rows.append((u, z))
    return rows


def _real_as_rows():
    rows = []
    for ci, (h, e, k) in enumerate(I1_CELLS):
        ed = f"e{int(round(e * 10)):02d}"
        for sky in C_SKIES:
            tag = f"AS_{ed}_h{hkey(h)}_k{k}_s{sky}"
            z = np.load(REAL_AS_LEDGER.format(tag=tag), allow_pickle=True)
            rows.append((ci, sky, z))
    return rows


def mode_analysis(args):
    rng = np.random.default_rng(GEN.BOOT_SEED + 7)
    print("=== COMPASS analysis (attribution ratio, constants, floor law) ===",
          flush=True)
    iso = _unit_rows("I1")
    real = _real_as_rows()
    out = {}
    print("\n-- ARM 1: the attribution table (pre-registered: rho = sd_sky iso/real) --",
          flush=True)
    for ci, (h, e, k) in enumerate(I1_CELLS):
        ic = [z for (u, z) in iso if u["ci"] == ci]
        rc = [z for (c, s, z) in real if c == ci]
        if len(ic) < 8 or len(rc) < 8:
            print(f"  cell {ci} INCOMPLETE ({len(ic)}/8 iso, {len(rc)}/8 real)",
                  flush=True)
            continue
        cin = np.array([float(z["corr"]) for z in ic])
        crn = np.array([float(z["corr"]) for z in rc])
        cin_p = np.array([float(z["corr_pess"]) for z in ic])
        crn_p = np.array([float(z["corr_pess"]) for z in rc])
        bs_i = np.array([float(z["corr_bs_sd"]) for z in ic])
        bs_r = np.array([float(z["corr_bs_sd"]) for z in rc])
        fi = np.array([float(z["floor_adopted"]) for z in ic])
        fr = np.array([float(z["floor_adopted"]) for z in rc])
        si, sr = np.std(cin, ddof=1), np.std(crn, ddof=1)
        rho = si / sr
        bs = np.array([np.std(cin[rng.integers(0, 8, 8)], ddof=1)
                       / max(np.std(crn[rng.integers(0, 8, 8)], ddof=1), 1e-9)
                       for _ in range(GEN.NBOOT)])
        lo, hi = np.percentile(bs, [16, 84])
        v_i = max(si ** 2 - np.mean(bs_i ** 2), 0.0)
        v_r = max(sr ** 2 - np.mean(bs_r ** 2), 1e-12)
        rho_dec = np.sqrt(v_i / v_r)
        verdict = ("COLLAPSES" if rho < 0.5 else
                   "PERSISTS" if rho >= 0.75 else "MARGINAL (band)")
        conf_i = int(np.sum(cin_p > 1.0)); conf_r = int(np.sum(crn_p > 1.0))
        print(f"  cell (h={h}, e={e}, {k}+{16-k}): iso {np.mean(cin):.2f} +/- {si:.2f} "
              f"(range {cin.min():.2f}-{cin.max():.2f}, {conf_i}/8 CONF) | real "
              f"{np.mean(crn):.2f} +/- {sr:.2f} ({conf_r}/8 CONF) | rho = {rho:.2f} "
              f"[{lo:.2f}, {hi:.2f}] (decomposed {rho_dec:.2f}) -> {verdict}", flush=True)
        print(f"      floors iso {fi.min():.1f}-{fi.max():.1f} (x{fi.max()/max(fi.min(),1e-9):.1f})"
              f" | real {fr.min():.1f}-{fr.max():.1f} (x{fr.max()/max(fr.min(),1e-9):.1f})",
              flush=True)
        out[f"cell{ci}"] = dict(h=h, e=e, k=k, iso_counts=cin, real_counts=crn,
                                iso_counts_pess=cin_p, real_counts_pess=crn_p,
                                rho=rho, rho_lo=lo, rho_hi=hi, rho_dec=rho_dec,
                                iso_floors=fi, real_floors=fr, verdict=verdict)
    print("\n-- ARM 2a: the channel-budget contrast (pooled protocol) --", flush=True)
    RECUT = {(1, H_CENSUS, 3): (0.70, "below"), (2, H_CENSUS, 3): (2.77, "CONFIRMED"),
             (3, H_CENSUS, 3): (2.50, "CONFIRMED")}
    for (u, z) in _unit_rows("I2"):
        key = (u["n_ecc"], u["h"], u["k"])
        ref = RECUT.get(key)
        na = (16 - u["n_ecc"]) + 8 * u["n_ecc"]     # S7.6.4c channel-budget convention
        print(f"  {u['tag']} (channels {na}): floor {float(z['floor_adopted']):.2f} "
              f"+/- {float(z['floor_adopted_sd']):.2f} [{z['floor_estimator']}, "
              f"zf {float(z['floor_zero_frac']):.2f}] count {float(z['corr']):.2f} "
              f"[{float(z['corr_pess']):.2f}] -> {z['grade']}"
              + (f" | real re-cut {ref[0]} {ref[1]}" if ref else " | real: I1 survivor"),
              flush=True)
    print("\n-- ARM 2b: union-ceiling analogue (same-protocol both sides) --", flush=True)
    uni = {}
    for side, files in (("iso", sorted(glob.glob(f"{OUT}/comp_sig_I1_*.npz"))),
                        ("real", sorted(glob.glob(f"{GEN_OUT}/gen_sig_AS_*.npz")))):
        led = {u["tag"]: float(np.load(comp_ledger_path(u))["floor_adopted"])
               for u in i1_units() if os.path.exists(comp_ledger_path(u))} \
            if side == "iso" else \
              {f"AS_{t}": float(np.load(REAL_AS_LEDGER.format(tag=f'AS_{t}'))["floor_adopted"])
               for t in [f"e{int(round(e*10)):02d}_h{hkey(h)}_k{k}_s{s}"
                         for (h, e, k) in I1_CELLS for s in C_SKIES]}
        freq = {}
        nrl = 0
        for f in files:
            d = np.load(f, allow_pickle=True)
            tag = "_".join(os.path.basename(f).split("_")[2:-2])
            fl = led.get(tag)
            if fl is None:
                continue
            nrl += 1
            dl = np.nan_to_num(np.asarray(d["dlnL_det"], float), posinf=1e30)
            fin = dl < 1e29
            cert = fin & (dl > np.maximum(np.asarray(d["lnK"], float), fl)) \
                & (np.asarray(d["qmax"]) > GEN.QBAR) & np.asarray(d["on_true"], bool)
            for nm in np.asarray(d["names"])[cert]:
                freq[str(nm)] = freq.get(str(nm), 0) + 1
        top = sorted(freq.items(), key=lambda kv: -kv[1])[:8]
        print(f"  {side}: union = {len(freq)} pulsars over {nrl} signal reals; "
              f"leaders: " + ", ".join(f"{n} ({c})" for n, c in top), flush=True)
        uni[side] = freq
    print("\n-- ARM 2c: the floor-law column (mean per-sky floor vs h; real slope "
          "benchmark 1.66) --", flush=True)
    col = [(-13.0, "e05", 3), (-12.75, "e03", 3), (-12.5, "e03", 3)]
    for side in ("iso", "real"):
        hs, fm = [], []
        for (h, ed, k) in col:
            if side == "iso":
                fls = [float(z["floor_adopted"]) for (u, z) in iso
                       if (u["h"], u["edist"], u["k"]) == (h, ed, k)]
            else:
                fls = [float(z["floor_adopted"]) for (ci, s, z) in real
                       if I1_CELLS[ci] == (h, float(ed[1:]) / 10.0, k)]
            if len(fls) == 8:
                hs.append(h); fm.append(np.mean(fls))
                print(f"  {side} h={h} ({ed},{k}+{16-k}): floor mean {np.mean(fls):.1f} "
                      f"(min-max {np.min(fls):.1f}-{np.max(fls):.1f})", flush=True)
        if len(hs) >= 2:
            slope = np.polyfit(np.array(hs) * np.log(10.0), np.log(fm), 1)[0]
            print(f"  {side}: floor ~ h^{slope:.2f}", flush=True)
            out[f"floor_slope_{side}"] = slope
    print("\n-- S1: shuffled-assignment robustness (survivor cell) --", flush=True)
    s1 = _unit_rows("S1")
    if len(s1) == 8 and "cell1" in out:
        cs = np.array([float(z["corr"]) for (u, z) in s1])
        ci_ = out["cell1"]["iso_counts"]
        print(f"  shuf counts {np.mean(cs):.2f} +/- {np.std(cs, ddof=1):.2f} vs main "
              f"{np.mean(ci_):.2f} +/- {np.std(ci_, ddof=1):.2f} "
              f"(|d mean| = {abs(np.mean(cs)-np.mean(ci_)):.2f}, in sky-sd units "
              f"{abs(np.mean(cs)-np.mean(ci_))/max(np.std(ci_, ddof=1),1e-9):.2f})",
              flush=True)
        out["s1_counts"] = cs
    flat = {}
    for k_, v in out.items():
        if isinstance(v, dict):
            for k2, v2 in v.items():
                flat[f"{k_}_{k2}"] = v2
        else:
            flat[k_] = v
    np.savez(f"{OUT}/comp_analysis.npz", **flat,
             union_iso=np.array(sorted(uni.get("iso", {}).items())) if uni.get("iso") else np.zeros((0, 2)),
             union_real=np.array(sorted(uni.get("real", {}).items())) if uni.get("real") else np.zeros((0, 2)))
    print(f"\n  saved {OUT}/comp_analysis.npz", flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "run", "ledgers", "status", "analysis"])
    ap.add_argument("--arm", choices=["I1", "I2", "S1"], default=None)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    a = ap.parse_args()
    if a.mode == "run" and a.arm is None:
        sys.exit("run mode needs --arm")
    sys.exit({"gate": mode_gate, "run": mode_run, "ledgers": mode_ledgers,
              "status": mode_status, "analysis": mode_analysis}[a.mode](a))
