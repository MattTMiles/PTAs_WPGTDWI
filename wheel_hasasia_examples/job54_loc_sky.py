"""JOB 54 measure (b): expected (noiseless) sky likelihood, other parameters at truth.

ln L(n; set) = sum_a inc_a(n) + sum_{a in phased set} [coh_a(n) - inc_a(n)]      (job54_loc_core.parts_fn)
Grid: rotated lon/lat (u, v) about the truth, uniform step h, cell solid angle h^2 cos v.
 1. whole sky at 0.5 deg, n_ph = 0  -> 99.9% set -> bounding box
 2. 512 x 512 over the box, n_ph = 0 -> 99.9% set, dilated by one cell = mask
 3. fine grid (step h) on masked cells, all cases at once (nested sets = cumulative sums)
 4. zoom: the 12 x 12 base cells around the truth are replaced by cells 64x finer
 5. 90% area from area/mass histograms in ln L (bin 0.02); islands by 8-connected labelling
"""
import numpy as np, jax, jax.numpy as jnp, os, sys, time, json
from scipy import ndimage
import job54_loc_core as C

HI, BW, NB = 3.0, 0.02, 3500
CH = 65536
ZR, ZW = 64, 6
ARR = C.build_all(); SRC = C.draw_sources()
ON, OR_ = jnp.asarray(ARR["order_near"]), jnp.asarray(ARR["order_rand"])
NPH = jnp.asarray(C.NPH_LIST)
SIGLS = (0.1, 0.3, 1.0, 3.0, 10.0); NPART = (4, 10, 116)
SKYDIR = os.path.join(C.DATA, "sky"); os.makedirs(SKYDIR, exist_ok=True)
S8 = np.ones((3, 3), int)


def cases_main(u, v, par):
    """rows 0..10: nearest-first n_ph list; rows 11..21: random order."""
    coh, inc, *_ = C.parts_fn(u, v, par)
    d = coh - inc; t0 = inc.sum(1)
    outs = []
    for order in (ON, OR_):
        cs = jnp.concatenate([jnp.zeros((d.shape[0], 1)), jnp.cumsum(d[:, order], 1)], 1)
        outs.append(t0[None, :] + cs[:, NPH].T)
    return jnp.concatenate(outs, 0)


def cases_partial(u, v, par, M=24):
    """row 0: n_ph=0; then for sigma_L in SIGLS, n_ph in NPART (nearest-first): distance marginalised
    over its Gaussian prior (phase-only; Jacobi-Anger series, Miller recurrence for I_m)."""
    coh, inc, a, alpha, sphi1, base = C.parts_fn(u, v, par)
    t0 = inc.sum(1)
    a_s = jnp.clip(a, 1e-2, None)
    Ip1 = jnp.zeros_like(a_s); Ic = jnp.full_like(a_s, 1e-250)
    norm = jnp.zeros_like(a_s); accs = [jnp.zeros_like(a_s) for _ in SIGLS]
    for m in range(M, 0, -1):
        norm = norm + 2 * Ic
        cm = 2 * Ic * jnp.cos(m * alpha)
        accs = [acc + cm * jnp.exp(-0.5 * (m * sl * sphi1) ** 2) for acc, sl in zip(accs, SIGLS)]
        Ic, Ip1 = (2.0 * m / a_s) * Ic + Ip1, Ic
    norm = norm + Ic
    rows = [t0]
    for acc in accs:
        pl = base + a_s + jnp.log(jnp.clip((acc + Ic) / norm, 1e-300, None))
        cs = jnp.cumsum((pl - inc)[:, ON], 1)
        for n in NPART:
            rows.append(t0 + cs[:, n - 1])
    return jnp.stack(rows, 0)


def make_chunk_fn(case_fn):
    @jax.jit
    def f(u, v, w, par, lt):
        rel = case_fn(u, v, par) - lt[:, None]
        b = jnp.clip(jnp.floor((HI - rel) / BW), 0, NB - 1).astype(jnp.int32)
        m = w[None, :] * jnp.exp(jnp.clip(rel, -700.0, 50.0))
        ha = jax.vmap(lambda bb: jax.ops.segment_sum(w, bb, NB))(b)
        hm = jax.vmap(lambda bb, mm: jax.ops.segment_sum(mm, bb, NB))(b, m)
        return jnp.clip(rel, -6e4, 6e4).astype(jnp.float16), ha, hm
    return f


CHUNK = {"main": make_chunk_fn(cases_main), "partial": make_chunk_fn(cases_partial)}
KCASE = {"main": 22, "partial": 1 + len(SIGLS) * len(NPART)}


def eval_points(kind, u, v, w, par, lt, keep=True):
    n = len(u); K = KCASE[kind]
    R = np.empty((K, n), np.float16) if keep else None
    ha = np.zeros((K, NB)); hm = np.zeros((K, NB))
    for i0 in range(0, n, CH):
        uu, vv, ww = u[i0:i0 + CH], v[i0:i0 + CH], w[i0:i0 + CH]
        m = len(uu)
        if m < CH:
            uu = np.pad(uu, (0, CH - m)); vv = np.pad(vv, (0, CH - m)); ww = np.pad(ww, (0, CH - m))
        r, a, b = CHUNK[kind](jnp.asarray(uu), jnp.asarray(vv), jnp.asarray(ww), par, lt)
        if keep:
            R[:, i0:i0 + m] = np.asarray(r)[:, :m]
        ha += np.asarray(a); hm += np.asarray(b)
    return R, ha, hm


def level(ha, hm, frac):
    cm = np.cumsum(hm); target = frac * cm[-1]
    j = int(np.searchsorted(cm, target)); prev = cm[j - 1] if j > 0 else 0.0
    fr = (target - prev) / max(hm[j], 1e-300)
    return ha[:j].sum() + fr * ha[j], HI - BW * (j + fr)


def top_set(rel, w, frac):
    o = np.argsort(-rel.ravel()); m = (np.exp(rel.ravel()[o] - rel.max()) * w.ravel()[o]).cumsum()
    k = int(np.searchsorted(m, frac * m[-1])) + 1
    msk = np.zeros(rel.size, bool); msk[o[:k]] = True
    return msk.reshape(rel.shape)


def run_pass(isrc, f, snr, kind="main", n_budget=2e7, mc=1e9, h_scale=1.0, save_maps=False, verbose=True):
    t0 = time.time()
    s = C.src_at(SRC, isrc)
    h_amp = snr * C.snr_unit_h(ARR, s, f, mc)
    par = C.sky_par(ARR, s, f, h_amp, mc)
    K = KCASE[kind]
    lt = jnp.asarray(CHUNK[kind](jnp.zeros(CH), jnp.zeros(CH), jnp.zeros(CH), par, jnp.zeros(K))[0][:, 0].astype(jnp.float64))
    # exact truth totals (float64) by a tiny direct call
    lt = (cases_main if kind == "main" else cases_partial)(jnp.zeros(1), jnp.zeros(1), par)[:, 0]
    # --- 1. whole sky, 0.5 deg
    nu, nv = 720, 360
    uc = (np.arange(nu) + 0.5) * (2 * np.pi / nu) - np.pi; vc = (np.arange(nv) + 0.5) * (np.pi / nv) - np.pi / 2
    U, V = np.meshgrid(uc, vc); W = np.cos(V) * (2 * np.pi / nu) * (np.pi / nv)
    R, _, _ = eval_points(kind, U.ravel(), V.ravel(), W.ravel(), par, lt)
    r0 = R[0].astype(float).reshape(nv, nu)
    m1 = ndimage.binary_dilation(top_set(r0, W, 0.999), S8)
    jv, ju = np.where(m1)
    ub = (uc[ju.min()] - np.pi / nu, uc[ju.max()] + np.pi / nu); vb = (vc[jv.min()] - 0.5 * np.pi / nv, vc[jv.max()] + 0.5 * np.pi / nv)
    sky_frac_coarse = float(W[m1].sum() / (4 * np.pi))
    # --- 2. 512 x 512 over the box
    n2 = 512
    du2, dv2 = (ub[1] - ub[0]) / n2, (vb[1] - vb[0]) / n2
    u2 = ub[0] + (np.arange(n2) + 0.5) * du2; v2 = vb[0] + (np.arange(n2) + 0.5) * dv2
    U, V = np.meshgrid(u2, v2); W = np.cos(V) * du2 * dv2
    R, _, _ = eval_points(kind, U.ravel(), V.ravel(), W.ravel(), par, lt)
    r0 = R[0].astype(float).reshape(n2, n2)
    m2 = ndimage.binary_dilation(top_set(r0, W, 0.999), S8)
    jv, ju = np.where(m2)
    ulo, uhi = ub[0] + ju.min() * du2, ub[0] + (ju.max() + 1) * du2
    vlo, vhi = vb[0] + jv.min() * dv2, vb[0] + (jv.max() + 1) * dv2
    a_mask = float(W[m2].sum()); a_box = (uhi - ulo) * (vhi - vlo)
    # --- 3. fine grid, truth at a cell centre
    h = max(np.sqrt(a_mask / n_budget), np.sqrt(a_box / 6e7)) * h_scale
    nxL, nxR = int(np.ceil(-ulo / h)), int(np.ceil(uhi / h)); nyL, nyR = int(np.ceil(-vlo / h)), int(np.ceil(vhi / h))
    nx, ny = nxL + nxR + 1, nyL + nyR + 1
    ucf = (np.arange(nx) - nxL) * h; vcf = (np.arange(ny) - nyL) * h
    iu = np.clip(np.floor((ucf - ub[0]) / du2).astype(int), 0, n2 - 1); iv = np.clip(np.floor((vcf - vb[0]) / dv2).astype(int), 0, n2 - 1)
    mf = m2[np.ix_(iv, iu)]
    wy0, wx0 = nyL - ZW, nxL - ZW
    mf[wy0:wy0 + 2 * ZW, wx0:wx0 + 2 * ZW] = True
    idx = np.flatnonzero(mf); del mf
    iy, ix = idx // nx, idx % nx
    wcell = h * h * np.cos(vcf[iy])
    R, ha, hm = eval_points(kind, ucf[ix], vcf[iy], wcell, par, lt)
    t_fine = time.time() - t0
    # window base cells (to subtract) and zoom cells
    inwin = (iy >= wy0) & (iy < wy0 + 2 * ZW) & (ix >= wx0) & (ix < wx0 + 2 * ZW)
    _, haw, hmw = eval_points(kind, ucf[ix[inwin]], vcf[iy[inwin]], wcell[inwin], par, lt, keep=False)
    nz = 2 * ZW * ZR; hz = h / ZR
    uz = (wx0 - nxL - 0.5) * h + (np.arange(nz) + 0.5) * hz; vz = (wy0 - nyL - 0.5) * h + (np.arange(nz) + 0.5) * hz
    Uz, Vz = np.meshgrid(uz, vz); Wz = hz * hz * np.cos(Vz)
    Rz, haz, hmz = eval_points(kind, Uz.ravel(), Vz.ravel(), Wz.ravel(), par, lt)
    Ha = np.clip(ha - haw, 0, None) + haz; Hm = np.clip(hm - hmw, 0, None) + hmz
    # fringe periods at truth [rad]: lambda_p / (L sin mu)
    n0 = np.asarray(par["n0"]); cmu = ARR["p_hats"] @ n0
    lam_pc = 299792458.0 / f / (C.p2.KPC_M / 1000.0)
    period = lam_pc / (ARR["L_pc"] * np.sqrt(np.clip(1 - cmu ** 2, 1e-12, None)))
    out = dict(isrc=isrc, f=f, snr=snr, h_rad=h, h_deg=np.degrees(h), nx=nx, ny=ny, n_cells=len(idx), a_mask_deg2=a_mask * C.DEG2,
               a_box_deg2=a_box * C.DEG2, sky_frac_coarse=sky_frac_coarse, box_deg=np.degrees([ulo, uhi, vlo, vhi]),
               period_rad=period, max_rel=np.maximum(R.max(1), Rz.reshape(K, -1).max(1)).astype(float))
    res = {k: np.full(K, np.nan) for k in ("A90", "A50", "A99", "thr90", "n_isl", "A_truth_isl", "truth_in", "zoom_frac",
                                           "n_members", "n_base_members_out", "pmin_over_h")}
    zt = ZW * ZR + ZR // 2
    base_area_map_w = wcell
    for k in range(K):
        A90, thr = level(Ha[k], Hm[k], 0.9)
        res["A90"][k] = A90 * C.DEG2; res["thr90"][k] = thr
        res["A50"][k] = level(Ha[k], Hm[k], 0.5)[0] * C.DEG2; res["A99"][k] = level(Ha[k], Hm[k], 0.99)[0] * C.DEG2
        res["zoom_frac"][k] = hmz[k].sum() / Hm[k].sum()
        Zm = Rz[k].reshape(nz, nz) > thr
        memb = (R[k] > thr) & ~inwin
        nout = int(memb.sum()); res["n_base_members_out"][k] = nout; res["n_members"][k] = nout + int(Zm.sum())
        res["truth_in"][k] = float(Zm[zt, zt])
        if nout == 0:
            lab, nl = ndimage.label(Zm, S8); res["n_isl"][k] = nl
            res["A_truth_isl"][k] = float(Wz[lab == lab[zt, zt]].sum() * C.DEG2) if Zm[zt, zt] else 0.0
        else:
            y0, y1, x0, x1 = min(iy[memb].min(), wy0), max(iy[memb].max(), wy0 + 2 * ZW - 1), min(ix[memb].min(), wx0), max(ix[memb].max(), wx0 + 2 * ZW - 1)
            am = np.zeros((y1 - y0 + 1, x1 - x0 + 1))
            am[iy[memb] - y0, ix[memb] - x0] = wcell[memb]
            am[wy0 - y0:wy0 - y0 + 2 * ZW, wx0 - x0:wx0 - x0 + 2 * ZW] = (Zm * Wz).reshape(2 * ZW, ZR, 2 * ZW, ZR).sum((1, 3))
            lab, nl = ndimage.label(am > 0, S8); res["n_isl"][k] = nl
            lt_ = lab[nyL - y0, nxL - x0]
            res["A_truth_isl"][k] = float(am[lab == lt_].sum() * C.DEG2) if lt_ > 0 else 0.0
    # resolution: finest fringe of each phased set / h
    if kind == "main":
        pm = []
        for order in (ARR["order_near"], ARR["order_rand"]):
            for n in C.NPH_LIST:
                pm.append(np.inf if n == 0 else period[order[:n]].min() / h)
        res["pmin_over_h"] = np.array(pm)
    else:
        pm = [np.inf] + [period[ARR["order_near"][:n]].min() / h for _ in SIGLS for n in NPART]
        res["pmin_over_h"] = np.array(pm)
    out.update(res); out["t_fine_s"] = t_fine
    if save_maps:
        ND = 256
        dj = (np.arange(nx) * ND) // nx; di = (np.arange(ny) * ND) // ny
        disp = np.full((11, ND, ND), -6e4, np.float32)
        for k in range(11):
            np.maximum.at(disp[k], (di[iy], dj[ix]), R[k].astype(np.float32))
            zmax = Rz[k].reshape(2 * ZW, ZR, 2 * ZW, ZR).max((1, 3)).astype(np.float32)
            wy = np.arange(wy0, wy0 + 2 * ZW); wx = np.arange(wx0, wx0 + 2 * ZW)
            np.maximum.at(disp[k], (di[wy][:, None].repeat(2 * ZW, 1), dj[wx][None, :].repeat(2 * ZW, 0)), zmax)
        ue = np.array([(ucf[dj == j].mean() if (dj == j).any() else np.nan) for j in range(ND)])
        ve = np.array([(vcf[di == i].mean() if (di == i).any() else np.nan) for i in range(ND)])
        out["disp"] = disp.astype(np.float16); out["disp_u_deg"] = np.degrees(ue); out["disp_v_deg"] = np.degrees(ve)
        # zoomed set: window 30x smaller (linear), 256 x 256, block max of an 8x finer evaluation
        sub = 8; nzm = ND * sub
        hu, hv = nx * h / 30.0, ny * h / 30.0
        uzz = -hu / 2 + (np.arange(nzm) + 0.5) * hu / nzm; vzz = -hv / 2 + (np.arange(nzm) + 0.5) * hv / nzm
        Uq, Vq = np.meshgrid(uzz, vzz)
        Rq, _, _ = eval_points(kind, Uq.ravel(), Vq.ravel(), np.ones(Uq.size), par, lt)
        zm = Rq[:11].reshape(11, ND, sub, ND, sub).max((2, 4))
        zm[:, ND // 2 - 1:ND // 2 + 1, ND // 2 - 1:ND // 2 + 1] = np.maximum(zm[:, ND // 2 - 1:ND // 2 + 1, ND // 2 - 1:ND // 2 + 1], 0)
        out["zoom"] = zm.astype(np.float16)
        out["zoom_u_deg"] = np.degrees(uzz.reshape(ND, sub).mean(1)); out["zoom_v_deg"] = np.degrees(vzz.reshape(ND, sub).mean(1))
        out["zoom_fine_step_deg"] = np.degrees(hu / nzm)
        for nm, (ua, va) in (("disp", (ue, ve)), ("zoom", (np.radians(out["zoom_u_deg"]), np.radians(out["zoom_v_deg"])))):
            Ug, Vg = np.meshgrid(np.nan_to_num(ua), np.nan_to_num(va))
            nvec = np.asarray(C.uv2n(jnp.asarray(Ug.ravel()), jnp.asarray(Vg.ravel()), par))
            out[nm + "_ra_deg"] = (np.degrees(np.arctan2(nvec[:, 1], nvec[:, 0])) % 360).reshape(ND, ND).astype(np.float32)
            out[nm + "_dec_deg"] = np.degrees(np.arcsin(np.clip(nvec[:, 2], -1, 1))).reshape(ND, ND).astype(np.float32)
        out["true_ra_deg"] = np.degrees(s["phi"]) % 360; out["true_dec_deg"] = 90 - np.degrees(s["theta"])
    out["t_total_s"] = time.time() - t0
    if verbose:
        print(f"src {isrc} f {f*1e9:.0f} snr {snr:.0f} {kind}: h {np.degrees(h):.2e} deg, cells {len(idx):.2e}, mask {a_mask*C.DEG2:.0f} deg2, "
              f"t {out['t_total_s']:.0f}s (fine {t_fine:.0f}s) maxrel {out['max_rel'].max():.2f}", flush=True)
        print("   A90 near:", " ".join(f"{x:.3g}" for x in out["A90"][:11 if kind == 'main' else K]), flush=True)
        print("   isl near:", " ".join(f"{x:.0f}" for x in out["n_isl"][:11 if kind == 'main' else K]), flush=True)
        if kind == "main":
            print("   A90 rand:", " ".join(f"{x:.3g}" for x in out["A90"][11:]), flush=True)
            print("   zoomfrac near:", " ".join(f"{x:.2f}" for x in out["zoom_frac"][:11]), " pmin/h:", " ".join(f"{x:.1f}" for x in out["pmin_over_h"][:11]), flush=True)
    return out


def fname(isrc, f, snr, kind, tag=""):
    return os.path.join(SKYDIR, f"{kind}_s{isrc:03d}_f{f*1e9:.0f}_snr{snr:.0f}{tag}.npz")


if __name__ == "__main__":
    # usage: job54_loc_sky.py kind f_nHz snr src_list budget [tag] [h_scale] [maps]
    kind = sys.argv[1]; f = float(sys.argv[2]) * 1e-9; snr = float(sys.argv[3])
    srcs = [int(x) for x in sys.argv[4].split(",")]; budget = float(sys.argv[5])
    tag = sys.argv[6] if len(sys.argv) > 6 else ""; hs = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    maps = len(sys.argv) > 8 and sys.argv[8] == "maps"
    for i in srcs:
        fn = fname(i, f, snr, kind, tag)
        if os.path.exists(fn):
            continue
        out = run_pass(i, f, snr, kind, budget, h_scale=hs, save_maps=maps)
        np.savez(fn, **out)
