# ================= JOB 53 Part 2 -- the heatmap ==================================================
# Controlled experiment (no skies).  N_s sources, each alone in its own bin with share beta:
#     z_k = e^{i kappa_k L_true} (1 + C_k),   C_k complex Gaussian, E|C_k|^2 = (1-beta)/beta
# Decoder: the joint fit.  With one unit-amplitude source per bin the joint-fit statistic
#     S(L) = sum_k [2 Re(z_k conj t_k(L)) - |t_k(L)|^2]/P - L^2/(2 sigma_d^2),  t_k = e^{i kappa_k L}
# has |t_k| = 1, so it equals the matched filter of Job 49 (D-mf) exactly; here it is evaluated on a
# fine grid (step <= 0.05 pc and <= one sixth of the shortest tooth spacing).
# Frequencies and angles as Job 47c: the N_s highest-ranked combs (by nats) of the machine's booked
# list on population t mod 50 (constructed populations; only the FREQUENCIES are used), angle
# mu = arccos U(-0.9, 0.9) per source per realisation, kappa = 2 pi f (1 - cos mu)/c.
# Modes: repro | grid | rows | overlay <beta...> | plot
import os, sys, json, time, csv
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np
import multiprocessing as mp
BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
sys.path.insert(0, BASE); os.chdir(BASE); os.makedirs("job53_data", exist_ok=True)
TW = 0.25
XS = np.geomspace(0.1, 300.0, 20)
_lg = np.linspace(np.log(0.02/0.98), np.log(0.98/0.02), 20)
BS = 1.0/(1.0+np.exp(-_lg))                     # 20 values of beta, equal steps in ln(beta/(1-beta))
NSL = (2, 4, 8, 12)
NREAL = 500
mode = sys.argv[1] if len(sys.argv) > 1 else "grid"
NW = int(os.environ.get("JOB53_NW", "4"))


def cross_x(row, xs, level=0.9):
    for a in range(len(xs)-1):
        if row[a] >= level >= row[a+1]:
            f = (row[a]-level)/max(row[a]-row[a+1], 1e-9)
            return float(10**(np.log10(xs[a])+f*(np.log10(xs[a+1])-np.log10(xs[a]))))
    return np.nan


if mode != "plot":
    # Job 49 D1 preamble (loader, comb sets, contributor shapes, purity_cell_sw), exec'd from the file
    src = open("job49_D1.py").read(); cut = src.index("XS = np.geomspace(0.1, 300.0, 20)")
    NS49 = {"__name__": "job49_D1_preamble"}
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(src[:cut], "<job49_D1 preamble>", "exec"), NS49)
    SETS = NS49["SETS"]; SHAPE = NS49["SHAPE"]; C_PC_S = NS49["C_PC_S"]; purity_cell_sw = NS49["purity_cell_sw"]
    nset = [0 if s is None else len(s) for s in SETS]


def fset(t, Ns):
    jj = t % 50
    while nset[jj] < Ns: jj = (jj+1) % 50
    return SETS[jj][:Ns]


def decode_fine(z, kap, sigd, P, Lt):
    dL = min(0.05, (2*np.pi/float(kap.max()))/6.0)
    M = int(np.floor(6*sigd/dL))+1; Lg = -3*sigd+dL*np.arange(M)
    S = np.zeros(M, np.float32)
    for c0 in range(0, M, 8192):
        L = Lg[c0:c0+8192]
        S[c0:c0+8192] = (2.0/P)*np.real(z[:, None]*np.exp(-1j*kap[:, None]*L[None, :])).sum(0)-L**2/(2*sigd**2)
    return abs(float(Lg[int(np.argmax(S))])-Lt) <= TW


def cell(args):
    """fraction decoded in one cell: N_s equal-share sources, share beta, prior width sigd"""
    Ns, sigd, beta, nreal, tag = args
    P = (1-beta)/beta; ok = 0
    for t in range(nreal):
        fk = fset(t, Ns)
        rg = np.random.default_rng([53, 2, Ns, int(round(1e6*sigd)), int(round(1e6*beta)), t, tag])
        mu = np.arccos(rg.uniform(-0.9, 0.9, Ns)); kap = 2*np.pi*fk*(1-np.cos(mu))/C_PC_S
        Lt = float(np.clip(rg.normal(0, sigd), -2.9*sigd, 2.9*sigd))
        C = (rg.normal(0, 1, Ns)+1j*rg.normal(0, 1, Ns))*np.sqrt(P/2)
        z = np.exp(1j*kap*Lt)*(1+C)
        ok += int(decode_fine(z, kap, sigd, P, Lt))
    return ok/nreal


def cell49(args):
    """Job 49 D1 physical D-mf cell (same streams, same data), decoded on its 601-point grid AND on the fine grid"""
    sigd, sigphi, nreal = args; P = sigphi**2; ok = np.zeros(2, int); NCOMB = 4
    for t in range(nreal):
        jj = t % 50; fk = SETS[jj]; shp = SHAPE[jj]
        if fk is None: fk = SETS[(jj+1) % 50]; shp = SHAPE[(jj+1) % 50]
        fk = fk[:NCOMB]
        rg = np.random.default_rng(470_000_000+t*7919+int(1e6*sigd) % 99991+int(1e6*sigphi) % 9973)
        mu = np.arccos(rg.uniform(-0.9, 0.9, len(fk))); kap = 2*np.pi*fk*(1-np.cos(mu))/C_PC_S
        Lt = float(np.clip(rg.normal(0, sigd), -2.9*sigd, 2.9*sigd))
        Lg = np.linspace(-3*sigd, 3*sigd, 601)
        rph = np.random.default_rng([49, 47, t, int(1e6*sigd), int(1e6*sigphi)])
        u = np.ones(len(fk), complex)
        for c in range(len(fk)):
            top, rem = shp[c]
            th = rph.uniform(0, 2*np.pi, len(top)); u[c] += np.sum(np.sqrt(P*top)*np.exp(1j*th))
            g = rph.normal(0, 1, 2)*np.sqrt(P*rem/2); u[c] += g[0]+1j*g[1]
        z = np.exp(1j*kap*Lt)*u
        S = (2.0/P)*np.real(z[:, None]*np.exp(-1j*kap[:, None]*Lg[None, :])).sum(0)-Lg**2/(2*sigd**2)
        ok[0] += int(abs(float(Lg[int(np.argmax(S))])-Lt) <= TW)
        ok[1] += int(decode_fine(z, kap, sigd, P, Lt))
    return ok/nreal


if __name__ == "__main__" and mode != "plot":
    t0 = time.time(); pool = mp.get_context("fork").Pool(NW)
    print(f"[{mode}] comb sets: {sum(n > 0 for n in nset)}/50 populations; sizes min {min(n for n in nset if n)}, max {max(nset)}; with >= 12: {sum(n >= 12 for n in nset)}", flush=True)
    if mode == "repro":
        YS = np.linspace(0.05, 3.0, 20)
        ref = np.load("job49_data/D1_grid.npz")
        R = np.array(pool.map(cell49, [(x, sp, 250) for sp in YS for x in XS])).reshape(20, 20, 2)
        d = float(np.max(np.abs(R[:, :, 0]-ref["G250"][2])))
        print(f"REPRO Job 49 D-mf panel (job49_data/D1_grid.npz, 250 realisations, 601-point grid): max |difference| over 400 cells = {d:.2e}")
        assert d == 0.0
        D = R[:, :, 1]-R[:, :, 0]; BE = 1/(1+YS**2)
        print(f"fine grid minus coarse grid, same data: mean {D.mean():+.4f}; cells changed {int((D != 0).sum())}/400; largest gain {D.max():+.3f}, largest loss {D.min():+.3f}")
        for lo, hi in ((0, 0.6), (0.6, 3), (3, 20), (20, 400)):
            m = (XS >= lo) & (XS < hi); print(f"   sigma_d in [{lo}, {hi}) pc: mean change {D[:, m].mean():+.4f}, max {D[:, m].max():+.3f}; coarse-grid step there {6*XS[m].min()/600:.4f} to {6*XS[m].max()/600:.3f} pc")
        print("   mean change by column (sigma_d pc: change): " + "  ".join(f"{x:.3g}: {D[:, ix].mean():+.3f}" for ix, x in enumerate(XS)))
        for b in (0.8, 0.5, 0.3, 0.2):
            iy = int(np.argmin(np.abs(BE-b))); print(f"   row beta = {BE[iy]:.3f}: 0.9 crossing coarse {cross_x(R[iy, :, 0], XS):.2f} pc, fine {cross_x(R[iy, :, 1], XS):.2f} pc")
        np.savez("job53_data/heat_repro.npz", XS=XS, YS=YS, coarse=R[:, :, 0], fine=R[:, :, 1])
    elif mode == "grid":
        G = {}
        for Ns in NSL:
            G[Ns] = np.array(pool.map(cell, [(Ns, x, b, NREAL, 0) for b in BS for x in XS], chunksize=4)).reshape(20, 20)
            np.savez("job53_data/heat_grid.npz", XS=XS, BS=BS, **{f"G{k}": v for k, v in G.items()})
            print(f"  N_s = {Ns}: grid done, {time.time()-t0:.0f} s", flush=True)
        # tooth counts: K = 6 sigma_d kappa / 2 pi, per pc of sigma_d
        rg = np.random.default_rng(5302); KK = {}
        for Ns in NSL:
            lo = []; hi = []
            for t in range(2000):
                fk = fset(t, Ns); kap = 2*np.pi*fk*(1-np.cos(np.arccos(rg.uniform(-0.9, 0.9, Ns))))/C_PC_S
                lo.append(6*kap.min()/(2*np.pi)); hi.append(6*kap.max()/(2*np.pi))
            KK[Ns] = (float(np.median(lo)), float(np.median(hi)))
            print(f"  N_s = {Ns}: K per pc of sigma_d, median over realisations: coarsest tooth {KK[Ns][0]:.2f}, finest tooth {KK[Ns][1]:.2f}; frequencies used 10/50/90%: {np.percentile(np.concatenate([fset(t, Ns) for t in range(50)])*1e9, [10, 50, 90]).round(1)} nHz")
        json.dump({str(k): v for k, v in KK.items()}, open("job53_data/heat_K.json", "w"))
    elif mode == "rows":
        XR = np.geomspace(0.1, 300.0, 40); out = {}
        for Ns in NSL:
            for b in (0.95, 0.8, 0.5, 0.3, 0.2, 0.1):
                r = np.array(pool.map(cell, [(Ns, x, b, 1000, 1) for x in XR], chunksize=2))
                out[f"{Ns}|{b}"] = dict(frac=r.tolist(), cross90=cross_x(r, XR, 0.9), cross50=cross_x(r, XR, 0.5), first=float(r[0]), last=float(r[-1]))
                print(f"  N_s = {Ns:2d}, beta = {b:.2f}: 0.9 crossing at sigma_d = {out[f'{Ns}|{b}']['cross90']:.2f} pc; 0.5 crossing {out[f'{Ns}|{b}']['cross50']:.2f} pc; fraction at 0.1 pc {r[0]:.3f}, at 300 pc {r[-1]:.3f}  ({time.time()-t0:.0f} s)", flush=True)
                json.dump(dict(XR=XR.tolist(), rows=out), open("job53_data/heat_rows.json", "w"), indent=1)
    elif mode == "overlay":
        betas = [float(x) for x in sys.argv[2:]]; X7 = (0.3, 1, 3, 10, 30, 100, 300); out = {}
        for Ns in (4, 12):
            for b in betas:
                r = np.array(pool.map(cell, [(Ns, float(x), b, 4000, 2) for x in X7], chunksize=1))
                out[f"{Ns}|{b}"] = r.tolist()
                print(f"  N_s = {Ns:2d}, beta = {b:.4f}: " + "  ".join(f"{x} pc: {v:.3f}" for x, v in zip(X7, r)) + f"   ({time.time()-t0:.0f} s)", flush=True)
        json.dump(dict(X=X7, betas=betas, rows=out), open("job53_data/heat_overlay.json", "w"), indent=1)
    pool.close(); pool.join()
    print(f"TOTAL {time.time()-t0:.0f} s")

if __name__ == "__main__" and mode == "plot":
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.ndimage import uniform_filter
    from scipy.special import erf
    Z = np.load("job53_data/heat_grid.npz"); KK = json.load(open("job53_data/heat_K.json")); MJ = json.load(open("job53_data/main.json"))
    b50, b90 = MJ["bun"][1], MJ["bun"][2]
    lg = np.log(BS/(1-BS)); dl = lg[1]-lg[0]; ye = np.concatenate([[lg[0]-dl/2], lg+dl/2])
    lx = np.log10(XS); dx = lx[1]-lx[0]; xe = 10**np.concatenate([[lx[0]-dx/2], lx+dx/2])
    yt = np.array([0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98]); lt = np.log(yt/(1-yt))
    # prior alone gives 90%: erf(0.25/(sqrt2 s)) = 0.9
    from scipy.optimize import brentq
    s_pr = brentq(lambda s: erf(0.25/(np.sqrt(2)*s))-0.9, 0.01, 10)

    def panel(ax, Ns, overlay, big=False):
        G = Z[f"G{Ns}"]; SM = uniform_filter(G, size=3, mode="nearest")
        pc = ax.pcolormesh(xe, ye, G, cmap="Blues", vmin=0, vmax=1, shading="flat")
        cs = ax.contour(XS, lg, SM, levels=[0.5, 0.9], colors=["0.35", "k"], linewidths=[1.0, 1.4], linestyles=["--", "-"])
        ax.clabel(cs, fmt={0.5: "50%", 0.9: "90%"}, fontsize=7)
        ax.set_xscale("log"); ax.set_yticks(lt); ax.set_yticklabels([f"{v:g}" for v in yt]); ax.set_ylim(ye[0], ye[-1]); ax.set_xlim(xe[0], xe[-1])
        ax.axvline(s_pr, color="#2ca02c", lw=1.2, ls=":")
        ax.text(s_pr*0.93, lt[1], "prior alone ≥ 90%", rotation=90, fontsize=6.5, color="#2ca02c", ha="right", va="bottom")
        if overlay:
            for b, lab in ((b50, f"median known pulsar term, natural skies (β = {b50:.3f})"), (b90, f"90th percentile (β = {b90:.3f})")):
                ax.axhline(np.log(b/(1-b)), color="#d62728", lw=1.1, ls="-" if b == b90 else "--")
                ax.text(xe[-1]*0.9, np.log(b/(1-b))+0.08, lab, fontsize=6.5 if not big else 8, color="#d62728", ha="right", va="bottom")
            for x, lab in ((31, "census 10th pct\n31 pc"), (242, "census median\n242 pc")):
                ax.plot([x], [ye[0]], marker="^", color="#ff7f0e", ms=9, clip_on=False, zorder=5)
                ax.text(x, ye[0]+0.15, lab, fontsize=6.5 if not big else 8, color="#ff7f0e", ha="center", va="bottom")
        kc, kf = KK[str(Ns)]
        sx = ax.secondary_xaxis("top", functions=(lambda x: x*kc, lambda k: k/kc))
        sx.set_xlabel(f"K, coarsest tooth  (finest tooth: K × {kf/kc:.1f})", fontsize=8 if not big else 10)
        ax.set_title(f"$N_s$ = {Ns} sources", fontsize=10 if not big else 12, pad=28)
        return pc

    fig, axs = plt.subplots(2, 2, figsize=(11.5, 9.6))
    for ax, Ns in zip(axs.ravel(), NSL):
        pc = panel(ax, Ns, Ns in (4, 12))
    for ax in axs[1]: ax.set_xlabel(r"prior width $\sigma_d$ [pc]")
    for ax in axs[:, 0]: ax.set_ylabel(r"$\beta$ = each source's share of its bin")
    fig.tight_layout(rect=(0, 0, 0.92, 1)); cax = fig.add_axes([0.935, 0.15, 0.015, 0.7]); fig.colorbar(pc, cax=cax, label="fraction decoded correctly")
    fig.savefig("job53_fig1_heatmaps.png", dpi=160); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.4, 5.2)); pc = panel(ax, 4, True, big=True)
    ax.set_xlabel(r"prior width $\sigma_d$ [pc]"); ax.set_ylabel(r"$\beta$ = each source's share of its bin"); fig.colorbar(pc, ax=ax, label="fraction decoded correctly")
    fig.tight_layout(); fig.savefig("job53_fig2_heatmap_Ns4.png", dpi=200); plt.close(fig)
    print(f"wrote job53_fig1_heatmaps.png, job53_fig2_heatmap_Ns4.png; overlays beta {b50:.4f}, {b90:.4f}; prior-alone-90% boundary {s_pr:.3f} pc; beta grid: {np.round(BS, 3).tolist()}")
    for Ns in NSL:
        G = Z[f"G{Ns}"]; SM = uniform_filter(G, size=3, mode="nearest")
        print(f"  N_s = {Ns}: grid max {G.max():.3f} min {G.min():.3f}; raw-cell binomial s.e. <= {0.5/np.sqrt(NREAL):.3f}")
