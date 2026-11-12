# ================= JOB 50 PART A -- numerical tests behind the facts table ==========
# A1: tau = 0 toggle per path, 10 skies.  A3: the 0.303.  A5: delta-kappa lever arms.
import os, sys, json, pickle, time
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job49_common as J
os.chdir(J.BASE); os.makedirs("job50_data", exist_ok=True)
NS = J.load(verbose=False)
OUT = {}
q = lambda x: [float(v) for v in np.quantile(x, [0.1, 0.5, 0.9])]
f_src, FP0 = NS["f_src"], NS["FP_ALL"].copy()
print("namespace has:", {k: (k in NS) for k in ("SRC_AMP", "SRC_MC", "SRC_Z", "mc_src", "SRC_P", "DET_A")})

# ---------------- (ii)+(iii)+(vii) harness: inventory share, decode share, p2 share ----
def harness(nsky=10):
    rec = J.collect36m_sw(nsky=nsky, fuzzy=True, full_dmu=True, marg=True)
    s_dec = np.concatenate([1/(1+r["sig"]**2) for r in rec])
    s_p2 = np.concatenate([1/(1+2*r["Xc"]+2/r["rho"]**2) for r in rec])
    inv = []
    rgA = np.random.default_rng(555)
    for j in range(nsky):
        A = 6.4e-15*np.exp(rgA.normal(0, 0.32)); rt = float(np.clip((A/6.4e-15)**2, 0.2, 3.5))
        sd, f2, h2, P0, sea = NS["_pop18"](j % 50, 0.35)
        r = NS["inv_rate"](f2, h2, P0, sea, 45, rate=rt)
        if r is not None:
            inv.append(1/(1+np.maximum(r[1]**2-1/r[3]**2, 0)))
    wr = np.array([[r["wrong_inc"], r["wrong_p2_mf"]] for r in rec])
    return dict(n=len(rec), inv=np.concatenate(inv), dec=s_dec, p2=s_p2, wr=wr)

base = harness()
# toggle 1: FP_ALL := f_src (tau = 0) WITHOUT re-deriving anything (does any booked fn read it?)
NS["FP_ALL"] = np.tile(f_src, (FP0.shape[0], 1))
t1 = harness()
# toggle 2: also re-derive depro's f_E table from the tau = 0 ledger (re-exec cell 17:1-127)
nb = json.load(open("active_wheel_demo_v4.ipynb"))
c17 = "\n".join("".join(nb["cells"][17]["source"]).splitlines()[:127])
import io, contextlib
fE_before = [NS["fE_of"](f) for f in (3e-9, 1e-8, 2.5e-8)]
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(c17, "<cell17 tau0>", "exec"), NS)
fE_after = [NS["fE_of"](f) for f in (3e-9, 1e-8, 2.5e-8)]
t2 = harness()
NS["FP_ALL"] = FP0
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(c17, "<cell17 restore>", "exec"), NS)
chk = harness()
assert np.array_equal(chk["dec"], base["dec"]), "restore failed"
print("\n== A1 harness paths, 10 headline skies (marg): share 10/50/90 pct  [per (pulsar,source) decode pair unless said]")
for lab, key in (("(ii) inv_rate share, per sky-source (pre-depro)", "inv"),
                 ("(iii) headline decode share 1/(1+sig^2)", "dec"),
                 ("(vii) p2 true two-term share", "p2")):
    print(f"  {lab}:")
    for nm, d in (("booked", base), ("FP_ALL:=f_E only", t1), ("FP_ALL:=f_E + f_E table re-derived", t2)):
        print(f"     {nm:>36}: n={len(d[key]):6d}  {np.round(q(d[key]), 4)}  identical to booked: {np.array_equal(d[key], base[key])}")
    OUT[key] = dict(booked=q(base[key]), tau0_only=q(t1[key]), tau0_fE=q(t2[key]),
                    ident_only=bool(np.array_equal(t1[key], base[key])),
                    ident_fE=bool(len(t2[key]) == len(base[key]) and np.array_equal(t2[key], base[key])))
print(f"  f_E table at 3/10/25 nHz: booked {np.round(fE_before, 3)} -> tau=0 {np.round(fE_after, 3)}")
print(f"  ungated purity inc/p2_mf: booked {1-base['wr'].mean(0)}  tau0+fE {1-t2['wr'].mean(0)} (n {base['n']} vs {t2['n']})")
OUT["fE"] = dict(before=fE_before, after=fE_after)

# ---------------- (iv) heatmap comb sets --------------------------------------------
def sets():
    out = []
    for j in range(10):
        sd, f2, h2, P0, sea = NS["_pop18"](j, 0.35)
        r = NS["inv_rate"](f2, h2, P0, sea, 45.0); out.append(None if r is None else r[2][:4])
    return out
s0 = sets(); NS["FP_ALL"] = np.tile(f_src, (FP0.shape[0], 1)); s1 = sets(); NS["FP_ALL"] = FP0
print(f"\n  (iv) heatmap comb sets (10 skies, top-4 f_k) identical under tau=0: "
      f"{all((a is None and b is None) or np.array_equal(a, b) for a, b in zip(s0, s1))};"
      f" sigma_phi is set by hand (job47c:132), so the share axis cannot move.")

# ---------------- (v) job46 share-vs-frequency figure (1 sky, pulsar 41) -------------
h_src, T = NS["h_src"], NS["T_OBS"]; edges = NS["bin_edges"]; sea = NS["S_FAINT"]; nbn = len(edges)-1
def shares(f, h):
    P = h**2*f*T; b = np.digitize(f, edges)-1; ok = (b >= 0) & (b < nbn)
    tot = np.bincount(b[ok], weights=P[ok], minlength=nbn)+sea
    s = np.full(len(f), np.nan); s[ok] = P[ok]/tot[b[ok]]; return s
fP = FP0[41]; hP = h_src*(fP/f_src)**(2/3)
sON, sOFF = shares(fP, hP), shares(f_src, h_src)
inb = lambda f: (f >= 15e-9) & (f < 40e-9)
top = lambda s, m: np.sort(s[m & np.isfinite(s)])[::-1][:20]
print(f"\n  (v) job46 figure, presentation sky, pulsar 41, 15-40 nHz: n terms ON {inb(fP).sum()} / OFF {inb(f_src).sum()};"
      f" share 10/50/90: ON {np.round(q(sON[inb(fP)]), 5)} OFF {np.round(q(sOFF[inb(f_src)]), 5)};"
      f" median of top-20 shares ON {np.median(top(sON, inb(fP))):.3f} OFF {np.median(top(sOFF, inb(f_src))):.3f}")
OUT["job46"] = dict(on=q(sON[inb(fP)]), off=q(sOFF[inb(f_src)]),
                    top20_on=float(np.median(top(sON, inb(fP)))), top20_off=float(np.median(top(sOFF, inb(f_src)))))

# ---------------- (i) notebook V4o audited share: target at retarded f, others E + P ----
# sigma_audited's ledger (v4 cell 14:1598-1617) with every source treated as undetected
# (the detected-set bookkeeping needs the hasasia curve; the retardation logic is the same)
P_A = h_src**2*f_src*T
loud = np.argsort(P_A)[::-1][:60]
def share_aud(i, ks, FP):
    out = []
    for k in ks:
        fpk = FP[i, k]
        w_e = np.sinc((f_src-fpk)*T)**2; w_p = np.sinc((FP[i]-fpk)*T)**2
        m = np.ones(len(f_src), bool); m[k] = False
        b = int(np.clip(np.digitize(fpk, edges)-1, 0, nbn-1))
        res = np.sum(w_e[m]*P_A[m])+np.sum(w_p[m]*P_A[m])+sea[b]
        out.append(P_A[k]/(P_A[k]+res))
    return np.array(out)
FPz = np.tile(f_src, (FP0.shape[0], 1))
a_on = np.concatenate([share_aud(i, loud, FP0) for i in (3, 41, 77)])
a_off = np.concatenate([share_aud(i, loud, FPz) for i in (3, 41, 77)])
print(f"\n  (i) notebook V4o audited ledger (60 loudest sources x pulsars 3,41,77; all-undetected variant):"
      f" share 10/50/90 retarded {np.round(q(a_on), 4)}  tau=0 {np.round(q(a_off), 4)}")
OUT["v4o"] = dict(on=q(a_on), off=q(a_off))
s_old = P_A[loud]/np.maximum(NS["hc2_bin"][np.clip(np.digitize(f_src[loud], edges)-1, 0, nbn-1)], 1e-300)
print(f"      notebook sigma_phi() bin share (cell 3:262-276, no tau input): {np.round(q(s_old), 4)}")

# ---------------- (vi) Job 41 displaced excess, 10 skies ------------------------------
L41 = open("job41_displaced_excess.py").read().splitlines()
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile("\n".join(L41[193-1:222]), "<job41:193-222>", "exec"), NS)
sig = np.mean([NS["excess_spectrum"](j) for j in range(10)], axis=0)
nul = np.mean([NS["excess_spectrum"](j, null=True) for j in range(10)], axis=0)
print(f"\n  (vi) Job 41 in-band/above-50 excess ratio, 10 skies: retarded {NS['ratio_of'](sig):.2f}  tau=0 {NS['ratio_of'](nul):.2f}"
      f"  (booked, 20 skies: 18.7 / 1.66)")
OUT["job41"] = dict(on=float(NS["ratio_of"](sig)), off=float(NS["ratio_of"](nul)))

# ---------------- A3: the 0.303 ------------------------------------------------------
R = pickle.load(open("job49_data/rec_marg.pkl", "rb"))["rec"]
sig_ = np.concatenate([r["sig"] for r in R]); Xc = np.concatenate([r["Xc"] for r in R])
rho = np.concatenate([r["rho"] for r in R]); sfz = np.concatenate([r["sfz"] for r in R])
s_m = 1/(1+sig_**2); s_nf = 1/(1+Xc+1/rho**2); s_c = 1/(1+Xc); b2 = 1/(1+2*Xc+2/rho**2)
print(f"\n== A3 (200 headline skies, {len(sig_)} decode pairs)")
print(f"  median machine share 1/(1+sigma^2), sigma^2 = Xc + 1/rho^2 + sfz^2 : {np.median(s_m):.4f}  -> s/(2-s) = {np.median(s_m)/(2-np.median(s_m)):.4f}")
print(f"  median 1/(1+Xc+1/rho^2) (no template term)                        : {np.median(s_nf):.4f}  -> s/(2-s) = {np.median(s_nf)/(2-np.median(s_nf)):.4f}")
print(f"  median 1/(1+Xc)  (confusion only)                                 : {np.median(s_c):.4f}  -> s/(2-s) = {np.median(s_c)/(2-np.median(s_c)):.4f}")
print(f"  median Job-49 'true two-term share' 1/(1+2Xc+2/rho^2)             : {np.median(b2):.4f}")
print(f"  variance budget medians: Xc {np.median(Xc):.3f}, 1/rho^2 {np.median(1/rho**2):.4f}, sfz^2 {np.median(sfz**2):.4f};"
      f" mean share of sigma^2 from the template term {np.mean(sfz**2/sig_**2):.3f}")
OUT["A3"] = dict(s_machine=float(np.median(s_m)), s_nofuzzy=float(np.median(s_nf)), s_conf=float(np.median(s_c)),
                 b2=float(np.median(b2)), tmpl_frac=float(np.mean(sfz**2/sig_**2)))

# ---------------- A5: delta-kappa lever arms -----------------------------------------
# replicate collect36m's rng stream (job43c:712-756) to recover Dk (true distance, kpc)
def lever(nsky=200, T=45, seed0=8000):
    rows = []; rgA = np.random.default_rng(555); n = 0
    for j in range(nsky):
        A = 6.4e-15*np.exp(rgA.normal(0, 0.32)); rt = float(np.clip((A/6.4e-15)**2, 0.2, 3.5))
        sd, f2, h2, P0, sea = NS["_pop18"](j % 50, 0.35)
        rng = np.random.default_rng(seed0+sd+977*j); Np = 120
        beta = np.arcsin(rng.uniform(-1, 1, Np)); Dk = rng.uniform(0.5, 3.0, Np)
        px = 2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/T); pb = 3.0*(T/30.0)**(-2.5)
        pri = np.maximum(np.where(rng.random(Np) < 0.7, np.minimum(px, pb), px), 0.05)
        Arn, _g = NS["rn_draw"](Np); rncl = Arn == 0.0
        r = NS["inv_rate"](f2, h2, P0, sea, T, rate=rt)
        if r is None: continue
        ks, sig, fk, rp, arr = r
        for i in range(Np):
            if not rncl[i]: continue
            sg = np.array([NS["depro_sigma"](np.sqrt(max(sig[a]**2-1/max(rp[a], 1e-9)**2, 0.0)),
                                             rp[a], fk[a], i, "loo") for a in range(len(ks))])
            ok = sg < np.pi
            if ok.sum() < 2: continue
            sgu = sg[ok]; fku = fk[ok]; arru = arr[ok]
            mu = np.arccos(rng.uniform(-0.9, 0.9, len(sgu)))
            kap = 2*np.pi*fku*(1-np.cos(mu))/NS["C_PC_S"]
            md = 1.0/np.maximum(arru, 1e-6); dmu = rng.rayleigh(md*np.sqrt(2/np.pi))
            dfq = 1.0/(T*NS["YR_S"]*np.maximum(arru, 1e-6))
            tf = (2*np.pi*(1-np.cos(mu))/NS["C_PC_S"])*dfq; tm = (2*np.pi*fku*np.sin(mu)/NS["C_PC_S"])*dmu
            dkap = np.sqrt(tf**2+tm**2)
            sgu2 = np.sqrt(sgu**2+(dkap*pri[i])**2); o2 = sgu2 < np.pi
            if o2.sum() < 2: continue
            Lt = float(np.clip(rng.normal(0, pri[i]), -2.9*pri[i], 2.9*pri[i])); rng.normal(0, sgu2[o2])
            if 1.0/np.sqrt(max(float(np.sum(kap[o2]**2/sgu2[o2]**2)), 1e-30)) >= NS["TOOTH_W"]: continue
            assert np.allclose(sgu2[o2], R[n]["sig"], rtol=1e-9); n += 1
            for c in np.where(o2)[0]:
                rows.append((dkap[c]*pri[i], dkap[c]*Dk[i]*1e3, tf[c]*Dk[i]*1e3, tm[c]*Dk[i]*1e3, arru[c], pri[i], Dk[i]))
    assert n == len(R)
    return np.array(rows)
LV = lever()
print(f"\n== A5 (200 headline skies, {len(LV)} decode pairs; stream check vs rec_marg.pkl passed)")
for lab, c in (("dkappa*sigma_d [rad]", 0), ("dkappa*L [rad]", 1), ("  frequency part 2pi(1-cos mu) df L/c", 2),
               ("  sky-position part 2pi f sin(mu) dmu L/c", 3), ("array SNR of the source", 4)):
    print(f"  {lab:>42}: 10/50/90 pct = {np.array2string(np.quantile(LV[:, c], [.1, .5, .9]), precision=3)}")
print(f"  fraction of decode pairs with dkappa*L < 1 rad: {np.mean(LV[:, 1] < 1):.4f}; < pi: {np.mean(LV[:, 1] < np.pi):.4f}")
OUT["A5"] = {k: q(LV[:, c]) for k, c in (("dk_sigd", 0), ("dk_L", 1), ("dk_L_freq", 2), ("dk_L_sky", 3), ("arr", 4))}
OUT["A5"]["frac_lt1"] = float(np.mean(LV[:, 1] < 1))
np.save("job50_data/A5_lever.npy", LV)
json.dump(OUT, open("job50_data/partA.json", "w"), indent=1)
print("\nwrote job50_data/partA.json")
