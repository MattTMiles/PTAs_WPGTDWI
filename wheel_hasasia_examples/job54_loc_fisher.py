"""JOB 54 measure (a): Fisher blocks per pulsar by jacfwd; areas by Schur complement.

Parameters g = [theta, phi, ln f, ln h, phi0, iota, psi, ln Mc] + one distance per pulsar [pc].
Per pulsar a (white noise sigma):  A = Jg^T Jg/s^2, b = Jg^T jL/s^2, c = jL^T jL/s^2.
Marginalising L_a with Gaussian prior sigma_La:  M_a = A - b b^T/(c + 1/sigma_La^2)
   = A_perp + b b^T * p/(c (c+p)),  p = 1/sigma_La^2,  A_perp = Jperp^T Jperp/s^2 (stable form).
F8 = sum_a M_a ;  Sigma = F8^-1 ;  A90 = 2 pi ln(10) sin(theta) sqrt(det Sigma[0:2,0:2]) (sr).
"""
import numpy as np, jax, jax.numpy as jnp, os, time
import job54_loc_core as C


def res_g(g, L, toas, p_hat, h):
    return C.residual(toas, p_hat, L, g[0], g[1], g[6], jnp.cos(g[5]), jnp.exp(g[2]),
                      h * jnp.exp(g[3]), g[4], jnp.exp(g[7]))


@jax.jit
def blocks(g, h, toas, p_hats, L_pc):
    def one(p_hat, L):
        Jg = jax.jacfwd(res_g, 0)(g, L, toas, p_hat, h)       # (nt, 8)
        jL = jax.jacfwd(res_g, 1)(g, L, toas, p_hat, h)       # (nt,)
        r = res_g(g, L, toas, p_hat, h)
        c = jL @ jL
        b = Jg.T @ jL
        Jp = Jg - jnp.outer(jL, b / c)
        return Jp.T @ Jp, b, c, Jg.T @ Jg, r @ r
    return jax.vmap(one)(p_hats, L_pc)


def compute_blocks(freqs=(1e-8, 2e-8), mc=1e9):
    arr = C.build_all(); src = C.draw_sources()
    toas, P, L = jnp.asarray(arr["toas"]), jnp.asarray(arr["p_hats"]), jnp.asarray(arr["L_pc"])
    out = {}
    for f in freqs:
        Ap, B, Cc, Af, h1 = [], [], [], [], []
        for i in range(len(src["theta"])):
            s = C.src_at(src, i)
            g = jnp.array([s["theta"], s["phi"], np.log(f), 0.0, s["phase0"], np.arccos(s["cinc"]), s["psi"], np.log(mc)])
            a, b, c, af, rr = blocks(g, 1.0, toas, P, L)
            hu = C.SIGMA / np.sqrt(float(jnp.sum(rr)))          # h for SNR = 1
            sc = hu ** 2 / C.SIGMA ** 2                          # blocks at SNR = 1
            Ap.append(np.asarray(a) * sc); B.append(np.asarray(b) * sc); Cc.append(np.asarray(c) * sc)
            Af.append(np.asarray(af) * sc); h1.append(hu)
        tag = f"f{f*1e9:.0f}"
        out[tag + "_Aperp"] = np.array(Ap); out[tag + "_b"] = np.array(B); out[tag + "_c"] = np.array(Cc)
        out[tag + "_Afull"] = np.array(Af); out[tag + "_h_snr1"] = np.array(h1)
    np.savez(os.path.join(C.DATA, "fisher_blocks.npz"), **out, **{k: v for k, v in arr.items()},
             **{"src_" + k: v for k, v in src.items()})
    return out


def F8(Ap, b, c, sig, snr):
    """Ap [P,8,8], b [P,8], c [P] at SNR=1; sig [P] prior sigma [pc]."""
    s2 = snr ** 2
    p = 1.0 / sig ** 2
    w = p / (s2 * c * (s2 * c + p)) * s2 ** 2      # (s2 b)(s2 b)^T * p/(s2 c (s2 c + p))
    return s2 * Ap.sum(0) + np.einsum("p,pi,pj->ij", w, b, b)


def area_from_F(F, theta, rcond=None, cond_only=False):
    if cond_only:
        S = np.linalg.inv(F[:2, :2])
    else:
        d = 1.0 / np.sqrt(np.diag(F))
        Fn = F * np.outer(d, d)
        Si = np.linalg.inv(Fn) if rcond is None else np.linalg.pinv(Fn, rcond=rcond, hermitian=True)
        S = (Si * np.outer(d, d))[:2, :2]
    det = S[0, 0] * S[1, 1] - S[0, 1] ** 2
    return 2 * np.pi * np.log(10) * np.sin(theta) * np.sqrt(max(det, 0.0)) * C.DEG2


def sig_vec(arr, order, n_ph, sigL, unphased="free"):
    """Prior sigma per pulsar [pc]. unphased='free': unphased pulsar-term phase is a free nuisance
    (sigma = inf; the honest choice when sigma_L >> lambda_gw). unphased='census': literal Gaussian
    census prior (linearised; over-optimistic because it cannot represent phase wrapping)."""
    sig = np.full(len(arr["sigL_pc"]), np.inf) if unphased == "free" else arr["sigL_pc"].astype(float).copy()
    sig[order[:n_ph]] = sigL
    return sig


NPH_F = sorted(set(C.NPH_LIST + [30]))

if __name__ == "__main__":
    t0 = time.time()
    fb = os.path.join(C.DATA, "fisher_blocks.npz")
    if not os.path.exists(fb):
        compute_blocks()
    print("blocks done", time.time() - t0, flush=True)
    d = np.load(fb, allow_pickle=True)
    arr = {k: d[k] for k in ("L_pc", "sigL_pc", "order_near", "order_rand")}
    th = d["src_theta"]; res = {}; ns = len(th)
    for tag in ("f10", "f20"):
        Ap, b, c = d[tag + "_Aperp"], d[tag + "_b"], d[tag + "_c"]
        for snr in (10.0, 20.0):
            for oname in ("near", "rand"):
                order = arr["order_" + oname]
                for sl in (0.01, 0.1, 0.3, 1.0, 3.0, 10.0):
                    for n in NPH_F:
                        for mode in ("free", "census"):
                            sig = sig_vec(arr, order, n, sl, mode)
                            Fs = [F8(Ap[i], b[i], c[i], sig, snr) for i in range(ns)]
                            key = f"{tag}_snr{snr:.0f}_{oname}_sl{sl}_n{n}_{mode}"
                            res[key] = np.array([area_from_F(Fs[i], th[i]) for i in range(ns)])
                            res[key + "_cond"] = np.array([area_from_F(Fs[i], th[i], cond_only=True) for i in range(ns)])
    np.savez(os.path.join(C.DATA, "fisher_areas.npz"), **res)
    chk = []
    for tag in ("f10", "f20"):
        Ap, b, c = d[tag + "_Aperp"], d[tag + "_b"], d[tag + "_c"]
        for n in (0, 4, 116):
            sig = sig_vec(arr, arr["order_near"], n, 0.01)
            for i in range(ns):
                F = F8(Ap[i], b[i], c[i], sig, 10.0)
                a0 = area_from_F(F, th[i]); a1 = area_from_F(F, th[i], rcond=1e-13); a2 = area_from_F(F, th[i], rcond=1e-10)
                a3 = area_from_F(F8(Ap[i], b[i] * 1e3, c[i] * 1e6, sig / 1e3, 10.0), th[i])
                p = 1 / sig ** 2
                Fn = 100 * d[tag + "_Afull"][i].sum(0) - np.einsum("p,pi,pj->ij", 1e4 / (100 * c[i] + p), b[i], b[i])
                a4 = area_from_F(Fn, th[i])
                chk.append([n, a1 / a0, a2 / a0, a3 / a0, a4 / a0])
    chk = np.array(chk); np.save(os.path.join(C.DATA, "fisher_stability.npy"), chk)
    for n in (0, 4, 116):
        m = chk[:, 0] == n
        print("stability n_ph", n, "max |ratio-1| / median |ratio-1|: pinv1e-13 %.2e/%.2e pinv1e-10 %.2e/%.2e kpc-units %.2e/%.2e naive-schur %.2e/%.2e" %
              tuple(np.ravel(np.stack([np.nanmax(np.abs(chk[m, 1:] - 1), 0), np.nanmedian(np.abs(chk[m, 1:] - 1), 0)], 1))))
    for mode in ("free", "census"):
        for tag in ("f10", "f20"):
            for snr in (10, 20):
                for n in (0, 116):
                    A = res[f"{tag}_snr{snr}_near_sl0.01_n{n}_{mode}"]; Ac = res[f"{tag}_snr{snr}_near_sl0.01_n{n}_{mode}_cond"]
                    print(mode, tag, snr, "n_ph", n, "A90 marg med [16,84] = %.3g [%.3g, %.3g] deg2 ; cond med %.3g" %
                          (np.median(A), *np.percentile(A, [16, 84]), np.median(Ac)))
    print("total", time.time() - t0)
