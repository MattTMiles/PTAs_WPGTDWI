import numpy as np
rng = np.random.default_rng(7)

def wn_logpdf(x, sig):
    k = np.arange(-6, 7)[:, None, None, None]
    return np.log(np.exp(-(x[None] + 2 * np.pi * k) ** 2 / (2 * sig ** 2)).sum(0) + 1e-300)

def run(N, K, beta, ntr=40000):
    P = (1 - beta) / beta                        # power of everything else / source power
    sig = np.sqrt(P)                             # the project's sigma_phi
    D = rng.uniform(0, 2 * np.pi, (ntr, N, K)); D[:, :, 0] = 0.0   # candidate 0 is the truth
    # physical: source phasor (amplitude 1) + complex Gaussian of total power P
    z = 1 + rng.normal(0, np.sqrt(P / 2), (ntr, N)) + 1j * rng.normal(0, np.sqrt(P / 2), (ntr, N))
    S_mf = np.real(z[:, :, None] * np.exp(-1j * D)).sum(1)          # matched filter
    S_ph = np.cos(np.angle(z)[:, :, None] - D).sum(1)               # phase-only cosine stack
    # phase error drawn as Gaussian with rms sigma_phi
    eps = rng.normal(0, sig, (ntr, N))
    S_wn = wn_logpdf(((eps[:, :, None] - D + np.pi) % (2 * np.pi)) - np.pi, sig).sum(1)
    # phase error drawn with the small-angle variance, sigma_phi^2 / 2
    eps2 = rng.normal(0, sig / np.sqrt(2), (ntr, N))
    S_w2 = wn_logpdf(((eps2[:, :, None] - D + np.pi) % (2 * np.pi)) - np.pi, sig / np.sqrt(2)).sum(1)
    f = lambda S: (S.argmax(1) == 0).mean()
    return f(S_mf), f(S_ph), f(S_wn), f(S_w2), sig

print("N  K  beta  sigma_phi | physical MF, physical phase-only | drawn sigma_phi | drawn sigma_phi/sqrt2")
for N, K, b in [(4, 11, 0.25), (4, 11, 0.45), (4, 11, 0.65), (8, 11, 1 / (1 + 1.81 ** 2)), (8, 11, 0.5),
                (4, 6, 0.2), (4, 6, 0.5), (4, 6, 0.8), (4, 60, 0.2), (4, 60, 0.5), (4, 60, 0.8),
                (4, 600, 0.2), (4, 600, 0.5), (4, 600, 0.8), (4, 600, 0.95)]:
    a, p, w, w2, s = run(N, K, b, 40000 if K <= 60 else 6000)
    print(f"{N}  {K:4d}  {b:4.2f}  {s:4.2f} | {a:5.3f}  {p:5.3f} | {w:5.3f} | {w2:5.3f}")
