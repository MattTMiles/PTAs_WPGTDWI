import numpy as np
rng = np.random.default_rng(12)

def run_weak(K, M, ntr, rot=0.0, chunk=50):
    # 10 moderate sources (true two-term share 0.261) + M weak sources (true share 0.01).
    # Physical data: unit phasor + complex Gaussian with variance Rs per quadrature. Matched filter.
    # The true candidate predicts phase 0 for every source; every other candidate a random phase.
    # rot = rms of an unknown template phase error on the WEAK sources only (rad).
    Rs = np.r_[np.full(10, 1.415), np.full(M, 49.5)]; N = len(Rs); ok = 0
    for c in range(0, ntr, chunk):
        n = min(chunk, ntr - c)
        D = rng.uniform(0, 2*np.pi, (n, N, K)); D[:, :, 0] = 0
        z = 1 + (rng.normal(0, 1, (n, N)) + 1j*rng.normal(0, 1, (n, N))) * np.sqrt(Rs)[None, :]
        if rot > 0:
            eta = np.zeros((n, N)); eta[:, 10:] = rng.normal(0, rot, (n, M)); z = z*np.exp(1j*eta)
        w = np.exp(-(np.r_[np.zeros(10), np.full(M, rot)])**2/2) / (2*Rs)
        S = np.real((z*w[None, :])[:, :, None] * np.exp(-1j*D)).sum(1)
        ok += (S.argmax(1) == 0).sum()
    return ok/ntr

for K, ntr in ((60, 3000), (600, 800)):
    print(K, [round(run_weak(K, M, ntr), 3) for M in (0, 100, 400, 1000)],
          [round(run_weak(K, 400, ntr, r), 3) for r in (1.0, 2.0, 3.0)])
