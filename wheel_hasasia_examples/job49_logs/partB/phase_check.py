import numpy as np
from scipy.special import i0e, i1e
rng = np.random.default_rng(1)
n = 2_000_000
print("beta  rms_phase  sqrt((1-b)/b)  sqrt((1-b)/(2b))  <cos>sim  exp(-s1^2/2)  exp(-s2^2/2)  <cos>Rice")
for b in (0.99, 0.95, 0.9, 0.7, 0.5, 0.3, 0.2, 0.1):
    A = 1.0
    N = A * A * (1 - b) / b                      # total power of everything else in the bin
    z = rng.normal(0, np.sqrt(N / 2), n) + 1j * rng.normal(0, np.sqrt(N / 2), n)
    ph = np.angle(A + z)
    s1 = np.sqrt((1 - b) / b)
    s2 = np.sqrt((1 - b) / (2 * b))
    g = A * A / N
    rice = np.sqrt(np.pi) / 2 * np.sqrt(g) * (i0e(g / 2) + i1e(g / 2))   # i0e, i1e carry the e^{-g/2}
    print(f"{b:4.2f}  {ph.std():7.3f}  {s1:7.3f}  {s2:7.3f}  {np.cos(ph).mean():6.3f}  "
          f"{np.exp(-s1**2/2):6.3f}  {np.exp(-s2**2/2):6.3f}  {rice:6.3f}")
print("rms of a uniform phase:", np.pi / np.sqrt(3))
