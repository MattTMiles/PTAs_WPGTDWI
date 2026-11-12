# ================= JOB 49 PART C -- unit test of the machine's phase-error generator =
# The machine has no generator FUNCTION: the phase error is produced inline.  This
# script execs the BOOKED SOURCE LINES themselves, read from the booked files at
# run time (not retyped), with measurement noise off (rho -> 1e12) and kappa = 0,
# L_true = 0, so eps IS the phase error:
#   share -> sigma : job43c_rebooking.py:505   (inv_rate; the headline path)
#                    job47c_heatmap.py s <-> sigma_phi (header l.10; axis l.231)
#   generator      : job43c_rebooking.py:749   (collect36m; the headline path)
#                    job47c_heatmap.py:124      (purity_cell; the heatmap path)
# References: exact Rice (physical complex sum), wrapped normal var (1-b)/b,
# wrapped normal var (1-b)/(2b).
import os, csv, json
import numpy as np
from scipy.special import i0e, i1e, erf
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
os.chdir(BASE)
L43 = open("job43c_rebooking.py").read().splitlines()
L47 = open("job47c_heatmap.py").read().splitlines()
SRC_SIG = L43[505-1].strip()        # sg[a]=np.sqrt(max(res,0)/max(Pt,1e-300)+1/max(rp[k],1e-9)**2)
SRC_GEN43 = L43[749-1].strip()      # eps=kap*Lt+rng.normal(0,sgu)
SRC_GEN47 = L47[124-1].strip()      # eps = kap*Lt + rg.normal(0, sigphi, len(kap))
print("booked lines exec'd:\n  ", SRC_SIG, "\n  ", SRC_GEN43, "\n  ", SRC_GEN47)
assert "rng.normal(0,sgu)" in SRC_GEN43 and "rg.normal(0, sigphi" in SRC_GEN47

N = 200_000
BETAS = (0.99, 0.9, 0.7, 0.5, 0.3, 0.2, 0.1)
wrap = lambda x: (x + np.pi) % (2*np.pi) - np.pi


def machine_sigma_43(beta):
    ns = dict(np=np, res=(1-beta)/beta, Pt=1.0, rp={0: 1e12}, k=0, sg=np.zeros(1), a=0)
    exec(SRC_SIG, ns)
    return float(ns["sg"][0])


def machine_eps_43(beta, seed):
    sg = machine_sigma_43(beta)
    ns = dict(np=np, kap=np.zeros(N), Lt=0.0, sgu=np.full(N, sg),
              rng=np.random.default_rng(seed))
    exec(SRC_GEN43, ns)
    return ns["eps"], sg


def machine_eps_47(beta, seed):
    sigphi = np.sqrt((1-beta)/beta)          # 47c: s = 1/(1+sigma_phi^2)
    ns = dict(np=np, kap=np.zeros(N), Lt=0.0, sigphi=sigphi,
              rg=np.random.default_rng(seed))
    exec(SRC_GEN47, ns)
    return ns["eps"], sigphi


def rice_cos(beta):
    g = beta/(1-beta)
    return np.sqrt(np.pi)/2*np.sqrt(g)*(i0e(g/2)+i1e(g/2))


def rice_pdf(phi, beta):
    g = beta/(1-beta); c = np.cos(phi)
    return np.exp(-g)/(2*np.pi)*(1+np.sqrt(np.pi*g)*c*np.exp(g*c**2)*(1+erf(np.sqrt(g)*c)))


def wn_pdf(phi, s):
    k = np.arange(-8, 9)[:, None]
    return np.exp(-(phi[None]+2*np.pi*k)**2/(2*s*s)).sum(0)/np.sqrt(2*np.pi*s*s)


PH = np.linspace(-np.pi, np.pi, 20001)
def rms_of(pdf):
    p = pdf/np.trapezoid(pdf, PH)
    return float(np.sqrt(np.trapezoid(PH**2*p, PH)))


rows = []
print(f"\n{'beta':>5} | {'machine43 rms(wrapped)':>22} {'<cos>':>13} | {'machine47 <cos>':>15} |"
      f" {'Rice rms':>8} {'<cos>':>6} | {'WN(1-b)/b rms':>13} {'<cos>':>6} | {'WN(1-b)/2b rms':>14} {'<cos>':>6}")
for n_, b in enumerate(BETAS):
    e43, s43 = machine_eps_43(b, 4900+n_)
    e47, s47 = machine_eps_47(b, 4950+n_)
    assert abs(s43-np.sqrt((1-b)/b)) < 1e-9 and abs(s47-np.sqrt((1-b)/b)) < 1e-12
    w43 = wrap(e43)
    c43 = np.cos(e43); c47 = np.cos(e47)
    s1 = np.sqrt((1-b)/b); s2 = np.sqrt((1-b)/(2*b))
    row = dict(beta=b, sigma_machine=s43,
               m43_rms_raw=float(e43.std()), m43_rms_wrapped=float(w43.std()),
               m43_cos=float(c43.mean()), m43_cos_se=float(c43.std()/np.sqrt(N)),
               m47_rms_wrapped=float(wrap(e47).std()), m47_cos=float(c47.mean()),
               m47_cos_se=float(c47.std()/np.sqrt(N)),
               rice_rms=rms_of(rice_pdf(PH, b)), rice_cos=float(rice_cos(b)),
               wn1_rms=rms_of(wn_pdf(PH, s1)), wn1_cos=float(np.exp(-s1**2/2)),
               wn2_rms=rms_of(wn_pdf(PH, s2)), wn2_cos=float(np.exp(-s2**2/2)))
    # which reference does the machine match? |delta<cos>| in units of MC s.e.
    for ref in ("rice", "wn1", "wn2"):
        row[f"z_{ref}"] = (row["m43_cos"]-row[f"{ref}_cos"])/row["m43_cos_se"]
    rows.append(row)
    print(f"{b:5.2f} | {row['m43_rms_wrapped']:10.3f} (raw {row['m43_rms_raw']:5.3f})"
          f" {row['m43_cos']:6.3f}+-{row['m43_cos_se']:.3f} | {row['m47_cos']:8.3f}+-{row['m47_cos_se']:.3f} |"
          f" {row['rice_rms']:8.3f} {row['rice_cos']:6.3f} | {row['wn1_rms']:13.3f} {row['wn1_cos']:6.3f} |"
          f" {row['wn2_rms']:14.3f} {row['wn2_cos']:6.3f}")
print("\n  z = (machine <cos> - reference)/MC s.e.:")
for r in rows:
    print(f"   beta {r['beta']:4.2f}:  vs Rice {r['z_rice']:+9.1f}   vs WN(1-b)/b {r['z_wn1']:+6.2f}"
          f"   vs WN(1-b)/2b {r['z_wn2']:+9.1f}")
match = all(abs(r["z_wn1"]) < 4 for r in rows)
match_rice = any(abs(r["z_rice"]) < 4 for r in rows if r["beta"] < 0.95)
print(f"\n  VERDICT: machine matches wrapped normal with var (1-b)/b at every beta: {match};"
      f" matches Rice anywhere below b=0.95: {match_rice}")

with open("job49_partC_table.csv", "w", newline="") as fh:
    cw = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); cw.writeheader()
    for r in rows: cw.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
print("  wrote job49_partC_table.csv")

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for a, b in zip(ax, (0.5, 0.2)):
    e43, _ = machine_eps_43(b, 4900+BETAS.index(b))
    a.hist(wrap(e43), bins=90, range=(-np.pi, np.pi), density=True, histtype="stepfilled",
           color="0.8", label="machine (job43c:749), wrapped")
    a.plot(PH, rice_pdf(PH, b), color="C0", lw=1.6, label="physical sum (Rice, exact)")
    a.plot(PH, wn_pdf(PH, np.sqrt((1-b)/b)), color="C3", lw=1.3, ls="--",
           label=r"wrapped normal, var $(1-\beta)/\beta$")
    a.plot(PH, wn_pdf(PH, np.sqrt((1-b)/(2*b))), color="C2", lw=1.3, ls=":",
           label=r"wrapped normal, var $(1-\beta)/(2\beta)$")
    a.set_title(rf"$\beta$ = {b}"); a.set_xlabel("phase error [rad]"); a.set_xlim(-np.pi, np.pi)
    a.grid(alpha=.25, lw=.5)
ax[0].set_ylabel("density"); ax[0].legend(fontsize=7.5, loc="upper left")
plt.tight_layout(); plt.savefig("job49_partC_hist.png", dpi=150); plt.close()
print("  wrote job49_partC_hist.png")
json.dump(dict(verdict_wn1=bool(match), verdict_rice=bool(match_rice)),
          open("job49_partC_verdict.json", "w"))
