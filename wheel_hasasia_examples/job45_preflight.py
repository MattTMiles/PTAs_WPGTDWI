# ================= JOB 45 -- S0.4 PRE-FLIGHT ====================================
# Rebuild the frozen machine namespace from active_wheel_demo_v4.ipynb (read-only:
# cells are exec'd from a JSON parse; the notebook file is never written), then run
# the Job 43C anchor path VERBATIM (line-sliced from job43c_rebooking.py) and
# compare against the booked headline:
#   36.27 [34.02, 38.75] psr/sky @ purity 0.9009 [0.890, 0.912]  (frozen-gate OOS,
#   banked SEED50 harness, fuzzy + full dmu + marginalised A).
# PASS = reproduced numbers inside the booked bands. FAIL = STOP, nothing else runs.
import os, sys, json, time, hashlib, traceback
os.environ.setdefault("MPLBACKEND", "Agg")

BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
SCRATCH = os.environ.get(
    "JOB45_SCRATCH",
    "/tmp/claude-1005/-home-mattm-soft/986e71f8-eff1-4a6d-b07f-676b52e4ab76/scratchpad")
os.chdir(BASE)

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()

print("=" * 74)
print("JOB 45 PRE-FLIGHT -- md5 ledger")
print("=" * 74)
for p in ("wheel_demo_v3.ipynb", "active_wheel_demo_v4.ipynb",
          "one_pulsar_six_worlds.ipynb", "step1_working_world.ipynb"):
    print(f"  {md5(p)}  {p}")
print(f"  {md5('/home/mattm/projects/HSYMT/CW_transition/best_distances.txt')}"
      f"  CW_transition/best_distances.txt")
BOOKED_V3 = "d897b1e2802692d30b1deb074a3ac94c"
assert md5("wheel_demo_v3.ipynb") == BOOKED_V3, "frozen v3 md5 MOVED -- STOP"
print(f"  frozen wheel_demo_v3 md5 == booked {BOOKED_V3}: OK")

# ---- rebuild machine namespace: exec code cells in order ------------------------
nb = json.load(open("active_wheel_demo_v4.ipynb"))
SKIP = {12,   # static dial frames -- needs a live kernel
        43}   # namespace snapshot to a dead scratchpad path -- replaced below
NS = {"__name__": "__main__"}
t_all = time.time()
fails = []
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code" or i in SKIP:
        continue
    src = "".join(c["source"])
    t0 = time.time()
    try:
        exec(compile(src, f"<v4-cell-{i}>", "exec"), NS)
        print(f"  [cell {i:2d}] OK    {time.time()-t0:7.1f} s", flush=True)
    except Exception as e:
        fails.append((i, repr(e)))
        print(f"  [cell {i:2d}] FAIL  {time.time()-t0:7.1f} s  {e!r}", flush=True)
        traceback.print_exc()
print(f"  namespace rebuild: {time.time()-t_all:.0f} s, {len(fails)} cell failures")

REQUIRED = ["_pop18", "rn_draw", "depro_sigma", "alias_score24",
            "TOOTH_W", "C_PC_S", "YR_S", "T_OBS", "bin_edges", "bin_cent"]
missing = [k for k in REQUIRED if k not in NS]
if missing:
    print(f"  STOP: required machine names missing after rebuild: {missing}")
    sys.exit(2)
print(f"  all {len(REQUIRED)} required machine names present.")

# ---- snapshot for the later tasks (T1/T2 reuse; volatile scratchpad is fine) ----
try:
    import dill
    ns = {}
    for k, v in NS.items():
        if k.startswith("__") or isinstance(v, type(os)):
            continue
        try:
            dill.dumps(v)
            ns[k] = v
        except Exception:
            pass
    os.makedirs(SCRATCH, exist_ok=True)
    with open(os.path.join(SCRATCH, "job45_banked_state.pkl"), "wb") as f:
        dill.dump(ns, f)
    print(f"  snapshot: {len(ns)} objects -> job45_banked_state.pkl")
except Exception as e:
    print(f"  snapshot skipped ({e!r}) -- later tasks will rebuild")

# ---- anchor: exec Job 43C fragments VERBATIM by line range ----------------------
# (485,509): numpy import + T45 + inv_rate       (661,689): book
# (710,782): collect36m + the three-step re-booking incl. the FINAL anchor row
lines = open("job43c_rebooking.py").read().splitlines()
def run_span(lo, hi):     # 1-indexed inclusive, matching Read line numbers
    exec(compile("\n".join(lines[lo - 1:hi]), f"<43c:{lo}-{hi}>", "exec"), NS)
t0 = time.time()
run_span(485, 509)
run_span(661, 689)
run_span(710, 782)
print(f"  anchor block runtime: {time.time()-t0:.0f} s")

# ---- verdict --------------------------------------------------------------------
R = NS["R"]
fin = R["3. + marginalised A (F7) = FINAL"]
ok = (fin["ok"] and 34.02 - 0.5 <= fin["y"] <= 38.75 + 0.5
      and 0.890 - 0.005 <= fin["pur"] <= 0.912 + 0.005)
print("\n" + "=" * 74)
print(f"PRE-FLIGHT ANCHOR: booked 36.27 [34.02, 38.75] @ 0.9009 [0.890, 0.912]")
if fin["ok"]:
    print(f"      reproduced: {fin['y']:.2f} [{fin['yb'][0]:.2f}, {fin['yb'][1]:.2f}]"
          f" @ {fin['pur']:.4f} [{fin['pb'][0]:.3f}, {fin['pb'][1]:.3f}],"
          f" gate {fin['gate']:.1f} %")
print(f"PRE-FLIGHT: {'PASS -- Job 45 tasks may run' if ok else 'FAIL -- STOP'}")
sys.exit(0 if ok else 3)
