# ================= JOB 49 D2/D3 -- headline + amplitude limit, all arms ============
# usage: python job49_D2D3.py marg            (D2: the booked FINAL row, 36.27 @ 0.9009)
#        python job49_D2D3.py rate 0.33 ...   (D3: T3 rows at fixed count-scale)
# Protocol per config:
#   (1) REPRODUCE with untouched booked code: book(collect36m(...)) for marg;
#       book(_fixed()) exec'd from job43c_rebooking.py:790-829 for a rate.
#       Compare to the booked numbers (report / fig43_corner_vs_A.csv). If they
#       differ the config is flagged REPRO_FAIL and its arms are not booked.
#   (2) collect36m_sw on the same seeds; assert records bit-identical to (1).
#   (3) every arm: gate re-frozen on train skies (<100) by book()'s own procedure,
#       applied blind to held-out skies; plus the INCUMBENT threshold applied
#       unchanged (labelled).  Alias scores are data-independent -> shared.
import os, sys, json, time, pickle
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job49_common as J

BASE = J.BASE
os.chdir(BASE); os.makedirs("job49_data", exist_ok=True)
BOOKED = {"marg": dict(y=36.270, pur=0.9009, gate=58.59),
          0.33: dict(y=42.350, pur=0.9974, gate=58.84),
          0.60: dict(y=40.990, pur=0.9709, gate=58.77),
          1.0: dict(y=31.570, pur=0.9035, gate=50.32),
          1.7: None, 2.7: None}                  # None = booked "no 0.90 crossing"
A_OF = {0.33: 3.7, 0.60: 5.0, 1.0: 6.4, 1.7: 8.3, 2.7: 10.6}

t_all = time.time()
NS = J.load(verbose=False)
print(f"[load] {time.time()-t_all:.0f} s", flush=True)
L43 = open("job43c_rebooking.py").read().splitlines()


def jsonable(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, np.ndarray): out[k] = v.tolist()
        elif isinstance(v, (np.floating, np.integer)): out[k] = v.item()
        else: out[k] = v
    return out


def run(cfg):
    t0 = time.time()
    tag = "marg" if cfg == "marg" else f"rate{cfg:.2f}"
    # ---- (1) untouched booked code ------------------------------------------------
    if cfg == "marg":
        orig = NS["collect36m"](fuzzy=True, full_dmu=True, marg=True)
    else:
        NS["rate"] = cfg
        src = "\n".join(l[4:] for l in L43[790-1:829])       # def _fixed(...) body
        exec(compile(src, "<43c:790-829 _fixed>", "exec"), NS)
        orig = NS["_fixed"]()
    rb = NS["book"](orig)
    bk = BOOKED[cfg]
    if bk is None:
        repro_ok = not rb["ok"]
        repro_txt = f"booked: no 0.90 crossing; reproduced: {'no crossing' if not rb['ok'] else 'CROSSING at %.2f @ %.4f' % (rb['y'], rb['pur'])}"
    else:
        repro_ok = rb["ok"] and abs(rb["y"]-bk["y"]) < 0.01 and abs(rb["pur"]-bk["pur"]) < 0.0001
        repro_txt = (f"booked {bk['y']:.2f} @ {bk['pur']:.4f} gate {bk['gate']:.2f}; reproduced "
                     + (f"{rb['y']:.3f} [{rb['yb'][0]:.2f},{rb['yb'][1]:.2f}] @ {rb['pur']:.4f}"
                        f" [{rb['pb'][0]:.3f},{rb['pb'][1]:.3f}] gate {rb['gate']:.2f}" if rb["ok"] else "NO CROSSING"))
    print(f"[{tag}] REPRODUCTION {'OK' if repro_ok else 'FAIL'}: {repro_txt}", flush=True)
    # ---- (2) switch, incumbent default: bit-identical ---------------------------
    rec = J.collect36m_sw(fuzzy=True, full_dmu=True, marg=(cfg == "marg"),
                          rate=None if cfg == "marg" else cfg, keep_detail=True)
    same = len(rec) == len(orig) and all(
        a["sky"] == b["sky"] and a["wrong"] == b["wrong"] and np.array_equal(a["sig"], b["sig"])
        and np.array_equal(a["kap"], b["kap"]) and a["Ltrue"] == b["Ltrue"] for a, b in zip(orig, rec))
    print(f"[{tag}] switch(incumbent) bit-identical to booked collector: {same} (n={len(rec)})", flush=True)
    al = np.array([NS["alias_score24"](r) for r in rec])
    res = dict(cfg=tag, repro_ok=bool(repro_ok), repro_txt=repro_txt, identical=bool(same),
               n=len(rec), arms={})
    thr_inc = None
    for arm in J.ARMS:
        J.set_arm(rec, arm)
        b = J.book_sw(rec, alias=al)
        if arm == "inc":
            thr_inc = b["thr"] if b["ok"] else None
            assert (not b["ok"] and not rb["ok"]) or (abs(b["y"]-rb["y"]) < 1e-12
                                                      and abs(b["pur"]-rb["pur"]) < 1e-12)
        bo = J.book_sw(rec, alias=al, thr_fixed=thr_inc) if thr_inc is not None else None
        ung = 1-np.mean([r["wrong"] for r in rec])
        res["arms"][arm] = dict(refrozen=jsonable(b), oldgate=jsonable(bo) if bo else None,
                                ungated_purity=float(ung))
        s = (f"{b['y']:6.2f} [{b['yb'][0]:.2f},{b['yb'][1]:.2f}] @ {b['pur']:.4f} "
             f"[{b['pb'][0]:.3f},{b['pb'][1]:.3f}] gate {b['gate']:.1f}%" if b["ok"]
             else f"no 0.90 crossing (max train {b['maxpur']:.4f})")
        s2 = (f" | old thr: {bo['y']:6.2f} @ {bo['pur']:.4f}" if bo else "")
        print(f"[{tag}] {arm:>8}: ungated {ung:.4f} | {s}{s2}", flush=True)
    J.set_arm(rec, "inc")
    with open(f"job49_data/rec_{tag}.pkl", "wb") as fh:
        pickle.dump(dict(rec=rec, alias=al), fh)
    json.dump(res, open(f"job49_data/D2D3_{tag}.json", "w"), indent=1)
    print(f"[{tag}] done {time.time()-t0:.0f} s", flush=True)


if sys.argv[1] == "marg":
    run("marg")
else:
    for x in sys.argv[2:]:
        run(float(x))
print(f"TOTAL {time.time()-t_all:.0f} s")
