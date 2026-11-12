#!/usr/bin/env python
"""Run a measurement script against the banked namespace. No notebook re-run.
Restored notebook functions resolve globals against __main__, so the banked
namespace is injected THERE rather than into a local dict."""
import sys, os, dill, numpy as np, matplotlib
matplotlib.use("Agg")
for _p in ("/home/mattm/projects/HSYMT/hpc_harbor/vernier",
           "/home/mattm/projects/HSYMT/wheel_hasasia_examples",
           "/home/mattm/projects/HSYMT"):
    if _p not in sys.path:
        sys.path.insert(0, _p)
SNAP = "/tmp/claude-1005/-home-mattm-soft/861befe2-5c21-43c4-97b2-391ae87f6ed3/scratchpad/banked_state.pkl"
script = os.path.abspath(sys.argv[1])
code = open(script).read()
with open(SNAP, "rb") as f:
    ns = dill.load(f)
print(f"[fastrun] loaded {ns.get('_SNAP_STAMP',{}).get('n_keys','?')} banked objects")
import __main__
__main__.__dict__.update(ns)
__main__.__dict__["np"] = np
import scipy.stats, scipy.special
__main__.__dict__.setdefault("norm", scipy.stats.norm)
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")
exec(compile(code, script, "exec"), __main__.__dict__)
