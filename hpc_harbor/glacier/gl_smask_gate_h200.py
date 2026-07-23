#!/usr/bin/env python3
"""FORGE-G2 smask battery ON THE H200 LANE -- reference-absent mode, declared.

The baked reference (CW_transition/b1_smask_baked_ref.npz) is a CRONUS (RTX 4090) capture.
Cross-host ABSOLUTE comparisons differ at ~1e-12 by reduction order -- the gl_g0.sbatch
scope line: only WITHIN-MACHINE contrasts are gated, and quoting arithmetic as a STOP
would be mistaking a gate for a goal. SG1/SG4-split/SG6-vs-incumbent are exactly such
cross-host compares with ==0.0 gates, so on this host they would FAIL spuriously.

Pointing REF at a nonexistent path takes the script's own designed SKIP branch (it prints
the skip loudly). What remains GATED is the host-portable subset, which is the item-2
authorisation spec verbatim:
  SG2  freshly-baked closure == runtime-arg call, 0.0 (bit-exact pattern flip, THIS host)
  SG3  warm pattern flip < 10 ms median, jit cache sizes constant (no recompile)
  SG5  grads through the double-where with a carried non-evaluable source
  SG6  flip stability 0.0 + batched-evaluator cache constancy + the one-time None->array
       structure retrace paid HERE (t_first printed), not in the fan
The vs-incumbent bit-level record stands from cronus: reports/FORGE_G2_runtime_smask.md,
ALL PASS, logs in CW_transition/logs/.
"""
import sys
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_b1_smask_gate as SG

SG.REF = SG.REF + ".CRONUS_CAPTURE_NOT_BIT_PORTABLE"
print(f"[H200 wrapper] reference-absent mode: REF -> {SG.REF}", flush=True)
print(f"jax {SG.jax.__version__}, {SG.jax.devices()}", flush=True)
sys.exit(SG.main())
