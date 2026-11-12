#!/usr/bin/env python
"""GPU worker for the wheel notebook's gate permutation surfaces.

Runs in the jug-gpu env (jax + CUDA, float64); the notebook kernel stays in
discotech and talks to this process over a length-prefixed pickle pipe. All
random draws happen notebook-side -- this worker is a pure deterministic
evaluator, so backends change speed, never numbers. Protocol:
  worker -> {"device": ...} handshake
  notebook -> {"DLG": ..., "TOOTH_W": ...} config
  then per call: notebook -> job dict, worker -> certified bool array (P,).
"""
import sys
import struct
import pickle

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp


def rd():
    hdr = sys.stdin.buffer.read(8)
    if len(hdr) < 8:
        sys.exit(0)
    n = struct.unpack("<Q", hdr)[0]
    buf = sys.stdin.buffer.read(n)
    return pickle.loads(buf)


def wr(obj):
    b = pickle.dumps(obj, protocol=5)
    sys.stdout.buffer.write(struct.pack("<Q", len(b)))
    sys.stdout.buffer.write(b)
    sys.stdout.buffer.flush()


wr({"device": str(jax.devices())})
cfg = rd()
DLG = jnp.asarray(cfg["DLG"])
DLG_SPAN = float(cfg["DLG"][-1] - cfg["DLG"][0])
TOOTH_W = float(cfg["TOOTH_W"])

while True:
    job = rd()
    kap = job["kap"]; rho_u = job["rho_u"]; eps = job["eps"]
    n_use = job["n_use"]; sig_prior = job["sig_prior"]; tshift = job["tshift"]
    pct = job["pct"]; marg = job["marg"]
    P, n = kap.shape
    M1 = eps.shape[1]                                   # 1 + N_SCR_W
    X = DLG[None, None, :]
    S = jnp.zeros((P, M1, DLG.shape[0]))
    for j in range(n):
        S = S - (jnp.asarray(rho_u[:, j])[:, None, None]**2/2) * \
            (1 - jnp.cos(jnp.asarray(kap[:, j])[:, None, None]*X
                         - jnp.asarray(eps[:, :, j])[:, :, None]))
    best = jnp.argmax(S, axis=2)
    smax = S.max(axis=2)
    rival = jnp.where(jnp.abs(DLG[None, None, :] - DLG[best][:, :, None]) > TOOTH_W,
                      S, -jnp.inf).max(axis=2)
    contr = smax - rival
    floor = jnp.percentile(contr[:, 1:], pct, axis=1) + marg + jnp.asarray(tshift)
    post = jnp.exp(S[:, 0, :] - S[:, 0, :].max(axis=1, keepdims=True)
                   - 0.5*(DLG[None, :]/jnp.asarray(sig_prior)[:, None])**2)
    on_tooth = jnp.abs(DLG[None, :] - DLG[best[:, 0]][:, None]) < TOOTH_W
    qmax = (post*on_tooth).sum(axis=1)/post.sum(axis=1)
    if job.get("want") == "stats":
        # raw statistics for rule comparisons (FDR vs FWER on the SAME draws):
        # per-pulsar contrast matrix (truth + nulls), winning tooth, posterior mass
        wr((np.asarray(contr), np.asarray(DLG[best[:, 0]]), np.asarray(qmax)))
        continue
    cert = (jnp.asarray(n_use) >= 3) & (contr[:, 0] > floor) & (qmax > 0.99)
    wrong = cert & (jnp.abs(DLG[best[:, 0]]) > TOOTH_W)   # certified off the true tooth
    wr((np.asarray(cert), np.asarray(wrong)))
