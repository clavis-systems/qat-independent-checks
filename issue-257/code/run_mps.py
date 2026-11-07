"""Run the full 60-rung circuit with the own symmetric MPS at a given chi."""
import argparse
import hashlib
import json
import os
import pickle
import platform
import resource
import sys
import time

import numpy as np
import scipy

import lsh


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def prepare(path, bits_override=None, site_offset=0, sites=None):
    circ = lsh.load_circuit(path)
    bits = circ["bits"] if bits_override is None else bits_override
    nq = circ["nq"]
    L = nq // 2
    states = [0] * L
    for w in range(nq):
        r, c = lsh.site_chain(w)
        if bits[w]:
            states[r] |= (2 if c == 0 else 1)
    segs = []
    prev = 0
    for cut in circ["cuts"]:
        glist = circ["logical"][prev:cut]
        sops = lsh.site_ops_from_logical(glist, lambda u: lsh.site_chain(u)[0], lambda u: lsh.site_chain(u)[1])
        segs.append(sops)
        prev = cut
    assert prev == len(circ["logical"]), "gates after last cut"
    return circ, states, segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qasm", required=True)
    ap.add_argument("--chi", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--steps", type=int, default=20)
    ap.add_argument("--ckpt", default=None)
    args = ap.parse_args()
    t0 = time.time()
    info = dict(qasm=args.qasm, sha256=sha256(args.qasm), chi=args.chi, python=sys.version,
                numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform())
    print(json.dumps(info), flush=True)
    circ, states, segs = prepare(args.qasm)
    L = len(states)
    worst = max(lsh.check_conservation(s) for s in segs)
    n2 = [sum(1 for u, U in s if len(u) == 2) for s in segs]
    print("prep %.1fs; L=%d; two-site ops per step=%s; conservation/unitarity worst=%.2e"
          % (time.time() - t0, L, n2, worst), flush=True)
    assert worst < 1e-12
    # final permutation consistency: output wire 2r,2r+1 <-> site r
    fp = circ["final_perm"]
    ok = all(lsh.site_chain(fp[w])[0] == w // 2 for w in range(circ["nq"]))
    print("final layout: output wires (2r,2r+1) hold logical site r:", ok, flush=True)
    mps = lsh.SymMPS(states)
    n0, _ = mps.measure()
    N0, Q0 = lsh.observables(n0)
    print("step 0 N=%.8f Q=%.8f" % (N0, Q0), flush=True)
    rows = []
    stats = {}
    start = 0
    prior = 0.0
    if args.ckpt and os.path.exists(args.ckpt):
        with open(args.ckpt, "rb") as f:
            ck = pickle.load(f)
        assert ck["chi"] == args.chi and ck["sha256"] == info["sha256"]
        mps.A, mps.center = ck["A"], ck["center"]
        mps.discarded, mps.max_trunc, mps.n_svd = ck["discarded"], ck["max_trunc"], ck["n_svd"]
        rows, start, prior = ck["rows"], ck["step"], ck["elapsed"]
        print("resumed from checkpoint after step %d (prior compute %.1fs)" % (start, prior), flush=True)
    for k in range(start, args.steps):
        ts = time.time()
        layers = lsh.schedule(segs[k], L)
        stats["max_chi"] = 0
        stats["max_S"] = 0.0
        lsh.run_layers(mps, layers, args.chi, stats)
        tm = time.time()
        n, normdev = mps.measure()
        tm = time.time() - tm
        N, Q = lsh.observables(n)
        bd = mps.bond_dims()
        row = dict(step=k + 1, N=N, Q=Q, n=n.tolist(), discarded=mps.discarded, max_trunc=mps.max_trunc,
                   max_bond=max(bd), max_entropy=stats["max_S"], normdev=normdev, meas_time=tm, bond_dims=bd,
                   step_time=time.time() - ts, elapsed=prior + time.time() - t0,
                   peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
                   mps_mb=mps.nbytes() / 2 ** 20)
        rows.append(row)
        print("step %2d N=%.8f Q=%.12f disc=%.3e maxtr=%.2e maxbond=%d Smax=%.3f t=%.1fs tot=%.1fs peakRSS=%.0fMB mps=%.0fMB"
              % (k + 1, N, Q, mps.discarded, mps.max_trunc, max(bd), stats["max_S"], row["step_time"],
                 row["elapsed"], row["peak_rss_mb"], row["mps_mb"]), flush=True)
        with open(args.out, "w") as f:
            json.dump(dict(info=info, rows=rows, N0=N0, Q0=Q0, svds=mps.n_svd), f)
        if args.ckpt:
            with open(args.ckpt + ".tmp", "wb") as f:
                pickle.dump(dict(chi=args.chi, sha256=info["sha256"], A=mps.A, center=mps.center,
                                 discarded=mps.discarded, max_trunc=mps.max_trunc, n_svd=mps.n_svd,
                                 rows=rows, step=k + 1, elapsed=row["elapsed"]), f, protocol=5)
            os.replace(args.ckpt + ".tmp", args.ckpt)
    print("done; this process %.1fs; total compute incl. resumed parts %.1fs" % (time.time() - t0, prior + time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
