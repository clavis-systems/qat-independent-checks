"""Test A: windowed sub-circuit (logical gates fully inside sites [a, a+W)),
dense gate-by-gate statevector (unfused) vs fused symmetric MPS without truncation."""
import sys
import time

import numpy as np

import lsh


def dense_window(circ, a, W):
    n = 2 * W
    pos = lambda L: 2 * (lsh.site_chain(L)[0] - a) + lsh.site_chain(L)[1]
    inwin = lambda L: a <= lsh.site_chain(L)[0] < a + W
    psi = np.zeros(2 ** n, dtype=complex)
    idx = 0
    for L in range(circ["nq"]):
        if inwin(L) and circ["bits"][L]:
            idx |= 1 << (n - 1 - pos(L))
    psi[idx] = 1.0
    psi = psi.reshape([2] * n)
    res = []
    prev = 0
    segs = []
    for cut in circ["cuts"]:
        seg = [g for g in circ["logical"][prev:cut] if all(inwin(L) for L in g[2])]
        segs.append(seg)
        for name, par, qs in seg:
            U = lsh.gate_matrix(name, par)
            if len(qs) == 1:
                p = pos(qs[0])
                psi = np.moveaxis(np.tensordot(U, psi, axes=(1, p)), 0, p)
            else:
                p, q = pos(qs[0]), pos(qs[1])
                psi = np.moveaxis(np.tensordot(U.reshape(2, 2, 2, 2), psi, axes=([2, 3], [p, q])), [0, 1], [p, q])
        prob = np.abs(psi) ** 2
        nq = np.array([prob.sum(axis=tuple(i for i in range(n) if i != j))[1] for j in range(n)])
        res.append(nq.reshape(W, 2))
        prev = cut
    return res, segs


def mps_window(circ, a, W, segs, chi):
    states = [0] * W
    for L in range(circ["nq"]):
        r, c = lsh.site_chain(L)
        if a <= r < a + W and circ["bits"][L]:
            states[r - a] |= (2 if c == 0 else 1)
    mps = lsh.SymMPS(states)
    out = []
    stats = {}
    for seg in segs:
        sops = lsh.site_ops_from_logical(seg, lambda u: lsh.site_chain(u)[0] - a, lambda u: lsh.site_chain(u)[1])
        assert lsh.check_conservation(sops) < 1e-12
        lsh.run_layers(mps, lsh.schedule(sops, W), chi, stats)
        n, _ = mps.measure()
        out.append(n)
    return out, mps.discarded


if __name__ == "__main__":
    path = sys.argv[1]
    circ = lsh.load_circuit(path)
    for spec in sys.argv[2:]:
        a, W = map(int, spec.split(":"))
        t0 = time.time()
        dres, segs = dense_window(circ, a, W)
        t1 = time.time()
        mres, disc = mps_window(circ, a, W, segs, chi=4 ** (W // 2 + 1))
        t2 = time.time()
        dev = [float(np.abs(d - m).max()) for d, m in zip(dres, mres)]
        ngates = sum(len(s) for s in segs)
        print("window a=%d W=%d (%d qubits, %d gates): max |n_dense - n_mps| per step: %s" % (a, W, 2 * W, ngates, " ".join("%.1e" % x for x in dev)))
        print("   overall max dev %.2e; discarded %.1e; dense %.1fs mps %.1fs" % (max(dev), disc, t1 - t0, t2 - t1))
        nd = dres[-1]
        print("   step-20 window occupations (dense) first rows:", np.round(nd[:3], 6).tolist(), flush=True)
