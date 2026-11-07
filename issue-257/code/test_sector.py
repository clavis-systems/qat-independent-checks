"""Test B: full 120-qubit raw circuit (SWAPs applied physically on wires), exact in a
few-particle sector, vs the MPS pipeline (SWAPs tracked as permutation, rung sites).
Initial states are custom few-particle bitstrings (not the benchmark states)."""
import itertools
import sys
import time

import numpy as np

import lsh
import run_mps


class Sector:
    def __init__(self, nq, wires):
        k = len(wires)
        self.k = k
        self.occ = np.array(list(itertools.combinations(range(nq), k)), dtype=np.int64)
        self.base = np.array([128 ** (k - 1 - j) for j in range(k)], dtype=np.int64)
        self.keys = self.occ @ self.base
        assert np.all(np.diff(self.keys) > 0)
        self.psi = np.zeros(len(self.occ), dtype=complex)
        self.psi[np.searchsorted(self.keys, np.array(sorted(wires)) @ self.base)] = 1.0
        self.has = [np.nonzero((self.occ == w).any(1))[0] for w in range(nq)]
        self.nst = len(self.occ)

    def apply(self, wires, U):
        if len(wires) == 1:
            assert abs(U[0, 1]) < 1e-12 and abs(U[1, 0]) < 1e-12
            self.psi *= U[0, 0]
            self.psi[self.has[wires[0]]] *= U[1, 1] / U[0, 0]
            return
        a, b = wires
        pc = np.array([0, 1, 1, 2])
        assert np.abs(U[pc[:, None] != pc[None, :]]).max() < 1e-12
        ia, ib = self.has[a], self.has[b]
        mb = np.zeros(self.nst, bool)
        mb[ib] = True
        ma = np.zeros(self.nst, bool)
        ma[ia] = True
        both = ia[mb[ia]]
        only_a = ia[~mb[ia]]
        only_b = ib[~ma[ib]]
        occ = self.occ[only_a].copy()
        occ[occ == a] = b
        occ.sort(axis=1)
        pa = np.searchsorted(self.keys, occ @ self.base)
        assert len(pa) == len(only_b) and np.array_equal(np.sort(pa), np.sort(only_b))
        new = self.psi * U[0, 0]
        new[both] = U[3, 3] * self.psi[both]
        new[only_a] = U[2, 2] * self.psi[only_a] + U[2, 1] * self.psi[pa]
        new[pa] = U[1, 1] * self.psi[pa] + U[1, 2] * self.psi[only_a]
        self.psi = new

    def n(self, w):
        return float(np.sum(np.abs(self.psi[self.has[w]]) ** 2))


def fuse_raw(seg):
    """Split the raw gate list into maximal CONTIGUOUS runs acting on at most two wires
    and multiply each run exactly (no reordering at all). In these QASM files every
    hop / electric / mass pattern is contiguous, so each run conserves particle number
    (checked in Sector.apply). Independent of lsh.fuse."""
    runs = []
    cur_w, cur = None, []
    for g in seg:
        ws = set(g[2])
        if cur and len(cur_w | ws) <= 2:
            cur_w |= ws
            cur.append(g)
        else:
            if cur:
                runs.append((cur_w, cur))
            cur_w, cur = set(ws), [g]
    if cur:
        runs.append((cur_w, cur))
    out = []
    for ws, gl in runs:
        wl = sorted(ws)
        M = np.eye(2 ** len(wl), dtype=complex)
        for name, par, qs in gl:
            U = lsh.gate_matrix(name, par)
            M = lsh.embed(U, tuple(wl.index(q) for q in qs), len(wl), 2) @ M
        out.append((tuple(wl), M))
    return out


def run(path, wires, chi=256):
    circ = lsh.load_circuit(path)
    nq = circ["nq"]
    bits = [0] * nq
    for w in wires:
        bits[w] = 1
    t0 = time.time()
    sec = Sector(nq, wires)
    raw = circ["raw"]
    start = circ["prefix"]
    sec_res = []
    bounds = list(circ["raw_cuts"])
    bounds[-1] = len(raw)  # include trailing gates (if any) in the last step
    for end in bounds:
        seg = raw[start:end]
        blocks = fuse_raw(seg)
        for wq, U in blocks:
            sec.apply(wq, U)
        sec_res.append(np.array([sec.n(w) for w in range(nq)]))
        start = end
    t1 = time.time()
    _, states, segs = run_mps.prepare(path, bits_override=bits)
    mps = lsh.SymMPS(states)
    stats = {}
    devs = []
    perms = list(circ["perms"])
    perms[-1] = circ["final_perm"]
    for k, sops in enumerate(segs):
        lsh.run_layers(mps, lsh.schedule(sops, len(states)), chi, stats)
        n, _ = mps.measure()
        nw = np.array([n[lsh.site_chain(perms[k][w])[0], lsh.site_chain(perms[k][w])[1]] for w in range(nq)])
        devs.append(float(np.abs(nw - sec_res[k]).max()))
    t2 = time.time()
    print("sector k=%d wires=%s dim=%d: max |n_sector - n_mps| per step: %s" % (len(wires), wires, sec.nst, " ".join("%.1e" % d for d in devs)))
    print("   overall %.2e; norm %.15f; discarded %.1e; sector %.1fs mps %.1fs; final total n=%.12f"
          % (max(devs), float(np.sum(np.abs(sec.psi) ** 2)), mps.discarded, t1 - t0, t2 - t1, float(sec_res[-1].sum())), flush=True)
    return max(devs)


if __name__ == "__main__":
    path = sys.argv[1]
    for spec in sys.argv[2:]:
        run(path, [int(x) for x in spec.split(",")])
