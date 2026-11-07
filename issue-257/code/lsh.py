"""Own implementation: QASM -> logical gates (SWAPs tracked as a permutation) ->
fused rung-level gates -> U(1)xU(1) block-sparse MPS with exact gates (TEBD-like).

Written for the independent check of tracker issues #254/#257. numpy/scipy only.
"""
import math
import time
from collections import defaultdict

import numpy as np
import scipy.linalg as sla

from qasm_parse import parse_qasm

# ---------------------------------------------------------------- gates
_H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2.0)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_CX = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
_SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)


def gate_matrix(name, params):
    """Matrix of a gate; 2-qubit gates on (a, b) use index 2*bit_a + bit_b."""
    if name == "x":
        return _X
    if name == "h":
        return _H
    if name == "rz":
        t = params[0]
        return np.diag([np.exp(-0.5j * t), np.exp(0.5j * t)])
    if name == "cx":
        return _CX
    if name == "swap":
        return _SWAP
    raise ValueError("gate not supported: %s" % name)


def embed(U, pos, n, d):
    """Embed U (acting on subsystems `pos`, in that order) into n subsystems of dim d."""
    k = len(pos)
    Ut = U.reshape([d] * (2 * k))
    full = np.eye(d ** n, dtype=complex).reshape([d] * (2 * n))
    res = np.tensordot(Ut, full, axes=(list(range(k, 2 * k)), list(pos)))
    rem = [i for i in range(n) if i not in pos]
    order_now = list(pos) + rem
    perm = [order_now.index(i) for i in range(n)] + list(range(n, 2 * n))
    return res.transpose(perm).reshape(d ** n, d ** n)


# ---------------------------------------------------------------- layout
def site_chain(L):
    """Initial wire L -> (site, chain). Deduced from the QASM (see NOTE.md)."""
    k, j = divmod(L, 4)
    return 2 * k + (j % 2), j // 2


def load_circuit(path, mass_angle=0.03, per_step_mass=120):
    """Parse QASM; return initial bits, logical gate list, step cuts, permutations."""
    nq, gates = parse_qasm(path)
    k = 0
    bits = [0] * nq
    while k < len(gates) and gates[k][0] == "x":
        bits[gates[k][2][0]] ^= 1
        k += 1
    perm = list(range(nq))  # wire -> logical qubit (= initial wire label)
    logical = []
    cuts, raw_cuts, perms = [], [], []
    nmass = 0
    for idx in range(k, len(gates)):
        name, par, qs = gates[idx]
        if name == "swap":
            a, b = qs
            perm[a], perm[b] = perm[b], perm[a]
            continue
        logical.append((name, par, tuple(perm[q] for q in qs)))
        if name == "rz" and abs(abs(par[0]) - mass_angle) < 1e-12:
            nmass += 1
            if nmass % per_step_mass == 0:
                cuts.append(len(logical))
                raw_cuts.append(idx + 1)
                perms.append(list(perm))
    return dict(nq=nq, bits=bits, prefix=k, raw=gates, logical=logical, cuts=cuts,
                raw_cuts=raw_cuts, perms=perms, final_perm=list(perm))


# ---------------------------------------------------------------- fusion
def fuse(ops, d):
    """Greedy exact fusion of a sequence of 1- and 2-unit ops into 2-unit blocks.

    ops: list of (units tuple, matrix in that unit order). Returns list of (units, matrix);
    single-unit leftovers (never followed by a 2-unit op) are appended at the end.
    """
    blocks = []
    last = {}
    pending = {}
    for units, U in ops:
        if len(units) == 1:
            u = units[0]
            if u in last:
                bu, BM = blocks[last[u]]
                p = bu.index(u)
                blocks[last[u]][1] = embed(U, (p,), len(bu), d) @ BM
            else:
                pending[u] = U @ pending.get(u, np.eye(d, dtype=complex))
        else:
            a, b = units
            la, lb = last.get(a), last.get(b)
            if la is not None and la == lb and set(blocks[la][0]) == {a, b}:
                bu, BM = blocks[la]
                blocks[la][1] = embed(U, (bu.index(a), bu.index(b)), 2, d) @ BM
            else:
                Pa = pending.pop(a, np.eye(d, dtype=complex))
                Pb = pending.pop(b, np.eye(d, dtype=complex))
                M = U @ np.kron(Pa, Pb)
                blocks.append([(a, b), M])
                last[a] = last[b] = len(blocks) - 1
    for u, P in pending.items():
        blocks.append([(u,), P])
    return [(tuple(b[0]), b[1]) for b in blocks]


CH = [(0, 0), (0, 1), (1, 0), (1, 1)]  # local rung state s = 2*b_chain0 + b_chain1


def site_ops_from_logical(glist, site_of, chain_of):
    """Logical gates -> fused qubit blocks -> rung (site) level ops (4x4 or 16x16)."""
    qops = [(g[2], gate_matrix(g[0], g[1])) for g in glist]
    qblocks = fuse(qops, 2)
    sops = []
    for units, U in qblocks:
        sites = [site_of(u) for u in units]
        r0 = min(sites)
        if max(sites) - r0 > 1:
            raise ValueError("non-local block %s" % (units,))
        if max(sites) == r0:
            pos = tuple(chain_of(u) for u in units)
            sops.append(((r0,), embed(U, pos, 2, 2)))
        else:
            pos = tuple(2 * (site_of(u) - r0) + chain_of(u) for u in units)
            sops.append(((r0, r0 + 1), embed(U, pos, 4, 2)))
    return fuse(sops, 4)


def check_conservation(sops, tol=1e-12):
    worst = 0.0
    for units, U in sops:
        if len(units) == 1:
            ch = np.array([CH[s] for s in range(4)])
        else:
            ch = np.array([(CH[a][0] + CH[b][0], CH[a][1] + CH[b][1]) for a in range(4) for b in range(4)])
        mism = (ch[:, None, :] != ch[None, :, :]).any(-1)
        worst = max(worst, float(np.abs(U[mism]).max()) if mism.any() else 0.0)
        un = np.abs(U.conj().T @ U - np.eye(U.shape[0])).max()
        worst = max(worst, 0.0 if un < tol else un)
    return worst


def schedule(sops, L):
    frontier = [0] * L
    layers = defaultdict(list)
    for units, U in sops:
        sites = units if len(units) == 1 else (units[0], units[1])
        lev = max(frontier[s] for s in sites)
        layers[lev].append((units, U))
        for s in sites:
            frontier[s] = lev + 1
    return [layers[i] for i in sorted(layers)]


# ---------------------------------------------------------------- symmetric MPS
def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


_PAIRS_BY_C = defaultdict(list)
for _s1 in range(4):
    for _s2 in range(4):
        _PAIRS_BY_C[_add(CH[_s1], CH[_s2])].append((_s1, _s2))


def _svd(M):
    try:
        return sla.svd(M, full_matrices=False, lapack_driver="gesdd", check_finite=False)
    except np.linalg.LinAlgError:
        return sla.svd(M, full_matrices=False, lapack_driver="gesvd", check_finite=False)


class SymMPS:
    """MPS on L rung sites (dim 4) with U(1)xU(1) charges (N_chain0, N_chain1) on bonds.

    A[r] is a dict {(qL, s): block (d_qL x d_qR)} with qR = qL + CH[s].
    """

    def __init__(self, states):
        self.L = len(states)
        self.A = []
        q = (0, 0)
        for s in states:
            self.A.append({(q, s): np.ones((1, 1), dtype=complex)})
            q = _add(q, CH[s])
        self.center = 0
        self.discarded = 0.0
        self.max_trunc = 0.0
        self.n_svd = 0

    # -- canonical form moves
    def _move_right(self, r):
        A = self.A[r]
        by_qR = defaultdict(list)
        for (qL, s), B in A.items():
            by_qR[_add(qL, CH[s])].append(((qL, s), B))
        newA, Rs = {}, {}
        for qR, lst in by_qR.items():
            M = np.vstack([B for _, B in lst])
            Q, R = np.linalg.qr(M)
            off = 0
            for key, B in lst:
                newA[key] = np.ascontiguousarray(Q[off:off + B.shape[0]])
                off += B.shape[0]
            Rs[qR] = R
        self.A[r] = newA
        self.A[r + 1] = {(qm, s): Rs[qm] @ B for (qm, s), B in self.A[r + 1].items() if qm in Rs}
        self.center = r + 1

    def _move_left(self, r):
        """center r+1 -> r"""
        A = self.A[r + 1]
        by_qL = defaultdict(list)
        for (qL, s), B in A.items():
            by_qL[qL].append(((qL, s), B))
        newA, Ls = {}, {}
        for qL, lst in by_qL.items():
            M = np.hstack([B for _, B in lst])
            Q, R = np.linalg.qr(M.conj().T)
            Qt = Q.conj().T
            Ls[qL] = R.conj().T
            off = 0
            for key, B in lst:
                newA[key] = np.ascontiguousarray(Qt[:, off:off + B.shape[1]])
                off += B.shape[1]
        self.A[r + 1] = newA
        newA0 = {}
        for (qL, s), B in self.A[r].items():
            qR = _add(qL, CH[s])
            if qR in Ls:
                newA0[(qL, s)] = B @ Ls[qR]
        self.A[r] = newA0
        self.center = r

    def move_to(self, t):
        while self.center < t:
            self._move_right(self.center)
        while self.center > t:
            self._move_left(self.center - 1)

    # -- gates
    def apply1(self, r, U):
        ph = np.diag(U)
        self.A[r] = {(qL, s): B * ph[s] for (qL, s), B in self.A[r].items()}

    def apply2(self, r, G, chi, right, cutoff=1e-13):
        if self.center not in (r, r + 1):
            self.move_to(r if right else r + 1)
        A1, A2 = self.A[r], self.A[r + 1]
        A2_by = defaultdict(dict)
        for (qm, s2), B in A2.items():
            A2_by[qm][s2] = B
        theta = defaultdict(dict)  # (qL, C) -> {(s1,s2): block}
        for (qL, s1), B1 in A1.items():
            qm = _add(qL, CH[s1])
            for s2, B2 in A2_by.get(qm, {}).items():
                C = _add(CH[s1], CH[s2])
                theta[(qL, C)][(s1, s2)] = B1 @ B2
        bym = defaultdict(dict)  # qM -> {(qL, s1, s2): block}
        G4 = G.reshape(4, 4, 4, 4)
        for (qL, C), d in theta.items():
            pairs = _PAIRS_BY_C[C]
            shp = next(iter(d.values())).shape
            stack = np.zeros((len(pairs),) + shp, dtype=complex)
            for i, p in enumerate(pairs):
                if p in d:
                    stack[i] = d[p]
            Gs = np.array([[G4[o[0], o[1], i[0], i[1]] for i in pairs] for o in pairs])
            out = np.tensordot(Gs, stack, axes=(1, 0))
            for i, (s1, s2) in enumerate(pairs):
                if not out[i].any():
                    continue
                bym[_add(qL, CH[s1])][(qL, s1, s2)] = out[i]
        svds = {}
        allS = []
        for qM, d in bym.items():
            rows = sorted({(k[0], k[1]) for k in d})
            cols = sorted({(k[2], _add(qM, CH[k[2]])) for k in d})
            rdim, cdim = {}, {}
            for (qL, s1, s2), T in d.items():
                rdim[(qL, s1)] = T.shape[0]
                cdim[(s2, _add(qM, CH[s2]))] = T.shape[1]
            roff, o = {}, 0
            for k in rows:
                roff[k] = o
                o += rdim[k]
            nr = o
            coff, o = {}, 0
            for k in cols:
                coff[k] = o
                o += cdim[k]
            nc = o
            M = np.zeros((nr, nc), dtype=complex)
            for (qL, s1, s2), T in d.items():
                ro, co = roff[(qL, s1)], coff[(s2, _add(qM, CH[s2]))]
                M[ro:ro + T.shape[0], co:co + T.shape[1]] = T
            U, S, Vh = _svd(M)
            self.n_svd += 1
            svds[qM] = (U, S, Vh, rows, roff, rdim, cols, coff, cdim)
            allS.append(S)
        qlist = list(svds.keys())
        Sall = np.concatenate(allS)
        sec = np.concatenate([np.full(len(svds[q][1]), i) for i, q in enumerate(qlist)])
        tot = float(np.sum(Sall ** 2))
        order = np.argsort(-Sall, kind="stable")
        keep_n = min(chi, int(np.sum(Sall > cutoff * Sall[order[0]])))
        kcount = np.bincount(sec[order[:keep_n]], minlength=len(qlist))
        kept_sum = 0.0
        newA1, newA2 = {}, {}
        for iq, qM in enumerate(qlist):
            U, S, Vh, rows, roff, rdim, cols, coff, cdim = svds[qM]
            k = int(kcount[iq])
            if k <= 0:
                continue
            Sk = S[:k]
            kept_sum += float(np.sum(Sk ** 2))
            for key in rows:
                blk = U[roff[key]:roff[key] + rdim[key], :k]
                newA1[key] = np.ascontiguousarray(blk) if right else blk * Sk[None, :]
            for key in cols:
                s2 = key[0]
                blk = Vh[:k, coff[key]:coff[key] + cdim[key]]
                newA2[(qM, s2)] = (Sk[:, None] * blk) if right else np.ascontiguousarray(blk)
        trunc = max(0.0, (tot - kept_sum) / tot)
        self.discarded += trunc
        self.max_trunc = max(self.max_trunc, trunc)
        norm = math.sqrt(kept_sum)
        if right:
            newA2 = {k: v / norm for k, v in newA2.items()}
        else:
            newA1 = {k: v / norm for k, v in newA1.items()}
        p = Sall ** 2 / tot
        p = p[p > 0]
        self.last_entropy = float(-np.sum(p * np.log(p)))
        self.A[r], self.A[r + 1] = newA1, newA2
        self.center = r + 1 if right else r
        return keep_n

    # -- observables
    def bond_dims(self):
        out = []
        for r in range(self.L - 1):
            dims = {}
            for (qL, s), B in self.A[r].items():
                dims[_add(qL, CH[s])] = B.shape[1]
            out.append(sum(dims.values()))
        return out

    def measure(self):
        """<n_chain0(r)>, <n_chain1(r)> for all sites; leaves center at L-1."""
        self.move_to(0)
        res = np.zeros((self.L, 2))
        norms = []
        for r in range(self.L):
            if r > 0:
                self._move_right(r - 1)
            p = np.zeros(4)
            for (qL, s), B in self.A[r].items():
                p[s] += float(np.vdot(B, B).real)
            nrm = p.sum()
            norms.append(nrm)
            res[r, 0] = (p[2] + p[3]) / nrm
            res[r, 1] = (p[1] + p[3]) / nrm
        return res, float(np.max(np.abs(np.array(norms) - 1.0)))

    def nbytes(self):
        return sum(B.nbytes for A in self.A for B in A.values())


def run_layers(mps, layers, chi, stats):
    for layer in layers:
        two = sorted([op for op in layer if len(op[0]) == 2], key=lambda o: o[0][0])
        for units, U in layer:
            if len(units) == 1:
                mps.apply1(units[0], U)
        if not two:
            continue
        lo, hi = two[0][0][0], two[-1][0][0]
        right = abs(mps.center - lo) <= abs(mps.center - (hi + 1))
        seq = two if right else two[::-1]
        for units, U in seq:
            kept = mps.apply2(units[0], U, chi, right)
            stats["max_chi"] = max(stats.get("max_chi", 0), kept)
            stats["max_S"] = max(stats.get("max_S", 0.0), mps.last_entropy)


def observables(n):
    """n: (L, 2) occupations per site. Returns N (sum of staggered n_f, #254 def), Q."""
    L = n.shape[0]
    nsum = n[:, 0] + n[:, 1]
    nf = np.where(np.arange(L) % 2 == 0, nsum, 2.0 - nsum)
    return float(nf.sum()), float(nsum.sum())
