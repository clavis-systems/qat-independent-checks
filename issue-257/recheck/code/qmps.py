"""Qubit-level MPS (one site per qubit) for nearest-neighbour circuits with gates
cx,h,rz,x,swap.  SWAP gates are applied physically.  Own code, independent of the
rung-grouped method."""
import re, time, numpy as np

# ---------------- circuit I/O ----------------
def parse_qasm(path):
    g = []
    for l in open(path):
        l = l.strip().rstrip(';')
        if not l or l.startswith(('OPENQASM', 'include', 'qreg')): continue
        m = re.match(r'(\w+)(?:\(([^)]*)\))?\s+(.*)', l)
        n, a, q = m.groups()
        q = [int(x) for x in re.findall(r'q\[(\d+)\]', q)]
        g.append((n, float(a) if a is not None else None, tuple(q)))
    return g

def gen_circuit(N, theta, steps, mass=0.03, hop=0.15):
    """Gate list of my reconstruction of the LSH circuit for N sites (2N qubits, N even)."""
    g = []
    phi = theta / 2.0
    def hopb(a, b):
        return [('cx', None, (b, a)), ('h', None, (b,)), ('rz', -hop, (b,)), ('cx', None, (a, b)),
                ('rz', hop, (b,)), ('cx', None, (a, b)), ('h', None, (b,)), ('cx', None, (b, a))]
    nq = 2 * N
    for _ in range(steps):
        for k in range(N): g += hopb(2 * k, 2 * k + 1)
        for j in range(N - 1): g.append(('swap', None, (2 * j + 1, 2 * j + 2)))
        for k in range(1, N - 1): g += hopb(2 * k, 2 * k + 1)
        for m in range(N // 2 - 1): g.append(('swap', None, (4 * m + 3, 4 * m + 4)))
        for k in range(N):
            a, b = 2 * k, 2 * k + 1
            g += [('x', None, (a,)), ('rz', -phi, (a,)), ('cx', None, (a, b)), ('rz', phi, (b,)),
                  ('cx', None, (a, b)), ('rz', -phi, (b,)), ('x', None, (a,))]
        for par in (0, 1):
            for k in range(par, N, 2):
                s = -mass if par == 0 else mass
                g += [('rz', s, (2 * k,)), ('rz', s, (2 * k + 1,))]
        for m in range(N // 2): g.append(('swap', None, (4 * m + 1, 4 * m + 2)))
    return g

def init_bits(N, meson_flips=()):
    """grouped layout: q[4m..4m+3] = (i_2m,i_2m+1,o_2m,o_2m+1); odd q occupied"""
    b = np.zeros(2 * N, int)
    b[1::2] = 1
    return b

def qubit_labels(N):
    lab = []
    for m in range(N // 2): lab += [('i', 2 * m), ('i', 2 * m + 1), ('o', 2 * m), ('o', 2 * m + 1)]
    return lab

# ---------------- gate matrices ----------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
H = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
def RZ(a): return np.diag([np.exp(-0.5j * a), np.exp(0.5j * a)])
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)
CX01 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], complex)  # control = first
def mat1(n, a):
    return {'x': X, 'h': H}[n] if n != 'rz' else RZ(a)
def gate_on_pair(n, a, q, lo):
    """4x4 matrix of gate on qubits (lo, lo+1)."""
    if len(q) == 1:
        m = mat1(n, a)
        return np.kron(m, I2) if q[0] == lo else np.kron(I2, m)
    if n == 'swap': return SWAP
    assert n == 'cx'
    if q[0] == lo: return CX01
    return SWAP @ CX01 @ SWAP

# ---------------- fusing ----------------
def fuse(gates):
    """Merge consecutive gates into 2-qubit ops (lo, M4) or 1-qubit ops (q, M2).
    Order-preserving: a gate is merged into the LAST op touching its qubits if that op
    covers all of them."""
    ops = []  # [kind, lo or q, matrix]
    last = {}  # qubit -> op index
    for n, a, q in gates:
        for qq in q: assert True
        j = max(last.get(qq, -1) for qq in q)
        if len(q) == 2: assert abs(q[0] - q[1]) == 1, q
        if j >= 0:
            op = ops[j]
            cov = {op[1]} if op[0] == 1 else {op[1], op[1] + 1}
            if set(q) <= cov:
                if op[0] == 1:
                    op[2] = mat1(n, a) @ op[2]
                else:
                    op[2] = gate_on_pair(n, a, q, op[1]) @ op[2]
                continue
        if len(q) == 1:
            # attach to existing op on that qubit? only if j covers; else new 1q op
            ops.append([1, q[0], mat1(n, a)]); last[q[0]] = len(ops) - 1
        else:
            lo = min(q)
            # absorb earlier 1q/2q ops? keep simple: new op, but absorb pending 1q op on same qubit
            M = gate_on_pair(n, a, q, lo)
            for qq in q:
                jj = last.get(qq, -1)
                if jj >= 0 and ops[jj][0] == 1 and jj == max(last.get(x, -1) for x in q if x == qq) and ops[jj][1] == qq:
                    # 1q op on qq is latest on qq; it can be absorbed iff no later op touches the other qubit before it?
                    # it is later-than-nothing on qq; but ordering vs other qubit irrelevant (different qubits commute)
                    P = ops[jj][2]
                    Pk = np.kron(P, I2) if qq == lo else np.kron(I2, P)
                    M = M @ Pk
                    ops[jj] = [0, 0, None]  # dead
            ops.append([2, lo, M])
            last[lo] = last[lo + 1] = len(ops) - 1
    return [o for o in ops if o[0] != 0]

# ---------------- dense reference ----------------
def dense_run(nq, bits, gates, obs_every=None):
    """statevector; qubit 0 = most significant axis."""
    psi = np.zeros([2] * nq, complex)
    psi[tuple(bits)] = 1.0
    for n, a, q in gates:
        if len(q) == 1:
            psi = np.moveaxis(np.tensordot(mat1(n, a), psi, axes=([1], [q[0]])), 0, q[0])
        else:
            if n == 'swap': M = SWAP
            else: M = CX01
            M = M.reshape(2, 2, 2, 2)
            psi = np.tensordot(M, psi, axes=([2, 3], [q[0], q[1]]))
            psi = np.moveaxis(psi, [0, 1], [q[0], q[1]])
    return psi

def dense_z(psi):
    nq = psi.ndim
    p = (np.abs(psi) ** 2)
    out = []
    for k in range(nq):
        pk = p.sum(axis=tuple(a for a in range(nq) if a != k))
        out.append(pk[0] - pk[1])
    return np.array(out)

# ---------------- MPS ----------------
class MPS:
    def __init__(self, bits, chi, cutoff=1e-14):
        self.n = len(bits)
        self.A = []
        for b in bits:
            t = np.zeros((1, 2, 1), complex); t[0, b, 0] = 1; self.A.append(t)
        self.c = 0
        self.chi = chi; self.cutoff = cutoff
        self.disc = 0.0
        self.maxbond = 1

    def _move(self, to):
        while self.c < to:
            k = self.c; a = self.A[k]; l, d, r = a.shape
            Q, R = np.linalg.qr(a.reshape(l * d, r))
            self.A[k] = Q.reshape(l, d, -1)
            self.A[k + 1] = np.tensordot(R, self.A[k + 1], axes=([1], [0]))
            self.c += 1
        while self.c > to:
            k = self.c; a = self.A[k]; l, d, r = a.shape
            Q, R = np.linalg.qr(a.reshape(l, d * r).T)
            self.A[k] = Q.T.reshape(-1, d, r)
            self.A[k - 1] = np.tensordot(self.A[k - 1], R.T, axes=([2], [0]))
            self.c -= 1

    def apply1(self, q, M):
        self.A[q] = np.tensordot(M, self.A[q], axes=([1], [1])).transpose(1, 0, 2)

    def apply2(self, lo, M):
        self._move(lo)
        a, b = self.A[lo], self.A[lo + 1]
        l = a.shape[0]; r = b.shape[2]
        th = np.tensordot(a, b, axes=([2], [0]))           # l,s,t,r
        th = np.tensordot(M.reshape(2, 2, 2, 2), th, axes=([2, 3], [1, 2]))  # s',t',l,r
        th = th.transpose(2, 0, 1, 3).reshape(l * 2, 2 * r)
        U, S, Vh = np.linalg.svd(th, full_matrices=False)
        w = S ** 2; tot = w.sum()
        # keep minimal number so that discarded weight <= cutoff*tot, capped at chi
        cum = np.cumsum(w[::-1])[::-1]     # cum[k] = weight of S[k:]
        k = len(S)
        ok = np.nonzero(cum <= self.cutoff * tot)[0]
        if len(ok): k = max(1, ok[0])
        k = min(k, self.chi)
        d = cum[k] if k < len(S) else 0.0
        self.disc += d / tot
        S = S[:k]; S = S / np.sqrt((S ** 2).sum())
        self.A[lo] = U[:, :k].reshape(l, 2, k)
        self.A[lo + 1] = (S[:, None] * Vh[:k]).reshape(k, 2, r)
        self.c = lo + 1
        self.maxbond = max(self.maxbond, k)

    def zexp(self):
        self._move(0)
        out = np.zeros(self.n)
        for p in range(self.n):
            self._move(p)
            a = self.A[p]
            w = (np.abs(a) ** 2).sum(axis=(0, 2))
            out[p] = w[0] - w[1]
        return out, None

def run(gates, ops_fn, bits, chi, step_ends, cutoff=1e-14, log=None):
    """gates already fused as ops list; step_ends: op indices (exclusive end) per step."""
    m = MPS(bits, chi, cutoff)
    res = []
    prev = 0
    t0 = time.time()
    for s, e in enumerate(step_ends):
        for op in ops_fn[prev:e]:
            if op[0] == 1: m.apply1(op[1], op[2])
            else: m.apply2(op[1], op[2])
        prev = e
        z, _ = m.zexp()
        res.append((z.copy(), m.disc, m.maxbond))
        if log: log(f'step {s+1} disc_cum={m.disc:.3e} maxbond={m.maxbond} t={time.time()-t0:.0f}s')
    return res

def perm_track(nq, gates):
    """position -> original qubit index whose content is there"""
    pos = list(range(nq))
    for n, a, q in gates:
        if n == 'swap': pos[q[0]], pos[q[1]] = pos[q[1]], pos[q[0]]
    return pos

def N_from_z(z, labels, perm):
    """z: Z expectation by position; perm[p]= initial qubit index now at position p.
    labels: initial qubit index -> (type, site).  N = 60 - 1/2 sum (-1)^r Z."""
    N = len(z) // 2
    s = 0.0
    for p, zz in enumerate(z):
        r = labels[perm[p]][1]
        s += (-1) ** r * zz
    return N - 0.5 * s, N - 0.5 * s  # placeholder second value

def Nobs(z, labels, perm):
    N = len(z) // 2
    s = sum((-1) ** labels[perm[p]][1] * z[p] for p in range(len(z)))
    return N - 0.5 * s

def Qobs(z):
    return len(z) // 2 - 0.5 * z.sum()

def step_boundaries(gates, ops_index_of_gate):
    pass
