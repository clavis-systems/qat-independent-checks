"""Own OLE pipeline.  f = 2^-n Tr(O C^dag O C),  C = U^dag V U  (file = A (U), M (V), B = A^dag).
Equivalent: f = 2^-n Tr(W V^dag W V),  W = U O U^dag."""
import numpy as np, collections
from parse import parse

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1]).astype(complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.diag([1, 1j])
CZ = np.diag([1, 1, 1, -1]).astype(complex).reshape(2, 2, 2, 2)  # [o1,o2,i1,i2]


def mat1(name, a=None):
    if name == "rx": return np.cos(a / 2) * I2 - 1j * np.sin(a / 2) * X
    if name == "rz": return np.cos(a / 2) * I2 - 1j * np.sin(a / 2) * Z
    if name == "s": return S
    if name == "sdg": return S.conj().T
    if name == "sx": return 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
    if name == "sxdg":  # as defined in the file header: s; h; s  (time order) -> matrix S H S
        return S @ H @ S
    raise ValueError(name)


def load(path):
    g = [x for x in parse(path) if x[0] != "barrier"]
    first = next(i for i, x in enumerate(g) if x[0] == "sdg")
    last = max(i for i, x in enumerate(g) if x[0] in ("s", "sx"))
    return g[:first], g[first:last + 1], g[last + 1:]


def dagger_gate(x):
    n, a, q = x
    if n in ("rx", "rz"): return (n, -a, q)
    if n == "cz": return x
    raise ValueError(n)


def dag_list(L):
    return [dagger_gate(x) for x in reversed(L)]


def extract_V(M, nq):
    """Single-qubit content of M: per qubit matrix product in time order; return dict q-> 2x2 (None if identity)."""
    per = {}
    for n, a, q in M:
        t = q[0]
        m = mat1(n, a)
        per[t] = m @ per.get(t, I2)
    return per


def forward_cone(A, Oq):
    """keep gates (time order) that touch the growing support of U O U^dag. exact."""
    sup = set(Oq); keep = []
    for x in A:
        if sup & set(x[2]):
            keep.append(x)
            if len(x[2]) == 2: sup |= set(x[2])
    return keep, sup


def backward_prune(G, Pq):
    """From the V side: drop outermost gates (last in time order) acting only outside supp of V' (exact)."""
    q = set(Pq); keep_rev = []
    for x in reversed(G):
        if q & set(x[2]):
            keep_rev.append(x)
            if len(x[2]) == 2: q |= set(x[2])
    return list(reversed(keep_rev)), q


# ---------------- dense reference ----------------
def dense_apply(state, x, nq):
    """state: array shape (2,)*nq + (batch,) ; apply gate x (qubits already relabelled 0..nq-1)"""
    n, a, q = x
    if n == "cz":
        idx = [slice(None)] * state.ndim
        idx[q[0]] = 1; idx[q[1]] = 1
        state[tuple(idx)] *= -1
        return state
    m = mat1(n, a)
    st = np.tensordot(m, state, axes=([1], [q[0]]))
    return np.moveaxis(st, 0, q[0])


def dense_f(A, Vmats, Oq, nq):
    """f = 2^-n sum_z sigma_z <z|C^dag O C|z>, C = U^dag V U, relabelled qubits 0..nq-1; Vmats: dict q->2x2."""
    d = 2 ** nq
    st = np.eye(d, dtype=complex).reshape((2,) * nq + (d,))  # column z = basis state
    for x in A: st = dense_apply(st, x, nq)
    for t, m in Vmats.items():
        st = np.moveaxis(np.tensordot(m, st, axes=([1], [t])), 0, t)
    for x in dag_list(A): st = dense_apply(st, x, nq)
    # O diag
    bits = np.indices((2,) * nq).reshape(nq, d)
    sig_state = np.ones(d)
    for t in Oq: sig_state *= 1 - 2 * bits[t]
    stf = st.reshape(d, d)  # [out, z]
    # <z|C^dag O C|z> = sum_out sig(out) |C[out,z]|^2 ; sigma_z=sig(z)
    val = (sig_state[:, None] * np.abs(stf) ** 2).sum(axis=0)
    return (sig_state * val).sum() / d
