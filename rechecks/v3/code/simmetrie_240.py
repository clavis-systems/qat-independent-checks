"""V3: #240, simmetrie di Pauli (commutante) e termini tra gruppi. Algebra lineare su GF(2), codice di V3."""
import sys, numpy as np, itertools
arr = np.load(sys.argv[1]); n = 49
terms = []
for p, c in arr:
    o = {n - 1 - k: ch for k, ch in enumerate(p) if ch != 'I'}
    if o: terms.append((o, complex(c).real))
def vec(o):
    v = np.zeros(2 * n, dtype=np.uint8)
    for q, ch in o.items():
        if ch in 'XY': v[q] = 1
        if ch in 'ZY': v[n + q] = 1
    return v
T = np.array([vec(o) for o, c in terms])
# P commuta con t  <=>  x_P.z_t + z_P.x_t = 0 (mod 2): matrice A con righe [z_t, x_t]
A = np.concatenate([T[:, n:], T[:, :n]], axis=1) % 2
def nullspace_gf2(A):
    A = A.copy() % 2; m, k = A.shape; piv = []; r = 0
    for c in range(k):
        rows = [i for i in range(r, m) if A[i, c]]
        if not rows: continue
        A[[r, rows[0]]] = A[[rows[0], r]]
        for i in range(m):
            if i != r and A[i, c]: A[i] ^= A[r]
        piv.append(c); r += 1
        if r == m: break
    free = [c for c in range(k) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(k, dtype=np.uint8); v[f] = 1
        for i, c in enumerate(piv):
            if A[i, f]: v[c] = 1
        basis.append(v)
    return basis
B = nullspace_gf2(A)
def name(v):
    s = []
    for q in range(n):
        x, z = v[q], v[n + q]
        if x or z: s.append(f"{'X' if x and not z else 'Z' if z and not x else 'Y'}{q}")
    return ' '.join(s)
print("dimensione del commutante (stringhe di Pauli conservate indipendenti):", len(B))
for v in B: print("  ", name(v))
