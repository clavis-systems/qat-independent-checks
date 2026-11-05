"""Issue #240: which Pauli operators commute with every term of the 49-qubit Hamiltonian?

Written by AI agents (Claude) directed by Emanuele Rizzan, 4 Nov 2026. Apache-2.0.
P = (x', z') commutes with a term (x, z) iff x'.z + z'.x = 0 (mod 2); the script computes the null space
over GF(2) of the 419 x 98 system and prints a basis of it (expected: only X on qubit 48).

Usage:  python -I symmetries_240.py 49Q_v1_hamiltonian.npy
"""
import sys
import numpy as np
arr = np.load(sys.argv[1], allow_pickle=False)
n = 49
rows = []
for p, c in arr:
    p = str(p)
    x = np.zeros(n, np.uint8); z = np.zeros(n, np.uint8)
    for k, ch in enumerate(reversed(p)):
        if ch in "XY": x[k] = 1
        if ch in "ZY": z[k] = 1
    rows.append(np.concatenate([z, x]))  # P=(x',z') commutes iff x'.z + z'.x = 0 -> row [z, x] against unknown [x', z']
A = np.array(rows, np.uint8)
# null space over GF(2)
M = A.copy(); m, N = M.shape; piv = []; r = 0
for col in range(N):
    sel = [i for i in range(r, m) if M[i, col]]
    if not sel: continue
    M[[r, sel[0]]] = M[[sel[0], r]]
    for i in range(m):
        if i != r and M[i, col]: M[i] ^= M[r]
    piv.append(col); r += 1
    if r == m: break
free = [c for c in range(N) if c not in piv]
print("rank", r, "| dimension of the commutant (Pauli operators commuting with every term):", len(free))
for fcol in free:
    v = np.zeros(N, np.uint8); v[fcol] = 1
    for i, pc in enumerate(piv):
        v[pc] = M[i, fcol]
    x, z = v[:n], v[n:]
    s = "".join("Y" if x[k] and z[k] else "X" if x[k] else "Z" if z[k] else "I" for k in range(n))
    print({k: s[k] for k in range(n) if s[k] != "I"})
