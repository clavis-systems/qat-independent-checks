"""V4: tutti gli operatori di Pauli che commutano con tutti i termini (nucleo su GF(2) della forma simplettica)."""
import sys, numpy as np
a = np.load(sys.argv[1], allow_pickle=False); n = len(a[0, 0])
rows = []
for lab, _ in a:
    x = [1 if ch in "XY" else 0 for ch in lab]; z = [1 if ch in "ZY" else 0 for ch in lab]
    rows.append(z + x)   # P commuta con T  <=>  x_T.z_P + z_T.x_P = 0  -> riga (z_T, x_T) contro (x_P, z_P)
A = np.array(rows, dtype=np.uint8) % 2
# eliminazione di Gauss su GF(2), nucleo di A (vettori v=(x_P, z_P))
M = A.copy(); m, k = M.shape; piv = []; r = 0
for c in range(k):
    p = next((i for i in range(r, m) if M[i, c]), None)
    if p is None: continue
    M[[r, p]] = M[[p, r]]
    for i in range(m):
        if i != r and M[i, c]: M[i] ^= M[r]
    piv.append(c); r += 1
free = [c for c in range(k) if c not in piv]
print("rango", r, "dimensione del nucleo", len(free))
for f in free:
    v = np.zeros(k, dtype=np.uint8); v[f] = 1
    for i, c in enumerate(piv): v[c] = M[i, f]
    s = "".join({(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}[(v[j], v[n + j])] for j in range(n))
    print("generatore:", s)
