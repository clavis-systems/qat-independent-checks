import numpy as np
from ole import *
A, M, B = load("/data/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm")
V = extract_V(M, 49)
Pq = sorted(x[2][0] for x in M if x[0] == "rz")
print("P =", Pq, len(Pq))
for t in (Pq[0], 7):
    print(t, np.round(V[t], 6).tolist())
# Check V_u = exp(-i 0.15 X) up to phase
v = V[Pq[0]]
ref = np.cos(0.15) * I2 - 1j * np.sin(0.15) * X
ph = v[0, 0] / ref[0, 0]
print("V == exp(-i 0.15 X) up to phase:", np.allclose(v, ph * ref), "phase", ph)
nonP = [t for t in V if t not in Pq]
print("non-P qubits net matrix proportional to identity:", all(np.allclose(V[t], V[t][0, 0] * I2) for t in nonP))
O = [52, 59, 72]
G, sup = forward_cone(A, O)
print("A gates", len(A), "cz", sum(x[0] == "cz" for x in A), "-> forward cone gates", len(G), "cz", sum(x[0] == "cz" for x in G), "support", len(sup))
G2, q = backward_prune(G, Pq)
print("after backward prune from P:", len(G2), "cz", sum(x[0] == "cz" for x in G2), "qubits", len(q))
print("all qubits touched by A:", len({t for x in A for t in x[2]}))
