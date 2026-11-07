import sys, collections
from parse import parse
g = parse("/data/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm")
g = [x for x in g if x[0] != "barrier"]
first = next(i for i, x in enumerate(g) if x[0] == "sdg")
last = max(i for i, x in enumerate(g) if x[0] in ("s", "sx") )
print("A=[0,%d) M=[%d,%d] B=(%d,%d)" % (first, first, last, last, len(g)))
A = g[:first]; M = g[first:last + 1]; B = g[last + 1:]
print(len(A), len(M), len(B))
print("M gates:", collections.Counter((x[0], x[1]) for x in M))
def inv(x):
    n, a, q = x
    if n in ("rx", "rz"): return (n, -a, q)
    if n == "cz": return x
    raise ValueError(n)
Adag = [inv(x) for x in reversed(A)]
def canon(L):  # order within commuting gates on different qubits may differ; compare per-qubit sequences
    per = collections.defaultdict(list)
    for n, a, q in L:
        for t in q:
            per[t].append((n, None if a is None else round(a, 12), q))
    return per
def same(L1, L2):
    p1, p2 = canon(L1), canon(L2)
    return p1 == p2
print("B == A (per-qubit sequences)?", same(B, A))
print("B == A^dagger?", same(B, Adag))
print("A qubits", len({t for x in A for t in x[2]}), "B qubits", len({t for x in B for t in x[2]}))
print("A rx angles", collections.Counter(round(x[1], 6) for x in A if x[0] == "rx"))
print("B rx angles", collections.Counter(round(x[1], 6) for x in B if x[0] == "rx"))
# split A into halves by cz count: first 9 layers vs last 9
def halves(L):
    nc = 0; out = []
    for x in L:
        out.append(x)
    return out
# compare A second half with inverse of first half etc.
cz_idx = [i for i, x in enumerate(A) if x[0] == "cz"]
print("A cz count", len(cz_idx))
# find the step boundary: first rx after the 9th-layer's cz block ... use half by cz count 162
cnt = 0
for i, x in enumerate(A):
    if x[0] == "cz":
        cnt += 1
        if cnt == 162:
            cut = i + 1; break
# extend cut to end of that layer's trailing 1q gates? just print first few gates after cut
print("after cut:", A[cut:cut + 6])
A1 = A[:cut]; A2 = A[cut:]
print("A rx in A1:", collections.Counter(round(x[1], 6) for x in A1 if x[0] == "rx"))
print("A rx in A2:", collections.Counter(round(x[1], 6) for x in A2 if x[0] == "rx"))
print("O qubits 52,59,72 in M? M touches", sorted({t for x in M for t in x[2]})[:60])
print("rz(0.3) qubits:", sorted(x[2][0] for x in M if x[0] == "rz"))
