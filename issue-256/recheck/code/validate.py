import numpy as np, random, time, collections, sys
from ole import *
from tnbuild import *
from peep import *

A, M, B = load("/data/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm")
V = extract_V(M, 49)
Pq = sorted(x[2][0] for x in M if x[0] == "rz")
adj = collections.defaultdict(set)
for n, a, q in A:
    if n == "cz": adj[q[0]].add(q[1]); adj[q[1]].add(q[0])

def ball(seed, k, rng):
    s = [seed]; fr = [seed]
    while len(s) < k:
        cand = sorted({y for x in s for y in adj[x]} - set(s))
        if not cand: break
        s.append(rng.choice(cand))
    return s

rng = random.Random(1)
worst = 0
cases = [(52, 8), (52, 10), (11, 9), (26, 10)]
for seed, k in cases:
    sub = ball(seed, k, rng)
    idx = {q: i for i, q in enumerate(sub)}
    nq = len(sub)
    Asub = [(n, a, tuple(idx[t] for t in q)) for n, a, q in A if all(t in idx for t in q)]
    Vm = {idx[t]: V[t] for t in sub if t in V}
    Psub = [idx[t] for t in sub if t in Pq]
    # O: 3 qubits, preferably containing the seed
    Oq = [idx[seed]] + rng.sample([i for i in range(nq) if i != idx[seed]], 2)
    t0 = time.time()
    fd = dense_f(Asub, Vm, Oq, nq)
    td = time.time() - t0
    res = []
    for prune in (False, True):
        G, wires = prepare(Asub, Vm, Oq, Psub, prune=prune)
        tn = make_tn(G, {t: m for t, m in Vm.items() if t in wires}, Oq, wires)
        tn.full_simplify_(seq="ADCRS", atol=1e-12)
        val = f_from_tn(tn, len(wires), optimize="auto")
        res.append(val)
    res.append(f_v2(Asub, Vm, Oq, Psub)); print("  v2 variant", res[-1])
    err = max(abs(r - fd) for r in res)
    worst = max(worst, err)
    print("sub seed %d nq=%d gates=%d cz=%d O=%s dense f=%.15f%+.3gj  TN(noprune)=%.15f TN(prune)=%.15f  maxerr=%.2e  (dense %.1fs)" % (
        seed, nq, len(Asub), sum(x[0] == "cz" for x in Asub), Oq, fd.real, fd.imag, res[0].real, res[1].real, err, td), flush=True)
print("WORST abs error", worst)
