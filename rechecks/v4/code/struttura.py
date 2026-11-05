"""V4: struttura della Hamiltoniana 49Q_v1 (codice indipendente, non derivato da C3)."""
import sys, hashlib, collections
import numpy as np

path = sys.argv[1]
print("sha256", hashlib.sha256(open(path, "rb").read()).hexdigest())
a = np.load(path, allow_pickle=False)
n = len(a[0, 0])
terms = []
for lab, c in a:
    z = complex(c)
    assert abs(z.imag) == 0.0
    terms.append((lab, z.real))
print("n qubit", n, "termini", len(terms), "etichette uniche", len({t[0] for t in terms}))
# posizione nella stringa = etichetta del qubit (convenzione irrilevante per il minimo)
w = collections.Counter(sum(ch != "I" for ch in lab) for lab, _ in terms)
print("pesi", dict(w))
kinds = collections.Counter("".join(sorted(ch for ch in lab if ch != "I")) for lab, _ in terms)
print("tipi", dict(kinds))
# grafo dei termini a 2+ qubit
adj = collections.defaultdict(set)
for lab, c in terms:
    q = [i for i, ch in enumerate(lab) if ch != "I"]
    for i in q:
        for j in q:
            if i != j:
                adj[i].add(j)
def comps(adjd, nodes):
    seen, out = set(), []
    for s in nodes:
        if s in seen: continue
        st, comp = [s], []
        seen.add(s)
        while st:
            u = st.pop(); comp.append(u)
            for v in adjd[u]:
                if v not in seen: seen.add(v); st.append(v)
        out.append(sorted(comp))
    return out
print("componenti con tutti i termini:", [len(c) for c in comps(adj, range(n))])
# elenco dei termini a 2 qubit per tipo e |c|
two = [(lab, c) for lab, c in terms if sum(ch != "I" for ch in lab) == 2]
byk = collections.Counter(("".join(lab[i] for i in range(n) if lab[i] != "I"), round(abs(c), 6)) for lab, c in two)
for k, v in sorted(byk.items()): print("  2q", k, v)
# gruppi: componenti del grafo dei soli termini XX/YY
adj2 = collections.defaultdict(set)
for lab, c in two:
    q = [i for i, ch in enumerate(lab) if ch != "I"]
    if {lab[q[0]], lab[q[1]]} <= {"X", "Y"}:
        adj2[q[0]].add(q[1]); adj2[q[1]].add(q[0])
cc = comps(adj2, range(n))
print("componenti XX/YY:", [len(c) for c in cc])
gid = {q: k for k, c in enumerate(cc) for q in c}
# termini XZ: dentro o tra componenti?
inside = between = 0
for lab, c in two:
    q = [i for i, ch in enumerate(lab) if ch != "I"]
    if "Z" in (lab[q[0]], lab[q[1]]):
        if gid[q[0]] == gid[q[1]]: inside += 1
        else: between += 1
print("termini con Z: dentro", inside, "tra componenti", between)
