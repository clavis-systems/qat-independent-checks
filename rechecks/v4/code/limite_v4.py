"""V4: limite inferiore rigoroso per 49Q_v1 (#240), codice scritto da zero.

Idea: per ogni termine tra gruppi  c P_a Q_b  (P, Q Pauli su qubit di gruppi diversi)
  c P_a Q_b = |c| (I + s P_a)(I + t Q_b) - |c| (I + s P_a + t Q_b),   con s*t = segno(c).
Il primo pezzo e' >= 0 (prodotto di due operatori >= 0 che commutano). Quindi
  H >= sum_g H_g(s) + costante, e  min H >= sum_g min H_g(s) + costante.
Scelta dei segni: ricerca ESAUSTIVA su tutte le 2^10 scelte (ogni scelta da' un limite valido).
"""
import sys, json, time, hashlib, itertools, collections
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

path = sys.argv[1]
t0 = time.time()
sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
arr = np.load(path, allow_pickle=False)
n = len(arr[0, 0])
terms = []
for lab, c in arr:
    z = complex(c); assert z.imag == 0.0
    terms.append((lab, z.real))

def support(lab): return [i for i, ch in enumerate(lab) if ch != "I"]

# gruppi = componenti dei termini XX/YY (qubit senza termini XX/YY = gruppo da solo)
par = list(range(n))
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for lab, c in terms:
    q = support(lab)
    if len(q) == 2 and lab[q[0]] in "XY" and lab[q[1]] in "XY":
        par[f(q[0])] = f(q[1])
groups = collections.defaultdict(list)
for i in range(n): groups[f(i)].append(i)
groups = sorted(groups.values(), key=lambda g: g[0])
gid = {q: k for k, g in enumerate(groups) for q in g}

local = [[] for _ in groups]   # termini interni a ogni gruppo
cross = []                    # termini tra gruppi
const = 0.0
for lab, c in terms:
    q = support(lab)
    if not q: const += c; continue
    gs = {gid[i] for i in q}
    if len(gs) == 1: local[gs.pop()].append((lab, c))
    else:
        assert len(q) == 2
        cross.append((lab, c))
print("gruppi", [len(g) for g in groups], "termini tra gruppi", len(cross), "costante", const)

def one_q(i, P):
    return "".join(P if k == i else "I" for k in range(n))

def group_matrix(g, tlist):
    """matrice sparsa reale del gruppo g (qubit in ordine di g; bit k = qubit g[k])."""
    m = len(g); dim = 1 << m; pos = {q: k for k, q in enumerate(g)}
    x = np.arange(dim, dtype=np.int64)
    acc = collections.defaultdict(lambda: np.zeros(dim))
    for lab, c in tlist:
        xm = zm = 0; ny = 0
        for i in support(lab):
            b = 1 << pos[i]; ch = lab[i]
            if ch == "X": xm |= b
            elif ch == "Z": zm |= b
            elif ch == "Y": xm |= b; zm |= b; ny += 1
        assert ny % 2 == 0
        ph = (-1) ** (ny // 2)          # i^ny reale
        sign = 1 - 2 * (np.bitwise_count(x & zm).astype(np.int64) & 1)
        acc[xm] = acc[xm] + c * ph * sign
    rows, cols, vals = [], [], []
    for xm, v in acc.items():
        nz = v != 0
        rows.append((x ^ xm)[nz]); cols.append(x[nz]); vals.append(v[nz])
    M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim))
    return M

def gmin(M, seed=0):
    if M.shape[0] <= 64:
        return float(np.linalg.eigvalsh(M.toarray())[0]), None
    v0 = np.random.default_rng(seed).standard_normal(M.shape[0])
    w, v = sla.eigsh(M, k=3, which="SA", v0=v0, tol=1e-14, ncv=40)
    i = np.argsort(w); return float(w[i[0]]), (w[i], v[:, i[0]])

# per ogni gruppo: quali termini tra gruppi lo toccano e su quale qubit
touch = [[] for _ in groups]
for e, (lab, c) in enumerate(cross):
    qa, qb = support(lab)
    touch[gid[qa]].append((e, qa)); touch[gid[qb]].append((e, qb))

def extra_terms(e, s):
    """termini a un qubit e costante prodotti dal termine e con segno s su qa."""
    lab, c = cross[e]; qa, qb = support(lab)
    t = s * (1 if c > 0 else -1)
    return {qa: (lab[qa], -abs(c) * s), qb: (lab[qb], -abs(c) * t)}, -abs(c)

# tabella: per gruppo, per ogni scelta dei segni dei termini che lo toccano, il minimo
E = len(cross)
cache = []
for k, g in enumerate(groups):
    es = sorted({e for e, _ in touch[k]})
    tab = {}
    for signs in itertools.product((1, -1), repeat=len(es)):
        tl = list(local[k])
        for e, s in zip(es, signs):
            ex, _ = extra_terms(e, s)
            for q, (P, cc) in ex.items():
                if gid[q] == k: tl.append((one_q(q, P), cc))
        tab[signs] = gmin(group_matrix(g, tl))[0]
    cache.append((es, tab))
    print(f"gruppo {k} ({len(g)} qubit): {len(tab)} scelte, {time.time()-t0:.0f} s", flush=True)

best = None
for signs in itertools.product((1, -1), repeat=E):
    val = const - sum(abs(c) for _, c in cross)
    for es, tab in cache: val += tab[tuple(signs[e] for e in es)]
    if best is None or val > best[0]: best = (val, signs)
print("migliore limite (ARPACK)", repr(best[0]), "segni", best[1])
json.dump({"sha256": sha, "groups": groups, "cross": cross, "best_bound": best[0], "signs": best[1]},
          open("v4_best.json", "w"), indent=1)
print("tempo", round(time.time() - t0, 1), "s")
