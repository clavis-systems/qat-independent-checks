"""V3: #240, energia del miglior stato prodotto tra i gruppi XX/YY (3 x 16 qubit + 1), codice scritto da V3.
Ciclo autoconsistente (campo medio tra gruppi) con matrici sparse reali; piu' partenze casuali;
residuo ||(H-E)psi|| calcolato con una decomposizione in pezzi ortogonali (niente cancellazioni);
controllo dell'energia ricalcolata termine per termine dalla lista di Pauli."""
import sys, json, hashlib, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
path = sys.argv[1]
print("sha256", hashlib.sha256(open(path, 'rb').read()).hexdigest())
arr = np.load(path); n = len(arr[0, 0])
terms = []
for p, c in arr:
    c = complex(c); assert abs(c.imag) < 1e-15
    terms.append(({n - 1 - k: ch for k, ch in enumerate(p) if ch != 'I'}, c.real))
const = sum(c for o, c in terms if not o)
# gruppi dai soli termini XX/YY
par = list(range(n))
def f(a):
    while par[a] != a: a = par[a]
    return a
for o, c in terms:
    if len(o) == 2 and set(o.values()) <= {'X', 'Y'}:
        a, b = list(o); par[f(b)] = f(a)
G = {}
for q in range(n): G.setdefault(f(q), []).append(q)
groups = sorted(G.values(), key=len, reverse=True)
gid = {q: i for i, g in enumerate(groups) for q in g}
print("gruppi:", [len(g) for g in groups])
inner = [[] for _ in groups]; cross = []
for o, c in terms:
    if not o: continue
    gs = {gid[q] for q in o}
    if len(gs) == 1: inner[gs.pop()].append((o, c))
    else:
        assert len(o) == 2 and len(gs) == 2; cross.append((o, c))
print("termini interni per gruppo:", [len(t) for t in inner], "termini tra gruppi:", len(cross))

def op_matrix(g, o):
    loc = {q: i for i, q in enumerate(g)}; dim = 2 ** len(g); idx = np.arange(dim)
    xmask = 0; ph = np.ones(dim, dtype=complex)
    for q, ch in o.items():
        i = loc[q]; bit = (idx >> i) & 1
        if ch in 'XY': xmask |= 1 << i
        if ch == 'Z': ph *= 1 - 2 * bit
        if ch == 'Y': ph *= np.where(bit == 0, 1j, -1j)
    return sp.csr_matrix((ph, (idx ^ xmask, idx)), shape=(dim, dim))

Hg = []
for gi, g in enumerate(groups):
    M = sum(c * op_matrix(g, o) for o, c in inner[gi])
    assert abs(M.imag).max() < 1e-14
    M = sp.csr_matrix(M.real); assert abs(M - M.T).max() < 1e-13
    Hg.append(M)
# operatori locali dei termini tra gruppi: (gruppo, matrice) per i due lati
X = []
for o, c in cross:
    sides = []
    for q, ch in o.items():
        m = op_matrix(groups[gid[q]], {q: ch}); assert abs(m.imag).max() == 0
        sides.append((gid[q], sp.csr_matrix(m.real)))
    X.append((c, sides))

def ground(M, v0=None):
    if M.shape[0] <= 64:
        w, v = np.linalg.eigh(M.toarray()); return w[0], v[:, 0]
    w, v = sla.eigsh(M, k=1, which='SA', tol=1e-14, v0=v0, ncv=40)
    return w[0], v[:, 0]

def energy(psi):
    e = const + sum(psi[g] @ (Hg[g] @ psi[g]) for g in range(len(groups)))
    for c, ((ga, A), (gb, B)) in X:
        e += c * (psi[ga] @ (A @ psi[ga])) * (psi[gb] @ (B @ psi[gb]))
    return e

def scf(psi, maxit=40):
    E_old = 0
    for it in range(maxit):
        for g in range(len(groups)):
            heff = Hg[g].copy()
            for c, ((ga, A), (gb, B)) in X:
                if ga == g: heff = heff + c * (psi[gb] @ (B @ psi[gb])) * A
                if gb == g: heff = heff + c * (psi[ga] @ (A @ psi[ga])) * B
            e, v = ground(sp.csr_matrix(heff), psi[g])
            psi[g] = v / np.linalg.norm(v)
        E = energy(psi)
        if abs(E - E_old) < 1e-12 and it > 2: break
        E_old = E
    return E, psi, it

def residual(psi):
    """||(H-E)psi||^2 = somma di pezzi: r_g (deviazione in un gruppo) e f_k (deviazioni in due gruppi)."""
    ng = len(groups); ex = {}
    vecs = []  # ogni vettore: dict gruppo -> vettore locale (gli altri gruppi = psi)
    for g in range(ng):
        heff = Hg[g].copy()
        for c, ((ga, A), (gb, B)) in X:
            if ga == g: heff = heff + c * (psi[gb] @ (B @ psi[gb])) * A
            if gb == g: heff = heff + c * (psi[ga] @ (A @ psi[ga])) * B
        hv = heff @ psi[g]; e = psi[g] @ hv
        vecs.append({g: hv - e * psi[g]})
    for c, ((ga, A), (gb, B)) in X:
        da = A @ psi[ga]; da = da - (psi[ga] @ da) * psi[ga]
        db = B @ psi[gb]; db = db - (psi[gb] @ db) * psi[gb]
        vecs.append({ga: c * da, gb: db})
    tot = 0.0
    for u in vecs:
        for w in vecs:
            if set(u) != set(w): continue  # altrimenti in un gruppo c'e' <psi|delta> = 0
            p = 1.0
            for g in u: p *= u[g] @ w[g]
            tot += p
    return np.sqrt(max(tot, 0)), tot

def energy_from_list(psi):
    """ricalcolo indipendente: <P> = prodotto dei valori attesi locali di ogni lato di ogni termine"""
    e = 0.0
    for o, c in terms:
        if not o: e += c; continue
        byg = {}
        for q, ch in o.items(): byg.setdefault(gid[q], {})[q] = ch
        val = 1.0
        for g, oo in byg.items():
            m = op_matrix(groups[g], oo); val *= (psi[g].conj() @ (m @ psi[g])).real
        e += c * val
    return e

out = {"runs": []}
rng = np.random.default_rng(1)
starts = [("stati fondamentali dei gruppi isolati", [ground(Hg[g])[1] for g in range(len(groups))])]
for s in range(2):
    starts.append((f"casuale {s}", [rng.standard_normal(Hg[g].shape[0]) for g in range(len(groups))]))
best = None
for name, psi in starts:
    psi = [v / np.linalg.norm(v) for v in psi]
    E, psi, it = scf(psi)
    r, r2 = residual(psi)
    print(f"partenza {name}: E = {E:.14f} dopo {it+1} giri; residuo ||(H-E)psi|| = {r:.2e}", flush=True)
    out["runs"].append({"partenza": name, "E": E, "giri": it + 1, "residuo": r})
    if best is None or E < best[0]: best = (E, psi)
E, psi = best
E2 = energy_from_list(psi)
print(f"MIGLIORE: {E:.14f}; ricalcolo termine per termine: {E2:.14f}; diff {E-E2:.1e}")
for ref, nome in [(-85.58662506986461, "C2"), (-85.586625069864652, "issue #240")]:
    print(f"differenza da {nome}: {E-ref:.2e}")
gs_iso = [ground(Hg[g])[0] for g in range(len(groups))]
print("somma dei minimi dei gruppi isolati (senza termini tra gruppi):", const + sum(gs_iso))
out.update({"E_best": E, "E_ricalcolo": E2, "diff_C2": E + 85.58662506986461, "diff_issue": E + 85.586625069864652,
            "gruppi": [len(g) for g in groups], "n_cross": len(cross)})
json.dump(out, open(sys.argv[2], "w"), indent=1)
