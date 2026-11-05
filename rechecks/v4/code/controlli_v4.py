"""V4: controlli sul limite con i segni migliori (da v4_best.json).
1) identita' H = sum_g H_g(s) + C + sum_e B_e, ricostruita sui coefficienti di Pauli
2) ogni B_e >= 0 (autovalori della matrice 4x4)
3) minimi dei gruppi con due metodi: ARPACK (Lanczos) e fattorizzazione LDL^T con inerzia (Sylvester)
4) limite superiore: energia vera dello stato prodotto dei minimi dei gruppi
"""
import sys, os, json, time, collections, runpy, functools
print = functools.partial(print, flush=True)
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, scipy.linalg as la

path = sys.argv[1]
sys.argv = [sys.argv[0], path]
# riuso le funzioni del mio script (stesso file, stessa sessione V4), senza rifare la ricerca
src = open(__file__.replace("controlli_v4.py", "limite_v4.py")).read().split("# tabella: per gruppo")[0]
ns = {"__file__": "limite_v4.py"}; exec(src, ns)
groups, cross, local, terms, n, gid = ns["groups"], ns["cross"], ns["local"], ns["terms"], ns["n"], ns["gid"]
support, one_q, group_matrix, extra_terms = ns["support"], ns["one_q"], ns["group_matrix"], ns["extra_terms"]
best = json.load(open("v4_best.json")); signs = best["signs"]
t0 = time.time()

# --- 1) identita' sui coefficienti
P1 = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0.]]), "Z": np.diag([1., -1])}
rec = collections.defaultdict(float)
glist = [list(local[k]) for k in range(len(groups))]
C = 0.0
for e, s in enumerate(signs):
    ex, cst = extra_terms(e, s); C += cst
    for q, (P, cc) in ex.items(): glist[gid[q]].append((one_q(q, P), cc))
for k in range(len(groups)):
    for lab, c in glist[k]: rec[lab] += c
rec["I" * n] += C
for e, s in enumerate(signs):
    lab, c = cross[e]; qa, qb = support(lab); t = s * (1 if c > 0 else -1)
    # |c|(I + s Pa)(I + t Qb) = |c|(I + s Pa + t Qb + s t Pa Qb)
    rec["I" * n] += abs(c); rec[one_q(qa, lab[qa])] += abs(c) * s; rec[one_q(qb, lab[qb])] += abs(c) * t
    rec[lab] += abs(c) * s * t
orig = collections.defaultdict(float)
for lab, c in terms: orig[lab] += c
keys = set(rec) | set(orig)
err = max(abs(rec[k] - orig[k]) for k in keys)
print("1) identita': errore massimo sui coefficienti", err)

# --- 2) B_e >= 0
mins = []
for e, s in enumerate(signs):
    lab, c = cross[e]; qa, qb = support(lab); t = s * (1 if c > 0 else -1)
    B = abs(c) * np.kron(np.eye(2) + s * P1[lab[qa]], np.eye(2) + t * P1[lab[qb]])
    mins.append(np.linalg.eigvalsh(B)[0])
print("2) autovalore minimo dei pezzi B: min", min(mins), "max", max(mins))

# --- 3) minimi dei gruppi con due metodi
tot_arp, tot_ldl, states = C, C, []
for k, g in enumerate(groups):
    M = group_matrix(g, glist[k])
    if M.shape[0] <= 64:
        w, v = np.linalg.eigh(M.toarray()); lam, psi = w[0], v[:, 0]
        tot_arp += lam; tot_ldl += lam; states.append(psi)
        print(f"   gruppo {k}: denso {lam!r}"); continue
    lams = []
    for seed in (1, 2, 3):
        w, v = sla.eigsh(M, k=4, which="SA", v0=np.random.default_rng(seed).standard_normal(M.shape[0]), tol=1e-14, ncv=40)
        i = np.argsort(w); lams.append(w[i[0]]); psi = v[:, i[0]]; gap = w[i[1]] - w[i[0]]
    lam = min(lams); res = np.linalg.norm(M @ psi - lams[-1] * psi)
    # inerzia: M - (lam - eps) I deve essere definita positiva (Sylvester): LU senza scambi di righe
    out = {1e-9: "non fatto", -1e-9: "non fatto"}
    for eps in ((1e-9, -1e-9) if os.environ.get("V4_INERZIA") == str(k) else ()):
        A = (M - (lam - eps) * sp.identity(M.shape[0])).tocsc()
        lu = sla.splu(A, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                      options={"SymmetricMode": True})
        same = np.array_equal(lu.perm_r, lu.perm_c)
        d = lu.U.diagonal(); out[eps] = (int((d < 0).sum()), same)
    print(f"   gruppo {k}: ARPACK 3 partenze {[repr(x) for x in lams]}, residuo {res:.1e}, salto {gap:.4f}, "
          f"pivot negativi con shift lam-1e-9: {out[1e-9]}, con lam+1e-9: {out[-1e-9]}  (n.neg, perm simmetrica)")
    tot_arp += lam; states.append(psi)
print("3) limite inferiore", repr(tot_arp))

# --- 4) limite superiore: <H> sullo stato prodotto dei minimi dei gruppi
E = 0.0
expv = {}
def exp1(k, q, P):
    g = groups[k]; pos = g.index(q); psi = states[k]; x = np.arange(len(psi))
    if P == "Z": return float(np.sum(psi**2 * (1 - 2 * ((x >> pos) & 1))))
    if P == "X": return float(psi @ psi[x ^ (1 << pos)])
for k in range(len(groups)):
    M = group_matrix(groups[k], local[k]); psi = states[k]; E += float(psi @ (M @ psi))
for lab, c in cross:
    qa, qb = support(lab); E += c * exp1(gid[qa], qa, lab[qa]) * exp1(gid[qb], qb, lab[qb])
print("4) energia dello stato prodotto (limite superiore)", repr(E), " differenza con il limite", E - tot_arp)
print("tempo", round(time.time() - t0, 1), "s")
