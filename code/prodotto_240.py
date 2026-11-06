"""Limite superiore variazionale per la #240 con uno stato prodotto tra i gruppi di qubit (scritto per C2).

Uso: python -I prodotto_240.py HAMILTONIANA.npy OUT.json
1. Gruppi = componenti del grafo dei termini XX/YY (16 + 16 + 16 + 1 qubit).
2. Per ogni gruppo da 16 qubit: Hamiltoniana con tutti i termini interni al gruppo, matrice sparsa 2^16,
   stato fondamentale con eigsh (nessuna ipotesi sul supporto).
3. Stato di prova = prodotto di stati dei gruppi, ottimizzato in modo autoconsistente (ogni gruppo vede i
   termini tra gruppi come campo medio degli altri; si ripete fino a convergenza). L'energia dello stato
   prodotto è esatta da calcolare (termini tra gruppi = prodotto dei valori medi) ed è un limite superiore
   all'energia fondamentale. Primo tentativo (solo stati fondamentali dei gruppi, senza ciclo): -85.1906,
   0,396 sopra; quel log non è incluso. Se coincide con il valore della #240, quel valore è
   raggiunto da uno stato esplicito; il limite inferiore NON è rifatto qui.
"""
import json
import sys
import time
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

t0 = time.time()
arr = np.load(sys.argv[1], allow_pickle=False)
n = len(str(arr[0][0]))
termini = []
for p, c in arr:
    p = str(p)
    c = complex(str(c).replace(" ", "")).real
    ops = {k: ch for k, ch in enumerate(reversed(p)) if ch != "I"}
    termini.append((ops, c))
par = list(range(n))


def f(a):
    while par[a] != a:
        par[a] = par[par[a]]
        a = par[a]
    return a


for ops, c in termini:
    if sorted(ops.values()) in (["X", "X"], ["Y", "Y"]):
        a, b = ops
        par[f(a)] = f(b)
gruppi = defaultdict(list)
for a in range(n):
    gruppi[f(a)].append(a)
gruppi = sorted(gruppi.values(), key=len, reverse=True)
gruppo_di = {q: g for g, qs in enumerate(gruppi) for q in qs}


def matrice(qs, lista):
    loc = {q: i for i, q in enumerate(qs)}
    m = len(qs)
    s = np.arange(2 ** m, dtype=np.int64)
    righe, col, val = [], [], []
    for ops, c in lista:
        x = z = ny = 0
        for q, ch in ops.items():
            if ch in "XY":
                x |= 1 << loc[q]
            if ch in "ZY":
                z |= 1 << loc[q]
            ny += ch == "Y"
        segno = 1 - 2 * (np.bitwise_count(s & z).astype(np.int64) % 2)  # cast: bitwise_count dà uint8
        righe.append(s ^ x)
        col.append(s)
        val.append((c * (1j) ** ny) * segno)
    H = sp.csr_matrix((np.concatenate(val), (np.concatenate(righe), np.concatenate(col))), shape=(2 ** m, 2 ** m))
    return H


def media_singolo(v, m, i, ch):
    s = np.arange(2 ** m, dtype=np.int64)
    if ch == "Z":
        return float(np.real(np.vdot(v, v * (1 - 2 * ((s >> i) & 1)))))
    if ch == "X":
        return float(np.real(np.vdot(v, v[s ^ (1 << i)])))
    raise ValueError(ch)


iso = [qs[0] for qs in gruppi if len(qs) == 1]
grandi = [g for g, qs in enumerate(gruppi) if len(qs) > 1]
interni = {g: [(o, c) for o, c in termini if all(gruppo_di[q] == g for q in o)] for g in grandi}
tra = [(o, c) for o, c in termini if not (len({gruppo_di[q] for q in o}) == 1 and len(gruppi[gruppo_di[next(iter(o))]]) > 1)]
H0 = {g: matrice(gruppi[g], interni[g]) for g in grandi}
medie = {}  # (qubit, lettera) -> valore medio nello stato prodotto corrente
stato_iso = {q: np.array([1.0, 0.0]) for q in iso}
X1 = np.array([[0.0, 1.0], [1.0, 0.0]])
Z1 = np.diag([1.0, -1.0])
for q in iso:
    medie[(q, "Z")], medie[(q, "X")] = 1.0, 0.0
for g in grandi:  # partenza: campo nullo (= stati fondamentali dei gruppi da soli al primo passo)
    for q in gruppi[g]:
        medie[(q, "Z")], medie[(q, "X")] = 0.0, 0.0
vettori = {}
storia = []


def energia_prodotto():
    """<Psi|H|Psi> esatta per lo stato prodotto corrente."""
    E = 0.0
    for g in grandi:
        v = vettori[g]
        E += float(np.real(np.vdot(v, H0[g] @ v)))
    for o, c in tra:
        t = c
        for q, ch in o.items():
            t *= medie[(q, ch)]
        E += t
    return E


def campo(qs_set):
    """Termini tra gruppi visti da un gruppo: coefficiente * medie degli altri fattori."""
    eff = []
    for o, c in tra:
        dentro = {q: ch for q, ch in o.items() if q in qs_set}
        if not dentro:
            continue
        fuori = {q: ch for q, ch in o.items() if q not in qs_set}
        k = c
        for q, ch in fuori.items():
            k *= medie[(q, ch)]
        eff.append((dentro, k))
    return eff


for giro in range(40):
    for g in grandi:
        qs = gruppi[g]
        H = H0[g] + (matrice(qs, campo(set(qs))) if campo(set(qs)) else 0)
        w, V = spla.eigsh(H, k=2, which="SA", tol=1e-13, v0=vettori.get(g))
        v = V[:, np.argmin(w)]
        vettori[g] = v
        for i, q in enumerate(qs):
            for ch in "XZ":
                medie[(q, ch)] = media_singolo(v, len(qs), i, ch)
    for q in iso:
        h = np.zeros((2, 2))
        for dentro, k in campo({q}):
            h += k * (X1 if dentro[q] == "X" else Z1)
        w, V = np.linalg.eigh(h)
        u = V[:, 0]
        stato_iso[q] = u
        medie[(q, "Z")] = float(u @ Z1 @ u)
        medie[(q, "X")] = float(u @ X1 @ u)
    E = energia_prodotto()
    storia.append(E)
    print("giro", giro, "E stato prodotto", repr(E), "tempo", round(time.time() - t0, 1), flush=True)
    if giro > 2 and abs(storia[-1] - storia[-2]) < 1e-13:
        break

rif = -85.586625069864652
E = storia[-1]

# residuo ||H Psi - E Psi|| per lo stato prodotto, senza formare vettori a 49 qubit:
# H Psi = somma di stati prodotto; la norma si ottiene dalle sovrapposizioni gruppo per gruppo.
P1 = {"X": X1, "Z": Z1, "Y": np.array([[0, -1j], [1j, 0]])}
base = {("g", g): vettori[g].astype(complex) for g in grandi}
base.update({("q", q): stato_iso[q].astype(complex) for q in iso})
componenti = []  # (coefficiente, {chiave: vettore modificato})
for g in grandi:
    componenti.append((1.0, {("g", g): H0[g] @ vettori[g]}))
for o, c in tra:
    mod = {}
    per_gruppo = defaultdict(dict)
    for q, ch in o.items():
        if q in stato_iso:
            mod[("q", q)] = P1[ch] @ stato_iso[q]
        else:
            per_gruppo[gruppo_di[q]][q] = ch
    for g, sub in per_gruppo.items():
        mod[("g", g)] = matrice(gruppi[g], [(sub, 1.0)]) @ vettori[g]
    componenti.append((c, mod))
componenti.append((-E, {}))
chiavi = list(base)
norma2 = 0.0
for ci, mi in componenti:
    for cj, mj in componenti:
        sov = 1.0 + 0j
        for k in chiavi:
            a = mi.get(k, base[k])
            b = mj.get(k, base[k])
            sov *= np.vdot(a, b)
        norma2 += ci * cj * sov.real
residuo = float(np.sqrt(max(norma2, 0.0)))
print("residuo ||H Psi - E Psi|| dello stato prodotto:", residuo, "(norma^2 grezza", norma2, ")", flush=True)
print("ENERGIA stato prodotto autoconsistente", repr(E), "valore #240", rif, "differenza", E - rif,
      "tempo", time.time() - t0, flush=True)
json.dump(dict(gruppi=[len(qs) for qs in gruppi], qubit_isolati=iso, termini_tra_gruppi=len(tra),
               storia_energia=storia, E_stato_prodotto=E, residuo_HPsi_meno_EPsi=residuo, valore_issue_240=rif, differenza=E - rif,
               secondi=time.time() - t0), open(sys.argv[2], "w"), indent=1)
