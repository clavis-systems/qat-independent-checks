"""Controlli di struttura per la #240 (scritto per C2, senza usare il codice di chi ha inviato la richiesta).

Uso: python -I struttura_240.py HAMILTONIANA.npy
Verifica tre affermazioni del README della #240 leggendo solo la lista di Pauli ufficiale:
 1. i termini XX/YY formano un grafo con tre componenti da 16 qubit più un qubit isolato;
 2. su ogni lato |coeff XX| = |coeff YY|;
 3. non ci sono termini ZZ.
"""
import sys
from collections import Counter, defaultdict

import numpy as np

arr = np.load(sys.argv[1], allow_pickle=False)
n = len(str(arr[0][0]))
xx, yy, tipi = {}, {}, Counter()
for p, c in arr:
    p = str(p)
    c = complex(str(c).replace(" ", "")).real
    q = tuple(k for k, ch in enumerate(reversed(p)) if ch != "I")
    t = "".join(ch for ch in reversed(p) if ch != "I")
    tipi[t] += 1
    if t == "XX":
        xx[q] = c
    if t == "YY":
        yy[q] = c
par = list(range(n))


def f(a):
    while par[a] != a:
        par[a] = par[par[a]]
        a = par[a]
    return a


for (a, b) in list(xx) + list(yy):
    par[f(a)] = f(b)
comp = defaultdict(list)
for a in range(n):
    comp[f(a)].append(a)
print("qubit", n, "| lati XX", len(xx), "lati YY", len(yy), "| stessi lati:", set(xx) == set(yy))
print("componenti del grafo XX/YY (dimensioni):", sorted((len(v) for v in comp.values()), reverse=True))
print("max ||XX| - |YY|| sui lati:", max(abs(abs(xx[e]) - abs(yy[e])) for e in xx))
print("lati con segni XX e YY opposti:", sum(1 for e in xx if np.sign(xx[e]) != np.sign(yy[e])))
print("tipi di termini:", dict(tipi), "| termini ZZ:", tipi.get("ZZ", 0))
