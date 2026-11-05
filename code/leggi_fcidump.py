"""Lettore FCIDUMP scritto per C2 (non usa il lettore di pyscf né quello di block2).

Restituisce n_orb, n_elec, ms2, h1 (n x n), g2 (n x n x n x n, notazione chimica (pq|rs)), e_core.
Applica la simmetria a 8 vie degli integrali reali.
"""
import re
import sys
import hashlib
import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blocco in iter(lambda: f.read(1 << 20), b""):
            h.update(blocco)
    return h.hexdigest()


def leggi(path):
    with open(path) as f:
        testo = f.read()
    fine = re.search(r"&END|/\s*\n", testo, re.IGNORECASE)
    intest, corpo = testo[: fine.start()], testo[fine.end():]
    norb = int(re.search(r"NORB\s*=\s*(\d+)", intest, re.I).group(1))
    nelec = int(re.search(r"NELEC\s*=\s*(\d+)", intest, re.I).group(1))
    ms2 = int(re.search(r"MS2\s*=\s*(-?\d+)", intest, re.I).group(1))
    n = norb
    h1 = np.zeros((n, n))
    g2 = np.zeros((n, n, n, n))
    ecore = 0.0
    righe = 0
    for riga in corpo.splitlines():
        p = riga.split()
        if len(p) != 5:
            continue
        v = float(p[0].replace("D", "E").replace("d", "e"))
        i, j, k, l = (int(x) for x in p[1:])
        righe += 1
        if i == j == k == l == 0:
            ecore += v
        elif k == 0 and l == 0:
            a, b = i - 1, j - 1
            h1[a, b] = h1[b, a] = v
        else:
            a, b, c, d = i - 1, j - 1, k - 1, l - 1
            for (w, x, y, z) in ((a, b, c, d), (b, a, c, d), (a, b, d, c), (b, a, d, c),
                                 (c, d, a, b), (d, c, a, b), (c, d, b, a), (d, c, b, a)):
                g2[w, x, y, z] = v
    return dict(norb=norb, nelec=nelec, ms2=ms2, h1=h1, g2=g2, ecore=ecore, righe=righe)


if __name__ == "__main__":
    d = leggi(sys.argv[1])
    print("sha256", sha256(sys.argv[1]))
    print("norb", d["norb"], "nelec", d["nelec"], "ms2", d["ms2"], "righe", d["righe"], "ecore", repr(d["ecore"]))
    print("h1 simmetrica:", np.allclose(d["h1"], d["h1"].T))
    g = d["g2"]
    print("g2 simmetria (pq|rs)=(rs|pq):", np.allclose(g, g.transpose(2, 3, 0, 1)))
