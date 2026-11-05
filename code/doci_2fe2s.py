"""DOCI (spazio a seniority zero) per [2Fe-2S] (30e, 20o), scritto per C2: controllo del primo numero della #248.

Uso: python -I doci_2fe2s.py FCIDUMP OUT.json

Nello spazio a seniority zero ogni orbitale è vuoto o doppiamente occupato: con 15 coppie in 20 orbitali
ci sono C(20,15) = 15.504 configurazioni.
  diagonale:      E(S) = 2 sum_{i in S} h_ii + sum_{i,j in S} (2 (ii|jj) - (ij|ji))
  fuori diagonale: spostare una coppia da i (occupato) ad a (vuoto) vale (ia|ia)
Controllo interno: la stessa matrice costruita in modo denso con un secondo codice (somma esplicita sulle
coppie) su un sottoproblema piccolo (8 orbitali, 4 coppie) contro la diagonalizzazione completa.
"""
import itertools
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from leggi_fcidump import leggi, sha256  # noqa: E402


def matrice_doci(h1, g2, npairs, ecore=0.0):
    n = h1.shape[0]
    conf = list(itertools.combinations(range(n), npairs))
    pos = {c: k for k, c in enumerate(conf)}
    J = np.einsum("iijj->ij", g2)
    K = np.einsum("ijji->ij", g2)
    P = np.einsum("iaia->ia", g2)  # (ia|ia)
    righe, col, val = [], [], []
    diag = np.empty(len(conf))
    for k, c in enumerate(conf):
        s = np.array(c)
        diag[k] = ecore + 2 * h1[s, s].sum() + (2 * J[np.ix_(s, s)] - K[np.ix_(s, s)]).sum()
        occ = set(c)
        vuoti = [a for a in range(n) if a not in occ]
        for i in c:
            for a in vuoti:
                nuovo = tuple(sorted((occ - {i}) | {a}))
                righe.append(k)
                col.append(pos[nuovo])
                val.append(P[i, a])
    H = sp.csr_matrix((val, (righe, col)), shape=(len(conf), len(conf)))
    H = H + sp.diags(diag)
    return H, conf


def controllo_piccolo(h1, g2):
    """Sottoproblema 8 orbitali / 4 coppie: DOCI qui contro un secondo codice basato su determinanti
    (operatori di seconda quantizzazione espliciti, solo determinanti con ogni orbitale vuoto o pieno)."""
    n, npairs = 8, 4
    h, g = h1[:n, :n], g2[:n, :n, :n, :n]
    H, conf = matrice_doci(h, g, npairs)
    e1 = np.linalg.eigvalsh(H.toarray())[0]
    # secondo codice: elementi di matrice <S'|H|S> con H = sum h_pq E_pq + 1/2 sum (pq|rs)(E_pq E_rs - d_qr E_ps)
    # applicato a stati a coppie rappresentati come insiemi di orbitali doppi; calcolo per spin-orbitali.
    def det(c):  # occupazioni spin-orbitali, ordine (0a,0b,1a,1b,...)
        return tuple(sorted([2 * i for i in c] + [2 * i + 1 for i in c]))

    def applica(ops, occ):
        occ = list(occ)
        segno = 1
        for tipo, so in reversed(ops):
            if tipo == "a":
                if so not in occ:
                    return None, 0
                k = occ.index(so)
                segno *= (-1) ** k
                occ.pop(k)
            else:
                if so in occ:
                    return None, 0
                k = sum(1 for x in occ if x < so)
                segno *= (-1) ** k
                occ.insert(k, so)
        return tuple(occ), segno

    dets = [det(c) for c in conf]
    idx = {dd: k for k, dd in enumerate(dets)}
    M = np.zeros((len(dets), len(dets)))
    so = range(2 * n)
    for k, dd in enumerate(dets):
        for p in so:
            for q in so:
                if p % 2 != q % 2:
                    continue
                v = h[p // 2, q // 2]
                if v == 0:
                    continue
                r, s_ = applica([("c", p), ("a", q)], dd)
                if r in idx:
                    M[idx[r], k] += s_ * v
        for p in so:
            for q in so:
                for r_ in so:
                    for s in so:
                        if p % 2 != s % 2 or q % 2 != r_ % 2:
                            continue
                        v = g[p // 2, s // 2, q // 2, r_ // 2]
                        if v == 0:
                            continue
                        r, s_ = applica([("c", p), ("c", q), ("a", r_), ("a", s)], dd)
                        if r in idx:
                            M[idx[r], k] += 0.5 * s_ * v
    e2 = np.linalg.eigvalsh(M)[0]
    return float(e1), float(e2), float(np.abs(H.toarray() - M).max())


if __name__ == "__main__":
    fcidump, out = sys.argv[1:3]
    t0 = time.time()
    d = leggi(fcidump)
    h1, g2, ecore = d["h1"], d["g2"], d["ecore"]
    npairs = d["nelec"] // 2
    e_a, e_b, dmax = controllo_piccolo(h1, g2)
    print("controllo 8 orbitali/4 coppie: E_doci", repr(e_a), "E_secondo_codice", repr(e_b),
          "max diff elementi", dmax, flush=True)
    H, conf = matrice_doci(h1, g2, npairs, ecore)
    print("dimensione", H.shape[0], "non zero", H.nnz, "tempo costruzione", time.time() - t0, flush=True)
    w, v = spla.eigsh(H, k=1, which="SA", tol=1e-12)
    e0 = float(w[0])
    res = float(np.linalg.norm(H @ v[:, 0] - e0 * v[:, 0]))
    print("E_DOCI", repr(e0), "residuo", res, "tempo", time.time() - t0, flush=True)
    json.dump(dict(fcidump_sha256=sha256(fcidump), dimensione=H.shape[0], E_DOCI=e0, residuo=res,
                   controllo_piccolo=dict(E_doci=e_a, E_secondo_codice=e_b, max_diff_elementi=dmax),
                   valore_issue_248=-116.24197244, differenza_con_248=e0 - (-116.24197244),
                   secondi=time.time() - t0), open(out, "w"), indent=1)
