"""Localizzazione di Edmiston-Ruedenberg dei 20 orbitali attivi, solo dagli integrali (scritto per C2).

Uso: python -I localizza_er.py FCIDUMP OUT_U.npy
Massimizza sum_i (ii|ii) con rotazioni di Jacobi tra coppie di orbitali. In un CAS completo l'energia esatta
non cambia ruotando gli orbitali attivi tra loro; cambia solo la difficoltà per il DMRG.
Per una coppia (i, j) ruotata di theta, (i'i'|i'i') + (j'j'|j'j') = a0 + a1 cos 4θ + b1 sin 4θ: tre valori
bastano per trovare il massimo; il codice controlla che la somma non scenda mai.
U: colonne = nuovi orbitali scritti nei vecchi (phi'_p = sum_q U[q, p] phi_q).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from leggi_fcidump import leggi  # noqa: E402


def ruota_coppia(g, h, i, j, th):
    c, s = np.cos(th), np.sin(th)
    R = np.array([[c, s], [-s, c]])
    ij = [i, j]
    h[ij, :] = R @ h[ij, :]
    h[:, ij] = h[:, ij] @ R.T
    for ax in range(4):
        gm = np.moveaxis(g, ax, 0)
        gm[ij] = np.tensordot(R, gm[ij], axes=1)
    return R


def f_coppia(g, i, j, th):
    c, s = np.cos(th), np.sin(th)
    R = np.array([[c, s], [-s, c]])
    T = g[np.ix_([i, j], [i, j], [i, j], [i, j])]
    T2 = np.einsum("ap,bq,cr,ds,pqrs->abcd", R, R, R, R, T)
    return T2[0, 0, 0, 0] + T2[1, 1, 1, 1]


def localizza(h1, g2, giri=30, tol=1e-10):
    n = h1.shape[0]
    h, g = h1.copy(), g2.copy()
    U = np.eye(n)
    storia = [float(np.einsum("iiii->", g))]
    for giro in range(giri):
        for i in range(n):
            for j in range(i + 1, n):
                f0 = f_coppia(g, i, j, 0.0)
                f1 = f_coppia(g, i, j, np.pi / 8)
                f2 = f_coppia(g, i, j, np.pi / 4)
                a0, a1 = (f0 + f2) / 2, (f0 - f2) / 2
                b1 = f1 - a0
                th = np.arctan2(b1, a1) / 4
                if f_coppia(g, i, j, th) > f0 + 1e-14:
                    R = ruota_coppia(g, h, i, j, th)
                    G = np.eye(n)
                    G[np.ix_([i, j], [i, j])] = R.T
                    U = U @ G
        storia.append(float(np.einsum("iiii->", g)))
        assert storia[-1] >= storia[-2] - 1e-12
        if storia[-1] - storia[-2] < tol:
            break
    return U, h, g, storia


if __name__ == "__main__":
    d = leggi(sys.argv[1])
    U, h, g, storia = localizza(d["h1"], d["g2"])
    # controllo: le integrali ruotate con U da zero coincidono con quelle aggiornate passo per passo
    h_chk = U.T @ d["h1"] @ U
    g_chk = np.einsum("pi,qj,rk,sl,pqrs->ijkl", U, U, U, U, d["g2"], optimize=True)
    print("giri", len(storia) - 1, "sum (ii|ii) prima", storia[0], "dopo", storia[-1])
    print("U ortogonale:", np.abs(U.T @ U - np.eye(len(U))).max(),
          "| scarto h", np.abs(h_chk - h).max(), "| scarto g", np.abs(g_chk - g).max())
    J = np.einsum("iijj->ij", g)
    np.set_printoptions(precision=3, suppress=True, linewidth=200)
    print("J_ii dopo:", np.diag(J))
    np.save(sys.argv[2], U)
