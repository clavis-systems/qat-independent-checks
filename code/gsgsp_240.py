"""Controllo della #240 (guided_sparse_ground_state_problem_49Q_v1, "Exact Diagonalization"), scritto per C2.

Uso: python -I gsgsp_240.py HAMILTONIANA.npy INDICE_STATO_GUIDA.npy OUT.json [LIMITE_STATI]

Idea (nostra, non presa dal codice di chi ha inviato la richiesta): ogni termine di Pauli P = (fase) X^x Z^z
manda uno stato di base |b> in un solo stato |b xor x>. Partendo dallo stato guida si visitano tutti gli stati
raggiungibili (ricerca in ampiezza sulle maschere x). Se l'insieme raggiungibile è piccolo, l'Hamiltoniana
ristretta a quell'insieme è esatta (il sottospazio è invariante) e si diagonalizza.
Controllo interno: su un'Hamiltoniana casuale a 10 qubit lo stesso codice deve dare l'autovalore più basso
del blocco raggiungibile calcolato con la matrice densa costruita con prodotti di Kronecker.
"""
import json
import sys
import time
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def parse(arr):
    termini = []
    for riga in arr:
        p, c = str(riga[0]), complex(str(riga[1]).replace(" ", ""))
        termini.append((p, c))
    return termini


def maschere(p):
    """Ordine Qiskit: l'ultimo carattere è il qubit 0."""
    n = len(p)
    x = z = 0
    ny = 0
    for k, ch in enumerate(reversed(p)):
        if ch in "XY":
            x |= 1 << k
        if ch in "ZY":
            z |= 1 << k
        if ch == "Y":
            ny += 1
    return x, z, ny, n


def prepara(termini):
    gruppi = defaultdict(list)  # x -> lista di (z, coeff * i^ny)
    n = None
    for p, c in termini:
        x, z, ny, n = maschere(p)
        gruppi[x].append((z, c * (1j) ** ny))
    return gruppi, n


def azione(gruppi, b):
    """<b xor x | P | b> per ogni x: Y = i X Z, quindi P|b> = c i^ny (-1)^{popcount(b & z)} |b xor x>."""
    out = {}
    for x, lst in gruppi.items():
        v = 0j
        for z, c in lst:
            v += c * (-1) ** bin(b & z).count("1")
        if v != 0:
            out[b ^ x] = out.get(b ^ x, 0) + v
    return out


def raggiungibili(gruppi, b0, limite):
    visti = {b0: 0}
    coda = [b0]
    righe, col, val = [], [], []
    k = 0
    while k < len(coda):
        b = coda[k]
        for b2, v in azione(gruppi, b).items():
            if b2 not in visti:
                if len(visti) >= limite:
                    return None
                visti[b2] = len(visti)
                coda.append(b2)
            righe.append(visti[b2])
            col.append(visti[b])
            val.append(v)
        k += 1
    H = sp.csr_matrix((val, (righe, col)), shape=(len(visti), len(visti)))
    return H, coda


def kron_denso(termini, n):
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1.0, -1.0])
    M = {"I": I, "X": X, "Y": Y, "Z": Z}
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for p, c in termini:
        A = np.array([[1.0]])
        for ch in p:  # primo carattere = qubit più alto: kron in quest'ordine dà l'ordine Qiskit
            A = np.kron(A, M[ch])
        H += c * A
    return H


def controllo_interno(seed=1):
    rng = np.random.default_rng(seed)
    n = 10
    termini = []
    for _ in range(25):
        p = "".join(rng.choice(list("IIXZZY"), size=n))
        termini.append((p, complex(rng.normal())))
    termini = [(p, c if any(ch in "XY" for ch in p) and p.count("Y") % 2 == 0 else c.real) for p, c in termini]
    # rendo hermitiana: coefficienti reali bastano per stringhe di Pauli
    termini = [(p, complex(c.real)) for p, c in termini]
    gruppi, _ = prepara(termini)
    b0 = 5
    H, stati = raggiungibili(gruppi, b0, 10 ** 6)
    e_nostro = np.linalg.eigvalsh(H.toarray())[0]
    Hd = kron_denso(termini, n)
    idx = np.array(stati)
    e_denso_blocco = np.linalg.eigvalsh(Hd[np.ix_(idx, idx)])[0]
    fuori = np.abs(Hd[np.ix_(np.setdiff1d(np.arange(2 ** n), idx), idx)]).max() if len(idx) < 2 ** n else 0.0
    return float(e_nostro), float(e_denso_blocco), float(fuori), len(stati)


if __name__ == "__main__":
    fh, fi, out = sys.argv[1:4]
    limite = int(sys.argv[4]) if len(sys.argv) > 4 else 2_000_000
    t0 = time.time()
    e_a, e_b, fuori, ns = controllo_interno()
    print("controllo interno 10 qubit: E nostro", e_a, "E denso", e_b, "max elemento fuori dal blocco", fuori,
          "stati", ns, flush=True)
    arr = np.load(fh, allow_pickle=False)
    b0 = int(np.load(fi, allow_pickle=False))
    termini = parse(arr)
    print("termini", len(termini), "forma file", arr.shape, "stato guida", b0, flush=True)
    imm = max(abs(c.imag) for _, c in termini)
    gruppi, n = prepara(termini)
    print("qubit", n, "maschere x distinte", len(gruppi), "max parte immaginaria coeff", imm, flush=True)
    r = raggiungibili(gruppi, b0, limite)
    if r is None:
        print("insieme raggiungibile oltre il limite", limite, "tempo", time.time() - t0, flush=True)
        json.dump(dict(esito="oltre il limite", limite=limite), open(out, "w"))
        sys.exit(0)
    H, stati = r
    herm = float(abs(H - H.getH()).max()) if H.nnz else 0.0
    print("stati raggiungibili", len(stati), "non zero", H.nnz, "errore hermitianità", herm,
          "tempo", time.time() - t0, flush=True)
    if len(stati) <= 3000:
        w = np.linalg.eigvalsh(H.toarray())
        e0 = float(w[0])
        gap = float(w[1] - w[0]) if len(w) > 1 else None
    else:
        w, _ = spla.eigsh(H, k=2, which="SA", tol=1e-14)
        e0, gap = float(w[0]), float(w[1] - w[0])
    print("E0 nel sottospazio invariante", repr(e0), "gap", gap, "tempo", time.time() - t0, flush=True)
    json.dump(dict(n_qubit=n, termini=len(termini), stato_guida=b0, stati_raggiungibili=len(stati),
                   E0=e0, gap_nel_sottospazio=gap, errore_hermitianita=herm,
                   valore_issue_240=-85.586625069864652, differenza_con_240=e0 - (-85.586625069864652),
                   controllo_interno=dict(E_nostro=e_a, E_denso_blocco=e_b, max_fuori_blocco=fuori, stati=ns),
                   secondi=time.time() - t0), open(out, "w"), indent=1)
