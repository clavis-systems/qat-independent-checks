"""V3: DOCI (seniority zero) per la #248, codice scritto da V3.
Configurazioni = insiemi di orbitali doppiamente occupati (bitmask). Formule:
 diag = sum_i n_i (2 h_ii + (ii|ii)) + sum_{i<j} n_i n_j (4 (ii|jj) - 2 (ij|ij))   [notazione chimica, orbitali reali]
 fuori diagonale (coppia i -> a): (ia|ia)
Controllo interno: la stessa matrice su un sottoproblema piccolo confrontata con il blocco a seniority zero della
matrice FCI completa costruita con pyscf (fci.direct_spin1.contract_2e su vettori unitari)."""
import sys, json, time, itertools, hashlib, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, scipy.linalg as sl
from pyscf.tools import fcidump
from pyscf import ao2mo, fci
from pyscf.fci import cistring

def doci_matrix(h, g, npair):
    n = len(h)
    confs = [sum(1 << k for k in c) for c in itertools.combinations(range(n), npair)]
    pos = {m: k for k, m in enumerate(confs)}
    J = np.einsum('iijj->ij', g); K = np.einsum('ijij->ij', g); hd = np.diag(h)
    rows, cols, vals = [], [], []
    diag = np.zeros(len(confs))
    for k, m in enumerate(confs):
        occ = [i for i in range(n) if m >> i & 1]; vir = [a for a in range(n) if not m >> a & 1]
        o = np.array(occ)
        diag[k] = np.sum(2 * hd[o] + J[o, o]) + (np.sum(4 * J[np.ix_(o, o)] - 2 * K[np.ix_(o, o)]) - np.sum(4 * J[o, o] - 2 * K[o, o])) / 2
        for i in occ:
            for a in vir:
                rows.append(pos[m ^ (1 << i) ^ (1 << a)]); cols.append(k); vals.append(K[i, a])
    Hm = sp.csr_matrix((vals, (rows, cols)), shape=(len(confs), len(confs))) + sp.diags(diag)
    return Hm, confs

def check_small(h, g, orbs, npair):
    hs = h[np.ix_(orbs, orbs)]; gs = g[np.ix_(orbs, orbs, orbs, orbs)]; n = len(orbs)
    Hm, confs = doci_matrix(hs, gs, npair)
    na = cistring.num_strings(n, npair)
    h2 = fci.direct_spin1.absorb_h1e(hs, gs, n, (npair, npair), .5)
    full = np.zeros((na * na, na * na))
    for c in range(na * na):
        v = np.zeros(na * na); v[c] = 1
        full[:, c] = fci.direct_spin1.contract_2e(h2, v.reshape(na, na), n, (npair, npair)).ravel()
    addr = [cistring.str2addr(n, npair, m) for m in confs]
    sel = [a * na + a for a in addr]
    blk = full[np.ix_(sel, sel)]
    diff = np.abs(blk - Hm.toarray()).max()
    e_doci = np.linalg.eigvalsh(Hm.toarray())[0]; e_fci = np.linalg.eigvalsh(full)[0]
    return diff, e_doci, e_fci

if __name__ == "__main__":
    t0 = time.time()
    path = sys.argv[1]
    print("sha256", hashlib.sha256(open(path, 'rb').read()).hexdigest())
    d = fcidump.read(path); n = d['NORB']; ne = d['NELEC']; ecore = d['ECORE']
    h = np.asarray(d['H1']); g = ao2mo.restore(1, d['H2'], n)
    out = {}
    for orbs, npair in [(list(range(10, 17)), 3), ([0, 3, 7, 12, 15, 18], 3), (list(range(12, 20)), 4)]:
        diff, ed, ef = check_small(h, g, orbs, npair)
        print(f"controllo pyscf su orbitali {orbs}, {npair} coppie: max|diff matrice| = {diff:.2e}; E_DOCI={ed:.10f} >= E_FCI={ef:.10f}", flush=True)
        out.setdefault("controlli", []).append({"orbitali": orbs, "coppie": npair, "diff_max": diff, "E_doci": ed, "E_fci": ef})
    Hm, confs = doci_matrix(h, g, ne // 2)
    print("dimensione", len(confs), "non zeri", Hm.nnz, "simmetria", abs(Hm - Hm.T).max(), flush=True)
    w, v = sla.eigsh(Hm, k=3, which='SA', tol=1e-14)
    E = w[0] + ecore; r = np.linalg.norm(Hm @ v[:, 0] - w[0] * v[:, 0])
    print(f"Lanczos: E_DOCI = {E:.12f}  (residuo {r:.1e}); successivi {w[1]+ecore:.8f} {w[2]+ecore:.8f}", flush=True)
    wd = sl.eigh(Hm.toarray(), eigvals_only=True, subset_by_index=[0, 1])
    print(f"densa (LAPACK): E_DOCI = {wd[0]+ecore:.12f}", flush=True)
    out.update({"E_DOCI_lanczos": float(E), "residuo": float(r), "E_DOCI_densa": float(wd[0] + ecore), "dim": len(confs),
                "confronto_C2": float(E - (-116.2419724396)), "confronto_issue": float(E - (-116.24197244)), "tempo_s": time.time() - t0})
    print(json.dumps(out, indent=1))
    json.dump(out, open(sys.argv[2], "w"), indent=1)
