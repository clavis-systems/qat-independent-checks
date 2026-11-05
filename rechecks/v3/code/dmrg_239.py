"""V3: un punto DMRG per la #239 con scelte diverse da C2.
- integrali letti con pyscf.tools.fcidump (non il lettore di C2);
- localizzazione Edmiston-Ruedenberg scritta da V3 (Jacobi a 2x2 con fit esatto c0+c1cos4t+c2sin4t);
- ordine degli orbitali: algoritmo genetico di block2 (gaopt), non Fiedler;
- simmetria SU(2), singoletto; M fino a 750; poi sweep a M decrescente per estrapolare."""
import sys, json, time, numpy as np
from pyscf.tools import fcidump
from pyscf import ao2mo
t0 = time.time()
d = fcidump.read(sys.argv[1])
n = d['NORB']; ne = d['NELEC']; ecore = d['ECORE']
h1 = np.asarray(d['H1']); g2 = ao2mo.restore(1, d['H2'], n)
print("norb", n, "nelec", ne, "ecore", ecore, flush=True)

def er_obj(g):
    return sum(g[k, k, k, k] for k in range(len(g)))

def rot_pair(h, g, U, i, j, t):
    c, s = np.cos(t), np.sin(t)
    R = np.eye(len(h)); R[i, i] = c; R[j, j] = c; R[i, j] = -s; R[j, i] = s  # nuovo_i = c*i + s*j ; nuovo_j = -s*i + c*j
    # colonna k della matrice R = coefficienti dell'orbitale nuovo k sui vecchi
    h = R.T @ h @ R
    g = np.einsum('pqrs,pa->aqrs', g, R, optimize=True)
    g = np.einsum('aqrs,qb->abrs', g, R, optimize=True)
    g = np.einsum('abrs,rc->abcs', g, R, optimize=True)
    g = np.einsum('abcs,sd->abcd', g, R, optimize=True)
    return h, g, U @ R

def gval(g, i, j, t):
    c, s = np.cos(t), np.sin(t)
    v1 = np.zeros(len(g)); v1[i] = c; v1[j] = s
    v2 = np.zeros(len(g)); v2[i] = -s; v2[j] = c
    idx = [i, j]
    sub = g[np.ix_(idx, idx, idx, idx)]
    a = v1[idx]; b = v2[idx]
    return np.einsum('pqrs,p,q,r,s', sub, a, a, a, a) + np.einsum('pqrs,p,q,r,s', sub, b, b, b, b)

U = np.eye(n); h, g = h1.copy(), g2.copy()
obj0 = er_obj(g)
for sweep in range(30):
    dmax = 0
    for i in range(n):
        for j in range(i + 1, n):
            g0 = gval(g, i, j, 0); g1 = gval(g, i, j, np.pi / 8); g2_ = gval(g, i, j, np.pi / 4)
            c0 = (g0 + g2_) / 2; c1 = (g0 - g2_) / 2; c2 = g1 - c0
            t = np.arctan2(c2, c1) / 4
            gain = np.hypot(c1, c2) - c1
            if gain > 1e-12:
                h, g, U = rot_pair(h, g, U, i, j, t); dmax = max(dmax, gain)
    print(f"ER sweep {sweep}: obiettivo {er_obj(g):.10f} guadagno max {dmax:.2e}", flush=True)
    if dmax < 1e-9:
        break
print("ER: obiettivo iniziale", obj0, "finale", er_obj(g), "ortogonalita'", np.abs(U.T @ U - np.eye(n)).max(), flush=True)
# controllo: integrali ruotati da zero
gchk = ao2mo.incore.full(ao2mo.restore(8, g2, n), U, compact=False).reshape(n, n, n, n)
print("controllo trasformazione da zero:", np.abs(gchk - g).max(), np.abs(U.T @ h1 @ U - h).max(), flush=True)
np.save(sys.argv[2] + "_U.npy", U)

from pyblock2.driver.core import DMRGDriver, SymmetryTypes
driver = DMRGDriver(scratch=sys.argv[3], symm_type=SymmetryTypes.SU2, n_threads=4, stack_mem=int(9e9))
idx = driver.orbital_reordering(h, g, method='gaopt')
print("ordine gaopt:", list(idx), flush=True)
h = h[idx][:, idx]; g = g[idx][:, idx][:, :, idx][:, :, :, idx]
driver.initialize_system(n_sites=n, n_elec=ne, spin=0, orb_sym=[0] * n)
mpo = driver.get_qc_mpo(h1e=h, g2e=g, ecore=ecore, iprint=0)
ket = driver.get_random_mps(tag="KET", bond_dim=250, nroots=1)
bd = [250] * 6 + [500] * 4 + [750] * 4 + [650, 550, 450, 350]
nz = [1e-5] * 6 + [1e-6] * 4 + [1e-7, 1e-7, 0, 0] + [0] * 4
th = [1e-7] * 10 + [1e-9] * 8
E = driver.dmrg(mpo, ket, n_sweeps=len(bd), bond_dims=bd, noises=nz, thrds=th, tol=0, iprint=1)
dm = driver._dmrg
res = {"bond_dims": [int(x) for x in dm.sweep_bond_dims], "energies": [float(np.array(e).ravel()[0]) for e in dm.sweep_energies],
       "discarded": [float(x) for x in dm.sweep_discarded_weights], "ordine": [int(x) for x in idx], "tempo_s": time.time() - t0}
print(json.dumps(res, indent=1), flush=True)
json.dump(res, open(sys.argv[2] + ".json", "w"), indent=1)
