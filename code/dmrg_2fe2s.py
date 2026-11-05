"""DMRG per [2Fe-2S] (30e, 20o) del Quantum Advantage Tracker, scritto per C2.

Uso: python -I dmrg_2fe2s.py FCIDUMP SIMMETRIA(SU2|SZ) ORDINE(fiedler|file) THREADS SCRATCH OUT.json [MMAX] [2S]
(2S = due volte lo spin totale, solo per SU2; 0 = singoletto, 2 = tripletto)

- Integrali letti con il nostro lettore (leggi_fcidump.py), non con quello di block2.
- Fase 1: sweep a M crescente (10 a M=250, 6 a M=500, con rumore), poi 4 sweep senza rumore a M massimo.
- Fase 2: un sweep per ogni M decrescente, senza rumore, per estrapolare l'energia rispetto al
  peso scartato (discarded weight).
- Controllo interno: energia ricalcolata contraendo 1-RDM e 2-RDM dell'MPS finale con i nostri integrali.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from leggi_fcidump import leggi, sha256  # noqa: E402

from pyblock2.driver.core import DMRGDriver, SymmetryTypes  # noqa: E402

fcidump, simm, ordine, nthr, scratch, out = sys.argv[1:7]
mmax = int(sys.argv[7]) if len(sys.argv) > 7 else 1500
due_s = int(sys.argv[8]) if len(sys.argv) > 8 else 0
nthr = int(nthr)
t0 = time.time()

d = leggi(fcidump)
n, ne, ms2 = d["norb"], d["nelec"], d["ms2"]
h1, g2, ecore = d["h1"], d["g2"], d["ecore"]
h1_orig, g2_orig = h1.copy(), g2.copy()
# rotazione facoltativa degli orbitali attivi (variabile d'ambiente DMRG_ROT = file .npy con U, colonne = nuovi
# orbitali nei vecchi). Nel CAS completo l'energia esatta non cambia; il controllo RDM finale si fa nella base
# ORIGINALE del file, con le RDM ruotate all'indietro.
rot = os.environ.get("DMRG_ROT")
U = np.load(rot, allow_pickle=False) if rot else np.eye(n)
if rot:
    h1 = U.T @ h1 @ U
    g2 = np.einsum("pi,qj,rk,sl,pqrs->ijkl", U, U, U, U, g2, optimize=True)
    print("rotazione applicata da", rot, "| sum (ii|ii) =", float(np.einsum("iiii->", g2)), flush=True)

st = SymmetryTypes.SU2 if simm == "SU2" else SymmetryTypes.SZ
drv = DMRGDriver(scratch=scratch, symm_type=st, n_threads=nthr, stack_mem=int(3e9))
spin = due_s if simm == "SU2" else ms2
drv.initialize_system(n_sites=n, n_elec=ne, spin=spin, orb_sym=[0] * n)

if ordine == "fiedler":
    idx = np.array(drv.orbital_reordering(h1, g2), dtype=int)
else:
    idx = np.arange(n)
h1r = h1[idx][:, idx]
g2r = g2[idx][:, idx][:, :, idx][:, :, :, idx]
print("ordine orbitali:", idx.tolist(), flush=True)

mpo = drv.get_qc_mpo(h1e=h1r, g2e=g2r, ecore=ecore, iprint=0)
ket = drv.get_random_mps(tag="KET", bond_dim=min(250, mmax), nroots=1)

# schema rivisto alle 19:40: molti sweep a M piccolo (economici) prima di salire; al M massimo niente rumore
# (con il rumore a M=1000 il processo superava i 13 GB ed e' stato ucciso; anche senza rumore a M=1000: 13,3 GB)
schema = {250: 8, 500: 6, 750: 4, 900: 4, 1000: 4}
passi = [m for m in (250, 500, 900 if mmax == 900 else 750, 1000) if m <= mmax] or [mmax]
bd, noise, thr = [], [], []
for m in passi[:-1]:
    k = schema.get(m, 4)
    bd += [m] * k
    noise += [1e-4] * (k // 2) + [1e-5] * (k - k // 2 - 1) + [0.0]
    thr += [1e-6] * (k // 2) + [1e-7] * (k - k // 2)
k = schema.get(passi[-1], 4)
bd += [passi[-1]] * k
noise += [0.0] * k
thr += [1e-8] * (k - 2) + [1e-9] * 2
e_cresc = drv.dmrg(mpo, ket, n_sweeps=len(bd), bond_dims=bd, noises=noise, thrds=thr,
                   iprint=1, twosite_to_onesite=None)
bds, dws, ens = drv.get_dmrg_results()
fase1 = [dict(M=int(b), dw=float(np.max(w)), E=float(np.min(e))) for b, w, e in zip(bds, dws, ens)]
print("fase1 energia", repr(float(e_cresc)), "tempo", time.time() - t0, flush=True)


def salva(extra):
    r = dict(fcidump_sha256=sha256(fcidump), simmetria=simm, due_S=due_s, rotazione=rot, ordine=ordine, ordine_orbitali=idx.tolist(),
             threads=nthr, Mmax=passi[-1], energia_Mmax=float(e_cresc), sweep_fase1=fase1,
             secondi=time.time() - t0)
    r.update(extra)
    with open(out, "w") as f:
        json.dump(r, f, indent=1)


salva({})

# fase 2: un sweep per M decrescente senza rumore
lista_giu = (800, 700, 600, 500, 400) if mmax == 900 else (1800, 1500, 1200, 1000, 800, 650, 550, 450, 350, 250, 40, 30, 20)
giu = [m for m in lista_giu if m < passi[-1]][:5]
punti = []
for m in giu:
    e = drv.dmrg(mpo, ket, n_sweeps=1, bond_dims=[m], noises=[0.0], thrds=[1e-9],
                 iprint=0, twosite_to_onesite=None)
    r = drv.get_dmrg_results()
    dw = float(np.max(r[1][-1]))
    punti.append(dict(M=m, energia=float(e), dw=dw))
    print("giu M", m, "E", repr(float(e)), "dw", dw, "tempo", time.time() - t0, flush=True)
    salva(dict(sweep_giu=punti))

x = np.array([p["dw"] for p in punti])
y = np.array([p["energia"] for p in punti])
estr = None
if len(punti) >= 3:
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    estr = dict(pendenza=float(coef[0]), E_dw0=float(coef[1]),
                errore_stimato_regola_1_5=float(abs(coef[1] - float(e_cresc)) / 5.0))
    print("estrapolazione E(dw->0)", repr(float(coef[1])), flush=True)
salva(dict(sweep_giu=punti, estrapolazione=estr))

# controllo interno: energia dell'MPS finale (M piu' piccolo) da RDM e nostri integrali
e_mps = float(drv.expectation(ket, mpo, ket) / drv.expectation(ket, drv.get_identity_mpo(), ket))
dm1 = drv.get_1pdm(ket)
dm2 = drv.get_2pdm(ket)
# convenzione verificata dalle tracce: dm2[p,q,r,s] = sum_{st} <a+_{p s} a+_{q t} a_{r t} a_{s s}>
# E = ecore + sum h_pq D_pq + 1/2 sum (pq|rs) dm2[p,r,s,q]
if simm == "SU2":
    D1, e2 = dm1, 0.5 * np.einsum("pqrs,prsq", g2r, dm2)
    tr2 = float(np.einsum("pqqp", dm2))
else:
    D1 = dm1[0] + dm1[1]
    e2 = 0.5 * (np.einsum("pqrs,prsq", g2r, dm2[0]) + 2 * np.einsum("pqrs,prsq", g2r, dm2[1])
                + np.einsum("pqrs,prsq", g2r, dm2[2]))
    tr2 = float(np.einsum("pqqp", dm2[0]) + 2 * np.einsum("pqqp", dm2[1]) + np.einsum("pqqp", dm2[2]))
e_rdm = float(ecore + np.einsum("pq,pq", h1r, D1) + e2)
# stessa energia nella base originale del file: W = matrice dai siti DMRG agli orbitali originali
W = U[:, idx]
D1o = W @ D1 @ W.T
if simm == "SU2":
    D2o = np.einsum("pa,qb,rc,sd,abcd->pqrs", W, W, W, W, dm2, optimize=True)
else:
    D2o = sum(f * np.einsum("pa,qb,rc,sd,abcd->pqrs", W, W, W, W, dm2[k], optimize=True)
              for k, f in ((0, 1.0), (1, 2.0), (2, 1.0)))
e_rdm_orig = float(ecore + np.einsum("pq,pq", h1_orig, D1o) + 0.5 * np.einsum("pqrs,prsq", g2_orig, D2o))
print("controllo base originale: E da RDM ruotate all'indietro", repr(e_rdm_orig), flush=True)
print("controllo: E(MPS finale)", repr(e_mps), "E da RDM", repr(e_rdm), "diff", e_rdm - e_mps,
      "traccia 1RDM", float(np.trace(D1)), "traccia 2RDM", tr2, flush=True)
salva(dict(sweep_giu=punti, estrapolazione=estr, controllo_rdm=dict(
    E_mps_finale=e_mps, E_da_rdm=e_rdm, differenza=e_rdm - e_mps, E_da_rdm_base_originale=e_rdm_orig,
    traccia_1rdm=float(np.trace(D1)), traccia_2rdm=tr2)))
print("FINE", time.time() - t0, flush=True)
