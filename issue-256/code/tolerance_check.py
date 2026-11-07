"""Robustness of the numerical thresholds used by the exact simplifications (tn.is_diag 1e-14,
tn.is_id_phase 1e-13): report the gap between 1q products treated as diagonal/identity and the
rest, over every merge performed during the L=3 peephole (instrumented copies of the tests)."""
import numpy as np
import structure, tn

orig_diag, orig_id = tn.is_diag, tn.is_id_phase
rec_off, rec_id = [], []


def is_diag(M, tol=1e-14):
    rec_off.append(max(abs(M[0, 1]), abs(M[1, 0])))
    return orig_diag(M, tol)


def is_id_phase(M, tol=1e-13):
    if max(abs(M[0, 1]), abs(M[1, 0])) < 1e-14:
        rec_id.append(abs(M[0, 0] - M[1, 1]))
    return orig_id(M, tol)


tn.is_diag, tn.is_id_phase = is_diag, is_id_phase
info = structure.analyse(verbose=False)
ops = tn.peephole(info["all_gates"])
ops = tn.lightcone(ops, structure.OBS)
off = np.array(rec_off); idd = np.array(rec_id)
print("off-diagonal magnitudes tested:", len(off))
print("  largest value treated as diagonal (<1e-14):", off[off < 1e-14].max() if (off < 1e-14).any() else None)
print("  smallest value treated as non-diagonal   :", off[off >= 1e-14].min())
print("identity tests on diagonal products:", len(idd))
print("  largest |d0-d1| treated as identity (<1e-13):", idd[idd < 1e-13].max() if (idd < 1e-13).any() else None)
print("  smallest |d0-d1| treated as non-identity  :", idd[idd >= 1e-13].min() if (idd >= 1e-13).any() else None)
print("final ops", len(ops))
