"""Extra check: recompute n_f(20) from the stored per-site occupations with the
#250 wire-based definition (output wire w, site r on wires 2r, 2r+1, sign (-1)^r),
using the final SWAP permutation, and compare with the #254 site-based value."""
import json
import sys

import lsh

qasm_scv, qasm_mid, res_scv, res_mid = sys.argv[1:5]
out = {}
for tag, q, r in (("SCV", qasm_scv, res_scv), ("MID", qasm_mid, res_mid)):
    circ = lsh.load_circuit(q)
    perm = circ["final_perm"]
    with open(r) as f:
        rows = json.load(f)["rows"]
    n = rows[-1]["n"]  # n[site][chain] of logical qubits at step 20
    nw = [n[lsh.site_chain(perm[w])[0]][lsh.site_chain(perm[w])[1]] for w in range(circ["nq"])]
    out[tag] = (nw, rows[-1]["N"])
nf_wire = sum((-1) ** r * ((out["MID"][0][2 * r] + out["MID"][0][2 * r + 1]) - (out["SCV"][0][2 * r] + out["SCV"][0][2 * r + 1]))
              for r in range(60))
nf_site = out["MID"][1] - out["SCV"][1]
print("n_f(20) wire-based (#250 def) = %.9f ; site-based (#254 def) = %.9f ; diff %.1e" % (nf_wire, nf_site, nf_wire - nf_site))
print("Q_SCV = %.12f, Q_MID = %.12f" % (sum(out["SCV"][0]), sum(out["MID"][0])))
