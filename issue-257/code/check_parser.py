"""Cross-check of the own QASM parser against qiskit.qasm2 (independent parser),
and of the logical layout: every non-SWAP 2-qubit gate must be same-rung or
nearest-neighbour same-chain after undoing the SWAP network."""
import collections
import sys

import qiskit
import qiskit.qasm2

import lsh
from qasm_parse import parse_qasm

for path in sys.argv[1:]:
    nq, mine = parse_qasm(path)
    qc = qiskit.qasm2.load(path, custom_instructions=qiskit.qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    theirs = []
    for inst in qc.data:
        theirs.append((inst.operation.name, tuple(float(p) for p in inst.operation.params),
                       tuple(qc.find_bit(q).index for q in inst.qubits)))
    same = len(mine) == len(theirs) and all(
        a[0] == b[0] and a[2] == b[2] and len(a[1]) == len(b[1]) and all(abs(x - y) < 1e-15 for x, y in zip(a[1], b[1]))
        for a, b in zip(mine, theirs))
    print(path, "qiskit", qiskit.__version__, "ops", len(theirs), "own parser identical:", same)
    circ = lsh.load_circuit(path)
    kinds = collections.Counter()
    for name, par, qs in circ["logical"]:
        if len(qs) == 2:
            (ra, ca), (rb, cb) = lsh.site_chain(qs[0]), lsh.site_chain(qs[1])
            kinds[("same-rung" if ra == rb else "dsite=%d" % abs(ra - rb)) + (" same-chain" if ca == cb else " cross-chain")] += 1
    print("  2q logical gate kinds:", dict(kinds))
    print("  steps found:", len(circ["cuts"]), "logical gates:", len(circ["logical"]), "last cut:", circ["cuts"][-1],
          "raw ops after last cut:", len(circ["raw"]) - circ["raw_cuts"][-1])
    fp = circ["final_perm"]
    print("  final perm: output wire 2r -> chain0 of site r, 2r+1 -> chain1:",
          all(lsh.site_chain(fp[w]) == (w // 2, w % 2) for w in range(nq)))
    init = [(lsh.site_chain(w)) for w in range(nq) if circ["bits"][w]]
    occ = collections.defaultdict(lambda: [0, 0])
    for r, c in init:
        occ[r][c] = 1
    print("  initial occupations (site: n_chain0 n_chain1) for occupied sites:",
          " ".join("%d:%d%d" % (r, occ[r][0], occ[r][1]) for r in sorted(occ)))
    # per-step gate statistics
    prev = 0
    hops = []
    for cut in circ["cuts"]:
        seg = circ["logical"][prev:cut]
        hops.append(sum(1 for g in seg if g[0] == "rz" and abs(abs(g[1][0]) - 0.15) < 1e-12))
        prev = cut
    print("  rz(+-0.15) per step:", hops)
