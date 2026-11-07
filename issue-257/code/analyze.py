"""Structure analysis of the LSH circuits: gate counts, SWAP tracking, layout."""
import collections
import sys

from qasm_parse import parse_qasm


def analyze(path):
    nq, gates = parse_qasm(path)
    print("file", path, "qubits", nq, "gates", len(gates))
    print("counts", dict(collections.Counter(g[0] for g in gates)))
    # leading x gates
    k = 0
    while k < len(gates) and gates[k][0] == "x":
        k += 1
    xs = sorted(g[2][0] for g in gates[:k])
    print("leading x gates:", k, "wires:", xs)
    rest = gates[k:]
    print("x gates after prefix:", sum(1 for g in rest if g[0] == "x"))
    # angle statistics
    angs = collections.Counter(round(g[1][0], 9) for g in rest if g[0] == "rz")
    print("rz angles:", sorted(angs.items()))
    # SWAP tracking: wire -> logical (initial wire label)
    perm = list(range(nq))
    pair_counter = collections.Counter()
    nswap = 0
    for name, par, qs in rest:
        if name == "swap":
            a, b = qs
            perm[a], perm[b] = perm[b], perm[a]
            nswap += 1
            continue
        if len(qs) == 2:
            la, lb = perm[qs[0]], perm[qs[1]]
            pair_counter[(min(la, lb), max(la, lb))] += 1
    print("swaps", nswap)
    # logical labels: initial wire w -> chain i if w<60 else o, site w%60  (hypothesis)
    def lab(l):
        return ("i" if l < 60 else "o", l % 60)
    kinds = collections.Counter()
    for (la, lb), c in pair_counter.items():
        ca, ra = lab(la)
        cb, rb = lab(lb)
        kinds[(ca + cb, abs(ra - rb))] += c
    print("2q gate kinds (chains, |dsite|) under hypothesis w<60 -> i:", dict(kinds))
    print("final permutation (wire -> logical label):")
    print([lab(p) for p in perm])
    return nq, gates, perm


if __name__ == "__main__":
    for p in sys.argv[1:]:
        analyze(p)
        print()
