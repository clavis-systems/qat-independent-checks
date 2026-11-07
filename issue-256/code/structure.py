"""Structure of the OLE circuit: lattice, site classes, halves, Trotter steps, symmetry checks.

Also builds the circuits for other depths L (1, 2, 3) from the Trotter step found in the file,
and checks that the rebuilt L=3 circuit is gate-for-gate identical to the file.
"""
import math
import collections
import numpy as np
from qasm_parse import parse_qasm, mat1

CIRC = "/circ/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm"
OBS = (52, 59, 72)


def key(g):
    return (g["name"], g["qubits"], tuple(round(p, 12) for p in g["params"]))


def inverse_gate(g):
    """Gate-level inverse, expressed with the same gate vocabulary."""
    n = g["name"]
    if n in ("rx", "rz"):
        p = (-g["params"][0],)
        return dict(name=n, qubits=g["qubits"], params=p, mat=mat1(n, p))
    if n == "cz":
        return dict(g)
    inv = {"s": "sdg", "sdg": "s", "sx": "sxdg", "sxdg": "sx"}[n]
    return dict(name=inv, qubits=g["qubits"], params=(), mat=g["mat"].conj().T)


def split(c):
    gates = c["gates"]
    names = [g["name"] for g in gates]
    i0 = names.index("sdg")                       # start of V block
    i1 = max(i for i, n in enumerate(names) if n == "s") + 1   # end of V block
    return gates[:i0], gates[i0:i1], gates[i1:]


def analyse(verbose=True):
    c = parse_qasm(CIRC)
    U, Vb, U2 = split(c)
    out = {}
    used = sorted({q for g in c["gates"] for q in g["qubits"]})
    edges = sorted({tuple(sorted(g["qubits"])) for g in c["gates"] if g["name"] == "cz"})
    deg = collections.Counter(q for e in edges for q in e)
    # site classes from first-half rx angles
    rxang = collections.defaultdict(set)
    for g in U:
        if g["name"] == "rx":
            rxang[g["qubits"][0]].add(round(abs(g["params"][0]), 12))
    S = sorted(q for q, a in rxang.items() if a == {0.5})
    F = sorted(q for q, a in rxang.items() if a == {round(3 * math.pi / 8, 12)})
    P = sorted(g["qubits"][0] for g in Vb if g["name"] == "rz")
    Vang = sorted({g["params"][0] for g in Vb if g["name"] == "rz"})
    out.update(used=used, edges=edges, S=S, F=F, P=P)
    # (a) second half == gate-wise inverse of first half
    invU = [inverse_gate(g) for g in reversed(U)]
    a_ok = len(invU) == len(U2) and all(key(x) == key(y) for x, y in zip(invU, U2))
    # (b) U = A then B with B == gate-wise inverse of A with rx sign flipped on S sites
    h = len(U) // 2
    A, B = U[:h], U[h:]
    def flipS(g):
        if g["name"] == "rx" and g["qubits"][0] in S:
            p = (-g["params"][0],)
            return dict(g, params=p, mat=mat1("rx", p))
        return g
    invA_S = [flipS(inverse_gate(g)) for g in reversed(A)]
    b_ok = len(A) == len(B) and all(key(x) == key(y) for x, y in zip(invA_S, B))
    # (c) A = 3 identical Trotter steps
    t = len(A) // 3
    steps = [A[i * t:(i + 1) * t] for i in range(3)]
    c_ok = len(A) % 3 == 0 and all(key(x) == key(y) for x, y in zip(steps[0], steps[1])) \
        and all(key(x) == key(y) for x, y in zip(steps[0], steps[2]))
    # V block: V = exp(-i delta X) on P?
    Vmat_ok = True
    for q in used:
        m = np.eye(2, dtype=complex)
        for g in Vb:
            if g["qubits"][0] == q:
                m = g["mat"] @ m
        ref = np.eye(2, dtype=complex)
        if q in P:
            d = 0.15
            ref = np.array([[math.cos(d), -1j * math.sin(d)], [-1j * math.sin(d), math.cos(d)]])
        # compare up to global phase
        ph = np.vdot(ref.ravel(), m.ravel()); ph /= abs(ph)
        if np.abs(m - ph * ref).max() > 1e-12:
            Vmat_ok = False
    if verbose:
        print("used qubits", len(used), used)
        print("edges", len(edges), "max degree", max(deg.values()))
        print("S (|rx|=0.5):", S)
        print("F (|rx|=3pi/8):", len(F))
        print("P (rz in V block):", len(P), P, "angles", Vang)
        print("P & S empty:", not set(P) & set(S))
        print("gates: first half", len(U), "V block", len(Vb), "second half", len(U2))
        print("(a) second half == gate-wise inverse of first half:", a_ok)
        print("(b) first half = A then B, B == inverse(A) with rx sign flipped on S:", b_ok)
        print("(c) A = 3 identical Trotter steps of", t, "gates:", c_ok)
        print("V block == exp(-i 0.15 X) on P, identity elsewhere (up to phase):", Vmat_ok)
        cz_step = sum(1 for g in steps[0] if g["name"] == "cz")
        print("CZ per Trotter step", cz_step)
    out.update(a_ok=a_ok, b_ok=b_ok, c_ok=c_ok, Vmat_ok=Vmat_ok, step=steps[0], Vb=Vb, U=U, U2=U2,
               all_gates=c["gates"])
    return out


def build_L(info, L):
    """C(L) = U(L)^dag V U(L), U(L) = [step^L] then [inverse(step^L) with S-flip]."""
    S = set(info["S"])
    A = [g for _ in range(L) for g in info["step"]]
    def flipS(g):
        if g["name"] == "rx" and g["qubits"][0] in S:
            p = (-g["params"][0],)
            return dict(g, params=p, mat=mat1("rx", p))
        return g
    B = [flipS(inverse_gate(g)) for g in reversed(A)]
    U = A + B
    U2 = [inverse_gate(g) for g in reversed(U)]
    return U + list(info["Vb"]) + U2


if __name__ == "__main__":
    info = analyse()
    C3 = build_L(info, 3)
    same = len(C3) == len(info["all_gates"]) and all(key(x) == key(y) for x, y in zip(C3, info["all_gates"]))
    print("rebuilt L=3 circuit identical to file:", same)
    for L in (1, 2):
        print("L =", L, "gates", len(build_L(info, L)))
