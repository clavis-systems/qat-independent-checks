"""Lower bound for the 49-qubit guided sparse problem (tracker issue #240).

Written by AI agents (Claude) directed by Emanuele Rizzan, 4 Nov 2026. Apache-2.0.
It reads only the official Pauli list of the tracker. It does not use the submitter's code.

Idea (operator inequality, no assumption on the ground state):
  every term that couples two different patches has the form  c * Z_a X_b  (a, b in different patches).
  For a sign s in {+1,-1} and kappa >= 1 the two-qubit operator
      B = |c| (I + s Z_a) (kappa I + t X_b),   t = s * sign(c)
  is positive semidefinite (product of two PSD factors on different qubits), and
      c Z_a X_b = B - |c| kappa I - |c| t X_b - |c| s kappa Z_a .
  So H = sum_p R_p + C + sum B  with R_p acting on one patch only, C a constant, B >= 0, hence
      lambda_min(H) >= sum_p lambda_min(R_p) + C          (rigorous lower bound).
  Each R_p acts on 16 qubits (or on the single qubit 48) and is diagonalised in its own 2^16 space.
  The bound is attained if the product of the patch ground states is annihilated by every B,
  which the script checks (sum of <B> over that product state, computed exactly from local expectations).

Checks done by the script:
  - the identity H = sum_p R_p + C + sum B is verified term by term on the Pauli expansion;
  - every B is checked to be PSD (4x4 eigenvalues);
  - the patch minima are computed twice (ARPACK Lanczos and LOBPCG, different random starts);
  - self-test (--selftest): on random 10-qubit instances with the same structure, the bound never exceeds the exact
    minimum from dense diagonalisation of the whole Hamiltonian.

Usage:  python -I lower_bound_240.py 49Q_v1_hamiltonian.npy OUT.json [kappa]
        python -I lower_bound_240.py --selftest
"""
import hashlib
import json
import resource
import sys
import time
from collections import defaultdict

import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as sla

E_ISSUE = -85.586625069864652


def parse(path):
    arr = np.load(path, allow_pickle=False)
    n = len(str(arr[0][0]))
    terms = []
    for p, c in arr:
        p = str(p)
        cc = complex(str(c).replace(" ", ""))
        assert abs(cc.imag) == 0.0
        # Qiskit order: rightmost character is qubit 0
        sup = {k: ch for k, ch in enumerate(reversed(p)) if ch != "I"}
        terms.append((sup, cc.real))
    return n, terms


def patches(n, terms):
    par = list(range(n))

    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a

    for sup, _ in terms:
        if len(sup) == 2 and set(sup.values()) in ({"X"}, {"Y"}):
            a, b = sorted(sup)
            par[f(a)] = f(b)
    comp = defaultdict(list)
    for a in range(n):
        comp[f(a)].append(a)
    return sorted((sorted(v) for v in comp.values()), key=lambda v: (-len(v), v[0]))


def pauli_sparse(local_terms, qubits):
    """Real sparse matrix of sum_k c_k P_k on the listed qubits (bit j of the index = qubit qubits[j])."""
    m = len(qubits)
    pos = {q: j for j, q in enumerate(qubits)}
    dim = 1 << m
    idx = np.arange(dim, dtype=np.int64)
    diag = np.zeros(dim)
    rows, cols, vals = [], [], []
    for sup, c in local_terms:
        xm = zm = 0
        ny = 0
        for q, ch in sup.items():
            bit = 1 << pos[q]
            if ch in "XY":
                xm |= bit
            if ch in "ZY":
                zm |= bit
            ny += ch == "Y"
        # P|i> = i^ny (-1)^popcount(i & zm) |i ^ xm>
        ph = (1j) ** ny
        assert abs(ph.imag) < 1e-15, "odd number of Y: not real"
        sign = 1.0 - 2.0 * (np.bitwise_count(idx & zm) & 1)
        v = c * ph.real * sign
        if xm == 0:
            diag += v
        else:
            rows.append(idx ^ xm)
            cols.append(idx)
            vals.append(v)
    M = sp.diags(diag, format="csr")
    if rows:
        M = M + sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim))
    return M


def lowest(M, seed=0, k=2):
    """Lowest eigenpairs with ARPACK (Lanczos); returns values, vectors and residual norms."""
    if M.shape[0] <= 64:
        w, v = np.linalg.eigh(M.toarray())
        return w[:k], v[:, :k], np.zeros(k)
    rng = np.random.default_rng(seed)
    w, v = sla.eigsh(M, k=k, which="SA", v0=rng.standard_normal(M.shape[0]), tol=1e-14, ncv=40)
    o = np.argsort(w)
    w, v = w[o], v[:, o]
    res = np.array([np.linalg.norm(M @ v[:, j] - w[j] * v[:, j]) for j in range(k)])
    return w, v, res


def lowest_lobpcg(M, seed=1):
    """Second, different solver for the patch minimum (LOBPCG, random block of 4 vectors)."""
    if M.shape[0] <= 64:
        return float(np.linalg.eigvalsh(M.toarray())[0])
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((M.shape[0], 4))
    w, _ = sla.lobpcg(M, X, largest=False, tol=1e-12, maxiter=2000)
    return float(np.min(w))


def local_expect(vec, qubits, q, ch):
    """<vec| P_q |vec> for a single-qubit Pauli P on qubit q of the patch."""
    M = pauli_sparse([({q: ch}, 1.0)], qubits)
    return float(vec @ (M @ vec))


def build(n, terms, groups, signs, kappa):
    gid = {q: g for g, qs in enumerate(groups) for q in qs}
    per = defaultdict(list)
    inter = []
    for sup, c in terms:
        gs = {gid[q] for q in sup}
        if len(gs) == 1:
            per[gs.pop()].append((sup, c))
        else:
            inter.append((sup, c))
    const = 0.0
    bonds = []
    for e, (sup, c) in enumerate(inter):
        assert len(sup) == 2 and sorted(sup.values()) == ["X", "Z"], sup
        a = [q for q in sup if sup[q] == "Z"][0]
        b = [q for q in sup if sup[q] == "X"][0]
        s = signs[e]
        t = s * np.sign(c)
        const -= abs(c) * kappa
        per[gid[b]].append(({b: "X"}, -abs(c) * t))
        per[gid[a]].append(({a: "Z"}, -abs(c) * s * kappa))
        bonds.append(dict(a=a, b=b, c=c, s=int(s), t=int(t)))
    return per, inter, bonds, const


def check_identity(terms, per, bonds, const, kappa):
    """Pauli expansion of sum_p R_p + C + sum B must equal the original list, term by term."""
    def key(sup):
        return tuple(sorted(sup.items()))

    orig, rebuilt = defaultdict(float), defaultdict(float)
    for sup, c in terms:
        orig[key(sup)] += c
    for lst in per.values():
        for sup, c in lst:
            rebuilt[key(sup)] += c
    rebuilt[()] += const
    for bd in bonds:
        a, b, c, s, t = bd["a"], bd["b"], bd["c"], bd["s"], bd["t"]
        w = abs(c)
        rebuilt[()] += w * kappa
        rebuilt[key({b: "X"})] += w * t
        rebuilt[key({a: "Z"})] += w * s * kappa
        rebuilt[key({a: "Z", b: "X"})] += w * s * t
    keys = set(orig) | set(rebuilt)
    return max(abs(orig.get(k, 0.0) - rebuilt.get(k, 0.0)) for k in keys)


def bond_min_eig(bd, kappa):
    Z = np.diag([1.0, -1.0])
    X = np.array([[0.0, 1.0], [1.0, 0.0]])
    Bm = abs(bd["c"]) * np.kron(np.eye(2) + bd["s"] * Z, kappa * np.eye(2) + bd["t"] * X)
    return float(np.linalg.eigvalsh(Bm)[0])


def bound(n, terms, groups, signs, kappa, want_vectors=False):
    per, inter, bonds, const = build(n, terms, groups, signs, kappa)
    lam, vecs, res, gap = [], [], [], []
    for g, qs in enumerate(groups):
        M = pauli_sparse(per[g], qs)
        w, v, r = lowest(M)
        lam.append(float(w[0]))
        gap.append(float(w[1] - w[0]) if len(w) > 1 else None)
        res.append(float(r[0]))
        vecs.append(v[:, 0])
    L = float(sum(lam) + const)
    out = dict(L=L, lam=lam, gap=gap, res=res, const=const, bonds=bonds)
    if want_vectors:
        out["vecs"] = vecs
    return out


def ascent(n, terms, groups, kappa, verbose=True):
    """Any sign choice gives a valid bound; coordinate ascent only makes it tighter."""
    gid = {q: g for g, qs in enumerate(groups) for q in qs}
    ninter = sum(1 for sup, _ in terms if len({gid[q] for q in sup}) > 1)
    signs = [1] * ninter
    best = bound(n, terms, groups, signs, kappa)
    if verbose:
        print("start signs", signs, "L", repr(best["L"]))
    improved, sweep = True, 0
    while improved and sweep < 4:
        improved = False
        sweep += 1
        for e in range(ninter):
            trial = list(signs)
            trial[e] = -trial[e]
            r = bound(n, terms, groups, trial, kappa)
            if r["L"] > best["L"] + 1e-12:
                signs, best, improved = trial, r, True
                if verbose:
                    print(f"  sweep {sweep} flip {e} -> L {r['L']!r}")
    return signs


def selftest(ntests=20, seed=7):
    """Random instances: 3 patches of 3 qubits + 1 isolated qubit, same term types as 49Q_v1."""
    rng = np.random.default_rng(seed)
    worst = -np.inf
    for it in range(ntests):
        groups_true = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9]]
        terms = []
        for g in groups_true[:3]:
            for i in range(3):
                for j in range(i + 1, 3):
                    c = rng.normal()
                    terms.append(({g[i]: "X", g[j]: "X"}, c))
                    terms.append(({g[i]: "Y", g[j]: "Y"}, c * rng.choice([-1, 1])))
                    if rng.random() < 0.5:
                        terms.append(({g[i]: "X", g[j]: "Z"}, rng.normal()))
            for q in g:
                terms.append(({q: "Z"}, rng.normal()))
                terms.append(({q: "X"}, rng.normal()))
        terms.append(({9: "X"}, rng.normal()))
        for (a, b) in [(1, 4), (2, 7), (5, 8), (0, 9), (6, 9)]:
            terms.append(({a: "Z", b: "X"}, rng.normal()))
            if b != 9:
                terms.append(({a: "X", b: "Z"}, rng.normal()))
        n = 10
        groups = patches(n, terms)
        signs = ascent(n, terms, groups, 1.0, verbose=False)
        L = bound(n, terms, groups, signs, 1.0)["L"]
        exact = float(np.linalg.eigvalsh(pauli_sparse(terms, list(range(n))).toarray())[0])
        worst = max(worst, L - exact)
        assert L <= exact + 1e-9, (it, L, exact)
    print(f"selftest: {ntests} random 10-qubit instances, max(bound - exact) = {worst:.3e} (must be <= 0)")
    return worst


def main():
    if sys.argv[1] == "--selftest":
        selftest()
        return
    path, outp = sys.argv[1], sys.argv[2]
    kappa = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    t0 = time.time()
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    n, terms = parse(path)
    groups = patches(n, terms)
    gid = {q: g for g, qs in enumerate(groups) for q in qs}
    inter = [(sup, c) for sup, c in terms if len({gid[q] for q in sup}) > 1]
    print("sha256", sha)
    print("qubits", n, "terms", len(terms), "patch sizes", [len(g) for g in groups], "inter-patch terms", len(inter))
    for sup, c in inter:
        print("  inter", {q: sup[q] for q in sorted(sup)}, "patches", sorted({gid[q] for q in sup}), c)

    # Step 1: choose the signs s (start from +1 everywhere, then coordinate ascent on the bound).
    signs = ascent(n, terms, groups, kappa)
    final = bound(n, terms, groups, signs, kappa, want_vectors=True)
    L = final["L"]

    # Step 2: algebra checks.
    per, _, bonds, const = build(n, terms, groups, signs, kappa)
    id_err = check_identity(terms, per, bonds, const, kappa)
    bmin = min(bond_min_eig(bd, kappa) for bd in bonds)

    # Step 3: is the bound attained? Product of the patch ground states: <B> from local expectations.
    vecs = final["vecs"]
    sumB = 0.0
    for bd in final["bonds"]:
        ga, gb = gid[bd["a"]], gid[bd["b"]]
        za = local_expect(vecs[ga], groups[ga], bd["a"], "Z")
        xb = local_expect(vecs[gb], groups[gb], bd["b"], "X")
        Bexp = abs(bd["c"]) * (1 + bd["s"] * za) * (kappa + bd["t"] * xb)
        bd.update(Za=za, Xb=xb, B_expect=Bexp)
        sumB += Bexp
    E_prod = L + sumB  # exact energy of that product state: sum lambda_p + C + sum <B>

    # Step 4: second solver on the final patch operators (LOBPCG, other random start).
    lam2 = [lowest_lobpcg(pauli_sparse(per[g], qs), seed=10 + g) for g, qs in enumerate(groups)]
    L2 = float(sum(lam2) + const)

    out = dict(
        file_sha256=sha,
        kappa=kappa,
        patch_sizes=[len(g) for g in groups],
        patches=groups,
        signs=signs,
        lower_bound=L,
        lower_bound_lobpcg=L2,
        patch_minima_eigsh=final["lam"],
        patch_minima_lobpcg=lam2,
        patch_gap_to_second=final["gap"],
        patch_residuals=final["res"],
        constant=final["const"],
        identity_max_coefficient_error=id_err,
        min_eigenvalue_of_bonds=bmin,
        bonds=final["bonds"],
        sum_B_on_product_state=sumB,
        product_state_energy=E_prod,
        issue_240_energy=E_ISSUE,
        lower_bound_minus_issue=L - E_ISSUE,
        seconds=time.time() - t0,
        max_rss_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        versions=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__),
    )
    print("identity H = sum R_p + C + sum B: max coefficient error", id_err, "| min eigenvalue of the bonds", bmin)
    print("patch minima (eigsh) ", [repr(x) for x in final["lam"]])
    print("patch minima (lobpcg)", [repr(x) for x in lam2])
    print("gap to 2nd level per patch", final["gap"], "| residuals", final["res"])
    print("constant", repr(final["const"]))
    print("LOWER BOUND", repr(L), "| with LOBPCG", repr(L2))
    print("sum <B> on product of patch ground states", repr(sumB), "-> energy of that state", repr(E_prod))
    print("issue #240", repr(E_ISSUE), "| lower bound - issue", L - E_ISSUE,
          "| seconds", round(out["seconds"], 1), "| max RSS MB", round(out["max_rss_MB"]))
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
