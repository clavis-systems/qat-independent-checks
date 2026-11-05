"""V4 autotest: Hamiltoniane casuali piccole con gli stessi tipi di termini (Z, X, XX+YY dentro i gruppi,
ZX dentro e tra gruppi). Il limite di limite_v4.py non deve mai superare il minimo esatto (diagonalizzazione densa)."""
import sys, os, subprocess, json, tempfile
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(7)
P = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0.]]), "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1., -1])}
def lab(n, d):
    return "".join(d.get(i, "I") for i in range(n))
def dense(terms, n):
    H = 0
    for l, c in terms:
        M = np.array([[1.]])
        for ch in l: M = np.kron(M, P[ch])
        H = H + c * M
    return np.real(H)
worst = []
for trial in range(12):
    groups = [[0], [1, 2, 3, 4], [5, 6, 7, 8]]; n = 9; terms = []
    for g in groups:
        for q in g:
            terms.append((lab(n, {q: "X"}), rng.normal()))
            if len(g) > 1: terms.append((lab(n, {q: "Z"}), rng.normal()))
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                c = rng.normal() * 0.5
                terms += [(lab(n, {g[i]: "X", g[j]: "X"}), c), (lab(n, {g[i]: "Y", g[j]: "Y"}), c)]
                if rng.random() < 0.5: terms.append((lab(n, {g[i]: "Z", g[j]: "X"}), rng.normal()))
    for a, b in [(1, 5), (2, 7), (8, 3), (4, 0), (6, 0)]:
        terms.append((lab(n, {a: "Z", b: "X"}), rng.choice([-1, 1]) * rng.uniform(0.2, 1.0)))
    d = tempfile.mkdtemp()
    np.save(os.path.join(d, "h.npy"), np.array([[l, str(complex(c))] for l, c in terms]))
    subprocess.run([sys.executable, "-I", os.path.join(here, "limite_v4.py"), os.path.join(d, "h.npy")],
                   cwd=d, check=True, capture_output=True)
    bound = json.load(open(os.path.join(d, "v4_best.json")))["best_bound"]
    exact = np.linalg.eigvalsh(dense(terms, n))[0]
    worst.append(exact - bound)
    print(f"prova {trial}: limite {bound:.10f}  minimo esatto {exact:.10f}  margine {exact - bound:.3e}")
print("margine minimo", min(worst), "-> OK" if min(worst) > -1e-10 else "-> ERRORE")
