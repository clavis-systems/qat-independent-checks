"""Validation on reduced circuits small enough for dense unitaries.

1. qiskit cross-parse: gate list from qiskit.qasm3 == own parser's gate list (names, qubits, angles).
2. For several qubit subsets Q (k <= 12; all gates of the full L=3 circuit acting inside Q),
   compare f from: own dense simulator (raw gates), qiskit Operator (independent parse+sim),
   own TN with and without the peephole/light-cone passes, different paths and forced slicing.
"""
import sys, time, math, resource, platform
from collections import deque
import numpy as np
import qasm_parse, structure, dense, pipeline

CIRC = structure.CIRC
OBS = structure.OBS


def ball(edges, centers, k):
    adj = {}
    for u, v in edges:
        adj.setdefault(u, set()).add(v); adj.setdefault(v, set()).add(u)
    seen, dq = [], deque(centers)
    while dq and len(seen) < k:
        q = dq.popleft()
        if q in seen:
            continue
        seen.append(q)
        for w in sorted(adj[q]):
            if w not in seen:
                dq.append(w)
    return sorted(seen)


def main():
    import cotengra, opt_einsum, qiskit
    print("python", platform.python_version(), "numpy", np.__version__, "cotengra", cotengra.__version__,
          "opt_einsum", opt_einsum.__version__, "qiskit", qiskit.__version__)
    info = structure.analyse(verbose=False)
    gates = info["all_gates"]
    # 1. cross-parse
    ql = dense.qiskit_gate_list(CIRC)
    mine = [(g["name"], g["qubits"], g["params"]) for g in gates]
    same = len(ql) == len(mine) and all(a[0] == b[0] and a[1] == b[1] and
                                          len(a[2]) == len(b[2]) and
                                          all(abs(x - y) < 1e-15 for x, y in zip(a[2], b[2]))
                                          for a, b in zip(ql, mine))
    print(f"qiskit cross-parse: {len(ql)} gates, identical to own parser: {same}")
    edges = info["edges"]
    S, P = info["S"], info["P"]
    # hand-picked connected patches mixing O sites, S sites (rx +-0.5) and P sites (perturbed)
    cases = [
        ("ball(O),k=10 [trivial]", ball(edges, list(OBS), 10)),
        # (k=12 patches were too slow for the dense references under the 4-CPU limit; k=10/11 used)
        ("O52+S39,50+P", sorted([33, 38, 39, 48, 49, 50, 51, 52, 53, 54])),
        ("O52,72+S50,69,70+P", sorted([49, 50, 51, 52, 58, 69, 70, 71, 72, 73])),
        ("O52,59+S39,50+P", sorted([33, 34, 39, 50, 51, 52, 53, 54, 55, 59])),
        ("O72+S57,67-70+P", sorted([47, 48, 57, 58, 67, 68, 69, 70, 71, 72])),
        ("S7-10,17+P, obs=[9]", sorted([7, 8, 9, 10, 11, 12, 17, 18, 26, 27])),
        ("O52+S39,50+P, k=11", sorted([29, 33, 38, 39, 48, 49, 50, 51, 52, 53, 54])),
    ]
    maxdev = 0.0
    if len(sys.argv) > 1:            # optional: run only the listed case numbers (1-based)
        keep = {int(x) for x in sys.argv[1].split(",")}
        cases = [c for i, c in enumerate(cases, 1) if i in keep]
    for name, Q in cases:
        obsQ = [q for q in OBS if q in Q]
        if not obsQ:
            obsQ = [9]           # patch without O sites: use Z_9 as observable
        sub = [g for g in gates if all(q in Q for q in g["qubits"])]
        print(f"\n[{name}] Q={Q} obs={obsQ} gates={len(sub)}", flush=True)
        t0 = time.time(); f_own = dense.own_dense_f(sub, Q, obsQ); t_own = time.time() - t0
        print(f"  own dense      f = {f_own:.15f}   ({t_own:.1f}s)", flush=True)
        t0 = time.time(); f_qk, nqk = dense.qiskit_f(CIRC, Q, obsQ); t_qk = time.time() - t0
        print(f"  qiskit         f = {f_qk:.15f}   ({t_qk:.1f}s, {nqk} gates)", flush=True)
        devs = {"qiskit": abs(f_qk - f_own)}
        variants = [
            ("tn_full", dict(peep=True, lc=True, target_log2=24)),
            ("tn_sliced", dict(peep=True, lc=True, target_log2=-4, seed=1, force_slices=3,
                               max_sliced=5)),
            ("tn_oe", dict(peep=True, lc=True, target_log2=24, method="oe-random-greedy",
                           max_repeats=8, max_time=10)),
            ("tn_lc_only", dict(peep=False, lc=True, target_log2=24, ctg_methods=("kahypar", "greedy"),
                                max_time=20, max_log2_flops=33)),
            ("tn_peep_only", dict(peep=True, lc=False, target_log2=24, ctg_methods=("kahypar", "greedy"),
                                  max_time=20, max_log2_flops=33)),
            ("tn_raw", dict(peep=False, lc=False, target_log2=24, ctg_methods=("kahypar", "greedy"),
                            max_time=20, max_log2_flops=33)),
        ]
        for vname, kw in variants:
            t0 = time.time()
            v, st = pipeline.value(sub, obsQ, **kw)
            dt = time.time() - t0
            if v != v:      # nan: skipped by cost guard
                print(f"  {vname:14s} skipped: path cost 2^{st.get('log2_flops_total')} > 2^33", flush=True)
                continue
            devs[vname] = abs(v - f_own)
            print(f"  {vname:14s} f = {v.real:.15f} {v.imag:+.1e}j  |dev|={abs(v - f_own):.1e}  "
                  f"({dt:.1f}s; ops {st['n_lightcone']}, tensors {st['n_tensors']}, "
                  f"slices {st.get('n_sliced')}, log2 flops {st.get('log2_flops_total')})", flush=True)
        m = max(devs.values()); maxdev = max(maxdev, m)
        print(f"  max |dev| this case {m:.2e}", flush=True)
    print(f"\nMAX |deviation| over all cases and methods: {maxdev:.3e}")
    print("peak RSS MB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)


if __name__ == "__main__":
    main()
