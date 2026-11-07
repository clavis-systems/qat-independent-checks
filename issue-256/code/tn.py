"""Own network builder, exact circuit simplifications, slicing and contraction executor.

f = 2^-k Tr(O C^dag O C) = 2^-k sum_{a,b} s_a s_b C_ba conj(C_ba),  O = prod_{q in obs} Z_q.

Pipeline
 1. ops = merged gate list: ('1q', q, M) with M a 2x2 matrix (product of consecutive 1q gates;
    diagonal products are carried forward through CZs, which they commute with), ('cz', u, v).
 2. peephole(ops): exact cancellation of CZ pairs on the same edge separated only by diagonal
    operations on both wires (CZ^2 = I, diagonal gates commute); re-merge; repeat to fixpoint.
    1q products equal to the identity up to a global phase are dropped (f ignores global phase).
 3. lightcone(ops, obs): drop operations that commute with the Heisenberg-evolved observable
    (backward pass) and with the Schroedinger-evolved observable (forward pass, Tr(C O C^dag O)).
 4. build_network: doubled network (C and conj(C)), diagonal operations as hyperedge tensors
    (a qubit-wire index is shared by every diagonal operation on that wire segment).
 5. path: cotengra (path only), then own greedy slicing and own pairwise contraction executor
    (batch-matmul, hyperindices kept as batch dimensions until their last use).
"""
import math
import itertools
import numpy as np

CZT = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.complex128)   # CZ as 2-index (hyper) tensor


def is_diag(M, tol=1e-14):
    return abs(M[0, 1]) < tol and abs(M[1, 0]) < tol


def is_id_phase(M, tol=1e-13):
    if not is_diag(M, tol):
        return False
    return abs(M[0, 0] - M[1, 1]) < tol


# ---------------------------------------------------------------- 1. merge
def merge(ops_or_gates):
    """Accepts parsed gate dicts or ops; returns merged ops list."""
    ops = []
    for g in ops_or_gates:
        if isinstance(g, dict):
            if g["name"] == "cz":
                ops.append(("cz",) + tuple(g["qubits"]))
            else:
                ops.append(("1q", g["qubits"][0], g["mat"]))
        else:
            ops.append(g)
    pend = {}
    out = []

    def flush(q):
        M = pend.pop(q, None)
        if M is not None and not is_id_phase(M):
            out.append(("1q", q, M))

    for op in ops:
        if op[0] == "1q":
            q, M = op[1], op[2]
            pend[q] = M @ pend[q] if q in pend else M.copy()
        else:
            for q in op[1:]:
                if q in pend and not is_diag(pend[q]):
                    flush(q)           # non-diagonal must stay before the CZ
            out.append(op)             # diagonal pending products commute with CZ
    for q in list(pend):
        flush(q)
    return out


# ---------------------------------------------------------------- 2. peephole
def _cancel_pass(ops):
    wires = {}
    for i, op in enumerate(ops):
        for q in (op[1:] if op[0] == "cz" else (op[1],)):
            wires.setdefault(q, []).append(i)
    pos = {}
    for q, lst in wires.items():
        for k, i in enumerate(lst):
            pos[(q, i)] = k
    removed = set()
    for i, op in enumerate(ops):
        if op[0] != "cz" or i in removed:
            continue
        u, v = op[1], op[2]
        edge = frozenset((u, v))
        partner = []
        for q in (u, v):
            lst = wires[q]
            found = None
            for j in lst[pos[(q, i)] + 1:]:
                if j in removed:
                    continue
                o = ops[j]
                if o[0] == "cz":
                    if frozenset(o[1:]) == edge:
                        found = j
                        break
                    continue           # other CZ: diagonal, commutes
                if is_diag(o[2]):
                    continue           # diagonal 1q: commutes
                break                  # non-diagonal 1q blocks
            partner.append(found)
        if partner[0] is not None and partner[0] == partner[1]:
            removed.add(i)
            removed.add(partner[0])
    return [op for i, op in enumerate(ops) if i not in removed], len(removed) // 2


def peephole(ops, verbose=False):
    ops = merge(ops)
    it = 0
    while True:
        ops, n = _cancel_pass(ops)
        ops = merge(ops)
        it += 1
        if verbose:
            print(f"  peephole iter {it}: cancelled {n} CZ pairs, ops now {len(ops)}")
        if n == 0:
            return ops


# ---------------------------------------------------------------- 3. light cone
def _lc_pass(ops, obs):
    supp = set(obs)                  # support of evolved observable
    nd = set()                       # qubits where it may fail to commute with Z_q
    keep = []
    for op in ops:
        if op[0] == "cz":
            qs = op[1:]
            if all(q not in nd for q in qs):
                continue             # diagonal op commuting with X
            supp.update(qs)
            keep.append(op)
        else:
            q, M = op[1], op[2]
            if is_diag(M):
                if q not in nd:
                    continue
                keep.append(op)
                continue
            if q not in supp:
                continue
            nd.add(q)
            keep.append(op)
    return keep


def lightcone(ops, obs):
    while True:
        n0 = len(ops)
        ops = _lc_pass(ops[::-1], obs)[::-1]     # backward: X = C^dag O C
        ops = _lc_pass(ops, obs)                 # forward:  Y = C O C^dag
        if len(ops) == n0:
            return ops


# ---------------------------------------------------------------- 4. network
def build_network(ops, obs):
    """Return (arrays, inputs, size_dict, k). Contraction of all = Tr(O C^dag O C)."""
    qubits = sorted({q for op in ops for q in (op[1:] if op[0] == "cz" else (op[1],))} | set(obs))
    cnt = itertools.count()
    arrays, inputs = [], []
    start = {q: next(cnt) for q in qubits}
    sgn = np.array([1.0, -1.0], dtype=np.complex128)
    for q in obs:
        arrays.append(sgn); inputs.append((start[q],))
    ends = []
    for conj in (False, True):
        cur = dict(start)
        for op in ops:
            if op[0] == "cz":
                arrays.append(CZT); inputs.append((cur[op[1]], cur[op[2]]))
            else:
                q, M = op[1], (op[2].conj() if conj else op[2])
                if is_diag(M):
                    arrays.append(np.diag(M).copy()); inputs.append((cur[q],))
                else:
                    n = next(cnt)
                    arrays.append(M.copy()); inputs.append((n, cur[q]))   # M[out, in]
                    cur[q] = n
        ends.append(cur)
    # close: ket output index == bra output index (sum over b), with sign on obs
    ren = {}
    for q in qubits:
        ren[ends[1][q]] = ends[0][q]
    inputs = [tuple(ren.get(i, i) for i in t) for t in inputs]
    for q in obs:
        arrays.append(sgn); inputs.append((ends[0][q],))
    size = {i: 2 for t in inputs for i in t}
    return arrays, inputs, size, len(qubits)


def absorb_vectors(arrays, inputs):
    """Multiply every 1-index tensor into another tensor that carries the same index (exact)."""
    arrays, inputs = list(arrays), list(inputs)
    scalar = 1.0 + 0j
    while True:
        vec = next((i for i, t in enumerate(inputs) if len(t) == 1), None)
        if vec is None:
            return arrays, inputs, scalar
        ix = inputs[vec][0]
        tgt = next((j for j, t in enumerate(inputs) if j != vec and ix in t), None)
        if tgt is None:                      # index used by this vector only: sum it out
            scalar *= complex(arrays[vec].sum())
            del arrays[vec]; del inputs[vec]
            continue
        t = inputs[tgt]
        ax = t.index(ix)
        shape = [1] * len(t); shape[ax] = 2
        arrays[tgt] = arrays[tgt] * arrays[vec].reshape(shape)
        del arrays[vec]; del inputs[vec]


# ---------------------------------------------------------------- 5. cost, slicing, executor
def ssa_from_linear(path, n):
    ids = list(range(n))
    nxt = n
    ssa = []
    for a, b in path:
        a, b = sorted((a, b))
        x, y = ids[a], ids[b]
        del ids[b]; del ids[a]
        ssa.append((x, y))
        ids.append(nxt); nxt += 1
    return ssa


def tree_sets(inputs, ssa, sliced=frozenset()):
    """Index sets of every intermediate and per-step info (own bookkeeping)."""
    sets = [frozenset(t) - sliced for t in inputs]
    from collections import Counter
    count = Counter(i for s in sets for i in s)
    steps = []
    for x, y in ssa:
        A, B = sets[x], sets[y]
        for i in A: count[i] -= 1
        for i in B: count[i] -= 1
        U = A | B
        keep = frozenset(i for i in U if count[i] > 0)
        for i in keep: count[i] += 1
        steps.append((x, y, U, keep))
        sets.append(keep)
    return sets, steps


def cost(inputs, ssa, sliced=frozenset()):
    _, steps = tree_sets(inputs, ssa, sliced)
    flops = sum(2.0 ** len(U) for _, _, U, _ in steps)
    maxsz = max(2.0 ** len(K) for _, _, _, K in steps)
    return flops, maxsz


def choose_slices(inputs, ssa, target_log2, max_slices=30):
    sliced = frozenset()
    allix = {i for t in inputs for i in t}
    while True:
        fl, mx = cost(inputs, ssa, sliced)
        if mx <= 2.0 ** target_log2 or len(sliced) >= max_slices:
            return sliced, fl, mx
        best = None
        _, steps = tree_sets(inputs, ssa, sliced)
        cands = set()
        for _, _, U, K in steps:
            if 2.0 ** len(K) >= mx:
                cands |= K
        for i in cands:
            f2, m2 = cost(inputs, ssa, sliced | {i})
            score = (m2, f2)
            if best is None or score < best[0]:
                best = (score, i)
        sliced = sliced | {best[1]}


def pair_contract(A, ia, B, ib, keep):
    """Contract A[ia] B[ib]; indices in `keep` survive (batch if shared). Own implementation."""
    sa, sb = set(ia), set(ib)
    # sum out indices private to one operand and not needed later
    for arr_name in ("A", "B"):
        arr, ix, other = (A, ia, sb) if arr_name == "A" else (B, ib, sa)
        drop = [k for k, i in enumerate(ix) if i not in other and i not in keep]
        if drop:
            arr = arr.sum(axis=tuple(drop))
            ix = tuple(i for i in ix if i in other or i in keep)
        if arr_name == "A":
            A, ia = arr, ix
        else:
            B, ib = arr, ix
    sa, sb = set(ia), set(ib)
    shared = [i for i in ia if i in sb]
    batch = [i for i in shared if i in keep]
    summ = [i for i in shared if i not in keep]
    fa = [i for i in ia if i not in sb]
    fb = [i for i in ib if i not in sa]
    pa = [ia.index(i) for i in batch + fa + summ]
    pb = [ib.index(i) for i in batch + summ + fb]
    nb, nfa, ns, nfb = 2 ** len(batch), 2 ** len(fa), 2 ** len(summ), 2 ** len(fb)
    Am = np.ascontiguousarray(A.transpose(pa)).reshape(nb, nfa, ns)
    Bm = np.ascontiguousarray(B.transpose(pb)).reshape(nb, ns, nfb)
    C = np.matmul(Am, Bm)
    out_ix = tuple(batch + fa + fb)
    return C.reshape((2,) * len(out_ix)), out_ix


def contract_slice(arrays, inputs, ssa, sliced_vals):
    """sliced_vals: dict index -> 0/1. Returns scalar."""
    arrs, ixs = [], []
    for a, t in zip(arrays, inputs):
        sel = tuple(sliced_vals[i] if i in sliced_vals else slice(None) for i in t)
        arrs.append(a[sel])
        ixs.append(tuple(i for i in t if i not in sliced_vals))
    sliced = frozenset(sliced_vals)
    sets, steps = tree_sets(inputs, ssa, sliced)
    for (x, y, U, K) in steps:
        C, cix = pair_contract(arrs[x], ixs[x], arrs[y], ixs[y], K)
        assert set(cix) == set(K)
        arrs[x] = arrs[y] = None
        arrs.append(C); ixs.append(cix)
    res = arrs[-1]
    assert ixs[-1] == ()
    return complex(res)


def contract(arrays, inputs, ssa, sliced, log=None):
    sliced = sorted(sliced)
    total = 0j
    for n, bits in enumerate(itertools.product((0, 1), repeat=len(sliced))):
        total += contract_slice(arrays, inputs, ssa, dict(zip(sliced, bits)))
        if log is not None and (n + 1) % max(1, 2 ** len(sliced) // 16) == 0:
            log(n + 1, 2 ** len(sliced), total)
    return total
