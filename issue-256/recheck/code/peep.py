"""Own peephole: merge 1q gates, drop identities (tracking phase), cancel CZ pairs separated only by diagonal ops.
Ops: ('1', q, mat) or ('cz', a, b).  Exact."""
import numpy as np
from ole import *
from tnbuild import build_ops

TOL = 1e-13

def to_ops(seq):
    ops = []
    for kind, x in seq:
        if kind == "V":
            for t, m in x.items(): ops.append(["1", t, m])
        elif kind == "O":
            for t in x: ops.append(["1", t, Z.copy()])
        else:
            n, a, q = x
            if n == "cz": ops.append(["cz", q[0], q[1]])
            else: ops.append(["1", q[0], mat1(n, a)])
    return ops

def is_diag(m): return abs(m[0, 1]) < TOL and abs(m[1, 0]) < TOL
def is_ident(m): return is_diag(m) and abs(m[0, 0] - m[1, 1]) < TOL

def wires_of(op): return (op[1],) if op[0] == "1" else (op[1], op[2])

def peephole(ops):
    ops = [list(o) for o in ops]
    scalar = 1.0 + 0j
    while True:
        changed = False
        # merge 1q on same wire; a pending DIAGONAL op is carried forward through CZs (commutes)
        pend = {}      # wire -> index in out of pending 1q op (still movable)
        out = []
        for o in ops:
            if o[0] == "1":
                t = o[1]
                if t in pend:
                    p = out[pend.pop(t)]
                    m = o[2] @ p[2]
                    p[0] = "x"            # mark dead
                    out.append(["1", t, m]); changed = True
                else:
                    out.append(o)
                pend[t] = len(out) - 1
            else:
                for w in (o[1], o[2]):
                    if w in pend and not is_diag(out[pend[w]][2]): pend.pop(w)
                out.append(o)
        ops = [o for o in out if o[0] != "x"]
        # drop identities
        out = []
        for o in ops:
            if o[0] == "1" and is_ident(o[2]):
                scalar *= o[2][0, 0]; changed = True
            else: out.append(o)
        ops = out
        # cancel CZ pairs
        dead = set()
        perwire = {}
        for i, o in enumerate(ops):
            for w in wires_of(o): perwire.setdefault(w, []).append(i)
        for i, o in enumerate(ops):
            if o[0] != "cz" or i in dead: continue
            a, b = o[1], o[2]
            # scan forward on wire a to find partner
            la = perwire[a]; lb = perwire[b]
            ia = la.index(i); ib = lb.index(i)
            j = None
            k = ia + 1
            while k < len(la):
                oj = ops[la[k]]
                if la[k] in dead: k += 1; continue
                if oj[0] == "cz":
                    if set((oj[1], oj[2])) == {a, b}: j = la[k]; break
                    k += 1; continue
                if is_diag(oj[2]): k += 1; continue
                break
            if j is None: continue
            # check wire b between i and j
            ok = True
            k = ib + 1
            while lb[k] != j:
                oj = ops[lb[k]]
                if lb[k] in dead: k += 1; continue
                if oj[0] == "1" and not is_diag(oj[2]): ok = False; break
                k += 1
            if ok:
                dead.add(i); dead.add(j); changed = True
        ops = [o for i, o in enumerate(ops) if i not in dead]
        if not changed: break
    return ops, scalar


def tn_from_ops(ops, wires, dtype="complex128"):
    import quimb.tensor as qtn
    cnt = [0]
    def new():
        cnt[0] += 1; return "i%d" % cnt[0]
    start = {t: new() for t in wires}; cur = dict(start); tensors = []
    for o in ops:
        if o[0] == "1":
            t = o[1]; n = new()
            tensors.append(qtn.Tensor(o[2].astype(dtype), inds=(n, cur[t]))); cur[t] = n
        else:
            a, b = o[1], o[2]; n0, n1 = new(), new()
            tensors.append(qtn.Tensor(CZ.astype(dtype), inds=(n0, n1, cur[a], cur[b]))); cur[a], cur[b] = n0, n1
    tn = qtn.TensorNetwork(tensors)
    for t in wires:
        if cur[t] == start[t]: raise RuntimeError("bare wire %d" % t)
        tn.reindex_({cur[t]: start[t]})
    return tn


def build_peep(G, Vmats, Oq, wires):
    ops, sc = peephole(to_ops(build_ops(G, {t: m for t, m in Vmats.items() if t in wires}, Oq)))
    return ops, sc


def f_peep(G, Vmats, Oq, wires, optimize="auto"):
    ops, sc = build_peep(G, Vmats, Oq, wires)
    used = sorted({w for o in ops for w in wires_of(o)})
    if not ops: return complex(sc)
    tn = tn_from_ops(ops, used)
    tn2 = tn.copy(); tn2.full_simplify_(seq="ADCRS", atol=1e-12)
    if tn2.num_tensors == 0: tn2 = tn
    val = tn2.contract(all, optimize=optimize, output_inds=())
    return complex(val) * sc / 2.0 ** len(used)


# ---- v2: peephole on U alone first, then cone on the reduced op list, then peephole on the whole sequence ----
def ops_dagger(ops):
    return [[o[0], o[1], o[2].conj().T] if o[0] == "1" else list(o) for o in reversed(ops)]

def ops_forward_cone(ops, Oq):
    sup = set(Oq); keep = []
    for o in ops:
        if sup & set(wires_of(o)):
            keep.append(o)
            if o[0] == "cz": sup |= {o[1], o[2]}
    return keep, sup

def ops_backward_prune(ops, Pq):
    q = set(Pq); rev = []
    for o in reversed(ops):
        if q & set(wires_of(o)):
            rev.append(o)
            if o[0] == "cz": q |= {o[1], o[2]}
    return list(reversed(rev)), q

def build_v2(A, Vmats, Oq, Pq):
    Aops, _ = peephole(to_ops([("g", x) for x in A]))
    G, _ = ops_forward_cone_nd(Aops, Oq)
    G, q = ops_backward_prune(G, Pq)
    wires = sorted(q | set(Oq))
    Gd = ops_dagger(G)
    Vm = {t: m for t, m in Vmats.items() if t in wires}
    Vd = {t: m.conj().T for t, m in Vm.items()}
    seq = ([["1", t, m] for t, m in Vm.items()] + Gd + [["1", t, Z.copy()] for t in Oq] + [list(o) for o in G]
           + [["1", t, m] for t, m in Vd.items()] + ops_dagger(G) + [["1", t, Z.copy()] for t in Oq] + [list(o) for o in G])
    ops, sc = peephole(seq)
    return ops, sc, wires, len(Aops), len(G)

def f_v2(A, Vmats, Oq, Pq, optimize="auto"):
    ops, sc, wires, _, _ = build_v2(A, Vmats, Oq, Pq)
    used = sorted({w for o in ops for w in wires_of(o)})
    if not ops: return complex(sc)
    tn = tn_from_ops(ops, used)
    tn2 = tn.copy(); tn2.full_simplify_(seq="ADCRS", atol=1e-12)
    if tn2.num_tensors == 0: tn2 = tn
    return complex(tn2.contract(all, optimize=optimize, output_inds=())) * sc / 2.0 ** len(used)


def ops_forward_cone_nd(ops, Oq):
    """Finer exact forward cone for W = U O U^dag with O diagonal (Z-type): tracks support and 'non-diagonal' wires."""
    supp = set(Oq); nd = set(); keep = []
    for o in ops:
        if o[0] == "cz":
            if o[1] not in nd and o[2] not in nd: continue      # commutes with a Z-type operator there
            supp |= {o[1], o[2]}; keep.append(o)
        else:
            t, m = o[1], o[2]
            if is_diag(m):
                if t not in nd: continue
                keep.append(o)
            else:
                if t not in supp: continue
                nd.add(t); keep.append(o)
    return keep, supp
