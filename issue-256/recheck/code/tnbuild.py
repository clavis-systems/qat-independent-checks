"""Build the tensor network for f = 2^-n Tr(W V^dag W V), W = U O U^dag, with exact light-cone pruning."""
import numpy as np, quimb.tensor as qtn
from ole import *


def build_ops(G, Vmats, Oq):
    """op sequence in time order (first applied first)."""
    Vd = {t: m.conj().T for t, m in Vmats.items()}
    ops = []
    ops.append(("V", Vmats))
    ops += [("g", x) for x in dag_list(G)]
    ops.append(("O", Oq))
    ops += [("g", x) for x in G]
    ops.append(("V", Vd))
    ops += [("g", x) for x in dag_list(G)]
    ops.append(("O", Oq))
    ops += [("g", x) for x in G]
    return ops


def make_tn(G, Vmats, Oq, wires, dtype="complex128"):
    ops = build_ops(G, Vmats, Oq)
    cnt = [0]
    def new():
        cnt[0] += 1; return "i%d" % cnt[0]
    start = {t: new() for t in wires}
    cur = dict(start)
    pend = {}
    tensors = []
    def flush(t):
        if t in pend:
            o = new()
            tensors.append(qtn.Tensor(pend.pop(t).astype(dtype), inds=(o, cur[t])))
            cur[t] = o
    for kind, x in ops:
        if kind == "V":
            for t, m in x.items():
                if t in cur: pend[t] = m @ pend.get(t, I2)
        elif kind == "O":
            for t in x: pend[t] = Z @ pend.get(t, I2)
        else:
            n, a, q = x
            if n == "cz":
                for t in q: flush(t)
                o0, o1 = new(), new()
                tensors.append(qtn.Tensor(CZ.astype(dtype), inds=(o0, o1, cur[q[0]], cur[q[1]])))
                cur[q[0]], cur[q[1]] = o0, o1
            else:
                pend[q[0]] = mat1(n, a) @ pend.get(q[0], I2)
    for t in list(pend): flush(t)
    # close the trace: rename final index of each wire to the start index
    tn = qtn.TensorNetwork(tensors)
    for t in wires:
        tn.reindex_({cur[t]: start[t]})
    return tn


def prepare(A, Vmats, Oq, Pq, prune=True):
    """returns G (pruned gate list), wires"""
    if prune:
        G, _ = forward_cone(A, Oq)
        G, q = backward_prune(G, Pq)
        wires = sorted(set(q) | set(Oq))
    else:
        G = list(A)
        wires = sorted({t for x in A for t in x[2]} | set(Oq) | set(Pq))
    return G, wires


def f_from_tn(tn, nwires, optimize="auto-hq"):
    val = tn.contract(all, optimize=optimize, output_inds=())
    return complex(val) / 2 ** nwires
