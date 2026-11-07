"""Own batch-matmul executor for a cotengra (sliced) tree on a hyper-index tensor network.
usage: mmexec.py <tree.pkl> <ckpt.jsonl> <n_slices_to_do|all>   (network rebuilt exactly as in run_full.py)"""
import sys, time, json, os, pickle, itertools, resource, numpy as np
from collections import Counter
from ole import *
from tnbuild import *
from peep import *

def build_network():
    A, M, B = load("/data/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm")
    V = extract_V(M, 49)
    Pq = sorted(x[2][0] for x in M if x[0] == "rz")
    Oq = [52, 59, 72]
    ops, sc, wires, nA, nG = build_v2(A, V, Oq, Pq)
    used = sorted({w for o in ops for w in wires_of(o)})
    tn = tn_from_ops(ops, used)
    tn.full_simplify_(seq="ADCRS", atol=1e-12)
    inputs, output, size_dict = tn.get_inputs_output_size_dict(output_inds=())
    return tn, [tuple(i) for i in inputs], sc, len(used)

def pair(Aa, ia, Bb, ib, keep):
    sa, sb = set(ia), set(ib)
    d = [k for k, i in enumerate(ia) if i not in sb and i not in keep]
    if d: Aa = Aa.sum(axis=tuple(d)); ia = tuple(i for i in ia if i in sb or i in keep)
    d = [k for k, i in enumerate(ib) if i not in set(ia) and i not in keep]
    if d: Bb = Bb.sum(axis=tuple(d)); ib = tuple(i for i in ib if i in set(ia) or i in keep)
    sa, sb = set(ia), set(ib)
    shared = [i for i in ia if i in sb]
    batch = [i for i in shared if i in keep]; summ = [i for i in shared if i not in keep]
    fa = [i for i in ia if i not in sb]; fb = [i for i in ib if i not in sa]
    pa = [ia.index(i) for i in batch + fa + summ]; pb = [ib.index(i) for i in batch + summ + fb]
    Am = np.ascontiguousarray(Aa.transpose(pa)).reshape(2 ** len(batch), 2 ** len(fa), 2 ** len(summ))
    Bm = np.ascontiguousarray(Bb.transpose(pb)).reshape(2 ** len(batch), 2 ** len(summ), 2 ** len(fb))
    out = tuple(batch + fa + fb)
    return np.matmul(Am, Bm).reshape((2,) * len(out)), out

STRIDE = int(sys.argv[4]) if len(sys.argv) > 4 else 1
OFF = int(sys.argv[5]) if len(sys.argv) > 5 else 0
def main():
    tree = pickle.load(open(sys.argv[1], "rb")); ck = sys.argv[2]
    tn, inputs, sc, nused = build_network()
    assert len(inputs) == len(tree.inputs)
    arrays = [t.data for t in tn.tensors]
    ssa = tree.get_ssa_path()
    sliced = sorted(tree.sliced_inds)
    print("tensors", len(inputs), "sliced", len(sliced), "steps", len(ssa), flush=True)
    # static bookkeeping: for each step which indices are kept
    base = [frozenset(t) - set(sliced) for t in inputs]
    cnt = Counter(i for s in base for i in s); sets = list(base); steps = []
    for x, y in ssa:
        a, b = sets[x], sets[y]
        for i in a: cnt[i] -= 1
        for i in b: cnt[i] -= 1
        keep = frozenset(i for i in a | b if cnt[i] > 0)
        for i in keep: cnt[i] += 1
        steps.append((x, y, keep)); sets.append(keep)
    done = {}
    if os.path.exists(ck):
        for l in open(ck):
            r = json.loads(l); done[r["i"]] = complex(r["re"], r["im"])
    nsl = 2 ** len(sliced); lim = nsl if sys.argv[3] == "all" else int(sys.argv[3])
    t0 = time.time(); tot = 0j
    for n, bits in enumerate(itertools.product((0, 1), repeat=len(sliced))):
        if n % STRIDE != OFF: continue
        if n in done: tot += done[n]; continue
        if n >= lim: break
        ts = time.time(); sv = dict(zip(sliced, bits))
        arrs, ixs = [], []
        for a, t in zip(arrays, inputs):
            arrs.append(a[tuple(sv[i] if i in sv else slice(None) for i in t)]); ixs.append(tuple(i for i in t if i not in sv))
        for x, y, keep in steps:
            C, cix = pair(arrs[x], ixs[x], arrs[y], ixs[y], keep)
            arrs[x] = arrs[y] = None; arrs.append(C); ixs.append(cix)
        v = complex(arrs[-1]); assert ixs[-1] == ()
        open(ck, "a").write(json.dumps(dict(i=n, re=v.real, im=v.imag, s=round(time.time() - ts, 2))) + chr(10))
        tot += v
        if n % 8 == 0: print("slice", n, "of", nsl, "%.1fs" % (time.time() - ts), "v=", v, flush=True)
    else:
        f = tot * sc / 2.0 ** nused
        print("RESULT f = %.16f %+.3ej  (sum of %d slices) wall %.0fs peak RSS %.0f MB" % (f.real, f.imag, nsl, time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), flush=True)
        return
    print("partial: did slices up to", lim, "wall %.0fs" % (time.time() - t0))
main()
