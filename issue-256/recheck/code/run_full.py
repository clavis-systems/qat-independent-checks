"""Full L=3 run: cone prune + own peephole + quimb simplify + cotengra path (+slicing) + contraction.
usage: run_full.py <optimizer: greedy|labels|oe> <seed> <search_seconds> <log2 slice target>"""
import sys, time, random, json, resource, numpy as np, cotengra as ctg
from ole import *
from tnbuild import *
from peep import *

opt_name = sys.argv[1]; seed = int(sys.argv[2]); budget = float(sys.argv[3]); target = int(sys.argv[4])
random.seed(seed); np.random.seed(seed)
t0 = time.time()
A, M, B = load("/data/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm")
V = extract_V(M, 49)
Pq = sorted(x[2][0] for x in M if x[0] == "rz")
Oq = [52, 59, 72]
G, wires = prepare(A, V, Oq, Pq, prune=True)
ops, sc, wires, nA, nG = build_v2(A, V, Oq, Pq)
print("peepholed U ops", nA, "cone ops", nG)
used = sorted({w for o in ops for w in wires_of(o)})
print("after prune: wires %d; peephole ops %d (cz %d), used wires %d, scalar %s" % (len(wires), len(ops), sum(o[0] == "cz" for o in ops), len(used), sc), flush=True)
tn = tn_from_ops(ops, used)
tn.full_simplify_(seq="ADCRS", atol=1e-12)
print("simplified TN tensors %d inds %d (%.1fs)" % (tn.num_tensors, tn.num_indices, time.time() - t0), flush=True)
inputs, output, size_dict = tn.get_inputs_output_size_dict(output_inds=())
inputs = [tuple(i) for i in inputs]
opt = ctg.HyperOptimizer(methods=opt_name.split(","), max_time=budget, max_repeats=100000, parallel=False, optlib="random",
                         minimize="flops", progbar=False)
import pickle, os
pf = "/w/logs/tree_%s_%d.pkl" % (sys.argv[1].replace(",", "+"), seed)
if os.path.exists(pf): tree = pickle.load(open(pf, "rb")); print("loaded saved tree", flush=True)
else: tree = opt.search(inputs, output, size_dict)
print("unsliced: log2 width %.1f log10 flops %.2f" % (tree.contraction_width(), np.log10(float(tree.contraction_cost()))), flush=True)
if target and not os.path.exists(pf): tree.slice_(target_size=2 ** target); pickle.dump(tree, open(pf, "wb"))
print("path: optimizer %s seed %d search %.0fs; log2 width %.1f; log10 flops(total) %.2f; nslices %d" % (
    sys.argv[1], seed, time.time() - t0, tree.contraction_width(), np.log10(float(tree.contraction_cost())), tree.nslices), flush=True)
arrays = [t.data for t in tn.tensors]
t1 = time.time()
import os
ck = "/w/logs/ckpt_%s_%d.jsonl" % (sys.argv[1].replace(",", "+"), seed)
done = {}
if os.path.exists(ck):
    for l in open(ck):
        r = json.loads(l); done[r["i"]] = complex(r["re"], r["im"])
tot = 0j
for i in range(tree.multiplicity):
    if i in done: tot += done[i]; continue
    ts = time.time(); v = complex(tree.contract_slice(arrays, i, backend="numpy"))
    open(ck, "a").write(json.dumps(dict(i=i, re=v.real, im=v.imag, s=round(time.time() - ts, 2))) + chr(10))
    tot += v
    if i % 16 == 0: print("slice", i, "of", tree.multiplicity, "%.1fs" % (time.time() - ts), flush=True)
val = tot
val = complex(val)
f = val * sc * 10.0 ** tn.exponent / 2.0 ** len(used) if hasattr(tn, "exponent") else val * sc / 2.0 ** len(used)
print("RESULT optimizer=%s seed=%d f = %.16f %+.3ej  contraction %.0fs wall total %.0fs peak RSS %.0f MB" % (
    sys.argv[1], seed, f.real, f.imag, time.time() - t1, time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), flush=True)
