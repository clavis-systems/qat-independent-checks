"""Build the full TN, simplify, find contraction trees with two optimizers, report cost/width. No contraction."""
import sys, time, json, pickle, numpy as np, cotengra as ctg
from ole import *
from tnbuild import *

import random; opt_name = sys.argv[1]; seed = int(sys.argv[2]); random.seed(seed); np.random.seed(seed); budget = float(sys.argv[3]); target = int(sys.argv[4]) if len(sys.argv) > 4 else None
A, M, B = load("/data/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm")
V = extract_V(M, 49)
Pq = sorted(x[2][0] for x in M if x[0] == "rz")
Oq = [52, 59, 72]
G, wires = prepare(A, V, Oq, Pq, prune=True)
print("G", len(G), "cz", sum(x[0] == "cz" for x in G), "wires", len(wires), flush=True)
t0 = time.time()
tn = make_tn(G, {t: m for t, m in V.items() if t in wires}, Oq, wires)
print("raw TN: tensors", tn.num_tensors, "inds", tn.num_indices, flush=True)
tn.full_simplify_(seq="ADCRS", atol=1e-12)
print("simplified TN: tensors", tn.num_tensors, "inds", tn.num_indices, "t=%.1fs" % (time.time() - t0), flush=True)
pickle.dump(tn, open("/w/tn_simplified.pkl", "wb"))
inputs, output, size_dict = tn.get_inputs_output_size_dict(output_inds=())
inputs = [tuple(i) for i in inputs]
if opt_name == "greedy":
    opt = ctg.HyperOptimizer(methods=["greedy"], max_time=budget, max_repeats=100000, parallel=False, optlib="random",
                             minimize="combo-64", slicing_opts=dict(target_size=2 ** target) if target else None, progbar=False)
elif opt_name == "kahypar":
    opt = ctg.HyperOptimizer(methods=["kahypar"], max_time=budget, parallel=False, optlib="random", minimize="combo-64",
                             slicing_opts=dict(target_size=2 ** target) if target else None)
else:
    opt = ctg.HyperOptimizer(methods=[opt_name], max_time=budget, max_repeats=100000, parallel=False, optlib="random",
                             minimize="combo-64", slicing_opts=dict(target_size=2 ** target) if target else None, progbar=False)
tree = opt.search(inputs, output, size_dict)
print("opt", opt_name, "seed", seed, "log2 width %.1f" % tree.contraction_width(), "log10 flops %.2f" % np.log10(float(tree.contraction_cost())),
      "nslices", getattr(tree, "nslices", 0), "t=%.1fs" % (time.time() - t0), flush=True)
if target:
    print("sliced inds", len(tree.sliced_inds), "slices", tree.nslices, "log2 slice width", tree.contraction_width(), flush=True)
pickle.dump(tree, open("/w/tree_%s_%d.pkl" % (opt_name, seed), "wb"))
