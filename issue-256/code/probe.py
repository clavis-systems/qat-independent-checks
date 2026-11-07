"""Network sizes after each exact simplification, and quick greedy path costs, for L = 1, 2, 3."""
import time, math, sys
import structure, pipeline, tn

info = structure.analyse(verbose=False)
for L in (1, 2, 3):
    gates = info["all_gates"] if L == 3 else structure.build_L(info, L)
    for peep, lc in ((False, False), (False, True), (True, True)):
        t0 = time.time()
        arrays, inputs, size, k, st = pipeline.prepare(gates, structure.OBS, peep, lc)
        line = (f"L={L} peephole={peep} lightcone={lc}: gates {st['n_gates']} merged {st['n_merged']} "
                f"after_peephole {st['n_peephole']} after_lc {st['n_lightcone']} (cz {st['n_cz']}) "
                f"qubits {k} tensors {st['n_tensors']} indices {st['n_indices']}")
        if st["n_tensors"] > 1:
            ssa, tp = pipeline.find_path(inputs, size, "cotengra", seed=0, max_repeats=8, max_time=30)
            fl, mx = tn.cost(inputs, ssa)
            line += f" | greedy path: log2 flops {math.log2(fl):.1f}, log2 max {math.log2(mx):.1f} ({tp:.0f}s)"
        print(line, flush=True)
        if peep and lc:
            print("   qubits in final network:", st["qubits"], flush=True)
