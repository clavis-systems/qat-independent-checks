"""Main exact runs. Usage: python run_main.py L peep lc method seed target_log2 max_time [ctg_methods]

Writes one JSON line with the value, the run statistics, versions, wall time and peak memory.
"""
import sys, time, json, math, resource, platform, hashlib
import numpy as np
import structure, pipeline


def main():
    L = int(sys.argv[1]); peep = sys.argv[2] == "1"; lc = sys.argv[3] == "1"
    method = sys.argv[4]; seed = int(sys.argv[5]); target = int(sys.argv[6])
    max_time = float(sys.argv[7])
    ctg_methods = tuple(sys.argv[8].split(",")) if len(sys.argv) > 8 else ("greedy",)
    ckpt = sys.argv[9] if len(sys.argv) > 9 else None      # e.g. /w/logs/ckpt_L3_runA
    import cotengra, opt_einsum
    t0 = time.time()
    info = structure.analyse(verbose=False)
    gates = info["all_gates"] if L == 3 else structure.build_L(info, L)
    f, st = pipeline.value(gates, structure.OBS, peep=peep, lc=lc, target_log2=target,
                           method=method, seed=seed, max_repeats=10 ** 6, max_time=max_time,
                           ctg_methods=ctg_methods, ckpt=ckpt,
                           log=lambda m: print("#", time.strftime("%H:%M:%S"), m, flush=True))
    rec = dict(L=L, peephole=peep, lightcone=lc, ctg_methods=ctg_methods,
               f_real=repr(f.real), f_imag=repr(f.imag),
               diff_from_claim=(f.real - 0.822488318196) if L == 3 else None,
               wall_s=round(time.time() - t0, 1),
               peak_rss_MB=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
               python=platform.python_version(), numpy=np.__version__,
               cotengra=cotengra.__version__, opt_einsum=opt_einsum.__version__,
               circuit_sha256=hashlib.sha256(open(structure.CIRC, "rb").read()).hexdigest())
    st.pop("qubits", None)
    rec.update({k: (v if not isinstance(v, complex) else repr(v)) for k, v in st.items()})
    print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
