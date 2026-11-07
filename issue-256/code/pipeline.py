"""Glue: gates -> ops -> network -> path (cotengra or opt_einsum) -> own slicing -> own contraction."""
import time
import math
import itertools
import numpy as np
import tn


def prepare(gates, obs, peep=True, lc=True, verbose=False):
    ops = tn.merge(gates)
    n_merge = len(ops)
    if peep:
        ops = tn.peephole(ops, verbose=verbose)
    n_peep = len(ops)
    if lc:
        ops = tn.lightcone(ops, obs)
    stats = dict(n_gates=len(gates), n_merged=n_merge, n_peephole=n_peep, n_lightcone=len(ops),
                 n_cz=sum(1 for o in ops if o[0] == "cz"),
                 qubits=sorted({q for o in ops for q in (o[1:] if o[0] == "cz" else (o[1],))}))
    arrays, inputs, size, k = tn.build_network(ops, obs)
    arrays, inputs, scalar = tn.absorb_vectors(arrays, inputs)
    size = {i: 2 for t in inputs for i in t}
    stats.update(k=k, n_tensors=len(arrays), n_indices=len(size), scalar=scalar)
    return arrays, inputs, size, k, stats


def find_path(inputs, size, method="cotengra", seed=0, max_repeats=32, max_time=60,
              ctg_methods=("greedy",), minimize="flops"):
    t0 = time.time()
    if method == "cotengra":
        import cotengra as ctg
        opt = ctg.HyperOptimizer(methods=list(ctg_methods), max_repeats=max_repeats,
                                 max_time=max_time, minimize=minimize, parallel=False,
                                 progbar=False, seed=seed)
        tree = opt.search(inputs, (), size)
        path = tree.get_path()
    elif method == "oe-random-greedy":
        import opt_einsum as oe
        syms = {}
        def sym(i):
            if i not in syms:
                syms[i] = oe.get_symbol(len(syms))
            return syms[i]
        eq = ",".join("".join(sym(i) for i in t) for t in inputs) + "->"
        shapes = [tuple(size[i] for i in t) for t in inputs]
        opt = oe.RandomGreedy(max_repeats=max_repeats, max_time=max_time, minimize="flops")
        path, _ = oe.contract_path(eq, *shapes, shapes=True, optimize=opt)
    else:
        raise ValueError(method)
    ssa = tn.ssa_from_linear(path, len(inputs))
    return ssa, time.time() - t0


def value(gates, obs, peep=True, lc=True, target_log2=24, method="cotengra", seed=0,
          max_repeats=32, max_time=60, ctg_methods=("greedy",), force_slices=0, log=print,
          ckpt=None, max_log2_flops=None, max_sliced=30):
    """ckpt: file prefix; if given, the path is saved to <ckpt>.path.json and every slice's
    partial value is appended to <ckpt>.slices.jsonl, so an interrupted run resumes."""
    import json, os, hashlib
    t0 = time.time()
    arrays, inputs, size, k, st = prepare(gates, obs, peep, lc)
    if len(arrays) == 0:
        st.update(total_time_s=round(time.time() - t0, 2), note="network reduced to a scalar")
        return st["scalar"] / 2 ** k, st
    net_hash = hashlib.sha256(repr(inputs).encode()).hexdigest()[:16]
    pfile = ckpt + ".path.json" if ckpt else None
    if pfile and os.path.exists(pfile):
        d = json.load(open(pfile))
        assert d["net_hash"] == net_hash, "network changed since checkpoint"
        ssa, tpath, sliced = [tuple(x) for x in d["ssa"]], d["path_time_s"], frozenset(d["sliced"])
        fl0, mx0 = tn.cost(inputs, ssa)
        fl, mx = tn.cost(inputs, ssa, sliced)
        st["resumed_path"] = True
    else:
        ssa, tpath = find_path(inputs, size, method, seed, max_repeats, max_time, ctg_methods)
        fl0, mx0 = tn.cost(inputs, ssa)
        sliced, fl, mx = tn.choose_slices(inputs, ssa, target_log2)
        if pfile:
            json.dump(dict(net_hash=net_hash, ssa=ssa, sliced=sorted(sliced), path_time_s=tpath,
                           method=method, seed=seed, ctg_methods=list(ctg_methods)),
                      open(pfile, "w"))
    st["net_hash"] = net_hash
    if ckpt:
        log(f"path ready: log2 flops unsliced {math.log2(fl0):.2f}, max {math.log2(mx0):.1f}; "
            f"{len(sliced)} sliced, per-slice log2 flops {math.log2(fl):.2f}, max {math.log2(mx):.1f}")
        st.update(path_method=method, path_seed=seed, path_time_s=round(tpath, 1),
                  log2_flops_unsliced=round(math.log2(fl0), 2),
                  log2_max_unsliced=round(math.log2(mx0), 2), n_sliced=len(sliced),
                  log2_flops_per_slice=round(math.log2(fl), 2), log2_max_sliced=round(math.log2(mx), 2),
                  log2_flops_total=round(math.log2(fl * 2 ** len(sliced)), 2))
        lim = float(os.environ.get("MAX_LOG2_FLOPS", "999"))
        if math.log2(fl * 2 ** len(sliced)) > lim:
            log(f"total cost above MAX_LOG2_FLOPS={lim}: stopping before contraction")
            st["stopped"] = "cost limit"
            return complex("nan"), st
        sfile = ckpt + ".slices.jsonl"
        done = {}
        if os.path.exists(sfile):
            for line in open(sfile):
                if line.strip():
                    r = json.loads(line)
                    done[r["n"]] = complex(r["re"], r["im"])
        sl = sorted(sliced)
        t1 = time.time()
        total = 0j
        nsl = 2 ** len(sl)
        for n, bits in enumerate(itertools.product((0, 1), repeat=len(sl))):
            if n in done:
                total += done[n]
                continue
            ts = time.time()
            v = tn.contract_slice(arrays, inputs, ssa, dict(zip(sl, bits)))
            with open(sfile, "a") as fh:
                fh.write(json.dumps(dict(n=n, re=v.real, im=v.imag, s=round(time.time() - ts, 2))) + "\n")
            total += v
            log(f"slice {n + 1}/{nsl} done in {time.time() - ts:.1f}s")
        st["contract_time_s_this_session"] = round(time.time() - t1, 2)
        st["slices_resumed"] = len(done)
        st["total_time_s"] = round(time.time() - t0, 2)
        return st["scalar"] * total / 2 ** k, st
    if target_log2 < 0:          # relative target: unsliced max intermediate minus |target|
        target_log2 = math.log2(mx0) + target_log2
    sliced, fl, mx = tn.choose_slices(inputs, ssa, target_log2, max_slices=max_sliced)
    if force_slices and len(sliced) < force_slices:
        # extra slices on the most frequent indices, only to exercise the slicing code
        extra = [i for i in sorted(size) if i not in sliced][:force_slices - len(sliced)]
        sliced = frozenset(sliced) | set(extra)
        fl, mx = tn.cost(inputs, ssa, frozenset(sliced))
    st.update(path_method=method, path_seed=seed, path_time_s=round(tpath, 1),
              log2_flops_unsliced=round(math.log2(fl0), 2), log2_max_unsliced=round(math.log2(mx0), 2),
              n_sliced=len(sliced), log2_flops_per_slice=round(math.log2(fl), 2),
              log2_max_sliced=round(math.log2(mx), 2),
              log2_flops_total=round(math.log2(fl * 2 ** len(sliced)), 2))
    if max_log2_flops is not None and st["log2_flops_total"] > max_log2_flops:
        st["stopped"] = "cost limit"
        return complex("nan"), st
    t1 = time.time()
    Z = st["scalar"] * tn.contract(arrays, inputs, ssa, sliced, log=None)
    st["contract_time_s"] = round(time.time() - t1, 2)
    st["total_time_s"] = round(time.time() - t0, 2)
    f = Z / 2 ** k
    return f, st
