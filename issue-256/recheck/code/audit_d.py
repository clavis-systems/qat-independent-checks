import json, math
for name, ck, res in (("A", "ckpt_L3_runA.slices.jsonl", "06_L3_runA_ctg_greedy_labels.jsonl"), ("B", "ckpt_L3_runB.slices.jsonl", "08_L3_runB_oe_randomgreedy.jsonl")):
    rows = [json.loads(l) for l in open("/s/logs/" + ck) if l.strip()]
    ns = sorted(r["n"] for r in rows)
    tot = complex(math.fsum(r["re"] for r in rows), math.fsum(r["im"] for r in rows))
    k = 44
    f = tot / 2 ** k
    last = [json.loads(l) for l in open("/s/logs/" + res) if l.startswith("{")][-1]
    print(name, "slices", len(rows), "unique", len(set(ns)), "range", ns[0], ns[-1], "sum/2^44 = %r %r" % (f.real, f.imag), "| logged", last["f_real"], last["f_imag"], "k", last["k"], "scalar", last["scalar"])
