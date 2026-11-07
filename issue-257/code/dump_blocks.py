"""Print the raw circuit compressed into runs of gates on the same wire pair."""
import sys

from qasm_parse import parse_qasm

nq, gates = parse_qasm(sys.argv[1])
start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
stop = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
k = 0
while gates[k][0] == "x":
    k += 1
runs = []
cur = None
for idx, (name, par, qs) in enumerate(gates[k:], start=k):
    if len(qs) == 2:
        key = tuple(sorted(qs))
    else:
        key = None
    if cur is not None and (key == cur[0] or (key is None and qs[0] in cur[0])):
        cur[1].append((name, par, qs))
        continue
    if key is None:
        runs.append(((qs[0],), [(name, par, qs)], idx))
        cur = None
        continue
    cur = [key, [(name, par, qs)], idx]
    runs.append(cur)
for r in runs[start:stop]:
    key, gl, idx = r[0], r[1], r[2]
    desc = " ".join("%s%s%s" % (g[0], "(%.6g)" % g[1][0] if g[1] else "", "" if len(g[2]) == 1 and len(key) == 1 else str(list(g[2]))) for g in gl)
    print(idx, key, desc)
