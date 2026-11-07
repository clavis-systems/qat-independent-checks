import sys, json, time, numpy as np
from common import *
f, chi, tag = sys.argv[1], int(sys.argv[2]), sys.argv[3]
nsteps = int(sys.argv[4]) if len(sys.argv) > 4 else 20
g = parse_qasm(f)
init, steps = split_file(g)
nq = 120
bits = bits_from_init(nq, init)
lab = qubit_labels(60)
ops_steps = [fuse(s) for s in steps]
def log(s): print(time.strftime('%H:%M:%S'), s, flush=True)
print('occupied qubits', [i for i in range(nq) if bits[i]], flush=True)
log(f'{f} chi={chi} ops/step={[len(s) for s in ops_steps][:3]}')
m = MPS(bits, chi, 1e-14)
out = []; pg = []; t0 = time.time()
for s in range(nsteps):
    for op in ops_steps[s]:
        if op[0] == 1: m.apply1(op[1], op[2])
        else: m.apply2(op[1], op[2])
    pg += steps[s]
    z, _ = m.zexp()
    perm = perm_track(nq, pg)
    out.append(dict(step=s+1, N=float(Nobs(z, lab, perm)), Q=float(Qobs(z)), disc_cum=float(m.disc), maxbond=int(m.maxbond), z=z.tolist(), perm=[int(x) for x in perm]))
    json.dump(out, open(f'/w/logs/{tag}.json', 'w'))
    log('step %d N=%.6f Q=%.6f disc_cum=%.3e maxbond=%d t=%.0fs' % (s+1, out[-1]['N'], out[-1]['Q'], m.disc, m.maxbond, time.time()-t0))
