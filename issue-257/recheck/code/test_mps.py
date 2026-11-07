import numpy as np, sys
from common import *
rng = np.random.default_rng(1)
# 1) gate-by-gate dense vs MPS on random NN circuits (all gate types)
for nq in (6, 8):
    gl = []
    for _ in range(300):
        t = rng.integers(0, 6)
        if t == 0: gl.append(('x', None, (int(rng.integers(nq)),)))
        elif t == 1: gl.append(('h', None, (int(rng.integers(nq)),)))
        elif t == 2: gl.append(('rz', float(rng.normal()), (int(rng.integers(nq)),)))
        else:
            a = int(rng.integers(nq-1)); q = (a, a+1) if rng.random()<.5 else (a+1, a)
            gl.append((['cx','swap','cx'][t-3], None, q))
    bits = rng.integers(0, 2, nq)
    psi = dense_run(nq, bits, gl)
    ops = fuse(gl)
    r = run(gl, ops, bits, 64, [len(ops)])
    print('random nq', nq, 'ops', len(ops), 'max|dZ|', np.abs(r[0][0]-dense_z(psi)).max(), 'disc', r[0][1])
# 2) generator vs real files
for f, th in (('SU2_T20_SCV.qasm', 0.226125), ('SU2_T20_MID.qasm', 0.226125)):
    g = parse_qasm('/data/'+f)
    init, steps = split_file(g)
    gen = gen_circuit(60, th, 20)
    flat = [x for s in steps for x in s]
    ok = all(a[0]==b[0] and a[2]==b[2] and ((a[1] is None and b[1] is None) or abs(a[1]-b[1])<1e-12) for a,b in zip(flat,gen)) and len(gen)-len(flat)==30 and all(x[0]=="swap" for x in gen[len(flat):])
    print(f, 'steps', len(steps), 'gen==file', ok, len(flat), len(gen), 'init x count', len(init))
# 3) small LSH circuit (N sites) dense vs MPS, per-step
for N in (4, 6):
    nq = 2*N
    gen = gen_circuit(N, 0.226125, 3)
    bits = init_bits(N); bits[2*N-3] ^= 1  # arbitrary extra flip
    lab = qubit_labels(N)
    per = len(gen)//3
    ops_steps = [fuse(gen[s*per:(s+1)*per]) for s in range(3)]
    ops = [o for s in ops_steps for o in s]; ends = np.cumsum([len(s) for s in ops_steps])
    r = run(gen, ops, bits, 256, list(ends))
    for s in range(3):
        psi = dense_run(nq, bits, gen[:(s+1)*per]); z = dense_z(psi); perm = perm_track(nq, gen[:(s+1)*per])
        print('N',N,'step',s+1,'|dZ|',np.abs(r[s][0]-z).max(),'Ndense',Nobs(z,lab,perm),'Nmps',Nobs(r[s][0],lab,perm),'Q',Qobs(z))
