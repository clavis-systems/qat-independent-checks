import numpy as np, sys
from common import *
lab=qubit_labels(60)
for f in ['/data/SU2_T20_SCV.qasm','/data/SU2_T20_MID.qasm','/t/x_100_SCV.qasm','/t/x_100_meson.qasm']:
    g=parse_qasm(f); init,steps=split_file(g); b=bits_from_init(120,init)
    occ={}
    for q in range(120):
        t,r=lab[q]; occ.setdefault(r,{})[t]=int(b[q])
    ref={r:((1,1) if r%2 else (0,0)) for r in range(60)}
    diff=[(r,ref[r],(occ[r]['i'],occ[r]['o'])) for r in range(60) if (occ[r]['i'],occ[r]['o'])!=ref[r]]
    print(f.split('/')[-1],'x-count in prefix',len(init),'sites differing from SCV (site, SCV(i,o), this(i,o)):',diff,'Q=',int(b.sum()))
    pg=[x for s in steps for x in s]; perm=perm_track(120,pg)
    ok=all(lab[perm[p]][1]==p//2 for p in range(120)); print('  final layout: position p holds site p//2 ->',ok, ' i/o order ok:', all(lab[perm[2*r]][0]=='i' and lab[perm[2*r+1]][0]=='o' for r in range(60)))
