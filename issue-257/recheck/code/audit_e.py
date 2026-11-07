import json,numpy as np
R='/v/results/'
def nf(inst,m,chi):
    S=json.load(open(R+f'{inst}_SCV_chi{chi}.json'))['rows']; M=json.load(open(R+f'{inst}_{m}_chi{chi}.json'))['rows']
    return np.array([b['N']-a['N'] for a,b in zip(S,M)]), S, M
a,_,_=nf('l12','meson',256); b,_,_=nf('l12','meson',2048)
print('l12 chi256 vs chi2048 steps1-11 max', np.abs(a-b)[:11].max(), 'steps', np.abs(a-b)[:11])
a,_,_=nf('l12','meson',1024); print('l12 chi1024 vs 2048 steps 1-15 max', np.abs(a-b)[:15].max(), ' s19,s20', (a-b)[18:])
a,_,_=nf('l300','MID',1024); b,S,M=nf('l300','MID',3072); print('l300 chi1024 vs 3072 steps<=16 max', np.abs(a-b)[:16].max(), 'steps17-20', (a-b)[16:])
a2,_,_=nf('l300','MID',2048); print('l300 2048 vs 3072 per step', np.abs(a2-b))
print('disc last', S[-1]['discarded'], M[-1]['discarded'], 'monotone', all(S[i]['discarded']<=S[i+1]['discarded'] for i in range(19)))
print('normdev max', max(r['normdev'] for r in S+M))
