import json,glob,re,numpy as np
R='/v/results/'
def rows(f): return json.load(open(f))['rows']
print(open(R+'l300_MID_chi3072.json').read()[:0])
d=json.load(open(R+'l300_MID_chi3072.json')); print(list(d.keys()), list(d['rows'][0].keys()), len(d['rows']))
def Nfrom(n):
    n=np.array(n); tot=0
    for r in range(len(n)):
        s=n[r].sum(); tot+= s if r%2==0 else 2-s
    return tot
for inst in ('l300','l12'):
    for chi in (256,512,1024,2048,3072):
        c={}
        for ck in ('SCV','MID' if inst=='l300' else 'meson'):
            f=R+f'{inst}_{ck}_chi{chi}.json'
            try: rr=rows(f)
            except Exception as e: break
            last=rr[-1]
            c[ck]=(last['N'],Nfrom(last['n']),last.get('disc'),last['step'],len(rr))
        if len(c)==2:
            a,b=list(c.values())
            print(inst,chi,'N_scv json %.6f recomputed %.6f | N_mid json %.6f recomp %.6f | nf %.6f | steps %d %d'%(a[0],a[1],b[0],b[1],b[0]-a[0],a[4],b[4]), 'disc',a[2],b[2])
# per-step nf chi3072 vs summary values
S=rows(R+'l300_SCV_chi3072.json'); M=rows(R+'l300_MID_chi3072.json')
print(' '.join('%.5f'%(m['N']-s['N']) for s,m in zip(S,M)))
print('Q dev max', max(abs(r['Q']-60) for r in S+M))
