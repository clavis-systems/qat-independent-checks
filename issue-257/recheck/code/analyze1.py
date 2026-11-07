import json,sys,numpy as np
def L(t):
    try: return json.load(open(f'/w/logs/{t}.json'))
    except Exception: return []
S256=L('snap1_scv_chi256'); M256=L('snap1_mid_chi256'); S128=L('scv_chi128'); M128=L('mid_chi128')
print('step | Nscv256 Nmid256 nf256 | disc_s disc_m | nf128 (chi128) | |nf256-nf128|')
for t in range(1,21):
    if t>len(S256) or t>len(M256): break
    a=S256[t-1];b=M256[t-1]
    s='%2d | %.6f %.6f %.6f | %.1e %.1e |'%(t,a['N'],b['N'],b['N']-a['N'],a['disc_cum'],b['disc_cum'])
    if t<=len(S128) and t<=len(M128):
        nf128=M128[t-1]['N']-S128[t-1]['N']; s+=' %.6f | %.2e'%(nf128,abs(nf128-(b['N']-a['N'])))
    print(s)
