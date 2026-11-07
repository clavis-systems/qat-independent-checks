import re,sys
L=[l.strip().rstrip(';') for l in open(sys.argv[1]) if l.strip()][3:]
g=[]
for l in L:
    m=re.match(r'(\w+)(?:\(([^)]*)\))?\s+(.*)',l)
    n,a,q=m.groups(); q=[int(x) for x in re.findall(r'q\[(\d+)\]',q)]
    g.append((n,float(a) if a else None,q))
lab=[]
for m in range(30): lab+= [('i',2*m),('i',2*m+1),('o',2*m),('o',2*m+1)]
# segment by kind: consecutive gates; print label info at "events": hop blocks (cx..h) and elec blocks
i=0;ev=[];lastkind=None
def flush(k,items):
    if items: print(k,len(items),items[:3],'...',items[-2:])
cur=None;items=[]
n_steps=int(sys.argv[2])
cnt=0
while i<len(g) and cnt<n_steps:
    n,a,q=g[i]
    if n=='x' and i<60: i+=1; continue
    if n=='swap':
        kind='swap'; item=(q[0],q[1]); lab[q[0]],lab[q[1]]=lab[q[1]],lab[q[0]]; i+=1
    elif n=='cx' and g[i+1][0]=='h':
        kind='hop%.3f'%abs(g[i+2][1]); item=(lab[q[1]],lab[q[0]]) ; i+=8 if g[i+7][0]=='cx' else 7
    elif n=='x' or n=='rz':
        # electric block: x a; rz a; cx a b; rz b; cx a b; rz b; x a  (or rz only)
        if n=='x':
            kind='el%.4f'%abs(g[i+1][1]); item=(g[i][2][0],g[i+2][2][1],lab[g[i][2][0]],lab[g[i+2][2][1]]); i+=7
        else:
            kind='rz%.3f'%g[i][1]; item=(q[0],lab[q[0]]); i+=1
    else:
        kind='?'+n; item=q; i+=1
    if kind!=cur:
        flush(cur,items); items=[]; cur=kind; cnt+=1
    items.append(item)
flush(cur,items)
print('final lab first 12',lab[:12]); print(lab[56:64])
