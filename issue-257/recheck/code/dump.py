import re,sys
L=[l.strip().rstrip(';') for l in open(sys.argv[1]) if l.strip()][3:]
g=[]
for l in L:
    m=re.match(r'(\w+)(?:\(([^)]*)\))?\s+(.*)',l)
    n,a,q=m.groups(); q=[int(x) for x in re.findall(r'q\[(\d+)\]',q)]
    g.append((n,float(a) if a else None,q))
# tokenise: collapse hopping block pattern cx(b,a) h b rz cx(a,b) rz cx(a,b) h cx(b,a)
out=[];i=0
while i<len(g):
    if g[i][0]=='cx' and i+7<=len(g) and [x[0] for x in g[i:i+7]]==['cx','h','rz','cx','rz','cx','h'] :
        out.append(('HOP',g[i][2][1],g[i][2][0],g[i+2][1],g[i+4][1],g[i+7][0] if i+7<len(g) else None)); i+=8 if g[i+7][0]=='cx' else 7
    else:
        out.append(g[i]);i+=1
for o in out[:int(sys.argv[2])]: print(o)
