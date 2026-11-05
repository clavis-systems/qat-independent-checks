"""V3: estrapolazione E(peso scartato) -> 0 dagli sweep a M decrescente (codice di V3)."""
import sys, json, numpy as np
r = json.load(open(sys.argv[1]))
M = np.array(r["bond_dims"]); E = np.array(r["energies"]); W = np.array(r["discarded"])
E_issue = -116.60560912042631
print("ultimo sweep a M=750:", E[13], " differenza dalla issue (µEh):", (E[13] - E_issue) * 1e6)
rev = list(range(14, len(M)))  # 650, 550, 450, 350
print("sweep a M decrescente:")
for i in rev: print(f"  M={M[i]}  E={E[i]:.10f}  peso scartato={W[i]:.3e}")
out = {}
for k in (2, 3, 4):
    idx = rev[:k]
    p = np.polyfit(W[idx], E[idx], 1)
    out[f"lineare_{k}_punti"] = p[1]
    print(f"lineare, {k} punti (M>={M[idx[-1]]}): E0={p[1]:.7f}  issue+{(p[1]-E_issue)*1e6:.1f} µEh")
p = np.polyfit(W[rev], E[rev], 2)
out["quadratica_4_punti"] = p[2]
print(f"quadratica, 4 punti: E0={p[2]:.7f}  issue+{(p[2]-E_issue)*1e6:.1f} µEh")
json.dump(out, open(sys.argv[2], "w"), indent=1)
