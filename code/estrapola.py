"""Estrapolazione dell'energia DMRG rispetto al peso scartato, su diversi sottoinsiemi di punti (scritto per C2).
Uso: python -I estrapola.py dmrg.json"""
import json
import sys

import numpy as np

r = json.load(open(sys.argv[1]))
E239 = -116.60560912042631
pts = sorted(r["sweep_giu"], key=lambda p: -p["M"])
out = {}
for k in (2, 3, 4, 5):
    x = np.array([p["dw"] for p in pts[:k]])
    y = np.array([p["energia"] for p in pts[:k]])
    a, b = np.polyfit(x, y, 1)
    out[f"lineare_{k}_punti_M>={pts[k-1]['M']}"] = dict(E0=float(b), diff_da_239_uEh=float((b - E239) * 1e6))
x = np.array([p["dw"] for p in pts[:4]])
y = np.array([p["energia"] for p in pts[:4]])
c = np.polyfit(x, y, 2)
out["quadratica_4_punti"] = dict(E0=float(c[-1]), diff_da_239_uEh=float((c[-1] - E239) * 1e6))
out["variazionale_Mmax"] = dict(M=r["Mmax"], E=r["energia_Mmax"], diff_da_239_uEh=(r["energia_Mmax"] - E239) * 1e6)
for k, v in out.items():
    print(k, v)
json.dump(out, open(sys.argv[1].replace(".json", "_estrapolazioni.json"), "w"), indent=1)
