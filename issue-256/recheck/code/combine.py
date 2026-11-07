import json, sys, glob
paths = sys.argv[2:]; nused = 44; sc = complex(sys.argv[1])
rows = {}
for p in paths:
    for l in open(p):
        r = json.loads(l); rows[r["i"]] = complex(r["re"], r["im"])
import math
print("slices found", len(rows), "unique idx", len(rows), "min/max", min(rows), max(rows))
tot = complex(math.fsum(v.real for v in rows.values()), math.fsum(v.imag for v in rows.values()))
f = tot * sc / 2.0 ** nused
print("f = %.16f %+.3ej" % (f.real, f.imag))
