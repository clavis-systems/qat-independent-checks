"""Summarize own MPS runs: n_f(t) = N_MID - N_SCV per chi, plus comparison columns
with numbers quoted in the tracker issues (#255 for l_i=12, #257 for l_i=300)."""
import glob
import json
import os
import sys

RES = sys.argv[1] if len(sys.argv) > 1 else "/res"

# quoted from issue texts (for comparison only; not used in any computation)
REF255_NF = [1.825336, 1.347301, 0.690350, 0.021282, -0.498568, -0.758322, -0.726334, -0.454570, -0.058119,
             0.323064, 0.567767, 0.608895, 0.450079, 0.159266, -0.157519, -0.392936, -0.473363, -0.381026,
             -0.156713, 0.116863]
REF255_STAG = [-54.845406, -40.722478, -21.265902, -1.358549, 14.254973, 22.278867, 21.713781, 14.002762,
               2.455765, -8.846881, -16.318932, -17.908604, -13.607112, -5.295247, 3.989935, 11.107720,
               13.835759, 11.537197, 5.276803, -2.623850]
REF257_SCV = [5.15459, 18.90818, 36.93115, 54.16146, 66.66919, 72.65165, 72.42626, 67.79642, 61.21971, 55.02823,
              50.82861, 49.20024, 49.77799, 51.65888, 53.90227, 55.86418, 57.25892, 58.03763, 58.24810, 57.96720]
REF257_MID = [8.80527, 21.62338, 38.41130, 54.44792, 66.07923, 71.64365, 71.45278, 67.18609, 61.11227, 55.37390,
              51.44685, 49.87159, 50.32907, 52.00758, 54.05500, 55.88232, 57.22023, 58.00948, 58.27527, 58.07005]


def load(lab, chi):
    p = os.path.join(RES, "%s_chi%d.json" % (lab, chi))
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def chis(lab):
    out = []
    for p in glob.glob(os.path.join(RES, lab + "_chi*.json")):
        out.append(int(p.rsplit("_chi", 1)[1].split(".")[0]))
    return sorted(out)


for inst, a, b in (("l12", "l12_SCV", "l12_meson"), ("l300", "l300_SCV", "l300_MID")):
    print("=" * 100)
    print("instance", inst)
    for chi in chis(a):
        A, B = load(a, chi), load(b, chi)
        if A is None or B is None:
            continue
        ra, rb = A["rows"], B["rows"]
        n = min(len(ra), len(rb))
        last = n
        sa, sb = ra[n - 1], rb[n - 1]
        print("chi=%5d steps=%2d | N_SCV=%.6f N_MID=%.6f n_f=%.6f | Q_SCV-60=%.1e Q_MID-60=%.1e | disc SCV=%.2e MID=%.2e"
              " | maxbond %d/%d | Smax %.3f/%.3f | time %.0f/%.0f s | peakRSS %.0f/%.0f MB"
              % (chi, last, sa["N"], sb["N"], sb["N"] - sa["N"], sa["Q"] - 60, sb["Q"] - 60, sa["discarded"],
                 sb["discarded"], sa["max_bond"], sb["max_bond"], sa["max_entropy"], sb["max_entropy"],
                 sa["elapsed"], sb["elapsed"], sa["peak_rss_mb"], sb["peak_rss_mb"]))
    # convergence: successive differences and a linear fit of n_f vs mean discarded weight (last 3 chi)
    pts = []
    for chi in chis(a):
        A, B = load(a, chi), load(b, chi)
        if A and B and len(A["rows"]) == 20 and len(B["rows"]) == 20:
            ra, rb = A["rows"][-1], B["rows"][-1]
            pts.append((chi, rb["N"] - ra["N"], 0.5 * (ra["discarded"] + rb["discarded"]), ra["N"], rb["N"]))
    for i in range(1, len(pts)):
        print("  chi %d -> %d: d n_f = %+.6f, d N_SCV = %+.6f, d N_MID = %+.6f"
              % (pts[i - 1][0], pts[i][0], pts[i][1] - pts[i - 1][1], pts[i][3] - pts[i - 1][3], pts[i][4] - pts[i - 1][4]))
    if len(pts) >= 3:
        import numpy as np
        x = np.array([p[2] for p in pts[-3:]])
        for j, name in ((1, "n_f"), (3, "N_SCV"), (4, "N_MID")):
            y = np.array([p[j] for p in pts[-3:]])
            c = np.polyfit(x, y, 1)
            print("  linear fit of %s vs mean discarded weight (chi=%s): intercept %.6f, slope %.4f"
                  % (name, [p[0] for p in pts[-3:]], c[1], c[0]))
    # per-step table for the largest complete chi
    full = [c for c in chis(a) if load(a, c) and load(b, c) and len(load(a, c)["rows"]) == 20 and len(load(b, c)["rows"]) == 20]
    if not full:
        continue
    c = full[-1]
    A, B = load(a, c)["rows"], load(b, c)["rows"]
    print("per-step at chi=%d (own) vs quoted refs" % c)
    for k in range(20):
        nf = B[k]["N"] - A[k]["N"]
        if inst == "l12":
            print("  step %2d n_f=%.6f  [#255 PT2: %.6f, diff %.1e]  stag_SCV=%.6f [#255: %.6f, diff %.1e]  disc %.1e/%.1e"
                  % (k + 1, nf, REF255_NF[k], nf - REF255_NF[k], A[k]["N"] - 60, REF255_STAG[k],
                     A[k]["N"] - 60 - REF255_STAG[k], A[k]["discarded"], B[k]["discarded"]))
        else:
            print("  step %2d N_SCV=%.5f [#257 %.5f d=%.1e] N_MID=%.5f [#257 %.5f d=%.1e] n_f=%.5f [#257 %.5f]  disc %.1e/%.1e"
                  % (k + 1, A[k]["N"], REF257_SCV[k], A[k]["N"] - REF257_SCV[k], B[k]["N"], REF257_MID[k],
                     B[k]["N"] - REF257_MID[k], nf, REF257_MID[k] - REF257_SCV[k], A[k]["discarded"], B[k]["discarded"]))
