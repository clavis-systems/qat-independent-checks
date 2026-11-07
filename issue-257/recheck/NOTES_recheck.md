# ricontrollo_257 - independent check of n_f(20), SU(2) LSH, l_i=300 (7/11/2026, written 10:55 local)

## Phase 1 (own code, done blind; logged 10:50:21 local in logs/PHASE1_RESULT.txt, snapshots logs/snap1_*.json)
Method: qubit-level MPS, 120 sites (one per qubit), SWAPs applied physically, exact gates fused only
into 4x4 blocks on the same qubit pair (code/qmps.py), SVD truncation chi (256; chi=128 for comparison),
no symmetry. Observables from Z expectations; N = 60 - 1/2 sum_r (-1)^r (Zi+Zo), site r from my own label tracking.
Validation: random NN circuits (cx,h,rz,x,swap) dense vs MPS 4e-15; reconstructed LSH circuit N=4,6 sites,
3 steps dense vs MPS <=4e-12; my generator == QASM gate-for-gate (both l300 files, all 20 steps, except the
last SWAP layer which the file omits); l=12 steps 1-6 n_f = 1.82534 1.34730 0.69035 0.02128 -0.49857 -0.75832.
Result l=300, chi=256 (disc = cumulative discarded weight): steps 1-9 disc<1e-7; 10: 1.3e-6; 11: 1.6e-5; 12: 1.3e-4.
n_f(t), t=1..12: 3.650671 2.715205 1.480150 0.286461 -0.589960 -1.007995 -0.973485 -0.610328 -0.107441
0.345667 0.618244 0.671375 (steps 13+ not trusted; run stopped at step 13, 10:52). Cost ~25 min per circuit to step 12
on 8 threads (chi=256, ~140 s/step at the end). chi=128 vs 256 agree to 3e-7 (n_f) up to step 9.
Cannot reach step 20: the non-symmetric qubit MPS needs chi~2000+ there.

## Phase 2: see final report (files read: verifica_257/NOTE.md, logs/summary_final.txt, results/*.json, issue_257.json)
