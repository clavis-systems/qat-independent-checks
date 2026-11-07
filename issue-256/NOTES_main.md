# Verifica issue #256 (operator_loschmidt_echo_49x648) — running notes

Goal: independent exact value of f for `49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm`,
claimed f = 0.822488318196 (issue #256, dylanneve1). Own parser, own network builder,
own contraction/slicing; cotengra used only to find the contraction path.
The submitter's code (qsim-lab) is NOT read before we have our own number.

Everything runs in Docker container `v256` (python:3.13-slim, --memory 3g --cpus 4),
folder mounted at /w, circuits mounted read-only at /circ.

## Log
- 16:34 start. Issue text read. Circuit SHA-256
  4aae5089337253921cd70df470dc52b5243e0c8410b36a6bcbb43497db48db29.
  Tracker README (operator_loschmidt_echo/README.md): f = 2^-n Tr(U O U^† V^† U O U^† V);
  sum form 2^-n sum_z sigma_z <z|...|z>. Circuit file = C = U2 . V . U1 (time order U1, V, U2),
  V = exp(-i 0.15 X) on 24 P sites (sdg,sxdg,rz(0.3),sx,s). Definition used here:
  f = 2^-n Tr(O C^† O C), O = Z52 Z59 Z72 (= expectation of O after C, averaged with sign sigma_z
  over all basis inputs z). Note: README sum formula as printed omits the last U; the trace
  form and the circuit structure agree with Tr(O C^† O C).
- 16:45 container v256 up: numpy 2.5.4, cotengra 0.8.3, opt_einsum 3.4.0, qiskit 2.6.2.

- 16:55 code/qasm_parse.py (own parser): 4756 gates, 49 used qubits of 156 (logs/01_parse.log).
- 17:00 code/structure.py (logs/02_structure.log): 54 edges, 54 CZ per Trotter step;
  S (rx +-0.5) = 13 sites [7,8,9,10,17,30,39,50,57,67,68,69,70]; P = 24 sites, P & S empty.
  Checked gate-for-gate: second half == inverse(first half), i.e. C = U^dag V U;
  first half U = A.B with B == inverse(A) with rx sign flipped on S (submitter's U2 = Z_S U1 Z_S);
  A = 3 identical Trotter steps of 378 gates. Rebuilt L=3 circuit == file, so L=1, L=2 circuits
  built the same way (code/structure.py build_L).
- 17:10 code/tn.py, pipeline.py: exact simplifications (merge 1q, cancel CZ pairs separated only
  by diagonal ops, light cone of O in both directions), doubled hyperedge network, cotengra path,
  own slicing + own batch-matmul executor. Probe (logs/03_probe.log): L=1 network vanishes
  entirely after the peephole (f(L=1) = 1 exactly); L=2: 18 qubits, 220 tensors, ~2^13 flops.
  L=3 probe running.
- 17:12 validate_small.py running (logs/04_validate_small.log).

  (Clock note: times above were estimates; `date` at 16:50 says the whole setup took ~16 min.)
- 16:48 L=3 probe: after peephole+light cone 554 ops (214 CZ), 44 qubits, 1108 tensors;
  a quick greedy path gives 2^36.6 flops, max intermediate 2^31 (slicing needed).
- 16:49 f(L=1) = 1.0 exactly (network collapses to a scalar) (logs/05_L1_peephole_lc.jsonl);
  f(L=2) = 0.8972543929065162, Im -4e-17 (logs/05_L2_peephole_lc.jsonl). Submitter: 1 and
  0.897254392906.
- 16:50 qiskit cross-parse: identical gate list (4756 gates). Validation on patches running.
- 16:50 L=3 run A started (cotengra greedy+labels, 300 s search, own slicing to 2^23).

- 16:50-18:28 PAUSED by the user (container stopped, run A and validation lost, nothing finished).
- 18:31 resumed. Added checkpointing: pipeline.value(ckpt=prefix) saves the path to
  logs/<prefix>.path.json and each slice's value to logs/<prefix>.slices.jsonl; rerunning the same
  command resumes. Restarted: run A (logs/06_L3_runA_ctg_greedy_labels.jsonl, ckpt_L3_runA.*) and
  validation (logs/04_validate_small.log).
  Restart commands (from Git Bash, MSYS_NO_PATHCONV=1, after `docker start v256`):
  docker exec -e OPENBLAS_NUM_THREADS=2 -w /w/code v256 python -u run_main.py 3 1 1 cotengra 0 23 300 greedy,labels /w/logs/ckpt_L3_runA
  docker exec -e OPENBLAS_NUM_THREADS=2 -w /w/code v256 python -u validate_small.py

- 18:33 container v256 stopped again (exit 255, not OOM; at the same time Docker destroyed the
  other session's v257_* containers, likely a Docker restart). Nothing had been saved yet.
  18:35 restarted; run A (path search 180 s) and validation relaunched detached (docker exec -d),
  output written straight to logs/ inside the mount.

- 18:38 RESULT run A (logs/06_L3_runA_ctg_greedy_labels.jsonl, path logs/ckpt_L3_runA.path.json,
  slices logs/ckpt_L3_runA.slices.jsonl): f(L=3) = 0.8224883181962652, Im = +4.9e-16.
  Difference from claim 0.822488318196: +2.65e-13 (claim given to 12 decimals: agrees).
  Network: 44 qubits, 1108 tensors; path cotengra greedy+labels (205 s search): 2^33.0 flops,
  max 2^26 unsliced; own slicing 3 indices -> 8 slices, 2^30.5 flops and 2^23 max per slice;
  contraction 13.8 s with own executor; peak RSS 486 MB; wall 220 s.
- 18:39 launched run B (opt_einsum random-greedy path, seed 7, slicing target 2^21,
  logs/08_L3_runB_oe_randomgreedy.jsonl) and run C (L=2 without peephole,
  logs/07_L2_no_peephole_lc_only.jsonl).

- 18:46 RESULT run B (logs/08_L3_runB_oe_randomgreedy.jsonl, ckpt_L3_runB.*): path from
  opt_einsum RandomGreedy (seed 7, 181 s; no cotengra), 2^36.7 flops unsliced, own slicing of
  9 indices -> 512 slices, 2^29.9 flops / 2^21 max per slice, 2^38.9 total; contraction 311 s;
  peak RSS 177 MB. f = 0.8224883181962643, Im +1.3e-16. |A - B| = 9e-16.
- 18:44 run C (L=2 without peephole, cotengra greedy+labels) killed: path 2^41 flops, max 2^36
  -> 2^52 after slicing, infeasible. Relaunched 18:47 with path from cotengra+KaHyPar (kahypar 1.3.7
  pip wheel; only for this cross-check of the peephole), cost cap 2^38
  (logs/07_L2_no_peephole_lc_only.jsonl, ckpt_L2_nopeep_kahypar.*).
- Run D (L=3 without peephole, cotengra greedy+labels, cost cap 2^39.5) still searching.

- 18:47 run D (L=3, light cone only, no peephole): best cotengra greedy+labels path 2^86 flops
  -> infeasible, stopped by the cost cap (logs/09_L3_runD_no_peephole_lc_only.jsonl).
- 18:47 validation with k=12 patches too slow under the CPU limit (partial log kept as
  logs/04a_validate_small_k12_aborted.log); relaunched with k=10/11 patches.
- 18:51 L=2 without peephole, KaHyPar path: 2^36.5 unsliced, 512 slices, 2^38.0 total; path saved,
  contraction running on the saved path (logs/07b_L2_no_peephole_lc_only_contract.jsonl).
- 18:52 run E (L=3 peephole, no light cone, greedy+labels): 2^71 -> stopped
  (logs/10_L3_runE_peephole_no_lightcone.jsonl). 18:55 retry with KaHyPar path, cap 2^38.5
  (logs/10b_L3_runE_peephole_no_lightcone_kahypar.jsonl).

- 18:59 validation stuck >10 min in the 2nd patch (probably the raw network with a bad greedy
  path). Killed; validate_small.py now prints each method as it finishes, and the three
  'reduced simplification' variants (light cone only, peephole only, raw) use a KaHyPar path with
  a cost guard 2^33 (skipped if above). Relaunched (logs/04_validate_small.log).

- 19:04 RESULT L=2 WITHOUT peephole (light cone only; 38 qubits, 1072 tensors; KaHyPar path,
  own slicing 512 slices, 2^38.0 flops; logs/07b_L2_no_peephole_lc_only_contract.jsonl):
  f(L=2) = 0.8972543929064776 vs 0.8972543929065162 with peephole: |diff| = 3.9e-14.
  -> full-scale check that the peephole cancellation does not change the value (at L=2).
- 19:00 run E with KaHyPar (L=3 peephole, no light cone): 2^56 -> stopped, not feasible
  (logs/10b_L3_runE_peephole_no_lightcone_kahypar.jsonl). The light-cone pass is checked only
  on the small patches (variant tn_peep_only) and indirectly by L=1/L=2 agreement.

- 19:06 the real cause of slow validation: the 'tn_sliced' variant used an absolute target of
  2^4 elements -> tens of thousands of slices. Now relative (unsliced max / 2^4, at least 3
  sliced indices). Relaunched (logs/04_validate_small.log; partial logs 04a/04b kept).

- 19:12 threshold robustness (code/tolerance_check.py, logs/11_tolerance_check.log): over 21573
  diagonal tests in the L=3 simplification, largest off-diagonal treated as zero 2.4e-16,
  smallest treated as non-zero 0.149 -> classification is not near the 1e-14 threshold.
  (The identity-drop test is never triggered for this circuit: identity products are carried
  as diagonal and absorbed into the next gate.) Versions in logs/00_versions.txt.

- 19:22 validation cases 1-4 complete (logs/04c_validate_small_cases1-4.log; max |dev| per case
  2.8e-14, 1.5e-14, 1.3e-14, 1.6e-14). Case 5 stuck again in tn_sliced (relative target on a
  tiny network still gave too many slices). Added max_sliced cap (5 indices) for that variant and
  ran cases 5-6 separately (logs/04d_validate_small_cases5-6.log).
  Code note: pipeline.py was edited during the session (checkpointing, cost guard, relative
  slicing target, max_sliced); none of these touch the code path of runs A/B (ckpt branch,
  absolute targets 2^23/2^21). To recompute A/B from scratch delete logs/ckpt_L3_run*.slices.jsonl.

- 19:28 VALIDATION cases 1-6 done (04c: cases 1-4, 04d: cases 5-6). Patches of 10 qubits, all
  gates of the full L=3 circuit inside the patch (942-1016 gates). Values: 1 (trivial patch),
  0.944575, 0.983433, 0.964973, 0.999580, 0.967669. For each: own dense unitary (raw
  gates) vs qiskit Operator (qiskit.qasm3 parse) vs six TN variants (full simplification;
  forced slicing; opt_einsum path; light cone only; peephole only; raw). Max |deviation| over
  everything 2.8e-14 (trivial patch, dense rounding); 1.6e-14 on non-trivial patches.
  Case 7 (11 qubits) running (logs/04e_validate_small_case7_k11.log).

- 19:31 case 7 (11 qubits, 1 more P site): f = 0.936388503630, max |dev| 1.6e-14
  (logs/04e_validate_small_case7_k11.log). Validation complete: 7 patches, max |dev| 2.8e-14.
- 19:33 SHA256SUMS.txt written (code/, logs/, circuit). logs/ckpt_L2_nopeep.slices.jsonl is an
  orphan of the killed run C (its path file was removed); ignore it.

## FINAL SUMMARY (19:33)
- f(L=3) = 0.8224883181962652 (+4.9e-16 i) run A; 0.8224883181962643 (+1.3e-16 i) run B.
  Claim 0.822488318196 -> difference +2.6e-13, i.e. agreement to all 12 quoted decimals.
- f(L=2) = 0.8972543929065162 (peephole) / 0.8972543929064776 (no peephole); claim
  0.897254392906 (|diff| 5.2e-13 / 4.8e-13; the value sits on the 12th-decimal rounding edge).
- f(L=1) = 1 exactly (claim 1). Structure U2 = Z_S U1 Z_S confirmed gate for gate; P and S disjoint.
- Not done: L=3 without the peephole (all paths found 2^56-2^86 flops); light-cone pass checked
  only on patches; the submitter's 12-qubit sub-circuit (subset unspecified) not reproduced.
- Submitter's code NOT read at any point.
