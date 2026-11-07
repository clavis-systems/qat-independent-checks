# Recheck of issue #256 (operator_loschmidt_echo_49x648) - notes
Clock: PC = UTC+1; container clock UTC. Versions: logs/00_versions.txt (python 3.13.18, numpy 2.5.4, quimb 1.15.0, cotengra 0.8.3).

## Blindness / independence - PLAIN STATEMENT
- Phase 1 (blind) ended WITHOUT a value: logs/05_phase1_result.log (22:47 PC, 6/11). Without the CZ-pair peephole the network
  needs 2^34 width / ~1e13 flops (or 1e15.65 sliced), infeasible here.
- Afterwards I READ verifica_256/ (NOTE.md, tn.py, pipeline.py, structure.py, logs) and then added my own CZ-pair peephole
  (code/peep.py) because their NOTE showed it is what makes the problem tractable. From that point my final value is NOT blind
  and NOT independent of that simplification (the code is my own, the idea is theirs). I did not read the submitter's repo or
  the GitHub issue comments; I read issue_256.json body only in Phase 2.
- I did not look at their numeric result before my contraction finished (I had seen NOTE.md's 0.82248831819626 earlier, though).

## Result (own code: parse.py, ole.py, tnbuild.py, peep.py, mmexec.py)
f = 2^-n Tr(O C^dag O C), C = U^dag V U (file = U, V, U^dag: second half is the exact gate-wise inverse of the first half,
checked in code/structure.py -> logs). V = exp(-i 0.15 X) on 24 sites (checked numerically). O = Z52 Z59 Z72.
Pipeline: peephole U alone (2268 -> 685 ops), forward light cone with Z-tracking (288 ops), backward prune from V sites,
build W V^dag W V loop, second peephole, quimb full_simplify (hyperedges), cotengra greedy path, slice 8 indices (256 slices, max 2^23).
- Path A (cotengra greedy, seed 5): f = 0.8224883181962831 - 1.6e-15 i   (logs/11_resultA.log)
- Path B (cotengra greedy, seed 6): f = 0.8224883181962837 - 1.7e-15 i   (logs/11_resultB.log)
- Previous agent: ...962652 / ...962643; difference to mine 1.8e-14 (abs). Submitter 0.822488318196: agrees to all 12 digits.
Cost: 1098 tensors, 44 qubits, 256 slices x ~4.5 s on 1 core each, 4 processes, ~5 min per path, peak RSS ~770 MB each;
path search ~220 s each. Slice-sum with math.fsum. Slice 0 of my executor equals quimb's contract_slice for the same tree to 1e-15 rel.
Validation: dense vs TN (unpruned, pruned, peephole v2) on real-gate subcircuits of 8-10 qubits: max abs err 2.4e-14
(logs/02_validate_part1.log, 02_validate_part2.log, 06_validate_peep.log).
Not done: 11-12 qubit dense checks (too slow), L=1/L=2 values.
Caveats: an early accidental empty `python3 -` call on the host (no code, no effect); a few containers were killed by pauses/Docker restarts.
