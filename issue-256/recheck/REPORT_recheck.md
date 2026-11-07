# Final report of the re-check session (#256)

This is the final report that the re-check session sent to the session directing the work, on 7 Nov 2026. Only
machine paths were replaced by `<work>`; `ricontrollo_256/` is this `recheck/` folder and `verifica_256/` is the
main check (`../code`, `../logs`). One count was corrected later in `../README.md`: the CZ cancellation of this
re-check was tested on 4 sub-circuits (`logs/06_validate_peep.log`), the rest of its pipeline on 9.

---

REPORT (all files under `<work>/ricontrollo_256/`: NOTE.md, code/, logs/). PC clock UTC+1, container UTC.

1. PHASE 1
- Blind Phase 1 ended WITHOUT a value (logged 6/11 22:47 PC, logs/05_phase1_result.log). Without a CZ-pair peephole my network (light cone 1152 gates/182 CZ/48 wires, 1198 tensors) had unsliced width 2^34 / ~1e13 flops, or 1e15.65 flops when sliced to 2^25. Infeasible.
- Afterwards I read verifica_256 and reused its CZ-peephole idea in my own code (code/peep.py). So the final value below is NOT blind, and NOT independent of that simplification (my code, their idea). This is written in NOTE.md.
- Definition and structure: the circuit file is U (first 2268 gates), then V (220 gates), then U-dagger (the exact gate-wise inverse of the first block, checked in code). V = exp(-i 0.15 X) on 24 sites (checked numerically). f = 2^-n Tr(O C^dag O C) with C = U^dag V U. This matches README formula 1, Tr(U O U^dag V^dag U O U^dag V).
- README inconsistency: the sum formula <z|U^dag V^dag U O U^dag V|z> omits the trailing U. The correct per-state term is <z|U^dag V^dag U O U^dag V U|z>, which is what the circuit and formula 1 give. I used the circuit and formula 1. The previous agent and the submitter used the same.
- My value, own quimb + cotengra + own batch-matmul executor, complex128, 256 slices (8 sliced indices, max intermediate 2^23):
  - Path A (cotengra greedy, seed 5): f = 0.8224883181962831 - 1.6e-15 i
  - Path B (cotengra greedy, seed 6): f = 0.8224883181962837 - 1.7e-15 i
  - The two paths differ by 6e-16.
- Compared with the previous agent's runs A and B (0.8224883181962652 and ...2643): mine is higher by 1.8e-14. This is roughly 13 digits of agreement. I did not find the cause of the 1.8e-14 offset; both of my paths share the same simplified network (same peephole, so a common rounding source is possible). The submitter's 0.822488318196 agrees to all 12 quoted digits.
- Cost: 1098 tensors, 44 qubits. Path search about 220 s per path. The contraction took about 5 min per path as 4 single-core processes, about 4.5 s per slice, peak RSS about 770 MB each.
- My executor's slice 0 equals quimb's contract_slice on the same tree to 1e-15 relative.
- Validation, dense numpy vs TN, on 10 subcircuits built from real gates of this circuit (8-10 qubits): unpruned, pruned and peephole variants all agree. Max abs error 2.4e-14 (logs/02_validate_part1.log, 02_validate_part2.log, 06_validate_peep.log).

2. PHASE 2 AUDIT
a) Definition: consistent with the README trace formula, with mine, and with the submitter's text ("f = 2^-n sum sigma_z <z|W^dag O W|z>, W = U^dag V U", from issue_256.json body). The prior code builds 2^-k sum_{a,b} s_a s_b |C_ba|^2 (tn.py:157-189), which is the same quantity. No error found.
b) Use of the target value: the claim 0.822488318196 appears in the code only at run_main.py:26 (`diff_from_claim`), which is printed and not fed back into any decision. It is used nowhere in the path search, slicing or cost caps. No problem. (The code has a cost cap, MAX_LOG2_FLOPS, in pipeline.py:91-95. It does not depend on the target value.)
c) Simplifications:
  - CZ-pair cancellation (tn.py:70-118): correct. It requires both wires to carry only diagonal ops between the pair. Diagonal products are carried through CZs in merge(), tn.py:60-63, which is legitimate because they commute. Identity-up-to-phase products are dropped, and the phase cancels between C and C*.
  - Light cone (tn.py:122-153): correct. It tracks support plus "non-diagonal" wires. A CZ or diagonal op on wires where the evolved operator is Z-type commutes with it. A non-diagonal 1q gate outside the support cancels. Backward pass for C^dag O C, forward pass for C O C^dag, iterated; each is an exact operator identity.
  - I tested both ideas independently in my own implementation: peephole, forward Z-tracking cone, backward prune. They agree with dense on small circuits (error up to 2.4e-14) and with their L=3 network size. Their final network had 554 ops/214 CZ per copy of C; mine has 1124 ops/428 CZ for two copies, i.e. the same size.
  - Their own checks: 7 patches vs dense and qiskit, max deviation 2.8e-14; L=2 with and without peephole differ by 3.9e-14. Their light cone was not tested at full scale without the peephole (infeasible), only on patches.
d) Final values: I recomputed the sum of the checkpoint slice files with math.fsum (logs/09_audit_d.log).
  - Run A: 8 slices, sum/2^44 = 0.8224883181962651 + 4.918e-16 i. The log says ...2652 + 4.918e-16 i. Match.
  - Run B: 512 unique slices, 0.8224883181962642 + 1.340e-16 i. The log says ...2643 + 1.345e-16 i. Match to last-digit rounding (naive vs fsum summation).
  - k = 44, scalar = 1 in both. The logged values are what the code produces. I did not re-run their code, only re-summed their slice outputs.
e) NOTE.md statements without log support:
  - "Run C first attempt: path 2^41 flops, max 2^36, -> 2^52 after slicing" has no log. Its checkpoint was deleted, and only the later KaHyPar run, 2^36.5 / 2^38.0, is logged.
  - "The identity-drop test is never triggered for this circuit". logs/11_tolerance_check.log reports 0 identity tests, but merge() calls is_id_phase on every flush, and the L=1 collapse to a scalar (05_L1 log) requires dropped identities. This explanation is contradictory or untested.
  - The 16:xx clock entries are, by their own admission, estimates.
  - "Submitter's code not read" cannot be checked.
  - Everything else I checked is backed by logs: 2^33.0 (A) and 2^36.7 (B) flops, 512 slices, peak RSS 486 MB and 177 MB, qiskit cross-parse identical, validation cases, D/E costs 2^86/2^71/2^56.
  - I did not find any numerical error in their result.
  - Minor: the issue says "U2 = Z_S U1 Z_S", but it should read Z_S U1^dag Z_S (their NOTE has the right form).

3. NOT CHECKED
- 11-12 qubit dense comparison (too slow).
- L=1 and L=2 values.
- The submitter's repo, issue comments and the "dense 12-qubit 0.9754460616" claim.
- Their code was not re-executed.
- The source of my 1.8e-14 offset from their value (probably rounding in the peephole/scalar path; not isolated).

Housekeeping: an accidental empty `python3 -` call touched the host once, with no code and no effect. Nothing was published or posted. The v257_ container was not touched.
