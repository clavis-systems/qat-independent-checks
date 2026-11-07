# Final report of the main check session (#256)

This is the final report that the agent session of the main check sent to the session directing the work, on
6 Nov 2026. Only the machine path was replaced by `<work>`; `<work>/verifica_256/` is this folder. Two points
were corrected later in `README.md`: the validation used 7 groups of 10 or 11 qubits with 942-1141 gates;
"U2 = Z_S U1 Z_S" needs a dagger; and "without them every path cost 2^56-2^86 flops" mixes two cases (without the
CZ cancellation 2^86, without the light cone 2^56-2^71).

---

**#256 verification: agreement confirmed (container v256 stopped, nothing posted)**

**1. Values**
- **L=3:** f = 0.8224883181962652, Im = +4.9e-16 (run A). A second run gave 0.8224883181962643, Im = +1.3e-16 (run B). The difference from the claimed 0.822488318196 is +2.6e-13, so the two agree to all 12 quoted decimals, and A and B differ by 9e-16.
- **Definition used:** f = 2^-n Tr(O C† O C) with C the QASM circuit and O = Z52 Z59 Z72. This is the tracker README trace formula, and it matches the submitter's formula. The README's sum formula as printed is missing the final U; the trace form is what agrees with the circuit file.
- **Method:** own parser and own doubled network, with diagonal gates as hyperedges. Two exact simplifications come before the contraction: CZ pairs separated only by diagonal gates cancel, and gates outside the light cone of O are dropped. The result is 44 qubits and 1108 tensors.
  - Run A: path from cotengra greedy+labels (no KaHyPar), 2^33.0 flops; my own slicing into 8 slices; my own numpy batch-matmul contraction.
  - Run B: path from opt_einsum RandomGreedy (no cotengra), 512 slices, 2^38.9 flops.
- **Other checks:**
  - f(L=1) = 1 exactly, as claimed.
  - f(L=2) = 0.8972543929065162 with the simplification and 0.8972543929064776 without it (KaHyPar path, 2^38). The claim is 0.897254392906; the difference is about 5e-13, and the value sits right on the 12th-decimal rounding edge.
  - Both parts of the submitter's structure claim hold gate for gate: U2 = Z_S U1 Z_S, and the circuit is C = U†VU. The perturbed sites and the b=0.25 sites (S) do not overlap.
  - The circuit has 4756 gates, as claimed.

**2. Validation**
- qiskit.qasm3 parses exactly the same 4756-gate list as my parser.
- I used 7 patches of 10–11 qubits, each keeping every gate of the real L=3 circuit inside the patch (942–1016 gates). Values were 1, 0.9446, 0.9834, 0.9650, 0.9996, 0.9677 and 0.9364.
- Each patch was computed by my own dense unitary, by qiskit's Operator, and by six network variants (full simplification, forced slicing, opt_einsum path, light cone only, peephole only, raw network). Max deviation was 2.8e-14.
- The tolerance check shows a gap from 2.4e-16 to 0.149 around the 1e-14 "is this zero?" threshold, so no gate classification is borderline.

**3. Versions, time, memory**
- Python 3.13.18, numpy 2.5.4 (OpenBLAS 0.3.34), cotengra 0.8.3, opt_einsum 3.4.0, qiskit 2.6.2, qiskit-qasm3-import 0.6.0, kahypar 1.3.7 (only for the L=2 cross-check and two infeasible attempts).
- All runs in Docker, 4 CPUs, 3 GB, no GPU. Circuit SHA-256: 4aae5089…db48db29.
- Run A: 205 s path search + 14 s contraction, peak 486 MB. Run B: 181 s + 311 s, peak 177 MB.
- Files for a write-up, all in `<work>/verifica_256/`:
  - `NOTE.md` and `SHA256SUMS.txt`
  - `code/` (qasm_parse, structure, tn, pipeline, dense, validate_small, run_main, probe, tolerance_check, hashes)
  - `logs/` 00–03, 04c/04d/04e (validation), 05 (L=1, L=2), 06 (run A), 08 (run B), 07b (L=2 without simplification), 11 (tolerance), plus `ckpt_*` (paths and per-slice values)
  - `logs/ckpt_L2_nopeep.slices.jsonl` is an orphan from a killed run and should be ignored.

**4. Limits**
- At L=3 the value depends on the two simplifications. Without them every path I found cost 2^56–2^86 flops, so the CZ cancellation is checked at full size only at L=2. The light-cone step is checked only on the small patches.
- The patches are sub-circuits, not the real circuit, and k=12 patches were too slow under the CPU limit.
- There is no rigorous error bound; the observed spread is about 1e-15, with about 1e-14 on the dense checks.
- The submitter's 12-qubit check (0.9754460616) could not be reproduced because the subset is not specified.
- Independence: the parser, network, simplifications, slicing and contraction are all mine, and the submitter's code was never read. What we share is numpy/BLAS, the circuit file, the definition of f and the doubled-network formulation.
- One deviation from the rules: early on I used host Python once, only to string-edit two source files; no computation ran on the host.
- The session was paused 16:50–18:28 and Docker stopped once more at 18:33. All results come from runs after 18:35.
