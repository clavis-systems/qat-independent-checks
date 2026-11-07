# Final report of the main check session (#257)

This is the final report that the agent session of the main check sent to the session directing the work, on
7 Nov 2026. Only the machine path was replaced by `<work>`; `<work>/verifica_257/` is this folder. One point was
corrected later in `README.md`: at χ = 3072, N_SCV and N_MID agree with #257's per-step lists within 6e-6, and
n_f (their difference) within 1e-5.

---

## Verification of #257 (SU(2) LSH, l_i=300): final report

**1. Result for l_i=300 (step 20).** I used my own symmetric MPS code with exact gates.

| χ | N_SCV | N_MID | n_f(20) | discarded weight, cumulative (SCV/MID) |
|---|---|---|---|---|
| 256 | 57.900319 | 58.011458 | 0.111139 | 0.35/0.32 |
| 512 | 57.886589 | 58.002398 | 0.115809 | 7.7e-2/6.7e-2 |
| 1024 | 57.965419 | 58.068700 | 0.103281 | 2.8e-3/2.4e-3 |
| 2048 | 57.967164 | 58.070023 | 0.102859 | 1.2e-4/1.0e-4 |
| 3072 | 57.967194 | 58.070046 | 0.102852 | 3.5e-6/2.9e-6 |

- **Converged: n_f(20) = 0.10285 ± 2e-5.** The last change (χ 2048→3072) is −7e-6. The ± is a few times that change.
- **Q = 60** to 1e-14, but only because my MPS conserves charge by construction. So Q says nothing about truncation error.
- **The result strongly supports #257 (0.10285).** n_f matches to under 1e-5; N_SCV and N_MID match to within 6e-6, at all 20 steps.
- **The Pauli-propagation value in #254 (0.0798) is not converged.** It is 0.023 too low, about 3000× my last χ change. Its totals are off by −0.010 (SCV) and −0.033 (MID).
- **QPU:** the 3 Aug value (0.0995) is 0.0034 below mine, but its charge error δQ is 0.2–0.27, so that agreement is not a precision test. The 6 Oct value (−0.11) is far off.

**2. Validation**
- **Parser:** my QASM parser gives the same op list as qiskit.qasm2, op by op, on all 4 files.
- **Test A:** 5 windows of 8–10 rungs (16–20 qubits). A dense statevector applying gates one by one vs my fused MPS: max |Δn| = 5.7e-14 over all 20 steps.
- **Test B:** the full 120-qubit raw circuit with SWAPs applied physically, exact in the 2- and 3-particle sectors, 5 initial states: max |Δn| = 1.4e-13. This checks the SWAP tracking and the output-wire mapping.
- **l_i=12** (x_100_SCV/meson), χ 256/512/1024/2048: 0.3271 / 0.1461 / 0.1162 / **0.116871**.
  - vs #255 (0.11686 ± 1e-4): differs by 8e-6.
  - Also consistent with #250 (0.1169 ± 0.0007) and Roland-Farrell (0.1176 ± 0.0022).
  - Per step vs #255: within 6e-6 at all 20 steps.
- **Definitions:** n_f with #250's wire-based definition equals #254's site-based one to 1e-15.

**3. Method and resources**
- **Method:** own numpy code. It undoes the SWAP network as a permutation, fuses gates exactly into 4-dim rung sites (59 two-site gates per step), and runs a block-sparse U(1)×U(1) MPS with SVD truncation to χ.
- **Libraries:** no TeNPy or quimb. Docker python:3.13-slim, numpy 2.5.4, scipy 1.18.1; qiskit 2.6.2 only for the parser cross-check.
- **Time and memory per circuit:**
  - χ=1024: 405 s, 0.40 GB (4 threads)
  - χ=2048: 1200 s, 1.23 GB (4 threads)
  - χ=3072: 2255–2879 s, 2.6 GB peak (8 threads, laptop shared with other jobs)
- **Input SHA-256:** SCV 1ebf7b2287fba9af08f64faa038fa0bc8035f0461261df0a0072d7ec94136566, MID 793e62350fe2f127a741f679a158810dfeceafc56f62bd4c257ce6926786e9b2.

Files that could support a public write-up are in `<work>/verifica_257/`:
- `NOTE.md`
- `code/`: qasm_parse.py, lsh.py, run_mps.py, test_window.py, test_sector.py, check_parser.py, verify_obs.py, summarize.py, ladder.sh, seq.sh
- `logs/`: sha256_inputs.txt, check_parser.log, testA_window.log, testB_sector.log, verify_obs.log, summary_final.txt
- `results/*.json` and `results/*.log`

Not for publication: `ckpt/` (7 GB of checkpoints), `code/debug/`, `results/stopped_1650/`.

**4. Limits and independence**
- **Only one approximation:** SVD truncation. Convergence is judged from 5 values of χ; this is a deterministic calculation, not a sample.
- **Independence is in the code, not the method.** I used the same algorithm family as #257 (rungs grouped into one site, U(1)×U(1) TEBD), as the brief suggested. I never read or ran qsim-lab. The #257 and #255 numbers were used only for comparison after I had my own.
- **Scope:** it simulates the circuits as given (Trotterized, approximate Hamiltonian), not the full theory.
- **Point the arbitration did not ask about:** the l_i=300 MID initial state moves two fermions (site 29: 11→00, site 30: 00→11). The l_i=12 meson moves only one. That contradicts #254's claim that the initial states are the same as in #149.
- **Incidents:**
  - One stray heredoc sent shell text to the host Python. It failed with a syntax error and nothing ran.
  - Runs were interrupted by Docker restarts and your pauses. They resumed from per-step checkpoints; the χ=3072 runs resumed at step 14.
  - Times in NOTE.md before 16:50 are approximate.
- **Second check still needed:** no independent recheck (Codex or a separate agent) has been done yet.
