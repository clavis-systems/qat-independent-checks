# Check of issue #256: exact value of the operator Loschmidt echo 49x648

Issue: [#256](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/256)
(dylanneve1, qsim-lab), circuit instance `operator_loschmidt_echo_49x648`.

This is not an official verification; see the main [README](../README.md). The code was written by AI agents
(Claude) that the author directs. The submitter's code was not read or run.

## Result

- Claim in the issue: f = **0.822488318196**, exact tensor-network contraction (quimb + cotengra/KaHyPar).
- Done here: two exact contractions written separately by two agent sessions (the second one is described under
  "Re-check" below), each with two different contraction paths.

| Who | Code | f |
|---|---|---|
| issue #256 | quimb + cotengra/KaHyPar | 0.822488318196 |
| main check (`code/`) | own network and contraction, path from cotengra (run A) or opt_einsum (run B) | 0.8224883181962652 / 0.8224883181962643 |
| re-check (`recheck/`) | own code with quimb, two cotengra greedy paths | 0.8224883181962831 / 0.8224883181962837 |

- **The issue value is reproduced to all 12 digits it gives.** The two checks here differ from each other by
  1.8-1.9e-14; the cause of this offset was not isolated. The imaginary parts are below 2e-15.
- Also checked: f(L=1) = 1 (issue: 1) and f(L=2) = 0.8972543929065162, which is 5.2e-13 from the issue's
  0.897254392906 (rounded to 12 digits, the value here would end in …907), on circuits for L = 1 and 2 rebuilt from
  the Trotter step found in the L = 3 file; 4756 gates; and the structure the issue describes: the file
  is C = U† V U, with V = exp(−i 0.15 X) on 24 sites, and the second half of U is the inverse of the first half with
  the sign of the `rx` angles flipped on the 13 sites with b = 0.25.

## Definition used

f = 2^-n Tr(O C† O C), with C the whole circuit in the file and O = Z52 Z59 Z72. This is the trace formula of the
tracker README for this instance (`data/observable-estimations/circuit-models/operator_loschmidt_echo/README.md`,
line 15) and the formula written in the issue.

Two small text points, neither of which changes the number:

- In the tracker README, the sum formula (line 27), ⟨z|U† V† U O U† V|z⟩, lacks the final U; the trace formula
  (line 15) and the circuit file agree with ⟨z|U† V† U O U† V U|z⟩.
- The issue writes "U2 = Z_S U1 Z_S" without defining U1 and U2. Reading them as U = U2 U1 (U1 applied first, as
  in the tracker README's U = (U_{b−α}†)^L U_b^L), its own conclusion W = Zt V Zt with Zt = U1† Z_S U1 needs
  U2 = Z_S U1† Z_S, so it looks like a missing dagger. The check of the circuit (`logs/02_structure.log`, point
  (b): the second half of U is the inverse of the first half with the `rx` signs flipped on S) agrees with this
  form. `NOTES_main.md` (lines 26 and 138) repeats the issue's wording.

## Main check (`code/`, `logs/`)

- Own OpenQASM 3 parser (`qasm_parse.py`); qiskit.qasm3 parses exactly the same 4756-gate list.
- Two exact simplifications before the contraction (`tn.py`): pairs of CZ gates separated only by diagonal gates
  cancel, and gates outside the light cone of O are dropped. The result is a doubled network on 44 qubits with
  1108 tensors. Own slicing and own contraction (numpy batched matrix products).
  - Run A (`logs/06_L3_runA_ctg_greedy_labels.jsonl`): path from cotengra (greedy + labels), 2^33.5 flops in total,
    8 slices; 205 s path search + 14 s contraction, peak memory 486 MB.
  - Run B (`logs/08_L3_runB_oe_randomgreedy.jsonl`): path from opt_einsum RandomGreedy, 2^38.9 flops in total,
    512 slices; 181 s + 311 s, peak memory 177 MB.
- Validation (`logs/04c`, `04d`, `04e`): 7 groups of qubits (six of 10 qubits, one of 11), each keeping every gate
  of the real L = 3 circuit that acts inside the group (942-1141 gates; in one group the network reduces to f = 1).
  For each group, f from an own dense unitary, from qiskit's `Operator` and from six network variants (with and
  without each simplification, forced slicing, another path) agree within 2.8e-14.
- The CZ cancellation was also checked at full size at L = 2: 0.8972543929065162 with it and 0.8972543929064776
  without it (`logs/07b_L2_no_peephole_lc_only_contract.jsonl`).
- Tolerance (`logs/11_tolerance_check.log`): the largest value treated as zero when deciding whether a gate is
  diagonal is 2.4e-16, the smallest treated as non-zero is 0.149, so no decision is borderline.
- Infeasible attempts are kept: at L = 3, without the CZ cancellation (light cone only) the path found cost 2^86
  flops unsliced (`logs/09_*`); without the light cone (CZ cancellation only) 2^71 and 2^56 (`logs/10_*`,
  `logs/10b_*`). `NOTES_main.md` (line 139) puts both cases under "without the peephole".

## Re-check (`recheck/`)

A separate agent session with its own code (quimb for building and simplifying the network, cotengra for the path,
own executor for the slices).

- It was asked to compute f first without reading the expected value. **That part ended without a value**
  (`recheck/logs/05_phase1_result.log`): without the CZ cancellation its network was too costly for the machine.
  The best unsliced path needed about 2^43.5 flops but intermediate tensors of 2^34 entries (about 275 GB); sliced
  to 2^25 entries it needed about 2^52 flops (`recheck/logs/03_*`, `05_*`).
- It then read the main check, re-implemented the CZ cancellation in its own code (`recheck/code/peep.py`) and ran
  the full contraction: 1098 tensors, 44 qubits, 256 slices, about 1e11.4 flops per path, about 5 minutes per path
  as 4 processes, peak memory about 770 MB each (`recheck/logs/08_*`, `10_*`, `11_*`). Each of the 4 processes
  prints the sum of its own 64 slices as "RESULT f" in `10_mm*` (a partial sum); `11_result*` is the total.
- So the re-check is **not blind** (it had seen the main check's notes, which contain the value) and **not
  independent of the CZ cancellation**: its code is its own, the idea is the main check's. The cancellation is an
  exact identity; the re-check reasoned through it and tested its version against dense calculations on 4
  sub-circuits of 8-10 qubits built from real gates (largest error 2.4e-14, `recheck/logs/06_validate_peep.log`);
  the rest of its pipeline was tested on 9 such sub-circuits (largest error 2.8e-14, `recheck/logs/02_validate_*`).
- Audit of the main check (from the re-check's final report, `recheck/REPORT_recheck.md`; only point (d) has its
  own log):
  - the definition of f is the same as the README trace formula and the issue's;
  - the expected value appears in the main code only to print a difference (`code/run_main.py`, `diff_from_claim`)
    and is not used for any decision;
  - the two simplifications are exact operator identities, and the code applies them correctly;
  - (d) re-adding the slice values of runs A and B gives the logged values (`recheck/logs/09_audit_d.log`);
  - in `NOTES_main.md`, a first "run C" attempt (2^41 flops) has no log and the clock times before 16:50 are
    estimates. The remark that the identity-drop test "is never triggered" agrees with
    `logs/11_tolerance_check.log` (0 identity tests), but the re-check found it hard to reconcile with the L = 1
    network reducing to a scalar; this was not resolved. None of these points concerns the L = 3 result.

## Limits

- At L = 3 both values here rely on the two simplifications above. The CZ cancellation was checked at full size only
  at L = 2; the light-cone step only on the small groups.
- The two simplifications are exact identities, so the small-group tests check their implementation, not the
  idea. The small groups are sub-circuits, not the real circuit; a run that was to include groups of 12 qubits was
  stopped as too slow before reaching them (`logs/04a_*`).
- There is no rigorous error bound, and the floating-point round-off over the 256-512 slices was not bounded; the
  spread between paths is about 1e-15, between the two codes 1.8-1.9e-14.
- Both agent sessions knew the claimed value. All sessions are Claude agents directed by the same person.
- The 12-qubit check in the issue (0.9754460616) was not reproduced, because the qubit subset is not given.
- Two slips, for the record (from the sessions' reports, `REPORT_main.md` and `recheck/REPORT_recheck.md`): the
  main session once used the host's Python to edit the text of two source files (no computation), and the
  re-check session once started the host's Python with no code. Everything else ran in Docker.

## Data

- Circuit: `data/observable-estimations/circuit-models/operator_loschmidt_echo/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm`
  in the tracker repository at commit `4d27f88` (the same file as in the issue).
  - SHA-256 `4aae5089337253921cd70df470dc52b5243e0c8410b36a6bcbb43497db48db29`

## How to re-run

All runs were made in Docker on Linux containers (`python:3.13-slim`), CPU only.

- Main check: Python 3.13, numpy 2.5.4, cotengra 0.8.3, opt_einsum 3.4.0, qiskit 2.6.2, qiskit-qasm3-import 0.6.0
  (full list in `logs/00_versions.txt`; kahypar 1.3.7 only for the L = 2 cross-check). The scripts read the circuit
  from `/circ/` and write under the current folder; the main runs used
  `python run_main.py L peep lc method seed target_log2 max_time [ctg_methods] [ckpt]` (see the first lines of
  `code/run_main.py`; the path method, seed and slicing of runs A and B are recorded in their JSON lines, and
  their paths are in `logs/ckpt_L3_run*.path.json`).
- Re-check: `recheck/Dockerfile` (numpy 2.5.4, quimb 1.15.0, cotengra 0.8.3 as recorded in
  `recheck/logs/00_versions.txt`; the Dockerfile itself does not pin versions). The scripts read the circuit from
  `/data/`; `recheck/run.sh` shows the mounts, with the machine paths replaced by `<work>`
  (`<work>/ricontrollo_256` is this `recheck/` folder, `<work>/tracker` the tracker clone).
  `run_full.py <greedy|labels|oe> <seed> <search_seconds> <log2 slice target>` builds the network, finds the
  path and saves it as a pickle; `mmexec.py <tree.pkl> <ckpt.jsonl> <n_slices|all>` contracts the slices (the
  final runs used 4 processes of 64 slices each) and `combine.py` adds the slice values. The path pickles are not
  included, so a re-run needs a new path search, which may give a different path. The audit scripts read the main
  check's folder from `/s/`.

## Files

- `code/`, `logs/`, `NOTES_main.md` (the session's running notes), `REPORT_main.md` (its final report): main
  check. `logs/ckpt_*` hold the contraction paths and the value of every slice.
- `recheck/`: the re-check (code, logs, `NOTES_recheck.md`, `REPORT_recheck.md`). In its notes, `verifica_256/`
  means the main check (this folder's `code/` and `logs/`).
- Not included: binary pickle files of the networks and paths (caches, and pickle files are unsafe to load from
  others), one orphan log of a killed run, and files named in the notes that are not needed here:
  `SHA256SUMS.txt` (hashes of the working folder), `logs/04_validate_small.log` (identical to `04c`), and
  `NOTE.md`, which is `NOTES_main.md` here.
