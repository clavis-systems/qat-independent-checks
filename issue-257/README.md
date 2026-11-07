# Check of issue #257: SU(2) hadron dynamics with l_i = 300 (problem of issue #254)

Issue: [#257](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/257)
(dylanneve1, qsim-lab). The problem is the one posed in
[#254](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/254)
(IBM Research India): 120 qubits, 60 lattice sites, 20 Trotter steps, observable n_f(t) = N_MID(t) − N_SCV(t).

This is not an official verification; see the main [README](../README.md). The code was written by AI agents
(Claude) that the author directs. The submitter's code was not read or run.

## Result

- Claim in #257: **n_f(20) = 0.10285 ± 0.00001** (N_SCV = 57.96720, N_MID = 58.07005), from an MPS (TEBD) with
  bond dimension χ = 2800.
- Value in #254 from Pauli propagation: 0.0798 (atol = 1e-5); IBM hardware: 0.0995 (3 Aug 2026), −0.11 (6 Oct 2026).
- Done here: own MPS code at five bond dimensions χ (`code/`, `results/`, `logs/summary_final.txt`).

| χ | N_SCV(20) | N_MID(20) | n_f(20) | discarded weight at the end (SCV / MID) |
|---|---|---|---|---|
| 256 | 57.900319 | 58.011458 | 0.111139 | 0.35 / 0.32 |
| 512 | 57.886589 | 58.002398 | 0.115809 | 7.7e-2 / 6.7e-2 |
| 1024 | 57.965419 | 58.068700 | 0.103281 | 2.8e-3 / 2.4e-3 |
| 2048 | 57.967164 | 58.070023 | 0.102859 | 1.2e-4 / 1.0e-4 |
| 3072 | 57.967194 | 58.070046 | 0.102852 | 3.5e-6 / 2.9e-6 |

- **n_f(20) = 0.10285, with an estimated uncertainty of 0.00002, which reproduces #257.** From χ = 2048 to 3072
  the value changes by −7e-6, and a linear fit in the discarded weight gives 0.102846; the uncertainty is a
  judgement from these two numbers, a few times the last change. The sequence is not monotone at small χ; only
  the last steps show convergence.
- At χ = 3072, N_SCV and N_MID agree with the per-step lists in #257 at all 20 steps within 6e-6 (#257 gives 5
  decimals); n_f, their difference, within 1e-5.
- The Pauli-propagation value of #254 (0.0798) is 0.023 lower. #257 already says that it is not converged; this
  check reproduces #257's numbers with separately written code. The 3 Aug hardware value is 0.0033 below; the
  hardware charge error δQ given in #254 is 0.2-0.27, so that closeness is not a precision test.
- The total charge Q = 60 holds within 1.4e-14 at every χ, but only because the MPS conserves it by construction;
  it says nothing about the truncation error.

## Method and tests (main check)

- Own numpy code (`code/lsh.py`, `code/run_mps.py`): the SWAP network is undone as a permutation of wires, the gates
  of each step are fused exactly into 4-dimensional rung sites, and a block-sparse MPS with the two conserved
  charges (U(1)×U(1)) is evolved with SVD truncation to χ. No TeNPy or quimb.
- Own QASM parser; on all four circuit files it gives the same operation list as qiskit.qasm2
  (`logs/check_parser.log`).
- Test A (`logs/testA_window.log`): 5 windows of 8-10 rungs (16-20 qubits), dense state vector applying the gates
  one by one versus the fused MPS: largest difference 5.7e-14 over all 20 steps.
- Test B (`logs/testB_sector.log`): the full 120-qubit circuit with the SWAPs applied physically, solved exactly in
  sectors with 2 and 3 particles (5 initial states, not the benchmark states): largest difference 1.4e-13. This tests
  the SWAP tracking and the output-wire mapping.
- Definitions (`logs/verify_obs.log`): n_f(20) computed with the wire-based definition of
  [#250](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/250) equals the
  site-based definition of #254 to 3e-15.
- Known case: on the older instance l_i = 12 (circuits of #149), χ = 256 / 512 / 1024 / 2048 give 0.327065 /
  0.146149 / 0.116162 / **0.116871**, versus 0.11686 ± 0.0001 in
  [#255](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/255) and
  0.1169 ± 0.0007 in #250. This instance is less converged here (the change from χ = 1024 to 2048 is +7.1e-4).
- Time and memory per circuit: χ = 1024 about 405 s and 0.40 GB (4 threads); χ = 2048 about 1200 s and 1.23 GB;
  χ = 3072 2255 s (SCV) and 2879 s (MID) and 2.6 GB (8 threads, on a laptop shared with other jobs).

## Re-check (`recheck/`)

A separate agent session with its own code and a differently organised method: an MPS with one site per qubit
(120 sites, no symmetry), the SWAP gates applied physically, χ = 256 (χ = 128 for comparison).

- It was asked to compute first without reading the other numbers, and logged its values at 10:50:21 local time on
  7 Nov 2026 (`recheck/logs/PHASE1_RESULT.txt`); its audit scripts were all written after that time.
- Its n_f for steps 1-12: 3.650671, 2.715205, 1.480150, 0.286461, −0.589960, −1.007995, −0.973485, −0.610328,
  −0.107441, 0.345667, 0.618244, 0.671375.
- Its own discarded weight says how far these can be trusted: below 1e-8 up to step 8, 8e-8 at step 9, 1.3e-6 at
  step 10, 1.6e-5 and 1.3e-4 at steps 11-12. It stopped at step 13; this kind of MPS would need a much larger χ for the later
  steps.
- N_SCV and N_MID agree with the main check and with #257 within 1e-5 at steps 1-10; at steps 11-12 the gaps are of
  the size of its own discarded weight.
- So the later steps, and n_f(20) itself, rest on the main check and on #257, which use the same kind of method
  (rung sites, U(1)×U(1) MPS) with separately written code.
- Audit of the main check (from the re-check's final report, `recheck/REPORT_recheck.md`; its audit scripts in
  `recheck/code/` print their results and wrote no log): its layout of wires and sites, and its definitions of
  N_SCV, N_MID and n_f, match #254; the table in `logs/summary_final.txt` is what `results/*.json` contain; the
  convergence claim is justified; no statement of `NOTES_main.md` without a log was found.

## The initial state MID differs from the `meson` state of #149

#254 says that the state labelled MID is the one called `meson` in the file names of #149, and that the initial
states are "the same as in #149". The circuit files do not agree with this. Sites are counted from 0; in the
initial layout, qubit w belongs to site 2·(w div 4) + (w mod 4) mod 2 (`NOTES_main.md`, line 13; the re-check found
the same map), so qubits 57 and 59 are the two occupations (n_i, n_o) of site 29, and qubits 60 and 62 those of
site 30. The SCV files set every odd qubit.

- `SU2_T20_MID.qasm` (l_i = 300): a second `x` on qubits 57 and 59 undoes the vacuum, and qubits 60 and 62 get an
  `x`. Site 29 goes from 11 to 00 and site 30 from 00 to 11.
- `x_100_meson.qasm` (l_i = 12): only qubit 59 is undone and qubit 60 gets an `x`. Site 29 goes from 11 to 10 and
  site 30 from 00 to 10.
- With #254's definition of n_f(r) (n_i + n_o on even sites, 2 − n_i − n_o on odd sites), n_f at t = 0 is 4 for
  the l_i = 300 files and 2 for the l_i = 12 files (`results/l300_MID_chi3072.log` and
  `results/l12_meson_chi2048.log`, "step 0"; N_SCV(0) = 0). Both states have Q = 60. Consistently, n_f after the
  first step is exactly twice as large for l_i = 300 (3.650671) as for l_i = 12 (1.825336).

#254 also says that the exact initial strings are fixed by the QASM files, so the files define the problem and
#257 computed the right quantity. It is not known here whether the different state is intended; if it is, the
text of #254 should say so, and the two instances are not directly comparable.

## Limits

- The only approximation in the main check is the SVD truncation; convergence is judged from five values of χ.
  This is a deterministic calculation, not a sample.
- The main check uses the same family of method as #257 (one site per rung, U(1)×U(1) TEBD), by design.
  Independence is in the code, not in the method; the re-check uses a different organisation but only reaches
  step 12.
- The circuits are simulated as given (Trotterized). Nothing here is about the continuum theory.
- Both agent sessions are Claude agents directed by the same person; the main one knew #257's value.
- Slips, for the record (`NOTES_main.md`, `REPORT_main.md`): once a piece of shell text reached the host's Python
  by mistake and stopped with a syntax error, without running anything. Runs were interrupted by pauses and restarted from per-step checkpoints (the
  χ = 3072 runs at step 14). Clock times in `NOTES_main.md` before 16:50 on 6 Nov are estimates.

## Data

- l_i = 300: `qasm_T20.zip` attached to #254 (https://github.com/user-attachments/files/33163937/qasm_T20.zip,
  SHA-256 `6edc750533b352a677129a1edbb5773ba387500a591d1c53a2093e186d9f303f`). Not copied here.
  - `SU2_T20_SCV.qasm` `1ebf7b2287fba9af08f64faa038fa0bc8035f0461261df0a0072d7ec94136566`
  - `SU2_T20_MID.qasm` `793e62350fe2f127a741f679a158810dfeceafc56f62bd4c257ce6926786e9b2`
- l_i = 12: `data/observable-estimations/circuit-models/su2_hadron_dynamics_lsh/x_100_SCV.qasm` and
  `x_100_meson.qasm` in the tracker repository at commit `4d27f88`.
  - `x_100_SCV.qasm` `af269dd7aeacfe6c9fef00bdab3bd1d0a641f5ee99be67c56ebf73a15791b6da`
  - `x_100_meson.qasm` `28b3d82a2572dac60110c94a0f58a945a117d414e94446866e0d87cc5121b4f9`

## How to re-run

All runs were made in Docker (`python:3.13-slim`), CPU only.

- Main check: numpy 2.5.4, scipy 1.18.1, psutil; qiskit 2.6.2 only for `check_parser.py`.
  `python run_mps.py --qasm FILE --chi CHI --out OUT.json [--steps 20] [--ckpt FILE]`;
  `code/ladder.sh` and `code/seq.sh` run several χ or several circuits in a row and expect the folders `/code`,
  `/res` and `/ckpt` in the container. `code/summarize.py` builds `logs/summary_final.txt` from `results/`.
- Re-check: numpy 2.5.4 only. `recheck/code/dr.sh` shows the container and its mounts, with the machine paths
  replaced by `<work>` (`<work>/ricontrollo_257` is this `recheck/` folder, `<work>/allegato_254/estratto` the
  unpacked attachment of #254, `<work>/tracker` the tracker clone).
  `python run_real.py FILE CHI TAG [NSTEPS]` runs one circuit.

## Files

- `code/`, `logs/`, `results/` (one JSON file and one log per circuit and χ), `NOTES_main.md` (the session's
  running notes), `REPORT_main.md` (its final report): main check.
- `recheck/`: the re-check (code, logs, `NOTES_recheck.md`, `REPORT_recheck.md`). In the notes, `verifica_257/`
  means the main check.
- Not included: the per-step checkpoints (about 6 GB), debugging scripts, copies of issue texts, intermediate
  summaries (`summary_1915.txt`, `summary_1935.txt`, superseded by `summary_final.txt`), a dump of circuit blocks,
  and the files of runs stopped on 6 Nov at 16:50. `NOTE.md` in the reports is `NOTES_main.md` here.
