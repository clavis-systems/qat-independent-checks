# Independent checks of three Quantum Advantage Tracker submissions (#239, #240, #248)

Author: Emanuele Rizzan. The code was written by AI agents (Claude) that the author directs.
License: Apache-2.0 (see `LICENSE` and `NOTICE`).

This folder re-checks three open submissions of the
[Quantum Advantage Tracker](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io)
using only the tracker's official data files. The submitters' code was not run.

## What is checked

| Issue | Problem | Claim in the issue | What was done here | Result |
|---|---|---|---|---|
| [#240](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/240) | `guided_sparse_ground_state_problem_49Q_v1` (49 qubits) | exact ground energy −85.586625069864652 | structure of the Pauli list; product state over the three 16-qubit patches (upper side); lower bound by splitting every inter-patch term into a positive semidefinite piece plus one-patch terms (lower side) | **confirmed from both sides**: lower bound −85.5866250698646 (LOBPCG) / −85.58662506986435 (ARPACK), explicit state at −85.58662506986461 and −85.58662506986437 (two codes); agreement within about 3e-13 |
| [#248](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/248) | [2Fe-2S] (30e, 20o) | first step: DOCI energy −116.24197244 | DOCI (seniority-zero) diagonalization, 15,504 configurations | **first step confirmed**: −116.2419724395527 (residual 7e-11). The final value of the submission (−116.46391144) was not re-run |
| [#239](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/239) | [2Fe-2S] (30e, 20o) | −116.60560912042631 (full-space Davidson, flagged "not converged") | DMRG (block2, SU(2) singlet) on Edmiston–Ruedenberg localized orbitals, M up to 750, extrapolation in the discarded weight | **consistent within a few µEh**: variational −116.6055954692 (13.7 µEh above), extrapolated −116.6056052 to −116.6056071 (2–4 µEh above, uncertainty about 2 µEh). Nothing lower than the issue value was found. Digits below the µEh are not checked |

### #240 lower bound in one paragraph

Every term that couples two patches has the form `c Z_a X_b`. For a sign `s` and `kappa = 1`,
`B = |c| (I + s Z_a)(kappa I + t X_b)` with `t = s * sign(c)` is positive semidefinite, and
`c Z_a X_b = B − |c| kappa − |c| t X_b − |c| s kappa Z_a`. So `H = sum_p R_p + C + sum B`, where each `R_p` acts on one
patch only, and `lambda_min(H) >= sum_p lambda_min(R_p) + C`. Each `R_p` is diagonalized in its own 2^16 space with two
different solvers. The product of the patch ground states has `sum <B> = −8e-15`, so it attains the bound.
The script checks the identity term by term on the Pauli expansion (maximum coefficient error 0.0), checks that every
`B` has minimum eigenvalue 0, and tests the bound on 20 random 10-qubit Hamiltonians against dense diagonalization.
Second decomposition (same signs s, different kappa, so different pieces B and R_p): with `kappa = 1.01` the bound is −85.58662506986457 (ARPACK) / −85.58662506986461 (LOBPCG),
again attained (`logs/lower_bound_240_kappa1.01.*`).
Limit: the bound is exact up to floating-point rounding and up to the iterative solvers finding the lowest eigenvalue
of each patch (two solvers with different random starts agree to 1.3e-13; Lanczos residuals ≤ 1.5e-13).

## Independent re-checks (separate agent sessions, own code)

Two more agent sessions redid the main numbers with code they wrote themselves, before reading the code above.
Their code and logs are in `rechecks/` (comments in these scripts are in Italian).

| Check | Folder | Result |
|---|---|---|
| #240 lower bound | `rechecks/v4` | same kind of positive pieces, but all 2^10 = 1,024 sign choices tried: best bound −85.5866250698647; the product of the patch ground states has energy −85.58662506986457 (difference 1.3e-13). Only Pauli symmetry: X_48 (kernel of dimension 1 over GF(2)). Self-test on 12 random 9-qubit Hamiltonians: the bound never exceeds the exact minimum |
| #240 energy | `rechecks/v3` | product state at −85.58662506986462 from three different starts; residual ‖(H−E)ψ‖ = 4.7e-14, so it is an exact eigenstate |
| #248 DOCI | `rechecks/v3` | −116.2419724395529 (Lanczos and dense LAPACK); on three small sub-problems the DOCI matrix matches the seniority-zero block of a full FCI matrix built with PySCF to 2.1e-14 |
| #239 DMRG | `rechecks/v3` | own Edmiston–Ruedenberg localization, orbital order from the block2 genetic algorithm: variational −116.6055920667 at M = 750, extrapolated 2–5 µEh above the issue value; nothing lower |

Limits of the re-checks:

- `rechecks/v3/code/dmrg_239.py` ran all sweeps, then failed while writing its JSON file (`sweep_bond_dims` is not an
  attribute of the block2 DMRG object; see the end of `rechecks/v3/logs/dmrg_239_gaopt.log`). The numbers in
  `dmrg_239_gaopt.json` were read from the printed log. The script has not been fixed.
- The LDLᵀ inertia check in `rechecks/v4/code/controlli_v4.py` (a non-iterative check of the patch minima) did not
  finish in the time available, so the v4 bound has the same weak point as the main one: it relies on iterative
  solvers finding the lowest eigenvalue of each patch.
- All sessions are Claude agents directed by the same person. "Independent" here means separate sessions with
  separately written code, not separate people.
- To re-run: each script takes the path of the tracker data file as its first argument (see its first lines).
  `limite_v4.py` writes `v4_best.json` in the current directory and `controlli_v4.py` reads it from there.

## Data (tracker commit `4d27f88`, 3 Nov 2026)

| File in the tracker repository | SHA-256 |
|---|---|
| `data/variational-problems/hamiltonians/2fe_2s/2fe_2s_30e_20o.fcidump` | `bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7` (same as written in #239) |
| `data/variational-problems/hamiltonians/guided_sparse_ground_state_problem/49Q_v1/49Q_v1_hamiltonian.npy` | `c15ad8d74b5a41508a7cf99c389ad8a45ac50206a0eb7a992a72a3e11bf07566` |
| `.../49Q_v1/49Q_v1_initial_state_index.npy` (not needed by the checks) | `dcdb99c48ad035351f3ba0361d02e3fc8c6f6a02b47200796b293c47f1f3ad4b` |

The data files are not copied here. Get them from the tracker:

```bash
git clone https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io tracker
```

## How to re-run (one command)

```bash
pip install -r requirements.txt
```

```bash
bash run_all.sh tracker
```

```bash
bash run_all.sh tracker --dmrg
```

The first command runs the #240 and #248 checks; `--dmrg` adds the #239 DMRG. Results go to `out/`; compare them
with `logs/`. Set `PYTHON=...` to choose the interpreter.

Tested on 4 Nov 2026 on a fresh copy of this folder (`logs/run_all_dmrg_test_2026-11-04.log`; that copy did not
yet have the `kappa = 1.01` line): all #240 and #248 numbers came out the same as in `logs/`; the DMRG gave
−116.6055954648 at M = 750 (logged run: −116.6055954692; the DMRG starts from a random state) and extrapolations
within 0.004 µEh of the logged ones; the orbital rotation was identical to `logs/U_er.npy`.

## Versions, time and memory

Python 3.13; numpy 2.5.3 and scipy 1.18.1 for the #240 lower bound run; block2 from pip on 4 Nov 2026 for the DMRG
(exact version not recorded; 0.5.4 was the latest on PyPI that day). All runs on 4 CPU cores, about 15 GB of RAM.

| Step | Script | Time | Peak memory |
|---|---|---|---|
| #240 structure, symmetries | `struttura_240.py`, `symmetries_240.py` | < 2 s | small |
| #240 product state | `prodotto_240.py` | about 32 s | not measured |
| #240 lower bound (with sign search and second solver), each kappa | `lower_bound_240.py` | about 120 s | 620 MB |
| lower bound self-test | `lower_bound_240.py --selftest` | about 21 s | small |
| #248 DOCI | `doci_2fe2s.py` | about 3 s | small |
| #239 orbital localization | `localizza_er.py` | seconds | small |
| #239 DMRG, M up to 750 | `dmrg_2fe2s.py` | 461 s (518 s in the re-run) | 7.7 GB |

## Files

- `code/` — all scripts. `leggi_fcidump.py` is the FCIDUMP reader used by the #239 and #248 scripts.
  `gsgsp_240.py` is a first #240 attempt that was a dead end (the X masks have rank 49, so the subspace reached
  from the guide state is the full space); it is kept for completeness and not run by `run_all.sh`.
- `logs/` — outputs of the runs reported above. `U_er.npy` is the orbital rotation used by the DMRG.
- `rechecks/` — the two separate re-checks described above (code and logs).
- `run_all.sh`, `requirements.txt`, `LICENSE`, `NOTICE`.

Notes: comments and printed messages of the scripts for #239, #248 and the first #240 scripts are in Italian
(the language of the first session); `lower_bound_240.py` and `symmetries_240.py` are in English. In the logs,
absolute paths of the machine were shortened to `logs/U_er.npy`.
In code comments, C2 and C3 are the sessions that wrote the main code, V3 and V4 the re-check sessions.

## What was not done

- #239: no second full-space method (the Davidson run in the issue needed 9 hours on 64 threads and 512 GiB).
  The singlet–triplet gap used for a Temple-type estimate (about 1.4 mEh, `dmrg_SU2_ER_M500_tripletto.*`) is itself a
  DMRG number, so that estimate is not a proof.
- #248: the 10 selection rounds after the DOCI step were not re-run.
- #240: the energy was not in dispute (SKQD, Truncated Arnoldi and DMRG entries in the table give the same value);
  this check confirms the value and the lower bound, not the runtime claim.
