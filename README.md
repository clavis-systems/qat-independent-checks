# Independent checks of three Quantum Advantage Tracker submissions (#239, #240, #248)

Author: Emanuele Rizzan. The code was written by AI agents (Claude) that the author directs.
License: Apache-2.0 (see `LICENSE` and `NOTICE`).

This folder re-checks three open submissions of the
[Quantum Advantage Tracker](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io)
using only the tracker's official data files. The submitters' code was not run.

This is not an official verification. The tracker's maintainers verify submissions themselves, and the author is
not affiliated with the tracker or with the submitters.

## What is checked

| Issue | Problem | Result |
|---|---|---|
| [#240](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/240) | 49 qubits | **confirmed from both sides** within about 3e-13, assuming the iterative eigensolvers found each patch minimum |
| [#248](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/248) | [2Fe-2S] (30e, 20o) | **first step confirmed**; the final value was not re-run |
| [#239](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/239) | [2Fe-2S] (30e, 20o) | **consistent at the level of about 10 µEh**; nothing lower than the issue value was found |

Details:

- **#240**, problem `guided_sparse_ground_state_problem_49Q_v1`.
  - Claim in the issue: exact ground energy −85.586625069864652, error bound ±5.7e-14.
  - Done here: structure of the Pauli list; product state over the three 16-qubit patches (upper side); lower bound
    by splitting every inter-patch term into a positive semidefinite piece plus one-patch terms (lower side).
  - Result: **confirmed from both sides**: lower bound −85.5866250698646 (LOBPCG) / −85.58662506986435 (ARPACK);
    an explicit product state has energy −85.58662506986461 (`prodotto_240.py`, energy computed directly from the
    state). For the state built in `lower_bound_240.py`, the bound plus `sum <B>` gives −85.58662506986437; this
    uses the same eigensolver output as the bound, so it is not a separate check. All these numbers and the issue
    value agree within about 3e-13. This is the resolution of this check: it does not confirm the issue's ±5.7e-14.
- **#248**, [2Fe-2S] (30e, 20o).
  - Claim in the issue: first step, DOCI energy −116.24197244.
  - Done here: DOCI (seniority-zero) diagonalization, 15,504 configurations.
  - Result: **first step confirmed**: −116.2419724395527 (residual 7e-11). The final value of the submission
    (−116.46391144) was not re-run.
- **#239**, [2Fe-2S] (30e, 20o).
  - Claim in the issue: −116.60560912042631 (full-space Davidson; the issue reports the solver's convergence flag
    as `false`).
  - Done here: DMRG (block2, SU(2) singlet) on Edmiston–Ruedenberg localized orbitals, M up to 750, extrapolation in
    the discarded weight.
  - Result: **consistent at the level of about 10 µEh**. The variational DMRG energy at M = 750 is −116.6055954692,
    13.7 µEh above the issue value; this DMRG was itself not fully converged (the last sweep still lowered the energy
    by 4.0e-7, and block2 printed "DMRG is not converged to desired tolerance of 1e-8"). Linear extrapolations in the
    discarded weight give between −116.6056014 (all 5 points, M ≥ 250: the fit the script reports, 7.7 µEh above)
    and −116.6056071 (4 points, M ≥ 350: 2.0 µEh above); a quadratic fit gives −116.6056052 (3.9 µEh above). All
    lie above the issue value, and the spread between fits is of the same size. So this DMRG finds nothing lower
    than the issue value and agrees with it at the level of about 10 µEh; it does not check smaller differences.

### #240 lower bound in one paragraph

Every term that couples two patches has the form `c Z_a X_b`. For a sign `s` and `kappa = 1`,
`B = |c| (I + s Z_a)(kappa I + t X_b)` with `t = s * sign(c)` is positive semidefinite, and
`c Z_a X_b = B − |c| kappa − |c| t X_b − |c| s kappa Z_a`. So `H = sum_p R_p + C + sum B`, where each `R_p` acts on one
patch only, and `lambda_min(H) >= sum_p lambda_min(R_p) + C`. Each `R_p` is diagonalized in its own 2^16 space with two
different solvers. The product of the patch ground states has `sum <B> = −8e-15`, so it attains the bound.
The script checks the identity term by term on the Pauli expansion (maximum coefficient error 0.0), checks that every
`B` has minimum eigenvalue 0, and tests the bound on 20 random 10-qubit Hamiltonians against dense diagonalization.
Second decomposition (same signs s, different kappa, so different pieces B and R_p): with `kappa = 1.01` the bound
is −85.58662506986457 (ARPACK) / −85.58662506986461 (LOBPCG), again attained (`logs/lower_bound_240_kappa1.01.*`).
Limit: the bound is exact up to floating-point rounding and up to the iterative solvers finding the lowest eigenvalue
of each patch (two solvers with different random starts agree to 1.3e-13; Lanczos residuals ≤ 1.5e-13). The upper
side (an explicit state) does not depend on this, but its patch states come from the same kind of solver, so the
agreement of the two sides cannot reveal a minimum missed by both solvers. The self-tests use small patches, which
are diagonalized densely, so they test the algebra of the bound but not the iterative solvers.

## Independent re-checks (separate agent sessions, own code)

Two more agent sessions redid the main numbers with code they wrote themselves. They knew the target values (and,
for v4, the method in words) and were told to write their code before reading the main code; this order cannot be
checked from the files. Their code and logs are in `rechecks/` (comments in these scripts are in Italian).

| Check | Folder | Result |
|---|---|---|
| #240 lower bound | `rechecks/v4` | same kind of positive pieces, but all 2^10 = 1,024 sign choices tried: best bound −85.58662506986451 (`ricerca_segni.log`); recomputed with three random starts per patch: bound −85.5866250698647, and the product of the patch ground states has energy −85.58662506986457 (difference 1.3e-13, `controlli.log`). Self-test on 12 random 9-qubit Hamiltonians: the bound never exceeds the exact minimum |
| #240 energy | `rechecks/v3` | product state at −85.58662506986462 from three different starts; residual ‖(H−E)ψ‖ = 4.7e-14, so it is an eigenstate to that accuracy |
| #248 DOCI | `rechecks/v3` | −116.2419724395529 (Lanczos and dense LAPACK); on three small sub-problems the DOCI matrix matches the seniority-zero block of a full FCI matrix built with PySCF to 2.1e-14 |
| #239 DMRG | `rechecks/v3` | a repeat with separately written code, not a different method: same program (block2), own Edmiston–Ruedenberg localization (it gives almost the same orbitals), orbital order from the block2 genetic algorithm: variational −116.6055920667 at M = 750, extrapolated 2–5 µEh above the issue value; nothing lower |

Limits of the re-checks:

- `rechecks/v3/code/dmrg_239.py` ran all sweeps, then failed while writing its JSON file (`sweep_bond_dims` is not an
  attribute of the block2 DMRG object; see the end of `rechecks/v3/logs/dmrg_239_gaopt.log`). The numbers in
  `dmrg_239_gaopt.json` were read from the printed log. The script has not been fixed.
- The LDLᵀ inertia check in `rechecks/v4/code/controlli_v4.py` (a non-iterative check of the patch minima) was not
  run: it runs only when the environment variable `V4_INERZIA` is set to a group number, and `controlli.log` says "non fatto" (not
  done). So the v4 bound has the same weak point as the main one: it relies on iterative solvers finding the lowest
  eigenvalue of each patch.
- `rechecks/v4/logs/run_all_copia_pulita.log` is a run of the main `run_all.sh` from a clean copy of this folder, made
  by the v4 session; it is not a run of the v4 code.
- All sessions are Claude agents directed by the same person. "Independent" here means separate sessions with
  separately written code, not separate people.
- To re-run: the scripts that read tracker data take the path of the data file as their first argument (see their
  first lines); `estrapola_239.py` reads the JSON written for the v3 DMRG, and `autotest_v4.py` needs no data.
  `limite_v4.py` writes `v4_best.json` in the current directory and `controlli_v4.py` reads it from there.

## Data (tracker commit `4d27f88`, 3 Nov 2026)

Files in the tracker repository and their SHA-256:

- `data/variational-problems/hamiltonians/2fe_2s/2fe_2s_30e_20o.fcidump`
  - `bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7` (same as written in #239)
- `data/variational-problems/hamiltonians/guided_sparse_ground_state_problem/49Q_v1/49Q_v1_hamiltonian.npy`
  - `c15ad8d74b5a41508a7cf99c389ad8a45ac50206a0eb7a992a72a3e11bf07566`
- `.../49Q_v1/49Q_v1_initial_state_index.npy` (not needed by the checks)
  - `dcdb99c48ad035351f3ba0361d02e3fc8c6f6a02b47200796b293c47f1f3ad4b`

The data files are not copied here. Get them from the tracker, at the commit used here:

```bash
git clone https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io tracker
```

```bash
git -C tracker checkout 4d27f88
```

## How to re-run

The scripts were run on Linux: `run_all.sh` uses `sha256sum`, and `lower_bound_240.py` uses the Unix-only
`resource` module. For the fast checks (#240 and #248) only numpy and scipy are needed:

```bash
pip install numpy==2.5.3 scipy==1.18.1
```

```bash
bash run_all.sh tracker
```

For everything, including block2 for the #239 DMRG and pyscf for `rechecks/v3`:

```bash
pip install -r requirements.txt
```

```bash
bash run_all.sh tracker --dmrg
```

`run_all.sh` first checks the SHA-256 of the two data files and stops if they differ. Results go to `out/`; compare
them with `logs/`. Set `PYTHON=...` to choose the interpreter. The singlet–triplet run
(`logs/dmrg_SU2_ER_M500_tripletto.*`) is not in `run_all.sh`; it used
`DMRG_ROT=U_er.npy python3 -I code/dmrg_2fe2s.py FCIDUMP SU2 fiedler 4 SCRATCH_DIR OUT.json 500 2`.

Re-runs:

- 4 Nov 2026, on a fresh copy of this folder (`logs/run_all_dmrg_test_2026-11-04.log`; that copy did not yet have
  the `kappa = 1.01` line): the #240 numbers matched `logs/`; DOCI gave −116.24197243955263 (logged
  −116.2419724395527) with residual 6.7e-12 (logged 6.7e-11); the DMRG gave −116.6055954648 at M = 750 (logged run:
  −116.6055954692; the DMRG starts from a random state) and extrapolations within 0.004 µEh of the logged ones; the
  orbital localization reached the same value of `sum (ii|ii)` (12.6833172432544) as `logs/localizza_er.log`.
- 6 Nov 2026, in a clean `python:3.13-slim` Docker container with only numpy 2.5.3 and scipy 1.18.1 installed,
  fast checks only (`logs/run_all_docker_2026-11-06.log`): the data hashes matched; #240 lower bound
  −85.58662506986465 (ARPACK) / −85.58662506986462 (LOBPCG), with `kappa = 1.01` −85.58662506986445 /
  −85.58662506986462; product state −85.58662506986461; DOCI −116.24197243955291 (residual 2.6e-11). All are within
  3e-13 of the values in `logs/`. The run took 12 minutes on a 16-core desktop.

## Versions, time and memory

Python 3.13, numpy 2.5.3 and scipy 1.18.1 for the #240 lower bound run (recorded in its JSON). For the other runs
the versions were not recorded, except that block2 came from pip on 4 Nov 2026 (0.5.4 was the latest on PyPI that
day). The sessions ran on a cloud machine with 4 CPU cores and about 15 GB of RAM (as reported by the sessions; not
in the logs). Times are approximate; only some are in the logs.

| Step | Script | Time | Peak memory |
|---|---|---|---|
| #240 structure, symmetries | `struttura_240.py`, `symmetries_240.py` | < 2 s | small |
| #240 product state | `prodotto_240.py` | about 32 s | not measured |
| #240 lower bound (with sign search and second solver), each kappa | `lower_bound_240.py` | about 120 s | 620 MB |
| lower bound self-test | `lower_bound_240.py --selftest` | about 21 s | small |
| #248 DOCI | `doci_2fe2s.py` | about 3 s | small |
| #239 orbital localization | `localizza_er.py` | seconds | small |
| #239 DMRG, M up to 750 | `dmrg_2fe2s.py` | 461 s (518 s in the re-run) | not logged for this run (see below) |

DMRG memory: `logs/memoria_dmrg.txt` is a partial memory trace of a DMRG process in the same session (it is not
clear which run); its peak is 5.0 GB. Runs at M = 1000 used more than 13 GB and were stopped.

## Files

- `code/` — all scripts. `leggi_fcidump.py` is the FCIDUMP reader used by the #239 and #248 scripts.
  `gsgsp_240.py` is a first #240 attempt that was a dead end (the X masks have rank 49, so the subspace reached
  from the guide state is the full space); it is kept for completeness and not run by `run_all.sh`.
- `logs/` — outputs of the runs reported above. `U_er.npy` is the orbital rotation used by the DMRG.
- `rechecks/` — the two separate re-checks described above (code and logs).
- `run_all.sh`, `requirements.txt`, `LICENSE`, `NOTICE`.

Notes:

- Comments and printed messages of the scripts for #239, #248 and the first #240 scripts are in Italian (the
  language of the first session); `lower_bound_240.py` and `symmetries_240.py` are in English. In code comments, C2
  and C3 are the sessions that wrote the main code, V3 and V4 the re-check sessions.
- In the logs, paths of the machines were replaced by short placeholders: `tracker/`, `out/`, `./out/`,
  `logs/U_er.npy`, `<scratch>/` and `rechecks/v3/code/`.
- `prodotto_240.py` prints a residual of 9.5e-7 for the product state. It is computed as the square root of a
  difference of large numbers (raw value 9.1e-13), so it only shows rounding; the v3 re-check computes the residual
  directly and gets 4.7e-14.
- The RDM energy check at the end of `logs/dmrg_SU2_ER_M750.log` (−116.6052214) refers to the final MPS after the
  reverse sweeps, truncated to M = 250, not to the M = 750 state. It checks the integrals and the RDM code, not the
  M = 750 energy.

## What was not done

- #239: no second full-space method (the Davidson run in the issue needed 9 hours on 64 threads and 512 GiB), and
  no DMRG beyond M = 750 (see the memory note above). The singlet–triplet gap used for a Temple-type estimate (about
  1.4 mEh, `dmrg_SU2_ER_M500_tripletto.*`) is itself a DMRG number, so that estimate is not a proof.
- #248: the 10 selection rounds after the DOCI step were not re-run.
- #240: the energy was not in dispute (SKQD, Truncated Arnoldi and DMRG entries in the table give the same value);
  this check confirms the value and the lower bound, not the runtime claim.
