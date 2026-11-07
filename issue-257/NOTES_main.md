# verifica_257 — independent check of n_f(20), SU(2) LSH, l_i = 300

Goal: own exact-gate MPS (numpy, own code, U(1)xU(1) block-sparse), validate, then
n_f(20) for l_i=300 at increasing chi. Everything runs in Docker (image verifica257:py313,
built from python:3.13-slim + pip numpy 2.5.4, scipy 1.18.1, qiskit 2.6.2, psutil).

## Log (6/11/2026; times before 16:50 are approximate: I did not run `date` at each step,
## the file timestamps show everything up to the chi=256 runs was done by 16:50)
- 16:34 start. Read #254, #257, #250 (+comments), #255. SHA-256 of inputs in logs/sha256_inputs.txt.
- ~16:38 structure (logs/analyze.log): 32070 ops; x 2460 (60 prefix + 120/step), cx 11840, h 4720,
  rz 10720, swap 2330. rz angles l_i=300: +-0.15 (hop), +-0.1130625 (= theta/2, electric), +-0.03 (mass).
  l_i=12 identical except +-0.005 electric.
- Layout deduced: initial wire w -> site 2*(w//4) + (w%4)%2, chain (w%4)//2 ; final SWAP permutation
  puts site r on output wires 2r, 2r+1 (matches #250 definition). SCV = all odd wires = (00,11,00,11...).
- NOTE: l_i=300 MID initial state differs from l_i=12 meson: MID clears wires 57,59 and sets 60,62
  (site 29: 11->00, site 30: 00->11); meson clears 59 and sets 60 (one fermion moved).
- ~16:42 code written: code/qasm_parse.py (own parser), code/lsh.py (SWAP permutation tracking,
  gate fusion, U(1)xU(1) block-sparse MPS, exact gates, global SVD truncation to chi),
  code/run_mps.py (driver), tests code/check_parser.py, code/test_window.py, code/test_sector.py.
- check_parser (logs/check_parser.log): own parser == qiskit.qasm2 (2.6.2) op by op, all 4 files;
  after SWAP removal all 2q gates are same-rung cross-chain (2400) or nearest-neighbour same-chain (9440);
  20 steps found (cut after every 120 mass rz); final perm = #250 output layout (wires 2r,2r+1 = site r).
- Test A (logs/testA_window.log): windows of 8/10 rungs (16/20 qubits), gate-by-gate dense
  statevector (unfused) vs fused symmetric MPS: max |dn| <= 5.7e-14 over all 20 steps (l300 MID 3 windows,
  l300 SCV, l12 meson).
- Test B (logs/testB_sector.log): FULL 120-qubit raw circuit with SWAPs applied physically on wires,
  exact in the 2- and 3-particle sectors (dim 7140 / 280840), vs MPS pipeline: max |dn| <= 1.4e-13 at every
  step, 5 initial states (l300 MID file x4, l12 meson x1). (First version failed its conservation assert
  because a reordering fusion split X-conjugations; replaced by contiguous runs, no reordering.)
- Slip, for the record: one Bash call accidentally fed a few shell lines to the HOST python via a stray
  heredoc; it failed with a syntax error and executed nothing. All computations ran in Docker.
- ~16:48 l300 SCV/MID chi=256 done (results/l300_*_chi256.log, 35 s each): n_f(20)=0.1111, disc 0.35 (not converged).
- 16:50 user paused everything; ladder containers stopped; partial logs moved to results/stopped_1650/.
- 18:28 resume. Added per-step checkpoint (--ckpt, ckpt/*.pkl) and code/ladder.sh (restartable).
- 18:30 started l12 SCV and meson ladders chi 256, 512, 1024 (4 CPUs each).
- 18:32 containers destroyed without a kill/die event (Docker daemon restart, it seems); resumed 18:35 from
  checkpoints with --restart unless-stopped (ladder.sh sleeps when finished; containers stopped by hand).
- l12 so far: chi=256 n_f(20)=0.3271 (disc 0.74), chi=512 0.1461 (disc 0.16); steps 1-11 match #255 PT2
  per-step values to <=2e-6 (n_f) already at chi=256. l_i=12 is MORE entangled than l_i=300 (Smax 5.48 vs 4.67).
- 18:39 lsh.py: kept SVD/QR blocks are now contiguous copies (views kept the whole U/Vh alive -> memory).
  Re-ran one Test A window + one Test B state after the change: 5.5e-14 / 7.6e-14 (logs/retest_after_contiguous.log).
  run_mps.py now also records meas_time and bond_dims per step.
- Background launcher: when l12 chi=1024 is done -> l300 SCV/MID ladders 512, 1024, 2048 (4 CPUs, 2.2 GB each).

- 18:45 l12 chi=1024 done (444 s/circuit, 4 threads, peak RSS ~0.6 GB): n_f(20)=0.116162 (disc 7.5e-3/5.9e-3),
  N_SCV=57.379204, N_MID=57.495366. Refs: 0.1169+-0.0007 (#250 PP), 0.11686+-1e-4 (#255 PT2), 0.1176+-0.0022
  (#227 MPS). Per step vs #255: |d n_f| <= 7e-7 for steps 1-15 (where own disc < 1e-6); deviation grows with
  discarded weight (step 19: -5.7e-4, step 20: -7.0e-4). Sequence chi 256/512/1024: 0.3271/0.1461/0.1162 -> not
  yet converged at 1e-3; chi=2048 needed for a firm l12 statement.
- 18:45 l300 ladder started (containers v257_l300_SCV / v257_l300_MID).

- 18:54 l300 results (n_f(20), N_SCV, N_MID, mean disc): chi=256 0.111139 57.900319 58.011458 (0.33);
  chi=512 0.115809 57.886589 58.002398 (0.072); chi=1024 0.103281 57.965419 58.068700 (2.6e-3), 405 s/circuit,
  peak RSS 0.40 GB. Per step vs #257 (chi 2800): <=5e-6 through step 16, then deviation grows with own disc.
- 18:55 l300 chi=2048 running; launcher (bg) then starts l12 chi=2048.

- 19:15 l300 chi=2048 done (1203 s/circuit, 4 threads, peak RSS 1.23 GB): N_SCV=57.967164, N_MID=58.070023,
  n_f(20)=0.102859, disc 1.24e-4/1.00e-4. Change 1024->2048: d n_f=-4.2e-4, d N_SCV=+1.7e-3, d N_MID=+1.3e-3.
- 19:15 l12 chi=2048 started (containers v257_l12b_*); after it, background launcher runs l300 chi=3072
  SCV then MID in one container (8 CPUs, 4.4 GB), code/seq.sh.

- 19:32 l12 chi=2048 done (~1010 s/circuit, peak 1.21 GB): n_f(20)=0.116871 (d from chi 1024: +7.1e-4),
  N_SCV=57.376122 (stag_SCV=-2.623878), N_MID=57.492993, disc 1.8e-4/1.4e-4. vs #255 0.11686+-1e-4 (diff
  +8e-6), stag_SCV -2.62385+-2e-4 (diff -2.8e-5); #250 0.1169+-0.0007; #227 0.1176+-0.0022. Per step vs #255:
  |d n_f| <= 6e-6 at all 20 steps. -> l12 validation passed. Summary: logs/summary_1935.txt.
- 19:32 l300 chi=3072 started: SCV then MID in container v257_l300_3072 (8 CPUs, 4.4 GB).

- 19:34 extra check (logs/verify_obs.log): n_f(20) recomputed with the #250 wire-based definition through the
  final SWAP permutation equals the site-based value to 3e-15 (l300 and l12, chi=2048).
- 19:51 l300 SCV chi=3072 at step 14; steps slowing (150 -> 238 s; another job r256val shares the CPU);
  peak RSS 2.5 GB, MPS 1.34 GB.

- 19:52 user paused everything; Docker Desktop stopped (container exit 137).
- 22:38 resume (new cap ~01:38). 22:39 restarted v257_l300_3072 (seq.sh): SCV chi=3072 resumed from the
  step-14 checkpoint (ckpt/l300_SCV_chi3072.pkl, written 19:51), then MID chi=3072.

- 22:57 l300 SCV chi=3072 done: N_SCV(20)=57.967194, disc 3.5e-6 (chi 2048 -> 3072: +3.0e-5); compute 2265 s
  total (8 threads, shared machine), peak RSS 2.6 GB, MPS 1.36 GB. MID chi=3072 started.

- 6/11 23:20 user paused everything (MID chi=3072 had reached step 14, checkpoint 23:17); Docker stopped.
- 7/11 09:34 resume (cap ~11:30): MID chi=3072 resumed from the step-14 checkpoint.

- 7/11 10:02 l300 MID chi=3072 done: N_MID(20)=58.070046, disc 2.9e-6. => chi=3072: n_f(20)=0.102852,
  N_SCV=57.967194, N_MID=58.070046. Change 2048->3072: d n_f=-7e-6, d N_SCV=+3.0e-5, d N_MID=+2.3e-5.
  Wire-based recheck at chi=3072 identical (1e-16). Final table: logs/summary_final.txt. Container stopped.
  (chi=3072 wall times 2255/2879 s are on a shared, partly throttled laptop, resumed across two pauses.)

## Result
- l_i=300: n_f(20) = 0.10285 +- ~2e-5 (converged in chi; error = a few x the 2048->3072 change).
  Supports #257 (0.10285). PP #254 (0.0798) is 0.023 lower, ~3000x the last chi change. QPU 3 Aug 0.0995 is
  0.0034 lower (its dQ is 0.2-0.27). ckpt/ (~7 GB pickles) is scratch, not for publication.
