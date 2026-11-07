# Final report of the re-check session (#257)

This is the final report that the re-check session sent to the session directing the work, on 7 Nov 2026. Only
the machine path was replaced by `<work>`; `ricontrollo_257/` is this `recheck/` folder and `verifica_257/` is
the main check. Its audit scripts (`code/audit_b.py`, `code/audit_e.py`, `code/check_init_final.py`) printed their
results to the session and did not write log files, so for those points this report is the record.

---

I could not reproduce n_f(20) myself, but my independent numbers agree with #257 and the other agent for every step I could trust. Their n_f(20) = 0.10285 stands as far as I can tell. Nothing was published or posted. Everything ran in Docker, and I did not touch verifica_257/, ricontrollo_256/ or any .gpu_lock.

## Phase 1: my computation (done blind, logged 10:50:21 local)

**Method.**
- Qubit-level MPS with 120 sites, one per qubit, no symmetry. SWAPs are applied physically.
- Gates are exact. I only fuse consecutive gates on the same qubit pair into 4x4 blocks. SVD truncation at chi = 256, with chi = 128 for comparison.
- I deduced the layout myself from the QASM. Each group of 4 wires is (i_2m, i_2m+1, o_2m, o_2m+1), and odd wires are occupied, so the SCV is (00,11,...). The electric block then acts on (i_r, o_r), and the mass sign matches H_M (rz(-0.03) on even r).
- At the end of the circuit, wires 2r and 2r+1 hold site r. The file omits the last SWAP layer.
- N = 60 - ½ Σ(-1)^r (Zi+Zo). Code is in `ricontrollo_257/code/` (`qmps.py`, `common.py`, `run_real.py`).

**Validation.**
- Random nearest-neighbour circuits with cx/h/rz/x/swap: dense statevector vs MPS agree to 4e-15.
- My reconstructed LSH circuit with N = 4 and 6 sites, 3 steps: dense vs MPS agree to ≤ 4e-12, and Q is conserved.
- My circuit generator matches both l=300 files gate for gate over all 20 steps, apart from the omitted last SWAP layer.
- On l=12, steps 1-6, my n_f is 1.82534, 1.34730, 0.69035, 0.02128, -0.49857, -0.75832. That matches the PP/TN table of #149 (1.8253, 1.3479 (PP GPU), ...).
- Q = 60 to about 1e-14.

**Results for l = 300** (chi = 256; "disc" is cumulative discarded weight, maximum over the two circuits):

| Step | n_f | disc |
|---|---|---|
| 1 | 3.650671 | 0 |
| 2 | 2.715205 | 2e-12 |
| 3 | 1.480150 | 4e-12 |
| 4 | 0.286461 | 6e-12 |
| 5 | -0.589960 | 9e-12 |
| 6 | -1.007995 | 1e-11 |
| 7 | -0.973485 | 7e-11 |
| 8 | -0.610328 | 3.1e-9 |
| 9 | -0.107441 | 8e-8 |
| 10 | 0.345667 | 1.3e-6 |
| 11 | 0.618244 | 1.6e-5 |
| 12 | 0.671375 | 1.3e-4 |

- N_SCV and N_MID per step are in `logs/PHASE1_RESULT.txt` and `logs/snap1_*.json`.
- **Trustworthy:** steps 1-8 for sure, steps 9-10 to about 1e-6, steps 11-12 only to about 1e-4 to 1e-3.
- I stopped at step 13 (10:52) and did not run steps 13-20. A non-symmetric qubit MPS would need chi of about 2000 or more there.
- chi = 128 vs 256 give the same n_f to ≤ 3.4e-7 up to step 9. A chi = 128 run on MID reached step 9.
- Cost: about 25 min per circuit on 8 threads to step 12, roughly 140 s per step at the end.
- I did not do l = 12 beyond step 6, at chi = 128.

**Comparison with the other agent and #257.**
- N_SCV and N_MID agree with the other agent's chi = 3072 values and with #257's table to ≤ 1e-5 at steps 1-10. Their 5-decimal rounding is 5e-6 or less.
- At steps 11-12 the gaps are 1e-5 and up to 1.7e-4, which is the size of my own discarded weight.
- n_f for steps 1-12 matches theirs to 5 decimals. The one gap that shows is step 12: mine 0.671375 vs 0.67135.
- All of this is with their grouped-rung, U(1)xU(1) method versus my plain qubit MPS, so it is a genuinely different route.

## Phase 2: audit

**a) Circuit reading and definitions: they match #254.**
- Their layout (`NOTE.md`: wire w -> site 2*(w//4)+(w%4)%2, chain (w%4)//2) is identical to mine, including the final wires 2r and 2r+1 for site r.
- Their `results/l300_*_chi*.json` have per-site occupations. I recomputed N from them with #254's n_f(r) formula (even r: n_i+n_o; odd r: 2-n_i-n_o). It equals their stored N to 1e-6 for every run.
- n_f = N_MID - N_SCV is as in #254. Q = 60 holds to 1e-14.

**b) Final table vs JSON: consistent.**
- I recomputed N_SCV, N_MID and n_f for all chi (l300: 256, 512, 1024, 2048, 3072; l12: 256-2048) from `results/*.json`.
- Every row of `logs/summary_final.txt` is reproduced, including the 20 per-step n_f at chi = 3072.
- The final discarded weight (3.52e-6 and 2.94e-6) matches the JSON, and the cumulative discarded weight is monotone.

**c) Convergence claim: the final number is justified, but "chi 256..3072" needs a caveat.**
- n_f(20) across chi is 0.1111, 0.1158, 0.1033, 0.102859, 0.102852.
- The sequence is not monotone and is not converged below chi = 1024. Only the last step (2048 to 3072, -7e-6, discarded weight 1.2e-4 down to 3e-6) shows convergence.
- A fit against discarded weight gives an intercept of 0.102846, 6e-6 from the chi = 3072 value.
- In my own runs the error in n_f is about 0.1 to 0.2 times the discarded weight, which would put chi = 3072 at less than 1e-6. So ±2e-5 is cautious.
- Per-step, chi = 2048 vs 3072 differ by 7e-6 at step 20 and 2e-6 at step 19; chi = 1024 vs 3072 differ by 4.3e-4 at step 20.
- Evidence of convergence below step 20 is solid, and my independent MPS matches theirs through step 10.

**d) The MID initial state: a real inconsistency in #254's text, not a misreading by #257.**
- From the QASM x gates, `SU2_T20_MID.qasm` has extra x gates on q57, q59, q60 and q62. The first two are doubled on top of the SCV prefix.
  - Site 29: (i,o) goes from 11 to 00.
  - Site 30: (i,o) goes from 00 to 11.
  - Two fermions are moved.
- `x_100_meson.qasm` (l=12) adds only x q59 and x q60. Site 29 becomes 10 and site 30 becomes 10, so one fermion is moved.
- Both states have Q = 60.
- #254's text says the initial states are "the same as in #149", but the files differ.
- Direct evidence: my l=300 n_f(1) = 3.650671 and the l=12 n_f(1) = 1.825336. The ratio is 2.0000.
- #254 also says the qubit layout and 120-bit initial strings are fixed by the QASM files. So the files define the benchmark, and #257 computed the right thing.
- I cannot tell whether the different state is intentional, because #149's own text is not in my folder (only the README).

**e) NOTE.md claims without support: I found none.**
- I checked these against logs:
  - the timings and discarded-weight figures;
  - the chi = 256/512/1024/2048/3072 values;
  - the cross-check of wire-based vs site-based n_f (`logs/verify_obs.log`);
  - `logs/testA_window.log` and `logs/testB_sector.log`, with max deviations of 1e-13 or less.
- "Steps 1-11 match #255 at chi = 256 to ≤ 2e-6": I only checked chi 256 vs 2048 on l12, which gives ≤ 1.5e-6. I did not check #255's values.
- Two weak points:
  - Test B covers only the 2- and 3-particle sectors, not the SCV/MID states themselves. Test A (windows) and my own run cover the rest.
  - The remark that l12 is "more entangled" compares chi = 512 with chi = 256 in one place. At equal chi, and at the high-chi values, it still holds: 5.34 vs 4.67, and 5.51 vs 4.81.

**Other observations:**
- #257's text says its value is "within 0.003" of the 3 Aug hardware value (0.0995), but the gap is 0.0033.
- #257's per-step N lists match mine and the other agent's to rounding.

## What I could not check
- n_f for steps 13-20 independently. Beyond my step 12, #257 and the other agent agree with each other (both use the same grouped-rung idea), but I have no third route.
- #255's per-step values, the l=12 validation of the other agent at chi = 2048 (I only did l=12 steps 1-6), the "g = 0 control", and #257's code (GitHub).
- Whether the l=300 MID state was meant to differ from #149's meson state.

Files are in `<work>/ricontrollo_257/`:
- `NOTE.md`
- `code/` (`qmps.py`, `common.py`, `run_real.py`, `test_mps.py`, `check_init_final.py`, audit scripts)
- `logs/PHASE1_RESULT.txt`, `logs/snap1_scv_chi256.json`, `logs/snap1_mid_chi256.json`
- `logs/*chi128*` for the chi = 128 comparison and the l=12 runs
- `logs/scv_chi256.*`, `logs/mid_chi256.*` for the full chi = 256 runs, including the step 13 of MID that ran after the snapshot (stopped there)
