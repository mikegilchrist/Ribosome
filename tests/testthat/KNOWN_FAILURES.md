# KNOWN_FAILURES -- findings from Part A (2026-10-09, host obed, commit f70c00a code)

Nothing under R/, Text/, DESCRIPTION or NAMESPACE was modified. Each entry is
either a skipped test (the assertion that would fail on the original code) or
a passing test that pins the observed behaviour. Evidence: the numbers below
come from the full 10 mu x 25 kappa grid, imax in {4, 39, 194}, both
which.interfere settings, compared with the base-R oracle
(Claude.Code_2026-10-09_reference.solver_check.share.R).

## 1. _newtau set: AGREES (no failure; recorded because the brief asked)
- RESOLVED 2026-10-09, Part C.3: the _newtau set was deleted from
  R/Polysome_functions_RAUC.R on the strength of the numbers below, and its
  two tests were removed with it.
- CalcUnmarkedClass_newtau / CalcMarkedClass_newtau / Classes_newtau equal the
  first set and the oracle. Max abs diff 2.1e-10, max rel diff 1.1e-13 over
  the grid (abs diff is large only where values are ~4e3 = 1/mu).
- They are textually identical to the first set apart from print(tau) and
  print(kappa) in CalcUnmarkedClass_newtau (tested: output is non-empty).
  CalcMarkedClass_newtau falls back to CalcUnmarkedClass (tests pass either way
  because the two are identical). Safe-to-merge candidate per the numbers.

## 2. *NP set (samples.R) DISAGREES with the first set and oracle at the same inputs
- Tests skipped: test-equivalence.R "NP set equals the reference at the SAME
  input vector" and "... first set at the SAME input vector".
- Why: NP uses the rescaled convention tau_i = i, kappa/tau, mu/tau, so the
  parameter "tau" means something different (no 1/(9*imax) factor). It also
  always applies kappa interference, regardless of which.interfere.
- Evidence: max abs diff 4.4e3, max rel diff 1.4e3 (kappa+tau, imax 39);
  min over settings 3.9e3 abs / 396 rel (imax 4).
- Relationship (passing test): NP(p) == Classes(p with tau <- 9*imax*tau,
  which.interfere = c("kappa","tau")) to 1e-10 for unmarked, marked, Classes,
  imax 4/39/194 (observed difference exactly 0).
  Also NP's CalcUnmarkedClassNP/CalcMarkedClassNP index parameters by
  POSITION (parameters[1..5]), not by name, and take no which.interfere.

## 3. Ribo_sd vector branch errors
- RESOLVED 2026-10-09, Part C.5: the two names were corrected to the
  variables the branch assigns; the skip and the error-pinning test were
  replaced by value tests (sd = 1 on the toy input, vector == 1-column matrix).
- Test skipped: test-summary-stats.R "Ribo_sd: vector branch".
- Why: it references mean.tmp and square.diff.tmp, never defined there
  (it assigns `mean` and `square.diff`). Evidence: Ribo_sd(c(1,2,3,4)) ->
  Error: object 'mean.tmp' not found. A passing test pins this error message.
  The matrix branch is correct (hand-computed values match).

## 4. CHECK_STEADY_STATE uses the wrong tau convention for the first-set solver
- Test skipped: test-steady-state-check.R "CHECK_STEADY_STATE on first-set output".
- Why: it builds tau_i = i (the NP convention), not tau*i/(9*imax), so it
  reports "Failed" (c(1,1)) for correct first-set output (grid row 100: kappa 0.5, mu 0.00385;
  class-0 residuals 2.04 unmarked, 347.3 marked). It cannot be
  reused as the A.4 residual test; the suite rebuilds the tridiagonal matrix
  from the reference instead. A passing test pins the c(1,1) result.
  It passes (c(0,0)) on NP output, which matches its convention.

## 5. CHECK_STEADY_STATE only ever checks class 0
- No skipped test; a passing test demonstrates it.
- Why: the `return(...)` is inside the for loop (R/Polysome_functions_RAUC.R
  ~lines 175-227), so the function exits after i == 1. Evidence: multiplying
  unmarked class 5 by 10 still returns c(0,0), and post=TRUE shows zeros for
  every class beyond 0.

## Other observations (not failures)
- CalcUnmarkedClass returns an (imax+1) x 1 matrix (Solve.tridiag);
  CalcMarkedClass returns a plain vector; Classes returns an (imax+1) x 1 matrix.
  Ribo_mean returns a 1x1 (vector input) or ncol x 1 (matrix input) matrix.
- The brief's yeast mu list has a duplicated 3rd/4th entry; kept as given.
- Ribo_frac_load prints length(Denominator) to stdout.
