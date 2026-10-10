# Code and equation notes for Ricardo, with suggested patches

Written 2026-10-10 on obed, Claude Code session ef5348d2, at MG's request.
Found during the documentation review of PR mikegilchrist/Ribosome#1. None
of these is applied: each is a suggestion for Ricardo to review and
implement (or reject). Patches are in `Claude.Code_2026-10-10_v1.0_patches/`
and were made against commit d2093ad (code) and polysomeModel 8eb28e6
(manuscript). Figure-level questions are separate, in
`Claude.Code_2026-10-09_v1.0_review.notes_figures.md`.

Apply a code patch from the repo root with `git apply <patch>`; the
manuscript patch from the polysomeModel root. All four pass
`git apply --check`.

## 1. Notes only (no patch)

1. **`rescale` does not change results.** With `rescale = TRUE`,
   `CalcUnmarkedClass` and `CalcMarkedClass` divide every rate *and*
   lambda by `parameters[rescale.val]`, so the steady state is unchanged
   up to rounding (max abs difference 3.6e-15 at tau = 7, imax = 39). If a
   nondimensionalisation was intended (e.g. reporting abundances in units
   of lambda/tau), the code does not do it. The comment
   `#Nondimensionalized by tau` above `CalcUnmarkedClass` is therefore
   misleading; consider rewording or removing it. The figures are not
   affected: they set tau = 1.

## 2. Code bugs (patches 01-03)

2. **`01_mu.interference_decapped.flux.patch` -- wrong decapped classes
   when decapping is density dependent.** With `"mu"` in
   `which.interfere`, `CalcMarkedClass` computes
   `m*_i = mu_i * sum_{j>=i} m_j / tau_i`. The decapped steady state
   requires `tau_i m*_i = sum_{j>=i} mu_j m_j` (i >= 1) and
   `delta m*_0 = sum_j mu_j m_j`: decapped class i collects the decapping
   flux of every capped class j >= i, each at its own rate mu_j. Example,
   imax = 5, kappa = 0.05, mu = 0.01, tau = 1: class 1 is 36.6 (code) vs
   29.7 (flux balance). The patch changes one line and adds
   `tests/testthat/test-mu-interference.R` (9 cases: fail before, pass
   after). With constant mu the change is algebraically a no-op; all
   existing tests still pass (FAIL 0 | SKIP 3 | PASS 7660). No figure
   uses `"mu"`.
3. **`02_check.dir_messages.patch` -- `Check_dir` messages.** The hint
   says `new.dir=T` but the argument that creates the directory is
   `add=T`; the "created" and "does not exist" messages omit the `/`
   between `dir` and `new.dir`. Messages only.
4. **`03_create.par.matrix_num.pars.patch` -- `CreateParMatrix`
   single-range branch.** `rep(parameter.range, num.pars=1)` passes
   `num.pars` as an argument `rep` ignores; the `fixed` branch two lines
   up uses `num.pars-1`, which is presumably what was meant. Effect for a
   length-1 `parameter.range` with no `fixed`: with `par.names` of length
   6 it currently errors and with the patch returns the full grid
   (64 x 6 for `list(1:2)`); with no `par.names` the result goes from
   2 columns to 1. The Rmd always passes a 6-element list, so this branch
   is not used there. Also fixes the "Paramter" typo in a `stop()`
   message. All tests pass.

## 3. Manuscript equations (patch 04)

5. **`04_manuscript_equation.typos.patch`** (for
   `polysomeModel/mRNA_Populations_PLOS_GENETICS.tex`). Checked against
   the paper's own ODEs (eq:kappa_i, eq:Capped_ODE, eq:Decapped_ODE):
   - eq:Capped_ODE, m_i line: initiation from class i-1 is written
     `\kappa0 ((i-1)/imax)`; eq:kappa_i gives `\kappa_0 (1 - (i-1)/imax)`.
   - eq:Decapped_ODE, last line: `\mu m*_imax - \tau m_imax` has m and m*
     swapped; should be `\mu m_imax - \tau m*_imax`.
   - eq:decapped_abundance: with termination `\tau i/imax`, the balance
     gives `m*_i = (imax/i)(mu/tau) sum_{j>=i} m_j`; the m*_1, m*_i and
     m*_imax lines drop or misplace the imax factor.
   For reference: the code's termination rate is `tau_c * i / (9 imax)`,
   so the paper's tau in these equations is the code's tau / 9; the paper
   already explains the tau_c conversion (tex L304-310).
