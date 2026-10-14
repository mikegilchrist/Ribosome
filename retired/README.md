# retired: superseded code, kept for review before removal

Moved here on 2026-10-09 (Claude Code session ef5348d2, MG's decision).
Everything in this folder has a newer counterpart in the package or the
figure Rmd. It is kept in git so it can be reviewed, and removing the
folder is a future todo. It is excluded from the built package
(`.Rbuildignore`) and nothing in `R/` or `Text/` uses it.

| File | Superseded by | Evidence |
|---|---|---|
| `Manuscript_Figures.R` | `Text/Paper final figures.Rmd` | earlier version of the Rmd; its 113 outputs not in the Rmd are draft plots, none in the paper (code map, Section 1) |
| `Model_Figures.R` | the Rmd and `R/Polysome_functions_RAUC.R` | older exploratory model figures; local copies of `Ribo_mean`, `Ribo_percent` |
| `samples.R` | `CalcUnmarkedClass`, `CalcMarkedClass`, `Classes` | its `*NP` solvers equal the main solver with tau = 9 * imax * tau_NP (exact; `tests/testthat/test-equivalence.R`) |
| `legacy_helpers.R` | nothing needed it | functions cut verbatim from `R/Polysome_functions_RAUC.R`: `CHECK_STEADY_STATE` (broken, KNOWN_FAILURES #4, #5; replaced by the residual test), `plot.series`, `plot.bar`, `plot.bars`, `plot.bars.mat`, `series.mat.plot`, `s2_approx` (only callers were the two scripts above) |

The scripts call helpers that now live in `legacy_helpers.R`; to run one,
source `R/Polysome_functions_RAUC.R` and `legacy_helpers.R` first. They
were never runnable top to bottom in a clean session anyway (code map,
Section 3).

Tests that exercise retired code (`*NP` equivalence, `CHECK_STEADY_STATE`
findings) read it from here when run from the repo root and skip, with
the reason stated, under `R CMD check`.

`Model_Figures.R` was `Model Figures.R` until Part C.5c (spaces removed);
the code map uses the old names and `R/` paths.
