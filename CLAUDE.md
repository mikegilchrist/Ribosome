# Ribosome (R package behind the polysomeModel manuscript)

Handoff Parts A-D are done and merged upstream (rurquidi/Ribosome#1,
2026-10-14); follow-up rurquidi/Ribosome#2 restores plotting helpers (open).
Open science questions: `Claude.Code_2026-10-09_v1.0_review.notes_figures.md`
-- never change a figure from it without MG and Ricardo.

- Remotes: `origin` = mikegilchrist/Ribosome (push here), `upstream` =
  rurquidi/Ribosome (Ricardo's; PRs only with MG's say-so). Never force-push.
- Real figure source is `Text/Paper final figures.Rmd` (sources only
  `R/Polysome_functions_RAUC.R`). Parameter grid is `expand.grid` with kappa
  fastest: each 5000-row block is one mu; select with `mu_block_rows(k)`.
- Input data lives outside the repo (`../Data`, `../../Data`); not pushed,
  ask Ricardo.
- Tests: `Rscript -e 'testthat::test_dir("tests/testthat")'` from repo root
  (limSolve in ~/R/lib/4.6 on obed). Oracle:
  `Claude.Code_2026-10-09_reference.solver_check.share.R`.
- R files, the Rmd and DESCRIPTION are CRLF: preserve it (Python: newline='').
- R CMD build in the sandbox chokes on device placeholders (.vscode, .idea):
  build from a copy of `git ls-files`.
- `retired/` = superseded code awaiting review/removal; `inst/scripts/` =
  data scripts with no replacement. Neither is used by the Rmd.
- Before moving code to `retired/`, check callers across the whole repo
  (including other retired scripts), not just the figure Rmd. Never delete
  outright; retire it. Ricardo's plotting helpers (`plot.*`, `s2_approx`) live
  in `R/Polysome_functions_RAUC.R`.
- `gh`, `git push` and `.git/config` writes fail in the sandbox (token
  unreadable; 401 / "could not lock config"): run them with the sandbox lifted.
  Auto mode may block `gh pr merge`; MG runs that one himself.
- Branches: local `master` tracks `origin/master`; upstream `dev` exists but
  Ricardo merges to `master`. Finished branches are renamed `retired/<name>`,
  never deleted. Upstream PR descriptions are shown to MG before sending.
- Do not read the 1800-line Rmd or retired/Manuscript_Figures.R whole; grep
  or delegate to a haiku subagent.
