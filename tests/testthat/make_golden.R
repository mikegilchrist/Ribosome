# make_golden.R -- generate golden/*.rds characterisation files.
# RUN ONLY ON UNMODIFIED CODE AT f70c00a (the original solvers in
# R/Polysome_functions_RAUC.R). Regenerating after any change to the solvers
# defeats the purpose: the files are the net that protects the refactor.
# Usage (from tests/testthat):  Rscript make_golden.R
# One file per (imax, which.interfere): golden/imax<I>_<setname>.rds holding
# list(grid, imax, which.interfere, unmarked, marked, classes); each matrix is
# (imax+1) x nrow(grid), columns in grid order (kappa fastest, one mu per block).
source("helper-params.R")
dir.create("golden", showWarnings = FALSE)
g <- make_grid()
for (nm in names(INTERFERE_SETS)) for (im in IMAX_SET) {
  wi <- INTERFERE_SETS[[nm]]
  out <- list(grid = g, imax = im, which.interfere = wi,
              unmarked = solve_grid(g, im, wi, "unmarked"),
              marked   = solve_grid(g, im, wi, "marked"),
              classes  = solve_grid(g, im, wi, "classes"))
  saveRDS(out, file.path("golden", sprintf("imax%d_%s.rds", im, nm)), compress = "xz")
}
cat("golden files written:\n"); print(file.info(list.files("golden", full.names = TRUE))[, "size", drop = FALSE])
