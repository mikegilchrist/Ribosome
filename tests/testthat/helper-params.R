# helper-params.R -- shared setup for the Ribosome characterisation suite.
# Sources the solver files WITHOUT modifying them. testthat sets the working
# directory to tests/testthat, so paths are relative to that.
suppressPackageStartupMessages({
  library(limSolve)
  library(data.table)
})

.repo_R <- normalizePath(file.path("..", "..", "R"), mustWork = TRUE)
.repo   <- normalizePath(file.path("..", ".."), mustWork = TRUE)

source(file.path(.repo_R, "Polysome_functions_RAUC.R"), local = FALSE)

# Pull only named function definitions out of a file that has side effects
# when sourced whole (samples.R builds a 14641-row grid and runs solvers).
.extract_functions <- function(file, fnames, env = globalenv()) {
  exprs <- parse(file, keep.source = FALSE)
  got <- character(0)
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% fnames) {
      eval(e, envir = env)
      got <- c(got, as.character(e[[2]]))
    }
  }
  missing <- setdiff(fnames, got)
  if (length(missing)) stop("not found in ", file, ": ", paste(missing, collapse = ", "))
  invisible(got)
}
.extract_functions(file.path(.repo_R, "samples.R"),
                   c("CalcUnmarkedClassNP", "CalcMarkedClassNP", "ClassesNP"))
# The oracle file guards its demo with sys.nframe()==0; extract the function only.
.extract_functions(file.path(.repo, "Claude.Code_2026-10-09_reference.solver_check.share.R"),
                   "solve_model")

# ---- canonical parameter sets (brief A.2) -------------------------------
IMAX_SET   <- c(4, 39, 194)
IMAX_MAIN  <- 39
TAU_C      <- 1
LAMBDA     <- 1
DELTA      <- 1e5
S2         <- 1
# The ten yeast values as listed in the brief (3rd and 4th are duplicates).
MU_YEAST <- c(0.01, 0.0057762265, 0.0038508177, 0.0038508177, 0.0023104906,
              0.0016503504, 0.0008251752, 0.0004443251, 0.0003300701, 0.0002221626)
KAPPA_GRID <- 10^seq(-4, log10(0.5), length.out = 25)
INTERFERE_SETS <- list(tau = c("tau"), kappa_tau = c("kappa", "tau"))

# Grid with kappa fastest (as CreateParMatrix / the Rmd); one mu per block.
make_grid <- function(mu = MU_YEAST, kappa = KAPPA_GRID) {
  g <- expand.grid(kappa = kappa, mu = mu)
  g$lambda <- LAMBDA; g$tau <- TAU_C; g$delta <- DELTA; g$s2 <- S2
  g[, c("kappa", "mu", "lambda", "tau", "delta", "s2")]
}
par_vec <- function(row) unlist(row)[c("kappa", "mu", "lambda", "tau", "delta", "s2")]

# Run a package solver over a grid; returns matrix (imax+1) x nrow(grid).
solve_grid <- function(grid, imax, which.interfere, what = c("unmarked", "marked", "classes")) {
  what <- match.arg(what)
  fn <- switch(what, unmarked = CalcUnmarkedClass, marked = CalcMarkedClass, classes = Classes)
  sapply(seq_len(nrow(grid)), function(r)
    unname(fn(par_vec(grid[r, ]), imax, which.interfere = which.interfere)))
}

# Reference solver generalised for the c("tau") interference setting, where the
# coded kappa is constant for classes 0..imax-1 and 0 for class imax. For
# which.interfere containing "kappa" this is exactly the oracle (checked in
# test-reference.R).
ref_solve <- function(kappa, mu, tau_c, imax, lambda = 1, delta = 1e5,
                      kappa_interf = TRUE) {
  n <- imax + 1
  taui <- tau_c / (9 * imax) * (0:imax)
  kapi <- if (kappa_interf) kappa * seq(1, 0, length.out = n) else c(rep(kappa, imax), 0)
  A <- matrix(0, n, n)
  for (i in 1:n) {
    A[i, i] <- -(kapi[i] + taui[i] + mu)
    if (i < n) { A[i + 1, i] <- kapi[i]; A[i, i + 1] <- taui[i + 1] }
  }
  m <- solve(A, c(-lambda, rep(0, imax)))
  S <- rev(cumsum(rev(m)))
  mstar <- c(mu * S[1] / delta, mu * S[-1] / taui[-1])
  list(m = m, mstar = mstar, A = A, b = c(lambda, rep(0, imax)))
}
ref_for <- function(row, imax, which.interfere)
  ref_solve(row$kappa, row$mu, row$tau, imax, row$lambda, row$delta,
            kappa_interf = "kappa" %in% which.interfere)
