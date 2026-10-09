

# Model Solutions ---------------------------------------------------------

#' Build a parameter grid
#'
#' Crosses the supplied parameter values with [expand.grid()], so the first
#' column varies fastest. In the manuscript figures kappa is the first
#' column, so each block of `length(kappa)` consecutive rows has one value of
#' the next parameter (mu).
#'
#' @param parameter.range A list of numeric vectors. If it has length 1 it is
#'   replicated (see `fixed`); otherwise it must have one element per entry of
#'   `par.names` (the function stops if the lengths differ). A non-list or an
#'   empty list stops with an error.
#' @param fixed Indices (or names) of list elements to set to `fixed.val`
#'   after replication; only used when `parameter.range` has length 1. If
#'   `NA` (default) the single element is simply repeated once, giving two
#'   columns.
#' @param fixed.val Value assigned to the elements chosen by `fixed`; default 1.
#' @param par.names Character vector of column names. If given with a length-1
#'   `parameter.range`, `num.pars` is taken as `length(par.names)`. Default
#'   `NA` leaves the columns unnamed.
#' @param num.pars Number of copies of a length-1 `parameter.range` to make
#'   when `fixed` is given and `par.names` is `NA`; default 1.
#' @return A data frame with one row per combination of parameter values and
#'   one column per parameter, named by `par.names` when supplied.
#' @examples
#' CreateParMatrix(list(kappa = c(0.1, 0.2), mu = c(0.01, 0.02), lambda = 1),
#'                 par.names = c("kappa", "mu", "lambda"))
#' @export
CreateParMatrix <- function(parameter.range,fixed = NA,fixed.val=1,par.names=NA,num.pars = 1){
  if(typeof(parameter.range)!="list"){
    stop("Parameter range provided is not a list.")
  }
  if(length(parameter.range)==0){
    stop("Paramter range provided is empty.")
  }

  if(length(parameter.range)==1){
    if(all(!is.na(par.names))){
      num.pars <- length(par.names)
    }
    if(all(!is.na(fixed))){
      parameter.range <- c(parameter.range,rep(parameter.range,num.pars-1))
      parameter.range[fixed] <- fixed.val
    }else{
      parameter.range <- c(parameter.range,rep(parameter.range,num.pars=1))
    }
  }else if(length(parameter.range)!=length(par.names)){
    stop("To specify different values for each parameter the number of items in the
         list must be equal to the number of parameters in par.names.")
  }

  par.range <- parameter.range
  par.mat <- expand.grid(par.range)
  if(all(!is.na(par.names))){
    colnames(par.mat) <- par.names
  }
  return(par.mat)
}



#Nondimensionalized by tau
#' Capped (unmarked) steady-state distribution over polysome classes
#'
#' Solves for the steady-state amount of capped mRNA in each polysome class
#' 0 to `imax`. The balance equations form a tridiagonal linear system,
#' solved with [limSolve::Solve.tridiag()]. Capped mRNA enters class 0 at
#' rate lambda, moves from class `i` to `i + 1` by initiation, from `i` to
#' `i - 1` by termination, and leaves each class by decapping at rate mu.
#' Initiation out of class `imax` is 0 and termination out of class 0 is 0.
#'
#' With the default `which.interfere`, termination from class `i` is
#' `tau * i / (9 * imax)` and initiation from class `i` is
#' `kappa * (1 - i / imax)`.
#'
#' @param parameters Numeric vector of model parameters in the order
#'   `kappa, mu, lambda, tau, delta, s2`: kappa is the initiation rate, mu the
#'   decapping rate, lambda the rate at which capped mRNA is produced, tau the
#'   elongation rate in codons per second, delta the degradation rate of
#'   ribosome-free decapped mRNA, and s2 is unused. If the vector has no names,
#'   `par.names` is assigned as its names; the code looks values up by name.
#' @param imax Integer, the maximum number of ribosomes (the highest polysome
#'   class). Classes run from 0 to `imax`; class `i` is at index `i + 1` of
#'   the result.
#' @param rescale Logical. If `TRUE`, the rates are divided by
#'   `parameters[rescale.val]` before solving (the default `FALSE` leaves them
#'   unscaled).
#' @param rescale.val Name of the parameter used to rescale when
#'   `rescale = TRUE`; default `"tau"`.
#' @param par.names Names assigned to `parameters` when it is unnamed.
#' @param interference Logical. If `TRUE` (default), the rates named in
#'   `which.interfere` depend on the polysome class; if `FALSE`, initiation is
#'   `kappa` from classes `0` to `imax - 1` and termination from class `i >= 1`
#'   is `tau / (9 * imax)`, both independent of `i`.
#' @param which.interfere Character vector naming the rates that depend on the
#'   class when `interference = TRUE`; any of `"kappa"`, `"tau"`, `"mu"`.
#'   Default `c("kappa", "tau")`. With `"kappa"`, initiation from class `i`
#'   into class `i + 1` is `kappa * (1 - i / imax)`; with `"tau"`, termination
#'   from class `i` is `tau * i / (9 * imax)`; with `"mu"`, the decapping rate
#'   from class `i` is `mu * (1 - i / imax)`.
#' @return A matrix with `imax + 1` rows and one column; row `i + 1` is the
#'   amount in polysome class `i`.
#' @examples
#' p <- c(kappa = 0.1, mu = 0.01, lambda = 1, tau = 1, delta = 1e5, s2 = 0)
#' CalcUnmarkedClass(p, imax = 4)
#' @seealso [CalcMarkedClass()], [Classes()]
#' @export
CalcUnmarkedClass <- function(parameters,imax,rescale = F ,rescale.val = "tau",
                              par.names=c("kappa","mu","lambda","tau","delta","s2"),
                              interference = T, which.interfere = c("kappa","tau")){
  if(is.null(names(parameters))){
    names(parameters) <- par.names
  }
  #load parameters
  tau <- parameters["tau"]
  tau <- c(0,rep(tau/(imax*9),imax))
  kappa <- c(rep(parameters["kappa"],imax),0)
  lambda <-  parameters["lambda"]
  mu <- rep(parameters["mu"],imax+1)
  #add interference
  if(interference==T){
    if(any(which.interfere=="kappa")){
      kappa <- rep(parameters["kappa"],imax+1)*seq(1,0,-1/imax)
    }
    if(any(which.interfere=="tau")){
      tau <- rep(parameters["tau"]/(imax*9),imax+1)*seq(0,imax,1)
    }
    if(any(which.interfere=="mu")){
      mu <- parameters["mu"]*seq(1,0,-1/imax)
    }
  }
  if(rescale==T){
    tau <- tau/parameters[rescale.val]
    kappa <- kappa/parameters[rescale.val]
    lambda <- lambda/parameters[rescale.val]
    mu <- mu/parameters[rescale.val]
  }
#  print(tau)
#  print(kappa)
  upper <- tau[-1]
  lower <- kappa[-(imax+1)]
  main <- -1*(kappa+tau+mu)
  B <-  c(-lambda,rep(0,imax))
  Inv <- limSolve::Solve.tridiag(diam1 = lower,dia = main,diap1 = upper,B = B)

  return(Inv)
}

#' Decapped (marked) steady-state distribution over polysome classes
#'
#' Computes the steady-state amount of decapped mRNA in each polysome class
#' from the capped distribution. For class `i >= 1` the value is
#' `mu_i * sum(unmarked[j], j >= i) / tau_i`, where `tau_i` is the termination
#' rate from class `i` (`tau * i / (9 * imax)` with interference in tau); for
#' class 0 it is `mu_0 * sum(unmarked) / delta`. Here `mu_i` equals `mu`
#' unless `"mu"` is in `which.interfere`. With `rescale = TRUE`, delta is
#' divided by `parameters[rescale.val]` along with the other rates.
#'
#' @param parameters Numeric vector of model parameters in the order
#'   `kappa, mu, lambda, tau, delta, s2`: kappa is the initiation rate, mu the
#'   decapping rate, lambda the rate at which capped mRNA is produced, tau the
#'   elongation rate in codons per second, delta the degradation rate of
#'   ribosome-free decapped mRNA, and s2 is unused. If the vector has no names,
#'   `par.names` is assigned as its names; the code looks values up by name.
#' @param imax Integer, the maximum number of ribosomes (the highest polysome
#'   class). Classes run from 0 to `imax`; class `i` is at index `i + 1` of
#'   the result.
#' @param rescale Logical. If `TRUE`, the rates are divided by
#'   `parameters[rescale.val]` before solving (the default `FALSE` leaves them
#'   unscaled).
#' @param rescale.val Name of the parameter used to rescale when
#'   `rescale = TRUE`; default `"tau"`.
#' @param par.names Names assigned to `parameters` when it is unnamed.
#' @param interference Logical. If `TRUE` (default), the rates named in
#'   `which.interfere` depend on the polysome class; if `FALSE`, initiation is
#'   `kappa` from classes `0` to `imax - 1` and termination from class `i >= 1`
#'   is `tau / (9 * imax)`, both independent of `i`.
#' @param which.interfere Character vector naming the rates that depend on the
#'   class when `interference = TRUE`; any of `"kappa"`, `"tau"`, `"mu"`.
#'   Default `c("kappa", "tau")`. With `"kappa"`, initiation from class `i`
#'   into class `i + 1` is `kappa * (1 - i / imax)`; with `"tau"`, termination
#'   from class `i` is `tau * i / (9 * imax)`; with `"mu"`, the decapping rate
#'   from class `i` is `mu * (1 - i / imax)`.
#' @param unmarked Optional capped distribution, as returned by
#'   [CalcUnmarkedClass()]. If `NULL` (default) it is computed with the same
#'   arguments.
#' @return A numeric vector of length `imax + 1` (no dimensions); element
#'   `i + 1` is the amount in polysome class `i`.
#' @examples
#' p <- c(kappa = 0.1, mu = 0.01, lambda = 1, tau = 1, delta = 1e5, s2 = 0)
#' CalcMarkedClass(p, imax = 4)
#' @seealso [CalcUnmarkedClass()], [Classes()]
#' @export
CalcMarkedClass <- function(parameters,imax,unmarked=NULL,rescale = F, rescale.val = "tau",
                            par.names=c("kappa","mu","lambda","tau","delta","s2"),
                            interference = T, which.interfere = c("kappa","tau")){

  if(is.null(names(parameters))){
    names(parameters) <- par.names
  }
  if(is.null(unmarked)){
    unmarked <- CalcUnmarkedClass(parameters = parameters,imax = imax, rescale = rescale,
                                  rescale.val = rescale.val, par.names = par.names,
                                  interference = interference, which.interfere = which.interfere)
  }
  #load parameters
  tau <- parameters["tau"]
  tau <- c(0,rep(tau/(imax*9),imax))
  kappa <- c(rep(parameters["kappa"],imax),0)
  lambda <-  parameters["lambda"]
  mu <- rep(parameters["mu"],imax+1)
  delta = parameters["delta"]
  #add interference
  if(interference==T){
    if(any(which.interfere=="kappa")){
      kappa <- rep(parameters["kappa"],imax+1)*seq(1,0,-1/imax)
    }
    if(any(which.interfere=="tau")){
      tau <- rep(parameters["tau"]/(imax*9),imax+1)*seq(0,imax,1)
    }
    if(any(which.interfere=="mu")){
      mu <- parameters["mu"]*seq(1,0,-1/imax)
    }
  }
  if(rescale==T){
    tau <- tau/parameters[rescale.val]
    kappa <- kappa/parameters[rescale.val]
    lambda <- lambda/parameters[rescale.val]
    mu <- mu/parameters[rescale.val]
    delta <- delta/parameters[rescale.val]
  }

  n <- imax+1
  gamma.func = function(i) {
    if(i == 0){
      delta
    } else{
      tau[i+1]
    }
  }
  marked <- vapply(0:imax,FUN.VALUE = n,FUN = function(i){mu[i+1]*sum(unmarked[(i+1):(n)])/gamma.func(i)})
  return(marked)
}


#' Total (capped plus decapped) distribution over polysome classes
#'
#' Adds the decapped and capped steady-state distributions class by class.
#' Either may be supplied to avoid recomputation; otherwise they are
#' computed with [CalcUnmarkedClass()] and [CalcMarkedClass()].
#'
#' @param parameters Numeric vector of model parameters in the order
#'   `kappa, mu, lambda, tau, delta, s2`: kappa is the initiation rate, mu the
#'   decapping rate, lambda the rate at which capped mRNA is produced, tau the
#'   elongation rate in codons per second, delta the degradation rate of
#'   ribosome-free decapped mRNA, and s2 is unused. If the vector has no names,
#'   `par.names` is assigned as its names; the code looks values up by name.
#' @param imax Integer, the maximum number of ribosomes (the highest polysome
#'   class). Classes run from 0 to `imax`; class `i` is at index `i + 1` of
#'   the result.
#' @param rescale Logical. If `TRUE`, the rates are divided by
#'   `parameters[rescale.val]` before solving (the default `FALSE` leaves them
#'   unscaled).
#' @param rescale.val Name of the parameter used to rescale when
#'   `rescale = TRUE`; default `"tau"`.
#' @param par.names Names assigned to `parameters` when it is unnamed.
#' @param interference Logical. If `TRUE` (default), the rates named in
#'   `which.interfere` depend on the polysome class; if `FALSE`, initiation is
#'   `kappa` from classes `0` to `imax - 1` and termination from class `i >= 1`
#'   is `tau / (9 * imax)`, both independent of `i`.
#' @param which.interfere Character vector naming the rates that depend on the
#'   class when `interference = TRUE`; any of `"kappa"`, `"tau"`, `"mu"`.
#'   Default `c("kappa", "tau")`. With `"kappa"`, initiation from class `i`
#'   into class `i + 1` is `kappa * (1 - i / imax)`; with `"tau"`, termination
#'   from class `i` is `tau * i / (9 * imax)`; with `"mu"`, the decapping rate
#'   from class `i` is `mu * (1 - i / imax)`.
#' @param unmarked Optional capped distribution from [CalcUnmarkedClass()].
#' @param marked Optional decapped distribution from [CalcMarkedClass()].
#' @return A matrix with `imax + 1` rows and one column (the sum of the
#'   capped matrix and the decapped vector); row `i + 1` is polysome class `i`.
#' @examples
#' p <- c(kappa = 0.1, mu = 0.01, lambda = 1, tau = 1, delta = 1e5, s2 = 0)
#' Classes(p, imax = 4)
#' @seealso [CalcUnmarkedClass()], [CalcMarkedClass()]
#' @export
Classes <- function(parameters,imax,unmarked=NULL, marked = NULL,
                    rescale = F, rescale.val = "tau",
                    par.names=c("kappa","mu","lambda","tau","delta","s2"),
                    interference = T, which.interfere = c("kappa","tau")){
  if(is.null(names(parameters))){
    names(parameters) <- par.names
  }
  if(is.null(unmarked)){
    unmarked <- CalcUnmarkedClass(parameters = parameters,imax = imax, rescale = rescale,
                                     rescale.val = rescale.val, par.names = par.names,
                                     interference = interference, which.interfere = which.interfere)
  }
  if(is.null(marked)){
    marked <- CalcMarkedClass(parameters = parameters,imax = imax, unmarked = unmarked, rescale = rescale,
                                 rescale.val = rescale.val, par.names = par.names,
                                 interference = interference, which.interfere = which.interfere)
  }
  Classes <- (marked + unmarked)
  return(Classes)
}




# Summary Statistics ------------------------------------------------------

#' Class probabilities
#'
#' Normalises a distribution over polysome classes so it sums to one.
#'
#' @param Data A numeric vector (one distribution) or a matrix with polysome
#'   classes 0 to `imax` in rows and one distribution per column. A
#'   one-column matrix is treated as a matrix.
#' @return For a vector, a vector of the same length; for a matrix, a matrix
#'   of the same shape with each column summing to one.
#' @examples
#' Ribo_percent(c(2, 1, 1))
#' @export
Ribo_percent <- function(Data){
  if(is.null(dim(Data))&length(unlist(Data))>1){
    tmp <- Data/sum(Data)
  } else {
    tmp <- apply(Data,2,function(i) i/sum(i))
  }
  return(tmp)
}

#' Mean number of ribosomes per mRNA
#'
#' Mean of the class index (0 to `imax`) under the normalised distribution.
#'
#' @inheritParams Ribo_percent
#' @return A one-column matrix with one row per column of `Data` (1 x 1 for a
#'   vector).
#' @examples
#' Ribo_mean(c(2, 1, 1))
#' @export
Ribo_mean <- function(Data){
  if(is.null(dim(Data))&length(unlist(Data))>1){
    mean.tmp <- unlist((Data/sum(Data))%*%c(0:(length(Data)-1)))
  }else{
    Data_pct <- apply(Data,2,function(i) i/sum(i))
    mean.tmp <- unlist(t(Data_pct)%*%c(0:(nrow(Data)-1)))
  }
  return(mean.tmp)
}


#' Standard deviation of the number of ribosomes per mRNA
#'
#' Standard deviation of the class index under the normalised distribution.
#'
#' @inheritParams Ribo_percent
#' @return A numeric vector with one value per column of `Data` (a 1 x 1
#'   matrix for a vector input).
#' @examples
#' Ribo_sd(c(2, 1, 1))
#' @export
Ribo_sd <- function(Data){
  if(is.null(dim(Data))&length(unlist(Data))>1){
    mean <- Ribo_mean(Data)
    square.diff <- ((0:(length(Data)-1))-as.vector(mean))^2
    data.pct <- Data/sum(Data)
    var <- t(data.pct)%*%square.diff
  }else{
    Data_pct <- apply(Data,2,function(i) i/sum(i))
    mean <- Ribo_mean(Data)
    classes <- matrix((0:(nrow(Data)-1)),nrow = nrow(Data),ncol = ncol(Data))
    square.diff <- lapply(1:ncol(classes),function(i) (classes[,i]-mean[i])^2)
    var <- unlist(lapply(1:length(square.diff),function(i) Data_pct[,i]%*%square.diff[[i]]))
  }
  return(sqrt(var))
}

#' Modal polysome class
#'
#' Class index (0 to `imax`) with the largest probability. If several classes
#' tie for the maximum, all are returned.
#'
#' @inheritParams Ribo_percent
#' @return The modal class index: one value per column of `Data` for a matrix.
#' @examples
#' Ribo_mode(c(2, 1, 1))
#' @export
Ribo_mode <- function(Data){
  if(is.null(dim(Data))&length(unlist(Data))>1){
    pct.tmp <- Ribo_percent(Data)
    mode <- which(pct.tmp==max(pct.tmp))
  }else{
    pct.tmp <- Ribo_percent(Data)
    mode <- unlist(apply(pct.tmp,2,function(i) which(i==max(i))))
  }
  return(mode-1)
}


#' Ratio of column totals
#'
#' For each column, the sum of `Numerator` divided by the sum of
#' `Denominator`, e.g. the fraction of a total that is in one state.
#'
#' @param Numerator,Denominator Matrices with classes in rows and one
#'   distribution per column. Only `ncol(Denominator)` columns are used, so
#'   `Numerator` must have at least that many.
#' @return A numeric vector with one ratio per column.
#' @examples
#' Ribo_frac(matrix(1:4, 2), matrix(2:5, 2))
#' @export
Ribo_frac <- function(Numerator,Denominator){
  return(unlist(lapply(1:ncol(Denominator),function(i) sum(Numerator[,i])/sum(Denominator[,i]))))
}


#' Fraction-weighted load ratio
#'
#' Elementwise `Fraction * Numerator / Denominator`, for `Fraction`,
#' `Numerator` and `Denominator` given per parameter set (for example the
#' output of [Ribo_frac()] and [Ribo_mean()]). Prints `length(Denominator)`
#' as a side effect.
#'
#' @param Fraction,Numerator,Denominator Numeric vectors (or one-column
#'   matrices) indexed by parameter set; the loop runs over
#'   `length(Denominator)`.
#' @return A numeric vector of the same length as `Denominator`.
#' @examples
#' Ribo_frac_load(c(0.5, 0.25), c(2, 4), c(4, 4))
#' @export
Ribo_frac_load <- function(Fraction,Numerator,Denominator){
  print(length(Denominator))
  return(unlist(lapply(1:length(Denominator),function(i) (Fraction[i]*Numerator[i])/Denominator[i])))
}



# Utilities ---------------------------------------------------------------

#' Check for (and optionally create) an output directory
#'
#' Looks in `dir` for an entry named `new.dir`. If it exists, returns a message
#' saying so; if not and `add = TRUE`, creates it; otherwise prints that it does
#' not exist.
#'
#' @param dir Parent directory.
#' @param new.dir Name of the subdirectory; default `NA` uses today's date as
#'   `YYYYMMDD`.
#' @param add Logical; create the directory when missing.
#' @return A character message, or the result of `cat()` (`NULL`) when the
#'   directory is missing and `add = FALSE`.
#' @keywords internal
#' @noRd
Check_dir <- function(dir,new.dir=NA,add=F){
  if(is.na(new.dir)){
    new.dir=Sys.Date()
    new.dir <- gsub(new.dir,pattern = "-",replacement = "")
  }
  tmp.dir <- list.files(dir)
  if(any(tmp.dir==as.character(new.dir))){
    return("Directory already exists")
  }else if(add==T) {
    dir.create(paste(dir,as.character(new.dir),sep="/"))
    return(paste("Directory ",dir,as.character(new.dir)," created",sep=""))
  } else{
    return(cat(paste("Directory ",dir,as.character(new.dir)," does not exist."
                     ,sep=""),"Run Check_dir with new.dir=T to create this directory.",sep="\n"))
  }
}



#' Filled contour plot without a colour key
#'
#' A copy of [graphics::filled.contour()] with the colour-key code commented
#' out, so only the contour panel is drawn. Arguments follow
#' [graphics::filled.contour()]; `...` is passed to [graphics::title()].
#'
#' @param x,y Locations of grid lines (increasing).
#' @param z Matrix of values to contour.
#' @param ... Passed to [graphics::title()].
#' @keywords internal
#' @noRd
filled.contour2 <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1,
                                                                          length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
                             ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
                             levels = pretty(zlim, nlevels), nlevels = 20, color.palette = function(n) hcl.colors(n,
                                                                                                                  "YlOrRd", rev = TRUE), col = color.palette(length(levels) -
                                                                                                                                                               1), plot.title, plot.axes, key.title, key.axes, asp = NA,
                             xaxs = "i", yaxs = "i", las = 1, axes = TRUE,
                             frame.plot = axes, ...)
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0))
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las",
                                 "mfrow")))$mar
  on.exit(par(par.orig))
  # w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  # layout(matrix(c(1, 1), ncol = 2L), widths = c(1, lcm(w)))
  # par(las = las)
  # mar <- mar.orig
  # mar[4L] <- mar[2L]
  # mar[2L] <- 1
  # par(mar = mar)
  # plot.new()
  # plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i",
  #             yaxs = "i")
  # # rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  # # if (missing(key.axes)) {
  # #   if (axes)
  # #     axis(4)
  # # }
  # # else key.axes
  # # box()
  # # if (!missing(key.title))
  # #   key.title
  # mar <- mar.orig
  # mar[4L] <- 1
  # par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs,
              asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot)
    box()
  if (missing(plot.title))
    title(...)
  else plot.title
  invisible()
}

