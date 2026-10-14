

# Model Solutions ---------------------------------------------------------

#' Build a parameter grid
#'
#' Crosses the parameter values with [expand.grid()], so the first column
#' varies fastest. With kappa first, each block of `length(kappa)` rows has
#' one value of the next parameter (mu).
#'
#' @param parameter.range A list with one numeric vector per entry of
#'   `par.names`, or a list of length 1 that is replicated (see `fixed`).
#'   Anything else stops with an error.
#' @param fixed Indices of the elements set to `fixed.val` after a length-1
#'   `parameter.range` is replicated to `num.pars` elements. If `NA`
#'   (default) the single element is repeated once, giving two columns.
#' @param fixed.val Value for the elements chosen by `fixed`.
#' @param par.names Column names; required when `parameter.range` has more
#'   than one element. With a length-1 `parameter.range`, its length
#'   overrides `num.pars`.
#' @param num.pars Number of elements made from a length-1
#'   `parameter.range` when `fixed` is given.
#' @return A data frame with one row per combination of values and one
#'   column per parameter.
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
#' Solves the steady state of the capped ODEs, a tridiagonal linear system,
#' with [limSolve::Solve.tridiag()]. Capped mRNA enters class 0 at rate
#' lambda, moves from class `i` to `i + 1` by initiation and from `i` to
#' `i - 1` by termination, and leaves every class by decapping at rate mu.
#' With the defaults, initiation from class `i` is `kappa * (1 - i / imax)`
#' and termination is `tau * i / (9 * imax)`.
#'
#' @param parameters Numeric vector `c(kappa, mu, lambda, tau, delta, s2)`:
#'   kappa the translation initiation rate (kappa_0 in the paper), mu the
#'   decapping rate, lambda the production rate of capped mRNA into class 0,
#'   tau the elongation rate in codons per second (the paper's tau is
#'   tau / 9), delta the removal rate of decapped mRNA from class 0, s2
#'   unused. Values are looked up by name; an unnamed vector is given
#'   `par.names`.
#' @param imax Maximum number of ribosomes (the highest polysome class).
#'   Class `i` is at index `i + 1` of the result.
#' @param rescale Logical; if `TRUE`, the rates and lambda are divided by
#'   `parameters[rescale.val]`. This leaves the result unchanged up to
#'   rounding.
#' @param rescale.val Name of the parameter used by `rescale`.
#' @param par.names Names given to an unnamed `parameters`.
#' @param interference Logical; if `FALSE`, `which.interfere` is ignored,
#'   initiation is `kappa` from classes 0 to `imax - 1` and termination is
#'   `tau / (9 * imax)` from classes 1 to `imax`.
#' @param which.interfere Rates that depend on class `i` when
#'   `interference = TRUE`: `"kappa"` makes initiation
#'   `kappa * (1 - i / imax)`, `"tau"` makes termination
#'   `tau * i / (9 * imax)`, `"mu"` makes decapping `mu * (1 - i / imax)`.
#'   Rates not named are as for `interference = FALSE`.
#' @return An `(imax + 1)` x 1 matrix; row `i + 1` is class `i`.
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
#' Computes the decapped steady state from the capped one, `m`: class
#' `i >= 1` is `mu_i * sum(m_j, j >= i) / tau_i`, where `tau_i` is the
#' termination rate from class `i`, and class 0 is
#' `mu_0 * sum(m_j) / delta`. `mu_i` is `mu` unless `"mu"` is in
#' `which.interfere`.
#'
#' @inheritParams CalcUnmarkedClass
#' @param unmarked Capped distribution from [CalcUnmarkedClass()]; if `NULL`
#'   (default) it is computed from the other arguments.
#' @return A numeric vector of length `imax + 1`; element `i + 1` is class
#'   `i`.
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
#' Class-by-class sum of [CalcUnmarkedClass()] and [CalcMarkedClass()];
#' either can be supplied to avoid recomputing it.
#'
#' @inheritParams CalcUnmarkedClass
#' @inheritParams CalcMarkedClass
#' @param marked Decapped distribution from [CalcMarkedClass()]; if `NULL`
#'   (default) it is computed.
#' @return An `(imax + 1)` x 1 matrix; row `i + 1` is class `i`.
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
#' Normalises a distribution over polysome classes to sum to one.
#'
#' @param Data A numeric vector of length > 1, or a matrix with classes 0 to
#'   `imax` in rows and one distribution per column.
#' @return An object of the same shape as `Data`, each column (or the
#'   vector) summing to one.
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
#' `sum(i * p_i)` over classes 0 to `imax`, for each distribution.
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
#' Standard deviation of the class index `i` for each distribution.
#'
#' @inheritParams Ribo_percent
#' @return A numeric vector with one value per column of `Data` (a 1 x 1
#'   matrix for a vector).
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
#' Class index (0 to `imax`) with the largest probability.
#'
#' @inheritParams Ribo_percent
#' @return The modal class index, one per column of `Data` for a matrix. Tied
#'   classes are all returned, concatenated across columns.
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
#' For each column, `sum(Numerator[, i]) / sum(Denominator[, i])`; e.g.
#' `Ribo_frac(marked, classes)` is the fraction of mRNA that is decapped.
#'
#' @param Numerator,Denominator Matrices with classes in rows and one
#'   distribution per column. Only the first `ncol(Denominator)` columns of
#'   `Numerator` are used.
#' @return A numeric vector with one ratio per column.
#' @examples
#' Ribo_frac(matrix(1:4, 2), matrix(2:5, 2))
#' @export
Ribo_frac <- function(Numerator,Denominator){
  return(unlist(lapply(1:ncol(Denominator),function(i) sum(Numerator[,i])/sum(Denominator[,i]))))
}


#' Fraction-weighted load ratio
#'
#' Elementwise `Fraction * Numerator / Denominator`. The Rmd calls
#' `Ribo_frac_load(Ribo_frac(M, F), Ribo_mean(M), Ribo_mean(F))`, which is
#' the fraction of bound ribosomes that are on the mRNA in `M`. Prints
#' `length(Denominator)`.
#'
#' @param Fraction,Numerator,Denominator Numeric vectors (or one-column
#'   matrices), one value per parameter set.
#' @return A numeric vector of the same length as `Denominator`.
#' @examples
#' Ribo_frac_load(c(0.5, 0.25), c(2, 4), c(4, 4))
#' @export
Ribo_frac_load <- function(Fraction,Numerator,Denominator){
  print(length(Denominator))
  return(unlist(lapply(1:length(Denominator),function(i) (Fraction[i]*Numerator[i])/Denominator[i])))
}



# Utilities ---------------------------------------------------------------

#' Report whether `new.dir` (default today's date, YYYYMMDD) exists in `dir`;
#' create it if `add = TRUE`.
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



#' [graphics::filled.contour()] with the colour key commented out, so only
#' the contour panel is drawn; `...` goes to [graphics::title()].
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

