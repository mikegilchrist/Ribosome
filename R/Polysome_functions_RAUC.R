

# Model Solutions ---------------------------------------------------------

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

Ribo_percent <- function(Data){
  if(is.null(dim(Data))&length(unlist(Data))>1){
    tmp <- Data/sum(Data)
  } else {
    tmp <- apply(Data,2,function(i) i/sum(i))
  }
  return(tmp)
}

Ribo_mean <- function(Data){
  if(is.null(dim(Data))&length(unlist(Data))>1){
    mean.tmp <- unlist((Data/sum(Data))%*%c(0:(length(Data)-1)))
  }else{
    Data_pct <- apply(Data,2,function(i) i/sum(i))
    mean.tmp <- unlist(t(Data_pct)%*%c(0:(nrow(Data)-1)))
  }
  return(mean.tmp)
}


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


Ribo_frac <- function(Numerator,Denominator){
  return(unlist(lapply(1:ncol(Denominator),function(i) sum(Numerator[,i])/sum(Denominator[,i]))))
}


Ribo_frac_load <- function(Fraction,Numerator,Denominator){
  print(length(Denominator))
  return(unlist(lapply(1:length(Denominator),function(i) (Fraction[i]*Numerator[i])/Denominator[i])))
}



# Utilities ---------------------------------------------------------------

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

