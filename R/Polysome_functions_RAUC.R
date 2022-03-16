library(limSolve)
library(data.table)


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
  Inv <- Solve.tridiag(diam1 = lower,dia = main,diap1 = upper,B = B)
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




CHECK_STEADY_STATE <- function(unmarked, marked,parameters, imax=NA,post=F,verbose=T,tol=10){
  if(verbose==T){
    print("Steady State Check")
    print("Parameters:")
    print(parameters)
  }
  parameters <- unlist(parameters)
  if(is.na(imax)){
    imax <- length(unmarked)-1
  }
  tau <- rep(parameters[4],imax+1)*seq(0,imax,1)/parameters[4]
  kappa <- rep(parameters[1],imax+1)*seq(1,0,-1/imax)/parameters[4]
  mu <- rep(parameters[2]/parameters[4],imax+1)
  lambda <-  c(parameters[3]/parameters[4],rep(0,imax))
  delta = c(parameters[5]/parameters[4],rep(0,imax))
  marked_check <- numeric(imax+1)
  unmarked_check <- numeric(imax+1)
  for(i in 1:(imax+1)){
    if(i==1){
      unmarked_check[i] <- (-mu[i]-kappa[i])*unmarked[i] + tau[i+1]*unmarked[i+1] + lambda[i]
      marked_check[i] <- (mu[i])*unmarked[i] + tau[i+1]*marked[i+1] - delta[i]*marked[i]
    }else if(i<(imax+1)){
      unmarked_check[i] <- unmarked[i]*(-kappa[i]-tau[i]-mu[i]) + kappa[i-1]*unmarked[i-1]+tau[i+1]*unmarked[i+1]
      marked_check[i] <- (mu[i])*unmarked[i] + tau[i+1]*marked[i+1] - tau[i]*marked[i]
    }else {
      unmarked_check[i] <- unmarked[i]*(-tau[i]-mu[i]) + kappa[i-1]*unmarked[i-1]
      marked_check[i] <- (mu[i])*unmarked[i] - tau[i]*marked[i]
    }
    unmarked_check <- round(unmarked_check,tol)
    marked_check <- round(marked_check,tol)
    test <- numeric(2)
    if(verbose==T){
      print("Unmarked Class Check...")
      if(all(unmarked_check==0)){
        print("Passed")
        test[1] <- 0
      }else{
        print("Failed")
        test[1] <- 1
      }
      print("Marked Class Check...")
      if(all(marked_check==0)){
        print("Passed")
        test[2] <- 0
      }else{
        print("Failed")
        test[2] <- 1
      }
    } else {
      if(all(unmarked_check==0)){
        test[1] <- 0
      }else{
        test[1] <- 1
      }
      if(all(marked_check==0)){
        test[2] <- 0
      }else{
        test[2] <- 1
      }
    }
    if(post==T){
      out.tmp <- rbind(unmarked_check,marked_check)
      colnames(out.tmp) <- seq(0,imax,1)
      rownames(out.tmp) <- c("unmarked","marked")
      return(list(test,out.tmp))
    }else{
      return(test)
    }

  }
}



# Plotting Functions ------------------------------------------------------




plot.series <- function(data,series,par,main){
  plot(0:(nrow(data)-1),data[,series[1]],type="p",pch=20,col=1,ylim=c(0,max(data)),xlab="Class",ylab="",main=main)
  lines(0:(nrow(data)-1),data[,series[1]],col=1,lty=2)
  for(i in 2:length(series[-1])){
    points(0:(nrow(data)-1),data[,series[i]],col=i,pch=20)
    lines(0:(nrow(data)-1),data[,series[i]],col=i,lty=i)
  }
  leg.titles <- paste(c("k","m","l","t","d")[c(1,2,5)],"=")
  leg.text <- as.character(unlist(apply(t(apply(par[series,c(1,2,5)],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 0), sep=" "))),1,function(j) paste(j,collapse = ", ") )))
  legend("topright",legend = leg.text,pch = 20 ,col = 1:length(series))

}


plot.bar <- function(data,series,pars,main,beside=F){
  bars.tmp <- barplot(data[,series[1]],ylim=c((max(data[,series])*0.01),max(data[,series])),beside=beside)
  mtext(text = c(0:(nrow(data)-1))[seq(1,nrow(data),2)],side = 1,
        at = bars.tmp[seq(1,nrow(data),2)],line = 0.5)
  mtext(text = c(0:(nrow(data)-1))[seq(2,nrow(data),2)],side = 1,
        at = bars.tmp[seq(2,nrow(data),2)],line = 1.5)
  mtext(text="Classes",side = 1,line = 2.5,font=1,cex=1.25)
  leg.titles <- paste(c("k","m","l","t","d")[c(1,2,5)],"=")
  leg.text <- as.character(unlist(apply(t(apply(pars[series,c(1,2,5)],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 0), sep=" "))),1,function(j) paste(j,collapse = ", ") )))
  mtext(text=paste(main,leg.text,sep=", "),side=3,line=1,font=2,cex=1.25)
  # lines(0:(nrow(data)-1),data[,series[1]],col=1,lty=2)
  # for(i in 2:length(series[-1])){
  #   points(0:(nrow(data)-1),data[,series[i]],col=i,pch=20)
  #   lines(0:(nrow(data)-1),data[,series[i]],col=i,lty=i)
  # }

}


plot.bars <- function(unmarked,marked,series,pars,main,
                      yscale=0.1,xaxt.plot=T,main.plot=T,
                      log=F,legend.plot=T,pars.leg=F,
                      global.ymax=NA,beside=F,xaxt.step=1,
                      set.params=NA){
  if(is.na(global.ymax)){
    ymax <- max(unmarked[,series[1]]+marked[,series[1]])
  } else {
    ymax <- global.ymax
  }
  ymax <- ymax + ymax*yscale
  if(log==T){
    ymax <- round(ymax)
    bars.tmp <- barplot(height = rbind(unmarked[,series[1]],marked[,series[1]]),
                        ylim=c(1,ymax),col = c("grey80","grey40"),log="y",
                        beside=beside,xpd =F)
    #axis(2,at = seq(0,ymax,ymax/4),labels = c("",round(log(seq(0,ymax,ymax/4)[2:5],10),2)))
  }else{
    bars.tmp <- barplot(height = rbind(unmarked[,series[1]],marked[,series[1]]),
                        ylim=c(ymax*0.01,ymax),col = c("grey80","grey40"),
                        beside=beside,xpd = F)
  }
  if(xaxt.plot ==T){
    if(beside==F){
      mtext(text = c(0:(nrow(unmarked)-1))[seq(1,nrow(unmarked),xaxt.step)],side = 1,
            at = bars.tmp[seq(1,nrow(unmarked),xaxt.step)],line = 0.5)
      if(xaxt.step==2){
        mtext(text = c(0:(nrow(unmarked)-1))[seq(2,nrow(unmarked),xaxt.step)],side = 1,
              at = bars.tmp[seq(2,nrow(unmarked),xaxt.step)],line = 1.5)
      }
    } else {
      mtext(text = c(0:(nrow(unmarked)-1))[seq(1,nrow(unmarked),xaxt.step)],side = 1,
            at = colMeans(bars.tmp[,seq(1,nrow(unmarked),xaxt.step)]),line = 0.5)
      if(xaxt.step==2){
        mtext(text = c(0:(nrow(unmarked)-1))[seq(2,nrow(unmarked),xaxt.step)],side = 1,
              at = colMeans(bars.tmp[,seq(2,nrow(unmarked),xaxt.step)]),line = 1.5)
      }
    }
    mtext(text="Classes",side = 1,line = 2.5,font=1,cex=1.25)
  }
  if(any(is.na(set.params))){
    leg.titles <- paste(c("k","m","l","t","d")[c(1,2,5)],"=")
    leg.text <- as.character(unlist(apply(t(apply(pars[series,c(1,2,5)],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 2), sep=" "))),1,function(j) paste(j,collapse = ", ") )))
  } else {
    leg.titles <- set.params
    leg.text <- as.character(unlist(apply(t(apply(pars[series,..set.params],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 2), sep=" "))),1,function(j) paste(j,collapse = ", ") )))
  }

  if(main.plot == T ){
    mtext(text=paste(main,leg.text,sep=", "),side=3,line=1,font=2,cex=1.25)
  }
  if(legend.plot==T){
    legend("topright",legend = c("Unmarked","Marked"),fill = c("grey80","grey40"))
  }
  if(pars.leg == T){
    legend("topright",legend = leg.text,bg = rgb(1,1,1,0),box.col = rgb(1,1,1,0))
  }
  box()
}


plot.bars.mat <- function(unmarked,marked,series,pars,main,
                          yscale=0.1,xaxt.plot=T,main.plot=T,
                          log=F,legend.plot=T,pars.leg=F,beside=T,
                          xaxt.step=1,global.ymax=F,set.params=NA){
  if(global.ymax==T){
    global.ymax <- max(unmarked[,series]+marked[,series])
  } else {
    global.ymax <- NA
  }

  #par(mfrow=c(4,3))
  count <- 1

  layout.matrix <- t(matrix(c(1:12,0,0,0), nrow = 3, ncol = 5))

  layout(mat = layout.matrix,
         heights = c(2, 2, 2 ,2,1), # Heights of the two rows
         widths = c(2, 2, 2)) # Widths of the two columns


  for(i in series){
    par(mar=c(0,2,1,1))
    if(count<9){
      plot.bars(unmarked = unmarked,marked = marked,series=i,
                pars=pars,main=main,yscale=yscale,
                xaxt.plot = F,main.plot = F,legend.plot = F,
                pars.leg = T,log = log,global.ymax = global.ymax,
                beside=beside,xaxt.step = xaxt.step,set.params=set.params)
    }else {
      par(mar=c(0,2,1,1))
      plot.bars(unmarked = unmarked,marked = marked,series=i,
                pars=pars,main=main,yscale=yscale,
                xaxt.plot = T,main.plot = F,legend.plot = F,
                pars.leg = T,log = log,global.ymax = global.ymax,
                beside=beside,xaxt.step = xaxt.step,set.params=set.params)
    }
    count <- count+1
  }
  plot(1:nrow(unmarked),rep(0,nrow(unmarked)),axes=F,col=rgb(0,0,0,0),ylab="",xlab="")
  legend("center",legend = c("Unmarked","Marked"),fill = c("grey80","grey40"),cex = 1.25)
}

series.mat.plot <- function(parameters,Unmarked,Marked,group.by ="kappa",set.params = c("mu","tau"),file.loc="",log=F,beside=T,xaxt.step=5,global.ymax=F){
  if(is.character(group.by)){
    groupsize <-nrow(unique(parameters[,..group.by]))
  }else if(is.numeric(group.by)){
    groupsize <- group.by
  }

  series.mat <- t(matrix(1:nrow(parameters),nrow=groupsize))
  for(i in 1:(nrow(parameters)/groupsize)){
    print(i)
    series <- unlist(series.mat[i,])
    pars=parameters
    leg.titles <- set.params
    leg.text <- as.character(unlist(apply(t(apply(pars[series[1],..set.params],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 0), sep="_"))),1,function(j) paste(j,collapse = "_") )))
    png(paste(file.loc,"/",group.by,"series","_",leg.text,"_globalymax",global.ymax,"_log",log,"_imax",nrow(Unmarked)-1,".png",sep=""),width = 8,height = 6,res = 600,units = "in")
    plot.bars.mat(unmarked = Unmarked,marked = Marked,series = series,
                  pars=parameters,main="Sum classes. l=1, d=10^5",log=log,beside=beside,xaxt.step = xaxt.step,global.ymax = global.ymax,set.params=set.params)
    dev.off()
  }

}



s2_approx <- function(parameters,imax,mean){
  if(is.null(names(parameters))){
    names(parameters) <- par.names
  }
  #load parameters
  tau <- parameters["tau"]
  kappa <- parameters["kappa"]
  mu <- parameters["mu"]
  approx <- 0.5*(kappa+(mean/imax)*(kappa-tau))/((kappa+tau)/imax-mu)
}


# New tau formula ---------------------------------------------------------


CalcUnmarkedClass_newtau <- function(parameters,imax,rescale = F ,rescale.val = "tau",
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
  print(tau)
  print(kappa)
  upper <- tau[-1]
  lower <- kappa[-(imax+1)]
  main <- -1*(kappa+tau+mu)
  B <-  c(-lambda,rep(0,imax))
  Inv <- Solve.tridiag(diam1 = lower,dia = main,diap1 = upper,B = B)
  return(Inv)
}




CalcMarkedClass_newtau <- function(parameters,imax,unmarked=NULL,rescale = F, rescale.val = "tau",
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


Classes_newtau <- function(parameters,imax,unmarked=NULL, marked = NULL,
                    rescale = F, rescale.val = "tau",
                    par.names=c("kappa","mu","lambda","tau","delta","s2"),
                    interference = T, which.interfere = c("kappa","tau")){
  if(is.null(names(parameters))){
    names(parameters) <- par.names
  }
  if(is.null(unmarked)){
    unmarked <- CalcUnmarkedClass_newtau(parameters = parameters,imax = imax, rescale = rescale,
                                     rescale.val = rescale.val, par.names = par.names,
                                     interference = interference, which.interfere = which.interfere)
  }
  if(is.null(marked)){
    marked <- CalcMarkedClass_newtau(parameters = parameters,imax = imax, unmarked = unmarked, rescale = rescale,
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
    square.diff <- ((0:(length(Data)-1))-as.vector(mean.tmp))^2
    data.pct <- Data/sum(Data)
    var <- t(data.pct)%*%square.diff.tmp
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




