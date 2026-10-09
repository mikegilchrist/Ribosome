# legacy_helpers.R -- functions retired from R/Polysome_functions_RAUC.R
# on 2026-10-09 (Claude Code session ef5348d2), moved verbatim (source lines
# 156-404 at commit 0399bd4). Their only callers were the retired scripts
# Manuscript_Figures.R and Model_Figures.R; the figure Rmd never used them.
# CHECK_STEADY_STATE is also broken (tests/testthat/KNOWN_FAILURES.md #4, #5);
# the test suite's residual test replaces it. Kept for review; see README.md.

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
