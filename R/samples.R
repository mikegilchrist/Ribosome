library(limSolve)
library(data.table)

CalcUnmarkedClassNP <- function(parameters,imax,rescale.fn = NULL,recscale.val = NULL){
  tau <- rep(parameters[4],imax+1)*seq(0,imax,1)/parameters[4]
  kappa <- rep(parameters[1],imax+1)*seq(1,0,-1/imax)/parameters[4]
  lambda <-  parameters[3]/parameters[4]
  mu <- parameters[2]/parameters[4]
  upper <- tau[-1]
  lower <- kappa[-(imax+1)]
  main <- -1*(kappa+tau+mu)
  B <-  c(-lambda,rep(0,imax))
  Inv <- Solve.tridiag(diam1 = lower,dia = main,diap1 = upper,B = B)
  return(Inv)
}



CalcMarkedClassNP <- function(parameters,imax,unmarked=NULL){
  if(is.null(unmarked)){
    unmarked <- CalcUnmarkedClassNP(parameters,imax)
  }
  
  kappa = parameters[1]/parameters[4]
  mu = parameters[2]/parameters[4]
  lambda = parameters[3]/parameters[4]
  delta = parameters[5]/parameters[4]
  tau = parameters[4]/parameters[4]
  
  n <- imax+1
  gamma.func = function(i) {
    if(i == 0) delta/tau else i
  }
  marked <- mu/tau*vapply(0:imax,FUN.VALUE = n,FUN = function(i){sum(unmarked[(i+1):(n)])/gamma.func(i)})
  return(marked) 
}


ClassesNP <- function(parameters,imax,unmarked=NULL,marked=NULL){
  if(is.null(unmarked)){
    unmarked <- CalcUnmarkedClassNP(parameters = parameters,imax=imax)    
  }
  if(is.null(marked)){
    marked <- CalcMarkedClassNP(parameters = parameters,imax = imax,unmarked = unmarked)    
  }
  Classes <- (marked + unmarked)
  return(Classes)
}


par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.range <- c(10^-5,10^-4,10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3,10^4,10^5)
par.mat <- expand.grid(par.range,par.range,par.range,par.range)
par.mat[,5:6] <- 1
par.mat[,4:5] <- par.mat[c(5,4)]
par.names <- c("kappa","mu","lambda","tau","delta","s2")
colnames(par.mat) <- par.names
par.mat.dt <- as.data.table(par.mat)
par.mat.dt[,"index":=1:nrow(par.mat.dt)]
par.mat.lambda1 <- par.mat.dt[lambda==1,]
par.mat.lambda1[,"index2":=(1:nrow(par.mat.lambda1))]

Unmarked_imax30 <- apply(par.mat.lambda1,1,function(i) CalcUnmarkedClassNP(parameters = i,imax = 30))
Marked_imax30 <- apply(par.mat.lambda1,1,function(i) CalcMarkedClassNP(parameters = i,imax = 30))
Classes_imax30 <- Unmarked_imax30+Marked_imax30
ncol(Classes_imax30)

Unmarked_imax30_I <- apply(par.mat.lambda1,1,function(i) CalcUnmarkedClassNP_I(parameters = i,imax = 30))
Marked_imax30_I <- apply(par.mat.lambda1,1,function(i) CalcMarkedClassNP(parameters = i,imax = 30,unmarked = Unmarked_imax30_I))
Classes_imax30_I <- Unmarked_imax30_I+Marked_imax30_I
ncol(Classes_imax30_I)

Unmarked_imax30_NI <- apply(par.mat.lambda1,1,function(i) CalcUnmarkedClassNP_NI(parameters = i,imax = 30))
Marked_imax30_NI <- apply(par.mat.lambda1,1,function(i) CalcMarkedClassNP(parameters = i,imax = 30,unmarked = Unmarked_imax30_NI))
Classes_imax30_NI <- Unmarked_imax30_NI+Marked_imax30_NI
ncol(Classes_imax30_NI)


max(Classes_imax30)

par.mat.lambda1

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

plot.series(Classes_imax30,1:5,par=par.mat.lambda1,main="Sum classes. l=1, t =1")
plot.series(Classes_imax30,6:11,par=par.mat.lambda1,main="Sum classes. l=1, t =1")


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
plot.bar(Classes_imax30,8,par=par.mat.lambda1,main="Sum classes. l=1, t =1")


plot.bars <- function(unmarked,marked,series,pars,main,
                      yscale=0.1,xaxt.plot=T,main.plot=T,
                      log=F,legend.plot=T,pars.leg=F,
                      global.ymax=NA,beside=F,xaxt.step=1){
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
  leg.titles <- paste(c("k","m","l","t","d")[c(1,2,5)],"=")
  leg.text <- as.character(unlist(apply(t(apply(pars[series,c(1,2,5)],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 0), sep=" "))),1,function(j) paste(j,collapse = ", ") )))
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


plot.bars(unmarked = Unmarked_imax30,marked = Marked_imax30,series = 789,pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,xaxt.plot = T,xaxt.step = 2)
plot.bars(Unmarked_imax30,Marked_imax30,series = 1221,pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,beside = T,xaxt.plot = T,xaxt.step = 2)
plot.bars(Unmarked_imax30_I,Marked_imax30_I,series = 1221,pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,beside = T,xaxt.plot = T,xaxt.step = 2)
plot.bars(Unmarked_imax30_NI,Marked_imax30_NI,series = 1218,pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,beside = T,xaxt.plot = T,xaxt.step = 2)


plot.bars.mat <- function(unmarked,marked,series,pars,main,
                          yscale=0.1,xaxt.plot=T,main.plot=T,
                          log=F,legend.plot=T,pars.leg=F,beside=T,xaxt.step=1,global.ymax=F){
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
                beside=beside,xaxt.step = xaxt.step)
    }else {
      par(mar=c(0,2,1,1))
      plot.bars(unmarked = unmarked,marked = marked,series=i,
                pars=pars,main=main,yscale=yscale,
                xaxt.plot = T,main.plot = F,legend.plot = F,
                pars.leg = T,log = log,global.ymax = global.ymax,
                beside=beside,xaxt.step = xaxt.step)
    }
    count <- count+1
  }
  plot(1:nrow(unmarked),rep(0,nrow(unmarked)),axes=F,col=rgb(0,0,0,0),ylab="",xlab="")
  legend("center",legend = c("Unmarked","Marked"),fill = c("grey80","grey40"),cex = 1.25)
}

plot.bars.mat(unmarked = Unmarked_imax30,marked = Marked_imax30,series = 1211:1221,
              pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,beside = F,xaxt.step = 5)

  
  
  
series <- 1247:1257

png("../../../Results/20210928/test_model_example.png",width = 8,height = 6,res = 600,units = "in")
plot.bars.mat(unmarked = Unmarked_imax30,marked = Marked_imax30,series = 1211:1221,
              pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,beside = T,xaxt.step = 5)
dev.off()


plot.bars.mat(unmarked = Unmarked_imax30,marked = Marked_imax30,series = 1321:1331,
              pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,beside = T,xaxt.step = 5)


series.mat <- t(matrix(1:nrow(par.mat.lambda1),nrow=11))
for(i in 1:(nrow(par.mat.lambda1)/11)){
  print(i)
  series <- unlist(series.mat[i,])
  pars=par.mat.lambda1
  leg.titles <- c("k","m","l","t","d")[c(2,5)]
  leg.text <- as.character(unlist(apply(t(apply(pars[series[1],c(2,5)],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 0), sep="_"))),1,function(j) paste(j,collapse = "_") )))
  png(paste("../../../Results/20211112/kseries_",leg.text,"_unlog_NOGLOBAL.png",sep=""),width = 8,height = 6,res = 600,units = "in")
  plot.bars.mat(unmarked = Unmarked_imax30,marked = Marked_imax30,series = series,
                pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,beside=T,xaxt.step = 5,global.ymax = F)
  dev.off()
}

for(i in 1:(nrow(par.mat.lambda1)/11)){
  print(i)
  series <- unlist(series.mat[i,])
  pars=par.mat.lambda1
  leg.titles <- c("k","m","l","t","d")[c(2,5)]
  leg.text <- as.character(unlist(apply(t(apply(pars[series[1],c(2,5)],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 0), sep="_"))),1,function(j) paste(j,collapse = "_") )))
  png(paste("../../../Results/20211112/kseries_",leg.text,"_I_unlog_NOGLOBAL.png",sep=""),width = 8,height = 6,res = 600,units = "in")
  plot.bars.mat(unmarked = Unmarked_imax30_I,marked = Marked_imax30_I,series = series,
                pars=par.mat.lambda1,main="Sum classes I. l=1, t =1",log=F,beside=T,xaxt.step = 5,global.ymax = F)
  dev.off()
}

for(i in 1:(nrow(par.mat.lambda1)/11)){
  print(i)
  series <- unlist(series.mat[i,])
  pars=par.mat.lambda1
  leg.titles <- c("k","m","l","t","d")[c(2,5)]
  leg.text <- as.character(unlist(apply(t(apply(pars[series[1],c(2,5)],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 0), sep="_"))),1,function(j) paste(j,collapse = "_") )))
  png(paste("../../../Results/20211112/kseries_",leg.text,"_NI_unlog_NOGLOBAL.png",sep=""),width = 8,height = 6,res = 600,units = "in")
  plot.bars.mat(unmarked = Unmarked_imax30_NI,marked = Marked_imax30_NI,series = series,
                pars=par.mat.lambda1,main="Sum classes NI_K. l=1, t =1",log=F,beside=T,xaxt.step = 5,global.ymax = F)
  dev.off()
}




U_M_ratio <- Unmarked_imax30/Marked_imax30

t(U_M_ratio)[1:10,]*100

U_percent <- Unmarked_imax30/(Marked_imax30+Unmarked_imax30)*100

M_percent <- Marked_imax30/(Marked_imax30+Unmarked_imax30)*100

t(U_percent)[1:10,]

U_sum <- colSums(Unmarked_imax30)
M_sum <- colSums(Marked_imax30)
#unmarked
plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),U_sum)
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,3]),10),U_sum)
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),U_sum)
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),U_sum)
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,2]),10),U_sum)

plot(log(unique(U_sum),10))

plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),log(U_sum,10))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,3]),10),log(U_sum,10))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log(U_sum,10))
abline(a=c(0,-1),col="red")
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),log(U_sum,10),col=rep(rainbow(11),121))
lines(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),log(U_sum,10))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,2]),10),log(U_sum,10))

plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),log(U_sum,10),col=rep(rainbow(11),121))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log(U_sum,10),col=rep(rainbow(11),121))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,4]),10),log(U_sum,10),col=rep(rainbow(11),121))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,5]),10),log(U_sum,10),col=rep(rainbow(11),121))


#Marked
plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),M_sum)
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,3]),10),M_sum)
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),M_sum)
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),M_sum)
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,2]),10),M_sum)

plot(log(unique(M_sum),10))

plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),log(M_sum,10))
plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,2]),10),log(M_sum,10))
plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,1]),10),log(M_sum,10))

plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),log(Marked_imax30[1,],10))
plot(log(unlist(par.mat.lambda1[,3]/par.mat.lambda1[,5]),10),log(Marked_imax30[1,],10))
plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,2]),10),log(Marked_imax30[1,],10))
plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,1]),10),log(Marked_imax30[1,],10))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log(Marked_imax30[1,],10))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,3]),10),log(Marked_imax30[1,],10))

plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),log(colSums(Marked_imax30[-1,]),10))
plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,2]),10),log(colSums(Marked_imax30[-1,]),10))
plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,1]),10),log(colSums(Marked_imax30[-1,]),10))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log(colSums(Marked_imax30[-1,]),10))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,3]),10),log(colSums(Marked_imax30[-1,]),10))


plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,3]),10),log(M_sum,10))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log(M_sum,10))
abline(a=c(0,-1),col="red")
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),log(M_sum,10),col=rep(rainbow(11),121))
lines(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),log(M_sum,10))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,2]),10),log(M_sum,10))
#combined

plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),(U_sum+M_sum))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,3]),10),(U_sum+M_sum))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),(U_sum+M_sum))
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),(U_sum+M_sum))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,2]),10),(U_sum+M_sum))

plot(log(unique((U_sum+M_sum)),10))

plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),log((U_sum+M_sum),10))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,3]),10),log((U_sum+M_sum),10))
abline(h=0,col="red")
abline(h=1,col="red")
abline(h=5,col="red")
abline(h=log(200000,10),col="red")
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log((U_sum+M_sum),10))
abline(a=c(0,-1),col="red")
plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),log((U_sum+M_sum),10),col=rep(rainbow(11),121))
lines(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,1]),10),log((U_sum+M_sum),10))
plot(log(unlist(par.mat.lambda1[,1]/par.mat.lambda1[,2]),10),log((U_sum+M_sum),10))


plot(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),log(U_sum,10),pch=3,col="red")
points(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),log(M_sum,10),pch=4,col="blue")
points(log(unlist(par.mat.lambda1[,5]/par.mat.lambda1[,3]),10),log((U_sum+M_sum),10),pch=20)

plot(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log(U_sum,10),pch=3,col="red")
points(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log(M_sum,10),pch=4,col="blue")
points(log(unlist(par.mat.lambda1[,2]/par.mat.lambda1[,3]),10),log((U_sum+M_sum),10),pch=20)

plot(log(U_sum,10),((log(M_sum,10)-log(U_sum,10))))
points(log(U_sum,10),log(M_sum,10),pch=20,col="red")

U_M_sum_ratio <- U_sum/M_sum
U_sum_percent <- U_sum/(U_sum+M_sum)*100
M_sum_percent <- M_sum/(U_sum+M_sum)*100

dev.off()
plot(sort(U_sum_percent))

cor(U_sum_percent,par.mat.lambda1[,kappa/mu])
cor(U_sum_percent,par.mat.lambda1[,1/mu])

hist(U_sum_percent)
hist(M_sum_percent)

plot(M_sum_percent,log(M_sum,10)) 

plot(log10(par.mat.lambda1[order(U_sum_percent),kappa/mu]),sort(U_sum_percent))

x1<- log10(par.mat.lambda1[delta==max(delta),kappa/mu])
y1 <- U_sum_percent[which(par.mat.lambda1[,delta]==max(par.mat.lambda1[,delta]))]
plot(x1[1:11],y1[1:11])
for(i in 1:11){
  lines(x1[((11*i)-10):(11*i)],y1[((11*i)-10):(11*i)],col=i)
}
legend("bottomright",legend =unique(par.mat.lambda1[,mu]),col = 1:11,lty = 1)



order(U_sum_percent)
par.mat.lambda1[,kappa/mu]




plot.bars.mat(unmarked = Unmarked_imax30,marked = Marked_imax30,series = 1211:1221,
              pars=par.mat.lambda1,main="Sum classes. l=1, t =1",log=F,beside = F,xaxt.step = 5)

plot.bars.mat(unmarked = Unmarked_imax30_I,marked = Marked_imax30_I,series = 1211:1221,
              pars=par.mat.lambda1,main="Sum classes I. l=1, t =1",log=F,beside = F,xaxt.step = 5)

plot.bars.mat(unmarked = Unmarked_imax30_NI,marked = Marked_imax30_NI,series = 1211:1221,
              pars=par.mat.lambda1,main="Sum classes NI. l=1, t =1",log=F,beside = F,xaxt.step = 5)


t(Unmarked_imax30/colSums(Unmarked_imax30))%*%(0:30)

hist((t(Unmarked_imax30)/colSums(Unmarked_imax30))*100,breaks = seq(0,100,1))

hist((t(Unmarked_imax30)/colSums(Unmarked_imax30))%*%(0:30),breaks = 0:30)





persp(x = log(par.mat.lambda1[,kappa],10),y = log(par.mat.lambda1[,mu],10),z = U_sum)

par.mat.lambda1[order(par.mat.lambda1[,kappa]),.(log(kappa,10),log(mu,10),index2)]

matrix(U_sum[par.mat.lambda1[order(par.mat.lambda1[,kappa]),index2]],nrow=11)[,1:11]

Unmarked_imax30[,par.mat.lambda1[mu==0.00001&delta==100000,index2]]
z<-Unmarked_imax30[,par.mat.lambda1[mu==0.00001&delta==100000,index2]]
z[which(z<1)] <- 1
persp(x=-5:5,y=0:30,z=t(log(z,10)),xlab = "kappa",ylab="ribosomes",zlab="mRNA",
      theta = 30,phi=30,xlim=c(-5,5),ylim=c(0,30),zlim=c(0,5),expand = 0.5, col = "lightblue")

min(log(z,10))


library(gplots)

heatmap.2(Unmarked_imax30,dendrogram = "none",Rowv = "none",Colv = "none",scale = "none",trace = "none")


x <- seq(-10, 10, length= 30)
y <- x
f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
z <- outer(x, y, f)
z[is.na(z)] <- 1
op <- par(bg = "white")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      ltheta = 120, shade = 0.75, ticktype = "detailed",
      xlab = "X", ylab = "Y", zlab = "Sinc( r )"
) -> res
round(res, 3)




# Decay rates -------------------------------------------------------------

#Test decay rate

MakeUnmarkedMatrix <- function(parameters,imax,rescale.fn = NULL,recscale.val = NULL){
  tau <- rep(parameters[4],imax+1)*seq(0,imax,1)/parameters[4]
  kappa <- rep(parameters[1],imax+1)*seq(1,0,-1/imax)/parameters[4]
  lambda <-  parameters[3]/parameters[4]
  mu <- parameters[2]/parameters[4]
  upper <- tau
  lower <- kappa
  main <- -1*(kappa+tau+mu)
  MAT <- matrix(0,ncol=imax+1,nrow=imax+1)
  for(i in 1:(imax+1)){
    if(i==1){
      MAT[1:2,i]  <- c(main[i],lower[i])
    }else if(i==(imax+1)){
      MAT[imax:(imax+1),i] <- c(upper[i],main[i])
    }else{
      MAT[(i-1):(i+1),i] <- c(upper[i],main[i],lower[i])
    }
  }
  return(MAT)
}

UMAT <- MakeUnmarkedMatrix(parameters = unlist(par.mat.lambda1[1096,]),imax = 30)

which(Unmarked_imax30[15,]==max(Unmarked_imax30[15,]))
plot(Unmarked_imax30[,1096])
SS <- Unmarked_imax30[,1096]
SS_M <- Marked_imax30[,1096]
UMAT%*%SS

SIM_UMAT <- function(UMAT,T0,imax,parameters,iter=10,no_input=F){
  counter <- 0
  if(no_input==T){
    lambda=0
  }else{
    lambda <- c(unlist(par.mat.lambda1[1096,3]),rep(0,imax))    
  }
  OUT_MAT <- matrix(0,nrow=imax+1,ncol=iter+1)
  D_MAT <- matrix(0,nrow=imax+1,ncol=iter+1)
  for(i in 1:(iter+1)){
    if(i==1){
      OUT_MAT[,i] <- T0
    }else{
      D_MAT[,i] <- UMAT%*%OUT_MAT[,i-1]
      INT <- OUT_MAT[,i-1] +  D_MAT[,i] + lambda
      index <- which(INT<0)
      if(length(index)>0){
        INT[index] <- 0
        counter <- counter + 1
      }
      OUT_MAT[,i] <- INT 
    }

  }
  print(counter)
  return(list(D_MAT,OUT_MAT))
}

SS_sim <- SIM_UMAT(UMAT,T0 = SS,imax=30,unlist(par.mat.lambda1[1096,]),iter = 10,no_input = F)

colSums(SS_sim)
apply(SS_sim,1,diff)

SS_sim <- SIM_UMAT(UMAT/20,T0 = SS,imax=30,unlist(par.mat.lambda1[1096,]),iter = 1500000,no_input = T)
cs <- colSums(SS_sim[[2]])

smp <- sample(colSums(SS_sim[[2]]),1000)
plot(sort(smp,decreasing = T))
plot(colSums(SS_sim[[2]]))
View(round(SS_sim[[1]],4))
View(SS_sim[[2]])

SS_new <- SS + UMAT%*%SS
sum(SS_new)
SS_new <- SS_new + UMAT%*%SS_new
sum(SS_new)


plot(1:31,SS_sim[,1],type="l")
points(1:31,SS_sim[,1])
for(i in 2:10){
  lines(1:31,SS_sim[,i],col=i)
  points(1:31,SS_sim[,i],col=i)
}

smp <- sample(1:1500000,1000)
plot(sort(smp/20),cs[sort(smp)])
lines(x=seq(0,1500000,1)/20,y=100000*exp(-0.00001*seq(0,1500000,10000)/20),col="red")

smp_ss_sim <- SS_sim[[2]][,seq(1,1500000,100000)]

plot(0:30,smp_ss_sim[,1],type="l")
points(0:30,smp_ss_sim[,1])
for(i in 2:10){
  lines(0:30,smp_ss_sim[,i],col=i)
  points(0:30,smp_ss_sim[,i],col=i)
}


pct_smp_ss_sim <-t(t(smp_ss_sim)/colSums(smp_ss_sim))

colSums(pct_smp_ss_sim)

apply(pct_smp_ss_sim,1,diff)

plot(0:30,pct_smp_ss_sim[,1],type="l")
points(0:30,pct_smp_ss_sim[,1])
for(i in 2:10){
  lines(0:30,pct_smp_ss_sim[,i],col=i)
  points(0:30,pct_smp_ss_sim[,i],col=i)
}


MakeMarkedMatrix <- function(parameters,imax,rescale.fn = NULL,recscale.val = NULL){
  tau <- rep(parameters[4],imax+1)*seq(0,imax,1)/parameters[4]
  kappa <- rep(parameters[1],imax+1)*seq(1,0,-1/imax)/parameters[4]
  lambda <-  parameters[3]/parameters[4]
  mu <- parameters[2]/parameters[4]
  delta <- parameters[5]/parameters[4]
  MATL <- diag(rep(mu,imax+1))
  MATR <- diag(-1*c(delta,tau[-1]))
  for(i in 1:(imax)){
      MATR[i,(i+1)] <- tau[i+1]
  }
  return(cbind(MATL,MATR))
}

MMAT <- MakeMarkedMatrix(parameters = unlist(par.mat.lambda1[1096,]),imax = 30)
c(SS,SS_M)
MMAT%*%c(SS,SS_M)


SIM_MMAT <- function(MMAT,T0_U,T0_M,imax,parameters,iter=10,no_input=F){
  counter <- 0
  OUT_MAT <- matrix(0,nrow=imax+1,ncol=iter+1)
  D_MAT <- matrix(0,nrow=imax+1,ncol=iter+1)
  if(no_input==T){
    MMAT <- MMAT[,(imax+2):ncol(MMAT)]
  }
  for(i in 1:(iter+1)){
    if(i==1){
      OUT_MAT[,i] <- T0_M
    }else{
      if(no_input==T){
        D_MAT[,i] <- MMAT%*%OUT_MAT[,i-1]
      }else{
        D_MAT[,i] <- MMAT%*%c(T0_U[,i-1],OUT_MAT[,i-1])        
      }

      INT <- OUT_MAT[,i-1] +  D_MAT[,i] 
      index <- which(INT<0)
      if(length(index)>0){
        #print(i)
        INT[index] <- 0   
        counter <- counter +1
      }
      OUT_MAT[,i] <- INT 
    }
    
  }
  print(counter)
  return(list(D_MAT,OUT_MAT))
}

SS_sim <- SIM_UMAT(UMAT/10000,T0 = SS,imax=30,unlist(par.mat.lambda1[1096,]),iter = 150000,no_input = T)
SS_M_sim <- SIM_MMAT(MMAT = MMAT/10000,T0_U = SS_sim[[2]],T0_M = SS_M,imax=30,parameters = unlist(par.mat.lambda1[1096,]),iter = ncol(SS_sim[[2]])-1)

SS_M_sim <- SIM_MMAT(MMAT = MMAT/10000,T0_U = SS_sim[[2]],T0_M = SS_M,imax=30,parameters = unlist(par.mat.lambda1[1096,]),iter = 250000,no_input = T)


cs_M <- colSums(SS_M_sim[[2]])
#plot((0:150000)/10000,cs_M)
smp <- sample(1:250000,10000)
plot(sort(smp/10000),cs_M[sort(smp)])
lines(x=seq(0,250000,10000)/10000,y=sum(SS_M)*exp(-1*tau_bar*seq(0,250000,10000)/10000),col="red",lwd=2)


smp_SS_M_sim <- SS_M_sim[[2]][,seq(0,250000,10000)]

plot(1:30,smp_SS_M_sim[-1,1],type="l")
points(1:30,smp_SS_M_sim[-1,1])
for(i in 2:10){
  lines(1:30,smp_SS_M_sim[-1,i],col=i)
  points(1:30,smp_SS_M_sim[-1,i],col=i)
}

pct_smp_SS_M_sim <-t(t(smp_SS_M_sim)/colSums(smp_SS_M_sim))

colSums(pct_smp_SS_M_sim)

round(apply(pct_smp_SS_M_sim,1,diff),10)

plot(1:30,pct_smp_SS_M_sim[-1,1],type="l",ylim=c(0,1.1))
points(1:30,pct_smp_SS_M_sim[-1,1])
for(i in 2:ncol(pct_smp_SS_M_sim[-1,])){
  lines(1:30,pct_smp_SS_M_sim[-1,i],col=i)
  points(1:30,pct_smp_SS_M_sim[-1,i],col=i)
}

#half-life based off only initial distribution
#Calculate half life for marked
w <- SS_M[-1]/sum(SS_M[-1])

t0.5 <- unlist(lapply(1:30,function(i) sum(1/1:i)*log(2)/unlist(par.mat.lambda1[1096,])[4]))
tbar <- t0.5%*%w
tau_bar <- as.vector(unlist(log(2)/tbar))


plot(sort(smp/10000),cs_M[sort(smp)])
lines(x=seq(0,250000,10000)/10000,y=sum(SS_M)*exp(-1*tau_bar*seq(0,250000,10000)/10000),col="green",lwd=2,lty=1)


plot(sort(smp/10000),log(cs_M[sort(smp)]))
lines(x=seq(0,250000,10000)/10000,y=log(sum(SS_M)*exp(-1*tau_bar*seq(0,250000,10000)/10000)),col="green",lwd=2,lty=1)


#half-life based off changing distribution
#class based half lives do not change as they are determined by the model
t0.5 <- unlist(lapply(1:30,function(i) sum(1/1:i)*log(2)/unlist(par.mat.lambda1[1096,])[4]))
#here the half life is estimated once per simutation timepoint
w <-t(t(SS_M_sim[[2]][-1,])/colSums(SS_M_sim[[2]][-1,]))
tbar <- t(w)%*%as.vector(t0.5)
tau_bar <- as.vector(unlist(log(2)/tbar))

m_decay <- lapply(seq(0,250000,10000)/10000,function(i) sum(SS_M)*exp(-1*mean(tau_bar[i+1])*i))

#plot(sort(smp/10000),cs_M[sort(smp)])
lines(x=seq(0,250000,10000)/10000,y=m_decay,col="red",lwd=2,lty=3)



#check assumption for single loaded class 5
SS_M_5 <- rep(0,31)
SS_M_5[6] <- 1
SS_M_sim <- SIM_MMAT(MMAT = MMAT/10000,T0_U = SS_sim[[2]],T0_M = SS_M_5,imax=30,parameters = unlist(par.mat.lambda1[1096,]),iter = 250000,no_input = T)

plot(seq(0,250000,10000)/10000,colSums(SS_M_sim[[2]])[seq(1,250001,10000)])

smp_SS_M_sim <- SS_M_sim[[2]][,seq(0,250000,10000)]

plot(1:30,smp_SS_M_sim[-1,1],type="l")
points(1:30,smp_SS_M_sim[-1,1])
for(i in 2:10){
  lines(1:30,smp_SS_M_sim[-1,i],col=i)
  points(1:30,smp_SS_M_sim[-1,i],col=i)
}

pct_smp_SS_M_sim <-t(t(smp_SS_M_sim)/colSums(smp_SS_M_sim))

colSums(pct_smp_SS_M_sim)

round(apply(pct_smp_SS_M_sim,1,diff),10)

plot(1:30,pct_smp_SS_M_sim[-1,1],type="l",ylim=c(0,1.1))
points(1:30,pct_smp_SS_M_sim[-1,1])
for(i in 2:ncol(pct_smp_SS_M_sim[-1,])){
  lines(1:30,pct_smp_SS_M_sim[-1,i],col=i)
  points(1:30,pct_smp_SS_M_sim[-1,i],col=i)
}


library("RConics")
# par.names <- c("kappa","mu","lambda","tau","delta","imax")
parameters <- c(0.66,0.0292,5.4,0.1,0.1,30)
UMAT <- MakeUnmarkedMatrix(parameters = parameters,imax = 30)
x<-adjoint(UMAT)%*%c(-parameters[3],rep(0,30))
y <- x/det(UMAT)
z <- CalcUnmarkedClassNP(parameters,imax = 30)
round(y,13)==round(z/10,13)
cbind(y,z/10)
cbind(y,z,a[2:32])
