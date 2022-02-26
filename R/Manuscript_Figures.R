source("R/Polysome_functions_RAUC.R")
library("viridis")

# Focused parameters ------------------------------------------------------

#Tau will be 0.1aa/s, 1aa/s, 5aa/s, 10aa/s, 20aa/s
#Kappa will be 0.1%, 1%, 25%, 50%, 75%, 100% and 200% of Each respective tau
#Kappa will be 0.1%, 1%, 25%, 50%, 75%, 100% and 200% of Tau = 1aa/s for each tau >=1aa/s
#Mu will be based on the following half lives from sorenson
#12min, 23min, 59min, 109min, 220min, 830min, 4018min The 1%, 5%, 25%, 50%, 75%, 95%, and 99% quantiles of arabidopsis
#decay constant mu = ln(2)/half life in seconds
#imax will be 6, 10 ,24, 40, 60, 112, 181 The 1%, 5%, 25%, 50%, 75%, 95%, and 99% quantiles of arabidopsis

HL_for_mu <- c(12,23,59,109,220,830,4018)*60
mus <- log(2)/HL_for_mu
par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.ranges <- list(c(0.0001,0.0005,0.001,0.005,0.01,0.25,0.5,0.75,1,2),c(mus),c(1),c(0.1,1,5,10,20),c(100000),c(1))

Focused.par.setk <- CreateParMatrix(parameter.range = par.ranges,fixed = F,par.names = par.names)
Focused.par.setk <- as.data.table(Focused.par.setk)
Focused.par.freek <- CreateParMatrix(parameter.range = par.ranges,fixed = F,par.names = par.names)
Focused.par.freek <- as.data.table(Focused.par.freek)
Focused.par.freek[,"kappa":=kappa*tau]

Unmarked_imax30_fk <- apply(Focused.par.freek,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax30_fk <- apply(Focused.par.freek,1,function(i) CalcMarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

Unmarked_imax30_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax30_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))


Unmarked_imax6_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 6,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax6_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 6,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

Unmarked_imax10_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 10,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax10_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 10,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

Unmarked_imax24_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 24,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax24_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 24,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

Unmarked_imax40_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 40,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax40_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 40,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

Unmarked_imax60_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 60,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax60_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 60,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

Unmarked_imax112_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 112,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax112_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 112,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

Unmarked_imax181_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax181_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))



library(RConics)


CalcUnmarkedAdj_noM <- function(parameters,imax,rescale = F ,rescale.val = "tau",
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
  main <- -1*(kappa+tau)
  MAT <- matrix(0,ncol=imax+1,nrow=imax+1)
  for(i in 1:(imax+1)){
    if(i==1){
      MAT[1:2,i]  <- c(main[i],lower[i])
    }else if(i==(imax+1)){
      MAT[imax:(imax+1),i] <- c(upper[i-1],main[i])
    }else{
      MAT[(i-1):(i+1),i] <- c(upper[i-1],main[i],lower[i])
    }
  }
  adj <- adjoint(MAT)%*%c(1,rep(0,imax))
  return(adj)
}




tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,0.11111,0.0001)),0.0000001,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,4,0.005)),0.0001,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,4,0.005)),1,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,0.1,0.0001)),0.01,1,1,100000,1),par.names = par.names)


# Manuscript figures ------------------------------------------------------

# Figure 1 ----------------------------------------------------------------

#Define parameter set. Since determinant method is being used I chose a very small mu. Prior exploration showed that a
#small mu behaves much like mu=0. The reason I am doing this is that calculating the adjugate matrix is slower than finding
#the full solution.
par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.ranges <- list(c(seq(0.000001,0.5,0.0001)),0.00000001,c(1),c(0.1,0.25,0.5,1,2.5,5,10,20),c(100000),c(1))
plot1.params <- CreateParMatrix(parameter.range = par.ranges,fixed = F,par.names = par.names)
plot1.params <- as.data.table(plot1.params)
plot1.params[,"index":=1:nrow(plot1.params)]
#calculate unmarked class
Unmarked_imax30_plot1 <- apply(plot1.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

#Calculate summary statistics
U_pct_imax30_plot1 <- Ribo_percent(Unmarked_imax30_plot1)
U_mean_imax30_plot1 <- Ribo_mean(Unmarked_imax30_plot1)
U_sd_imax30_plot1 <- Ribo_sd(Unmarked_imax30_plot1)
U_mode_imax30_plot1 <- Ribo_mode(Unmarked_imax30_plot1)

#plot curves for figure 1 example
peaks_mean <- round(U_mean_imax30_plot1[plot1.params[tau==1,index]])
peaks_mode <- U_mode_imax30_plot1[plot1.params[tau==1,index]]
peaks_mean_index <- unlist(lapply(c(0,5,10,15,25,30), function(i) which(peaks_mean==i)[1]))
peaks_mode_index <- unlist(lapply(c(0,5,10,15,25,30), function(i) which(peaks_mode==i)[1]))
peaks_mean_index[1] <- 10

png("../Results/20220218/Fig1_Ribosome_loading_no_int.png",width = 5,height = 5,units = "in",res=300)
plot(0:30,Ribo_percent(Unmarked_imax30_plot1[,plot1.params[tau==1,index]][,peaks_mean_index[1]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,
     main = bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
for(i in 2:length(peaks_mean_index)){
  lines(0:30,Ribo_percent(Unmarked_imax30_plot1[,plot1.params[tau==1,index]][,peaks_mean_index[i]]),col=i,lwd=2)
}
dev.off()

png("../Results/20220218/Fig1_Ribosome_loading_no_int_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
kappas <- c(round(plot1.params[tau==1,kappa][peaks_mean_index[1]],4),round(plot1.params[tau==1,kappa][peaks_mean_index[2:6]],3))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(kappas, function(x) bquote(kappa==.(x)))),col=1:6,lty=1)
dev.off()

png("../Results/20220218/Fig1_Ribosome_loading_no_int_bw.png",width = 5,height = 5,units = "in",res=300)
plot(0:30,Ribo_percent(Unmarked_imax30_plot1[,plot1.params[tau==1,index]][,peaks_mean_index[1]]),type="l",ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,
     main = bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
for(i in 2:length(peaks_mean_index)){
  lines(0:30,Ribo_percent(Unmarked_imax30_plot1[,plot1.params[tau==1,index]][,peaks_mean_index[i]]),lty=i,lwd=2)
}
dev.off()

png("../Results/20220218/Fig1_Ribosome_loading_no_int_legend_bw.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
kappas <- c(round(plot1.params[tau==1,kappa][peaks_mean_index[1]],4),round(plot1.params[tau==1,kappa][peaks_mean_index[2:6]],3))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(kappas, function(x) bquote(kappa==.(x)))),col=1,lty=1:6)
dev.off()



png("../Results/20220218/Fig1_Ribosome_loading_no_int_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
kappas <- c(round(plot1.params[tau==1,kappa][peaks_mean_index[1]],4),round(plot1.params[tau==1,kappa][peaks_mean_index[2:6]],3))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(kappas, function(x) bquote(kappa==.(x)))),col=1:6,lty=1)
dev.off()


bp.data <- Ribo_percent(Unmarked_imax30_plot1[,plot1.params[tau==1,index]][,peaks_mean_index])
bp.data.trim <- bp.data
bp.data.trim[which(bp.data.trim>0.3)] <- bp.data.trim[which(bp.data.trim>0.3)] -0.2
#bp.data.trim[1,1] <- bp.data.trim[1,1]-0.15

png("../Results/20220218/Fig1_Ribosome_loading_no_int_bp_cut.png",width = 5,height = 5,units = "in",res=300)
bp.tmp <-barplot(t(bp.data.trim),beside = T,col=c("black","grey80","grey40","grey60","grey20","white"),ylim=c(0,0.65),width = 2,axes=F,main=bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
mtext(side = 1,at = colMeans(bp.tmp)[seq(1,31,2)],text = seq(0,30,2))
mtext(side = 1,at = colMeans(bp.tmp)[seq(2,31,2)],text = seq(1,30,2),line=1)
mtext(side = 1,at=mean(colMeans(bp.tmp)),text = "Ribosomal class",line=2)
mtext(side = 2,at=0.3,text = "Percent Ribosomal Load",line=2.5)
axis(2,at = c(seq(0,0.3,0.1),0.4,0.5,0.6),labels = c(0,0.1,0.2,0.3,0.6,0.7,0.8),las=2)
par(xpd=T)
polygon(x=c(-20,max(bp.tmp)+10,max(bp.tmp)+10,-20,-20),y=c(0.3,0.3,0.395,0.395,0.3),col = "white",border = rgb(0,0,0,0))
lines(x = c(-20,-10),y=-0.001*(c(-20,-10)+20)+0.402,lwd=2)
lines(x = c(-20,-10),y=-0.001*(c(-20,-10)+20)+0.302,lwd=2)
box()
dev.off()

png("../Results/20220218/Fig1_Ribosome_loading_no_int_legend_bp.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
kappas <- c(round(plot1.params[tau==1,kappa][peaks_mean_index[1]],4),round(plot1.params[tau==1,kappa][peaks_mean_index[2:6]],3))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(kappas, function(x) bquote(kappa==.(x)))),fill=c("black","grey80","grey40","grey60","grey20","white"))
dev.off()

png("../Results/20220218/Fig1_Ribosome_loading_no_int_bp_no_cut.png",width = 5,height = 5,units = "in",res=300)
bp.data <- Ribo_percent(Unmarked_imax30_plot1[,plot1.params[tau==1,index]][,peaks_mean_index])
bp.tmp <-barplot(t(bp.data),beside = T,col=c("black","grey80","grey40","grey60","grey20","white"),ylim=c(0,0.8),width = 2,main=bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
mtext(side = 1,at = colMeans(bp.tmp)[seq(1,31,2)],text = seq(0,30,2))
mtext(side = 1,at = colMeans(bp.tmp)[seq(2,31,2)],text = seq(1,30,2),line=1)
mtext(side = 1,at=mean(colMeans(bp.tmp)),text = "Ribosomal class",line=2)
mtext(side = 2,at=0.3,text = "Percent Ribosomal Load",line=2.5)
box()
dev.off()

#SD plot
png("../Results/20220218/Fig1_SD_tau_range.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[tau==0.1,kappa],U_sd_imax30_plot1[plot1.params[tau==0.1,index]],type="l",ylim=c(0,max(U_sd_imax30_plot1)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
counter <- 2
counter <- 2
for(i in par.ranges[[4]][1:4]){
  lines(plot1.params[tau==i,kappa],U_sd_imax30_plot1[plot1.params[tau==i,index]],lty=counter,lwd=2)
  counter <- counter + 1
}
counter <- 1
for(i in par.ranges[[4]][5:8]){
  lines(plot1.params[tau==i,kappa],U_sd_imax30_plot1[plot1.params[tau==i,index]],lty=counter,col="grey50",lwd=2)
  counter <- counter + 1
}
dev.off()
#Mean plot
png("../Results/20220218/Fig1_mean_tau_range.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[tau==0.1,kappa],U_mean_imax30_plot1[plot1.params[tau==0.1,index]],type="l",ylim=c(0,max(U_mean_imax30_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
     main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in par.ranges[[4]][1:4]){
  lines(plot1.params[tau==i,kappa],U_mean_imax30_plot1[plot1.params[tau==i,index]],lty=counter,lwd=2)
  counter <- counter + 1
}
counter <- 1
for(i in par.ranges[[4]][5:8]){
  lines(plot1.params[tau==i,kappa],U_mean_imax30_plot1[plot1.params[tau==i,index]],lty=counter,col="grey50",lwd=2)
  counter <- counter + 1
}
dev.off()


plot(plot1.params[tau==0.1,kappa]/plot1.params[tau==0.1,tau],U_mean_imax30_plot1[plot1.params[tau==0.1,index]],type="l",ylim=c(0,max(U_mean_imax30_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
     main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in par.ranges[[4]][1:4]){
  lines(plot1.params[tau==i,kappa]/plot1.params[tau==i,tau],U_mean_imax30_plot1[plot1.params[tau==i,index]],lty=counter,lwd=2,col=counter)
  counter <- counter + 1
}
counter <- 1
for(i in par.ranges[[4]][5:8]){
  lines(plot1.params[tau==i,kappa]/plot1.params[tau==i,tau],U_mean_imax30_plot1[plot1.params[tau==i,index]],lty=counter,lwd=2,col=counter)
  counter <- counter + 1
}

#Mean and SD legend
png("../Results/20220218/Fig1_Ribosome_mean_sd_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))
dev.off()
#--------------------------------------------------------
#imax range plot, compare to tau=1 range of params with imax 30
plot1.params.imaxrange <- plot1.params[tau==1,-c(7)]
#imax will be 6, 10 ,24, 40, 60, 112, 181 The 1%, 5%, 25%, 50%, 75%, 95%, and 99% quantiles of arabidopsis
imax.range <- c(6,10,24,30,40,60,112,181)
imax.percentiles <- c("1%","5%","25%","34%","50%","75%","95%","99%")
Unmarked_imax_range_plot1 <- lapply(imax.range,function(j) apply(plot1.params.imaxrange,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))
names(Unmarked_imax_range_plot1) <- imax.percentiles

#Calculate summary statistics
U_pct_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_percent)
U_mean_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_mean)
U_sd_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_sd)

#profile plots imax range
peak_positions <- c(0.05,0.25,0.5,0.75,0.95)

peaks_mean_imax_range <- lapply(1:length(imax.range),function(i) round(U_mean_imax_range_plot1[[i]]/imax.range[i],2))
peaks_mean_index_imax_range <- lapply(1:length(imax.range), function(i) unlist(lapply(peak_positions,function(j) which(peaks_mean_imax_range[[i]]==j)[1])))

#plot curves for figure 1 imax range example



png("../Results/20220224/Fig1_Ribosome_loading_no_int_imax_range_halfway_loading.png",width = 5,height = 5,units = "in",res=300)
 plot((0:imax.range[4])/imax.range[4],Ribo_percent(Unmarked_imax_range_plot1[[4]][,peaks_mean_index_imax_range[[4]][3]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,ylim=c(0,0.3),
     main = bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in c(1:3,5:8)){
  lines((0:imax.range[i])/imax.range[i],Ribo_percent(Unmarked_imax_range_plot1[[i]][,peaks_mean_index_imax_range[[i]][3]]),col=counter,lwd=2)
  counter <- counter+1
}
abline(v=peak_positions[3],col="grey50",lty=2)
dev.off()

png("../Results/20220224/Fig1_Ribosome_loading_no_int_legend_imax_range.png",width = 1.75,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = paste("imax",imax.range,sep=" = "),col=c(2:4,1,5:8),lty=1)
dev.off()



png("../Results/20220224/Fig1_Ribosome_loading_no_int_imax_range_kappa=0.05.png",width = 5,height = 5,units = "in",res=300)
plot((0:imax.range[4]),Ribo_percent(Unmarked_imax_range_plot1[[4]][,501]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,ylim=c(0,0.3),xlim=c(0,105),
     main = bquote("No interference, imax=30, "*kappa*"="*.(round(plot1.params.imaxrange[501,kappa],2))*", "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),cex.main=0.95)
counter <- 2
for(i in c(1:3,5:8)){
  lines((0:imax.range[i]),Ribo_percent(Unmarked_imax_range_plot1[[i]][,501]),col=counter,lwd=2)
  counter <- counter+1
}
dev.off()



#SD plot imax range
png("../Results/20220218/Fig1_SD_imax_range.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params.imaxrange[,kappa],U_sd_imax30_plot1[plot1.params[tau==0.1,index]],type="l",ylim=c(0,max(U_sd_imax30_plot1)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
counter <- 2
counter <- 2
for(i in par.ranges[[4]][1:4]){
  lines(plot1.params[tau==i,kappa],U_sd_imax30_plot1[plot1.params[tau==i,index]],lty=counter,lwd=2)
  counter <- counter + 1
}
counter <- 1
for(i in par.ranges[[4]][5:8]){
  lines(plot1.params[tau==i,kappa],U_sd_imax30_plot1[plot1.params[tau==i,index]],lty=counter,col="grey50",lwd=2)
  counter <- counter + 1
}
dev.off()

#Mean plot
png("../Results/20220218/Fig1_mean_imax_range_tau=1.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[4]][1:2500]/imax.range[4])*100,type="l",ylim=c(0,100),
     las=1,ylab="Mean Ribsome Load Percent",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
     main=bquote("No interference, imax range, "*tau*"="*1*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in c(1:3)){
  lines(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[i]][1:2500]/imax.range[i])*100,lty=counter,lwd=2)
  counter <- counter + 1
}
counter <- 1
for(i in 5:8){
  lines(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[i]][1:2500]/imax.range[i])*100,col="grey50",lty=counter,lwd=2)
  counter <- counter + 1
}
dev.off()

png("../Results/20220224/Fig1_mean_imax_range_tau=1.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[4]][1:2500]/imax.range[4])*100,type="l",ylim=c(0,100),
     las=1,ylab="Mean Ribsome Load Percent",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
     main=bquote("No interference, imax range, "*tau*"="*1*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in c(1:3,5:8)){
  lines(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[i]][1:2500]/imax.range[i])*100,col=counter,lwd=2)
  counter <- counter + 1
}
dev.off()


plot(log(plot1.params.imaxrange[500:2500,kappa],10),(U_mean_imax_range_plot1[[4]][500:2500]/imax.range[4])*100,type="l",ylim=c(0,100),
     las=1,ylab="Mean Ribsome Load Percent",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
     main=bquote("No interference, imax range, "*tau*"="*1*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in c(1:3,5:8)){
  lines(log(plot1.params.imaxrange[500:2500,kappa],10),(U_mean_imax_range_plot1[[i]][500:2500]/imax.range[i])*100,col=counter,lwd=2)
  counter <- counter + 1
}

png("../Results/20220224/Fig1_Ribosome_loading_no_int_imax_range_means.png.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params.imaxrange[,kappa],U_mean_imax_range_plot1[[4]],type="l",ylim=c(0,181),
     las=1,ylab="Mean Ribsome Load Percent",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
     main=bquote("No interference, imax range, "*tau*"="*1*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in c(1:3,5:8)){
  lines(plot1.params.imaxrange[,kappa],U_mean_imax_range_plot1[[i]],col=counter,lwd=2)
  counter <- counter + 1
}
dev.off()

#Mean and SD legend
png("../Results/20220218/Fig1_Ribosome_mean_sd_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))
dev.off()







# Figure 2 ----------------------------------------------------------------
#




