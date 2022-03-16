source("R/Polysome_functions_RAUC.R")
library("viridis")
library(GenomicFeatures)

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
#histogram of Protein lengths in arabidopsis.
#BiocManager::install("GenomicFeatures")
rep.models <- fread("../Results/20220301/At_Representative_Gene_models.txt",sep="\t",header = T,quote="",stringsAsFactors = F)
setkeyv(rep.models,"Locus Identifier")
GTF <- rtracklayer::import.gff("../../RNA_SEQ_AVA/genome/Araport11_GFF3_genes_transposons.201606.gtf")
At.Tx.Db <- makeTxDbFromGRanges(GTF)
CDS <- as.data.table(select(x = At.Tx.Db,keys=sort(unlist(rep.models[,`Gene Model Name`])),columns = c("CDSSTART","CDSEND"),keytype = "TXNAME"))
CDS[,"CDS_LENGTH":= (CDSEND-CDSSTART)+1]
CDS_LENGTHS <- CDS[,sum(CDS_LENGTH),TXNAME]
setkeyv(CDS_LENGTHS,"TXNAME")
cbind(CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),],CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),2]%%3)
rep.models[as.character(gsub(x=unlist(CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),1]),replacement = "",pattern = "\\.\\d")),c(2,5)]
CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),V1:=V1-CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),2]%%3]
desc <- rep.models[as.character(gsub(x=unlist(CDS_LENGTHS[,1]),replacement = "",pattern = "\\.\\d")),c(1,3)]

#Remove genes marked as hypothetical (too short)
hypothetical <- as.numeric(unlist(lapply(desc[,2],function(i)grep(x = i,pattern = "^hypothetical protein"))))
desc <- desc[-unlist(hypothetical),]

#Remove genes marked as generic transmembrane (too short)
transmembrane <- as.numeric(unlist(lapply(desc[,2],function(i)grep(x = i,pattern = "^transmembrane"))))
desc <- desc[-unlist(transmembrane),]

#This gene encodes a small protein and has evolotuionary similarity to something (uncertain existance)
evo <- as.numeric(unlist(lapply(desc[,2],function(i)grep(x = i,pattern = "^This gene encodes a small protein and has"))))

hist(unlist(CDS_LENGTHS[-hypothetical,][-transmembrane,2]/3),col="red",breaks = seq(0,18000/3,100))
hist(unlist(CDS_LENGTHS[-hypothetical,][-transmembrane,][evo,2]/3),breaks = seq(0,18000/3,100),add=T,col="white")
quantile(sort(unlist(CDS_LENGTHS[-unlist(hypothetical),][-transmembrane,][-evo,2]/3)),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
protein_quant <- quantile(unlist(CDS_LENGTHS[-hypothetical,][-transmembrane,][-evo,2]),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)/3
hist(log(unlist(CDS_LENGTHS[-hypothetical,][-transmembrane,][-evo,2]/3),10),breaks = seq(0,log(18000/3,2),1/10),col="white")


comp_and_exp_w_annotation <- CDS_LENGTHS[-hypothetical,][-transmembrane,][-evo,]

old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
# Protein length distribution ---------------------------------------------
png("Figure1_At_Protein_length_histogram.png",width = 5,height = 5,units = "in",res = 300)
prot.length.hist <- hist(unlist(comp_trimmed[,2])/3,breaks = seq(0,1500,20),ylim=c(0,1100),col="white",axes=F,xlab="",main="Arabidopsis Protein Lengths")
abline(v=quantile(comp_trimmed[,2]/3,c(0.1,0.5,0.9),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,10^max(nchar(prot.length.hist$counts)),250),labels = seq(0,10^max(nchar(prot.length.hist$counts)),250))
axis(1,at=seq(0,1500,500),labels = c("0","500","1000",">1500"))
mtext(1,line=2.25,text = "Protein Length (a.a.)")
legend("topright",legend = c("Small","10% (140aa)","Medium","50% (368aa)","Large","90% (780aa)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
box()
dev.off()
setwd(old.dir)

# Half Life Distribution --------------------------------------------------
HL_At <- fread("../../Project_1/Results/20220117/Half life Sorenson.csv",stringsAsFactors = F,header=T)
hist(log(unlist(HL_At[which(unlist(HL_At[,3])<100000),3]),10))
median(unlist(HL_At[,3]),na.rm=T)
quantile(unlist(HL_At[which(unlist(HL_At[,3])<100000),3]),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
quantile(unlist(HL_At[,3]),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))

hist(unlist(HL_At[,3]))

HL_At <- fread("../Results/20220301/Sorenson_results_WT.csv",stringsAsFactors = F,header=T)

HL_minutes <- unlist(log(2)/HL_At[-which(HL_At[,2]==0),2])

median(sort(HL_minutes)[-c(17900:18585)])
HL_minutes <- sort(HL_minutes)[-c(17900:18585)]

quantile(HL_minutes,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
HL_minutes[(which(HL_minutes>1500))]<-1500

ylim=c(0,10)
old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Figure1_At_Transcript_half_life_histogram.png",width = 5,height = 5,units = "in",res = 300)
half.life.hist <- hist(HL_minutes,col="white",breaks = seq(0,1500,25),axes=F,xlab="",main="Arabidopsis Half-life")
abline(v=quantile(HL_minutes,c(0.1,0.5,0.9),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,10^max(nchar(half.life.hist$counts)),500),labels = seq(0,10^max(nchar(half.life.hist$counts)),500),las=2)
axis(1,at=seq(0,1500,500),labels = c("0","500","1000",">1500"))
mtext(1,line=2.25,text = "Transcript Half-life (min.)")
legend("topright",legend = c("Short","10% (33 min.)","Medium","50% (107 min.)","Long","90% (427 min.)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
box()
dev.off()
setwd(old.dir)

#Define parameter set. Since determinant method is being used I chose a very small mu. Prior exploration showed that a
#small mu behaves much like mu=0. The reason I am doing this is that calculating the adjugate matrix is slower than finding
#the full solution.
#mu values 1% 9.994377e-04 5% 5.068461e-04 10% 3.500904e-04 25% 1.973310e-04 50% 1.081415e-04
#75% 5.500915e-05 90% 2.702400e-05 95%  1.763636e-05  99% 8.043137e-06
par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.ranges <- list(c(seq(0.0001,0.5,0.0001)),c(0.01,9.994377*10^-4,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.081415*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,8.043137*10^-6) ,c(1),c(1),c(100000),c(1))
plot1.params <- CreateParMatrix(parameter.range = par.ranges,fixed = F,par.names = par.names)
plot1.params <- as.data.table(plot1.params)
plot1.params[,"index":=1:nrow(plot1.params)]
#calculate unmarked class
Unmarked_imax_median_plot1 <- apply(plot1.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = 41,rescale = T,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

#Calculate summary statistics
U_pct_imax_median_plot1 <- Ribo_percent(Unmarked_imax_median_plot1)
U_mean_imax_median_plot1 <- Ribo_mean(Unmarked_imax_median_plot1)
U_sd_imax_median_plot1 <- Ribo_sd(Unmarked_imax_median_plot1)
U_mode_imax_median_plot1 <- Ribo_mode(Unmarked_imax_median_plot1)

#plot curves for figure 1 example
peaks_mean <- round(U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]])
peaks_mode <- U_mode_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]]
#peaks_mean_index <- unlist(lapply(c(0,5,10,15,25,30), function(i) which(peaks_mean==i)[1]))
peaks_mean_index <- unlist(lapply(c(2,15,30), function(i) which(peaks_mean==i)[1]))
#peaks_mode_index <- unlist(lapply(c(0,5,10,15,25,30), function(i) which(peaks_mode==i)[1]))
#peaks_mean_index[1] <- 10
Initiation_rates <- unlist(plot1.params[peaks_mean_index,1])
names(Initiation_rates) <- c("Low","Mid","High")

old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Fig1_Ribosome_loading_no_int.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=1.1,col=col[1],ylim=c(0,35),
     main ="No interference, median protein size, Long Half-life")
points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[i]]),col=col[i],lwd=2)
  points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 21,36),y=c(25,10,8),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
kappas <- c(round(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index[1]],4),round(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index[2:6]],3))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
text(x = 0.5,y=0.8,"Scaled Initiation",font=2)
legend("center", legend = c("Low","Medium","High"),fill=col)
dev.off()

# png("Fig1_Ribosome_loading_no_int_bw.png",width = 5,height = 5,units = "in",res=300)
# plot(0:30,Ribo_percent(Unmarked_imax30_plot1[,plot1.params[tau==1,index]][,peaks_mean_index[1]]),type="l",ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,
#      main = bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
# for(i in 2:length(peaks_mean_index)){
#   lines(0:30,Ribo_percent(Unmarked_imax30_plot1[,plot1.params[tau==1,index]][,peaks_mean_index[i]]),lty=i,lwd=2)
# }
# dev.off()
#
# png("Fig1_Ribosome_loading_no_int_legend_bw.png",width = 1.5,height = 2,units = "in",res=300)
# par(mar=c(0.1,0.1,0.1,0.1))
# kappas <- c(round(plot1.params[tau==1,kappa][peaks_mean_index[1]],4),round(plot1.params[tau==1,kappa][peaks_mean_index[2:6]],3))
# plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
# legend("left", legend = as.expression(lapply(kappas, function(x) bquote(kappa==.(x)))),col=1,lty=1:6)
# dev.off()
#
#
#
# png("Fig1_Ribosome_loading_no_int_legend.png",width = 1.5,height = 2,units = "in",res=300)
# par(mar=c(0.1,0.1,0.1,0.1))
# kappas <- c(round(plot1.params[tau==1,kappa][peaks_mean_index[1]],4),round(plot1.params[tau==1,kappa][peaks_mean_index[2:6]],3))
# plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
# legend("left", legend = as.expression(lapply(kappas, function(x) bquote(kappa==.(x)))),col=1:6,lty=1)
# dev.off()

#
# bp.data <- Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[tau==1,index]][,peaks_mean_index])
# bp.data.trim <- bp.data
# bp.data.trim[which(bp.data.trim>0.3)] <- bp.data.trim[which(bp.data.trim>0.3)] -0.2
#bp.data.trim[1,1] <- bp.data.trim[1,1]-0.15

# png("Fig1_Ribosome_loading_no_int_bp_cut.png",width = 5,height = 5,units = "in",res=300)
# bp.tmp <-barplot(t(bp.data.trim),beside = T,col=c("black","grey80","grey40","grey60","grey20","white"),ylim=c(0,0.65),width = 2,axes=F,main=bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
# mtext(side = 1,at = colMeans(bp.tmp)[seq(1,31,2)],text = seq(0,30,2))
# mtext(side = 1,at = colMeans(bp.tmp)[seq(2,31,2)],text = seq(1,30,2),line=1)
# mtext(side = 1,at=mean(colMeans(bp.tmp)),text = "Ribosomal class",line=2)
# mtext(side = 2,at=0.3,text = "Percent Ribosomal Load",line=2.5)
# axis(2,at = c(seq(0,0.3,0.1),0.4,0.5,0.6),labels = c(0,0.1,0.2,0.3,0.6,0.7,0.8),las=2)
# par(xpd=T)
# polygon(x=c(-20,max(bp.tmp)+10,max(bp.tmp)+10,-20,-20),y=c(0.3,0.3,0.395,0.395,0.3),col = "white",border = rgb(0,0,0,0))
# lines(x = c(-20,-10),y=-0.001*(c(-20,-10)+20)+0.402,lwd=2)
# lines(x = c(-20,-10),y=-0.001*(c(-20,-10)+20)+0.302,lwd=2)
# box()
# dev.off()
#
# png("Fig1_Ribosome_loading_no_int_legend_bp.png",width = 1.5,height = 2,units = "in",res=300)
# par(mar=c(0.1,0.1,0.1,0.1))
# kappas <- c(round(plot1.params[tau==1,kappa][peaks_mean_index[1]],4),round(plot1.params[tau==1,kappa][peaks_mean_index[2:6]],3))
# plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
# legend("left", legend = as.expression(lapply(kappas, function(x) bquote(kappa==.(x)))),fill=c("black","grey80","grey40","grey60","grey20","white"))
# dev.off()
#
# png("Fig1_Ribosome_loading_no_int_bp_no_cut.png",width = 5,height = 5,units = "in",res=300)
# bp.data <- Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[tau==1,index]][,peaks_mean_index])
# bp.tmp <-barplot(t(bp.data),beside = T,col=c("black","grey80","grey40","grey60","grey20","white"),ylim=c(0,0.8),width = 2,main=bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
# mtext(side = 1,at = colMeans(bp.tmp)[seq(1,31,2)],text = seq(0,30,2))
# mtext(side = 1,at = colMeans(bp.tmp)[seq(2,31,2)],text = seq(1,30,2),line=1)
# mtext(side = 1,at=mean(colMeans(bp.tmp)),text = "Ribosomal class",line=2)
# mtext(side = 2,at=0.3,text = "Percent Ribosomal Load",line=2.5)
# box()
# dev.off()

#SD plot
png("Fig1_SD.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
points(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_SD_logx.png",width = 5,height = 5,units = "in",res=300)
plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
points(log(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],10),U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(-2.6,-1.8,-1.35),y=c(3,5,7),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_SD_mu_range.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
counter <- 2
for(i in par.ranges[[2]][-1]){
  lines(plot1.params[mu==i,kappa],U_sd_imax_median_plot1[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
  counter <- counter + 1
}
dev.off()

png("Fig1_SD_mu_range_logx.png",width = 5,height = 5,units = "in",res=300)
plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
counter <- 2
for(i in par.ranges[[2]][-1]){
  lines(log(plot1.params[mu==i,kappa],10),U_sd_imax_median_plot1[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
  counter <- counter + 1
}
dev.off()

#Mean plot
png("Fig1_mean.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Long Half-life")
counter <- 2
points(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_mean_logx.png",width = 5,height = 5,units = "in",res=300)
plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote(log[10]*" Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Long Half-life")
points(log(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(-2.6,-1.8,-1.35),y=c(5,16,31),labels = c("Low","Medium","High"),font=2)
dev.off()

#half life quantiles and mu ranges
quantile(HL_minutes,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
log(2)/(quantile(HL_minutes,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)*60)

png("Fig1_mean_mu_range.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Long Half-life",col=cols2)
counter <- 2
points(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,4,2,1)]){
  lines(plot1.params[mu==i,kappa],U_mean_imax_median_plot1[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
  points(plot1.params[mu==i,kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
  counter <- counter + 1
}
legend("bottomright",fill=cols2,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
dev.off()

png("Fig1_mean_mu_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Long Half-life",lty=1)
counter <- 2
points(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,4,2,1)]){
  lines(plot1.params[mu==i,kappa],U_mean_imax_median_plot1[plot1.params[mu==i,index]],lty=counter,lwd=2)
  points(plot1.params[mu==i,kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
  counter <- counter + 1
}
legend("bottomright",lty = 1:5,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
dev.off()

png("Fig1_mean_mu_range_logx.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[mu==par.ranges[[2]][8],kappa],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Half-life range",col=cols2[1])
counter <- 2
points(log(plot1.params[mu==par.ranges[[2]][8],kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]][peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,4,2,1)]){
  lines(log(plot1.params[mu==i,kappa],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
  points(log(plot1.params[mu==i,kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
#  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",fill = cols2,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
dev.off()

png("Fig1_mean_mu_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[mu==par.ranges[[2]][8],kappa],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Half-life range",lty=1)
counter <- 2
points(log(plot1.params[mu==par.ranges[[2]][8],kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]][peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,4,2,1)]){
  lines(log(plot1.params[mu==i,kappa],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]],lty=counter,lwd=2)
  points(log(plot1.params[mu==i,kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
dev.off()


#Mean and SD legend
png("Fig1_Ribosome_mean_sd_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))
dev.off()











#--------------------------------------------------------
#imax range plot, compare to tau=1 range of params with imax 30
protein_quant
protein_quant/9
#imax will be 8,15,41,87,164 The 1%, 10%, 50%, 90%, and 99% quantiles of arabidopsis
imax.range <- c(8,15,41,87,164)
imax.percentiles <- c("1%","10%","50%","90%","99%")
Unmarked_imax_range_plot1 <- lapply(imax.range,function(j) apply(plot1.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))
names(Unmarked_imax_range_plot1) <- imax.percentiles

#Calculate summary statistics
U_pct_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_percent)
U_mean_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_mean)
U_sd_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_sd)

#profile plots imax range
peak_positions <- c(0.01,0.10,0.5,0.90,0.99)

peaks_mean_imax_range <- lapply(1:length(imax.range),function(i) round(U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]]/imax.range[i],2))
peaks_mean_index_imax_range <- lapply(1:length(imax.range), function(i) unlist(lapply(peak_positions,function(j) which(peaks_mean_imax_range[[i]]==j)[1])))














#plot curves for figure 1 imax range example

png("Fig1_Ribosome_loading_no_int_imax_range_halfway_loading.png",width = 5,height = 5,units = "in",res=300)
col3 = cividis(5)
plot((0:imax.range[3])/imax.range[3],Ribo_percent(Unmarked_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][10],index],][,peaks_mean_index_imax_range[[3]][3]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,ylim=c(0,0.3),col=col3[3],
     main = bquote("No interference, imax=30, "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in c(1,2,4,5)){
  lines((0:imax.range[i])/imax.range[i],Ribo_percent(Unmarked_imax_range_plot1[[i]][,peaks_mean_index_imax_range[[i]][3]]),col=col3[counter],lwd=2)
  counter <- counter+1
}
abline(v=peak_positions[3],col="grey50",lty=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_legend_imax_range.png",width = 1.75,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = paste("imax",imax.range,sep=" = "),col=c(2:4,1,5:8),lty=1)
dev.off()



png("Fig1_Ribosome_loading_no_int_imax_range_kappa=0.05.png",width = 5,height = 5,units = "in",res=300)
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
png("Fig1_SD_imax_range.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params.imaxrange[,kappa],U_sd_imax_median_plot1[plot1.params[tau==0.1,index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
counter <- 2
counter <- 2
for(i in par.ranges[[4]][1:4]){
  lines(plot1.params[tau==i,kappa],U_sd_imax_median_plot1[plot1.params[tau==i,index]],lty=counter,lwd=2)
  counter <- counter + 1
}
counter <- 1
for(i in par.ranges[[4]][5:8]){
  lines(plot1.params[tau==i,kappa],U_sd_imax_median_plot1[plot1.params[tau==i,index]],lty=counter,col="grey50",lwd=2)
  counter <- counter + 1
}
dev.off()

#Mean plot
png("Fig1_mean_imax_range_tau=1.png",width = 5,height = 5,units = "in",res=300)
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

png("Fig1_mean_imax_range_tau=1.png",width = 5,height = 5,units = "in",res=300)
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

png("Fig1_Ribosome_loading_no_int_imax_range_means.png",width = 5,height = 5,units = "in",res=300)
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
png("Fig1_Ribosome_mean_sd_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))
dev.off()




setwd(old.dir)


# Figure 2 ----------------------------------------------------------------
#




