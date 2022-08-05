source("R/Polysome_functions_RAUC.R")
library("viridis")
library(GenomicFeatures)

# Manuscript figures ------------------------------------------------------

# Figure 1 ----------------------------------------------------------------
#histogram of Protein lengths in arabidopsis.
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
comp_trimmed <- comp_and_exp_w_annotation
comp_trimmed[which(comp_trimmed[,2]>4500),2] <- 4500

old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
# Protein length distribution ---------------------------------------------
png("Figure1_At_Protein_length_histogram.png",width = 5,height = 5,units = "in",res = 300)
prot.length.hist <- hist(unlist(comp_trimmed[,2])/3,breaks = seq(0,1500,20),ylim=c(0,1100),col="white",axes=F,xlab="",main="Arabidopsis Protein Lengths")
abline(v=quantile(comp_trimmed[,2]/3,c(0.1,0.5,0.95),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,10^max(nchar(prot.length.hist$counts)),250),labels = seq(0,10^max(nchar(prot.length.hist$counts)),250))
axis(1,at=seq(0,1500,500),labels = c("0","500","1000",">1500"))
mtext(1,line=2.25,text = "Protein Length (a.a.)")
legend("topright",legend = c("Small","10% (140aa)","Medium","50% (368aa)","Large","95% (973aa)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
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
abline(v=quantile(HL_minutes,c(0.1,0.5,0.99),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,10^max(nchar(half.life.hist$counts)),500),labels = seq(0,10^max(nchar(half.life.hist$counts)),500),las=2)
axis(1,at=seq(0,1500,500),labels = c("0","500","1000",">1500"))
mtext(1,line=2.25,text = "Transcript Half-life (min.)")
legend("topright",legend = c("Short","10% (33 min.)","Medium","50% (107 min.)","Long","95% (655 min.)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
box()
dev.off()
setwd(old.dir)


# Figure 1 no int ---------------------------------------------------------


#Define parameter set. Since determinant method is being used I chose a very small mu. Prior exploration showed that a
#small mu behaves much like mu=0. The reason I am doing this is that calculating the adjugate matrix is slower than finding
#the full solution.
#mu values Hypothetical fast 1e-2  1% 9.994377e-04 5% 5.068461e-04 10% 3.500904e-04 25% 1.973310e-04 50% 1.081415e-04
#75% 5.500915e-05 90% 2.702400e-05 95%  1.763636e-05  99% 8.043137e-06
#Late March. Mike suggested keeping it to more round numbers ~99% 1500, ~50% 115, ~1% 11.5, sim 1.15
par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.ranges <- list(c(seq(0.0001,0.5,0.0001)),c(0.01,1.004561*10^-3,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.004561*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,7.701635*10^-6) ,c(1),c(1),c(100000),c(1))
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
peaks_mean_index <- unlist(lapply(c(2,21,41), function(i) which(peaks_mean==i)[1]))
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
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,70),
     main ="Density Independent Initiation, median protein size, Long Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[i]]),col=col[i],lwd=2)
 # points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
kappas <- c(round(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index[1]],4),round(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index[2:6]],3))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
text(x = 0.5,y=0.8,"Scaled Initiation",font=2)
legend("center", legend = c("Low","Medium","High"),fill=col)
dev.off()

#SD plots
png("Fig1_SD.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main="Density Independent Initiation, median protein size, Long Half-life",lwd=2,cex.lab=1.2,font=2,cex.main=0.85)
points(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(0.05,0.08,0.21),y=c(1.2,6.5,3),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_SD_logx.png",width = 5,height = 5,units = "in",res=300)
plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
     xlab=bquote("log"[10]*(kappa)),ylab="Ribsomal class sd",main="Density Independent Initiation, median protein size, Long Half-life",lwd=2,cex.lab=1.2,font=2,cex.main=0.85)
points(log(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],10),U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(-2.6,-1.75,-0.45),y=c(2,5,2.5),labels = c("Low","Medium","High"),font=2)
dev.off()

cols2 <- cividis(3)

# png("Fig1_SD_mu_range.png",width = 5,height = 5,units = "in",res=300)
# plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
#      xlab=bquote(kappa),ylab="Ribsomal class sd",main="Density Independent Initiation, median protein size, long half-life",cex.main=0.85,lwd=2,cex.lab=1.2,font=2)
# counter <- 2
# for(i in par.ranges[[2]][c(2,1)]){
#   lines(plot1.params[mu==i,kappa],U_sd_imax_median_plot1[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
#   counter <- counter + 1
# }
# dev.off()
#
# png("Fig1_SD_mu_range_logx.png",width = 5,height = 5,units = "in",res=300)
# plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
#      xlab=bquote(kappa),ylab="Ribsomal class sd",main="Density Independent Initiation, median protein size, long half-life",lwd=2,cex.lab=1.2,font=2,cex.main=0.85)
# counter <- 2
# for(i in par.ranges[[2]][c(2,1)]){
#   lines(log(plot1.params[mu==i,kappa],10),U_sd_imax_median_plot1[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
#   counter <- counter + 1
# }
# dev.off()

#Mean plot
png("Fig1_mean.png",width = 5,height = 5,units = "in",res=300)
cols <- cividis(3)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=0.85,font=2,
     main="Density Independent Initiation, median protein size, Long Half-life")
counter <- 2
points(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(0.04,0.115,0.365),y=c(1.5,20.5,38.5),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_mean_logx.png",width = 5,height = 5,units = "in",res=300)
plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),lwd=2,cex.main=0.85,font=2,
     main="Density Independent Initiation, median protein size, Long Half-life")
points(log(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(-2.6,-1.7,-1.1),y=c(4,21,39.5),labels = c("Low","Medium","High"),font=2)
dev.off()

#half life quantiles and mu ranges
quantile(HL_minutes,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
log(2)/(quantile(HL_minutes,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)*60)

# png("Fig1_mean_mu_range.png",width = 5,height = 5,units = "in",res=300)
# cols2 <- cividis(4)
# plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
#      las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
#      main="No interference, median protein size, Long Half-life",col=cols2)
# counter <- 2
# points(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
# #text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
# for(i in par.ranges[[2]][c(6,2,1)]){
#   lines(plot1.params[mu==i,kappa],U_mean_imax_median_plot1[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
#   points(plot1.params[mu==i,kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
#   counter <- counter + 1
# }
# legend("bottomright",fill=cols2,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
# dev.off()

png("Fig1_mean_mu_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=0.85,font=2,
     main="Density Independent Initiation, median protein size",lty=1)
counter <- 2
points(plot1.params[mu==par.ranges[[2]][10],kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,2,1)]){
  lines(plot1.params[mu==i,kappa],U_mean_imax_median_plot1[plot1.params[mu==i,index]],lty=counter,lwd=2)
  points(plot1.params[mu==i,kappa][peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
  counter <- counter + 1
}
legend("bottomright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

# png("Fig1_mean_mu_range_logx.png",width = 5,height = 5,units = "in",res=300)
# cols2 <- cividis(5)
# plot(log(plot1.params[mu==par.ranges[[2]][8],kappa],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
#      las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
#      main="No interference, median protein size, Half-life range",col=cols2[1])
# counter <- 2
# points(log(plot1.params[mu==par.ranges[[2]][8],kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]][peaks_mean_index],pch=20,cex=2,col=col)
# #text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
# for(i in par.ranges[[2]][c(6,4,2,1)]){
#   lines(log(plot1.params[mu==i,kappa],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
#   points(log(plot1.params[mu==i,kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
# #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
#   counter <- counter + 1
# }
# legend("topleft",fill = cols2,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
# dev.off()

png("Fig1_mean_mu_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[mu==par.ranges[[2]][8],kappa],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=0.85,font=2,
     main="Density Independent Initiation, median protein size",lty=1)
counter <- 2
points(log(plot1.params[mu==par.ranges[[2]][8],kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]][peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,2,1)]){
  lines(log(plot1.params[mu==i,kappa],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]],lty=counter,lwd=2)
  points(log(plot1.params[mu==i,kappa][peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()


#Mean and SD legend
png("Fig1_Ribosome_mu_range_legend.png",width = 2,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()



png("Fig1_Ribosome_loading_no_int_mu_range_low.png",width = 5,height = 5,units = "in",res=300)
col <- c("grey80","grey60","grey40","grey20")
xrange <- 1:11
plot(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[1]])[xrange],type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,ylim=c(0,80),
     main ="Density independent initiation, median protein size")#,col=col[1]
#points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[1]])[xrange],pch=20)
mean.ribo.load <- U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index][peaks_mean_index[1]]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[1]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points(mean.ribo.load,interpolation,pch=20)
lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation))
counter <- 2
for(i in c(6,2,1)){
  lines(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[1]])[xrange],lty=counter,lwd=2)#col=col[counter]
#  points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[1]])[xrange],pch=20)
  mean.ribo.load <- U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][i],index][peaks_mean_index[1]]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[1]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points(mean.ribo.load,interpolation,pch=20)
  lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=counter)
   counter <- counter+1
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
#text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
legend("topright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_mu_range_medium.png",width = 5,height = 5,units = "in",res=300)
col <- c("grey80","grey60","grey40","grey20")
xrange <- 1:42
plot(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[2]])[xrange],type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,ylim=c(0,20),
     main ="Density independent initiation, median protein size")#,col=col[1]
#points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[2]])[xrange],pch=20,cex=0.7)
mean.ribo.load <- U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index][peaks_mean_index[2]]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[2]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points(mean.ribo.load,interpolation,pch=20)
lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=counter)
counter <- 2
for(i in c(6,2,1)){
  lines(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[2]])[xrange],lty=counter,lwd=2)#col=col[counter]
#  points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[2]])[xrange],pch=20,cex=0.7)
  mean.ribo.load <- U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][i],index][peaks_mean_index[2]]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[2]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points(mean.ribo.load,interpolation,pch=20)
  lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
#text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
legend("topright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_mu_range_high.png",width = 5,height = 5,units = "in",res=300)
col <- c("grey80","grey60","grey40","grey20")
xrange <- 1:42
plot(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[3]])[xrange],type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,ylim=c(0,20),
     main ="Density independent initiation, median protein size")#,col=col[1]
#points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[2]])[xrange],pch=20,cex=0.7)
mean.ribo.load <- U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index][peaks_mean_index[3]]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[3]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points(mean.ribo.load,interpolation,pch=20)
lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(6,2,1)){
  lines(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[3]])[xrange],lty=counter,lwd=2)#col=col[counter]
  #  points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[2]])[xrange],pch=20,cex=0.7)
  mean.ribo.load <- U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][i],index][peaks_mean_index[3]]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][i],index]][,peaks_mean_index[3]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points(mean.ribo.load,interpolation,pch=20)
  lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
#text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)

dev.off()


png("Fig1_Ribosome_loading_no_int_99_mu.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,80),
     main ="Density independent initiation, median protein size, Long Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[i]]),col=col[i],lwd=2)
# points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_50_mu.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][6],index]][,peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,80),
     main ="Density independent initiation, median protein size, Median Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][6],index]][,peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][6],index]][,peaks_mean_index[i]]),col=col[i],lwd=2)
#  points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][6],index]][,peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_1_mu.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][2],index]][,peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,80),
     main ="Density independent initiation, median protein size, Short Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][2],index]][,peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][2],index]][,peaks_mean_index[i]]),col=col[i],lwd=2)
# points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][2],index]][,peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_min_mu.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][1],index]][,peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,80),
     main ="Density independent initiation, median protein size, Min Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][1],index]][,peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][1],index]][,peaks_mean_index[i]]),col=col[i],lwd=2)
#  points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][1],index]][,peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(3,10,36),y=c(25,10,4),labels = c("Low","Medium","High"),font=2)
dev.off()





# Figure 1 imax range--------------------------------------------------------
#imax range plot, compare to tau=1 range of params with imax 30
protein_quant
protein_quant/9
#imax will be 8,15,41,87,164 The 1%, 10%, 50%, 90%, and 99% quantiles of arabidopsis
imax.range <- c(8,15,41,87,164)
imax.percentiles <- c("1%","10%","50%","90%","99%")
names(imax.range) <- imax.percentiles
Unmarked_imax_range_plot1 <- lapply(imax.range,function(j) apply(plot1.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))
names(Unmarked_imax_range_plot1) <- imax.percentiles

#Calculate summary statistics
U_pct_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_percent)
U_mean_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_mean)
U_mode_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_mode)
U_sd_imax_range_plot1 <- lapply(Unmarked_imax_range_plot1,Ribo_sd)

U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]

#profile plots imax range
peak_positions <- c(0.01,0.10,0.5,0.90,0.99)

peaks_mean_imax_range <- lapply(1:length(imax.range),function(i) round(U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]]/imax.range[i],2))
peaks_mean_index_imax_range <- lapply(1:length(imax.range), function(i) unlist(lapply(peak_positions,function(j) which(peaks_mean_imax_range[[i]]==j)[1])))

#plot curves for figure 1 imax range example

# Imax range no int norm length -------------------------------------------------------

png("Fig1_Ribosome_loading_no_int_imax_range_low_initiation_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.main=0.85,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density independent initiation, Low initiation, long half-life",lty=1)
print(U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1

}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_imax_range_medium_initiation_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density independent initiation, medium initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_imax_range_high_initiation_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density independent initiation, high initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1

}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topleft",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_imax_range_low_initiation_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density independent initiation, low initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_imax_range_medium_initiation_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density independent initiation, medium initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig1_Ribosome_loading_no_int_imax_range_high_initiation_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density independent initiation, high initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()



# bw intiation plots imax range ------------------------------------------------------

png("Fig1_mean_imax_range_logx_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/41,type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, median protein size",lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/41,pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(1,5)){
  lines(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[i]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[i],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[i],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",legend = c("Short (70 a.a. ~1%)","Median (370 a.a. ~50%)","Long (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig1_mean_imax_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,164),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, median protein size",lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(1,5)){
  lines(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[i]][plot1.params[(mu==par.ranges[[2]][10]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",legend = c("Short (70 a.a. ~1%)","Median (370 a.a. ~50%)","Long (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig1_mean_imax_range_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/41,type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, median protein size",lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/41,pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(1,5)){
  lines(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[i]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[i],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[i],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("right",legend = c("Short (70 a.a. ~1%)","Median (370 a.a. ~50%)","Long (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig1_mean_imax_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,164),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, median protein size",lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(1,5)){
  lines(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[i]][plot1.params[(mu==par.ranges[[2]][10]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("right",legend = c("Short (70 a.a. ~1%)","Median (370 a.a. ~50%)","Long (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()


#Mean and SD legend
png("Fig1_Ribosome_mean_sd_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))
dev.off()


#imax range and mu range

png("Fig1_mean_imax_long_mu_range_logx_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[5]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[5],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, Long protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[5]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[5],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[5]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[5],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[5]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[5],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_mean_imax_short_mu_range_logx_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[1]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[1],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, short protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[1]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[1],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[1]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[1],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[1]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[1],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_mean_imax_medium_mu_range_logx_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[3],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, medium protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[3],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[3],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[3],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

####
png("Fig1_mean_imax_long_mu_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[5]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,164),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, Long protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[5]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[5]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[5]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_mean_imax_short_mu_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[1]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,10),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, short protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[1]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[1]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[1]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_mean_imax_medium_mu_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, medium protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

#############

png("Fig1_mean_imax_long_mu_range_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[5]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[5],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, Long protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[5]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[5],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[5]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[5],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[5]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[5],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("right",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_mean_imax_short_mu_range_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[1]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[1],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, short protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[1]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[1],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[1]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[1],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[1]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[1],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("right",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_mean_imax_medium_mu_range_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[3],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, medium protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[3],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[3],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[3],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("bottomright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

####
png("Fig1_mean_imax_long_mu_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[5]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,164),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, Long protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[5]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[5]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[5]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("right",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_mean_imax_short_mu_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[1]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,10),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, short protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[1]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[1]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[1]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("bottomright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig1_mean_imax_medium_mu_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density independent initiation, medium protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1[[3]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index])
  counter <- counter + 1
}
legend("bottomright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()


#
#
# png("Fig1_Ribosome_loading_no_int_imax_range_kappa=0.05.png",width = 5,height = 5,units = "in",res=300)
# plot((0:imax.range[4]),Ribo_percent(Unmarked_imax_range_plot1[[4]][,501]),type="l",
#      ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,ylim=c(0,0.3),xlim=c(0,105),
#      main = bquote("No interference, imax=30, "*kappa*"="*.(round(plot1.params.imaxrange[501,kappa],2))*", "*tau*"=1"*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),cex.main=0.95)
# counter <- 2
# for(i in c(1:3,5:8)){
#   lines((0:imax.range[i]),Ribo_percent(Unmarked_imax_range_plot1[[i]][,501]),col=counter,lwd=2)
#   counter <- counter+1
# }
# dev.off()
#
#
#
# #SD plot imax range
# png("Fig1_SD_imax_range.png",width = 5,height = 5,units = "in",res=300)
# plot(plot1.params.imaxrange[,kappa],U_sd_imax_median_plot1[plot1.params[tau==0.1,index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1)),
#      xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
# counter <- 2
# counter <- 2
# for(i in par.ranges[[4]][1:4]){
#   lines(plot1.params[tau==i,kappa],U_sd_imax_median_plot1[plot1.params[tau==i,index]],lty=counter,lwd=2)
#   counter <- counter + 1
# }
# counter <- 1
# for(i in par.ranges[[4]][5:8]){
#   lines(plot1.params[tau==i,kappa],U_sd_imax_median_plot1[plot1.params[tau==i,index]],lty=counter,col="grey50",lwd=2)
#   counter <- counter + 1
# }
# dev.off()
#
# #Mean plot
# png("Fig1_mean_imax_range_tau=1.png",width = 5,height = 5,units = "in",res=300)
# plot(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[4]][1:2500]/imax.range[4])*100,type="l",ylim=c(0,100),
#      las=1,ylab="Mean Ribsome Load Percent",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
#      main=bquote("No interference, imax range, "*tau*"="*1*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
# counter <- 2
# for(i in c(1:3)){
#   lines(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[i]][1:2500]/imax.range[i])*100,lty=counter,lwd=2)
#   counter <- counter + 1
# }
# counter <- 1
# for(i in 5:8){
#   lines(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[i]][1:2500]/imax.range[i])*100,col="grey50",lty=counter,lwd=2)
#   counter <- counter + 1
# }
# dev.off()
#
# png("Fig1_mean_imax_range_tau=1.png",width = 5,height = 5,units = "in",res=300)
# plot(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[4]][1:2500]/imax.range[4])*100,type="l",ylim=c(0,100),
#      las=1,ylab="Mean Ribsome Load Percent",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
#      main=bquote("No interference, imax range, "*tau*"="*1*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
# counter <- 2
# for(i in c(1:3,5:8)){
#   lines(plot1.params.imaxrange[1:2500,kappa],(U_mean_imax_range_plot1[[i]][1:2500]/imax.range[i])*100,col=counter,lwd=2)
#   counter <- counter + 1
# }
# dev.off()
#
#
# plot(log(plot1.params.imaxrange[500:2500,kappa],10),(U_mean_imax_range_plot1[[4]][500:2500]/imax.range[4])*100,type="l",ylim=c(0,100),
#      las=1,ylab="Mean Ribsome Load Percent",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
#      main=bquote("No interference, imax range, "*tau*"="*1*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
# counter <- 2
# for(i in c(1:3,5:8)){
#   lines(log(plot1.params.imaxrange[500:2500,kappa],10),(U_mean_imax_range_plot1[[i]][500:2500]/imax.range[i])*100,col=counter,lwd=2)
#   counter <- counter + 1
# }
#
# png("Fig1_Ribosome_loading_no_int_imax_range_means.png",width = 5,height = 5,units = "in",res=300)
# plot(plot1.params.imaxrange[,kappa],U_mean_imax_range_plot1[[4]],type="l",ylim=c(0,181),
#      las=1,ylab="Mean Ribsome Load Percent",xlab= bquote(kappa),lwd=2,cex.lab=1.2,font=2,
#      main=bquote("No interference, imax range, "*tau*"="*1*", "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
# counter <- 2
# for(i in c(1:3,5:8)){
#   lines(plot1.params.imaxrange[,kappa],U_mean_imax_range_plot1[[i]],col=counter,lwd=2)
#   counter <- counter + 1
# }
# dev.off()
#
# #Mean and SD legend
# png("Fig1_Ribosome_mean_sd_legend.png",width = 1.5,height = 2,units = "in",res=300)
# par(mar=c(0.1,0.1,0.1,0.1))
# plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
# legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))
# dev.off()
#
#


setwd(old.dir)


# Figure 2 ----------------------------------------------------------------
#
# Figure 2 w/int ---------------------------------------------------------


#Define parameter set. Since determinant method is being used I chose a very small mu. Prior exploration showed that a
#small mu behaves much like mu=0. The reason I am doing this is that calculating the adjugate matrix is slower than finding
#the full solution.
#mu values Hypothetical fast 1e-2  1% 9.994377e-04 5% 5.068461e-04 10% 3.500904e-04 25% 1.973310e-04 50% 1.081415e-04
#75% 5.500915e-05 90% 2.702400e-05 95%  1.763636e-05  99% 8.043137e-06
par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.ranges <- list(c(seq(0.0001,0.5,0.0001)),c(0.01,1.004561*10^-3,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.004561*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,7.701635*10^-6) ,c(1),c(1),c(100000),c(1))
plot1.params <- CreateParMatrix(parameter.range = par.ranges,fixed = F,par.names = par.names)
plot1.params <- as.data.table(plot1.params)
plot1.params[,"index":=1:nrow(plot1.params)]
#calculate unmarked class
Unmarked_imax_median_plot1_int <- apply(plot1.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = 41,rescale = T,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

#Calculate summary statistics
U_pct_imax_median_plot1_int <- Ribo_percent(Unmarked_imax_median_plot1_int)
U_mean_imax_median_plot1_int <- Ribo_mean(Unmarked_imax_median_plot1_int)
U_sd_imax_median_plot1_int <- Ribo_sd(Unmarked_imax_median_plot1_int)
U_mode_imax_median_plot1_int <- Ribo_mode(Unmarked_imax_median_plot1_int)

#plot curves for figure 1 example
int_peaks_mean <- round(U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]])
int_peaks_mode <- U_mode_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]]
#int_peaks_mean_index <- unlist(lapply(c(0,5,10,15,25,30), function(i) which(int_peaks_mean==i)[1]))
int_peaks_mean_index <- which(unlist(plot1.params[(mu==par.ranges[[2]][10]),kappa])%in%Initiation_rates)
#int_peaks_mode_index <- unlist(lapply(c(0,5,10,15,25,30), function(i) which(int_peaks_mode==i)[1]))
#int_peaks_mean_index[1] <- 10
#Initiation_rates <- unlist(plot1.params[int_peaks_mean_index,1])
#names(Initiation_rates) <- c("Low","Mid","High")

old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Fig2_Ribosome_loading_w_int.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85 ,col=col[1],ylim=c(0,40),
     main ="Density Dependent Initiation, median protein size, Long Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(int_peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[i]]),col=col[i],lwd=2)
  #points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 14,31),y=c(25,16.5,16.5),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_Ribosome_loading_comp.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.8,col=col[1],ylim=c(0,70),
     main ="Density dependent vs independent, median protein size, Long Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(int_peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[i]]),col=col[i],lwd=2)
#  points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[i]]),col=col[i],pch=20)
}
for(i in 1:length(peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[i]]),col=col[i],lwd=2,lty=2)
#  points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
#text(x = c(5,14,25),y=c(25,15,15),labels = c("Low","Medium","High"),font=2)
legend("top",legend = c("Low","Medium","High","Density Dependent","Density Independent"),fill = c(col,rep(rgb(0,0,0,0),2)),lty = c(1,1,1,1,2),col = c(rep(rgb(0,0,0,0),3),rep("black",2)),border = F,)
abline(h=0,col="black",lwd=2)
dev.off()


png("Fig2_Ribosome_loading_no_int_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
kappas <- c(round(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index[1]],4),round(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index[2:6]],3))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
text(x = 0.5,y=0.8,"Scaled Initiation",font=2)
legend("center", legend = c("Low","Medium","High"),fill=col)
dev.off()

#SD plots
png("Fig2_SD_int.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_sd_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1_int)),
     xlab=bquote(kappa),ylab="Ribsomal class sd",main="No interference, median protein size, Long Half-life",lwd=2,cex.lab=1.2,font=2,cex.main=1.1)
points(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],U_sd_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(0.04,0.09,0.355),y=c(1.2,3.75,3.4),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_SD_int_logx.png",width = 5,height = 5,units = "in",res=300)
plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_sd_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1_int)),
     xlab=bquote("log"[10]*(kappa)),ylab="Ribsomal class sd",main="No interference, median protein size, Long Half-life",lwd=2,cex.lab=1.2,font=2,cex.main=1.1)
points(log(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],10),U_sd_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(-2.55,-1.5,-0.6),y=c(1.85,3.75,3.5),labels = c("Low","Medium","High"),font=2)
dev.off()

cols2 <- cividis(3)

# png("Fig2_SD_mu_range_int.png",width = 5,height = 5,units = "in",res=300)
# plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_sd_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1_int)),
#      xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("Density Dependent initiation, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
# counter <- 2
# for(i in par.ranges[[2]][c(2,1)]){
#   lines(plot1.params[mu==i,kappa],U_sd_imax_median_plot1_int[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
#   counter <- counter + 1
# }
# dev.off()

# png("Fig2_SD_mu_range_int_logx.png",width = 5,height = 5,units = "in",res=300)
# plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_sd_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_sd_imax_median_plot1_int)),
#      xlab=bquote(kappa),ylab="Ribsomal class sd",main=bquote("Density Dependent Initiation, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5),lwd=2,cex.lab=1.2,font=2)
# counter <- 2
# for(i in par.ranges[[2]][c(2,1)]){
#   lines(log(plot1.params[mu==i,kappa],10),U_sd_imax_median_plot1_int[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
#   counter <- counter + 1
# }
# dev.off()

#Mean plot
png("Fig2_mean_int.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1_int)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=0.85,font=2,
     main="Density Dependent, median protein size, Long Half-life")
counter <- 2
points(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(0.0475,0.115,.3),y=c(1.5,13.65,32),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_mean_comp.png",width = 5,height = 5,units = "in",res=300)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=0.8,font=2,col="grey50",
     main="Density Dependent vs Independent, median protein size, Long Half-life")
counter <- 2
lines(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],lwd=2)
points(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
points(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.0475,0.115,.3),y=c(1.5,13.65,32),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_mean_int_logx.png",width = 5,height = 5,units = "in",res=300)
plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote(log[10]*" Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=0.85,font=2,
     main="Density Dependent, median protein size, Long Half-life")
points(log(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
text(x = c(-2.5,-1.7,-0.8),y=c(3.5,14,30.75),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_mean_comp_logx.png",width = 5,height = 5,units = "in",res=300)
plot(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote(log[10]*" Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=0.8,font=2,
     main="Density Dependent vs Independent, median protein size, Long Half-life",col="grey50")
lines(log(plot1.params[mu==par.ranges[[2]][10],kappa],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1_int)),lwd=2)
points(log(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
points(log(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(-2.5,-1.7,-0.8),y=c(3.5,14,30.75),labels = c("Low","Medium","High"),font=2)
dev.off()


#half life quantiles and mu ranges
quantile(HL_minutes,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
log(2)/(quantile(HL_minutes,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)*60)

# png("Fig2_mean_mu_range.png",width = 5,height = 5,units = "in",res=300)
# cols2 <- cividis(4)
# plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1_int)),
#      las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
#      main="No interference, median protein size, Long Half-life",col=cols2)
# counter <- 2
# points(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
# #text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
# for(i in par.ranges[[2]][c(6,2,1)]){
#   lines(plot1.params[mu==i,kappa],U_mean_imax_median_plot1_int[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
#   points(plot1.params[mu==i,kappa][int_peaks_mean_index],U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index],pch=20,cex=2,col=col)
#   counter <- counter + 1
# }
# legend("bottomright",fill=cols2,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
# dev.off()

png("Fig2_mean_mu_range_int_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Long Half-life",lty=1)
counter <- 2
points(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,2,1)]){
  lines(plot1.params[mu==i,kappa],U_mean_imax_median_plot1_int[plot1.params[mu==i,index]],lty=counter,lwd=2)
  points(plot1.params[mu==i,kappa][int_peaks_mean_index],U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index],pch=20,cex=2,col=col)
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,436 (99%)","107    (50%)","11.5   (  1%)","1.15"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()


png("Fig2_mean_mu_range_comp_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,70),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Long Half-life",lty=1)
counter <- 2
points(plot1.params[mu==par.ranges[[2]][10],kappa][int_peaks_mean_index],U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
counter <- 1
for(i in par.ranges[[2]][c(6,2,1)]){
  lines(plot1.params[mu==i,kappa],U_mean_imax_median_plot1_int[plot1.params[mu==i,index]],lty=counter,lwd=2)
  points(plot1.params[mu==i,kappa][int_peaks_mean_index],U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index],pch=20,cex=2,col=col)
  counter <- counter + 1
}

for(i in par.ranges[[2]][c(10,6,2,1)]){
  lines(plot1.params[mu==i,kappa],U_mean_imax_median_plot1[plot1.params[mu==i,index]],lty=counter,lwd=2,col="grey50")
  points(plot1.params[mu==i,kappa][int_peaks_mean_index],U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,436 (99%)","107    (50%)","11.5   (  1%)","1.15"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

mu_range_indvsdep_comp <- 100*c( U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][6],index]][int_peaks_mean_index]/U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index]
                                ,U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][6],index]][peaks_mean_index]/U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index]
                                ,U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][2],index]][int_peaks_mean_index]/U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index]
                                ,U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][2],index]][peaks_mean_index]/U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index]
                                ,U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][1],index]][int_peaks_mean_index]/U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index]][int_peaks_mean_index]
                                ,U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][1],index]][peaks_mean_index]/U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][peaks_mean_index])

mu_range_indvsdep_comp <- matrix(mu_range_indvsdep_comp,nrow=3)

colnames(mu_range_indvsdep_comp) <- c("Median Dep","Median Ind","1% Dep","1% Ind","Min Dep","Min Ind")
rownames(mu_range_indvsdep_comp) <- c("Low","Mid","High")

png("Fig2_mean_mu_range_comp_barplot_grouped.png",width = 5,height = 5,units = "in",res=300)
bar.tmp <- barplot(mu_range_indvsdep_comp,beside = T,las=2,ylim=c(0,100),ylab="Percent Ribosomal load of Long Half-life",
        legend.text = c("Low","Mid","High"),main="Effect of Dependence and marking rate on ribosomal load",
        cex.main=0.85,font=2,cex.lab=1.2,names.arg = rep("",6))
text(x = colMeans(bar.tmp),y=-3,labels =  c("Median Dep","Median Ind","1% Dep","1% Ind","Min Dep","Min Ind"),font = 2,srt=45,xpd=T,cex=1,adj = 1)
box()
dev.off()

mu_range_indvsdep_comp_rearrange <- rbind(mu_range_indvsdep_comp[,1:2],mu_range_indvsdep_comp[,3:4],mu_range_indvsdep_comp[,5:6])

png("Fig2_mean_mu_range_comp_barplot.png",width = 5,height = 5,units = "in",res=300)
bar.tmp <-barplot(t(mu_range_indvsdep_comp_rearrange),beside = T,ylim=c(0,100),ylab="Percent Ribosomal load of Long Half-life",
        legend.text = c("Dependent","Independent"),names.arg =rep("",9),args.legend = c(cex=0.8),
        main="Effect of Dependence and marking rate on ribosomal load",cex.main=0.85,font=2,cex.lab=1.2)
text(x = colMeans(bar.tmp),y=-3,labels =  c("Median Low","Median Mid","Median High","1% Low","1% Mid","1% High","Min Low","Min Mid","Min High"),font = 2,srt=45,xpd=T,cex=1,adj = 1)
box()
dev.off()

# png("Fig2_mean_mu_range_logx.png",width = 5,height = 5,units = "in",res=300)
# cols2 <- cividis(5)
# plot(log(plot1.params[mu==par.ranges[[2]][8],kappa],10),U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][8],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1_int)),
#      las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
#      main="No interference, median protein size, Half-life range",col=cols2[1])
# counter <- 2
# points(log(plot1.params[mu==par.ranges[[2]][8],kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][8],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
# #text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
# for(i in par.ranges[[2]][c(6,4,2,1)]){
#   lines(log(plot1.params[mu==i,kappa],10),U_mean_imax_median_plot1_int[plot1.params[mu==i,index]],col=cols2[counter],lwd=2)
#   points(log(plot1.params[mu==i,kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index],pch=20,cex=2,col=col)
# #  abline(h=U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
#   counter <- counter + 1
# }
# legend("topleft",fill = cols2,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
# dev.off()

png("Fig2_mean_mu_range_int_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[mu==par.ranges[[2]][8],kappa],10),U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][8],index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Half-life range",lty=1)
counter <- 2
points(log(plot1.params[mu==par.ranges[[2]][8],kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][8],index]][int_peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,2,1)]){
  lines(log(plot1.params[mu==i,kappa],10),U_mean_imax_median_plot1_int[plot1.params[mu==i,index]],lty=counter,lwd=2)
  points(log(plot1.params[mu==i,kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index],pch=20,cex=2,col=col)
  #  abline(h=U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
dev.off()

png("Fig2_mean_mu_range_comp_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[mu==par.ranges[[2]][8],kappa],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="No interference, median protein size, Half-life range",lty=1,col="gray50")
counter <- 2
points(log(plot1.params[mu==par.ranges[[2]][8],kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][8],index]][peaks_mean_index],pch=20,cex=2,col=col)
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in par.ranges[[2]][c(6,2,1)]){
  lines(log(plot1.params[mu==i,kappa],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]],lty=counter,lwd=2,col="gray50")
  points(log(plot1.params[mu==i,kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1[plot1.params[mu==i,index]][peaks_mean_index],pch=20,cex=2,col=col)
  #  abline(h=U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
for(i in par.ranges[[2]][c(10,6,2,1)]){
  lines(log(plot1.params[mu==i,kappa],10),U_mean_imax_median_plot1_int[plot1.params[mu==i,index]],lty=counter,lwd=2)
  points(log(plot1.params[mu==i,kappa][int_peaks_mean_index],10),U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index],pch=20,cex=2,col=col)
  #  abline(h=U_mean_imax_median_plot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}

legend("topleft",lty = 1:5,legend = c("1,436 (99%)","107 (50%)","33 (10%)","11.5 (1%)","1.15"),title = "Half-life (min.)",xjust = 1)
dev.off()


#Mean and SD legend
png("Fig2_Ribosome_mean_sd_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))
dev.off()



png("Fig2_Ribosome_loading_int_mu_range_low.png",width = 5,height = 5,units = "in",res=300)
col <- c("grey80","grey60","grey40","grey20")
xrange <- 1:11
plot(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]])[xrange],type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,ylim=c(0,80),
     main ="Density dependent low intitiation, median protein size")#,col=col[1]
#points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]])[xrange],pch=20)
mean.ribo.load <- U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index][int_peaks_mean_index[1]]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points(mean.ribo.load,interpolation,pch=20)
lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation))
counter <- 2
for(i in c(6,2,1)){
  lines(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[1]])[xrange],lty=counter,lwd=2)#col=col[counter]
  #  points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[1]])[xrange],pch=20)
  mean.ribo.load <- U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][i],index][int_peaks_mean_index[1]]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[1]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points(mean.ribo.load,interpolation,pch=20)
  lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
#text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_Ribosome_loading_int_mu_range_medium.png",width = 5,height = 5,units = "in",res=300)
col <- c("grey80","grey60","grey40","grey20")
xrange <- 1:42
plot(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[2]])[xrange],type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,ylim=c(0,20),
     main ="Density Dependent medium initiation, median protein size")#,col=col[1]
#points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[2]])[xrange],pch=20,cex=0.7)
mean.ribo.load <- U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index][int_peaks_mean_index[2]]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[2]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points(mean.ribo.load,interpolation,pch=20)
lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=counter)
counter <- 2
for(i in c(6,2,1)){
  lines(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[2]])[xrange],lty=counter,lwd=2)#col=col[counter]
  #  points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[2]])[xrange],pch=20,cex=0.7)
  mean.ribo.load <- U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][i],index][int_peaks_mean_index[2]]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[2]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points(mean.ribo.load,interpolation,pch=20)
  lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
#text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_Ribosome_loading_int_mu_range_high.png",width = 5,height = 5,units = "in",res=300)
col <- c("grey80","grey60","grey40","grey20")
xrange <- 1:42
plot(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[3]])[xrange],type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,ylim=c(0,20),
     main ="Density Dependent high initiation, median protein size")#,col=col[1]
#points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[2]])[xrange],pch=20,cex=0.7)
mean.ribo.load <- U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][10],index][int_peaks_mean_index[3]]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[3]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points(mean.ribo.load,interpolation,pch=20)
lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(6,2,1)){
  lines(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[3]])[xrange],lty=counter,lwd=2)#col=col[counter]
  #  points(xrange-1,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[2]])[xrange],pch=20,cex=0.7)
  mean.ribo.load <- U_mean_imax_median_plot1_int[plot1.params[mu==par.ranges[[2]][i],index][int_peaks_mean_index[3]]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][i],index]][,int_peaks_mean_index[3]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points(mean.ribo.load,interpolation,pch=20)
  lines(x=rep(mean.ribo.load,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
#text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()


png("Fig2_Ribosome_loading_int_99_mu.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,80),
     main ="Density Dependent, median protein size, Long Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(int_peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[i]]),col=col[i],lwd=2)
 # points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]][,int_peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_Ribosome_loading_int_50_mu.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][6],index]][,int_peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,80),
     main ="Density Dependent, median protein size, median Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][6],index]][,int_peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(int_peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][6],index]][,int_peaks_mean_index[i]]),col=col[i],lwd=2)
#  points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][6],index]][,int_peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_Ribosome_loading_int_1_mu.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][2],index]][,int_peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,80),
     main ="Density Dependent, median protein size, short Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][2],index]][,int_peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(int_peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][2],index]][,int_peaks_mean_index[i]]),col=col[i],lwd=2)
#  points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][2],index]][,int_peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(5, 20,36),y=c(25,11,25),labels = c("Low","Medium","High"),font=2)
dev.off()

png("Fig2_Ribosome_loading_int_min_mu.png",width = 5,height = 5,units = "in",res=300)
col <- cividis(3)
plot(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][1],index]][,int_peaks_mean_index[1]]),type="l",
     ylab="",xlab="",lwd=2,font=2,cex.main=0.85,col=col[1],ylim=c(0,80),
     main ="Density Dependent, median protein size, min Half-life")
#points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][1],index]][,int_peaks_mean_index[1]]),pch=20,col=col[1])
for(i in 2:length(int_peaks_mean_index)){
  lines(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][1],index]][,int_peaks_mean_index[i]]),col=col[i],lwd=2)
#  points(0:41,100*Ribo_percent(Unmarked_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][1],index]][,int_peaks_mean_index[i]]),col=col[i],pch=20)
}
mtext(1,text = "Ribosomal Class",line=2.5,font=2)
mtext(2,text = "Percent Ribosomal Load",line=2.5,font=2)
text(x = c(3,10,36),y=c(25,10,4),labels = c("Low","Medium","High"),font=2)
dev.off()


# Figure 2 imax range--------------------------------------------------------
#imax range plot, compare to tau=1 range of params with imax 30
protein_quant
protein_quant/9
#imax will be 8,15,41,87,164 The 1%, 10%, 50%, 90%, and 99% quantiles of arabidopsis
imax.range <- c(8,15,41,87,164)
imax.percentiles <- c("1%","10%","50%","90%","99%")
names(imax.range) <- imax.percentiles
Unmarked_imax_range_plot1_int <- lapply(imax.range,function(j) apply(plot1.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau"))))
names(Unmarked_imax_range_plot1_int) <- imax.percentiles

#Calculate summary statistics
U_pct_imax_range_plot1_int <- lapply(Unmarked_imax_range_plot1_int,Ribo_percent)
U_mean_imax_range_plot1_int <- lapply(Unmarked_imax_range_plot1_int,Ribo_mean)
U_mode_imax_range_plot1_int <- lapply(Unmarked_imax_range_plot1_int,Ribo_mode)
U_sd_imax_range_plot1_int <- lapply(Unmarked_imax_range_plot1_int,Ribo_sd)

U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]

#profile plots imax range
int_peak_positions <- c(0.01,0.10,0.5,0.90,0.99)

int_peaks_mean_imax_range <- lapply(1:length(imax.range),function(i) round(U_mean_imax_range_plot1_int[[i]][plot1.params[mu==par.ranges[[2]][10],index]]/imax.range[i],2))
int_peaks_mean_index_imax_range <- lapply(1:length(imax.range), function(i) unlist(lapply(int_peak_positions,function(j) which(int_peaks_mean_imax_range[[i]]==j)[1])))

#plot curves for figure 2 imax range example

# Fig2 Imax range int norm length -------------------------------------------------------

png("Fig2_Ribosome_loading_int_imax_range_low_initiation_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.main=0.85,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, Low initiation, long half-life",lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1

}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()


png("Fig2_Ribosome_loading_int_imax_range_low_comp_initiation_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.main=0.85,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, Low initiation, long half-life",lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1

}
counter <- 1
ltype <- c(2,1,3)
for(i in c(1,3,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),lwd=2,lty=ltype[counter],col="gray50")
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=ltype[counter],col="gray50")
  counter <- counter+1
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)","Dependent","Independent"),
       lty=c(2,1,3,1,1),col=c(rep("black",3),rep(rgb(0,0,0,0),2)),lwd=2,fill=c(rep(rgb(0,0,0,0),3),"black","gray50"),border=NA)

dev.off()

png("Fig2_Ribosome_loading_int_imax_range_medium_initiation_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, medium initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()


png("Fig2_Ribosome_loading_int_imax_range_medium_comp_initiation_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, medium initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1
}
counter <- 1
ltype <- c(2,1,3)
for(i in c(1,3,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),col="gray50",lwd=2,lty=ltype[counter])
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=ltype[counter],col="gray50")
  counter <- counter+1
}

#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)","Dependent","Independent"),
       lty=c(2,1,3,1,1),col=c(rep("black",3),rep(rgb(0,0,0,0),2)),lwd=2,fill=c(rep(rgb(0,0,0,0),3),"black","gray50"),border=NA)
dev.off()

png("Fig2_Ribosome_loading_int_imax_range_high_initiation_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, high initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1

}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topleft",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig2_Ribosome_loading_int_imax_range_high_initiation_comp_longHL_lengthnorm.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3])/imax.range[3]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Number of Occupied Sites (%)",lwd=2,font=2,cex.lab=1.2,ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, high initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[3]]))
mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)

counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1_int[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),lty=counter)
  counter <- counter+1

}
counter <- 1
ltype=c(2,1,3)
for(i in c(1,3,5)){
  lines((0:imax.range[i])/imax.range[i]*100,100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),col="gray50",lwd=2,lty=ltype[counter])
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]/nrow(Unmarked_imax_range_plot1[[i]]))
  mean.ribo.load <- U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  points((mean.ribo.load/imax.range[i])*100,interpolation,pch=20)
  lines(x=rep((mean.ribo.load/imax.range[i])*100,2),y=c(-5,interpolation),col="gray50",lty=ltype[counter])
  counter <- counter+1

}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topleft",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)","Dependent","Independent"),
       lty=c(2,1,3,1,1),col=c(rep("black",3),rep(rgb(0,0,0,0),2)),lwd=2,fill=c(rep(rgb(0,0,0,0),3),"black","gray50"),border=NA)
dev.off()


png("Fig2_Ribosome_loading_int_imax_range_low_initiation_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, low initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig2_Ribosome_loading_int_imax_range_low_initiation_comp_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, low initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
counter <- 1
ltype=c(2,1,3)
for(i in c(1,3,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]]),col="gray50",lwd=2,lty=ltype[counter])
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[1])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)","Dependent","Independent"),
       lty=c(2,1,3,1,1),col=c(rep("black",3),rep(rgb(0,0,0,0),2)),lwd=2,fill=c(rep(rgb(0,0,0,0),3),"black","gray50"),border=NA)
dev.off()

png("Fig2_Ribosome_loading_int_imax_range_medium_initiation_comp_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, medium initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
counter <- 1
ltype=c(2,1,3)
for(i in c(1,3,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),col="gray50",lwd=2,lty=ltype[counter])
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)","Dependent","Independent"),
       lty=c(2,1,3,1,1),col=c(rep("black",3),rep(rgb(0,0,0,0),2)),lwd=2,fill=c(rep(rgb(0,0,0,0),3),"black","gray50"),border=NA)
dev.off()

png("Fig2_Ribosome_loading_int_imax_range_medium_initiation_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, medium initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[2])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig2_Ribosome_loading_int_imax_range_high_initiation_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, high initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()


png("Fig2_Ribosome_loading_int_imax_range_high_initiation_comp_longHL.png",width = 5,height = 5,units = "in",res=300)
col3 = rep("black",3)
plot((0:imax.range[3]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),type="l",
     ylab="Percent Ribosomal Load",xlab="Ribosomal Class",lwd=2,font=2,cex.lab=1.2,xlim=c(0,164),ylim=c(0,80),col=col3[3],
     main = "Density dependent initiation, high initiation, long half-life",cex.main=0.85,lty=1)
print(U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])
# mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
# mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
# interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
# points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
# lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
counter <- 2
for(i in c(1,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1_int[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),col=col3[counter],lwd=2,lty=counter)
  print(U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
counter <- 1
ltype=c(2,1,3)
for(i in c(1,3,5)){
  lines((0:imax.range[i]),100*Ribo_percent(Unmarked_imax_range_plot1[[i]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]),col="gray50",lwd=2,lty=ltype[counter])
  print(U_mean_imax_range_plot1[[i]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]],lty=counter)
  counter <- counter+1
  # mean.ribo.load <- U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]]
  # mean.class.load <- 100*Ribo_percent(Unmarked_imax_range_plot1_int[[3]][,plot1.params[(kappa==Initiation_rates[3])&(mu==par.ranges[[2]][10]),index]])[(floor(mean.ribo.load)+1):(ceiling(mean.ribo.load)+1)]
  # interpolation <- mean.class.load[1]+(mean.ribo.load-floor(mean.ribo.load))*(mean.class.load[2]-mean.class.load[1])
  # points((mean.ribo.load/imax.range[3])*100,interpolation,pch=20)
  # lines(x=rep((mean.ribo.load/imax.range[3])*100,2),y=c(-5,interpolation),lty=1)
}
#abline(v=peak_positions[3],col="grey50",lty=2)
legend("topright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)","Dependent","Independent"),
       lty=c(2,1,3,1,1),col=c(rep("black",3),rep(rgb(0,0,0,0),2)),lwd=2,fill=c(rep(rgb(0,0,0,0),3),"black","gray50"),border=NA)
dev.off()



#Fig 2  bw intiation plots imax range ------------------------------------------------------

png("Fig2_mean_imax_range_logx_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/41,type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, median protein size",lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/41,pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(1,5)){
  lines(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[i]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[i],lty=counter,lwd=2,col="green")
  # for(j in 1:3){
  #   points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[i],pch=20,cex=2,col=col[j])
  # }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",legend = c("Short (70 a.a. ~1%)","Median (370 a.a. ~50%)","Long (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig2_mean_imax_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,164),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, median protein size",lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(1,5)){
  lines(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[i]][plot1.params[(mu==par.ranges[[2]][10]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",legend = c("Short (70 a.a. ~1%)","Median (370 a.a. ~50%)","Long (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig2_mean_imax_range_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/41,type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, median protein size",lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/41,pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(1,5)){
  lines(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[i]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[i],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[i],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("right",legend = c("Short (70 a.a. ~1%)","Median (370 a.a. ~50%)","Long (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()

png("Fig2_mean_imax_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,164),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, median protein size",lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(1,5)){
  lines(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[i]][plot1.params[(mu==par.ranges[[2]][10]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[i]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("right",legend = c("Short (70 a.a. ~1%)","Median (370 a.a. ~50%)","Long (1470 a.a. ~99%)"),lty=c(2,1,3),lwd=2)
dev.off()


#Mean and SD legend
png("Fig2_Ribosome_mean_sd_legend.png",width = 1.5,height = 2,units = "in",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))
dev.off()


#imax range and mu range

png("Fig2_mean_imax_long_mu_range_logx_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[5]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[5],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, Long protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[5]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[5],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[5]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[5],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[5]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[5],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig2_mean_imax_short_mu_range_logx_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[1]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[1],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, short protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[1]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[1],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[1]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[1],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[1]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[1],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig2_mean_imax_medium_mu_range_logx_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[3],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, medium protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[3],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[3],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[3],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

####
png("Fig2_mean_imax_long_mu_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[5]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,164),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, Long protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[5]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[5]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[5]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig2_mean_imax_short_mu_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[1]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,10),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, short protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[1]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[1]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[1]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig2_mean_imax_medium_mu_range_logx_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(log(plot1.params[(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, medium protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(log(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(log(plot1.params[(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(log(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],10),U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("topleft",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

#############

png("Fig2_mean_imax_long_mu_range_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[5]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[5],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, Long protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[5]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[5],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[5]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[5],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[5]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[5],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("right",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig2_mean_imax_short_mu_range_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[1]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[1],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, short protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[1]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[1],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[1]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[1],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[1]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[1],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("right",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig2_mean_imax_medium_mu_range_percent_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]]*100/imax.range[3],type="l",ylim=c(0,100),
     las=1,ylab="Percent Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, medium protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]]*100/imax.range[3],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][i]),index]]*100/imax.range[3],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]]*100/imax.range[3],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("bottomright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

####
png("Fig2_mean_imax_long_mu_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[5]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,164),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, Long protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[5]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[5]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[5]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("right",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig2_mean_imax_short_mu_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[1]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,10),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, short protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[1]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[1]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[1]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("bottomright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()

png("Fig2_mean_imax_medium_mu_range_bw.png",width = 5,height = 5,units = "in",res=300)
cols2 <- cividis(5)
plot(plot1.params[(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][10]),index]],type="l",ylim=c(0,45),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("log"[10]*" [ Scaled Initiation ("*kappa/tau*") ]"),lwd=2,cex.main=1,font=2,
     main="Density dependent initiation, medium protein size",cex.main=0.85,lty=1)
counter <- 2
for(i in 1:3){
  points(plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[i])&(mu==par.ranges[[2]][10]),index]],pch=20,cex=2,col=col[i])
}
#text(x = c(0.04,0.1,0.12),y=c(2,14.5,29),labels = c("Low","Medium","High"),font=2)
for(i in c(6,2,1)){
  lines(plot1.params[(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(mu==par.ranges[[2]][i]),index]],lty=counter,lwd=2)
  for(j in 1:3){
    points(plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),kappa],U_mean_imax_range_plot1_int[[3]][plot1.params[(kappa==Initiation_rates[j])&(mu==par.ranges[[2]][i]),index]],pch=20,cex=2,col=col[j])
  }
  #  abline(h=U_mean_imax_medianplot1_int[plot1.params[mu==i,index]][int_peaks_mean_index])
  counter <- counter + 1
}
legend("bottomright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
       title = "Half-life (min.)",xjust = 1,lwd=2)
dev.off()


setwd(old.dir)



cub1 <- function(x){
    return(x+2*x^2+3*x^3+4*x^4+5*x^5+6*x^6)
}
cub1(1)


cub2 <- function(x){
  return(1+x+x^2+x^3+x^4+x^5+x^6)
}


plot(seq(0.001,0.5,0.001),cub1(seq(0.001,0.5,0.001)),type="l")
plot(seq(0.001,0.5,0.001),cub2(seq(0.001,0.5,0.001)),type="l")
plot(seq(0.001,0.5,0.001),cub1(seq(0.001,0.5,0.001))/cub2(seq(0.001,0.5,0.001)),type="l")
abline(a=0,b=1,col="red")

plot(seq(0.001,1,0.001),cub1(seq(0.001,1,0.001))/cub1(seq(0.001,1,0.001)),type="l")



pm_generate <- function(parameters,imax,rescale = T ,rescale.val = "tau",
                        par.names=c("kappa","mu","lambda","tau","delta","s2")){
  if(is.null(names(parameters))){
    names(parameters) <- par.names
  }
  #load parameters
  tau <- parameters["tau"]
  kappa <- parameters["kappa"]
  lambda <-  parameters["lambda"]
  mu <- parameters["mu"]
  #add interference
  if(rescale==T){
    tau <- tau/parameters[rescale.val]
    kappa <- kappa/parameters[rescale.val]
    lambda <- lambda/parameters[rescale.val]
    mu <- mu/parameters[rescale.val]
  }

  pm_vec <- numeric(imax+1)
  for(i in seq(0,imax,1)){
    if(i==0){
      pm_vec[i+1] <- prod(1:imax)*tau^imax
    }else{

      pm_vec[i+1] <- prod(1:imax)/prod(1:i)*(tau^(imax-i))*kappa^(i)*prod(seq((imax+1-i),imax,1))/(imax^i)

    }

  }
  return(lambda*pm_vec)

}

plot(plot1.params[mu==par.ranges[[2]][10],kappa],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]],type="l",ylim=c(0,max(U_mean_imax_median_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote("Scaled Initiation ("*kappa/tau*")"),lwd=2,cex.main=0.85,font=2,
     main="Density Independent Initiation, median protein size, Long Half-life")

plot1.params
pm_vec <- pm_generate(unlist(plot1.params[mu==par.ranges[[2]][10],][1,]),imax=41)

m<- 100*pm_vec/sum(pm_vec)
Ribo_percent(Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index][1]])
m

0:41%*%m
0:41%*%Unmarked_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index][1]]/100

plot(plot1.params[mu==par.ranges[[2]][10],kappa][1:750],U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                     "y"=U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)

fit.tmp2 <-lm(formula = y~x-1,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                             "y"=U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp2,col="blue",lty=3,lwd=2)

plot1.params[mu==par.ranges[[2]][10],kappa][1:750]*367.5


plot(plot1.params[mu==par.ranges[[2]][10],kappa][1:750],U_mean_imax_range_plot1[[2]][plot1.params[mu==par.ranges[[2]][10],index]][1:750],type="l",lwd=2)
fit.tmp2 <-lm(formula = y~x-1,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                                "y"=U_mean_imax_range_plot1[[2]][plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp2,col="blue",lty=3,lwd=2)

library(truncnorm)

U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][1:750]

plot.ind <- 1500
tn <- truncnorm::dtruncnorm(seq(0,41,1),a = 0,b = 41,mean = 367.7*plot1.params[mu==par.ranges[[2]][10],kappa][plot.ind],U_sd_imax_median_plot1[plot1.params[mu==par.ranges[[2]][10],index]][plot.ind])
plot(seq(0,41,1),tn/sum(tn),type="l",lwd=2,ylim=c(0,1.2))
lines(seq(0,41,1),U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index][plot.ind]],col="red",lty=2,lwd=2)


#Is the slope of the linear part of the mean ribosome load plot related to # of codons

#Imax 8 1%
plot(plot1.params[mu==par.ranges[[2]][10],kappa][1:750],U_mean_imax_range_plot1[[1]][plot1.params[mu==par.ranges[[2]][10],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[1]][plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)

#Imax 15 10%
plot(plot1.params[mu==par.ranges[[2]][10],kappa][1:750],U_mean_imax_range_plot1[[2]][plot1.params[mu==par.ranges[[2]][10],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[2]][plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)

#Imax 41 50%
plot(plot1.params[mu==par.ranges[[2]][10],kappa][1:750],U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][10],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)


#Imax 87 90%
plot(plot1.params[mu==par.ranges[[2]][10],kappa][1:750],U_mean_imax_range_plot1[[4]][plot1.params[mu==par.ranges[[2]][10],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[4]][plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)

#Imax 164 99%
plot(plot1.params[mu==par.ranges[[2]][10],kappa][1:750],U_mean_imax_range_plot1[[5]][plot1.params[mu==par.ranges[[2]][10],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[5]][plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)


#marking rate effect on slope
#Imax 41 50% mu is small 99% 1500 min entry 10
plot(plot1.params[mu==par.ranges[[2]][10],kappa][1:750],U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][10],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][10],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][10],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)

#slope is 367.5

#Imax 41 50% mu is median 50% 115 min entry 6
plot(plot1.params[mu==par.ranges[[2]][6],kappa][1:750],U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][6],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][6],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][6],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)
#slope is 355.4

#Imax 41 50% mu is large 1% 11.5 min entry 2
plot(plot1.params[mu==par.ranges[[2]][2],kappa][1:750],U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][2],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][2],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][2],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)
#slope is 269.1

#Imax 41 50% mu is very large sim 1.15 min entry 1
plot(plot1.params[mu==par.ranges[[2]][1],kappa][1:750],U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][1],index]][1:750],type="l",lwd=2)
fit.tmp <-lm(formula = y~x,data = data.frame("x"=plot1.params[mu==par.ranges[[2]][1],kappa][1:750],
                                             "y"=U_mean_imax_range_plot1[[3]][plot1.params[mu==par.ranges[[2]][1],index]][1:750]))
abline(fit.tmp,col="red",lty=2,lwd=2)

#slope is 78.678

plot(par.ranges[[2]][c(10,6,2,1)],1-abs((c(367.5,355.4,269.1,78.7)/367.5)))

plot(par.ranges[[2]][c(10,6,2,1)],(par.ranges[[2]][c(10,6,2,1)])^0.35)


fact_exp <- numeric(30)
for(i in 1:30){
  fact_exp[i] <- prod(1:i)/i^i
}

plot(1:30,fact_exp)


x1<-c(seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(2,100,1))
unm<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.0000000000001,1,1,1,1),6,rescale = F,interference = T,which.interfere = c("tau")))

plot(x1,lapply(unm,Ribo_mean))


calc_adjugate <- function(parameters,imax,density_dependent = F,rescale =F,rescale.val,
                          par.names=c("kappa","mu","lambda","tau","delta","s2")){

  if(is.null(names(parameters))){
    names(parameters) <- par.names
  }
  #load parameters
  tau <- parameters["tau"]
  tau <- rep(tau/(imax*9),imax+1)
  kappa <- rep(parameters["kappa"],imax+1)
  lambda <-  parameters["lambda"]
  mu <- rep(parameters["mu"],imax+1)

  if(rescale==T){
    tau <- tau/parameters[rescale.val]
    kappa <- kappa/parameters[rescale.val]
    lambda <- lambda/parameters[rescale.val]
    mu <- mu/parameters[rescale.val]
  }

  if(density_dependent==F){
    adjugate <- numeric(imax+1)
    for(i in 0:imax){
      if(i==0){
        adjugate[i+1] <- prod(1:imax) * tau[i+1]^imax

      }else if (i<imax&&i>0){
        adjugate[i+1] <- prod((i+1):imax) * tau[i+1]^(imax-i) * kappa[i+1]^(i)

      }else{
        adjugate[i+1] <- kappa[i+1]^(i)

      }
    }

  }else {
    adjugate <- numeric(imax+1)
    for(i in 0:imax){
      if(i==0){
        adjugate[i+1] <- prod(1:imax) * tau[i+1]^imax

      }else if (i<imax&&i>0){
        adjugate[i+1] <- prod((i+1):imax) * tau[i+1]^(imax-i) * prod(1-(0:(i-1))/imax) * kappa[i+1]^(i)

      }else{
        adjugate[i+1] <- prod(1-(0:(i-1))/imax) * kappa[i+1]^(i)

      }
    }

  }

  return(adjugate)

}


Adj <- calc_adjugate(c(x1[1],0.0000000000001,1,1,1,1),imax = 6,density_dependent = F)

Adj <- lapply(x1,function(i) calc_adjugate(c(i,0.0000000000001,1,1,1,1),imax = 6,density_dependent = F))
#Adj/sum(Adj)*100

A_pct <- lapply(Adj,Ribo_percent)
A_ribomean <- lapply(A_pct,Ribo_mean)

plot(x1[1:131],unlist(lapply(unm,Ribo_mean))[1:131])
lines(x1[1:131],unlist(A_ribomean[1:131]),col="red")
plot(x1[1:131],A_ribomean[1:131],col="red",type="l")

unm1<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.000001,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("kappa","tau")))
unm2<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.001,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("kappa","tau")))
unm3<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.01,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("kappa","tau")))
Adj <- lapply(x1,function(i) calc_adjugate(c(i,0.000001,1,1,1,1),imax = 41,density_dependent = T))
Adj2 <- lapply(x1,function(i) calc_adjugate(c(i,0.000001,1,1,1,1),imax = 41,density_dependent = F))
prod_sol <- unlist(lapply(0:41,function(i) prod(1-(0:i)/41)))

approx_density <- (41*9*x1)*(1)
# approx_density <- approx_density[1:70]*c(rep(prod_sol[1],5),rep(prod_sol[2],5),rep(prod_sol[3],5),rep(prod_sol[4],5),
#                               rep(prod_sol[5],5),rep(prod_sol[6],5),rep(prod_sol[7],5),rep(prod_sol[8],5),
#                               rep(prod_sol[9],5),rep(prod_sol[11],5),rep(prod_sol[11],5),rep(prod_sol[12],5),
#                               rep(prod_sol[13],5),rep(prod_sol[4],5))
plot(prod_sol)
plot(log(x1[1:41],10),(unlist(lapply(Adj,Ribo_mean))/unlist(lapply(Adj2,Ribo_mean)))[1:41])

plot(x1[1:41],unlist(lapply(unm1,Ribo_mean))[1:41])
lines(x1[1:41],unlist(lapply(Adj,Ribo_mean))[1:41])
lines(x1[1:41],approx_density[1:41])

#



unm1<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.000001,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("tau")))
unm2<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.001,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("tau")))
unm3<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.01,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("tau")))
Adj <- lapply(x1,function(i) calc_adjugate(c(i,0.000001,1,1,1,1),imax = 41,density_dependent = F))

A_pct <- lapply(Adj,Ribo_percent)
A_ribomean <- lapply(A_pct,Ribo_mean)

A_ribomean_0.001 <- unlist(lapply(unlist(A_ribomean)/x1,function(i) i/(1+i*0.001)))*x1
A_ribomean_0.01 <- unlist(lapply(unlist(A_ribomean)/x1,function(i) i/(1+i*0.01)))*x1

A_ribomean_0.001v2 <- unlist(lapply(1:length(x1),function(i) unlist(A_ribomean)[i]/x1[i]/(1+unlist(A_ribomean)[i]/x1[i]*0.001/exp(x1[i]/2))))*x1
A_ribomean_0.01v2 <- unlist(lapply(1:length(x1),function(i) unlist(A_ribomean)[i]/x1[i]/(1+unlist(A_ribomean)[i]/x1[i]*0.01/exp(x1[i]/2))))*x1

plot(x1,(unlist(A_ribomean)/x1)/369)
plot(x1,(unlist(A_ribomean)/x1)*0.001)
plot(x1,(unlist(A_ribomean)/x1)*0.001^(unlist(A_ribomean)/x1/369))
plot(x1[1:81],unlist(lapply(unm2,Ribo_mean))[1:81]-unlist(A_ribomean_0.001)[1:81],pch=20,cex=0.8)
plot(x1[1:81],unlist(lapply(unm3,Ribo_mean))[1:81]-unlist(A_ribomean_0.01)[1:81],pch=20,cex=0.8)

plot(log(x1[1:81],10),unlist(lapply(unm2,Ribo_mean))[1:81]-unlist(A_ribomean_0.001)[1:81],pch=20,cex=0.8,ylim=c(-2,2))
lines(log(x1[1:81],10),unlist(lapply(unm2,Ribo_mean))[1:81]-unlist(A_ribomean_0.001v2)[1:81],col="red")

plot(log(x1[1:81],10),unlist(lapply(unm3,Ribo_mean))[1:81]-unlist(A_ribomean_0.01)[1:81],pch=20,cex=0.8,ylim=c(-2,5))
lines(log(x1[1:81],10),unlist(lapply(unm3,Ribo_mean))[1:81]-unlist(A_ribomean_0.01v2)[1:81],col="red")


lines(x1,unlist(lapply(unm2,Ribo_mean))-unlist(A_ribomean_0.001))
plot(x1,unlist(lapply(unm2,Ribo_mean))-unlist(A_ribomean_0.001),type='l')

plot(x1,(0.001^(1+unlist(A_ribomean)/x1)))

plot(x1,1-(x1/unlist(A_ribomean))/max(x1/unlist(A_ribomean)))

plot(x1,unlist(A_ribomean)/x1)

plot(x1[1:41],unlist(lapply(unm1,Ribo_mean))[1:41],cex=0.7)
lines(x1[1:41],unlist(A_ribomean[1:41]),col="red")
points(x1[1:41],unlist(lapply(unm2,Ribo_mean))[1:41],pch=3,cex=0.7)
lines(x1[1:41],unlist(A_ribomean_0.001[1:41]),col="blue")
lines(x1[1:41],unlist(A_ribomean_0.001v2[1:41]),col="darkblue",lty=2)
points(x1[1:41],unlist(lapply(unm3,Ribo_mean))[1:41],pch=6,cex=0.7)
lines(x1[1:41],unlist(A_ribomean_0.01[1:41]),col="green")
lines(x1[1:41],unlist(A_ribomean_0.01v2[1:41]),col="darkgreen",lty=2)

plot(x1[1:131],A_ribomean[1:131],col="red",type="l")


plot(x1,unlist(lapply(unm2,Ribo_mean))-A_ribomean_0.001)
plot(x1,unlist(lapply(unm3,Ribo_mean))-A_ribomean_0.01)

plot(x1[1:111],(unlist(lapply(unm2,Ribo_mean))-A_ribomean_0.001)[1:111])
plot(log(x1[1:111],10),(unlist(lapply(unm2,Ribo_mean))-A_ribomean_0.001)[1:111])
plot(x1[1:111],(unlist(lapply(unm3,Ribo_mean))-A_ribomean_0.01)[1:111])
plot(log(x1[1:111],10),(unlist(lapply(unm3,Ribo_mean))-A_ribomean_0.01)[1:111])


Adj <- calc_adjugate(c(0.0001,0.0000000000001,1,1,1,1),imax = 6,density_dependent = F)
Adj/sum(Adj)
(Adj/sum(Adj))%*%0:6

Adj <- calc_adjugate(c(1,0.0000000000001,1,1,1,1),imax = 6,density_dependent = F)
Adj/sum(Adj)
(Adj/sum(Adj))%*%0:6

A_scaled <- Adj*0:6

A_scaled[2]/Adj[1]
A_scaled[3]/Adj[2]
A_scaled[4]/Adj[3]
A_scaled[5]/Adj[4]
A_scaled[6]/Adj[5]
A_scaled[7]/Adj[6]


54^6/720 + 54^5/120 + 54^4/24 + 54^3/6 + 54^2/2 + 54 + 1
54^6/720 + 54^5/24  + 54^4/6 + 54^3/2 + 54^2 + 54

(54^6/120 + 54^5/24  + 54^4/6 + 54^3/2 + 54^2 + 54)/(54^6/720 + 54^5/120 + 54^4/24 + 54^3/6 + 54^2/2 + 54 + 1)


720/121

calc_adjugate_w_mu <- function(parameters,imax,density_dependent = F,rescale =F,rescale.val,
                          par.names=c("kappa","mu","lambda","tau","delta","s2")){

  if(is.null(names(parameters))){
    names(parameters) <- par.names
  }
  #load parameters
  tau <- parameters["tau"]
  tau <- rep(tau/(imax*9),imax+1)
  kappa <- rep(parameters["kappa"],imax+1)
  lambda <-  parameters["lambda"]
  mu <- rep(parameters["mu"],imax+1)

  if(rescale==T){
    tau <- tau/parameters[rescale.val]
    kappa <- kappa/parameters[rescale.val]
    lambda <- lambda/parameters[rescale.val]
    mu <- mu/parameters[rescale.val]
  }

  if(density_dependent==F){
    adjugate <- numeric(imax+1)
    for(i in 0:imax){
      if(i==0){
        adjugate[i+1] <- prod(1:imax) * tau[i+1]^imax + mu[1]

      }else if (i<imax&&i>0){
        adjugate[i+1] <- prod((i+1):imax) * tau[i+1]^(imax-i) * kappa[i+1]^(i) + mu[1]*kappa[i+1]^(i)

      }else{
        adjugate[i+1] <- kappa[i+1]^(i) + mu[1]*kappa[i+1]^(i)

      }
    }

  }else {
    adjugate <- numeric(imax+1)
    for(i in 0:imax){
      if(i==0){
        adjugate[i+1] <- prod(1:imax) * tau[i+1]^imax

      }else if (i<imax&&i>0){
        adjugate[i+1] <- prod((i+1):imax) * tau[i+1]^(imax-i) * prod(1-(0:(i-1))/imax) * kappa[i+1]^(i)

      }else{
        adjugate[i+1] <- prod(1-(0:(i-1))/imax) * kappa[i+1]^(i)

      }
    }

  }

  return(adjugate)

}



unm1<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.000001,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("tau")))
unm2<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.001,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("tau")))
unm3<- lapply(x1,function(i) CalcUnmarkedClass(c(i,0.01,1,1,1,1),41,rescale = F,interference = T,which.interfere = c("tau")))
Adj1 <- lapply(x1,function(i) calc_adjugate_w_mu(c(i,0.000001,1,1,1,1),imax = 41,density_dependent = F))
Adj2 <- lapply(x1,function(i) calc_adjugate_w_mu(c(i,0.001,1,1,1,1),imax = 41,density_dependent = F))
Adj3 <- lapply(x1,function(i) calc_adjugate_w_mu(c(i,0.01,1,1,1,1),imax = 41,density_dependent = F))

A_pct1 <- lapply(Adj1,Ribo_percent)
A_ribomean1 <- lapply(A_pct1,Ribo_mean)

A_pct2 <- lapply(Adj2,Ribo_percent)
A_ribomean2 <- lapply(A_pct2,Ribo_mean)

A_pct3 <- lapply(Adj3,Ribo_percent)
A_ribomean3 <- lapply(A_pct3,Ribo_mean)


plot(x1[1:111],unlist(lapply(unm1,Ribo_mean))[1:111],cex=0.7)
lines(x1[1:111],unlist(A_ribomean1[1:111]),col="red")
points(x1[1:111],unlist(lapply(unm2,Ribo_mean))[1:111],pch=3,cex=0.7)
lines(x1[1:111],unlist(A_ribomean2[1:111]),col="blue")
points(x1[1:111],unlist(lapply(unm3,Ribo_mean))[1:111],pch=6,cex=0.7)
lines(x1[1:111],unlist(A_ribomean3[1:111]),col="green")





Adj <- lapply(x1,function(i) calc_adjugate(c(i,0.000001,1,1,1,1),imax = 41,density_dependent = F))

Adj_M <- matrix(unlist(Adj),nrow=length(Adj[[1]]))

Adj_M
Peaks_Adj_M <- matrix(NA,ncol=ncol(Adj_M),nrow=nrow(Adj_M))
for(i in 1:ncol(Adj_M)){
  Peaks_Adj_M[,i] <- order(Adj_M[,i],decreasing = F)
}

cols <- rainbow(7)
plot(0,0,axes=F,col=rgb(0,0,0,0),xlab="",ylab="",
     xlim=log(c(min(x1[1:40]),max(x1[1:40])),10),ylim=c(1,8),
     main="Density independent")
for(i in 1:nrow(Peaks_Adj_M)){
  lines(log(x1[1:40],10),Peaks_Adj_M[i,1:40],col=cols[i],lwd=2)
}
axis(1,at = log(c(0.0001,0.001,0.01,0.1,1),10),labels = c(0.0001,0.001,0.01,0.1,1),las=2)
axis(2,at= 1:7,labels = 7:1,las=2)
legend("topleft",fill = cols,legend = c("p0","p1","p2","p3","p4","p5","p6"))
box()



cols <- rainbow(nrow(Peaks_Adj_M))
plot(0,0,axes=F,col=rgb(0,0,0,0),xlab="",ylab="",
     xlim=c(min(x1[1:31]),max(x1[1:31])),ylim=c(1,nrow(Peaks_Adj_M)+1),
     main="Density independent")
for(i in 1:nrow(Peaks_Adj_M)){
  lines(x1[1:31],Peaks_Adj_M[i,1:31],col=cols[i],lwd=2)
}
axis(1,at = c(0.0001,0.001,0.01,0.1,1),labels = c(0.0001,0.001,0.01,0.1,1),las=2)
axis(2,at= 1:nrow(Peaks_Adj_M),labels = (nrow(Peaks_Adj_M)):1,las=2)
#legend("topleft",fill = cols,legend = paste("p",0:(nrow(Peaks_Adj_M)-1),sep=""))
lines(x1[1:31],apply(Adj_M[,1:31],2,Ribo_mean)+1,lwd=2)
axis(4,at = 1:(nrow(Peaks_Adj_M)),labels = 0:(nrow(Peaks_Adj_M)-1))
box()



prod(1:6)/54^6

c(x1[1]^5+#green
    (1*(x1[1]^4)+6*(x1[1]^4))/54+#blue
    (2*(x1[1]^3)+6*(x1[1]^3)+30*(x1[1]^3))/54^2+#pink
    (6*(x1[1]^2)+12*(x1[1]^2)+30*(x1[1]^2)+120*(x1[1]^4))/54^3+ #red
    (24*(x1[1])+36*(x1[1])+ 60*(x1[1])+120*(x1[1])+360*(x1[1]))/54^4+#grey
    (prod(1:5)+prod(c(2,3,4,6))+prod(c(2,3,5,6))+prod(c(2,4,5,6))+prod(c(3,4,5,6))+prod(c(2,3,4,5,6)))/54^5 #yellow
    )

Ribo_mean(unm3[[1]])/x1[1]

54/(1+0.01*54)


# Result 1 ----------------------------------------------------------------

#Define parameter set. Since determinant method is being used I chose a very small mu. Prior exploration showed that a
#small mu behaves much like mu=0. The reason I am doing this is that calculating the adjugate matrix is slower than finding
#the full solution.
#mu values Hypothetical fast 1e-2  1% 9.994377e-04 5% 5.068461e-04 10% 3.500904e-04 25% 1.973310e-04 50% 1.081415e-04
#75% 5.500915e-05 90% 2.702400e-05 95%  1.763636e-05  99% 8.043137e-06
#Late March. Mike suggested keeping it to more round numbers ~99% 1500, ~50% 115, ~1% 11.5, sim 1.15
par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.ranges <- list(c(seq(0.0001,0.5,0.0001)),c(0.01,1.004561*10^-3,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.004561*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,7.701635*10^-6) ,c(1),c(1),c(100000),c(1))
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
peaks_mean_index <- unlist(lapply(c(2,21,41), function(i) which(peaks_mean==i)[1]))
#peaks_mode_index <- unlist(lapply(c(0,5,10,15,25,30), function(i) which(peaks_mode==i)[1]))
#peaks_mean_index[1] <- 10
Initiation_rates <- unlist(plot1.params[peaks_mean_index,1])
names(Initiation_rates) <- c("Low","Mid","High")


log((peaks_mean/41)/(1-(peaks_mean/41)),10)
min((1-peaks_mean)/41)

max((peaks_mean/41)/(1-(peaks_mean/41)))

max(peaks_mean)

min((peaks_mean/41))

(peaks_mean/41)



for(i in 1:length(par.ranges[[2]])){
  lods <- log((U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][i],index]]/41)/(1-(U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][i],index]]/41)))
  if(i==1){
    plot((U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][i],index]])/41,lods,xlim=c(0,1),ylim=c(-8,8))
  }else {
    lines((U_mean_imax_median_plot1[plot1.params[mu==par.ranges[[2]][i],index]])/41,lods,col=i)
  }
}

for(i in 1:5){
  lods <- log((U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]]/imax.range[i])/(1-(U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]]/imax.range[i])))
  if(i==1){
    plot((U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]])/imax.range[i],lods,xlim=c(0,1),ylim=c(-8,8))
  }else {
    lines((U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]])/imax.range[i],lods,col=i)
  }
}


for(i in 1:5){
  lods <- log((U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]]/imax.range[i])/(1-(U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]]/imax.range[i])))
  if(i==1){
    plot((U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]])/imax.range[i],lods,xlim=c(0,1),ylim=c(-8,8))
  }else {
    lines((U_mean_imax_range_plot1[[i]][plot1.params[mu==par.ranges[[2]][10],index]])/imax.range[i],lods,col=i)
  }
}


U_mean_imax_range_plot1
par.ranges






percent_var <- function(initial.val,percentages=1){
  output <- numeric(length(percentages)+1)
  for(i in 1: length(output)){
    if(i == 1){
      output[i] <- initial.val
    }else {
      output[i] <- initial.val*percentages[i-1]
    }
  }
  return(sort(unique(output),decreasing = F))
}

percentage_variances <-c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1,1.01,1.05,1.1,1.25,1.5,1.75,1.9,1.95,1.99)
percent_var(0.05,percentage_variances)


c(0.005,0.025,0.05,0.25,0.5)
par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.low.ranges <- list(c(percent_var(0.005,percentage_variances)),c(0.01,1.004561*10^-3,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.004561*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,7.701635*10^-6) ,c(1),c(1),c(100000),c(1))
plot.low.params <- CreateParMatrix(parameter.range = par.low.ranges,fixed = F,par.names = par.names)
plot.low.params <- as.data.table(plot.low.params)
plot.low.params[,"index":=1:nrow(plot.low.params)]

par.midlow.ranges <- list(c(percent_var(0.025,percentage_variances)),c(0.01,1.004561*10^-3,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.004561*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,7.701635*10^-6) ,c(1),c(1),c(100000),c(1))
plot.midlow.params <- CreateParMatrix(parameter.range = par.midlow.ranges,fixed = F,par.names = par.names)
plot.midlow.params <- as.data.table(plot.midlow.params)
plot.midlow.params[,"index":=1:nrow(plot.midlow.params)]

par.mid.ranges <- list(c(percent_var(0.05,percentage_variances)),c(0.01,1.004561*10^-3,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.004561*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,7.701635*10^-6) ,c(1),c(1),c(100000),c(1))
plot.mid.params <- CreateParMatrix(parameter.range = par.mid.ranges,fixed = F,par.names = par.names)
plot.mid.params <- as.data.table(plot.mid.params)
plot.mid.params[,"index":=1:nrow(plot.mid.params)]

par.midhigh.ranges <- list(c(percent_var(0.25,percentage_variances)),c(0.01,1.004561*10^-3,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.004561*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,7.701635*10^-6) ,c(1),c(1),c(100000),c(1))
plot.midhigh.params <- CreateParMatrix(parameter.range = par.midhigh.ranges,fixed = F,par.names = par.names)
plot.midhigh.params <- as.data.table(plot.midhigh.params)
plot.midhigh.params[,"index":=1:nrow(plot.midhigh.params)]

par.high.ranges <- list(c(percent_var(0.5,percentage_variances)),c(0.01,1.004561*10^-3,5.068461*10^-4,3.500904*10^-4,1.973310*10^-4,1.004561*10^-4,5.500915*10^-5,2.702400*10^-5,1.763636*10^-5,7.701635*10^-6) ,c(1),c(1),c(100000),c(1))
plot.high.params <- CreateParMatrix(parameter.range = par.high.ranges,fixed = F,par.names = par.names)
plot.high.params <- as.data.table(plot.high.params)
plot.high.params[,"index":=1:nrow(plot.high.params)]




imax.range <- c(8,15,41,87,164)
imax.percentiles <- c("1%","10%","50%","90%","99%")
names(imax.range) <- imax.percentiles
Unmarked_imax_range_plot_low <- lapply(imax.range,function(j) apply(plot.low.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))
names(Unmarked_imax_range_plot_low) <- imax.percentiles

Unmarked_imax_range_plot_midlow <- lapply(imax.range,function(j) apply(plot.midlow.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))
names(Unmarked_imax_range_plot_midlow) <- imax.percentiles

Unmarked_imax_range_plot_mid <- lapply(imax.range,function(j) apply(plot.mid.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))
names(Unmarked_imax_range_plot_mid) <- imax.percentiles

Unmarked_imax_range_plot_midhigh <- lapply(imax.range,function(j) apply(plot.midhigh.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))
names(Unmarked_imax_range_plot_midhigh) <- imax.percentiles

Unmarked_imax_range_plot_high <- lapply(imax.range,function(j) apply(plot.high.params,1,function(i) CalcUnmarkedClass(parameters = i,imax = j,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))
names(Unmarked_imax_range_plot_high) <- imax.percentiles



#Calculate summary statistics
U_pct_imax_range_plot_low <- lapply(Unmarked_imax_range_plot_low,Ribo_percent)
U_mean_imax_range_plot_low <- lapply(Unmarked_imax_range_plot_low,Ribo_mean)
U_mode_imax_range_plot_low <- lapply(Unmarked_imax_range_plot_low,Ribo_mode)
U_sd_imax_range_plot_low <- lapply(Unmarked_imax_range_plot_low,Ribo_sd)

U_pct_imax_range_plot_midlow <- lapply(Unmarked_imax_range_plot_midlow,Ribo_percent)
U_mean_imax_range_plot_midlow <- lapply(Unmarked_imax_range_plot_midlow,Ribo_mean)
U_mode_imax_range_plot_midlow <- lapply(Unmarked_imax_range_plot_midlow,Ribo_mode)
U_sd_imax_range_plot_midlow <- lapply(Unmarked_imax_range_plot_midlow,Ribo_sd)

U_pct_imax_range_plot_mid <- lapply(Unmarked_imax_range_plot_mid,Ribo_percent)
U_mean_imax_range_plot_mid <- lapply(Unmarked_imax_range_plot_mid,Ribo_mean)
U_mode_imax_range_plot_mid <- lapply(Unmarked_imax_range_plot_mid,Ribo_mode)
U_sd_imax_range_plot_mid <- lapply(Unmarked_imax_range_plot_mid,Ribo_sd)

U_pct_imax_range_plot_midhigh <- lapply(Unmarked_imax_range_plot_midhigh,Ribo_percent)
U_mean_imax_range_plot_midhigh <- lapply(Unmarked_imax_range_plot_midhigh,Ribo_mean)
U_mode_imax_range_plot_midhigh <- lapply(Unmarked_imax_range_plot_midhigh,Ribo_mode)
U_sd_imax_range_plot_midhigh <- lapply(Unmarked_imax_range_plot_midhigh,Ribo_sd)

U_pct_imax_range_plot_high <- lapply(Unmarked_imax_range_plot_high,Ribo_percent)
U_mean_imax_range_plot_high <- lapply(Unmarked_imax_range_plot_high,Ribo_mean)
U_mode_imax_range_plot_high <- lapply(Unmarked_imax_range_plot_high,Ribo_mode)
U_sd_imax_range_plot_high <- lapply(Unmarked_imax_range_plot_high,Ribo_sd)




low_vshort <- U_mean_imax_range_plot_low[[1]][plot.low.params[mu==par.low.ranges[[2]][10],index]]
low_short <- U_mean_imax_range_plot_low[[2]][plot.low.params[mu==par.low.ranges[[2]][10],index]]
low_median <- U_mean_imax_range_plot_low[[3]][plot.low.params[mu==par.low.ranges[[2]][10],index]]
low_long <- U_mean_imax_range_plot_low[[4]][plot.low.params[mu==par.low.ranges[[2]][10],index]]
low_vlong <- U_mean_imax_range_plot_low[[5]][plot.low.params[mu==par.low.ranges[[2]][10],index]]

midlow_vshort <- U_mean_imax_range_plot_midlow[[1]][plot.midlow.params[mu==par.midlow.ranges[[2]][10],index]]
midlow_short <- U_mean_imax_range_plot_midlow[[2]][plot.midlow.params[mu==par.midlow.ranges[[2]][10],index]]
midlow_median <- U_mean_imax_range_plot_midlow[[3]][plot.midlow.params[mu==par.midlow.ranges[[2]][10],index]]
midlow_long <- U_mean_imax_range_plot_midlow[[4]][plot.midlow.params[mu==par.midlow.ranges[[2]][10],index]]
midlow_vlong <- U_mean_imax_range_plot_midlow[[5]][plot.midlow.params[mu==par.midlow.ranges[[2]][10],index]]

mid_vshort <- U_mean_imax_range_plot_mid[[1]][plot.mid.params[mu==par.mid.ranges[[2]][10],index]]
mid_short <- U_mean_imax_range_plot_mid[[2]][plot.mid.params[mu==par.mid.ranges[[2]][10],index]]
mid_median <- U_mean_imax_range_plot_mid[[3]][plot.mid.params[mu==par.mid.ranges[[2]][10],index]]
mid_long <- U_mean_imax_range_plot_mid[[4]][plot.mid.params[mu==par.mid.ranges[[2]][10],index]]
mid_vlong <- U_mean_imax_range_plot_mid[[5]][plot.mid.params[mu==par.mid.ranges[[2]][10],index]]

midhigh_vshort <- U_mean_imax_range_plot_midhigh[[1]][plot.midhigh.params[mu==par.midhigh.ranges[[2]][10],index]]
midhigh_short <- U_mean_imax_range_plot_midhigh[[2]][plot.midhigh.params[mu==par.midhigh.ranges[[2]][10],index]]
midhigh_median <- U_mean_imax_range_plot_midhigh[[3]][plot.midhigh.params[mu==par.midhigh.ranges[[2]][10],index]]
midhigh_long <- U_mean_imax_range_plot_midhigh[[4]][plot.midhigh.params[mu==par.midhigh.ranges[[2]][10],index]]
midhigh_vlong <- U_mean_imax_range_plot_midhigh[[5]][plot.midhigh.params[mu==par.midhigh.ranges[[2]][10],index]]

high_vshort <- U_mean_imax_range_plot_high[[1]][plot.high.params[mu==par.high.ranges[[2]][10],index]]
high_short <- U_mean_imax_range_plot_high[[2]][plot.high.params[mu==par.high.ranges[[2]][10],index]]
high_median <- U_mean_imax_range_plot_high[[3]][plot.high.params[mu==par.high.ranges[[2]][10],index]]
high_long <- U_mean_imax_range_plot_high[[4]][plot.high.params[mu==par.high.ranges[[2]][10],index]]
high_vlong <- U_mean_imax_range_plot_high[[5]][plot.high.params[mu==par.high.ranges[[2]][10],index]]



jitter(rep(1,19))

plot(jitter(rep(1,19)),low_vshort-low_vshort[10],xlim=c(0,6),ylim=c(-10,10),pch=20)
points(jitter(rep(2,19)),low_short-low_short[10],pch=20)
points(jitter(rep(3,19)),low_median-low_median[10],pch=20)
points(jitter(rep(4,19)),low_long-low_long[10],pch=20)
points(jitter(rep(5,19)),low_vlong-low_vlong[10],pch=20)

cols <- cividis(5)
plot(percentage_variances,low_vshort-low_vshort[10],ylim=c(-10,10),pch = 20, xlab = bquote("Percent change in "*kappa/tau), ylab = bquote(Delta*" Mean Ribosomal Load"),col=cols[1])
points(percentage_variances,low_short-low_short[10],pch=20,col=cols[2])
points(percentage_variances,low_median-low_median[10],pch=20,col=cols[3])
points(percentage_variances,low_long-low_long[10],pch=20,col=cols[4])
points(percentage_variances,low_vlong-low_vlong[10],pch=20,col=cols[5])

png("Variation_in_Ribosomal_loading_low_initiation_imax_range.png",width = 6,height = 6,units = "in",res = 300)
plot(percentage_variances,low_vshort-low_vshort[10],ylim=c(-10,10),
     pch = 20, xlab = bquote("Percent change in "*kappa/tau),
     ylab = bquote(Delta*" Mean Ribosomal Load"),col=cols[1],
     main="Variation in ribosomal loading, D.I.I, low initiation",cex=1.4)
lines(percentage_variances,low_vshort-low_vshort[10],lwd=2)
# points(percentage_variances,low_short-low_short[10],pch=20,col=cols[2])
# lines(percentage_variances,low_short-low_short[10],pch=20,col=cols[2])
points(percentage_variances,low_median-low_median[10],pch=20,col=cols[3],cex=1.4)
lines(percentage_variances,low_median-low_median[10],pch=20,col=cols[3],lwd=2)
# points(percentage_variances,low_long-low_long[10],pch=20,col=cols[4])
# lines(percentage_variances,low_long-low_long[10],pch=20,col=cols[4])
points(percentage_variances,low_vlong-low_vlong[10],pch=20,col=cols[5],cex=1.4)
lines(percentage_variances,low_vlong-low_vlong[10],pch=20,col=cols[5],lwd=2)
legend("topleft",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),fill = cols[c(1,3,5)])
dev.off()


png("Variation_in_Ribosomal_loading_midlow_initiation_imax_range.png",width = 6,height = 6,units = "in",res = 300)
plot(percentage_variances,midlow_vshort-midlow_vshort[10],ylim=c(-20,20),
     pch = 20, xlab = bquote("Percent change in "*kappa/tau),
     ylab = bquote(Delta*" Mean Ribosomal Load"),col=cols[1],
     main="Variation in ribosomal loading, D.I.I, midlow initiation",cex=1.4)
lines(percentage_variances,midlow_vshort-midlow_vshort[10],lwd=2)
# points(percentage_variances,midlow_short-midlow_short[10],pch=20,col=cols[2])
# lines(percentage_variances,midlow_short-midlow_short[10],pch=20,col=cols[2])
points(percentage_variances,midlow_median-midlow_median[10],pch=20,col=cols[3],cex=1.4)
lines(percentage_variances,midlow_median-midlow_median[10],pch=20,col=cols[3],lwd=2)
# points(percentage_variances,midlow_long-midlow_long[10],pch=20,col=cols[4])
# lines(percentage_variances,midlow_long-midlow_long[10],pch=20,col=cols[4])
points(percentage_variances,midlow_vlong-midlow_vlong[10],pch=20,col=cols[5],cex=1.4)
lines(percentage_variances,midlow_vlong-midlow_vlong[10],pch=20,col=cols[5],lwd=2)
legend("topleft",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),fill = cols[c(1,3,5)])
dev.off()



png("Variation_in_Ribosomal_loading_mid_initiation_imax_range.png",width = 6,height = 6,units = "in",res = 300)
plot(percentage_variances,mid_vshort-mid_vshort[10],ylim=c(-40,40),
     pch = 20, xlab = bquote("Percent change in "*kappa/tau),
     ylab = bquote(Delta*" Mean Ribosomal Load"),col=cols[1],
     main="Variation in ribosomal loading, D.I.I, mid initiation",cex=1.4)
lines(percentage_variances,mid_vshort-mid_vshort[10],lwd=2)
# points(percentage_variances,mid_short-mid_short[10],pch=20,col=cols[2])
# lines(percentage_variances,mid_short-mid_short[10],pch=20,col=cols[2])
points(percentage_variances,mid_median-mid_median[10],pch=20,col=cols[3],cex=1.4)
lines(percentage_variances,mid_median-mid_median[10],pch=20,col=cols[3],lwd=2)
# points(percentage_variances,mid_long-mid_long[10],pch=20,col=cols[4])
# lines(percentage_variances,mid_long-mid_long[10],pch=20,col=cols[4])
points(percentage_variances,mid_vlong-mid_vlong[10],pch=20,col=cols[5],cex=1.4)
lines(percentage_variances,mid_vlong-mid_vlong[10],pch=20,col=cols[5],lwd=2)
legend("topleft",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),fill = cols[c(1,3,5)])
dev.off()

png("Variation_in_Ribosomal_loading_midhigh_initiation_imax_range.png",width = 6,height = 6,units = "in",res = 300)
plot(percentage_variances,midhigh_vshort-midhigh_vshort[10],ylim=c(-80,10),
     pch = 20, xlab = bquote("Percent change in "*kappa/tau),
     ylab = bquote(Delta*" Mean Ribosomal Load"),col=cols[1],
     main="Variation in ribosomal loading, D.I.I, midhigh initiation",cex=1.4)
lines(percentage_variances,midhigh_vshort-midhigh_vshort[10],lwd=2)
# points(percentage_variances,midhigh_short-midhigh_short[10],pch=20,col=cols[2])
# lines(percentage_variances,midhigh_short-midhigh_short[10],pch=20,col=cols[2])
points(percentage_variances,midhigh_median-midhigh_median[10],pch=20,col=cols[3],cex=1.4)
lines(percentage_variances,midhigh_median-midhigh_median[10],pch=20,col=cols[3],lwd=2)
# points(percentage_variances,midhigh_long-midhigh_long[10],pch=20,col=cols[4])
# lines(percentage_variances,midhigh_long-midhigh_long[10],pch=20,col=cols[4])
points(percentage_variances,midhigh_vlong-midhigh_vlong[10],pch=20,col=cols[5],cex=1.4)
lines(percentage_variances,midhigh_vlong-midhigh_vlong[10],pch=20,col=cols[5],lwd=2)
legend("bottomright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),fill = cols[c(1,3,5)])
dev.off()

png("Variation_in_Ribosomal_loading_high_initiation_imax_range.png",width = 6,height = 6,units = "in",res = 300)
plot(percentage_variances,high_vshort-high_vshort[10],ylim=c(-80,10),
     pch = 20, xlab = bquote("Percent change in "*kappa/tau),
     ylab = bquote(Delta*" Mean Ribosomal Load"),col=cols[1],
     main="Variation in ribosomal loading, D.I.I, high initiation",cex=1.4)
lines(percentage_variances,high_vshort-high_vshort[10],lwd=2)
# points(percentage_variances,high_short-high_short[10],pch=20,col=cols[2])
# lines(percentage_variances,high_short-high_short[10],pch=20,col=cols[2])
points(percentage_variances,high_median-high_median[10],pch=20,col=cols[3],cex=1.4)
lines(percentage_variances,high_median-high_median[10],pch=20,col=cols[3],lwd=2)
# points(percentage_variances,high_long-high_long[10],pch=20,col=cols[4])
# lines(percentage_variances,high_long-high_long[10],pch=20,col=cols[4])
points(percentage_variances,high_vlong-high_vlong[10],pch=20,col=cols[5],cex=1.4)
lines(percentage_variances,high_vlong-high_vlong[10],pch=20,col=cols[5],lwd=2)
legend("bottomright",legend = c("Short    (70 a.a.     ~1%)","Median (370 a.a.   ~50%)","Long     (1470 a.a. ~99%)"),fill = cols[c(1,3,5)])
dev.off()


# Contour plot fig 1 ------------------------------------------------------
#DII
png(paste(Sys.Date(),"_contour_DII_medianlength_low_marking.png",sep=""),width = 6,height = 8,units = "in",res=300)
filled.contour2(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]])),
               col = c(grey.colors(n = 10,start = 1,end = 0,alpha = 1),rep(rgb(0,0,0,1),10)),nlevels = 20,plot.axes = c(axis(1,at = seq(0,40,5),labels = seq(0,40,5)),axis(2,at = seq(0,-5,-1),labels = seq(0,-5,-1))),
               key.axes = c(axis(4,at = seq(0,1,0.2),labels = seq(0,100,20))),key.title = "Ribosomal Density")
#contour.default(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]])))
mtext(1,text = "Ribosomal Class",line=2,at = 20,cex=1.2)
mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=2,cex=1.2)
#mtext(4,text = "Ribosomal Density",line=1,at = -2,cex=1.2)
dev.off()

plot(0,0,col=rgb(0,0,0,0),xlim=c(0,41),ylim=c(min(log(par.ranges[[1]],10)),max(log(par.ranges[[1]],10))))

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


#DDI
png(paste(Sys.Date(),"_contour_DDI_medianlength_low_marking.png",sep=""),width = 6,height = 8,units = "in",res=300)
filled.contour(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]])),
               color.palette = cividis,plot.axes = c(axis(1,at = seq(0,40,5),labels = seq(0,40,5)),axis(2,at = seq(0,-5,-1),labels = seq(0,-5,-1))),
               key.axes = c(axis(4,at = seq(0,1,0.2),labels = seq(0,100,20))),key.title = "Ribosomal Density")
mtext(1,text = "Ribosomal Class",line=2,at = 15,cex=1.2)
mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=2,cex=1.2)
mtext(4,text = "Ribosomal Density",line=1,at = -2,cex=1.2)
dev.off()




#Figure 1

png(paste(Sys.Date(),"_Figure1_DIIvsDDI_medianlength_low_marking.png",sep=""),width = 8,height = 6,units = "in",res=300)

layout(matrix(c(1,2,3,4,4,4,5,5,5,6,6,6,7,8,9),nrow = 3,ncol=5),widths = c(2,3,1,3,2),heights = c(1,1,1))

for(i in 3:1){
  if(i==1){
    par(mar=c(3.1,4.1,1.1,0.1))
  }else if(i==3){
    par(mar=c(1.1,4.1,1.1,0.1))
  }else{
    par(mar=c(1.1,4.1,1.1,0.1))
  }
  plot(0:41,U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],][kappa==Initiation_rates[i],index]],pch=20,xlab="",ylab="",ylim=c(0,0.7),axes=F)
  lines(0:41,U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],][kappa==Initiation_rates[i],index]],pch=20)
  axis(1,at = seq(0,40,10),labels = seq(0,40,10),las=2)
  axis(2,at = seq(0,0.7,0.1),labels = seq(0,0.7,0.1),las=2)
  box()
  if(i == 1) {
    mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
    text(x=10,y=0.6,"low",font=2,cex=1.2)
  } else if(i == 2){
    mtext(2,text = "Ribosomal Location Probability",line=2,cex=1,xpd=F)
    text(x=10,y=0.6,"mid",font=2,cex=1.2)
  } else {
    text(x=10,y=0.6,"high",font=2,cex=1.2)
  }

}

#plot(0,0,col=rgb(0,0,0,0),axes=F,xlim=c(0,41),ylim=c(min(log(par.ranges[[1]],10)),max(log(par.ranges[[1]],10))))
par(mar=c(3.1,3.1,1.1,0.1))

image(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]])),
      col = c(rgb(1,1,1,1),grey.colors(n = 12,start = 1,end = 0,alpha = 1)[4:12],rep(rgb(0,0,0,1),10)),axes=F,xlab="",ylab="")
abline(h=log(Initiation_rates,10),col="black",lwd=1)
text(x=35,y=log(Initiation_rates,10)[1]*0.95,"low",font=2,cex=1.2)
text(x=35,y=log(Initiation_rates,10)[2]*0.9,"mid",font=2,cex=1.2)
text(x=10,y=log(Initiation_rates,10)[3]*0.8,"high",font=2,cex=1.2)
for(i in c(seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,0.5,0.1))){
  if(grepl(as.numeric(i),pattern = "5$|^5e")){
    print(i)
    lines(x=c(-2,-0.5),y=log(c(i,i),10),xpd=T)
  }else{
    lines(x=c(-1.5,-0.5),y=log(c(i,i),10),xpd=T)
  }
}
axis(1,at = seq(0,40,10),labels = seq(0,40,10))
axis(2,at = c(-4,-3,-2,-1),labels = c(-4,-3,-2,-1),las=2)
mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=1.5,cex=1,xpd=T)
box()


# filled.contour2(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]])),
#                 col = c(grey.colors(n = 10,start = 1,end = 0,alpha = 1),rep(rgb(0,0,0,1),10)),nlevels = 20,plot.axes = c(axis(1,at = seq(0,40,5),labels = seq(0,40,5)),axis(2,at = seq(0,-5,-1),labels = seq(0,-5,-1))),
#                 key.axes = c(axis(4,at = seq(0,1,0.2),labels = seq(0,100,20))),key.title = "Ribosomal Density")
#contour.default(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]])))
#mtext(1,text = "Ribosomal Class",line=0,cex=1.2,xpd=F)
#mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=0,cex=1.2,xpd=T)
#mtext(4,text = "Ribosomal Density",line=1,at = -2,cex=1.2)

#colorbar
par(mar=c(1.1,1.1,1.1,1.1))

plot(0,0,col=rgb(0,0,0,0),axes=F,xlim=c(0,2),ylim=c(0,55),xlab="",ylab="")
for(i in 0:9){
  cols <- c(grey.colors(n = 12,start = 1,end = 0,alpha = 1)[4:12],rep(rgb(0,0,0,1),1))
  polygon(x = c(0.5,1.5,1.5,0.5,0.5), y = c(i*5,i*5,i*5+5,i*5+5,i*5),col=cols[i])
}
# for(j in 0:5){
#   lines(x=c(1.5,1.6),y=rep(j*10,2))
# }
mtext(3,text = "Ribosomal Location Probability",font=2)
text(x=2,y=seq(0,50,10),label= c(seq(0,45,10),">50"),las=2,xpd=T)
#dev.off()

par(mar=c(3.1,3.1,1.1,0.1))
#plot(0,0,col=rgb(0,0,0,0),axes=F,xlim=c(0,41),ylim=c(min(log(par.ranges[[1]],10)),max(log(par.ranges[[1]],10))))
image(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]])),
      col = c(rgb(1,1,1,1),grey.colors(n = 12,start = 1,end = 0,alpha = 1)[4:12],rep(rgb(0,0,0,1),10)),axes=F,xlab="",ylab="")
abline(h=log(Initiation_rates,10),col="black",lwd=1)
text(x=35,y=log(Initiation_rates,10)[1]*0.95,"low",font=2,cex=1.2)
text(x=35,y=log(Initiation_rates,10)[2]*0.9,"mid",font=2,cex=1.2)
text(x=10,y=log(Initiation_rates,10)[3]*0.8,"high",font=2,cex=1.2)
for(i in c(seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,0.5,0.1))){
  if(grepl(as.numeric(i),pattern = "5$|^5e")){
    print(i)
    lines(x=c(-2,-0.5),y=log(c(i,i),10),xpd=T)
  }else{
    lines(x=c(-1.5,-0.5),y=log(c(i,i),10),xpd=T)
  }
}

axis(1,at = seq(0,40,10),labels = seq(0,40,10))
axis(2,at = c(-4,-3,-2,-1),labels = c(-4,-3,-2,-1),las=2)
mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=1.5,cex=1,xpd=T)
box()

# filled.contour2(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],index]])),
#                col = c(grey.colors(n = 10,start = 1,end = 0,alpha = 1),rep(rgb(0,0,0,1),10)),nlevels = 20,plot.axes = c(axis(1,at = seq(0,40,5),labels = seq(0,40,5)),axis(2,at = seq(0,-5,-1),labels = seq(0,-5,-1))),
#                key.axes = c(axis(4,at = seq(0,1,0.2),labels = seq(0,100,20))),key.title = "Ribosomal Density")
#contour.default(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]])))
#mtext(1,text = "Ribosomal Class",line=0,cex=1.2,xpd=F)
#mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=0,cex=1.2,xpd=T)
#dev.off()


for(i in 3:1){
  if(i==1){
    par(mar=c(3.1,3.1,1.1,1.1))
  }else if(i==3){
    par(mar=c(1.1,3.1,1.1,1.1))
  }else{
    par(mar=c(1.1,3.1,1.1,1.1))
  }

  plot(0:41,U_pct_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],][kappa==Initiation_rates[i],index]],pch=20,xlab="",ylab="",ylim=c(0,0.7),axes=F)
  lines(0:41,U_pct_imax_median_plot1_int[,plot1.params[mu==par.ranges[[2]][10],][kappa==Initiation_rates[i],index]],pch=20)
  axis(1,at = seq(0,40,10),labels = seq(0,40,10),las=2)
  axis(2,at = seq(0,0.7,0.1),labels = seq(0,0.7,0.1),las=2)
  box()
  if(i == 1) {
    mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
    text(x=10,y=0.6,"low",font=2,cex=1.2)
  } else if(i == 2){
    mtext(2,text = "Ribosomal Location Probability",line=2,cex=1,xpd=F)
    text(x=10,y=0.6,"mid",font=2,cex=1.2)
  } else if(i==3){

    text(x=10,y=0.6,"high",font=2,cex=1.2)

  }

}


dev.off()



# Figure 2 Marking rate ---------------------------------------------------


png(paste(Sys.Date(),"_Figure2_Marking_Rate_range_medianlength.png",sep=""),width = 8,height = 6,units = "in",res=300)

layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,8,9),nrow = 3),widths = c(1,2,2,2,2,1,3),heights = c(1,1,1))
par(mar=c(3.1,0.1,0.1,0.1))

plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")

#low marking
par(mar=c(3.1,3.1,0.1,0.1))
image(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],index]])),
      col = c(rgb(1,1,1,1),grey.colors(n = 12,start = 1,end = 0,alpha = 1)[4:12],rep(rgb(0,0,0,1),10)),axes=F,xlab="",ylab="")
abline(h=log(Initiation_rates,10),col="black",lwd=1)
text(x=35,y=log(Initiation_rates,10)[1]*0.95,"low",font=2,cex=1.2)
text(x=35,y=log(Initiation_rates,10)[2]*0.9,"mid",font=2,cex=1.2)
text(x=10,y=log(Initiation_rates,10)[3]*0.8,"high",font=2,cex=1.2)
for(i in c(seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,0.5,0.1))){
  if(grepl(as.numeric(i),pattern = "5$|^5e")){
    print(i)
    lines(x=c(-2,-0.5),y=log(c(i,i),10),xpd=T)
  }else{
    lines(x=c(-1.5,-0.5),y=log(c(i,i),10),xpd=T)
  }
}
axis(1,at = seq(0,40,10),labels = seq(0,40,10))
axis(2,at = c(-4,-3,-2,-1),labels = c(-4,-3,-2,-1),las=2)
mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=1.5,cex=1,xpd=T)
box()

#mid marking
par(mar=c(3.1,0.1,1.1,0.1))
image(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][6],index]])),
      col = c(rgb(1,1,1,1),grey.colors(n = 12,start = 1,end = 0,alpha = 1)[4:12],rep(rgb(0,0,0,1),10)),axes=F,xlab="",ylab="")
abline(h=log(Initiation_rates,10),col="black",lwd=1)
text(x=35,y=log(Initiation_rates,10)[1]*0.95,"low",font=2,cex=1.2)
text(x=35,y=log(Initiation_rates,10)[2]*0.9,"mid",font=2,cex=1.2)
text(x=10,y=log(Initiation_rates,10)[3]*0.8,"high",font=2,cex=1.2)
# for(i in c(seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,0.5,0.1))){
#   if(grepl(as.numeric(i),pattern = "5$|^5e")){
#     print(i)
#     lines(x=c(-2,-0.5),y=log(c(i,i),10),xpd=T)
#   }else{
#     lines(x=c(-1.5,-0.5),y=log(c(i,i),10),xpd=T)
#   }
# }
axis(1,at = seq(0,40,10),labels = seq(0,40,10))
# axis(2,at = c(-4,-3,-2,-1),labels = c(-4,-3,-2,-1),las=2)
mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
# mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=1.5,cex=1,xpd=T)
box()

#High marking
par(mar=c(3.1,0.1,1.1,0.1))
image(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][2],index]])),
      col = c(rgb(1,1,1,1),grey.colors(n = 39,start = 1,end = 0,alpha = 1),rep(rgb(0,0,0,1),10)),axes=F,xlab="",ylab="")
abline(h=log(Initiation_rates,10),col="black",lwd=1)
text(x=35,y=log(Initiation_rates,10)[1]*0.95,"low",font=2,cex=1.2)
text(x=35,y=log(Initiation_rates,10)[2]*0.9,"mid",font=2,cex=1.2)
text(x=10,y=log(Initiation_rates,10)[3]*0.8,"high",font=2,cex=1.2)
# for(i in c(seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,0.5,0.1))){
#   if(grepl(as.numeric(i),pattern = "5$|^5e")){
#     print(i)
#     lines(x=c(-2,-0.5),y=log(c(i,i),10),xpd=T)
#   }else{
#     lines(x=c(-1.5,-0.5),y=log(c(i,i),10),xpd=T)
#   }
# }
axis(1,at = seq(0,40,10),labels = seq(0,40,10))
# axis(2,at = c(-4,-3,-2,-1),labels = c(-4,-3,-2,-1),las=2)
mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
# mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=1.5,cex=1,xpd=T)
box()

#sim
par(mar=c(3.1,0.1,1.1,0.1))
image(x=0:41,y=log(par.ranges[[1]],10),z=as.matrix((U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][1],index]])),
      col = c(rgb(1,1,1,1),grey.colors(n = 19,start = 1,end = 0,alpha = 1)[4:12],rep(rgb(0,0,0,1),70)),axes=F,xlab="",ylab="")
abline(h=log(Initiation_rates,10),col="black",lwd=1)
text(x=35,y=log(Initiation_rates,10)[1]*0.95,"low",font=2,cex=1.2)
text(x=35,y=log(Initiation_rates,10)[2]*0.9,"mid",font=2,cex=1.2)
text(x=10,y=log(Initiation_rates,10)[3]*0.8,"high",font=2,cex=1.2)
# for(i in c(seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,0.5,0.1))){
#   if(grepl(as.numeric(i),pattern = "5$|^5e")){
#     print(i)
#     lines(x=c(-2,-0.5),y=log(c(i,i),10),xpd=T)
#   }else{
#     lines(x=c(-1.5,-0.5),y=log(c(i,i),10),xpd=T)
#   }
# }
axis(1,at = seq(0,40,10),labels = seq(0,40,10))
# axis(2,at = c(-4,-3,-2,-1),labels = c(-4,-3,-2,-1),las=2)
mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
# mtext(2,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=1.5,cex=1,xpd=T)
box()


#colorbar
par(mar=c(1.1,1.1,1.1,1.1))

plot(0,0,col=rgb(0,0,0,0),axes=F,xlim=c(0,2),ylim=c(0,55),xlab="",ylab="")
for(i in 0:9){
  cols <- c(grey.colors(n = 12,start = 1,end = 0,alpha = 1)[4:12],rep(rgb(0,0,0,1),1))
  polygon(x = c(0.5,1.5,1.5,0.5,0.5), y = c(i*5,i*5,i*5+5,i*5+5,i*5),col=cols[i])
}
# for(j in 0:5){
#   lines(x=c(1.5,1.6),y=rep(j*10,2))
# }
mtext(3,text = "Ribosomal Location Probability",font=2)
text(x=2,y=seq(0,50,10),label= c(seq(0,45,10),">50"),las=2,xpd=T)
#dev.off()



for(i in 3:1){
  if(i==1){
    par(mar=c(3.1,4.1,1.1,0.1))
  }else if(i==3){
    par(mar=c(1.1,4.1,1.1,0.1))
  }else{
    par(mar=c(1.1,4.1,1.1,0.1))
  }
  plot(0:41,U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][10],][kappa==Initiation_rates[i],index]],type="l",lty=1,xlab="",ylab="",ylim=c(0,0.7),axes=F)
  lines(0:41,U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][6],][kappa==Initiation_rates[i],index]],lty=2)
  lines(0:41,U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][2],][kappa==Initiation_rates[i],index]],lty=3)
  lines(0:41,U_pct_imax_median_plot1[,plot1.params[mu==par.ranges[[2]][1],][kappa==Initiation_rates[i],index]],lty=4)
  axis(1,at = seq(0,40,10),labels = seq(0,40,10),las=2)
  axis(2,at = seq(0,0.7,0.1),labels = seq(0,0.7,0.1),las=2)
  box()
  if(i == 1) {
    mtext(1,text = "Ribosomal Class",line=2,cex=1,xpd=F)
    text(x=10,y=0.6,"low",font=2,cex=1.2)
  } else if(i == 2){
    mtext(2,text = "Ribosomal Location Probability",line=2,cex=1,xpd=F)
    text(x=10,y=0.6,"mid",font=2,cex=1.2)
    legend("topright",lty = 1:5,legend = c("1,500 (~99%)","115    (~50%)","11.5   (~1%)","1.15   (sim)"),
           title = "Half-life (min.)",xjust = 1,lwd=2)

  } else {
    text(x=10,y=0.6,"high",font=2,cex=1.2)
  }

}


dev.off()





















