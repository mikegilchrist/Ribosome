source("R/Polysome_functions_RAUC.R")
library(Biostrings)
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
#histogram of Protein lengths in arabidopsis.
#BiocManager::install("GenomicFeatures")
#library(GenomicFeatures)
At <- Biostrings::readDNAStringSet("../../Project_1/Results/20220117/CDS_arabidopsis.txt")
hist(width(At)/3)
At_lengths <- sort(width(At))
quantile(At_lengths,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))/3
hist(At_lengths[which(At_lengths<3040)]/3,breaks = seq(1,1030,10),axes=F,xlab="Protein length (aa)")
axes(1,at=c())
boxplot(log(At_lengths/3,10))

which(width(At)<300)
which(width(At)>300&width(At)<402)
At[which(width(At)==9)]

allnames <- unlist(lapply(names(At[which(width(At)>300&width(At)<402)]),function(i) base::strsplit(x = i,split = "\\|")[[1]][2]))

names(At_lengths) <- allnames

allAGI <- unlist(lapply(names(At),function(i) base::strsplit(x = i,split = "\\|")[[1]][1]))
allreps <- unlist(lapply(names(At),function(i) base::strsplit(x = i,split = "\\|")[[1]][2]))
#write.table(sort(unique(allAGI)),"../Results/20220301/20220301_All_At_AGI.txt",quote=F,col.names = F,row.names = F,sep="\t")
rep.models <- fread("../Results/20220301/At_Representative_Gene_models.txt",sep="\t",header = T,quote="",stringsAsFactors = F)
setkeyv(rep.models,"Locus Identifier")

GTF <- rtracklayer::import.gff("D:/Research/Araport11_GFF3_genes_transposons.201606.gtf")
At.Tx.Db <- makeTxDbFromGRanges(GTF)
seqlevels(At.Tx.Db)
columns(At.Tx.Db)
keytypes(At.Tx.Db)
keys(At.Tx.Db,"TXNAME")
CDS <- as.data.table(select(x = At.Tx.Db,keys=sort(unlist(rep.models[,`Gene Model Name`])),columns = c("CDSSTART","CDSEND"),keytype = "TXNAME"))
CDS[,"CDS_LENGTH":= (CDSEND-CDSSTART)+1]
CDS_LENGTHS <- CDS[,sum(CDS_LENGTH),TXNAME]
setkeyv(CDS_LENGTHS,"TXNAME")
cbind(CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),],CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),2]%%3)
rep.models[as.character(gsub(x=unlist(CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),1]),replacement = "",pattern = "\\.\\d")),c(2,5)]
CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),V1:=V1-CDS_LENGTHS[which((CDS_LENGTHS[,2]%%3)!=0),2]%%3]


hist(CDS_LENGTHS[,V1]/3)
quantile(CDS_LENGTHS[,V1],c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)/27
quantile(CDS_LENGTHS[,V1],c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)/3
hist(At_lengths[which(At_lengths<3000)]/3,breaks = seq(1,1010,10))

sort(CDS_LENGTHS[,V1])[1:100]
CDS_LENGTHS[order(CDS_LENGTHS[,V1]),][1:100]
rep.models[as.character(gsub(x=unlist(CDS_LENGTHS[order(CDS_LENGTHS[,V1]),1][1:100]),replacement = "",pattern = "\\.\\d")),c(1,3)]
rep.models["AT1G24660"]
desc <- rep.models[as.character(gsub(x=unlist(CDS_LENGTHS[,1]),replacement = "",pattern = "\\.\\d")),c(1,3)]

hypothetical <- as.numeric(unlist(lapply(desc[,2],function(i)grep(x = i,pattern = "^hypothetical protein"))))
plot(unlist(hypothetical)/length(hypothetical))
hist(unlist(hypothetical),breaks = seq(0,28000,500))


hist(unlist(CDS_LENGTHS[,2]/3),col="red",breaks = seq(0,18000/3,100))
hist(unlist(CDS_LENGTHS[unlist(hypothetical),2]/3),breaks = seq(0,18000/3,100),add=T,col="white")
hist(unlist(CDS_LENGTHS[-unlist(hypothetical),2])/3,breaks = seq(0,18000,100)/3)
quantile(sort(unlist(CDS_LENGTHS[-unlist(hypothetical),2])),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)/3
quantile(sort(unlist(CDS_LENGTHS[,2])),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)/3


#transmembrane
desc <- desc[-unlist(hypothetical),]
transmembrane <- as.numeric(unlist(lapply(desc[,2],function(i)grep(x = i,pattern = "^transmembrane"))))

hist(unlist(CDS_LENGTHS[-hypothetical,2]/3),col="red",breaks = seq(0,18000/3,100))
hist(unlist(CDS_LENGTHS[-hypothetical,][unlist(transmembrane),2]/3),breaks = seq(0,18000/3,100),add=T,col="white")
hist(unlist(CDS_LENGTHS[-unlist(hypothetical),2])/3,breaks = seq(0,18000,100)/3)
quantile(sort(unlist(CDS_LENGTHS[-unlist(hypothetical),][-unlist(transmembrane),2]/3)),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
quantile(sort(unlist(CDS_LENGTHS[,2])),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)/3

#This gene encodes a small protein and has
desc <- desc[-unlist(transmembrane),]
evo <- as.numeric(unlist(lapply(desc[,2],function(i)grep(x = i,pattern = "^This gene encodes a small protein and has"))))

hist(unlist(CDS_LENGTHS[-hypothetical,][-transmembrane,2]/3),col="red",breaks = seq(0,18000/3,100))
hist(unlist(CDS_LENGTHS[-hypothetical,][-transmembrane,][evo,2]/3),breaks = seq(0,18000/3,100),add=T,col="white")
hist(unlist(CDS_LENGTHS[-unlist(hypothetical),2])/3,breaks = seq(0,18000,100)/3)
quantile(sort(unlist(CDS_LENGTHS[-unlist(hypothetical),][-transmembrane,][-evo,2]/3)),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
protein_quant <- quantile(sort(unlist(CDS_LENGTHS[,2])),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)/3
hist(log(unlist(CDS_LENGTHS[-hypothetical,][-transmembrane,][-evo,2]/3),10),breaks = seq(0,log(18000/3,2),1/10),col="white")


comp_and_exp_w_annotation <- CDS_LENGTHS[-hypothetical,][-transmembrane,][-evo,]

hist(unlist(comp_and_exp_w_annotation[,2])/3,breaks = seq(0,18000/3,100),col="white")
comp_trimmed <- CDS_LENGTHS[-hypothetical,][-transmembrane,][-evo,]
comp_trimmed[V1>4500,2] <- 4500


# Protein length distribution ---------------------------------------------

png("../Results/20220301/Figure1_At_Protein_length_histogram.png",width = 5,height = 5,units = "in",res = 300)
prot.length.hist <- hist(unlist(comp_trimmed[,2])/3,breaks = seq(0,1500,20),ylim=c(0,1100),col="white",axes=F,xlab="",main="Arabidopsis Protein Lengths")
abline(v=quantile(comp_trimmed[,2]/3,c(0.1,0.5,0.9),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,10^max(nchar(prot.length.hist$counts)),250),labels = seq(0,10^max(nchar(prot.length.hist$counts)),250))
axis(1,at=seq(0,1500,500),labels = c("0","500","1000",">1500"))
mtext(1,line=2.25,text = "Protein Length (a.a.)")
legend("topright",legend = c("Small","10% (140aa)","Medium","50% (368aa)","Large","90% (780aa)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
box()
dev.off()


# Half Life Distribution --------------------------------------------------
HL_At <- fread("../../Project_1/Results/20220117/Half life Sorenson.csv",stringsAsFactors = F,header=T)
hist(log(unlist(HL_At[which(unlist(HL_At[,3])<100000),3]),10))
median(unlist(HL_At[,3]),na.rm=T)
quantile(unlist(HL_At[which(unlist(HL_At[,3])<100000),3]),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
quantile(unlist(HL_At[,3]),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))

hist(unlist(HL_At[,3]))

HL_At <- fread("../Results/20220301/Sorenson_results_WT.csv",stringsAsFactors = F,header=T)

Sorenson_model <- function(time,alpha,beta){
  tmp.res <- numeric(length(time))
  if(beta!=0){
    for(i in 1:length(tmp.res)){
      tmp.res[i] <- exp(-1*(alpha/beta)*(1-exp(-1*beta*time[i])))
    }
  }else{
    for(i in 1:length(tmp.res)){
      tmp.res[i] <- exp(-1*(alpha*time[i]))
    }
  }

  return(tmp.res)
}
time = c(seq(0,10,0.1),seq(11,5000,1))
Sorenson_model(time,alpha = as.numeric(HL_At[1,2]),beta = as.numeric(HL_At[1,3]))[3000]

Sorenson_model_results <- apply(as.data.frame(HL_At),1,function(i) Sorenson_model(time = c(seq(0,10,0.1),seq(11,5000,1)),alpha = as.numeric(i[2]),beta = as.numeric(i[3])))
dim(Sorenson_model_results)

Sorenson_model_results[,1]

HL_CALC <- apply(Sorenson_model_results,2,function(i) which(i<0.5))
HL_index <- lapply(HL_CALC,function(i) which(length(i)==0))
length(which(HL_index==1))
halves <- Sorenson_model_results[,-which(HL_index==1)]
halves_index <- lapply(HL_CALC,function(i) i[1])

for(i in 1:10){#ncol(Sorenson_model_results)){
  if(is.na(halves_index[[i]])){
    next
  }
  print(i)
  print(Sorenson_model_results[(halves_index[[i]]-3):((halves_index[[i]]+3)),i])
  print(time[(i-1):i])
}

counter=1
times.tmp <- numeric(length = (ncol(Sorenson_model_results)-length(which(is.na(halves_index)))))
names.tmp <- character(length = (ncol(Sorenson_model_results)-length(which(is.na(halves_index)))))

for(i in 1:ncol(Sorenson_model_results)){
#  print(i)
#  print(counter)
  if(is.na(halves_index[[i]])){
#    print("next")
    next
  }
  times.tmp[counter] <- time[i]
  names.tmp[counter] <- HL_At[i,1]
  counter = counter + 1
}
names(times.tmp) <- names.tmp

hist(times.tmp)
median(times.tmp,na.rm = T)
mean(times.tmp,na.rm = T)
max(times.tmp,na.rm = T)

HL_minutes <- unlist(log(2)/HL_At[-which(HL_At[,2]==0),2])

HL_minutes[(which(HL_minutes>10000))]<-10000

png("../Results/20220301/Figure1_At_Protein_length_histogram.png",width = 5,height = 5,units = "in",res = 300)
half.life.hist <- hist(HL_minutes,col="white",breaks = seq(0,10000,100),ylim=c(0,10),axes=F,xlab="",main="Arabidopsis Protein Lengths")
abline(v=quantile(comp_trimmed[,2]/3,c(0.1,0.5,0.9),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,10^max(nchar(prot.length.hist$counts)),250),labels = seq(0,10^max(nchar(prot.length.hist$counts)),250))
axis(1,at=seq(0,1500,500),labels = c("0","500","1000",">1500"))
mtext(1,line=2.25,text = "Protein Length (a.a.)")
legend("topright",legend = c("Small","10% (140aa)","Medium","50% (368aa)","Large","90% (780aa)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
box()
dev.off()




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




