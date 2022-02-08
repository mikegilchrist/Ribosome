source("https://bioconductor.org/biocLite.R")
install.packages("RColorBrewer")
install.packages("cluster")
biocLite("AnnotationDbi")
biocLite("limma")
biocLite("statmod")
biocLite("affy")
#unloadNamespace("affy")
biocLite("GO.db")
BiocManager::install("gcrma")
require(AnnotationDbi)
require(GO.db)
require(affy)
require(ath1121501.db)
require(limma)
require(statmod)
require(gplots)
library(RColorBrewer)
library(gcrma)
require(cluster)

setwd("~/utk/Research/Project_1/Results/2021_06_reanalysis/")


Dir <- "../../data/cel/"

Files <- lapply(list.dirs(Dir)[2:4],function(i) as.vector(unlist(list.celfiles(i,full.names=T))))

names(Files) <- c("Rep 1","Rep 2","Rep 3")

unlist(Files)


Data <- ReadAffy(filenames = unname(unlist(Files)))

Data$sample<- lapply(unlist(Files), function(i) strsplit(x = strsplit(i,split = "Col_")[[1]][2],split = "\\.")[[1]][1])

Data_TX_col <- Data[,grep("T",colnames(Data))]
Data_TL_col <- Data[,-grep("T",colnames(Data))]

id_data <<- AnnotationDbi::select(ath1121501.db, keys = keys(ath1121501.db),
                                  columns = c("TAIR","SYMBOL","PROBEID","CHRLOC","CHRLOCEND"),
                                  keytype = "PROBEID")
id_data[, "LENGTH"] = abs(id_data$CHRLOCEND - id_data$CHRLOC) + 1
id_data$LENGTH[(is.null(id_data$LENGTH))] = NA

affy.2.tair = function(data, replace.x = F) {
  if(length(data) > 0) {
    my_list = lapply(data, function(x) id_data$TAIR[(id_data$PROBEID %in% x)])
    names(my_list) = data
    if(replace.x == T) {
      replace_vector = sapply(1:length(my_list), function(x) {if (length(my_list[[x]]) >= 1 && !is.na(my_list[[x]])) {
        return(my_list[[x]][(!is.na(my_list[[x]]))][1])} else {return(data[x])}}) 
      return(replace_vector)
    }
    else {
      return(my_list)
      
    }
  }
  else {
    return(NULL)
  }
}


#Plot Raw Data boxplots
png(filename = "Raw_Data_boxplot.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(log(pm(Data),2),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample,cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


png(filename = "Raw_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(log(pm(Data)[,c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],2),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()

#Plot Raw data correlations
old.dir <- getwd()
setwd(paste(old.dir,"Results/20180524_rerun/TX/",sep="/"))


cplot <- function(Sample_names,Data,plot='Log',label=""){
  for(i in 1:(length(Sample_names)-1)){
    for(j in (i+1):length(Sample_names)){
      title1 <- paste(c(Sample_names[i],"vs",Sample_names[j],label),collapse=" ")
      title2 <- paste(c(Sample_names[i],"vs",Sample_names[j],"log",label),collapse=" ")
      if(plot=="Raw"|plot=="Both"){
        png(filename =paste(title1,".png",sep="") ,width = 2500,height = 2000,units = "px",res=300)  
        plot(pm(Data[,i]),pm(Data[,j]),main=title1,xlab=Sample_names[i],ylab=Sample_names[j],xlim=c(0,42000),ylim=c(0,42000))
        abline(a = c(0,1),col="red")
        box()
        dev.off()
      }
      if(plot=="Log"|plot=="Both"){
        png(filename =paste(title2,".png",sep="") ,width = 2500,height = 2000,units = "px",res=300)
        plot(log(pm(Data[,i]),2),log(pm(Data[,j]),2),main=title2,xlab=Sample_names[i],ylab=Sample_names[j],xlim=c(0,20),ylim=c(0,20))
        abline(a = c(0,1),col="red")
        box()
        dev.off()
      }
    }
  }
}  

cplot2 <- function(Sample_names,Data,label=""){
  for(i in 1:(length(Sample_names)-1)){
    for(j in (i+1):length(Sample_names)){
      title1 <- paste(c(Sample_names[i],"vs",Sample_names[j],label),collapse=" ")
      png(filename =paste(title1,".png",sep="") ,width = 2500,height = 2000,units = "px",res=300)  
      plot((Data[,i]),(Data[,j]),main=title1,xlab=Sample_names[i],ylab=Sample_names[j],xlim=c(0,18),ylim=c(0,18))
      abline(a = c(0,1),col="red")
      box()
      dev.off()
    }
  }
}  

T_samples <- grep(x = Data$sample,pattern = "^T")
NP <- grep(x = Data$sample,pattern = "^NP")
F5 <- grep(x = Data$sample,pattern = "^5")
F6 <- grep(x = Data$sample,pattern = "^6")
F7 <- grep(x = Data$sample,pattern = "^7")
F8 <- grep(x = Data$sample,pattern = "^8")
F9 <- grep(x = Data$sample,pattern = "^9")
F10 <- grep(x = Data$sample,pattern = "^10")
F11 <- grep(x = Data$sample,pattern = "^11")
F12 <- grep(x = Data$sample,pattern = "^12")

cplot(Sample_names = Data$sample[T_samples],Data = Data[,T_samples],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[NP],Data = Data[,NP],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[F5],Data = Data[,F5],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[F6],Data = Data[,F6],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[F7],Data = Data[,F7],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[F8],Data = Data[,F8],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[F9],Data = Data[,F9],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[F10],Data = Data[,F10],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[F11],Data = Data[,F11],plot = "Both",label = "Raw")
cplot(Sample_names = Data$sample[F12],Data = Data[,F12],plot = "Both",label = "Raw")

setwd(old.dir)
count =0
for(i in seq(47,1,-1)){
  count = count + i
}


raw_cor <- cor(pm(Data),method = "spearman")
raw_cor_log <- cor(log(pm(Data),2),method = "spearman")
colnames(raw_cor) <- Data$sample
rownames(raw_cor) <- Data$sample
colnames(raw_cor_log) <- Data$sample
rownames(raw_cor_log) <- Data$sample
write.csv(raw_cor,"Raw cor.csv")
write.csv(raw_cor_log,"Raw cor log.csv")


#Mas5 TX
m5calls <- mas5calls(Data)
m5pvals <-  assayData(m5calls)[["se.exprs"]]
m5calls <- as.matrix(exprs(m5calls))

m5stats <- apply(m5calls,2, function(i) table(i))

View(head(m5calls))

plot(seq(0.25,.75,0.25),m5stats[3,c(1,11,21)],ylim=c(0,18000),xlim=c(0,11))
lines(seq(0.25,.75,0.25),m5stats[3,c(1,11,21)],ylim=c(0,18000),xlim=c(0,11))
lines(x=c(0.25,0.75),rep(mean(m5stats[3,c(1,11,21)]),2),col="red",lwd=3)
for(i in 2:10){
  points(seq(0.25,.75,0.25)+i,m5stats[3,c(i,10+i,20+i)])
  lines(seq(0.25,.75,0.25)+i,m5stats[3,c(i,10+i,20+i)])
  lines(x=c(0.25,0.75)+i,rep(mean(m5stats[3,c(i,10+i,20+i)]),2),col="red",lwd=3)
}
for(i in 1:10){
  points(seq(0.25,.75,0.25)+i,m5stats[2,c(i,10+i,20+i)],pch=2)
  lines(seq(0.25,.75,0.25)+i,m5stats[2,c(i,10+i,20+i)])
  lines(x=c(0.25,0.75)+i,rep(mean(m5stats[2,c(i,10+i,20+i)]),2),col="red",lwd=3)
}
for(i in 1:10){
  points(seq(0.25,.75,0.25)+i,m5stats[1,c(i,10+i,20+i)],pch=0)
  lines(seq(0.25,.75,0.25)+i,m5stats[1,c(i,10+i,20+i)])
  lines(x=c(0.25,0.75)+i,rep(mean(m5stats[1,c(i,10+i,20+i)]),2),col="red",lwd=3)
}

colSums(m5stats)
Mmatall <- matrix(ncol=10,nrow=nrow(m5calls))
for(i in 1:10){
  Mmatall[,i] <- apply(m5calls[,c(i,10+i,20+i)],1,function(j) all(j=="M"))
  
}
Pmatall <- matrix(ncol=10,nrow=nrow(m5calls))
for(i in 1:10){
  Pmatall[,i] <- apply(m5calls[,c(i,10+i,20+i)],1,function(j) all(j=="P"))
  
}
Amatall <- matrix(ncol=10,nrow=nrow(m5calls))
for(i in 1:10){
  Amatall[,i] <- apply(m5calls[,c(i,10+i,20+i)],1,function(j) all(j=="A"))
  
}
Mixedmatany <- matrix(ncol=10,nrow=nrow(m5calls))
for(i in 1:10){
  Mixedmatany[,i] <- apply(m5calls[,c(i,10+i,20+i)],1,function(j) (any(j=="M")&(length(which(j=="M"))<3))|(any(j=="A")&(length(which(j=="A"))<3)))
  
}

a<-apply(Pmatall,2,function(i) length(which(i)))
b<-apply(Amatall,2,function(i) length(which(i)))
c<-apply(Mmatall,2,function(i) length(which(i)))
d<-apply(Mixedmatany,2,function(i) length(which(i)))
a+b+c+d

plot(1:10,a,ylim=c(0,25000),pch=1)
lines(1:10,a)
points(1:10,b,pch=0)
lines(1:10,b)
points(1:10,c,pch=2)
lines(1:10,c)
points(1:10,d,pch=3)
lines(1:10,d)
points(1:10,a+b+c+d,pch=4)
lines(1:10,a+b+c+d)
legend("top",horiz = T,legend = c("Pcalls","Acalls","Mcalls","Mixed","All"),pch = c(1,0,2,3,4))

Pall_index <-which(Pmatall[,1])
length(Pall_index)
#This is what was done in Ben's data
Pall_TX_TL_index <- which(apply(Pmatall[Pall_index,2:10],1,function(i) length(which(i))>=1))
Filtered_data_rn <- rownames(m5calls[Pall_index[Pall_TX_TL_index],])
Filtered_data_rn <- Filtered_data_rn[-grep(x = Filtered_data_rn,pattern = "^AFFX")]
length(grep(x = Filtered_data_rn,pattern = "_s_"))
length(grep(x = Filtered_data_rn,pattern = "^AFFX"))


Filtered_data_TAIR <- affy.2.tair(Filtered_data_rn)
Filtered_data_TAIR <- lapply(Filtered_data_TAIR,unique)
GDA_TAIR <- affy.2.tair(GDA[,1])
GDA_TAIR <- lapply(GDA_TAIR,unique)
length(which(unlist(lapply(Filtered_data_TAIR,function(i) length(i)>1))))
Filtered_data_rn[which(unlist(lapply(Filtered_data_TAIR,function(i) length(i)>1)))]
Filtered_data_TAIR[which(unlist(lapply(Filtered_data_TAIR,function(i) length(i)>1)))]
length(Filtered_data_TAIR[grep(x = Filtered_data_TAIR,pattern = "ATM")])
Filtered_data_rn <- Filtered_data_rn[-grep(x = Filtered_data_TAIR,pattern = "ATM")]
Filtered_data_TAIR <- Filtered_data_TAIR[-grep(x = Filtered_data_TAIR,pattern = "ATM")]
length(Filtered_data_TAIR[grep(x = Filtered_data_TAIR,pattern = "ATC")])
Filtered_data_rn <- Filtered_data_rn[-grep(x = Filtered_data_TAIR,pattern = "ATC")]
Filtered_data_TAIR <- Filtered_data_TAIR[-grep(x = Filtered_data_TAIR,pattern = "ATC")]
length(Filtered_data_TAIR[which(is.na(Filtered_data_TAIR))])
length(Filtered_data_TAIR)
length(Filtered_data_rn)


which((GDA[,1]%in%Filtered_data_rn)==F)
GDA_TAIR[which((GDA_TAIR%in%Filtered_data_TAIR)==F)]
Filtered_data_TAIR[which((Filtered_data_TAIR%in%GDA_TAIR)==F)]


Filtered_data_rn[which((Filtered_data_rn%in%GDA[,1])==F)]
GDA[which((GDA[,1]%in%Filtered_data_rn)==F),1]
Filtered_data_TAIR[which((Filtered_data_rn%in%GDA[,1])==F)]


length(which(apply(Pmatall[Pall_index,2:10],1,function(i) length(which(i))>=2)))

length(which(apply(Pmatall[,2:10],1,function(i) length(which(i))>=1)))
length(which(apply(Pmatall[,2:10],1,function(i) length(which(i))>=2)))
length(which(apply(Pmatall[,2:10],1,function(i) length(which(i))>=3)))

dim(m5calls[-which(rownames(m5calls)%in%names(Filtered_data_TAIR)==F),])





old.dir <- getwd()
setwd(paste(old.dir,"TX/",sep="/"))
colnames(col_m5calls_TX) <- sub(".CEL$","",sub("_AVA\\d{3}_","_",colnames(col_m5calls_TX)))
write.csv(col_m5calls_TX,"col_M5callsTX.csv")
setwd(old.dir)

TX_T <- apply(col_m5calls_TX,1,function(i) all(i=="P")==T)

index_TX <- which(TX_T==T)

#rma All together
All_data_rma <- rma(object = Data[,])
All_data_rma_expr <- exprs(All_data_rma) 

png(filename = "RMA_Data_boxplot_rep_order_all_together.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(All_data_rma_expr,las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


#rma seperate
T_data_rma <- rma(object = Data[,grep(Data$sample,pattern = "T")])
T_data_rma_expr <- exprs(T_data_rma) 

TL_data_rma <- rma(object = Data[,-grep(Data$sample,pattern = "T")])
TL_data_rma_expr <- exprs(TL_data_rma) 

png(filename = "RMA_Data_boxplot_rep_order_separate.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(T_data_rma_expr,TL_data_rma_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


#gcrma together
All_data_gcrma <- gcrma(object = Data)
All_data_gcrma_expr <- exprs(All_data_gcrma) 

#gcrma seperate
T_data_gcrma <- gcrma(object = Data[,grep(Data$sample,pattern = "T")])
T_data_gcrma_expr <- exprs(T_data_gcrma) 

TL_data_gcrma <- gcrma(object = Data[,-grep(Data$sample,pattern = "T")])
TL_data_gcrma_expr <- exprs(TL_data_gcrma) 

#rma All separate
T_data_rma <- rma(object = Data[,T_samples])
T_data_rma_expr <- exprs(T_data_rma) 

NP_data_rma <- rma(object = Data[,NP])
NP_data_rma_expr <- exprs(NP_data_rma) 

F5_data_rma <- rma(object = Data[,F5])
F5_data_rma_expr <- exprs(F5_data_rma) 

F6_data_rma <- rma(object = Data[,F6])
F6_data_rma_expr <- exprs(F6_data_rma) 

F7_data_rma <- rma(object = Data[,F7])
F7_data_rma_expr <- exprs(F7_data_rma) 

F8_data_rma <- rma(object = Data[,F8])
F8_data_rma_expr <- exprs(F8_data_rma) 

F9_data_rma <- rma(object = Data[,F9])
F9_data_rma_expr <- exprs(F9_data_rma) 

F10_data_rma <- rma(object = Data[,F10])
F10_data_rma_expr <- exprs(F10_data_rma) 

F11_data_rma <- rma(object = Data[,F11])
F11_data_rma_expr <- exprs(F11_data_rma) 

F12_data_rma <- rma(object = Data[,F12])
F12_data_rma_expr <- exprs(F12_data_rma) 

png(filename = "RMA_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(NP_data_rma_expr,F5_data_rma_expr,F6_data_rma_expr,F7_data_rma_expr,F8_data_rma_expr,
              F9_data_rma_expr,F10_data_rma_expr,F11_data_rma_expr,F12_data_rma_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()



#gcrma all separate
T_data_gcgcrma <- gcrma(object = Data[,T_samples])
T_data_gcrma_expr <- exprs(T_data_gcrma) 

NP_data_gcrma <- gcrma(object = Data[,NP])
NP_data_gcrma_expr <- exprs(NP_data_gcrma) 

F5_data_gcrma <- gcrma(object = Data[,F5])
F5_data_gcrma_expr <- exprs(F5_data_gcrma) 

F6_data_gcrma <- gcrma(object = Data[,F6])
F6_data_gcrma_expr <- exprs(F6_data_gcrma) 

F7_data_gcrma <- gcrma(object = Data[,F7])
F7_data_gcrma_expr <- exprs(F7_data_gcrma) 

F8_data_gcrma <- gcrma(object = Data[,F8])
F8_data_gcrma_expr <- exprs(F8_data_gcrma) 

F9_data_gcrma <- gcrma(object = Data[,F9])
F9_data_gcrma_expr <- exprs(F9_data_gcrma) 

F10_data_gcrma <- gcrma(object = Data[,F10])
F10_data_gcrma_expr <- exprs(F10_data_gcrma) 

F11_data_gcrma <- gcrma(object = Data[,F11])
F11_data_gcrma_expr <- exprs(F11_data_gcrma) 

F12_data_gcrma <- gcrma(object = Data[,F12])
F12_data_gcrma_expr <- exprs(F12_data_gcrma) 

png(filename = "Gcrma_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(NP_data_gcrma_expr,F5_data_gcrma_expr,F6_data_gcrma_expr,F7_data_gcrma_expr,F8_data_gcrma_expr,
              F9_data_gcrma_expr,F10_data_gcrma_expr,F11_data_gcrma_expr,F12_data_gcrma_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


# RMA substraction only ---------------------------------------------------
T_data_rma_sub <- bg.correct(object = Data[,T_samples],method = "rma")
T_data_rma_sub_expr <- exprs(T_data_rma_sub) 

NP_data_rma_sub <- bg.correct(object = Data[,NP],method = "rma")
NP_data_rma_sub_expr <- exprs(NP_data_rma_sub) 

F5_data_rma_sub <- bg.correct(object = Data[,F5],method = "rma")
F5_data_rma_sub_expr <- exprs(F5_data_rma_sub) 

F6_data_rma_sub <- bg.correct(object = Data[,F6],method = "rma")
F6_data_rma_sub_expr <- exprs(F6_data_rma_sub) 

F7_data_rma_sub <- bg.correct(object = Data[,F7],method = "rma")
F7_data_rma_sub_expr <- exprs(F7_data_rma_sub) 

F8_data_rma_sub <- bg.correct(object = Data[,F8],method = "rma")
F8_data_rma_sub_expr <- exprs(F8_data_rma_sub) 

F9_data_rma_sub <- bg.correct(object = Data[,F9],method = "rma")
F9_data_rma_sub_expr <- exprs(F9_data_rma_sub) 

F10_data_rma_sub <- bg.correct(object = Data[,F10],method = "rma")
F10_data_rma_sub_expr <- exprs(F10_data_rma_sub) 

F11_data_rma_sub <- bg.correct(object = Data[,F11],method = "rma")
F11_data_rma_sub_expr <- exprs(F11_data_rma_sub) 

F12_data_rma_sub <- bg.correct(object = Data[,F12],method = "rma")
F12_data_rma_sub_expr <- exprs(F12_data_rma_sub) 

png(filename = "rma_sub_only_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(log2(cbind(T_data_rma_sub_expr,NP_data_rma_sub_expr,F5_data_rma_sub_expr,F6_data_rma_sub_expr,F7_data_rma_sub_expr,F8_data_rma_sub_expr,
              F9_data_rma_sub_expr,F10_data_rma_sub_expr,F11_data_rma_sub_expr,F12_data_rma_sub_expr)),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


# RMA summarization only --------------------------------------------------
T_data_rma_sum <- rma(object = Data[,T_samples],normalize=F,background = F)
T_data_rma_sum_expr <- exprs(T_data_rma_sum) 

NP_data_rma_sum <- rma(object = Data[,NP],normalize=F,background = F)
NP_data_rma_sum_expr <- exprs(NP_data_rma_sum) 

F5_data_rma_sum <- rma(object = Data[,F5],normalize=F,background = F)
F5_data_rma_sum_expr <- exprs(F5_data_rma_sum) 

F6_data_rma_sum <- rma(object = Data[,F6],normalize=F,background = F)
F6_data_rma_sum_expr <- exprs(F6_data_rma_sum) 

F7_data_rma_sum <- rma(object = Data[,F7],normalize=F,background = F)
F7_data_rma_sum_expr <- exprs(F7_data_rma_sum) 

F8_data_rma_sum <- rma(object = Data[,F8],normalize=F,background = F)
F8_data_rma_sum_expr <- exprs(F8_data_rma_sum) 

F9_data_rma_sum <- rma(object = Data[,F9],normalize=F,background = F)
F9_data_rma_sum_expr <- exprs(F9_data_rma_sum) 

F10_data_rma_sum <- rma(object = Data[,F10],normalize=F,background = F)
F10_data_rma_sum_expr <- exprs(F10_data_rma_sum) 

F11_data_rma_sum <- rma(object = Data[,F11],normalize=F,background = F)
F11_data_rma_sum_expr <- exprs(F11_data_rma_sum) 

F12_data_rma_sum <- rma(object = Data[,F12],normalize=F,background = F)
F12_data_rma_sum_expr <- exprs(F12_data_rma_sum) 

png(filename = "rma_sum_only_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(T_data_rma_sum_expr,NP_data_rma_sum_expr,F5_data_rma_sum_expr,F6_data_rma_sum_expr,F7_data_rma_sum_expr,F8_data_rma_sum_expr,
                   F9_data_rma_sum_expr,F10_data_rma_sum_expr,F11_data_rma_sum_expr,F12_data_rma_sum_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


# RMA sub and summarization only ------------------------------------------

T_data_rma_sub_sum <- rma(object = Data[,T_samples],normalize=F,background = T)
T_data_rma_sub_sum_expr <- exprs(T_data_rma_sub_sum) 

NP_data_rma_sub_sum <- rma(object = Data[,NP],normalize=F,background = T)
NP_data_rma_sub_sum_expr <- exprs(NP_data_rma_sub_sum) 

F5_data_rma_sub_sum <- rma(object = Data[,F5],normalize=F,background = T)
F5_data_rma_sub_sum_expr <- exprs(F5_data_rma_sub_sum) 

F6_data_rma_sub_sum <- rma(object = Data[,F6],normalize=F,background = T)
F6_data_rma_sub_sum_expr <- exprs(F6_data_rma_sub_sum) 

F7_data_rma_sub_sum <- rma(object = Data[,F7],normalize=F,background = T)
F7_data_rma_sub_sum_expr <- exprs(F7_data_rma_sub_sum) 

F8_data_rma_sub_sum <- rma(object = Data[,F8],normalize=F,background = T)
F8_data_rma_sub_sum_expr <- exprs(F8_data_rma_sub_sum) 

F9_data_rma_sub_sum <- rma(object = Data[,F9],normalize=F,background = T)
F9_data_rma_sub_sum_expr <- exprs(F9_data_rma_sub_sum) 

F10_data_rma_sub_sum <- rma(object = Data[,F10],normalize=F,background = T)
F10_data_rma_sub_sum_expr <- exprs(F10_data_rma_sub_sum) 

F11_data_rma_sub_sum <- rma(object = Data[,F11],normalize=F,background = T)
F11_data_rma_sub_sum_expr <- exprs(F11_data_rma_sub_sum) 

F12_data_rma_sub_sum <- rma(object = Data[,F12],normalize=F,background = T)
F12_data_rma_sub_sum_expr <- exprs(F12_data_rma_sub_sum) 

png(filename = "rma_sub_sum_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(T_data_rma_sub_sum_expr,NP_data_rma_sub_sum_expr,F5_data_rma_sub_sum_expr,F6_data_rma_sub_sum_expr,F7_data_rma_sub_sum_expr,F8_data_rma_sub_sum_expr,
              F9_data_rma_sub_sum_expr,F10_data_rma_sub_sum_expr,F11_data_rma_sub_sum_expr,F12_data_rma_sub_sum_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


# RMA sub and summarization with constant norm ------------------------------------------

T_data_rma_sub_sum_norm <- expresso(Data[,T_samples],normalize=T,normalize.method="constant",
                           bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                           pmcorrect.method = "pmonly")
T_data_rma_sub_sum_norm_expr <- exprs(T_data_rma_sub_sum_norm) 

NP_data_rma_sub_sum_norm <- expresso(Data[,NP],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                 pmcorrect.method = "pmonly")
NP_data_rma_sub_sum_norm_expr <- exprs(NP_data_rma_sub_sum_norm) 

F5_data_rma_sub_sum_norm <- expresso(Data[,F5],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                 pmcorrect.method = "pmonly")
F5_data_rma_sub_sum_norm_expr <- exprs(F5_data_rma_sub_sum_norm) 

F6_data_rma_sub_sum_norm <- expresso(Data[,F6],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                 pmcorrect.method = "pmonly")
F6_data_rma_sub_sum_norm_expr <- exprs(F6_data_rma_sub_sum_norm) 

F7_data_rma_sub_sum_norm <- expresso(Data[,F7],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                 pmcorrect.method = "pmonly")
F7_data_rma_sub_sum_norm_expr <- exprs(F7_data_rma_sub_sum_norm) 

F8_data_rma_sub_sum_norm <- expresso(Data[,F8],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                 pmcorrect.method = "pmonly")
F8_data_rma_sub_sum_norm_expr <- exprs(F8_data_rma_sub_sum_norm) 

F9_data_rma_sub_sum_norm <- expresso(Data[,F9],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                 pmcorrect.method = "pmonly")
F9_data_rma_sub_sum_norm_expr <- exprs(F9_data_rma_sub_sum_norm) 

F10_data_rma_sub_sum_norm <- expresso(Data[,F10],normalize=T,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                  pmcorrect.method = "pmonly")
F10_data_rma_sub_sum_norm_expr <- exprs(F10_data_rma_sub_sum_norm) 

F11_data_rma_sub_sum_norm <- expresso(Data[,F11],normalize=T,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                  pmcorrect.method = "pmonly")
F11_data_rma_sub_sum_norm_expr <- exprs(F11_data_rma_sub_sum_norm) 

F12_data_rma_sub_sum_norm <- expresso(Data[,F12],normalize=T,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "rma",summary.method = "medianpolish",
                                  pmcorrect.method = "pmonly")
F12_data_rma_sub_sum_norm_expr <- exprs(F12_data_rma_sub_sum_norm) 


png(filename = "rma_sub_sum_norm_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(T_data_rma_sub_sum_norm_expr,NP_data_rma_sub_sum_norm_expr,F5_data_rma_sub_sum_norm_expr,F6_data_rma_sub_sum_norm_expr,F7_data_rma_sub_sum_norm_expr,F8_data_rma_sub_sum_norm_expr,
              F9_data_rma_sub_sum_norm_expr,F10_data_rma_sub_sum_norm_expr,F11_data_rma_sub_sum_norm_expr,F12_data_rma_sub_sum_norm_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


# MAS5 background sub only ------------------------------------------------
T_data_mas_sub <- bg.correct(object = Data[,T_samples],method = "mas")
T_data_mas_sub_expr <- exprs(T_data_mas_sub) 

NP_data_mas_sub <- bg.correct(object = Data[,NP],method = "mas")
NP_data_mas_sub_expr <- exprs(NP_data_mas_sub) 

F5_data_mas_sub <- bg.correct(object = Data[,F5],method = "mas")
F5_data_mas_sub_expr <- exprs(F5_data_mas_sub) 

F6_data_mas_sub <- bg.correct(object = Data[,F6],method = "mas")
F6_data_mas_sub_expr <- exprs(F6_data_mas_sub) 

F7_data_mas_sub <- bg.correct(object = Data[,F7],method = "mas")
F7_data_mas_sub_expr <- exprs(F7_data_mas_sub) 

F8_data_mas_sub <- bg.correct(object = Data[,F8],method = "mas")
F8_data_mas_sub_expr <- exprs(F8_data_mas_sub) 

F9_data_mas_sub <- bg.correct(object = Data[,F9],method = "mas")
F9_data_mas_sub_expr <- exprs(F9_data_mas_sub) 

F10_data_mas_sub <- bg.correct(object = Data[,F10],method = "mas")
F10_data_mas_sub_expr <- exprs(F10_data_mas_sub) 

F11_data_mas_sub <- bg.correct(object = Data[,F11],method = "mas")
F11_data_mas_sub_expr <- exprs(F11_data_mas_sub) 

F12_data_mas_sub <- bg.correct(object = Data[,F12],method = "mas")
F12_data_mas_sub_expr <- exprs(F12_data_mas_sub) 

png(filename = "mas_sub_only_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(log2(cbind(T_data_mas_sub_expr,NP_data_mas_sub_expr,F5_data_mas_sub_expr,F6_data_mas_sub_expr,F7_data_mas_sub_expr,F8_data_mas_sub_expr,
                   F9_data_mas_sub_expr,F10_data_mas_sub_expr,F11_data_mas_sub_expr,F12_data_mas_sub_expr)),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


# MAS5 sum only -----------------------------------------------------------
T_data_mas_sum <- expresso(Data[,T_samples],normalize=F,normalize.method="constant",
                                bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                pmcorrect.method = "pmonly")
T_data_mas_sum_expr <- exprs(T_data_mas_sum) 

NP_data_mas_sum <- expresso(Data[,NP],normalize=F,normalize.method="constant",
                                 bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
NP_data_mas_sum_expr <- exprs(NP_data_mas_sum) 

F5_data_mas_sum <- expresso(Data[,F5],normalize=F,normalize.method="constant",
                                 bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F5_data_mas_sum_expr <- exprs(F5_data_mas_sum) 

F6_data_mas_sum <- expresso(Data[,F6],normalize=F,normalize.method="constant",
                                 bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F6_data_mas_sum_expr <- exprs(F6_data_mas_sum) 

F7_data_mas_sum <- expresso(Data[,F7],normalize=F,normalize.method="constant",
                                 bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F7_data_mas_sum_expr <- exprs(F7_data_mas_sum) 

F8_data_mas_sum <- expresso(Data[,F8],normalize=F,normalize.method="constant",
                                 bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F8_data_mas_sum_expr <- exprs(F8_data_mas_sum) 

F9_data_mas_sum <- expresso(Data[,F9],normalize=F,normalize.method="constant",
                                 bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F9_data_mas_sum_expr <- exprs(F9_data_mas_sum) 

F10_data_mas_sum <- expresso(Data[,F10],normalize=F,normalize.method="constant",
                                  bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F10_data_mas_sum_expr <- exprs(F10_data_mas_sum) 

F11_data_mas_sum <- expresso(Data[,F11],normalize=F,normalize.method="constant",
                                  bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F11_data_mas_sum_expr <- exprs(F11_data_mas_sum) 

F12_data_mas_sum <- expresso(Data[,F12],normalize=F,normalize.method="constant",
                                  bg.correct = F,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F12_data_mas_sum_expr <- exprs(F12_data_mas_sum) 

png(filename = "mas_sum_only_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(T_data_mas_sum_expr,NP_data_mas_sum_expr,F5_data_mas_sum_expr,F6_data_mas_sum_expr,F7_data_mas_sum_expr,F8_data_mas_sum_expr,
              F9_data_mas_sum_expr,F10_data_mas_sum_expr,F11_data_mas_sum_expr,F12_data_mas_sum_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()



# MAS5 sub and sum --------------------------------------------------------
T_data_mas_sub_sum <- expresso(Data[,T_samples],normalize=F,normalize.method="constant",
                                bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                pmcorrect.method = "pmonly")
T_data_mas_sub_sum_expr <- exprs(T_data_mas_sub_sum) 

NP_data_mas_sub_sum <- expresso(Data[,NP],normalize=F,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
NP_data_mas_sub_sum_expr <- exprs(NP_data_mas_sub_sum) 

F5_data_mas_sub_sum <- expresso(Data[,F5],normalize=F,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F5_data_mas_sub_sum_expr <- exprs(F5_data_mas_sub_sum) 

F6_data_mas_sub_sum <- expresso(Data[,F6],normalize=F,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F6_data_mas_sub_sum_expr <- exprs(F6_data_mas_sub_sum) 

F7_data_mas_sub_sum <- expresso(Data[,F7],normalize=F,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F7_data_mas_sub_sum_expr <- exprs(F7_data_mas_sub_sum) 

F8_data_mas_sub_sum <- expresso(Data[,F8],normalize=F,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F8_data_mas_sub_sum_expr <- exprs(F8_data_mas_sub_sum) 

F9_data_mas_sub_sum <- expresso(Data[,F9],normalize=F,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F9_data_mas_sub_sum_expr <- exprs(F9_data_mas_sub_sum) 

F10_data_mas_sub_sum <- expresso(Data[,F10],normalize=F,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F10_data_mas_sub_sum_expr <- exprs(F10_data_mas_sub_sum) 

F11_data_mas_sub_sum <- expresso(Data[,F11],normalize=F,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F11_data_mas_sub_sum_expr <- exprs(F11_data_mas_sub_sum) 

F12_data_mas_sub_sum <- expresso(Data[,F12],normalize=F,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F12_data_mas_sub_sum_expr <- exprs(F12_data_mas_sub_sum) 

png(filename = "mas_sub_sum_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(T_data_mas_sub_sum_expr,NP_data_mas_sub_sum_expr,F5_data_mas_sub_sum_expr,F6_data_mas_sub_sum_expr,F7_data_mas_sub_sum_expr,F8_data_mas_sub_sum_expr,
              F9_data_mas_sub_sum_expr,F10_data_mas_sub_sum_expr,F11_data_mas_sub_sum_expr,F12_data_mas_sub_sum_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()


# MAS5 sub and sum with constant norm pre sum intra fraction norm only -----------------------------
T_data_mas_sub_sum_norm <- expresso(Data[,T_samples],normalize=T,normalize.method="constant",
                                bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                pmcorrect.method = "pmonly")
T_data_mas_sub_sum_norm_expr <- exprs(T_data_mas_sub_sum_norm) 

NP_data_mas_sub_sum_norm <- expresso(Data[,NP],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
NP_data_mas_sub_sum_norm_expr <- exprs(NP_data_mas_sub_sum_norm) 

F5_data_mas_sub_sum_norm <- expresso(Data[,F5],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F5_data_mas_sub_sum_norm_expr <- exprs(F5_data_mas_sub_sum_norm) 

F6_data_mas_sub_sum_norm <- expresso(Data[,F6],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F6_data_mas_sub_sum_norm_expr <- exprs(F6_data_mas_sub_sum_norm) 

F7_data_mas_sub_sum_norm <- expresso(Data[,F7],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F7_data_mas_sub_sum_norm_expr <- exprs(F7_data_mas_sub_sum_norm) 

F8_data_mas_sub_sum_norm <- expresso(Data[,F8],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F8_data_mas_sub_sum_norm_expr <- exprs(F8_data_mas_sub_sum_norm) 

F9_data_mas_sub_sum_norm <- expresso(Data[,F9],normalize=T,normalize.method="constant",
                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                 pmcorrect.method = "pmonly")
F9_data_mas_sub_sum_norm_expr <- exprs(F9_data_mas_sub_sum_norm) 

F10_data_mas_sub_sum_norm <- expresso(Data[,F10],normalize=T,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F10_data_mas_sub_sum_norm_expr <- exprs(F10_data_mas_sub_sum_norm) 

F11_data_mas_sub_sum_norm <- expresso(Data[,F11],normalize=T,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F11_data_mas_sub_sum_norm_expr <- exprs(F11_data_mas_sub_sum_norm) 

F12_data_mas_sub_sum_norm <- expresso(Data[,F12],normalize=T,normalize.method="constant",
                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
                                  pmcorrect.method = "pmonly")
F12_data_mas_sub_sum_norm_expr <- exprs(F12_data_mas_sub_sum_norm) 

png(filename = "mas_sub_sum_norm_Data_boxplot_rep_order.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(cbind(T_data_mas_sub_sum_norm_expr,NP_data_mas_sub_sum_norm_expr,F5_data_mas_sub_sum_norm_expr,F6_data_mas_sub_sum_norm_expr,F7_data_mas_sub_sum_norm_expr,F8_data_mas_sub_sum_norm_expr,
              F9_data_mas_sub_sum_norm_expr,F10_data_mas_sub_sum_norm_expr,F11_data_mas_sub_sum_norm_expr,F12_data_mas_sub_sum_norm_expr),las=2,xaxt='n',yaxt='n',ann=F, main="Raw Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data))),labels = Data$sample[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)],cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()



# MAS5 sub and sum with constant norm post sum intra fraction norm only  ---------------------------
# T_data_mas_sum_norm <- expresso(Data[,T_samples],normalize=T,normalize.method="constant",
#                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                 pmcorrect.method = "pmonly")
# T_data_mas_sum_norm_expr <- exprs(T_data_mas_sum_norm) 
# 
# NP_data_mas_sum_norm <- expresso(Data[,NP],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# NP_data_mas_sum_norm_expr <- exprs(NP_data_mas_sum_norm) 
# 
# F5_data_mas_sum_norm <- expresso(Data[,F5],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F5_data_mas_sum_norm_expr <- exprs(F5_data_mas_sum_norm) 
# 
# F6_data_mas_sum_norm <- expresso(Data[,F6],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F6_data_mas_sum_norm_expr <- exprs(F6_data_mas_sum_norm) 
# 
# F7_data_mas_sum_norm <- expresso(Data[,F7],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F7_data_mas_sum_norm_expr <- exprs(F7_data_mas_sum_norm) 
# 
# F8_data_mas_sum_norm <- expresso(Data[,F8],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F8_data_mas_sum_norm_expr <- exprs(F8_data_mas_sum_norm) 
# 
# F9_data_mas_sum_norm <- expresso(Data[,F9],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F9_data_mas_sum_norm_expr <- exprs(F9_data_mas_sum_norm) 
# 
# F10_data_mas_sum_norm <- expresso(Data[,F10],normalize=T,normalize.method="constant",
#                                   bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                   pmcorrect.method = "pmonly")
# F10_data_mas_sum_norm_expr <- exprs(F10_data_mas_sum_norm) 
# 
# F11_data_mas_sum_norm <- expresso(Data[,F11],normalize=T,normalize.method="constant",
#                                   bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                   pmcorrect.method = "pmonly")
# F11_data_mas_sum_norm_expr <- exprs(F11_data_mas_sum_norm) 
# 
# F12_data_mas_sum_norm <- expresso(Data[,F12],normalize=T,normalize.method="constant",
#                                   bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                   pmcorrect.method = "pmonly")
# F12_data_mas_sum_norm_expr <- exprs(F12_data_mas_sum_norm) 
# 
# 
# # MAS5 sub and sum with constant norm pre sum and post sum inter fraction norm only ----------------
# T_data_mas_sum_norm <- expresso(Data[,T_samples],normalize=T,normalize.method="constant",
#                                 bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                 pmcorrect.method = "pmonly")
# T_data_mas_sum_norm_expr <- exprs(T_data_mas_sum_norm) 
# 
# NP_data_mas_sum_norm <- expresso(Data[,NP],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# NP_data_mas_sum_norm_expr <- exprs(NP_data_mas_sum_norm) 
# 
# F5_data_mas_sum_norm <- expresso(Data[,F5],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F5_data_mas_sum_norm_expr <- exprs(F5_data_mas_sum_norm) 
# 
# F6_data_mas_sum_norm <- expresso(Data[,F6],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F6_data_mas_sum_norm_expr <- exprs(F6_data_mas_sum_norm) 
# 
# F7_data_mas_sum_norm <- expresso(Data[,F7],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F7_data_mas_sum_norm_expr <- exprs(F7_data_mas_sum_norm) 
# 
# F8_data_mas_sum_norm <- expresso(Data[,F8],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F8_data_mas_sum_norm_expr <- exprs(F8_data_mas_sum_norm) 
# 
# F9_data_mas_sum_norm <- expresso(Data[,F9],normalize=T,normalize.method="constant",
#                                  bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                  pmcorrect.method = "pmonly")
# F9_data_mas_sum_norm_expr <- exprs(F9_data_mas_sum_norm) 
# 
# F10_data_mas_sum_norm <- expresso(Data[,F10],normalize=T,normalize.method="constant",
#                                   bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                   pmcorrect.method = "pmonly")
# F10_data_mas_sum_norm_expr <- exprs(F10_data_mas_sum_norm) 
# 
# F11_data_mas_sum_norm <- expresso(Data[,F11],normalize=T,normalize.method="constant",
#                                   bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                   pmcorrect.method = "pmonly")
# F11_data_mas_sum_norm_expr <- exprs(F11_data_mas_sum_norm) 
# 
# F12_data_mas_sum_norm <- expresso(Data[,F12],normalize=T,normalize.method="constant",
#                                   bg.correct = T,bgcorrect.method = "mas",summary.method = "mas",
#                                   pmcorrect.method = "pmonly")
# F12_data_mas_sum_norm_expr <- exprs(F12_data_mas_sum_norm) 
# 
# 






#Comparing RMA methods
All_data_rma_expr[1:10,grep(colnames(All_data_rma_expr),pattern = "T")]
T_data_rma_expr[1:10,]

#comparing rma to GCRMA

T_data_gcrma_expr[1:10,]
T_data_rma_expr[1:10,]

T_data_gcrma_expr[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),][1:10,]
T_data_rma_expr[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),][1:10,]

T_rma_sd <- hist(apply(T_data_rma_expr,1,sd),breaks=seq(0,3.5,0.1))
T_gcrma_sd <- hist(apply(T_data_gcrma_expr,1,sd),breaks=seq(0,3.5,0.1))

barplot(T_gcrma_sd$counts)
barplot(T_rma_sd$counts,add=T,col=rgb(1,0,0,0.5))


#Reconstituting TX from TL

All_data_gcrma_expr_R1 <- All_data_gcrma_expr[,grep(colnames(All_data_gcrma_expr),pattern = "R1")]
All_data_gcrma_expr_R2 <- All_data_gcrma_expr[,grep(colnames(All_data_gcrma_expr),pattern = "R2")]
All_data_gcrma_expr_R3 <- All_data_gcrma_expr[,grep(colnames(All_data_gcrma_expr),pattern = "R3")]

TL_data_gcrma_expr_R1 <- TL_data_gcrma_expr[,grep(colnames(TL_data_gcrma_expr),pattern = "R1")]
TL_data_gcrma_expr_R2 <- TL_data_gcrma_expr[,grep(colnames(TL_data_gcrma_expr),pattern = "R2")]
TL_data_gcrma_expr_R3 <- TL_data_gcrma_expr[,grep(colnames(TL_data_gcrma_expr),pattern = "R3")]

SEP_data_gcrma_expr_R1 <- cbind(NP_data_gcrma_expr[,1],F5_data_gcrma_expr[,1],F6_data_gcrma_expr[,1],
                                F7_data_gcrma_expr[,1],F8_data_gcrma_expr[,1],F9_data_gcrma_expr[,1],
                                F10_data_gcrma_expr[,1],F11_data_gcrma_expr[,1],F12_data_gcrma_expr[,1])
SEP_data_gcrma_expr_R2 <- cbind(NP_data_gcrma_expr[,2],F5_data_gcrma_expr[,2],F6_data_gcrma_expr[,2],
                                F7_data_gcrma_expr[,2],F8_data_gcrma_expr[,2],F9_data_gcrma_expr[,2],
                                F10_data_gcrma_expr[,2],F11_data_gcrma_expr[,2],F12_data_gcrma_expr[,2])
SEP_data_gcrma_expr_R3 <- cbind(NP_data_gcrma_expr[,3],F5_data_gcrma_expr[,3],F6_data_gcrma_expr[,3],
                                F7_data_gcrma_expr[,3],F8_data_gcrma_expr[,3],F9_data_gcrma_expr[,3],
                                F10_data_gcrma_expr[,3],F11_data_gcrma_expr[,3],F12_data_gcrma_expr[,3])


rep1_dilution_factors <- c(9.25,10.4,8.3,8.46,7.82,5.35,1.51,1.22,1)
rep1_dilution_factors <- c(.41,.46,.37,.67,.35,.24,.07,.05,0.04)

dilution_factors_wfullNP <- c(0.312400421, 0.171593299, 0.13446095, 0.124347128, 0.101602326, 0.074285634, 0.037868024, 0.025067032, 0.018375186)
dilution_factors_wNP_2 <- c(0.185233544, 0.203442948, 0.159383444, 0.147252893, 0.120144984, 0.087913686, 0.045044533, 0.029776716, 0.021807253)


#expand.grid(1,2,3)
#?expand.grid
library(gtools)
df_perm <- permutations(n=9,r=9,v=dilution_factors_wfullNP)
df_perm2 <- df_perm[sample(1:nrow(df_perm),size = 10000,replace = F),]
#df_perm2 <- rbind(rev(dilution_factors_wfullNP),dilution_factors_wfullNP)

out.cor.p <- numeric(nrow(df_perm2))
out.cor.s <- numeric(nrow(df_perm2))
for(i in 1:nrow(df_perm2)){
  if(i%%1000==0){
    print(i)
  }
  x1 <- All_data_gcrma_expr_R1[-which(rownames(All_data_gcrma_expr_R1)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(df_perm2[i,],2),length(Filtered_data_TAIR)),nrow=9))
  x1 <- log(rowSums(2^x1),2)
  x2 <- All_data_gcrma_expr_R2[-which(rownames(All_data_gcrma_expr_R2)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(df_perm2[i,],2),length(Filtered_data_TAIR)),nrow=9))
  x2 <- log(rowSums(2^x2),2)
  x3 <- All_data_gcrma_expr_R3[-which(rownames(All_data_gcrma_expr_R3)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(df_perm2[i,],2),length(Filtered_data_TAIR)),nrow=9))
  x3 <- log(rowSums(2^x3),2)
  xm <- matrix(c(x1,x2,x3),ncol=3)
  rownames(xm) <- names(x3)
  x <- log(rowMeans(2^xm),2)
  
  y <- cbind(All_data_gcrma_expr_R1[,1],All_data_gcrma_expr_R2[,1],All_data_gcrma_expr_R3[,1])
  y <- log(rowMeans(2^y),2)
  y <- y[-which(names(y)%in%names(Filtered_data_TAIR)==F)]
  
  #z <- log(y-x,2)
  z <- y-x
  #z<-z[-which(names(z)%in%names(Filtered_data_TAIR)==F)]
  #zul <- 2^y-2^x
  
  #hist(zul/2^y)
  
  out.cor.p[i] <- cor(sort(y),unlist(x[order(y)]),method = "pearson")
  out.cor.s[i] <- cor(sort(y),unlist(x[order(y)]),method = "spearman")
  #plot(sort(y),z[order(y)])
  #plot(sort(y),x[order(y)],main="with full NP")
  #text(x =3 ,y=15,labels = paste("cor= ",round(cor(sort(y),unlist(x[order(y)])),4),sep=""))
  #abline(a=c(0,1),col="red")
  
}

hist(out.cor.p)
hist(out.cor.s)

x <- All_data_gcrma_expr_R1[,2:10]+t(matrix(rep(log(dilution_factors_wfullNP,2),nrow(All_data_gcrma)),nrow=9))
x <- log(rowSums(2^x),2)
x1 <- All_data_gcrma_expr_R1[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x1 <- log(rowSums(2^x1),2)
x2 <- All_data_gcrma_expr_R2[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x2 <- log(rowSums(2^x2),2)
x3 <- All_data_gcrma_expr_R3[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x3 <- log(rowSums(2^x3),2)
xm <- matrix(c(x1,x2,x3),ncol=3)
rownames(xm) <- names(x3)
x <- log(rowMeans(2^xm),2)

x<- rowSums(2^(All_data_gcrma_expr_R1[,2:10])*rep1_dilution_factors)/2.77
x<- rowSums(All_data_gcrma_expr_R1[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]*rep1_dilution_factors)/2.77
x<- rowSums(TL_data_gcrma_expr_R1*rep1_dilution_factors)/2.77

y <- All_data_gcrma_expr_R1[,1]
y <- All_data_gcrma_expr_R1[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),1]
y <- T_data_gcrma_expr[,1]
y <- cbind(All_data_gcrma_expr_R1[,1],All_data_gcrma_expr_R2[,1],All_data_gcrma_expr_R3[,1])
y <- log(rowMeans(2^y),2)
y <- y[-which(names(y)%in%names(Filtered_data_TAIR)==F)]

#z <- log(y-x,2)
z <- y-x
#z<-z[-which(names(z)%in%names(Filtered_data_TAIR)==F)]
zul <- 2^y-2^x

hist(zul/2^y)

plot(sort(y),z[order(y)])
plot(sort(y),x[order(y)],main="with full NP")
text(x =3 ,y=15,labels = paste("cor= ",round(cor(sort(y),unlist(x[order(y)])),4),sep=""))
abline(a=c(0,1),col="red")
hist(out.cor.p)
abline(v=cor(sort(y),unlist(x[order(y)]),method="pearson"),col="red")
hist(out.cor.s)
abline(v=cor(sort(y),unlist(x[order(y)]),method="spearman"),col="red")

plot(sort(y),zul[order(y)])
plot(sort(y),zul[order(y)]/2^sort(y))

a<- hist(z,breaks=seq(-5,5,0.1),main="with full NP")
ac <-a$counts
names(ac) <- a$mids

length(which((z)>1))
length(which((z)>2))
length(which((z)>3))
length(which((z)< -1))
length(which((z)< -2))

#using sep
x1 <- SEP_data_gcrma_expr_R1[-which(rownames(SEP_data_gcrma_expr_R1)%in%names(Filtered_data_TAIR)==F),]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x1 <- log(rowSums(2^x1),2)
x2 <- SEP_data_gcrma_expr_R2[-which(rownames(SEP_data_gcrma_expr_R2)%in%names(Filtered_data_TAIR)==F),]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x2 <- log(rowSums(2^x2),2)
x3 <- SEP_data_gcrma_expr_R3[-which(rownames(SEP_data_gcrma_expr_R3)%in%names(Filtered_data_TAIR)==F),]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x3 <- log(rowSums(2^x3),2)
xm <- matrix(c(x1,x2,x3),ncol=3)
rownames(xm) <- names(x3)
x <- log(rowMeans(2^xm),2)

y <- T_data_gcrma_expr
y <- log(rowMeans(2^y),2)
y <- y[-which(names(y)%in%names(Filtered_data_TAIR)==F)]

#z <- log(y-x,2)
z <- y-x
#z<-z[-which(names(z)%in%names(Filtered_data_TAIR)==F)]
#zul <- 2^y-2^x

#hist(zul/2^y)

#out.cor[i] <- cor(sort(y),unlist(x[order(y)]))
#plot(sort(y),z[order(y)])
#plot(sort(y),x[order(y)],main="with full NP")
#text(x =3 ,y=15,labels = paste("cor= ",round(cor(sort(y),unlist(x[order(y)])),4),sep=""))
#abline(a=c(0,1),col="red")

a<- hist(z,breaks=seq(-10,10,0.1),main="with half NP")
ac <-a$counts
names(ac) <- a$mids

length(which((z)>1))
length(which((z)>2))
length(which((z)>3))
length(which((z)< -1))
length(which((z)< -2))




# using half np
x <- All_data_gcrma_expr_R1[,2:10]+t(matrix(rep(log(dilution_factors_wfullNP,2),nrow(All_data_gcrma)),nrow=9))
x <- log(rowSums(2^x),2)
x1 <- All_data_gcrma_expr_R1[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(dilution_factors_wNP_2,2),length(Filtered_data_TAIR)),nrow=9))
x1 <- log(rowSums(2^x1),2)
x2 <- All_data_gcrma_expr_R2[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(dilution_factors_wNP_2,2),length(Filtered_data_TAIR)),nrow=9))
x2 <- log(rowSums(2^x2),2)
x3 <- All_data_gcrma_expr_R3[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+t(matrix(rep(log(dilution_factors_wNP_2,2),length(Filtered_data_TAIR)),nrow=9))
x3 <- log(rowSums(2^x3),2)
xm <- matrix(c(x1,x2,x3),ncol=3)
rownames(xm) <- names(x3)
x <- log(rowMeans(2^xm),2)

y <- All_data_gcrma_expr_R1[,1]
y <- All_data_gcrma_expr_R1[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),1]
y <- T_data_gcrma_expr[,1]
y <- cbind(All_data_gcrma_expr_R1[,1],All_data_gcrma_expr_R2[,1],All_data_gcrma_expr_R3[,1])
y <- log(rowMeans(2^y),2)
y <- y[-which(names(y)%in%names(Filtered_data_TAIR)==F)]

#z <- log(y-x,2)
z <- y-x
#z<-z[-which(names(z)%in%names(Filtered_data_TAIR)==F)]
zul <- 2^y-2^x

hist(zul/2^y)

plot(sort(y),z[order(y)])
plot(sort(y),x[order(y)],main="with half NP")
text(x =3 ,y=15,labels = paste("cor= ",round(cor(sort(y),unlist(x[order(y)])),4),sep=""))
abline(a=c(0,1),col="red")

plot(sort(y),zul[order(y)])
plot(sort(y),zul[order(y)]/2^sort(y))

a<- hist(z,breaks=seq(-5,5,0.1),main="with half NP")
ac <-a$counts
names(ac) <- a$mids

length(which((z)>1))
length(which((z)>2))
length(which((z)>3))
length(which((z)< -1))
length(which((z)< -2))












#Using separate TX and TL

x <- TL_data_gcrma_expr_R1[,2:10]+t(matrix(rep(log(dilution_factors_wfullNP,2),nrow(TL_data_gcrma)),nrow=9))
x <- log(rowSums(2^x),2)
x1 <- TL_data_gcrma_expr_R1[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x1 <- log(rowSums(2^x1),2)
x2 <- TL_data_gcrma_expr_R2[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x2 <- log(rowSums(2^x2),2)
x3 <- TL_data_gcrma_expr_R3[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),]+t(matrix(rep(log(dilution_factors_wfullNP,2),length(Filtered_data_TAIR)),nrow=9))
x3 <- log(rowSums(2^x3),2)
xm <- matrix(c(x1,x2,x3),ncol=3)
rownames(xm) <- names(x3)
x <- log(rowMeans(2^xm),2)

x<- rowSums(2^(TL_data_gcrma_expr_R1[,2:10])*rep1_dilution_factors)/2.77
x<- rowSums(TL_data_gcrma_expr_R1[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]*rep1_dilution_factors)/2.77
x<- rowSums(TL_data_gcrma_expr_R1*rep1_dilution_factors)/2.77

y <- T_data_gcrma_expr[,1]
y <- T_data_gcrma_expr
y <- log(rowMeans(2^y),2)
y <- y[-which(names(y)%in%names(Filtered_data_TAIR)==F)]

#z <- log(y-x,2)
z <- y-x
#z<-z[-which(names(z)%in%names(Filtered_data_TAIR)==F)]
zul <- 2^y-2^x
plot(2^z,zul)

hist(zul/2^y)

plot(sort(y),z[order(y)])
plot(sort(y),x[order(y)])

plot(sort(y),zul[order(y)])
plot(sort(y),zul[order(y)]/2^sort(y))

a<- hist(z,breaks=seq(-10,10,0.1))
ac <-a$counts
names(ac) <- a$mids

length(which((z)>1))
length(which((z)>2))
length(which((z)>3))
length(which((z)< -1))
length(which((z)< -2))




j=10001
plot(All_data_gcrma_expr_R1[j,2:10]*rep1_dilution_factors,ylim=c(0,15))
lines(All_data_gcrma_expr_R1[j,2:10]*rep1_dilution_factors)
for(i in j:(j+9)){
  points(All_data_gcrma_expr_R1[i,2:10]*rep1_dilution_factors,pch=1:10)
  lines(All_data_gcrma_expr_R1[i,2:10]*rep1_dilution_factors,lty=1:10)
}



names(which((y-x)>2))










z <- t(apply(All_data_gcrma_expr_R1[,2:10]*rep1_dilution_factors,1,order,decreasing=T))

All_data_gcrma_expr_R1[1:10,2:10]*rep1_dilution_factors

rep3_dilution_factors <- c(7.9*2,8.5,6.5,5.5,3.8,3,2.5,1.5,1)
rep3_dilution_factors <- c(0.88,0.95,.73,.61,.42,.34,.28,.16,.11)

x<- rowSums(All_data_gcrma_expr_R3[,2:10]*rep3_dilution_factors)/5.37

y <- All_data_gcrma_expr_R3[,1]

a<- hist(y-x,breaks=seq(-15,15,0.1))
ac <-a$counts
names(ac) <- a$mids

length(which((y-x)>1))
length(which((y-x)>2))
length(which((y-x)>3))
length(which((y-x)< -1))
length(which((y-x)< -2))

j=1
plot(All_data_gcrma_expr_R3[j,2:10]+log(dilution_factors_wfullNP,2),ylim=c(-15,15))
lines(All_data_gcrma_expr_R3[j,2:10]+log(dilution_factors_wfullNP,2))
lines(dilution_factors_wfullNP*10,col="red")
for(i in j:(j+9)){
  points(All_data_gcrma_expr_R3[i,2:10]+log(dilution_factors_wfullNP,2),pch=i)
  lines(All_data_gcrma_expr_R3[i,2:10]+log(dilution_factors_wfullNP,2),lty=i)
}

plot(dilution_factors_wfullNP*10,col="red")


m <- cbind(All_data_gcrma_expr_R1[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+log(dilution_factors_wfullNP,2),
           All_data_gcrma_expr_R2[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+log(dilution_factors_wfullNP,2),
           All_data_gcrma_expr_R3[-which(rownames(Data)%in%names(Filtered_data_TAIR)==F),2:10]+log(dilution_factors_wfullNP,2))


heatmap.2(m,Rowv = "none",Colv = "none",dendrogram = "none",scale = "none",tra)



j=1
plot(All_data_gcrma_expr_R3[j,2:10]*dilution_factors_wNP_2,ylim=c(0,15))
lines(All_data_gcrma_expr_R3[j,2:10]*dilution_factors_wNP_2)
for(i in j:(j+9)){
  points(All_data_gcrma_expr_R3[i,2:10]*dilution_factors_wNP_2,pch=i)
  lines(All_data_gcrma_expr_R3[i,2:10]*dilution_factors_wNP_2,lty=i)
}



z <- t(apply(All_data_gcrma_expr_R3[,2:10]*rep3_dilution_factors,1,order,decreasing=T))

All_data_gcrma_expr_R3[1:10,2:10]*rep3_dilution_factors







#Mas5 TL
col_m5calls_TL <- mas5calls(Data_TL_col)
col_m5calls_TL <- as.matrix(exprs(col_m5calls_TL))

old.dir <- getwd()
setwd(paste(old.dir,"TX/",sep="/"))
colnames(col_m5calls_TL) <- sub(".CEL$","",sub("_AVA\\d{3}_","_",colnames(col_m5calls_TL)))
write.csv(col_m5calls_TL[,c(IND_NP,IND_5,IND_6,IND_7,IND_8,IND_9,IND_10,IND_11,IND_12)],"col_M5callsTX.csv")
setwd(old.dir)


TL_NP <- apply(col_m5calls_TL[,1:3],1,function(i) all(i=="P")==T)

index_TL <- which(TL_NP==T)



#RMA TX

col_tx_rma <- rma(Data_TX_col)
col_tx_rma <- exprs(col_tx_rma)

old.dir <- getwd()
setwd(paste(old.dir,"TX/rma",sep="/"))
png(filename = "RMA TX Data boxplot.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(col_tx_rma,las=2,xaxt='n',yaxt='n',ann=F, main="RMA TX Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data_TX_col))),labels = sub(".CEL$","",sub("AVA\\d{3}_","",colnames(Data_TX_col))),cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()

setwd(old.dir)

old.dir <- getwd()
setwd(paste(old.dir,"TX/rma",sep="/"))
png(filename = "TX RMA MDS.png" ,width = 2500,height = 2000,units = "px",res=300)  
plotMDS(col_tx_rma[index_TX,],col=c(rep("blue",1),rep("green3",1),rep("red",1)),main = "TX RMA MDS")
dev.off()
setwd(old.dir)



#GCRMA TX
col_tx_gcrma <- gcrma(Data_TX_col)
col_tx_gcrma <- exprs(col_tx_gcrma)

setwd(old.dir)

old.dir <- getwd()
setwd(paste(old.dir,"TX/gcrma",sep="/"))
png(filename = "TX GCRMA MDS.png" ,width = 2500,height = 2000,units = "px",res=300)  
plotMDS(col_tx_gcrma[index_TX,],col=c(rep("blue",1),rep("green3",1),rep("red",1)),main = "TX GCRMA MDS")
dev.off()
setwd(old.dir)



#RMA TL
col_tl_rma <- rma(Data_TL_col)
col_tl_rma <- exprs(col_tl_rma)


colnames(col_tl_rma) <- sub(".CEL$","",sub("AVA\\d{3}_","",colnames(col_tl_rma)))


old.dir <- getwd()
setwd(paste(old.dir,"TL/rma",sep="/"))
png(filename = "RMA TL Data boxplot.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(col_tl_rma,las=2,xaxt='n',yaxt='n',ann=F, main="RMA TL Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data_TL_col))),labels = sub(".CEL$","",sub("AVA\\d{3}_","",colnames(Data_TL_col))),cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()
setwd(old.dir)


col_tl_rma <- col_tl_rma[,c(IND_NP,IND_5,IND_6,IND_7,IND_8,IND_9,IND_10,IND_11,IND_12)]

colnames(col_tl_rma)

old.dir <- getwd()
setwd(paste(old.dir,"TL/rma",sep="/"))
png(filename = "TL RMA MDS.png" ,width = 2500,height = 2000,units = "px",res=300)  
plotMDS(col_tl_rma[index_TL,],col=c(rep("blue",3),rep("green3",3),rep("red",3),rep("lightcyan3",3),rep("lightcoral",3),rep("black",3),rep("forestgreen",3),rep("aquamarine4",3),rep("bisque4",3),rep("orange",3)),main = "TL RMA MDS",top=10000)
dev.off()
setwd(old.dir)

#Remove Batch effects from data using limma

batch.factor <- rep(c("R1","R2","R3"),9)
batch.factor <- factor(batch.factor,levels=unique(batch.factor))
condition.factor <-  c(rep("NP",3),rep("5",3),rep("6",3),rep("7",3),rep("8",3),rep("9",3),rep("10",3),rep("11",3),rep("12",3))
condition.factor <- factor(condition.factor,levels=unique(condition.factor))

design.1 <- model.matrix(~ 0 + condition.factor)
rownames(design.1) <- colnames(col_tl_rma)
col.tl.rma.rm.batch <- removeBatchEffect(col_tl_rma, batch.factor, design=design.1)


old.dir <- getwd()
setwd(paste(old.dir,"TL/rma",sep="/"))
png(filename = "TL_rma_batch_removed_MDS.png" ,width = 2500,height = 2000,units = "px",res=300)  
plotMDS(col.tl.rma.rm.batch,col=c(rep("blue",3),rep("green3",3),rep("red",3),rep("lightcyan3",3),rep("lightcoral",3),rep("black",3),rep("forestgreen",3),rep("aquamarine4",3),rep("bisque4",3),rep("orange",3)),main = "TL RMA MDS",top=10000)
dev.off()
setwd(old.dir)



#GCRMA TL
col_tl_gcrma <- gcrma(Data_TL_col)
col_tl_gcrma <- exprs(col_tl_gcrma)

colnames(col_tl_gcrma) <- sub(".CEL$","",sub("AVA\\d{3}_","",colnames(col_tl_gcrma)))

old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma",sep="/"))
png(filename = "GCRMA TL Data boxplot.png" ,width = 2500,height = 2000,units = "px",res=300)  
par(mar=c(8.1, 5.1, 4.1, 2.1))
boxplot(col_tl_gcrma,las=2,xaxt='n',yaxt='n',ann=F, main="Raw TL Data box plot",ylab=expression(paste('log'["2"]," expression")),cex.lab=1.2,font.lab=2)
axis(side = 1,tick = T,las=2,at = c(1:length(colnames(Data_TL_col))),labels = sub(".CEL$","",sub("AVA\\d{3}_","",colnames(Data_TL_col))),cex=2,font=2)
axis(side = 2,tick = T,las=2,at = seq(2,14,2),labels = seq(2,14,2),cex=2,font=2)
dev.off()

setwd(old.dir)

col_tl_gcrma <- col_tl_gcrma[,c(IND_NP,IND_5,IND_6,IND_7,IND_8,IND_9,IND_10,IND_11,IND_12)]

colnames(col_tl_gcrma)

old.dir <- getwd()
setwd(paste(old.dir,"TL/rma",sep="/"))
png(filename = "TL RMA MDS.png" ,width = 2500,height = 2000,units = "px",res=300)  
plotMDS(col_tl_rma[index_TL,],col=c(rep("blue",3),rep("green3",3),rep("red",3),rep("lightcyan3",3),rep("lightcoral",3),rep("black",3),rep("forestgreen",3),rep("aquamarine4",3),rep("bisque4",3),rep("orange",3)),main = "TL RMA MDS")
dev.off()
setwd(old.dir)


#Remove Batch effects from data using limma

batch.factor <- rep(c("R1","R2","R3"),9)
batch.factor <- factor(batch.factor,levels=unique(batch.factor))
condition.factor <-  c(rep("NP",3),rep("5",3),rep("6",3),rep("7",3),rep("8",3),rep("9",3),rep("10",3),rep("11",3),rep("12",3))
condition.factor <- factor(condition.factor,levels=unique(condition.factor))

design.1 <- model.matrix(~ 0 + condition.factor)
rownames(design.1) <- colnames(col_tl_gcrma)
col.tl.gcrma.rm.batch <- removeBatchEffect(col_tl_gcrma, batch.factor, design=design.1)


old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma",sep="/"))
png(filename = "TL_gcrma_batch_removed_MDS.png" ,width = 2500,height = 2000,units = "px",res=300)  
plotMDS(col.tl.gcrma.rm.batch,col=c(rep("blue",3),rep("green3",3),rep("red",3),rep("lightcyan3",3),rep("lightcoral",3),rep("black",3),rep("forestgreen",3),rep("aquamarine4",3),rep("bisque4",3),rep("orange",3)),main = "TL GCRMA MDS")
dev.off()
setwd(old.dir)


#HEATMAPS
#RMA
cor_col_TL_rma_col <- cor(col_tl_rma)
cor_col_TL_rma_row <- cor(t(col_tl_rma))


which(is.na(cor_TX_expr_gcrma_row[1,]))
x <- which(apply(cor_TX_expr_gcrma_row,1,function(i) any(is.na(i))))

TL_rma_hc_col <- hclust(as.dist(1-cor_col_TL_rma_col))
TL_rma_hc_row <- hclust(as.dist(1-cor_col_TL_rma_row))

plot(TX_expr_gcrma_hc_col)
plot(TX_expr_gcrma_hc_row)

old.dir <- getwd()
setwd(paste(old.dir,"TL/rma/",sep="/"))
png(filename = "TL_rma_row_all_heatmap.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_rma,Rowv = as.dendrogram(TL_rma_hc_row), Colv = as.dendrogram(TL_rma_hc_col),trace="none",scale="row")
dev.off()
setwd(old.dir)


old.dir <- getwd()
setwd(paste(old.dir,"TL/rma/",sep="/"))
png(filename = "TL_rma_row_all_heatmap_no_col_dend.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_rma,Rowv = as.dendrogram(TL_rma_hc_row), Colv = 'none',trace="none",scale="row")
dev.off()
setwd(old.dir)


rm(cor_col_TL_rma_col)
rm(cor_col_TL_rma_row)
rm(TL_rma_hc_col)
rm(TL_rma_hc_row)


#ave heatmap

col_tl_rma_mean <- matrix(numeric(),ncol=ncol(col_tl_rma)/3,nrow=nrow(col_tl_rma))

for(i in 1:9){
  
  col_tl_rma_mean[,i] <- apply(col_tl_rma[,seq((3*i-2),3*i,1)],1,mean)
  
}

colnames(col_tl_rma_mean) <- c("NP","F5","F6","F7","F8","F9","F10","F11","F12")


cor_col_TL_rma_mean_col <- cor(col_tl_rma_mean)
cor_col_TL_rma_mean_row <- cor(t(col_tl_rma_mean))


TL_rma_mean_hc_col <- hclust(as.dist(1-cor_col_TL_rma_mean_col))
TL_rma_mean_hc_row <- hclust(as.dist(1-cor_col_TL_rma_mean_row))

plot(TX_expr_gcrma_hc_col)
plot(TX_expr_gcrma_hc_row)

old.dir <- getwd()
setwd(paste(old.dir,"TL/rma/",sep="/"))
png(filename = "TL_rma_mean_row_all_heatmap.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_rma_mean,Rowv = as.dendrogram(TL_rma_mean_hc_row), Colv = as.dendrogram(TL_rma_mean_hc_col),trace="none",scale="row")
dev.off()
setwd(old.dir)


old.dir <- getwd()
setwd(paste(old.dir,"TL/rma/",sep="/"))
png(filename = "TL_rma_mean_row_all_heatmap_no_col_dend.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_rma_mean,Rowv = as.dendrogram(TL_rma_mean_hc_row), Colv = 'none',trace="none",scale="row")
dev.off()
setwd(old.dir)


rm(cor_col_TL_rma_mean_col)
rm(cor_col_TL_rma_mean_row)
rm(TL_rma_mean_hc_col)
rm(TL_rma_mean_hc_row)


#RMA batch removed

index_TL_RMA_sd <- apply(col.tl.gcrma.rm.batch,1,sd)
index_TL_RMA_sd <- which(index_TL_RMA_sd==0)


cor_col_TL_rma_batch_col <- cor(col.tl.gcrma.rm.batch[-index_TL_RMA_sd,])
cor_col_TL_rma_batch_row <- cor(t(col.tl.gcrma.rm.batch[-index_TL_RMA_sd,]))

TL_rma_batch_hc_col <- hclust(as.dist(1-cor_col_TL_rma_batch_col))
TL_rma_batch_hc_row <- hclust(as.dist(1-cor_col_TL_rma_batch_row))

plot(TX_expr_gcrma_hc_col)
plot(TX_expr_gcrma_hc_row)

old.dir <- getwd()
setwd(paste(old.dir,"TL/rma/",sep="/"))
png(filename = "TL_rma_batch_row_all_heatmap.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col.tl.gcrma.rm.batch[-index_TL_RMA_sd,],Rowv = as.dendrogram(TL_rma_batch_hc_row), Colv = as.dendrogram(TL_rma_batch_hc_col),trace="none",scale="row")
dev.off()
setwd(old.dir)


old.dir <- getwd()
setwd(paste(old.dir,"TL/rma/",sep="/"))
png(filename = "TL_rma_batch_row_all_heatmap_no_col_dend.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col.tl.gcrma.rm.batch[-index_TL_RMA_sd,],Rowv = as.dendrogram(TL_rma_batch_hc_row), Colv = 'none',trace="none",scale="row")
dev.off()
setwd(old.dir)


rm(cor_col_TL_rma_batch_col)
rm(cor_col_TL_rma_batch_row)
rm(TL_rma_batch_hc_col)
rm(TL_rma_batch_hc_row)


#ave heatmap

col_tl_rma_batch_mean <- matrix(numeric(),ncol=ncol(col.tl.gcrma.rm.batch)/3,nrow=nrow(col.tl.gcrma.rm.batch))

for(i in 1:9){
  
  col_tl_rma_batch_mean[,i] <- apply(col.tl.gcrma.rm.batch[,seq((3*i-2),3*i,1)],1,mean)
  
}

colnames(col_tl_rma_batch_mean) <- c("NP","F5","F6","F7","F8","F9","F10","F11","F12")


index_TL_RMA_mean_sd <- apply(col_tl_rma_batch_mean,1,sd)
index_TL_RMA_mean_sd <- which(index_TL_RMA_mean_sd==0)

cor_col_TL_rma_batch_mean_col <- cor(col_tl_rma_batch_mean[-index_TL_RMA_mean_sd,])
cor_col_TL_rma_batch_mean_row <- cor(t(col_tl_rma_batch_mean[-index_TL_RMA_mean_sd,]))


TL_rma_batch_mean_hc_col <- hclust(as.dist(1-cor_col_TL_rma_batch_mean_col))
TL_rma_batch_mean_hc_row <- hclust(as.dist(1-cor_col_TL_rma_batch_mean_row))

plot(TX_expr_gcrma_hc_col)
plot(TX_expr_gcrma_hc_row)

old.dir <- getwd()
setwd(paste(old.dir,"TL/rma/",sep="/"))
png(filename = "TL_rma_batch_mean_row_all_heatmap.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_rma_batch_mean[-index_TL_RMA_mean_sd,],Rowv = as.dendrogram(TL_rma_batch_mean_hc_row), Colv = as.dendrogram(TL_rma_batch_mean_hc_col),trace="none",scale="row")
dev.off()
setwd(old.dir)


old.dir <- getwd()
setwd(paste(old.dir,"TL/rma/",sep="/"))
png(filename = "TL_rma_batch_mean_row_all_heatmap_no_col_dend.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_rma_batch_mean[-index_TL_RMA_mean_sd,],Rowv = as.dendrogram(TL_rma_batch_mean_hc_row), Colv = 'none',trace="none",scale="row")
dev.off()
setwd(old.dir)


rm(cor_col_TL_rma_batch_mean_col)
rm(cor_col_TL_rma_batch_mean_row)
rm(TL_rma_batch_mean_hc_col)
rm(TL_rma_batch_mean_hc_row)



#GCRMA
index_TL_GCRMA_sd <- apply(col_tl_gcrma,1,sd)
index_TL_GCRMA_sd <- which(index_TL_GCRMA_sd==0)

cor_col_TL_gcrma_col <- cor(col_tl_gcrma[-index_TL_GCRMA_sd,])
cor_col_TL_gcrma_row <- cor(t(col_tl_gcrma[-index_TL_GCRMA_sd,]))

TL_gcrma_hc_col <- hclust(as.dist(1-cor_col_TL_gcrma_col))
TL_gcrma_hc_row <- hclust(as.dist(1-cor_col_TL_gcrma_row))

old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma/",sep="/"))
png(filename = "TL_gcrma_row_all_heatmap.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_gcrma[-index_TL_GCRMA_sd,],Rowv = as.dendrogram(TL_gcrma_hc_row), Colv = as.dendrogram(TL_gcrma_hc_col),trace="none",scale="row")
dev.off()
setwd(old.dir)


old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma/",sep="/"))
png(filename = "TL_gcrma_row_all_heatmap_no_col_dend.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_gcrma[-index_TL_GCRMA_sd,],Rowv = as.dendrogram(TL_gcrma_hc_row), Colv = 'none',trace="none",scale="row")
dev.off()
setwd(old.dir)


rm(cor_col_TL_gcrma_col)
rm(cor_col_TL_gcrma_row)
rm(TL_gcrma_hc_col)
rm(TL_gcrma_hc_row)


#ave heatmap

col_tl_gcrma_mean <- matrix(numeric(),ncol=ncol(col_tl_gcrma)/3,nrow=nrow(col_tl_gcrma))

for(i in 1:9){
  
  col_tl_gcrma_mean[,i] <- apply(col_tl_gcrma[,seq((3*i-2),3*i,1)],1,mean)
  
}

colnames(col_tl_gcrma_mean) <- c("NP","F5","F6","F7","F8","F9","F10","F11","F12")

index_TL_GCRMA_mean_sd <- apply(col_tl_gcrma_mean,1,sd)
index_TL_GCRMA_mean_sd <- which(index_TL_GCRMA_mean_sd==0)

cor_col_TL_gcrma_mean_col <- cor(col_tl_gcrma_mean[-index_TL_GCRMA_mean_sd,])
cor_col_TL_gcrma_mean_row <- cor(t(col_tl_gcrma_mean[-index_TL_GCRMA_mean_sd,]))


TL_gcrma_mean_hc_col <- hclust(as.dist(1-cor_col_TL_gcrma_mean_col))
TL_gcrma_mean_hc_row <- hclust(as.dist(1-cor_col_TL_gcrma_mean_row))


old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma/",sep="/"))
png(filename = "TL_gcrma_mean_row_all_heatmap.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_gcrma_mean[-index_TL_GCRMA_mean_sd,],Rowv = as.dendrogram(TL_gcrma_mean_hc_row), Colv = as.dendrogram(TL_gcrma_mean_hc_col),trace="none",scale="row")
dev.off()
setwd(old.dir)


old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma/",sep="/"))
png(filename = "TL_gcrma_mean_row_all_heatmap_no_col_dend.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_gcrma_mean[-index_TL_GCRMA_mean_sd,],Rowv = as.dendrogram(TL_gcrma_mean_hc_row), Colv = 'none',trace="none",scale="row")
dev.off()
setwd(old.dir)


rm(cor_col_TL_gcrma_mean_col)
rm(cor_col_TL_gcrma_mean_row)
rm(TL_gcrma_mean_hc_col)
rm(TL_gcrma_mean_hc_row)


#gcrma batch removed

index_TL_gcrma_sd <- apply(col.tl.gcrma.rm.batch,1,sd)
index_TL_gcrma_sd <- which(index_TL_gcrma_sd==0)


cor_col_TL_gcrma_batch_col <- cor(col.tl.gcrma.rm.batch[-index_TL_gcrma_sd,])
cor_col_TL_gcrma_batch_row <- cor(t(col.tl.gcrma.rm.batch[-index_TL_gcrma_sd,]))

TL_gcrma_batch_hc_col <- hclust(as.dist(1-cor_col_TL_gcrma_batch_col))
TL_gcrma_batch_hc_row <- hclust(as.dist(1-cor_col_TL_gcrma_batch_row))


old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma/",sep="/"))
png(filename = "TL_gcrma_batch_row_all_heatmap.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col.tl.gcrma.rm.batch[-index_TL_gcrma_sd,],Rowv = as.dendrogram(TL_gcrma_batch_hc_row), Colv = as.dendrogram(TL_gcrma_batch_hc_col),trace="none",scale="row")
dev.off()
setwd(old.dir)


old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma/",sep="/"))
png(filename = "TL_gcrma_batch_row_all_heatmap_no_col_dend.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col.tl.gcrma.rm.batch[-index_TL_gcrma_sd,],Rowv = as.dendrogram(TL_gcrma_batch_hc_row), Colv = 'none',trace="none",scale="row")
dev.off()
setwd(old.dir)


rm(cor_col_TL_gcrma_batch_col)
rm(cor_col_TL_gcrma_batch_row)
rm(TL_gcrma_batch_hc_col)
rm(TL_gcrma_batch_hc_row)


#ave heatmap

col_tl_gcrma_batch_mean <- matrix(numeric(),ncol=ncol(col.tl.gcrma.rm.batch)/3,nrow=nrow(col.tl.gcrma.rm.batch))

for(i in 1:9){
  
  col_tl_gcrma_batch_mean[,i] <- apply(col.tl.gcrma.rm.batch[,seq((3*i-2),3*i,1)],1,mean)
  
}

colnames(col_tl_gcrma_batch_mean) <- c("NP","F5","F6","F7","F8","F9","F10","F11","F12")


index_TL_gcrma_mean_sd <- apply(col_tl_gcrma_batch_mean,1,sd)
index_TL_gcrma_mean_sd <- which(index_TL_gcrma_mean_sd==0)

cor_col_TL_gcrma_batch_mean_col <- cor(col_tl_gcrma_batch_mean[-index_TL_gcrma_mean_sd,])
cor_col_TL_gcrma_batch_mean_row <- cor(t(col_tl_gcrma_batch_mean[-index_TL_gcrma_mean_sd,]))


TL_gcrma_batch_mean_hc_col <- hclust(as.dist(1-cor_col_TL_gcrma_batch_mean_col))
TL_gcrma_batch_mean_hc_row <- hclust(as.dist(1-cor_col_TL_gcrma_batch_mean_row))


old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma/",sep="/"))
png(filename = "TL_gcrma_batch_mean_row_all_heatmap.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_gcrma_batch_mean[-index_TL_gcrma_mean_sd,],Rowv = as.dendrogram(TL_gcrma_batch_mean_hc_row), Colv = as.dendrogram(TL_gcrma_batch_mean_hc_col),trace="none",scale="row")
dev.off()
setwd(old.dir)


old.dir <- getwd()
setwd(paste(old.dir,"TL/gcrma/",sep="/"))
png(filename = "TL_gcrma_batch_mean_row_all_heatmap_no_col_dend.png" ,width = 2500,height = 2000,units = "px",res=300)  
heatmap.2(col_tl_gcrma_batch_mean[-index_TL_gcrma_mean_sd,],Rowv = as.dendrogram(TL_gcrma_batch_mean_hc_row), Colv = 'none',trace="none",scale="row")
dev.off()
setwd(old.dir)


rm(cor_col_TL_gcrma_batch_mean_col)
rm(cor_col_TL_gcrma_batch_mean_row)
rm(TL_gcrma_batch_mean_hc_col)
rm(TL_gcrma_batch_mean_hc_row)

















#GCRMA





































