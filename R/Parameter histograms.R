library(Biostrings)
library(data.table)

sc.proteins <- Biostrings::readAAStringSet("../Data/yeast_proteins.FASTA")

hist(width(sc.proteins),breaks = seq(0,5000,20))
quantile(x = width(sc.proteins),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
quantile(x = width(sc.proteins),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))/9

hs.proteins <- Biostrings::readAAStringSet("../Data/human_proteins.FASTA")

tmp <- hist(width(hs.proteins),breaks = c(seq(0,5000,20),36000))
quantile(x = width(hs.proteins),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))



#yeast half lives
Chan_2018 <- fread("../Data/Sc_Halflife_Chan_2018elife-32536-fig1-data2-v4.txt")
Duttagupta_2005 <- fread("../Data/SC_halflives_Duttagupta_2005_Table_S1_abc.txt")
Presnyak_2015 <- fread("../Data/SC_halflives_Presnyak_2015.txt")
Shalem_2008 <- fread("../Data/SC_halflives_Shalem_2008_abc.txt")


quantile(Chan_2018[,3],c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
quantile(Presnyak_2015[,4],c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
quantile(as.numeric(unlist(Shalem_2008[-2,24])),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)

mean(unlist(Chan_2018[,3]),na.rm=T)
mean(unlist(Presnyak_2015[,4]),na.rm=T)
mean(unlist(as.numeric(unlist(Shalem_2008[-2,24]))),na.rm=T)


median(unlist(Chan_2018[,3]),na.rm=T)
median(unlist(Presnyak_2015[,4]),na.rm=T)
median(unlist(as.numeric(unlist(Shalem_2008[-2,24]))),na.rm=T)


marking_rate_Presnyak <- log(2)/(round(quantile(Presnyak_2015[,4],c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T),0)*60)

marking_rate_Chan <- log(2)/(quantile(Chan_2018[,3],c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)*60)






hist(as.numeric(unlist(Presnyak_2015[,4])),breaks=seq(0,60,1),ylim=c(0,0.25),probability = T,col=rgb(0.7,0.7,0.7,0.5))
hist(as.numeric(unlist(Chan_2018[,3])),add=T,breaks=seq(0,60,1),probability = T,col=rgb(0.5,0.5,0.5,0.5))

Presnyak_2015[,c(1,4)]
Chan_2018[,c(1,3)]
setkeyv(Presnyak_2015,"Gene ID")
setkeyv(Chan_2018,"gene_id")
index <- unlist(Presnyak_2015[,c(1)][which(unlist(Presnyak_2015[,c(1)])%in%unlist(Chan_2018[,c(1)])),])

Presnyak_2015[index,c(1,4)]
Chan_2018[index,c(1,3)]
all(Presnyak_2015[index,c(1)]==Chan_2018[index,c(1)])
which(is.na(Chan_2018[index,c(3)]))


cor(as.numeric(unlist(Presnyak_2015[index,c(4)][-which(is.na(Chan_2018[index,c(3)]))])),as.numeric(unlist(Chan_2018[index,c(3)][-which(is.na(Chan_2018[index,c(3)]))])))

plot(as.numeric(unlist(Presnyak_2015[index,c(4)][-which(is.na(Chan_2018[index,c(3)]))]))/max(as.numeric(unlist(Presnyak_2015[index,c(4)][-which(is.na(Chan_2018[index,c(3)]))]))),
     as.numeric(unlist(Chan_2018[index,c(3)][-which(is.na(Chan_2018[index,c(3)]))]))/max(as.numeric(unlist(Chan_2018[index,c(3)][-which(is.na(Chan_2018[index,c(3)]))]))))



old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
# Protein length distribution ---------------------------------------------
png("Figure1_Sc_Protein_length_histogram.png",width = 5,height = 5,units = "in",res = 300)
sc.proteins.length <- width(sc.proteins)
sc.proteins.trim <- sc.proteins.length
sc.proteins.trim[which(sc.proteins.length>1500)] <- 1500
prot.length.hist <- hist(sc.proteins.trim,breaks = seq(0,1500,20),ylim=c(0,1100),col="white",axes=F,xlab="",main="Yeast Protein Lengths")
abline(v=quantile(sc.proteins.trim,c(0.01,0.5,0.99),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,10^max(nchar(prot.length.hist$counts)),250),labels = seq(0,10^max(nchar(prot.length.hist$counts)),250))
axis(1,at=seq(0,1500,500),labels = c("0","500","1000",">1500"))
mtext(1,line=2.25,text = "Protein Length (a.a.)")
legend("topright",legend = c("Small","1% (38aa)","Medium","50% (354aa)","Large","99% (1750aa)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
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

HL_At[-which(HL_At[,2]==0),"HL":=unlist(log(2)/HL_At[-which(HL_At[,2]==0),2])]

median(sort(HL_minutes)[-c(17900:18585)])
HL_minutes <- sort(HL_minutes)[-c(17900:18585)]

quantile(HL_minutes,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99),na.rm=T)
HL_minutes[(which(HL_minutes>1500))]<-1500

ylim=c(0,10)
old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Figure1_Sc_Transcript_half_life_histogram.png",width = 5,height = 5,units = "in",res = 300)
half.life.hist <- hist(unlist(Presnyak_2015[,4]),col="white",breaks = seq(0,60,1),axes=F,xlab="",main="Yeast Half-life",ylim=c(0,500))
abline(v=quantile(unlist(Presnyak_2015[,4]),c(0.015,0.475,0.989),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,500,100),labels = seq(0,500,100),las=2)
axis(1,at=seq(0,60,10),labels =seq(0,60,10))
mtext(1,line=2.25,text = "Transcript Half-life (min.)")
legend("topright",legend = c("Short","~1% (2 min.)","Medium","~50% (7 min.)","Long","~99% (52 min.)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
box()
dev.off()
setwd(old.dir)



#comp with duc

sc.proteins
#Read in intiation rates from Duc 2018
Duc_Initiation <- fread("../../Data/Duc2018/Initiation_rates_DUC.txt", sep = "\t", stringsAsFactors = F,header=T)

#Read in elongation rates from Duc 2018
Duc_Elongation <- fread("../../Data/Duc2018/Elongation_rates_DUC.txt", sep = "\t", stringsAsFactors = F,header=F)

View(Duc_Elongation)

Duc_Genes <- Duc_Elongation[seq(1,nrow(Duc_Elongation),7),1]
Duc_Sequence <- Duc_Elongation[seq(2,nrow(Duc_Elongation),7),]
Duc_Experimental <- Duc_Elongation[seq(3,nrow(Duc_Elongation),7),]
Duc_Infered <- Duc_Elongation[seq(4,nrow(Duc_Elongation),7),]
Duc_Sequence_lengths <- unlist(apply(Duc_Sequence,1,function(i) grep(x = unlist(i),pattern = "TGA|TAG|TAA")))


Arith_mean_inferred <- apply(Duc_Infered,1,function(i) mean(as.numeric(unlist(i)),na.rm=T))
names(Arith_mean_inferred) <- unlist(Duc_Genes)



Arith_mean_inferred
Duc.DT <- cbind(Duc_Initiation[,c(1,3)],Arith_mean_inferred)
setkeyv(Duc.DT,"ORF")

plot(unlist(Duc.DT[,2]),unlist(Duc.DT[,3]))


Presnyak_2015_Duc <- Presnyak_2015[which(unlist(Presnyak_2015[,c(1)])%in%unlist(Duc.DT[,1])),]

Duc.DT[unlist(Presnyak_2015_Duc[,1]),]

plot(unlist(Presnyak_2015_Duc[,4]),unlist(Duc.DT[unlist(Presnyak_2015_Duc[,1]),2]))
plot(unlist(Presnyak_2015_Duc[,4]),unlist(Duc.DT[unlist(Presnyak_2015_Duc[,1]),3]))
plot(unlist(Presnyak_2015_Duc[,4]),log(unlist(Duc.DT[unlist(Presnyak_2015_Duc[,1]),Inferred_Initiation_rate/Arith_mean_inferred]),10))

cor.test(unlist(Presnyak_2015_Duc[,4]),unlist(Duc.DT[unlist(Presnyak_2015_Duc[,1]),2]))
cor.test(unlist(Presnyak_2015_Duc[,4]),unlist(Duc.DT[unlist(Presnyak_2015_Duc[,1]),3]))

sc.names <- unlist(lapply(names(sc.proteins),function(i) strsplit(x=i,split="\\|")[[1]][1]))

sc.proteins_Duc <- width(sc.proteins)[which(sc.names%in%unlist(Duc.DT[,1]))]
names(sc.proteins_Duc) <- sc.names[which(sc.names%in%unlist(Duc.DT[,1]))]

plot(unlist(sc.proteins_Duc),unlist(Duc.DT[unlist(names(sc.proteins_Duc)),Inferred_Initiation_rate]))
plot(unlist(sc.proteins_Duc),unlist(Duc.DT[unlist(names(sc.proteins_Duc)),Arith_mean_inferred]))
plot(unlist(sc.proteins_Duc),log(unlist(Duc.DT[unlist(names(sc.proteins_Duc)),Inferred_Initiation_rate/Arith_mean_inferred]),10))

cor.test(unlist(sc.proteins_Duc),unlist(Duc.DT[unlist(names(sc.proteins_Duc)),Inferred_Initiation_rate]))
cor.test(unlist(sc.proteins_Duc),unlist(Duc.DT[unlist(names(sc.proteins_Duc)),Arith_mean_inferred]))
cor.test(unlist(sc.proteins_Duc),log(unlist(Duc.DT[unlist(names(sc.proteins_Duc)),Inferred_Initiation_rate/Arith_mean_inferred]),10))


Weinberg_mRNA_RPKM <- fread("../Data/GSE75897_RiboZero_RPKMs.txt",sep="\t",stringsAsFactors = F,header = F)
setkey(Weinberg_mRNA_RPKM,"V1")

Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),]
unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2])

plot(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2]),10),unlist(Duc.DT[,Inferred_Initiation_rate]))
plot(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2]),10),unlist(Duc.DT[,Arith_mean_inferred]))
plot(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2]),10),log(unlist(Duc.DT[,Inferred_Initiation_rate/Arith_mean_inferred]),10))
plot(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[names(unlist(sc.proteins_Duc)),1]),2]),10),unlist(sc.proteins_Duc))
plot(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[unlist(Presnyak_2015_Duc[,1]),1]),2]),10),unlist(Presnyak_2015_Duc[,4]))


cor.test(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2]),10),unlist(Duc.DT[,Inferred_Initiation_rate]))
cor.test(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2]),10),unlist(Duc.DT[,Arith_mean_inferred]))
cor.test(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2]),10),log(unlist(Duc.DT[,Inferred_Initiation_rate/Arith_mean_inferred]),10))
cor.test(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[names(unlist(sc.proteins_Duc)),1]),2]),10),unlist(sc.proteins_Duc))
cor.test(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[unlist(Presnyak_2015_Duc[,1]),1]),2]),10),unlist(Presnyak_2015_Duc[,4]))


Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),]

hist(unlist(Weinberg_mRNA_RPKM[,2]))


old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Figure1_Sc_Transcript_RPKM_histogram.png",width = 5,height = 5,units = "in",res = 300)
half.life.hist <- hist(log(unlist(Weinberg_mRNA_RPKM[,2])+1,2),col="white",breaks = seq(0,15,0.5),axes=F,xlab="",main="Yeast RPKM",ylim=c(0,750),xlim=c(0,20))
abline(v=quantile(log(unlist(Weinberg_mRNA_RPKM[,2])+1,2),c(0.01,0.5,0.99),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,750,75),labels = seq(0,750,75),las=2)
axis(1,at=seq(-5,15,5),labels =seq(-5,15,5))
mtext(1,line=2.25,text = bquote(log[2]*" [Transcript abundance (RPKM)]"))
legend("topright",legend = c("Low","~1% (1.29 RPKM)","Median","~50% (47 RPKM)","High","~99% (3,283 RPKM)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
box()
dev.off()
setwd(old.dir)


old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Figure1_Sc_Transcript_log2FC_histogram.png",width = 5,height = 5,units = "in",res = 300)
half.life.hist <- hist(log(unlist(Weinberg_mRNA_RPKM[,2])+1,2)-log(47.29,2),col="white",breaks = seq(-6,9,0.5),axes=F,xlab="",main="Yeast RPKM",ylim=c(0,750),xlim=c(-6,15))
abline(v=quantile(log(unlist(Weinberg_mRNA_RPKM[,2])+1,2)-log(47.29,2),c(0.01,0.5,0.99),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,750,75),labels = seq(0,750,75),las=2)
axis(1,at=seq(-6,9,2),labels =seq(-6,9,2))
mtext(1,line=2.25,text = bquote(log[2]*" Fold Change[Transcript abundance (RPKM)]"))
legend("topright",legend = c("Low","~1% (0.29 RPKM)","Median","~50% (46 RPKM)","High","~99% (3,283 RPKM)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5)
box()
dev.off()
setwd(old.dir)

old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Figure1_Sc_Transcript_log10FC_histogram.png",width = 5,height = 5,units = "in",res = 300)
half.life.hist <- hist(log(unlist(Weinberg_mRNA_RPKM[,2])+1,10)-log(47.29,10),col="white",breaks = seq(-2,3,0.1),axes=F,xlab="",main="Yeast RPKM",ylim=c(0,600),xlim=c(-2,))
abline(v=quantile(log(unlist(Weinberg_mRNA_RPKM[,2])+1,10)-log(47.29,10),c(0.01,0.5,0.99),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,600,100),labels = seq(0,600,100),las=2)
axis(1,at=seq(-2,3,1),labels =seq(-2,3,1))
mtext(1,line=2.25,text = bquote(log[10]*" Fold Change[Transcript abundance (RPKM)]"))
legend("topright",legend = c("Low","~1% (0.29 RPKM)","Median","~50% (46 RPKM)","High","~99% (3,283 RPKM)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5)
box()
dev.off()
setwd(old.dir)

#just duc genes

old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Figure1_Sc_Transcript_RPKM_histogram.png",width = 5,height = 5,units = "in",res = 300)
half.life.hist <- hist(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2])+1,2),col="white",breaks = seq(0,15,0.5),axes=F,xlab="",main="Yeast RPKM",ylim=c(0,750),xlim=c(0,20))
abline(v=quantile(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2])+1,2),c(0.01,0.5,0.99),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,750,75),labels = seq(0,750,75),las=2)
axis(1,at=seq(-5,15,5),labels =seq(-5,15,5))
mtext(1,line=2.25,text = bquote(log[2]*" [Transcript abundance (RPKM)]"))
legend("topright",legend = c("Low","~1% (1.29 RPKM)","Median","~50% (47 RPKM)","High","~99% (3,283 RPKM)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5,)
box()
dev.off()
setwd(old.dir)



old.dir <- getwd()
Check_dir("../Results/",add=T)
setwd(paste(old.dir,"../Results/",gsub(Sys.Date(),pattern = "-",replacement = ""),sep="/"))
png("Figure1_Sc_Transcript_log10FC_histogram.png",width = 5,height = 5,units = "in",res = 300)
half.life.hist <- hist(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2])+1,10)-log(47.29,10),col="white",breaks = seq(-2,3,0.1),axes=F,xlab="",main="Yeast RPKM",ylim=c(0,600),xlim=c(-2,))
abline(v=quantile(log(unlist(Weinberg_mRNA_RPKM[unlist(Duc.DT[,1]),2])+1,10)-log(47.29,10),c(0.01,0.5,0.99),na.rm=T),col="black",lty=c(3,2,4),lwd=2)
axis(2,at=seq(0,600,100),labels = seq(0,600,100),las=2)
axis(1,at=seq(-2,3,1),labels =seq(-2,3,1))
mtext(1,line=2.25,text = bquote(log[10]*" Fold Change[Transcript abundance (RPKM)]"))
legend("topright",legend = c("Low","~1% (0.29 RPKM)","Median","~50% (46 RPKM)","High","~99% (3,283 RPKM)"),lty=c(3,3,2,2,4,4),col=rep(c("black","white"),3),lwd=2,xjust = 0.5)
box()
dev.off()
setwd(old.dir)
