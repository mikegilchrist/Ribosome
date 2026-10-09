library(data.table)

elong <- fread("../Data/pgen.1007166.s023_elong.txt")

mean_elong_per_transcript <- apply(elong[seq(4,nrow(elong),7),],1,function(i) mean(as.numeric(unlist(i)),na.rm=T))

init <- fread("../Data/pgen.1007166.s023_init.txt",header=T,sep="\t")

mean(unlist(init[,3]))


unlist(init[,3])
length(unlist(init[,3]))
length(mean_elong_per_transcript)

hist(unlist(init[,3])/mean_elong_per_transcript)
hist(log(unlist(init[,3])/mean_elong_per_transcript,10),breaks=seq(-3,-1,0.1))
hist(mean_elong_per_transcript)

png("../Results/20221026/Scaled_initiation_histogram_DUC_pgen_2018.png",width = 6,height = 6,units = "in",res = 300)
hist(log(unlist(init[,3])/mean_elong_per_transcript,10),breaks=seq(-3,-1,0.1),main="",axes=F,xlab="",ylab="")
axis(2,at=seq(0,200,50),labels = seq(0,200,50),cex.axis=1.25)
axis(1,at=seq(-3,-1,0.5),labels = seq(-3,-1,0.5),cex.axis=1.25)
mtext(2,text =  "Frequency",font=1,line=2.5,cex=1.25,xpd=T)
mtext(1,text =  bquote(log[10]*" [Scaled Initiation ("*kappa/tau*")]"),font=2,line=2.5,cex=1.25,xpd=T)
box()
dev.off()



