source("R/Polysome_functions_RAUC.R")

#Create parameter matrix with kappa, mu, lambda and delta ranging from 10^-5 to 10^5
#tau is set to 1 and sd = 1 (inherited from fitting functions)
par.names <- c("kappa","mu","lambda","tau","delta","s2")
par.range <- list(c(10^-5,10^-4,10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3,10^4,10^5))
par.mat <- CreateParMatrix(parameter.range = par.range, fixed = c(4,6),fixed.val = 1,par.names = par.names )
#store as data.table for quick access
par.mat.dt <- as.data.table(par.mat)
#index each parameter entry
par.mat.dt[,"index":=1:nrow(par.mat.dt)]


#Lambda is simply a scaling factor for the whole system. Therefore I will limit the search
#space to lambda = 1
par.mat.lambda1 <- par.mat.dt[lambda==1,]
par.mat.lambda1[,"index2":=(1:nrow(par.mat.lambda1))]
setkeyv(par.mat.lambda1,"index2")

# Nondimensionalized results ----------------------------------------------
# All parameters are scaled by tau

Unmarked_imax30 <- apply(par.mat.lambda1,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax30 <- apply(par.mat.lambda1,1,function(i) CalcMarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

U_sum <- colSums(Unmarked_imax30)
U_percent <- apply(Unmarked_imax30,1,function(i) (i/U_sum))
M_percent <- Marked_imax30/(Marked_imax30+Unmarked_imax30)*100
U_Mean_Ribosomes <- unlist(U_percent%*%c(0:30))

plot(sort(unique(U_percent%*%c(0:30))))

which(unlist(U_percent%*%c(0:30))%in%unlist(unique(U_percent%*%c(0:30))))

x <- lapply(unique(U_Mean_Ribosomes),function(i) grep(pattern = i,U_Mean_Ribosomes))

par.mat.dt[x[[2]],]
par.mat.dt[x[[2]],kappa/mu]
U_Mean_Ribosomes[x[[2]]]

y<-lapply(1:length(x),function(i) unique(par.mat.dt[x[[i]],kappa/mu]))

par.mat.lambda1.delta100000.tau1 <- par.mat.lambda1[which(par.mat.lambda1[,delta]==100000),]
par.mat.lambda1.delta100000.tau1[,"k.m":=kappa/mu]
par.mat.l1.d1e5.t1.unique <- par.mat.lambda1.delta100000.tau1[-which(duplicated(round(par.mat.lambda1.delta100000.tau1[,k.m],6))),]
par.mat.l1.d1e5.t1.unique <- par.mat.l1.d1e5.t1.unique[order(unique(round(par.mat.l1.d1e5.t1.unique[,k.m],6))),]
format(sort(unique(round(par.mat.lambda1.delta100000.tau1[,k.m],6))),scientific=FALSE)

Unmarked_imax30_NPM <- apply(par.mat.l1.d1e5.t1.unique,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax30_NPM <- apply(par.mat.l1.d1e5.t1.unique,1,function(i) CalcMarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))


Unmarked_imax30_NPM
plot(0:30,log(Unmarked_imax30_NPM[,1],10),type="l",ylim=c(0,log(max(Unmarked_imax30_NPM),10)))
for(i in 2:ncol(Unmarked_imax30_NPM)){
  lines(0:30, log(Unmarked_imax30_NPM[,i],10),col=i)
}

U_sum_NPM <- colSums(Unmarked_imax30_NPM)
U_percent_NPM <- t(apply(Unmarked_imax30_NPM,1,function(i) (i/U_sum_NPM)))
M_percent_NPM <- Marked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM)%*%c(0:30))

plot(0:30,U_percent_NPM[,1],type="l",ylim=c(0,1))
for(i in 2:ncol(U_percent_NPM)){
  lines(0:30,U_percent_NPM[,i],col=i)
}


#filled.contour(x=0:30,y=log(par.mat.l1.d1e5.t1.unique[,k.m],10),z=U_percent_NPM,)



par.mat.lambda1.delta100000.tau1.mu1 <- par.mat.lambda1[which(par.mat.lambda1[,delta]==100000&par.mat.lambda1[,mu]==1),]
par.mat.lambda1.delta100000.tau1.mu1[,"k.m":=kappa/mu]

Unmarked_imax30_NPM <- apply(par.mat.lambda1.delta100000.tau1.mu1,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax30_NPM <- apply(par.mat.lambda1.delta100000.tau1.mu1,1,function(i) CalcMarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))


Unmarked_imax30_NPM
plot(0:30,Unmarked_imax30_NPM[,1],type="l",ylim=c(0,1))
for(i in 2:ncol(Unmarked_imax30_NPM)){
  lines(0:30, Unmarked_imax30_NPM[,i],col=i)
}

U_sum_NPM <- colSums(Unmarked_imax30_NPM)
U_percent_NPM <- t(apply(Unmarked_imax30_NPM,1,function(i) (i/U_sum_NPM)))
M_percent_NPM <- Marked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM)%*%c(0:30))

plot(0:30,U_percent_NPM[,1],type="l",ylim=c(0,1))
for(i in 2:ncol(U_percent_NPM)){
  lines(0:30,U_percent_NPM[,i],col=i)
}

filled.contour(x=0:30,y=-5:5,z=U_percent_NPM)
# k/m series --------------------------------------------------------------

krange <- list(c(seq(0.00001,0.00009,0.00001),seq(0.0001,0.0009,0.0001),seq(0.001,0.009,0.001),
                 seq(0.01,0.99,0.01),seq(1,9.9,0.1),seq(10,100000,10)))
par.mat.k.only <- CreateParMatrix(parameter.range = krange, fixed = c(2:6),
                           fixed.val = 1,par.names = par.names )
#store as data.table for quick access
par.mat.k.only.dt <- as.data.table(par.mat.k.only)
#index each parameter entry
par.mat.k.only.dt[,"index":=1:nrow(par.mat.k.only.dt)]


Unmarked_imax30_NPM <- apply(par.mat.k.only.dt,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax30_NPM <- apply(par.mat.k.only.dt,1,function(i) CalcMarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

#Unmarked_imax30_NPM
plot(0:30,log(Unmarked_imax30_NPM[,1],10),type="l",ylim=c(0,1))
for(i in 2:ncol(Unmarked_imax30_NPM)){
  lines(0:30, Unmarked_imax30_NPM[,i],col=i)
}

U_sum_NPM <- colSums(Unmarked_imax30_NPM)
U_percent_NPM <- t(apply(Unmarked_imax30_NPM,1,function(i) (i/U_sum_NPM)))
M_percent_NPM <- Marked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM)%*%c(0:30))
difference <- apply(matrix(rep(c(0:30),ncol(Unmarked_imax30_NPM)),nrow=nrow(Unmarked_imax30_NPM)),1,function(i) i-U_Mean_Ribosomes_NPM)
U_SD_Ribosomes_NPM <- diag(sqrt(unlist(t(U_percent_NPM)%*%t((difference)^2))))


plot(0:30,U_percent_NPM[,1],type="l",ylim=c(0,1))
for(i in 2:ncol(U_percent_NPM)){
  lines(0:30,U_percent_NPM[,i],col=i)
}


approx_s2 <- unlist(lapply(1:length(U_Mean_Ribosomes_NPM),function(i) s2_approx(parameters = unlist(par.mat.k.only.dt[i,]),imax = 30,mean = U_Mean_Ribosomes_NPM[i])))


png("../Results/20220210/contourl1m1d1t1_k1e-5_1e5.png",width = 6,height = 8,units = "in",res=600)
filled.contour(x=0:30,y=log(par.mat.k.only.dt[,kappa],10),z=U_percent_NPM,color.palette = rainbow)
dev.off()
contour(x=0:30,y=log(par.mat.k.only.dt[,kappa],10),z=U_percent_NPM)
.filled.contour(x=0:30,y=log(par.mat.k.only.dt[,kappa],10),z=U_percent_NPM,col = rainbow(21),levels = seq(0,1,0.05))
contour(x=0:30,y=log(par.mat.k.only.dt[,kappa],10),z=U_percent_NPM,add=T,levels = seq(0,1,0.05))
#lines((U_Mean_Ribosomes_NPM)/1.65-7,log(par.mat.k.only.dt[,kappa],10)) approximate do not use. find actual way to add
lines(U_Mean_Ribosomes_NPM,log(par.mat.k.only.dt[,kappa],10),col="green")

plot(log(par.mat.k.only.dt[,kappa],10),U_Mean_Ribosomes_NPM,type="l",ylim=c(0,30))
plot(U_Mean_Ribosomes_NPM,U_SD_Ribosomes_NPM,type="l")
plot(log(par.mat.k.only.dt[,kappa],10),U_SD_Ribosomes_NPM,type="l")
abline(v=mean(U_SD_Ribosomes_NPM))
mean(U_SD_Ribosomes_NPM)
sd(U_SD_Ribosomes_NPM)
plot(dnorm(log(par.mat.k.only.dt[,kappa],10),mean = mean(U_SD_Ribosomes_NPM),sd = sd(U_SD_Ribosomes_NPM)))
lines(log(par.mat.k.only.dt[,kappa],10),approx_s2,col="red")

# k/m series no interference --------------------------------------------------------------

krange <- list(c(seq(0.00001,0.00009,0.00001),seq(0.0001,0.0009,0.0001),seq(0.001,0.009,0.001),
                 seq(0.01,0.99,0.01),seq(1,9.9,0.1),seq(10,100000,10)))
par.mat.k.only <- CreateParMatrix(parameter.range = krange, fixed = c(2:6),
                                  fixed.val = 1,par.names = par.names )
#store as data.table for quick access
par.mat.k.only.dt <- as.data.table(par.mat.k.only)
#index each parameter entry
par.mat.k.only.dt[,"index":=1:nrow(par.mat.k.only.dt)]
par.mat.k.only.dt[,"mu":=1/36000]

Unmarked_imax30_NPM_NI <- apply(par.mat.k.only.dt,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("tau")))
Marked_imax30_NPM_NI <- apply(par.mat.k.only.dt,1,function(i) CalcMarkedClass(parameters = i,imax = 30,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("tau")))

#Unmarked_imax30_NPM_NI
plot(0:30,log(Unmarked_imax30_NPM_NI[,1],10),type="l",ylim=c(0,1))
for(i in 2:ncol(Unmarked_imax30_NPM_NI)){
  lines(0:30, Unmarked_imax30_NPM_NI[,i],col=i)
}

U_sum_NPM_NI <- colSums(Unmarked_imax30_NPM_NI)
U_percent_NPM_NI <- t(apply(Unmarked_imax30_NPM_NI,1,function(i) (i/U_sum_NPM_NI)))
M_percent_NPM_NI <- Marked_imax30_NPM_NI/(Marked_imax30_NPM_NI+Unmarked_imax30_NPM_NI)*100
U_Mean_Ribosomes_NPM_NI <- unlist(t(U_percent_NPM_NI)%*%c(0:30))
difference <- apply(matrix(rep(c(0:30),ncol(Unmarked_imax30_NPM_NI)),nrow=nrow(Unmarked_imax30_NPM_NI)),1,function(i) i-U_Mean_Ribosomes_NPM_NI)
U_SD_Ribosomes_NPM_NI <- diag(sqrt(unlist(t(U_percent_NPM_NI)%*%t((difference)^2))))


plot(0:30,U_percent_NPM_NI[,1],type="l",ylim=c(0,1))
for(i in 2:ncol(U_percent_NPM_NI)){
  lines(0:30,U_percent_NPM_NI[,i],col=i)
}


approx_s2 <- unlist(lapply(1:length(U_Mean_Ribosomes_NPM_NI),function(i) s2_approx(parameters = unlist(par.mat.k.only.dt[i,]),imax = 30,mean = U_Mean_Ribosomes_NPM_NI[i])))


png("../Results/20220210/contourl1m1d1t1_k1e-5_1e5_noint.png",width = 6,height = 8,units = "in",res=600)
filled.contour(x=0:30,y=log(par.mat.k.only.dt[,kappa],10),z=U_percent_NPM_NI,color.palette = rainbow)
dev.off()
contour(x=0:30,y=log(par.mat.k.only.dt[,kappa],10),z=U_percent_NPM_NI)
.filled.contour(x=0:30,y=log(par.mat.k.only.dt[,kappa],10),z=U_percent_NPM_NI,col = rainbow(21),levels = seq(0,1,0.05))
contour(x=0:30,y=log(par.mat.k.only.dt[,kappa],10),z=U_percent_NPM_NI,add=T,levels = seq(0,1,0.05))
#lines((U_Mean_Ribosomes_NPM_NI)/1.65-7,log(par.mat.k.only.dt[,kappa],10)) approximate do not use. find actual way to add
lines(U_Mean_Ribosomes_NPM_NI,log(par.mat.k.only.dt[,kappa],10))
lines(U_Mean_Ribosomes_NPM,log(par.mat.k.only.dt[,kappa],10),col="green")

plot(log(par.mat.k.only.dt[,kappa],10),U_Mean_Ribosomes_NPM_NI,type="l",ylim=c(0,30))
lines(log(par.mat.k.only.dt[,kappa],10),U_Mean_Ribosomes_NPM,col='red')

plot(U_Mean_Ribosomes_NPM_NI,log(par.mat.k.only.dt[,kappa],10),type="l",ylim=c(-5,5))
lines(U_Mean_Ribosomes_NPM,log(par.mat.k.only.dt[,kappa],10),col='red')

plot(log(par.mat.k.only.dt[,kappa],10),U_Mean_Ribosomes_NPM_NI-U_Mean_Ribosomes_NPM,type="l",ylim=c(0,30))
plot(log(par.mat.k.only.dt[,kappa],10),U_Mean_Ribosomes_NPM_NI/U_Mean_Ribosomes_NPM,type="l",ylim=c(0,2))

plot(U_Mean_Ribosomes_NPM_NI,U_SD_Ribosomes_NPM_NI,type="l")
lines(U_Mean_Ribosomes_NPM_NI,U_SD_Ribosomes_NPM,col="red")

plot(log(par.mat.k.only.dt[,kappa],10),U_SD_Ribosomes_NPM_NI,type="l")
lines(log(par.mat.k.only.dt[,kappa],10),U_SD_Ribosomes_NPM,col="red")


# k/m series --------------------------------------------------------------

krange <- list(c(seq(0.01,0.9,0.05),seq(1,20,1),seq(20,1000,40),seq(1000,10000,5000)))
par.mat.k.tau.mu <- CreateParMatrix(parameter.range = krange, fixed = c(3,5,6),
                                  fixed.val = 1,par.names = par.names )
#store as data.table for quick access
par.mat.k.tau.mu.dt <- as.data.table(par.mat.k.tau.mu)
#index each parameter entry
par.mat.k.tau.mu.dt[,"index":=1:nrow(par.mat.k.tau.mu.dt)]

#tau range
pars <- par.mat.k.tau.mu.dt[which(par.mat.k.tau.mu.dt[,mu]==1&par.mat.k.tau.mu.dt[,kappa]==100),]
pars <- pars[which(pars[,tau]<=20),]
Unmarked_imax30_NPM <- apply(pars,1,function(i) CalcUnmarkedClassNPM(parameters = i,imax = 30))
Marked_imax30_NPM <- apply(pars,1,function(i) CalcMarkedClassNPM(parameters = i,imax = 30))
plot(0:30,Unmarked_imax30_NPM[,1],type="l",ylim=c(0,max(Unmarked_imax30_NPM)))
for(i in 2:ncol(Unmarked_imax30_NPM)){
  lines(0:30,Unmarked_imax30_NPM[,i],col=i)
}

count <- 1
png("../../../Results/20220105/contourl1m1d1t0.01_t20_k1e-5_1e5_mean.png",width = 6,height = 8,units = "in",res=600)
for(i in c(0.01,0.1,1,10,100,1000,10000)){
  pars[,kappa:=..i]
  Unmarked_imax30_NPM <- apply(pars,1,function(i) CalcUnmarkedClassNPM(parameters = i,imax = 30))
  Marked_imax30_NPM <- apply(pars,1,function(i) CalcMarkedClassNPM(parameters = i,imax = 30))
  if(count==1){
    plot(pars[,tau],t(Unmarked_imax30_NPM)%*%0:30,ylim=c(0,30),xlim=c(0,40),pch=20,
         xlab="tau",ylab="mean i")
    lines(pars[,tau],t(Unmarked_imax30_NPM)%*%0:30)
  }else{
    points(pars[,tau],t(Unmarked_imax30_NPM)%*%0:30,pch=20,col=count)
    lines(pars[,tau],t(Unmarked_imax30_NPM)%*%0:30,col=count)
  }
  count <- count + 1

}
legend("right",legend = paste("k/m = ",c(0.01,0.1,1,10,100,1000,10000),sep=""),fill = 1:count)
dev.off()


count <- 1
png("../../../Results/20220105/contourl1m1d1t0.01_t20_k1e-5_1e5_sd.png",width = 6,height = 8,units = "in",res=600)
for(i in c(0.01,0.1,1,10,100,1000,10000)){
  pars[,kappa:=..i]
  Unmarked_imax30_NPM <- apply(pars,1,function(i) CalcUnmarkedClassNPM(parameters = i,imax = 30))
  Marked_imax30_NPM <- apply(pars,1,function(i) CalcMarkedClassNPM(parameters = i,imax = 30))
  Unmarked_imax30_NPM_mean <- t(Unmarked_imax30_NPM)%*%0:30
  var<-((t(matrix(0:30,ncol = nrow(Unmarked_imax30_NPM_mean),nrow=31))-matrix(Unmarked_imax30_NPM_mean,nrow=nrow(Unmarked_imax30_NPM_mean),ncol=31))^2)%*%0:30
  print(apply(Unmarked_imax30_NPM,2,sd))
  print(sqrt(var))
  if(count==1){
    plot(pars[,tau],log(sqrt(var),10),ylim=c(2,2.8),xlim=c(0,40),pch=20,
         xlab = "tau",ylab="sd i")
    lines(pars[,tau],log(sqrt(var),10))
  }else{
    points(pars[,tau],log(sqrt(var),10),pch=20,col=count)
    lines(pars[,tau],log(sqrt(var),10),col=count)
  }
  # if(count==1){
  #   plot(pars[,tau],apply(Unmarked_imax30_NPM,2,sd),ylim=c(0,0.2),xlim=c(0,40),pch=20)
  #   lines(pars[,tau],apply(Unmarked_imax30_NPM,2,sd))
  # }else{
  #   points(pars[,tau],apply(Unmarked_imax30_NPM,2,sd),pch=20,col=count)
  #   lines(pars[,tau],apply(Unmarked_imax30_NPM,2,sd),col=count)
  # }
  #
  count <- count + 1

}
legend("right",legend = paste("k/m = ",c(0.01,0.1,1,10,100,1000,10000),sep=""),fill = 1:count)
dev.off()



U_sum_NPM <- colSums(Unmarked_imax30_NPM)
U_percent_NPM <- t(apply(Unmarked_imax30_NPM,1,function(i) (i/U_sum_NPM)))
M_percent_NPM <- Marked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM)%*%c(0:30))

plot(0:30,U_percent_NPM[,1],type="l",ylim=c(0,1))
for(i in 2:ncol(U_percent_NPM)){
  lines(0:30,U_percent_NPM[,i],col=i)
}


png("../../../Results/20220105/contourl1m1d1t1_k1e-5_1e5.png",width = 6,height = 8,units = "in",res=600)
filled.contour(x=0:30,y=log(par.mat.k.tau.mu.dt[,kappa],10),z=U_percent_NPM,color.palette = rainbow)
dev.off()



# imax range based on k range ---------------------------------------------

krange <- list(c(seq(0.00001,0.00009,0.00001),seq(0.0001,0.0009,0.0001),seq(0.001,0.009,0.001),
                 seq(0.01,0.99,0.01),seq(1,9.9,0.1),seq(10,100000,10)))
par.mat.k.only <- CreateParMatrix(parameter.range = krange, fixed = c(2:6),
                                  fixed.val = 1,par.names = par.names )
#store as data.table for quick access
par.mat.k.only.dt <- as.data.table(par.mat.k.only)
#index each parameter entry
par.mat.k.only.dt[,"index":=1:nrow(par.mat.k.only.dt)]


Imax_Range_tau1_mean_Ribosomes_U <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau1_sd_Ribosomes_U <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau1_mean_Ribosomes_M <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau1_sd_Ribosomes_M <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])

#tau=1
count <- 1
for(k in c(2:50,100,150,200,250,300)){
print(k)
Unmarked_imax30_NPM <- apply(par.mat.k.only.dt,1,function(i) CalcUnmarkedClass(parameters = i,imax = k,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax30_NPM <- apply(par.mat.k.only.dt,1,function(i) CalcMarkedClass(parameters = i,imax = k,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

U_sum_NPM <- colSums(Unmarked_imax30_NPM)
M_sum_NPM <- colSums(Marked_imax30_NPM)
U_percent_NPM_unmarked <- t(apply(Unmarked_imax30_NPM,1,function(i) (i/U_sum_NPM)))*100
M_percent_NPM_marked <- t(apply(Marked_imax30_NPM,1,function(i) (i/M_sum_NPM)))*100
U_percent_NPM_all <- Unmarked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
M_percent_NPM_all <- Marked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM_unmarked)%*%c(0:k))
M_Mean_Ribosomes_NPM <- unlist(t(M_percent_NPM_marked)%*%c(0:k))
difference <- apply(matrix(rep(c(0:k),ncol(Unmarked_imax30_NPM)),nrow=nrow(Unmarked_imax30_NPM)),1,function(i) i-U_Mean_Ribosomes_NPM/100)
U_SD_Ribosomes_NPM <- diag(sqrt(unlist(t(U_percent_NPM_unmarked)%*%t((difference)^2))))
difference <- apply(matrix(rep(c(0:k),ncol(Marked_imax30_NPM)),nrow=nrow(Marked_imax30_NPM)),1,function(i) i-M_Mean_Ribosomes_NPM/100)
M_SD_Ribosomes_NPM <- diag(sqrt(unlist(t(M_percent_NPM_marked)%*%t((difference)^2))))

Imax_Range_tau1_mean_Ribosomes_U[count,] <- U_Mean_Ribosomes_NPM
Imax_Range_tau1_sd_Ribosomes_U[count,] <- U_SD_Ribosomes_NPM
Imax_Range_tau1_mean_Ribosomes_M[count,] <- M_Mean_Ribosomes_NPM
Imax_Range_tau1_sd_Ribosomes_M[count,] <- M_SD_Ribosomes_NPM

count <- count + 1
#U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM)%*%c(0:k))


# png(paste("../../../Results/20220117/contour_norm_unmarked_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_unmarked,2),color.palette = rainbow,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
# png(paste("../../../Results/20220117/contour_norm_marked_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_marked,2),color.palette = rainbow,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
# png(paste("../../../Results/20220117/contour_norm_all_unmarked_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_all,2),color.palette = rainbow,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
# png(paste("../../../Results/20220117/contour_norm_all_marked_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_all,2),color.palette = rainbow,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
#
# png(paste("../../../Results/20220117/contour_sum_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(Unmarked_imax30_NPM+Marked_imax30_NPM,2),color.palette = rainbow,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
#

# png(paste("../../../Results/20220117/contour_grey_norm_unmarked_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_unmarked,2),color.palette = grey.colors,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
# png(paste("../../../Results/20220117/contour_grey_norm_marked_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_marked,2),color.palette = grey.colors,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
# png(paste("../../../Results/20220117/contour_grey_norm_unmarked_all_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_all,2),color.palette = grey.colors,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
# png(paste("../../../Results/20220117/contour_grey_norm_marked_all_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_all,2),color.palette = grey.colors,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()
#
# png(paste("../../../Results/20220117/contour_grey_sum_l1m1d1t1_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
# filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(Unmarked_imax30_NPM+Marked_imax30_NPM,2),color.palette = grey.colors,xlab="class",
#                ylab="Log 10 k/m",key.title = "% Ribosome Load")
# dev.off()


}


Imax_Range_tau10_mean_Ribosomes_U <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau10_sd_Ribosomes_U <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau10_mean_Ribosomes_M <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau10_sd_Ribosomes_M <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])

#tau=10
count <- 1
par.mat.k.only.dt[,tau:=10]
for(k in c(2:50,100,150,200,250,300)){
  print(k)
  Unmarked_imax30_NPM <- apply(par.mat.k.only.dt,1,function(i) CalcUnmarkedClass(parameters = i,imax = k,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
  Marked_imax30_NPM <- apply(par.mat.k.only.dt,1,function(i) CalcMarkedClass(parameters = i,imax = k,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

  U_sum_NPM <- colSums(Unmarked_imax30_NPM)
  M_sum_NPM <- colSums(Marked_imax30_NPM)
  U_percent_NPM_unmarked <- t(apply(Unmarked_imax30_NPM,1,function(i) (i/U_sum_NPM)))*100
  M_percent_NPM_marked <- t(apply(Marked_imax30_NPM,1,function(i) (i/M_sum_NPM)))*100
  U_percent_NPM_all <- Unmarked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
  M_percent_NPM_all <- Marked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
  U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM_unmarked)%*%c(0:k))
  M_Mean_Ribosomes_NPM <- unlist(t(M_percent_NPM_marked)%*%c(0:k))
  difference <- apply(matrix(rep(c(0:k),ncol(Unmarked_imax30_NPM)),nrow=nrow(Unmarked_imax30_NPM)),1,function(i) i-U_Mean_Ribosomes_NPM/100)
  U_SD_Ribosomes_NPM <- diag(sqrt(unlist(t(U_percent_NPM_unmarked)%*%t((difference)^2))))
  difference <- apply(matrix(rep(c(0:k),ncol(Marked_imax30_NPM)),nrow=nrow(Marked_imax30_NPM)),1,function(i) i-M_Mean_Ribosomes_NPM/100)
  M_SD_Ribosomes_NPM <- diag(sqrt(unlist(t(M_percent_NPM_marked)%*%t((difference)^2))))


  Imax_Range_tau10_mean_Ribosomes_U[count,] <- U_Mean_Ribosomes_NPM
  Imax_Range_tau10_sd_Ribosomes_U[count,] <- U_SD_Ribosomes_NPM
  Imax_Range_tau10_mean_Ribosomes_M[count,] <- M_Mean_Ribosomes_NPM
  Imax_Range_tau10_sd_Ribosomes_M[count,] <- M_SD_Ribosomes_NPM

  count <- count + 1

#  U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM)%*%c(0:k))

  # png(paste("../../../Results/20220117/contour_norm_unmarked_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_unmarked,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_norm_marked_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_marked,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_norm_unmarked_all_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_all,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_norm_marked_all_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_all,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  #
  # png(paste("../../../Results/20220117/contour_sum_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(Unmarked_imax30_NPM+Marked_imax30_NPM,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  #
  #
  # png(paste("../../../Results/20220117/contour_grey_norm_unmarked_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_unmarked,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_grey_norm_marked_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_marked,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_grey_norm_unmarked_all_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_all,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_grey_norm_marked_all_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_all,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  #
  # png(paste("../../../Results/20220117/contour_grey_sum_l1m1d1t10_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(Unmarked_imax30_NPM+Marked_imax30_NPM,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()


}

Imax_Range_tau20_mean_Ribosomes_U <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau20_sd_Ribosomes_U <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau20_mean_Ribosomes_M <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])
Imax_Range_tau20_sd_Ribosomes_M <- matrix(0,nrow=length(c(2:50,100,150,200,250,300)),ncol=par.mat.k.only.dt[,.N])

count <- 1
#tau=20
par.mat.k.only.dt[,tau:=20]
for(k in c(2:50,100,150,200,250,300)){
  print(k)
  Unmarked_imax30_NPM <- apply(par.mat.k.only.dt,1,function(i) CalcUnmarkedClass(parameters = i,imax = k,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
  Marked_imax30_NPM <- apply(par.mat.k.only.dt,1,function(i) CalcMarkedClass(parameters = i,imax = k,rescale = T,rescale.val = "mu",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))

  U_sum_NPM <- colSums(Unmarked_imax30_NPM)
  M_sum_NPM <- colSums(Marked_imax30_NPM)
  U_percent_NPM_unmarked <- t(apply(Unmarked_imax30_NPM,1,function(i) (i/U_sum_NPM)))*100
  M_percent_NPM_marked <- t(apply(Marked_imax30_NPM,1,function(i) (i/M_sum_NPM)))*100
  U_percent_NPM_all <- Unmarked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
  M_percent_NPM_all <- Marked_imax30_NPM/(Marked_imax30_NPM+Unmarked_imax30_NPM)*100
  U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM_unmarked)%*%c(0:k))
  M_Mean_Ribosomes_NPM <- unlist(t(M_percent_NPM_marked)%*%c(0:k))
  difference <- apply(matrix(rep(c(0:k),ncol(Unmarked_imax30_NPM)),nrow=nrow(Unmarked_imax30_NPM)),1,function(i) i-U_Mean_Ribosomes_NPM/100)
  U_SD_Ribosomes_NPM <- diag(sqrt(unlist(t(U_percent_NPM_unmarked)%*%t((difference)^2))))
  difference <- apply(matrix(rep(c(0:k),ncol(Marked_imax30_NPM)),nrow=nrow(Marked_imax30_NPM)),1,function(i) i-M_Mean_Ribosomes_NPM/100)
  M_SD_Ribosomes_NPM <- diag(sqrt(unlist(t(M_percent_NPM_marked)%*%t((difference)^2))))


  Imax_Range_tau20_mean_Ribosomes_U[count,] <- U_Mean_Ribosomes_NPM
  Imax_Range_tau20_sd_Ribosomes_U[count,] <- U_SD_Ribosomes_NPM
  Imax_Range_tau20_mean_Ribosomes_M[count,] <- M_Mean_Ribosomes_NPM
  Imax_Range_tau20_sd_Ribosomes_M[count,] <- M_SD_Ribosomes_NPM

  count <- count + 1

#  U_Mean_Ribosomes_NPM <- unlist(t(U_percent_NPM)%*%c(0:k))

  # png(paste("../../../Results/20220117/contour_norm_unmarked_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_unmarked,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_norm_marked_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_marked,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_norm_unmarked_all_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_all,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_norm_marked_all_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_all,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  #
  # png(paste("../../../Results/20220117/contour_sum_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(Unmarked_imax30_NPM+Marked_imax30_NPM,2),color.palette = rainbow,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  #
  #
  # png(paste("../../../Results/20220117/contour_grey_norm_unmarked_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_unmarked,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_grey_norm_marked_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_marked,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_grey_norm_unmarked_all_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(U_percent_NPM_all,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  # png(paste("../../../Results/20220117/contour_grey_norm_marked_all_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(M_percent_NPM_all,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()
  #
  # png(paste("../../../Results/20220117/contour_grey_sum_l1m1d1t20_k1e-5_1e5_imax_",k,".png",sep=""),width = 6,height = 8,units = "in",res=600)
  # filled.contour(x=0:k,y=log(par.mat.k.only.dt[,kappa],10),z=round(Unmarked_imax30_NPM+Marked_imax30_NPM,2),color.palette = grey.colors,xlab="class",
  #                ylab="Log 10 k/m",key.title = "% Ribosome Load")
  # dev.off()


}




plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_mean_Ribosomes_U[1,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes)]/100,type="l",ylim=c(0,300))
for(i in 2:nrow(Imax_Range_tau1_mean_Ribosomes_U)){
  lines(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_mean_Ribosomes_U[i,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes)]/100,col=i)
}


plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_mean_Ribosomes_U[11,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes)]/100,type="l",ylim=c(0,50))
for(i in 12:49){
  lines(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_mean_Ribosomes_U[i,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes)]/100,col=i)
}

plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau10_mean_Ribosomes_U[11,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes)]/100,type="l",ylim=c(0,50))
for(i in 12:49){
  lines(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau10_mean_Ribosomes_U[i,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes)]/100,col=i)
}

plot(log(par.mat.k.only.dt[kappa>0.9,kappa],20),Imax_Range_tau10_mean_Ribosomes_U[11,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes)]/100,type="l",ylim=c(0,50))
for(i in 12:49){
  lines(log(par.mat.k.only.dt[kappa>0.9,kappa],20),Imax_Range_tau10_mean_Ribosomes_U[i,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes)]/100,col=i)
}


plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes_U)]/100,type="l",ylim=c(0,50))
lines(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau10_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau10_mean_Ribosomes_U)]/100,col="red")
lines(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau20_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau20_mean_Ribosomes_U)]/100,col="blue")
plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes_U)]/100-Imax_Range_tau10_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau10_mean_Ribosomes_U)]/100,type="l",ylab="")
plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes_U)]/100-Imax_Range_tau20_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau20_mean_Ribosomes_U)]/100,type="l",ylab="")
plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau10_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_mean_Ribosomes_U)]/100-Imax_Range_tau20_mean_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau20_mean_Ribosomes_U)]/100,type="l",ylab="")


plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_sd_Ribosomes_U)]/100,type="l",ylim=c(0,500))
lines(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau10_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau10_sd_Ribosomes_U)]/100,col="red")
lines(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau20_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau20_sd_Ribosomes_U)]/100,col="blue")
plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_sd_Ribosomes_U)]/100-Imax_Range_tau10_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau10_sd_Ribosomes_U)]/100,type="l",ylab="")
plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_sd_Ribosomes_U)]/100-Imax_Range_tau20_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau20_sd_Ribosomes_U)]/100,type="l",ylab="")
plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau10_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_sd_Ribosomes_U)]/100-Imax_Range_tau20_sd_Ribosomes_U[29,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau20_sd_Ribosomes_U)]/100,type="l",ylab="")

plot(log(par.mat.k.only.dt[kappa>0.9,kappa],10),Imax_Range_tau1_sd_Ribosomes_U[1,par.mat.k.only.dt[kappa>0.9,index][1]:ncol(Imax_Range_tau1_sd_Ribosomes_U)]/100,type="l",ylim=c(0,500))


Imax_Range_tau20_sd_Ribosomes_U
Imax_Range_tau20_mean_Ribosomes_M
Imax_Range_tau20_sd_Ribosomes_M


Imax_Range_tau20_mean_Ribosomes_U[,600]/100
imax_range <- c(2:50,100,150,200,250,300)



# focused examples --------------------------------------------------------
library(Biostrings)
At <- Biostrings::readDNAStringSet("../../Project_1/Results/20220117/CDS_arabidopsis.txt")
hist(width(At)/3)
At_lengths <- sort(width(At))
quantile(At_lengths,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))/27
hist(At_lengths[which(At_lengths<3000)]/3,breaks = seq(1,1010,10))

Hs <- Biostrings::readDNAStringSet("../../../Results/20220117/CDS_human(1).txt")
Hs_lengths <- sort(width(Hs))
quantile(Hs_lengths[-which(Hs_lengths<=9)],c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))/3
hist(Hs_lengths[which(Hs_lengths<3000)]/3,10,breaks = seq(1,1010,10))

Sc <- Biostrings::readDNAStringSet("../../../Results/20220117/CDS_yeast.txt")
Sc_lengths <- sort(width(Sc))
quantile(Sc_lengths,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))/3
hist(Sc_lengths[which(Sc_lengths<3000)]/3,breaks = seq(1,1010,10))


Dm <- Biostrings::readDNAStringSet("../../../Results/20220117/CDS_drosophila.txt")
Dm_lengths <- sort(width(Dm))
quantile(Dm_lengths,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))/3
hist(Dm_lengths[which(Dm_lengths<3000)]/3,breaks = seq(1,1010,10))


HL_At <- fread("../Data/Half life Sorenson.csv",stringsAsFactors = F,header=T)
hist(log(unlist(HL_At[which(unlist(HL_At[,3])<100000),3]),10))
median(unlist(HL_At[,3]),na.rm=T)
quantile(unlist(HL_At[which(unlist(HL_At[,3])<100000),3]),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
log(2)/(quantile(unlist(HL_At[which(unlist(HL_At[,3])<100000),3]),c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))*60)



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


U_sum <- colSums(Unmarked_imax30_fk)
M_sum <- colSums(Marked_imax30_fk)
U_percent_unmarked <- t(apply(Unmarked_imax30_fk,1,function(i) (i/U_sum)))*100
M_percent_marked <- t(apply(Marked_imax30_fk,1,function(i) (i/M_sum)))*100
U_percent_all <- Unmarked_imax30_fk/(Marked_imax30_fk+Unmarked_imax30_fk)*100
M_percent_all <- Marked_imax30_fk/(Marked_imax30_fk+Unmarked_imax30_fk)*100
U_Mean_Ribosomes <- unlist(t(U_percent_unmarked)%*%c(0:30))
M_Mean_Ribosomes <- unlist(t(M_percent_marked)%*%c(0:30))
difference <- apply(matrix(rep(c(0:30),ncol(Unmarked_imax30_fk)),nrow=nrow(Unmarked_imax30_fk)),1,function(i) i-U_Mean_Ribosomes/100)
U_SD_Ribosomes <- diag(sqrt(unlist(t(U_percent_unmarked)%*%t((difference)^2))))
difference <- apply(matrix(rep(c(0:30),ncol(Marked_imax30_fk)),nrow=nrow(Marked_imax30_fk)),1,function(i) i-M_Mean_Ribosomes/100)
M_SD_Ribosomes <- diag(sqrt(unlist(t(M_percent_marked)%*%t((difference)^2))))
peak_unmarked <- apply(Unmarked_imax30_fk,2,function(j) which(j==max(j)))

U_Mean_Ribosomes_mat <- matrix(U_Mean_Ribosomes,nrow=7)/100
rownames(U_Mean_Ribosomes_mat) <- unique(unlist(Focused.par.setk[1:7,kappa]))
colnames(U_Mean_Ribosomes_mat) <- unique(paste("m=",unlist(Focused.par.freek[,2]),", t=",unlist(Focused.par.freek[,4]),sep=""))

APPROX <- matrix(numeric(),nrow=nrow(Unmarked_imax30_fk),ncol = ncol(Unmarked_imax30_fk))
for(i in 1:nrow(Focused.par.freek)){
  kappa <- Focused.par.freek[i,kappa]
  tau <- Focused.par.freek[i,tau]
  mu <- Focused.par.freek[i,mu]
  kappa <- rep(kappa,31)*(1-0:30/30)
  tau <- rep(tau,31)*(0:30)

    for(j in 1:31){
    if(j==1){
      APPROX[j,i] <- sum(log(tau[(j+1):31],10))

    } else if(j==31){
      APPROX[j,i] <- sum(log(kappa[1:(j-1)],10))

    } else {
      APPROX[j,i] <- sum(log(kappa[1:(j-1)],10))+sum(log(tau[(j+1):30],10))

    }

  }

}


# result <- apply(10^APPROX,2,function(i) i/sum(i))
# result <- lapply(1:ncol(result),function(i) result[,i]/unlist(Focused.par.freek[,mu])[i])
# result <- matrix(unlist(result),nrow=31)
#
# plot(result[,8],type='l')
# lines(result[,9],col=2)
# lines(result[,9],col=3)
# lines(result[,9],col=4)
# lines(result[,9],col=2)
# lines(result[,10],col=3)
# lines(result[,11],col=4)
# lines(result[,12],col=5)
# lines(result[,13],col=6)
# lines(result[,14],col=7)
# plot(Unmarked_imax30_fk[,8],type='l')
# lines(Unmarked_imax30_fk[,9],col=2)
# lines(Unmarked_imax30_fk[,9],col=3)
# lines(Unmarked_imax30_fk[,9],col=4)
# lines(Unmarked_imax30_fk[,9],col=2)
# lines(Unmarked_imax30_fk[,10],col=3)
# lines(Unmarked_imax30_fk[,11],col=4)
# lines(Unmarked_imax30_fk[,12],col=5)
# lines(Unmarked_imax30_fk[,13],col=6)
# lines(Unmarked_imax30_fk[,14],col=7)
#
# plot((colSums((result- U_percent_unmarked/100)^2)),pch=20)
# abline(v=seq(0.5,245.5,7),col="red")
# plot(result[,245],type="l")
# lines(Unmarked_imax30_fk[,245],col="red")


series.mat.plot(parameters=Focused.par.freek,Unmarked=Unmarked_imax30_fk,Marked = Marked_imax30_fk,
                group.by =length(unique(Focused.par.setk[,kappa])),set.params = c("kappa","mu","tau"),file.loc="../Results/20220210",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.freek,Unmarked=Unmarked_imax30_fk,Marked = Marked_imax30_fk,
                group.by =length(unique(Focused.par.setk[,kappa])),set.params = c("kappa","mu","tau"),file.loc="../Results/20220210",log=F,beside=T,xaxt.step=5,global.ymax=T)


Unmarked_imax30_sk <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))
Marked_imax30_sk <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))



U_sum <- colSums(Unmarked_imax30_sk)
M_sum <- colSums(Marked_imax30_sk)
U_percent_unmarked <- t(apply(Unmarked_imax30_sk,1,function(i) (i/U_sum)))*100
M_percent_marked <- t(apply(Marked_imax30_sk,1,function(i) (i/M_sum)))*100
U_percent_all <- Unmarked_imax30_sk/(Marked_imax30_sk+Unmarked_imax30_sk)*100
M_percent_all <- Marked_imax30_sk/(Marked_imax30_sk+Unmarked_imax30_sk)*100
U_Mean_Ribosomes <- unlist(t(U_percent_unmarked)%*%c(0:30))
M_Mean_Ribosomes <- unlist(t(M_percent_marked)%*%c(0:30))
difference <- apply(matrix(rep(c(0:30),ncol(Unmarked_imax30_sk)),nrow=nrow(Unmarked_imax30_sk)),1,function(i) i-U_Mean_Ribosomes/100)
U_SD_Ribosomes <- diag(sqrt(unlist(t(U_percent_unmarked)%*%t((difference)^2))))
difference <- apply(matrix(rep(c(0:30),ncol(Marked_imax30_sk)),nrow=nrow(Marked_imax30_sk)),1,function(i) i-M_Mean_Ribosomes/100)
M_SD_Ribosomes <- diag(sqrt(unlist(t(M_percent_marked)%*%t((difference)^2))))
peak_unmarked <- apply(Unmarked_imax30_sk,2,function(j) which(j==max(j)))

U_Mean_Ribosomes_mat <- matrix(U_Mean_Ribosomes,nrow=7)/100
rownames(U_Mean_Ribosomes_mat) <- unique(unlist(Focused.par.setk[1:7,kappa]))
colnames(U_Mean_Ribosomes_mat) <- unique(paste("m=",unlist(Focused.par.freek[,2]),", t=",unlist(Focused.par.freek[,4]),sep=""))

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax30_sk,Marked = Marked_imax30_sk,
                group.by =length(unique(Focused.par.setk[,kappa])),set.params = c("kappa","mu","tau"),file.loc="../Results/20220210",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax30_sk,Marked = Marked_imax30_sk,
                group.by =length(unique(Focused.par.setk[,kappa])),set.params = c("kappa","mu","tau"),file.loc="../Results/20220210",log=F,beside=T,xaxt.step=5,global.ymax=T)


# Focused parameters with imax --------------------------------------------
#Utilizing set k parameter set
#Tau will be 0.1aa/s, 1aa/s, 5aa/s, 10aa/s, 20aa/s
#Kappa will be 0.1%, 1%, 25%, 50%, 75%, 100% and 200% of Tau = 1aa/s for each tau
#Mu will be based on the following half lives from sorenson
#12min, 23min, 59min, 109min, 220min, 830min, 4018min The 1%, 5%, 25%, 50%, 75%, 95%, and 99% quantiles of arabidopsis
#decay constant mu = ln(2)/half life in seconds
#imax will be 6, 10 ,24, 40, 60, 112, 181 The 1%, 5%, 25%, 50%, 75%, 95%, and 99% quantiles of arabidopsis


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



Ribo_percent <- function(Data){
  tmp <- apply(Data,2,function(i) i/sum(i))
  return(tmp)
}
col <- rainbow(50)
col <- c(col,rep(col[50],51))
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,2/9,0.0001)),0.0000001,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,4,0.005)),0.0001,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,4,0.005)),1,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,0.1,0.0001)),0.0001,1,1,100000,1),par.names = par.names)
filled.contour(x=0:10,y=log(tmp.par[,"kappa"],10),z=Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 10,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))),color.palette = rainbow)
filled.contour(x=0:30,y=log(tmp.par[,"kappa"],10),z=Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))),color.palette = rainbow)
filled.contour(x=0:181,y=log(tmp.par[,"kappa"],10),z=Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau")))),nlevels = 101,col = col,
               plot.axes = {axis(1,at=c(0,seq(10,180,10)),labels = c(0,seq(10,180,10)),las=2)
                 axis(2,at=seq(-4,1,1),labels = seq(-4,1,1))},key.title = "Percent Ribosome load",xlab="Ribosome Class",ylab=expression(paste(log[10],kappa,sep="")))

filled.contour(x=0:10,y=log(tmp.par[,"kappa"],10),z=Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau"))))[1:11,],nlevels = 101,col = col,
               plot.axes = {axis(1,at=c(0:10,seq(10,180,10)),labels = c(0:10,seq(10,180,10)),las=2)
                 axis(2,at=seq(-4,1,1),labels = seq(-4,1,1))},key.title = "Percent Ribosome load",xlab="Ribosome Class",ylab=expression(paste(log[10],kappa,sep="")))

#no int
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,0.11111,0.0001)),0.0000001,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,4,0.005)),0.0001,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,4,0.005)),1,1,1,100000,1),par.names = par.names)
tmp.par <- CreateParMatrix(parameter.range = list(c(seq(0.0001,0.1,0.0001)),0.01,1,1,100000,1),par.names = par.names)
filled.contour(x=0:10,y=log(tmp.par[,"kappa"],10),z=Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 10,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))),color.palette = rainbow)
filled.contour(x=0:30,y=log(tmp.par[,"kappa"],10),z=Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))),color.palette = rainbow)
filled.contour(x=0:181,y=log(tmp.par[,"kappa"],10),z=Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))),nlevels = 101,col = col,
               plot.axes = {axis(1,at=c(0,seq(10,180,10)),labels = c(0,seq(10,180,10)),las=2)
                 axis(2,at=seq(-4,1,1),labels = seq(-4,1,1))},key.title = "Percent Ribosome load",xlab="Ribosome Class",ylab=expression(paste(log[10],kappa,sep="")))

filled.contour(x=0:10,y=log(tmp.par[,"kappa"],10),z=Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))[1:11,],nlevels = 101,col = col,
               plot.axes = {axis(1,at=c(0:10,seq(10,180,10)),labels = c(0:10,seq(10,180,10)),las=2)
                 axis(2,at=seq(-4,1,1),labels = seq(-4,1,1))},key.title = "Percent Ribosome load",xlab="Ribosome Class",ylab=expression(paste(log[10],kappa,sep="")))
plot(tmp.par[,"kappa"],colSums(Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))[-1,])*100,ylab="%loaded transcripts",xlab="kappa")


Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("kappa","tau"))))


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


tmp<-CalcUnmarkedAdj_noM(parameters = unlist(tmp.par[1,]),imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = F,which.interfere = c("tau"))
#tmp <- tmp%*%c(1,rep(0,10))
tmp <- tmp/sum(tmp)*100

Unmarked_no_M <- apply(tmp.par,1,function(i) CalcUnmarkedAdj_noM(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = F,which.interfere = c("tau")))

plot(tmp.par[,"kappa"][1:63],colSums(Ribo_percent(Unmarked_no_M)[-1,1:63])*100,ylab="%loaded transcripts",xlab="kappa")
plot(tmp.par[,"kappa"][1:63],apply(Unmarked_no_M,2,Ribo_mean)[1:63],ylab="mean ribosome load",xlab="kappa")
abline(v=1/(30*9))
abline(h=15)
plot(0:30,Ribo_percent(Unmarked_no_M[,1:2])[,1])

Unmarked_no_M <- apply(tmp.par,1,function(i) CalcUnmarkedAdj_noM(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

plot(tmp.par[,"kappa"][1:63],colSums(Ribo_percent(Unmarked_no_M)[-1,1:63])*100,ylab="%loaded transcripts",xlab="kappa")
plot(tmp.par[,"kappa"][1:63],apply(Unmarked_no_M,2,Ribo_mean)[1:63],ylab="mean ribosome load",xlab="kappa")
abline(v=1/(30*9))
abline(h=15)
plot(0:30,Ribo_percent(Unmarked_no_M[,1:63])[,37])

Unmarked_with_M <- apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

plot(tmp.par[,"kappa"][1:63],colSums(Ribo_percent(Unmarked_no_M)[-1,1:63])*100,ylab="%loaded transcripts",xlab="kappa")
lines(tmp.par[,"kappa"][1:63],colSums(Ribo_percent(Unmarked_with_M)[-1,1:63])*100,col="red")
plot(tmp.par[,"kappa"][1:63],apply(Unmarked_no_M,2,Ribo_mean)[1:63],ylab="mean ribosome load",xlab="kappa")
lines(tmp.par[,"kappa"][1:63],apply(Unmarked_no_M,2,Ribo_mean)[1:63],col="red")
abline(v=1/(30*9))
abline(h=15)
plot(0:30,Ribo_percent(Unmarked_no_M[,1:63])[,37])

tmp.par[,"mu"]<-0.0000001
plot(tmp.par[,"kappa"],apply(Unmarked_no_M,2,Ribo_mean),type="l")
lines(tmp.par[,"kappa"],tmp.par[,"kappa"]/(1/270),col="red")
abline(a=0,b=270,col="red")
lines(tmp.par[,"kappa"],apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean),type="l",col= "red",lty=2)
tmp.dat <- apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean)
lm(formula = y~x-1,data = data.frame("x"=tmp.par[1:500,"kappa"],"y"=tmp.dat[1:500]))

tmp.par[,"mu"]<-0.000001
lines(tmp.par[,"kappa"],apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean),type="l",col= "blue",lty=2)
tmp.dat <- apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean)
lm(formula = y~x-1,data = data.frame("x"=tmp.par[1:500,"kappa"],"y"=tmp.dat[1:500]))

tmp.par[,"mu"]<-0.00001
lines(tmp.par[,"kappa"],apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean),type="l",col= "green",lty=2)
tmp.dat <- apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean)
lm(formula = y~x-1,data = data.frame("x"=tmp.par[1:500,"kappa"],"y"=tmp.dat[1:500]))

tmp.par[,"mu"]<-0.0001
lines(tmp.par[,"kappa"],apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean),type="l",col= "orange",lty=2)
tmp.dat <- apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean)
lm(formula = y~x-1,data = data.frame("x"=tmp.par[1:500,"kappa"],"y"=tmp.dat[1:500]))

tmp.par[,"mu"]<-0.001
lines(tmp.par[,"kappa"],apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean),type="l",col= "purple",lty=2)
tmp.dat <- apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean)
lm(formula = y~x-1,data = data.frame("x"=tmp.par[1:500,"kappa"],"y"=tmp.dat[1:500]))
abline(a=0,b=212.6)

tmp.par[,"mu"]<-0.01
lines(tmp.par[,"kappa"],apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean),type="l",col= "grey",lty=2)
tmp.dat <- apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean)
lm(formula = y~x-1,data = data.frame("x"=tmp.par[1:500,"kappa"],"y"=tmp.dat[1:500]))

tmp.par[,"mu"]<-0.1
lines(tmp.par[,"kappa"],apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean),type="l",col= "magenta",lty=2)
tmp.dat <- apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean)
lm(formula = y~x-1,data = data.frame("x"=tmp.par[1:500,"kappa"],"y"=tmp.dat[1:500]))

tmp.par[,"mu"]<-1
lines(tmp.par[,"kappa"],apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean),type="l",col= "black",lty=3)
tmp.dat <- apply(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))),2,Ribo_mean)
lm(formula = y~x-1,data = data.frame("x"=tmp.par[1:500,"kappa"],"y"=tmp.dat[1:500]))

plot(c(0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1),c(270,269.9,269.3,262.9,212.6,72.97,9.643,0.9963))


plot(tmp.par[,"kappa"],colSums(Ribo_percent(apply(tmp.par,1,function(i) CalcUnmarkedClass(parameters = i,imax = 30,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))))[-1,]),type="l",col= "black",lty=3,ylab="%transcripts loaded")
lines()

Ribo_mean <- function(Data){
  if(is.null(dim(Data))&length(unlist(Data))>1){
    mean.tmp <- unlist((Data/sum(Data))%*%c(0:(length(Data)-1)))
  }else{
    Data_pct <- apply(Data,2,function(i) i/sum(i))
    mean.tmp <- unlist(t(Data_pct)%*%c(0:(nrow(Data)-1)))
  }
  return(mean.tmp)
}

U_Mean_Ribosomes_imax6 <- Ribo_mean(Unmarked_imax6_sk)
U_Mean_Ribosomes_imax10 <- Ribo_mean(Unmarked_imax10_sk)
U_Mean_Ribosomes_imax24 <- Ribo_mean(Unmarked_imax24_sk)
U_Mean_Ribosomes_imax40 <- Ribo_mean(Unmarked_imax40_sk)
U_Mean_Ribosomes_imax60 <- Ribo_mean(Unmarked_imax60_sk)
U_Mean_Ribosomes_imax112 <- Ribo_mean(Unmarked_imax112_sk)
U_Mean_Ribosomes_imax181 <- Ribo_mean(Unmarked_imax181_sk)


max(c(U_Mean_Ribosomes_imax6))
max(c(U_Mean_Ribosomes_imax10))
max(c(U_Mean_Ribosomes_imax24))
max(c(U_Mean_Ribosomes_imax40))
max(c(U_Mean_Ribosomes_imax60))
max(c(U_Mean_Ribosomes_imax112))
max(c(U_Mean_Ribosomes_imax181))

mean(c(U_Mean_Ribosomes_imax6))
mean(c(U_Mean_Ribosomes_imax10))
mean(c(U_Mean_Ribosomes_imax24))
mean(c(U_Mean_Ribosomes_imax40))
mean(c(U_Mean_Ribosomes_imax60))
mean(c(U_Mean_Ribosomes_imax112))
mean(c(U_Mean_Ribosomes_imax181))

diff(c(mean(c(U_Mean_Ribosomes_imax6)),mean(c(U_Mean_Ribosomes_imax10)),
       mean(c(U_Mean_Ribosomes_imax24)),mean(c(U_Mean_Ribosomes_imax40)),
       mean(c(U_Mean_Ribosomes_imax60)),mean(c(U_Mean_Ribosomes_imax112)),
       mean(c(U_Mean_Ribosomes_imax181))))

plot(c(6,10,24,40,60,112,181),c(mean(c(U_Mean_Ribosomes_imax6)),mean(c(U_Mean_Ribosomes_imax10)),
                                      mean(c(U_Mean_Ribosomes_imax24)),mean(c(U_Mean_Ribosomes_imax40)),
                                      mean(c(U_Mean_Ribosomes_imax60)),mean(c(U_Mean_Ribosomes_imax112)),
                                      mean(c(U_Mean_Ribosomes_imax181))))



series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax6_sk,Marked = Marked_imax6_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax10_sk,Marked = Marked_imax10_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax24_sk,Marked = Marked_imax24_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax40_sk,Marked = Marked_imax40_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax60_sk,Marked = Marked_imax60_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax112_sk,Marked = Marked_imax112_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax181_sk,Marked = Marked_imax181_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)




plot(Unmarked_imax6_sk[,7]/sum(Unmarked_imax6_sk[,7]),ylim = c(0,1),xlim=c(0,40),type="l")
lines(Unmarked_imax10_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=2)
lines(Unmarked_imax24_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=3)
lines(Unmarked_imax40_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=4)
lines(Unmarked_imax60_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=5)
lines(Unmarked_imax112_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=6)
lines(Unmarked_imax181_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=7)

for(i in seq(1,7,1)){
  plot(Unmarked_imax6_sk[,i]/sum(Unmarked_imax6_sk[,i]),ylim = c(0,1),xlim=c(0,40),type="l",main=i)
  lines(Unmarked_imax10_sk[,i]/sum(Unmarked_imax10_sk[,i]),col=2)
  lines(Unmarked_imax24_sk[,i]/sum(Unmarked_imax24_sk[,i]),col=3)
  lines(Unmarked_imax40_sk[,i]/sum(Unmarked_imax40_sk[,i]),col=4)
  lines(Unmarked_imax60_sk[,i]/sum(Unmarked_imax60_sk[,i]),col=5)
  lines(Unmarked_imax112_sk[,i]/sum(Unmarked_imax112_sk[,i]),col=6)
  lines(Unmarked_imax181_sk[,i]/sum(Unmarked_imax181_sk[,i]),col=7)
  legend("topright",legend = c(paste("kappa=",Focused.par.freek[i,kappa],", tau=",Focused.par.freek[i,tau],sep=""),
                               "imax=6","imax=10","imax=24","imax=40","imax=60","imax=112","imax=181"),fill = 0:7,border = F)
}





for(i in seq(1,7,1)){
  plot(Unmarked_imax6_sk[,i]/sum(Unmarked_imax6_sk[,i]),ylim = c(0,1),xlim=c(0,40),type="l",main=i)
  lines(Unmarked_imax10_sk[,i*2]/sum(Unmarked_imax10_sk[,i*2]),col=2)
  lines(Unmarked_imax24_sk[,i*3]/sum(Unmarked_imax24_sk[,i*3]),col=3)
  lines(Unmarked_imax40_sk[,i*4]/sum(Unmarked_imax40_sk[,i*4]),col=4)
  lines(Unmarked_imax60_sk[,i*5]/sum(Unmarked_imax60_sk[,i*5]),col=5)
  lines(Unmarked_imax112_sk[,i*6]/sum(Unmarked_imax112_sk[,i*6]),col=6)
  lines(Unmarked_imax181_sk[,i*7]/sum(Unmarked_imax181_sk[,i*7]),col=7)
  legend("topright",legend = c(paste("kappa=",Focused.par.freek[i,kappa],", tau=",Focused.par.freek[i,tau],sep=""),
                               "imax=6","imax=10","imax=24","imax=40","imax=60","imax=112","imax=181"),fill = 0:7,border = F)
}



for(i in 1){
  plot(Unmarked_imax6_sk[,i]/sum(Unmarked_imax6_sk[,i]),ylim = c(0,1),xlim=c(0,40),type="l",main=i)
  lines(Unmarked_imax6_sk[,i+7]/sum(Unmarked_imax6_sk[,i*2]),col=2)
  lines(Unmarked_imax6_sk[,i+14]/sum(Unmarked_imax6_sk[,i*3]),col=3)
  lines(Unmarked_imax6_sk[,i+21]/sum(Unmarked_imax6_sk[,i*4]),col=4)
  lines(Unmarked_imax6_sk[,i]/sum(Unmarked_imax6_sk[,i*5]),col=5)
  lines(Unmarked_imax6_sk[,i*6]/sum(Unmarked_imax6_sk[,i*6]),col=6)
  lines(Unmarked_imax6_sk[,i*7]/sum(Unmarked_imax6_sk[,i*7]),col=7)
  legend("topright",legend = c(paste("kappa=",Focused.par.freek[i,kappa],", tau=",Focused.par.freek[i,tau],sep=""),
                               "imax=6","imax=10","imax=24","imax=40","imax=60","imax=112","imax=181"),fill = 0:7,border = F)
}

c(6,10,24,40,60,112,181)




# Focused parameters no interference --------------------------------------


Unmarked_imax6_NI <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 6,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))
Marked_imax6_NI <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 6,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

Unmarked_imax10_NI <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 10,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))
Marked_imax10_NI <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 10,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

Unmarked_imax24_NI <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 24,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))
Marked_imax24_NI <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 24,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

Unmarked_imax40_NI <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 40,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))
Marked_imax40_NI <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 40,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

Unmarked_imax60_NI <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 60,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))
Marked_imax60_NI <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 60,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

Unmarked_imax112_NI <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 112,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))
Marked_imax112_NI <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 112,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))

Unmarked_imax181_NI <- apply(Focused.par.setk,1,function(i) CalcUnmarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))
Marked_imax181_NI <- apply(Focused.par.setk,1,function(i) CalcMarkedClass(parameters = i,imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau")))



U_Mean_Ribosomes_imax6_NI <- Ribo_mean(Unmarked_imax6_NI)
U_Mean_Ribosomes_imax10_NI <- Ribo_mean(Unmarked_imax10_NI)
U_Mean_Ribosomes_imax24_NI <- Ribo_mean(Unmarked_imax24_NI)
U_Mean_Ribosomes_imax40_NI <- Ribo_mean(Unmarked_imax40_NI)
U_Mean_Ribosomes_imax60_NI <- Ribo_mean(Unmarked_imax60_NI)
U_Mean_Ribosomes_imax112_NI <- Ribo_mean(Unmarked_imax112_NI)
U_Mean_Ribosomes_imax181_NI <- Ribo_mean(Unmarked_imax181_NI)


max(c(U_Mean_Ribosomes_imax6))
max(c(U_Mean_Ribosomes_imax10))
max(c(U_Mean_Ribosomes_imax24))
max(c(U_Mean_Ribosomes_imax40))
max(c(U_Mean_Ribosomes_imax60))
max(c(U_Mean_Ribosomes_imax112))
max(c(U_Mean_Ribosomes_imax181))

mean(c(U_Mean_Ribosomes_imax6))
mean(c(U_Mean_Ribosomes_imax10))
mean(c(U_Mean_Ribosomes_imax24))
mean(c(U_Mean_Ribosomes_imax40))
mean(c(U_Mean_Ribosomes_imax60))
mean(c(U_Mean_Ribosomes_imax112))
mean(c(U_Mean_Ribosomes_imax181))

diff(c(mean(c(U_Mean_Ribosomes_imax6)),mean(c(U_Mean_Ribosomes_imax10)),
       mean(c(U_Mean_Ribosomes_imax24)),mean(c(U_Mean_Ribosomes_imax40)),
       mean(c(U_Mean_Ribosomes_imax60)),mean(c(U_Mean_Ribosomes_imax112)),
       mean(c(U_Mean_Ribosomes_imax181))))

plot(c(6,10,24,40,60,112,181),c(mean(c(U_Mean_Ribosomes_imax6)),mean(c(U_Mean_Ribosomes_imax10)),
                                mean(c(U_Mean_Ribosomes_imax24)),mean(c(U_Mean_Ribosomes_imax40)),
                                mean(c(U_Mean_Ribosomes_imax60)),mean(c(U_Mean_Ribosomes_imax112)),
                                mean(c(U_Mean_Ribosomes_imax181))))


series.mat.plot <- function(parameters,Unmarked,Marked,group.by ="kappa",set.params = c("mu","tau"),file.loc="",log=F,beside=T,xaxt.step=5,global.ymax=F){
  series.mat <- t(matrix(1:nrow(parameters),nrow=nrow(unique(parameters[,..group.by]))))
  for(i in 1:(nrow(parameters)/(nrow(unique(parameters[,..group.by]))))){
    print(i)
    series <- unlist(series.mat[i,])
    pars=parameters
    leg.titles <- set.params
    leg.text <- as.character(unlist(apply(t(apply(pars[series[1],..set.params],1,function(i) paste(leg.titles,formatC(unlist(i), format = "e", digits = 0), sep="_"))),1,function(j) paste(j,collapse = "_") )))
    png(paste(file.loc,"/",group.by,"series","_",leg.text,"_globalymax",global.ymax,"_log",log,"_imax",nrow(Unmarked)-1,".png",sep=""),width = 8,height = 6,res = 600,units = "in")
    plot.bars.mat(unmarked = Unmarked,marked = Marked,series = series,
                  pars=parameters,main="Sum classes. l=1, d=10^5",log=log,beside=beside,xaxt.step = xaxt.step,global.ymax = global.ymax)
    dev.off()
  }

}

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax6_sk,Marked = Marked_imax6_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax10_sk,Marked = Marked_imax10_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax24_sk,Marked = Marked_imax24_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax40_sk,Marked = Marked_imax40_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax60_sk,Marked = Marked_imax60_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax112_sk,Marked = Marked_imax112_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)

series.mat.plot(parameters=Focused.par.setk,Unmarked=Unmarked_imax181_sk,Marked = Marked_imax181_sk,
                group.by ="kappa",set.params = c("mu","tau"),file.loc="../../../Results/20220125",log=F,beside=T,xaxt.step=5,global.ymax=F)




plot(Unmarked_imax6_sk[,7]/sum(Unmarked_imax6_sk[,7]),ylim = c(0,1),xlim=c(0,40),type="l")
lines(Unmarked_imax10_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=2)
lines(Unmarked_imax24_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=3)
lines(Unmarked_imax40_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=4)
lines(Unmarked_imax60_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=5)
lines(Unmarked_imax112_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=6)
lines(Unmarked_imax181_sk[,7]/sum(Unmarked_imax6_sk[,7]),col=7)

for(i in seq(1,7,1)){
  plot(Unmarked_imax6_sk[,i]/sum(Unmarked_imax6_sk[,i]),ylim = c(0,1),xlim=c(0,40),type="l",main=i)
  lines(Unmarked_imax10_sk[,i]/sum(Unmarked_imax10_sk[,i]),col=2)
  lines(Unmarked_imax24_sk[,i]/sum(Unmarked_imax24_sk[,i]),col=3)
  lines(Unmarked_imax40_sk[,i]/sum(Unmarked_imax40_sk[,i]),col=4)
  lines(Unmarked_imax60_sk[,i]/sum(Unmarked_imax60_sk[,i]),col=5)
  lines(Unmarked_imax112_sk[,i]/sum(Unmarked_imax112_sk[,i]),col=6)
  lines(Unmarked_imax181_sk[,i]/sum(Unmarked_imax181_sk[,i]),col=7)
  legend("topright",legend = c(paste("kappa=",Focused.par.freek[i,kappa],", tau=",Focused.par.freek[i,tau],sep=""),
                               "imax=6","imax=10","imax=24","imax=40","imax=60","imax=112","imax=181"),fill = 0:7,border = F)
}





for(i in seq(1,7,1)){
  plot(Unmarked_imax6_sk[,i]/sum(Unmarked_imax6_sk[,i]),ylim = c(0,1),xlim=c(0,40),type="l",main=i)
  lines(Unmarked_imax10_sk[,i*2]/sum(Unmarked_imax10_sk[,i*2]),col=2)
  lines(Unmarked_imax24_sk[,i*3]/sum(Unmarked_imax24_sk[,i*3]),col=3)
  lines(Unmarked_imax40_sk[,i*4]/sum(Unmarked_imax40_sk[,i*4]),col=4)
  lines(Unmarked_imax60_sk[,i*5]/sum(Unmarked_imax60_sk[,i*5]),col=5)
  lines(Unmarked_imax112_sk[,i*6]/sum(Unmarked_imax112_sk[,i*6]),col=6)
  lines(Unmarked_imax181_sk[,i*7]/sum(Unmarked_imax181_sk[,i*7]),col=7)
  legend("topright",legend = c(paste("kappa=",Focused.par.freek[i,kappa],", tau=",Focused.par.freek[i,tau],sep=""),
                               "imax=6","imax=10","imax=24","imax=40","imax=60","imax=112","imax=181"),fill = 0:7,border = F)
}


c(6,10,24,40,60,112,181)






for(i in 1:7){
plot(Unmarked_imax6_NI[,i]/sum(Unmarked_imax6_NI[,i]),type="l",ylim=c(0,1),xlim=c(0,40),main=i)
lines(Unmarked_imax6_sk[,i]/sum(Unmarked_imax6_sk[,i]),col=2,lty=2)
lines(Unmarked_imax10_NI[,i]/sum(Unmarked_imax10_NI[,i]),col=3)
lines(Unmarked_imax10_sk[,i]/sum(Unmarked_imax10_sk[,i]),col=4,lty=2)
lines(Unmarked_imax24_NI[,i]/sum(Unmarked_imax24_NI[,i]),col=5)
lines(Unmarked_imax24_sk[,i]/sum(Unmarked_imax24_sk[,i]),col=6,lty=2)
lines(Unmarked_imax40_NI[,i]/sum(Unmarked_imax40_NI[,i]),col=7)
lines(Unmarked_imax40_sk[,i]/sum(Unmarked_imax40_sk[,i]),col=8,lty=2)
lines(Unmarked_imax60_NI[,i]/sum(Unmarked_imax60_NI[,i]),col=9)
lines(Unmarked_imax60_sk[,i]/sum(Unmarked_imax60_sk[,i]),col=10,lty=2)
lines(Unmarked_imax112_NI[,i]/sum(Unmarked_imax112_NI[,i]),col=11)
lines(Unmarked_imax112_sk[,i]/sum(Unmarked_imax112_sk[,i]),col=12,lty=2)
lines(Unmarked_imax181_NI[,i]/sum(Unmarked_imax181_NI[,i]),col=13)
lines(Unmarked_imax181_sk[,i]/sum(Unmarked_imax181_sk[,i]),col=14,lty=2)

}

param.test <- Focused.par.setk[1,]

param.test[,kappa:=100000000]

xt <- CalcUnmarkedClass(parameters = unlist(param.test[1,]),imax = 181,rescale = F,rescale.val = "tau",par.names = par.names,interference = T,which.interfere = c("tau"))

t(xt/sum(xt))%*%(0:181)



plot(x=Focused.par.setk[,`k/t`],U_Mean_Ribosomes_imax181_NI)
abline(a=0,b=1)
plot(x=Focused.par.setk[,`k/t`],U_Mean_Ribosomes_imax181)
abline(a=0,b=1,col="red")
x=sort(unlist(Focused.par.setk[,`k/t`]))
y=U_Mean_Ribosomes_imax181[order(unlist(Focused.par.setk[,`k/t`]))]
m<-nls(y~a*x+b,start = list("a"=1,"b"=-1))
cor(y,m$m$predict(sort(x)))
plot(x,y)
lines(sort(x),m$m$predict(sort(x)))


n<-predict(m)
plot(U_Mean_Ribosomes_imax181[order(unlist(Focused.par.setk[,`k/t`]))]-sort(unlist(Focused.par.setk[,`k/t`])))


log(x)

lines(x=seq(0,20,0.1),y=)
max(Focused.par.setk[,`k/t`])
min(Focused.par.setk[,`k/t`])

