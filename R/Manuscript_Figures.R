source("R/Polysome_functions_RAUC.R")


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
U_mode_imax30_plot1

#plot

plot(plot1.params[tau==0.1,kappa],U_sd_imax30_plot1[plot1.params[tau==0.1,index]],type="l",ylim=c(0,max(U_sd_imax30_plot1)))
counter <- 2
counter <- 2
for(i in par.ranges[[4]][1:4]){
  lines(plot1.params[tau==i,kappa],U_sd_imax30_plot1[plot1.params[tau==i,index]],lty=counter)
  counter <- counter + 1
}
counter <- 1
for(i in par.ranges[[4]][5:8]){
  lines(plot1.params[tau==i,kappa],U_sd_imax30_plot1[plot1.params[tau==i,index]],lty=counter,col="grey50")
  counter <- counter + 1
}


plot(plot1.params[tau==0.1,kappa],U_mean_imax30_plot1[plot1.params[tau==0.1,index]],type="l",ylim=c(0,max(U_mean_imax30_plot1)),
     las=1,ylab="Mean Ribsome Load",xlab= bquote(kappa),
     main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
counter <- 2
for(i in par.ranges[[4]][1:4]){
  lines(plot1.params[tau==i,kappa],U_mean_imax30_plot1[plot1.params[tau==i,index]],lty=counter)
  counter <- counter + 1
}
counter <- 1
for(i in par.ranges[[4]][5:8]){
  lines(plot1.params[tau==i,kappa],U_mean_imax30_plot1[plot1.params[tau==i,index]],lty=counter,col="grey50")
  counter <- counter + 1
}


# plot(plot1.params[tau==0.1,kappa][1:1112],U_mean_imax30_plot1[plot1.params[tau==0.1,index]][1:1112],type="l",ylim=c(0,max(U_mean_imax30_plot1)),
#      las=1,ylab="Mean Ribsome Load",xlab= bquote(kappa),
#      main=bquote("No interference, imax=30, "*lambda*"="*1*", "*mu*"="*0*", "*delta*"="*10^5))
# counter <- 2
# for(i in par.ranges[[4]][1:4]){
#   lines(plot1.params[tau==i,kappa][1:1112],U_mean_imax30_plot1[plot1.params[tau==i,index]][1:1112],lty=counter)
#   counter <- counter + 1
# }
# counter <- 1
# for(i in par.ranges[[4]][5:8]){
#   lines(plot1.params[tau==i,kappa][1:1112],U_mean_imax30_plot1[plot1.params[tau==i,index]][1:1112],lty=counter,col="grey50")
#   counter <- counter + 1
# }

abline(v=.1/9,col="red")
points(plot1.params[3334,kappa],U_mean_imax30_plot1[plot1.params[tau==0.1,index]][3334],pch=20,col="red")
points(plot1.params[324,kappa],U_mean_imax30_plot1[plot1.params[tau==0.1,index]][324],pch=20,col="red")
points(plot1.params[112,kappa],U_mean_imax30_plot1[plot1.params[tau==0.1,index]][112],pch=20,col="red")

which(U_mean_imax30_plot1<15.1&U_mean_imax30_plot1>14.9)



plot(x=0,y=0,ylim=c(0,1),xlim=c(0,1),axes=F,ylab="",xlab="",col=rgb(0,0,0,0))
#legend("left",legend = lapply(par.ranges[[4]],function(j) as.expression(substitute(tau == p,list("p"=as.name(j))))),lty = c(1:length(par.ranges[[4]])))
legend("left", legend = as.expression(lapply(par.ranges[[4]], function(x) bquote(tau==.(x)))), lty = c(1:4,1:4),col=c(rep("black",4),rep("grey50",4)))

plot(0:30,U_pct_imax_plot1[,1],ylim=c(0,1))

points(0:30,U_pct_imax_plot1[,2])
points(0:30,U_pct_imax_plot1[,11])

points(0:30,U_pct_imax_plot1[,31])

points(0:30,U_pct_imax_plot1[,51])
points(0:30,U_pct_imax_plot1[,21])
points(0:30,U_pct_imax_plot1[,41])
png("../Results/20220216/No_int_example.png",width = 4,height = 4,res = 300,units = "in")


dev.off()







