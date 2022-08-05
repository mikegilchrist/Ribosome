#Continuous model
#install.packages("deSolve")
library(deSolve)
#install.packages("bvpSolve")
library(bvpSolve)

# par.names <- c("kappa","mu","lambda","tau","delta","imax")
parameters <- c(1,1,1,1,1,1)
parameters <- c(0.66,0.0292,5.4,0.1,0.1,30)
#Eqtn is second order ode
#eqtn is 0 = a*x + b*dx + c*ddx
#split into two first order ode
#dx=y (then ddx=dy)
#dy = (a*x + b*y)/c


cont <- function(x,y,pars){

    # A <- pars[1]/pars[6]-pars[2]+pars[4]
    # B <- ((x-1)*pars[1])/pars[6]-pars[1]+(x+1)*pars[4]
    # C <- 0.5*(pars[1]+x*pars[4]-(x*pars[1])/pars[6])

    list(c(y[2],
           ((pars[1]/pars[6]-pars[2]+pars[4])*y[1] + (((x-1)*pars[1])/pars[6]-pars[1]+(x+1)*pars[4])*y[2])/(-0.5*(pars[1]+x*pars[4]-(x*pars[1])/pars[6]))))

}
step.size <- 1
times <- seq(0,30,step.size)


unmarked <- bvptwp(yini = c(1,NA),yend = c(0,NA),func = cont,x = times,parms = parameters)
unmarked <- unmarked[which(unmarked[,1]%%(1)==0),]
unmarked[,2] <- unmarked[,2]*parameters[3]/(parameters[2]*sum(unmarked[,2])*step.size)

plot(unmarked[,1],unmarked[,2],type="l")

marked <- cbind(unmarked[,1],rep(NA,nrow(unmarked)))

marked.cont <- function(unmarked,par){
  marked <- cbind(unmarked[,1],rep(NA,nrow(unmarked)))
  for(i in 1:(nrow(unmarked)-1)){
    marked[(i+1),2]<- (par[2]/(par[4]*i))*sum(unmarked[-1*(1:i),2])*step.size
  }
  marked[1,2] <- (par[4]*marked[2,2]+par[2]*unmarked[1,2])/par[5]
  return(marked)
}

marked <- marked.cont(unmarked = unmarked,par = parameters)


plot(unmarked[,1],unmarked[,2],type="l",ylim=c(0,70),xlim=c(0,30))
lines(marked[,1],marked[,2],col="red")
lines(unmarked[,1],unmarked[,2]+marked[,2],col="blue")




# Attempt 2 ---------------------------------------------------------------


# par.names <- c("kappa","mu","lambda","tau","delta","imax")

unmarked.full <- function(t,x,params){

  kappa <- params[1]
  mu <- params[2]
  lambda <- params[3]
  tau <- params[4]
  delta <- params[5]
  imax <- params[6]
  kappa_vec <- rep(kappa,imax+1)*(1-seq(0,imax,1)/imax)
  tau_vec <- rep(tau,imax+1)*seq(0,imax,1)


  dm0dt <- lambda + tau_vec[2]*x[2] -(kappa_vec[1]+mu)*x[1]
  dm1dt <- kappa_vec[1]*x[1] + tau_vec[3]*x[3] -(kappa_vec[2]+mu + tau_vec[2])*x[2]
  dm2dt <- kappa_vec[2]*x[2] + tau_vec[4]*x[4] -(kappa_vec[3]+mu + tau_vec[3])*x[3]
  dm3dt <- kappa_vec[3]*x[3] + tau_vec[5]*x[5] -(kappa_vec[4]+mu + tau_vec[4])*x[4]
  dm4dt <- kappa_vec[4]*x[4] + tau_vec[6]*x[6] -(kappa_vec[5]+mu + tau_vec[5])*x[5]
  dm5dt <- kappa_vec[5]*x[5] + tau_vec[7]*x[7] -(kappa_vec[6]+mu + tau_vec[6])*x[6]
  dm6dt <- kappa_vec[6]*x[6] + tau_vec[8]*x[8] -(kappa_vec[7]+mu + tau_vec[7])*x[7]
  dm7dt <- kappa_vec[7]*x[7] + tau_vec[9]*x[9] -(kappa_vec[8]+mu + tau_vec[8])*x[8]
  dm8dt <- kappa_vec[8]*x[8] + tau_vec[10]*x[10] -(kappa_vec[9]+mu + tau_vec[9])*x[9]
  dm9dt <- kappa_vec[9]*x[9] + tau_vec[11]*x[11] -(kappa_vec[10]+mu + tau_vec[10])*x[10]
  dm10dt <- kappa_vec[10]*x[10] + tau_vec[12]*x[12] -(kappa_vec[11]+mu + tau_vec[11])*x[11]
  dm11dt <- kappa_vec[11]*x[11] + tau_vec[13]*x[13] -(kappa_vec[12]+mu + tau_vec[12])*x[12]
  dm12dt <- kappa_vec[12]*x[12] + tau_vec[14]*x[14] -(kappa_vec[13]+mu + tau_vec[13])*x[13]
  dm13dt <- kappa_vec[13]*x[13] + tau_vec[15]*x[15] -(kappa_vec[14]+mu + tau_vec[14])*x[14]
  dm14dt <- kappa_vec[14]*x[14] + tau_vec[16]*x[16] -(kappa_vec[15]+mu + tau_vec[15])*x[15]
  dm15dt <- kappa_vec[15]*x[15] + tau_vec[17]*x[17] -(kappa_vec[16]+mu + tau_vec[16])*x[16]
  dm16dt <- kappa_vec[16]*x[16] + tau_vec[18]*x[18] -(kappa_vec[17]+mu + tau_vec[17])*x[17]
  dm17dt <- kappa_vec[17]*x[17] + tau_vec[19]*x[19] -(kappa_vec[18]+mu + tau_vec[18])*x[18]
  dm18dt <- kappa_vec[18]*x[18] + tau_vec[20]*x[20] -(kappa_vec[19]+mu + tau_vec[19])*x[19]
  dm19dt <- kappa_vec[19]*x[19] + tau_vec[21]*x[21] -(kappa_vec[20]+mu + tau_vec[20])*x[20]
  dm20dt <- kappa_vec[20]*x[20] + tau_vec[22]*x[22] -(kappa_vec[21]+mu + tau_vec[21])*x[21]
  dm21dt <- kappa_vec[21]*x[21] + tau_vec[23]*x[23] -(kappa_vec[22]+mu + tau_vec[22])*x[22]
  dm22dt <- kappa_vec[22]*x[22] + tau_vec[24]*x[24] -(kappa_vec[23]+mu + tau_vec[23])*x[23]
  dm23dt <- kappa_vec[23]*x[23] + tau_vec[25]*x[25] -(kappa_vec[24]+mu + tau_vec[24])*x[24]
  dm24dt <- kappa_vec[24]*x[24] + tau_vec[26]*x[26] -(kappa_vec[25]+mu + tau_vec[25])*x[25]
  dm25dt <- kappa_vec[25]*x[25] + tau_vec[27]*x[27] -(kappa_vec[26]+mu + tau_vec[26])*x[26]
  dm26dt <- kappa_vec[26]*x[26] + tau_vec[28]*x[28] -(kappa_vec[27]+mu + tau_vec[27])*x[27]
  dm27dt <- kappa_vec[27]*x[27] + tau_vec[29]*x[29] -(kappa_vec[28]+mu + tau_vec[28])*x[28]
  dm28dt <- kappa_vec[28]*x[28] + tau_vec[30]*x[30] -(kappa_vec[29]+mu + tau_vec[29])*x[29]
  dm29dt <- kappa_vec[29]*x[29] + tau_vec[31]*x[31] -(kappa_vec[30]+mu + tau_vec[30])*x[30]
  dm30dt <- kappa_vec[30]*x[30] -(mu + tau_vec[31])*x[31]


  list(c(dm0dt,dm1dt,dm2dt,dm3dt,dm4dt,dm5dt,dm6dt,dm7dt,dm8dt,dm9dt,
         dm10dt,dm11dt,dm12dt,dm13dt,dm14dt,dm15dt,dm16dt,dm17dt,dm18dt,dm19dt,
         dm20dt,dm21dt,dm22dt,dm23dt,dm24dt,dm25dt,dm26dt,dm27dt,dm28dt,dm29dt,dm30dt))

}

xstart <- rep(0,31)
# par.names <- c("kappa","mu","lambda","tau","delta","imax")
parameters <- c(1,1,1,1,1,30)
parameters <- c(0.66,0.0292,5.4,0.1,0.1,30)

res <- ode(func=full.sys, y=xstart, times=seq(from=0,to=500,by=1/20),parms=parameters)

#plot(res)

dim(res)

head(res)
View(res)
5.4/0.0292
max(rowSums(res[,-1]))

plot(res[,1],res[,2],type="l")
plot(rowSums(res[,-1]))
barplot(res[5000,-1])


marked.full <- function(unmarked,xstart,params){
  kappa <- params[1]
  mu <- params[2]
  lambda <- params[3]
  tau <- params[4]
  delta <- params[5]
  imax <- params[6]
  kappa_vec <- rep(kappa,imax+1)*(1-seq(0,imax,1)/imax)
  tau_vec <- rep(tau,imax+1)*seq(0,imax,1)

  marked <- matrix(0,nrow=nrow(unmarked),ncol=ncol(unmarked))
  marked[1,] <- xstart

  for(i in 2:(nrow(unmarked))){
    marked[i,1]<- mu * unmarked[i,1] + tau_vec[2]*marked[(i-1),2] - delta * marked[(i-1),1]
    for(j in 2:imax){
      marked[i,j] <- mu * unmarked[i,j] + tau_vec[j+1]*marked[(i-1),j+1] - tau_vec[j]*marked[(i-1),j]
    }
    marked[i,imax] <- mu * unmarked[i,imax] - tau_vec[imax]*marked[(i-1),imax]
    if(any(marked[i,]<0)){
      marked[i,which(marked[i,]<0)] <- 0
    }
  }
  return(marked)
}

marked <- marked.full(unmarked = res[,-1],xstart = rep(0,31),params = parameters)


head(marked)

plot(res[,1],marked[,1],type="l")
plot(rowSums(marked[,-1]))
barplot(marked[5002,-1])

barplot(colMeans(marked[5000:5500,-1]))

test <- res[5000,-1]+ marked[5000,]

plot(0:30,test,type = "l")
lines(0:30,res[5000,-1],col="red")
lines(0:30,marked[5000,],col="blue")

Unmarked_imax30 <- CalcUnmarkedClassNP(parameters = parameters,imax = 30)
Marked_imax30 <- CalcMarkedClassNP(parameters = parameters,imax = 30)


plot(0:14,test[-c(15:30)],type = "l",lwd=2)
lines(0:14,res[5000,-c(1,15:30)],col="red",lwd=2)
lines(0:14,marked[5000,-c(15:30)],col="blue",lwd=2)
lines(0:14,Unmarked_imax30[-c(15:30)],col="orange",lty=2,lwd=2)
lines(0:14,Marked_imax30[-c(15:30)],col="green",lwd=2,lty=2)

plot(0:14,test[-c(15:30)],type = "l",lwd=2)
lines(0:14,res[5000,-c(1,15:30)],col="red",lwd=2)
lines(0:14,colMeans(marked[5000:5500,-1])[1:15],col="blue",lwd=2)
lines(0:14,Unmarked_imax30[-c(15:30)],col="orange",lty=2,lwd=2)



full.sys <- function(t,x,params){

  kappa <- params[1]
  mu <- params[2]
  lambda <- params[3]
  tau <- params[4]
  delta <- params[5]
  imax <- params[6]
  kappa_vec <- rep(kappa,imax+1)*(1-seq(0,imax,1)/imax)
  tau_vec <- rep(tau/(imax*9),imax+1)*seq(0,imax,1)


  du0dt <- lambda + tau_vec[2]*x[2] -(kappa_vec[1]+mu)*x[1]
  du1dt <- kappa_vec[1]*x[1] + tau_vec[3]*x[3] -(kappa_vec[2]+mu + tau_vec[2])*x[2]
  du2dt <- kappa_vec[2]*x[2] + tau_vec[4]*x[4] -(kappa_vec[3]+mu + tau_vec[3])*x[3]
  du3dt <- kappa_vec[3]*x[3] + tau_vec[5]*x[5] -(kappa_vec[4]+mu + tau_vec[4])*x[4]
  du4dt <- kappa_vec[4]*x[4] + tau_vec[6]*x[6] -(kappa_vec[5]+mu + tau_vec[5])*x[5]
  du5dt <- kappa_vec[5]*x[5] + tau_vec[7]*x[7] -(kappa_vec[6]+mu + tau_vec[6])*x[6]
  du6dt <- kappa_vec[6]*x[6] + tau_vec[8]*x[8] -(kappa_vec[7]+mu + tau_vec[7])*x[7]
  du7dt <- kappa_vec[7]*x[7] + tau_vec[9]*x[9] -(kappa_vec[8]+mu + tau_vec[8])*x[8]
  du8dt <- kappa_vec[8]*x[8] + tau_vec[10]*x[10] -(kappa_vec[9]+mu + tau_vec[9])*x[9]
  du9dt <- kappa_vec[9]*x[9] + tau_vec[11]*x[11] -(kappa_vec[10]+mu + tau_vec[10])*x[10]
  du10dt <- kappa_vec[10]*x[10] + tau_vec[12]*x[12] -(kappa_vec[11]+mu + tau_vec[11])*x[11]
  du11dt <- kappa_vec[11]*x[11] + tau_vec[13]*x[13] -(kappa_vec[12]+mu + tau_vec[12])*x[12]
  du12dt <- kappa_vec[12]*x[12] + tau_vec[14]*x[14] -(kappa_vec[13]+mu + tau_vec[13])*x[13]
  du13dt <- kappa_vec[13]*x[13] + tau_vec[15]*x[15] -(kappa_vec[14]+mu + tau_vec[14])*x[14]
  du14dt <- kappa_vec[14]*x[14] + tau_vec[16]*x[16] -(kappa_vec[15]+mu + tau_vec[15])*x[15]
  du15dt <- kappa_vec[15]*x[15] + tau_vec[17]*x[17] -(kappa_vec[16]+mu + tau_vec[16])*x[16]
  du16dt <- kappa_vec[16]*x[16] + tau_vec[18]*x[18] -(kappa_vec[17]+mu + tau_vec[17])*x[17]
  du17dt <- kappa_vec[17]*x[17] + tau_vec[19]*x[19] -(kappa_vec[18]+mu + tau_vec[18])*x[18]
  du18dt <- kappa_vec[18]*x[18] + tau_vec[20]*x[20] -(kappa_vec[19]+mu + tau_vec[19])*x[19]
  du19dt <- kappa_vec[19]*x[19] + tau_vec[21]*x[21] -(kappa_vec[20]+mu + tau_vec[20])*x[20]
  du20dt <- kappa_vec[20]*x[20] + tau_vec[22]*x[22] -(kappa_vec[21]+mu + tau_vec[21])*x[21]
  du21dt <- kappa_vec[21]*x[21] + tau_vec[23]*x[23] -(kappa_vec[22]+mu + tau_vec[22])*x[22]
  du22dt <- kappa_vec[22]*x[22] + tau_vec[24]*x[24] -(kappa_vec[23]+mu + tau_vec[23])*x[23]
  du23dt <- kappa_vec[23]*x[23] + tau_vec[25]*x[25] -(kappa_vec[24]+mu + tau_vec[24])*x[24]
  du24dt <- kappa_vec[24]*x[24] + tau_vec[26]*x[26] -(kappa_vec[25]+mu + tau_vec[25])*x[25]
  du25dt <- kappa_vec[25]*x[25] + tau_vec[27]*x[27] -(kappa_vec[26]+mu + tau_vec[26])*x[26]
  du26dt <- kappa_vec[26]*x[26] + tau_vec[28]*x[28] -(kappa_vec[27]+mu + tau_vec[27])*x[27]
  du27dt <- kappa_vec[27]*x[27] + tau_vec[29]*x[29] -(kappa_vec[28]+mu + tau_vec[28])*x[28]
  du28dt <- kappa_vec[28]*x[28] + tau_vec[30]*x[30] -(kappa_vec[29]+mu + tau_vec[29])*x[29]
  du29dt <- kappa_vec[29]*x[29] + tau_vec[31]*x[31] -(kappa_vec[30]+mu + tau_vec[30])*x[30]
  du30dt <- kappa_vec[30]*x[30] -(mu + tau_vec[31])*x[31]
  dm0dt <- mu*x[1] + tau_vec[2]*x[33] - delta*x[32]
  dm1dt <- mu*x[2] + tau_vec[3]*x[34] - tau_vec[2]*x[33]
  dm2dt <- mu*x[3] + tau_vec[4]*x[35] - tau_vec[3]*x[34]
  dm3dt <- mu*x[4] + tau_vec[5]*x[36] - tau_vec[4]*x[35]
  dm4dt <- mu*x[5] + tau_vec[6]*x[37] - tau_vec[5]*x[36]
  dm5dt <- mu*x[6] + tau_vec[7]*x[38] - tau_vec[6]*x[37]
  dm6dt <- mu*x[7] + tau_vec[8]*x[39] - tau_vec[7]*x[38]
  dm7dt <- mu*x[8] + tau_vec[9]*x[40] - tau_vec[8]*x[39]
  dm8dt <- mu*x[9] + tau_vec[10]*x[41] - tau_vec[9]*x[40]
  dm9dt <- mu*x[10] + tau_vec[11]*x[42] - tau_vec[10]*x[41]
  dm10dt <- mu*x[11] + tau_vec[12]*x[43] - tau_vec[11]*x[42]
  dm11dt <- mu*x[12] + tau_vec[13]*x[44] - tau_vec[12]*x[43]
  dm12dt <- mu*x[13] + tau_vec[14]*x[45] - tau_vec[13]*x[44]
  dm13dt <- mu*x[14] + tau_vec[15]*x[46] - tau_vec[14]*x[45]
  dm14dt <- mu*x[15] + tau_vec[16]*x[47] - tau_vec[15]*x[46]
  dm15dt <- mu*x[16] + tau_vec[17]*x[48] - tau_vec[16]*x[47]
  dm16dt <- mu*x[17] + tau_vec[18]*x[49] - tau_vec[17]*x[48]
  dm17dt <- mu*x[18] + tau_vec[19]*x[50] - tau_vec[18]*x[49]
  dm18dt <- mu*x[19] + tau_vec[20]*x[51] - tau_vec[19]*x[50]
  dm19dt <- mu*x[20] + tau_vec[21]*x[52] - tau_vec[20]*x[51]
  dm20dt <- mu*x[21] + tau_vec[22]*x[53] - tau_vec[21]*x[52]
  dm21dt <- mu*x[22] + tau_vec[23]*x[54] - tau_vec[22]*x[53]
  dm22dt <- mu*x[23] + tau_vec[24]*x[55] - tau_vec[23]*x[54]
  dm23dt <- mu*x[24] + tau_vec[25]*x[56] - tau_vec[24]*x[55]
  dm24dt <- mu*x[25] + tau_vec[26]*x[57] - tau_vec[25]*x[56]
  dm25dt <- mu*x[26] + tau_vec[27]*x[58] - tau_vec[26]*x[57]
  dm26dt <- mu*x[27] + tau_vec[28]*x[59] - tau_vec[27]*x[58]
  dm27dt <- mu*x[28] + tau_vec[29]*x[60] - tau_vec[28]*x[59]
  dm28dt <- mu*x[29] + tau_vec[30]*x[61] - tau_vec[29]*x[60]
  dm29dt <- mu*x[30] + tau_vec[31]*x[62] - tau_vec[30]*x[61]
  dm30dt <- mu*x[31] - tau_vec[31]*x[62]

  list(c(du0dt,du1dt,du2dt,du3dt,du4dt,du5dt,du6dt,du7dt,du8dt,du9dt,
         du10dt,du11dt,du12dt,du13dt,du14dt,du15dt,du16dt,du17dt,du18dt,du19dt,
         du20dt,du21dt,du22dt,du23dt,du24dt,du25dt,du26dt,du27dt,du28dt,du29dt,du30dt,
         dm0dt,dm1dt,dm2dt,dm3dt,dm4dt,dm5dt,dm6dt,dm7dt,dm8dt,dm9dt,
         dm10dt,dm11dt,dm12dt,dm13dt,dm14dt,dm15dt,dm16dt,dm17dt,dm18dt,dm19dt,
         dm20dt,dm21dt,dm22dt,dm23dt,dm24dt,dm25dt,dm26dt,dm27dt,dm28dt,dm29dt,dm30dt))

}
xstart <- rep(0,62)

parameters <- c(1/17,0.0001,1,3.5,100000,30)
res <- ode(func=full.sys, y=xstart, times=seq(from=0,to=200000,by=10),parms=parameters)

Unmarked_imax30_full <- res[,1:32]
Marked_imax30_full <- res[,c(1,33:63)]

test <- rowSums(res[,-1])
plot(Marked_imax30_full[5000,-1])

plot(0:14,Unmarked_imax30_full[17000,c(2:16)]+Marked_imax30_full[10000,c(2:16)],type = "l",lwd=2)
lines(0:14,Unmarked_imax30_full[17000,2:16],col="red",lwd=2)
lines(0:14,Marked_imax30_full[17000,2:16],col="blue",lwd=2)
lines(0:14,Unmarked_imax30[-c(15:30)],col="green",lty=2,lwd=2)
lines(0:14,Marked_imax30[1:15],col="orange",lwd=2,lty=2)
lines(0:14,newtau_test[1:15],col="green",lty=2,lwd=2)
newtau_test <- CalcUnmarkedClass_newtau(parameters = c(1/17,0.0001,1,3.5,1000000,1),imax = 30,rescale = F,interference = T)

Unmarked_imax30_full[5000,-1]


s2 <- function(params,i){
  kappa <- params[1]
  mu <- params[2]
  lambda <- params[3]
  tau <- params[4]
  delta <- params[5]
  imax <- params[6]

  return(0.5*((kappa+i/imax)/(((kappa+tau)/imax)-mu)))

}

library(truncnorm)

which(unlist(Unmarked_imax30_full[5000,2:16])%in%max(Unmarked_imax30_full[5000,2:16]))
plot(0:30,Unmarked_imax30_full[5000,-1],col="red",lwd=2,type="l")
lines(0:30,trnorm*sum(Unmarked_imax30_full[5000,-1]),lty=2)

dnorm(0:30,mean=5,sd = sqrt(s2(parameters,5)))


i_mu <- sum((Unmarked_imax30_full[5000,-1]/sum(Unmarked_imax30_full[5000,-1]))*(0:30))
i_var <- sum((Unmarked_imax30_full[5000,-1]/sum(Unmarked_imax30_full[5000,-1]))*(0:30-i_mu)^2)

trnorm <- truncnorm::dtruncnorm(0:30,mean=i_mu,sd = sqrt(i_var),a = 0,b = 30)


mean(Unmarked_imax30_full[5000,-1]%*%t(0:31))

# Half life ---------------------------------------------------------------

xstart <- res[5000,-1]

res <- ode(func=full.sys, y=xstart, times=seq(from=0,to=500,by=1/20),parms=c( 0.6600,0.0292,0,0.1000,0.1000,30.0000))

Unmarked_imax30_full <- res[,1:32]
Marked_imax30_full <- res[,c(1,33:63)]

plot(Unmarked_imax30_full[,1],rowSums(Unmarked_imax30_full[,-1]),pch=20)
lines(Unmarked_imax30_full[,1],sum(xstart[1:31])*exp(-0.0292*Unmarked_imax30_full[,1]),col="red")

plot(Marked_imax30_full[,1],rowSums(Marked_imax30_full[,-c(1:2)]),pch=20)
plot(Marked_imax30_full[,1],Marked_imax30_full[,2],pch=20)

xstart_hl <- rep(0,62)
xstart_hl[33] <- 100

res <- ode(func=full.sys, y=xstart_hl, times=seq(from=0,to=125,by=0.01),parms=c( 0.6600,0.0292,0,0.1000,0.1000,30.0000))


Unmarked_imax30_full <- res[,1:32]
Marked_imax30_full <- res[,c(1,33:63)]


plot(Unmarked_imax30_full)

plot(Marked_imax30_full[,1],rowSums(Marked_imax30_full[,-c(1:2)]),pch=20)
lines(Marked_imax30_full[,1],sum(xstart[1:31])*exp(-0.1*1.35*Marked_imax30_full[,1]),col="red")
plot(Marked_imax30_full[,1],Marked_imax30_full[,2],pch=20)

plot(Marked_imax30_full[,1],rowSums(Marked_imax30_full[,-c(1:2)]),type="l")
which(rowSums(Marked_imax30_full[,-c(1,2)])<50.1&rowSums(Marked_imax30_full[,-c(1,2)])>49.9)
rowSums(Marked_imax30_full[1679:1680,-1])
Marked_imax30_full[694:695,1]
lines(Marked_imax30_full[,1],sum(xstart[1:31])*exp(-1*(log(2)/6.94)*Marked_imax30_full[,1]),col="red",lty=2)
lines(Marked_imax30_full[,1],sum(xstart[1:31])*exp(-1*0.1*Marked_imax30_full[,1]),col="blue",lty=2)
abline(v=6.94,lty=2)
abline(h=50,lty=2)

y<-lapply(1:nrow(Marked_imax30_full),function(i) sum(Marked_imax30_full[i,-c(1,2)])*exp(-1*0.1*Marked_imax30_full[i,1]))
lines(Marked_imax30_full[,1],y,col="blue",lty=2)


plot(Marked_imax30_full[,1],Marked_imax30_full[,3])
lines(Marked_imax30_full[,1],Marked_imax30_full[,2],col="red")
abline(h=(1+log(0.5))*100)

plot(Marked_imax30_full[,1],log(rowSums(Marked_imax30_full[,-c(1:2)]),10),type="l")

plot(res[,1],rowSums(res[,-1]))
plot(0:50,10^8*(exp(-1*(log(2)/120)*0:50)),ylim=c(10^6,10^9))
10^8*(exp(-1*(log(2)/120)*0:50))











