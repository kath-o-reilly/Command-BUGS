library(nimble)
library(coda)
data <- read.csv("data.csv")

N <- nrow(data)
nT <- max(data$t)
n <- max(data$cohort)
Ni <- data$Ni
Nf <- data$Nf
Si <- data$Si
Sf <- data$Sf
Ii <- data$Ii
If <- data$If
C <- data$C
Ri <- data$Ri
In <- data$In
Rn <- data$Rn
cohort <- data$cohort
Rf <- data$Rf
dat <- as.matrix(data[,2:3],rownames.force=NA)

# the model
nbModel <- nimbleCode({ 
  for(j in 1:N){
    # N = total number of rows in the dataset
    #binomial part for passage from susceptibles to infectious
    # C = number of new infectious in each time step for each cohort
    # Si = number of susceptible animals at the beginning of each time step
    # Ni = total number of animals at the beginning of each time step
    # Ii = number of infectious animals at the beginning of each time step
    C[j]~dbin(p[j],Si[j])
    Ioc[j]~dbin(Psen,Ni[j])
    I[j]<- Ii[j]+Ioc[j]
    dumm[j]<-step(0-I[j])
    p[j]<-1-q[j]
    log(q[j])<--B[j]
    # dat is a matrix which includes for each row of the dataset the group number and the time step number (group is the first 
    # column and time is the second column)
    log(B[j])<-int1+log(I[j]+dumm[j])+(-log(Ni[j])*(1-dumm[j]))+alpha[dat[j,1],dat[j,2]]  
    #poisson part for passage from infectious to carrier
    #Rn = number of new carriers
    Rn[j]~dbin(p1[j],I[j])
    p1[j]<-1-q1[j]
    log(q1[j])<--A[j]
    log(A[j])<- delta[dat[j,1]]+int2
    #poisson part for passage from carrier to infectious
    # In = number of new infectious animals that in the time step before were carriers
    In[j]~dpois(phi[j])
    dumm2[j]<-step(0-Ri[j])
    log(phi[j])<-beta[dat[j,1]]+log(Ri[j]+dumm2[j])+int3
  }
  for(m in 1:n){
    # n = total number of groups/cohorts
    alpha[m,1]~dnorm(0,tau5)
    delta[m]~dnorm(0,tau2)
    beta[m]~dnorm(0,tau3)
    for(r in 2:nT){
      # t= total number of time steps
      alpha[m,r]~dnorm(alpha[m,r-1],tau5)
    }
  } 
  senC~dbeta(44.8,36.7)
  senE~dbeta(47.4,31.6)
  Psen <- (1-senC)*(1-senE)
  tau2~dgamma(0.5,0.005)
  tau3~dgamma(0.5,0.005) 
  tau5~dgamma(0.5,0.005)
  int1~dnorm(0,0.01)
  int2~dnorm(0,0.01)
  int3~dnorm(0,0.01)
})

#binomial from s to i and from i to r
nimbleData <- list(Ni=Ni,Si=Si,Ii=Ii,C=C,dat=dat,In=In,Rn=Rn,Ri=Ri)
Consts <- list(N=N,n=n,nT=nT)
inits <- list(tau2=1,tau3=1,tau5=1,int1=-1,int2=0,int3=0)

nimbleNBModel <- nimbleModel(code = nbModel, name = 'nimbleNBModel', 
                             constants = Consts, data = nimbleData, inits = inits)
MCMCconfig <- configureMCMC(nimbleNBModel,monitors=c("tau5","alpha","int1","int2","int3","delta","beta","tau2","tau3","Psen")) 

nbMCMC <- buildMCMC(MCMCconfig)
Cnb <- compileNimble(nimbleNBModel)
CnbMCMC <- compileNimble(nbMCMC, project = nimbleNBModel)

results <- runMCMC(CnbMCMC,niter=300000,nburnin=100000,nchains=3,inits=inits,thin=10,samplesAsCodaMCMC = T)
results <- as.mcmc.list(results)
plot(results,ask=T)
gelman.diag(results)

# int1 is log(beta) in the paper
b <- unlist(results[,"int1"])
# int2 is log(alpha) in the paper
a <- unlist(results[,"int2"])
# int3 is log(nu) in the paper
nu <- unlist(results[,"int3"])
# R_0

# Plots
x11(width=24,height=12)
par(mfrow=c(2,3),mar = c(4, 4, 1, 1),cex=1.3,lwd=3)
plot(density(exp(b)),ylab="density",xlab=expression(beta),main="")
plot(density(exp(a)),ylab="density",xlab=expression(alpha),main="")
plot(density(exp(nu)),ylab="density",xlab=expression(nu),main="")
plot(density(exp(b)/exp(a)),ylab="density",xlab=expression(R[0]),main="")
RanEf_Coh1 <- results[,c("alpha[1, 1]","alpha[1, 2]","alpha[1, 3]","alpha[1, 4]","alpha[1, 5]","alpha[1, 6]","alpha[1, 7]","alpha[1, 8]","alpha[1, 9]")]
RanEf_Coh1 <- rbind(RanEf_Coh1[[1]],RanEf_Coh1[[2]],RanEf_Coh1[[3]])
plot(1:9,apply(RanEf_Coh1,2,mean),xlab="time (bi-weekly steps)",ylab="cohort 1 temporal effect",ylim=c(-5,5),pch=20)
points(1:9,apply(RanEf_Coh1,2,quantile,probs=0.025),pch=4)
points(1:9,apply(RanEf_Coh1,2,quantile,probs=0.975),pch=4)
plot(density(unlist(results[,"Psen"])),ylab="density",xlab="pND",main="")
dev.print(pdf,"figures.pdf")
