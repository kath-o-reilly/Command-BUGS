# simulate data from OPV example...
library(R2jags)
library(runjags)      # this one required JAGS installation
library(mcmcplots)

# step 1: generate sample data
prov <- 11  # DRC has 11 provinces
year<- 5    # 5 years of data
N <- 100    # assume 100 samples within each province
tmp <- data.frame(expand.grid(prov=c(1:prov),year=c(1:5)))
tmp$dd <- 1
tmp$ypd <- ((tmp$year-1)*11)+tmp$prov
tmp$sia.d <- rbinom(n=dim(tmp)[1],size=3,p=0.5)+1  # max 4
tmp$sr.d <- ((tmp$prov-1)*4)+tmp$sia.d
tmp$st.d <- ((tmp$year-1)*4)+tmp$sia.d
tmp$ad <- round(rlnorm(dim(tmp)[1],log(12),0.5)) # age
head(tmp)
# parameters are treated as data for the simulation step
#data<-list(N=N,x=x,alpha=alpha,beta=beta,tau=tau)
data <- list(D=dim(tmp)[1],
          R=prov,Y=year,YP=prov*year,S=4,SR=prov*4,ST=year*4,
          ad=tmp$ad,
          rd=tmp$prov,
          dd=tmp$dd,
          td=tmp$year,
          ypd=tmp$ypd,
          sia.d=tmp$sia.d,
          sr.d=tmp$sr.d,
          st.d=tmp$st.d,
          r=seq(-2,2,length.out = prov),   #province
          t=rnorm(year,0,2),   # year of vaccination
          yp=rnorm(prov*year,0,2),   #year:province interaction
          s=rnorm(4,0,2),     #sia count (assume just 4(5 inc baseline) categories)
          sr=rnorm(4*prov,0,2),     #sia count:province interaction
          st=rnorm(4*year,0,2),      #sia count:year interaction
          sigma_r=1,sigma_d=1,sigma_t=1,sigma_yp=1,sigma_s=1,sigma_sr=1,sigma_st=1)

# run jags
out <- run.jags(model="vacc.jag", data = data,
                monitor=c("xd"),
                sample=1, 
                n.chains=1, summarise=FALSE)
# reformat the outputs
Simulated <- coda::as.mcmc(out)
#Simulated
dim(Simulated)
tmp$xd.sim = as.vector(Simulated)

tapply(Simulated/tmp$sia.d,tmp$prov,mean)
tapply(Simulated/tmp$sia.d,tmp$ad,mean)

write.csv(tmp,"cov_sim_data_9March19.csv",row.names = F)

#
