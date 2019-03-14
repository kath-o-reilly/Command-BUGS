# Estimate parameters from the AFP and DHA data
rm(list = ls())

library(R2jags)
library(runjags)      # this one required JAGS installation
library(mcmcplots)
# load in data
dd <- read.csv("covDHS_sim_data_9March19.csv")
da <- read.csv("covAFP_sim_data_9March19.csv")
prov <- 11
year <- 5

parameters <- c("province", "year", "yp","d1.OR", "d.avg.OR", "sia", "sr", "st",
                "sigma_r", "sigma_d", "sigma_t", "sigma_yp", "sigma_s", "sigma_sr", "sigma_st")

inits <- function(){list(r = rnorm(prov, 0, 0.1), #province
                         t = rnorm(year, 0, 0.1), #year of vacc
                         yp = rnorm(prov*year, 0, 0.1), #year:province interaction
                         s = rnorm(4, 0, 0.1), #sia count 
                         sr = rnorm(4*prov, 0, 0.1), #sia count:province interaction
                         st = rnorm(4*year, 0, 0.1))} #sia count:year interaction

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

#MODEL DATA
jags.data <- list(D = dim(dd)[1],
                  ad = dd$ad,
                  xd = dd$xd.sim,
                  dd = rep(1,dim(dd)[1]),
                  rd = dd$prov,
                  td = dd$year,
                  ypd = dd$ypa,
                  sia.d = dd$sia.a, #need to shift SIA terms so there is no zero node (sia=0 --> s[1])
                  sr.d = dd$sr.a,
                  st.d = dd$st.a,
                  A = dim(da)[1],
                  xa = da$xd.sim,
                  aa = da$ad,
                  ra = da$prov,
                  ta = da$year,
                  ypa = da$ypd,
                  sia.a = da$sia.d,
                  sr.a = da$sr.d,
                  st.a = da$st.d,
                  R = prov,
                  Y = year,
                  YP = prov*year,
                  S = 4,
                  SR = 4*prov,
                  ST = 4*year)

parameters <- c("province", "year")

out <- jags(jags.data, inits, parameters, model.file="vacc_est.jag", 
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
                    working.directory = getwd())

print(out, digits = 3)


# end
