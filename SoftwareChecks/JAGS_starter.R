# JAGS example
# COMMAND lecture
# Author: Kath O'Reilly

# installing relevant libraries
library(R2jags)
library(runjags)      # this one requires JAGS installation
library(mcmcplots)

rm(list = ls())

# load data
dat <- read.csv("dat100_starter.csv")

jags.data <- list(N = dim(dat)[1], 
                  y = dat$y, 
                  x = dat$x,
                  a=0,b=5,
                  m.alpha=0,p.alpha=1,
                  m.beta=0,p.beta=1)

# initial values for the priors
inits <- function(){list(alpha = rnorm(1), 
                         beta = rnorm(1), 
                         log.sigma = runif(1,0,1))}  
# parameters monitored
parameters <- c("alpha", "beta", "tau")

# MCMC settings
ni <- 100000   # number of iterations
nt <- 10       # thinning of mcmc chain
nb <- 5000     # burn-in
nc <- 3        # number of chains

# call JAGS from R

system.time(res <- jags(jags.data, inits, parameters, model.file="regression01.jag", 
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
                        working.directory = getwd()))

# summarize posteriors
print(res, digits = 3)   
