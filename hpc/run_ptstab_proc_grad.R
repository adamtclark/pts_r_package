#!/usr/bin/env Rscript

if(FALSE) {
  error
  rm(list=ls())
  setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
}

## Load functions:
commArgin<-commandArgs(1)
if(length(commArgin)==0) {
  commArgin<-round(runif(1)*1e6)
  commArg_ps<-1
} else {
  commArg_ps<-as.numeric(commArgin)
}

#process error to use
procuse<-seq(1,2,length=20)[commArg_ps]
print(procuse)

require(BayesianTools)
require(rEDM)

source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")


## Simulate data
pars<-list(obs=c(log(1e-2), log(0.1)),
           proc=c(-2, log(procuse)),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(log(3),log(1)))

datout<-makedynamics(n = 100, obs = pars$obs, proc = pars$proc, r = pars$det[1],
                     K = pars$det[2], pcol = pars$pcol)
y<-datout$obs

## Run filter
N = 1e3

## Run optimizer
pars$proc[2]<-log(1.5) #set prior to center of distribution
param<-unlist(pars)[1:6]

#create priors
density_fun_USE<-function(param) density_fun0(param = param, pars = pars)
sampler_fun_USE<-function(x) sampler_fun0(n = 1, pars = pars)
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = rep(-10,6), upper = rep(10,6))

#number of MCMC iterations - increase for more accurate results
#note - runtime will be long for EDM example
niter<-5001
nburn<-1002

#with detfun0
likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0)
bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior)
out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0,
                       settings = list(iterations=niter, burnin=nburn))

#with EDM
likelihood_EDM<-function(x) likelihood0(param = x, y=y, parseparam = parseparam0,
                                        detfun = EDMfun0, edmdat = list(E=2))
bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior)
out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM,
                   settings = list(iterations=niter, burnin=nburn))

#save outputs
save(list = c("out_detfun0", "out_EDM", "procuse"), file = paste("datout/mcmcout_", commArgin, ".rda", sep=""))



