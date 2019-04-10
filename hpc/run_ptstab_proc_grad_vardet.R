#!/usr/bin/env Rscript

if(FALSE) {
  error
  rm(list=ls())
  setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
}

########################################
## Load functions:
########################################

commArgin<-commandArgs(1)
if(length(commArgin)==0) {
  commArgin<-round(runif(1)*1e6)
  commArg_ps<-1
} else {
  commArg_ps<-as.numeric(commArgin)
}
print(commArg_ps)

require(BayesianTools)
require(rEDM)

source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")


########################################
#NEW FUNCTIONS
########################################

#NEW DETERMINISTIC FUNCTION
detfun_new<-function(sdet, xt) {
  xt = xt*exp(exp(sdet[1])*(1-xt/exp(sdet[2])))
  return(xt)
}

#NEW PROCESS NOISE FUNCTION
procfun_new<-function(sp, xt) {
  sm<-length(xt[xt>0])

  #process noise, following classic Ricker function
  #N(t+1) = N(t)*exp(r*(1-N(t)/K)+rnorm(1, 0, exp(sp[1])))
  xt[xt>0] = exp(log(xt[xt>0])+rnorm(sm, 0, exp(sp[1])))

  #add in mortality following random binomial distribution
  #with pmor = ilogit(sp[2]+sp[3]*xt)
  xt[xt>0] = rbinom(length(xt>0), 1, 1-ilogit(sp[2]+sp[3]*xt[xt>0]))*xt[xt>0]

  return(xt)
}

#NEW PARSING FUNCTIONS
parseparam_NEW<-function(param, detparam=c(log(3),log(1))) {
  if(length(param)==7) {
    pars<-list(obs=c(param[1], param[2]),
               proc=c(param[3], param[4], param[5]),
               pcol=c(param[6], param[7]),
               det=detparam)
  } else if(length(param)==9) {
    pars<-list(obs=c(param[1], param[2]),
               proc=c(param[3], param[4], param[5]),
               pcol=c(param[6], param[7]),
               det=c(param[8], param[9]))
  } else {
    return("error: param must be either length 7 or 9")
  }

  return(pars)
}

sampler_fun_NEW<-function(n=1, pars=pars, priorsd=c(0.5, 0.5, 0.5, 0.5, 2, 0.5, 0.5)){
  d1 = rnorm(n, mean = lognormal_imode(pars$obs[1], priorsd[1]), sd = priorsd[1])
  d2 = rnorm(n, mean = lognormal_imode(pars$obs[2], priorsd[2]), sd = priorsd[2])
  d3 = rnorm(n, mean = lognormal_imode(pars$proc[1], priorsd[3]), sd = priorsd[3])
  d4 = rnorm(n, mean = logitnormal_imode(pars$proc[2], priorsd[4]), sd = priorsd[4])
  d5 = rnorm(n, mean = pars$proc[3], sd=priorsd[5])
  d6 = rnorm(n, mean = logitnormal_imode(pars$pcol[1], priorsd[6]), sd = priorsd[6])
  d7 = rnorm(n, mean = lognormal_imode(pars$pcol[2], priorsd[7]), sd=priorsd[7])

  if(length(priorsd)==9) {
    d8 = rnorm(n, mean = lognormal_imode(pars$det[1], priorsd[8]), sd=priorsd[8])
    d9 = rnorm(n, mean = lognormal_imode(pars$det[2], priorsd[9]), sd=priorsd[9])
    return(cbind(d1,d2,d3,d4,d5,d6,d7,d8,d9))
  } else {
    return(cbind(d1,d2,d3,d4,d5,d6,d7))
  }
}

inv_fun_NEW<-function(x) {
  if(ncol(x)==7) {
    cbind(exp(x[,1]), exp(x[,2]), exp(x[,3]), ilogit(x[,4]), x[,5], ilogit(x[,6]), exp(x[,7]))
  } else {
    cbind(exp(x[,1]), exp(x[,2]), exp(x[,3]), ilogit(x[,4]), x[,5], ilogit(x[,6]), exp(x[,7]), exp(x[,8]), exp(x[,9]))
  }
}

density_fun_NEW<-function(param, pars=pars, priorsd=c(0.5, 0.5, 0.5, 0.5, 2, 0.5, 0.5)){
  dsum = dnorm(param[1], mean = lognormal_imode(pars$obs[1], priorsd[1]), sd =  priorsd[1], log = TRUE)
  dsum = dsum+dnorm(param[2], mean = lognormal_imode(pars$obs[2], priorsd[2]), sd =  priorsd[2], log = TRUE)

  dsum = dsum+dnorm(param[3], mean = lognormal_imode(pars$proc[1], priorsd[3]), sd = priorsd[3], log = TRUE)
  dsum = dsum+dnorm(param[4], mean = logitnormal_imode(pars$proc[2], priorsd[4]), sd = priorsd[4], log = TRUE)
  dsum = dsum+dnorm(param[5], mean = pars$proc[3], sd = priorsd[5], log = TRUE)

  dsum = dsum+dnorm(param[6], mean = logitnormal_imode(pars$pcol[1], priorsd[6]), sd = priorsd[6], log = TRUE)
  dsum = dsum+dnorm(param[7], mean = lognormal_imode(pars$pcol[2], priorsd[7]), sd = priorsd[7], log = TRUE)

  if(length(param)==9) {
    dsum = dsum+dnorm(param[8], mean = lognormal_imode(pars$det[1], priorsd[8]), sd = priorsd[8], log = TRUE)
    dsum = dsum+dnorm(param[9], mean = lognormal_imode(pars$det[2], priorsd[9]), sd = priorsd[9], log = TRUE)
  }

  return(dsum)
}


########################################
## Simulate data
########################################

#sample from priors
pars<-list(obs=c(log(1e-2), log(0.1)),
           proc=c(log(0.1), logit(0.1), -1),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(log(3),log(1)))

pars_sim<-parseparam_NEW(sampler_fun_NEW(n=1, pars = pars))

datout<-makedynamics_general(n=100, n0=0.1,
                     pdet=pars_sim$det, proc=pars_sim$proc,
                     obs=pars_sim$obs, pcol=pars_sim$pcol,
                     detfun = detfun_new, procfun = procfun_new, obsfun=obsfun0, colfun=colfun0,
                     doplot=FALSE)
#plot(datout$true, type="l")
y<-datout$obs


########################################
## Run filter
########################################

N = 1e3

## Run optimizer
param<-unlist(pars)[1:7]

#create priors
density_fun_USE<-function(param) density_fun_NEW(param = param, pars = pars)
sampler_fun_USE<-function(x) sampler_fun_NEW(n = 1, pars = pars)
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = rep(-10,7), upper = rep(10,7))

#number of MCMC iterations - increase for more accurate results
#note - runtime will be long for EDM example
niter<-5001
nburn<-1002

#with detfun0
likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam_NEW, detfun = detfun_new)
bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior)
out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0,
                       settings = list(iterations=niter, burnin=nburn))

#with EDM
likelihood_EDM<-function(x) likelihood0(param = x, y=y, parseparam = parseparam_NEW,
                                        detfun = EDMfun0, edmdat = list(E=2))
bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior)
out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM,
                   settings = list(iterations=niter, burnin=nburn))

#run PF's
smp_detfun0_untr<-(getSample(out_detfun0))
smp_EDM_untr<-(getSample(out_EDM))

pfdet<-particleFilterLL(y=datout$obs, pars=parseparam_NEW(colMeans(smp_detfun0_untr)), detfun = detfun_new, dotraceback = TRUE)
pfedm<-particleFilterLL(y=datout$obs, pars=parseparam_NEW(colMeans(smp_EDM_untr)), detfun = EDMfun0, edmdat = list(E=2), dotraceback = TRUE)
pftrue<-particleFilterLL(y=datout$obs, pars=pars, detfun = detfun_new, dotraceback = TRUE)

truecol<-sum(datout$true[-1]>0 & datout$true[-length(datout$true)]==0)/sum(datout$true[-length(datout$true)]==0)
truemor<-sum(datout$true[-1]==0 & datout$true[-length(datout$true)]>0)/sum(datout$true[-length(datout$true)]>0)

demdat<-list(pfdet=pfdet, pfedm=pfedm, pftrue=pftrue, truecol=truecol, truemor=truemor)

#save outputs
save(list = c("out_detfun0", "out_EDM", "pars_sim", "datout", "demdat"), file = paste("datout/mcmcout_", commArgin, "_varfxn.rda", sep=""))



