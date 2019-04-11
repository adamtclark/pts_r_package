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
print(commArg_ps)

require(BayesianTools)
require(rEDM)

source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")


## Simulate data
#sample from priors
pars<-list(obs=c(log(1e-2), log(0.1)),
           proc=c(-2, log(1.5)),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(log(3),log(1)))

pars_sim<-parseparam0(sampler_fun0(n=1, pars = pars))


datout<-makedynamics(n = 100, obs = pars_sim$obs, proc = pars_sim$proc, r = pars_sim$det[1],
                     K = pars_sim$det[2], pcol = pars_sim$pcol)
y<-datout$obs

## Run filter
N = 1e3

## Run optimizer
param<-unlist(pars)[1:6]

#create priors
density_fun_USE<-function(param) density_fun0(param = param, pars = pars)
sampler_fun_USE<-function(x) sampler_fun0(n = 1, pars = pars)
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = rep(-10,6), upper = rep(10,6))

#number of MCMC iterations - increase for more accurate results
#note - runtime will be long for EDM example
niter<-10002
nburn<-2001

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

#run PF's
smp_detfun0_untr<-(getSample(out_detfun0))
smp_EDM_untr<-(getSample(out_EDM))
#smp_detfun0_untr<-smp_EDM_untr<-t(unlist(pars_sim)[1:6])

pfdet<-particleFilterLL(y=datout$obs, pars=parseparam0(colMeans(smp_detfun0_untr)), detfun = detfun0, dotraceback = TRUE)
pfedm<-particleFilterLL(y=datout$obs, pars=parseparam0(colMeans(smp_EDM_untr)), detfun = EDMfun0, edmdat = list(E=2), dotraceback = TRUE)
pftrue<-particleFilterLL(y=datout$obs, pars=pars, detfun = detfun0, dotraceback = TRUE)


#calculate demographic rates
datout_long<-makedynamics(n = 2e4, obs = pars_sim$obs, proc = pars_sim$proc,
                          r = pars_sim$det[1], K = pars_sim$det[2], pcol = pars_sim$pcol)

etdfilter_det<-extend_particleFilter(pfout=pfdet, pars=parseparam0(colMeans(smp_detfun0_untr)),
                                     Next = 1e3, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL)
etdfilter_edm<-extend_particleFilter(pfout=pfedm, pars=parseparam0(colMeans(smp_EDM_untr)),
                                     Next = 1e3, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=2))
etdfilter_true<-extend_particleFilter(pfout=pftrue, pars=pars_sim,
                                     Next = 1e3, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL)
#mean(1/etdfilter_det$demdat$text,na.rm=T)
#mean(1/etdfilter_edm$demdat$text,na.rm=T)
#mean(1/etdfilter_true$demdat$text,na.rm=T)
#getcm(datout_long$true)$pm

demdat_true<-getcm(datout_long$true)
demdat_short<-getcm(datout$true)

demdat<-list(demdat_true=demdat_true, demdat_short=demdat_short, demdat_det=etdfilter_det$demdat, demdat_edm=etdfilter_edm$demdat)

#save outputs
save(list = c("out_detfun0", "out_EDM", "pars_sim", "datout", "demdat"), file = paste("datout/mcmcout_", commArgin, "_full.rda", sep=""))



