#Potentially, large increase in N plus variable parameter inputs
#will allow direct implementation of ABC?



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

pars_sim<-parseparam0(sampler_fun0(n=1, pars = pars, priorsd = c(0.5, 0.5, 2, 0.5, 0.5, 0.5)))

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
                     lower = c(rep(-6.9,2),-29.9,rep(-6.9,3)), upper = c(rep(2.9,2),29.9,rep(2.9,3)))

#number of MCMC iterations - increase for more accurate results
#note - runtime will be long for EDM example
niter<-10002
nburn<-2001

#with detfun0
#likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0)
likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0)

bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior)
#out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0,
#                       settings = list(iterations=niter, burnin=nburn))
out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0,
                       settings = list(initialParticles=1000, iterations=10), sampler="SMC")

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
#par(mfrow=c(3,2)); for(i in 1:6) {hist(smp_detfun0_untr[,i]); abline(v=unlist(pars_sim)[i]); abline(v=unlist(pars)[i], col=2)}


pfdet<-particleFilterLL(y=datout$obs, pars=parseparam0(colMeans(smp_detfun0_untr)), detfun = detfun0, dotraceback = TRUE)
pfedm<-particleFilterLL(y=datout$obs, pars=parseparam0(colMeans(smp_EDM_untr)), detfun = EDMfun0, edmdat = list(E=2), dotraceback = TRUE)
pftrue<-particleFilterLL(y=datout$obs, pars=pars, detfun = detfun0, dotraceback = TRUE)


#calculate demographic rates
datout_long<-makedynamics(n = 2e4, obs = pars_sim$obs, proc = pars_sim$proc,
                          r = pars_sim$det[1], K = pars_sim$det[2], pcol = pars_sim$pcol)

#Hmm... think of sampling from full distribution of parameter values...

etdfilter_det<-extend_particleFilter(pfout=pfdet, pars=parseparam0(colMeans(smp_detfun0_untr)),
                                     Next = 2e4, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL)
etdfilter_edm<-extend_particleFilter(pfout=pfedm, pars=parseparam0(colMeans(smp_EDM_untr)),
                                     Next = 2e4, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=2))
etdfilter_true<-extend_particleFilter(pfout=pftrue, pars=pars_sim,
                                     Next = 2e4, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL)
#mean(etdfilter_det$demdat$text,na.rm=T)
#mean(etdfilter_edm$demdat$text,na.rm=T)
#mean(etdfilter_true$demdat$text,na.rm=T)
#1/getcm(datout_long$true)$pm
#1/getcm(datout$true)$pm

demdat_true<-getcm(datout_long$true)
demdat_short<-getcm(datout$true)

demdat<-list(demdat_true=demdat_true, demdat_short=demdat_short, demdat_det=etdfilter_det$demdat, demdat_edm=etdfilter_edm$demdat, demdat_filt_true=etdfilter_true$demdat)

#save outputs
save(list = c("out_detfun0", "out_EDM", "pars_sim", "datout", "demdat"), file = paste("datout/mcmcout_", commArgin, "_full.rda", sep=""))



