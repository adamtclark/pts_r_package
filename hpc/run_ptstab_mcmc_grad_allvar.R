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

#require(BayesianTools)
require(mvtnorm)
require(rEDM)
require(BayesianTools)

source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")

## Simulate data
#sample from priors
pars0<-list(obs=(-2.5),
                proc=(-2.5),
                pcol=c((-1.5), NA),
                det=c(log(3),log(1)))
pars0$pcol[2]<-pars0$proc

search_rng<-3
prior_sd<-1

#create priors
density_fun_USE<-function(param) density_fun0(param = param, pars = pars0, priorsd = rep(prior_sd, 3))
sampler_fun_USE<-function(x) sampler_fun0(n = 1, pars = pars0, priorsd = rep(prior_sd, 3),
                                          minv = c(pars0$obs, pars0$proc, pars0$pcol[1])-search_rng,
                                          maxv = c(pars0$obs, pars0$proc, pars0$pcol[1])+search_rng)
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = c(pars0$obs, pars0$proc, pars0$pcol[1])-search_rng,
                     upper = c(pars0$obs, pars0$proc, pars0$pcol[1])+search_rng)

y<-0
while(sum(y>0)<=(length(y)/20)) { # want at least 5% nonzero values
  pars_sim<-parseparam0(unname(sampler_fun_USE()))
  #parseparam0(sampler_fun0(n=1, pars = pars0, priorsd = c(1, 1, 1)))

  datout<-makedynamics_general(n = 100, n0 = exp(rnorm(1,0,0.1)), pdet=pars_sim$det,
                               proc = pars_sim$proc, obs = pars_sim$obs, pcol = pars_sim$pcol,
                               detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
  y<-datout$obs
}

## Run filter
p0<-unname(unlist(pars0)[1:3])
ptrue<-unname(unlist(pars_sim)[1:3])

#set number of iterations
niter<-6e3
N<-2e3
neff_use<-10

#set up likelihoods
likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0, N = N, neff = neff_use)
bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior)

likelihood_EDM<-function(x) likelihood0(param = x, y=y, parseparam = parseparam0,
                                        detfun = EDMfun0, edmdat = list(E=2), N = N, neff = neff_use)
bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior)

#run MCMC chains
out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0, settings = list(iterations=niter, consoleUpdates=10))
out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=10))

## extract parameters
smp_detfun0<-getSample(out_detfun0, start = 500)
smp_EDM<-getSample(out_EDM, start=500)

parsest_det<-cbind(colMeans(smp_detfun0), apply(smp_detfun0, 2, sd))
parsest_edm<-cbind(colMeans(smp_EDM), apply(smp_EDM, 2, sd))

#based on detful0
filterout_det<-particleFilterLL(y, pars=parseparam0(parsest_det[,1]), detfun = detfun0,
                                dotraceback = TRUE)
#based on EDM
filterout_edm<-particleFilterLL(y, pars=parseparam0(parsest_edm[,1]), detfun = EDMfun0, edmdat = list(E=2),
                                dotraceback = TRUE)
#based on true values
filterout_true<-particleFilterLL(y, pars=parseparam0(ptrue), detfun = detfun0,
                                dotraceback = TRUE)
#based on EDM, with correct values
filterout_edm_true<-particleFilterLL(y, pars=parseparam0(ptrue), detfun = EDMfun0, edmdat = list(E=2),
                                dotraceback = TRUE)


## save outputs
#from simulated time series
cm_true<-getcm(datout$true)
cm_obs<-getcm(datout$obs)

#from extended time series
datout_long<-makedynamics_general(n = 2e4, n0 = exp(rnorm(1,0,0.1)), pdet=pars_sim$det,
                                  proc = pars_sim$proc, obs = pars_sim$obs, pcol = pars_sim$pcol,
                                  detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
cm_long<-getcm(datout_long$true)

#demographics
demdat<-list(cm_true=unlist(cm_true), cm_obs=unlist(cm_obs), cm_long=unlist(cm_long),
             demdat_det=unlist(filterout_det$dem[c("mucol", "mumor")]), demdat_edm=unlist(filterout_edm$dem[c("mucol", "mumor")]), demdat_true=unlist(filterout_true$dem[c("mucol", "mumor")]), demdat_edm_true=unlist(filterout_edm_true$dem[c("mucol", "mumor")]))

#paramters
pars_sim_realized<-log(c(sd(datout$obs-datout$true), sd(datout$noproc-datout$true)))
pars_sim_realized_long<-log(c(sd(datout_long$obs-datout_long$true), sd(datout_long$noproc-datout_long$true)))

parslst<-list(pars0=pars0, pars_sim=pars_sim, p0=p0, ptrue=ptrue,
              pars_sim_realized=pars_sim_realized, pars_sim_realized_long=pars_sim_realized_long,
              parsest_det=parsest_det, parsest_edm=parsest_edm)

#fits
cordat<-list(obs=cor(datout$true, datout$obs),
             det=cor(datout$true, filterout_det$rN),
             edm=cor(datout$true, filterout_edm$rN),
             true=cor(datout$true, filterout_true$rN))

#filter output
filterdat<-list(filterout_det=filterout_det, filterout_edm=filterout_edm, filterout_true=filterout_true, filterout_edm_true=filterout_edm_true)

#optimizer outputs
optdat<-list(optout_det=out_detfun0, optout_edm=out_EDM)

#simulation outputs
simdat<-list(datout=datout, datout_long=datout_long)

save(list = c("simdat", "parslst", "optdat", "filterdat", "demdat"), file = paste("datout/mcmcout_", commArgin, "_full.rda", sep=""), version=2)



