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

source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")


## Simulate data
#sample from priors
pars0<-list(obs=log(0.2),
                proc=log(0.1),
                pcol=c(logit(0.2), log(0.1)),
                det=c(log(3),log(1)))

y<-0
while(sum(y>0)<=(length(y)/5)) { # want at least 20% nonzero values
  pars_sim<-parseparam0(sampler_fun0(n=1, pars = pars0, priorsd = c(2, 2, 2)))

  datout<-makedynamics_general(n = 100, n0 = exp(rnorm(1,0,0.1)), pdet=pars_sim$det,
                               proc = pars_sim$proc, obs = pars_sim$obs, pcol = pars_sim$pcol,
                               detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
  y<-datout$obs
}

## Run filter
likelihoodEDM<-function(param, y, parseparam, N) {likelihood0(param, y, parseparam, N, detfun = EDMfun0, edmdat = list(E=2))}
p0<-unname(unlist(pars0)[1:3])
ptrue<-unname(unlist(pars_sim)[1:3])

#run optimizer
optout_det<-run_ABC_optim(y, sd0 = c(2,2,2), p0, likelihood = likelihood0, fretain = 0.5, niter_optim = 100, silent = TRUE)
optout_edm<-run_ABC_optim(y, sd0 = c(2,2,2), p0, likelihood = likelihoodEDM, fretain = 0.5, niter_optim = 100, silent = TRUE)

#process outputs
dens_out_det<-abc_densities(optout = optout_det, param0 = p0, param_true = ptrue, fretain = 0.75, enp.target = 4, nbootstrap = 100, nobs = length(y), doplot = FALSE)
dens_out_edm<-abc_densities(optout = optout_edm, param0 = p0, param_true = ptrue, fretain = 0.75, enp.target = 4, nbootstrap = 100, nobs = length(y), doplot = FALSE)

## calculate demographic rates
parsest_det<-cbind(dens_out_det$muest, dens_out_det$sdest)
parsest_edm<-cbind(dens_out_edm$muest, dens_out_edm$sdest)

#based on detful0
filterout_det<-particleFilterLL(y, pars=parseparam0(dens_out_det$muest), detfun = detfun0,
                                dotraceback = TRUE)
#based on EDM
filterout_edm<-particleFilterLL(y, pars=parseparam0(dens_out_edm$muest), detfun = EDMfun0, edmdat = list(E=2),
                                dotraceback = TRUE)
#based on true values
filterout_true<-particleFilterLL(y, pars=parseparam0(ptrue), detfun = detfun0,
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
demdat<-list(cm_true=cm_true, cm_obs=cm_obs, cm_long=cm_long,
             demdat_det=filterout_det$dem[c("mucol", "mumor")], demdat_edm=filterout_edm$dem[c("mucol", "mumor")], demdat_true=filterout_true$dem[c("mucol", "mumor")])

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
filterdat<-list(filterout_det=filterout_det, filterout_edm=filterout_edm, filterout_true=filterout_true)

#optimizer outputs
optdat<-list(optout_det=optout_det, optout_edm=optout_edm)

#simulation outputs
simdat<-list(datout=datout, datout_long=datout_long)

#density function outputs
densout<-list(dens_out_det=dens_out_det, dens_out_edm=dens_out_edm)

save(list = c("simdat", "parslst", "optdat", "densout", "filterdat", "demdat"), file = paste("datout/abcout_", commArgin, "_full.rda", sep=""), version=2)



