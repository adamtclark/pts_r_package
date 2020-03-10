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
pars0<-pars_true<-list(obs=log(0.2),
                       proc=c(log(0.1), log(1.2)),
                       pcol=c(logit(0.2), log(0.1)),
                       det=c(log(2),log(1)))

detfun0_sin<-function(sdet, xt, time=NULL) {
  K<-((sin(time/(pi*2))+exp(sdet[2]))*0.5+0.5)
  xt = xt*exp(exp(sdet[1])*(1-xt/K))
  return(xt)
}

#create priors
p0<-list(c(log(0.001), log(0.5)), c(log(0.001), log(0.1)), c(log(0.001), log(3)))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(log(0.001), log(0.5)), c(log(0.001), log(0.1)), c(log(0.001), log(3)))#, c(-5, 2))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

density_fun_USE<-function(param) density_fun0(param = param, minv = minvUSE, maxv=maxvUSE)
sampler_fun_USE<-function(x) sampler_fun0(n = 1, minv = minvUSE, maxv=maxvUSE)
prior_USE <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = minvUSE, upper = maxvUSE)

density_fun_USE_edm<-function(param) density_fun0(param = param, minv = minvUSE_edm, maxv=maxvUSE_edm)
sampler_fun_USE_edm<-function(x) sampler_fun0(n = 1, minv = minvUSE_edm, maxv=maxvUSE_edm)
prior_edm <- createPrior(density = density_fun_USE_edm, sampler = sampler_fun_USE_edm,
                         lower = minvUSE_edm, upper = maxvUSE_edm)
y<-0
while(sum(y>0)<=(length(y)/5)) { # want at least 5% nonzero values
  pars_sim<-parseparam0(unname(sampler_fun_USE()))
  #parseparam0(sampler_fun0(n=1, pars = pars0, priorsd = c(1, 1, 1)))

  datout<-makedynamics_general(n = 100, n0 = exp(rnorm(1,0,0.1)), pdet=pars_sim$det,
                               proc = pars_sim$proc, obs = pars_sim$obs, pcol = pars_sim$pcol,
                               detfun = detfun0_sin, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
  y<-datout$obs
}

## Run filter
ptrue<-unname(unlist(pars_sim)[1:3])

if(FALSE) {
  exp(ptrue)
  plot(y, type="l"); abline(h=0, lty=3)
  sd(datout$true-datout$noproc)
  sd(datout$true-datout$obs)
}

#set number of iterations
niter<-10000
N<-2e3
Euse<-2

#get theta
tmp<-s_map(y, E=Euse, silent = TRUE)
thuse<-tmp$theta[which.max(tmp$rho)]

#set up likelihoods
likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0, N = N, detfun = detfun0_sin)
bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior_USE)

likelihood_EDM<-function(x) {
  xuse<-x#x[1:(length(x)-1)]
  tuse_edm<-thuse#exp(x[length(x)])
  likelihood0(param = xuse, y=y, parseparam = parseparam0,
              detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse_edm), N = N)
}

bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior_edm)

#run MCMC chains
out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0, settings = list(iterations=niter, consoleUpdates=200))
out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=200))

## extract parameters
smp_detfun0<-getSample(out_detfun0, start = 2000)
smp_EDM<-getSample(out_EDM, start=2000)

parsest_det<-cbind(colMeans(smp_detfun0), apply(smp_detfun0, 2, sd))
parsest_edm<-cbind(colMeans(smp_EDM), apply(smp_EDM, 2, sd))

#based on detful0
filterout_det<-particleFilterLL(y, pars=parseparam0(parsest_det[,1]), detfun = detfun0_sin,
                                dotraceback = TRUE)
#based on EDM
filterout_edm<-particleFilterLL(y, pars=parseparam0(parsest_edm[1:length(ptrue),1]), detfun = EDMfun0, edmdat = list(E=Euse, theta=thuse),
                                dotraceback = TRUE)
#based on true values
filterout_true<-particleFilterLL(y, pars=parseparam0(ptrue), detfun = detfun0_sin,
                                dotraceback = TRUE)
#based on EDM, with correct values
sout<-s_map(y, E=Euse,silent = TRUE)
filterout_edm_true<-particleFilterLL(y, pars=parseparam0(ptrue), detfun = EDMfun0, edmdat = list(E=Euse, theta=thuse),
                                dotraceback = TRUE)


## save outputs
parslst<-list(ptrue=ptrue,
              parsest_det=parsest_det, parsest_edm=parsest_edm)

#fits
cordat<-list(obs=cor(datout$true, datout$obs),
             det=cor(datout$true, filterout_det$Nest),
             edm=cor(datout$true, filterout_edm$Nest),
             true=cor(datout$true, filterout_true$Nest),
             edm_true=cor(datout$true, filterout_edm_true$Nest))

#filter output
filterdat<-list(filterout_det=filterout_det, filterout_edm=filterout_edm, filterout_true=filterout_true, filterout_edm_true=filterout_edm_true)

#optimizer outputs
optdat<-list(optout_det=out_detfun0, optout_edm=out_EDM)

#simulation outputs
simdat<-list(datout=datout)

save(list = c("simdat", "parslst", "optdat", "filterdat"), file = paste("datout/mcmcout_", commArgin, "_full_taylor.rda", sep=""), version=2)



