#!/usr/bin/env Rscript

if(FALSE) {
  error
  rm(list=ls())
  setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
}

## Load functions:
nclus <- 32
nrep <- 10
commArgin<-commandArgs(1)
if(length(commArgin)==0) {
  commArgin<-round(runif(1)*1e6)
  commArg_ps<-1
} else {
  commArg_ps<-as.numeric(commArgin)
}
print(commArg_ps)

if(length(dir("/cl_tmp/clarka/Rpkg/"))>0) {
  .libPaths("/cl_tmp/clarka/Rpkg/")
}
require(mvtnorm)
require(rEDM)
require(BayesianTools)


commArgin_kp = commArg_ps
tmp = NULL

simlength = 1000
filterlength = 150

for(iclu in 1:nrep) {
  commArgin = commArgin_kp + (iclu - 1)*nclus

  source("../pttstability/R/bayesfun.R")
  source("../pttstability/R/fake_data.R")
  source("../pttstability/R/logit_funs.R")
  source("../pttstability/R/particlefilter.R")

  ## Simulate data
  #sample from priors
  pars0<-pars_true<-list(obs=c(log(0.2)),
                         proc=c(log(0.2), log(1.2)),
                         pcol=c(logit(0.2), log(0.1)),
                         det=c(log(1.2),log(1)))

  #create priors
  p0<-list(c(log(0.01), log(0.5)), c(log(0.01), log(1)))
  minvUSE<-unlist(lapply(p0, function(x) x[1]))
  maxvUSE<-unlist(lapply(p0, function(x) x[2]))

  p0_edm<-list(c(log(0.01), log(0.5)), c(log(0.01), log(1)))
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
  while(sum(y>0)<=(length(y)/3)) { # want at least 5% nonzero values
    pars_sim<-parseparam0(c(unname(sampler_fun_USE())))
    pars_sim$proc=c(pars_sim$proc, pars_sim$det[1])
    #exp(unlist(pars_sim)[1:3])

    datout<-makedynamics_general(n = simlength, n0 = (sin(1/2)+1+0.5)/2, pdet=pars_sim$det,
                                 proc = pars_sim$proc, obs = pars_sim$obs, pcol = pars_sim$pcol,
                                 detfun = detfun0_sin, procfun = procfun_ct, obsfun=obsfun0, colfun=colfun0)
    y<-datout$obs[1:filterlength]
    #plot(y, type="l"); abline(h=0, lty=3)
  }

  ## Run filter
  ptrue<-unname(unlist(pars_sim)[1:2])
  sd_abs<-sdproc_abstract(exp(pars_sim$proc[1]), exp(pars_sim$proc[2]))

  #set number of iterations
  niter<-6e3
  N<-1e3

  sout<-NULL
  for(E in 2:4) {
    sout<-rbind(sout, s_map(y, E=E, silent = TRUE))
  }

  tuse<-sout$theta[which.min(sout$rmse)]
  Euse<-sout$E[which.min(sout$rmse)]
  #plot(sout[sout[,"E"]==Euse,"theta"], sout[sout[,"E"]==Euse,"rho"], type = "b")
  ssave<-s_map(y, E=Euse, theta = tuse, silent = TRUE)


  #set up likelihoods
  likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0, N = N, detfun = detfun0_sin)
  bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior_USE)

  likelihood_EDM<-function(x) {
    likelihood0(param = x, y=y, parseparam = parseparam0,
                detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse, ytot=y), N = N)
  }

  bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior_edm)

  if(FALSE) {
    #likelihoods at "correct" parameters
    tmp<-s_map(y, E = Euse, theta = tuse, silent = TRUE)
    tmp$rho
    likelihood_detfun0(ptrue)
    likelihood_EDM(ptrue)
  }

  #run MCMC chains
  out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0, settings = list(iterations=niter, consoleUpdates=200))
  out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=200))

  ## extract parameters
  smp_detfun0<-getSample(out_detfun0, start = 1000)
  smp_EDM<-getSample(out_EDM, start=1000)

  if(FALSE) {
    #plot results, with a 200-step burn-in
    plot(out_detfun0, start = 1000)
    plot(out_EDM, start = 1000)

    par(mfrow=c(2,2))
    hist(exp(smp_detfun0[,1]), main="det. function", xlab="obs", breaks = 20); abline(v=exp(ptrue[1]), col=2)
    hist(exp(smp_detfun0[,2]), main="det. function", xlab="proc", breaks = 20); abline(v=sd_abs, col=2)

    hist(exp(smp_EDM[,1]), main="EDM function", xlab="obs", breaks = 20); abline(v=exp(ptrue[1]), col=2)
    hist(exp(smp_EDM[,2]), main="EDM function", xlab="proc", breaks = 20); abline(v=sd_abs, col=2)
  }


  parsest_det<-cbind(colMeans(smp_detfun0), apply(smp_detfun0, 2, sd))
  parsest_edm<-cbind(colMeans(smp_EDM), apply(smp_EDM, 2, sd))

  ptrue_abs = c(ptrue[1], log(sd_abs))

  #based on detful0
  filterout_det<-particleFilterLL(y, pars=parseparam0(parsest_det[,1]), detfun = detfun0_sin,
                                  dotraceback = TRUE)
  #based on EDM
  filterout_edm<-particleFilterLL(y, pars=parseparam0(parsest_edm[,1]), detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse, ytot=y),
                                  dotraceback = TRUE)
  #based on true values
  filterout_true<-particleFilterLL(y, pars=parseparam0(ptrue_abs), detfun = detfun0_sin,
                                  dotraceback = TRUE)
  #based on EDM, with correct values
  filterout_edm_true<-particleFilterLL(y, pars=parseparam0(ptrue_abs), detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse, ytot=y),
                                  dotraceback = TRUE)

  ## save outputs
  parslst<-list(ptrue=ptrue,
                parsest_det=parsest_det, parsest_edm=parsest_edm,
                Euse=Euse, tuse=tuse,
                pars_sim=pars_sim,
                sd_abs=sd_abs)

  #fits
  cordat<-list(obs=cor(datout$true[1:filterlength], datout$obs[1:filterlength]),
               det=cor(datout$true[1:filterlength], filterout_det$Nest),
               edm=cor(datout$true[1:filterlength], filterout_edm$Nest),
               true=cor(datout$true[1:filterlength], filterout_true$Nest),
               edm_true=cor(datout$true[1:filterlength], filterout_edm_true$Nest))

  #filter output
  filterdat<-list(filterout_det=filterout_det, filterout_edm=filterout_edm, filterout_true=filterout_true, filterout_edm_true=filterout_edm_true)

  #optimizer outputs
  optdat<-list(optout_det=out_detfun0, optout_edm=out_EDM)

  #simulation outputs
  simdat<-list(datout=datout, ssave=ssave)

  save(list = c("simdat", "parslst", "optdat", "filterdat", "cordat"), file = paste("datout/mcmcout_", commArgin, "_211005.rda", sep=""), version=2)
}
