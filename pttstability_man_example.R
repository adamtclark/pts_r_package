#!/usr/bin/env Rscript

#error
#rm(list=ls())

## Load functions:
#setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/pttstability/")
require(BayesianTools)
require(mvtnorm)
require(rEDM)
source("R/bayesfun.R")
source("R/fake_data.R")
source("R/logit_funs.R")
source("R/particlefilter.R")

#Set seed
set.seed(2110)

## Simulate data
pars_true<-list(obs=log(0.2),
            proc=c(log(0.1), log(0.9)),
            pcol=c(logit(0.2), log(0.1)),
            det=c(log(2),log(0.2)))

#generate random parameter values
datout<-makedynamics_general(n = 100, n0 = exp(pars_true$det[2]), pdet=pars_true$det,
                             proc = pars_true$proc, obs = pars_true$obs, pcol = pars_true$pcol,
                             detfun = detfun0, procfun = procfun_ct, obsfun=obsfun0, colfun=colfun0)
y<-datout$obs
plot(y, type = "l", xlab="time", ylab="observed abundance")

#get theta paramter for s-mapping
s<-s_map(y, E=2, silent = TRUE)
tuse<-s$theta[which.min(s$rmse)]
plot(s$theta, s$rho, type="b")

## Run filter
N = 1e3
#based on detful0
filterout_det<-particleFilterLL(y, pars=pars_true, N, detfun = detfun0, procfun = procfun_ct,
                                dotraceback = TRUE, fulltraceback = TRUE)
#based on EDM
filterout_edm<-particleFilterLL(y, pars=pars_true, N, detfun = EDMfun0, edmdat = list(E=2, theta=tuse), procfun = procfun_ct,
                                dotraceback = TRUE, fulltraceback = TRUE)

#plot filter output
par(mar=c(4,4,2,2), mfrow=c(3,1))
#plot 30 of the 1000 particles to show general trend
matplot(1:length(y), filterout_det$fulltracemat[,1:30], col=adjustcolor(1,alpha.f = 0.5), lty=3,
        type="l", xlab="time", ylab="abund", main="detfun0")
lines(1:length(y), y, col=2, lwd=1.5) #observations
lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2) #true values
lines(1:length(y), filterout_det$Nest, col=3, lwd=1.5) #mean estimate from filter

matplot(1:length(y), filterout_edm$fulltracemat[,1:30], col=adjustcolor(1,alpha.f = 0.5), lty=3,
        type="l", xlab="time", ylab="abund", main="EDM")
lines(1:length(y), y, col=2, lwd=1.5)
lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2)
lines(1:length(y), filterout_edm$Nest, col=3, lwd=1.5)

plot(filterout_det$Nest, datout$true, xlim=range(c(filterout_det$Nest, filterout_edm$Nest)),
     xlab="predicted", ylab="true value", col=4)
points(filterout_edm$Nest, datout$true, col=2)
points(y, datout$true, col=3)
abline(a=0, b=1, lty=2)
legend("topleft", c("detfun0", "EDM", "obs"), pch=1, col=c(4,2,3), bty="n")

#note improvement in fit, for both filters
cor(datout$true, datout$obs)^2 #observations
cor(datout$true, filterout_det$Nest)^2 #deterministic filter
cor(datout$true, filterout_edm$Nest)^2 #EDM filter

## Run optimizers
#\dontrun{
#create priors
minvUSE<-c(-4, -4, -4) #minimum interval for obs, proc, and rgr
maxvUSE<-c(0, 0, 1) #maximum interval for obs, proc, and rgr

minvUSE_edm<-c(-4, -4, -4, -5) #minimum interval for obs, proc, rgr, and theta
maxvUSE_edm<-c(0, 0, 1, 2) #maximum interval for obs, proc, rgr, and theta

#density, sampler, and prior functions for deterministic function
density_fun_USE<-function(param) density_fun0(param = param, minv = minvUSE, maxv=maxvUSE)
sampler_fun_USE<-function(x) sampler_fun0(n = 1, minv = minvUSE, maxv=maxvUSE)
prior_USE <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                         lower = minvUSE, upper = maxvUSE)

#density, sampler, and prior functions for EDM function
density_fun_USE_edm<-function(param) density_fun0(param = param, minv = minvUSE_edm, maxv=maxvUSE_edm)
sampler_fun_USE_edm<-function(x) sampler_fun0(n = 1, minv = minvUSE_edm, maxv=maxvUSE_edm)
prior_edm <- createPrior(density = density_fun_USE_edm, sampler = sampler_fun_USE_edm,
                         lower = minvUSE_edm, upper = maxvUSE_edm)
## Run filter
niter<-2000 #number of steps for the MCMC sampler
N<-1e3 #number of particles
Euse<-2 #number of embedding dimensions

#likelihood and bayesian set-ups for deterministic functions
likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0, procfun = procfun_ct, N = N)
bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior_USE)

#likelihood and bayesian set-ups for EDM functions
likelihood_EDM<-function(x) {
  xuse<-x[1:3]
  tuse_edm<-exp(x[4])
  likelihood0(param = xuse, y=y, parseparam = parseparam0, procfun = procfun_ct,
              detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse_edm), N = N)
}

bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior_edm)

#run MCMC optimization
out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0, settings = list(iterations=niter, consoleUpdates=20))
out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=20))

#plot results, with a 200-step burn-in
plot(out_detfun0, start = 200)
plot(out_EDM, start = 200)

## extract and plot parameter distributions
smp_detfun0<-getSample(out_detfun0, start = 200)
smp_EDM<-getSample(out_EDM, start=200)

par(mfrow=c(2,3))
hist(smp_detfun0[,1], main="det. function", xlab="obs", breaks = 20); abline(v=pars_true$obs, col=2)
hist(smp_detfun0[,2], main="det. function", xlab="proc", breaks = 20); abline(v=pars_true$proc[1], col=2)
hist(smp_detfun0[,3], main="det. function", xlab="rgr", breaks = 20); abline(v=pars_true$proc[2], col=2)

hist(smp_EDM[,1], main="EDM function", xlab="obs", breaks = 20); abline(v=pars_true$obs, col=2)
hist(smp_EDM[,2], main="EDM function", xlab="proc", breaks = 20); abline(v=pars_true$proc[1], col=2)
hist(smp_EDM[,3], main="EDM function", xlab="rgr", breaks = 20); abline(v=pars_true$proc[2], col=2)

hist(smp_EDM[,4], main="EDM function", xlab="theta", breaks = 20)
#}
