#!/usr/bin/env Rscript

#error
#rm(list=ls())

## Load functions:
#setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/pttstability/")
require(BayesianTools)
require(rEDM)
source("R/bayesfun.R")
source("R/fake_data.R")
source("R/logit_funs.R")
source("R/particlefilter.R")

#Set seed
set.seed(2110)

## Simulate data
pars_true<-list(obs=log(0.2),
            proc=c(log(0.25), log(1.2)),
            pcol=c(logit(0.2), log(0.1)),
            det=c(log(1.2),log(1)))
#abstracted process noise standard dev., based on formula for
#a linear process: see ?sdproc_abstract for details.
sd_abs = sdproc_abstract(exp(pars_true$proc[1]), exp(pars_true$proc[2]))
#parameters for the filter
pars_filter<-pars_true; pars_filter$proc = pars_true$proc[1]

#generate random parameter values
datout<-makedynamics_general(n = 100, n0 = exp(pars_true$det[2]), pdet=pars_true$det,
                             proc = pars_true$proc, obs = pars_true$obs, pcol = pars_true$pcol,
                             detfun = detfun0_sin, procfun = procfun_ct, obsfun=obsfun0, colfun=colfun0)
y<-datout$obs
plot(y, type = "l", xlab="time", ylab="observed abundance")

#get theta paramter for s-mapping
# assume best E = 4
# alternatively, we could run e.g. s_map(y, E=(2:5), silent = TRUE)
# to test a range of potential E values
Euse = 4
#run leave-one-out cross validation
s<-s_map(y, E=Euse, silent = TRUE)
tuse<-s$theta[which.min(s$rmse)] # retain best theta
plot(s$theta, s$rho, type="b")

## Run filter with "correct" parameter values
N = 1e3 # number of particles
#based on detful0
filterout_det<-particleFilterLL(y, pars=pars_filter, N, detfun = detfun0_sin, procfun = procfun0,
                                dotraceback = TRUE, fulltraceback = TRUE)
#based on EDM
filterout_edm<-particleFilterLL(y, pars=pars_filter, N, detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse), procfun = procfun0,
                                dotraceback = TRUE, fulltraceback = TRUE)

#plot filter output
par(mar=c(4,4,2,2), mfrow=c(3,1))
#plot 30 of the 1000 particles to show general trend
# correct deterministic function
matplot(1:length(y), filterout_det$fulltracemat[,1:30], col=adjustcolor(1,alpha.f = 0.5), lty=3,
        type="l", xlab="time", ylab="abund", main="detfun0")
lines(1:length(y), y, col=2, lwd=1.5) #observations
lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2) #true values
lines(1:length(y), filterout_det$Nest, col=3, lwd=1.5) #mean estimate from filter

# EDM function
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
minvUSE<-c(-4, -4) #minimum interval for obs and proc
maxvUSE<-c(0, 0) #maximum interval for obs and proc

minvUSE_edm<-c(-4, -4) #minimum interval for obs and proc
maxvUSE_edm<-c(0, 0) #maximum interval for obs and proc

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
niter<-5000 #number of steps for the MCMC sampler
N<-1e3 #number of particles
Euse<-Euse #number of embedding dimensions

#likelihood and bayesian set-ups for deterministic functions
likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0, detfun = detfun0_sin, procfun = procfun0, N = N)
bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior_USE)

#likelihood and bayesian set-ups for EDM functions
likelihood_EDM<-function(x) {
  likelihood0(param = x, y=y, parseparam = parseparam0, procfun = procfun0,
              detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse), N = N)
}

bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior_edm)

#run MCMC optimization
out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0, settings = list(iterations=niter, consoleUpdates=20))
out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=20))

#plot results, with a 1000-step burn-in
plot(out_detfun0, start = 1000, thin = 2)
plot(out_EDM, start = 1000, thin = 2)

## extract and plot parameter distributions
smp_detfun0<-getSample(out_detfun0, start = 1000, thin = 2)
smp_EDM<-getSample(out_EDM, start=1000, thin = 2)

par(mfrow=c(2,2))
hist(exp(smp_detfun0[,1]), main="det. function", xlab="obs", breaks = 20); abline(v=exp(pars_true$obs), col=2)
hist(exp(smp_detfun0[,2]), main="det. function", xlab="proc", breaks = 20); abline(v=sd_abs, col=2)

hist(exp(smp_EDM[,1]), main="EDM function", xlab="obs", breaks = 20); abline(v=exp(pars_true$obs), col=2)
hist(exp(smp_EDM[,2]), main="EDM function", xlab="proc", breaks = 20); abline(v=sd_abs, col=2)

## compare total EDM coefficient prediction error to total model error
s_full<-s_map(y, E=Euse, theta=tuse, silent = TRUE)

# over-estimation of stochastic variance by EDM
(mean(exp(smp_EDM[,1])^2 + exp(smp_EDM[,2])^2)-(sd_abs^2+exp(pars_true$obs[1])^2))
# rmse due to EDM fitting error
s_full$rmse[[1]]^2-(sd_abs^2+exp(pars_true$obs[1])^2)
#}
