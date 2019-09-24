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
set.seed(190908)

## Simulate data
pars_true<-list(obs=log(0.2),
            proc=log(0.1),
            pcol=c(logit(0.2), log(0.1)),
            det=c(log(3),log(1)))

#generate random parameter values
datout<-makedynamics_general(n = 100, n0 = exp(rnorm(1,0,0.1)), pdet=pars_true$det,
                             proc = pars_true$proc, obs = pars_true$obs, pcol = pars_true$pcol,
                             detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
y<-datout$obs

## Run filter
N = 1e3
#based on detful0
filterout_det<-particleFilterLL(y, pars=pars_true, N, detfun = detfun0,
                                dotraceback = TRUE)
#based on EDM
filterout_edm<-particleFilterLL(y, pars=pars_true, N, detfun = EDMfun0, edmdat = list(E=2),
                                dotraceback = TRUE)

#plot filter output
par(mar=c(4,4,2,2), mfrow=c(3,1))
matplot(1:length(y), t(filterout_det$x)[,1:30], col=adjustcolor(1,alpha.f = 0.5), lty=3,
        type="l", xlab="time", ylab="abund", main="detfun0")
lines(1:length(y), y, col=2, lwd=1.5)
lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2)
lines(1:length(y), filterout_det$rN, col=3, lwd=1.5)

matplot(1:length(y), t(filterout_edm$x)[,1:30], col=adjustcolor(1,alpha.f = 0.5), lty=3,
        type="l", xlab="time", ylab="abund", main="EDM")
lines(1:length(y), y, col=2, lwd=1.5)
lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2)
lines(1:length(y), filterout_edm$rN, col=3, lwd=1.5)

plot(filterout_det$rN, datout$true, xlim=range(c(filterout_det$rN, filterout_edm$rN)),
     xlab="predicted", ylab="true value", col=4)
points(filterout_edm$rN, datout$true, col=2)
points(y, datout$true, col=3)
abline(a=0, b=1, lty=2)
legend("topleft", c("detfun0", "EDM", "obs"), pch=1, col=c(4,2,3), bty="n")

#note improvement in fit, for both filters
cor(datout$true, datout$obs)^2
cor(datout$true, filterout_det$rN)^2
cor(datout$true, filterout_edm$rN)^2

## estimate extinction rate
demdat_true<-getcm(datout$true)
demdat_obs<-getcm(datout$obs)

demdat_det<-rbinom(n = 1e4, size = length(y), filterout_det$dem$mumor)/length(y)
demdat_edm<-rbinom(n = 1e4, size = length(y), filterout_edm$dem$mumor)/length(y)

#plot rates
par(mar=c(4,4,2,2), mfrow=c(1,1))
ddet<-density(demdat_det, from = 0, bw = 0.01)
dedm<-density(demdat_edm, from = 0, bw = 0.01)

xlm<-range(c(ddet$x, dedm$x, demdat_true$pm, demdat_obs$pm), na.rm=T)
ylm<-range(c(ddet$y, dedm$y), na.rm=T)

plot(ddet$x, ddet$y, xlim=xlm, ylim=ylm, type="l", lwd=2,
     col=adjustcolor("blue", alpha.f = 0.5),
     xlab="extinction rate", ylab="density", main="")
lines(dedm$x, dedm$y, xlim=xlm, lwd=2,
     col=adjustcolor("red", alpha.f = 0.5))
abline(v=demdat_true$pm, lwd=2, lty=2, col=1)
abline(v=demdat_obs$pm, lwd=2, lty=2, col=3)
abline(h=0, v=0, lty=3)

legend("topright", c("det. est.", "EDM est.", "true value", "raw est."),
       fill = adjustcolor(c("blue", "red", NA, NA), alpha.f = 0.5),
       lty=c(NA, NA, 2, 2), lwd=c(NA, NA, 2, 2), col=c(NA, NA, 1, 3), border = c(1,1,NA,NA))

## Run optimizers
#\dontrun{

#prior parameters for optimizer - note that these differ from pars_true
pars0<-list(obs=pars_true$obs-1,
            proc=pars_true$proc+1,
            pcol=c(pars_true$pcol[1]-1, pars_true$pcol[2]),
            det=pars_true$det)
ptrue<-unlist(pars_true)[1:3]
p0<-unname(unlist(pars0)[1:3])

#likelihood function for EDM
likelihoodEDM<-function(param, y, parseparam, N) {likelihood0(param, y, parseparam, N, detfun = EDMfun0, edmdat = list(E=2))}

## ABC optimizer
#runs very slowly - load pre-compiled results instead.
#set to TRUE if you want to run the optimizer yourself
if(FALSE) {
  #run optimizer for deterministic model
  optout_det<-run_ABC_optim(y, sd0 = c(2,2,2), p0,likelihood = likelihood0, fretain = 0.5)

  #extend run for another 20 iterations
  optout_det_ext<-run_ABC_optim(oldrun = optout_det, niter_optim = 20)

  #run optimizer for EDM model (will run ~10x slower than deterministic model)
  optout_edm<-run_ABC_optim(y, sd0 = c(2,2,2), p0, likelihood = likelihoodEDM, fretain = 0.5)
  #extend run for another 20 iterations
  optout_edm_ext<-run_ABC_optim(oldrun = optout_edm, niter_optim = 20)
} else {
  load("data/abcout_rv2.rda") #functions run slowly - load pre-run results
}

## plotting
#plotting likelihoods, deterministic function
par(mfrow=c(2,3), mar=c(4,4,2,2))
plot_abc_likelihood(optout_det, logx = TRUE)
plot_abc_rmse(optout_det, ptrue)
#plot results for extended run
plot_abc_likelihood(optout_det_ext, logx = TRUE)
plot_abc_rmse(optout_det_ext, ptrue)

#EDM
plot_abc_likelihood(optout_edm, logx = TRUE)
plot_abc_rmse(optout_edm, ptrue)
#EDM extended
plot_abc_likelihood(optout_edm_ext, logx = TRUE)
plot_abc_rmse(optout_edm_ext, ptrue)

#plotting parameter estimates
par(mfrow=c(2,3), mar=c(4,4,2,2))
#deterministic run
plot_abc_params(optout_det, param0 = p0, param_true = ptrue)
#extended deterministic run
plot_abc_params(optout_det_ext, param0 = p0, param_true = ptrue)
#EDM run
plot_abc_params(optout_edm, param0 = p0, param_true = ptrue)
#EDM run extended
plot_abc_params(optout_edm_ext, param0 = p0, param_true = ptrue)
#true values are in red, priors are in blue

#extracting estimates from likelihood surfaces
par(mfrow=c(2,3), mar=c(4,4,2,2))
dens_out_det<-abc_densities(optout = optout_det, param0 = p0, param_true = ptrue, fretain = 0.75, enp.target = 4, nbootstrap = 100, nobs = length(y))
dens_out_det_ext<-abc_densities(optout = optout_det_ext, param0 = p0, param_true = ptrue, fretain = 0.75, enp.target = 4, nbootstrap = 100, nobs = length(y))
dens_out_edm<-abc_densities(optout = optout_edm, param0 = p0, param_true = ptrue, fretain = 0.75, enp.target = 4, nbootstrap = 100, nobs = length(y))
dens_out_edm_ext<-abc_densities(optout = optout_edm_ext, param0 = p0, param_true = ptrue, fretain = 0.75, enp.target = 4, nbootstrap = 100, nobs = length(y))
#true values are in red, priors are in blue

#show estimates
dens_out_det_ext
dens_out_edm_ext
ptrue

## MCMC optimizer
#create priors
density_fun_USE<-function(param) density_fun0(param = param, pars = pars0)#/length(y)
sampler_fun_USE<-function(x) sampler_fun0(n = 1, pars = pars0, priorsd = c(1,1,1), minv = rep(-4, 3), maxv = rep(2,3))
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = rep(-4,3), upper = rep(2,3))

#set number of iterations
niter<-1e4

#set up likelihoods
likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0, N = N)*length(y)
bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior)

likelihood_EDM<-function(x) likelihood0(param = x, y=y, parseparam = parseparam0,
                                        detfun = EDMfun0, edmdat = list(E=2), N = N)*length(y)
bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior)

#run MCMC chains
#runs very slowly - load pre-compiled results instead.
#set to TRUE if you want to run the optimizer yourself
if(FALSE) {
  out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0, settings = list(iterations=niter, consoleUpdates=10))
  out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=10))
} else {
  load("data/mcmcout_rv2.rda")
}

plot(out_detfun0, start=500)
plot(out_EDM, start=500)
gelmanDiagnostics(out_detfun0, plot = FALSE, start=500)
gelmanDiagnostics(out_EDM, plot = FALSE, start=500)
#note, convergence is worse for EDM - needs to run for longer

## Summarize outputs
smp_detfun0<-getSample(out_detfun0, start = 500)
smp_EDM<-getSample(out_EDM, start=500)

#check for correlations among estimates
correlationPlot(out_detfun0)
correlationPlot(out_EDM)

#plot priors and posteriors
marginalPlot(out_detfun0, prior = TRUE)
marginalPlot(out_EDM, prior = TRUE)

#plot posteriors vs. true values
par(mar=c(4,4,2,2), mfrow=c(2,3))
for(i in 1:3) {
  xrng<-range(c(smp_detfun0[,i], ptrue[i], p0[i]))
  hist(smp_detfun0[,i],breaks = 20, probability = TRUE, main="", xlim=xrng);
  abline(v=ptrue[i], col=c(1), lty=2)
  abline(v=ptrue[i], col=c(3), lty=2)
}

for(i in 1:3) {
  xrng<-range(c(smp_EDM[,i], ptrue[i], p0[i]))
  hist(smp_EDM[,i],breaks = 20, probability = TRUE, main="", xlim=xrng);
  abline(v=ptrue[i], col=c(1), lty=2)
  abline(v=p0[i], col=c(3), lty=2)
}

#}
