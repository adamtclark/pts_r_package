#REMOVE
if(TRUE) {
  error
  rm(list=ls())
  setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/r_package/hpc/")
}

## Load functions:
commArgin<-commandArgs(1)
if(length(commArgin)==0) {
  commArgin<-round(runif(1)*1e6)
}

require(BayesianTools)
require(rEDM)
source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")


## Simulate data
pars<-list(obs=c(log(1e-2), log(0.1)),
           proc=c(-2, log(1.2)),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(log(3),log(1)))

datout<-makedynamics(n = 100, obs = pars$obs, proc = pars$proc, r = pars$det[1],
                     K = pars$det[2], pcol = pars$pcol)
y<-datout$obs


## Run filter
N = 1e3
#based on detful0
filterout_det<-particleFilterLL(y, pars=pars, N, detfun = detfun0)
#based on EDM
filterout_edm<-particleFilterLL(y, pars=pars, N, detfun = EDMfun0, edmdat = list(E=2))

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

plot(filterout_det$rN, y, xlim=range(c(filterout_det$rN, filterout_edm$rN)),
     xlab="predicted", ylab="observed", col=4); points(filterout_edm$rN, y, col=2)
abline(a=0, b=1, lty=2)
legend("topleft", c("detfun0", "EDM"), pch=1, col=c(4,2), bty="n")


## Run optimizer
param<-unlist(pars)[1:6]

#create priors
density_fun_USE<-function(param) density_fun0(param = param, pars = pars)
sampler_fun_USE<-function(x) sampler_fun0(n = 1, pars = pars)
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = rep(-10,6), upper = rep(10,6))

likelihood0(param, y = y)+density_fun_USE(param)
likelihood0(param, y=y, detfun = EDMfun0, edmdat = list(E=2))+density_fun_USE(param)

#plot priors
tmp<-inv_fun0(sampler_fun0(1e4,pars = pars))
hist(tmp[,1], breaks = 30); abline(v=exp(pars$obs[1]), col=2, lwd=2, lty=2)
hist(tmp[,2], breaks = 30); abline(v=exp(pars$obs[2]), col=2, lwd=2, lty=2)
hist(tmp[,3], breaks = 30); abline(v=(pars$proc[1]), col=2, lwd=2, lty=2)
hist(tmp[,4], breaks = 30); abline(v=exp(pars$proc[2]), col=2, lwd=2, lty=2)
hist(tmp[,5], breaks = 30); abline(v=ilogit(pars$pcol[1]), col=2, lwd=2, lty=2)
hist(tmp[,6], breaks = 30); abline(v=exp(pars$pcol[2]), col=2, lwd=2, lty=2)

#number of MCMC iterations - increase for more accurate results
#note - runtime will be long for EDM example
niter<-501
nburn<-102

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


## Summarize outputs
smp_detfun0<-inv_fun0(getSample(out_detfun0))
smp_EDM<-inv_fun0(getSample(out_EDM))

plot(data.frame(smp_detfun0))
plot(data.frame(smp_EDM))

#plot posteriors
truepars_transformed<-inv_fun0(t(as.matrix(unlist(pars))))

par(mar=c(4,4,2,2), mfrow=c(3,2))
for(i in 1:6) {
  hist(smp_detfun0[,i],breaks = 30); abline(v=truepars_transformed[i], col=c(2), lty=2)
}

for(i in 1:6) {
  hist(smp_EDM[,i],breaks = 30); abline(v=truepars_transformed[i], col=c(2), lty=2)
}

## Plot demographic rates
pfout_detfun0<-particleFilterLL(y, pars=pars, N, detfun = detfun0)
pfout_EDM<-particleFilterLL(y, pars=pars, N, detfun = EDMfun0, edmdat = list(E=2))

par(mar=c(4,4,2,2), mfrow=c(2,1))
matplot(1:nrow(pfout_detfun0$dem$col),
        cbind(pfout_detfun0$dem$col[,1], pfout_EDM$dem$col[,1]),
        pch=1, col=c(2,4), xlab="time", ylab="pr[col]")
matplot(1:nrow(pfout_detfun0$dem$col),
        cbind(pfout_detfun0$dem$mor[,1], pfout_EDM$dem$mor[,1]),
        pch=1, col=c(2,4), xlab="time", ylab="pr[ext]")
