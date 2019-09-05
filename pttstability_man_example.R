error
rm(list=ls())


## Load functions:
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/pttstability/")
require(BayesianTools)
require(rEDM)
source("R/bayesfun.R")
source("R/fake_data.R")
source("R/logit_funs.R")
source("R/particlefilter.R")

#Set seed
set.seed(41019)

## Simulate data
pars0<-list(obs=log(0.1),
            proc=log(0.1),
            pcol=c(logit(0.2), log(1e-2)),
            det=c(log(3),log(1)))

#generate random parameter values
pars_true<-parseparam0(sampler_fun0(n=1, pars = pars0))

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

plot(filterout_det$rN, y, xlim=range(c(filterout_det$rN, filterout_edm$rN)),
     xlab="predicted", ylab="observed", col=4); points(filterout_edm$rN, y, col=2)
abline(a=0, b=1, lty=2)
legend("topleft", c("detfun0", "EDM"), pch=1, col=c(4,2), bty="n")

#note improvement in fit, for both filters
cor(datout$true, datout$obs)^2
cor(datout$true, filterout_det$rN)^2
cor(datout$true, filterout_edm$rN)^2


#estimate demographic rates from extended timeseries
etdfilter_det<-extend_particleFilter(pfout=filterout_det, pars=pars_true,
                                     Next = 1e3, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL)
etdfilter_edm<-extend_particleFilter(pfout=filterout_edm, pars=pars_true,
                                     Next = 1e3, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=2))
#compare to results from a long simulation
datout_long<-makedynamics_general(n = 2e4, n0 = exp(rnorm(1,0,0.1)), pdet=pars_true$det,
                          proc = pars_true$proc, obs = pars_true$obs, pcol = pars_true$pcol,
                          detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
demdat_true<-getcm(datout_long$true)

#plot rates
par(mar=c(4,4,2,2), mfrow=c(2,1))
hist(etdfilter_det$demdat$text, xlim=c(40,160),
     col=adjustcolor("blue", alpha.f = 0.5), breaks = 30, xlab="time to extinction", main="")
par(new=TRUE)
hist(etdfilter_edm$demdat$text,xlim=c(40,160),
     col=adjustcolor("red", alpha.f = 0.5), breaks = 30, axes=FALSE, xlab="", ylab="", main="")
abline(v=1/demdat_true$pm, lwd=2, lty=2)

hist(etdfilter_det$demdat$tcol, xlim=c(5, 20),
     col=adjustcolor("blue", alpha.f = 0.5), breaks = 30, xlab="time to colonization", main="")
par(new=TRUE)
hist(etdfilter_edm$demdat$tcol, xlim=c(5, 12),
     col=adjustcolor("red", alpha.f = 0.5), breaks = 30, axes=FALSE, xlab="", ylab="", main="")
abline(v=1/demdat_true$pc, lwd=2, lty=2)

legend("topright", c("deterministic est.", "EDM est.", "true value"),
       fill = adjustcolor(c("blue", "red", NA), alpha.f = 0.5),
       lty=c(NA, NA, 2), lwd=c(NA, NA, 2), col=c(NA, NA, 1), border = c(1,1,NA))


## Run optimizer
#\dontrun{
param<-unlist(pars0)[1:3]

#create priors
density_fun_USE<-function(param) density_fun0(param = param, pars = pars0)
sampler_fun_USE<-function(x) sampler_fun0(n = 1, pars = pars0)
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = rep(-10,3), upper = rep(10,3))

likelihood0(param, y = y)+density_fun_USE(param)
likelihood0(param, y=y, detfun = EDMfun0, edmdat = list(E=2))+density_fun_USE(param)

#plot priors
par(mar=c(4,4,2,2), mfrow=c(3,1))
tmp<-inv_fun0(sampler_fun0(1e4,pars = pars0))
hist(tmp[,1], breaks = 30); abline(v=exp(pars0$obs[1]), col=2, lwd=2, lty=2)
hist(tmp[,2], breaks = 30); abline(v=exp(pars0$proc[1]), col=2, lwd=2, lty=2)
hist(tmp[,3], breaks = 30); abline(v=ilogit(pars0$pcol[1]), col=2, lwd=2, lty=2)

#number of MCMC iterations - increase for more accurate results
#note - runtime will be long for EDM example
niter<-1000
nburn<-500

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

#check for correlations among estimates
plot(data.frame(smp_detfun0))
plot(data.frame(smp_EDM))

#plot posteriors
truepars_transformed<-inv_fun0(t(as.matrix(unlist(pars_true))))
priorpars_transformed<-inv_fun0(t(as.matrix(unlist(pars0))))

par(mar=c(4,4,2,2), mfrow=c(1,3))
for(i in 1:3) {
  hist(smp_detfun0[,i],breaks = 30)
  abline(v=truepars_transformed[i], col=c(2), lty=2)
  abline(v=truepars_transformed[i], col=c(2), lty=2)
  abline(v=priorpars_transformed[i], col=c(4), lty=2)
}

for(i in 1:6) {
  hist(smp_EDM[,i],breaks = 30); abline(v=truepars_transformed[i], col=c(2), lty=2)
}

## Plot demographic rates
pfout_detfun0<-particleFilterLL(y, pars=pars, N, detfun = detfun0, dotraceback = TRUE)
pfout_EDM<-particleFilterLL(y, pars=pars, N, detfun = EDMfun0, edmdat = list(E=2),
                            dotraceback = TRUE)

par(mar=c(4,4,2,2), mfrow=c(2,1))
matplot(1:nrow(pfout_detfun0$dem$col),
        cbind(pfout_detfun0$dem$col[,1]/pfout_detfun0$dem$col[,2],
              pfout_EDM$dem$col[,1]/pfout_EDM$dem$col[,2]),
        pch=1, col=c(2,4), xlab="time", ylab="pr[col]")
matplot(1:nrow(pfout_detfun0$dem$col),
        cbind(pfout_detfun0$dem$mor[,1]/pfout_detfun0$dem$mor[,2],
              pfout_EDM$dem$mor[,1]/pfout_EDM$dem$mor[,2]),
        pch=1, col=c(2,4), xlab="time", ylab="pr[ext]")
#}
