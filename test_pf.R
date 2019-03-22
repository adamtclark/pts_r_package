error
rm(list=ls())

#LOAD:
source("~/Dropbox/Rfunctions/logit_funs.R")
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/pttstability/")
source("R/fake_data.R"); source("R/particlefilter.R")
library(BayesianTools)

########################
## Run simulation
########################
par(mar=c(4,4,2,2), mfrow=c(1,1))

#set.seed(1403)
pars0<-list(obs=c(log(1e-2), log(0.1)),
           proc=c(-2, log(1.2)),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(3,1))

datout<-makedynamics(n = 100, obs = pars0$obs, proc = pars0$proc, r = pars0$det[1], K = pars0$det[2], pcol = pars0$pcol)
y<-datout$obs

N<-1e3
#set.seed(1403)
system.time(tmp<-particleFilterLL(y, pars=pars0, N, detfun = detfun0, edmdat = NULL))#, obsfun, procfun, detfun))
system.time(tmp2<-particleFilterLL(y, pars=pars0, N, detfun = EDMfun0, edmdat = list(E=2)))#, obsfun, procfun, detfun))

tmp$LL
tmp2$LL

par(mar=c(4,4,2,2), mfrow=c(2,1))
matplot(1:length(y), t(tmp$x)[,1:100], col=adjustcolor(1,alpha.f = 0.5), lty=3, type="l", xlab="time", ylab="abund", ylim=c(0, max(c(y, tmp$x[,1:100], tmp$x2[,1:100]))))
lines(1:length(y), y, col=2, lwd=1.5)
lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2)
lines(1:length(y), tmp$rN, col=3, lwd=1.5)

matplot(1:length(y), t(tmp2$x)[,1:100], col=adjustcolor(1,alpha.f = 0.5), lty=3, type="l", xlab="time", ylab="abund", ylim=c(0, max(c(y, tmp$x[,1:100], tmp$x2[,1:100]))))
lines(1:length(y), y, col=2, lwd=1.5)
lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2)
lines(1:length(y), tmp2$rN, col=3, lwd=1.5)

par(mar=c(4,4,2,2), mfrow=c(1,1))
plot(tmp$rN, y, xlim=range(c(tmp$rN, tmp2$rN))); points(tmp2$rN, y, col=2); abline(a=0, b=1, lty=2)

########################
## Run optimizer
########################

#New parameters for fitting
pars<-list(obs=c(log(1e-2), log(0.1)),
           proc=c(-2, log(1.2)),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(3,1))

#Likelihood
likelihood <- function(param) {
  pars<-list(obs=c(param[1], param[2]),
             proc=c(param[3], param[4]),
             pcol=c(param[5], param[6]),
             det=c(3,1))

  #tmp<-particleFilterLL(y, pars, N=1e3)
  tmp<-particleFilterLL(y, pars, N=1e3, detfun = EDMfun0, edmdat = list(E=2))
  tmp$LL
}


param<-c(pars$obs, pars$proc, pars$pcol)
param0<-c(pars0$obs, pars0$proc, pars0$pcol)

niter<-1
x<-0; for(i in 1:niter) {x<-x+likelihood(param)};
param2<-param*c(1,1,1,1,1,1)
x2<-0; for(i in 1:niter) {x2<-x2+likelihood(param2)};
(x/niter-x2/niter) #positive means correct values are better
x/niter


#Prior
density_fun = function(param){
  priorsd<-0.5
  d1 = dnorm(param[1], mean = pars$obs[1]+priorsd^2, sd =  priorsd, log = TRUE)
  d2 = dnorm(param[2], mean = pars$obs[2]+priorsd^2, sd =  priorsd, log = TRUE)
  d3 = dnorm(param[3], mean = pars$proc[1], sd = 2*priorsd, log = TRUE)
  d4 = dnorm(param[4], mean = pars$proc[2]+priorsd^2, sd = priorsd, log = TRUE)
  d5 = dnorm(param[5], mean = pars$pcol[1] - (priorsd)^2*(2*ilogit(pars$pcol[1])-1), sd = priorsd, log = TRUE)
  d6 = dnorm(param[6], mean = pars$pcol[2]+priorsd^2, sd = priorsd, log = TRUE)

  return(d1 + d2 + d3 + d4 + d5 + d6)
}

sampler_fun = function(n=1){
  priorsd<-0.5
  d1 = rnorm(n, mean = pars$obs[1]+priorsd^2, sd = priorsd)
  d2 = rnorm(n, mean = pars$obs[2]+priorsd^2, sd = priorsd)
  d3 = rnorm(n, mean = pars$proc[1], sd = 2*priorsd)
  d4 = rnorm(n, mean = pars$proc[2]+priorsd^2, sd=priorsd)
  d5 = rnorm(n, mean = pars$pcol[1] - (priorsd)^2*(2*ilogit(pars$pcol[1])-1), sd = priorsd)
  d6 = rnorm(n, mean = pars$pcol[2]+priorsd^2, sd=priorsd)
  return(cbind(d1,d2,d3,d4,d5,d6))
}

inv_fun<-function(x) {
  cbind(exp(x[,1]), exp(x[,2]), x[,3], exp(x[,4]), ilogit(x[,5]), exp(x[,6]))
}

prior <- createPrior(density = density_fun, sampler = sampler_fun,
                     lower = rep(-10,6), upper = rep(10,6), best = NULL)

likelihood(param)+density_fun(param)
likelihood(param0)+density_fun(param0)

#plot priors
tmp<-inv_fun(sampler_fun(1e4))
hist(tmp[,1], breaks = 30); abline(v=exp(pars0$obs[1]), col=2, lwd=2, lty=2)
hist(tmp[,2], breaks = 30); abline(v=exp(pars0$obs[2]), col=2, lwd=2, lty=2)
hist(tmp[,3], breaks = 30); abline(v=(pars0$proc[1]), col=2, lwd=2, lty=2)
hist(tmp[,4], breaks = 30); abline(v=exp(pars0$proc[2]), col=2, lwd=2, lty=2)
hist(tmp[,5], breaks = 30); abline(v=ilogit(pars0$pcol[1]), col=2, lwd=2, lty=2)
hist(tmp[,6], breaks = 30); abline(v=exp(pars0$pcol[2]), col=2, lwd=2, lty=2)


param<-c(pars$obs[1], pars$obs[2], pars$proc[1], pars$proc[2], pars$pcol[1], pars$pcol[2])
likelihood(param)

bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior)

#settings <- list(iterations = 4002, burnin=1002, adapt = T, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = F, nrChains = 1)
#out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = list(iterations=4002, burnin=1002))


########################
## Summarize outputs
########################

## plot parameter values

#plot(out)
summary(out)

smp<-inv_fun(getSample(out))[-c(1:10),]
plot(data.frame(smp))
cor(smp)

par(mar=c(4,4,2,2), mfrow=c(3,2))
plot(density(smp[,1],bw = 0.001, from = 0, to = 0.05)); abline(v=c(exp(pars0$obs[1]), exp(pars$obs[1])), col=c(2,4), lty=2)
plot(density(smp[,2],bw = 0.005, from = 0, to = 0.25)); abline(v=c(exp(pars0$obs[2]), exp(pars$obs[2])), col=c(2,4), lty=2)
plot(density(smp[,3],bw = 0.05, from = -3, to = 0)); abline(v=c(pars0$proc[1], pars$proc[1]), col=c(2,4), lty=2)
plot(density(smp[,4],bw = 0.05, from = 0, to = 2)); abline(v=c(exp(pars0$proc[2]), exp(pars$proc[2])), col=c(2,4), lty=2)
plot(density(smp[,5],bw = 0.01, from = 0, to = 1)); abline(v=c(ilogit(pars0$pcol[1]), ilogit(pars$pcol[1])), col=c(2,4), lty=2)
plot(density(smp[,6],bw = 0.001, from = 0, to = 0.05)); abline(v=c(exp(pars0$pcol[2]), exp(pars$pcol[2])), col=c(2,4), lty=2)

par(mar=c(4,4,2,2), mfrow=c(3,2))
hist(smp[,1],breaks = 30); abline(v=c(exp(pars0$obs[1]), exp(pars$obs[1])), col=c(2,4), lty=2)
hist(smp[,2],breaks = 30); abline(v=c(exp(pars0$obs[2]), exp(pars$obs[2])), col=c(2,4), lty=2)
hist(smp[,3],breaks = 30); abline(v=c(pars0$proc[1], pars$proc[1]), col=c(2,4), lty=2)
hist(smp[,4],breaks = 30); abline(v=c(exp(pars0$proc[2]), exp(pars$proc[2])), col=c(2,4), lty=2)
hist(smp[,5],breaks = 30); abline(v=c(ilogit(pars0$pcol[1]), ilogit(pars$pcol[1])), col=c(2,4), lty=2)
hist(smp[,6],breaks = 30); abline(v=c(exp(pars0$pcol[2]), exp(pars$pcol[2])), col=c(2,4), lty=2)

apply(smp, 2, function(x) quantile(x, c(0.025, pnorm(c(-1,0,1)), 0.975)))



## demographic rates
pfout_1<-particleFilterLL(y, pars=pars0, N, detfun = detfun0, edmdat = NULL)
pfout_2<-particleFilterLL(y, pars=pars0, N, detfun = EDMfun0, edmdat = list(E=2))


pfout_1$dem$col<-cbind(pfout_1$dem$col[,1:2], pfout_1$dem$col[,1]*(1-pfout_1$dem$col[,1])/sqrt(pfout_1$dem$col[,2]))
pfout_1$dem$mor<-cbind(pfout_1$dem$mor[,1:2], pfout_1$dem$mor[,1]*(1-pfout_1$dem$mor[,1])/sqrt(pfout_1$dem$mor[,2]))
pfout_2$dem$col<-cbind(pfout_2$dem$col[,1:2], pfout_2$dem$col[,1]*(1-pfout_2$dem$col[,1])/sqrt(pfout_2$dem$col[,2]))
pfout_2$dem$mor<-cbind(pfout_2$dem$mor[,1:2], pfout_2$dem$mor[,1]*(1-pfout_2$dem$mor[,1])/sqrt(pfout_2$dem$mor[,2]))

#colonization
clmat<-cbind(pfout_1$dem$col[,1], pfout_2$dem$col[,1])
clmat[clmat==0]<-NA
ps<-which(is.finite(rowSums(clmat)))

matplot(ps, clmat[ps,], pch=1, col=c(2,4), type="b", lty=3, ylim=c(0.1, 0.17), xlab="time", ylab="pr[col]")
segments(ps, clmat[ps,1]+pfout_1$dem$col[ps,3], ps, clmat[ps,1]-pfout_1$dem$col[ps,3], col=2)
segments(ps, clmat[ps,2]+pfout_2$dem$col[ps,3], ps, clmat[ps,2]-pfout_2$dem$col[ps,3], col=4)
abline(h=sum(clmat[ps,1]*(1/pfout_1$dem$col[ps,3]))/sum(1/pfout_1$dem$col[ps,3]), lty=2, col=2)
abline(h=sum(clmat[ps,2]*(1/pfout_2$dem$col[ps,3]))/sum(1/pfout_2$dem$col[ps,3]), lty=2, col=4)

#extinction
mrmat<-cbind(pfout_1$dem$mor[,1], pfout_2$dem$mor[,1])
mrmat[mrmat==0]<-NA
ps<-which(is.finite(rowSums(mrmat)))

matplot(ps, mrmat[ps,], pch=1, col=c(2,4), type="b", lty=3, ylim=c(0, 0.42), xlab="time", ylab="pr[mor]")
segments(ps, mrmat[ps,1]+pfout_1$dem$mor[ps,3], ps, mrmat[ps,1]-pfout_1$dem$mor[ps,3], col=2)
segments(ps, mrmat[ps,2]+pfout_2$dem$mor[ps,3], ps, mrmat[ps,2]-pfout_2$dem$mor[ps,3], col=4)
abline(h=sum(mrmat[ps,1]*sqrt(pfout_1$dem$mor[ps,2]))/sum(sqrt(pfout_1$dem$mor[ps,2])), lty=2, col=2)
abline(h=sum(mrmat[ps,2]*sqrt(pfout_2$dem$mor[ps,2]))/sum(sqrt(pfout_2$dem$mor[ps,2])), lty=2, col=4)

matplot(datout$true, mrmat, pch=1, col=c(2,4), xlab="abundance", ylab="pr[mor]")



#NOTES:
#2. Package up test script
#3. Run test of fitting across parameter ranges (with real function, and with EDM)
#4. Send to Benjamin - think about R package structure... (ideally with lmer notation for blocks etc.)


