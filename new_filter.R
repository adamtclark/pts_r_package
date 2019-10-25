#!/usr/bin/env Rscript

error
rm(list=ls())

## Load functions:
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/pttstability/")
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
                det=c(log(2),log(1)))

#obsfun0<-function(so, yt, xt=NULL, inverse=FALSE, N=NULL, time=NULL) {
#  if(inverse) {
#    pmax(0, rnorm(n = N, mean = yt, sd = exp(so[1])*yt))
#  } else {
#    std_tmp<-exp(so[1])*yt
#    #Tobit distribution:
#    if(yt==0) {
#      pnorm(0, mean=xt/std_tmp, log.p = TRUE)
#    } else {
#      dnorm((yt-xt)/std_tmp,log=TRUE)-log(std_tmp)
#    }
#  }
#}

#generate random parameter values
datout<-makedynamics_general(n = 100, n0 = exp(rnorm(1,0,0.1)), pdet=pars_true$det,
                             proc = pars_true$proc, obs = pars_true$obs, pcol = pars_true$pcol,
                             detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
y<-datout$obs
plot(y, type="b")

s<-simplex(y, silent = TRUE); plot(s$E, s$rho, type="l")
s$E[which.min(s$rmse)]
s<-s_map(y, E=2, silent = TRUE); plot(s$theta, s$rho)
tuse<-s$theta[which.min(s$rmse)]

#try simple filter
plot(y, type="l")

EDMfun0<-function(smp_cf, yp, x, time) {
  nD<-ncol(smp_cf)
  if(nD>2) {
    sum(smp_cf[time,-1]*c(yp[(time-(nD-2)):(time-1)], 1))+smp_cf[time,1]*x
  } else {
    smp_cf[time,-1]+smp_cf[time,1]*x
  }
}



particleFilterLL<-function(y, pars, N=1e3, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=FALSE) {
  LL<-rep(NA, length(y))
  
  if(dotraceback) {
    Nest<-rep(NA, length(y))
    Nsd<-rep(NA, length(y))
  } else {
    Nest<-NULL
    Nsd<-NULL
  }
  
  #pr<-rnorm(N, 0, exp(pars$proc[1]))
  
  tstart<-max(c(2, edmdat$E))
  for(i in 1:(tstart)) {
    #prd<-detfun(pars$det, rnorm(N, y[i], exp(pars$obs)))
    #prd<-rnorm(N, y[i], exp(pars$obs))
    prd<-obsfun(so = pars$obs, yt = y[i], xt = proc, time = i, inverse = TRUE, N = N)
    if(dotraceback) {
      Nest[i]<-mean(prd)
      Nsd[i]<-sd(prd)
    }
  }
  
  if(!is.null(edmdat)) {
    smp<-s_map(y, E=edmdat$E, theta = edmdat$theta, silent = TRUE, save_smap_coefficients = TRUE)
    smp_cf<-smp$smap_coefficients[[1]]
  }
  
  for(i in tstart:length(y)) {
    #prd is deterministic estimate for t=i
    
    #add process noise
    #proc<-pr+prd
    proc<-procfun(sp = pars$proc, xt = prd, time = i)
    
    #get likelihood of proc given y[i]
    #dobs<-dnorm(y[i], proc, exp(pars$obs[1]),log=TRUE)
    dobs<-obsfun(so = pars$obs, yt = y[i], xt = proc, time = i)
    mxd<-max(dobs, na.rm=T)
    wdobs<-exp(dobs-mxd)/sum(exp(dobs-mxd))
    #dobs<-exp(dobs)/sum(exp(dobs))
    
    #estimates of true state for y[i]
    post_smp<-sample(proc, N, rep=T, prob = wdobs)
    
    #estimates of deterministic state at t=i+1
    #prd<-sum(smp_cf[i,2:3]*c(y[i-1], 1))+smp_cf[i,1]*post_smp
    if(is.null(edmdat)) {
      prd<-detfun(sdet = pars$det, xt = post_smp, time = i)
    } else {
      prd<-detfun(smp_cf = smp_cf, yp = y, x = post_smp, time = i)
    }
    
    #save likelihood
    LL[i]<-log(mean(exp(dobs-mxd)))+mxd
    #LL[i]<-log(mean(exp(dobs)))
    
    #save state
    if(dotraceback) {
      Nest[i]<-mean(post_smp)
      Nsd[i]<-sd(post_smp)
    }
  }
  LLtot <- sum(LL[is.finite(LL)], na.rm=T)
  
  return(list(LL = LLtot, LLlst=LL, Nest=Nest, Nsd=Nsd))
}


N<-2e4
pars_est<-pars_true
pars_est$obs<-log(0.2)
pars_est$proc<-log(0.1)
system.time({pfout1<-particleFilterLL(y=y, pars=pars_est, N=N, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=TRUE)})
pfout1$LL

pars_est$obs<-log(0.2)
pars_est$proc<-log(0.1)
system.time({pfout2<-particleFilterLL(y=y, pars=pars_est, N=N, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=2, theta=tuse), dotraceback=TRUE)})
pfout2$LL

#problem - cov with obs and proc


cor(pfout1$Nest, datout$true, use="pairwise")^2
cor(pfout2$Nest, datout$true, use="pairwise")^2
cor(y, datout$true, use="pairwise")^2



#test fitting
p0<-c(-2, -2)

parseparam0<-function(param, detparam=c(log(2),log(1))) {
  pars<-list(obs=c(param[1]),
             proc=c(param[2]),
             det=detparam[c(1:2)])
  return(pars)
}

sampler_fun0<-function(n=1, pars=pars, priorsd=c(1, 1),
                       minv=c(rep(-6.9,2)),
                       maxv=c(rep(2.9,2))) {
  
  d1 = rnorm(n, mean = pars$obs[1], sd = priorsd[1])
  d2 = rnorm(n, mean = pars$proc[1], sd = priorsd[2])
  
  nex<-which((d1<=minv[1]) | (d1>=maxv[1]))
  while(length(nex)>0) {
    d1[nex]<-rnorm(length(nex), mean = pars$obs[1], sd = priorsd[1])
    nex<-which((d1<=minv[1]) | (d1>=maxv[1]))
  }
  
  nex<-which((d2<=minv[2]) | (d2>=maxv[2]))
  while(length(nex)>0) {
    d2[nex]<-rnorm(length(nex), mean = pars$proc[1], sd = priorsd[2])
    nex<-which((d2<=minv[2]) | (d2>=maxv[2]))
  }

  return(cbind(d1,d2))
}

density_fun0<-function(param, pars=pars, priorsd=c(1, 1, 1)){
  dsum = dnorm(param[1], mean = pars$obs[1], sd =  priorsd[1], log = TRUE)
  dsum = dsum+dnorm(param[2], mean = pars$proc[1], sd = priorsd[2], log = TRUE)
  
  return(dsum)
}

density_fun_USE<-function(param) density_fun0(param = param, pars = parseparam0(p0))
sampler_fun_USE<-function(x) sampler_fun0(n = 1, pars = parseparam0(p0), priorsd = c(1,1,1), minv = rep(-4, 2), maxv = rep(0, 2))
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = rep(-4,2), upper = rep(0,2))

#set number of iterations
niter<-5000
N<-2e3

likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0, N = N)
bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior)

likelihood_EDM<-function(x) likelihood0(param = x, y=y, parseparam = parseparam0,
                                        detfun = EDMfun0, edmdat = list(E=2, theta=0.75), N = N)
bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior)

#run MCMC chains
out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0, settings = list(iterations=niter, consoleUpdates=10))
out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=10))
  


#plot outputs
sp<-1000

plot(out_detfun0, start=sp)
plot(out_EDM, start=sp)
gelmanDiagnostics(out_detfun0, plot = FALSE, start=sp)
gelmanDiagnostics(out_EDM, plot = FALSE, start=sp)

## Summarize outputs
smp_detfun0<-getSample(out_detfun0, start = sp)
smp_EDM<-getSample(out_EDM, start=sp)

#check for correlations among estimates
correlationPlot(out_detfun0, start = sp)
correlationPlot(out_EDM, start = sp)

#plot priors and posteriors
marginalPlot(out_detfun0, prior = TRUE, start = sp)
marginalPlot(out_EDM, prior = TRUE, start = sp)

#plot posteriors vs. true values
ptrue<-unlist(pars_true[1:2])

par(mar=c(4,4,2,2), mfrow=c(2,2))
for(i in 1:2) {
  xrng<-range(c(smp_detfun0[,i], ptrue[i], p0[i]))
  hist(smp_detfun0[,i],breaks = 20, probability = TRUE, main="", xlim=xrng);
  abline(v=ptrue[i], col=c(1), lty=2)
  abline(v=p0[i], col=c(3), lty=2)
}
for(i in 1:2) {
  xrng<-range(c(smp_EDM[,i], ptrue[i], p0[i]))
  hist(smp_EDM[,i],breaks = 20, probability = TRUE, main="", xlim=xrng);
  abline(v=ptrue[i], col=c(1), lty=2)
  abline(v=p0[i], col=c(3), lty=2)
}







#gradient
niter<-101

pgrad<-matrix(nrow=niter, ncol=2)
pgrad[,1]<-seq(ptrue[1]-2, ptrue[1]+2, length=niter)
pgrad[,2]<-seq(ptrue[2]-2, ptrue[2]+2, length=niter)

LLmat1<-matrix(nrow=niter, ncol=2)
LLmat2<-matrix(nrow=niter, ncol=2)

for(i in 1:niter) {
  LLmat1[i,1]<-likelihood_detfun0(c(ptrue[1], pgrad[i,2]))
  LLmat1[i,2]<-likelihood_EDM(c(ptrue[1], pgrad[i,2]))
  
  LLmat2[i,1]<-likelihood_detfun0(c(pgrad[i,1], ptrue[2]))
  LLmat2[i,2]<-likelihood_EDM(c(pgrad[i,1], ptrue[2]))
  
  print(round(i/niter,3))
}

par(mar=c(4,4,2,2), mfrow=c(2,2))
plot(exp(pgrad[,1]), exp(LLmat1[,1]), xlab="error", ylab="likelihood", col=4, log="x"); abline(v=exp(ptrue[1]), lty=3)
plot(exp(pgrad[,1]), exp(LLmat1[,2]), xlab="error", ylab="likelihood", col=2, log="x"); abline(v=exp(ptrue[1]), lty=3)

plot(exp(pgrad[,2]), exp(LLmat2[,1]), xlab="error", ylab="likelihood", col=4, log="x"); abline(v=exp(ptrue[2]), lty=3)
plot(exp(pgrad[,2]), exp(LLmat2[,2]), xlab="error", ylab="likelihood", col=2, log="x"); abline(v=exp(ptrue[2]), lty=3)

#estimates
exp(pgrad[which.max(LLmat1[,1]),1])
exp(pgrad[which.max(LLmat1[,2]),1])

exp(pgrad[which.max(LLmat2[,1]),2])
exp(pgrad[which.max(LLmat2[,2]),2])
