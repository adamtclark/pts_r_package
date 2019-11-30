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

#new detfun
detfun0<-function(sdet, xt, time=NULL) {
  K<-(sin(time/2)+exp(sdet[2])+0.5)/2
  xt = xt*exp(exp(sdet[1])*(1-xt/K))
  return(xt)
}

estmorfun<-function(x, psd) {
  tmp<-mean(pnorm(0, x$Nest[pnorm(0, x$Nest, x$Nsd)<0.05], psd),na.rm=T)
  return(tmp)
}

N<-2e3
Euse<-6
tuse<-4
minv<-1e-4

## Simulate data
niter<-1000
if(FALSE) {
  morout<-matrix(nrow=niter, ncol=6)
  colnames(morout)<-c("true", "edm", "det", "true_short", "obs_short", "puse")

  for(itr in 1:niter) {
    lpuse<-runif(1, min = log(0.035), max = log(0.2))
    pars_sim<-pars_true<-list(obs=log(0.1),
                           proc=c(lpuse),
                           pcol=c(logit(0.2), log(0.1)),
                           det=c(log(2),log(1)))

    datout_long<-makedynamics_general(n = 1/minv, n0 = (sin(1/2)+1+0.5)/2, pdet=pars_sim$det,
                                   proc = pars_sim$proc, obs = pars_sim$obs, pcol = pars_sim$pcol,
                                   detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
    datout<-datout_long[1:100,]
    y<-datout$obs

    ## Run filter
    #based on detful0
    filterout_det<-particleFilterLL(y, pars=pars_sim, detfun = detfun0, N = N,
                                    dotraceback = TRUE)
    #based on EDM
    filterout_edm<-particleFilterLL(y, pars=pars_sim, detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse), N = N,
                                    dotraceback = TRUE)

    #estimate mortality rates
    if(FALSE) {
      matplot(filterout_edm$Nest+cbind(-filterout_edm$Nsd*2, 0, filterout_edm$Nsd*2), type="l", lty=c(2,1,2), col=2, xlab="time", ylab="state")
      points(datout$true, col=1, cex=0.5)
      points(y, col=3, cex=0.5)
    }

    morout[itr,]<-c(getcm(datout_long$true)$pm,
                 estmorfun(filterout_edm, exp(pars_sim$proc)),
                 estmorfun(filterout_det, exp(pars_sim$proc)),
                 getcm(datout$true)$pm,
                 getcm(datout$obs)$pm,
                 exp(lpuse))

    if(itr/10==floor(itr/10)) {
      print(round(itr/niter,2))
    }
  }

  write.csv(morout, "datout/morsimout.csv", row.names=F)
  morout<-data.frame(morout)
} else {
  morout<-read.csv("datout/morsimout.csv")
}

collst<-adjustcolor(c("red", "blue", "green", "purple"), alpha.f = 0.8)

plot(morout[,6], morout[,1]+minv, xlab="proc", ylab="mortality", log="xy")
morout$true_plus<-morout$true+minv
pmod<-loess(log(true_plus)~log(puse), morout)
morout$true_pred<-exp(predict(pmod))
lines(sort(morout$puse), morout$true_pred[order(morout$puse)], col=2)

par(mar=c(5.5,5.5,2,2))
axsq<-c(minv, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05)
matplot(morout$true_pred, morout[,2:5]+minv, col=collst, pch=16,
        xlab="", ylab="", log="xy", cex=0.3, axes=F)
axis(1, at=axsq, labels = c(paste("<", axsq[1]), axsq[-1]), las=2)
axis(2, at=axsq, labels = c(paste("<", axsq[1]), axsq[-1]), las=2)
box()
mtext("true", 1, line=4)
mtext("predicted", 2, line=4)

abline(a=0, b=1, lty=3)
psq<-seq(min(morout$true_pred), max(morout$true_pred), length=1000)

for(i in c(2,3,5)) {
  mod<-loess(log(morout[,i]+minv)~log(true_pred), morout, enp.target = 10)
  pred<-predict(mod, newdata=data.frame(true_pred=psq), se=TRUE)
  polygon(c(psq, rev(psq)), exp(c(pred$fit+qt(0.95,pred$df)*pred$se.fit, rev(pred$fit-qt(0.95,pred$df)*pred$se.fit))), col=collst[i-1], border=NA)
}
abline(v=0.01, lty=2)

