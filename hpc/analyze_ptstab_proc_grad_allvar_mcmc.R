#!/usr/bin/env Rscript

#error
#rm(list=ls())
#setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")

#load packages and functions
require(BayesianTools)
require(rEDM)
#require(mvtnorm)
require(msir)
source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")

#load data
pars0<-list(obs=(-2.5),
            proc=(-2.5),
            pcol=c((-1.5), NA),
            det=c(log(3),log(1)))
pars0$pcol[2]<-pars0$proc

search_rng<-2.5
prior_sd<-1

#create priors
density_fun_USE<-function(param) density_fun0(param = param, pars = pars0, priorsd = rep(prior_sd, 3))
sampler_fun_USE<-function(x) sampler_fun0(n = 1, pars = pars0, priorsd = rep(prior_sd, 3),
                                          minv = c(pars0$obs, pars0$proc, pars0$pcol[1])-search_rng,
                                          maxv = c(pars0$obs, pars0$proc, pars0$pcol[1])+search_rng)
prior <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
                     lower = c(pars0$obs, pars0$proc, pars0$pcol[1])-search_rng,
                     upper = c(pars0$obs, pars0$proc, pars0$pcol[1])+search_rng)

flst<-dir("datout")
rmps<-grep("_full.rda", flst)
if(length(rmps)>0) {
  flst<-flst[rmps]
}

nvar<-3

qtllst<-c(pnorm(-2:2))
simdatsum<-array(dim=c(length(qtllst),nvar,length(flst),2))
trueparsum<-matrix(nrow=length(flst),ncol=nvar)
colrates<-matrix(nrow=length(flst), ncol=7)
morrates<-matrix(nrow=length(flst), ncol=7)
#simplexdat<-matrix(nrow=length(flst), ncol=3)
rhatdat<-array(dim=c(length(flst),nvar,2))
#colnames(simplexdat)<-c("rho", "mae", "rmse")
rhodat<-matrix(nrow=length(flst), ncol=5)

#if(FALSE) {
  for(ifl in 1:length(flst)) {
    flnm<-paste("datout/", flst[ifl], sep="")
    load(flnm)

    #summarize outputs
    pars<-parslst$ptrue

    truepars_transformed<-inv_fun0(t(as.matrix(unlist(pars))))
    p0_transformed<-inv_fun0(t(as.matrix(unlist(parslst$p0))))

    smp_detfun0_untr<-(getSample(optdat$optout_det))
    smp_EDM_untr<-(getSample(optdat$optout_edm))

    smp_detfun0<-inv_fun0(smp_detfun0_untr)
    smp_EDM<-inv_fun0(smp_EDM_untr)

    simdatsum[,,ifl,1]<-apply(smp_detfun0, 2, function(x) quantile(x, qtllst, na.rm=T))
    simdatsum[,,ifl,2]<-apply(smp_EDM, 2, function(x) quantile(x, qtllst, na.rm=T))

    try(rhatdat[ifl,,1]<-gelmanDiagnostics(optdat$optout_det)$psrf[,1])
    try(rhatdat[ifl,,2]<-gelmanDiagnostics(optdat$optout_edm)$psrf[,1])

    trueparsum[ifl,]<-truepars_transformed[1:nvar]

    #get col/mor
    colrates[ifl,1]<-demdat$cm_obs["pc"]
    colrates[ifl,2]<-demdat$cm_true["pc"]
    colrates[ifl,3]<-demdat$cm_long["pc"]
    colrates[ifl,4]<-demdat$demdat_det["mucol"]
    colrates[ifl,5]<-demdat$demdat_edm["mucol"]
    colrates[ifl,6]<-demdat$demdat_true["mucol"]
    colrates[ifl,7]<-demdat$demdat_edm_true["mucol"]

    morrates[ifl,1]<-demdat$cm_obs["pm"]
    morrates[ifl,2]<-demdat$cm_true["pm"]
    morrates[ifl,3]<-demdat$cm_long["pm"]
    morrates[ifl,4]<-demdat$demdat_det["mumor"]
    morrates[ifl,5]<-demdat$demdat_edm["mumor"]
    morrates[ifl,6]<-demdat$demdat_edm["mumor"]
    morrates[ifl,7]<-demdat$demdat_edm_true["mumor"]

    #smp_out<-unlist(simplex(datout$obs, E = 2, silent = TRUE))
    #simplexdat[ifl,]<-smp_out[6:8]

    #save prediction strengths
    rhodat[ifl,1]<-1-sum((simdat$datout$true-simdat$datout$obs)^2)/sum((simdat$datout$true-mean(simdat$datout$true))^2)#cor(simdat$datout$true, simdat$datout$obs)
    rhodat[ifl,2]<-1-sum((simdat$datout$true-filterdat$filterout_det$rN)^2)/sum((simdat$datout$true-mean(simdat$datout$true))^2)#cor(simdat$datout$true, filterdat$filterout_det$rN)
    rhodat[ifl,3]<-1-sum((simdat$datout$true-filterdat$filterout_edm$rN)^2)/sum((simdat$datout$true-mean(simdat$datout$true))^2)#cor(simdat$datout$true, filterdat$filterout_edm$rN)
    rhodat[ifl,4]<-1-sum((simdat$datout$true-filterdat$filterout_true$rN)^2)/sum((simdat$datout$true-mean(simdat$datout$true))^2)#cor(simdat$datout$true, filterdat$filterout_true$rN)
    #rhodat[ifl,5]<-1-sum((simdat$datout$true-filterdat$filterout_edm_true$rN)^2)/sum((simdat$datout$true-mean(simdat$datout$true))^2)cor(simdat$datout$true, filterdat$filterout_edm_true$rN)

    if(ifl/20 == floor(ifl/20)) {
      print(round(ifl/length(flst),3))
    }
  }

  save.image("summarydata/saved_summary_full_mcmc.rda")
#} else {
#  load("summarydata/saved_summary_full_mcmc.rda")
#}

#check r-hat
hist(rhatdat[,,1])
hist(rhatdat[,,2])


#plot parameters
pltnames<-c("obs", "proc", "col")
collst<-adjustcolor(c("black", "darkgreen", "cornflowerblue", "coral3"), alpha.f = 0.5)
dxl<-c(-0.01, 0.01)

pdf("plotout/plot_pstab_proc_grad_full_mcmc.pdf", width=6, height=12, colormodel = "cmyk", useDingbats = FALSE)
  ##############################
  #plot pars
  ##############################
  uselog<-"xy"
  pssbs<-pssbs<-1:nrow(trueparsum)

  m<-rbind(c(1),
           c(2),
           c(3))
  layout(m)

  par(mar=c(4,4,2,2))
  for(i in 1:length(pltnames)) {
    rngp<-range(c(trueparsum[pssbs,i],simdatsum[2:4,i,pssbs,]),na.rm=T)
    matplot(trueparsum[pssbs,i], simdatsum[3,i,pssbs,],
            col=collst[c(3:4)], type="p", pch=1, cex=0.5,
            ylim=rngp,xlim=rngp,
            xlab="true", ylab="estimated", main=pltnames[i], log=uselog)#, xaxs="i", yaxs="i")
    abline(a=0, b=1, lty=3)
    segments(trueparsum[pssbs,i], simdatsum[2,i,pssbs,1], trueparsum[pssbs,i], simdatsum[4,i,pssbs,1], col=collst[3], lend=2)
    segments(trueparsum[pssbs,i], simdatsum[2,i,pssbs,2], trueparsum[pssbs,i], simdatsum[4,i,pssbs,2], col=collst[4], lend=2)

    x<-trueparsum[pssbs,i]
    y1<-simdatsum[3,i,pssbs,1]
    y2<-simdatsum[3,i,pssbs,2]
    sd1<-simdatsum[3,i,pssbs,1]-simdatsum[2,i,pssbs,1]
    sd2<-simdatsum[3,i,pssbs,2]-simdatsum[2,i,pssbs,2]

    sbsp<-which(is.finite(x) & is.finite(y1) & is.finite(y2) & is.finite(sd1) & is.finite(sd2))
    if(uselog=="xy") {
      x<-log(x[sbsp]); y1<-log(y1[sbsp]); y2<-log(y2[sbsp])
    } else {
      x<-x[sbsp]; y1<-y1[sbsp]; y2<-y2[sbsp]
    }
    sd1<-sd1[sbsp]; sd2<-sd2[sbsp]

    lom1<-loess.sd(y1~x, weights = 1/sd1, enp.target=2, nsigma = 1)
    lom2<-loess.sd(y2~x, weights = 1/sd2, enp.target=2, nsigma = 1)

    px1<-c((sort(lom1$x)), rev((sort(lom1$x))))
    py1<-c((lom1$upper[order(lom1$x)]), rev((lom1$lower[order(lom1$x)])))
    px2<-c((sort(lom2$x)), rev((sort(lom2$x))))
    py2<-c((lom2$upper[order(lom2$x)]), rev((lom2$lower[order(lom2$x)])))

    if(uselog=="xy") {
      px1<-exp(px1)
      py1<-exp(py1)
      px2<-exp(px2)
      py2<-exp(py2)
    }

    polygon(px1,py1,
            col = collst[3])
    polygon(px2,py2,
            col = collst[4])


    abline(h=p0_transformed[i], col=1, lty=3)

    #Add R2
    rssobs<-c(sum((lom1$x-lom1$y)^2*(1/sd1))/sum(1/sd1),
              sum((lom2$x-lom2$y)^2*(1/sd2))/sum(1/sd2))
    rsstot<-mean((lom1$x-mean(lom1$x))^2)
    r2est<-(1-rssobs/rsstot)
    legend("topleft", legend = round(r2est,2), fill = collst[3:4], bty="n", title = expression(paste("E"[2])))

    par(new=TRUE)

    if(uselog=="xy") {
      hist(log(trueparsum[pssbs,i]), xlim=log(rngp), breaks=20, xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(trueparsum)))
    } else {
      hist(trueparsum[pssbs,i], xlim=rngp, breaks=20, xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(trueparsum)))
    }
  }

  ##############################
  #plot rates
  ##############################
  m<-rbind(c(1,4),
           c(2,5),
           c(3,6))
  layout(m)

  collst<-adjustcolor(c("black", "darkgreen", "purple", "cornflowerblue", "coral3", "cornflowerblue", "coral3"), alpha.f = 0.5)

  #col
  x<-trueparsum[,3]
  ptlrng<-range(c(colrates[log(is.finite(colrates))], (trueparsum[,3])),na.rm=T)

  par(mar=c(4,4,2,2))
  for(i in 1:dim(colrates)[2]) {
    sbs<-which(is.finite(log(colrates[,i])))
    if(i%in%c(1,4,6)) {
      plot(ptlrng, ptlrng,
      type="n", xlab="true tcol", ylab="est tcol", log="xy")
      abline(a=0, b=1, lty=3)
    }
    points(trueparsum[sbs,3], colrates[sbs,i], col=collst[i], cex=0.5)

    y<-colrates[sbs,i]
    x1<-x[sbs]
    lss<-loess.sd(log(y)~log(x1), nsigma = 1)
    pd1y<-exp(c(lss$upper, rev(lss$lower)))
    polygon(exp(c(lss$x, rev(lss$x))), pd1y, col = adjustcolor(collst[i]))
  }

  #mor
  x<-morrates[,3]
  ptlrng<-range(c(morrates[is.finite(log(morrates))]),na.rm=T)

  for(i in 1:dim(morrates)[2]) {
    sbs<-which(is.finite(log(morrates[,i])) & is.finite(log(morrates[,3])))
    if(i%in%c(1,4,6)) {
      plot(ptlrng, ptlrng,
           type="n", xlab="true tmor", ylab="est tmor", log="xy")
      abline(a=0, b=1, lty=3)
      abline(v=1/length(simdat$datout$obs), lty=3)
    }
    points(morrates[sbs,3], morrates[sbs,i], col=collst[i], cex=0.5)

    y<-morrates[sbs,i]
    x1<-x[sbs]
    lss<-loess.sd(log(y)~log(x1), nsigma = 1)
    pd1y<-exp(c(lss$upper, rev(lss$lower)))
    polygon(exp(c(lss$x, rev(lss$x))), pd1y, col = adjustcolor(collst[i]))
  }

  par(mfrow=c(1,1))
  matplot(trueparsum[,1], rhodat^2, xlab="obs", ylab="rho", col=collst[1:5], type="p", pch=1)
  x<-trueparsum[,1]
  for(i in 1:ncol(rhodat)) {
    y<-rhodat[,i]^2
    lss<-loess.sd((y)~log(x), nsigma = 1)
    pd1y<-lss$y#(c(lss$y, rev(lss$y)))
    #polygon(exp(c(lss$x, rev(lss$x))), pd1y, col = adjustcolor(collst[i]))
    lines(exp(lss$x), pd1y, col=collst[i], lwd=2)
  }
  legend("bottomleft", c("obs", "det", "edm", "det_true", "edm_true"), lty=1, lwd=2, col=collst, bty="n")
dev.off()


