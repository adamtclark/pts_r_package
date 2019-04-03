error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")

#load packages and functions
require(BayesianTools); require(rEDM)
source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")

#load data
pars<-list(obs=c(log(1e-2), log(0.1)),
           proc=c(-2, NA),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(log(3),log(1)))

flst<-dir("datout")
rmps<-grep("full", flst)
if(length(rmps)>0) {
  flst<-flst[rmps]
}

qtllst<-c(0.025, pnorm(-1:1), 0.975)
simdatsum<-array(dim=c(length(qtllst),6,length(flst),2))
trueparsum<-matrix(nrow=length(flst),ncol=6)
colrates<-matrix(nrow=length(flst), ncol=4)
morrates<-matrix(nrow=length(flst), ncol=4)

if(FALSE) {
  for(ifl in 1:length(flst)) {
    flnm<-paste("datout/", flst[ifl], sep="")
    load(flnm)

    #summarize outputs
    pars<-pars_sim

    truepars_transformed<-inv_fun0(t(as.matrix(unlist(pars))))

    smp_detfun0_untr<-(getSample(out_detfun0))
    smp_EDM_untr<-(getSample(out_EDM))

    smp_detfun0<-inv_fun0(smp_detfun0_untr)
    smp_EDM<-inv_fun0(smp_EDM_untr)

    simdatsum[,,ifl,1]<-apply(smp_detfun0, 2, function(x) quantile(x, qtllst, na.rm=T))
    simdatsum[,,ifl,2]<-apply(smp_EDM, 2, function(x) quantile(x, qtllst, na.rm=T))

    trueparsum[ifl,]<-truepars_transformed[1:6]

    #get col/mor
    pfdet<-particleFilterLL(y=datout$obs, pars=parseparam0(colMeans(smp_detfun0_untr)), detfun = detfun0, dotraceback = TRUE)
    pfedm<-particleFilterLL(y=datout$obs, pars=parseparam0(colMeans(smp_EDM_untr)), detfun = EDMfun0, edmdat = list(E=2), dotraceback = TRUE)
    pftrue<-particleFilterLL(y=datout$obs, pars=pars, detfun = detfun0, dotraceback = TRUE)

    truecol<-sum(datout$true[-1]>0 & datout$true[-length(datout$true)]==0)/sum(datout$true[-length(datout$true)]==0)
    truemor<-sum(datout$true[-1]==0 & datout$true[-length(datout$true)]>0)/sum(datout$true[-length(datout$true)]>0)

    colrates[ifl,1]<-truecol
    colrates[ifl,2]<-pftrue$dem$mucol
    colrates[ifl,3]<-pfdet$dem$mucol
    colrates[ifl,4]<-pfedm$dem$mucol

    morrates[ifl,1]<-truemor
    morrates[ifl,2]<-pftrue$dem$mumor
    morrates[ifl,3]<-pfdet$dem$mumor
    morrates[ifl,4]<-pfedm$dem$mumor

    if(ifl/20 == floor(ifl/20)) {
      print(round(ifl/length(flst),3))
    }
  }
  save.image("summarydata/saved_summary_full.rda")
} else {
  load("summarydata/saved_summary_full.rda")
}


#plot parameters
pltnames<-c("obs b0", "obs b1", "proc b0", "proc b1", "colp", "col0")
collst<-adjustcolor(c(4,2),alpha.f = 0.3)
dxl<-c(-0.01, 0.01)

pdf("plotout/plot_pstab_proc_grad_full.pdf", width=12, height=8, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,3), mar=c(4,4,2,2))
  for(i in 1:length(pltnames)) {
    if(min(simdatsum[3,i,,])>0 & min(trueparsum[,i])>0) {
      uselog<-"xy"
    } else {
      uselog<-""
    }

    if(i %in% c(5,6)) {
      pssbs<-which(is.finite(colrates[,1]) & is.finite(morrates[,1]))
    } else {
      pssbs<-1:nrow(morrates)
    }

    matplot(trueparsum[pssbs,i], simdatsum[3,i,pssbs,],
            col=collst, type="p", pch=1, cex=0.5,
            ylim=range(simdatsum[3,i,,]),
            xlim=range(simdatsum[3,i,,]),
            xlab="true", ylab="estimated", main=pltnames[i], log=uselog)
    abline(a=0, b=1, lty=3)

    x<-trueparsum[pssbs,i]
    y1<-simdatsum[3,i,pssbs,1]
    y2<-simdatsum[3,i,pssbs,2]
    sd1<-simdatsum[3,i,pssbs,1]-simdatsum[2,i,pssbs,1]
    sd2<-simdatsum[3,i,pssbs,2]-simdatsum[2,i,pssbs,2]

    if(uselog=="xy") {
      x<-log(x); y1<-log(y1); y2<-log(y2)
    }

    lom1<-lm(y1~x, weights = 1/sd1)
    lom2<-lm(y2~x, weights = 1/sd1)

    xsq<-seq(min(x), max(x), length=1000)
    pd1<-predict(lom1, newdata = data.frame(x=xsq), se.fit = TRUE, interval="confidence")
    pd2<-predict(lom2, newdata = data.frame(x=xsq), se.fit = TRUE, interval="confidence")

    pd1y<-c(pd1$fit[,2], rev(pd1$fit[,3]))
    pd2y<-c(pd2$fit[,2], rev(pd2$fit[,3]))

    if(uselog=="xy") {
      pd1y<-exp(pd1y)
      pd2y<-exp(pd2y)
      xsq<-exp(xsq)
    }

    polygon(c(xsq, rev(xsq)), pd1y, border = NA, col = collst[1])
    polygon(c(xsq, rev(xsq)), pd2y, border = NA, col = collst[2])
  }



  #plot rates
  par(mfrow=c(1,2))
  plot(range(colrates,na.rm=T), range(colrates,na.rm=T), type="n", xlab="true col. rate", ylab="est col. rate")
  points(colrates[,1], colrates[,2], col=adjustcolor("black", 0.5), cex=0.5)
  points(colrates[,1], colrates[,3], col=collst[1], cex=0.5)
  points(colrates[,1], colrates[,4], col=collst[2], cex=0.5)
  abline(a=0, b=1, lty=3)

  for(i in 2:4) {
    sbs<-which(is.finite(rowSums(colrates)))
    x<-colrates[sbs,1]
    y1<-colrates[sbs,i]

    lom1<-lm(y1~x)

    xsq<-seq(min(x), max(x), length=1000)
    pd1<-predict(lom1, newdata = data.frame(x=xsq), se.fit = TRUE, interval="confidence")

    pd1y<-c(pd1$fit[,2], rev(pd1$fit[,3]))

    polygon(c(xsq, rev(xsq)), pd1y, border = NA, col = c(adjustcolor(1, alpha.f = 0.3), collst)[i-1])
  }



  plot(range(morrates,na.rm=T), range(morrates,na.rm=T), type="n", xlab="true mor. rate", ylab="est mor. rate")
  points(morrates[,1], morrates[,2], col=adjustcolor("black", 0.3), cex=0.5)
  points(morrates[,1], morrates[,3], col=collst[1], cex=0.5)
  points(morrates[,1], morrates[,4], col=collst[2], cex=0.5)
  abline(a=0, b=1, lty=3)


  for(i in 2:4) {
    sbs<-which(is.finite(rowSums(morrates)))
    x<-morrates[sbs,1]
    y1<-morrates[sbs,i]

    lom1<-lm(y1~x)

    xsq<-seq(min(x), max(x), length=1000)
    pd1<-predict(lom1, newdata = data.frame(x=xsq), se.fit = TRUE, interval="confidence")

    pd1y<-c(pd1$fit[,2], rev(pd1$fit[,3]))

    polygon(c(xsq, rev(xsq)), pd1y, border = NA, col = c(adjustcolor(1, alpha.f = 0.3), collst)[i-1])
  }
dev.off()



