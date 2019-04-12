error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")

#load packages and functions
require(BayesianTools); require(rEDM); require(RColorBrewer)
require(msir)
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
rmps<-grep("_full.rda", flst)
if(length(rmps)>0) {
  flst<-flst[rmps]
}

qtllst<-c(0.025, pnorm(-1:1), 0.975)
simdatsum<-array(dim=c(length(qtllst),6,length(flst),2))
trueparsum<-matrix(nrow=length(flst),ncol=6)
colrates<-matrix(nrow=length(flst), ncol=2)
morrates<-matrix(nrow=length(flst), ncol=2)
colrates_q<-array(dim=c(length(flst), 3, 5))
morrates_q<-array(dim=c(length(flst), 3, 5))
simplexdat<-matrix(nrow=length(flst), ncol=3)
rhatdat<-array(dim=c(length(flst),6,2))
colnames(simplexdat)<-c("rho", "mae", "rmse")

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

    try(rhatdat[ifl,,1]<-gelmanDiagnostics(out_detfun0)$psrf[,2]<=1.1)
    try(rhatdat[ifl,,2]<-gelmanDiagnostics(out_EDM)$psrf[,2]<=1.1)

    trueparsum[ifl,]<-truepars_transformed[1:6]

    #get col/mor
    colrates[ifl,1]<-1/demdat$demdat_true$pc
    colrates[ifl,2]<-1/demdat$demdat_short$pc

    colrates_q[ifl,1,]<-quantile(demdat$demdat_filt_true$tcol, qtllst, na.rm=T)
    colrates_q[ifl,2,]<-quantile(demdat$demdat_det$tcol, qtllst, na.rm=T)
    colrates_q[ifl,3,]<-quantile(demdat$demdat_edm$tcol, qtllst, na.rm=T)

    morrates[ifl,1]<-1/demdat$demdat_true$pm
    morrates[ifl,2]<-1/demdat$demdat_short$pm

    morrates_q[ifl,1,]<-quantile(demdat$demdat_filt_true$text, qtllst, na.rm=T)
    morrates_q[ifl,2,]<-quantile(demdat$demdat_det$text, qtllst, na.rm=T)
    morrates_q[ifl,3,]<-quantile(demdat$demdat_edm$text, qtllst, na.rm=T)

    smp_out<-unlist(simplex(datout$obs, E = 2, silent = TRUE))
    simplexdat[ifl,]<-smp_out[6:8]

    if(ifl/20 == floor(ifl/20)) {
      print(round(ifl/length(flst),3))
    }
  }

  save.image("summarydata/saved_summary_full.rda")
} else {
  load("summarydata/saved_summary_full.rda")
}

#mask values that did not converge
for(i in 1:6) {
  simdatsum[,i,!rhatdat[,i,1],1]<-NA
  simdatsum[,i,!rhatdat[,i,2],2]<-NA
}
morrates_q[!rhatdat[,i,1],2,]<-NA
colrates_q[!rhatdat[,i,1],2,]<-NA
morrates_q[!rhatdat[,i,2],3,]<-NA
colrates_q[!rhatdat[,i,2],3,]<-NA

#plot parameters
pltnames<-c("obs b0", "obs b1", "proc b0", "proc b1", "colp", "col0")
collst<-adjustcolor(c("black", "darkgreen", "cornflowerblue", "coral3"), alpha.f = 0.2)
dxl<-c(-0.01, 0.01)

pdf("plotout/plot_pstab_proc_grad_full.pdf", width=12, height=6, colormodel = "cmyk", useDingbats = FALSE)
  m<-rbind(c(1,2,3),
           c(4,5,6))
  layout(m)

  par(mar=c(4,4,2,2))
  for(i in 1:length(pltnames)) {
    if(min(simdatsum[3,i,,],na.rm=T)>0 & min(trueparsum[,i],na.rm=T)>0) {
      uselog<-"xy"
    } else {
      uselog<-""
    }

    if(i %in% c(5,6)) {
      pssbs<-which(is.finite(colrates[,1]) & is.finite(morrates[,1]))
    } else {
      pssbs<-1:nrow(morrates)
    }

    rngp<-range(c(trueparsum[pssbs,i],simdatsum[2:4,i,pssbs,]),na.rm=T)
    matplot(trueparsum[pssbs,i], simdatsum[3,i,pssbs,],
            col=collst[c(3:4)], type="p", pch=1, cex=0.5,
            ylim=rngp,xlim=rngp,
            xlab="true", ylab="estimated", main=pltnames[i], log=uselog, xaxs="i", yaxs="i")
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

    lom1<-loess.sd(y1~x, weights = 1/sd1, nsigma = 1)
    lom2<-loess.sd(y2~x, weights = 1/sd1, nsigma = 1)

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

    #Add R2
    if(uselog=="xy") {
      r2est<-c(1-mean((log(simdatsum[3,i,pssbs,1])-log(trueparsum[pssbs,i]))^2,na.rm=T)/mean((mean(log(trueparsum[pssbs,i]),na.rm=T)-log(trueparsum[pssbs,i]))^2,na.rm=T),
               1-mean((log(simdatsum[3,i,pssbs,2])-log(trueparsum[pssbs,i]))^2,na.rm=T)/mean((mean(log(trueparsum[pssbs,i]),na.rm=T)-log(trueparsum[pssbs,i]))^2,na.rm=T))
      legend("topleft", legend = round(r2est,2), fill = collst[3:4], bty="n", title = expression(paste("R"^2)))
    } else {
      r2est<-c(1-mean(((simdatsum[3,i,pssbs,1])-trueparsum[pssbs,i])^2,na.rm=T)/mean((mean(trueparsum[pssbs,i],na.rm=T)-trueparsum[pssbs,i])^2,na.rm=T),
               1-mean(((simdatsum[3,i,pssbs,2])-trueparsum[pssbs,i])^2,na.rm=T)/mean((mean(trueparsum[pssbs,i],na.rm=T)-trueparsum[pssbs,i])^2,na.rm=T))
      legend("topleft", legend = round(r2est,2), fill = collst[3:4], bty="n", title = expression(paste("R"^2)))
    }


    par(new=TRUE)

    if(uselog=="xy") {
      hist(log(trueparsum[pssbs,i]), xlim=log(rngp), breaks=20, xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(trueparsum)))
    } else {
      hist(trueparsum[pssbs,i], xlim=rngp, breaks=20, xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(trueparsum)))
    }
  }

  #plot rates

  m<-rbind(c(1,2),
           c(3,4))
  layout(m)

  #col
  ptlrng<-range(c(colrates[is.finite(colrates)], colrates_q[,,3][colrates_q[,,3]]),na.rm=T)
  plot(ptlrng, ptlrng,
       type="n", xlab="true tcol", ylab="est tcol", log="xy")
  points(colrates[,1], colrates[,2], col=collst[1], cex=0.5)
  points(colrates[,1], colrates_q[,1,3], col=collst[2], cex=0.5)
  abline(a=0, b=1, lty=3)

  sbs<-which(is.finite(colrates[,2]) & is.finite(colrates[,1]))
  x<-colrates[sbs,1]
  y1<-colrates[sbs,2]

  lom1<-loess.sd(log(y1)~log(x), nsigma = 1)
  xsq<-exp(sort(lom1$x))
  pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))
  polygon(c(xsq, rev(xsq)), pd1y, col = c(adjustcolor(1, alpha.f = 0.3), collst)[1])

  sbs<-which(is.finite(colrates_q[,1,3]) & is.finite(colrates[,1]))
  x<-colrates[sbs,1]
  y1<-colrates_q[sbs,1,3]
  lom1<-loess.sd(log(y1)~log(x), nsigma = 1)
  xsq<-exp(sort(lom1$x))
  pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))
  polygon(c(xsq, rev(xsq)), pd1y, col = collst[2])

  ps<-is.finite(log(colrates[,1])) & is.finite(log(colrates[,2])) & is.finite(log(colrates_q[,1,3]))
  rsstot<-sum((mean(log(colrates[ps,1]),na.rm=T)-log(colrates[ps,1]))^2,na.rm=T)
  r2est<-c(1-sum((log(colrates[ps,2])-log(colrates[ps,1]))^2,na.rm=T)/rsstot,
           1-sum((log(colrates_q[ps,1,3])-log(colrates[ps,1]))^2,na.rm=T)/rsstot)
  legend("topleft", legend = round(r2est,2), fill = collst[1:2], bty="n", title = expression(paste("R"^2)))

  par(new=TRUE)
  hist(log(colrates[,1]), breaks = 20, probability = FALSE, xlim=log(ptlrng), xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(colrates)))


  plot(ptlrng, ptlrng,
       type="n", xlab="true tcol", ylab="est tcol", log="xy")
  points(colrates[,1], colrates_q[,2,3], col=collst[3], cex=0.5)
  points(colrates[,1], colrates_q[,3,3], col=collst[4], cex=0.5)
  abline(a=0, b=1, lty=3)

  for(i in 2:3) {
    sbs<-which(is.finite(colrates_q[,i,3]) & is.finite(colrates[,1]))
    x<-colrates[sbs,1]
    y1<-colrates_q[sbs,i,3]

    lom1<-loess.sd(log(y1)~log(x), nsigma = 1)

    xsq<-exp(sort(lom1$x))
    pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))

    polygon(c(xsq, rev(xsq)), pd1y, col = collst[i+1])
  }

  ps<-is.finite(log(colrates[,1])) & is.finite(log(colrates_q[,2,3])) & is.finite(log(colrates_q[,3,3]))
  rsstot<-sum((mean(log(colrates[ps,1]),na.rm=T)-log(colrates[ps,1]))^2,na.rm=T)
  r2est<-c(1-sum((log(colrates_q[ps,2,3])-log(colrates[ps,1]))^2,na.rm=T)/rsstot,
           1-sum((log(colrates_q[ps,3,3])-log(colrates[ps,1]))^2,na.rm=T)/rsstot)
  legend("topleft", legend = round(r2est,2), fill = collst[3:4], bty="n", title = expression(paste("R"^2)))

  par(new=TRUE)
  hist(log(colrates[,1]), breaks = 20, probability = FALSE, xlim=log(ptlrng), xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(colrates)))


  #mor
  ptlrng<-range(c(morrates[is.finite(morrates)], morrates_q[,,3][morrates_q[,,3]]),na.rm=T)
  ptlrng[2]<-quantile(c(morrates[is.finite(morrates)], morrates_q[,,3][morrates_q[,,3]]), 0.999, na.rm=T)
  plot(ptlrng, ptlrng,
       type="n", xlab="true tmor", ylab="est tmor", log="xy")
  points(morrates[,1], morrates[,2], col=collst[1], cex=0.5)
  points(morrates[,1], morrates_q[,1,3], col=collst[2], cex=0.5)
  abline(a=0, b=1, lty=3)

  sbs<-which(is.finite(morrates[,2]) & is.finite(morrates[,1]))
  x<-morrates[sbs,1]
  y1<-morrates[sbs,2]

  lom1<-loess.sd(log(y1)~log(x), nsigma = 1)
  xsq<-exp(sort(lom1$x))
  pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))
  polygon(c(xsq, rev(xsq)), pd1y, col = c(adjustcolor(1, alpha.f = 0.3), collst)[1])

  sbs<-which(is.finite(morrates_q[,1,3]) & is.finite(morrates[,1]))
  x<-morrates[sbs,1]
  y1<-morrates_q[sbs,1,3]
  lom1<-loess.sd(log(y1)~log(x), nsigma = 1)
  xsq<-exp(sort(lom1$x))
  pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))
  polygon(c(xsq, rev(xsq)), pd1y, col = collst[2])
  abline(v=100, h=100, lty=2)

  ps<-which(is.finite(log(morrates[,1])) & is.finite(log(morrates[,2])) & is.finite(log(morrates_q[,1,3])))
  rsstot<-sum((mean(log(morrates[ps,1]),na.rm=T)-log(morrates[ps,1]))^2,na.rm=T)
  r2est<-c(1-sum((log(morrates[ps,2])-log(morrates[ps,1]))^2,na.rm=T)/rsstot,
           1-sum((log(morrates_q[ps,1,3])-log(morrates[ps,1]))^2,na.rm=T)/rsstot)
  legend("topleft", legend = round(r2est,2), fill = collst[1:2], bty="n", title = expression(paste("R"^2)))

  par(new=TRUE)
  hist(log(morrates[,1]), breaks = 20, probability = FALSE, xlim=log(ptlrng), xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(morrates)))


  plot(ptlrng, ptlrng,
       type="n", xlab="true tmor", ylab="est tmor", log="xy")
  points(morrates[,1], morrates_q[,2,3], col=collst[3], cex=0.5)
  points(morrates[,1], morrates_q[,3,3], col=collst[4], cex=0.5)
  abline(a=0, b=1, lty=3)

  for(i in 2:3) {
    sbs<-which(is.finite(morrates_q[,i,3]) & is.finite(morrates[,1]))
    x<-morrates[sbs,1]
    y1<-morrates_q[sbs,i,3]

    lom1<-loess.sd(log(y1)~log(x), nsigma = 1)

    xsq<-exp(sort(lom1$x))
    pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))

    polygon(c(xsq, rev(xsq)), pd1y, col = collst[i+1])
  }
  abline(v=1000, h=1000, lty=2)

  ps<-which(is.finite(log(morrates[,1])) & is.finite(log(morrates_q[,2,3])) & is.finite(log(morrates_q[,3,3])))
  rsstot<-sum((mean(log(morrates[ps,1]),na.rm=T)-log(morrates[ps,1]))^2,na.rm=T)
  r2est<-c(1-sum((log(morrates_q[ps,2,3])-log(morrates[ps,1]))^2,na.rm=T)/rsstot,
           1-sum((log(morrates_q[ps,3,3])-log(morrates[ps,1]))^2,na.rm=T)/rsstot)
  legend("topleft", legend = round(r2est,2), fill = collst[3:4], bty="n", title = expression(paste("R"^2)))

  par(new=TRUE)
  hist(log(morrates[,1]), breaks = 20, probability = FALSE, xlim=log(ptlrng), xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(morrates)))


dev.off()



  #ERROR PLOTS


  par(mfrow=c(1,3), mar=c(4,4,2,2))
  for(i in 1:6) {
    x<-simplexdat[,3]; y<-sqrt(((t(simdatsum[3,i,,2])-trueparsum[,i])^2))
    plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rho", ylab=paste("rmse", pltnames[i]))
    abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
    lm1<-loess.sd(log(y)~log(x), nsigma = 1)
    matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
    par(new=TRUE)
    hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")
  }


  #Plot simplex error
  par(mfrow=c(1,3), mar=c(4,4,2,2))
  x<-simplexdat[,3]; y<-sqrt(rowSums((t(simdatsum[3,,,2])-trueparsum[,])^2))
  plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rho", ylab="rmse, parameters")
  abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
  lm1<-loess.sd(log(y)~log(x), nsigma = 1)
  matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
  par(new=TRUE)
  hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")

  x<-simplexdat[,3]; y<-sqrt((colrates[,1]-colrates[,4])^2)
  sbs<-which(is.finite(x) & is.finite(y) & y>0) ; x<-x[sbs]; y<-y[sbs]
  plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rho", ylab="rmse, col.")
  abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
  lm1<-loess.sd(log(y)~log(x), nsigma = 1)
  matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
  par(new=TRUE)
  hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")

  x<-simplexdat[,3]; y<-sqrt((morrates[,1]-morrates[,4])^2)
  sbs<-which(is.finite(x) & is.finite(y) & y>0) ; x<-x[sbs]; y<-y[sbs]
  plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rho", ylab="rmse, mor.")
  abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
  lm1<-loess.sd(log(y)~log(x), nsigma = 1)
  matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
  par(new=TRUE)
  hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")

dev.off()



