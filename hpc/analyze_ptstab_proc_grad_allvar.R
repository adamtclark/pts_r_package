#NOTE - need to somehow mark cases were mortality or colonization were never observed


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

    try(rhatdat[ifl,,1]<-gelmanDiagnostics(out_detfun0)$psrf[,1])#<=1.1)
    try(rhatdat[ifl,,2]<-gelmanDiagnostics(out_EDM)$psrf[,1])#<=1.1)

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
  simdatsum[,i,rhatdat[,i,1]>1.1,1]<-NA
  simdatsum[,i,rhatdat[,i,2]>1.1,2]<-NA
}
#morrates_q[!rhatdat[,i,1],2,]<-NA
#colrates_q[!rhatdat[,i,1],2,]<-NA
#morrates_q[!rhatdat[,i,2],3,]<-NA
#colrates_q[!rhatdat[,i,2],3,]<-NA

morrates_q[which(apply(rhatdat[,3:4,1]<=1.1,1,prod)==0),2,]<-NA
colrates_q[which(apply(rhatdat[,5:6,1]<=1.1,1,prod)==0),2,]<-NA
morrates_q[which(apply(rhatdat[,3:4,2]<=1.1,1,prod)==0),3,]<-NA
colrates_q[which(apply(rhatdat[,5:6,2]<=1.1,1,prod)==0),3,]<-NA

#simdatsum[,,trueparsum[,4]>3,]<-NA
#colrates_q[trueparsum[,4]>3,,]<-NA
#colrates_q[trueparsum[,4]>3,,]<-NA


#plot parameters
pltnames<-c("obs b0", "obs b1", "proc b0", "proc b1", "colp", "col0")
collst<-adjustcolor(c("black", "darkgreen", "cornflowerblue", "coral3"), alpha.f = 0.2)
dxl<-c(-0.01, 0.01)

pdf("plotout/plot_pstab_proc_grad_full.pdf", width=12, height=6, colormodel = "cmyk", useDingbats = FALSE)
  ##############################
  #plot pars
  ##############################

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
    #rngp<-range(c(trueparsum[pssbs,i],simdatsum[3,i,pssbs,]),na.rm=T)
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

    lom1<-loess.sd(y1~x, weights = 1/sd1, enp.target=2, nsigma = qnorm(0.975))
    lom2<-loess.sd(y2~x, weights = 1/sd2, enp.target=2, nsigma = qnorm(0.975))

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
    rssobs<-c(sum((lom1$x-lom1$y)^2*(1/sd1))/sum(1/sd1),
              sum((lom2$x-lom2$y)^2*(1/sd2))/sum(1/sd2))
    rsstot<-mean((lom1$x-mean(lom1$x))^2)
    r2est<-(1-rssobs/rsstot)
    legend("topleft", legend = round(r2est,2), fill = collst[3:4], bty="n", title = expression(paste("R"^2)))

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
  m<-rbind(c(1,2),
           c(3,4))
  layout(m)
  cv_cutoff<-0.2


  #col
  sbs<-which(is.finite(colrates[,2]) & is.finite(colrates[,1]) & is.finite(colrates_q[,1,3]) & ((colrates_q[,1,3]-colrates_q[,1,2])/colrates_q[,1,3])<cv_cutoff)

  ptlrng<-range(c(colrates[sbs,][is.finite(colrates[sbs,])], colrates_q[,,3][sbs,][is.finite(colrates_q[,,3][sbs,])]),na.rm=T)
  plot(ptlrng, ptlrng,
       type="n", xlab="true tcol", ylab="est tcol", log="xy")
  points(colrates[sbs,1], colrates[sbs,2], col=collst[1], cex=0.5)
  points(colrates[sbs,1], colrates_q[sbs,1,3], col=collst[2], cex=0.5)
  abline(a=0, b=1, lty=3)

  segments(colrates[sbs,1], colrates_q[sbs,1,2], colrates[sbs,1], colrates_q[sbs,1,4], col=collst[2], lend=2)

  x<-colrates[sbs,1]
  y1<-colrates[sbs,2]

  lom1<-loess.sd(log(y1)~log(x), nsigma = qnorm(0.975), epn.target=2)
  xsq<-exp(sort(lom1$x))
  pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))
  polygon(c(xsq, rev(xsq)), pd1y, col = c(adjustcolor(1, alpha.f = 0.3), collst)[1])

  x<-colrates[sbs,1]
  y1<-colrates_q[sbs,1,3]
  sd2<-(colrates_q[,1,3]-colrates_q[,1,2])[sbs]

  lom2<-loess.sd(log(y1)~log(x), weights=1/sd2, nsigma = qnorm(0.975), epn.target=2)
  xsq<-exp(sort(lom2$x))
  pd1y<-exp(c(lom2$upper[order(lom2$x)], rev(lom2$lower[order(lom2$x)])))
  polygon(c(xsq, rev(xsq)), pd1y, col = collst[2])

  ps<-is.finite(log(colrates[,1])) & is.finite(log(colrates[,2])) & is.finite(log(colrates_q[,1,3]))

  rssobs<-c(mean((lom1$x-lom1$y)^2),
            sum((lom2$x-lom2$y)^2*(1/sd2))/sum(1/sd2))
  rsstot<-mean((lom1$x-mean(lom1$x))^2)
  r2est<-(1-rssobs/rsstot)
  legend("topleft", legend = round(r2est,2), fill = collst[1:2], bty="n", title = expression(paste("R"^2)))

  par(new=TRUE)
  hist(log(colrates[,1]), breaks = 20, probability = FALSE, xlim=log(ptlrng), xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(colrates)))


  sbs<-which(is.finite(colrates_q[,2,3]) & is.finite(colrates_q[,3,3]) & is.finite(colrates[,1]) &
               ((colrates_q[,2,3]-colrates_q[,2,2])/colrates_q[,2,3])<cv_cutoff & ((colrates_q[,3,3]-colrates_q[,3,2])/colrates_q[,2,3])<cv_cutoff)
  plot(ptlrng, ptlrng,
       type="n", xlab="true tcol", ylab="est tcol", log="xy")
  points(colrates[sbs,1], colrates_q[sbs,2,3], col=collst[3], cex=0.5)
  points(colrates[sbs,1], colrates_q[sbs,3,3], col=collst[4], cex=0.5)
  abline(a=0, b=1, lty=3)

  segments(colrates[sbs,1], colrates_q[sbs,2,2], colrates[sbs,1], colrates_q[sbs,2,4], col=collst[3], lend=2)
  segments(colrates[sbs,1], colrates_q[sbs,3,2], colrates[sbs,1], colrates_q[sbs,3,4], col=collst[4], lend=2)

  n<-1
  rs<-numeric(2)
  for(i in 2:3) {
    x<-colrates[sbs,1]
    y1<-colrates_q[sbs,i,3]
    sd1<-colrates_q[sbs,i,3]-colrates_q[sbs,i,2]

    lom1<-loess.sd(log(y1[sd1>0])~log(x[sd1>0]), nsigma = qnorm(0.975), weights=1/sd1[sd1>0], enp.target=2)

    xsq<-exp(sort(lom1$x))
    pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))

    polygon(c(xsq, rev(xsq)), pd1y, col = collst[i+1])
    rs[n]<-sum((lom1$x-lom1$y)^2*(1/sd1[sd1>0]))/sum(1/sd1[sd1>0])
    n<-n+1
  }
  rsstot<-mean((lom1$x-mean(lom1$x))^2)

  r2est<-(1-rs/rsstot)
  legend("topleft", legend = round(r2est,2), fill = collst[3:4], bty="n", title = expression(paste("R"^2)))

  par(new=TRUE)
  hist(log(colrates[,1]), breaks = 20, probability = FALSE, xlim=log(ptlrng), xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(colrates)))


  #mor
  sbs<-which(is.finite(morrates[,2]) & is.finite(morrates[,1]) & is.finite(morrates_q[,1,3]) & ((morrates_q[,1,3]-morrates_q[,1,2])/morrates_q[,1,3])<cv_cutoff &
               morrates[,1]<3500)

  ptlrng<-range(c(morrates[sbs,][is.finite(morrates[sbs,])], morrates_q[,,3][sbs,][morrates_q[,,3][sbs,]]),na.rm=T)
  #ptlrng[2]<-quantile(c(morrates[sbs,][is.finite(morrates[sbs,])], morrates_q[,,3][sbs,][morrates_q[,,3][sbs,]]), 0.999, na.rm=T)
  plot(ptlrng, ptlrng,
       type="n", xlab="true tmor", ylab="est tmor", log="xy")
  points(morrates[sbs,1], morrates[sbs,2], col=collst[1], cex=0.5)
  points(morrates[sbs,1], morrates_q[sbs,1,3], col=collst[2], cex=0.5)
  abline(a=0, b=1, lty=3)

  segments(morrates[sbs,1], morrates_q[sbs,1,2], morrates[sbs,1], morrates_q[sbs,1,4], col=collst[2], lend=2)

  x<-morrates[sbs,1]
  y1<-morrates[sbs,2]
  y1<-y1[order(x)]
  x<-sort(x)

  lom1<-loess.sd(log(y1)~log(x), nsigma = qnorm(0.975), edm.target=2)
  ps2<-(lom1$upper-lom1$lower)>0
  xsq<-exp(sort(lom1$x[ps2]))
  pd1y<-exp(c(lom1$upper[ps2][order(lom1$x[ps2])], rev(lom1$lower[ps2][order(lom1$x[ps2])])))
  pd1y[pd1y>100]<-100
  polygon(c(xsq, rev(xsq)), pd1y, col = c(adjustcolor(1, alpha.f = 0.3), collst)[1])

  sd2<-morrates_q[sbs,1,3]-morrates_q[sbs,1,2]
  x<-morrates[sbs,1]
  y1<-morrates_q[sbs,1,3]
  y1<-y1[order(x)]
  x<-sort(x)

  lom2<-loess.sd(log(y1)~log(x), nsigma = qnorm(0.975), edm.target=2, wts=1/sd2)
  ps2<-(lom2$upper-lom2$lower)>0
  xsq<-exp(sort(lom2$x[ps2]))
  pd1y<-exp(c(lom2$upper[ps2][order(lom2$x[ps2])], rev(lom2$lower[ps2][order(lom2$x[ps2])])))
  polygon(c(xsq, rev(xsq)), pd1y, col = collst[2])
  abline(v=100, h=100, lty=2)

  rssobs<-c(mean((lom1$x-lom1$y)^2),
            sum((lom2$x-lom2$y)^2*(1/sd2))/sum(1/sd2))
  rsstot<-mean((lom1$x-mean(lom1$x))^2)
  r2est<-(1-rssobs/rsstot)
  legend("topleft", legend = round(r2est,2), fill = collst[1:2], bty="n", title = expression(paste("R"^2)))

  par(new=TRUE)
  hist(log(morrates[,1]), breaks = 20, probability = FALSE, xlim=log(ptlrng), xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(morrates)))


  sbs<-which(is.finite(morrates_q[,2,3]) & is.finite(morrates_q[,3,3]) & is.finite(morrates[,1]) &
               ((morrates_q[,2,3]-morrates_q[,2,2])/morrates_q[,2,3])<cv_cutoff & ((morrates_q[,3,3]-morrates_q[,3,2])/morrates_q[,2,3])<cv_cutoff &
               morrates[,1]<2500)

  plot(ptlrng, ptlrng,
       type="n", xlab="true tmor", ylab="est tmor", log="xy")
  points(morrates[sbs,1], morrates_q[sbs,2,3], col=collst[3], cex=0.5)
  points(morrates[sbs,1], morrates_q[sbs,3,3], col=collst[4], cex=0.5)
  abline(a=0, b=1, lty=3)

  segments(morrates[sbs,1], morrates_q[sbs,2,2], morrates[sbs,1], morrates_q[sbs,2,4], col=collst[3], lend=2)
  segments(morrates[sbs,1], morrates_q[sbs,3,2], morrates[sbs,1], morrates_q[sbs,3,4], col=collst[4], lend=2)

  n<-1
  rs<-numeric(2)
  for(i in 2:3) {
    x<-morrates[sbs,1]
    y1<-morrates_q[sbs,i,3]
    sd1<-morrates_q[sbs,i,3]-morrates_q[sbs,i,2]

    lom1<-loess.sd(log(y1[sd1>0])~log(x[sd1>0]), nsigma = qnorm(0.975), weights=1/sd1[sd1>0], enp.target=2)

    xsq<-exp(sort(lom1$x))
    pd1y<-exp(c(lom1$upper[order(lom1$x)], rev(lom1$lower[order(lom1$x)])))

    polygon(c(xsq, rev(xsq)), pd1y, col = collst[i+1])

    rs[n]<-sum((lom1$x-lom1$y)^2*(1/sd1[sd1>0]))/sum(1/sd1[sd1>0])
    n<-n+1
  }
  rsstot<-mean((lom1$x-mean(lom1$x))^2)

  r2est<-(1-rs/rsstot)
  legend("topleft", legend = round(r2est,2), fill = collst[3:4], bty="n", title = expression(paste("R"^2)))

  par(new=TRUE)
  hist(log(morrates[,1]), breaks = 20, probability = FALSE, xlim=log(ptlrng), xlab="", main="", axes=F, ylab="", ylim=c(0, nrow(morrates)))
dev.off()



#ERROR PLOTS
pdf("plotout/plot_pstab_proc_grad_full_SIMPLEXERROR.pdf", width=8, height=6, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,3), mar=c(4,4,2,2))
  for(i in 1:6) {
    x<-c(simplexdat[,3]); y<-c(sqrt(((t(simdatsum[3,i,,2])-trueparsum[,i])^2)))
    plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rmse, EDM", ylab=paste("rmse", pltnames[i]))
    abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
    ps<-is.finite(log(x)) & is.finite(log(y))
    y<-y[ps]; x<-x[ps]
    lm1<-loess.sd(y = log(y), x = log(x), nsigma = 1)
    matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
    par(new=TRUE)
    hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")
  }


  #Plot simplex error
  par(mfrow=c(3,3), mar=c(4,4,2,2))
  x<-simplexdat[,3]; y<-sqrt(rowMeans((t(simdatsum[3,,,2])-trueparsum[,])^2,na.rm = TRUE))
  plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rmse, EDM", ylab="rmse, parameters")
  abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
  ps<-is.finite(log(x)) & is.finite(log(y))
  y<-y[ps]; x<-x[ps]
  lm1<-loess.sd(log(y)~log(x), nsigma = 1)
  matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
  par(new=TRUE)
  hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")

  x<-simplexdat[,3]; y<-sqrt((colrates[,1]-colrates[,2])^2)
  sbs<-which(is.finite(x) & is.finite(y) & y>0) ; x<-x[sbs]; y<-y[sbs]
  plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rmse_EDM", ylab="rmse, col. short")
  abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
  ps<-is.finite(log(x)) & is.finite(log(y))
  y<-y[ps]; x<-x[ps]
  lm1<-loess.sd(log(y)~log(x), nsigma = 1)
  matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
  par(new=TRUE)
  hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")

  x<-simplexdat[,3]; y<-sqrt((morrates[,1]-morrates[,2])^2)
  sbs<-which(is.finite(x) & is.finite(y) & y>0) ; x<-x[sbs]; y<-y[sbs]
  plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rmse_EDM", ylab="rmse, mor. short")
  abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
  ps<-is.finite(log(x)) & is.finite(log(y))
  y<-y[ps]; x<-x[ps]
  lm1<-loess.sd(log(y)~log(x), nsigma = 1)
  matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
  par(new=TRUE)
  hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")

  plttype = c("true filter", "det. filter", "EDM filter")
  for(i in 1:3) {
    x<-simplexdat[,3]; y<-sqrt((colrates_q[,i,3]-colrates[,1])^2)
    sbs<-which(is.finite(x) & is.finite(y) & y>0) ; x<-x[sbs]; y<-y[sbs]
    plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rmse, EDM", ylab=paste("rmse, col.", plttype[i]))
    abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
    ps<-is.finite(log(x)) & is.finite(log(y))
    y<-y[ps]; x<-x[ps]
    lm1<-loess.sd(log(y)~log(x), nsigma = 1)
    matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
    par(new=TRUE)
    hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")
  }
  for(i in 1:3) {
    x<-simplexdat[,3]; y<-sqrt((morrates_q[,i,3]-morrates[,1])^2)
    sbs<-which(is.finite(x) & is.finite(y) & y>0) ; x<-x[sbs]; y<-y[sbs]
    plot(x, y, log="xy", col=adjustcolor(1, alpha.f = 0.5), xlab="rmse, EDM", ylab=paste("rmse, col.", plttype[i]))
    abline(h=10^seq(-12,5), v=c(0.2, 0.5, 1, 2, 5, 10), col="grey")
    ps<-is.finite(log(x)) & is.finite(log(y))
    y<-y[ps]; x<-x[ps]
    lm1<-loess.sd(log(y)~log(x), nsigma = 1)
    matlines(exp(sort(lm1$x)), exp(cbind(lm1$upper, lm1$lower)[order(lm1$x),]), lty=2, col=2)
    par(new=TRUE)
    hist(log10(simplexdat[,3]), axes=F, ylim=c(0, nrow(simplexdat)), main="", xlab="", ylab="")
  }
dev.off()



