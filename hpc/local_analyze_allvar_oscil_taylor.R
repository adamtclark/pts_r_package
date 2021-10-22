rm(list=ls())

set.seed(211005)

setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
require(BayesianTools)
require(mvtnorm)
require(rEDM)
require(msir)
require(pttstability)

#new detfun
pars0<-pars_true<-list(obs=c(log(0.2)),
                       proc=c(log(0.2), log(1.2)),
                       pcol=c(logit(0.2), log(0.1)),
                       det=c(log(1.2),log(1)))

#settings
N<-1e3
nobs<-150
collst<-adjustcolor(c("purple", "blue", "red", "black"),alpha.f = 0.3)

p0<-list(c(log(0.01), log(0.5)), c(log(0.01), log(1)))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(log(0.01), log(0.5)), c(log(0.01), log(1)))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

#load summary data
load("datout/summarydat_211005.rda")

#Plotting and summary functions
e2fun<-function(x,y,ybar=NULL) {
  if(is.null(ybar)) {
    ybar<-mean(y,na.rm=T)
  }
  1-mean((x-y)^2,na.rm=T)/mean((y-ybar)^2,na.rm=T)
}

eifun<-function(x,y,i=1,ybar=NULL) {
  if(is.null(ybar)) {
    ybar<-mean(y,na.rm=T)
  }
  1-mean(abs(x-y)^i,na.rm=T)/mean(abs(y-ybar)^i,na.rm=T)
}

rhokernel<-function(x, y, byvar, nsteps=20, niter=1e3,ei=2) {
  h<-1.06*sd(byvar,na.rm=T)*(length(byvar[is.finite(byvar)])^(-1/5))
  #Silverman, B.W. (1986) "rule of thumb"

  bylst<-seq(min(byvar, na.rm=T), max(byvar, na.rm=T), length=nsteps)
  rholst<-matrix(nrow=length(bylst), ncol=niter)
  eilst<-matrix(nrow=length(bylst), ncol=niter)
  ybar<-mean(y, na.rm=T)

  for(i in 1:length(bylst)) {
    wt<-dnorm((bylst[i]-byvar)/h)
    wt[byvar==bylst[i]]<-0
    wt[is.na(x) | is.na(y)]<-0
    wt[is.na(wt)]<-0
    wt<-wt/sum(wt,na.rm=T)

    for(j in 1:niter) {
      #smp<-sample(length(x), size = length(x)/nsteps, rep=FALSE, prob = wt)
      smp<-sample(length(x), rep=TRUE, prob = wt)
      rholst[i,j]<-cor(x[smp],y[smp])
      eilst[i,j]<-eifun(x[smp],y[smp],ybar,i = ei)
    }
  }

  rhoout<-t(apply(rholst, 1, function(x) quantile(x, pnorm(-1:1))))
  eiout<-t(apply(eilst, 1, function(x) quantile(x, pnorm(-1:1))))

  return(list(rhoout=rhoout, rholst=rholst, eiout=eiout, eilst=eilst, bylst=bylst))
}


rhokernel_2d<-function(x, y, byvar, nsteps=20, niter=1e3, ei=2) {
  byvarold = byvar
  for(i in 1:ncol(byvar)) {
    byvar[,i] = byvar[,i]-mean(byvar[,i], na.rm=T)
    byvar[,i] = byvar[,i]/sd(byvar[,i], na.rm=T)
  }


  h<-diag(ncol(byvar))
  d<-ncol(byvar)
  for(i in 1:d) {
    h[i,i]<-(4/(d+2))^(1/(d+4))*nrow(byvar)^(-1/(d+4))*sd(byvar[,i])
  }
  #Silverman, B.W. (1986) "rule of thumb"

  bylst<-apply(byvar, 2, function(x) seq(min(x, na.rm=T), max(x, na.rm=T), length=nsteps))
  rholst<-array(dim=c(nrow(bylst), nrow(bylst), niter))
  eilst<-array(dim=c(nrow(bylst), nrow(bylst), niter))
  ybar<-mean(y, na.rm=T)

  for(i in 1:nrow(bylst)) {
    for(j in 1:nrow(bylst)) {
      xdiff<-cbind(byvar[,1]-bylst[i,1],
                   byvar[,2]-bylst[j,2])

      wt<-dmvnorm(xdiff, sigma = h)
      wt[is.na(x) | is.na(y)]<-0
      wt[is.na(wt)]<-0
      wt<-wt/sum(wt,na.rm=T)

      for(k in 1:niter) {
        #smp<-sample(length(x), size = length(x)/nsteps, rep=FALSE, prob = wt)
        smp<-sample(length(x), rep=TRUE, prob = wt)

        if(FALSE) {
          eifun(x[smp],y[smp],ybar,i = ei)
          cor(x[smp],y[smp])
          hist(byvar[smp,1])
          hist(byvar[smp,2])
          plot(byvar[smp,], xlim = c(-2,4), ylim=c(-2,4))
          plot(x[smp],y[smp]); abline(a= 0, b = 1)
        }

        rholst[i,j,k]<-cor(x[smp],y[smp])
        eilst[i,j,k]<-eifun(x[smp],y[smp],ybar,i = ei)
      }
    }
  }

  rhoout<-apply(rholst, 1:2, function(x) quantile(x, pnorm(-1:1)))
  eiout<-apply(eilst, 1:2, function(x) quantile(x, pnorm(-1:1)))

  for(i in 1:ncol(byvarold)) {
    bylst[,i] = bylst[,i]*sd(byvarold[,i], na.rm=T)
    bylst[,i] = bylst[,i]+mean(byvarold[,i], na.rm=T)
  }

  return(list(rhoout=rhoout, rholst=rholst, eiout=eiout, eilst=eilst, bylst=bylst))
}


lf<-function(x) {x[is.finite(x) & x<=0]<-min(x[is.finite(x) & x>0], na.rm=T); log10(x)}

pf<-function(x,y,...) {
    points(x,y)
    abline(a=0, b=1, lty=2, col="blue", lwd=1.5)

    ps<-is.finite(x)&is.finite(y)
    mod<-loess.sd(x = x[ps], y = y[ps], nsigma = 1)
    polygon(c(mod$x, rev(mod$x)), c(mod$lower, rev(mod$upper)), col = adjustcolor(1, alpha.f = 0.2), adjustcolor(1,alpha.f = 0.5))
}

pf2<-function(x1,x2,y,category,labels=NULL,rngx=NULL,rngy=NULL,mnlst=NULL,ladj=0,wts1=NULL,wts2=NULL,x3=NULL,wts3=NULL,vline=NULL,...) {
  clevels<-sort(unique(category))
  if(is.null(rngx)) {
    rngx<-range(c(x1,x2,y),na.rm=T)
    if(!is.null(x3)) {
      rngx<-range(c(x1,x2,x3,y),na.rm=T)
    }
  }
  if(is.null(rngy)) {
    rngy<-range(c(x1,x2,y),na.rm=T)
  }

  if(is.null(wts1)) {
    wts1<-rep(1, length(x1))
  }
  if(is.null(wts2)) {
    wts2<-rep(1, length(x2))
  }
  if(is.null(wts3)) {
    wts3<-rep(1, length(x3))
  }

  for(i in 1:length(clevels)) {
    psl<-which(category==clevels[i])
    x1subs<-x1[psl]; x2subs<-x2[psl]; ysubs<-y[psl]
    if(!is.null(x3)) {
      x3subs<-x3[psl]
    }

    if(!is.null(mnlst)) {
      mn<-mnlst[i]
    } else {
      mn<-clevels[i]
    }

    plot(x1subs,ysubs, xlab="", ylab="", xlim=rngx, ylim=rngy, main=mn, type="n")
    points(x1subs,ysubs, col = adjustcolor("dodgerblue", alpha.f = 0.2), pch=16, cex=0.6)
    points(x2subs,ysubs, col = adjustcolor("firebrick", alpha.f = 0.2), pch=17, cex=0.6)

    if(!is.null(x3)) {
      points(x3subs,ysubs, col = adjustcolor("gold", alpha.f = 0.2), pch=18, cex=0.6)
    }

    abline(a=0, b=1, lty=2, col="black", lwd=1.5)

    mod1<-loess.sd(x = x1subs[is.finite(x1subs)], y = ysubs[is.finite(x1subs)], nsigma = 1, weights = 1/wts1[psl][is.finite(x1subs)], enp.target=4)
    mod2<-loess.sd(x = x2subs[is.finite(x2subs)], y = ysubs[is.finite(x2subs)], nsigma = 1, weights = 1/wts2[psl][is.finite(x2subs)], enp.target=4)

    lines(mod1$x, mod1$y, lwd=1.5, col="dodgerblue")
    lines(mod2$x, mod2$y, lwd=1.5, col="firebrick")

    polygon(c(mod1$x, rev(mod1$x)), c(mod1$lower, rev(mod1$upper)), col = adjustcolor("dodgerblue", alpha.f = 0.5), adjustcolor(1,alpha.f = 0.5))
    polygon(c(mod2$x, rev(mod2$x)), c(mod2$lower, rev(mod2$upper)), col = adjustcolor("firebrick", alpha.f = 0.5), adjustcolor(1,alpha.f = 0.5))

    if(!is.null(x3)) {
      mod3<-loess.sd(x = x3subs[is.finite(x3subs)], y = ysubs[is.finite(x3subs)], nsigma = 1, weights = 1/wts3[psl][is.finite(x3subs)], enp.target=4)
      lines(mod3$x, mod3$y, lwd=1.5, col="gold")
      polygon(c(mod3$x, rev(mod3$x)), c(mod3$lower, rev(mod3$upper)), col = adjustcolor("gold", alpha.f = 0.5), adjustcolor(1,alpha.f = 0.5))
    }

    if(!is.null(vline)) {
      abline(v=vline, lty=2, col="black")
    }


    title(paste(letters[i+ladj],".", sep=""), line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)
  }
}

pf3<-function(x,y,minv=-Inf,span=0.75,...) {
  ps<-is.finite(x)&is.finite(y)&(x>minv)&(y>minv)
  mod<-loess.sd(x = x[ps], y = y[ps], nsigma = 1, span=span)

  points(x,y,...)
  polygon(c(mod$x, rev(mod$x)), c(mod$lower, rev(mod$upper)), adjustcolor(1,alpha.f = 0.5), ...)
}

plot_log = function(x,y,minv = 1e-3, axsq = c(0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5), span=0.75, nobs = 150, vline = TRUE, ...) {
  x[x<minv] = minv
  y[y<minv] = minv

  plot(log(x), log(y), axes = FALSE, pch=16, cex=0.8,...)

  axsq[axsq==0] = minv
  axlbs = as.character(axsq)
  axlbs[axlbs==as.character(minv)] = paste("<",minv,sep="")

  axis(1, at = log(axsq), labels = axlbs)
  axis(2, at = log(axsq), labels = axlbs)
  abline(v=log(minv), h=log(minv), lty=3)
  box()

  ps<-is.finite(x)&is.finite(y)&x>minv
  mod<-loess.sd(x = log(x[ps]), y = log(y[ps]), nsigma = 1, span=span)

  polygon(c(mod$x, rev(mod$x)), c(mod$lower, rev(mod$upper)), adjustcolor(1,alpha.f = 0.5), ...)
  if(vline)
    abline(v=log(1/nobs),lty=2)
}

plot_lin = function(x,y,minv = 0, axsq = c(0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5), span=0.75, nobs = 150, vline = TRUE, linmod = FALSE, ...) {
  x[x<minv] = minv
  y[y<minv] = minv

  plot((x), (y), axes = FALSE, pch=16, cex=0.8,...)

  axsq[axsq==0] = minv
  axlbs = as.character(axsq)
  axlbs[axlbs==as.character(minv)] = paste("<",minv,sep="")

  axis(1, at = (axsq), labels = axlbs)
  axis(2, at = (axsq), labels = axlbs)
  abline(v=(minv), h=(minv), lty=3)
  box()

  ps<-is.finite(x)&is.finite(y)&x>minv
  if(linmod) {
    ym = y[ps]; xm = x[ps]
    mod<-lm(ym ~ xm)
    xsq = seq(min(xm), max(xm), length=100)
    pred = predict(mod, newdata=data.frame(xm = xsq), interval = "prediction")

    polygon(c(xsq, rev(xsq)), c(pred[,"lwr"], rev(pred[,"upr"])), adjustcolor(1,alpha.f = 0.5), ...)
  } else {
    mod<-loess.sd(x = (x[ps]), y = (y[ps]), nsigma = 1, span=span)
    polygon(c(mod$x, rev(mod$x)), c(mod$lower, rev(mod$upper)), adjustcolor(1,alpha.f = 0.5), ...)
  }

  if(vline)
    abline(v=(1/nobs),lty=2)
}

#total fit
pdf("plotout/totfit.pdf", width=4, height=5, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,1), mar=c(2,2,1,1), oma=c(2,2,0,0))
  ps<-!is.na(summarydat$gelmandet) & !is.na(summarydat$gelmanedm) & summarydat$gelmandet<1.1 & summarydat$gelmanedm<1.1

  plot(range(summarydat$obs0[ps]), c(0, 1), type="n", xlab="", ylab="")
  pf3(summarydat$obs0[ps], summarydat$cor0[ps], col=adjustcolor("gold", alpha.f = 0.5), pch=16, cex=0.8)
  pf3(summarydat$obs0[ps], summarydat$cor_det[ps], col=adjustcolor("dodgerblue", alpha.f = 0.5), pch=16, cex=0.8)
  pf3(summarydat$obs0[ps], summarydat$cor_edm[ps], col=adjustcolor("firebrick", alpha.f = 0.5), pch=16, cex=0.8)
  abline(h=c(0,1), v=c(0), lty=3)
  mtext(expression(paste("Pearson Correlation, ", rho)), 2, line=2.8)
  title("a.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  plot(range(summarydat$obs0[ps]), c(0, 1), type="n", xlab="", ylab="")
  pf3(summarydat$obs0[ps], (1-summarydat$rmse0/summarydat$proc_true_mu)[ps], col=adjustcolor("gold", alpha.f = 0.5), pch=16, cex=0.8)
  pf3(summarydat$obs0[ps], (1-summarydat$rmse_det/summarydat$proc_true_mu)[ps], col=adjustcolor("dodgerblue", alpha.f = 0.5), pch=16, cex=0.8)
  pf3(summarydat$obs0[ps], (1-summarydat$rmse_edm/summarydat$proc_true_mu)[ps], col=adjustcolor("firebrick", alpha.f = 0.5), pch=16, cex=0.8)
  abline(h=c(0,1), v=c(0), lty=3)
  mtext(expression(paste("Coefficient of Efficiency, ", E[2])), 2, line=2.8)
  title("b.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  mtext(expression(paste("Observation Error, ", sigma[italic(O)])), 1, line=2.8)
dev.off()

#obs error and proc noise
pdf("plotout/obsproc.pdf", width=4, height=7.5, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(3,1), mar=c(4,2,1,1), oma=c(0,3,0,0))

  ps = which(!is.na(summarydat$gelmandet) & summarydat$gelmandet<1.1)
  rng = range(c(summarydat$obs0[ps], summarydat$det_obs_mu0[ps], summarydat$edm_obs_mu0[ps]))
  plot(rng, rng, type="n", xlab="", ylab="")
  pf3(summarydat$det_obs_mu0[ps], summarydat$obs0[ps], col=adjustcolor("dodgerblue", alpha.f = 0.5), pch=16, cex=0.8)
  ps = which(!is.na(summarydat$gelmanedm) & summarydat$gelmanedm<1.1)
  pf3(summarydat$edm_obs_mu0[ps], summarydat$obs0[ps], col=adjustcolor("firebrick", alpha.f = 0.5), pch=16, cex=0.8)
  abline(h=c(0,1), v=c(0), lty=3)
  abline(a=0, b=1, lty=2)
  mtext(expression(paste("True Obs. Error, ", sigma[italic(O)])), 2, line=2.8)
  title("a.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)
  mtext(expression(paste("Estimated Obs. Error, ", hat(sigma)[italic(O)])), 1, line=2.8)

  ps<-which(!is.na(summarydat$gelmandet) & summarydat$gelmandet<1.1 & summarydat$obs0<0.1)
  rng = c(0, 0.8)#range(c(summarydat$proc0[ps], summarydat$det_proc_mu0[ps], summarydat$edm_proc_mu0[ps]))
  plot(rng, rng, type="n", xlab="", ylab="")
  pf3(summarydat$det_proc_mu0[ps], summarydat$proc0[ps], col=adjustcolor("dodgerblue", alpha.f = 0.5), pch=16, cex=0.8)
  ps<-which(!is.na(summarydat$gelmanedm) & summarydat$gelmanedm<1.1 & summarydat$obs0<0.1)
  pf3(summarydat$edm_proc_mu0[ps], summarydat$proc0[ps], col=adjustcolor("firebrick", alpha.f = 0.5), pch=16, cex=0.8)
  ps<-which(summarydat$obs0<0.1)
  pf3(summarydat$proc0_mu[ps], summarydat$proc0[ps], col=adjustcolor("gold", alpha.f = 0.5), pch=16, cex=0.8)

  abline(h=c(0,1), v=c(0), lty=3)
  abline(a=0, b=1, lty=2)
  mtext(expression(paste("True Proc. Noise, ", sigma[italic(P)], " = ", sigma^2,italic(lambda),"/",2,italic(r))), 2, line=0.8, adj = 0.33, outer = TRUE)
  title("b.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)
  mtext(expression(paste("Estimated Proc. Noise, ", hat(sigma)[italic(P)])), 1, line=2.8)
  text(rng[1]+diff(rng)*0.15, rng[2]*.92, expression(paste(sigma[italic(O)], " < 0.1")), cex=1.2)

  ps<-which(!is.na(summarydat$gelmandet) & summarydat$gelmandet<1.1 & summarydat$obs0>0.1)
  #rng = range(c(summarydat$proc0[ps], summarydat$det_proc_mu0[ps], summarydat$edm_proc_mu0[ps]))
  plot(rng, rng, type="n", xlab="", ylab="")
  pf3(summarydat$det_proc_mu0[ps], summarydat$proc0[ps], col=adjustcolor("dodgerblue", alpha.f = 0.5), pch=16, cex=0.8)
  ps<-which(!is.na(summarydat$gelmanedm) & summarydat$gelmanedm<1.1 & summarydat$obs0>0.1)
  pf3(summarydat$edm_proc_mu0[ps], summarydat$proc0[ps], col=adjustcolor("firebrick", alpha.f = 0.5), pch=16, cex=0.8)
  ps<-which(summarydat$obs0>0.1)
  pf3(summarydat$proc0_mu[ps], summarydat$proc0[ps], col=adjustcolor("gold", alpha.f = 0.5), pch=16, cex=0.8)

  abline(h=c(0,1), v=c(0), lty=3)
  abline(a=0, b=1, lty=2)
  #mtext(expression(paste("True Proc. Noise, ", sigma[italic(P)])), 2, line=2.8)
  title("c.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  mtext(expression(paste("Estimated Proc. Noise, ", hat(sigma)[italic(P)])), 1, line=2.8)
  text(rng[1]+diff(rng)*0.15, rng[2]*.92, expression(paste(sigma[italic(O)], " > 0.1")), cex=1.2)
dev.off()



#proc noise by obs error
rhokern_proc0_det<-rhokernel(x=summarydat[summarydat$gelmandet<=1.1,]$det_proc_mu0,
                             y=summarydat[summarydat$gelmandet<=1.1,]$proc0,
                             byvar=summarydat[summarydat$gelmandet<=1.1,]$obs0,
                             nsteps = 20)
rhokern_proc0_edm<-rhokernel(x=summarydat[summarydat$gelmanedm<=1.1,]$edm_proc_mu0,
                             y=summarydat[summarydat$gelmanedm<=1.1,]$proc0,
                             byvar=summarydat[summarydat$gelmanedm<=1.1,]$obs0,
                             nsteps = 20)
pdf("plotout/proc_by_obs.pdf", width=6, height=4, colormodel = "cmyk", useDingbats = FALSE)
  m<-cbind(c(1,1,1), c(1,1,1), c(2,3,4))
  layout(m)
  par(mar=c(2,4,2,1), oma=c(2,0.5,0,0))
  matplot(1,1, xlim=c(0.006, 0.5), ylim=c(0,1.05), xlab="", ylab="", type="n")
  mtext(expression(paste("Goodness of Fit, ", sigma[italic(P)])), 2, line=2.5)
  mtext(expression(paste("Observation Error, ", sigma[italic(O)])), 1, line=2.8)

  polygon(c(rhokern_proc0_det$bylst, rev(rhokern_proc0_det$bylst)),
          c(rhokern_proc0_det$rhoout[,1], rev(rhokern_proc0_det$rhoout[,3])),
          col=adjustcolor("dodgerblue", alpha.f = 0.5), border=NA)
  lines(rhokern_proc0_det$bylst, rhokern_proc0_det$rhoout[,2],
        col="dodgerblue", lwd=1.5, lty=1)
  polygon(c(rhokern_proc0_edm$bylst, rev(rhokern_proc0_edm$bylst)),
          c(rhokern_proc0_edm$rhoout[,1], rev(rhokern_proc0_edm$rhoout[,3])),
          col=adjustcolor("firebrick", alpha.f = 0.5), border=NA)
  lines(rhokern_proc0_edm$bylst, rhokern_proc0_edm$rhoout[,2],
        col="firebrick", lwd=1.5, lty=1)

  polygon(c(rhokern_proc0_det$bylst, rev(rhokern_proc0_det$bylst)),
          c(rhokern_proc0_det$eiout[,1], rev(rhokern_proc0_det$eiout[,3])),
          col=adjustcolor("dodgerblue", alpha.f = 0.5), border=NA)
  lines(rhokern_proc0_det$bylst, rhokern_proc0_det$eiout[,2],
        col="dodgerblue", lwd=1.5, lty=2)
  polygon(c(rhokern_proc0_edm$bylst, rev(rhokern_proc0_edm$bylst)),
          c(rhokern_proc0_edm$eiout[,1], rev(rhokern_proc0_edm$eiout[,3])),
          col=adjustcolor("firebrick", alpha.f = 0.5), border=NA)
  lines(rhokern_proc0_edm$bylst, rhokern_proc0_edm$eiout[,2],
        col="firebrick", lwd=1.5, lty=2)

  legend(0, -0.5,
         c("Analytical Function","EDM Estimate",
           expression(paste("Pearson Correlation, ", rho)),
           expression(paste("Coefficient of Efficiency, ", E[2]))),
         fill = c("dodgerblue", "firebrick", NA, NA), border = c(1, 1, NA, NA),
         lty=c(NA, NA, 1:2), lwd=c(NA, NA, 1.5,1.5), col=c(NA, NA, 1,1), bty="n")

  title("a.", line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)

  abline(h=0, lty=2)
  abline(h=c(-1,1), lty=3)

  ctlvlslin<-(c(0, 0.1, 0.3, 0.5))
  ps<-summarydat$gelmandet<1.1 & summarydat$gelmanedm<1.1
  pf2(x1 = (summarydat$det_proc_mu0[ps]),
      x2 = (summarydat$edm_proc_mu0[ps]),
      y = (summarydat$proc0[ps]),
      wts1=(summarydat$det_proc_qt0.3[ps]-summarydat$det_proc_qt0.2[ps]),
      wts2=(summarydat$edm_proc_qt0.3[ps]-summarydat$edm_proc_qt0.2[ps]),
      category = cut((summarydat[ps,]$obs0),ctlvlslin),
      rngx = (c(0,0.5)), rngy = (c(0,0.5)), ladj = 1,
      mnlst = c(expression(paste(sigma[italic(O)] %in% "(0,0.1]")),
                expression(paste(sigma[italic(O)] %in% "(0.1,0.3]")),
                expression(paste(sigma[italic(O)] %in% "(0.3,0.5]"))))

  mtext(expression(paste("True ", sigma[italic(P)])), 2, outer = TRUE, line=-32)
  mtext(expression(paste("Estimated ", sigma[italic(P)])), 1, line=2.8)
dev.off()









#mortality rate
pdf("plotout/mort.pdf", width=5, height=6.5, colormodel = "cmyk", useDingbats = FALSE)
  rng = c(log(0.001), log(0.5))
  par(mfcol=c(3,2), mar=c(2.2,2,1,1), oma=c(3,3.5,2,0))
  ps = which(!is.na(summarydat$gelmandet) & summarydat$gelmandet<1.1 & summarydat$obs0<0.1)
  plot_log(summarydat$pmdet_analy_noproc[ps], summarydat$pm_actual[ps],
           xlab = "det est", ylab = "true", col=adjustcolor("dodgerblue", alpha.f = 0.5),
           xlim=rng, ylim=rng); abline(a=0, b=1, lty=2)
  mtext(expression(paste(sigma[italic(O)], " > 0.1")), 3, line=0.8)
  title("a.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  ps = which(!is.na(summarydat$gelmanedm) & summarydat$gelmanedm<1.1 & summarydat$obs0<0.1)
  plot_log(summarydat$pmedm_analy_noproc[ps], summarydat$pm_actual[ps], xlab = "edm est", ylab = "true", col=adjustcolor("firebrick", alpha.f = 0.5),
           xlim=rng, ylim=rng); abline(a=0, b=1, lty=2)
  title("b.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  ps = which(summarydat$obs0<0.1)
  plot_log(summarydat$pmobs[ps], summarydat$pm_actual[ps], xlab = "obs est", ylab = "true", col=adjustcolor("gold", alpha.f = 0.5),
           xlim=rng, ylim=rng); abline(a=0, b=1, lty=2)
  title("c.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  ps = which(!is.na(summarydat$gelmandet) & summarydat$gelmandet<1.1 & summarydat$obs0>0.1)
  plot_log(summarydat$pmdet_analy_noproc[ps], summarydat$pm_actual[ps], xlab = "det est", ylab = "true", col=adjustcolor("dodgerblue", alpha.f = 0.5),
           xlim=rng, ylim=rng); abline(a=0, b=1, lty=2)
  mtext(expression(paste(sigma[italic(O)], " > 0.1")), 3, line=0.8)
  title("d.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  ps = which(!is.na(summarydat$gelmanedm) & summarydat$gelmanedm<1.1 & summarydat$obs0>0.1)
  plot_log(summarydat$pmedm_analy_noproc[ps], summarydat$pm_actual[ps], xlab = "edm est", ylab = "true", col=adjustcolor("firebrick", alpha.f = 0.5),
           xlim=rng, ylim=rng); abline(a=0, b=1, lty=2)
  title("e.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  ps = which(summarydat$obs0>0.1)
  plot_log(summarydat$pmobs[ps], summarydat$pm_actual[ps], xlab = "obs est", ylab = "true", col=adjustcolor("gold", alpha.f = 0.5),
           xlim=rng, ylim=rng); abline(a=0, b=1, lty=2)
  title("f.", line=0.2, xpd=NA, adj=0.02, cex.main=1.5)

  mtext(expression(paste("Estimated Mortality Probability, ", hat("Pr")[mor])), 1, line=1.5, outer = TRUE)
  mtext(expression(paste("True Mortality Probability, ", "Pr"[mor])), 2, line=1.5, outer = TRUE)
dev.off()




# total variance and biase
pdf("plotout/totvar.pdf", width=4, height=3.5, colormodel = "cmyk", useDingbats = FALSE)
  par(mfcol=c(1,1),mar=c(4,4,2,2))
  plot_lin((summarydat$edm_var^2-summarydat$obs0^2)[ps], (summarydat$edm_proc_mu0^2)[ps],
           xlim=c(0,0.55), ylim=c(0,0.55),
           axsq = c(seq(0,1,by=0.05)), minv = -1,
           xlab = "", ylab = "", col=adjustcolor("firebrick", alpha.f = 0.5),
           vline = FALSE); abline(a=0, b=1, lty=2)
  mtext(expression(paste("Est. Proc. Noise, ", hat(sigma)[italic(P)["EDM"]]^2)), 2, line=2)
  mtext(expression(paste("Resudual EDM Error, ", epsilon["EDM"]^2-sigma[italic(O)]^2)), 1, line=2.8)
dev.off()


