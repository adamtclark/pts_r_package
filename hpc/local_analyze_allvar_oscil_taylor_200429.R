rm(list=ls())

setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
require(BayesianTools)
require(mvtnorm)
require(rEDM)
require(msir)

source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")

#new detfun
pars0<-pars_true<-list(obs=c(log(0.2)),
                       proc=c(log(0.2), log(1.5)),
                       pcol=c(logit(0.2), log(0.1)),
                       det=c(log(1.2),log(1)))

detfun0_sin<-function(sdet, xt, time=NULL) {
  K<-(sin(time/2)+exp(sdet[2])+0.5)/2
  xt = xt*exp(exp(sdet[1])*(1-xt/K))
  return(xt)
}

#settings
sp<-2000
N<-2e3
collst<-adjustcolor(c("purple", "blue", "red", "black"),alpha.f = 0.3)

p0<-list(c(log(0.01), log(0.5)), c(log(0.01), log(0.5)), c(log(0.5), log(3)))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(log(0.01), log(0.5)), c(log(0.01), log(0.5)), c(log(0.5), log(3)))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

flst<-dir("datout")
#flst<-flst[grep("200429", flst)]
#flst<-flst[grep("200507", flst)]
flst<-flst[grep("200430", flst)]
flst<-flst[grep("rda", flst)]

if(FALSE) {
  summarydat<-data.frame(summed_obs_error=NA,
                         summed_proc_error=NA,
                     obs0=rep(NA, length(flst)),
                     obs1=NA,
                     proc0=NA,
                     proc1=NA,
                     det_obs_mu0=NA,
                     det_obs_mu1=NA,
                     det_proc_mu0=NA,
                     det_proc_mu1=NA,
                     edm_obs_mu0=NA,
                     edm_obs_mu1=NA,
                     edm_proc_mu0=NA,
                     edm_proc_mu1=NA,
                     proc_true_mu=NA,
                     proc0_mu=NA,
                     proc_det_mu=NA,
                     proc_edm_mu=NA,
                     cor0=NA,
                     cor_det=NA,
                     cor_edm=NA,
                     rmse0=NA,
                     rmse_det=NA,
                     rmse_edm=NA,
                     det_obs_qt0=matrix(nrow=length(flst), ncol=5),
                     det_obs_qt1=matrix(nrow=length(flst), ncol=5),
                     det_proc_qt0=matrix(nrow=length(flst), ncol=5),
                     det_proc_qt1=matrix(nrow=length(flst), ncol=5),
                     edm_obs_qt0=matrix(nrow=length(flst), ncol=5),
                     edm_obs_qt1=matrix(nrow=length(flst), ncol=5),
                     edm_proc_qt0=matrix(nrow=length(flst), ncol=5),
                     edm_proc_qt1=matrix(nrow=length(flst), ncol=5),
                     gelmandet=NA,
                     gelmanedm=NA,
                     edm_var=NA,
                     Euse=NA,
                     tuse=NA,
                     pcdet=NA,
                     pcedm=NA,
                     pctrue=NA,
                     pcobs=NA,
                     pmdet=NA,
                     pmedm=NA,
                     pmtrue=NA,
                     pmobs=NA,
                     pmdet_analy=NA,
                     pmedm_analy=NA,
                     pmtrue_analy=NA)


  for(ifl in 1:length(flst)) {
    #ifl<-sample(1:length(flst), 1)
    #ifl<-ifl+1
    #ifl<-6
    load(paste("datout/", flst[ifl], sep=""))
    #exp(parslst$ptrue)
    #exp(parslst$parsest_det[,1])
    #exp(parslst$parsest_edm[,1])


    datout<-simdat$datout
    y<-simdat$datout$obs
    #plot(y)

    out_detfun0<-optdat$optout_det
    out_EDM<-optdat$optout_edm
    pars_true<-parslst$ptrue

    #exp(pars_true)

    Euse<-parslst$Euse
    tuse<-parslst$tuse

    summarydat$summed_obs_error[ifl]<-sd(simdat$datout$true-simdat$datout$obs)
    summarydat$summed_proc_error[ifl]<-sd(simdat$datout$true-simdat$datout$noproc)

    summarydat$obs0[ifl]<-exp(pars_true[1])
    summarydat$obs1[ifl]<-NA
    summarydat$proc0[ifl]<-exp(pars_true[2])
    summarydat$proc1[ifl]<-exp(pars_true[3])

    summarydat$det_obs_mu0[ifl]<-exp(parslst$parsest_det[1,1])
    summarydat$det_obs_mu1[ifl]<-NA
    summarydat$det_proc_mu0[ifl]<-exp(parslst$parsest_det[2,1])
    summarydat$det_proc_mu1[ifl]<-exp(parslst$parsest_det[3,1])
    summarydat$edm_obs_mu0[ifl]<-exp(parslst$parsest_edm[1,1])
    summarydat$edm_obs_mu1[ifl]<-NA
    summarydat$edm_proc_mu0[ifl]<-exp(parslst$parsest_edm[2,1])
    summarydat$edm_proc_mu1[ifl]<-exp(parslst$parsest_edm[3,1])

    tmp<-try(gelmanDiagnostics(out_detfun0, plot = FALSE, start=sp)$psrf[,1], silent = TRUE)
    if(!is.character(tmp)) {
      summarydat$gelmandet[ifl]<-max(tmp)
    }
    tmp<-try(gelmanDiagnostics(out_EDM, plot = FALSE, start=sp)$psrf[,1], silent = TRUE)
    if(!is.character(tmp)) {
      summarydat$gelmanedm[ifl]<-max(tmp)
    }

    if(FALSE) {
      plot(out_detfun0, start=sp)
      plot(out_EDM, start=sp)
      gelmanDiagnostics(out_detfun0, plot = FALSE, start=sp)
      gelmanDiagnostics(out_EDM, plot = FALSE, start=sp)
    }

    ## Summarize outputs
    smp_detfun0<-getSample(out_detfun0, start = sp)
    smp_EDM<-getSample(out_EDM, start=sp)

    #check for correlations among estimates
    if(FALSE) {
      correlationPlot(out_detfun0, start = sp)
      correlationPlot(out_EDM, start = sp)
    }

    #plot priors and posteriors
    if(FALSE) {
      marginalPlot(out_detfun0, prior = TRUE, start = sp)
      marginalPlot(out_EDM, prior = TRUE, start = sp)
    }

    #plot posteriors vs. true values
    ptrue<-unlist(pars_true)

    if(FALSE) {
      par(mar=c(4,4,2,2), mfcol=c(3,2))
      for(i in 1:length(ptrue)) {
        xrng<-exp(range(c(smp_detfun0[,i], ptrue[i], p0[[i]])))
        hist(exp(smp_detfun0[,i]),breaks = 20, probability = TRUE, main="", xlim=xrng);
        abline(v=exp(p0[[i]]), col=c(3), lty=2)
        abline(v=exp(ptrue[i]), col=c(2), lty=2)
      }
      for(i in 1:length(ptrue)) {
        xrng<-exp(range(c(smp_EDM[,i], ptrue[i], p0[[i]])))
        hist(exp(smp_EDM[,i]),breaks = 20, probability = TRUE, main="", xlim=xrng);
        abline(v=exp(p0[[i]]), col=c(3), lty=2)
        abline(v=exp(ptrue[i]), col=c(2), lty=2)
      }
    }

    #run at optimum parameters
    pfout1_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_detfun0)), N=N, detfun=detfun0_sin, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=TRUE, fulltraceback = TRUE)
    pfout2_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_EDM)), N=N, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=Euse, theta=tuse, ytot=y), dotraceback=TRUE, fulltraceback = TRUE)

    summarydat$cor_det[ifl]<-cor(pfout1_opt$Nest, datout$true, use="pairwise")^2
    summarydat$cor_edm[ifl]<-cor(pfout2_opt$Nest, datout$true, use="pairwise")^2
    summarydat$cor0[ifl]<-cor(y, datout$true, use="pairwise")^2

    summarydat$rmse_det[ifl]<-sqrt(mean((pfout1_opt$Nest-datout$true)^2))
    summarydat$rmse_edm[ifl]<-sqrt(mean((pfout2_opt$Nest-datout$true)^2))
    summarydat$rmse0[ifl]<-sqrt(mean((y-datout$true)^2))

    tmp<-exp(apply(smp_detfun0, 2, function(x) quantile(x, pnorm(-2:2))))
    summarydat[ifl,grep("det_obs_qt0", colnames(summarydat))]<-unname(tmp[,1])
    summarydat[ifl,grep("det_proc_qt0", colnames(summarydat))]<-unname(tmp[,2])
    summarydat[ifl,grep("det_proc_qt1", colnames(summarydat))]<-unname(tmp[,3])

    tmp<-exp(apply(smp_EDM, 2, function(x) quantile(x, pnorm(-2:2))))
    summarydat[ifl,grep("edm_obs_qt0", colnames(summarydat))]<-unname(tmp[,1])
    summarydat[ifl,grep("edm_proc_qt0", colnames(summarydat))]<-unname(tmp[,2])
    summarydat[ifl,grep("edm_proc_qt1", colnames(summarydat))]<-unname(tmp[,3])


    summarydat$proc0_mu[ifl]<-sd(datout$obs)
    summarydat$proc_det_mu[ifl]<-sd(pfout1_opt$Nest)
    summarydat$proc_edm_mu[ifl]<-sd(pfout2_opt$Nest)
    summarydat$proc_true_mu[ifl]<-sd(datout$true)

    if(FALSE) {
      par(mfrow=c(1,1))
      plot(datout$true, datout$obs, col="red")
      points(datout$true, pfout2_opt$Nest, col="blue")
      points(datout$true, pfout1_opt$Nest, col="black")
      abline(a=0, b=1, lty=2)

      summarydat$cor_det[ifl]
      summarydat$cor_edm[ifl]
      summarydat$cor0[ifl]
    }

    stmp<-s_map(datout$obs,silent = TRUE,E=Euse, theta=tuse)
    summarydat$edm_var[ifl]<-stmp$rmse
    summarydat$Euse[ifl]<-Euse
    summarydat$tuse[ifl]<-tuse


    #extinction and colonization rates
    prtout_det<-indexsort(pfout1_opt$fulltracemat, pfout1_opt$fulltraceindex, nsmp=1)
    prtout_edm<-indexsort(pfout2_opt$fulltracemat, pfout2_opt$fulltraceindex, nsmp=1)

    prtout_det_noproc<-indexsort(pfout1_opt$fulltracemat_noproc, pfout1_opt$fulltraceindex, nsmp=1)
    prtout_edm_noproc<-indexsort(pfout2_opt$fulltracemat_noproc, pfout2_opt$fulltraceindex, nsmp=1)

    cmdet<-getcm(prtout_det)
    cmedm<-getcm(prtout_edm)
    cmtrue<-getcm(datout$true)
    cmobs<-getcm(datout$obs)

    summarydat$pcdet[ifl]<-cmdet$pc
    summarydat$pcedm[ifl]<-cmedm$pc
    summarydat$pctrue[ifl]<-cmtrue$pc
    summarydat$pcobs[ifl]<-cmobs$pc

    summarydat$pmdet[ifl]<-cmdet$pm
    summarydat$pmedm[ifl]<-cmedm$pm
    summarydat$pmtrue[ifl]<-cmtrue$pm
    summarydat$pmobs[ifl]<-cmobs$pm

    xt<-datout$noproc
    std<-sqrt(exp(parslst$ptrue[2])*xt^exp(parslst$ptrue[3]))
    summarydat$pmtrue_analy[ifl]<-sum(pnorm(0, xt[xt>0], std[xt>0]))/sum(xt>0)

    xtdet<-prtout_det_noproc
    stddet<-sqrt(exp(parslst$parsest_det[2])*xtdet^exp(parslst$parsest_det[3]))
    summarydat$pmdet_analy[ifl]<-sum(pnorm(0, xtdet[!is.na(xtdet) & xtdet>0], stddet[!is.na(xtdet) & xtdet>0]))/sum(!is.na(xtdet) & xtdet>0)

    xtedm<-prtout_edm_noproc
    stdedm<-sqrt(exp(parslst$parsest_edm[2])*xtedm^exp(parslst$parsest_edm[3]))
    summarydat$pmedm_analy[ifl]<-sum(pnorm(0, xtedm[!is.na(xtedm) & xtedm>0], stdedm[!is.na(xtedm) & xtedm>0]))/sum(!is.na(xtedm) & xtedm>0)

    if(ifl/100 == floor(ifl/100)) {
      print(round(ifl/length(flst),2))
    }
  }
  write.csv(summarydat, "datout/summarydat_allvar_oscil_taylor_200430.csv", row.names = FALSE)
} else {
  summarydat<-read.csv("datout/summarydat_allvar_oscil_taylor_200430.csv")
}





#Plot outputs
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

rhokernel<-function(x, y, byvar, nsteps=20, niter=1000,ei=1) {
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
      smp<-sample(length(x), rep=TRUE, prob = wt)
      rholst[i,j]<-cor(x[smp],y[smp])
      eilst[i,j]<-eifun(x[smp],y[smp],ybar,i = ei)
    }
  }

  rhoout<-t(apply(rholst, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
  eiout<-t(apply(eilst, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))

  return(list(rhoout=rhoout, rholst=rholst, eiout=eiout, eilst=eilst, bylst=bylst))
}


rhokernel_2d<-function(x, y, byvar, nsteps=20, niter=1000, ei=1) {
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
        smp<-sample(length(x), rep=TRUE, prob = wt)
        rholst[i,j,k]<-cor(x[smp],y[smp])
        eilst[i,j,k]<-eifun(x[smp],y[smp],ybar,i = ei)
      }
    }
  }

  rhoout<-apply(rholst, 1:2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
  eiout<-apply(eilst, 1:2, function(x) quantile(x, c(0.025, 0.5, 0.975)))

  return(list(rhoout=rhoout, rholst=rholst, eiout=eiout, eilst=eilst, bylst=bylst))
}


lf<-function(x) {x[is.finite(x) & x<=0]<-min(x[is.finite(x) & x>0], na.rm=T); log10(x)}

pf<-function(x,y,...) {
    points(x,y)
    abline(a=0, b=1, lty=2, col="blue", lwd=1.5)

    ps<-is.finite(x)&is.finite(y)
    mod<-loess.sd(x = x[ps], y = y[ps], nsigma = 1)
    polygon(c(mod$x, rev(mod$x)), c(mod$lower, rev(mod$upper)), col = adjustcolor(1, alpha.f = 0.2), border = NA)
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
      points(x3subs,ysubs, col = adjustcolor("goldenrod", alpha.f = 0.2), pch=18, cex=0.6)
    }

    abline(a=0, b=1, lty=2, col="black", lwd=1.5)

    mod1<-loess.sd(x = x1subs[is.finite(x1subs)], y = ysubs[is.finite(x1subs)], nsigma = 1, weights = 1/wts1[psl][is.finite(x1subs)], enp.target=4)
    mod2<-loess.sd(x = x2subs[is.finite(x2subs)], y = ysubs[is.finite(x2subs)], nsigma = 1, weights = 1/wts2[psl][is.finite(x2subs)], enp.target=4)

    lines(mod1$x, mod1$y, lwd=1.5, col="dodgerblue")
    lines(mod2$x, mod2$y, lwd=1.5, col="firebrick")

    polygon(c(mod1$x, rev(mod1$x)), c(mod1$lower, rev(mod1$upper)), col = adjustcolor("dodgerblue", alpha.f = 0.5), border = NA)
    polygon(c(mod2$x, rev(mod2$x)), c(mod2$lower, rev(mod2$upper)), col = adjustcolor("firebrick", alpha.f = 0.5), border = NA)

    if(!is.null(x3)) {
      mod3<-loess.sd(x = x3subs[is.finite(x3subs)], y = ysubs[is.finite(x3subs)], nsigma = 1, weights = 1/wts3[psl][is.finite(x3subs)], enp.target=4)
      lines(mod3$x, mod3$y, lwd=1.5, col="goldenrod")
      polygon(c(mod3$x, rev(mod3$x)), c(mod3$lower, rev(mod3$upper)), col = adjustcolor("goldenrod", alpha.f = 0.5), border = NA)
    }

    if(!is.null(vline)) {
      abline(v=vline, lty=2, col="black")
    }


    title(paste(letters[i+ladj],".", sep=""), line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)
  }
}


summarydat$lvl<-1
ctlvlslin<-(c(0, 0.1, 0.2, 0.5))
ctlvls<-lf(c(1e-6, 0.1, 0.2, 0.5))

### obs
rhokern_obs0_det<-rhokernel(x=summarydat[summarydat$gelmandet<=1.1,]$det_obs_mu0,
         y=summarydat[summarydat$gelmandet<=1.1,]$obs0,
         byvar=summarydat[summarydat$gelmandet<=1.1,]$summed_proc_error,
         nsteps = 20)
rhokern_obs0_edm<-rhokernel(x=summarydat[summarydat$gelmanedm<=1.1,]$edm_obs_mu0,
         y=summarydat[summarydat$gelmanedm<=1.1,]$obs0,
         byvar=summarydat[summarydat$gelmanedm<=1.1,]$summed_proc_error,
         nsteps = 20)

pdf("plotout/local_analyze_allvar_oscil_taylor_200429_obs0.pdf", width=6, height=4, colormodel = "cmyk", useDingbats = FALSE)
  m<-cbind(c(1,1,1), c(1,1,1), c(2,3,4))
  layout(m)
  par(mar=c(2,4,2,1), oma=c(2,0.5,0,0))
  matplot(1,1, xlim=c(0.06, 0.45), ylim=c(-1,1.05), xlab="", ylab="", type="n")
  mtext(expression(paste("Goodness of Fit, ", italic(beta)[obs])), 2, line=2.5)
  mtext(expression(paste("Process Noise, ", sigma[italic(P)[italic(tot)]])), 1, line=2.8)

  polygon(c(rhokern_obs0_det$bylst, rev(rhokern_obs0_det$bylst)),
          c(rhokern_obs0_det$rhoout[,1], rev(rhokern_obs0_det$rhoout[,3])),
          col=adjustcolor("dodgerblue", alpha.f = 0.5), border=NA)
  lines(rhokern_obs0_det$bylst, rhokern_obs0_det$rhoout[,2],
        col="dodgerblue", lwd=1.5, lty=1)
  polygon(c(rhokern_obs0_edm$bylst, rev(rhokern_obs0_edm$bylst)),
          c(rhokern_obs0_edm$rhoout[,1], rev(rhokern_obs0_edm$rhoout[,3])),
          col=adjustcolor("firebrick", alpha.f = 0.5), border=NA)
  lines(rhokern_obs0_edm$bylst, rhokern_obs0_edm$rhoout[,2],
        col="firebrick", lwd=1.5, lty=1)

  polygon(c(rhokern_obs0_det$bylst, rev(rhokern_obs0_det$bylst)),
          c(rhokern_obs0_det$eiout[,1], rev(rhokern_obs0_det$eiout[,3])),
          col=adjustcolor("dodgerblue", alpha.f = 0.5), border=NA)
  lines(rhokern_obs0_det$bylst, rhokern_obs0_det$eiout[,2],
        col="dodgerblue", lwd=1.5, lty=2)
  polygon(c(rhokern_obs0_edm$bylst, rev(rhokern_obs0_edm$bylst)),
          c(rhokern_obs0_edm$eiout[,1], rev(rhokern_obs0_edm$eiout[,3])),
          col=adjustcolor("firebrick", alpha.f = 0.5), border=NA)
  lines(rhokern_obs0_edm$bylst, rhokern_obs0_edm$eiout[,2],
        col="firebrick", lwd=1.5, lty=2)

  legend(0.228, -0.5,
         c("Analytical Function","EDM Estimate",
           expression(paste("Pearson Correlation, ", rho)),
           expression(paste("Coefficient of Efficiency, ", E[1]))),
         fill = c("dodgerblue", "firebrick", NA, NA), border = c(1, 1, NA, NA),
         lty=c(NA, NA, 1:2), lwd=c(NA, NA, 1.5,1.5), col=c(NA, NA, 1,1), bty="n")

  title("a.", line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)

  abline(h=0, lty=2)
  abline(h=c(-1,1), lty=3)


  ps<-summarydat$gelmandet<1.1 & summarydat$gelmanedm<1.1
  pf2(x1 = (summarydat$det_obs_mu0[ps]),
      x2 = (summarydat$edm_obs_mu0[ps]),
      y = (summarydat$obs0[ps]),
      wts1=(summarydat$det_obs_qt0.3[ps]-summarydat$det_obs_qt0.2[ps]),
      wts2=(summarydat$edm_obs_qt0.3[ps]-summarydat$edm_obs_qt0.2[ps]),
      category = cut((summarydat[ps,]$summed_proc_error),ctlvlslin),
      rngx = (c(0,0.5)), rngy = (c(0,0.5)), ladj = 1,
      mnlst = c(expression(paste(sigma[italic(P)[tot]] %in% "(0,0.1]")),
                 expression(paste(sigma[italic(P)[tot]] %in% "(0.1,0.2]")),
                 expression(paste(sigma[italic(P)[tot]] %in% "(0.2,0.5]"))))

  mtext(expression(paste("True ", italic(beta)[obs])), 2, outer = TRUE, line=-32)
  mtext(expression(paste("Predicted ", italic(beta)[obs])), 1, line=2.8)
dev.off()

### proc1
rhokern_proc0_det<-rhokernel(x=summarydat[summarydat$gelmandet<=1.1,]$det_proc_mu0,
                            y=summarydat[summarydat$gelmandet<=1.1,]$proc0,
                            byvar=summarydat[summarydat$gelmandet<=1.1,]$summed_obs_error,
                            nsteps = 20)
rhokern_proc0_edm<-rhokernel(x=summarydat[summarydat$gelmanedm<=1.1,]$edm_proc_mu0,
                            y=summarydat[summarydat$gelmanedm<=1.1,]$proc0,
                            byvar=summarydat[summarydat$gelmanedm<=1.1,]$summed_obs_error,
                            nsteps = 20)

pdf("plotout/local_analyze_allvar_oscil_taylor_200429_proc0.pdf", width=6, height=4, colormodel = "cmyk", useDingbats = FALSE)
  m<-cbind(c(1,1,1), c(1,1,1), c(2,3,4))
  layout(m)
  par(mar=c(2,4,2,1), oma=c(2,0.5,0,0))
  matplot(1,1, xlim=c(0.006, 0.45), ylim=c(-1,1.05), xlab="", ylab="", type="n")
  mtext(expression(paste("Goodness of Fit, ", italic(beta)[proc[0]])), 2, line=2.5)
  mtext(expression(paste("Observation Error, ", sigma[italic(O)[italic(tot)]])), 1, line=2.8)

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

  legend(0.196, -0.5,
         c("Analytical Function","EDM Estimate",
           expression(paste("Pearson Correlation, ", rho)),
           expression(paste("Coefficient of Efficiency, ", E[1]))),
         fill = c("dodgerblue", "firebrick", NA, NA), border = c(1, 1, NA, NA),
         lty=c(NA, NA, 1:2), lwd=c(NA, NA, 1.5,1.5), col=c(NA, NA, 1,1), bty="n")

  title("a.", line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)

  abline(h=0, lty=2)
  abline(h=c(-1,1), lty=3)


  ps<-summarydat$gelmandet<1.1 & summarydat$gelmanedm<1.1
  pf2(x1 = (summarydat$det_proc_mu0[ps]),
      x2 = (summarydat$edm_proc_mu0[ps]),
      y = (summarydat$proc0[ps]),
      wts1=(summarydat$det_proc_qt0.3[ps]-summarydat$det_proc_qt0.2[ps]),
      wts2=(summarydat$edm_proc_qt0.3[ps]-summarydat$edm_proc_qt0.2[ps]),
      category = cut((summarydat[ps,]$summed_obs_error),ctlvlslin),
      rngx = (c(0,0.5)), rngy = (c(0,0.5)), ladj = 1,
      mnlst = c(expression(paste(sigma[italic(O)[tot]] %in% "(0,0.1]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0.1,0.2]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0.2,0.5]"))))

  mtext(expression(paste("True ", italic(beta)[proc[0]])), 2, outer = TRUE, line=-32)
  mtext(expression(paste("Predicted ", italic(beta)[proc[0]])), 1, line=2.8)
dev.off()



#proc1
rhokern_proc1_det<-rhokernel_2d(x=summarydat[summarydat$gelmandet<=1.1,]$det_proc_mu1,
                             y=summarydat[summarydat$gelmandet<=1.1,]$proc1,
                             byvar=cbind(summarydat[summarydat$gelmandet<=1.1,]$summed_obs_error,
                                         summarydat[summarydat$gelmandet<=1.1,]$summed_proc_error))
rhokern_proc1_edm<-rhokernel_2d(x=summarydat[summarydat$gelmanedm<=1.1,]$edm_proc_mu1,
                             y=summarydat[summarydat$gelmanedm<=1.1,]$proc1,
                             byvar=cbind(summarydat[summarydat$gelmanedm<=1.1,]$summed_obs_error,
                                         summarydat[summarydat$gelmanedm<=1.1,]$summed_proc_error))


pdf("plotout/local_analyze_allvar_oscil_taylor_200429_proc1.pdf", width=7, height=4, colormodel = "cmyk", useDingbats = FALSE)
  m<-cbind(c(1,1,1,2,2,2), c(1,1,1,2,2,2),
           c(1,1,1,2,2,2), c(1,1,1,2,2,2),
           c(3,3,3,4,4,4),c(3,3,3,4,4,4),
           c(3,3,3,4,4,4),c(3,3,3,4,4,4),
           c(5,5,5,6,6,6),c(5,5,5,6,6,6),
           c(5,5,5,6,6,6),c(5,5,5,6,6,6))
  layout(m)
  par(mar=c(2,2,2,2), oma=c(2,2,0,0.5))

  contour(rhokern_proc1_edm$bylst[,1], rhokern_proc1_edm$bylst[,2], rhokern_proc1_edm$rhoout[2,,],levels = seq(0, 1, by=0.05), col="firebrick", method="flattest")
  contour(rhokern_proc1_det$bylst[,1], rhokern_proc1_det$bylst[,2], rhokern_proc1_det$rhoout[2,,],levels = seq(0, 1, by=0.05), method="edge", col="dodgerblue", add=TRUE)
  title("a.", line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)
  title(expression(paste("Pearson Correlation, ", rho)))

  contour(rhokern_proc1_edm$bylst[,1], rhokern_proc1_edm$bylst[,2], rhokern_proc1_edm$eiout[2,,],levels = seq(-2, 1, by=0.05), lty=2, col="firebrick", method="flattest")
  contour(rhokern_proc1_det$bylst[,1], rhokern_proc1_det$bylst[,2], rhokern_proc1_det$eiout[2,,],levels = seq(0, 1, by=0.05), lty=2, col="dodgerblue", method="edge", add=TRUE)
  title("b.", line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)
  title(expression(paste("Coefficient of Efficiency, ", E[1])))

  mtext(expression(paste("Process Noise, ", sigma[italic(P)[italic(tot)]])), 2, line=0, outer=TRUE)
  mtext(expression(paste("Observation Error, ", sigma[italic(O)[italic(tot)]])), 1, line=2.8)

  par(mar=c(2,2,2,1))
  ps<-summarydat$gelmandet<1.1 & summarydat$gelmanedm<1.1
  pf2(x1 = (summarydat$det_proc_mu1[ps]),
      x2 = (summarydat$edm_proc_mu1[ps]),
      y = (summarydat$proc1[ps]),
      wts1=(summarydat$det_proc_qt1.3[ps]-summarydat$det_proc_qt1.2[ps]),
      wts2=(summarydat$edm_proc_qt1.3[ps]-summarydat$edm_proc_qt1.2[ps]),
      category = paste(cut((summarydat[ps,]$summed_obs_error),ctlvlslin2), cut((summarydat[ps,]$summed_proc_error),ctlvlslin2)),
      rngx = (c(0.5,3)), rngy = (c(0.5,3)), ladj = 2,
      mnlst = c(expression(paste(sigma[italic(O)[tot]] %in% "(0,0.2], ", sigma[italic(P)[tot]] %in% "(0,0.2]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0,0.2], ", sigma[italic(P)[tot]] %in% "(0.2,0.5]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0.2,0.5], ", sigma[italic(P)[tot]] %in% "(0,0.2]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0.2,0.5], ", sigma[italic(P)[tot]] %in% "(0.2,0.5]"))))

  mtext(expression(paste("True ", italic(beta)[proc[1]])), 2, outer = TRUE, line=-17)
  mtext(expression(paste("Predicted ", italic(beta)[proc[1]])), 1, line=0.9, outer = TRUE, adj = 0.7)
dev.off()





#Analytical times to extinction
libl<-150
rhokern_mor_det<-rhokernel(x=lf(summarydat[summarydat$gelmandet<=1.1,]$pmdet_analy),
                             y=lf(summarydat[summarydat$gelmandet<=1.1,]$pmtrue_analy),
                             byvar=summarydat[summarydat$gelmandet<=1.1,]$summed_obs_error,
                             nsteps = 20)
rhokern_mor_edm<-rhokernel(x=lf(summarydat[summarydat$gelmanedm<=1.1,]$pmedm_analy),
                           y=lf(summarydat[summarydat$gelmanedm<=1.1,]$pmtrue_analy),
                           byvar=summarydat[summarydat$gelmanedm<=1.1,]$summed_obs_error,
                           nsteps = 20)
rhokern_mor_obs<-rhokernel(x=lf(summarydat$pmobs),
                           y=lf(summarydat$pmtrue_analy),
                           byvar=summarydat$summed_obs_error,
                           nsteps = 20)

pdf("plotout/local_analyze_allvar_oscil_taylor_200429_mor.pdf", width=6, height=4, colormodel = "cmyk", useDingbats = FALSE)
  m<-cbind(c(1,1,1), c(1,1,1), c(2,3,4))
  layout(m)

  par(mar=c(2,4,2,1), oma=c(2,0.5,0,0))
  matplot(1,1, xlim=c(0.006, 0.45), ylim=c(-1,1.05), xlab="", ylab="", type="n")
  mtext(expression(paste("Goodness of Fit, ", italic(beta)[proc[0]])), 2, line=2.5)
  mtext(expression(paste("Observation Error, ", sigma[italic(O)[italic(tot)]])), 1, line=2.8)

  polygon(c(rhokern_mor_det$bylst, rev(rhokern_mor_det$bylst)),
          c(rhokern_mor_det$rhoout[,1], rev(rhokern_mor_det$rhoout[,3])),
          col=adjustcolor("dodgerblue", alpha.f = 0.5), border=NA)
  lines(rhokern_mor_det$bylst, rhokern_mor_det$rhoout[,2],
        col="dodgerblue", lwd=1.5, lty=1)

  polygon(c(rhokern_mor_edm$bylst, rev(rhokern_mor_edm$bylst)),
          c(rhokern_mor_edm$rhoout[,1], rev(rhokern_mor_edm$rhoout[,3])),
          col=adjustcolor("firebrick", alpha.f = 0.5), border=NA)
  lines(rhokern_mor_edm$bylst, rhokern_mor_edm$rhoout[,2],
        col="firebrick", lwd=1.5, lty=1)

  polygon(c(rhokern_mor_obs$bylst, rev(rhokern_mor_obs$bylst)),
          c(rhokern_mor_obs$rhoout[,1], rev(rhokern_mor_obs$rhoout[,3])),
          col=adjustcolor("goldenrod", alpha.f = 0.5), border=NA)
  lines(rhokern_mor_obs$bylst, rhokern_mor_obs$rhoout[,2],
        col="goldenrod", lwd=1.5, lty=1)


  polygon(c(rhokern_mor_det$bylst, rev(rhokern_mor_det$bylst)),
          c(rhokern_mor_det$eiout[,1], rev(rhokern_mor_det$eiout[,3])),
          col=adjustcolor("dodgerblue", alpha.f = 0.5), border=NA)
  lines(rhokern_mor_det$bylst, rhokern_mor_det$eiout[,2],
        col="dodgerblue", lwd=1.5, lty=2)

  polygon(c(rhokern_mor_edm$bylst, rev(rhokern_mor_edm$bylst)),
          c(rhokern_mor_edm$eiout[,1], rev(rhokern_mor_edm$eiout[,3])),
          col=adjustcolor("firebrick", alpha.f = 0.5), border=NA)
  lines(rhokern_mor_edm$bylst, rhokern_mor_edm$eiout[,2],
        col="firebrick", lwd=1.5, lty=2)

  polygon(c(rhokern_mor_obs$bylst, rev(rhokern_mor_obs$bylst)),
          c(rhokern_mor_obs$eiout[,1], rev(rhokern_mor_obs$eiout[,3])),
          col=adjustcolor("goldenrod", alpha.f = 0.5), border=NA)
  lines(rhokern_mor_obs$bylst, rhokern_mor_obs$eiout[,2],
        col="goldenrod", lwd=1.5, lty=2)

  abline(h=0, lty=2)
  abline(h=c(-1,1), lty=3)

  legend(0.196, -0.5,
         c("Analytical Function","EDM Estimate", "Raw Observation",
           expression(paste("Pearson Correlation, ", rho)),
           expression(paste("Coefficient of Efficiency, ", E[1]))),
         fill = c("dodgerblue", "firebrick", "goldenrod", NA, NA), border = c(1, 1, 1, NA, NA),
         lty=c(NA, NA, NA, 1:2), lwd=c(NA, NA, NA, 1.5,1.5), col=c(NA, NA, NA, 1,1), bty="n")


  ps<-summarydat$gelmandet<1.1 & summarydat$gelmanedm<1.1
  pf2(x1 = lf(summarydat$pmdet_analy[ps]),
      x2 = lf(summarydat$pmedm_analy[ps]),
      x3 = lf(pmax(summarydat$pmobs[ps], 1/libl)),
      y = lf(summarydat$pmtrue_analy[ps]),
      category = cut((summarydat[ps,]$summed_obs_error),ctlvlslin),
      rngx = (c(-22,-0.5)), rngy = (c(-22,-0.5)), ladj = 2,vline = lf(1/libl),
      mnlst = c(expression(paste(sigma[italic(O)[tot]] %in% "(0,0.1]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0.1,0.2]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0.2,0.5]"))))

  mtext(expression(paste("True ", log[10], "(", Pr[mor], ")")), 2, outer = TRUE, line=-32)
  mtext(expression(paste("Predicted ", log[10], "(", Pr[mor], ")")), 1, line=2.8)
dev.off()





rhokern_col_det<-rhokernel(x=lf(summarydat[summarydat$gelmandet<=1.1 & !is.na(summarydat$pctrue) & summarydat$pctrue>0,]$pcdet),
                           y=lf(summarydat[summarydat$gelmandet<=1.1 & !is.na(summarydat$pctrue) & summarydat$pctrue>0,]$pctrue),
                           byvar=summarydat[summarydat$gelmandet<=1.1 & !is.na(summarydat$pctrue) & summarydat$pctrue>0,]$summed_obs_error,
                           nsteps = 20)
rhokern_col_edm<-rhokernel(x=lf(summarydat[summarydat$gelmanedm<=1.1 & !is.na(summarydat$pctrue) & summarydat$pctrue>0,]$pcedm),
                           y=lf(summarydat[summarydat$gelmanedm<=1.1 & !is.na(summarydat$pctrue) & summarydat$pctrue>0,]$pctrue),
                           byvar=summarydat[summarydat$gelmanedm<=1.1 & !is.na(summarydat$pctrue) & summarydat$pctrue>0,]$summed_obs_error,
                           nsteps = 20)
rhokern_col_obs<-rhokernel(x=lf(summarydat[!is.na(summarydat$pctrue) & summarydat$pctrue>0,]$pcobs),
                           y=lf(summarydat[!is.na(summarydat$pctrue) & summarydat$pctrue>0,]$pctrue),
                           byvar=summarydat[!is.na(summarydat$pctrue) & summarydat$pctrue>0,]$summed_obs_error,
                           nsteps = 20)


pdf("plotout/local_analyze_allvar_oscil_taylor_200429_col.pdf", width=6, height=4, colormodel = "cmyk", useDingbats = FALSE)
  m<-cbind(c(1,1,1), c(1,1,1), c(2,3,4))
  layout(m)

  par(mar=c(2,4,2,1), oma=c(2,0.5,0,0))
  matplot(1,1, xlim=c(0.006, 0.45), ylim=c(-1,1.05), xlab="", ylab="", type="n")
  mtext(expression(paste("Goodness of Fit, ", italic(beta)[proc[0]])), 2, line=2.5)
  mtext(expression(paste("Observation Error, ", sigma[italic(O)[italic(tot)]])), 1, line=2.8)

  polygon(c(rhokern_col_det$bylst, rev(rhokern_col_det$bylst)),
          c(rhokern_col_det$rhoout[,1], rev(rhokern_col_det$rhoout[,3])),
          col=adjustcolor("dodgerblue", alpha.f = 0.5), border=NA)
  lines(rhokern_col_det$bylst, rhokern_col_det$rhoout[,2],
        col="dodgerblue", lwd=1.5, lty=1)

  polygon(c(rhokern_col_edm$bylst, rev(rhokern_col_edm$bylst)),
          c(rhokern_col_edm$rhoout[,1], rev(rhokern_col_edm$rhoout[,3])),
          col=adjustcolor("firebrick", alpha.f = 0.5), border=NA)
  lines(rhokern_col_edm$bylst, rhokern_col_edm$rhoout[,2],
        col="firebrick", lwd=1.5, lty=1)

  polygon(c(rhokern_col_obs$bylst, rev(rhokern_col_obs$bylst)),
          c(rhokern_col_obs$rhoout[,1], rev(rhokern_col_obs$rhoout[,3])),
          col=adjustcolor("goldenrod", alpha.f = 0.5), border=NA)
  lines(rhokern_col_obs$bylst, rhokern_col_obs$rhoout[,2],
        col="goldenrod", lwd=1.5, lty=1)


  polygon(c(rhokern_col_det$bylst, rev(rhokern_col_det$bylst)),
          c(rhokern_col_det$eiout[,1], rev(rhokern_col_det$eiout[,3])),
          col=adjustcolor("dodgerblue", alpha.f = 0.5), border=NA)
  lines(rhokern_col_det$bylst, rhokern_col_det$eiout[,2],
        col="dodgerblue", lwd=1.5, lty=2)

  polygon(c(rhokern_col_edm$bylst, rev(rhokern_col_edm$bylst)),
          c(rhokern_col_edm$eiout[,1], rev(rhokern_col_edm$eiout[,3])),
          col=adjustcolor("firebrick", alpha.f = 0.5), border=NA)
  lines(rhokern_col_edm$bylst, rhokern_col_edm$eiout[,2],
        col="firebrick", lwd=1.5, lty=2)

  polygon(c(rhokern_col_obs$bylst, rev(rhokern_col_obs$bylst)),
          c(rhokern_col_obs$eiout[,1], rev(rhokern_col_obs$eiout[,3])),
          col=adjustcolor("goldenrod", alpha.f = 0.5), border=NA)
  lines(rhokern_col_obs$bylst, rhokern_col_obs$eiout[,2],
        col="goldenrod", lwd=1.5, lty=2)

  abline(h=0, lty=2)
  abline(h=c(-1,1), lty=3)

  legend(0, -0.5,
         c("Analytical Function","EDM Estimate", "Raw Observation",
           expression(paste("Pearson Correlation, ", rho)),
           expression(paste("Coefficient of Efficiency, ", E[1]))),
         fill = c("dodgerblue", "firebrick", "goldenrod", NA, NA), border = c(1, 1, 1, NA, NA),
         lty=c(NA, NA, NA, 1:2), lwd=c(NA, NA, NA, 1.5,1.5), col=c(NA, NA, NA, 1,1), bty="n")


  ps<-summarydat$gelmandet<1.1 & summarydat$gelmanedm<1.1 & !is.na(summarydat$pctrue) & summarydat$pctrue>0
  pf2(x1 = lf(summarydat$pcdet[ps]),
      x2 = lf(summarydat$pcedm[ps]),
      x3 = lf(pmax(summarydat$pcobs[ps], 1/libl)),
      y = lf(summarydat$pctrue[ps]),
      category = cut((summarydat[ps,]$summed_obs_error),ctlvlslin),
      rngx = (c(-1.5,-0)), rngy = (c(-1.5,-0)), ladj = 2, vline = lf(1/libl),
      mnlst = c(expression(paste(sigma[italic(O)[tot]] %in% "(0,0.1]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0.1,0.2]")),
                expression(paste(sigma[italic(O)[tot]] %in% "(0.2,0.5]"))))

  mtext(expression(paste("True ", log[10], "(", Pr[col], ")")), 2, outer = TRUE, line=-32)
  mtext(expression(paste("Predicted ", log[10], "(", Pr[col], ")")), 1, line=2.8)
dev.off()


