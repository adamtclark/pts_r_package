rm(list=ls())

setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
require(BayesianTools)
require(mvtnorm)
require(rEDM)

source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")
source("plotfun.R")

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
#p0<-list(c(log(0.01), log(1)), c(log(0.01), log(1)), c(log(0.01), log(3)))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(log(0.01), log(0.5)), c(log(0.01), log(0.5)), c(log(0.5), log(3)))
#p0_edm<-list(c(log(0.01), log(1)), c(log(0.01), log(1)), c(log(0.01), log(3)))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

flst<-dir("datout")
#flst<-flst[grep("200429", flst)]
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
    #if(FALSE) {
      pfout1_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_detfun0)), N=N, detfun=detfun0_sin, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=TRUE, fulltraceback = TRUE)
      pfout2_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_EDM)), N=N, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=Euse, theta=tuse, ytot=y), dotraceback=TRUE, fulltraceback = TRUE)
    #}

    #pfout1_opt<-filterdat$filterout_det
    #pfout2_opt<-filterdat$filterout_edm

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
require(msir)
require(mvtnorm)

e2fun<-function(x,y,ybar=NULL) {
  if(is.null(ybar)) {
    ybar<-mean(y,na.rm=T)
  }
  1-mean((x-y)^2,na.rm=T)/mean((y-ybar)^2,na.rm=T)
}


rhokernel<-function(x, y, byvar, nsteps=20, niter=1000) {
  h<-1.06*sd(byvar,na.rm=T)*(length(byvar[is.finite(byvar)])^(-1/5))
  #Silverman, B.W. (1986) "rule of thumb"

  bylst<-seq(min(byvar, na.rm=T), max(byvar, na.rm=T), length=nsteps)
  rholst<-matrix(nrow=length(bylst), ncol=niter)

  for(i in 1:length(bylst)) {
    wt<-dnorm((bylst[i]-byvar)/h)
    wt[byvar==bylst[i]]<-0
    wt[is.na(x) | is.na(y)]<-0
    wt[is.na(wt)]<-0
    wt<-wt/sum(wt,na.rm=T)

    for(j in 1:niter) {
      smp<-sample(length(x), rep=TRUE, prob = wt)
      rholst[i,j]<-cor(x[smp],y[smp])
    }
  }

  rhoout<-t(apply(rholst, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))

  return(list(rhoout=rhoout, rholst=rholst, bylst=bylst))
}


rhokernel_2d<-function(x, y, byvar, nsteps=20, niter=1000) {
  h<-diag(ncol(byvar))
  d<-ncol(byvar)
  for(i in 1:d) {
    h[i,i]<-(4/(d+2))^(1/(d+4))*nrow(byvar)^(-1/(d+4))*sd(byvar[,i])
  }
  #Silverman, B.W. (1986) "rule of thumb"

  bylst<-apply(byvar, 2, function(x) seq(min(x, na.rm=T), max(x, na.rm=T), length=nsteps))
  rholst<-array(dim=c(nrow(bylst), nrow(bylst), niter))

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
      }
    }
  }

  rhoout<-apply(rholst, 1:2, function(x) quantile(x, c(0.025, 0.5, 0.975)))

  return(list(rhoout=rhoout, rholst=rholst, bylst=bylst))
}


lf<-function(x) {log10(x)}
pf<-function(x,y,category,labels=NULL,rngx=NULL,rngy=NULL,mnlst=NULL,...) {
  clevels<-sort(unique(category))
  if(is.null(rngx)) {
    rngx<-range(c(x,y),na.rm=T)
  }
  if(is.null(rngy)) {
    rngy<-range(c(x,y),na.rm=T)
  }
  for(i in 1:length(clevels)) {
    psl<-which(category==clevels[i])
    xsubs<-x[psl]; ysubs<-y[psl]

    if(!is.null(mnlst)) {
      if(mnlst[1]=="none") {
        mn<-""
      } else {
        mn<-mnlst[i]
      }
    } else {
      mn<-clevels[i]
    }

    plot(xsubs,ysubs, xlab="", ylab="", xlim=rngx, ylim=rngy, main=mn)
    abline(a=0, b=1, lty=2, col="blue", lwd=1.5)

    mod<-loess.sd(x = xsubs, y = ysubs, nsigma = qnorm(0.975))
    polygon(c(mod$x, rev(mod$x)), c(mod$lower, rev(mod$upper)), col = adjustcolor(1, alpha.f = 0.2), border = NA)

    e2<-pmax(0, 1-mean((xsubs-ysubs)^2,na.rm=T)/mean((ysubs-mean(ysubs,na.rm=T))^2,na.rm=T))
    r2<-pmax(0, 1-mean((mod$x-mod$y)^2,na.rm=T)/mean((mod$y-mean(mod$y,na.rm=T))^2,na.rm=T))

    text(rngx[1]+diff(rngx)*0.2, rngy[2]-diff(rngy)*0.1, bquote(E[2] == .(round(e2,3))), ps=4)
    text(rngx[1]+diff(rngx)*0.2, rngy[2]-diff(rngy)*0.25, bquote(R^2 == .(round(r2,3))), ps=4)
  }
}

summarydat$lvl<-1
ctlvlslin<-(c(0, 0.1, 0.2, 0.5))
ctlvls<-lf(c(1e-6, 0.1, 0.2, 0.5))

#obs
par(mfcol=c(3,2), mar=c(4,4,2,2))
pf((summarydat[summarydat$gelmandet<=1.1,]$det_obs_mu0),
   (summarydat[summarydat$gelmandet<=1.1,]$obs0),
   cut((summarydat[summarydat$gelmandet<=1.1,]$summed_proc_error),ctlvlslin),
   rngx = c(0,0.5), rngy = c(0,0.5))

pf((summarydat[summarydat$gelmanedm<=1.1,]$edm_obs_mu0),
   (summarydat[summarydat$gelmanedm<=1.1,]$obs0),
   cut((summarydat[summarydat$gelmanedm<=1.1,]$summed_proc_error),ctlvlslin),
   rngx = c(0,0.5), rngy = c(0,0.5))

rhokern_obs0_det<-rhokernel(x=summarydat[summarydat$gelmandet<=1.1,]$det_obs_mu0,
         y=summarydat[summarydat$gelmandet<=1.1,]$obs0,
         byvar=summarydat[summarydat$gelmandet<=1.1,]$summed_proc_error)
rhokern_obs0_edm<-rhokernel(x=summarydat[summarydat$gelmanedm<=1.1,]$edm_obs_mu0,
                          y=summarydat[summarydat$gelmanedm<=1.1,]$obs0,
                          byvar=summarydat[summarydat$gelmanedm<=1.1,]$summed_proc_error)
matplot(rhokern_obs0_det$bylst, rhokern_obs0_det$rhoout, type="l",
        lty=c(2,1,2), col="blue", xlim=c(0, 0.5), ylim=c(0,1),
        xlab="Observation Error", ylab="Pearson Correlation")
matlines(rhokern_obs0_edm$bylst, rhokern_obs0_edm$rhoout,
        lty=c(2,1,2), col="red")
abline(h=0, lty=2)
abline(h=c(-1,1), lty=3)



#proc0
par(mfcol=c(3,2), mar=c(4,4,2,2))
pf((summarydat[summarydat$gelmandet<=1.1,]$det_proc_mu0),
   (summarydat[summarydat$gelmandet<=1.1,]$proc0),
   cut((summarydat[summarydat$gelmandet<=1.1,]$summed_obs_error),ctlvlslin),
   rngx = c(0,0.5), rngy = c(0,0.5))

pf((summarydat[summarydat$gelmanedm<=1.1,]$edm_proc_mu0),
   (summarydat[summarydat$gelmanedm<=1.1,]$proc0),
   cut((summarydat[summarydat$gelmanedm<=1.1,]$summed_obs_error),ctlvlslin),
   rngx = c(0,0.5), rngy = c(0,0.5))

rhokern_proc0_det<-rhokernel(x=summarydat[summarydat$gelmandet<=1.1,]$det_proc_mu0,
                          y=summarydat[summarydat$gelmandet<=1.1,]$proc0,
                          byvar=summarydat[summarydat$gelmandet<=1.1,]$summed_obs_error)
rhokern_proc0_edm<-rhokernel(x=summarydat[summarydat$gelmanedm<=1.1,]$edm_proc_mu0,
                          y=summarydat[summarydat$gelmanedm<=1.1,]$proc0,
                          byvar=summarydat[summarydat$gelmanedm<=1.1,]$summed_obs_error)
matplot(rhokern_proc0_det$bylst, rhokern_proc0_det$rhoout, type="l",
        lty=c(2,1,2), col="blue", xlim=c(0, 0.5), ylim=c(0,1),
        xlab="Process Noise", ylab="Pearson Correlation")
matlines(rhokern_proc0_edm$bylst, rhokern_proc0_edm$rhoout,
         lty=c(2,1,2), col="red")
abline(h=0, lty=2)
abline(h=c(-1,1), lty=3)

#proc1
ctlvlslin2<-(c(0, 0.2, 0.5))

par(mfcol=c(2,2), mar=c(4,4,2,2))
pf((summarydat[summarydat$gelmandet<=1.1,]$det_proc_mu1),
   (summarydat[summarydat$gelmandet<=1.1,]$proc1),
   paste(cut((summarydat[summarydat$gelmandet<=1.1,]$summed_obs_error),ctlvlslin2),
         cut((summarydat[summarydat$gelmandet<=1.1,]$summed_proc_error),ctlvlslin2)),
   rngx = c(0.5,3), rngy = c(0.5,3))

pf((summarydat[summarydat$gelmanedm<=1.1,]$edm_proc_mu1),
   (summarydat[summarydat$gelmanedm<=1.1,]$proc1),
   paste(cut((summarydat[summarydat$gelmanedm<=1.1,]$summed_obs_error),ctlvlslin2),
         cut((summarydat[summarydat$gelmanedm<=1.1,]$summed_proc_error),ctlvlslin2)),
   rngx = c(0.5,3), rngy = c(0.5,3))


rhokern_proc1_det<-rhokernel_2d(x=summarydat[summarydat$gelmandet<=1.1,]$det_proc_mu1,
                             y=summarydat[summarydat$gelmandet<=1.1,]$proc1,
                             byvar=cbind(summarydat[summarydat$gelmandet<=1.1,]$summed_obs_error,
                                         summarydat[summarydat$gelmandet<=1.1,]$summed_proc_error))
rhokern_proc1_edm<-rhokernel_2d(x=summarydat[summarydat$gelmanedm<=1.1,]$edm_proc_mu1,
                             y=summarydat[summarydat$gelmanedm<=1.1,]$proc1,
                             byvar=cbind(summarydat[summarydat$gelmanedm<=1.1,]$summed_obs_error,
                                         summarydat[summarydat$gelmanedm<=1.1,]$summed_proc_error))
par(mfrow=c(2,1), mar=c(4,4,2,2))
contour(rhokern_proc1_det$bylst[,1], rhokern_proc1_det$bylst[,2], rhokern_proc1_det$rhoout[1,,],levels = seq(0, 1, by=0.02))
contour(rhokern_proc1_edm$bylst[,1], rhokern_proc1_edm$bylst[,2], rhokern_proc1_edm$rhoout[1,,],levels = seq(0, 1, by=0.02))


#Check rates, mor
libl<-150
coplot(lf(pmax(pmtrue, 1/libl))~lf(pmax(pmdet, 1/libl))|lvl, summarydat[summarydat$gelmandet<=1.1,])
coplot(lf(pmax(pmtrue, 1/libl))~lf(pmax(pmdet, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmandet<=1.1,])

coplot(lf(pmax(pmtrue, 1/libl))~lf(pmax(pmedm, 1/libl))|lvl, summarydat[summarydat$gelmanedm<=1.1,])
coplot(lf(pmax(pmtrue, 1/libl))~lf(pmax(pmedm, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmanedm<=1.1,])

coplot(lf(pmax(pmtrue, 1/libl))~lf(pmax(pmobs, 1/libl))|lvl, summarydat, panel=pf)
coplot(lf(pmax(pmtrue, 1/libl))~lf(pmax(pmobs, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat, panel=pf)




#Check rates, col
coplot(lf(pmax(pctrue, 1/libl))~lf(pmax(pcdet, 1/libl))|lvl, summarydat[summarydat$gelmandet<=1.1 & summarydat$pmtrue>0,], panel=pf)
coplot(lf(pmax(pctrue, 1/libl))~lf(pmax(pcdet, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmandet<=1.1 & summarydat$pmtrue>0,], panel=pf)

coplot(lf(pmax(pctrue, 1/libl))~lf(pmax(pcedm, 1/libl))|lvl, summarydat[summarydat$gelmanedm<=1.1 & summarydat$pmtrue>0,], panel=pf)
coplot(lf(pmax(pctrue, 1/libl))~lf(pmax(pcedm, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmanedm<=1.1 & summarydat$pmtrue>0,], panel=pf)

coplot(lf(pmax(pctrue, 1/libl))~lf(pmax(pcobs, 1/libl))|lvl, summarydat[summarydat$pmtrue>0,], panel=pf)
coplot(lf(pmax(pctrue, 1/libl))~lf(pmax(pcobs, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$pmtrue>0,], panel=pf)




#Analytical times to extinction
coplot(lf(pmtrue_analy)~lf(pmdet_analy)|lvl, summarydat[summarydat$gelmandet<=1.1,], panel=pf)
coplot(lf(pmtrue_analy)~lf(pmdet_analy)|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmandet<=1.1,], panel=pf)

coplot(lf(pmtrue_analy)~lf(pmedm_analy)|lvl, summarydat[summarydat$gelmanedm<=1.1,], panel=pf)
coplot(lf(pmtrue_analy)~lf(pmedm_analy)|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmanedm<=1.1,], panel=pf)

coplot(lf(pmtrue_analy)~lf(pmax(pmobs, 1/libl))|lvl, summarydat, panel=pf)
coplot(lf(pmtrue_analy)~lf(pmax(pmobs, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat, panel=pf)

coplot(lf(pmtrue_analy)~lf(pmdet_analy)|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmandet<=1.1 & summarydat$pmtrue_analy<(1/libl),], panel=pf)
coplot(lf(pmtrue_analy)~lf(pmedm_analy)|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmanedm<=1.1 & summarydat$pmtrue_analy<(1/libl),], panel=pf)
coplot(lf(pmtrue_analy)~lf(pmax(pmobs, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$pmtrue_analy<(1/libl),], panel=pf)


coplot(lf(pmtrue_analy)~lf(pmdet_analy)|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmandet<=1.1 & summarydat$pmtrue_analy>(1/libl),], panel=pf)
coplot(lf(pmtrue_analy)~lf(pmedm_analy)|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$gelmanedm<=1.1 & summarydat$pmtrue_analy>(1/libl),], panel=pf)
coplot(lf(pmtrue_analy)~lf(pmax(pmobs, 1/libl))|cut(lf(summed_obs_error),ctlvls), summarydat[summarydat$pmtrue_analy>(1/libl),], panel=pf)


