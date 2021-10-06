#!/usr/bin/env Rscript
rm(list=ls())

if(FALSE) {
  setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
  source("../pttstability/R/bayesfun.R")
  source("../pttstability/R/fake_data.R")
  source("../pttstability/R/logit_funs.R")
  source("../pttstability/R/particlefilter.R")
}

if(length(dir("/cl_tmp/clarka/Rpkg/"))>0) {
  .libPaths("/cl_tmp/clarka/Rpkg/")
}

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
sp<-1000
N<-1e3
nobs<-150
collst<-adjustcolor(c("purple", "blue", "red", "black"),alpha.f = 0.3)

p0<-list(c(log(0.01), log(0.5)), c(log(0.01), log(1)))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(log(0.01), log(0.5)), c(log(0.01), log(1)))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

flst<-dir("datout")
flst<-flst[grep("211005", flst)]
flst<-flst[grep("rda", flst)]

zero_cutoff = 1e-3


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
                       pc_actual=NA,
                       pm_actual=NA,
                       pmdet_analy=NA,
                       pmedm_analy=NA,
                       pmtrue_analy=NA)


for(ifl in 1:length(flst)) {
  load(paste("datout/", flst[ifl], sep=""))
  #c(exp(parslst$ptrue[1]), parslst$sd_abs)
  #exp(parslst$parsest_det[,1])
  #exp(parslst$parsest_edm[,1])


  datout<-simdat$datout
  y<-simdat$datout$obs[1:nobs]
  #plot(y, type = "l")

  out_detfun0<-optdat$optout_det
  out_EDM<-optdat$optout_edm
  pars_true<-c(parslst$ptrue[1], log(parslst$sd_abs))

  #exp(pars_true)

  Euse<-parslst$Euse
  tuse<-parslst$tuse

  summarydat$summed_obs_error[ifl]<-sd(simdat$datout$true[1:nobs]-simdat$datout$obs[1:nobs])
  summarydat$summed_proc_error[ifl]<-sd(simdat$datout$true[1:nobs]-simdat$datout$noproc[1:nobs])

  summarydat$obs0[ifl]<-exp(pars_true[1])
  summarydat$obs1[ifl]<-NA
  summarydat$proc0[ifl]<-parslst$sd_abs
  summarydat$proc1[ifl]<-NA

  summarydat$det_obs_mu0[ifl]<-exp(parslst$parsest_det[1,1])
  summarydat$det_obs_mu1[ifl]<-NA
  summarydat$det_proc_mu0[ifl]<-exp(parslst$parsest_det[2,1])
  summarydat$det_proc_mu1[ifl]<-NA
  summarydat$edm_obs_mu0[ifl]<-exp(parslst$parsest_edm[1,1])
  summarydat$edm_obs_mu1[ifl]<-NA
  summarydat$edm_proc_mu0[ifl]<-exp(parslst$parsest_edm[2,1])
  summarydat$edm_proc_mu1[ifl]<-NA

  tmp<-try(gelmanDiagnostics(out_detfun0, plot = FALSE, start=sp)$psrf[,1], silent = TRUE)
  if(!is.character(tmp)) {
    summarydat$gelmandet[ifl]<-max(tmp)
  }
  tmp<-try(gelmanDiagnostics(out_EDM, plot = FALSE, start=sp)$psrf[,1], silent = TRUE)
  if(!is.character(tmp)) {
    summarydat$gelmanedm[ifl]<-max(tmp)
  }

  if(FALSE) {
    par(mar=c(4,4,2,2))
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
    par(mar=c(4,4,2,2), mfcol=c(2,2))
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
  pfout1_opt<-filterdat$filterout_det
  pfout2_opt<-filterdat$filterout_edm

  summarydat$cor_det[ifl]<-cor(pfout1_opt$Nest, datout$true[1:nobs], use="pairwise")^2
  summarydat$cor_edm[ifl]<-cor(pfout2_opt$Nest, datout$true[1:nobs], use="pairwise")^2
  summarydat$cor0[ifl]<-cor(y, datout$true[1:nobs], use="pairwise")^2

  summarydat$rmse_det[ifl]<-sqrt(mean((pfout1_opt$Nest-datout$true[1:nobs])^2))
  summarydat$rmse_edm[ifl]<-sqrt(mean((pfout2_opt$Nest-datout$true[1:nobs])^2))
  summarydat$rmse0[ifl]<-sqrt(mean((y-datout$true[1:nobs])^2))

  tmp<-exp(apply(smp_detfun0, 2, function(x) quantile(x, pnorm(-2:2))))
  summarydat[ifl,grep("det_obs_qt0", colnames(summarydat))]<-unname(tmp[,1])
  summarydat[ifl,grep("det_proc_qt0", colnames(summarydat))]<-unname(tmp[,2])

  tmp<-exp(apply(smp_EDM, 2, function(x) quantile(x, pnorm(-2:2))))
  summarydat[ifl,grep("edm_obs_qt0", colnames(summarydat))]<-unname(tmp[,1])
  summarydat[ifl,grep("edm_proc_qt0", colnames(summarydat))]<-unname(tmp[,2])

  summarydat$proc0_mu[ifl]<-sd(datout$obs[1:nobs])
  summarydat$proc_det_mu[ifl]<-sd(pfout1_opt$Nest)
  summarydat$proc_edm_mu[ifl]<-sd(pfout2_opt$Nest)
  summarydat$proc_true_mu[ifl]<-sd(datout$true[1:nobs])

  if(FALSE) {
    par(mfrow=c(1,1))
    plot(datout$true[1:nobs], datout$obs[1:nobs], col="red")
    points(datout$true[1:nobs], pfout2_opt$Nest[1:nobs], col="blue")
    points(datout$true[1:nobs], pfout1_opt$Nest[1:nobs], col="black")
    abline(a=0, b=1, lty=2)

    summarydat$cor_det[ifl]
    summarydat$cor_edm[ifl]
    summarydat$cor0[ifl]
  }

  stmp<-simdat$ssave
  summarydat$edm_var[ifl]<-stmp$rmse[[1]]
  summarydat$Euse[ifl]<-Euse
  summarydat$tuse[ifl]<-tuse


  #extinction and colonization rates
  #if(FALSE) {
    #based pfout1_opt detful0
    pfout1_opt<-particleFilterLL(y, pars=parseparam0(parslst$parsest_det[,1]), detfun = detfun0_sin,
                                    dotraceback = TRUE, fulltraceback = TRUE)
    #based on EDM
    pfout2_opt<-particleFilterLL(y, pars=parseparam0(parslst$parsest_edm[,1]), detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse, ytot=y),
                                    dotraceback = TRUE, fulltraceback = TRUE)
  #}

  prtout_det<-indexsort(pfout1_opt$fulltracemat, pfout1_opt$fulltraceindex, nsmp=1)
  prtout_edm<-indexsort(pfout2_opt$fulltracemat, pfout2_opt$fulltraceindex, nsmp=1)

  prtout_det_noproc<-indexsort(pfout1_opt$fulltracemat_noproc, pfout1_opt$fulltraceindex, nsmp=1)
  prtout_edm_noproc<-indexsort(pfout2_opt$fulltracemat_noproc, pfout2_opt$fulltraceindex, nsmp=1)

  cmdet<-getcm(prtout_det)
  cmedm<-getcm(prtout_edm)
  cmtrue<-getcm(datout$true)
  cmobs<-getcm(datout$obs[1:nobs])

  summarydat$pcdet[ifl]<-cmdet$pc
  summarydat$pcedm[ifl]<-cmedm$pc
  summarydat$pctrue[ifl]<-cmtrue$pc
  summarydat$pcobs[ifl]<-cmobs$pc

  summarydat$pmdet[ifl]<-cmdet$pm
  summarydat$pmedm[ifl]<-cmedm$pm
  summarydat$pmtrue[ifl]<-cmtrue$pm
  summarydat$pmobs[ifl]<-cmobs$pm

  summarydat$pc_actual[ifl]<-ilogit(pars0$pcol[1])
  summarydat$pm_actual[ifl]<-sum(datout$noproc>0 & datout$true==0)/sum(datout$noproc>0)

  xt<-datout$true[1:nobs]#datout$noproc[1:nobs]
  std<-parslst$sd_abs
  summarydat$pmtrue_analy[ifl]<-sum(pnorm(0, xt[xt>zero_cutoff], std))/sum(xt>zero_cutoff)

  xtdet<-prtout_det#prtout_det_noproc
  stddet<-exp(parslst$parsest_det[2,1])
  summarydat$pmdet_analy[ifl]<-sum(pnorm(0, xtdet[!is.na(xtdet) & xtdet>zero_cutoff], stddet))/sum(!is.na(xtdet) & xtdet>zero_cutoff)

  xtedm<-prtout_edm#prtout_edm_noproc
  stdedm<-exp(parslst$parsest_edm[2,1])
  summarydat$pmedm_analy[ifl]<-sum(pnorm(0, xtedm[!is.na(xtedm) & xtedm>zero_cutoff], stdedm))/sum(!is.na(xtedm) & xtedm>zero_cutoff)

  if(ifl/100 == floor(ifl/100)) {
    print(round(ifl/length(flst),2))
  }
  save(list = c("summarydat"), file = "datout/summarydat_211005.rda", version=2)
}
