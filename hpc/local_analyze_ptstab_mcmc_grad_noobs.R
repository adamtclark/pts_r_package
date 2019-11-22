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

#settings
sp<-1000
N<-2e3
collst<-adjustcolor(c("purple", "blue", "red", "black"),alpha.f = 0.3)

p0<-list(c(-5, 0), c(-5, 0))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(-5, 0), c(-5, 0), c(-5, 2))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

flst<-dir("datout")
flst<-flst[grep("mcmc", flst)]
flst<-flst[grep("noobs", flst)]

if(FALSE) {
  summarydat<-data.frame(obs=rep(NA, length(flst)),
                     proc=NA,
                     det_proc_mu=NA,
                     edm_proc_mu=NA,
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
                     det_proc_qt=matrix(nrow=length(flst), ncol=5),
                     edm_proc_qt=matrix(nrow=length(flst), ncol=5),
                     gelmandet=NA,
                     gelmanedm=NA)

  for(ifl in 1:length(flst)) {
    load(paste("datout/", flst[ifl], sep=""))

    datout<-simdat$datout
    y<-simdat$datout$obs
    out_detfun0<-optdat$optout_det
    out_EDM<-optdat$optout_edm
    pars_true<-parslst$ptrue

    summarydat$obs[ifl]<-exp(pars_true[1])
    summarydat$proc[ifl]<-exp(pars_true[2])

    summarydat$det_proc_mu[ifl]<-exp(parslst$parsest_det[1,1])
    summarydat$edm_proc_mu[ifl]<-exp(parslst$parsest_edm[1,1])


    if(FALSE) {
      plot(out_detfun0, start=sp)
      plot(out_EDM, start=sp)
    }
    summarydat$gelmandet[ifl]<-gelmanDiagnostics(out_detfun0, plot = FALSE, start=sp)$psrf[1]
    summarydat$gelmanedm[ifl]<-gelmanDiagnostics(out_EDM, plot = FALSE, start=sp)$psrf[1]

    ## Summarize outputs
    smp_detfun0<-getSample(out_detfun0, start = sp)
    smp_EDM<-getSample(out_EDM, start=sp)

    #plot posteriors vs. true values
    ptrue<-unlist(pars_true[1:2])

    pfout1_opt<-filterdat$filterout_det
    pfout2_opt<-filterdat$filterout_edm

    summarydat$cor_det[ifl]<-cor(pfout1_opt$Nest, datout$true, use="pairwise")^2
    summarydat$cor_edm[ifl]<-cor(pfout2_opt$Nest, datout$true, use="pairwise")^2
    summarydat$cor0[ifl]<-cor(y, datout$true, use="pairwise")^2

    summarydat$rmse_det[ifl]<-sqrt(mean((pfout1_opt$Nest-datout$true)^2))
    summarydat$rmse_edm[ifl]<-sqrt(mean((pfout2_opt$Nest-datout$true)^2))
    summarydat$rmse0[ifl]<-sqrt(mean((y-datout$true)^2))

    tmp<-exp(apply(smp_detfun0, 2, function(x) quantile(x, pnorm(-2:2))))
    summarydat[ifl,grep("det_proc_qt", colnames(summarydat))]<-unname(tmp[,1])

    tmp<-exp(apply(smp_EDM, 2, function(x) quantile(x, pnorm(-2:2))))
    summarydat[ifl,grep("edm_proc_qt", colnames(summarydat))]<-unname(tmp[,1])

    summarydat$proc0_mu[ifl]<-sd(datout$obs)
    summarydat$proc_det_mu[ifl]<-sd(pfout1_opt$Nest)
    summarydat$proc_edm_mu[ifl]<-sd(pfout2_opt$Nest)
    summarydat$proc_true_mu[ifl]<-sd(datout$true)
  }
  write.csv(summarydat, "datout/summarydat_noobs.csv", row.names = FALSE)
} else {
  summarydat<-read.csv("datout/summarydat_noobs.csv")
}

## set up for plotting
cutlst<-c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 1)
mrowpar<-c(2,3)

summarydat$obsccut<-cut(summarydat$obs, breaks = cutlst)
summarydat$proccut<-cut(summarydat$proc, breaks = cutlst)
obscutlst<-sort(unique(summarydat$obsccut))
proccutlst<-sort(unique(summarydat$proccut))

## plot proc error
plotfun(plotvar="proc", byvar="obs", summarydat=summarydat, cutlst=cutlst, mrowpar=mrowpar, collst=collst, xlim=c(0.02, 1), ylim=c(0.02,1), doci=FALSE)

## plot model fit
plotfun(plotvar="fit", byvar="proc", summarydat=summarydat, cutlst=cutlst, mrowpar=mrowpar, collst=collst, xlim=c(0.02, 1), ylim=c(0.01,1))

## plot variability of true dynamics
plotfun(plotvar="det", byvar="proc", summarydat=summarydat, cutlst=cutlst, mrowpar=mrowpar, collst=collst, xlim=c(0.02, 1), ylim=c(0.05,1))


