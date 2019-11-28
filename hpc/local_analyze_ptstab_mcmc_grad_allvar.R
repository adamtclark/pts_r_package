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
flst<-flst[grep("full", flst)]
flst<-flst[-grep("oscil", flst)]

if(FALSE) {
  summarydat<-data.frame(obs=rep(NA, length(flst)),
                     proc=NA,
                     det_obs_mu=NA,
                     det_proc_mu=NA,
                     edm_obs_mu=NA,
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
                     det_obs_qt=matrix(nrow=length(flst), ncol=5),
                     det_proc_qt=matrix(nrow=length(flst), ncol=5),
                     edm_obs_qt=matrix(nrow=length(flst), ncol=5),
                     edm_proc_qt=matrix(nrow=length(flst), ncol=5),
                     gelmandet=NA,
                     gelmanedm=NA,
                     edm_var=NA)


  for(ifl in 1:length(flst)) {
    load(paste("datout/", flst[ifl], sep=""))

    datout<-simdat$datout
    y<-simdat$datout$obs
    out_detfun0<-optdat$optout_det
    out_EDM<-optdat$optout_edm
    pars_true<-parslst$ptrue

    summarydat$obs[ifl]<-exp(pars_true[1])
    summarydat$proc[ifl]<-exp(pars_true[2])

    summarydat$det_obs_mu[ifl]<-exp(parslst$parsest_det[1,1])
    summarydat$det_proc_mu[ifl]<-exp(parslst$parsest_det[2,1])
    summarydat$edm_obs_mu[ifl]<-exp(parslst$parsest_edm[1,1])
    summarydat$edm_proc_mu[ifl]<-exp(parslst$parsest_edm[2,1])

    tmp<-try(gelmanDiagnostics(out_detfun0, plot = FALSE, start=sp)$psrf[,1], silent = TRUE)
    if(!is.character(tmp)) {
      summarydat$gelmandet[ifl]<-max(tmp)
    }
    tmp<-try(gelmanDiagnostics(out_EDM, plot = FALSE, start=sp)$psrf[,1], silent = TRUE)
    if(!is.character(tmp)) {
      summarydat$gelmanedm[ifl]<-max(tmp[1:2])
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
    ptrue<-unlist(pars_true[1:2])

    if(FALSE) {
      par(mar=c(4,4,2,2), mfrow=c(2,2))
      for(i in 1:2) {
        xrng<-exp(range(c(smp_detfun0[,i], ptrue[i], p0[[i]])))
        hist(exp(smp_detfun0[,i]),breaks = 20, probability = TRUE, main="", xlim=xrng);
        abline(v=exp(ptrue[i]), col=c(2), lty=2)
        abline(v=exp(p0[[i]]), col=c(3), lty=2)
      }
      for(i in 1:2) {
        xrng<-exp(range(c(smp_EDM[,i], ptrue[i], p0[[i]])))
        hist(exp(smp_EDM[,i]),breaks = 20, probability = TRUE, main="", xlim=xrng);
        abline(v=exp(ptrue[i]), col=c(2), lty=2)
        abline(v=exp(p0[[i]]), col=c(3), lty=2)
      }
    }

    #run at optimum parameters
    if(FALSE) {
      pfout1_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_detfun0)), N=N, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=TRUE)
      pfout2_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_EDM[,1:2])), N=N, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=2, theta=exp(mean(smp_EDM[,3]))), dotraceback=TRUE)
    }

    pfout1_opt<-filterdat$filterout_det
    pfout2_opt<-filterdat$filterout_edm

    summarydat$cor_det[ifl]<-cor(pfout1_opt$Nest, datout$true, use="pairwise")^2
    summarydat$cor_edm[ifl]<-cor(pfout2_opt$Nest, datout$true, use="pairwise")^2
    summarydat$cor0[ifl]<-cor(y, datout$true, use="pairwise")^2

    summarydat$rmse_det[ifl]<-sqrt(mean((pfout1_opt$Nest-datout$true)^2))
    summarydat$rmse_edm[ifl]<-sqrt(mean((pfout2_opt$Nest-datout$true)^2))
    summarydat$rmse0[ifl]<-sqrt(mean((y-datout$true)^2))

    tmp<-exp(apply(smp_detfun0, 2, function(x) quantile(x, pnorm(-2:2))))
    summarydat[ifl,grep("det_obs_qt", colnames(summarydat))]<-unname(tmp[,1])
    summarydat[ifl,grep("det_proc_qt", colnames(summarydat))]<-unname(tmp[,2])

    tmp<-exp(apply(smp_EDM, 2, function(x) quantile(x, pnorm(-2:2))))
    summarydat[ifl,grep("edm_obs_qt", colnames(summarydat))]<-unname(tmp[,1])
    summarydat[ifl,grep("edm_proc_qt", colnames(summarydat))]<-unname(tmp[,2])


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
    }

    stmp<-s_map(datout$obs,silent = TRUE,E=2)
    summarydat$edm_var[ifl]<-min(stmp$rmse,na.rm=T)

    if(ifl/100 == floor(ifl/100)) {
      print(round(ifl/length(flst),2))
    }
  }
  write.csv(summarydat, "datout/summarydat.csv", row.names = FALSE)
} else {
  summarydat<-read.csv("datout/summarydat.csv")
}

## set up for plotting
cutoff<-0.3
cutlst<-c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 1)
mrowpar<-c(2,3)

summarydat$obsccut<-cut(summarydat$obs, breaks = cutlst)
summarydat$proccut<-cut(summarydat$proc, breaks = cutlst)
obscutlst<-sort(unique(summarydat$obsccut))
proccutlst<-sort(unique(summarydat$proccut))

#show total error
plot(sqrt(summarydat$obs^2+summarydat$proc^2), summarydat$edm_var, ylab="est tot error", xlab="obs tot error", col=collst[3], log="xy"); abline(a=0, b=1, lty=3)
points(sqrt(summarydat$obs^2+summarydat$proc^2), summarydat$proc0_mu, col=collst[1])
abline(v=cutoff, h=cutoff, lty=3)

## plot proc error
plotfun(plotvar="proc", byvar="obs", summarydat=summarydat, cutlst=cutlst, mrowpar=mrowpar, collst=collst, xlim=c(0.02, 1), ylim=c(0.02,1), doci=FALSE, cutoff = cutoff)

## plot obs error
plotfun(plotvar="obs", byvar="proc", summarydat=summarydat, cutlst=cutlst, mrowpar=mrowpar, collst=collst, xlim=c(0.02, 1), ylim=c(0.02,1), doci=FALSE, cutoff = cutoff)

## plot model fit
plotfun(plotvar="fit", byvar="proc", summarydat=summarydat, cutlst=cutlst, mrowpar=mrowpar, collst=collst, xlim=c(0.02, 1), ylim=c(0.01,1), cutoff = cutoff)

## plot variability of true dynamics
plotfun(plotvar="det", byvar="proc", summarydat=summarydat, cutlst=cutlst, mrowpar=mrowpar, collst=collst, xlim=c(0.02, 1), ylim=c(0.05,1), cutoff = cutoff)

