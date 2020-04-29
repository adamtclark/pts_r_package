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

p0<-list(c(log(0.01), log(1)), c(log(0.01), log(1)), c(log(0.01), log(3)))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(log(0.01), log(1)), c(log(0.01), log(1)), c(log(0.01), log(3)))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

flst<-dir("datout")
flst<-flst[grep("200429", flst)]
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
                     pmobs=NA)


  for(ifl in 1:length(flst)) {
    #ifl<-sample(1:length(flst), 1)

    load(paste("datout/", flst[ifl], sep=""))

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

    if(ifl/100 == floor(ifl/100)) {
      print(round(ifl/length(flst),2))
    }
  }
  write.csv(summarydat, "datout/summarydat_allvar_oscil_taylor_200429.csv", row.names = FALSE)
} else {
  summarydat<-read.csv("datout/summarydat_allvar_oscil_taylor_200429.csv")
}



## set up for plotting
if(FALSE) {
cutoff<-0.2
cutlst<-c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 1)
mrowpar<-c(2,3)

summarydat$obsccut<-cut(summarydat$summed_obs_error, breaks = cutlst)
summarydat$proccut<-cut(summarydat$summed_proc_error, breaks = cutlst)
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
}





#tmp
cutoff<-0.2
cutoffo<-0.2
cutoffp<-0.2
ps<-sqrt(summarydat$summed_proc_error^2+summarydat$summed_obs_error^2)<cutoff
pso<-(summarydat$summed_proc_error)<cutoffo
psp<-(summarydat$summed_obs_error)<cutoffp
mean(ps); mean(pso); mean(psp)

hist(summarydat$summed_obs_error, breaks = 20); abline(v=cutoff, col=2, lty=3)
hist(summarydat$summed_proc_error, breaks = 20); abline(v=cutoff, col=2, lty=3)

lf<-function(x) {log10(x)}
#lf<-function(x) {(x)}
summary(modo0d<-lm(lf(obs0)~lf(det_obs_mu0), summarydat[pso & summarydat$gelmandet<=1.1,]))
summary(modp0d<-lm(lf(proc0)~lf(det_proc_mu0), summarydat[psp & summarydat$gelmandet<=1.1,]))
summary(modp1d<-lm(lf(proc1)~lf(det_proc_mu1), summarydat[psp & summarydat$gelmandet<=1.1,]))

summary(modo0e<-lm(lf(obs0)~lf(edm_obs_mu0), summarydat[pso & summarydat$gelmandet<=1.1,]))
summary(modp0e<-lm(lf(proc0)~lf(edm_proc_mu0), summarydat[psp & summarydat$gelmandet<=1.1,]))
summary(modp1e<-lm(lf(proc1)~lf(edm_proc_mu1), summarydat[psp & summarydat$gelmandet<=1.1,]))


#obs
plot(lf(obs0)~lf(det_obs_mu0), summarydat[pso & summarydat$gelmandet<=1.1,]); abline(a=0, b=1, lty=3)
abline(modo0d, col=2)
tmp<-summarydat[pso & summarydat$gelmandet<=1.1,]
plot(tapply(lf(tmp$det_obs_mu0), cut(lf(tmp$det_obs_mu0), 20), mean), tapply(lf(tmp$obs0), cut(lf(tmp$det_obs_mu0), 20), mean), xlab="pred", ylab="obs")
abline(a=0, b=1, lty=3)
abline(modo0d, col=2)

plot(lf(obs0)~lf(edm_obs_mu0), summarydat[pso & summarydat$gelmanedm<=1.1,]); abline(a=0, b=1, lty=3)
abline(modo0e, col=2)
tmp<-summarydat[pso & summarydat$gelmanedm<=1.1,]
plot(tapply(lf(tmp$edm_obs_mu0), cut(lf(tmp$edm_obs_mu0), 20), mean), tapply(lf(tmp$obs0), cut(lf(tmp$edm_obs_mu0), 20), mean), xlab="pred", ylab="obs")
abline(a=0, b=1, lty=3)
abline(modo0e, col=2)

#proc0
plot(lf(proc0)~lf(det_proc_mu0), summarydat[psp & summarydat$gelmandet<=1.1,]); abline(a=0, b=1, lty=3)
abline(modp0d, col=2)
tmp<-summarydat[psp & summarydat$gelmandet<=1.1,]
plot(tapply(lf(tmp$det_proc_mu0), cut(lf(tmp$det_proc_mu0), 20), mean), tapply(lf(tmp$proc0), cut(lf(tmp$det_proc_mu0), 20), mean), xlab="pred", ylab="obs")
abline(a=0, b=1, lty=3)
abline(modp0d, col=2)

plot(lf(proc0)~lf(edm_proc_mu0), summarydat[psp & summarydat$gelmanedm<=1.1,]); abline(a=0, b=1, lty=3)
abline(modp0e, col=2)
tmp<-summarydat[psp & summarydat$gelmanedm<=1.1,]
plot(tapply(lf(tmp$edm_proc_mu0), cut(lf(tmp$edm_proc_mu0), 20), mean), tapply(lf(tmp$proc0), cut(lf(tmp$edm_proc_mu0), 20), mean), xlab="pred", ylab="obs")
abline(a=0, b=1, lty=3)
abline(modp0e, col=2)

#proc1
plot(lf(proc1)~lf(det_proc_mu1), summarydat[psp & summarydat$gelmandet<=1.1,]); abline(a=0, b=1, lty=3)
abline(modp1d, col=2)
abline(v=lf(2), lty=2)
tmp<-summarydat[psp & summarydat$gelmandet<=1.1,]
plot(tapply(lf(tmp$det_proc_mu1), cut(lf(tmp$det_proc_mu1), 20), mean), tapply(lf(tmp$proc1), cut(lf(tmp$det_proc_mu1), 20), mean), xlab="pred", ylab="obs")
abline(a=0, b=1, lty=3)
abline(modp1d, col=2)

plot(lf(proc1)~lf(edm_proc_mu1), summarydat[psp & summarydat$gelmanedm<=1.1,]); abline(a=0, b=1, lty=3)
abline(modp1e, col=2)
abline(v=lf(2), lty=2)
tmp<-summarydat[psp & summarydat$gelmanedm<=1.1,]
plot(tapply(lf(tmp$edm_proc_mu1), cut(lf(tmp$edm_proc_mu1), 20), mean), tapply(lf(tmp$proc1), cut(lf(tmp$edm_proc_mu1), 20), mean), xlab="pred", ylab="obs")
abline(a=0, b=1, lty=3)
abline(modp1e, col=2)



#Check rates
plot(lf(pmtrue)~lf(pmdet), summarydat[psp & pso & summarydat$gelmandet<=1.1,]); abline(a=0, b=1, lty=3)
abline(modp1e, col=2)
abline(v=lf(2), lty=2)
tmp<-summarydat[psp & summarydat$gelmanedm<=1.1,]
plot(tapply(lf(tmp$edm_proc_mu1), cut(lf(tmp$edm_proc_mu1), 20), mean), tapply(lf(tmp$proc1), cut(lf(tmp$edm_proc_mu1), 20), mean), xlab="pred", ylab="obs")
abline(a=0, b=1, lty=3)
abline(modp1e, col=2)







coplot(lf(proc1)~lf(edm_proc_mu1)|lf(edm_proc_mu0), summarydat, panel=function(x,y,...) {points(x,y); abline(a=0, b=1, lty=3, col=2)})


pscor<-(summarydat$summed_obs_error>0.1)
plot(summarydat$cor0[pscor], summarydat$cor_det[pscor]); abline(a=0, b=1, lty=3)
plot(summarydat$cor0[pscor], summarydat$cor_edm[pscor]); abline(a=0, b=1, lty=3)
mean(summarydat$cor0[pscor]>summarydat$cor_det[pscor])
mean(summarydat$cor0[pscor]>summarydat$cor_edm[pscor])

#perhaps plot by category?
