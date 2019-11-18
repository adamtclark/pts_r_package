rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")

load("datout/mcmcout_120_full.rda")

sp<-1000
require(BayesianTools)
require(mvtnorm)
require(rEDM)
source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")
N<-2e3

p0<-list(c(-5, 0), c(-5, 0))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(-5, 0), c(-5, 0), c(-5, 2))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

datout<-simdat$datout
y<-simdat$datout$obs
out_detfun0<-optdat$optout_det
out_EDM<-optdat$optout_edm
pars_true<-parslst$ptrue

plot(out_detfun0, start=sp)
plot(out_EDM, start=sp)
gelmanDiagnostics(out_detfun0, plot = FALSE, start=sp)
gelmanDiagnostics(out_EDM, plot = FALSE, start=sp)

## Summarize outputs
smp_detfun0<-getSample(out_detfun0, start = sp)
smp_EDM<-getSample(out_EDM, start=sp)

#check for correlations among estimates
correlationPlot(out_detfun0, start = sp)
correlationPlot(out_EDM, start = sp)

#plot priors and posteriors
marginalPlot(out_detfun0, prior = TRUE, start = sp)
marginalPlot(out_EDM, prior = TRUE, start = sp)

#plot posteriors vs. true values
ptrue<-unlist(pars_true[1:2])

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

#run at optimum parameters
pfout1_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_detfun0)), N=N, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=TRUE)
pfout2_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_EDM[,1:2])), N=N, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=2, theta=exp(mean(smp_EDM[,3]))), dotraceback=TRUE)

cor(pfout1_opt$Nest, datout$true, use="pairwise")^2
cor(pfout2_opt$Nest, datout$true, use="pairwise")^2
cor(y, datout$true, use="pairwise")^2
mean((pfout1_opt$Nest-datout$true)^2)
mean((pfout2_opt$Nest-datout$true)^2)
mean((y-datout$true)^2)


exp(ptrue)
exp(apply(smp_detfun0, 2, function(x) quantile(x, pnorm(-2:2))))
exp(apply(smp_EDM, 2, function(x) quantile(x, pnorm(-2:2))))

par(mfrow=c(1,1))
plot(datout$true, datout$obs, col="red")
points(datout$true, pfout2_opt$Nest, col="blue")
points(datout$true, pfout1_opt$Nest, col="black")
abline(a=0, b=1, lty=2)


