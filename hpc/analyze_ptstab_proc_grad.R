error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")

#load packages and functions
require(BayesianTools); require(rEDM)
source("../pttstability/R/bayesfun.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/particlefilter.R")

#load data
pars<-list(obs=c(log(1e-2), log(0.1)),
           proc=c(-2, NA),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(log(3),log(1)))

simdatsum<-array(dim=c(4002,6,20,2))
trueparsum<-matrix(nrow=20,ncol=6)

for(ifl in 1:20) {
  flst<-dir("datout")
  
  flnm<-paste("datout/mcmcout_", ifl, ".rda", sep="")
  load(flnm)
  
  #summarize outputs
  procuse<-seq(1,2,length=20)[ifl]
  pars$proc[2]<-log(procuse)
  
  truepars_transformed<-inv_fun0(t(as.matrix(unlist(pars))))
  
  smp_detfun0<-inv_fun0(getSample(out_detfun0))
  smp_EDM<-inv_fun0(getSample(out_EDM))
  
  simdatsum[,,ifl,1]<-smp_detfun0
  simdatsum[,,ifl,2]<-smp_EDM
  
  trueparsum[ifl,]<-truepars_transformed[1:6]
}


#get quantiles
qtlf<-function(x) {quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=T)}
qtlout_detful0<-apply(simdatsum[,,,1],2:3,qtlf)
qtlout_EDM<-apply(simdatsum[,,,2],2:3,qtlf)

#plot
pltnames<-c("obs b0", "obs b1", "proc b0", "proc b1", "colp", "col0")
collst<-adjustcolor(c(4,2),alpha.f = 0.5)
dxl<-c(-0.01, 0.01)

pdf("plotout/plot_pstab_proc_grad.pdf", width=12, height=8, colormodel = "cmyk", useDingbats = FALSE)
par(mfrow=c(2,3), mar=c(4,4,2,2))
for(i in 1:length(pltnames)) {
  if(i==4) {
    matplot(cbind(trueparsum[,i]+dxl[1], trueparsum[,i]+dxl[2]), cbind(qtlout_detful0[3,i,], qtlout_EDM[3,i,]),
            col=collst, type="p", pch=16,
            ylim=range(c(qtlout_detful0[,i,], qtlout_EDM[,i,])),
            xlab="true", ylab="estimated", main=pltnames[i])
    abline(a=0, b=1, lty=3)
    
    segments(trueparsum[,4]+dxl[1], qtlout_detful0[2,4,], trueparsum[,4]+dxl[1], qtlout_detful0[4,4,], col=collst[1], lend=2, lwd=3)
    segments(trueparsum[,4]+dxl[1], qtlout_detful0[1,4,], trueparsum[,4]+dxl[1], qtlout_detful0[5,4,], col=collst[1], lend=2, lwd=1)
    
    segments(trueparsum[,4]+dxl[2], qtlout_EDM[2,4,], trueparsum[,4]+dxl[2], qtlout_EDM[4,4,], col=collst[2], lend=2, lwd=3)
    segments(trueparsum[,4]+dxl[2], qtlout_EDM[1,4,], trueparsum[,4]+dxl[2], qtlout_EDM[5,4,], col=collst[2], lend=2, lwd=1)
  } else {
    matplot(cbind(trueparsum[,4]+dxl[1], trueparsum[,4]+dxl[2]), cbind(qtlout_detful0[3,i,], qtlout_EDM[3,i,]),
            col=collst, type="p", pch=16,
            ylim=range(c(qtlout_detful0[,i,], qtlout_EDM[,i,])),
            xlab="true proc b0", ylab="estimated", main=pltnames[i])
    abline(h=mean(trueparsum[,i]), lty=3)
    
    segments(trueparsum[,4]+dxl[1], qtlout_detful0[2,i,], trueparsum[,4]+dxl[1], qtlout_detful0[4,i,], col=collst[1], lend=2, lwd=3)
    segments(trueparsum[,4]+dxl[1], qtlout_detful0[1,i,], trueparsum[,4]+dxl[1], qtlout_detful0[5,i,], col=collst[1], lend=2, lwd=1)
    
    segments(trueparsum[,4]+dxl[2], qtlout_EDM[2,i,], trueparsum[,4]+dxl[2], qtlout_EDM[4,i,], col=collst[2], lend=2, lwd=3)
    segments(trueparsum[,4]+dxl[2], qtlout_EDM[1,i,], trueparsum[,4]+dxl[2], qtlout_EDM[5,i,], col=collst[2], lend=2, lwd=1)
  }
}
dev.off()
