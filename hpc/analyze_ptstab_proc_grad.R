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

nsteps<-20; nitertot<-100
prcsq<-rep(seq(1,2, length=nsteps), each=nitertot)
flst<-dir("datout")

qtllst<-c(0.025, pnorm(-1:1), 0.975)
simdatsum<-array(dim=c(length(qtllst),6,length(prcsq),2))
trueparsum<-matrix(nrow=length(prcsq),ncol=6)

for(ifl in 1:length(prcsq)) {
  flnm<-paste("datout/mcmcout_", ifl, ".rda", sep="")
  load(flnm)
  
  #summarize outputs
  pars$proc[2]<-log(procuse)
  
  truepars_transformed<-inv_fun0(t(as.matrix(unlist(pars))))
  
  smp_detfun0<-inv_fun0(getSample(out_detfun0))
  smp_EDM<-inv_fun0(getSample(out_EDM))
  
  simdatsum[,,ifl,1]<-apply(smp_detfun0, 2, function(x) quantile(x, qtllst, na.rm=T))
  simdatsum[,,ifl,2]<-apply(smp_EDM, 2, function(x) quantile(x, qtllst, na.rm=T))
  
  trueparsum[ifl,]<-truepars_transformed[1:6]
  
  if(ifl/20 == floor(ifl/20)) {
    print(round(ifl/length(prcsq),3))
  }
}


#get quantiles
prclst<-sort(unique(trueparsum[,4]))
qtlout_detful0<-array(dim=c(5,6,length(prclst)))
qtlout_EDM<-array(dim=c(5,6,length(prclst)))

for(i in 1:length(prclst)) {
  sbs<-which(trueparsum[,4]==prclst[i])
  qtlout_detful0[,,i]<-apply(simdatsum[,,sbs,1], 1:2, mean)
  qtlout_EDM[,,i]<-apply(simdatsum[,,sbs,2], 1:2, mean)
}


#plot
pltnames<-c("obs b0", "obs b1", "proc b0", "proc b1", "colp", "col0")
collst<-adjustcolor(c(4,2),alpha.f = 0.5)
dxl<-c(-0.01, 0.01)

pdf("plotout/plot_pstab_proc_grad.pdf", width=12, height=8, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,3), mar=c(4,4,2,2))
  for(i in 1:length(pltnames)) {
    if(i==4) {
      matplot(cbind(prclst+dxl[1], prclst+dxl[2]), cbind(qtlout_detful0[3,i,], qtlout_EDM[3,i,]),
              col=collst, type="p", pch=16,
              ylim=range(c(qtlout_detful0[,i,], qtlout_EDM[,i,])),
              xlab="true", ylab="estimated", main=pltnames[i])
      abline(a=0, b=1, lty=3)
      
      segments(prclst+dxl[1], qtlout_detful0[2,4,], prclst+dxl[1], qtlout_detful0[4,4,], col=collst[1], lend=2, lwd=3)
      segments(prclst+dxl[1], qtlout_detful0[1,4,], prclst+dxl[1], qtlout_detful0[5,4,], col=collst[1], lend=2, lwd=1)
      
      segments(prclst+dxl[2], qtlout_EDM[2,4,], prclst+dxl[2], qtlout_EDM[4,4,], col=collst[2], lend=2, lwd=3)
      segments(prclst+dxl[2], qtlout_EDM[1,4,], prclst+dxl[2], qtlout_EDM[5,4,], col=collst[2], lend=2, lwd=1)
    } else {
      matplot(cbind(prclst+dxl[1], prclst+dxl[2]), cbind(qtlout_detful0[3,i,], qtlout_EDM[3,i,]),
              col=collst, type="p", pch=16,
              ylim=range(c(qtlout_detful0[,i,], qtlout_EDM[,i,])),
              xlab="true proc b0", ylab="estimated", main=pltnames[i])
      abline(h=mean(trueparsum[,i]), lty=3)
      
      segments(prclst+dxl[1], qtlout_detful0[2,i,], prclst+dxl[1], qtlout_detful0[4,i,], col=collst[1], lend=2, lwd=3)
      segments(prclst+dxl[1], qtlout_detful0[1,i,], prclst+dxl[1], qtlout_detful0[5,i,], col=collst[1], lend=2, lwd=1)
      
      segments(prclst+dxl[2], qtlout_EDM[2,i,], prclst+dxl[2], qtlout_EDM[4,i,], col=collst[2], lend=2, lwd=3)
      segments(prclst+dxl[2], qtlout_EDM[1,i,], prclst+dxl[2], qtlout_EDM[5,i,], col=collst[2], lend=2, lwd=1)
    }
  }
dev.off()
