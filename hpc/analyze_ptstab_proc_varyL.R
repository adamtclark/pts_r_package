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
           proc=c(-2, log(1.5)),
           pcol=c(logit(0.2), log(1e-2)),
           det=c(log(3),log(1)))


nsteps<-6; nitertot<-100
lsq<-rep(c(10, 20, 30, 50, 75, 100), each=nitertot)
flst<-dir("datout")
rmps<-grep("varyL", flst)
if(length(rmps)>0) {
  flst<-flst[rmps]
}

qtllst<-c(0.025, pnorm(-1:1), 0.975)
simdatsum<-array(dim=c(length(qtllst),6,length(lsq),2))
trueparsum<-matrix(nrow=length(lsq),ncol=6)
colrates<-matrix(nrow=length(lsq), ncol=4)
morrates<-matrix(nrow=length(lsq), ncol=4)
morrates_pest<-matrix(nrow=length(lsq), ncol=4)
libsq<-numeric(length(lsq))

if(FALSE) {
  for(ifl in 1:length(flst)) {
    flnm<-paste("datout/", flst[ifl], sep="")
    load(flnm)

    #summarize outputs
    truepars_transformed<-inv_fun0(t(as.matrix(unname(unlist(pars_sim)))))
    Luse<-length(datout$true)
    libsq[ifl]<-Luse
    
    smp_detfun0_untr<-(getSample(out_detfun0))
    smp_EDM_untr<-(getSample(out_EDM))

    smp_detfun0<-inv_fun0(smp_detfun0_untr)
    smp_EDM<-inv_fun0(smp_EDM_untr)

    simdatsum[,,ifl,1]<-apply(smp_detfun0, 2, function(x) quantile(x, qtllst, na.rm=T))
    simdatsum[,,ifl,2]<-apply(smp_EDM, 2, function(x) quantile(x, qtllst, na.rm=T))

    trueparsum[ifl,]<-truepars_transformed[1:6]

    #get col/mor
    pfdet<-demdat$pfdet
    pfedm<-demdat$pfedm
    pftrue<-demdat$pftrue
    
    morrates_pest[ifl,1]<-mean(pnorm(0, datout$true[datout$true>0], sqrt(exp(colSums(rbind(1, log(datout$true[datout$true>0]))*truepars_transformed[3:4])))))
    morrates_pest[ifl,2]<-mean(pnorm(0, pftrue$rN[pftrue$rN>0], sqrt(exp(colSums(rbind(1, log(pftrue$rN[pftrue$rN>0]))*truepars_transformed[3:4])))))
    morrates_pest[ifl,3]<-mean(pnorm(0, pfdet$rN[pfdet$rN>0], sqrt(exp(colSums(rbind(1, log(pfdet$rN[pfdet$rN>0]))*simdatsum[3,3:4,ifl,1])))))
    morrates_pest[ifl,4]<-mean(pnorm(0, pfedm$rN[pfedm$rN>0], sqrt(exp(colSums(rbind(1, log(pfedm$rN[pfedm$rN>0]))*simdatsum[3,3:4,ifl,2])))))
    
    truecol<-demdat$truecol
    truemor<-demdat$truemor

    colrates[ifl,1]<-truecol
    colrates[ifl,2]<-pftrue$dem$mucol
    colrates[ifl,3]<-pfdet$dem$mucol
    colrates[ifl,4]<-pfedm$dem$mucol

    morrates[ifl,1]<-truemor
    morrates[ifl,2]<-pftrue$dem$mumor
    morrates[ifl,3]<-pfdet$dem$mumor
    morrates[ifl,4]<-pfedm$dem$mumor

    if(ifl/20 == floor(ifl/20)) {
      print(round(ifl/length(flst),3))
    }
  }
  save.image("summarydata/saved_summary_varL.rda")
} else {
  load("summarydata/saved_summary_varL.rda")
}


#get quantiles
liblst<-sort(unique(lsq))
qtlout_detful0<-array(dim=c(5,6,length(liblst)))
qtlout_EDM<-array(dim=c(5,6,length(liblst)))
qtlout_col<-array(dim=c(5,4,length(liblst)))
qtlout_mor<-array(dim=c(5,4,length(liblst)))
qtlout_mor_p<-array(dim=c(5,4,length(liblst)))

for(i in 1:length(liblst)) {
  sbs<-which(libsq==liblst[i])
  qtlout_detful0[,,i]<-apply(simdatsum[,,sbs,1], 1:2, mean)
  qtlout_EDM[,,i]<-apply(simdatsum[,,sbs,2], 1:2, mean)

  qtlout_col[,,i]<-apply(colrates[sbs,], 2, function(x) quantile(x, qtllst, na.rm=T))
  qtlout_mor[,,i]<-apply(morrates[sbs,], 2, function(x) quantile(x, qtllst, na.rm=T))
  qtlout_mor_p[,,i]<-apply(morrates_pest[sbs,], 2, function(x) quantile(x, qtllst, na.rm=T))
}


#plot parameters
pltnames<-c("obs b0", "obs b1", "proc b0", "proc b1", "colp", "col0")
collst<-adjustcolor(c(4,2),alpha.f = 0.5)
dxl<-c(-1, 1)

pdf("plotout/plot_pstab_proc_varyL.pdf", width=12, height=8, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,3), mar=c(4,4,2,2))
  for(i in 1:length(pltnames)) {
    matplot(cbind(liblst+dxl[1], liblst+dxl[2]), cbind(qtlout_detful0[3,i,], qtlout_EDM[3,i,]),
            col=collst, type="p", pch=16,
            ylim=range(c(qtlout_detful0[,i,], qtlout_EDM[,i,])),
            xlab="L length", ylab="estimated", main=pltnames[i])
    abline(h=mean(trueparsum[,i]), lty=3)

    segments(liblst+dxl[1], qtlout_detful0[2,i,], liblst+dxl[1], qtlout_detful0[4,i,], col=collst[1], lend=2, lwd=3)
    segments(liblst+dxl[1], qtlout_detful0[1,i,], liblst+dxl[1], qtlout_detful0[5,i,], col=collst[1], lend=2, lwd=1)

    segments(liblst+dxl[2], qtlout_EDM[2,i,], liblst+dxl[2], qtlout_EDM[4,i,], col=collst[2], lend=2, lwd=3)
    segments(liblst+dxl[2], qtlout_EDM[1,i,], liblst+dxl[2], qtlout_EDM[5,i,], col=collst[2], lend=2, lwd=1)
  }



  dxlst<-seq(-3, 3, length=4)
  #plot rates
  par(mfrow=c(1,3))
  plot(range(liblst)+range(dxlst), range(qtlout_col[2:4,,]), type="n", xlab="L length", ylab="est col. rate")
  collst<-c("black", "blue", "red", "green")
  
  for(i in 1:4) {
    lines(liblst+dxlst[i], qtlout_col[3,i,], col=collst[i], type="b")
    segments(liblst+dxlst[i], qtlout_col[2,i,], liblst+dxlst[i], qtlout_col[4,i,], lend=3, col=collst[i], lwd=2)
    segments(liblst+dxlst[i], qtlout_col[1,i,], liblst+dxlst[i], qtlout_col[5,i,], lend=3, col=collst[i], lwd=1)
  }
  abline(h=c(0,1), lty=3)
  
  plot(range(liblst)+range(dxlst), range(qtlout_mor[2:4,,]), type="n", xlab="L length", ylab="est mor. rate")
  for(i in 1:4) {
    lines(liblst+dxlst[i], qtlout_mor[3,i,], col=collst[i], type="b")
    segments(liblst+dxlst[i], qtlout_mor[2,i,], liblst+dxlst[i], qtlout_mor[4,i,], lend=3, col=collst[i], lwd=2)
    segments(liblst+dxlst[i], qtlout_mor[1,i,], liblst+dxlst[i], qtlout_mor[5,i,], lend=3, col=collst[i], lwd=1)
  }
  abline(h=c(0,1), lty=3)
  
  
  plot(range(liblst)+range(dxlst), range(qtlout_mor_p[2:4,,]), type="n", xlab="L length", ylab="est mor. rate")
  for(i in 1:4) {
    lines(liblst+dxlst[i], qtlout_mor_p[3,i,], col=collst[i], type="b")
    segments(liblst+dxlst[i], qtlout_mor_p[2,i,], liblst+dxlst[i], qtlout_mor_p[4,i,], lend=3, col=collst[i], lwd=2)
    segments(liblst+dxlst[i], qtlout_mor_p[1,i,], liblst+dxlst[i], qtlout_mor_p[5,i,], lend=3, col=collst[i], lwd=1)
  }
  legend("topright", c("true", "true_pf", "det.", "edm"), ncol=2, col=c("black", "blue", "red", "green"), lty=1, lwd=2, bty="n")
  abline(h=c(0,1), lty=3)
  
  
dev.off()


