rm(list=ls())
require(BayesianTools)
require(mvtnorm)
require(rEDM)
require(msir)
require(viridis)
require(pttstability)
#if(FALSE) {
  setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
  source("../pttstability/R/bayesfun.R")
  source("../pttstability/R/fake_data.R")
  source("../pttstability/R/logit_funs.R")
  source("../pttstability/R/particlefilter.R")
#}

eifun<-function(x,y,i=1,ybar=NULL) {
  if(is.null(ybar)) {
    ybar<-mean(y,na.rm=T)
  }
  1-mean(abs(x-y)^i,na.rm=T)/mean(abs(y-ybar)^i,na.rm=T)
}

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

p0<-list(c(log(0.01), log(0.5)), c(log(0.01), log(0.5)))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(log(0.01), log(0.5)), c(log(0.01), log(0.5)))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

flst<-dir("datout")
flst<-flst[grep("211005", flst)]
flst<-flst[grep("rda", flst)]
if(length(grep("summarydat", flst))>0) {
  flst<-flst[-grep("summarydat", flst)]
}

#Plot an example run
load("datout/summarydat_211005.rda")
pspos = which(summarydat$obs0>0.1 & summarydat$obs0<0.15 & summarydat$proc0<0.2 & summarydat$proc0>0.1 & !is.na(summarydat$gelmandet) & !is.na(summarydat$gelmanedm) & summarydat$gelmandet<1.1 & summarydat$gelmanedm<1.1)
psuse<-559 #pspos[17]
load(paste("datout/", flst[psuse], sep=""))
(cdet<-summarydat$cor_det[psuse])
(cedm<-summarydat$cor_edm[psuse])
(c0<-summarydat$cor0[psuse])

datout<-simdat$datout
y<-simdat$datout$obs[1:nobs]

out_detfun0<-optdat$optout_det
out_EDM<-optdat$optout_edm
pars_true<-parslst$ptrue
pars_true[2] = log(parslst$sd_abs)

Euse<-parslst$Euse
tuse<-parslst$tuse

smp_detfun0<-getSample(out_detfun0, start = sp)
smp_EDM<-getSample(out_EDM, start=sp)
ptrue<-unlist(pars_true)

pfout1_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_detfun0)), N=N, detfun=detfun0_sin, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=TRUE, fulltraceback = TRUE)
pfout2_opt<-particleFilterLL(y=y, pars=parseparam0(colMeans(smp_EDM)), N=N, detfun=EDMfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=list(E=Euse, theta=tuse, ytot=y), dotraceback=TRUE, fulltraceback = TRUE)

qprt_det<-cbind(pmax(0, pfout1_opt$Nest-pfout1_opt$Nsd),
                pfout1_opt$Nest,
                pfout1_opt$Nest+pfout1_opt$Nsd)

qprt_edm<-cbind(pmax(0, pfout2_opt$Nest-pfout2_opt$Nsd),
                pfout2_opt$Nest,
                pfout2_opt$Nest+pfout2_opt$Nsd)

collst<-c("gold", "firebrick", "dodgerblue")

pdf("plotout/example_timeseries.pdf",
    width=8, height=5, colormodel = "cmyk", useDingbats = FALSE)
  m = cbind(c(1,1,1,2),c(1,1,1,2),c(1,1,1,2),
            c(1,1,1,2),c(1,1,1,2),c(1,1,1,2),
            c(1,1,1,2),
            c(1,1,1,3),c(1,1,1,3))
  layout(m)

  par(mar=c(3,2,1,1), oma=c(1,3.4,0.5,0))
  plot(1,1, type="n", ylim=c(0, 3), xlim=c(0, 50), xlab="", ylab="", xaxs="i")
  polygon(c(1:150, 150:1), c(qprt_det[,1], rev(qprt_det[,3])), col=adjustcolor(collst[3], alpha.f = 0.3), border=NA)
  polygon(c(1:150, 150:1), c(qprt_edm[,1], rev(qprt_edm[,3])), col=adjustcolor(collst[2], alpha.f = 0.3), border=NA)
  lines(y, col=collst[1], lwd=1.5)
  lines(datout$true, lwd=1.5, lty = 2)
  abline(h=0, lty=3)
  title("a.", line=0.3, xpd=NA, adj=0.0, cex.main=1.5)

  legend("topright", legend = c("Analytical Function", "EDM Estimate", "Raw Observation", "True Dynamics"),
           fill = c(collst[3], collst[2], NA, NA),
           border = c(1,1, NA, NA),
           lty=c(NA, NA, 1, 2), lwd=c(NA, NA, 1.5, 1.5),
           col=c(NA, NA, collst[1], 1),
           bty="n", cex=1.2)
  mtext("Time", side = 1, line=2.5)
  mtext("Abundance", side = 2, line=2.5)

  bardat = cbind(c(mean((exp(parslst$parsest_det[1])*pfout1_opt$Nest)^2),
                        exp(parslst$parsest_det[2])^2,
                        var(pfout1_opt$Nest_noproc,na.rm=TRUE)),
                      c(mean((exp(parslst$parsest_edm[1])*pfout2_opt$Nest)^2),
                        exp(parslst$parsest_edm[2])^2,
                        var(pfout2_opt$Nest_noproc,na.rm=TRUE)),
                 c(0,var(y),0),
                 c(mean((exp(pars_true[1])*simdat$datout$true)^2),
                   exp(pars_true[2])^2,
                   var(simdat$datout$noproc)))
  bardat[3,-3] = var(y)-colSums(bardat[-3,-3]) # standardize to constant variance

  colnames(bardat) = c("Analytical Fun.", "EDM Estimate", "Observation","True Dynamic")

  acol = function(col, a=c(0.9, 0.4, 0)) {
    c(adjustcolor(col, a[1]),
      adjustcolor(col, a[2]),
      adjustcolor(col, a[3]))
  }

  par(mar=c(3,5.5,1,1))
  tmp = barplot(bardat[,4:1], axes=F, names.arg = rep("", ncol(bardat)),
                col = cbind(acol("black")), horiz = TRUE)
  axis(1)
  axis(2, at = tmp, labels = rev(colnames(bardat)), cex.axis=1.2, las=2)
  title("b.", line=0.3, xpd=NA, adj=0.0, cex.main=1.5)
  mtext("Temporal Variance", 1, line = 2.4)
  #box()
  legend(0.275, 5.8, legend = c("Observation Error", "Process Noise", "Deterministic Variation"),
         fill = acol("black"),
         border = 1,
         bty="n", xpd = NA, cex=1.2)
dev.off()
