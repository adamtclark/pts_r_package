rm(list=ls())

setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")
require(BayesianTools)
require(mvtnorm)
require(rEDM)
require(msir)
require(viridis)
require(pttstability)

eifun<-function(x,y,i=1,ybar=NULL) {
  if(is.null(ybar)) {
    ybar<-mean(y,na.rm=T)
  }
  1-mean(abs(x-y)^i,na.rm=T)/mean(abs(y-ybar)^i,na.rm=T)
}

#new detfun
pars0<-pars_true<-list(obs=c(log(0.2)),
                       proc=c(log(0.2), log(1.5)),
                       pcol=c(logit(0.2), log(0.1)),
                       det=c(log(1.2),log(1)))

detfun0_sin<-function(sdet, xt, time=NULL) {
  K<-(sin(time/2)+exp(sdet[2])+0.5)*(2/3)
  xt = xt*exp(exp(sdet[1])*(1-xt/K))
  return(xt)
}

#settings
sp<-2000
N<-2e3
collst<-adjustcolor(c("purple", "blue", "red", "black"),alpha.f = 0.3)

p0<-list(c(log(0.01), log(0.5)), c(log(0.01), log(0.5)), c(log(0.5), log(3)))
minvUSE<-unlist(lapply(p0, function(x) x[1]))
maxvUSE<-unlist(lapply(p0, function(x) x[2]))

p0_edm<-list(c(log(0.01), log(0.5)), c(log(0.01), log(0.5)), c(log(0.5), log(3)))
minvUSE_edm<-unlist(lapply(p0_edm, function(x) x[1]))
maxvUSE_edm<-unlist(lapply(p0_edm, function(x) x[2]))

flst<-dir("datout")
flst<-flst[grep("200430", flst)]
flst<-flst[grep("rda", flst)]






#Plot an example run
summarydat<-read.csv("datout/summarydat_allvar_oscil_taylor_200430.csv")

#631, 735, 851, 985
psuse<-735
load(paste("datout/", flst[psuse], sep=""))
cdet<-summarydat$cor_det[psuse]
cedm<-summarydat$cor_edm[psuse]
c0<-summarydat$cor0[psuse]

datout<-simdat$datout
y<-simdat$datout$obs

out_detfun0<-optdat$optout_det
out_EDM<-optdat$optout_edm
pars_true<-parslst$ptrue

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

pdf("plotout/local_analyze_allvar_oscil_taylor_200429_example_timeseries.pdf",
    width=6, height=5, colormodel = "cmyk", useDingbats = FALSE)
  #m<-cbind(c(1,1,1), c(1,1,1), c(1,1,1), c(1,1,1), c(1,1,1), c(1,1,1), c(2,3,4), c(2,3,4), c(2,3,4), c(5,6,7), c(5,6,7),  c(5,6,7), c(8,8,8),c(8,8,8),c(8,8,8),c(8,8,8))
  m<-cbind(c(1,1,1,1,1,1,1,1,8,8,8,8,8,8,8), c(1,1,1,1,1,1,1,1,8,8,8,8,8,8,8), c(2,2,2,2,2,3,3,3,3,3,4,4,4,4,4), c(5,5,5,5,5,6,6,6,6,6,7,7,7,7,7))
  layout(m)

  par(mar=c(3,2,1,3), oma=c(1,2,0,0))
  plot(1,1, type="n", ylim=c(0, 3), xlim=c(0, 50), xlab="", ylab="", xaxs="i")
  polygon(c(1:150, 150:1), c(qprt_det[,1], rev(qprt_det[,3])), col=adjustcolor(collst[3], alpha.f = 0.3), border=NA)
  polygon(c(1:150, 150:1), c(qprt_edm[,1], rev(qprt_edm[,3])), col=adjustcolor(collst[2], alpha.f = 0.3), border=NA)
  lines(y, col=collst[1], lwd=1.5)
  lines(datout$true, lwd=1.5)
  abline(h=0, lty=3)
  title("a.", line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)

  legend("topright", legend = c("Analytical Function", "EDM Estimate", "Raw Observation", "True Dynamics"),
           fill = c(collst[3], collst[2], NA, NA),
           border = c(1,1, NA, NA),
           lty=c(NA, NA, 1, 1), lwd=c(NA, NA, 1.5, 1.5),
           col=c(NA, NA, collst[1], 1),
           bty="n")
  mtext("Time", side = 1, line=2.5)
  mtext("Abundance", side = 2, line=2.5)

  par(mar=c(3,1.5,1,1))
  parnames<-c(expression(paste("Observation Error Slope, ", beta[obs])),
              expression(paste("Process Noise Intercept, ", beta[proc[0]])),
              expression(paste("Process Noise Slope, ", beta[proc[1]])))
  lnnum<-c(-25, -12.5, -0.2)
  for(i in 1:length(ptrue)) {
    xrng<-exp(range(c(smp_detfun0[,i], ptrue[i], p0[[i]])))
    dns<-density(exp(smp_detfun0[,i]), from = xrng[1], to = xrng[2], bw = diff(range(xrng))*0.02)
    plot(dns$x, dns$y, xlim=xrng, col=collst[3], lwd=1.5, type="l",
         xlab="", ylab="")
    title(paste(letters[i+1], ".", sep=""), line=-1.05, xpd=NA, adj=0.08, cex.main=1.5)

    if(i==2) {
      mtext("Probability Density", side = 2, line=2.2)
    }

    abline(v=exp(p0[[i]]), h=0, col=c(1), lty=3)
    abline(v=exp(ptrue[i]), col=c(1), lty=2, lwd=1.5)

    mtext(parnames[i], side=1, adj=0.9, line=lnnum[i], outer=TRUE, xpd=NA)
  }
  for(i in 1:length(ptrue)) {
    xrng<-exp(range(c(smp_EDM[,i], ptrue[i], p0_edm[[i]])))
    dns<-density(exp(smp_EDM[,i]), from = xrng[1], to = xrng[2], bw = diff(range(xrng))*0.02)
    plot(dns$x, dns$y, xlim=xrng, col=collst[2], lwd=1.5, type="l",
         xlab="", ylab="")
    title(paste(letters[i+4], ".", sep=""), line=-1.05, xpd=NA, adj=0.08, cex.main=1.5)

    abline(v=exp(p0[[i]]), h=0, col=c(1), lty=3)
    abline(v=exp(ptrue[i]), col=c(1), lty=2, lwd=1.5)
  }

  par(mar=c(3,2,1,3))
  mxy<-max(c(datout$true, y, pfout1_opt$Nest, pfout2_opt$Nest))
  plot(y, datout$true, col=adjustcolor(collst[1], alpha.f = 0.5), xlim=c(0, mxy), ylim=c(0, mxy), cex=0.9, pch=18, xlab="", ylab="")
  points(pfout1_opt$Nest, datout$true, col=adjustcolor(collst[3], alpha.f = 0.5), cex=0.9, pch=16)
  points(pfout2_opt$Nest, datout$true, col=adjustcolor(collst[2], alpha.f = 0.5), cex=0.9, pch=17)
  abline(a=0, b=1, lty=2, lwd=1.5)
  title("h.", line=-1.05, xpd=NA, adj=0.02, cex.main=1.5)

  legend("bottomright",
         legend=c(as.expression(bquote(rho==.(paste(round(cdet,2), ";", sep="")) ~ E[1] ~ "=" ~ .(round(eifun(pfout1_opt$Nest, datout$true), 2)))),
                  as.expression(bquote(rho==.(paste(round(cedm,2), ";", sep="")) ~ E[1] ~ "=" ~ .(round(eifun(pfout2_opt$Nest, datout$true), 2)))),
                  as.expression(bquote(rho==.(paste(round(c0,2), ";", sep="")) ~ E[1] ~ "=" ~ .(round(eifun(y, datout$true), 2))))),
         pch=16:18, col=c(collst[3], collst[2], collst[1]), bty="n")
  mtext("Estimated Abundance", side = 1, line=2.5)
  mtext("True Abundance", side = 2, line=2.5)

dev.off()
