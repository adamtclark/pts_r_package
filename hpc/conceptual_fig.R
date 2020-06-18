error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")

#load scripts
require(rEDM)
require(scatterplot3d)
source("../pttstability/R/particlefilter.R")
source("../pttstability/R/fake_data.R")
source("../pttstability/R/logit_funs.R")
source("../pttstability/R/bayesfun.R")

#set up simulations
#seeduse<-round(runif(1)*1e5)
#69663
#37519
set.seed(37519)

stm<-1.7
detfun0_sin<-function(sdet, xt, time=NULL) {
  K<-(sin(time/stm)+exp(sdet[2])+0.5)/2
  xt = xt*exp(exp(sdet[1])*(1-xt/K))
  return(xt)
}


pars_sim<-pars_true<-list(obs=c(log(0.3)),
                          proc=c(log(0.1), log(1.5)),
                          pcol=c(logit(0.2), log(0.1)),
                          det=c(log(1.2),log(1)))


pdf("plotout/workflow.pdf", width=6, height=12, colormodel = "cmyk", useDingbats = FALSE)
layout(rbind(c(1,1),
             c(6,1),
             c(6,2),
             c(6,2),
             c(3,3),
             c(3,3),
             c(4,5),
             c(4,5)))
par(mar=c(3,3,2,2), oma=c(2,2,0,0))

#step 1: make and plot the time series:
n<-500
datout<-makedynamics_general(n = n, n0 = (sin(1/2)+1+0.5)/2, pdet=pars_sim$det,
                             proc = pars_sim$proc, obs = pars_sim$obs, pcol = pars_sim$pcol,
                             detfun = detfun0_sin, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
#plot(datout$true, type="l")

tm<-(101:150)-10
obstm<-seq(tm[1], tm[length(tm)], by=1)

x<-datout$true[tm]
#mnx<-mean(x)
mnx<-1
x<-x/mnx

obs<-datout$obs[tm]
obs<-obs/mnx

plot(tm-min(tm)+1, x, type="l", lwd=1.5, xlab="", ylab="",
     ylim=c(0, 2))
abline(h=0, lty=3)
points(tm-min(tm)+1, obs, col="gold", pch=18, cex=1.2)
segments(tm-min(tm)+1, x, tm-min(tm)+1, obs, col="gold", lwd=1.5)
mtext(expression(paste("Time")), 1, line = 2.4, cex=1.2)
mtext(expression(paste("Abundance")), 2, line = 2.4, cex=1.2)
title(main="a.", cex.main=1.6, adj=0.01, line=-1)

legend("topright", c(expression(paste("True Value, ", italic(N)[true])), expression(paste("Observation, ", italic(N)[obs]))), cex=1.4,
       pch=c(NA, 18), lty=c(1, 1), lwd=c(1.5, 1.5),
       col=c("black", "gold"), bty="n")

#step 2: plot the embedding
#nmx<-161:310
nmx<-100+(1:150)
xfull<-datout$obs[nmx]
mxf<-mean(xfull, na.rm=T)
xfull<-xfull/mnx
gfunfull<-splinefun(nmx, xfull)
nsq<-1000
tmsqfull<-seq(min(nmx), max(nmx), length=nsq)
trprdfull<-pmax(0, gfunfull(tmsqfull))
stepsize<-nsq/length(nmx)
xdet<-(sin(tmsqfull/stm)+exp(pars_sim$det[2])+0.5)/2
xdet<-xdet/mnx

ntm2<-trprdfull[1:(nsq-2*stepsize)]
ntm1<-trprdfull[(1+stepsize):(nsq-stepsize)]
nt<-trprdfull[(1+2*stepsize):(nsq)]

psuse<-1:150
etmp<-s_map((datout$obs/mnx)[psuse], E=2, theta=2, silent = TRUE, stats_only = FALSE)
plot(etmp$model_output[[1]]$pred, (datout$true/mnx)[psuse],
     xlim=c(0,2), ylim=c(0, 2))
abline(a=0, b=1, lty=2, lwd=1.5)

mtext(expression(paste("Takens Prediction, ", italic(N)[Takens])), 1, line = 3.2, cex=1.2)
mtext(expression(paste("True Abundance, ", italic(N)[true])), 2, line = 2.4, cex=1.2)
title(main="c.", cex.main=1.6, adj=0.02, line=-1)

#cor(etmp$model_output[[1]]$pred, (datout$true/mnx)[psuse], use = "pairwise")
#1-mean((etmp$model_output[[1]]$pred-(datout$true/mnx)[psuse])^2, na.rm=T)/mean((mean((datout$true/mnx)[psuse],na.rm=T)-(datout$true/mnx)[psuse])^2, na.rm=T)

legend("bottomright", c(expression(paste(rho, " = 0.83")),
                     expression(paste(E[2], " = 0.67"))), cex=1.4,
       col=c("black", "gold"), bty="n")

#step 3: run particle filter
#run filter
nparts<-30  # number of points to sample
tmps<-23:30 # positions to draw
pfout<-particleFilterLL(y = datout$obs[psuse], pars = pars_sim, N = 2e3, detfun = detfun0_sin, dotraceback = TRUE, fulltraceback = TRUE)

#simulate new point at position 8
newpoints<-pfout$fulltracemat[tm[tmps[length(tmps)]],1:nparts]

#plot trajectory
par(mar=c(3,3,6,2))
plot(c(1,15.5), c(0, 2), type="n", xlab="", ylab="", axes=F)
axis(2); axis(1, at=c(2,4,6,8), seq(tmps[1]+1, tmps[length(tmps)], by=2)); box()
sq<-seq(0, 3, by=0.5)
axis(3, at=c(sq*3+8.4), sq)

lines(1:8, datout$true[tm[tmps]], lwd=1.5)
lines(1:(length(tmps)-1), pfout$Nest[tm[tmps[-length(tmps)]]], col="dodgerblue", lwd=1.5)
matlines((length(tmps)-1):length(tmps), rbind(rep(pfout$Nest[tm[tmps[length(tmps)-1]]], nparts), newpoints), type="l",
         col=adjustcolor("dodgerblue", alpha.f = 0.2), lty=1)
points(rep(8, nparts), newpoints, col=adjustcolor("dodgerblue", 0.2), pch=16, cex=0.6)

points(1:8, datout$obs[tm[tmps]], col="gold", pch=18, cex=1.2)
segments(1:8, datout$true[tm[tmps]], 1:8, datout$obs[tm[tmps]], col="gold", lwd=1.5)

#density functions
pr_ntp1_g_proc<-density(pfout$fulltracemat[tm[tmps[length(tmps)]],], bw = 0.15, from = 0, to=2.5)
pr_ntp1_g_obs<-density(obsfun0(so=pars_sim$obs, yt = datout$obs[tm[tmps[length(tmps)]]], inverse=TRUE, N=ncol(pfout$fulltracemat)), bw = 0.15, from = 0, to=2.5)

lines(pr_ntp1_g_proc$y*3+8.4, pr_ntp1_g_proc$x, col="dodgerblue")
lines(pr_ntp1_g_obs$y*3+8.4, pr_ntp1_g_obs$x, col="gold")

conv<-pr_ntp1_g_proc$y*pr_ntp1_g_obs$y
conv<-conv/max(conv)
conv<-conv/sum(mean(diff(pr_ntp1_g_proc$x))*conv)
lines(conv*3+8.4, pr_ntp1_g_proc$x, col="grey40")
abline(h=0, lty=3)

text(9.8, 1.6, expression(paste(italic(N),"(", italic(t), ") ~ ", italic(sigma)[italic(P)])), cex=1.2, col="dodgerblue")
text(11.8, 0.12, expression(paste(italic(N),"(", italic(t), ") ~ ", italic(N)[obs])), cex=1.2, col="gold")
text(14.45, 0.25, expression(paste(italic(N),"(", italic(t), ") ~ {", italic(N)[obs], ", ", italic(sigma)[italic(P)], "}")), cex=1.2, col="grey40")
text(15, 1.0, expression(paste(italic(N)[true],"(", italic(t), ")")), cex=1.2, col="black")


mtext(expression(paste("Time")), 1, line = 2.4, cex=1.2, adj=1/4)
mtext(expression(paste("Abundance")), 2, line = 2.4, cex=1.2)
mtext(expression(paste("Probability Density")), 3, line = 1.8, cex=1.2, adj=0.9)

lines(c(8.2, 22), rep(datout$true[tm[tmps[length(tmps)]]],2), lwd=1.5, lty=2)

abline(v=8.2)

legend(1.2, 2.16, c(expression(paste("Mean Filter Estimate, ", widehat(italic(N)))),
                    expression(paste("Individual Particles"))), cex=1.4,
       pch=c(NA, 16), lty=c(1, 1), lwd=c(1.5, 1),
       col=c("dodgerblue", adjustcolor("dodgerblue", 0.4)), bty="n")
title(main="d.", cex.main=1.6, adj=0.01, line=-1.2)
title(main="e.", cex.main=1.6, adj=0.99, line=-1.2)


#step 4: calculate extinction probability
tmps2<-27:36

par(mar=c(3,3,2,2))
plot(c(1,length(tmps2)), c(0, 2), type="n", xlab="", ylab="", axes=F)
axis(2); axis(1, at=c(2,4,6,8,10), seq(tmps2[1]+1, tmps2[length(tmps2)], length=5))
box()

lines(1:10, datout$true[tm[tmps2]], lwd=1.5)
points(1:10, datout$noproc[tm[tmps2]], col="forestgreen", pch=17, cex=1.2)
segments(1:10, datout$true[tm[tmps2]], 1:10, datout$noproc[tm[tmps2]], col="forestgreen", lwd=1.5)

pnoproc<-indexsort(pfout$fulltracemat_noproc[tm[tmps2]-1,], pfout$fulltraceindex[tm[tmps2]-1,], nsmp = 2e3)

pqtl<-t(apply(pnoproc, 1, function(x) quantile(x, pnorm(-1:1), na.rm=T)))

polygon(c(1:10, 10:1),
        c(pqtl[,1], rev(pqtl[,2])),
        col=adjustcolor("dodgerblue", 0.5))
abline(h=0, lty=3)

legend(1, 2.15, c(expression(paste("Pre-Disturbance State")), expression(paste("Filter Estimate"))), cex=1.4,
       pch=c(19, NA), border = c(NA, 1), fill = c(NA, adjustcolor("dodgerblue", alpha.f = 0.5)), lty=c(1, NA), lwd=c(1.5, 1.5),
       col=c("forestgreen", "dodgerblue"), bty="n")
title(main="f.", cex.main=1.6, adj=0.015, line=-1.2)

mtext(expression(paste("Abundance")), 1, line = 2.4, cex=1.2, adj=1/4)
mtext(expression(paste("Time")), 2, line = 2.4, cex=1.2, adj=1/4)


#mu vs. sd
par(mar=c(3,3,2,4))
musq<-c(0, exp(seq(log(1e-6), log(2), length=1000)))
sdsq<-sqrt(exp(pars_sim$proc[1])*musq^exp(pars_sim$proc[2]))
pmor<-pnorm(0, musq, sdsq)

plot(musq,pmor, type="l", col="sienna", lwd=1.5, log="y", xlab="", ylab="",lty=2, xlim=c(-0.08, 2))

par(new=TRUE)
plot(musq,sdsq, type="l", lwd=1.5, xlab="", ylab="",  axes=F, col="forestgreen", xlim=c(-0.08, 2))
axis(4)

abline(h=0, v=0, lty=3)

legend(0.01, 0.57, c(expression(paste("Pr"[mor])), expression(paste(sigma[italic(P)[italic(tot)]]))), cex=1.4,
       lty=c(2,1), lwd=c(1.5, 1.5),
       col=c("sienna", "forestgreen"), bty="n")
title(main="g.", cex.main=1.6, adj=0.01, line=-1.2)


mtext(expression(paste("Pre-Disturbane Abundance")), 1, line = 2.4, cex=1.2, adj=1/4)
mtext(expression(paste("Process Noise, ", sigma[italic(P)[italic(tot)]])), 4, line = 2.9, cex=1.2, adj=1/4)
mtext(expression(paste("Mortality Probability, Pr"[mor])), 2, line = 2.4, cex=1.2, adj=1/4)



#step 2: plot the embedding (phase space)
sptmp<-scatterplot3d(ntm2,ntm1, nt, label.tick.marks = FALSE, type = "l",
                     xlab = "", ylab="", zlab="", box = FALSE, color=adjustcolor(1, alpha.f = 0.15),
                     mar=c(2,0,8,0),
                     asp = 1.5, scale.y = 0.8, angle = 50)
mtext(expression(paste(italic(N)[obs], "(", italic(t), "-2)")), 1, line = 0.2, adj=0.4)
mtext(expression(paste(italic(N)[obs], "(", italic(t), ")")), 2, line = -2.8, adj=0.35)
mtext(expression(paste(italic(N)[obs], "(", italic(t), "-1)")), 4, line = -3.8, adj=0.3)

tmp<-data.frame(xdet[1:(nsq-2*stepsize)],xdet[(1+stepsize):(nsq-stepsize)], xdet[(1+2*stepsize):(nsq)])
sptmp$points3d(tmp[,1], tmp[,2], tmp[,3], type = "l")
mtext(expression(paste("Phase Space Reconstruction")), 1, line = 2.4, cex=1.2)
title(main="b.", cex.main=1.6, adj=0.18, line=-6.1)


dev.off()
