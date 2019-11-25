plotfun<-function(plotvar, byvar, summarydat, cutlst, mrowpar, collst, cutoff=0.4, doci=TRUE, ...) {
  if(plotvar=="proc") {
    detv<-summarydat$det_proc_mu
    edmv<-summarydat$edm_proc_mu

    detq1<-summarydat$det_proc_qt.1
    detq2<-summarydat$det_proc_qt.2
    detq3<-summarydat$det_proc_qt.3
    detq4<-summarydat$det_proc_qt.4
    detq5<-summarydat$det_proc_qt.5

    edmq1<-summarydat$edm_proc_qt.1
    edmq2<-summarydat$edm_proc_qt.2
    edmq3<-summarydat$edm_proc_qt.3
    edmq4<-summarydat$edm_proc_qt.4
    edmq5<-summarydat$edm_proc_qt.5

    xv<-summarydat$proc
  } else if(plotvar=="obs") {
    detv<-summarydat$det_obs_mu
    edmv<-summarydat$edm_obs_mu

    detq1<-summarydat$det_obs_qt.1
    detq2<-summarydat$det_obs_qt.2
    detq3<-summarydat$det_obs_qt.3
    detq4<-summarydat$det_obs_qt.4
    detq5<-summarydat$det_obs_qt.5

    edmq1<-summarydat$edm_obs_qt.1
    edmq2<-summarydat$edm_obs_qt.2
    edmq3<-summarydat$edm_obs_qt.3
    edmq4<-summarydat$edm_obs_qt.4
    edmq5<-summarydat$edm_obs_qt.5

    xv<-summarydat$obs
  } else if(plotvar=="fit") {
    detv<-summarydat$rmse_det
    edmv<-summarydat$rmse_edm
    v0<-summarydat$rmse0

    detq1<-NA
    detq2<-NA
    detq3<-NA
    detq4<-NA
    detq5<-NA

    edmq1<-NA
    edmq2<-NA
    edmq3<-NA
    edmq4<-NA
    edmq5<-NA

    xv<-summarydat$obs
  } else if(plotvar=="det") {
    detv<-summarydat$proc_det_mu
    edmv<-summarydat$proc_edm_mu
    v0<-summarydat$proc0_mu

    detq1<-NA
    detq2<-NA
    detq3<-NA
    detq4<-NA
    detq5<-NA

    edmq1<-NA
    edmq2<-NA
    edmq3<-NA
    edmq4<-NA
    edmq5<-NA

    xv<-summarydat$obs
  }

  detv[summarydat$gelmandet>1.1]<-NA
  edmv[summarydat$gelmanedm>1.1]<-NA
  detq1[summarydat$gelmandet>1.1]<-NA
  edmq1[summarydat$gelmanedm>1.1]<-NA

  if(byvar=="obs") {
    bylist<-obscutlst
    byall<-summarydat$obsccut
  } else if(byvar=="proc") {
    bylist<-proccutlst
    byall<-summarydat$proccut
  }

  pchtype<-summarydat$edm_var<=cutoff

  par(mfrow=mrowpar, mar=c(2,2,2,2), oma=c(2,2,0,0))
  for(i in 1:length(bylist)) {
    ps<-which(byall==bylist[i])
    sq<-seq(0, 1, length=1000)

    if(plotvar%in%c("obs", "proc")) {
      pmat<-cbind(detv[ps],edmv[ps])
      cuse<-collst[-1]
    } else {
      pmat<-cbind(v0[ps], detv[ps],edmv[ps])
      cuse<-collst
    }
    matplot(xv[ps],
            pmat,
            xlab="", ylab="",
            type="n", main=bylist[i],
            log="xy", ...)
    for(ii in 1:ncol(pmat)) {
      points(xv[ps],
           pmat[,ii], col=cuse[ii], cex=0.8, pch=c(4,1)[pchtype[ps]+1])
    }
    if(plotvar%in%c("obs", "proc") & doci==TRUE) {
      segments(xv[ps], detq2[ps], xv[ps], detq4[ps], col=collst[2], lend=4)
      segments(xv[ps], edmq2[ps], xv[ps], edmq4[ps], col=collst[3], lend=4)
      wts1<-((detq3[ps]-detq2[ps])+(detq4[ps]-detq3[ps]))/2
      wts2<-((edmq3[ps]-edmq2[ps])+(edmq4[ps]-edmq3[ps]))/2
    } else {
      wts1<-rep(1, length(ps))
      wts2<-rep(1, length(ps))
    }
    if(plotvar!="det") {
      abline(a=0, b=1, lty=3)
      abline(v=cutoff, lty=3)
    } else {
      abline(v=c(0.1, cutoff), lty=3)
    }

    if(sum(summarydat$edm_var[ps]<cutoff)>0 & plotvar%in%c("obs", "proc")) {
      pscut<-which(summarydat$edm_var[ps]<cutoff)
      r2<-apply(pmat[pscut,], 2, function(x) cor(log(x), log(xv[ps][pscut]), use="pairwise")^2)
      legend("topleft", legend = round(r2, 2), col=adjustcolor(cuse,alpha.f = 2), pch=1, bty="n", title = expression(paste(R^2)))
    }

    ps1<-which(is.finite(wts1) & wts1>0)
    ps2<-which(is.finite(wts2) & wts2>0)
    xtmp<-xv[ps][ps1]
    ytmp<-detv[ps][ps1]
    mod0<-loess(ytmp~xtmp, enp.target = 2, weights = 1/wts1[ps1])
    xtmp<-xv[ps][ps2]
    ytmp<-edmv[ps][ps2]
    mod1<-loess(ytmp~xtmp, enp.target = 2, weights = 1/wts2[ps2])

    prd1<-predict(mod0, newdat=data.frame(xtmp=sq), se=TRUE)
    prd2<-predict(mod1, newdat=data.frame(xtmp=sq), se=TRUE)
    ps1<-is.finite(prd1$fit)
    ps2<-is.finite(prd2$fit)

    polygon(c(sq[ps1], rev(sq[ps1])), c(prd1$fit[ps1]-prd1$se.fit[ps1], rev(prd1$fit[ps1]+prd1$se.fit[ps1])), col=collst[2], border = NA)
    polygon(c(sq[ps2], rev(sq[ps2])), c(prd2$fit[ps2]-prd2$se.fit[ps2], rev(prd2$fit[ps2]+prd2$se.fit[ps2])), col=collst[3], border = NA)

    if(plotvar%in%c("fit", "det")) {
      xtmp<-xv[ps]
      ytmp<-v0[ps]
      mod3<-loess(ytmp~xtmp, enp.target = 2)

      prd3<-predict(mod3, newdat=data.frame(xtmp=sq), se=TRUE)
      ps3<-is.finite(prd1$fit)

      polygon(c(sq[ps3], rev(sq[ps3])), c(prd3$fit[ps3]-prd3$se.fit[ps3], rev(prd3$fit[ps3]+prd3$se.fit[ps3])), col=collst[1], border = NA)
    }
    if(plotvar=="det") {
      mutmp<-mean(summarydat$proc_true_mu[ps], na.rm=T)
      sdtmp<-sd(summarydat$proc_true_mu[ps], na.rm=T)

      polygon(c(min(sq[sq>0]), max(sq)*2, max(sq)*2, min(sq[sq>0])),
              c(mutmp+sdtmp, mutmp+sdtmp, mutmp-sdtmp, mutmp-sdtmp),
              col=collst[4], border=NA)
    }
  }

  if(plotvar=="fit") {
    mtext("rmse", 2, outer=TRUE, line=0.6)
    mtext("obs error", 1, outer=TRUE, line=0.6)
  } else if(plotvar=="det") {
    mtext("deterministic variation", 2, outer=TRUE, line=0.6)
    mtext("obs error", 1, outer=TRUE, line=0.6)
  } else {
    mtext(paste(plotvar, "error, predicted"), 2, outer=TRUE, line=0.6)
    mtext(paste(plotvar, "error, observed"), 1, outer=TRUE, line=0.6)
  }
}
