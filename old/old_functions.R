

#' Extend output from particle filter
#'
#' Function for extending the particles (i.e. extrapolating timeseries into the future), based on output from particleFilterLL function.
#' Function is adapted from the R code of Knape and Valpine (2012), Ecology 93:256-263.
#' @param pfout Output from a previous run of the particleFilterLL function.
#' @param pars A list of parameter values. Must include elements obs (observation error parameters), proc (process noise parameters), and pcol (colonization parameters), which are passed on the their respecive functions, described below. If edmdat=NULL, then element det (deterministic process parameters) must be included.
#' @param Next Length of extended filter - i.e. number of timesteps to predict forward. Defaults to 1e3.
#' @param detfun A function that simulates deterministic dynamics, which takes in arguments sdet (parameters for deterministic model, taken from pars$proc), and xt, observed abundances at time t. Returns estimated abundances at time t+1 based on deterministic function (either a parametric function or an EDM function). Defaults to detfun0.
#' @param procfun A function that simulates process noise, which takes in arguments sp (parameters for process noise function, taken from pars$proc) and xt (abundances prior to process noise). Returns abundances after process noise has occurred. Defaults to procfun0.
#' @param obsfun An observation function, which takes in up to five variables, including so (a vector of parameter values, inherited from pars$obs), yt (a number, showing observed abundance at time t), xt (predicted abundances), binary value "inverse", and number "N". If inverse = TRUE,
#' then function should simulate N draws from the observation function, centered around value yt. If inverse = FALSE, then function should return log probability denisty of observed value yt given predicted values in xt. Defaults to obsfun0.
#' @param colfun A function simulating colonization events, that takes in two arguments: co, a vector of parameter values taken from pars$pcol, and xt, a number or numeric vector of abundances at time t, before colonization has occurred. Returns predicted abundances after colonization has occurred. Defaults to colful0.
#' @param edmdat A list including arguments to be passed to block_lnlp from rEDM package - see block_lnlp help file for details. Can also include optional matrix "extra_columns", a matrix with length(y) rows including extra covariates for attractor reconstruction, which defaults to NULL (i.e. no additional columns).
#' Default for edmdat is NULL, which implies that EDM will not be applied - instead, a detfun and pars$det must be included.
#' @param minval Minimum value, below will observations will be assumed to be zero. Can be helpful for simplex projections, as rounding errors can result in slightly nonzero values. Defaults to 1e-6.
#' @param trimq Either NULL (no limits), or a vector of length 2, including the lower and upper quantiles to keep for extinction and colonization projections. Can be helpful for excluding degenerate particles. Defaults to c(0.01, 0.99).
#' @param randstart Randomize starting positions within the time vector for the new particles? If FALSE, then all particles start at the end of the observed time series. Otherwise particles start at randomly chosen points, which can be helpful for getting better coverage of the full dynamic manifold. Defaults to TRUE.
#' @param concatts Concatenate original and extended particles? Defaults to FALSE.
#' @source Adapted from Knape and Valpine (2012), Ecology 93:256-263.
#' @keywords particle filter, stability, time-series, Taylor power law
#' @return LL, P, rN, x, dem(col, mor)
#' @export

extend_particleFilter = function(pfout, pars, Next = 1e3, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, minval=1e-16, trimq=c(0.01, 0.99), randstart=TRUE, concatts=FALSE) {
  #Adapted from Knape and Valpine (2012), Ecology 93:256-263.

  #extract parameters
  so = pars$obs             # observation error
  sp = pars$proc            # process noise
  sdet = c(pars$det)        # deterministic function
  co<-pars$pcol             # colonization

  #matrices for storing simulations
  N = nrow(pfout$x)                     # Number of particles
  n = ncol(pfout$x)                     # Length of the time series

  xsort = matrix(rep(0, N * n), c(N, n))      # A matrix containing the past particles
  xnew = matrix(rep(0, N * Next), c(N, Next)) # A matrix containing the future estimates

  x = pfout$x
  ind = pfout$ind                       # Index for the previous position of each particle

  #sort particles
  itrace = ind[, n]
  xsort[,n] = x[itrace, n]
  for (t in (n - 1):1) {
    itrace = ind[itrace, t]
    xsort[,t] = x[itrace, t]
  }

  if(Next>0) {
    #set up initial timestep
    if(randstart) {
      nps<-sample(x = 1:n, size = N, replace=TRUE)
    } else {
      nps<-rep(n, N)
    }
    if(is.null(edmdat)) {
      xnew[,1]<-xsort[cbind(1:N, nps)]
      tstart<-2
    } else { #If using lagged embeddings, calculate first N positions
      #set defaults
      if(is.null(edmdat$E))
        edmdat$E<-2
      if(is.null(edmdat$lib))
        edmdat$lib<-c(1, n)
      if(is.null(edmdat$norm))
        edmdat$norm<-2
      if(is.null(edmdat$method))
        edmdat$method="simplex"
      if(is.null(edmdat$tp))
        edmdat$tp<-0
      if(is.null(edmdat$num_neighbors))
        edmdat$num_neighbors<-ifelse(edmdat$method=="simplex", "e+1", 0)

      nps[nps<edmdat$E]<-edmdat$E
      for(i in 1:edmdat$E) {
        xnew[, edmdat$E - (i-1)] = xsort[cbind(1:N, nps-(i-1))]
      }

      tstart<-edmdat$E+1

      #make y block
      yuse<-colMeans(xsort)
      yuse[yuse<minval]<-0
      yblock<-as.matrix(make_block(yuse, max_lag=edmdat$E+1, lib = edmdat$lib)[,-1])

      if(!is.null(edmdat$extra_columns)) { #additional predictors, if applicable
        yblock<-cbind(yblock, edmdat$extra_columns)
      }
    }

    #run for all timesteps
    for (t in tstart:Next) {
      # Project particles forward
      #colonization
      cps<-xnew[, t - 1]>0 #which are >0?
      xnew[!cps, t] = colfun(co=co, xt=xnew[, t - 1][!cps]) #colonize empty sites

      #deterministic
      if(sum(cps)>0) {
        if(is.null(edmdat$E)) {
          xnew[cps, t] = detfun(sdet = sdet, xt = xnew[, t - 1][cps]) #deterministic change in filled, non-colonized sites
        } else {
          smcps<-sum(cps,na.rm=T)
          xtedm<-cbind(NA, xnew[cps,(t-1):(t-edmdat$E)])

          if(!is.null(edmdat$extra_columns)) { #additional predictors for s-mapping
            xtedm<-cbind(xtedm, unname(edmdat$extra_columns[rep(t-1, nrow(xtedm)),]))
          }

          xnew[cps, t] = detfun(edmdat, xt=xtedm, yblock)
        }
      }
      xnew[, t][xnew[, t]<=minval]<-0

      #process noise
      cps[!cps]<-xnew[!cps, t]>0 #which are >0 now?
      xnew[cps, t] = procfun(sp = sp, xt = xnew[cps, t])
    }
  }
  if(concatts) {
    xnew<-cbind(xsort, xnew)
  }


  #Calculate demographic rates
  text<-1/(rowSums((xnew[,-1]<=minval & xnew[,-ncol(xnew)]>minval))/pmax(1, rowSums(xnew[,-ncol(xnew)]>minval)))
  tcol<-1/(rowSums((xnew[,-1]>minval & xnew[,-ncol(xnew)]<=minval))/pmax(1, rowSums(xnew[,-ncol(xnew)]<=minval)))
  text[!is.finite(text)]<-NA
  tcol[!is.finite(tcol)]<-NA

  #trim extreme values
  if(!is.null(trimq)) {
    tmp<-quantile(text, trimq, na.rm=T)
    text[text<tmp[1]]<-NA
    text[text>tmp[2]]<-NA

    tmp<-quantile(tcol, trimq, na.rm=T)
    tcol[tcol<tmp[1]]<-NA
    tcol[tcol>tmp[2]]<-NA
  }

  return(list(xnew=xnew, xsort=xsort, demdat=list(text=text, tcol=tcol)))
}

#' Abundance density kernel estimation
#'
#' Estimates the probability density function associated with an abundance distribution
#' @param x A vector of abundances
#' @param from Lower boundary for the kernel - defaults to 0.
#' @param bw Bandwidth for kernel - defauls to 0.1.
#' @param breaks Number of breaks for a comparative histogram - defaults to 20.
#' @param doplot Plot a comparative histogram and the density kernel? Defaults to TRUE.
#' @param excludezero Should zero values be ignored? Defaults to TRUE.
#' @return A fitted density kernel.
#' @export

abund_density<-function(x, from = 0, bw = 0.1, breaks = 20, doplot = TRUE, excludezero=TRUE, ...) {
  if(excludezero) {
    x<-x[x>0]
  }

  df<-density(c(x), from = 0, bw = bw, ...)

  if(doplot) {
    hist(x, breaks = breaks, xlab="abundance", probability = TRUE, main="")
    lines(df$x, df$y, lwd=2)
  }

  return(df)
}

#' Density extinction estimate.
#'
#' Estimate extinction rates based on the results from a call of  the abund_density function.
#' @param df Density kernel resulting from a call of the 'abund_density' function.
#' @param sp Parameters to be handed to the 'procfun' function.
#' @param tlength Length of timeseries to simulate.
#' @param niter Number of iterations from which to draw the estimate.
#' @param procfun Process noise function. Defaults to procfun0.
#' @return Mean extinction rate, standard deviation of rate, and list of niter estimates.
#' @export

est_mor<-function(df, sp, tlength = 100, niter = 1e3, procfun = procfun0) {
  mest<-numeric(niter)

  if(niter==1) {
    mest<-sum(procfun(sp, df$x, inverse = TRUE)*(df$y/sum(df$y)))
  }
  for(i in 1:niter) {
    smp <- sample(length(df$x), tlength, prob = df$y, replace = TRUE)
    mest[i] <- mean(procfun(sp, df$x[smp], inverse = TRUE))
  }

  return(list(mu = mean(mest), sd=sd(mest), mest=mest))
}






#' Approximate Bayesian Computation-based optimizer
#'
#' Attempts to find maximum likelihood parameter estimates for parameter set p0 and
#' associated models, given observed time-series y. Can be used to either run a new
#' optimization routine, or to extend an existing run to include more iterations.
#' @param y Observed values, against which model likelihood will be estimated. No default is included, though the value can be taken from 'oldrun' if provided.
#' @param p0 Initial mean values for each parameter to be estimated. Used to sample first sets of values for the optimizer, and to parameterize the prior. No default is included, though the value can be taken from 'oldrun' if provided.
#' @param sd0 Initial standard deviation for each parameter in 'p0'. Used to sample first sets of values for the optimizer, and to parameterize the prior. Defaults to 1 for every parameter.
#' @param N Number of particle to use for the particle filter - passed to 'likelihood' function. Defaults to 1e3.
#' @param niter_optim Number of optimizer iterations to perform. Defaults to 20. Note, iterations can be extended post-hoc.
#' @param niter_pars Number of parameter sets to test per iteration. Defaults to 100.
#' @param fretain Optimizer retains fraction 'fretain' parameter sets with the highest likelihoods in each iteration, and uses these to update the values passed to 'sampler_fun'. Defaults to 0.1.
#' @param burnin Number of iterations to burn-in, during which initial values p0 and sd0 are passed to 'the 'sampler_fun'. Defaults to 2.
#' @param priorLL Should the likelihood of the parameters under the prior be added to the likelihood used for optimization? Can be useful when fitting models with hyperparameters, or if there is high uncertainty. Defaults to FALSE.
#' @param oldrun The output of a previous run of the run_ABC_optim function. If included, then this funciton extends that run. Any parameters with user-defined values NULL will take their values from this previous run, with two expetions:
#' burnin defaults to zero, and niter_pars will always be taken from the previous run, regardless of user input.
#' @param sampler_fun A function that can sample new values for the optimizer. Must be able to take two types of input -
#' either a vector of means and standard deviations, or a vector of means and a covariance matrix. See sampler_fun0 (default) for details.
#' @param parseparam A function for parsing a vector of parameters into a form that is readable by the likelihood function. Defaults to parseparam0.
#' @param likelihood A function that returns log likelihood estimates of the observed data given some model (e.g. a particle filter). See likelihood0 (default) for examples.
#' @param density_fun A function that returns log likelihood estimates for parameter p0 given prior estimates of mean and standard deviation. Defaults to density_fun0. Is used to calculate likelihood under the prior, if requested.
#' @param bootstrap_subs Should the 'fretain' best observations be bootstrapped (with replacement) before moving to the next iteration? Can help prevent outliers from dominating optimization. Defaults to TRUE.
#' @param report_interval How frequently should progress be reported? Defaults to the minimum of every 10 iterations, or 1/10 of the total number of iterations.
#' @param silent If FALSE (default), then progress is reported.
#' @keywords ABC optimization
#' @return Returns a list, including mean parameter and likelihood values
#' for each iteration (parsmean and LLmean), an apporoximation of the
#' covariance matrix relating parameters (parscv), the full information for
#' likelihood and parameter estimates for each iteration (LLout and parout),
#' and a list of the parameters used to run the model (runstats).
#'
#' @import mvtnorm
#' @export

run_ABC_optim<-function(y=NULL,
                        p0=NULL, sd0=NULL, N=NULL,
                        niter_optim=NULL, niter_pars=NULL,
                        fretain=NULL, burnin=NULL, priorLL=NULL,
                        oldrun = NULL,
                        sampler_fun=NULL, parseparam=NULL, likelihood=NULL, density_fun=NULL,
                        bootstrap_subs=NULL,
                        report_interval=pmin(round(niter_optim/10), 10),
                        silent=FALSE) {

  ## get parameter values from old run, if applicable, or enter default values
  if(is.null(oldrun)) {
    if(is.null(y) | is.null(p0)) {
      return("error: y and p0 must be provided, either directly or via variable 'oldrun'")
    }
  } else {
    if(is.null(y)) {
      y<-oldrun$runstats$y
    }
    if(is.null(p0)) {
      p0<-oldrun$runstats$p0
    }
  }
  npars<-length(p0) #number of parameters

  if(is.null(sd0)) {
    if(is.null(oldrun)) {
      sd0<-rep(1, length(p0))
    } else {
      sd0<-oldrun$runstats$sd0
    }
  }
  if(is.null(N)) {
    if(is.null(oldrun)) {
      N<-1e3
    } else {
      N<-oldrun$runstats$N
    }
  }
  if(is.null(niter_optim)) {
    if(is.null(oldrun)) {
      niter_optim<-20
    } else {
      niter_optim<-oldrun$runstats$niter_optim
    }
  }
  if(is.null(oldrun)) {
    if(is.null(niter_pars)) {
      niter_pars<-100
    }
  } else {
    #niter_pars must match that from previous run
    niter_pars<-oldrun$runstats$niter_pars
  }
  if(is.null(fretain)) {
    if(is.null(oldrun)) {
      fretain<-0.1
    } else {
      fretain<-oldrun$runstats$fretain
    }
  }
  if(is.null(burnin)) {
    if(is.null(oldrun)) {
      burnin<-2
    } else {
      #defaults to no burnin if extending an old run
      burnin<-0
    }
  }
  if(is.null(priorLL)) {
    if(is.null(oldrun)) {
      priorLL<-FALSE
    } else {
      priorLL<-oldrun$runstats$priorLL
    }
  }
  if(is.null(sampler_fun)) {
    if(is.null(oldrun)) {
      sampler_fun<-sampler_fun0
    } else {
      sampler_fun<-oldrun$runstats$sampler_fun
    }
  }
  if(is.null(parseparam)) {
    if(is.null(oldrun)) {
      parseparam<-parseparam0
    } else {
      parseparam<-oldrun$runstats$parseparam
    }
  }
  if(is.null(likelihood)) {
    if(is.null(oldrun)) {
      likelihood<-likelihood0
    } else {
      likelihood<-oldrun$runstats$likelihood
    }
  }
  if(is.null(density_fun)) {
    if(is.null(oldrun)) {
      density_fun<-density_fun0
    } else {
      density_fun<-oldrun$runstats$density_fun
    }
  }
  if(is.null(bootstrap_subs)) {
    if(is.null(oldrun)) {
      bootstrap_subs<-TRUE
    } else {
      bootstrap_subs<-oldrun$runstats$bootstrap_subs
    }
  }

  ## set up output matrices
  #initial draw of parameters
  if(is.null(oldrun)) {
    prl<-(sampler_fun(n=niter_pars, pars=parseparam(p0), priorsd = sd0))
  } else {
    prl<-sampler_fun(niter_pars, oldrun$parsmean[nrow(oldrun$parsmean),], oldrun$parscv[nrow(oldrun$parsmean),,])
  }
  #likelihood
  LLout<-matrix(nrow=niter_optim, ncol=niter_pars)
  if(!is.null(oldrun)) {
    LLout<-rbind(oldrun$LLout, LLout)
  }
  #parameters
  if(!is.null(oldrun)) {
    parout<-array(dim=c(niter_optim+oldrun$runstats$niter_optim, niter_pars, npars))
    parout[1:oldrun$runstats$niter_optim,,]<-oldrun$parout
  } else {
    parout<-array(dim=c(niter_optim, niter_pars, npars))
  }
  #mean outputs
  LLmean<-numeric(niter_optim)
  if(!is.null(oldrun)) {
    LLmean<-c(oldrun$LLmean, LLmean)
  }
  parsmean<-matrix(nrow=niter_optim, ncol=npars)
  if(!is.null(oldrun)) {
    parsmean<-rbind(oldrun$parsmean, parsmean)
  }

  #covariance matrix
  if(!is.null(oldrun)) {
    parscv<-array(dim=c(niter_optim+oldrun$runstats$niter_optim, npars, npars))
    parscv[1:oldrun$runstats$niter_optim,,]<-oldrun$parscv
  } else {
    parscv<-array(dim=c(niter_optim,npars, npars))
  }

  ## find optimal parameters
  if(!is.null(oldrun)) {
    ifrom<-oldrun$runstats$niter_optim+1
    ito<-ifrom+niter_optim-1
    burnin_use<-burnin+ifrom-1
  } else {
    ifrom<-1
    ito<-niter_optim
    burnin_use<-burnin
  }

  for(i in ifrom:ito) {
    #run filter for each parameter set
    out<-apply(prl, 1, function(x) likelihood(x, y=y, parseparam = parseparam, N=1e3))

    #save outputs
    LLout[i,]<-out
    parout[i,,]<-prl

    #include prior in likelihood calculation?
    if(priorLL) {
      #likelihood of parameters under the prior
      LL_prior<-apply(prl, 1, function(x) density_fun(x, parseparam(p0), sd0))
      LLout[i,]<-LLout[i,]+LL_prior
    }

    #get weights for fraction "fretain" best parameter sets
    sbs<-which(LLout[1:i,]>quantile(LLout[1:i,], 1-fretain))
    if(bootstrap_subs) {
      sbs<-sample(sbs, length(sbs), replace = TRUE)
    }

    pd<-exp(LLout[1:i,][sbs]-max(LLout[1:i,][sbs]))/sum(exp(LLout[1:i,][sbs]-max(LLout[1:i,][sbs]))) #calculate weights

    #calculate weighted mean and variance/covariance
    mu_est<-numeric(npars)
    cvmat_est<-matrix(nrow=npars, ncol=npars)
    for(ii in 1:npars) {
      mu_est[ii]<-sum(parout[1:i,,ii][sbs]*pd)
      cvmat_est[ii,ii]<-sum(((parout[1:i,,ii][sbs]-mu_est[ii])^2)*pd*(1-sum(pd^2)/(sum(pd)^2)))
      for(jj in (ii+1):npars) {
        if(jj<=npars) {
          cvmat_est[ii,jj]<-cvmat_est[jj,ii]<-sum(((parout[1:i,,ii][sbs]-mu_est[ii])*
                                                     (parout[1:i,,jj][sbs]-mu_est[jj]))*pd)*
            (1/(1-sum(pd^2)))
        }
      }
    }

    #generate new parameter set
    if(i>=burnin_use) {
      prl<-sampler_fun(niter_pars, mu_est, cvmat_est)
    } else {
      prl<-sampler_fun(n=niter_pars, pars=parseparam(p0))
    }

    #save outputs
    parsmean[i,]<-mu_est
    parscv[i,,]<-cvmat_est
    LLmean[i]<-mean(LLout[1:i,][sbs])

    #report progress
    if(!silent) {
      if(i==ifrom) {
        cat("progress:")
      }
      if(i/report_interval == floor(i/report_interval)) {
        cat(paste("..", 100*round((i-ifrom+1)/(ito-ifrom+1),ceiling(log10(niter_optim))), "%", sep=""))
      }
    }
  }

  #save original burnin value
  if(!is.null(oldrun)) {
    burnin<-oldrun$runstats$burnin
  }

  #return output
  return(list(#summary data
    parsmean=parsmean, parscv=parscv, LLmean=LLmean,
    #full data
    LLout=LLout, parout=parout,
    #data needed for extending a run
    runstats=list(y=y,p0=p0, sd0=sd0,N=N,
                  niter_optim=niter_optim, niter_pars=niter_pars,
                  fretain=fretain, burnin=burnin, priorLL=priorLL,
                  sampler_fun=sampler_fun,
                  parseparam=parseparam,
                  likelihood=likelihood,
                  density_fun=density_fun,
                  bootstrap_subs=bootstrap_subs)))
}


#' Plot likelihoods from optimization
#'
#' Plots the change in log likelihood across simulations based on an output from the run_ABC_optim function.
#' @param optout Output from the run_ABC_optim function.
#' @param logx Should the x-axis (number of iterations) be plotted in log space? Defaults to FALSE.
#' @keywords ABC optimization
#' @return a plot
#' @export

plot_abc_likelihood<-function(optout, logx=FALSE) {
  if(logx) {
    tmp<-"x"
  } else {
    tmp<-""
  }
  plot(optout$LLmean, type="l", xlab="iteration", ylab="LL", log = tmp, lwd=2)
  abline(v=optout$runstats$burnin, lty=2)
}


#' Plot rmse of parameter estimates
#'
#' Plots the root mean square error of the optimizer over iterations, if true parameter values are known a prior (e.g. to test settings on simulated data).
#' @param optout Output from the run_ABC_optim function.
#' @param param_true A vector containing the "true" parameter values.
#' @param doLL Should rmse also be plotted alongside log likelihood? Defaults to TRUE.
#' @keywords ABC optimization
#' @return a plot
#' @export

plot_abc_rmse<-function(optout, param_true, doLL=TRUE) {
  rmse<-colSums((t(optout$parsmean)-param_true)^2)
  plot(rmse, type="l", xlab="iteration", ylab="rmse", lwd=2)

  if(doLL) {
    plot(rmse, optout$LLmean, xlab="rmse", ylab="LL", log="x")
    mod<-loess(optout$LLmean~log(rmse), enp.target = 5)
    sq<-seq(min(rmse), max(rmse), length=1000)
    lines(sq, predict(mod, newdata=data.frame(rmse=sq)), lwd=2)
  }
}


#' Plot parameter estimates from optimization
#'
#' Plots the change in parameter estimates over the course of optimization. Shows mean, confidence interval, and (optionally) positions of true values and priors.
#' @param optout Output from the run_ABC_optim function.
#' @param ifun A list of inverse functions to be applied to the parameters for plotting - defaults to no tranformation.
#' @param param0 An optional vector of prior values for the paramters, to be plotted in blue.
#' @param param_true An optional vector of true values for parameters, to be plotted in red.
#' @param alpha Confidence interval size to be plotted for the parameter estimates. Defaults to 0.95.
#' @keywords ABC optimization
#' @return a plot
#' @export

plot_abc_params<-function(optout, ifun=list(function(x) {x}), param0=NULL, param_true=NULL, alpha=0.95) {
  np<-length(optout$runstats$p0)
  if(length(ifun)==1 && np > 1) {
    for(i in 2:np) {
      ifun[[i]]<-ifun[[1]]
    }
  }

  for(i in 1:np) {
    xv<-optout$parsmean[,i]+qnorm(1-(1-alpha)/2)*cbind(-sqrt(optout$parscv[,i,i]), 0, sqrt(optout$parscv[,i,i]))
    xv[!is.finite(xv)]<-NA

    matplot(1:nrow(optout$parsmean), ifun[[i]](xv),
            type="l", col=1, lty=c(2,1,2),
            xlab="iterations", ylab=names(param0[i]),
            ylim=ifun[[i]](range(c(xv, param_true[i], param0[i]), na.rm=T)), lwd=c(1,2,1))
    if(!is.null(param_true)) {
      abline(h=ifun[[i]](param_true[i]), lty=3, col=2)
    }
    if(!is.null(param0)) {
      abline(h=ifun[[i]](param0[i]), lty=3, col=4)
    }
  }
}


#' Estimate parameter distribution from likelihood density function.
#'
#' Extracts (and optionally plots) information from the density function generated across runs in the run_ABC_optim function.
#' Calculates estimates of the mean and standard deviation of parameters from this distribution, based on weighted averaging
#' @param optout Output from the run_ABC_optim function.
#' @param param0 An optional vector of prior values for the paramters, to be plotted in blue.
#' @param param_true An optional vector of true values for parameters, to be plotted in red.
#' @param fretain Retain fraction 'fretain' parameter estiamtes with the highest likelihoods from the optimization process for subsequent calculations. Defaults to value used for optimizer call.
#' @param nbootstrap Number of bootstrapped iterations to use for estimating density kernel. Defaults to 1 (i.e. use full sample).
#' @param df_smooth Smoothing parameter, passed to smooth.spline function for estimating density function. Defaults to 5. Note - numbers closer to zero generally yield smoother, but less detailed, estimates.
#' @param sm_steps Number of steps to estimate from the smoother for calculating means and standard deviations. Defaults to 1000.
#' @param pltnames Names to be included in plots. Defaults to names(param0).
#' @param nobs Scaling factor for the likelihoods used to calculate the mean and standard deviation - defaults to 1.
#' @param doplot Should the distributions be plotted? Defaults to TRUE.
#' @keywords ABC optimization
#' @return A list including mean parameter values (muest) and standard deviations for parameter values (sdest) calculated from the density function. These
#' are probably the best estimates to report and use for subsequent analysis, as the covariance outputs from the optimizer are generally too small as a result
#' of the range of space that is sampled by 'sampler_fun'. Optionally also plots the density point estimates and kernel smoother -
#' note that smoothed function is plotted higher above the points than it actually appears, for better visualization.
#' Also returns full results from smoothing, including the mean smoothed estimates (mupred), and array with all smoothed estimates generated across bootstrapping
#' iterations (predarray), and an array of the values for which smoothed estimates were generated (sqarray).
#' @export

abc_densities<-function(optout, param0=NULL, param_true=NULL, fretain=NULL, nbootstrap=NULL, df_smooth=5, sm_steps=1000, pltnames=names(param0), nobs=1, doplot=TRUE) {
  if(is.null(fretain)) {
    fretain<-optout$runstats$fretain
  }
  if(is.null(nbootstrap)) {
    nbootstrap<-1
  }
  nparm<-length(optout$runstats$p0)
  muest<-numeric(nparm)
  sdest<-numeric(nparm)

  if(length(df_smooth)==1) {
    df_smooth<-rep(df_smooth, nparm)
  }

  mupred<-array(dim=c(nparm, sm_steps, 2))
  predarray<-array(dim=c(nparm, sm_steps, nbootstrap))
  sqarray<-array(dim=c(nparm, sm_steps))

  for(i in 1:nparm) {
    tmp<-optout$parout[,,i]
    ps<-which(optout$LLout>quantile(optout$LLout,1-fretain))

    if(doplot) {
      plot(tmp[ps], c(optout$LLout)[ps],
           main=pltnames[i],
           xlab="parameter", ylab="LL", cex=0.5,
           xlim=range(c(tmp[ps], param0[i], param_true[i])))
      if(!is.null(param0)) {
        abline(v=c(param0[i]), lty=2, lwd=2, col=4)
      }
      if(!is.null(param_true)) {
        abline(v=c(param_true[i]), lty=2, lwd=2, col=2)
      }
    }

    sqarray[i,]<-seq(min(tmp[ps]), max(tmp[ps]), length=sm_steps)

    for(j in 1:nbootstrap) {
      if(nbootstrap==1) {
        ps2<-1:length(ps)
      } else {
        ps2<-sample(length(ps), replace = TRUE)
      }

      xtmp<-tmp[ps[ps2]]
      ytmp<-optout$LLout[ps[ps2]]

      spl<-smooth.spline(xtmp, ytmp, df = df_smooth[i])
      predarray[i,,j]<-predict(spl, sqarray[i,])$y
    }

    if(nbootstrap>1) {
      mupred[i,,]<-t(apply(predarray[i,,],1,function(x) cbind(mean(x,na.rm=T), sd(x,na.rm=T))))
    } else {
      mupred[i,,1]<-predarray[i,,1]
      mupred[i,,2]<-0
    }

    if(doplot) {
      matlines(sqarray[i,], cbind(mupred[i,,1], mupred[i,,1]+mupred[i,,2], mupred[i,,1]-mupred[i,,2]),
               col=3, lwd=c(2,1,1), lty=c(1,2,2),
               xlim=range(c(tmp[ps], param0[i], param_true[i])),
               axes=F, xlab="", ylab="", type="l")
    }


    wts<-exp(mupred[i,,1]*nobs-max(mupred[i,,1]*nobs))/sum(exp(mupred[i,,1]*nobs-max(mupred[i,,1]*nobs)))

    muest[i]<-sum(seq(min(tmp[ps]), max(tmp[ps]), length=sm_steps)*wts)
    sdest[i]<-sqrt(sum(((seq(min(tmp[ps]), max(tmp[ps]), length=sm_steps)-muest[i])^2)*wts*(1-sum(wts^2)/(sum(wts)^2))))
  }

  return(list(muest=muest, sdest=sdest, mupred=mupred, predarray=predarray, sqarray=sqarray))
}


#' Coerce p-value to fall into lower tail.
#'
#' Returns either x, or 1-x, whichever is smaller.
#' @param x The value to be transformed
#' @keywords helper function
#' @return The lower-tail probability estimate
#' @export

lowertail<-function(x) {
  pmin(x, 1-x)
}


