#' default deterministic function
#'
#' Simulates deterministic component of Ricker model, of the form xt+1 = xt exp(exp(sdet[1])*(1-xt/exp(sdet[2])))
#' @param sdet a numeric vector of length two, specifying growth rate and carrying capacity
#' @param xt a number or numeric vector of abundances at time t
#' @param time the timestep - defaults to NULL (i.e. not used)
#' @keywords deterministic function, discrete-time model, time-series, fake data
#' @return a number or numeric vector of length xt, with predicted abundances at time t+1
#' @export

detfun0<-function(sdet, xt, time=NULL) {
  xt = xt*exp(exp(sdet[1])*(1-xt/exp(sdet[2])))
  return(xt)
}

#' REDM deterministic function
#'
#' Estimates future states of xt based on based benavior
#' @param edmdat A list of paramter values, which is passed on to the block_lnlp function. See block_lnlp for more details.
#' @param xt A matrix of dynamics for which abundances at time t+1 are to be predicted. Note that first column is assumed to represent t+1.
#' @param yblock A matrix of observed values, to be used as training data. Note that first column is assumed to represent t+1.
#' @param tm Time corresponding to position of xt in y, used to mask observations from training data. Defaults to NULL, meaning that no observations are masked.
#' @keywords applies EDM from the rEDM package to reconstruct deterministic time series dynamics
#' @return a number or numeric vector of length xt, with predicted abundances at time t+1
#' @import rEDM
#' @source @source Adapted from Ye, Sugihara, et al. (2015), PNAS 112:E1569-E1576.
#' @export

EDMfun0<-function(smp_cf, yp, x, minest=0, time) {
  nD<-ncol(smp_cf)
  if(nD>2) {
    out<-sum(smp_cf[time,-1]*c(yp[(time-(nD-2)):(time-1)], 1))+smp_cf[time,1]*x
  } else {
    out<-smp_cf[time,-1]+smp_cf[time,1]*x
  }
  out[out<minest]<-minest
  out
}

#' default process noise function
#'
#' Simulates effects of process noise following a Gaussian perturbation.
#' Note that process noise only influences positive abundances (i.e. process noise cannot contribute to colonization)
#' @param sp a numeric vector of length one, specifying the log-transformed standard deviation of the process noise function
#' @param xt a number or numeric vector of abundances at time t, before process noise has occurred
#' @param inverse a logical specifying whether the inverse (i.e. probability of drawing a value of zero given xt and sp) should be calcualted
#' @param time the timestep - defaults to NULL (i.e. not used)
#' @keywords process noise
#' @return a number or numeric vector of length xt, with predicted abundances after process noise has occurred
#' @import stats
#' @export

procfun0<-function(sp, xt, inverse = FALSE, time=NULL) {
  if(length(sp)==1) {
    std_tmp<-exp(sp[1])
  } else {
    std_tmp<-exp(sp[1])*xt^exp(sp[2])
  }

  if(!inverse) {
    sm<-length(xt)
    xt = pmax(0, xt + rnorm(sm, 0, std_tmp))
    return(xt)
  } else {
    pd <- pnorm(0, xt, std_tmp)
    return(pd)
  }
}


#' default observation noise function
#'
#' Two options: If inverse=FALSE, calculates the log probability density of observation yt based on true state xt and observation error.
#' Otherwise, simulates N random observations of yt.
#' Observation error follows a Gaussian distribution truncated at zero, using a Tobit distribution.
#' Note that probability density is calculated based on a Tobit distribution, with lower boundary zero.
#' @param so a numeric vector of length one, specifying log-transformed standard deviation of the observation error
#' @param yt a number, representing a potential observed value of xt
#' @param xt a number or numeric vector of "true" (or simulated) abundances at time t, from which the likelihood of yt will be calculated - defaults to NULL for inverse=TRUE
#' @param inverse a logical specifying whether inverse (i.e. random number generator) function should be implemented - defaults to FALSE
#' @param N number of draws from the random number generator, if inverse=TRUE - defaults to NULL
#' @param time the timestep - defaults to NULL (i.e. not used)
#' @keywords observation error
#' @return If inverse=FALSE, a number or numeric vector of length xt, with predicted log likelihoods of observation yt. If inverse = FALSE, returns N random draws from the observation function.
#' @import stats
#' @export

obsfun0<-function(so, yt, xt=NULL, inverse=FALSE, N=NULL, minsd=0.01, time=NULL) {
  if(inverse) {
    std_tmp<-exp(so[1])*yt
    std_tmp[std_tmp<minsd]<-minsd
    pmax(0, rnorm(n = N, mean = yt, sd = std_tmp))
  } else {
    std_tmp<-exp(so[1])*xt
    std_tmp[std_tmp<minsd]<-minsd

    ps<-(xt==0)
    LL<-numeric(length(xt))
    #Tobit distribution:
    LL[ps]<-pnorm(0, mean=yt/std_tmp[ps], log.p = TRUE)
    LL[!ps]<-dnorm((yt-xt[!ps])/std_tmp[!ps],log=TRUE)-log(std_tmp[!ps])
    LL
  }
}


#' default colonization function
#'
#' Simulates colonization events - events occur as a binomial random process with probability ilogit(p), and populations are seeded with abundance exp(A).
#' @param co a numeric vector of length two (p, A), specifying the logit-transformed colonization probability when abundance is zero, and the log-transformed abundance observed immediately after a colonization event
#' @param xt a number or numeric vector of abundances at time t, before colonization has occurred
#' @keywords colonization
#' @return a numeric, including number or numeric vector of length xt, with predicted abundances after colonization has occurred
#' @export

colfun0<-function(co, xt) {
  sm<-length(xt)
  xt <- xt+rbinom(sm, 1, ilogit(co[1]))*abs(rnorm(sm,0,exp(co[2])))

  return(xt)
}


#' particle filter
#'
#' General function for caluclating the log-likeihood of a stochastic discrete-time model,
#' based on a noisy observation of time-series y. Returns estimates of true values of y, as well as for process noise, observation error, colonization rates, and extinction rates.
#' Function is adapted from the R code of Knape and Valpine (2012), Ecology 93:256-263.
#' @param y A numeric vector of observed values, from which the likelihood of parameters and functions will be determined.
#' @param pars A list of parameter values. Must include elements obs (observation error parameters), proc (process noise parameters), and pcol (colonization parameters), which are passed on the their respecive functions, described below. If edmdat=NULL, then element det (deterministic process parameters) must be included.
#' @param N Number of particles to simulate. Defaults to 1e3.
#' @param detfun A function that simulates deterministic dynamics, which takes in arguments sdet (parameters for deterministic model, taken from pars$proc), and xt, observed abundances at time t. Returns estimated abundances at time t+1 based on deterministic function (either a parametric function or an EDM function). Defaults to detfun0.
#' @param procfun A function that simulates process noise, which takes in arguments sp (parameters for process noise function, taken from pars$proc) and xt (abundances prior to process noise). Returns abundances after process noise has occurred. Defaults to procfun0.
#' @param obsfun An observation function, which takes in up to five variables, including so (a vector of parameter values, inherited from pars$obs), yt (a number, showing observed abundance at time t), xt (predicted abundances), binary value "inverse", and number "N". If inverse = TRUE,
#' then function should simulate N draws from the observation function, centered around value yt. If inverse = FALSE, then function should return log probability denisty of observed value yt given predicted values in xt. Defaults to obsfun0.
#' @param colfun A function simulating colonization events, that takes in two arguments: co, a vector of parameter values taken from pars$pcol, and xt, a number or numeric vector of abundances at time t, before colonization has occurred. Returns predicted abundances after colonization has occurred. Defaults to colful0.
#' @param edmdat A list including arguments to be passed to block_lnlp from rEDM package - see block_lnlp help file for details. Can also include optional matrix "extra_columns", a matrix with length(y) rows including extra covariates for attractor reconstruction, which defaults to NULL (i.e. no additional columns).
#' Default for edmdat is NULL, which implies that EDM will not be applied - instead, a detfun and pars$det must be included.
#' @param dotraceback A logical, indicating whether estimated values and demographic rates should be reported - defaults to FALSE
#' @source Adapted from Knape and Valpine (2012), Ecology 93:256-263.
#' @keywords particle filter, stability, time-series, Taylor power law
#' @return LL, P, rN, x, ind, xsort, dem(col, mor, mucol, mumor)
#' @export

particleFilterLL<-function(y, pars, N=1e3, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=FALSE) {
  LL<-rep(NA, length(y))

  #set up matrices for storage
  if(dotraceback) {
    Nest<-rep(NA, length(y))
    Nsd<-rep(NA, length(y))

    Nest_noproc<-rep(NA, length(y))
    Nsd_noproc<-rep(NA, length(y))
  } else {
    Nest<-NULL
    Nsd<-NULL

    Nest_noproc<-NULL
    Nsd_noproc<-NULL
  }

  #initialize particles
  tstart<-max(c(2, edmdat$E))
  for(i in 1:(tstart)) {
    prd<-obsfun(so = pars$obs, yt = y[i], xt = proc, time = i, inverse = TRUE, N = N)
    if(dotraceback) {
      Nest[i]<-mean(prd)
      Nsd[i]<-sd(prd)
    }
  }

  if(!is.null(edmdat)) {
    smp<-s_map(y, E=edmdat$E, theta = edmdat$theta, silent = TRUE, save_smap_coefficients = TRUE)
    smp_cf<-smp$smap_coefficients[[1]]
  }

  for(i in tstart:length(y)) {
    #prd is deterministic estimate for t=i

    #add process noise
    proc<-procfun(sp = pars$proc, xt = prd, time = i)

    #get likelihood of proc given y[i]
    dobs<-obsfun(so = pars$obs, yt = y[i], xt = proc, time = i)
    mxd<-max(dobs, na.rm=T)
    wdobs<-exp(dobs-mxd)/sum(exp(dobs-mxd))

    #estimates of true state for y[i]
    post_smp<-sample(proc, N, rep=T, prob = wdobs)

    #estimates of deterministic state at t=i+1
    if(is.null(edmdat)) {
      prd<-detfun(sdet = pars$det, xt = post_smp, time = i)
    } else {
      prd<-detfun(smp_cf = smp_cf, yp = y, x = post_smp, time = i)
    }

    #add colonization

    #save likelihood
    LL[i]<-log(mean(exp(dobs-mxd)))+mxd

    #save state
    if(dotraceback) {
      Nest[i]<-mean(post_smp)
      Nsd[i]<-sd(post_smp)

      if(i<length(y)) {
        Nest_noproc[i+1]<-mean(prd)
        Nsd_noproc[i+1]<-sd(prd)
      }
    }
  }
  LLtot <- sum(LL[is.finite(LL)], na.rm=T)

  return(list(LL = LLtot, LLlst=LL, Nest=Nest, Nsd=Nsd, Nest_noproc=Nest_noproc, Nsd_noproc=Nsd_noproc))
}


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





