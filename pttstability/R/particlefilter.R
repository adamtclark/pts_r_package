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
#' @param smp_cf a matrix of s-map coefficients, taken from the s_map function. Columns correspond to intercept and time lags, rows to observations.
#' @param yp a matrix of covariates to be multiplied by the smp_cf (typically time lags). Should have one fewer column than smp_cf.
#' @param x observation at time-1, to be used to make the prediction.
#' @param minest minimum value to return for prediction - defaults to 0.
#' @param time the time step (i.e. position in smp_cf) for the desired prediction. Prediction will be made based on observation in preceding time point (i.e. time-1).
#' @keywords applies EDM from the rEDM package to reconstruct deterministic time series dynamics
#' @return a number or numeric vector of length xt, with predicted abundances at time t+1
#' @import rEDM
#' @source @source Adapted from Ye, Sugihara, et al. (2015), PNAS 112:E1569-E1576.
#' @export

EDMfun0<-function(smp_cf, yp, x, minest=0, time) {
  time_use<-time-1

  nD<-ncol(smp_cf)
  if(nD>2) {
    out<-sum(smp_cf[time_use,-1]*c(yp[rev((time_use-(nD-2)):(time_use-1))], 1))+smp_cf[time_use,1]*x
  } else {
    out<-smp_cf[time_use,-1]+smp_cf[time_use,1]*x
  }
  out[out<minest]<-minest
  out<-out*(x>0)
  out
}

#' default process noise function
#'
#' Simulates effects of process noise following a Gaussian perturbation.
#' Note that process noise only influences positive abundances (i.e. process noise cannot contribute to colonization)
#' @param sp a numeric vector of length one or two, specifying either the log-transformed standard deviation of the process noise function, or an intercept and slope for calculating process nosise based on a power function of x
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
    std_tmp<-sqrt(exp(sp[1])*xt^exp(sp[2]))
  }

  if(!inverse) {
    sm<-length(xt)
    xt = pmax(0, xt + rnorm(sm, 0, std_tmp))*(xt>0)
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
#' @param so a numeric vector of length one, specifying log-transformed standard deviation of the observation error, as a fraction of the observation
#' @param yt a number, representing a potential observed value of xt
#' @param xt a number or numeric vector of "true" (or simulated) abundances at time t, from which the likelihood of yt will be calculated - defaults to NULL for inverse=TRUE
#' @param inverse a logical specifying whether inverse (i.e. random number generator) function should be implemented - defaults to FALSE
#' @param N number of draws from the random number generator, if inverse=TRUE - defaults to NULL
#' @param minsd minimum observation error allowed (e.g. if observation = 0), to prevent log likelihoods of -infinity - defaults to 0.01
#' @param time the timestep - defaults to NULL (i.e. not used)
#' @keywords observation error
#' @return If inverse=FALSE, a number or numeric vector of length xt, with predicted log likelihoods of observation yt.
#' If inverse = FALSE, returns N random draws from the observation function.
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
  ps<-which(xt==0)

  if(length(ps)>0) {
    sm<-length(xt[ps])
    xt[ps] <- xt[ps]+rbinom(sm, 1, ilogit(co[1]))*abs(rnorm(sm,0,exp(co[2])))
  }

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
#' @param fulltraceback A logical, indicating whether full matrix of particles for all time steps should be returned.
#' @source Adapted from Knape and Valpine (2012), Ecology 93:256-263.
#' @keywords particle filter, stability, time-series, Taylor power law
#' @return LL (total log likelihood), LLlst (log likelihood for each time step), Nest (mean estimated state), Nsd (standard deviation of estimated state), Nest_noproc (mean estimated state without process error), Nsd_noproc (standard deviation of estimated state without process error), fulltracemat (full traceback of particle paths)
#' @export

particleFilterLL<-function(y, pars, N=1e3, detfun=detfun0, procfun=procfun0, obsfun=obsfun0, colfun=colfun0, edmdat=NULL, dotraceback=FALSE, fulltraceback=FALSE) {
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
    prd<-obsfun(so = pars$obs, yt = y[i], time = i, inverse = TRUE, N = N)
    if(dotraceback) {
      Nest[i]<-mean(prd)
      Nsd[i]<-sd(prd)
    }
  }

  if(!is.null(edmdat)) {
    if(is.null(edmdat$smp_cf)) {
      smp<-s_map(y, E=edmdat$E, theta = edmdat$theta, silent = TRUE, save_smap_coefficients = TRUE)
      smp_cf<-smp$smap_coefficients[[1]]
    } else {
      smp_cf<-edmdat$smp_cf
    }
  }

  if(fulltraceback) {
    fulltracemat<-matrix(nrow=length(y), ncol=N)
  } else {
    fulltracemat<-NULL
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
    post_smp<-sample(proc, N, replace=T, prob = wdobs)
    if(fulltraceback) {
      fulltracemat[i,]<-post_smp
    }

    #estimates of deterministic state at t=i+1
    if(is.null(edmdat)) {
      prd<-detfun(sdet = pars$det, xt = post_smp, time = i+1)
    } else {
      prd<-detfun(smp_cf = smp_cf, yp = y, x = post_smp, time = i+1)
    }

    #save state
    if(dotraceback) {
      Nest[i]<-mean(post_smp)
      Nsd[i]<-sd(post_smp)

      if(i<length(y)) {
        Nest_noproc[i+1]<-mean(prd)
        Nsd_noproc[i+1]<-sd(prd)
      }
    }

    #add colonization
    prd<-colfun(co=pars$pcol, xt=prd)

    #save likelihood
    LL[i]<-log(mean(exp(dobs-mxd)))+mxd
  }
  LLtot <- sum(LL[is.finite(LL)], na.rm=T)

  return(list(LL = LLtot, LLlst=LL, Nest=Nest, Nsd=Nsd, Nest_noproc=Nest_noproc, Nsd_noproc=Nsd_noproc, fulltracemat=fulltracemat))
}
