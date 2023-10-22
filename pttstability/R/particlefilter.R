#' default deterministic function
#'
#' Simulates deterministic component of Ricker model, of the form xt+1 = xt exp(exp(sdet[1])*(1-xt/exp(sdet[2])))
#' @param sdet a numeric vector of length two, specifying growth rate and carrying capacity
#' @param xt a number or numeric vector of abundances at time t
#' @param time the timestep - defaults to NULL (i.e. not used)
#' @param ... additional arguments, for compatability with other usages of the function - values are not used in this implementation
#' @keywords deterministic discrete-time time-series
#' @return a number or numeric vector of length xt, with predicted abundances at time t+1
#' @export

detfun0<-function(sdet, xt, time=NULL, ...) {
  xt = xt*exp(exp(sdet[1])*(1-xt/exp(sdet[2])))
  return(xt)
}

#' deterministic function with time-varying carrying capacity
#'
#' Simulates deterministic component of Ricker model, of the form xt+1 = xt exp(exp(sdet[1])*(1-xt/K))
#' where K varies with time as (sin(time/2)+exp(sdet[2])+0.5)*(2/3). Function is calibrated such that
#' for exp(sdet[2]) = 1, mean(K) = 1.
#' @param sdet a numeric vector of length two, specifying growth rate and carrying capacity
#' @param xt a number or numeric vector of abundances at time t
#' @param time the timestep - defaults to NULL (i.e. not used)
#' @param ... additional arguments, for compatability with other usages of the function - values are not used in this implementation
#' @keywords deterministic discrete-time time-series
#' @return a number or numeric vector of length xt, with predicted abundances at time t+1
#' @export
#'
detfun0_sin<-function(sdet, xt, time=NULL, ...) {
  K<-(sin(time/2)+exp(sdet[2])+0.5)*(2/3)
  xt = xt*exp(exp(sdet[1])*(1-xt/K))
  return(xt)
}

#' calculate estimated total variance
#'
#' Function for estimating stochastic variation in linar process x as a function of relative growth rate and disturbance regime standard deviation.
#' @param sd_proc standard deviation of the (Gaussian) disturbance process
#' @param rgr relative growth rate of the linear process
#' @param waiting_time average waiting time between (random exponentially distributed through time) disturbance events
#' @keywords stochastic linear system
#' @return standard deviation of stochastic variability in x
#' @export
#'
#'
sdproc_abstract<-function(sd_proc, rgr, waiting_time = 1) {
  sqrt((sd_proc^2/waiting_time)/(2*rgr))
}

#' REDM deterministic function
#'
#' Estimates future states of xt based on based behaviour
#' @param smp_cf a matrix of s-map coefficients. Columns correspond to intercept and time lags, rows to observations. Final column corresponds to intercept term.
#' @param yp a matrix of covariates to be multiplied by the smp_cf (typically time lags). Should have one fewer column than smp_cf.
#' @param x observation at time-1, to be used to make the prediction.
#' @param minest minimum value to return for prediction - defaults to 0.
#' @param maxest maximum value to return for prediction - defaults to NULL (no maximum)
#' @param time the time step (i.e. position in smp_cf) for the desired prediction. Prediction will be made based on observation in preceding time point (i.e. time-1).
#' @keywords EDM
#' @return a number or numeric vector of length xt, with predicted abundances at time t+1
#' @source Adapted from Ye, Sugihara, et al. (2015), PNAS 112:E1569-E1576.
#' @export

EDMfun0<-function(smp_cf, yp, x, minest=0, maxest=NULL, time) {
  time_use<-time-1

  nD<-ncol(smp_cf)
  if(nD>2) {
    out<-sum(smp_cf[time_use,-1]*c(yp[rev((time_use-(nD-2)):(time_use-1))], 1))+smp_cf[time_use,1]*x
  } else {
    out<-smp_cf[time_use,-1]+smp_cf[time_use,1]*x
  }
  out[out<minest]<-minest
  if(!is.null(maxest)) {
    out[out>maxest]<-maxest
  }

  out<-out*(x>0)
  out
}

#' Process s-mapping coefficients
#'
#' Processes s-mapping coefficients from rEDM into a matrix of form C1, C2, C3, ... C0, where C0 is the intercept,
#' C1 is the current time step t, C2 is timestep t-1, C3 is timestep t-2, and so on.
#' Rows correspond to the time step used to produce the prediction, e.g. row 4 is used to calculate
#' predicted value for time step 5. This is the format expected by the EDMfun0 function.
#' Note - the format produced by the rEDM package has changed substantially over time,
#' and so if you find that predictions are no longer working in this package, it is likely that this is
#' due to another change in reporting format, in which case you may need to update this function
#' accordingly (e.g. to re-align columns or rows).
#' @param smap_coefs a matrix of s-map coefficients, taken from the SMap function.
#' @return a matrix of s-mapping coefficients
#' @export

process_scof <- function(smap_coefs) {
  ps <- which(colnames(smap_coefs)=="Index")
  if(length(ps) > 0) {
    index = smap_coefs[,ps]
    smap_coefs = smap_coefs[,-ps]

    ps2 <- which(colnames(smap_coefs)=="C0")
    if(length(ps2) > 0) {
      smap_coefs = smap_coefs[,c(c(1:ncol(smap_coefs))[-ps2], ps2)]
    }
    colnames(smap_coefs) <- paste("C", c(1:(ncol(smap_coefs)-1), 0), sep="")

    if(min(index) > 2) {
      for(i in 1:(min(index)-2)) {
        smap_coefs = rbind(rep(NA, ncol(smap_coefs)), smap_coefs)
      }
    }
  }

  smap_coefs
}

#' default process noise function
#'
#' Simulates effects of process noise following a Gaussian perturbation.
#' Note that process noise only influences positive abundances (i.e. process noise cannot contribute to colonization)
#' @param sp a numeric vector of length one or two, specifying either the log-transformed standard deviation of the process noise function,
#' or an intercept and slope for calculating variance of process noise based on a power function of x, of the form var=exp(B0)*x^exp(B1)
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




#' continuous-time process noise function
#'
#' Simulates effects of process noise following a Gaussian perturbation.
#' Note that process noise only influences positive abundances (i.e. process noise cannot contribute to colonization)
#' @param sp a numeric vector of length two or three, where terms 1-2 specify either the log-transformed standard deviation of the process noise function,
#' or an intercept and slope for calculating variance of process noise based on a power function of x, of the form var=exp(B0)*x^exp(B1)
#' The final term in the vector represents the recovery rate - i.e. the continuous time rate at which abundances recover from perturbation
#' @param xt a number or numeric vector of abundances at time t, before process noise has occurred
#' @param waiting_time average time between disturbance events: defaults to 1
#' @param time the timestep - defaults to NULL (i.e. not used)
#' @keywords process noise
#' @return a number or numeric vector of length xt, with predicted abundances after process noise has occurred
#' @import stats
#' @export

procfun_ct<-function(sp, xt, waiting_time = 1, time=NULL) {
  if(length(sp)==2) {
    std_tmp<-exp(sp[1])
    ps2 = 1
  } else {
    std_tmp<-sqrt(exp(sp[1])*xt^exp(sp[2]))
  }
  rgr = exp(sp[length(sp)])

  # number of disturbances in this timestep
  ne = rpois(length(xt),waiting_time)

  maxne = max(ne)
  if(maxne > 0) {
    # disturbance times
    dist_times = matrix(nrow = length(xt), ncol = maxne, data = qexp(pexp(runif(length(xt)*maxne),1/waiting_time),1/waiting_time))

    # apply disturbances
    for(n in 1:maxne) {
      ps = which(ne == n & xt>0) # cycle trough disturbance events: disturbances cannot reverse extinction
      if(length(ps)>0) {
        #nt_dist = nt_det + dist*exp(-rt)
        if(length(std_tmp) == length(xt)) {
          ps2 = ps
        }
        xt[ps] = pmax(0, xt[ps] + rowSums(rnorm(length(ps)*n, 0, std_tmp[ps2])*exp(-rgr*(1-dist_times[ps,1:n,drop = FALSE]))))
      }
    }
  }

  # return output
  return(xt)
}


#' default observation noise function
#'
#' Two options: If inverse=FALSE, calculates the log probability density of observation yt based on true state xt and observation error.
#' Otherwise, simulates N random observations of yt.
#' Observation error follows a Gaussian distribution truncated at zero, using a Tobit distribution.
#' Note that probability density is calculated based on a Tobit distribution, with lower boundary zero.
#' @param so a numeric vector of length one, specifying either log-transformed standard deviation of the observation error as a fraction of the observation,
#' or two log-transformed parameters of the form sd=exp(B0)+exp(B1)*x.
#' @param yt a number, representing a potential observed value of xt
#' @param xt a number or numeric vector of "true" (or simulated) abundances at time t, from which the likelihood of yt will be calculated - defaults to NULL for inverse=TRUE
#' @param inverse a logical specifying whether inverse (i.e. random number generator) function should be implemented - defaults to FALSE
#' @param N number of draws from the random number generator, if inverse=TRUE - defaults to NULL
#' @param minsd minimum observation error allowed (e.g. if observation = 0), to prevent log likelihoods of -infinity - defaults to 0.01
#' @param time the timestep - defaults to NULL (i.e. not used)
#' @keywords observation error
#' @return If inverse=FALSE, returns a list including LL, a number or numeric vector of length xt, with predicted log likelihoods of observation yt,
#' and wts, a number or vector with weights corresponding to the relative likelihood of each observation (after accounting for variable continuous vs. discrete probability distributions).
#' If inverse = FALSE, returns N random draws from the observation function.
#' @import stats
#' @export

obsfun0<-function(so, yt, xt=NULL, inverse=FALSE, N=NULL, minsd=0.01, time=NULL) {
  if(inverse) {
    if(length(so)==1) {
      std_tmp<-exp(so[1])*yt
    } else {
      std_tmp<-exp(so[1])+exp(so[2])*yt
    }
    std_tmp[std_tmp<minsd]<-minsd
    pmax(0, rnorm(n = N, mean = yt, sd = std_tmp))
  } else {
    if(length(so)==1) {
      std_tmp<-exp(so[1])*xt
    } else {
      std_tmp<-exp(so[1])+exp(so[2])*xt
    }
    std_tmp[std_tmp<minsd]<-minsd

    #Tobit distribution (modified from script by J.Clark):
    if(yt == 0) {
      tiny <- 10e-6
      lo = -Inf; hi = 0 # lower and upper limits of normal quantiles for censored region in Tobit (i.e. values <= 0)

      q1 <- pnorm(lo,xt,std_tmp)
      q2 <- pnorm(hi,xt,std_tmp)

      # simulate potential values of y along a normal distribution given modeled mean and sd
      z <- runif(n = length(xt), min = q1, max = q2)
      z <- qnorm(z,xt,std_tmp)

      z[is.infinite(z) & sign(z)>0]  <- lo + tiny
      z[is.infinite(z) & sign(z)<0] <- hi - tiny

      LL = dnorm(z, xt, std_tmp, log = TRUE)
    } else {
      LL = dnorm(yt, xt, std_tmp, log = TRUE)
    }

    return(LL)
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
#' @param edmdat A list including arguments to be passed to SMap from rEDM package - see SMap help file for details. Alternatively, the user can provide a matrix of pre-computed S-map coefficients, in element "smp_cf".
#' Default for edmdat is NULL, which implies that EDM will not be applied - instead, a detfun and pars$det must be included.
#' @param dotraceback A logical, indicating whether estimated values and demographic rates should be reported - defaults to FALSE
#' @param fulltraceback A logical, indicating whether full matrix of particles for all time steps should be returned.
#' @source Adapted from Knape and Valpine (2012), Ecology 93:256-263.
#' @keywords particle filter stability time-series Taylor power law
#' @return LL (total log likelihood), LLlst (log likelihood for each time step), Nest (mean estimated state), Nsd (standard deviation of estimated state), Nest_noproc (mean estimated state at time t+1 without process error), Nsd_noproc (standard deviation of estimated state at time t+1 without process error),
#' fulltracemat (full traceback of particle paths), fulltracemat_noproc (full traceback of particle paths at time t+1 without process noise), and fulltraceindex (index positions for the particle traces over time)
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
      if (!requireNamespace("rEDM", quietly = TRUE)) {
        stop("Package \"rEDM\" needed for this function to work. Please either install it, or provide an \"smp_cf\" matrix in the \"edmdat\" list.",
             call. = FALSE)
      }
      ydf = data.frame(Index = 1:length(y), y = y)
      smp<-rEDM::SMap(dataFrame = ydf, E=edmdat$E, theta = edmdat$theta, lib = c(1, nrow(ydf)), pred = c(1, nrow(ydf)), columns = "y")
      smp_cf<-process_scof(smap_coefs = smp$coefficients)
    } else {
      smp_cf<-edmdat$smp_cf
    }

    if(!is.null(edmdat$ytot)) {
      mx<-max(edmdat$ytot, na.rm=T)
    } else {
      mx<-NULL
    }

    if(!is.null(edmdat$maxest)) {
      mx<-edmdat$maxest
    }
  }

  if(fulltraceback) {
    fulltracemat<-matrix(nrow=length(y), ncol=N)
    fulltraceindex<-matrix(nrow=length(y), ncol=N)
    fulltracemat_noproc<-matrix(nrow=length(y), ncol=N)

    #state before process noise
    if(tstart>1) {
      fulltracemat_noproc[tstart-1,]<-prd
    }
  } else {
    fulltracemat<-NULL
    fulltraceindex<-NULL
    fulltracemat_noproc<-NULL
  }

  for(i in tstart:length(y)) {
    #prd is deterministic estimate for t=i

    #add process noise
    proc<-procfun(sp = pars$proc, xt = prd, time = i)

    #get likelihood of proc given y[i]
    dobs<-obsfun(so = pars$obs, yt = y[i], xt = proc, time = i)
    #tmpobsout<-obsfun(so = pars$obs, yt = y[i], xt = proc, time = i)
    #dobs<-tmpobsout$LL
    mxd<-max(dobs, na.rm=T)
    wdobs<-exp(dobs-mxd)/sum(exp(dobs-mxd))
    #wdobs<-tmpobsout$wts

    #estimates of true state for y[i]
    if(fulltraceback) {
      smp_tmp<-sample(length(proc), N, replace=T, prob = wdobs)
      post_smp<-proc[smp_tmp]
    } else {
      post_smp<-sample(proc, N, replace=T, prob = wdobs)
    }

    #now, prd is estimate of deterministic state at t=i+1
    if(is.null(edmdat)) {
      prd<-detfun(sdet = pars$det, xt = post_smp, time = i+1)
    } else {
      if(is.null(mx)) {
        prd<-detfun(smp_cf = smp_cf, yp = y, x = post_smp, time = i+1)
      } else {
        prd<-detfun(smp_cf = smp_cf, yp = y, x = post_smp, time = i+1, maxest = mx)
      }
    }

    #add colonization
    prd<-colfun(co=pars$pcol, xt=prd)

    #save state
    if(dotraceback) {
      #state at time t
      Nest[i]<-mean(post_smp)
      Nsd[i]<-sd(post_smp)

      #state at time t+1, before process noise
      Nest_noproc[i]<-mean(prd)
      Nsd_noproc[i]<-sd(prd)
    }

    if(fulltraceback) {
      #state at time t
      fulltracemat[i,]<-post_smp
      fulltraceindex[i,]<-smp_tmp

      #state at time t+1, before process noise
      fulltracemat_noproc[i,]<-prd
    }

    #save likelihood
    LL[i]<-log(mean(exp(dobs-mxd)))+mxd
  }
  LLtot <- sum(LL[is.finite(LL)], na.rm=T)

  return(list(LL = LLtot, LLlst=LL, Nest=Nest, Nsd=Nsd, Nest_noproc=Nest_noproc, Nsd_noproc=Nsd_noproc, fulltracemat=fulltracemat, fulltracemat_noproc=fulltracemat_noproc, fulltraceindex=fulltraceindex))
}


#' Sort output of particle filter
#'
#' Sorts outputs of particle filter based on index - returns a sorted list of particles, based on the
#' sampling trajectory through time. This is a somewhat more accurate estiamte of the true posterior than
#' are the stepwise samples provided by the filter.
#' @param fulltracemat full output of particles from the particleFilterLL function
#' @param fulltraceindex full output of particle indices from the particleFilterLL function
#' @param nsmp number of particle paths to sample - defaults to NULL, which samples all paths
#' @keywords particle filter
#' @return an index-sorted matrix - each column shows the trajectory of a single particle
#' @export

indexsort<-function(fulltracemat, fulltraceindex, nsmp=NULL) {
  sortedmat<-matrix(nrow=nrow(fulltracemat), ncol=ncol(fulltracemat))
  startps<-min(which(apply(fulltraceindex, 1, function(x) mean(is.na(x)))!=1))


  if(is.null(nsmp)) {
    smpps<-1:ncol(sortedmat)
  } else {
    smpps<-sample(ncol(sortedmat), nsmp, replace=FALSE)
  }
  for(j in smpps) {
    indx<-j
    for(i in nrow(sortedmat):startps) {
      sortedmat[i,j]<-fulltracemat[i,indx]
      indx<-fulltraceindex[i,indx]
    }
  }
  sortedmat<-sortedmat[,colMeans(is.na(sortedmat))!=1]
  sortedmat
}

#' calculate likelihood for piecewise data
#'
#' Calculates likelihoods across several segments of data - e.g. multiple plots from a single experiment.
#' See documentation for particleFilterLL_piecewise for examples of use.
#' @param param parameters to be passed to likelihood0 function
#' @param y the time series to be analyzed
#' @param libuse_y a matrix with two columns, specifying the start end end positions of segments within vector y
#' @param smap_coefs a matrix of s-mapping coefficients
#' @param Euse embedding dimension for the s-mapping analysis
#' @param tuse theta for s-mapping analysis
#' @param N number of particles
#' @param colpar parameters to be passed to the colfun0 - defaults to c(logit(1e-6), log(0.1))
#' @keywords dewdrop regression particle filter
#' @return summed log likelihood across all segments
#' @export

likelihood_EDM_piecewise<-function(param, y, libuse_y, smap_coefs, Euse, tuse, N, colpar=c(logit(1e-6), log(0.1))) {
  LLtot<-0

  for(i in 1:nrow(libuse_y)) {
    ysegment<-y[libuse_y[i,1]:libuse_y[i,2]]
    smap_coefs_segment<-smap_coefs[libuse_y[i,1]:libuse_y[i,2],]

    LLtot<-LLtot+likelihood0(param = param, y=ysegment, parseparam = function(x) parseparam0(x, colparam=colpar),
                             detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse, smp_cf=smap_coefs_segment, ytot=y), N = N)
  }

  return(LLtot)
}


#' run particle filter across piecewise data
#'
#' Calculates likelihoods across several segments of data - e.g. multiple plots from a single experiment.
#' Requires several implicitely defined variables to run:
#' @param param parameters to be passed to parseparam0 function
#' @param N number of particles
#' @param y the time series to be analyzed
#' @param libuse_y a matrix with two columns, specifying the start end end positions of segments within vector y
#' @param smap_coefs a matrix of s-mapping coefficients
#' @param Euse embedding dimension for the s-mapping analysis
#' @param tuse theta for s-mapping analysis
#' @param colpar parameters to be passed to the colfun0 - defaults to c(logit(1e-6), log(0.1))
#' @param nsmp number of sample particle trajectories to return - defaults to 1
#' @param lowerbound minimum accepted likelihood - used to automatically select number of particles. Defaults to -999
#' @param maxNuse maximum number of particles to simulate - defaults to 512000
#' @keywords dewdrop regression particle filter
#' @return results from particle filter - including mean estimates (Nest) and standard deviations (Nsd), across particles,
#' and sample particle trajectories with (Nsmp) and without (Nsmp_noproc) process noise
#' @export
#' @examples
#' # load data
#' data(dat)
#' # sort by index
#' dat = dat[order(dat$treatment, dat$number, dat$time),]
#' 
#' \donttest{
#' # make list of starting and ending positions for each replicate in the dat list
#' libmat<-NULL
#' trtmat<-data.frame(trt=as.character(sort(unique(dat$treatment))))
#' datnum<-1:nrow(dat)
#'
#' for(i in 1:nrow(trtmat)) {
#'   ps1<-which(dat$treatment==trtmat$trt[i])
#'   replst<-sort(unique(dat$number[ps1]))
#'
#'   for(j in 1:length(replst)) {
#'     ps2<-which(dat$number[ps1]==replst[j])
#'     libmat<-rbind(libmat, data.frame(trt=trtmat$trt[i], rep=replst[j],
#'       start=min(datnum[ps1][ps2]), end=max(datnum[ps1][ps2])))
#'   }
#' }
#'
#' ## run particle filter
#' # select treatment to analyse: enter either "LSA" or "LSP"
#' trtuse<-"HSP"
#' # extract library positions for treatment
#' libuse<-as.matrix(libmat[libmat$trt==trtuse,c("start", "end")])
#' # save abundance data to variable y
#' yps<-which(dat$treatment==trtuse)
#' y<-dat[,"Chlamydomonas.terricola"][yps]
#' libuse_y<-libuse-min(libuse)+1 # translate positions in dat to positions in y vector
#' y<-y/sd(y) # standardize to mean of one
#' timesteps<-dat$time[yps]
#'
#' # get EDM parameters
#' require(rEDM) # load rEDM package
#' sout<-NULL
#' ydf = data.frame(Index = 1:length(y), y=y)
#' for(E in 2:4) {
#'   # note: ignore disjoint predictions warning
#'   sout<-rbind(sout, data.frame(E = E, PredictNonlinear(dataFrame = ydf, E = E,
#'                                                        columns = "y",
#'                                                        lib = c(t(libuse_y)), pred = c(t(libuse_y)),
#'                                                        showPlot = FALSE)))
#' }
#' tuse<-sout$Theta[which.max(sout$rho)] # find theta (nonlinerity) parameter
#' euse<-sout$E[which.max(sout$rho)] # find embedding dimension
#' spred<-SMap(dataFrame = ydf, E = euse,
#'             theta = tuse,
#'             columns = "y", lib = c(t(libuse_y)), pred = c(t(libuse_y)))
#'
#' # set priors (log-transformed Beta_obs, Beta_proc1, and Beta_proc2)
#' minvUSE_edm<-c(log(0.001), log(0.001))  # lower limits
#' maxvUSE_edm<-c(log(2), log(2)) # upper limits
#' }
#' 
#' \dontrun{
#' ## Run filter
#' # density, sampler, and prior functions for EDM function
#' # Commented-out code: Install BayesianTools package from GitHub if needed
#' #require(devtools)
#' #install_github("florianhartig/BayesianTools/BayesianTools")
#' # see BayesianTools documentation for details
#' require(BayesianTools)
#' density_fun_USE_edm<-function(param) density_fun0(param = param,
#'   minv = minvUSE_edm, maxv=maxvUSE_edm)
#' sampler_fun_USE_edm<-function(x) sampler_fun0(n = 1, minv = minvUSE_edm, maxv=maxvUSE_edm)
#' prior_edm <- createPrior(density = density_fun_USE_edm, sampler = sampler_fun_USE_edm,
#'                          lower = minvUSE_edm, upper = maxvUSE_edm)
#'   niter<-5e3 # number of steps for the MCMC sampler
#'   N<-1e3 # number of particles
#'   smap_coefs<-process_scof(spred$coefficients) # coefficients from s-mapping routine
#'
#'   # likelihood and bayesian set-ups for EDM functions
#'   likelihood_EDM_piecewise_use<-function(x) {
#'     # default values for filter - see likelihood_EDM_piecewise documentation for details
#'     # note that colpar are set near zero because we expect no colonisation into a closed microcosm.
#'     likelihood_EDM_piecewise(param=x, y, libuse_y, smap_coefs, euse, tuse, N,
#'                              colpar = c(logit(1e-06), log(0.1)))
#'   }
#'
#'   bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM_piecewise_use,
#'      prior = prior_edm)
#'
#'   # run MCMC optimization (will take ~ 15-20 min)
#'   out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM,
#'      settings = list(iterations=niter, consoleUpdates=20))
#'   burnin<-floor(niter/5) # burnin period
#'   plot(out_EDM, start=burnin) # plot MCMC chains
#'   gelmanDiagnostics(out_EDM, start=burnin) # calculate Gelman statistic
#'   summary(out_EDM, start=burnin) # coefficient summary
#'
#'   ## extract abundance estimate from particle filter
#'   # use final estimate from MCMC chain
#'   smp_EDM<-(getSample(out_EDM, start=floor(niter/5)))
#'   tmp<-particleFilterLL_piecewise(param = smp_EDM[nrow(smp_EDM),], N=N, y = y, libuse_y = libuse_y,
#'                                   smap_coefs = smap_coefs, Euse = euse, tuse = tuse)
#'   # mean estimated abundance
#'   simout<-tmp$Nest
#'   # sd estimated abundance
#'   sdout<-tmp$Nsd
#'   # sample from true particle trajectory
#'   simout_smp<-tmp$Nsmp
#'   # sample from true particle trajectory pre-process noise
#'   simout_smp_noproc<-tmp$Nsmp_noproc
#'
#'   plot(timesteps, simout, xlab="Time", ylab="Abundance")
#'   abline(h=0, lty=3)
#' }

particleFilterLL_piecewise<-function(param, N, y, libuse_y, smap_coefs, Euse, tuse, colpar=c(logit(1e-6), log(0.1)), nsmp=1, lowerbound = -999, maxNuse = 512000) {
  # piecewise particle filter function with automatic N selection for each time series chunk
  pars<-parseparam0(param, colparam=colpar)
  tuse_edm<-tuse

  pfout<-data.frame(rep=NA, Nest=rep(NA, length(y)), Nsd=NA, Nsmp=NA, Nsmp_noproc=NA)
  for(i in 1:nrow(libuse_y)) {
    ysegment<-y[libuse_y[i,1]:libuse_y[i,2]]
    smap_coefs_segment<-smap_coefs[libuse_y[i,1]:libuse_y[i,2],]

    Nuse = N
    LLtmp = -Inf
    while(LLtmp <=lowerbound & Nuse <= maxNuse) {
      tmp<-try(particleFilterLL(ysegment, pars, N=N, detfun = EDMfun0,
                            edmdat = list(E=Euse, theta=tuse_edm, smp_cf=smap_coefs_segment),
                            dotraceback = TRUE, fulltraceback = TRUE))
      if(is.character(tmp)) {
        LLtmp = -Inf
      } else {
        LLtmp = tmp$LL
      }
      Nuse = 2*Nuse
    }


    pfout$Nest[libuse_y[i,1]:libuse_y[i,2]]<-tmp$Nest
    pfout$Nsd[libuse_y[i,1]:libuse_y[i,2]]<-tmp$Nsd
    pfout$rep[libuse_y[i,1]:libuse_y[i,2]]<-i
    pfout$Nsmp[libuse_y[i,1]:libuse_y[i,2]]<-c(indexsort(tmp$fulltracemat, tmp$fulltraceindex, nsmp=1))
    pfout$Nsmp_noproc[libuse_y[i,1]:libuse_y[i,2]]<-c(indexsort(tmp$fulltracemat_noproc, tmp$fulltraceindex, nsmp=1))
  }

  return(pfout)
}
