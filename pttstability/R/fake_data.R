#' Get rates
#'
#' Calculates colonization rate, mortality rate, and expected mean occupancy time based on a time series
#' @param dat a numeric vector, including the timeseries
#' @keywords stability time-series
#' @return a list including colonization and mortality probability per time step (pc and pm, respectively), and pocc, the expected fraction of time that the species will be present
#' @export

getcm<-function(dat) {
  #time to colonization
  pc<-sum((!is.na(dat[-1]) & dat[-1]>0) & (!is.na(dat[-length(dat)]) & dat[-length(dat)]==0), na.rm=T)/sum((!is.na(dat[-length(dat)]) & dat[-length(dat)]==0), na.rm=T)

  #time to mortality
  pm<-sum((!is.na(dat[-1]) & dat[-1]==0) & (!is.na(dat[-length(dat)]) & dat[-length(dat)]>0), na.rm=T)/sum((!is.na(dat[-length(dat)]) & dat[-length(dat)]>0), na.rm=T)

  #mean occupancy
  pocc<-(1/pc*0+1/pm*1)/(1/pc+1/pm)

  #return output
  return(list(pc=pc, pm=pm, pocc=pocc))
}


#' Simulate general time series
#'
#' Simulates a time series following a user-defined deterministic function, observation function, process noise function, and colonization function.
#' @param n number of timesteps to simulate
#' @param n0 starting population size
#' @param pdet a numeric vector of parameters for the deterministic function
#' @param proc a numeric vector of parameters for the process noise function
#' @param obs a numeric vector of parameters for the observation error function
#' @param pcol a numeric vector of parameters for the colonization function
#' @param detfun A function that simulates deterministic dynamics, which takes in arguments sdet (parameters for deterministic model, taken from pars$proc), and xt, observed abundances at time t. Returns estimated abundances at time t+1 based on deterministic function (either a parametric function or an EDM function). Defaults to detfun0.
#' @param procfun A function that simulates process noise, which takes in arguments sp (parameters for process noise function, taken from pars$proc) and xt (abundances prior to process noise). Returns abundances after process noise has occurred. Defaults to procfun0.
#' @param obsfun An observation function, which takes in up to five variables, including so (a vector of parameter values, inherited from pars$obs), yt (a number, showing observed abundance at time t), xt (predicted abundances), binary value "inverse", and number "N". If inverse = TRUE,
#' then function should simulate N draws from the observation function, centered around value yt. If inverse = FALSE, then function should return log probability denisty of observed value yt given predicted values in xt. Defaults to obsfun0.
#' @param colfun A function simulating colonization events, that takes in two arguments: co, a vector of parameter values taken from pars$pcol, and xt, a number or numeric vector of abundances at time t, before colonization has occurred. Returns predicted abundances after colonization has occurred. Defaults to colful0.
#' @param doplot a logical specifying wether output should be plotted - defaults to FALSE
#' @keywords Taylor power law stability time-series
#' @return An n-by-3 dataframe of states, including obs (observed values), truth (true values), and noproc (values without process noise)
#' @import graphics
#' @import stats
#' @export
#' @examples
#'  #run function
#'  datout<-makedynamics_general(n=2e4, proc = c(-2,log(1.2)))
#'
#'  #show regression of variance vs. mean for binned data
#'  datout_ps<-datout[datout$true>0 & datout$noproc>0,]
#'  #bins
#'  sq<-seq(0, quantile(datout$true, 0.95), length=50)
#'  ctd<-cut(datout_ps$noproc, sq)
#'  #calculate mean and variance by bin
#'  tdat<-data.frame(mu=(sq[-1]+sq[-length(sq)])/2,
#'       var=tapply((datout_ps$true-datout_ps$noproc)^2, ctd, mean))
#'  #plot result
#'  plot(log(tdat$mu), log(tdat$var), xlab="mu", ylab="var")
#'  #show regression
#'  summary(mod<-lm(log(var)~log(mu), tdat)); abline(mod, col=2)

makedynamics_general<-function(n=1000, n0=0.1,
                       pdet=c(log(3), log(1)), proc=c(log(1)),
                       obs=c(log(1)), pcol=c(logit(0.2), log(1)),
                       detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0,
                       doplot=FALSE) {

  #Set up matrices
  Sdat_obs<-Sdat<-Sdat_noproc<-rep(NA,n)

  #Starting value
  Sdat[1]<-Sdat_noproc[1]<-n0

  for(j in 2:(n)) {
    n_tmp<-Sdat[j-1]

    if(n_tmp==0) {
      #Colonization
      n_tmp<-colfun(pcol, n_tmp)
    } else {
      #Deterministic dynamic
      n_tmp<-detfun(pdet, n_tmp, time=j)
    }

    #Process noise (following Taylor's law)
    if(n_tmp>0) {
      n_tmp_proc<-max(c(0, procfun(proc, n_tmp, time=j)))
    }

    Sdat_noproc[j]<-n_tmp
    Sdat[j]<-n_tmp_proc
  }

  #Observation error (linear function of true abundance)
  for(j in 1:(n)) {
    Sdat_obs[j]<-pmax(0, obsfun(obs, Sdat[j], inverse=TRUE, N=1, time=j))
  }
  Sdat_obs[Sdat==0]<-0 #no false positives

  #Output matrix
  datout<-data.frame(obs=Sdat_obs, true=Sdat, noproc=Sdat_noproc)

  if(doplot) {
    matplot(1:n, datout, type="l", xlab="t", ylab="N",
            col=c("blue", "black", "red"), lty=1)
    abline(h=0, lty=3)
    legend("topleft", c("obs", "true", "noproc"), col=c("blue", "black", "red"), lty=1, bty="n")
  }

  return(datout)
}






