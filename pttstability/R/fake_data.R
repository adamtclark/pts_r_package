#' Simulate fake data
#'
#' Simulates a time series with known process noise and observation error functions. Implemented May's chaotic logistic function, X(t+1) = X(t) exp(r (1-X(t)/K)).
#' @param n number of timesteps to simulate
#' @param n0 starting population size - defaults to NULL, which starts at K/10
#' @param obs a numeric vector of length two (b0, b1) specifying observation error, parameterizing the function sd(obs) = b0 + b1*X, where X is the system state
#' @param proc a numeric vector of length two (lc, z) specifying process noise, parameterizing the function var(proc) = exp(lc + z*log(X)), where X is the system state
#' @param pcol a numeric vector of length two, including probability of colonization per timestep given that abundance = 0, and abundance after colonization
#' @param r intrinsic growth rate for May's logistic function
#' @param K carrying capacity for May's logistic function
#' @param doplot a logical specifying wether output should be plotted - defaults to FALSE
#' @keywords Taylor power law, stability, time-series, fake data
#' @return An n-by-3 dataframe of states, including obs (observed values), truth (true values), and noproc (values without process noise)
#' @import graphics
#' @import stats
#' @export
#' @examples
#'  #run function
#'  datout<-makedynamics(n=2e4, proc = c(-2,1.2))
#'
#'  #show regression of variance vs. mean for binned data
#'  datout_ps<-datout[datout$true>0 & datout$noproc>0,]
#'  #bins
#'  sq<-seq(0, quantile(datout$true, 0.95), by=0.01)
#'  ctd<-cut(datout_ps$noproc, sq)
#'  #calculate mean and variance by bin
#'  tdat<-data.frame(mu=(sq[-1]+sq[-length(sq)])/2,
#'       var=tapply((datout_ps$true-datout_ps$noproc)^2, ctd, mean),
#'      unk=tapply((datout_ps$true-datout_ps$noproc)^2, ctd, sd))
#'  #plot result
#'  plot(log(tdat$mu), log(tdat$var), xlab="mu", ylab="var")
#'  #show regression
#'  summary(mod<-lm(log(var)~log(mu), tdat, weights = 1/tdat$unk)); abline(mod, col=2)

makedynamics<-function(n=1000, n0=NULL,
                       obs=c(log(0.01), log(0.1)), proc=c(-2, log(1.2)),
                       pcol=c(logit(0.2), log(1e-2)), r=log(3), K=log(1),
                       doplot=FALSE) {

  #transform deterministic variables
  r<-exp(r); K<-exp(K)

  #Set up matrices
  Sdat<-Sdat_noproc<-rep(NA,n)

  #Starting value
  if(is.null(n0)) {
    n0<-K*0.1
  }
  Sdat[1]<-Sdat_noproc[1]<-n0

  #Colonization events
  col_event<-rbinom(n, 1, ilogit(pcol[1]))*exp(pcol[2])

  for(j in 2:(n)) {
    n_tmp<-Sdat[j-1]

    if(n_tmp==0) {
      #Colonization
      n_tmp<-n_tmp+col_event[j]
    } else {
      #Deterministic dynamic
      n_tmp<-n_tmp*exp(r*(1-n_tmp/K))
    }

    #Process noise (following Taylor's law)
    if(n_tmp>0) {
      pvar<-exp(proc[1]+exp(proc[2])*log(n_tmp))
      n_tmp_proc<-max(c(0, n_tmp+rnorm(1,0,sqrt(pvar))))
    }

    Sdat_noproc[j]<-n_tmp
    Sdat[j]<-n_tmp_proc
  }

  #Observation error (linear function of true abundance)
  sd_obs<-exp(obs[1])+exp(obs[2])*Sdat
  Sdat_obs<-pmax(0, Sdat+rnorm(n, 0, sd_obs))
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


#' Get rates
#'
#' Calculates colonization rate, mortality rate, and expected mean occupancy time based on a time series
#' @param dat a numeric vector, including the timeseries
#' @keywords stability, time-series
#' @return a list including colonization and mortality probability per time step (pc and pm, respectively), and pocc, the expected fraction of time that the species will be present

getcm<-function(dat) {
  #time to colonization
  pc<-sum(dat[-1]>0 & dat[-length(dat)]==0)/sum(dat[-length(dat)]==0)

  #time to mortality
  pm<-sum(dat[-1]==0 & dat[-length(dat)]>0)/sum(dat[-length(dat)]>0)

  #mean occupancy
  pocc<-(1/pc*0+1/pm*1)/(1/pc+1/pm)

  #return output
  return(list(pc=pc, pm=pm, pocc=pocc))
}
