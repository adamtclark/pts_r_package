#' Logit
#'
#' Returns the logit transformation of x
#' @param x a number, vector, matrix, etc. to be transformed from (0, 1) to (-inf inf) by the logit transform
#' @param ... additional arguments to be passed to plogis
#' @keywords logit
#' @return transformed result - impossible values are replaced with NA, without warnings
#' @import stats
#' @export

logit<-function(x, ...) qlogis(x, ...)


#' Inverse logit
#'
#' Returns the inverse logit transformation of x
#' @param x a number, vector, matrix, etc. to be transformed from (-inf, inf) to (0 1) by the inverse logit transform
#' @param ... additional arguments to be passed to plogis
#' @keywords logit
#' @return transformed result
#' @import stats
#' @export

ilogit<-function(x, ...) plogis(x, ...)

#' Get inverse logit-normal mode
#'
#' Returns a mean for a logit normal such that the mode will be centered around mu
#' @param mu the value around which the mode should be centered (in logit space)
#' @param sd the standard deviation of the logit distribution (in logit space)
#' @keywords MCMC optimization
#' @return the proposed mean for the distribution
#' @export

logitnormal_imode = function(mu, sd){
  mode<-mu - (sd)^2*(2*ilogit(mu)-1)
  return(mode)
}


#' Get inverse log-normal mode
#'
#' Returns a mean for a lognormal such that the mode will be centered around mu
#' @param mu the value around which the mode should be centered (in log space)
#' @param sd the standard deviation of the lognormal distribution (in log space)
#' @keywords MCMC optimization
#' @return the proposed mean for the distribution
#' @export

lognormal_imode = function(mu, sd){
  mode<-mu+sd^2
  return(mode)
}


