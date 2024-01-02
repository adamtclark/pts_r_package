#' pttstability: Methods for Measuring Stability in Systems Without Static Equilibria
#'
#' The pttstability ("ParTicle-Takens Stability") package is a collection of functions
#' that can be used to estimate the parameters of a stochastic state space model (i.e.
#' a model where a time series is observed with error).
#'
#' @section Applications:
#'
#' The goal of this package is to estimate the variability around a deterministic process, both in terms of
#' observation error - i.e. variability due to imperfect observations that does not influence system state -
#' and in terms of process noise - i.e. stochastic variation in the actual state of the process.
#' Unlike classical methods for estmiating variability, this package does not necesarilly assume that
#' the deterministic state is fixed (i.e. a fixed-point equilibrium), meaning that variability around a
#' dynamic trajectory can be estimated (e.g. stochastic fluctuations during predator-prey dynamics).
#'
#' By combining information about both the estimated deterministic state of the system and the estimated
#' effects of process noise, this package can be used to compute a dynamic analog of various stability metrics
#'  - e.g. coefficient of variation (CV) or invariability. Estimated extinction rates and colonization rates
#'  can also be estimated.
#'
#' @section Contents:
#'
#' This package builds on three existing toolkits. First, it applies an updated version of the
#' "particleFilterLL" particle filter function of Knape and Valpine (2012) to calculate likelihoods of
#' observed time series given modeled dynamics.
#' Second, it applies empirical dynamic modeling (EDM) methods to estimate deterministic
#' dynamics even in cases where the underlying equations governing system behavior are not known.
#' These models are based on Takens delay-embedding theorem, from which this package takes its name.
#' Finally, it (optionally) uses the MCMC fitting methods from the BayesianTools package to estimate paramter values for
#' the observation error, process noise, and (optionally) deterministic functions underlying observed dynamics.
#'
#' The default observation error and process noise functions in this package (obsfun0 and procfun0)
#' take advantage of the Taylor Power law to separate noise components for relatively short time series.
#' Observation error is assumed to scale with the observed state as sd_obs(x) = x*exp(obs),
#' Process noise is either a constant (i.e. sd_proc(x) = exp(proc)), or, if two variables are given,
#' process noise scales as a power function of the observed value as sd_proc(x) = sqrt(exp(proc1)*x^exp(proc2))
#'
#' Note that although we include default functions in this package, users are able (and encouraged!) to write
#' their own (including for observation error, process noise, deterministic dynamics, priors, and likelihoods).
#'
#' @docType package
#' @source Adam T. Clark, Lina K. Mühlbauer, Helmut Hillebrand, and Canan Karakoç. (2022). Measuring Stability in Ecological Systems Without Static Equilibria. Ecosphere 13(12):e4328.
#' @source Knape, J., and Valpine, P. (2012). Fitting complex population models by combining particle filters with Markov chain Monte Carlo. Ecology 93:256-263.
#' @source Ye, H., Sugihara, G., et al. (2015). Equation-free ecosystem forecasting. PNAS 112:E1569-E1576.
#' @source Ye, H., et al. (2019). rEDM: Applications of Empirical Dynamic Modeling from Time Series. R package version 0.7.2.
#' @source Hartig, F., et al. (2019). BayesianTools: General-Purpose MCMC and SMC Samplers and Tools for Bayesian Statistics. R package version 0.1.6.
#' @name pttstability
#' @keywords internal
#' "_PACKAGE"
#' @examples
#' #Set seed
#' set.seed(5826)
#'
#' ## Simulate data
#' pars_true<-list(obs=log(0.15),
#'             proc=c(log(0.1)),
#'             pcol=c(logit(0.2), log(0.1)),
#'             det=c(log(1.2),log(1)))
#' #parameters for the filter
#' pars_filter<-pars_true
#'
#' #generate random parameter values
#' datout<-makedynamics_general(n = 100, n0 = exp(pars_true$det[2]),
#'                              pdet=pars_true$det, proc = pars_true$proc,
#'                              obs = pars_true$obs, pcol = pars_true$pcol,
#'                              detfun = detfun0_sin, procfun = procfun0,
#'                              obsfun=obsfun0, colfun=colfun0)
#' y<-datout$obs
#' plot(y, type = "l", xlab="time", ylab="observed abundance")
#'
#' #get parameters for S-mapping
#' sout<-data.frame(E = 2:4, theta = NA, RMSE = NA)
#' for(i in 1:nrow(sout)) {
#'   optout = optimize(f = function(x) {S_map_Sugihara1994(Y = y, E = sout$E[i],
#'     theta = x)$RMSE}, interval = c(0,10))
#'   sout$theta[i] = optout$minimum
#'   sout$RMSE[i] = optout$objective
#' }
#' tuse<-sout$theta[which.min(sout$RMSE)] # find theta (nonlinerity) parameter
#' Euse<-sout$E[which.min(sout$RMSE)] # find embedding dimension
#' spred<-S_map_Sugihara1994(Y = y, E = Euse, theta = tuse) # fit S-mapping for best paramter set
#' plot(spred$Y, spred$Y_hat); abline(a=0, b=1, lty=2) # observed vs. predicted
#' 
#' ## Run filter with "correct" parameter values
#' N = 1e3 # number of particles
#' #based on detful0
#' filterout_det<-particleFilterLL(y, pars=pars_filter, N,
#'                                 detfun = detfun0_sin, procfun = procfun0,
#'                                 dotraceback = TRUE, fulltraceback = TRUE)
#' #based on EDM
#' filterout_edm<-particleFilterLL(y, pars=pars_filter, N, detfun = EDMfun0,
#'                                 edmdat = list(E=Euse, theta=tuse),
#'                                 procfun = procfun0, dotraceback = TRUE,
#'                                 fulltraceback = TRUE)
#'
#' #plot filter output
#' op = par(mar=c(4,4,2,2), mfrow=c(3,1))
#' #plot 30 of the 1000 particles to show general trend
#' # correct deterministic function
#' matplot(1:length(y), filterout_det$fulltracemat[,1:30],
#'         col=adjustcolor(1,alpha.f = 0.5), lty=3,
#'         type="l", xlab="time", ylab="abund", main="detfun0")
#' lines(1:length(y), y, col=2, lwd=1.5) #observations
#' lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2) #true values
#' lines(1:length(y), filterout_det$Nest, col=3, lwd=1.5) #mean filter estimate
#'
#' # EDM function
#' matplot(1:length(y), filterout_edm$fulltracemat[,1:30],
#'         col=adjustcolor(1,alpha.f = 0.5), lty=3,
#'         type="l", xlab="time", ylab="abund", main="EDM")
#' lines(1:length(y), y, col=2, lwd=1.5)
#' lines(1:length(y), datout$true, col="blue", lwd=1.5, lty=2)
#' lines(1:length(y), filterout_edm$Nest, col=3, lwd=1.5)
#'
#' plot(filterout_det$Nest, datout$true, xlim=range(c(filterout_det$Nest,
#'      filterout_edm$Nest)),
#'      xlab="predicted", ylab="true value", col=4)
#' points(filterout_edm$Nest, datout$true, col=2)
#' points(y, datout$true, col=3)
#' abline(a=0, b=1, lty=2)
#' legend("topleft", c("detfun0", "EDM", "obs"), pch=1, col=c(4,2,3), bty="n")
#'
#' #note improvement in fit, for both filters
#' cor(datout$true, datout$obs)^2 #observations
#' cor(datout$true, filterout_det$Nest)^2 #deterministic filter
#' cor(datout$true, filterout_edm$Nest)^2 #EDM filter
#' par(op) # reset plotting parameters
#' 
#' \dontrun{
#' # Commented-out code: Install BayesianTools if needed
#' #install.packages("BayesianTools")
#' require(BayesianTools)
#' ## Run optimizers
#' #create priors
#' minvUSE<-c(-4, -4) #minimum interval for obs and proc
#' maxvUSE<-c(0, 0) #maximum interval for obs and proc
#'
#' minvUSE_edm<-c(-4, -4) #minimum interval for obs and proc
#' maxvUSE_edm<-c(0, 0) #maximum interval for obs and proc
#'
#' #density, sampler, and prior functions for deterministic function
#' density_fun_USE<-function(param) density_fun0(param = param, minv = minvUSE,
#'                          maxv=maxvUSE)
#' sampler_fun_USE<-function(x) sampler_fun0(n = 1, minv = minvUSE,
#'                          maxv=maxvUSE)
#' prior_USE <- createPrior(density = density_fun_USE,
#'                          sampler = sampler_fun_USE,
#'                          lower = minvUSE, upper = maxvUSE)
#'
#' #density, sampler, and prior functions for EDM function
#' density_fun_USE_edm<-function(param) density_fun0(param = param,
#'                          minv = minvUSE_edm, maxv=maxvUSE_edm)
#' sampler_fun_USE_edm<-function(x) sampler_fun0(n = 1, minv = minvUSE_edm,
#'                          maxv=maxvUSE_edm)
#' prior_edm <- createPrior(density = density_fun_USE_edm,
#'                          sampler = sampler_fun_USE_edm,
#'                          lower = minvUSE_edm, upper = maxvUSE_edm)
#' ## Run filter
#' niter<-5000 #number of steps for the MCMC sampler
#' N<-1e3 #number of particles
#' Euse<-Euse #number of embedding dimensions
#'
#' #likelihood and bayesian set-ups for deterministic functions
#' likelihood_detfun0<-function(x) likelihood0(param=x, y=y,
#'                          parseparam = parseparam0, detfun = detfun0_sin,
#'                          procfun = procfun0, N = N)
#' bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0,
#'                          prior = prior_USE)
#'
#' #likelihood and bayesian set-ups for EDM functions
#' likelihood_EDM<-function(x) {
#'   likelihood0(param = x, y=y, parseparam = parseparam0, procfun = procfun0,
#'               detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse), N = N)
#' }
#'
#' bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM,
#'                          prior = prior_edm)
#' #run MCMC optimization
#' out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0,
#'                    settings = list(iterations=niter, consoleUpdates=20))
#' out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM,
#'                    settings = list(iterations=niter, consoleUpdates=20))
#'
#' #plot results, with a 1000-step burn-in
#' plot(out_detfun0, start = 1000, thin = 2)
#' plot(out_EDM, start = 1000, thin = 2)
#'
#' ## extract and plot parameter distributions
#' smp_detfun0<-getSample(out_detfun0, start = 1000, thin = 2)
#' smp_EDM<-getSample(out_EDM, start=1000, thin = 2)
#'
#' op = par(mfrow=c(2,2))
#' hist(exp(smp_detfun0[,1]), xlim=c(exp(minvUSE[1]), exp(maxvUSE[1])),
#'                     main="det. function", xlab="obs", breaks = 20)
#' abline(v=exp(pars_true$obs), col=2) # true value
#' abline(v=c(exp(minvUSE[1]), exp(maxvUSE[1])), col=1, lty=2)# Priors
#'
#' hist(exp(smp_detfun0[,2]), xlim=c(exp(minvUSE[2]), exp(maxvUSE[2])),
#'                     main="det. function", xlab="proc", breaks = 20)
#' abline(v=exp(pars_true$proc), col=2) # true value
#' abline(v=c(exp(minvUSE[2]), exp(maxvUSE[2])), col=1, lty=2)# Priors
#'
#' hist(exp(smp_EDM[,1]), xlim=c(exp(minvUSE_edm[1]), exp(maxvUSE_edm[1])),
#'                     main="EDM function", xlab="obs", breaks = 20)
#' abline(v=exp(pars_true$obs), col=2) # true value
#' abline(v=c(exp(minvUSE_edm[1]), exp(maxvUSE_edm[1])), col=1, lty=2)# Priors
#'
#' hist(exp(smp_EDM[,2]), xlim=c(exp(minvUSE_edm[2]), exp(maxvUSE_edm[2])),
#'                     main="EDM function", xlab="proc", breaks = 20)
#' abline(v=exp(pars_true$proc), col=2) # true value
#' abline(v=c(exp(minvUSE_edm[2]), exp(maxvUSE_edm[2])), col=1, lty=2)# Priors
#' # note slight over-estimate of process noise due to EDM error
#' par(op) # reset plotting parameters
#' }
NULL

