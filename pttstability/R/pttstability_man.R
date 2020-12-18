#' pttstability: A Particle-Takens Filter for Measuring Stability of Nonstationary Systems.
#'
#' The pttstability (or "ParTicle-Takens Stability" package) is a collection of functions
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
#' Second, it applies empirical dynamic modeling (EDM) methods from the rEDM package to estimate deterministic
#' dynamics even in cases where the underlying equations governing system behavior are not known.
#' These models are based on Takens delay-embedding theorem, from which this package takes its name.
#' Finally, it uses the MCMC fitting methods from the BayesianTools package to estimate paramter values for
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
#' @source Knape, J., and Valpine, P. (2012). Fitting complex population models by combining particle filters with Markov chain Monte Carlo. Ecology 93:256-263.
#' @source Ye, H., Sugihara, G., et al. (2015). Equation-free ecosystem forecasting. PNAS 112:E1569-E1576.
#' @source Ye, H., et al. (2019). rEDM: Applications of Empirical Dynamic Modeling from Time Series. R package version 0.7.2.
#' @source Hartig, F., et al. (2019). BayesianTools: General-Purpose MCMC and SMC Samplers and Tools for Bayesian Statistics. R package version 0.1.6.
#' @name pttstability
#' @examples
#' #Load packages
#' require(rEDM)
#' require(BayesianTools)
#'
#' #Load packages
#' require(rEDM)
#' require(BayesianTools)
#'
#' #Set seed
#' set.seed(201116)
#'
#' ## Simulate data
#' pars_true<-list(obs=log(0.2),
#'                 proc=c(log(0.05), log(1.2)),
#'                 pcol=c(logit(0.2), log(0.1)),
#'                 det=c(log(1.2),log(1)))
#'
#' #generate dynamics
#' datout<-makedynamics_general(n = 100, n0 = 1, pdet=pars_true$det,
#'                              proc = pars_true$proc, obs = pars_true$obs, pcol = pars_true$pcol,
#'                              detfun = detfun0, procfun = procfun0, obsfun=obsfun0, colfun=colfun0)
#' y<-datout$obs # noisy observed abundances
#'
#' #get theta paramter for s-mapping
#' #note - rEDM must be installed to conduct
#' #nonparametric tests
#' s<-rEDM::s_map(y, E=2, silent = TRUE)
#' tuse<-s$theta[which.min(s$rmse)]
#'
#' ## Run filter
#' N = 2e3
#' #based on detful0
#' filterout_det<-particleFilterLL(y, pars=pars_true, N, detfun = detfun0,
#'                                 dotraceback = TRUE, fulltraceback = TRUE)
#' #based on EDM
#' filterout_edm<-particleFilterLL(y, pars=pars_true, N, detfun = EDMfun0,
#'                                 edmdat = list(E=2, theta=tuse),
#'                                 dotraceback = TRUE, fulltraceback = TRUE)
#'
#' #get sorted samples from the particle filters
#' sorted_filter_det<-indexsort(fulltracemat = filterout_det$fulltracemat,
#'                              fulltraceindex = filterout_det$fulltraceindex, nsmp = 100)
#' sorted_filter_edm<-indexsort(fulltracemat = filterout_edm$fulltracemat,
#'                              fulltraceindex = filterout_edm$fulltraceindex, nsmp = 100)
#'
#' #plot filter output
#' par(mar=c(4,4,2,2), mfrow=c(3,1))
#' #plot 100 particles to show general trend
#' matplot(1:length(y), sorted_filter_det, col=adjustcolor("dodgerblue",alpha.f = 0.5), lty=3,
#'         type="l", xlab="time", ylab="abund", main="detfun0")
#' lines(1:length(y), y, col="goldenrod", lwd=1.5) #observations
#' lines(1:length(y), datout$true, col="black", lwd=1.5, lty=2) #true values
#'
#' matplot(1:length(y), sorted_filter_edm, col=adjustcolor("firebrick",alpha.f = 0.5), lty=3,
#'         type="l", xlab="time", ylab="abund", main="EDM")
#' lines(1:length(y), y, col="goldenrod", lwd=1.5)
#' lines(1:length(y), datout$true, col="black", lwd=1.5, lty=2)
#'
#' plot(apply(sorted_filter_det, 1, mean), datout$true, xlim=c(0,2.1), ylim=c(0, 2.1),
#'      xlab="predicted", ylab="true value", col="dodgerblue")
#' points(apply(sorted_filter_edm, 1, mean), datout$true, col="firebrick")
#' points(y, datout$true, col="goldenrod")
#' abline(a=0, b=1, lty=2)
#' legend("topleft", c("detfun", "EDM", "obs"),
#'        pch=1, col=c("dodgerblue","firebrick","goldenrod"), bty="n")
#'
#' #note improvement in mean absolute error, for both filters
#' mean(abs(datout$true-datout$obs)) #observations
#' mean(abs(datout$true-filterout_det$Nest)) #deterministic filter
#' mean(abs(datout$true-filterout_edm$Nest)) #EDM filter
#'
#' ## Run optimizers
#' \dontrun{
#' #create priors
#' minvUSE<-c(log(0.01), log(0.01), log(0.5)) #minimum interval for obs and proc
#' maxvUSE<-c(log(0.5), log(0.5), log(3)) #maximum interval for obs and proc
#'
#' minvUSE_edm<-c(log(0.01), log(0.01), log(0.5)) #minimum interval for obs, proc
#' maxvUSE_edm<-c(log(0.5), log(0.5), log(3)) #maximum interval for obs, proc
#'
#' #density, sampler, and prior functions for deterministic function
#' density_fun_USE<-function(param) density_fun0(param = param, minv = minvUSE, maxv=maxvUSE)
#' sampler_fun_USE<-function(x) sampler_fun0(n = 1, minv = minvUSE, maxv=maxvUSE)
#' prior_USE <- createPrior(density = density_fun_USE, sampler = sampler_fun_USE,
#'                          lower = minvUSE, upper = maxvUSE)
#'
#' #density, sampler, and prior functions for EDM function
#' density_fun_USE_edm<-function(param) density_fun0(param = param,
#'                                 minv = minvUSE_edm, maxv=maxvUSE_edm)
#' sampler_fun_USE_edm<-function(x) sampler_fun0(n = 1, minv = minvUSE_edm, maxv=maxvUSE_edm)
#' prior_edm <- createPrior(density = density_fun_USE_edm, sampler = sampler_fun_USE_edm,
#'                          lower = minvUSE_edm, upper = maxvUSE_edm)
#' ## Run filter
#' niter<-5000 #number of steps for the MCMC sampler
#' N<-2e3 #number of particles
#' Euse<-2 #number of embedding dimensions
#'
#' #likelihood and bayesian set-ups for deterministic functions
#' likelihood_detfun0<-function(x) likelihood0(param=x, y=y, parseparam = parseparam0, N = N)
#' bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0, prior = prior_USE)
#'
#' #likelihood and bayesian set-ups for EDM functions
#' likelihood_EDM<-function(x) {
#'   likelihood0(param = x, y=y, parseparam = parseparam0,
#'               detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse), N = N)
#' }
#'
#' bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM, prior = prior_edm)
#'
#' #run MCMC optimization
#' #~5 min runtime on a regular laptop
#' out_detfun0 <- runMCMC(bayesianSetup = bayesianSetup_detfun0,
#'                                 settings = list(iterations=niter, consoleUpdates=20))
#' #~10 min runtime on a regular laptop
#' out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM,
#'                                 settings = list(iterations=niter, consoleUpdates=20))
#'
#' #plot results, with a 1000-step burn-in
#' plot(out_detfun0, start = 1000)
#' plot(out_EDM, start = 1000)
#'
#' ## extract and plot parameter distributions
#' smp_detfun0<-getSample(out_detfun0, start = 1000)
#' smp_EDM<-getSample(out_EDM, start=1000)
#'
#' par(mfcol=c(3,2))
#' hist(exp(smp_detfun0[,1]), main="det. function", xlab="obs0", breaks=20, xlim=c(0.01, 0.5))
#' abline(v=exp(pars_true$obs), col=2); abline(v=exp(mean(smp_detfun0[,1])), col=4)
#' hist(exp(smp_detfun0[,2]), main="det. function", xlab="proc0", breaks=20, xlim=c(0.01, 0.5))
#' abline(v=exp(pars_true$proc[1]), col=2); abline(v=exp(mean(smp_detfun0[,2])), col=4)
#' hist(exp(smp_detfun0[,3]), main="det. function", xlab="proc1", breaks=20, xlim=c(0.5, 3))
#' abline(v=exp(pars_true$proc[2]), col=2); abline(v=exp(mean(smp_detfun0[,3])), col=4)
#'
#' hist(exp(smp_EDM[,1]), main="EDM function", xlab="obs0", breaks=20, xlim=c(0.01, 0.5))
#' abline(v=exp(pars_true$obs), col=2); abline(v=exp(mean(smp_EDM[,1])), col=4)
#' hist(exp(smp_EDM[,2]), main="EDM function", xlab="proc0", breaks=20, xlim=c(0.01, 0.5))
#' abline(v=exp(pars_true$proc[1]), col=2); abline(v=exp(mean(smp_EDM[,2])), col=4)
#' hist(exp(smp_EDM[,3]), main="EDM function", xlab="proc1", breaks=20, xlim=c(0.5, 3))
#' abline(v=exp(pars_true$proc[2]), col=2); abline(v=exp(mean(smp_EDM[,3])), col=4)
#' }
NULL

