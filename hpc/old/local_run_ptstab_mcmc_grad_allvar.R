error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/hpc/")

for(i in 1:100) {
  system(paste("./run_ptstab_mcmc_grad_allvar.R", i))
}
