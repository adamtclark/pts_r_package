rm(list = ls())

require(devtools)
require(roxygen2)

#setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/")
#create("pttstability")

setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/pts_r_package/pttstability/")
document()
check(cran = TRUE)

build()
load_all()

