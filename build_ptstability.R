rm(list = ls())

require(devtools)
require(roxygen2)

#create("pttstability")

setwd("~/Dropbox/GitProjects/R_packages/pts_r_package/pttstability/")
document()
check(cran = TRUE)

build()
load_all()


# The code still includes two instances of \dontrun, and one
# instance of \donttest. The \dontrun examples (in the main help page,
# and the particleFilterLL_piecewise function) are both optional
# package extensions which use the recommended BayesianTools package. If this
# is not acceptable, then I can also simply remove these examples.
# The only usage of \donttest is for the "particleFilterLL_piecewise" function,
# for running a piece of very slow code - I'm afraid there is no short toy example
# for this segment, as it requires an n2 search over all possible sample size
# inputs - even the example shown here is for just a few lines of toy data.
# Since this example takes more than 2 minutes to run, I have wrapped it in
# \donttest.
