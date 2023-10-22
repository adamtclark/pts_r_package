rm(list = ls())

require(devtools)
require(roxygen2)

#create("pttstability")

setwd("pttstability/")
document()
check(cran = TRUE)

build()
load_all()



# I have updated the citation, reset the user's par settings to
# defaults, and have altered the \dontrun examples as requested.
# Apologies for the long delay, but I needed to wait for the DOI to be
# posted by the journal before I could add it to the description file.
# The code still includes two instances of \dontrun, and one
# instance of \donttest. The \dontrun examples (in the main help page,
# and the particleFilterLL_piecewise function) are both optional
# package extensions which use the BayesianTools package, which is
# currently not available on CRAN. Since this software is not available
# in the CRAN checks, I have wrapped these examples in \dontrun. If this
# is not acceptable, then I can also simply remove these examples.
# The only usage of \donttest is for the "particleFilterLL_piecewise" function,
# for running a piece of very slow code - I'm afraid there is no short toy example
# for this segment, as it requires an n2 search over all possible sample size
# inputs - even the example shown here is for just a few lines of toy data.
# Since this example takes more than 2 minutes to run, I have wrapped it in
# \donttest.
