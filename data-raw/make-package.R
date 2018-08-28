library(devtools)

wd <- '/home/pastor/projects/packages/skima'
setwd(wd)

devtools::create(wd)

devtools::use_data_raw()

devtools::document() # uses roxygen to compile the documentation

install.packages(devtools::build()) # produces a compressed version of the package

devtools::use_vignette('user_guide')

