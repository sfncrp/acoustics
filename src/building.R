
library("devtools")
library("roxygen2")

#create("acoustics")

getwd()

setwd("acoustics/")

build()

document()

build()

?lden
?pannoyed

pkg_path <- devtools::build(binary = TRUE, args = c('--preclean'))

#file.copy(pkg_path, )


setwd("../")
install.packages("acoustics_0.0.0.9000.tar.gz")

setwd("..")

devtools::install("acoustics")

devtools::install("acoustics_0.0.0.9000.tar.gz")

check("acoustics")
