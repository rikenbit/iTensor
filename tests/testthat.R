library("rTensor")
library("testthat")

options(testthat.use_colours = FALSE)

# source("../R/ICA.R") # comment out
# source("../R/GroupICA.R") # comment out
# source("../R/MultilinearICA.R") # comment out

test_file("testthat/test_ICA.R")
test_file("testthat/test_GroupICA.R")
test_file("testthat/test_MultilinearICA.R")
