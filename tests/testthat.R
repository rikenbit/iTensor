library(testthat)
library(iTensor)
library(rTensor)

options(testthat.use_colours = FALSE)

# Basic usage
test_file("testthat/test-ICA.R")
test_file("testthat/test-ICA2.R")
test_file("testthat/test-MICA.R")
test_file("testthat/test-GroupICA.R")
test_file("testthat/test-MultilinearICA.R")
