# Simulation Datasets
X <- array(runif(10*20), dim=c(10,20))
Y <- array(runif(15*20), dim=c(15,20))

## Perform MICA against Simulation Datasets
J <- 20
out <- MICA(X, Y, J=J)

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(out), TRUE)

### Test O-2: Object Names
expect_identical(
    names(out),
    c("U", "V", "A", "B", "J", "eta", "num.iter", "verbose", "ABChange"))

### Test O-3: U
expect_identical(is.matrix(out$U), TRUE)
expect_equal(dim(out$U), c(nrow(X), J))

### Test O-4: V
expect_identical(is.matrix(out$V), TRUE)
expect_equal(dim(out$V), c(nrow(Y), J))

### Test O-5: A
expect_identical(is.matrix(out$A), TRUE)
expect_equal(dim(out$A), c(ncol(X), J))

### Test O-6: B
expect_identical(is.matrix(out$B), TRUE)
expect_equal(dim(out$B), c(ncol(Y), J))

### Test O-7: J
expect_identical(out$J, J)

### Test 0-8: eta
expect_identical(out$eta, 1000*1e-4)

### Test O-9: verbose
expect_identical(out$verbose, formals(MICA)$verbose)

### Test O-10: RecError
expect_identical(is.vector(out$ABChange), TRUE)

## Test Orthgonalization
orthgonalized <- iTensor:::.orthgonalize(
     matrix(c(8, 1, 9, 8, 1, 4, 6, 0, 5), 3, 3))
expect_equal((orthgonalized[1, ] %*% orthgonalized[1, ])[1], 1, tolerance=1e-8)
expect_equal((orthgonalized[1, ] %*% orthgonalized[2, ])[1], 0, tolerance=1e-8)
expect_equal((orthgonalized[1, ] %*% orthgonalized[3, ])[1], 0, tolerance=1e-8)
expect_equal((orthgonalized[2, ] %*% orthgonalized[2, ])[1], 1, tolerance=1e-8)
expect_equal((orthgonalized[2, ] %*% orthgonalized[3, ])[1], 0, tolerance=1e-8)

## Test Error
### Test E-1: X
expect_error(MICA(X, J=J))
expect_error(MICA(Y, J=J))
expect_error(MICA(as.data.frame(X), Y, J=J))
expect_error(MICA(X, as.data.frame(Y), J=J))

### Test E-2: J
expect_error(MICA(X, Y, eta="5"), regexp = "J")
expect_error(MICA(X, Y, eta=c(2,4)), regexp = "J")
expect_error(MICA(X, Y, eta=10^10))

### Test E-3: eta
expect_error(MICA(X, Y, J=J, eta="5"), regexp = "is.numeric(eta)*")
expect_error(MICA(X, Y, J=J, eta=c(2,4)), regexp = "length")

### Test E-4: verbose
expect_error(MICA(X, Y, J=J, verbose="verbose"), regexp = "is.logical")
