# Simulation Dataset
X <- matrix(runif(100*200), nrow=100, ncol=200)

## Perform ICA against Simulation Dataset
J <- 5
out.FastICA <- ICA(X, J=J, algorithm="FastICA")
out.InfoMax <- ICA(X, J=J, algorithm="InfoMax")
out.ExtInfoMax <- ICA(X, J=J, algorithm="ExtInfoMax")

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(out.FastICA), TRUE)
expect_identical(is.list(out.InfoMax), TRUE)
expect_identical(is.list(out.ExtInfoMax), TRUE)

### Test O-2: Object Names
expect_identical(names(out.FastICA),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "WChange"))
expect_identical(names(out.InfoMax),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "WChange"))
expect_identical(names(out.ExtInfoMax),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "WChange"))

### Test 0-3: A
expect_identical(is.matrix(out.FastICA$A), TRUE)
expect_identical(is.matrix(out.InfoMax$A), TRUE)
expect_identical(is.matrix(out.ExtInfoMax$A), TRUE)

expect_equal(dim(out.FastICA$A), c(J, ncol(X)))
expect_equal(dim(out.InfoMax$A), c(J, ncol(X)))
expect_equal(dim(out.ExtInfoMax$A), c(J, ncol(X)))

### Test 0-4: S
expect_identical(is.matrix(out.FastICA$S), TRUE)
expect_identical(is.matrix(out.InfoMax$S), TRUE)
expect_identical(is.matrix(out.ExtInfoMax$S), TRUE)

expect_equal(dim(out.FastICA$S), c(nrow(X), J))
expect_equal(dim(out.InfoMax$S), c(nrow(X), J))
expect_equal(dim(out.ExtInfoMax$S), c(nrow(X), J))

### Test 0-5: J
expect_identical(out.FastICA$J, J)
expect_identical(out.InfoMax$J, J)
expect_identical(out.ExtInfoMax$J, J)

### Test 0-6: algorithm
expect_identical(out.FastICA$algorithm, "FastICA")
expect_identical(out.InfoMax$algorithm, "InfoMax")
expect_identical(out.ExtInfoMax$algorithm, "ExtInfoMax")

### Test 0-7: num.iter
expect_identical(out.FastICA$num.iter, formals(ICA)$num.iter)
expect_identical(out.InfoMax$num.iter, formals(ICA)$num.iter)
expect_identical(out.ExtInfoMax$num.iter, formals(ICA)$num.iter)

### Test 0-8: thr
expect_identical(out.FastICA$thr, formals(ICA)$thr)
expect_identical(out.InfoMax$thr, formals(ICA)$thr)
expect_identical(out.ExtInfoMax$thr, formals(ICA)$thr)

### Test 0-9: verbose
expect_identical(out.FastICA$verbose, formals(ICA)$verbose)
expect_identical(out.InfoMax$verbose, formals(ICA)$verbose)
expect_identical(out.ExtInfoMax$verbose, formals(ICA)$verbose)

### Test 0-10: RecError
expect_identical(is.vector(out.FastICA$WChange), TRUE)
expect_identical(is.vector(out.InfoMax$WChange), TRUE)
expect_identical(is.vector(out.ExtInfoMax$WChange), TRUE)

## Test Error
### Test E-1: X
expect_error(ICA(as.data.frame(X), J=J))

### Test E-2: J
expect_error(ICA(X, J="5"))
expect_error(ICA(X, J=c(2,4)))
expect_error(ICA(X, J=10^10))

### Test E-3: num.iter
expect_error(ICA(X, J=J, num.iter="100"))
expect_error(ICA(X, J=J, num.iter=-1))

### Test E-4: thr
expect_error(ICA(X, J=J, thr="0.1"))

### Test E-5: verbose
expect_error(ICA(X, J=J, verbose="verbose"))

## Test Orthgonalization

orthgonalized <- iTensor:::.orthgonalize(
    matrix(c(8, 1, 9, 8, 1, 4, 6, 0, 5), 3, 3))
expect_equal((orthgonalized[1, ] %*% orthgonalized[1, ])[1], 1, tolerance=1e-8)
expect_equal((orthgonalized[1, ] %*% orthgonalized[2, ])[1], 0, tolerance=1e-8)
expect_equal((orthgonalized[1, ] %*% orthgonalized[3, ])[1], 0, tolerance=1e-8)
expect_equal((orthgonalized[2, ] %*% orthgonalized[2, ])[1], 1, tolerance=1e-8)
expect_equal((orthgonalized[2, ] %*% orthgonalized[3, ])[1], 0, tolerance=1e-8)
