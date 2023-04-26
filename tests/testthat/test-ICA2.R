# Simulation Dataset
X <- matrix(runif(200 * 100), nrow = 200, ncol = 100)

## Perform ICA against Simulation Dataset
J <- 5
out.JADE <- ICA2(X, J = J, algorithm = "JADE")
out.AuxICA1 <- ICA2(X, J = J, algorithm = "AuxICA1")
out.AuxICA2 <- ICA2(X, J = J, algorithm = "AuxICA2")
out.IPCA <- ICA2(X, J = J, algorithm = "IPCA")
out.SIMBEC <- ICA2(X, J = J, algorithm = "SIMBEC")
out.AMUSE <- ICA2(X, J = J, algorithm = "AMUSE")
out.SOBI <- ICA2(X, J = J, algorithm = "SOBI")
out.FOBI <- ICA2(X, J = J, algorithm = "FOBI")
out.ProDenICA <- ICA2(X, J = J, algorithm = "ProDenICA", num.iter = 1)
out.RICA <- ICA2(X, J = J, algorithm = "RICA", num.iter = 1)

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(out.JADE), TRUE)
expect_identical(is.list(out.AuxICA1), TRUE)
expect_identical(is.list(out.AuxICA2), TRUE)
expect_identical(is.list(out.IPCA), TRUE)
expect_identical(is.list(out.SIMBEC), TRUE)
expect_identical(is.list(out.AMUSE), TRUE)
expect_identical(is.list(out.SOBI), TRUE)
expect_identical(is.list(out.FOBI), TRUE)
expect_identical(is.list(out.ProDenICA), TRUE)
expect_identical(is.list(out.RICA), TRUE)

### Test O-2: Object Names
expect_identical(
    sort(names(out.JADE)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.AuxICA1)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.AuxICA2)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.IPCA)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.SIMBEC)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.AMUSE)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.SOBI)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.FOBI)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.ProDenICA)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))
expect_identical(
    sort(names(out.RICA)),
    sort(c("A", "W", "S", "X", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange")))

### Test 0-3: A
expect_identical(is.matrix(out.JADE$A), TRUE)
expect_identical(is.matrix(out.AuxICA1$A), TRUE)
expect_identical(is.matrix(out.AuxICA2$A), TRUE)
expect_identical(is.matrix(out.IPCA$A), TRUE)
expect_identical(is.matrix(out.SIMBEC$A), TRUE)
expect_identical(is.matrix(out.AMUSE$A), TRUE)
expect_identical(is.matrix(out.SOBI$A), TRUE)
expect_identical(is.matrix(out.FOBI$A), TRUE)
expect_identical(is.matrix(out.ProDenICA$A), TRUE)
expect_identical(is.matrix(out.RICA$A), TRUE)

# Except for IPCA, A is a square matrix because ICA is performed after decorrelating.
expect_equal(dim(out.JADE$A), c(J, J))
expect_equal(dim(out.AuxICA1$A), c(J, J))
expect_equal(dim(out.AuxICA2$A), c(J, J))
expect_equal(dim(out.IPCA$A), c(J, ncol(X)))
expect_equal(dim(out.SIMBEC$A), c(J, J))
expect_equal(dim(out.AMUSE$A), c(J, J))
expect_equal(dim(out.SOBI$A), c(J, J))
expect_equal(dim(out.FOBI$A), c(J, J))
expect_equal(dim(out.ProDenICA$A), c(J, J))
expect_equal(dim(out.RICA$A), c(J, J))

### Test 0-4: S
expect_identical(is.matrix(out.JADE$S), TRUE)
expect_identical(is.matrix(out.AuxICA1$S), TRUE)
expect_identical(is.matrix(out.AuxICA2$S), TRUE)
expect_identical(is.matrix(out.IPCA$S), TRUE)
expect_identical(is.matrix(out.SIMBEC$S), TRUE)
expect_identical(is.matrix(out.AMUSE$S), TRUE)
expect_identical(is.matrix(out.SOBI$S), TRUE)
expect_identical(is.matrix(out.FOBI$S), TRUE)
expect_identical(is.matrix(out.ProDenICA$S), TRUE)
expect_identical(is.matrix(out.RICA$S), TRUE)

expect_equal(dim(out.JADE$S), c(nrow(X), J))
expect_equal(dim(out.AuxICA1$S), c(nrow(X), J))
expect_equal(dim(out.AuxICA2$S), c(nrow(X), J))
expect_equal(dim(out.IPCA$S), c(nrow(X), J))
expect_equal(dim(out.SIMBEC$S), c(nrow(X), J))
expect_equal(dim(out.AMUSE$S), c(nrow(X), J))
expect_equal(dim(out.SOBI$S), c(nrow(X), J))
expect_equal(dim(out.FOBI$S), c(nrow(X), J))
expect_equal(dim(out.ProDenICA$S), c(nrow(X), J))
expect_equal(dim(out.RICA$S), c(nrow(X), J))

### Test 0-5: J
expect_identical(out.JADE$J, J)
expect_identical(out.AuxICA1$J, J)
expect_identical(out.AuxICA2$J, J)
expect_identical(out.IPCA$J, J)
expect_identical(out.SIMBEC$J, J)
expect_identical(out.AMUSE$J, J)
expect_identical(out.SOBI$J, J)
expect_identical(out.FOBI$J, J)
expect_identical(out.ProDenICA$J, J)
expect_identical(out.RICA$J, J)

### Test 0-6: algorithm
expect_identical(out.JADE$algorithm, "JADE")
expect_identical(out.AuxICA1$algorithm, "AuxICA1")
expect_identical(out.AuxICA2$algorithm, "AuxICA2")
expect_identical(out.IPCA$algorithm, "IPCA")
expect_identical(out.SIMBEC$algorithm, "SIMBEC")
expect_identical(out.AMUSE$algorithm, "AMUSE")
expect_identical(out.SOBI$algorithm, "SOBI")
expect_identical(out.FOBI$algorithm, "FOBI")
expect_identical(out.ProDenICA$algorithm, "ProDenICA")
expect_identical(out.RICA$algorithm, "RICA")

### Test 0-7: num.iter
# For iterative algorithm, num.iter matches set value of num.iter.
expect_identical(out.AuxICA1$num.iter, 30)
expect_identical(out.AuxICA2$num.iter, 30)

# For ProDenICA and RICA, num.iter are configured to 1.
expect_identical(out.ProDenICA$num.iter, 1)
expect_identical(out.RICA$num.iter, 1)

# For non-iterative algorithm, num.iter are configured to 1.
expect_identical(out.JADE$num.iter, 1)
expect_identical(out.IPCA$num.iter, 1)
expect_identical(out.SIMBEC$num.iter, 1)
expect_identical(out.AMUSE$num.iter, 1)
expect_identical(out.SOBI$num.iter, 1)
expect_identical(out.FOBI$num.iter, 1)

### Test 0-8: thr
expect_identical(out.JADE$thr, NULL)
expect_identical(out.AuxICA1$thr, formals(ICA)$thr)
expect_identical(out.AuxICA2$thr, formals(ICA)$thr)
expect_identical(out.IPCA$thr, NULL)
expect_identical(out.SIMBEC$thr, NULL)
expect_identical(out.AMUSE$thr, NULL)
expect_identical(out.SOBI$thr, NULL)
expect_identical(out.FOBI$thr, NULL)
expect_identical(out.ProDenICA$thr, formals(ICA)$thr)
expect_identical(out.RICA$thr, formals(ICA)$thr)

### Test 0-9: verbose
expect_identical(out.JADE$verbose, formals(ICA)$verbose)
expect_identical(out.AuxICA1$verbose, formals(ICA)$verbose)
expect_identical(out.AuxICA2$verbose, formals(ICA)$verbose)
expect_identical(out.IPCA$verbose, formals(ICA)$verbose)
expect_identical(out.SIMBEC$verbose, formals(ICA)$verbose)
expect_identical(out.AMUSE$verbose, formals(ICA)$verbose)
expect_identical(out.SOBI$verbose, formals(ICA)$verbose)
expect_identical(out.FOBI$verbose, formals(ICA)$verbose)
expect_identical(out.ProDenICA$verbose, formals(ICA)$verbose)
expect_identical(out.RICA$verbose, formals(ICA)$verbose)

### Test 0-10: RecError
expect_identical(is.vector(out.JADE$RecError), TRUE)
expect_identical(is.vector(out.AuxICA1$RecError), TRUE)
expect_identical(is.vector(out.AuxICA2$RecError), TRUE)
expect_identical(is.vector(out.IPCA$RecError), TRUE)
expect_identical(is.vector(out.SIMBEC$RecError), TRUE)
expect_identical(is.vector(out.AMUSE$RecError), TRUE)
expect_identical(is.vector(out.SOBI$RecError), TRUE)
expect_identical(is.vector(out.FOBI$RecError), TRUE)
expect_identical(is.vector(out.ProDenICA$RecError), TRUE)
expect_identical(is.vector(out.RICA$RecError), TRUE)

### Test O-11: RelChange
expect_identical(is.vector(out.AuxICA1$RelChange), TRUE)
expect_identical(is.vector(out.AuxICA2$RelChange), TRUE)
expect_identical(is.vector(out.ProDenICA$RelChange), TRUE)
expect_identical(is.vector(out.RICA$RelChange), TRUE)

## Test Error
### Test E-1: X
expect_error(ICA2(as.data.frame(X), J=J))

### Test E-2: J
expect_error(ICA2(X, J="5"))
expect_error(ICA2(X, J=c(2,4)))
expect_error(ICA2(X, J=10^10))

### Test E-3: num.iter
expect_error(ICA2(X, J=J, num.iter="100"))
expect_error(ICA2(X, J=J, num.iter=-1))

### Test E-4: thr
expect_error(ICA2(X, J=J, thr="0.1"))

### Test E-5: verbose
expect_error(ICA2(X, J=J, verbose="verbose"))
