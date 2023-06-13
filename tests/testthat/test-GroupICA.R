# Simulation Datasets
d <- 6
X1 <- matrix(runif(20 * d), nrow = 20, ncol = d)
X2 <- matrix(runif(15 * d), nrow = 15, ncol = d)
X3 <- matrix(runif(25 * d), nrow = 25, ncol = d)

Xs <- list(X1 = X1, X2 = X2, X3 = X3)

## Perform GroupICA against Simulation Datasets
J1 <- d
set.seed(12345)
out.pooled <- GroupICA(Xs, J1=J1, algorithm="pooled", ica.algorithm = "FastICA", verbose=TRUE)
set.seed(12345)
out.Calhoun2009 <- GroupICA(Xs, J1=J1, algorithm="Calhoun2009", ica.algorithm = "FastICA", verbose=TRUE)
set.seed(12345)
out.Pfister2018 <- suppressWarnings(GroupICA(Xs, J1=J1, algorithm="Pfister2018", ica.algorithm = "FastICA", verbose=TRUE))

set.seed(12345)
out.pooled_2 <- GroupICA(Xs, J1=J1, algorithm="pooled", ica.algorithm = "JADE", verbose=TRUE)
set.seed(12345)
out.Calhoun2009_2 <- GroupICA(Xs, J1=J1, algorithm="Calhoun2009", ica.algorithm = "JADE", verbose=TRUE)
set.seed(12345)
out.Pfister2018_2 <- suppressWarnings(GroupICA(Xs, J1=J1, algorithm="Pfister2018", ica.algorithm = "JADE", verbose=TRUE))

# Test Hidden Functions
set.seed(12345)
out.pooled_ <- iTensor:::.pooled(Xs, J1, J2=J1,
    ica.algorithm="FastICA",
    num.iter=30, thr=1E-10, verbose=FALSE)
set.seed(12345)
out.Calhoun2009_ <- iTensor:::.Calhoun2009(Xs, J1, J2=J1,
    ica.algorithm="FastICA",
    num.iter=30, thr=1E-10, verbose=FALSE)
set.seed(12345)
out.Pfister2018_ <- suppressWarnings(iTensor:::.Pfister2018(Xs, J1, J2=J1,
    ica.algorithm="FastICA",
    num.iter=30, thr=1E-10, verbose=FALSE))
expect_identical(out.pooled, out.pooled_)
expect_identical(out.Calhoun2009, out.Calhoun2009_)
expect_identical(out.Pfister2018, out.Pfister2018_)
expect_true(is.numeric(iTensor:::.eachidx(J1, 1, 1)))
expect_true(is.numeric(iTensor:::.generate_subgroups(1, 1)))

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(out.pooled), TRUE)
expect_identical(is.list(out.Calhoun2009), TRUE)
expect_identical(is.list(out.Pfister2018), TRUE)
expect_identical(is.list(out.pooled_2), TRUE)
expect_identical(is.list(out.Calhoun2009_2), TRUE)
expect_identical(is.list(out.Pfister2018_2), TRUE)

### Test O-2: Object Names
expect_identical(names(out.pooled),
    c("A", "Ss", "RecError", "RelChange"))
expect_identical(names(out.Calhoun2009),
    c("A", "Ss", "RecError", "RelChange"))
expect_identical(names(out.Pfister2018),
    c("A", "Ss", "RecError", "RelChange"))
expect_identical(names(out.pooled_2),
    c("A", "Ss", "RecError", "RelChange"))
expect_identical(names(out.Calhoun2009_2),
    c("A", "Ss", "RecError", "RelChange"))
expect_identical(names(out.Pfister2018_2),
    c("A", "Ss", "RecError", "RelChange"))

### Test O-3: A
expect_equal(
    dim(out.pooled$A),
    c(ncol(Xs[[1]]), J1))
expect_equal(
    dim(out.Calhoun2009$A),
    c(ncol(Xs[[1]]), J1))
expect_equal(
    dim(out.Pfister2018$A),
    c(ncol(Xs[[1]]), J1))
expect_equal(
    dim(out.pooled_2$A),
    c(ncol(Xs[[1]]), J1))
expect_equal(
    dim(out.Calhoun2009_2$A),
    c(ncol(Xs[[1]]), J1))
expect_equal(
    dim(out.Pfister2018_2$A),
    c(ncol(Xs[[1]]), J1))

expect_true(is.matrix(out.pooled$A))
expect_true(is.matrix(out.Calhoun2009$A))
expect_true(is.matrix(out.Pfister2018$A))

expect_true(is.matrix(out.pooled_2$A))
expect_true(is.matrix(out.Calhoun2009_2$A))
expect_true(is.matrix(out.Pfister2018_2$A))

expect_equal(dim(out.pooled$A), c(ncol(Xs[[1]]), J1))
expect_equal(dim(out.Calhoun2009$A), c(ncol(Xs[[1]]), J1))
expect_equal(dim(out.Pfister2018$A), c(ncol(Xs[[1]]), J1))
expect_equal(dim(out.pooled_2$A), c(ncol(Xs[[1]]), J1))
expect_equal(dim(out.Calhoun2009_2$A), c(ncol(Xs[[1]]), J1))
expect_equal(dim(out.Pfister2018_2$A), c(ncol(Xs[[1]]), J1))

### Test O-4: Ss
expect_identical(is.matrix(out.pooled$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.pooled$Ss[[2]]), TRUE)
expect_identical(is.matrix(out.Calhoun2009$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.Calhoun2009$Ss[[2]]), TRUE)
expect_identical(is.matrix(out.Pfister2018$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.Pfister2018$Ss[[2]]), TRUE)
expect_identical(is.matrix(out.pooled_2$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.pooled_2$Ss[[2]]), TRUE)
expect_identical(is.matrix(out.Calhoun2009_2$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.Calhoun2009_2$Ss[[2]]), TRUE)
expect_identical(is.matrix(out.Pfister2018_2$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.Pfister2018_2$Ss[[2]]), TRUE)

expect_equal(dim(out.pooled$Ss[[1]]), c(nrow(Xs[[1]]), J1))
expect_equal(dim(out.pooled$Ss[[2]]), c(nrow(Xs[[2]]), J1))
expect_equal(dim(out.Calhoun2009$Ss[[1]]), c(nrow(Xs[[1]]), J1))
expect_equal(dim(out.Calhoun2009$Ss[[2]]), c(nrow(Xs[[2]]), J1))
expect_equal(dim(out.Pfister2018$Ss[[1]]), c(nrow(Xs[[1]]), J1))
expect_equal(dim(out.Pfister2018$Ss[[2]]), c(nrow(Xs[[2]]), J1))
expect_equal(dim(out.pooled_2$Ss[[1]]), c(nrow(Xs[[1]]), J1))
expect_equal(dim(out.pooled_2$Ss[[2]]), c(nrow(Xs[[2]]), J1))
expect_equal(dim(out.Calhoun2009_2$Ss[[1]]), c(nrow(Xs[[1]]), J1))
expect_equal(dim(out.Calhoun2009_2$Ss[[2]]), c(nrow(Xs[[2]]), J1))
expect_equal(dim(out.Pfister2018_2$Ss[[1]]), c(nrow(Xs[[1]]), J1))
expect_equal(dim(out.Pfister2018_2$Ss[[2]]), c(nrow(Xs[[2]]), J1))

### Test O-5: RecError
expect_identical(is.null(out.pooled$RecError), TRUE)
expect_identical(is.null(out.Calhoun2009$RecError), TRUE)
expect_identical(is.null(out.Pfister2018$RecError), TRUE)
expect_identical(is.null(out.pooled_2$RecError), FALSE)
expect_identical(is.null(out.Calhoun2009_2$RecError), TRUE)
expect_identical(is.null(out.Pfister2018_2$RecError), TRUE)

# JADE does not return relative change
expect_identical(is.null(out.pooled$RelChange), TRUE)
expect_identical(is.null(out.Calhoun2009$RelChange), TRUE)
expect_identical(is.null(out.Pfister2018$RelChange), TRUE)
expect_identical(is.null(out.pooled_2$RelChange), TRUE)
expect_identical(is.null(out.Calhoun2009_2$RelChange), TRUE)
expect_identical(is.null(out.Pfister2018_2$RelChange), TRUE)

## Test Error
### Test E-1: X
expect_error(GroupICA(X, J1=J1))

### Test E-2: J1
expect_error(GroupICA(Xs, J1="5"))
expect_error(GroupICA(Xs, J1=c(2,4)))
expect_error(GroupICA(Xs, J1=10^10))

### Test E-3: J2
expect_error(GroupICA(Xs, J2="5"))
expect_error(GroupICA(Xs, J2=c(2,4)))
expect_error(GroupICA(Xs, J2=10^10))

### Test E-4: algorithm
expect_error(GroupICA(Xs, algorithm="poooled"))

### Test E-5: ica.algorithm
expect_error(GroupICA(Xs, ica.algorithm="JAAE"))

### Test E-6: num.iter
expect_error(GroupICA(Xs, J=J, num_iter="100"))
expect_error(GroupICA(Xs, J=J, num_iter=-1))

### Test E-7: thr
expect_error(GroupICA(X, J=J, thr="0.1"))

### Test E-8: verbose
expect_error(GroupICA(X, J=J, verbose="verbose"))
