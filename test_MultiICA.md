# Test code of MultiICA Function
## Dependent Packages Installation
if (!requireNamespace("testthat", quietly = TRUE)){
    install.packages("testthat")
}
if (!requireNamespace("rTensor", quietly = TRUE)){
    install.packages("rTensor")
}
if (!requireNamespace("nnTensor", quietly = TRUE)){
    install.packages("nnTensor")
}
library("testthat")
library("rTensor")
library("nnTensor")

## MultiICA Function (Paste your MultiICA code here)

MultiICA <- function(X, J=NULL, modes=1:3, algorithm=c("FastICA", "InfoMax", "ExtInfoMax")){
    ######################################
    # Argument Check
    ######################################
    .checkMultiICA(X, J, modes)
    algorithm <- <- match.arg(algorithm)
    ######################################
    # Initialization
    ######################################
    int <- .initMultiICA(X, J, modes)
    X <- int$X
    A <- int$A
    S <- int$S
    RecError <- int$RecError
    for(m in modes){
        Xm <- unfold(X, m=m)
        A[m] <- ICA(X=Xm, J=J[m], algorithm=algorithm)$A
    }
    # After Update
    S <- .Projection(X, A, ms=modes)
    X_bar <- ttl(S, A, ms=modes)
    RecError <- .recError(X, X_bar)
	# Output
    list(A=A, S=S, J=J, algorithm=algorithm, RecError=RecError)
}

.ndim <- function(X){
    length(dim(X))
}

.checkMultiICA <- function(X, J, modes){
    stopifnot(is.array(X))
    stopifnot(is.vector(J))
    stopifnot(.ndim(X) == length(J))
    lapply(seq_along(J), function(i){
        stopifnot(dim(X)[i] >= J[i])
    })
    stopifnot(all(1 <= modes))
    stopifnot(all(.ndim(X) >= modes))
    stopifnot(all(.ndim(X) >= length(modes)))
}

.initMultiICA <- function(X, J, modes){
	# A
	nrs <- dim(X)
	A <- lapply(seq_along(nrs), function(x){
        matrix(runif(nrs[x]*J[x]), nrow=nrs[x], ncol=J[x])
    })
    # S
	S <- .Projection(X,  A, ms=modes)
	list(X=X, A=A, S=S)
}

.Projection <- function(X, A, ms=modes){
    # XにginvしたAnをかける
    ...
}

## Simulation Datasets
X <- array(runif(10*20*30), dim=c(10,20,30))

## Perform MultiICA against Simulation Datasets
J <- c(3,4,5)
out.FastICA <- MultiICAICA(X, J=J, algorithm="FastICA")
out.InfoMax <- MultiICAICA(X, J=J, algorithm="InfoMax")
out.ExtInfoMax <- MultiICAICA(X, J=J, algorithm="ExtInfoMax")

## Test Input object / type
### Test I-1: Object Names
expect_identical(names(formals(MultiICAICA)), c("X", "J", "modes", "algorithm"))
### Test I-2: X
expect_identical(as.character(formals(ICA)$X), "")
### Test I-3: J
expect_identical(as.character(formals(ICA)$J), "")
### Test I-4: modes
expect_identical(as.character(formals(ICA)$modes), 1:3)
### Test I-5: algorithm
expect_identical(as.character(formals(ICA)$algorithm), c("FastICA", "InfoMax", "ExtInfoMax"))

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(out.FastICA), TRUE)
expect_identical(is.list(out.InfoMax), TRUE)
expect_identical(is.list(out.ExtInfoMax), TRUE)

### Test O-2: Object Names
expect_identical(names(out.FastICA),
c("A", "S", "J", "algorithm", "RecError"))
expect_identical(names(out.InfoMax),
c("A", "S", "J", "algorithm", "RecError"))
expect_identical(names(out.ExtInfoMax),
c("A", "S", "J", "algorithm", "RecError"))

### Test 0-3: A
expect_identical(is.list(out.FastICA$A), TRUE)
expect_identical(is.list(out.InfoMax$A), TRUE)
expect_identical(is.list(out.ExtInfoMax$A), TRUE)

expect_identical(is.matrix(out.FastICA$A[[1]]), TRUE)
expect_identical(is.matrix(out.InfoMax$A[[1]]), TRUE)
expect_identical(is.matrix(out.ExtInfoMax$A[[1]]), TRUE)

expect_identical(is.matrix(out.FastICA$A[[2]]), TRUE)
expect_identical(is.matrix(out.InfoMax$A[[2]]), TRUE)
expect_identical(is.matrix(out.ExtInfoMax$A[[2]]), TRUE)

expect_identical(is.matrix(out.FastICA$A[[3]]), TRUE)
expect_identical(is.matrix(out.InfoMax$A[[3]]), TRUE)
expect_identical(is.matrix(out.ExtInfoMax$A[[3]]), TRUE)

expect_identical(dim(out.FastICA$A[[1]]), c(dim(X)[1], J[1]))
expect_identical(dim(out.InfoMax$A[[1]]), c(dim(X)[1], J[1]))
expect_identical(dim(out.ExtInfoMax$A[[1]]), c(dim(X)[1], J[1]))

expect_identical(dim(out.FastICA$A[[2]]), c(dim(X)[2], J[2]))
expect_identical(dim(out.InfoMax$A[[2]]), c(dim(X)[2], J[2]))
expect_identical(dim(out.ExtInfoMax$A[[2]]), c(dim(X)[2], J[2]))

expect_identical(dim(out.FastICA$A[[3]]), c(dim(X)[3], J[3]))
expect_identical(dim(out.InfoMax$A[[3]]), c(dim(X)[3], J[3]))
expect_identical(dim(out.ExtInfoMax$A[[3]]), c(dim(X)[3], J[3]))

### Test 0-4: S
expect_identical(is.array(out.FastICA$S), TRUE)
expect_identical(is.array(out.InfoMax$S), TRUE)
expect_identical(is.array(out.ExtInfoMax$S), TRUE)

expect_identical(dim(out.FastICA$S), J)
expect_identical(dim(out.InfoMax$S), J)
expect_identical(dim(out.ExtInfoMax$S), J)

### Test 0-5: J
expect_identical(out.FastICA$J, J)
expect_identical(out.InfoMax$J, J)
expect_identical(out.ExtInfoMax$J, J)

### Test 0-6: algorithm
expect_identical(out.FastICA$algorithm, "FastICA")
expect_identical(out.InfoMax$algorithm, "InfoMax")
expect_identical(out.ExtInfoMax$algorithm, "ExtInfoMax")

### Test 0-7: RecError
expect_identical(is.vector(out.FastICA$RecError), TRUE)
expect_identical(is.vector(out.InfoMax$RecError), TRUE)
expect_identical(is.vector(out.ExtInfoMax$RecError), TRUE)

## Test Error
### Test E-1: X
expect_error(MultiICA(as.array(X), J=J))
### Test E-2: J
expect_error(MultiICAICA(X, J="5"))
expect_error(MultiICAICA(X, J=c(2,4)))
expect_error(MultiICAICA(X, J=10^10)
### Test E-2: modes
expect_error(MultiICA(as.array(X), modes=1:4))
expect_error(MultiICA(as.array(X), modes=4))
expect_error(MultiICA(as.array(X), modes=-1:2))

## Session Information
sessionInfo()
