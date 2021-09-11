# Test code of ICA Function
## Dependent Packages Installation
if (!requireNamespace("testthat", quietly = TRUE)){
    install.packages("testthat")
}
library("testthat")

## ICA Function (Paste your ICA code here)

ICA <- function(X, J, algorithm=c("FastICA", "InfoMax", "ExtInfoMax"), num.iter=30, thr=1E-10, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkICA(X, J, num.iter, thr, verbose)
    algorithm <- match.arg(algorithm)
    ######################################
    # Initialization (e.g. Whiteniing)
    ######################################
    int <- .initICA(X, J)
    X <- int$X
    A <- int$A
    S <- int$S
    RecError <- int$RecError
    RelChange <- int$RelChange
    ######################################
    # Iteration
    ######################################
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
        if(algorithm="FastICA"){
        	A <- .FastICA(X, A, S, J)
        }
        if(algorithm="InfoMax"){
        	A <- .InfoMax(X, A, S, J)
        }
        if(algorithm="ExtInfoMax"){
        	A <- .ExtInfoMax(X, A, S, J)
        }
        # After Update
        S <- X %*% ginv(A)
        X_bar <- A %*% S
        iter <- iter + 1
        RecError[iter] <- .recError(X, X_bar)
        RelChange[iter] <- .relChange(iter, RecError)
        # Verbose
        if(verbose){
             cat(paste0(iter, " / ", num.iter,
                " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
    }
    # Output
	names(RecError) <- c("offset", 1:(iter - 1))
    names(RelChange) <- c("offset", 1:(iter - 1))
	# Output
    list(A=A, S=S, J=J, algorithm=algorithm, num.iter=num.iter,
    thr=thr, verbose=verbose, RecError=RecError, RelChange=RelChange)
}

.FastICA <- function(X, J){
	...
}

.InfoMax <- function(X, J){
	...
}

.ExtInfoMax <- function(X, J){
	...
}

.checkICA <- function(X, J, num.iter, thr, verbose){
	stopifnot(is.matrix(X))
	stopifnot(is.numeric(J))
	stopifnot(length(J) == 1)
	stopifnot(min(dim(X)) >= J)
	stopifnot(is.numeric(num.iter))
	stopifnot(num.iter > 0)
	stopifnot(is.numeric(thr))
	stopifnot(is.logical(verbose))
}

.initICA <- function(X, J){
	X <- .whitening(X)
	# A/S
	nr <- nrow(X)
	A <- matrix(runif(nr*J), nrow=nr, ncol=J)
	S <- X %*% ginv(A)
	X_bar <- A %*% S
	# Reconstruction Error
	RecError <- .recError(X, X_bar)
	# Relative Change
	RelChange <- thr * 10
	list(X=X, A=A, S=S, RecError=RecError, RelChange=RelChange)
}

.whitening <- function(X){
	...
}

.recError <- function(X, X_bar){
	...
}

.relChange <- function(iter, RecError){
	...
}


## Simulation Datasets
X <- matrix(runif(100*200), nrow=100, ncol=100)

## Perform ICA against Simulation Datasets
J <- 5
out.FastICA <- ICA(X, J=J, algorithm="FastICA")
out.InfoMax <- ICA(X, J=J, algorithm="InfoMax")
out.ExtInfoMax <- ICA(X, J=J, algorithm="ExtInfoMax")

## Test Input object / type
### Test I-1: Object Names
expect_identical(names(formals(ICA)), c("X", "J", "algorithm", "num.iter", "thr", "verbose"))
### Test I-2: X
expect_identical(as.character(formals(ICA)$X), "")
### Test I-3: J
expect_identical(as.character(formals(ICA)$J), "")
### Test I-4: algorithm
expect_identical(as.character(formals(ICA)$algorithm), c("FastICA", "InfoMax", "ExtInfoMax"))
### Test I-5: num.iter
expect_identical(as.character(formals(ICA)$num.iter), 30)
### Test I-6: thr
expect_identical(as.character(formals(ICA)$thr), 1E-10)
### Test I-7: verbose
expect_identical(as.character(formals(ICA)$verbose), FALSE)

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(out.FastICA), TRUE)
expect_identical(is.list(out.InfoMax), TRUE)
expect_identical(is.list(out.ExtInfoMax), TRUE)

### Test O-2: Object Names
expect_identical(names(out.FastICA), c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.InfoMax), c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.ExtInfoMax), c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))

### Test 0-3: A
expect_identical(is.matrix(out.FastICA$A), TRUE)
expect_identical(is.matrix(out.InfoMax$A), TRUE)
expect_identical(is.matrix(out.ExtInfoMax$A), TRUE)

expect_identical(dim(out.FastICA$A), c(nrow(X), J))
expect_identical(dim(out.InfoMax$A), c(nrow(X), J))
expect_identical(dim(out.ExtInfoMax$A), c(nrow(X), J))

### Test 0-4: S
expect_identical(is.matrix(out.FastICA$S), TRUE)
expect_identical(is.matrix(out.InfoMax$S), TRUE)
expect_identical(is.matrix(out.ExtInfoMax$S), TRUE)

expect_identical(dim(out.FastICA$S), c(J, ncol(X)))
expect_identical(dim(out.InfoMax$S), c(J, ncol(X)))
expect_identical(dim(out.ExtInfoMax$S), c(J, ncol(X)))

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
expect_identical(is.vector(out.FastICA$RecError), TRUE)
expect_identical(is.vector(out.InfoMax$RecError), TRUE)
expect_identical(is.vector(out.ExtInfoMax$RecError), TRUE)

### Test 0-11: RelChange
expect_identical(is.vector(out.FastICA$RelChange), TRUE)
expect_identical(is.vector(out.InfoMax$RelChange), TRUE)
expect_identical(is.vector(out.ExtInfoMax$RelChange), TRUE)

## Test Error
### Test E-1: X
expect_error(ICA(as.data.frame(X), J=J))
### Test E-2: J
expect_error(ICA(X, J="5"))
expect_error(ICA(X, J=c(2,4)))
expect_error(ICA(X, J=10^10)
### Test E-3: num.iter
expect_error(ICA(X, J=J, num.iter="100"))
expect_error(ICA(X, J=J, num.iter=-1))
### Test E-4: thr
expect_error(ICA(X, J=J, thr="0.1"))
### Test E-5: verbose
expect_error(ICA(X, J=J, verbose="verbose"))

## Test Decrease of Error
### Test D-1: RecError
.sampleRank <- function(x){
	rank(c(x[2], median(x), rev(x)[1]))
}
expect_identical(.sampleRank(out.FastICA$RecError), 3:1)
expect_identical(.sampleRank(out.InfoMax$RecError), 3:1)
expect_identical(.sampleRank(out.ExtInfoMax$RecError), 3:1)

### Test D-2: RelChange
expect_identical(.sampleRank(out.FastICA$RelChange), 3:1)
expect_identical(.sampleRank(out.InfoMax$RelChange), 3:1)
expect_identical(.sampleRank(out.ExtInfoMax$RelChange), 3:1)

## Session Information
sessionInfo()


