# Test code of MICA Function
## Dependent Packages Installation

## MICA Function (Paste your MICA code here)
MICA <- function(X, Y, J, eta=1E-10,
	num.iter=30, thr=1E-10, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkMICA(X, Y, J, eta, num.iter, thr, verbose)
    ######################################
    # Initialization (e.g. CCA)
    ######################################
    int <- .initMICA(X, Y, J)
    U <- int$U
    V <- int$V
    A <- int$A
    B <- int$B
    RecError <- int$RecError
    RelChange <- int$RelChange
    ######################################
    # Iteration
    ######################################
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
    	A <- A - eta * .gradA(U, V, A, B, J) # SGD
    	B <- B - eta * .gradB(U, V, A, B, J) # SGD
    	U <- X %*% ginv(A)
    	V <- Y %*% ginv(B)
        # After Update
        X_bar <- U %*% A
        Y_bar <- V %*% B
        iter <- iter + 1
        RecError[iter] <- .recError(X, X_bar) + .recError(Y, Y_bar)
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
    list(U=U, V=V, A=A, B=B, J=J, eta=eta, num.iter=num.iter,
    thr=thr, verbose=verbose, RecError=RecError, RelChange=RelChange)
}

.checkICA <- function(X, Y, J, eta, num.iter, thr, verbose){
	stopifnot(is.matrix(X))
	stopifnot(is.matrix(Y))
	stopifnot(nrow(X) == ncol(Y))
	stopifnot(is.numeric(J))
	stopifnot(length(J) == 1)
	stopifnot(is.numeric(eta))
	stopifnot(length(eta) == 1)
	stopifnot(min(ncol(X), ncol(Y)) >= J)
	stopifnot(is.numeric(num.iter))
	stopifnot(num.iter > 0)
	stopifnot(is.numeric(thr))
	stopifnot(is.logical(verbose))
}

.initICA <- function(X, Y, J){
	# CCA
	res.cca <- .CCA(X, Y, J)
	U <- res.cca$Xscore
	V <- res.cca$Yscore
	A <- res.cca$XLoading
	B <- res.cca$YLoading
	X_bar <- X %*% ginv(A)
	Y_bar <- Y %*% ginv(B)
	# Reconstruction Error
	RecError <- .recError(X, X_bar) + .recError(Y, Y_bar)
	# Relative Change
	RelChange <- thr * 10
	list(U=U, V=V, A=A, B=B,
		RecError=RecError, RelChange=RelChange)
}

.CCA <- function(X, Y, J){
   ...
}

## Simulation Datasets
X <- array(runif(10*15), dim=c(10,15))
Y <- array(runif(10*20), dim=c(10,20))

## Perform ICA against Simulation Datasets
J <- 5
out.MICA <- MICA(X, Y, J=J)

## Test Input object / type
### Test I-1: Object Names
expect_identical(names(formals(ICA)), c("X", "Y", "J", "eta", "num.iter", "thr", "verbose"))
### Test I-2: X
expect_identical(as.character(formals(ICA)$X), "")
### Test I-3: Y
expect_identical(as.character(formals(ICA)$Y), "")
### Test I-4: J
expect_identical(as.character(formals(ICA)$J), "")
### Test I-5: eta
expect_identical(as.character(formals(ICA)$eta), 1E-10)
### Test I-6: num.iter
expect_identical(as.character(formals(ICA)$num.iter), 30)
### Test I-7: thr
expect_identical(as.character(formals(ICA)$thr), 1E-10)
### Test I-8: verbose
expect_identical(as.character(formals(ICA)$verbose), FALSE)

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(out.MICA), TRUE)
### Test O-2: Object Names
expect_identical(names(out.MICA), c("U", "V", "A", "B", "J", "eta",
	"num.iter", "thr", "verbose", "RecError", "RelChange"))
### Test 0-3: U
expect_identical(is.matrix(out.MICA$U), TRUE)
expect_identical(dim(out.MICA$U), c(nrow(X), J))
### Test 0-4: V
expect_identical(is.matrix(out.MICA$V), TRUE)
expect_identical(dim(out.MICA$V), c(nrow(Y), J))
### Test 0-5: A
expect_identical(is.matrix(out.MICA$A), TRUE)
expect_identical(dim(out.MICA$A), c(ncol(X), J))
### Test 0-6: B
expect_identical(is.matrix(out.MICA$V), TRUE)
expect_identical(dim(out.MICA$V), c(ncol(Y), J))
### Test 0-7: J
expect_identical(out.MICA$J, J)
### Test 0-8: eta
expect_identical(out.MICA$eta, eta)
### Test 0-9: num.iter
expect_identical(out.MICA$num.iter, formals(ICA)$num.iter)
### Test 0-10: thr
expect_identical(out.MICA$thr, formals(ICA)$thr)
### Test 0-11: verbose
expect_identical(out.MICA$verbose, formals(ICA)$verbose)
### Test 0-12: RecError
expect_identical(is.vector(out.MICA$RecError), TRUE)
### Test 0-13: RelChange
expect_identical(is.vector(out.MICA$RelChange), TRUE)

## Test Error
### Test E-1: X
expect_error(ICA(X, J=J))
expect_error(ICA(Y, J=J))
expect_error(ICA(as.data.frame(X), Y, J=J))
expect_error(ICA(X, as.data.frame(Y), J=J))
### Test E-2: J
expect_error(ICA(X, Y, J="5"))
expect_error(ICA(X, Y, J=c(2,4)))
expect_error(ICA(X, Y, J=10^10)
### Test E-3: eta
expect_error(ICA(X, Y, eta="5"))
expect_error(ICA(X, Y, eta=c(2,4)))
expect_error(ICA(X, Y, eta=10^10)
### Test E-4: num.iter
expect_error(ICA(X, J=J, num.iter="100"))
expect_error(ICA(X, J=J, num.iter=-1))
### Test E-5: thr
expect_error(ICA(X, J=J, thr="0.1"))
### Test E-6: verbose
expect_error(ICA(X, J=J, verbose="verbose"))

## Test Decrease of Error
### Test D-1: RecError
.sampleRank <- function(x){
	rank(c(x[2], median(x), rev(x)[1]))
}
expect_identical(.sampleRank(out.MICA$RecError), 3:1)

### Test D-2: RelChange
expect_identical(.sampleRank(out.MICA$RelChange), 3:1)

## Session Information
sessionInfo()
