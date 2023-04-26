#' Independent Component Analysis (Classic Methods)
#'
#' The input data is assumed to be a matrix.
#' ICA decomposes the matrix and extract the components that
#' are statistically independent each other.
#' @param X A matrix
#' @param J Rank parameter to decompose
#' @param algorithm The decomposition algorithm (Default: "FastICA")
#' @param num.iter The number of iteration
#' @param thr The threshold to terminate the iteration (Default: 1E-10)
#' @param nonlinear_func The function used in FastICA (Default: "tanh")
#' @param learning_rate The learning rate used in InfoMax or ExtInfoMax
#' @param verbose Verbose option
#' @return A list containing the result of the decomposition
#' @examples
#' X <- matrix(runif(100*200), nrow=100, ncol=200)
#' J <- 5
#' out.FastICA <- ICA(X, J=J, algorithm="FastICA")
#' out.InfoMax <- ICA(X, J=J, algorithm="InfoMax")
#' out.ExtInfoMax <- ICA(X, J=J, algorithm="ExtInfoMax")
#' @importFrom MASS ginv
#' @importFrom stats cov
#' @export
ICA <- function(X, J,
    algorithm=c("FastICA", "InfoMax", "ExtInfoMax"),
    num.iter=100, thr=1E-10,
    nonlinear_func=c("tanh", "exp", "kurtosis"),
    learning_rate=1.0, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkICA(X, J, num.iter, thr, learning_rate, verbose)
    algorithm <- match.arg(algorithm)
    nonlinear_func <- match.arg(nonlinear_func)
    ######################################
    # Initialization (e.g. Whitening)
    ######################################
    int <- .initICA(X, J, thr=thr, algorithm=algorithm)
    X <- int$X
    W <- int$W
    WhiteningMatrix <- int$WhiteningMatrix
    S <- int$S
    WChange <- int$WChange
    update_func <- int$update_func
    ######################################
    # Iteration
    ######################################
    iter <- 1
    while ((WChange[iter] > thr) && (iter <= num.iter)){
        # Update W
        WNew <- update_func(X, W, S, J, nonlinear_func, learning_rate)
        # After Update
        S <- X %*% WNew
        A <- ginv(WhiteningMatrix %*% WNew)
        X_bar <- S %*% A
        iter <- iter + 1
        # The maximum value of the change in each component
        # of W is defined as the amount of change.
        WChange[iter] <- max(abs(WNew - W))
        W <- WNew
        # Verbose
        if (verbose){
            cat(paste0(iter, " / ", num.iter,
                "max(abs(WNew - W)) = ",
                WChange[iter], "\n"))
        }
    }
    # Output
    names(WChange) <- c("offset", seq_len(iter - 1))
    list(
        A=A, S=S, J=J, algorithm=algorithm, num.iter=num.iter,
        thr=thr, verbose=verbose, WChange=WChange)
}

.checkICA <- function(X, J, num.iter, thr, learning_rate, verbose){
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(J))
    stopifnot(length(J) == 1)
    stopifnot(min(dim(X)) >= J)
    stopifnot(is.numeric(num.iter))
    stopifnot(num.iter >= 0)
    stopifnot(is.numeric(thr))
    stopifnot(is.logical(verbose))
    stopifnot(is.numeric(learning_rate))
}

.initICA <- function(X, J, algorithm, thr = 0.1){
    WhiteningResult <- .whitening(X, J)
    X <- WhiteningResult$XWhiten
    J <- WhiteningResult$JNew
    WhiteningMatrix <- WhiteningResult$WhiteningMatrix
    # initialize recovering matrix
    W <- diag(J)
    # If the initial value of the restoration matrix is not a unit matrix,
    # the following calculation makes sense.
    S <- X %*% W
    X_bar <- S %*% ginv(W)
    # Relative Change
    WChange <- thr * 2
    # select update algorithm
    update_func <- .update[[algorithm]]
    list(X=X, W=W, S=S, WChange=WChange,
        WhiteningMatrix=WhiteningMatrix, update_func=update_func)
}

.whitening <- function(X, J){
    # apply centering
    XMeans <- colMeans(X)
    XCentered <- scale(X, scale=FALSE)
    # compute covariance matrix
    Sigma <- cov(XCentered)
    # compute eigenvectors
    EigenSigma <- eigen(Sigma)
    # check number of nonzero elements
    # TODO: includes machine epsilon
    NumNonzeroEigenvalues <- sum(EigenSigma$values > 0)
    if (NumNonzeroEigenvalues < J){
        J <- NumNonzeroEigenvalues
    }
    LambdaSqrtInv <- diag(1 / sqrt(EigenSigma$values[1:J]))
    WhiteningMatrix <- EigenSigma$vectors[, 1:J] %*% t(LambdaSqrtInv)
    XWhiten <- XCentered %*% WhiteningMatrix
    # Output
    return(list(XWhiten=XWhiten, JNew=J,
        WhiteningMatrix=WhiteningMatrix, XMeans=XMeans))
}

.orthgonalize <- function(W){
    WNewSvd <- svd(W)
    tcrossprod(WNewSvd$u, WNewSvd$v)
}

# learning_rate is ignored
.FastICA <- function(X, W, S, J, nonlinear_func, learning_rate){
    # Setting
    n <- nrow(X)
    J <- ncol(X)
    # In order to conform to the notation of the paper,
    # W and X is treated as transposed.
    W <- t(W)  # J * J
    X <- t(X)  # J * n
    # Select Non-linear Function
    tmp <- .nonlinear_transform_function_fastica(nonlinear_func)
    g <- tmp$g
    gPrime <- tmp$gPrime
    #########################################
    # Symmetric FastICA
    #########################################
    WX <- W %*% X                                # J * n
    gWX <- g(WX)                                 # J * n
    gPrimeWX <- gPrime(WX)                       # J * n
    EgPrimeWX <- apply(gPrimeWX, 1, FUN = mean)  # J
    EgPrimeWXW <- W %*% diag(EgPrimeWX)          # J * J
    XgWX <- gWX %*% t(X) / n                     # J * J
    WNew <- XgWX - EgPrimeWXW                    # J * J
    # Orthogonalization
    WNew <- .orthgonalize(WNew)
    # Update
    W <- WNew
    # Output
    return(t(W))
}

.nonlinear_transform_function_fastica <- function(nonlinear_func){
    if (nonlinear_func == "tanh"){
        alpha <- 1.0
        g <- function(x){ tanh(alpha* x) }
        gPrime <- function(x){alpha*(1 - tanh(alpha*x)^2)}
    } else if (nonlinear_func == "exp"){
        g <- function(x){x * exp(-x*2/2)}
        gPrime <- function(x){(1-x^2) * exp(-x^2/2)}
    } else if (nonlinear_func == "kurtosis"){
        g <- function(x){x^3}
        gPrime <- function(x){3*x^2}
    }
    list(g=g, gPrime=gPrime)
}

# nonlinear_func is ignored
.InfoMax <- function(X, W, S, J, nonlinear_func, learning_rate){
    ####################################################
    # Gradient method using inverse matrix
    ####################################################
    U <- X %*% W
    Y <- .nonlinear_transform_function_infomax(U)
    delta_W_without_rate <- .dydx_div_dy2dxdw(X, Y, W, nrow(X))
    delta_W <- learning_rate * delta_W_without_rate
    # if you use bias term, include these.
    # delta_W0_without_rate <- dydx_div_dy2dxdw0(X, Y, W)
    # delta_W0 <- learning_rate * delta_W0_without_rate
    # Update
    new_W <- W + delta_W
    # Orthgonalize
    new_W <- .orthgonalize(new_W)
    # Output
    return(new_W)
}

.nonlinear_transform_function_infomax <- function(u){1 / (1 + exp(-u))}

.dydx_div_dy2dxdw <- function(X, Y, W, num_observation){
    return(ginv(t(W)) + crossprod(X, (1-2*Y)) / num_observation)
}

# nonlinear_func is ignored
.ExtInfoMax <- function(X, W, S, J, nonlinear_func, learning_rate){
    ####################################################
    # Gradient method using Natural Gradient
    # It doesn't use inverse matrix.
    ####################################################
    U <- X %*% W
    # Automatically determines whether it is a sub-Gaussian or a super-Gaussian.
    K <- diag(sign(colMeans(.sech(U)^2)*colMeans(U^2) - colMeans(tanh(U)*U)))
    delta_W <- learning_rate * .delta_W_without_learning_rate(U, W, signs=K)
    new_W <- W + delta_W
    # Orthgonalize
    new_W <- .orthgonalize(new_W)
    # Output
    return(new_W)
}

# This is a nonlinear function for both sub-Gaussian and super-Gaussian.
.delta_W_without_learning_rate <- function(U, W, signs){
    return((1-signs %*% crossprod(tanh(U), U) - crossprod(U, U)) %*% W)
}

.sech <- function(x){ 1 / cosh(x) }

.update <- list(
    FastICA=.FastICA,
    InfoMax=.InfoMax,
    ExtInfoMax=.ExtInfoMax)
