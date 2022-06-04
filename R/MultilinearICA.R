#' Multilinear independent component analysis
#' The input object is assumed to be a Tensor object defined by rTensor package.
#' In MultilinearICA, ICA function is performed in each mode of the tensor.
#' @param X An rTensor object
#' @param Js A vector to specify the rank in each mode (Default: c(3,3,3))
#' @param modes A vector to specify which modes are decomposed (Default: 1:3)
#' @param algorithm The algorithm to decompose the input tensor in each mode (Default: "FastICA")
#' @return A list containing the result of the decomposition
#' @examples
#' library("rTensor")
#' arrX <- array(runif(10*20*30), dim=c(10,20,30))
#' X <- as.tensor(arrX)
#' Js <- c(2,3,4)
#' out <- MultilinearICA(X, Js=Js)
#' @importFrom rTensor cs_unfold ttl ttm as.tensor
#' @importFrom methods is
#' @export
MultilinearICA <- function(X, Js=c(3,3,3), modes=1:3,
    algorithm=c("FastICA", "InfoMax", "ExtInfoMax")){
    ######################################
    # Argument Check
    ######################################
    .checkMultilinearICA(X, Js, modes)
    algorithm <- match.arg(algorithm)
    ######################################
    # Initialization
    ######################################
    int <- .initMultilinearICA(X, Js, modes)
    As <- int$As
    S <- int$S
    ######################################
    # ICA in each mode
    ######################################
    nrs <- dim(X)
    for(m in seq_along(nrs)){
        if(m %in% modes){
            position <- which(modes == m)
            Xm <- cs_unfold(X, m=m)@data
            As[[m]] <- ICA(X=Xm, J=Js[position], algorithm=algorithm)$A
        }
    }
    # Core tensor
    S <- .projection(X, As, modes=modes)
	# Output
    list(As=As, S=S, Js=Js, algorithm=algorithm)
}

.checkMultilinearICA <- function(X, Js, modes){
    stopifnot(is(X)[1] == "Tensor")
    stopifnot(is.vector(Js))
    stopifnot(length(Js) == length(modes))
    for(m in seq_along(modes)){
        stopifnot(dim(X)[modes[m]] >= Js[m])
    }
    stopifnot(all(.ndim(X) >= length(Js)))
    stopifnot(all(1 <= modes))
    stopifnot(all(.ndim(X) >= length(modes)))
}

.ndim <- function(X){
    length(dim(X))
}

.initMultilinearICA <- function(X, Js, modes){
	# As
	nrs <- dim(X)
    As <- list()
    length(As) <- length(nrs)
    for(m in seq_along(nrs)){
        if(m %in% modes){
            position <- which(modes == m)
            As[[m]] <- matrix(0, nrow=Js[position], ncol=nrs[m])
        }else{
            As[[m]] <- diag(nrs[m])
        }
    }
    # S
	S <- as.tensor(array(0, dim=Js))
	list(As=As, S=S)
}

.projection <- function(X, As, modes=modes){
    out <- X
    for(m in seq_len(.ndim(X))){
        if(m %in% modes){
            out <- ttm(out, ginv(t(As[[m]])), m=m)
        }else{
            out <- ttm(out, t(As[[m]]), m=m)
        }
    }
    out
}
