#' Generalized Eigenvalue Decomposition
#'
#' Solves the generalized eigenvalue problem A x = lambda B x.
#' This is a pure-R replacement for the geigen package.
#' @param A A square matrix
#' @param B A square matrix of the same dimension as A
#' @param symmetric Logical, whether A and B are symmetric.
#'   Currently only symmetric = TRUE is supported.
#' @param only.values Logical, if TRUE only eigenvalues are computed
#'   (Default: FALSE)
#' @return A list with components:
#'   \item{values}{The generalized eigenvalues, in ascending order.}
#'   \item{vectors}{The corresponding eigenvectors (columns),
#'     or NULL if only.values = TRUE.}
#'   \item{alpha}{NULL (included for compatibility with geigen package).}
#'   \item{beta}{NULL (included for compatibility with geigen package).}
#' @examples
#' A <- matrix(c(2, 1, 1, 3), 2, 2)
#' B <- matrix(c(1, 0.5, 0.5, 2), 2, 2)
#' geigen(A, B, symmetric = TRUE)
#' @export
geigen <- function(A, B, symmetric, only.values = FALSE) {
    if (!symmetric) {
        stop("Only symmetric = TRUE is supported")
    }
    # Solve via B^{-1} A and reorder to ascending eigenvalues
    eig <- eigen(solve(B) %*% A, symmetric = FALSE)
    ord <- order(Re(eig$values))
    list(
        values = Re(eig$values[ord]),
        vectors = if (only.values) NULL
                  else Re(eig$vectors[, ord, drop = FALSE]),
        alpha = NULL,
        beta = NULL
    )
}
