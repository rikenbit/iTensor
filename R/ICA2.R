#' Independent Component Analysis (Modern Methods)
#'
#' The input data is assumed to be a matrix.
#' ICA decomposes the matrix and extract the components that
#' are statistically independent each other.
#' @param X A matrix
#' @param J Rank parameter to decompose
#' @param algorithm The decomposition algorithm (Default: "JADE")
#' @param num.iter The number of iteration
#' @param thr The threshold to terminate the iteration (Default: 1E-10)
#' @param r_list List of r-th order cumulants used in SIMBEC (Default: NULL)
#' @param omega_for_each_r Weight vector of r_list used in SIMBEC (Default: NULL)
#' @param a_r_for_each_r Parameter vector to specify the shape of partial activation function in SIMBEC (Default: NULL)
#' @param tau_list List of lags to consider the auto-correlation used in AMUSE and SOBI (Default: NULL)
#' @param num_bins Number of bins for histgram in ProDenICA (Default: NULL)
#' @param alpha Learning rate used for gradient descent in RICA (Default: NULL)
#' @param num_epoch Number of epoch used for gradient descent in RICA (Default: NULL)
#' @param verbose Verbose option
#' @return A list containing the result of the decomposition
#' @examples
#' ICA2
#' @importFrom MASS ginv
#' @importFrom graphics hist
#' @importFrom stats cov runif dnorm predict var sd
#' @importFrom mgcv gam
#' @importFrom geigen geigen
#' @importFrom einsum einsum
#' @export
ICA2 <- function(
    X, J,
    algorithm = c(
        "JADE", "AuxICA1", "AuxICA2", "IPCA", "SIMBEC", "AMUSE",
        "SOBI", "FOBI", "ProDenICA", "RICA"),
    num.iter = NULL, thr = 1E-10,
    r_list = NULL,
    omega_for_each_r = NULL,
    a_r_for_each_r = NULL,
    tau_list = NULL,
    num_bins = NULL,
    alpha = NULL,
    num_epoch = NULL,
    verbose = FALSE){
    ######################################
    # Argument Check
    ######################################
    algorithm <- match.arg(algorithm)
    .checkICA2(X, J, num.iter, thr, r_list, omega_for_each_r,
        a_r_for_each_r, tau_list, num_bins,
        alpha, num_epoch, verbose)
    ######################################
    # Iterative or Non-Iterative
    ######################################
    args <- .initICA2_1(X, J, algorithm, num.iter, thr, r_list,
        omega_for_each_r, a_r_for_each_r, tau_list,
        num_bins, alpha, num_epoch, verbose)
    if(algorithm %in% c("JADE", "IPCA", "SIMBEC", "AMUSE", "SOBI", "FOBI")){
        return(do.call(.non_iterative_ICA, args))
    }else if(algorithm %in% c("AuxICA1", "AuxICA2", "ProDenICA", "RICA")){
        return(do.call(.iterative_ICA, args))
    }
}

.non_iterative_ICA <- function(
        X, J,
        algorithm = c("JADE", "IPCA", "SIMBEC", "AMUSE", "SOBI", "FOBI"),
        ...){
    # Argument Check
    algorithm <- match.arg(algorithm)
    args <- list(...)
    verbose <- args[["verbose"]]
    # Performs preprocessing according to the algorithm.
    args[["X"]] <- .preprocess_non_iterative_ICA_input(X = X, J = J, algorithm = algorithm)
    args[["J"]] <- J
    # Runs the algorithm.
    ica_result <- do.call(.ica_function_list[[algorithm]], args)
    # Adds a common attribute to the output that does not depend on the algorithm.
    ica_result[["algorithm"]] <- algorithm
    ica_result[["verbose"]] <- verbose
    ica_result
}

.preprocess_non_iterative_ICA_input <- function(X, J, algorithm){
    # IPCA and JADE differ from the other algorithms in how X is initialized.
    if(algorithm == "IPCA"){
        # IPCA does not apply preprocessing to X.
    }else if(algorithm == "JADE"){
        # For JADE, no scaling is applied to X.
        if(nrow(X) < ncol(X)){
            stop("The number of samples must be larger than the number of dimensions.")
        }
        # Initialize X.
        X_initialized <- .init_X(X, J, scale_X = FALSE)
        X <- X_initialized$X
    }else{
        # For non-IPCA, nrow (X)> = ncol (X).
        if(nrow(X) < ncol(X)){
            stop("The number of samples must be larger than the number of dimensions.")
        }
        # Initialize X.
        X_initialized <- .init_X(X, J)
        X <- X_initialized$X
    }
    X
}

.iterative_ICA <- function(
        X, J,
        algorithm = c("AuxICA1", "AuxICA2", "ProDenICA", "RICA"),
        num.iter = 30, thr = 1E-10, verbose = FALSE, ...){
    ######################################
    # Argument Check
    ######################################
    algorithm <- match.arg(algorithm)
    args <- list(...)
    args[["X"]] <- X
    args[["J"]] <- J
    args[["algorithm"]] <- algorithm
    args[["num.iter"]] <- num.iter
    args[["thr"]] <- thr
    args[["verbose"]] <- verbose
    ######################################
    # Initialization (e.g. Whiteniing)
    ######################################
    if(algorithm == "RICA"){
        allow_overcomplete <- TRUE
    }else{
        allow_overcomplete <- FALSE
    }
    int <- .initICA2_2(X, J, thr, allow_overcomplete)
    X <- int$X
    A <- int$A
    W <- int$W
    S <- int$S
    RecError <- int$RecError
    RelChange <- int$RelChange
    ######################################
    # Iteration
    ######################################
    iter <- 1
    # TODO: Replace the following line with the next line.
    # while ((RelChange[iter] > thr) && (iter <= num.iter)){
    while (iter <= num.iter){
        # Runs the algorithm.
        args[["X"]] <- X
        args[["W"]] <- W
        args[["J"]] <- J
        ica_result <- do.call(.ica_function_list[[algorithm]], args)
        W <- ica_result$W
        # After Update
        A <- ginv(W)
        S <- X %*% W
        X_bar <- S %*% A
        iter <- iter + 1
        RecError[iter] <- .recError(X, X_bar)
        RelChange[iter] <- .relChange(iter, RecError)
        # Verbose
        if(verbose){
            cat(paste0(
                iter, " / ", num.iter,
                " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"
            ))
        }
    }
    # Output
    names(RecError) <- c("offset", 1:(iter - 1))
    names(RelChange) <- c("offset", 1:(iter - 1))
    list(
        A = A, S = S, J = J, W = W, X = X, algorithm = algorithm, num.iter = num.iter,
        thr = thr, verbose = verbose, RecError = RecError, RelChange = RelChange
    )
}

.JADE <- function(X, W = NULL, J = NULL, num.iter = 500, verbose = FALSE, ...){
    num_dimensions <- ncol(X)
    # Estimates the cumulant tensor.
    # NOTE: It is a good idea to estimate the covariance matrix
    #       in advance if you want to speed up the following process.
    Q <- array(0.0, dim = rep(num_dimensions, 4))
    for(i in 1:num_dimensions){
        for(j in 1:num_dimensions){
            for(k in 1:num_dimensions){
                for(l in 1:num_dimensions){
                    Q[i, j, k, l] <- mean(X[, i] * X[, j] * X[, k] * X[, l]) -
                        mean(X[, i] * X[, j]) * mean(X[, k] * X[, l]) -
                        mean(X[, i] * X[, k]) * mean(X[, j] * X[, l]) -
                        mean(X[, i] * X[, l]) * mean(X[, j] * X[, k])
                }
            }
        }
    }
    Q_unfolded <- matrix(as.numeric(Q), nrow = num_dimensions^2, ncol = num_dimensions^2)
    eigen_Q <- eigen(Q_unfolded)
    # jointly diagonalise the set \hat{N}^e = \{ \hat{lambda}_r \hat{M}_r | 1 \le r \le n}
    n_eigen_matrices_diagonalized <- num_dimensions
    num_components <- num_dimensions
    # Create a list of λ M in array format.
    # The last dimension is the number of the eigenmatrix to be diagonalized.
    lambda_M_list <- array(
        0.0,
        dim = c(num_dimensions, num_dimensions, n_eigen_matrices_diagonalized))
    order_of_abs_eigenvalues <- order(abs(eigen_Q$values), decreasing = TRUE)
    for(i in 1:n_eigen_matrices_diagonalized){
        lambda_M_list[, , i] <- eigen_Q$values[order_of_abs_eigenvalues[i]] *
            matrix(eigen_Q$vectors[, order_of_abs_eigenvalues[i]], num_dimensions, num_dimensions)
    }
    V <- diag(num_components)
    for(i_loop in 1:num.iter){
        for(p in 1:num_components){
            for(q in 1:num_components){
                if(p == q) next
                # Use Cyclic-by-Row Algorithm
                # Reference: Gene H. Golub, Charles F. Van Loan Matrix Computations p.480
                g_r <- matrix(NA, nrow = num_components, ncol = 3)  # An array of num_component * 3.
                for(r in 1:num_components){
                    g_r[r, ] <- c(
                        lambda_M_list[p, p, r] - lambda_M_list[q, q, r],
                        lambda_M_list[p, q, r],
                        lambda_M_list[q, p, r]
                    )
                }
                B <- matrix(c(1, 0, 0, 0, 1, 1, 0, -1i, 1i), nrow = 3, ncol = 3, byrow = TRUE)
                ReGHG <- Re(B %*% t(g_r) %*% g_r %*% t(B))
                ReGHG_eigen <- eigen(ReGHG)
                angles <- ReGHG_eigen$vectors[, 1]
                if(angles[1] < 0) angles <- -angles
                c_ <- sqrt(0.5 + angles[1] / 2)
                s_ <- 0.5 * (angles[2] - 1i * angles[3]) / c_
                # Creates a Givens rotation matrix.
                GivensRot <- .R_ijcs(p, q, c_, s_, num_components)
                V <- V %*% GivensRot
                for(r in 1:num_components){
                    lambda_M_list[, , r] <- t(GivensRot) %*% lambda_M_list[, , r] %*% GivensRot
                }
            }
        }
    }
    A <- V
    un_mixing_matrix <- solve(A)
    S <- X %*% V
    X_bar <- S %*% A
    # Output
    list(
        A = A,
        W = un_mixing_matrix,
        S = S,
        X = X_bar,
        J = J,
        num.iter = 1,
        thr = NULL,
        RecError = .recError(X, X_bar),  # Returns only the last reconstruction error.
        RelChange = NULL
    )
}

# Generate Givens rotation matrix.
.R_ijcs <- function(i, j, c, s, n){
    R <- diag(n)
    R[i, i] <- c
    R[i, j] <- -Conj(s)
    R[j, i] <- s
    R[j, j] <- c
    R
}

.AuxICA1 <- function(X, W = NULL, J = NULL, verbose = FALSE, ...){
    stopifnot(!is.null(W))
    # Within this function, we treat the rows as sensors and the columns as samples,
    # according to the notation in the paper.
    # X_t: J * n
    X_t <- t(X)
    n <- nrow(X)
    # W: J * J
    # Each row of W is a row vector of w_1, w_2, ... w_J.
    W <- t(W)
    G_prime <- function(x){
        # Let G (x) = log(cosh(abs(x))), where G'(x) = tanh(x).
        return(tanh(x))
    }
    for(k in 1:J){
        # Auxiliary variable updates
        w_k <- W[k, ]
        # n * J * J
        V_k_i_list <- array(0.0, dim = c(n, J, J))
        for(i in 1:n){
            x <- X_t[, i]
            # r_k is inner product of w_k and x (scalar)
            r_k <- abs(as.numeric(w_k %*% x))
            # J * J
            V_k_i <- G_prime(r_k) / r_k * outer(x, x)
            V_k_i_list[i, , ] <- V_k_i
        }
        # J * J
        V_k <- apply(V_k_i_list, c(2, 3), mean)
        # (J - 1) * J
        W_k_minus_1 <- t(W[-k, ])
        # J * (J - 1)
        P <- V_k %*% W_k_minus_1
        # Demixing matrix updates
        w_k <- w_k - P %*% solve(t(P) %*% P) %*% t(P) %*% w_k
        w_k <- w_k / sqrt(as.numeric(t(w_k) %*% V_k %*% w_k))
        W[k, ] <- w_k
    }
    # Output
    list(
        A = NULL,
        W = t(W),
        S = NULL,
        X = X,
        J = J,
        num.iter = NULL,
        thr = NULL,
        RecError = NULL,
        RelChange = NULL
    )
}

.AuxICA2 <- function(X, W = NULL, J = NULL, verbose = FALSE, ...){
    # Within this function, we treat the rows as sensors and the columns as samples,
    # according to the notation in the paper.
    # X_t: J * num_sample
    X_t <- t(X)
    num_sample <- nrow(X)
    # W: J * J
    # Each row of W is a row vector of w_1, w_2, ... w_J.
    W <- t(W)
    G_prime <- function(x){
        # Let G(x) = log(cosh(abs(x))), where G'(x) = tanh(abs(x)).
        return(tanh(x))
    }
    for(m in 1:(J - 1)){
        for(n in (m + 1):J){
            # Auxiliary variable updates
            # J * 1
            w_m <- W[m, ]
            # J * 1
            w_n <- W[n, ]
            # n * 2 * 2
            U_m_i_list <- array(0.0, dim = c(num_sample, 2, 2))
            # n * 2 * 2
            U_n_i_list <- array(0.0, dim = c(num_sample, 2, 2))
            for(i in 1:num_sample){
                x <- X_t[, i]
                w_m_h_x <- as.numeric(w_m %*% x)
                w_n_h_x <- as.numeric(w_n %*% x)
                u <- c(w_m_h_x, w_n_h_x)
                r_m <- abs(w_m_h_x)
                r_n <- abs(w_n_h_x)
                # 2 * 2
                U_m_i <- G_prime(r_m) / r_m * outer(u, u)
                U_n_i <- G_prime(r_n) / r_n * outer(u, u)

                U_m_i_list[i, , ] <- U_m_i
                U_n_i_list[i, , ] <- U_n_i
            }
            # J * J
            U_m <- apply(U_m_i_list, c(2, 3), mean)
            U_n <- apply(U_n_i_list, c(2, 3), mean)
            # Demixing matrix updates
            eig <- geigen(U_m, U_n, symmetric = TRUE)
            h_m <- eig$vectors[, 1]
            h_n <- eig$vectors[, 2]
            h_m <- h_m / sqrt(as.numeric(t(h_m) %*% U_m %*% h_m))
            h_n <- h_n / sqrt(as.numeric(t(h_n) %*% U_n %*% h_n))
            w_m_w_n <- cbind(w_m, w_n) %*% cbind(h_m, h_n)
            w_m <- w_m_w_n[, 1]
            w_n <- w_m_w_n[, 2]
            W[m, ] <- w_m
            W[n, ] <- w_n
        }
    }
    # Output
    list(
        A = NULL,
        W = t(W),
        S = NULL,
        X = X,
        J = J,
        num.iter = NULL,
        thr = NULL,
        RecError = NULL,
        RelChange = NULL
    )
}

.IPCA <- function(X, W = NULL, J = NULL, verbose = FALSE, num.iter = 200, ...){
    # IPCA assumes num_columns > num_rows.
    num_components <- J
    X_centered <- scale(X, scale = FALSE) # size: n * p
    X_svd <- svd(X_centered)
    V <- X_svd$v # p * n
    U <- X_svd$u
    V_scaled <- scale(V) # p * n
    V_t <- t(V_scaled)[1:num_components, ] # n * p -> num_components * p
    transposed_Vt_whitening_result <- .whitening2(t(V_t))
    # initial_W: num_components * num_components
    initial_W <- diag(num_components)
    # W_: num_components * num_components
    W_ <- initial_W
    for(i_iteration in 1:num.iter){
        W_ <- .FastICA2(
            transposed_Vt_whitening_result$XWhiten,
            W = W_, S = NULL, J = num_components)
    }
    W_ <- solve(W_)
    # S: num_components * p
    S <- W_ %*% V_t
    kurtosis <- function(x){
        X_scaled <- scale(x)
        mean(X_scaled^4) - 3
    }
    S_kurtosis <- apply(S, MARGIN = 1, FUN = kurtosis)
    order_S_kurtosis <- order(S_kurtosis, decreasing = TRUE)
    S <- S[order_S_kurtosis, ]
    U_tilde <- X_centered %*% t(S) # n * num_components
    X_bar <- U_tilde %*% ginv(t(S))
    # Output
    list(
        A = ginv(t(S)),
        W = t(S),
        S = U_tilde,
        X = X_bar,
        J = J,
        num.iter = 1,
        thr = NULL,
        RecError = .recError(X_centered, X_bar),
        RelChange = NULL
    )
}

.SIMBEC <- function(X, W = NULL, J = NULL, verbose = FALSE, num.iter = 1000,
                    r_list = c(4), omega_for_each_r = c(1.0),
                    a_r_for_each_r = c(1.0), ...){
    # Within this function, we treat the rows as sensors and the columns as samples,
    # according to the notation in the paper.
    X_centered <- t(X)
    n_r <- length(r_list)
    E <- J
    N <- nrow(X_centered)
    # U: E * N
    U <- matrix(0.0, nrow = E, ncol = N)
    for(i_component in 1:E){
        U[i_component, i_component] <- 1.0
    }
    for(i_iter in 1:num.iter){
        Y <- U %*% X_centered
        Y_centered <- scale(Y, scale = FALSE)
        # |C^4_{y_i}| ... forth_order_auto_cumulants
        # Reference: https://mathworld.wolfram.com/Cumulant.html
        forth_order_auto_cumulants_of_y <- rowMeans(Y_centered^4) - 3 * apply(Y_centered, 1, var)^2
        # \mu
        mu <- 2 / (3 * max(abs(forth_order_auto_cumulants_of_y)))
        added_term_without_mu <- matrix(0.0, nrow = E, ncol = N)
        for(i_r in 1:n_r){
            # Within this loop, we compute:
            # \omega _{r}( D_{2}^{r}C_{y,x}^{r-1,1}-C_{y,y}^{1,r-1}D_{z}^{r}U^{\left( k\right) }
            r <- r_list[i_r]
            omega_r <- omega_for_each_r[i_r]
            # |C^r_{y_i}| ... r-th_order_auto_cumulants
            # refer: https://mathworld.wolfram.com/Cumulant.html
            if(r == 4){
                r_th_order_auto_cumulants_of_y <- rowMeans(Y_centered^4) - 3 * apply(Y_centered, 1, var)^2
            }else if(r == 3){
                r_th_order_auto_cumulants_of_y <- rowMeans(Y_centered^3)
            }else{
                stop("r must be 3 or 4.")
            }
            # \alpha_r
            alpha_r <- a_r_for_each_r[i_r]
            # \boldmath{D}^r_y
            # Size is E * E
            # It is computed by the elementwise product.
            D_r_y <- diag(sign(r_th_order_auto_cumulants_of_y) * abs(r_th_order_auto_cumulants_of_y)^(alpha_r - 1))
            # \boldmath{C}_{y,x}^{r-1,1}
            # Size is E * N
            # E(xi^3 * xj) - E(xixj)E(xixi) - E(xixj)E(xixi) - E(xixj)E(xixi)
            # = E(xi^3 xj) - 3E(xixj)E(xixi)
            r_th_order_cross_cumulants_of_x_y <- matrix(0.0, nrow = E, ncol = N)
            for(i_row in 1:E){
                for(i_col in 1:N){
                    r_th_order_cross_cumulants_of_x_y[i_row, i_col] <- mean(
                        Y_centered[i_row, ]^(r - 1) * X_centered[i_col, ]
                    ) -
                        3 * mean(Y_centered[i_row, ] * X_centered[i_col, ]) *
                            mean(Y_centered[i_row, ] * Y_centered[i_row, ])
                }
            }
            # \boldmath{C}_{y,y}^{1,r-1}
            # Size is E * E
            # E(xi * xj^3) - E(xixj)E(xjxj) - E(xixj)E(xjxj) - E(xixj)E(xjxj)
            # = E(xi xj^3) - 3E(xixj)E(xjxj)
            r_th_order_cross_cumulants_of_y_y <- matrix(0.0, nrow = E, ncol = E)
            for(i_row in 1:E){
                for(i_col in 1:E){
                    r_th_order_cross_cumulants_of_y_y[i_row, i_col] <- mean(
                        Y_centered[i_row, ] * Y_centered[i_col, ]^(r - 1)
                    ) -
                        3 * mean(Y_centered[i_row, ] * Y_centered[i_col, ]) *
                            mean(Y_centered[i_col, ] * Y_centered[i_col, ])
                }
            }
            added_term_without_mu <- added_term_without_mu +
                omega_r * (
                    D_r_y %*% r_th_order_cross_cumulants_of_x_y -
                        r_th_order_cross_cumulants_of_y_y %*% D_r_y %*% U)
        }
        U <- U + mu * added_term_without_mu
        U <- .orthgonalize(U)
    }
    S <- t(U %*% X_centered)
    A <- ginv(t(U))
    X_bar <- S %*% A
    # Output
    list(
        A = A,
        W = t(U),
        S = S,
        X = X_bar,
        J = J,
        num.iter = 1,
        thr = NULL,
        RecError = .recError(t(X_centered), X_bar), # Returns only the estimate completion error.
        RelChange = NULL
    )
}

# Converts the given matrix to an orthonormal matrix
.orthgonalize <- function(W){
    W_svd <- svd(W)
    tcrossprod(W_svd$u, W_svd$v)
}

.AMUSE <- function(X, W = NULL, J = NULL, verbose = FALSE, tau_list = 1:50, ...){
    # This function treats rows as sensors and columns as samples (time).
    X_centered <- t(X)
    Y <- t(X_centered)
    max_eigen_values <- c()
    for(tau in tau_list){
        n <- nrow(Y)
        R_y_tau <- cov(Y[1:(n - tau), ], Y[(1 + tau):n, ])
        eig_R_y_tau <- eigen((R_y_tau + t(R_y_tau)) / 2)
        max_eigen_values <- c(max_eigen_values, max(eig_R_y_tau$values))
    }
    tau <- which.max(max_eigen_values)
    R_y_tau <- cov(Y[1:(n - tau), ], Y[(1 + tau):n, ])
    eig_R_y_tau <- eigen((R_y_tau + t(R_y_tau)) / 2)
    W <- eig_R_y_tau$vectors
    A <- ginv(W)
    S <- t(X_centered) %*% W
    X_bar <- S %*% A
    # Output
    list(
        A = ginv(W),
        W = W,
        S = S,
        X = X_bar,
        J = J,
        num.iter = 1,
        thr = NULL,
        # Return only the error at last estimate.
        RecError = .recError(t(X_centered), X_bar),
        RelChange = NULL
    )
}

.SOBI <- function(X, W = NULL, J = NULL, verbose = FALSE, num.iter = 100, tau_list = 1:10, ...){
    num_components <- J
    # Create Givens rotation matrix.
    R_ijcs <- function(i, j, c, s, n, use_complex = FALSE){
        R <- diag(n)
        R[i, i] <- c
        if(use_complex){
            R[i, j] <- -Conj(s)
        }else{
            R[i, j] <- -s
        }
        R[j, i] <- s
        R[j, j] <- c
        R
    }
    # This function treats rows as sensors and columns as samples (time).
    X_centered <- t(X)
    n_covariance_matrices <- length(tau_list)
    # Combines the covariance matrices into an array.
    covariance_matrices <- array(
        0.0,
        dim = c(num_components, num_components, n_covariance_matrices)
    )
    Y <- t(X_centered)
    n <- nrow(Y)
    for(i_tau in 1:n_covariance_matrices){
        tau <- tau_list[i_tau]
        R_tau <- cov(Y[1:(n - tau), ], Y[(1 + tau):n, ])
        covariance_matrices[, , i_tau] <- R_tau
    }
    # NOTE: You may be able to change to an iterative approach by tweaking this.
    V <- diag(num_components)
    for(i_loop in 1:num.iter){
        for(p in 1:num_components){
            for(q in 1:num_components){
                if(p == q) next
                # Use Cyclic-by-Row Algorithm
                # See Golub and Loan Matrix Computations p.430.
                # Size is num_components * 3
                g_r <- matrix(NA, nrow = num_components, ncol = 3)
                for(r in 1:num_components){
                    g_r[r, ] <- c(
                        covariance_matrices[p, p, r] - covariance_matrices[q, q, r],
                        covariance_matrices[p, q, r],
                        covariance_matrices[q, p, r]
                    )
                }
                ReGHG <- Re(t(g_r) %*% g_r)
                ReGHG_eigen <- eigen(ReGHG)
                angles <- ReGHG_eigen$vectors[, 1]
                c_ <- sqrt(0.5 + angles[1] / 2)
                use_complex <- FALSE
                if(use_complex){
                    s_ <- 0.5 * (angles[2] - 1i * angles[3]) / c_
                }else{
                    s_ <- 0.5 * angles[2] / c_
                }
                GivensRot <- R_ijcs(p, q, c_, s_, num_components, use_complex = use_complex)
                V <- V %*% GivensRot
                for(r in 1:num_components){
                    covariance_matrices[, , r] <- t(GivensRot) %*% covariance_matrices[, , r] %*% GivensRot
                }
            }
        }
    }
    W <- V
    A <- ginv(V)
    S <- t(X_centered) %*% W
    X_bar <- S %*% A
    # Output
    list(
        A = ginv(W),
        W = W,
        S = S,
        X = X_bar,
        J = J,
        num.iter = 1,
        thr = NULL,
        # Return only the error at last estimate.
        RecError = .recError(t(X_centered), X_bar),
        RelChange = NULL
    )
}

.FOBI <- function(X, W = NULL, J = NULL, verbose = FALSE, ...){
    X_t <- t(X)
    X_norm <- sqrt(colSums(X_t^2))
    X_mul_norm <- t(t(X_t) * X_norm)
    Omega <- X_mul_norm %*% t(X_mul_norm)
    eigen_Omega <- eigen(Omega)
    unmixing_mat <- t(eigen_Omega$vectors)
    W <- unmixing_mat
    A <- ginv(W)
    S <- t(X_t) %*% W
    X_bar <- S %*% A
    # Output
    list(
        A = A,
        W = W,
        S = S,
        X = X_bar,
        J = J,
        num.iter = 1,
        thr = NULL,
        RecError = .recError(t(X_t), X_bar),
        RelChange = NULL
    )
}

.ProDenICA <- function(
        X, W = NULL, J = NULL,
        num.iter = 100, num_bins = 20,
        verbose = FALSE, ...){
    n_sample <- nrow(X)
    for(i_iter in 1:num.iter){
        W <- .orthgonalize(W)
        S <- X %*% W
        Gs <- matrix(NA, nrow = n_sample, ncol = J)
        G_prime_primes <- matrix(NA, nrow = n_sample, ncol = J)
        for(j in 1:J){
            # Estiamte G
            s_vector <- scale(S[, j])
            bins <- seq(min(s_vector), max(s_vector),, num_bins + 1)
            count_on_grid <- hist(s_vector, breaks = bins, plot = FALSE)
            y <- count_on_grid$counts
            grid <- count_on_grid$mids
            gauss <- dnorm(grid)
            gam_model <- gam(
                y ~ s(grid, k = 6), family = "poisson", offset = log(gauss))
            # Estimate G'
            fine_grid <- seq(min(s_vector), max(s_vector), length.out = 1000)
            G_fine_grid <- predict(
                gam_model, newdata = data.frame(grid = fine_grid), type = "response")
            log_G_fine_grid_without_gauss <- log(G_fine_grid / dnorm(fine_grid))

            gam_model_deriv_log_G_fine_grid_without_gauss <- gam(
                y ~ s(fine_grid, k = 6), family = "gaussian",
                data = list(
                    y = diff(log_G_fine_grid_without_gauss) / diff(range(fine_grid)) * length(fine_grid),
                    fine_grid = fine_grid[2:length(fine_grid)]))
            # Estimate G''
            gam_model_deriv_deriv_log_G_fine_grid_without_gauss <- gam(
                y ~ s(fine_grid, k = 6), family = "gaussian",
                data = list(
                    y = diff(diff(log_G_fine_grid_without_gauss)) / diff(range(fine_grid)) * length(fine_grid),
                    fine_grid = fine_grid[3:length(fine_grid)]))
            Gs[, j] <- predict(
                gam_model_deriv_log_G_fine_grid_without_gauss,
                newdata = data.frame(fine_grid = s_vector), type = "response")
            G_prime_primes[, j] <- predict(
                gam_model_deriv_deriv_log_G_fine_grid_without_gauss,
                newdata = data.frame(fine_grid = s_vector), type = "response")
        }
        # Fixed point algorithm
        X_Gs_div_n_sample <- t(X) %*% Gs / n_sample
        W <- X_Gs_div_n_sample - scale(
            W, center = FALSE, scale = 1 / colMeans(G_prime_primes))
    }
    W <- .orthgonalize(W)
    # Output
    list(
        A = NULL,
        W = W,
        S = NULL,
        X = X,
        J = J,
        num.iter = NULL,
        thr = NULL,
        RecError = NULL,
        RelChange = NULL
    )
}

# Solve Rica with SGD.
.RICA <- function(X, W = NULL, J = NULL, verbose = FALSE, alpha = 0.001, num_epoch = 100, ...){
    # NOTE: RICA can also be used for overcomplete condition.
    K_component <- J
    n_input <- J
    # Within this function, we treat rows as sensors and columns as samples (time).
    num_sample <- nrow(X)
    X_whiten <- t(X)
    rica_loss_and_gradient_one_sample <- function(W, x, K_component, n_input){
        # Calculate with the following configuration:
        # 1. Reconstruction error term
        # 1-1. Forward
        # 1-2. Reverse Direction
        # 2. ICA term
        # 2-1. Forward
        # 2-2. Reverse Direction
        ## 1. Reconstruction error term
        ### 1-1. Forward
        y040 <- t(W) %*% W
        y030 <- c(y040 %*% x)
        y020 <- y030 - x
        y010 <- y020^2
        L <- sum(y010)
        ### 1.2. Reverse Direction
        dL_dy010 <- rep(1, n_input)
        dy010_dy020 <- 2 * y020
        dy020_dy030 <- rep(1, n_input)
        dy030_dy040 <- array(0, dim = c(n_input, n_input, n_input))
        for(i in 1:n_input){
            for(j in 1:n_input){
                for(k in 1:n_input){
                    dy030_dy040[i, j, k] <- (i == k) * x[j]
                }
            }
        }
        dy040_dW <- array(0, dim = c(K_component, n_input, n_input, n_input))
        for(i in 1:K_component){
            for(j in 1:n_input){
                for(k in 1:n_input){
                    for(l in 1:n_input){
                        dy040_dW[i, j, k, l] <- (j == k) * W[i, l] + (j == l) * W[i, k]
                    }
                }
            }
        }
        interim <- einsum("k,ijk->ij", dL_dy010 * dy010_dy020 * dy020_dy030, dy030_dy040)
        diff_f_W_estimate <- einsum("kl,ijkl->ij", interim, dy040_dW)
        ## 2. ICA term
        g <- function(x){
            1/ 2 * log(cosh(2 * x))
        }
        diff_g <- function(x){
            tanh(2 * x)
        }
        ### 2.1. Forward
        z020 <- c(W %*% x)
        z010 <- g(z020)
        Lica <- sum(z010)
        ### 2.2. Reverse Direction
        dLica_dz010 <- rep(1, K_component)
        dz010_dz020 <- diff_g(z020)
        dz020_dW <- array(0, dim = c(K_component, n_input, K_component))
        for(i in 1:K_component){
            for(j in 1:n_input){
                for(k in 1:K_component){
                    dz020_dW[i, j, k] <- (i == k) * x[j]
                }
            }
        }
        diff_W_ICA_estimate <- einsum(
            "k,ijk->ij",
            dLica_dz010 * dz010_dz020, dz020_dW
        )
        # Output
        return(list(
            loss = L + Lica,
            grad = diff_f_W_estimate + diff_W_ICA_estimate,
            L = L, Lica = Lica, diff_f_W_estimate = diff_f_W_estimate,
            diff_W_ICA_estimate = diff_W_ICA_estimate
        ))
    }
    W_ <- W
    for(i_epoch in 1:num_epoch){
        # decay learning rate
        alpha <- alpha * 0.99
        sum_loss <- 0
        sum_L <- 0
        sum_Lica <- 0
        accumulated_grad <- matrix(0, nrow = K_component, ncol = n_input)
        for(i_sample in 1:num_sample){
            xi <- X_whiten[, i_sample]
            result <- rica_loss_and_gradient_one_sample(W_, xi, K_component, n_input)
            accumulated_grad <- accumulated_grad + result$grad
            # W_ <- W_ - alpha * result$grad  # TODO: For SGD, use this line.
            sum_loss <- sum_loss + result$loss
            sum_L <- sum_L + result$L
            sum_Lica <- sum_Lica + result$Lica
        }
        W_ <- W_ - alpha * accumulated_grad / num_sample
        if(verbose){
            cat("epoch: ", i_epoch, " loss: ", sum_loss,
                " L: ", sum_L, " Lica: ", sum_Lica,
                " diff_f_W_estimate: ", result$diff_f_W_estimate,
                " diff_W_ICA_estimate: ", result$diff_W_ICA_estimate,
                " W_: ", W_, ", ",
                sep = "\n"
            )
        }
    }
    # Output
    list(
        A = NULL,
        W = t(W_),
        S = NULL,
        X = X,
        J = J,
        num.iter = NULL,
        thr = NULL,
        RecError = NULL,
        RelChange = NULL
    )
}

.checkICA2 <- function(X, J, num.iter, thr, r_list,
    omega_for_each_r, a_r_for_each_r, tau_list,
    num_bins, alpha, num_epoch, verbose){
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(J))
    stopifnot(length(J) == 1)
    stopifnot(min(dim(X)) >= J)
    if(!is.null(num.iter)){
        stopifnot(is.numeric(num.iter))
        stopifnot(num.iter >= 1)
    }
    stopifnot(is.numeric(thr))
    if(!is.null(r_list)){
        stopifnot(is.numeric(r_list))
    }
    if(!is.null(omega_for_each_r)){
        stopifnot(is.numeric(omega_for_each_r))
    }
    if(!is.null(a_r_for_each_r)){
        stopifnot(is.numeric(a_r_for_each_r))
    }
    if(!is.null(tau_list)){
        stopifnot(is.numeric(tau_list))
    }
    if(!is.null(num_bins)){
        stopifnot(is.numeric(num_bins))
    }
    if(!is.null(alpha)){
        stopifnot(is.numeric(alpha))
    }
    if(!is.null(num_epoch)){
        stopifnot(is.numeric(num_epoch))
    }
    stopifnot(is.logical(verbose))
}

.init_X <- function(X, J, scale_X = TRUE){
    # Gets the mean and standard deviation vectors of X.
    XMeans <- colMeans(X)
    XSds <- apply(X, 1, sd)
    # Whitens X.
    WhiteningResult <- .whitening2(X, scale = scale_X)
    X <- WhiteningResult$XWhiten
    WhiteningMatrix <- WhiteningResult$WhiteningMatrix
    # Truncate X according to the number of elements specified.
    # if J <ncol (X), truncate.
    # if ncol (X) <= J, leave as is.
    if(J < ncol(X)){
        X <- X[, 1:J]
    }
    # Output
    list(X = X,
        XMeans = XMeans,
        XSds = XSds,
        WhiteningMatrix = WhiteningMatrix
    )
}

.initICA2_1 <- function(X, J, algorithm, num.iter, thr, r_list,
    omega_for_each_r, a_r_for_each_r, tau_list,
    num_bins, alpha, num_epoch, verbose){
    # Algorithm's default parameters
    if(algorithm == "JADE"){
        num.iter <- if(is.null(num.iter)) 500 else num.iter
    }else if(algorithm == "AuxICA1"){
        num.iter <- if(is.null(num.iter)) 30 else num.iter
    }else if(algorithm == "AuxICA2"){
        num.iter <- if(is.null(num.iter)) 30 else num.iter
    }else if(algorithm == "IPCA"){
        num.iter <- if(is.null(num.iter)) 200 else num.iter
    }else if(algorithm == "SIMBEC"){
        num.iter <- if(is.null(num.iter)) 1000 else num.iter
        r_list <- if(is.null(r_list)) c(4) else r_list
        omega_for_each_r <- if(is.null(omega_for_each_r)) c(1.0) else omega_for_each_r
        a_r_for_each_r <- if(is.null(a_r_for_each_r)) c(1.0) else a_r_for_each_r
    }else if(algorithm == "AMUSE"){
        num.iter <- if(is.null(num.iter)) 1 else num.iter
        tau_list <- if(is.null(tau_list)) 1:50 else tau_list
    }else if(algorithm == "SOBI"){
        num.iter <- if(is.null(num.iter)) 100 else num.iter
        tau_list <- if(is.null(tau_list)) 1:10 else tau_list
    }else if(algorithm == "FOBI"){
        num.iter <- if(is.null(num.iter)) 100 else num.iter
    }else if(algorithm == "ProDenICA"){
        num.iter <- if(is.null(num.iter)) 100 else num.iter
        num_bins <- if(is.null(num_bins)) 20 else num_bins
    }else if(algorithm == "RICA"){
        num.iter <- if(is.null(num.iter)) 100 else num.iter
        alpha <- if(is.null(alpha)) 0.001 else alpha
        num_epoch <- if(is.null(num_epoch)) 100 else num_epoch
    }
    # Output
    list(X = X, J = J, algorithm = algorithm,
        num.iter = num.iter, thr = thr,
        r_list = r_list,
        omega_for_each_r = omega_for_each_r,
        a_r_for_each_r = a_r_for_each_r,
        tau_list = tau_list,
        num_bins = num_bins,
        alpha = alpha,
        num_epoch = num_epoch,
        verbose = verbose)
}

.initICA2_2 <- function(X, J, thr, allow_overcomplete){
    X_initialized <- .init_X(X, J)
    X <- X_initialized$X

    num_columns <- ncol(X)
    if(allow_overcomplete){
        W <- matrix(runif(num_columns * J), nrow = num_columns, ncol = J)
        A <- ginv(matrix(runif(J * num_columns), nrow = J, ncol = num_columns))
    }else{
        W <- matrix(runif(J * J), nrow = J, ncol = J)
        A <- ginv(W)
    }
    S <- X %*% W
    XBar <- S %*% A
    # Reconstruction Error
    RecError <- .recError(X, XBar)
    # Relative Change
    RelChange <- thr * 10
    # Output
    list(
        X = X, A = A, W = W, S = S,
        RecError = RecError, RelChange = RelChange,
        WhiteningMatrix = X_initialized$WhiteningMatrix,
        XMeans = X_initialized$XMeans, XSds = X_initialized$XSds
    )
}

.whitening2 <- function(X, scale = TRUE){
    num_sample <- nrow(X)
    X_scaled <- scale(X, scale = scale)
    eigen_X_scaled <- eigen(t(X_scaled) %*% X_scaled / num_sample)
    W <- eigen_X_scaled$vectors %*% diag(eigen_X_scaled$values^(-1 / 2)) %*% t(eigen_X_scaled$vectors)
    XWhiten <- X_scaled %*% W
    list(XWhiten = XWhiten, WhiteningMatrix = W)
}

.recError <- function(X, XBar){
    sum((X - XBar)^2) / sum(X^2)
}

.relChange <- function(iter, RecError){
    (RecError[iter - 1] - RecError[iter]) / RecError[iter - 1]
}

.FastICA2 <- function(X, W, S, J){
    n <- nrow(X)
    J <- ncol(X)
    # In order to conform to the notation of the paper,
    # W and X is treated as transposed.
    W <- t(W) # J * J
    X <- t(X) # J * n
    #########################################
    # Symmetric FastICA
    #########################################
    nonlinear_transform_function_name <- "tanh"
    # nonlinear_transform_function_name <- "exp"
    # nonlinear_transform_function_name <- "kurtosis"
    if(nonlinear_transform_function_name == "tanh"){
        alpha <- 1.0
        g <- function(x) tanh(alpha * x)
        gPrime <- function(x) alpha * (1 - tanh(alpha * x)^2)
    }else if(nonlinear_transform_function_name == "exp"){
        g <- function(x) x * exp(-x * 2 / 2)
        gPrime <- function(x) (1 - x^2) * exp(-x^2 / 2)
    }else if(nonlinear_transform_function_name == "kurtosis"){
        g <- function(x) x^3
        gPrime <- function(x) 3 * x^2
    }
    WX <- W %*% X # J * n
    gWX <- g(WX) # J * n
    gPrimeWX <- gPrime(WX) # J * n
    EgPrimeWX <- apply(gPrimeWX, 1, FUN = mean) # J
    EgPrimeWXW <- W %*% diag(EgPrimeWX) # J * J
    XgWX <- gWX %*% t(X) / n # J * J
    WNew <- XgWX - EgPrimeWXW # J * J
    # Orthgonalize
    WNewSvd <- svd(WNew)
    WNew <- tcrossprod(WNewSvd$u, WNewSvd$v)
    # Update
    t(WNew)
}

.ica_function_list <- list(
    "JADE" = .JADE,
    "AuxICA1" = .AuxICA1,
    "AuxICA2" = .AuxICA2,
    "IPCA" = .IPCA,
    "SIMBEC" = .SIMBEC,
    "AMUSE" = .AMUSE,
    "SOBI" = .SOBI,
    "FOBI" = .FOBI,
    "ProDenICA" = .ProDenICA,
    "RICA" = .RICA
)
