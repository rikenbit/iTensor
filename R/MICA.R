#' Multimodal independent component analysis
#'
#' The input datasets are assumed to be two matrices sharing the column space.
#' MICA decomposes the matrices simutanously
#' and extracts the components that maximizes the mutual information
#' between the components.
#' @param X A matrix sharing the column space with Y (??? x N)
#' @param Y A matrix sharing the column space with X (??? x N)
#' @param J The rank parameter to decompose the matrices
#' @param eta A learning rate parameter of stochastic gradient descent
#' @param verbose Verbose option
#' @param mu A learning rate parameter of stochastic gradient descent
#' @param gamma_ts Weighting factor for dependence on independence
#' @return A list containing the result of the decomposition
#' @examples
#' X <- array(runif(10*20), dim=c(10,20))
#' Y <- array(runif(15*20), dim=c(15,20))
#' J <- 20
#' out <- MICA(X, Y, J=J)
#' @importFrom stats rnorm
#' @export
MICA <- function(X, Y, J, eta=1000*1e-4, verbose=FALSE,
    mu=50.0*1e-4, gamma_ts=1.0){
    ###########
    # gamma_ts example:
    #   gamma_ts <- 0.2 / (c(1:nSample) + 1)  # used in simulation study of MICA paper
    #   gamma_ts <- rep(0, nSample)           # only independence is evaluated.
    #   gamma_ts <- rep(1, nSample)           # only dependence is evaluated.
    ######################################
    # Argument Check
    ######################################
    .checkMICA(X, Y, J, eta, verbose)
    ######################################
    # Initialization (e.g. CCA)
    ######################################
    int <- .initMICA(X, Y, J)
    A <- int$A
    B <- int$B
    X <- int$XScaled
    Y <- int$YScaled
    nSample <- nrow(X)
    N <- J  # Match the notation for MICA paper.
    ABChange <- 1
    # Time varying learning rate for online estimation of statistics
    if (length(mu) < 2){
        # if constant mu is given, convert it to timeseries.
        mu <- rep(mu, nSample)
    }
    # Weighting factor for dependence on independence
    if (length(gamma_ts) < 2){
        # if constant gamma_ts is given, convert it to timeseries.
        gamma_ts <- rep(gamma_ts, nSample)
    }
    # None of the diagonal elements must be 1.
    # because rho_bar_i is computed using the term R[i,i] / (1 - R[i,i]^2).
    R <- matrix(0, ncol = N, nrow = N)
    # initialize \beta_i^{k, l}
    # definition: \beta_i^{k, l} == beta_i_k_l__zero_origin[i, k+1, l+1]
    beta_i_k_l__zero_origin <- array(data = 1, dim = c(N, 5, 5))
    for (i in 1:N){
        for (k in 0:4){
            for (l in 0:4){
                # use joint cumulants of random normal variable
                beta_i_k_l__zero_origin[i, k+1, l+1] <- mean(
                    rnorm(1000)^k * rnorm(1000)^l) - .beta_zero_k_l(k, l)
            }
        }
    }
    ######################################
    # Iteration
    ######################################
    # use online algorithm
    for (iter in 1:nSample){
        u <- as.numeric(A %*% as.numeric(X[iter, ]))  # size = (N, )
        v <- as.numeric(B %*% as.numeric(Y[iter, ]))  # size = (N, )
        r <- rep(0, N)  # size = (N, )
        s <- rep(0, N)  # size = (N, )
        # Make u and v to uncorrelated variables r and s.
        for (i in 1:N){
            rho_i <- R[i, i]
            c_plus_i <- 1 / 2 * ((1 + rho_i)^(-1/2) + (1 - rho_i)^(-1/2))
            c_minus_i <- 1 / 2 * ((1 + rho_i)^(-1/2) - (1 - rho_i)^(-1/2))
            r[i] <- c_plus_i * u[i] +  c_minus_i * v[i]
            s[i] <- c_minus_i * u[i] +  c_plus_i * v[i]
        }
        # Online estimation of statistic beta_i^{k,l}
        for (i in 1:N){
            for (k in 0:4){
                for (l in 0:4){
                    beta_i_k_l__zero_origin[i, k+1, l+1] <-
                        beta_i_k_l__zero_origin[i, k+1, l+1] -
                            mu[iter] * (
                                beta_i_k_l__zero_origin[i, k+1, l+1] -
                                    r[i]^k * s[i]^l + .beta_zero_k_l(k, l))
                }
            }
        }
        # Online estimation of statistic R
        for (i in 1:N){
            for (j in 1:N){
                R[i, j] <- R[i, j] - mu[iter] * (R[i, j] - u[i]*v[j])
            }
        }
        # .gradA() means (-1) * ∂ε_tot/∂A A^TA
        A_new <- A - eta * .gradA(u, v, A, B, N, J, beta_i_k_l__zero_origin, R, r, s, gamma_ts[iter])
        B_new <- B - eta * .gradB(u, v, A, B, N, J, beta_i_k_l__zero_origin, R, r, s, gamma_ts[iter])
        # Orthgonalize
        A_new <- .orthgonalize(A_new)
        B_new <- .orthgonalize(B_new)
        # Record change of A and B
        # ABChange[iter] <- max(max(abs(A_new - A)), max(abs(B_new - B)))
        ABChange[iter] <- (mean(abs(A_new - A)) + mean(abs(B_new - B))) / 2
        # Verbose
        if (verbose){
             cat(paste0(iter,
                " |max(max(abs(A_new - A)), max(abs(B_new - B)))| = ",
                ABChange[iter], "\n"))
        }
        # Update
        A <- A_new
        B <- B_new
    }
    # Estimate original source using last A and B.
    # Note that this method is difficult in an online setting.
    U <- X %*% t(A)
    V <- Y %*% t(B)
    # Output
    list(U=U, V=V, A=A, B=B, J=J, eta=eta, num.iter=iter,
         verbose=verbose, ABChange=ABChange)
}

.checkMICA <- function(X, Y, J, eta, verbose){
    stopifnot(is.matrix(X))
    stopifnot(is.matrix(Y))
    stopifnot(is.numeric(J))
    stopifnot(length(J) == 1)
    stopifnot(is.numeric(eta))
    stopifnot(length(eta) == 1)
    stopifnot(min(ncol(X), ncol(Y)) >= J)
    stopifnot(is.logical(verbose))
}

.initMICA <- function(X, Y, J){
    XScaled <- scale(X)
    YScaled <- scale(Y)
    A <- diag(J)
    B <- diag(J)
    list(A=A, B=B, XScaled=XScaled, YScaled=YScaled)
}

# factorial function
.fact <- function(x){
    ifelse (x != floor(x) | x < 0, NA, ifelse(x == 0 , 1, prod(1:x)))
}

# phi() function are used in Section: 6  in
# Yang, Howard Hua, and Shun-ichi Amari.
# "Adaptive online learning algorithms for blind separation:
# maximum entropy and minimum mutual information." Neural computation 9.7 (1997): 1457-1482.
.phi <- function(y){2 * tanh(y)}
# phi <- function(y){ y^3 }
# phi <- function(y){ 3 / 4 * y^11 + 15 / 4 * y^9 - 14 / 3 * y^7 - 29 / 4 * y^5 + 29 / 4 * y^3 }

.beta_zero_k_l <- function(k, l){
    if ((k == 4) || (l == 4)){
        return(3)
    } else if ((k == 2) && (l == 2)){
        return(1)
    } else {
        return(0)
    }
}

.der_H_beta_i_der_beta_i_k_l <- function(i, k, l, beta_i_k_l__zero_origin){
    if ((k==3) && (l==0)){
        return(- 1 / (2 * .fact(3) ) * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else if ((k==2) && (l==1)){
        return(- 1 / (2 * .fact(3) ) * 3 * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else if ((k==1) && (l==2)){
        return(- 1 / (2 * .fact(3) ) * 3 * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else if ((k==0) && (l==3)){
        return(- 1 / (2 * .fact(3) ) * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else if ((k==4) && (l==0)){
        return(- 1 / (2 * .fact(4)) * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else if ((k==3) && (l==1)){
        return(- 1 / (2 * .fact(4)) * 4 * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else if ((k==2) && (l==2)){
        return(- 1 / (2 * .fact(4)) * 6 * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else if ((k==1) && (l==3)){
        return(- 1 / (2 * .fact(4)) * 4 * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else if ((k==0) && (l==4)){
        return(- 1 / (2 * .fact(4)) * 2 * .beta_i_k_l(i, k, l, beta_i_k_l__zero_origin))
    } else {
        stop("given term is not included in this approximation.")
    }
}

.beta_i_k_l <- function(i, k, l, beta_i_k_l__zero_origin){
    beta_i_k_l__zero_origin[i, k+1, l+1]
}

.der_riksil_der_risiT <- function(k, l, r_i, s_i){
    if (k > 0){
        der_riksil_der_ri <- k * r_i^(k-1) * s_i^l
    } else {
        der_riksil_der_ri <- 0
    }
    if (l > 0){
        der_riksil_der_si <- l * s_i^(l-1) * r_i^k
    } else {
        der_riksil_der_si <- 0
    }
    return(c(der_riksil_der_ri, der_riksil_der_si))
}

.gradA <- function(u, v, A, B, N, J, beta_i_k_l__zero_origin, R, r, s, gamma_t){

    Phi <- rep(0, N)      # size = (N, )
    Xi <- rep(0, N)       # size = (N, )
    Psi <- rep(0, N)      # size = (N, )
    Rho_bar <- rep(0, J)  # size = (N, )

    for (i in 1:N){

        rho_i <- R[i, i]
        c_plus_i <- 1 / 2 * ((1 + rho_i)^(-1/2) + (1 - rho_i)^(-1/2))
        c_minus_i <- 1 / 2 * ((1 + rho_i)^(-1/2) - (1 - rho_i)^(-1/2))

        f_i_g_i <- c(0, 0)
        for (k in 0:5){
            for (l in 0:5){
                # Only all subscript sets of k and l contained in equation (12) are considered.
                if ((k + l <= 2) || (k + l > 4)) next
                f_i_g_i <- f_i_g_i +
                    .der_H_beta_i_der_beta_i_k_l(i, k, l, beta_i_k_l__zero_origin) *
                    .der_riksil_der_risiT(k, l, r[i], s[i])
            }
        }
        f_i <- as.numeric(f_i_g_i[1])
        g_i <- as.numeric(f_i_g_i[2])

        xi_i <- c_plus_i * f_i + c_minus_i * g_i

        d_plus_i <- - 1 / 4 * ((1 + rho_i)^(-3/2) + (1 - rho_i)^(-3/2))
        d_minus_i <- - 1 / 4 * ((1 + rho_i)^(-3/2) - (1 - rho_i)^(-3/2))

        psi_i <- (d_plus_i * u[i] + d_minus_i * v[i]) * f_i + (d_minus_i * u[i] + d_plus_i * v[i]) * g_i
        rho_bar_i <- rho_i / (1-rho_i^2)
        # store them.
        Phi[i] <- .phi(u[i])
        Xi[i] <- xi_i
        Psi[i] <- psi_i
        Rho_bar[i] <- rho_bar_i
    }
    der_A_der_t_div_minus_eta_t <-
        (diag(N) + (1 - gamma_t) * Phi %*% t(u) - gamma_t * Xi %*% t(u) +
             gamma_t * diag(Rho_bar - Psi) %*% t(R)) %*% A
    return(- der_A_der_t_div_minus_eta_t)
}

.gradB <- function(u, v, A, B, N, J, beta_i_k_l__zero_origin, R, r, s, gamma_t){
    .gradA(v, u, B, A, N, J,
        aperm(beta_i_k_l__zero_origin,
            c(1,3,2)), t(R), s, r, gamma_t)
}
