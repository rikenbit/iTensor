.ICA_Type1 <- function(){
    # Generate independent components.
    num_sample <- 500
    S_true <- cbind(
        runif(num_sample, min = -1, max = 1),
        runif(num_sample, min = -1, max = 1),
        runif(num_sample, min = -1, max = 1))
    # Get the number of elements.
    num_components <- ncol(S_true)
    # Generates a random mixing matrix.
    A_true <- matrix(
        runif(num_components^2, min = -1, max = 1),
        ncol = num_components)
    # Generate observed data.
    X_observed <- S_true %*% A_true
    list(X_observed=X_observed, S_true=S_true, A_true=A_true)
}

.ICA_Type2 <- function(){
    # Generate independent components.
    num_sample <- 500
    S_true <- cbind(
            rt(num_sample, df = 3),
            rt(num_sample, df = 3),
            rt(num_sample, df = 3)
        )
    num_components <- ncol(S_true)
    # Generate a random mixing matrix.
    A_true <- matrix(
        runif(num_components^2, min = -1, max = 1),
        ncol = num_components)
    X_observed <- S_true %*% A_true
    list(X_observed=X_observed, S_true=S_true, A_true=A_true)
}

.ICA_Type3 <- function(){
    # Generate independent components.
    num_sample <- 500
    S_true <- cbind(
            rt(num_sample, df = 1),
            rnorm(num_sample),
            runif(num_sample, min = -1, max = 1)
        )
    # Generate a random mixing matrix.
    num_components <- ncol(S_true)
    A_true <- matrix(
        runif(num_components^2, min = -1, max = 1),
        ncol = num_components)
    X_observed <- S_true %*% A_true
    list(X_observed=X_observed, S_true=S_true, A_true=A_true)
}

.ICA_Type4 <- function(){
    num_sample <- 500
    # Moving average window size.
    n_smooth <- 10
    # Generate independent components.
    S_true <- cbind(
      sin((1:num_sample) / 3),
      sin((1:num_sample) / 4),
      cos((1:num_sample) / 5))
    #  Smooth by moving average.
    for (i_var in 1:ncol(S_true)) {
      S_true[, i_var] <- filter(S_true[, i_var], rep(1, n_smooth)/n_smooth)
    }
    S_true <- S_true[(n_smooth + 1):(nrow(S_true) - n_smooth), ]
    num_sample <- nrow(S_true)
    num_components <- ncol(S_true)
    # Generate a random mixing matrix.
    A_true <- matrix(
        runif(num_components^2, min = -1, max = 1),
        ncol = num_components)
    X_observed <- S_true %*% A_true
    list(X_observed=X_observed, S_true=S_true, A_true=A_true)
}

.ICA_Type5 <- function(){
    env <- new.env()
    data(liver.toxicity, package="mixOmics", envir=env)
    env$liver.toxicity
}

.MICA <- function(){
    # Generates an observation vector of three mixed signals.
    # Read `5.Simulation` section in MICA paper.
    t_series <- seq(from = 0.00, to = 1.000, by = 1e-4)
    u1 <- runif(n = length(t_series), min = -1, max = 1)
    u2 <- sin(2 * pi * 800 * t_series + 6 * cos(2 * pi * 60 * t_series))
    u3 <- sin(2 * pi * 90 * t_series)
    u <- cbind(u3, u2, u1)
    v1 <- abs(u1)
    v2 <- abs(u2)
    v3 <- abs(u3)
    v <- cbind(v3, v2, v1)
    A_true <- matrix(runif(n = 3 * 3, min = -1, max = 1), nrow = 3, ncol = 3)
    B_true <- matrix(runif(n = 3 * 3, min = -1, max = 1), nrow = 3, ncol = 3)
    X <- u %*% A_true
    Y <- v %*% B_true
    list(X=X, Y=Y)
}

.GroupICA <- function(){
    # Generate data from a block-wise variance model
    d <- 6  # modified from 2 to 6 because Calhoun2009 with JADE does not work with d=2
    m <- 10
    n <- 5000
    group_index <- rep(c(1,2), each=n)
    partition_index <- rep(rep(1:m, each=n/m), 2)
    S <- matrix(NA, 2*n, d)
    H <- matrix(NA, 2*n, d)
    for(i in unique(group_index)){
      varH <- abs(rnorm(d))/4
      H[group_index==i, ] <- matrix(rnorm(d*n)*rep(varH, each=n), n, d)
      for(j in unique(partition_index[group_index==i])){
        varS <- abs(rnorm(d))
        index <- partition_index==j & group_index==i
        S[index,] <- matrix(rnorm(d*n/m)*rep(varS, each=n/m), n/m, d)
      }
    }
    A <- matrix(rnorm(d^2), d, d)
    A <- A%*%t(A)
    X <- t(A%*%t(S + H))
    Xs <- list()
    for (i_group in unique(group_index)) {
        Xs[[i_group]] <- X[group_index==i_group, ]
    }
    list(X=Xs[[1]], Y=Xs[[2]])
}

.flist <- list(
    ICA_Type1 = .ICA_Type1,
    ICA_Type2 = .ICA_Type2,
    ICA_Type3 = .ICA_Type3,
    ICA_Type4 = .ICA_Type4,
    ICA_Type5 = .ICA_Type5,
    MICA = .MICA,
    GroupICA = .GroupICA
)

#' Toy model data for using ICA, MICA, and GroupICA
#' There are 7 types of simulation:
#' ICA_Type1: Time-independent sub-gaussian data
#' ICA_Type2: Time-independent super-gaussian data
#' ICA_Type3: Data mixed with signals having no time dependence and different kurtosis
#' ICA_Type4: Time-dependent data
#' ICA_Type5: Toydata to model IPCA in N < P systems
#' MICA: Two time-serices data to model MICA
#' GroupICA: Toydata to model GroupICA
#' @param model "ICA_Type1", "ICA_Type2", "ICA_Type3", "ICA_Type4",
#' and "ICA_Type5", "MICA", and "GrouICA" are available
#' (Default: "ICA_Type1").
#' @param seeds Random number for setting set.seeds in the function (Default: 123).
#' @return A list containing simulation data sets.
#' @examples
#' data1 <- toyModel("ICA_Type1")
#' data2 <- toyModel("ICA_Type2")
#' data3 <- toyModel("ICA_Type3")
#' data4 <- toyModel("ICA_Type4")
#' data5 <- toyModel("ICA_Type5")
#' data6 <- toyModel("MICA")
#' data7 <- toyModel("GroupICA")
#' @importFrom utils data
#' @importFrom groupICA groupICA
#' @importFrom stats runif rt filter
#' @import mixOmics
#' @export
toyModel <- function(model = "ICA_Type1", seeds=123){
    set.seed(seeds)
    out <- .flist[[model]]()
    set.seed(NULL)
    out
}
