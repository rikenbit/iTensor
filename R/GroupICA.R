#' Group Independent Component Analysis (GroupICA)
#'
#' The input data is assumed to be a list containing multiple matrices, which share common column.
#' @param Xs A list containing multiple matrices
#' @param J1 Rank parameter to decompose
#' @param J2 Rank parameter used in Calhoun2009
#' @param algorithm Pool algorithm to merge multiple ICA results (Default: pooled)
#' @param ica.algorithm The decomposition algorithm (Default: "FastICA")
#' @param num.iter The number of iterations
#' @param thr The threshold to terminate the iteration (Default: 1E-10)
#' @param verbose Verbose option
#' @return A list containing the result of the decomposition
#' @examples
#' X1 <- matrix(runif(100*200), nrow=100, ncol=200)
#' X2 <- matrix(runif(150*200), nrow=150, ncol=200)
#' Xs <- list(X1=X1, X2=X2)
#' out <- GroupICA(Xs, J1=5)
#' @importFrom jointDiag uwedge
#' @importFrom utils str
#' @importFrom stats prcomp
#' @export
GroupICA <- function(Xs, J1, J2 = J1,
    algorithm = c("pooled", "Calhoun2009", "Pfister2018"),
    ica.algorithm = c("FastICA", "InfoMax", "ExtInfoMax",
        "JADE", "AuxICA1", "AuxICA2", "IPCA", "SIMBEC",
        "AMUSE", "SOBI", "FOBI", "ProDenICA", "RICA"),
    num.iter = 30, thr = 1E-10, verbose = FALSE){
    ######################################
    # Argument Check
    ######################################
    algorithm <- match.arg(algorithm)
    ica.algorithm <- match.arg(ica.algorithm)
    .checkGroupICA(Xs, J1, J2, algorithm, num.iter, thr, verbose)
    l <- length(Xs)
    if(algorithm == "pooled"){
        out <- .pooled(Xs, J1, J2, ica.algorithm, num.iter, thr, verbose)
    }
    if(algorithm == "Calhoun2009"){
        out <- .Calhoun2009(Xs, J1, J2, ica.algorithm, num.iter, thr, verbose)
    }
    if(algorithm == "Pfister2018"){
        out <- .Pfister2018(Xs, J1, J2, ica.algorithm, num.iter, thr, verbose)
    }
    # Output
    list(A = out$A, Ss = out$Ss, RecError = out$RecError, RelChange = out$RelChange)
}

.pooled <- function(Xs, J1, J2, ica.algorithm, num.iter, thr, verbose){
    Xpooled <- do.call(rbind, Xs)
    if(ica.algorithm %in% c("FastICA", "InfoMax", "ExtInfoMax")){
        out <- ICA(Xpooled, J = J1, algorithm = ica.algorithm,
            num.iter = num.iter, thr = thr, verbose = verbose)
    }else{
    out <- ICA2(Xpooled, J = J1, algorithm = ica.algorithm,
        num.iter = num.iter, thr = thr, verbose = verbose)
    }
    A <- out$A
    group_index <- cumsum(as.numeric(unlist(lapply(
        Xs, function(x){
            l <- numeric(nrow(x))
            l[1] <- 1
            l
        }))))
    Ss <- list()
    for(i_group in seq_len(max(group_index))){
        Ss[[i_group]] <- out$S[group_index == i_group, ]
    }
    RecError <- out$RecError
    RelChange <- out$RelChange
    list(A = A, Ss = Ss, RecError = RecError, RelChange = RelChange)
}

.Calhoun2009 <- function(Xs, J1, J2, ica.algorithm, num.iter, thr, verbose){
    # The variable names are the same as those in the paper.
    # K is the number of samples in each group
    Ks <- lapply(Xs, function(x){
        nrow(x)
    })
    # L is the size of time dimension following reduction.
    L <- J2
    # V is the number of input components
    V <- ncol(Xs[[1]])
    # M is the number of groups
    M <- length(Xs)
    # N is the number of components to be estimated
    N <- J1
    if(verbose){
        print(paste("Ks:", Ks))
        print(paste("L:", L))
        print(paste("V:", V))
        print(paste("M:", M))
        print(paste("N:", N))
    }
    # each element of Ys is a K x V matrix
    Ys <- lapply(Xs, function(x){
        # scale adds unnecessary attributes, so we cut them off with `[,]`.
        scale(x, center=TRUE, scale=FALSE)[,]
    })
    if(verbose){
        print(paste("Ys:", lapply(Ys, dim)))
        print("contents of Ys:")
        print(str(Ys))
    }
    Y_prcomps <- lapply(Ys, function(Y){
        prcomp(t(Y), center=TRUE, scale=FALSE)
    })
    # each element of F_inv_Ys is a L x V matrix
    F_inv_Ys <- lapply(Y_prcomps, function(Y_prcomp){
        t(Y_prcomp$x[, seq_len(L)])
    })
    if(verbose){
        print(paste("F_inv_Ys:", lapply(F_inv_Ys, dim)))
        print("contents of F_inv_Ys:")
        print(str(F_inv_Ys))
    }
    # each element of F_invs is a L x K matrix
    F_invs <- lapply(Y_prcomps, function(Y_prcomp){
        t(Y_prcomp$rotation[, seq_len(L)])
    })
    if(verbose){
        print(paste("F_invs:", lapply(F_invs, dim)))
        print("contents of F_invs:")
        print(str(F_invs))
    }
    # F_inv_Y_concatenated ... (LM x V)
    F_inv_Y_concatenated <- do.call(rbind, F_inv_Ys)
    if(verbose){
        print("F_inv_Y_concatenated:")
        print(dim(F_inv_Y_concatenated))
    }
    # G_inv ... (N x LM)
    G_inv <- t(prcomp(t(F_inv_Y_concatenated), center=TRUE, scale=FALSE)$rotation[, seq_len(N)])
    if(verbose){
        print("G_inv:")
        print(dim(G_inv))
        print("contents of G_inv:")
        print(G_inv)
    }
    X <- t(prcomp(t(F_inv_Y_concatenated), center=TRUE, scale=FALSE)$x[, seq_len(N)])
    if(verbose){
        print("X:")
        print(dim(X))
        print("contents of X:")
        print(X)
    }
    if(ica.algorithm %in% c("FastICA", "InfoMax", "ExtInfoMax")){
        X.ica <- ICA(X, J = V, algorithm = ica.algorithm,
            num.iter = num.iter, thr = thr, verbose = verbose)
    }else{
        X.ica <- ICA2(X, J = V, algorithm = ica.algorithm,
            num.iter = num.iter, thr = thr, verbose = verbose)
    }
    # A ... N x N
    A <- X.ica$A
    if(verbose){
        print("A:")
        print(dim(A))
        print("contents of A:")
        print(A)
    }
    # S ... N x V
    S <- X.ica$S
    if(verbose){
        print("S:")
        print(dim(S))
        print("contents of S:")
        print(S)
    }
    # G ... LM x N
    G <- ginv(G_inv)
    if(verbose){
        print("G:")
        print(dim(G))
        print("contents of G:")
        print(G)
    }
    # each element of Gs is a L x N matrix
    Gs <- list()
    for(i_group in seq_len(M)){
        Gs[[i_group]] <- G[seq((i_group - 1) * L + 1, i_group * L), ]
    }
    if(verbose){
        print(paste("Gs:", lapply(Gs, dim)))
        print("contents of Gs:")
        print(str(Gs))
    }
    # each element of Ss is a N x V matrix
    Ss <- list()
    for(i_group in seq_len(M)){
        Ss[[i_group]] <- ginv(Gs[[i_group]] %*% A) %*% F_inv_Ys[[i_group]]
    }
    if(verbose){
        print(paste("Ss:", lapply(Ss, dim)))
        print("contents of Ss:")
        print(str(Ss))
    }
    # calculate ICA time courses
    # each element of FGAs is a K x N matrix
    FGAs <- lapply(seq_len(M), function(i_group){
        # F_invs[[i_group]] ... L x K
        # Gs[[i_group]] ... L x N
        # A ... N x N
        ginv(F_invs[[i_group]]) %*% Gs[[i_group]] %*% A
    })
    if(verbose){
        print(paste("FGAs:", lapply(FGAs, dim)))
        print("contents of FGAs:")
        print(str(FGAs))
    }
    # Because F_i G_i A is considered to be the original signal,
    # we output it as S_i.
    # Note that S_i is N x V.
    Ss <- FGAs
    RecError <- NULL
    RelChange <- NULL
    list(A = A, Ss = Ss, RecError = RecError, RelChange = RelChange)
}

.Pfister2018 <- function(Xs, J1, J2, ica.algorithm, num.iter, thr, verbose){
    ######################################
    # Initialization (e.g. Whiteniing)
    ######################################
    int <- .initGroupICA(Xs)
    X <- int$X
    g <- int$g    # group index
    Pg <- int$Pg  # group-wise partition
    M <- int$M    # empty list
    whitening_result <- .whitening2(X)
    X <- whitening_result$XWhiten
    X <- X[, seq_len(J1)]
    groups <- unique(g)
    for(g_i in groups){
        partitions <- unique(Pg[g == g_i])
        for(Pg_i in partitions){
            e <- which(g == g_i & Pg == Pg_i)
            # which(g != g_i | Pg != Pg_i) is not necessarily correct.
            # which(g != g_i & Pg == Pg_i) may be correct.
            # e_complement <- which(g != g_i | Pg != Pg_i)
            e_complement <- which(g != g_i | Pg == Pg_i)
            M[[length(M) + 1]] <- cov(X[e, ]) - cov(X[e_complement, ])
            M[[length(M) + 1]] <- cov(X[e_complement, ]) - cov(X[e, ])
        }
    }
    n_eigen_matrices_diagonalized <- length(M)
    M_stacked <- array(0.0, dim = c(J1, J1, n_eigen_matrices_diagonalized))
    for(i in seq_len(n_eigen_matrices_diagonalized)){
        M_stacked[, , i] <- M[[i]]
    }
    A <- .ApproximateJointDiagonalizer(M_stacked)
    # Output
    S <- X %*% A
    group_index <- cumsum(as.numeric(unlist(lapply(
        Xs, function(x){
            l <- numeric(nrow(x))
            l[1] <- 1
            l
        }))))
    Ss <- list()
    for(i_group in seq_len(max(group_index))){
        Ss[[i_group]] <- S[group_index == i_group, ]
    }
    RecError <- NULL
    RelChange <- NULL
    list(A = A, Ss = Ss, RecError = RecError, RelChange = RelChange)
}

.checkGroupICA <- function(Xs, J1, J2, algorithm, num.iter, thr, verbose){
    stopifnot(is.list(Xs))
    nr <- lapply(Xs, nrow)
    all.equal(length(unique(nr)), 1)
    stopifnot(is.numeric(J1))
    stopifnot(length(J1) == 1)
    lapply(Xs, function(x){
        stopifnot(min(dim(x)) >= J1)
    })
    if(algorithm == "Calhoun2009"){
        stopifnot(is.numeric(J2))
        stopifnot(length(J2) == 1)
        # Currently, J2 is treated as $L$ in the paper.
        for(i in seq_len(length(Xs))){
            stopifnot(nrow(Xs[[i]]) >= J2)
            stopifnot(ncol(Xs[[i]]) >= J2)
        }
    }
    stopifnot(is.numeric(num.iter))
    stopifnot(num.iter > 0)
    stopifnot(is.numeric(thr))
    stopifnot(is.logical(verbose))
}

.initGroupICA <- function(Xs, min_sample_per_group = 200){
    # Each X is considered to be each group.
    # For Pfister2018, Xs is a list of data separated by groups.
    # Divide each group equally into partitions.
    # min_sample_per_group is the minimum number of samples in each partition.
    n_group <- length(Xs)
    pooled_X <- do.call(rbind, Xs)
    # Generates a pooled group index according to the number of rows in each element of Xs.
    pooled_grouping <- cumsum(as.numeric(unlist(sapply(
        Xs, function(x){
            s <- numeric(nrow(x)) * 0
            s[1] <- 1
            s
            }))))
    # Divide each group into subgroups with at least min_sample_per_group elements.
    pooled_partition <- c()
    for(i_group in 1:n_group){
        current_X <- Xs[[i_group]]
        current_n <- nrow(current_X)
        min_sample_per_group <- min(floor(current_n / 2), min_sample_per_group)
        pooled_partition <- c(
            pooled_partition,
            .generate_subgroups(current_n, min_sample_per_group))
    }
    list(g = pooled_grouping, Pg = pooled_partition, X = pooled_X, M = list())
}

.eachidx <- function(J1, l, x){
    out <- 1:(J1 * l)
    start <- seq(from = 1, to = J1 * l, by = J1)[x]
    end <- seq(from = J1, to = J1 * l, by = J1)[x]
    out[start:end]
}

.ApproximateJointDiagonalizer <- function(M){
    t(jointDiag::uwedge(M)$B)
}

# Generate subgroups.
# Each subgroup has at least min_sample_per_subgroup samples.
# For example, if n = 9 and min_sample_per_subgroup = 2,
# the result will be 1, 1, 2, 2, 3, 3, 4, 4, 4.
.generate_subgroups <- function(n, min_sample_per_subgroup){
    min_sample_per_subgroup <- as.integer(min_sample_per_subgroup)
    if(n < min_sample_per_subgroup){
        stop("n must be greater than min_sample_per_subgroup")
    }
    # If n is not divisible by min_sample_per_subgroup,
    # each subgroup will have at least min_sample_per_subgroup samples.
    # For example, if n = 9 and min_sample_per_subgroup = 2,
    # the result will be 1, 1, 2, 2, 3, 3, 4, 4, 4.
    if(n %% min_sample_per_subgroup == 0){
        rep(1:(n / min_sample_per_subgroup), each = min_sample_per_subgroup)
    }else{
        # At first, calculate the number of samples that can be divided
        # by min_sample_per_subgroup.
        n_subgroup <- floor(n / min_sample_per_subgroup)
        # then, calculate the remainder of the number of samples divided by
        n_remain <- n %% min_sample_per_subgroup
        # Assign the number of samples that can be divided by min_sample_per_subgroup.
        res <- rep(1:n_subgroup, each = min_sample_per_subgroup)
        # Assign the remainder of the number of samples to the last subgroup.
        res <- c(res, rep(n_subgroup, n_remain))
        res
    }
}
