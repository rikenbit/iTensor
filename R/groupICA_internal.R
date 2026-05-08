##############################################
## Internal functions from groupICA package
## (Pfister, N. and Weichwald, S.)
## Original package: https://github.com/sweichwald/groupICA-R
## License: AGPL-3
## These functions are included here because the groupICA package
## is no longer maintained on CRAN.
##############################################

#' @importFrom MASS ginv
#' @importFrom stats cor

.groupICA_Pfister <- function(X,
                     group_index=NA,
                     partition_index=NA,
                     n_components=NA,
                     n_components_uwedge=NA,
                     rank_components=FALSE,
                     pairing='complement',
                     groupsize=1,
                     partitionsize=NA,
                     max_iter=1000,
                     tol=1e-12,
                     silent=TRUE){

  d <- dim(X)[2]
  n <- dim(X)[1]

  # generate group and partition indices as needed
  if(!is.numeric(group_index) & is.na(groupsize)){
    group_index <- rep(0, n)
  }
  else if(!is.numeric(group_index)){
    group_index <- .rigidgroup(n, groupsize)
  }
  if(!is.numeric(partition_index) & is.na(partitionsize)){
    smallest_group <- min(unique(group_index, return_counts=TRUE)$counts)
    partition_index <- .rigidpartition(group_index,
                                      max(c(d, floor(smallest_group/2))))
  }
  else if(!is.numeric(partition_index)){
    partition_index <- .rigidpartition(group_index, partitionsize)
  }


  no_groups <- length(unique(group_index))

  if(!silent){
    print('groupICA: computing covmats')
  }

  # estimate covariances (depending on pairing)
  if(pairing == 'complement'){
    no_pairs <- 0
    for(env in unique(group_index)){
      if(length(unique(partition_index[group_index == env]))>1){
        no_pairs <- no_pairs+length(unique(partition_index[group_index == env]))
      }
    }
    covmats <- vector("list", no_pairs)
    idx <- 1
    for(env in unique(group_index)){
      if(length(unique(partition_index[group_index == env]))==1){
        warning(paste("Removing group", toString(env),
                      "since the partition is trivial, i.e., contains only one set"))
      }
      else{
        for(subenv in unique(partition_index[group_index == env])){
          ind1 <- ((partition_index == subenv) &
                     (group_index == env))
          ind2 <- ((partition_index != subenv) &
                     (group_index == env))
          covmats[[idx]] <- cov(X[ind1,]) - cov(X[ind2,])
          idx <- idx + 1
        }
      }
    }
  }
  else if(pairing == 'allpairs'){
    no_pairs <-  0
    subvec <- rep(0, no_groups)
    for(i in 1:no_groups){
      env <- unique(group_index)[i]
      subvec[i] <- length(unique(partition_index[group_index == env]))
      no_pairs <- no_pairs + subvec[i]*(subvec[i]-1)/2
    }
    covmats <- vector("list", no_pairs)
    idx <- 1
    for(count in 1:no_groups){
      env <- unique(group_index)[count]
      unique_subs <- unique(partition_index[group_index == env])
      if(subvec[count] == 1){
        warning(paste("Removing group", toString(env),
                      "since the partition is trivial, i.e., contains only one set"))
      }
      else{
        for(i in 1:(subvec[count]-1)){
          for(j in (i+1):subvec[count]){
            ind1 <- ((partition_index == unique_subs[i]) &
                       (group_index == env))
            ind2 <- ((partition_index == unique_subs[j]) &
                       (group_index == env))
            covmats[[idx]] <- cov(X[ind1,]) - cov(X[ind2,])
            idx <- idx + 1
          }
        }
      }
    }
  }
  else{
    stop('no appropriate pairing specified')
  }

  # check if there are sufficiently many covariance matrices
  if(length(covmats)<=0){
    stop("Not sufficiently many covariance matrices.")
  }

  # add total observational covariance for normalization
  covmats <- append(list(cov(X)), covmats)

  if(!silent){
    print('groupICA: computed cov matrices')
  }

  # joint diagonalisation
  adj_res <- .uwedge_Pfister(covmats,
                    init=NA,
                    rm_x0=TRUE,
                    return_diag=FALSE,
                    tol=tol,
                    max_iter=max_iter,
                    n_components=n_components_uwedge,
                    silent=silent)
  V <- adj_res$V

  if(!silent){
    print('groupICA: finished uwedge ajd')
  }

  # rank components
  if(rank_components | !is.na(n_components)){
    A <- ginv(V)
    colcorrs <- rep(0, dim(V)[1])
    # running average
    for(k in 1:length(covmats)){
      colcorrs <- colcorrs + diag(abs(cor(A, covmats[[k]] %*% t(V)))) / length(covmats)
    }
    sorting <- order(colcorrs, decreasing=FALSE)
    V <- V[sorting,]
  }
  if(!is.na(n_components)){
    V <- V[1:n_components,]
  }

  # source estimation
  Shat <- t(V %*% t(X))


  return(list(V=V,
              Shat=Shat,
              converged=adj_res$converged,
              n_iter=adj_res$iteration,
              meanoffdiag=adj_res$meanoffdiag))
}


.rigidpartition <- function(group, nosamples){
  partition <- rep(0, length(group))
  for(e in unique(group)){
    partition[group == e] <- .rigidgroup(sum(group == e),
                                               nosamples) + max(partition)
  }
  return(partition)
}


.rigidgroup <- function(len, nosamples){
  groups <- floor(len / nosamples)
  changepoints <- round(seq(1, len, length.out=groups+1))
  index <- rep(0, len)
  for(i in 1:(length(changepoints)-1)){
    index[changepoints[i]:changepoints[i+1]] <- i
  }
  return(index)
}


.uwedge_Pfister <- function(Rx,
                   init=NA,
                   rm_x0=TRUE,
                   return_diag=FALSE,
                   tol=1e-10,
                   max_iter=1000,
                   n_components=NA,
                   silent=TRUE){


  # 0) Preprocessing

  # Remove and remember 0st matrix
  Rx0 <- Rx[[1]]
  if(rm_x0){
    Rx <- Rx[-1]
  }
  d <- dim(Rx0)[1]
  M <- length(Rx)

  if(is.na(n_components)){
    n_components <- d
  }

  # Initial guess
  if(!is.matrix(init) & n_components == d){
    EH <- eigen(Rx[[1]], symmetric=TRUE)
    V <- diag(1/sqrt(abs(EH$values))) %*% t(EH$vectors)
  }
  else if(!is.matrix(init)){
    EH <- eigen(Rx[[1]], symmetric=TRUE)
    mat <- matrix(0, n_components, d)
    diag(mat) <- 1/sqrt(abs(EH$values))[1:n_components]
    V <- mat %*% t(EH$vectors)
  }
  else{
    V <- init[1:n_components,]
  }

  V <- V / matrix(sqrt(rowSums(V^2)), n_components, d)

  converged <- FALSE
  for(iteration in 1:max_iter){

    # 1) Generate Rs
    Rs <- lapply(Rx, function(x) V %*% x %*% t(V))

    # 2) Use Rs to construct A, equation (24) in paper with W=Id
    # 3) Set A1=Id and substitute off-diagonals
    Rsdiag <- sapply(Rs, diag)
    Rsdiagprod <- Rsdiag %*% t(Rsdiag)
    denom_mat <- outer(diag(Rsdiagprod), diag(Rsdiagprod)) - Rsdiagprod^2
    Rkl_list = lapply(Rs, function(x) matrix(diag(x), n_components, n_components, byrow=T)*x)
    Rkl <- Reduce("+", Rkl_list)/M
    num_mat <- matrix(diag(Rsdiagprod), n_components, n_components)*Rkl - Rsdiagprod*t(Rkl)
    denom_mat[abs(denom_mat) < .Machine$double.tol] <- .Machine$double.tol
    diag(denom_mat) <- 1
    A <- num_mat / denom_mat
    diag(A) <- 1

    # 4) Set new V
    Vold <- V
    V <- qr.solve(A, Vold)

    # 5) Normalise V
    V <- V / matrix(sqrt(rowSums(V^2)), n_components, d)

    # 6) Check convergence
    changeinV = max(abs(V-Vold))
    if(changeinV < tol){
      converged = TRUE
      break
    }
  }

  # Rescale
  normaliser <- diag(V %*% Rx0 %*% t(V))
  V <- V / matrix((sign(normaliser)*sqrt(abs(normaliser))), n_components, d)
  Rxdiag <- lapply(Rx, function(x) V %*% x %*% t(V))
  entries_tot <- M*(n_components^2-n_components)
  meanoffdiag <- sqrt(sum(sapply(Rxdiag, function(x) sum(x^2)-sum(diag(x^2))))/entries_tot)

  # Return
  if(return_diag){
    res <- list(V=V,
                Rxdiag=Rxdiag,
                converged=converged,
                iteration=iteration,
                meanoffdiag=meanoffdiag)
  }
  else{
    res <- list(V=V,
                Rxdiag=NA,
                converged=converged,
                iteration=iteration,
                meanoffdiag=meanoffdiag)
  }

  return(res)
}
