#' Use clustering consensus to infer cluster number
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#' @param range a vector specifying the min and max
#'     number of clusters to iteratively test when building the
#'     consensus matrix
#' @param eigengap whether or not to use the max
#'     eigengap (upper bound) cluster count
#'
#' @return the number of clusters
#'
CountClusters <- function(data, tol = 0.01, range = 1:20, eigengap = TRUE){
  # compute the drop tolerance, enforcing parsimony of components
  solo_count <- GetComponents(data, tol = 0.01)$n_eigs # without consensus get a pre-estimate
  if(solo_count <= 5){
    tau = 0.3
  } else if(solo_count <= 10){
    tau = 0.4
  } else {
    tau = 0.5
  }
  cmatrix <- GetEnsemble(data, tol, tau, n_prcs = 3, range = range, method = 'kmeans')
  eigs <- GetComponents(cmatrix, tol = 0.01)

  # compute the largest eigengap
  gaps <- eigs$val[2:length(eigs$val)] - eigs$val[1:(length(eigs$val)-1)]
  upper_bound <- which(gaps == max(gaps))

  # compute the number of zero eigenvalues
  lower_bound <- length(eigs$val[which(eigs$val <= tol)])

  if(eigengap){
    return(upper_bound)
  } else {
    return(lwer_bound)
  }
}

#' Use the graph laplacian to get the number of graph components
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#'
#' @return the number of components
#'
GetComponents <- function(data,
                          tol = 0.01){
  n = nrow(data) # number of data points
  D = diag(as.vector(data %*% matrix(1, n, 1))) # degree matrix, diagonal with i,i = sum of w_i,j for all j
  Lsym = diag(1, n) - diag(diag(D)^(-1/2)) %*% data %*% diag(diag(D)^(-1/2)) # the symmetric normalized laplacian of the similarity matrix
  eigs <- abs(Re(eigen(Lsym, only.values = TRUE)$values))
  n_zeros <- length(eigs[which(eigs <= tol)])
  return(list(val = sort(eigs), n_eigs = n_zeros))
}

#' Produce a truncated ensemble consensus matrix
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#' @param prcs_dim  the number of pcs to use for clustering method
#' @param tau  the drop tolerance, controlling the sparsification
#'     (uncoupling) of the consensus matrix
#' @param range a vector specifying the min and max
#'     number of clusters to iteratively test when building the
#'     consensus matrix
#' @param method the clustering method for building consensus clusters
#'
#' @return a truncated ensemble consensus matrix
#'
GetEnsemble <- function(data,
                        tol,
                        n_prcs = 3,
                        tau,
                        range = 1:20,
                        method = 'kmeans'){
  no_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(no_cores)

  n_samples <- dim(data)[2]
  if(method == 'kmeans'){
    # reduce dimensionality
    prcs <- prcomp(t(data), center = TRUE)

    # perform kmeans clustering over the range specified
    cluster_assign <- sapply(range, function(x){
      kmeans(prcs$x[,1:n_prcs], centers = x)$cluster
    })


    # generate consensus matrices
    parallel::clusterExport(cl, c("cluster_assign", "GetConsensus"), envir = environment())
    consen_list <- parallel::parLapply(cl,
                             split(cluster_assign, col(cluster_assign)),
                             function(x){
                               GetConsensus(x)
                             })
    parallel::stopCluster(cl)
    # generate the ensemble
    consen <- Reduce("+", consen_list)

    # truncate the ensemble consensus matrix
    consen[which(consen <= length(range) * tau)] <- 0

    # normalize and make symmetric
    consen <- (consen + t(consen))/2
  }
}

#' Produce a consensus matrix
#'
#' @param clusters the cluster assignment of each cell
#'
#' @return a consensus matrix
#'
GetConsensus <- function(clusters){
  n_samples <- length(clusters)
  consen <- matrix(0, n_samples, n_samples)
  for(i in 1:n_samples){
    for(ii in 1:n_samples){
      if(clusters[i] == clusters[ii]){
        consen[i, ii] <- 1
      }
    }
  }
  return(consen)
}
