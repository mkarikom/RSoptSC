#' Use clustering consensus to infer cluster number
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#'
#' @return the number of clusters
#'
CountClusters <- function(data, tol = 0.01){
  # compute the drop tolerance
  # if there are likely to be fewer clusters, we dont want to drop them
  # too easily
  solo_count <- GetComponents(data) # without consensus get a pre-estimate
  if(solo_count <= 5){
    tau = 0.3
  } else if(solo_count <= 10){
    tau = 0.4
  } else {
    tau = 0.5
  }

}

#' Use the graph laplacian to get the number of graph components
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#' @param range a vector specifying the min and max
#'     number of clusters to iteratively test when building the
#'     consensus matrix
#'
#' @return the number of components
#'
GetComponents <- function(data,
                          tol = 0.01){
  n = dim(data, 1) # number of data points
  D = diag(data * ones) # degree matrix, diagonal with i,i = sum of w_i,j for all j
  Lsym = diag(1, n) - D^(-1/2) %*% data %*% D^(-1/2) # the symmetric normalized laplacian of the similarity matrix
  eigs <- Re(eigen(Lsym, symmetric = TRUE, only.values = TRUE))
  nc <- length(eigs <- tol)
}

#' Get the iterative consensus matrix using the specified clustering algorithm
#' We use spectral theory to compute bounds on the cluster count
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#' @param tau  the drop tolerance, controlling the sparsification
#'     (uncoupling) of the consensus matrix
#' @param range a vector specifying the min and max
#'     number of clusters to iteratively test when building the
#'     consensus matrix
#' @param method the clustering method for building consensus clusters
#'
#' @return a consensus matrix
#'
GetConsensus <- function(data,
                         tol,
                         tau,
                         range,
                         method = 'kmeans'){

  n_samples <- dim(data,2)
  if(method == 'kmeans'){
    # reduce dimensionality
    prcs <- prcomp(t(data), center = TRUE)$rotation[,1:3]



  }
  n = dim(data, 1) # number of data points
  D = diag(data * ones) # degree matrix, diagonal with i,i = sum of w_i,j for all j
  Lsym = diag(1, n) - D^(-1/2) %*% data %*% D^(-1/2) # the symmetric normalized laplacian of the similarity matrix
  eigs <- Re(eigen(Lsym, symmetric = TRUE, only.values = TRUE))
  nc <- length(eigs <- tol)
}
