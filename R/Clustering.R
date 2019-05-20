#' Cluster the Cells
#' 
#' Use NMF to cluster the cells.  Currently, the Lee and Seung algorithm is implemented via a call to renozao/NMF.  If cluster-number is inferred then output relevent diagnostics (eg the eigenvalue spaceing of the ensemble method). 
#'
#' @param similarityMatrix a symmetric nonnegative similarity matrix
#' @param n_clusters default is NULL, ensemble method will be called
#' @param n_comp the number of similarity components to use for ensemble clustering
#' @param ... additional parameters for NMF
#'
#' @return a list containing 
#'     \item{H}{the cluster weight matrix, the factorization of the similarity matrix}
#'     \item{labels}{the cluster labels}
#'
#' @importFrom NMF basis nmf seed
#'
#' @export
#'
ClusterCells <- function(similarityMatrix = NULL, 
                         n_clusters = NULL,
                         n_comp = 3,
                         ...){
  
  if(is.null(n_clusters)){
    clusters <- CountClusters(data = similarityMatrix, n_comp = n_comp)
    n_clusters <- clusters$upper_bound
  }
  
  output_NMF <- nmf(x = S$W,
                   rank = n_clusters,
                   method = 'lee',
                   seed = 'nndsvd',
                   ...)
  H <- basis(output_NMF)
  
  labels <- apply(H, 1, function(x){
    which(x == max(x))})
  
  return(list(H = H, 
              labels = labels))
}

#' Use clustering consensus to infer cluster number
#' 
#' Use clustering consensus to infer cluster number.  Plot the the eigenvalue spacing.
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#' @param n_comp  the number of similarity components to use for ensemble clustering
#' @param range a vector specifying the min and max
#'     number of clusters to iteratively test when building the
#'     consensus matrix
#' @param eigengap whether or not to use the max
#'     eigengap (upper bound) cluster count
#'
#' @return the number of clusters and the vector of eigenvalues
#'
#' @export
#'
CountClusters <- function(data, tol = 0.01, range = 2:20, eigengap = TRUE, n_comp){
  # compute the drop tolerance, enforcing parsimony of components
  solo_count <- GetComponents(data, tol = tol)$n_eigs # without consensus get a pre-estimate
  if(solo_count <= 5){
    tau = 0.3
  } else if(solo_count <= 10){
    tau = 0.4
  } else {
    tau = 0.5
  }
  cmatrix <- GetEnsemble(data, tol, tau, n_comp = n_comp, range = range)
  eigs <- GetComponents(cmatrix, tol = tol)
  # compute the largest eigengap
  gaps <- eigs$val[2:length(eigs$val)] - eigs$val[1:(length(eigs$val)-1)]
  upper_bound <- which(gaps == max(gaps))
  
  # compute the number of zero eigenvalues
  lower_bound <- length(eigs$val[which(eigs$val <= tol)])
  
  plot(c(1:20), 
       eigs$val[1:20],
       xlab = NA,
       ylab = 'eigenvalues',
       main = 'Eigenvalues of the Graph Laplacian')
  
  return(list(upper_bound = upper_bound, 
              lower_bound = lower_bound,
              eigs = eigs))
}

#' Get graph components
#' 
#' Use the graph laplacian to get the number of graph components.
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#'
#' @importFrom expm sqrtm
#'
#' @return the number of components
#'
GetComponents <- function(data,
                          tol = 0.01){
  n = nrow(data) # number of data points
  D = diag(as.vector(data %*% matrix(1, n, 1))) # degree matrix, diagonal with i,i = sum of w_i,j for all j
  DD <- solve(sqrtm(D))
  Lsym = diag(1, n) - DD %*% data %*% DD # the symmetric normalized laplacian of the similarity matrix
  eigs <- abs(Re(eigen(Lsym, only.values = TRUE)$values))
  n_zeros <- length(eigs[which(eigs <= tol)])
  return(list(val = sort(eigs), n_eigs = n_zeros))
}

#' Produce consensus matrix
#' 
#' Produce a truncated ensemble consensus matrix.
#'
#' @param data a symmetric nonnegative similarity matrix
#' @param tol  cutoff for lambda zero
#' @param n_comp  the number of similarity components to use for ensemble clustering
#' @param tau  the drop tolerance, controlling the sparsification
#'     (uncoupling) of the consensus matrix
#' @param range a vector specifying the min and max
#'     number of clusters to iteratively test when building the
#'     consensus matrix
#'
#' @importFrom cluster pam
#' @importFrom stats prcomp
#'
#' @return a truncated ensemble consensus matrix
#'
GetEnsemble <- function(data,
                        tol,
                        n_comp = 3,
                        tau,
                        range = 2:20){
  
  n_samples <- dim(data)[2]
  # reduce dimensionality
  prcs <- prcomp(t(data), center = TRUE)
  
  # perform kmeans clustering over the range specified
  cluster_assign <- sapply(range, function(x){
    pam(prcs$x[,1:n_comp], k = x )$clustering
  })
  
  
  # generate consensus matrices
  
  consen <- matrix(0, n_samples, n_samples)
  for(i in 1:ncol(cluster_assign)){
    cons <- GetConsensus(cluster_assign[,i])
    consen <- consen + cons
  }
  
  # truncate the ensemble consensus matrix
  consen[which(consen <= length(range) * tau)] <- 0
  
  # normalize and make symmetric
  consen <- (consen + t(consen))/2
}

#' Produce a consensus matrix
#' 
#' Produce a consensus matrix for ensemble analysis.
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
