#' Compute the similarity matrix
#'
#' Computes low dim embedding, constructs KNN graph on the embedding -> unweighted adjacency
#' Calls manifold learning algorithm which uses the normalized sample vectors and the
#' unweighted adjacency matrix to compute a low rank approximation of the data.
#'
#' @param data the expression data, where each column is treated as a normalized vector
#' @param lambda the balance term between the rank of Z and the error, default is 0.5
#' @param pre_embed_method how the initial non-linear embedding is performed, default is 'umap'
#' @param comps_knn number of components to use for knn, overrides eigengap-based inference
#' @param use_umap_indices use the knn indices computed during umap embedding to impose sparsity on L2R2, instead of recomputing based on the layout.
#' @param ... extra arguments passed to umap or Rtsne
#' 
#' @return a list containing 
#'     \item{W}{the similarity matrix}
#'     \item{E}{the error of the ADMM step}
#'     \item{nl_embedding}{the KNN sparsity constraint is based on this embedding}
#'
#' @importFrom Rtsne Rtsne
#' @importFrom Matrix as.matrix norm
#' @importFrom FNN get.knnx
#' @importFrom stats prcomp
#' @import umap
#'
#' @export
#'
SimilarityM <- function(lambda = 0.5, data, comps_knn = NULL, 
                        use_umap_indices = FALSE, pre_embed_method = 'umap', ...){
  set.seed(1)
  
  m = nrow(data)
  n = ncol(data)
  
  X <- data - min(data)
  X <- X / max(X)
  
  for(i in 1:n){
    X[,i] <- X[,i] / norm(X[,i], "2")
  }
  
  if (m >= 60){
    
    pca_data = prcomp(t(X))
    pca_eigvalue1 <- pca_data$sdev^2
    
  } else {
    pca_data = prcomp(t(X))
    pca_eigvalue1 <- pca_data$sdev^2
  }
  
  eigengaps = abs(pca_eigvalue1[2:(length(pca_eigvalue1)-1)] - pca_eigvalue1[-(1:2)])
  No_Comps1 = which(max(eigengaps)==eigengaps)
  if (No_Comps1>=1){
    No_Comps1 = 1
  }
  No_Comps1 = No_Comps1 + 2
  if(!is.null(comps_knn)){
     comps_knn = No_Comps1
  }
  
  cc = cumsum(pca_eigvalue1[-(1)])
  dd = cc[-(1)]/sum(pca_eigvalue1[-(1)])
  
  K1 = length(which(dd<=0.3))
  
  if (K1 <= 10){
    K = 10
  } else if (K1 >=30){
    K = 30
  } else {
    K = K1+1
  }
  if(pre_embed_method == 'tsne'){
    X2 <- Rtsne(X = t(as.matrix(X)), ...)
    D <- matrix(1, n, n)
    knn_object = get.knnx(X2$Y[,1:(No_Comps1)], X2$Y[,1:(No_Comps1)], k = K)
    IDX = knn_object$nn.index
  } else if (pre_embed_method == 'umap'){
    data.umap <- umap(d = t(as.matrix(data)), n_components = comps_knn, ...)
    X2 <- data.umap$layout
    
    D <- matrix(1, n, n)
    knn_object = get.knnx(X2[,1:No_Comps1], X2[,1:No_Comps1], k = K)
    IDX <- knn_object$nn.index
    if(use_umap_indices){
      IDX <- data.umap$knn$indexes
    }
  }
  
  for (jj in 1:n){
    D[jj,IDX[jj,]] = 0
  }
  Z <- computeM(D,X,lambda)
  Z$Z[Z$Z <= .Machine$double.eps] <- 0
  W <- 0.5*(abs(Z$Z)+abs(t(Z$Z)))
  return(list(W = W, E = Z$E, nl_embedding = X2))
}


#' Perform ADMM
#'
#' Compute cell-to-cell similarity matrix by solving the following
#' optimization problem via ADMM
#' \deqn{
#'   min_{Z,E}  ||Z||_* + lambda ||E||_{2,1}\\
#'   s.t.       X = XZ + E\\
#'              Z'1 = 1\\
#'               Z_{i,j} = 0 for (i,j) \in Omega}
#'
#' @param D the unweighted KNN adjacency matrix
#' @param X normalized sample vectors
#' @param lambda the balance term between the rank of Z and the error
#' @return a list containing the low rank approximation of X and manifold learning error

computeM <- function(D,X,lambda){
  ## ADMM iteration
  maxiter = 100
  Err <- matrix(0, maxiter, 2)
  rho <- 5        # 5
  mu <- 10^(-6)   # 10^(-6)
  mumax <- 10^(6)
  epsilon <- 10^(-5)
  m <- nrow(X)
  n <- ncol(X)
  Z <- matrix(0,n,n)
  E <- matrix(0, m, n)
  Y1 <- matrix(0, m, n)
  Y2 <- matrix(0, 1, n)
  Y3 <- matrix(0, n, n)
  iter <- 0
  XXZ <- X-X%*%Z
  
  while(TRUE){
    iter <- iter + 1
    if (iter >= maxiter){
      print("reached max iter")
      return(list(Z = Z, E = norm(XXZ,"f")))
    }
    
    
    # step 1: Update J
    mu <- min(rho*mu,mumax)
    svd_data <- svd(Z-Y3/mu, nu = nrow(Z-Y3/mu), nv = ncol(Z-Y3/mu))
    U <- svd_data$u
    S <- diag(svd_data$d, nrow=nrow(Z), ncol=ncol(Z))
    V <- svd_data$v
    
    R <- length(diag(S))
    Dmu <- matrix(0, R, R)
    
    MM <- pmax(diag(S)-(1/mu)*matrix(1, R, 1), matrix(0, R,1))
    for (i in 1:R){
      Dmu[i,i] <- MM[i]
    }
    
    J <- U%*%Dmu%*%t(V)
    
    if (iter >= 3){
      if (Err[iter-1,1] >= Err[iter-2,1] || norm(as.matrix(XXZ),"2") <= epsilon){
        print("convergence by error min")
        return(list(Z = Z, E = norm(XXZ,"f")))
      }
    }
    # Update E
    QQ <- XXZ+Y1/mu
    E <- apply(QQ, 2, function(x){
      normX <- norm(as.matrix(x), type = '2')
      as.vector(((normX - lambda/mu) / normX)) * as.matrix(x)
    })
    normsQQ <- apply(QQ, 2, function(x){
      norm(x, type = '2')
    })
    for(j in 1:n){
      if(normsQQ[j] <= lambda/mu){
        E[,j] <- matrix(0,length(QQ[,j]))
      }
    }
    
    eta <- norm(X, "2")^2 + norm(matrix(1,n,1),"2")^2+1
    H <- -t(X)%*%(XXZ-E + (1/mu)*Y1) - matrix(1,n,1)%*%(matrix(1,1,n)- matrix(1,1,n)%*%Z + (1/mu)*Y2) + (Z - J + (1/mu)*Y3)
    Z <- Z - (1/eta)*H
    Z[D>0] <- 0
    XXZ <- X-X%*%Z
    
    # Update Dual variable
    Y1 <- Y1 + mu*(XXZ-E)
    Y2 <- Y2 + mu*(matrix(1,1,n) - matrix(1,1,n)%*%Z)
    Y3 <- Y3 + mu*(Z - J)
    
    Err[iter,] <- c(norm(XXZ-E,"2"), norm(Z-J,"2"))
    print(paste0("Err = ", Err[iter,1]))
    if (max(Err[iter,]) <= epsilon){
      print("reached epsilon")
      return(list(Z = Z, E = norm(XXZ,"f")))
    }
  }
}
