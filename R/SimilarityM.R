#' Compute the similarity matrix
#'
#' Computes low dim embedding, constructs KNN graph on the embedding -> unweighted adjacency
#' Calls manifold learning algorithm which uses the normalized sample vectors and the
#' unweighted adjacency matrix to compute a low rank approximation of the data
#'
#' @param data the expression data, where each column is treated as a normalized vector
#' @param lambda the balance term between the rank of Z and the error
#' @return a list containing the symetric cell to cell similarity matrix and
#'     manifold learning error

SimilarityM <- function(lambda,data){

  m = nrow(data)
  n = ncol(data)

  X <- data - min(data)
  X = X / max(X)

  for(i in 1:n){
    X[,i] <- X[,i] / norm(X[,i], "2")
  }

  if (m >= 60){

    pca_data = prcomp(t(X))
    coeff1 <- pca_data$rotation[,1:60]
    X1 <- pca_data$x[,1:60]
    pca_eigvalue1 <- pca_data$sdev

  } else {
    pca_data = prcomp(t(X))
    coeff1 <- pca_data$rotation[,1:m]
    X1 <- pca_data$X1[,1:m]
    pca_eigvalue1 <- pca_data$sdev
  }

  eigengaps = abs(pca_eigvalue1[2:(length(pca_eigvalue1)-1)] - pca_eigvalue1[-(1:2)])
  No_Comps1 = which(max(eigengaps)==eigengaps)

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
  dim_init = 3
  # InitY_data = prcomp(t(data))
  # InitY = InitY_data$x[,1:dim_init]
  # Currently this is not "Standardize" as in MATLAB tsne('Standardize', TRUE)
  # Currently this is not "InitialY" as in MATLAB tsne('InitialY', InitY)

  #X2 = tsne::tsne(X = t(X), perplexity = 20, k = dim_init)
  #X2 = as.matrix(X2)
  X2 = Rtsne::Rtsne(X = t(X), perplexity = 20, pca_center = TRUE, pca_scale = TRUE, dims = 3)
  D = matrix(1, n, n)
  if (No_Comps1>=1){
    No_Comps1 = 1
  }
  knn_object = FNN::get.knnx(X2$Y[,1:(No_Comps1+2)], X2$Y[,1:(No_Comps1+2)], k = K)
  IDX = knn_object$nn.index

  for (jj in 1:n){
    D[jj,IDX[jj,]] = 0
  }

  Z <- computM(D,X,lambda)
  Z$Z[Z$Z <= .Machine$double.eps] <- 0
  W <- 0.5*(abs(Z$Z)+abs(t(Z$Z)))
  return(list(W = W, E = Z$E))
}


#' Perform ADMM on the
#'
#' Computing cell-to-cell similarity matrix by solving the following
#' optimization problem via ADMM
#'
#'   min_{Z,E}  ||Z||_* + lambda ||E||_{2,1}
#'   s.t.       X = XZ + E;
#'              Z'1 = 1;
#'               Z_{i,j} = 0 for (i,j)\in Omega
#'
#' @param D the unweighted KNN adjacency matrix
#' @param X normalized sample vectors
#' @param lambda the balance term between the rank of Z and the error
#' @return a list containing the low rank approximation of X and manifold learning error

computM <- function(D,X,lambda){
  ## ADMM iteration
  maxiter = 100
  Err = matrix(0, maxiter, 2)
  rho = 5        # 5
  mu = 10^(-6)   # 10^(-6)
  mumax = 10^(6)
  epsilon = 10^(-5)
  m = nrow(X)
  n = ncol(X)
  Z = matrix(0,n,n)
  E = matrix(0, m, n)
  Y1 = matrix(0, m, n)
  Y2 = matrix(0, 1, n)
  Y3 = matrix(0, n, n)
  iter = 0

  print('Iter  Err')

  while(TRUE){
    iter = iter + 1
    if (iter >= maxiter){
      print("reached max iter")
      return(list(Z = Z, E = norm(X-X%*%Z,"f")))
    }


    # step 1: Update J
    mu = min(rho*mu,mumax)
    svd_data = svd(Z-Y3/mu, nu = nrow(Z-Y3/mu), nv = ncol(Z-Y3/mu))
    U = svd_data$u
    S = diag(svd_data$d, nrow=nrow(Z), ncol=ncol(Z))
    V = svd_data$v

    R = length(diag(S))
    Dmu = matrix(0, R, R)

    MM = pmax(diag(S)-(1/mu)*matrix(1, R, 1), matrix(0, R,1))
    for (i in 1:R){
      Dmu[i,i] = MM[i]
    }

    J = U%*%Dmu%*%t(V)

    if (iter >= 3){
      if (Err[iter-1,1] >= Err[iter-2,1] || norm(as.matrix(X-X%*%Z),"2") <= epsilon){
        print("not reach epsilon")
        return(list(Z = Z, E = norm(X-X%*%Z,"f")))
      }
    }
    # Update E
    QQ = X-X%*%Z+Y1/mu
    for (j in 1:n){
      if (sqrt(sum(QQ[,j]^2)) > lambda/mu){
        E[,j] = (sqrt(sum(QQ[,j]^2))-lambda/mu)/(sqrt(sum(QQ[,j]^2)))*as.matrix(QQ[,j])
      } else {
        E[,j] = matrix(0,length(QQ[,j]))
      }
    }

    eta = norm(X, "2")^2 + norm(matrix(1,n,1),"2")^2+1
    H = -t(X)%*%(X - X%*%Z - E + (1/mu)*Y1) - matrix(1,n,1)%*%(matrix(1,1,n)- matrix(1,1,n)%*%Z + (1/mu)*Y2) + (Z - J + (1/mu)*Y3)
    Z = Z - (1/eta)*H
    Z[D>0] = 0
    # Update Dual variable
    Y1 = Y1 + mu*(X-X%*%Z-E)
    Y2 = Y2 + mu*(matrix(1,1,n) - matrix(1,1,n)%*%Z)
    Y3 = Y3 + mu*(Z - J)

    pracma::fprintf('%d, %8.6f\n',iter,norm(X-X%*%Z-E,"2"))
    Err[iter,] = c(norm(X-X%*%Z-E,"2"), norm(Z-J,"2"))

    if (max(Err[iter,]) <= epsilon){
      print("reached epsilon")
      return(list(Z = Z, E = norm(X-X%*%Z,"f")))
    }
 }
}
