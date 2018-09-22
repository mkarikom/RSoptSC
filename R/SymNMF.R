#' Initialize non-negative factorization of the similarity matrix
#'
#' @param A The similarity matrix.
#' @param k The rank of the output.
#' @return \code{W} and \code{H} such that \code{A} = \code{W} * \code{D} * \code{H}.


SymNMF <- function(A, nC, H, gamma = 0.000001, mu = 10^(-6), maxiter = 1000000){
  n = nrow(H)
  k = ncol(H)
  Hnew = vec(matrix(.Machine$integer.max, n, k))
  H = vec(H)
  Hinit = H
  Hold = H
  H = H = nnproj(H - gamma*vgrad(A, H, n))
  iter = 1
  olditer = 1
  initGradNorm = norm(matrix(vgrad(A, Hinit, n),nrow=n),"F")

  while(iter < maxiter && norm(matrix(vgrad(A, H, n),nrow=n),"F") > mu*norm(matrix(vgrad(A, Hinit, n), nrow=n),"F")){
    iter = iter + 1
    Hold = H
    H = nnproj(H - gamma*vgrad(A, H, n))
    if(objective(A, H, n) > objective(A, Hold, n)){
      return(vmat(Hold, n))
    }
    # if(objective(A, H, n) < 0.0038){
    #   return(vmat(Hold, n))
    # }
    if(iter == olditer + 1000){
      print(objective(A, H, n))
      print(iter)
      olditer = iter
    }
  }
  return(list(H = vmat(H, n), obj = objective(A, H, n), iter = iter, initGradNorm = initGradNorm, GradNorm = norm(matrix(vgrad(A, H, n),nrow=n),"F")))
}

objective <- function(A, H, n){
  H = vmat(H, n)
  obj = norm(A, 'f')^2 - 2 * psych::tr(t(H) %*% (A %*% H)) + psych::tr(t(H) %*% H %*% t(H) %*% H )
  H = vec(H)
  return(obj)
}

vgrad <- function(A, H, n){
  H = vmat(H, n)
  mGrad = 4*(H%*%t(H)-A)%*%H
  vGrad = vec(mGrad)
  H = vec(H)
  return(vGrad)
}

vec <- function(M){
  return(sapply(M, function(x) t(x)))
}

vmat <- function(M, rows){
  return(matrix(M, nrow = rows))
}

vnorm <- function(v){
  return(sqrt(sum(v^2)))
}

nnproj <- function(M){
  mat = M * (M>=0) + 0 * nproj(M)
  return(mat)
}

nproj <- function(M){
  return(M * (M<0))
}
