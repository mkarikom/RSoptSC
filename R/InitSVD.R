#' Initialize non-negative factorization of the similarity matrix
#'
#' @param A The similarity matrix.
#' @param k The rank of the output.
#' @return \code{W} and \code{H} such that \code{A} = \code{W} * \code{H}.

InitSVD <- function(A, k){
  A_svd = RSpectra::svds(A, k)
  W_init = matrix(0, nrow(A), k)
  H_init = matrix(0, k, ncol(A))

  W_init[,1] = sqrt(A_svd$d[1])*abs(A_svd$u[,1])
  H_init[1,] = sqrt(A_svd$d[1])*abs(t(A_svd$v[,1]))

  for(i in 2:k){
    u = as.matrix(A_svd$u[,i])
    v = as.matrix(A_svd$v[,i])
    pu = u * (u>=0)
    nu = -u * (u<0)
    pv = v * (v>=0)
    nv = -v * (v<0)

    punorm = sqrt(sum(pu^2))
    nunorm = sqrt(sum(nu^2))
    pvnorm = sqrt(sum(pv^2))
    nvnorm = sqrt(sum(nv^2))

    pos = punorm * pvnorm
    neg = nunorm * nvnorm

    if(pos >= neg){
      unorm = pu / punorm
      vnorm = t(pv) / pvnorm
      sigma = pos
    }else{
      unorm = nu / nunorm
      vnorm = t(nv) / nvnorm
      sigma = neg
    }

    W_init[,i] = sqrt(A_svd$d[i]*sigma)*unorm
    H_init[i,] = sqrt(A_svd$d[i]*sigma)*vnorm
  }

  W_init[W_init < 0.00000000001] = mean(A)
  H_init[H_init < 0.00000000001] = mean(A)


  return(list(W = W_init,H = H_init))
}
