#' Generate symmetric matrix for SymNMF testing
#'
#' @param rows The number of rows/columns in \code{A}
#' @param sparsity binomial p for the nonzero entries in \code{W}
#' @param min The minimum value of \code{W}
#' @param max The maximum value of \code{W}
#' @param k The rank of W.
#' @param nn A boolean value whether the values are integers, otherwise the round of runif is used
#' @param savecsv A boolean value whether or not to save the matrix as a csv
#' @return \code{W} and \code{A} such that \code{A} = \code{W} * \code{t(W)}.
#' @examples
#' makeSymmetricM(rows = 25, sparsity = 0.7, min = 1, max = 5, k = 3)
#'
#'  Generates A, W s.t. A = WW'.
#'  Use to generate sample matrices to test factorization.
#'  If \code{savecsv} is set, a csv is saved in inst/testdata
#'
#'  newtime <- as.integer(as.POSIXct(Sys.time()))
#'  test <- makeSymmetricM()
#'  newname <- paste0("mat3x25x", newtime)
#'  assign(newname, test)
#'  write.csv(newname, paste0("inst/testdata/", newname, ".csv"))
#'  devtools::use_data(newname, internal = TRUE)
makeSymmetricM <- function(rows = 25, sparsity = 0.7, min = 1,
                           max = 5, k = 3, nn = TRUE, savecsv = FALSE){
  if(nn){
    testW <- matrix(round(runif(n = rows*k, min = min, max = max)), nrow = rows, ncol = k)
  }else{
    testW <- matrix(runif(n = rows*k, min = min, max = max), nrow = rows, ncol = k)
  }
  testW_drop <- matrix(rbinom(n = rows*k, size = 1, prob = sparsity), nrow = rows, ncol = k)
  W <- testW*testW_drop
  A <- W%*%t(W)
  output <- list(A = A, W = W)
  if(savecsv == TRUE){
    newtime <- as.integer(as.POSIXct(Sys.time()))
    newtime <- str_sub(newtime, -4)
    Aname <- paste0("Sym_", newtime, "_", "A_", k, "x", rows)
    Wname <- paste0("Sym_", newtime, "_", "W_", k, "x", rows)
    write.csv(output$A, paste0("inst/testdata/", Aname, ".csv"))
    write.csv(output$W, paste0("inst/testdata/", Wname, ".csv"))
  }
  return(output)
}
