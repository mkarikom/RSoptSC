#' Generate a factored matrix for SymNMF testing
#'
#' @param rows The number of rows in \code{A}
#' @param columns The number of columns in \code{A}
#' @param sparsity binomial p for the nonzero entries in \code{W}/\code{H}
#' @param min The minimum value of \code{W}/\code{H}
#' @param max The maximum value of \code{W}/\code{H}
#' @param k The rank of \code{W}.
#' @param nn A boolean value whether the values are integers, otherwise the round of runif is used
#' @param savecsv A boolean value whether or not to save the matrix as a csv
#' @return \code{W}, \code{H}, and \code{A} such that \code{A} = \code{W} * \code{H}.
#' @examples
#' makeSymmetricM(rows = 25, sparsity = 0.7, min = 1, max = 5, k = 3)
#'
#'  Generates A, W s.t. A = WW'.
#'  Use to generate sample matrices to test factorization.
#'  If \code{savecsv} is set, a csv is saved in inst/testdata
#'  Load this data for testing:
#'  mydata <- read.csv(system.file("testdata","0152_A_3x50.csv",package="RSoptSC"))
#'  devtools::use_data(mydata, internal = TRUE)
makeM <- function(rows = 25, columns = 25,
                           sparsity = 0.7, min = 1, max = 5, k = 3,
                           nn = TRUE, savecsv = FALSE){
  if(nn){
    testW <- matrix(round(runif(n = rows*k, min = min, max = max)), nrow = rows, ncol = k)
    testH <- matrix(round(runif(n = columns*k, min = min, max = max)), nrow = k, ncol = columns)
  }else{
    testW <- matrix(runif(n = rows*k, min = min, max = max), nrow = rows, ncol = k)
    testH <- matrix(runif(n = columns*k, min = min, max = max), nrow = k, ncol = columns)
  }
  testW_drop <- matrix(rbinom(n = rows*k, size = 1, prob = sparsity), nrow = rows, ncol = k)
  testH_drop <- matrix(rbinom(n = columns*k, size = 1, prob = sparsity), nrow = k, ncol = columns)

  W <- testW*testW_drop
  H <- testH*testH_drop
  A <- W%*%H
  output <- list(A = A, W = W, H = H)
  if(savecsv == TRUE){
    newtime <- as.integer(as.POSIXct(Sys.time()))
    newtime <- str_sub(newtime, -4)
    Aname <- paste0(newtime, "_", "A_", k, "x", rows)
    Wname <- paste0(newtime, "_", "W_", k, "x", rows)
    Hname <- paste0(newtime, "_", "H_", k, "x", rows)
    write.csv(output$A, paste0("inst/testdata/", Aname, ".csv"))
    write.csv(output$W, paste0("inst/testdata/", Wname, ".csv"))
    write.csv(output$H, paste0("inst/testdata/", Hname, ".csv"))
  }
  return(output)
}
