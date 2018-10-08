# These functions are used mainly to generate matrices whose factorization
# is known.  Random matrices are produced according to some parameters, so
# while it may not be ideal to call these functions during testing, they
# can be used to generate matrices which are later placed in:
# 1) data/<object name>.rda (usually reserved for data supporting package functions)
# 2) inst/testdata/<object name>.rda (a popular choice)
# 3) R/sysdata.rda (where the bulk of test data is currently stored and called)

#' Read in a CSV file
#'
#' Call systemfile and load some data
#' @param filename the name of the csv file
#' @param foldername the parent directory of \code{filename}
#' @param thispackage this data will always be in RSoptSC/extdata
#' @return a matrix containing the data
#'
csv2matrix <- function(filename,
                       foldername = "testdata", thispackage = "RSoptSC" ){
  pathstring <- system.file(foldername,filename,thispackage)
  dframe <- read.csv(pathstring)
  dmatrix <- as.matrix(dframe[-1,])
  return(dmatrix)
}


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
#' @param labeltime A boolean value whether or not to label the output file with the last 4 of the unix time
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
                  nn = TRUE, savecsv = FALSE, labeltime = FALSE){
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
    if(labeltime){
      newtime <- as.integer(as.POSIXct(Sys.time()))
      newtime <- str_sub(newtime, -4)
    }else{
      newtime = c()
    }
    Aname <- paste0(newtime, "_", "A_", k, "x", rows)
    Wname <- paste0(newtime, "_", "W_", k, "x", rows)
    Hname <- paste0(newtime, "_", "H_", k, "x", rows)
    write.csv(output$A, paste0("inst/testdata/", Aname, ".csv"))
    write.csv(output$W, paste0("inst/testdata/", Wname, ".csv"))
    write.csv(output$H, paste0("inst/testdata/", Hname, ".csv"))
  }
  return(output)
}

#' Generate symmetric matrix for SymNMF testing
#'
#' @param rows The number of rows/columns in \code{A}
#' @param sparsity binomial p for the nonzero entries in \code{W}
#' @param min The minimum value of \code{W}
#' @param max The maximum value of \code{W}
#' @param k The rank of W.
#' @param nn A boolean value whether the values are integers, otherwise the round of runif is used
#' @param savecsv A boolean value whether or not to save the matrix as a csv
#' @param labeltime A boolean value whether or not to label the output file with the last 4 of the unix time
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
#'
#'  Load this data for testing:
#'  mydata <- read.csv(system.file("testdata","0152_A_3x50.csv",package="RSoptSC"))
#'  devtools::use_data(mydata, internal = TRUE)
makeSymmetricM <- function(rows = 25, sparsity = 0.7, min = 1,
                           max = 5, k = 3, nn = TRUE, savecsv = FALSE,
                           labeltime = FALSE){
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
    if(labeltime){
      newtime <- as.integer(as.POSIXct(Sys.time()))
      newtime <- str_sub(newtime, -4)
    }else{
      newtime = c()
    }
    Aname <- paste0("Sym_", newtime, "_", "A_", k, "x", rows)
    Wname <- paste0("Sym_", newtime, "_", "W_", k, "x", rows)
    write.csv(output$A, paste0("inst/testdata/", Aname, ".csv"))
    write.csv(output$W, paste0("inst/testdata/", Wname, ".csv"))
  }
  return(output)
}
