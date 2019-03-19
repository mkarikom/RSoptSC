context("ADMM result")
library(RSoptSC)

test_that("norm(X-X%*%Z,\"f\") is correct", {
  skip_on_travis()
  data <- Matrix::as.matrix(joostTest$data)
  X <- data - min(data)
  X <- X / max(X)
  for(i in 1:ncol(data)){
    X[,i] <- X[,i] / norm(X[,i], "2")
  }
  # computes the similarity matrix using a low rank representation of X subject to the adjacency constraint D
  M <- computeM(lambda = 0.5, X = X, D = Matrix::as.matrix(joostTest$D))
  expect_true(abs(M$E-16.7979) <= 1e-5)
})
