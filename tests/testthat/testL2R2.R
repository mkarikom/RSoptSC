context("similarity matrix based on L2R2 manifold embedding")
library(RSoptSC)

test_that("norm(X-X%*%Z,\"f\") is less than or equal to matlab", {
  
  W <- SimilarityM(lambda = JoostL2R2$lambda, 
                   data = JoostL2R2$data,
                   perplexity = 25, 
                   dims = 3, 
                   pca = TRUE,
                   theta = 0)
  expect_true(W$E <= JoostL2R2$FNorm_XXZ_error)
})
