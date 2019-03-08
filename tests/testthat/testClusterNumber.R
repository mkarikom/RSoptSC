context("cluster number inference on Joost data")
library(RSoptSC)

test_that("components based on Laplacian L0 work for W", {
  # test that # of zero eigs is same
  comps <- GetComponents(data = JoostCluster$W, tol = 0.01)
  expect_equal(comps$n_eigs, JoostCluster$No_cluster1)

  # test that the eigs are close enough
  deviation <- abs(comps$val - JoostCluster$ZZ)
  percentage <- deviation / JoostCluster$ZZ
  avg_deviation_percentage <- mean(percentage)
                                   
  expect_true(avg_deviation_percentage < 10^(-1))
})
