context("cluster number inference on Joost data")
library(RSoptSC)

test_that("components based on Laplacian L0 work for W", {
  browser()
  # test that # of zero eigs is same
  comps <- GetComponents(data = RSoptSC::JoostCluster$W, tol = 0.01)
  expect_equal(comps$n_eigs, RSoptSC::JoostCluster$No_cluster1)

  # test that the eigs are close enough
  error <- GetError(comps$eigs, RSoptSC::JoostCluster$ZZ)
  expect_true(error < 10^(-1))
})
