context("cluster number inference on Joost data")
library(RSoptSC)

test_that("initial components are correct", {
  comps <- GetComponents(data = joostTest$S, tol = 0.01)
  expect_equal(comps$n_eigs, 19)
})