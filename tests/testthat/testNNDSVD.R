context("decomposition values")
library(RSoptSC)

# test_that("NMF precision 3", {
#   data = RSoptSC::SymNMF(RSoptSC::symA3x25, nC = 3, H = RSoptSC::symWinit3x25)
#   product = prod(signif(data$H, digits = 2) - signif(RSoptSC::symW3x25, digits = 2) == 0)
#   print(product)
#   print(data)
#   print(RSoptSC::symW3x25)
#   expect_equal(product, 1)
# })

test_that("NMF grad < tol * initGrad", {
  skip("takes too long")
  tol = 10^(-6)
  data = RSoptSC::SymNMF(RSoptSC::A0146, nC = 3, H = RSoptSC::Winit0146)
  A = RSoptSC::A0146
  Hinit = RSoptSC::Winit0146
  initNorm = norm(4*(Hinit%*%t(Hinit)-A)%*%Hinit, "F")
  H = data$H
  finalNorm = norm(4*(H%*%t(H)-A)%*%H, "F")
  expect_true(finalNorm < tol * initNorm)
})
