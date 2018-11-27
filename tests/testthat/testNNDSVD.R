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

test_that("1st order NMF grad < tol * initGrad", {
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

test_that("Lee NMF generates labels whose ARI/NMI is within 2% of matlab", {
  res <- NMF::nmf(x = RSoptSC:::JoostCluster$W,
                  rank = n_clusters,
                  method = 'lee',
                  seed = 'nndsvd',
                  .options = 'nP');
  H <- NMF::basis(res)
  labels <- apply(H, 1, function(x){
    which(x == max(x))})
  true_labels <- as.numeric(as.factor(RSoptSC::GSE67602_Joost$annotation))
  matlab_labels <- as.numeric(RSoptSC:::JoostMarkers$cluster_label)

  ARIresClusterMatlab <- ClusterR::external_validation(true_labels = true_labels,
                                                       clusters = matlab_labels,
                                                       method = 'adjusted_rand_index')
  ARIresCluster <- ClusterR::external_validation(true_labels = true_labels,
                                                 clusters = as.numeric(labels),
                                                 method = 'adjusted_rand_index')
  NMIresClusterMatlab <- ClusterR::external_validation(true_labels = true_labels,
                                                       clusters = matlab_labels,
                                                       method = 'nmi')
  NMIresCluster <- ClusterR::external_validation(true_labels = true_labels,
                                                 clusters = as.numeric(labels),
                                                 method = 'nmi')
  ARItol <- ARIresClusterMatlab * 0.02
  NMItol <- NMIresClusterMatlab * 0.02

  expect_true(abs(ARIresCluster - ARIresClusterMatlab) < ARItol)
  expect_true(abs(NMIresCluster - NMIresClusterMatlab) < NMItol)
})
