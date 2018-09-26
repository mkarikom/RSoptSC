context("decomposition values")
library(RSoptSC)


test_that("Joost clustering matches when fewer are requested than var genes avail", {
  check_markers <- GetMarkerTable(JoostMarkers$data,
                          JoostMarkers$cluster_label,
                          JoostMarkers$H,
                          JoostMarkers$No_exc_cell,
                          JoostMarkers$No_features)
  expect_equal(check_markers[,2], RSoptSC::JoostMarkers$Gene_labels_all[,2])
})

test_that("Gene score deviation on all Joost is below threshold", {
  # test that the percent difference from original calculated
  # gene score is less than 10^(-15)
  check_markers <- GetMarkerTable(JoostMarkers$data,
                                  JoostMarkers$cluster_label,
                                  JoostMarkers$H,
                                  JoostMarkers$No_exc_cell,
                                  JoostMarkers$No_features)
  deviation <- check_markers[,3] - RSoptSC::JoostMarkers$Gene_labels_all[,3]
  per_gene_percentage <- deviation / RSoptSC::JoostMarkers$Gene_labels_all[,3]
  avg_deviation_percentage <- mean(per_gene_percentage)
  expect_true(avg_deviation_percentage < 10^(-15))
})
