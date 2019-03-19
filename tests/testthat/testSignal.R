context("signaling partners")
library(RSoptSC)
library(Matrix)

data <-as.matrix(joostTest$data)
gene_names <- joostTest$gene_names

test_that("signaling is accurate when only up targs are used", {
  lig_rec_path <- system.file("extdata", "tgfb_lig_rec.tsv", package = "RSoptSC")
  rec_target_path <- system.file("extdata", "tgfb_rec_target_up.tsv", package = "RSoptSC")
  pathway <- ImportPathway(lig_table_path = lig_rec_path,
                           rec_table_path = rec_target_path,
                           data = data,
                           gene_names = gene_names)
  Pmats <- GetSignalingPartners(data,
                                gene_names,
                                pathway$pathway_removed)
  expect_true(sum(abs(Pmats$P_agg - as.matrix(joostTest$signal_all_up))) <= 1e-15)
})

test_that("signaling is accurate when up and down targs are used", {
  lig_rec_path <- system.file("extdata", "tgfb_lig_rec.tsv", package = "RSoptSC")
  rec_target_path <- system.file("extdata", "tgfb_rec_target_both.tsv", package = "RSoptSC")
  pathway <- ImportPathway(lig_table_path = lig_rec_path,
                            rec_table_path = rec_target_path,
                            data = data,
                            gene_names = gene_names)
  Pmats <- GetSignalingPartners(data,
                                 gene_names,
                                 pathway$pathway_removed)
  print(dim(sparseJoostTest$signal_all_both))
  expect_true(sum(abs(Pmats$P_agg - as.matrix(joostTest$signal_all_both))) <= 1e-15)
})

test_that("signaling is accurate when no targets are used", {
  lig_rec_path <- system.file("extdata", "tgfb_lig_rec.tsv", package = "RSoptSC")
  pathway <- ImportPathway(lig_table_path = lig_rec_path,
                           data = data,
                           gene_names = gene_names)
  Pmats <- GetSignalingPartners(data,
                                gene_names,
                                pathway$pathway_removed)
  expect_true(sum(abs(Pmats$P_agg - as.matrix(joostTest$signal_all_notarg))) <= 1e-15)
})

