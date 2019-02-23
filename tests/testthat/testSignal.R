context("signaling partners")
library(RSoptSC)


test_that("signaling probability between cells is equal to matlab", {
  # a list of ligands,
  # each of which contains a list of receptor targets,
  # each of which contains both
  # a list of upregulated and a list of downregulated genes
  pathway <- structure(list(ligand = c("Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2"), receptor = c("Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2"), direction = c("up", "up", "up", "up", "up", "up", "up", "down", "down", "down", "down", "down", "down", "down", "down", "down", "up", "up", "up", "up", "up", "up", "up", "down", "down", "down", "down", "down", "down", "down", "down", "down", "up", "up", "up", "up", "up", "up", "up", "down", "down", "down", "down", "down", "down", "down", "down", "down", "up", "up", "up", "up", "up", "up", "up", "down", "down", "down", "down", "down", "down", "down", "down", "down"), target = structure(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L), .Label = c("Zeb2", "Smad2", "Wnt4", "Wnt11", "Bmp7", "Sox9", "Notch1", "Crebbp", "Fos", "Id1", "Jun", "Runx1", "Smad1", "Smad5", "Sox4", "Cdh1"), class = "factor")), class = "data.frame", row.names = c(NA, -64L))
    
    
    
  browser()  
  Pmats <- GetSignalingPartners(JoostSignal$data,
                                  JoostSignal$genes,
                                  pathway)
  expect_equal(Pmats$P, JoostSignal$P)
  expect_equal(Pmats$P_agg, JoostSignal$P_agg)
})

test_that("signaling probability between cells is equal to matlab", {
  # a list of ligands,
  # each of which contains a list of receptor targets,
  # each of which contains ONLY
  # a list of upregulated genes
  pathway <- structure(list(ligand = c("Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb1", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2", "Tgfb2"), receptor = c("Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr1", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2", "Tgfbr2"), direction = c("up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up", "up"), target = structure(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 1L, 2L, 3L, 4L, 5L, 6L, 7L), .Label = c("Zeb2", "Smad2", "Wnt4", "Wnt11", "Bmp7", "Sox9", "Notch1"), class = "factor")), class = "data.frame", row.names = c(NA, -28L))
  
  
  Pmats <- GetSignalingPartners(JoostSignal$data,
                                JoostSignal$genes,
                                pathway)
  expect_equal(Pmats$P, JoostSignal$P_ind_up_only)
  expect_equal(Pmats$P_agg, JoostSignal$P_agg_up_only)
})
