context("signaling partners")
library(RSoptSC)


test_that("signaling probability between cells is equal to matlab", {
  # a list of ligands,
  # each of which contains a list of receptor targets,
  # each of which contains both
  # a list of upregulated and a list of downregulated genes
  ligands <- list(Tgfb1 = list(Tgfbr1 = list(up = list('Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'),
                                             down = list('Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1')),
                               Tgfbr2 = list(up = list('Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'),
                                             down = list('Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1'))),
                  Tgfb2 = list(Tgfbr1 = list(up = list('Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'),
                                             down = list('Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1')),
                               Tgfbr2 = list(up = list('Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'),
                                             down = list('Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1'))))
  Pmats <- GetSignalingPartners(JoostSignal$data,
                                  JoostSignal$genes,
                                  ligands)
  expect_equal(Pmats$P, RSoptSC::JoostSignal$P)
  expect_equal(Pmats$P_agg, RSoptSC::JoostSignal$P_agg)
})