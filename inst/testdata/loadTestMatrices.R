### Load the test matrices into sysdata

A3x50 <- as.matrix(read.csv(system.file("testdata", "0152_A_3x50.csv",
                                        package = "RSoptSC"), header = TRUE)[,-1])
H3x50 <- as.matrix(read.csv(system.file("testdata", "0152_H_3x50.csv",
                                        package = "RSoptSC"), header = TRUE)[,-1])
W3x50 <- as.matrix(read.csv(system.file("testdata", "0152_W_3x50.csv",
                                        package = "RSoptSC"), header = TRUE)[,-1])
A0146 <- as.matrix(read.csv(system.file("testdata", "Sym_0146_A_3x25.csv",
                                        package = "RSoptSC"), header = TRUE)[,-1])
W0146 <- as.matrix(read.csv(system.file("testdata", "Sym_0146_W_3x25.csv",
                                        package = "RSoptSC"), header = TRUE)[,-1])
Winit0146 <- as.matrix(read.csv(system.file("testdata", "Sym_0146_W_init_3x25.csv",
                                            package = "RSoptSC"), header = FALSE))
Wnmf0146 <- as.matrix(read.csv(system.file("testdata", "Sym_0146_Wnmf_3x25.csv",
                                           package = "RSoptSC"), header = FALSE))
Dist0146 <- as.matrix(read.csv(system.file("testdata", "Sym_0146_distance_25x25.csv",
                                           package = "RSoptSC"), header = FALSE))
Latent0146 <- as.matrix(read.csv(system.file("testdata", "Sym_0146_latent_2x25.csv",
                                             package = "RSoptSC"), header = FALSE))
Labels0146 <- as.matrix(read.csv(system.file("testdata", "Sym_0146_labels_1x25.csv",
                                             package = "RSoptSC"), header = FALSE))
load(system.file("testdata", "GuoPtime.rda", package = "RSoptSC"))

devtools::use_data(A3x50, H3x50, W3x50, A0146, W0146, Winit0146,
                   Wnmf0146, Dist0146, Latent0146, Labels0146, GuoPtime, internal = TRUE, overwrite = TRUE)
